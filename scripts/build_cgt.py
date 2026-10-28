#!/usr/bin/env python3
"""Build the clinical ground truth (DESIGN v3 §10) from the curated mapping and the frozen sources.

Inputs (read-only): configs/cgt_curated_mapping.csv, configs/cgt_allele_rules.json, configs/cgt_corrections.csv,
configs/fusion_known_positives.csv, external/clinical/{civic,moalmanac,fda_cdx,ctgov,pubmed}, and molecular data only
(OmicsSomaticMutations, P1 cn/damaging/expression/hotspot/fusion caches, Model.csv, OmicsGlobalSignatures via
biomark.molecular). No gene-effect values, predictions or P1 test metrics are read (static scan in
tests/test_cgt_table.py).

Steps: (1) resolve every source record id (CIViC EID/AID, MOAlmanac assertion_id, FDA row_index|PMA, cached CTG record
with biomarker regex, cached PubMed esummary); rows without a resolvable id, T3 rows without a regex-matching CTG
record, and pmid_only rows outside `clinically-refuted` are dropped to results/clinical/cgt_dropped.csv;
(2) collapse drug duplicates per test; (3) compute feature carriers on dev (D) and test (R) lines; (4) testability flags;
(5) write configs/clinical_ground_truth.csv, results/cgt/summary.json, results/cgt/cgt_feature_carriers.csv.

Usage: python -X utf8 scripts/build_cgt.py
"""
from __future__ import annotations

import json
import os
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

os.environ.setdefault("LOKY_MAX_CPU_COUNT", "1")
warnings.filterwarnings("ignore", message="This pattern is interpreted as a regular expression")
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from biomark import genesets, molecular, p1  # noqa: E402  (molecular data only)

CFG = ROOT / "configs"
EXT = ROOT / "external" / "clinical"
MAPPING = CFG / "cgt_curated_mapping.csv"
RULES = CFG / "cgt_allele_rules.json"
CORRECTIONS = CFG / "cgt_corrections.csv"
KNOWN_POS = CFG / "fusion_known_positives.csv"
OUT_CGT = CFG / "clinical_ground_truth.csv"
OUT_DIR = ROOT / "results" / "cgt"
DROPPED = ROOT / "results" / "clinical" / "cgt_dropped.csv"
MUT_CACHE = ROOT / "cache" / "cgt" / "mutations_cgt.parquet"
MUT_COUNT_CACHE = ROOT / "cache" / "cgt" / "protein_altering_counts.parquet"
MUT_FILE = ROOT / "external" / "depmap_extra" / "OmicsSomaticMutations.csv"
MOA_JSON = EXT / "moalmanac" / "extracted" / "vanallenlab-moalmanac-db-6d3bd09" / "molecular-oncology-almanac.json"

STANDARD_FAMILIES = {"AMP", "DEL", "LOSS", "LOFBI", "HOT", "SIL", "EXPDOSE"}
SCORED_PARTITIONS = {"sensitivity", "bypass-resistance", "clinically-refuted"}
P2_HEADLINE_TIERS = {"T1", "T2", "T3"}  # PI: 'preclinical' (RB1->CDK4/6) is a secondary, never in the P2 headline
P2_KO_CLASSES = {"concordant", "concordant-family", "concordant-KO/combination-required"}
ALIASES = {"ABL1": ["ABL"], "ERBB2": ["HER2", "HER-2"], "ESR1": ["ESR1", "ER ", "ER+", "ER-positive", "estrogen", "Overexpression"],
           "KMT2A": ["KMT2A", "MLL"], "CDKN2A": ["CDKN2A", "CCND"], "BRCA1;BRCA2": ["BRCA"], "MSI": ["MSI", "Microsatellite", "MLH1", "mismatch"],
           "TP53": ["17p", "TP53"], "MTAP": ["MTAP"], "CD274": ["PD-L1"], "FOLR1": ["FOLR1"], "CLDN18": ["Claudin", "CLDN18"]}

INDICATION_MAP = CFG / "indication_map.csv"


def load_indications() -> dict[str, list[tuple[str, str]]]:
    """indication code -> [(match_field, match_value)] from configs/indication_map.csv."""
    im = pd.read_csv(INDICATION_MAP, dtype=str, keep_default_na=False)
    out: dict[str, list[tuple[str, str]]] = {}
    for _, r in im.iterrows():
        out.setdefault(r["indication_code"], []).append((r["match_field"], r["match_value"]))
    return out


IND = load_indications()


def log(*a):
    print("[build_cgt]", *a, flush=True)


# ------------------------------------------------------------------------------------------------- sources
def load_sources() -> dict:
    civ_e = pd.read_csv(EXT / "civic" / "nightly-ClinicalEvidenceSummaries.tsv", sep="\t", low_memory=False,
                        usecols=["evidence_id", "molecular_profile", "therapies", "evidence_type", "evidence_status",
                                 "evidence_direction", "significance"])
    civ_a = pd.read_csv(EXT / "civic" / "nightly-AssertionSummaries.tsv", sep="\t",
                        usecols=["assertion_id", "molecular_profile", "therapies", "assertion_type"])
    moa = {int(x["assertion_id"]): x for x in json.loads(MOA_JSON.read_text(encoding="utf-8"))}
    fda = pd.read_csv(EXT / "fda_cdx" / "fda_cdx_parsed.csv").set_index("row_index")
    ctg_man = json.loads((EXT / "ctgov" / "manifest.json").read_text(encoding="utf-8"))
    pm_man = json.loads((EXT / "pubmed" / "manifest.json").read_text(encoding="utf-8"))
    return {"eid": civ_e.set_index("evidence_id"), "aid": civ_a.set_index("assertion_id"), "moa": moa, "fda": fda,
            "ctg_man": ctg_man, "pm_man": pm_man}


def _moa_text(x: dict) -> str:
    return " ".join(str(x.get(k) or "") for k in ("gene", "gene1", "gene2", "protein_change", "description", "status", "therapy_name"))


def ctg_text(nct: str) -> str:
    js = json.loads((EXT / "ctgov" / f"{nct}.json").read_text(encoding="utf-8"))
    ps = js.get("protocolSection", {})
    parts = [ps.get("eligibilityModule", {}).get("eligibilityCriteria", ""),
             ps.get("identificationModule", {}).get("briefTitle", ""),
             ps.get("identificationModule", {}).get("officialTitle", ""),
             " ".join(ps.get("conditionsModule", {}).get("conditions", []) or []),
             " ".join(ps.get("conditionsModule", {}).get("keywords", []) or [])]
    return "\n".join(p for p in parts if p)


REFUTING = {("Does Not Support", "Sensitivity/Response"), ("Supports", "Resistance")}


def fda_kind(pma: str) -> str:
    return ("De Novo" if pma.startswith("DEN") else "510(k)" if pma.startswith("K") else "HDE" if pma.startswith("H")
            else "PMA" if pma.startswith("P") else "other")


def resolve_id(tok: str, src: dict, ctg_regex: str, refuted: bool = False) -> tuple[bool, str, str]:
    """(resolved, reason, record_text). For clinically-refuted rows a CIViC EID must be a refuting record
    (Does Not Support + Sensitivity/Response, or Supports + Resistance; AUD-2)."""
    try:
        if tok.startswith("CIViC:EID"):
            i = int(tok[9:])
            if i not in src["eid"].index:
                return False, "EID not in CIViC nightly", ""
            r = src["eid"].loc[i]
            if refuted and (r["evidence_direction"], r["significance"]) not in REFUTING:
                return False, f"not a refuting record ({r['evidence_direction']} / {r['significance']})", ""
            return True, "", f"{r['molecular_profile']} {r['therapies']}"
        if tok.startswith("CIViC:AID"):
            i = int(tok[9:])
            if i not in src["aid"].index:
                return False, "AID not in CIViC nightly", ""
            r = src["aid"].loc[i]
            return True, "", f"{r['molecular_profile']} {r['therapies']}"
        if tok.startswith("MOA:"):
            i = int(tok[4:])
            if i not in src["moa"]:
                return False, "assertion_id not in MOAlmanac", ""
            x = src["moa"][i]
            if x.get("_deprecated"):
                return False, "MOAlmanac assertion deprecated", ""
            return True, "", _moa_text(x)
        if tok.startswith("FDA:"):
            ri, pma = tok[4:].split("|", 1)
            ri = int(ri)
            if ri not in src["fda"].index:
                return False, "FDA row_index not in scrape", ""
            r = src["fda"].loc[ri]
            pmas = [p.strip() for p in str(r["pma_ids"]).split(";")]
            if pma not in pmas:
                return False, f"PMA {pma} not in FDA row {ri} ({r['pma_ids']})", ""
            return True, "", f"{r['Biomarker(s)']} {r['Biomarker(s) (Details)']} {r['Drug Trade Name (Generic) NDA / BLA']}"
        if tok.startswith("CTG:"):
            nct = tok[4:]
            f = EXT / "ctgov" / f"{nct}.json"
            if not f.exists() or nct not in src["ctg_man"]:
                return False, "CTG record not cached", ""
            txt = ctg_text(nct)
            if not ctg_regex or not re.search(ctg_regex, txt, flags=re.I):
                return False, f"CTG record lacks biomarker regex /{ctg_regex}/", ""
            return True, "", txt[:2000]
        if tok.startswith("PMID:"):
            pm = tok[5:]
            f = EXT / "pubmed" / f"{pm}.json"
            if not f.exists() or pm not in src["pm_man"]:
                return False, "PMID not verified (no cached esummary)", ""
            js = json.loads(f.read_text(encoding="utf-8"))
            rec = js.get("result", {}).get(pm, {})
            if not rec or "error" in rec:
                return False, "esummary error", ""
            return True, "", rec.get("title", "")
    except Exception as exc:  # malformed token
        return False, f"malformed id ({exc})", ""
    return False, "unknown id prefix", ""


def gene_tokens(feature: str, rules: dict, target_set: str) -> list[str]:
    fam, name = feature.split(":", 1)
    if feature in rules["features"]:
        g = rules["features"][feature].get("gene", "")
    elif fam == "SIG":
        g = "MSI"
    else:
        g = name
    toks = [t for t in re.split(r"[;]", g) if t]
    out = set(toks) | {name.split("_")[0]}
    for k, v in ALIASES.items():
        if g == k or k in toks:
            out |= set(v)
    out |= set(target_set.split(";"))  # a record naming the target drug class is also consistent
    return sorted(out)


# ------------------------------------------------------------------------------------------------- molecular data
def load_mutations(genes: set[str]) -> tuple[pd.DataFrame, pd.Index, pd.Series]:
    """Default-entry calls for the CGT genes, profiled ModelIDs, protein-altering counts per line."""
    gfile = MUT_CACHE.with_suffix(".genes.json")
    if MUT_CACHE.exists() and MUT_COUNT_CACHE.exists() and gfile.exists() and             MUT_CACHE.stat().st_mtime >= MUT_FILE.stat().st_mtime and genes <= set(json.loads(gfile.read_text())):
        m = pd.read_parquet(MUT_CACHE)
        cnt = pd.read_parquet(MUT_COUNT_CACHE)["n"]
        if "Hotspot" in m.columns:
            return m, cnt.index, cnt
    import pyarrow.csv as pc
    cols = ["ModelID", "IsDefaultEntryForModel", "HugoSymbol", "ProteinChange", "VariantType", "VariantInfo", "Exon",
            "Intron", "DNAChange", "AF", "LikelyLoF", "VepImpact", "Hotspot"]
    tb = pc.read_csv(MUT_FILE, convert_options=pc.ConvertOptions(include_columns=cols))
    m = tb.to_pandas()
    m = m[m["IsDefaultEntryForModel"] == "Yes"]
    pa = m["VepImpact"].isin(["HIGH", "MODERATE"])
    cnt = m[pa].groupby("ModelID").size().reindex(m["ModelID"].unique()).fillna(0).astype(int)
    m = m[m["HugoSymbol"].isin(genes)].drop(columns="IsDefaultEntryForModel").reset_index(drop=True)
    MUT_CACHE.parent.mkdir(parents=True, exist_ok=True)
    m.to_parquet(MUT_CACHE)
    cnt.rename("n").to_frame().to_parquet(MUT_COUNT_CACHE)
    gfile.write_text(json.dumps(sorted(genes)))
    return m, cnt.index, cnt


def mut_class(m: pd.DataFrame) -> pd.Series:
    pc = m["ProteinChange"].fillna("")
    vi = m["VariantInfo"].fillna("")
    cls = pd.Series("other", index=m.index)
    cls[pc.str.match(r"^p\.[A-Z]\d+[A-Z]$")] = "missense"
    cls[pc.str.contains(r"(del|ins|dup)", regex=True) & ~pc.str.contains("fs")] = "inframe"
    cls[pc.str.contains("fs")] = "frameshift"
    cls[(cls == "other") & vi.str.contains(r"splice_(donor|acceptor)", regex=True)] = "splice"
    return cls


def fusion_select(f: pd.DataFrame, r: dict) -> pd.Series:
    """Rows of the filtered fusion table (g5 = 5-prime, g3 = 3-prime, self-fusions already removed) matching a
    fusion rule; partners == "*" accepts any partner gene."""
    drv = r["driver"] if isinstance(r["driver"], list) else [r["driver"]]
    own, other = ("g3", "g5") if r["driver_end"] == "3p" else ("g5", "g3")
    sel = f[own].isin(drv) & ~f[other].isin(drv)
    if r["partners"] != "*":
        sel &= f[other].isin(r["partners"])
    if r.get("exclude_partner_regex"):
        sel &= ~f[other].str.contains(r["exclude_partner_regex"], regex=True)
    if r.get("exclude_partner_without_ensembl"):
        col = "Gene2" if other == "g3" else "Gene1"
        sel &= ~f[col].str.contains(r"\(\.\)", regex=True)
    if r.get("exclude_partner_within_bp"):
        co = _coords_by_symbol()
        w = int(r["exclude_partner_within_bp"])
        for dsym in drv:
            if dsym not in co.index:
                continue
            ch, st, en = co.loc[dsym, ["chrom", "start", "end"]]
            pc = co.reindex(f[other])
            near = (pc["chrom"].to_numpy() == ch) & (pc["end"].to_numpy() >= st - w) & (pc["start"].to_numpy() <= en + w)
            sel &= ~pd.Series(near, index=f.index)
    return sel


_COORDS = None


def _coords_by_symbol() -> pd.DataFrame:
    global _COORDS
    if _COORDS is None:
        _COORDS = genesets.gene_coords().drop_duplicates("symbol").set_index("symbol")[["chrom", "start", "end"]]
    return _COORDS


class Molecular:
    def __init__(self, rules: dict, pop: pd.DataFrame, mt: pd.DataFrame):
        self.rules, self.P, self.pop, self.mt = rules, rules["params"], pop, mt
        self.ids = pop.index
        genes = set()
        for r in rules["features"].values():
            genes |= {g for g in str(r.get("gene", "")).split(";") if g}
        genes |= {"CDKN2A", "CCND1", "TP53", "BRCA1", "BRCA2", "PTEN", "RB1", "ATM", "SMARCA4", "MTAP", "ERBB2", "ESR1", "CCNE1", "EGFR"}
        self.mut, self.mut_profiled, self.pa_count = load_mutations(genes)
        self.mut["cls"] = mut_class(self.mut)
        self.mut["codon"] = pd.to_numeric(self.mut["ProteinChange"].str.extract(r"^p\.[A-Z\*](\d+)")[0], errors="coerce")
        self.cache: dict[str, pd.Series] = {}
        self._fus = None
        self._msi = None

    # -- helpers --
    def _ez(self, sym):
        e = p1.entrez_of([sym]).iloc[0]
        return None if pd.isna(e) else str(e)

    def _mat(self, name, sym):
        ez = self._ez(sym)
        if ez is None:
            return pd.Series(np.nan, index=self.ids)
        df = p1.raw_matrix(name, columns=[ez], ids=self.ids, strict=False)
        if ez not in df.columns:
            return pd.Series(np.nan, index=self.ids)
        return df[ez].reindex(self.ids).astype(float)

    def cn(self, sym):
        return self._mat("cn", sym)

    def dmg(self, sym):
        return self._mat("damaging", sym)

    def expr(self, sym):
        return self._mat("expression", sym)

    def hot(self, sym):
        return self._mat("hotspot", sym)

    def _mprof(self, s: pd.Series) -> pd.Series:
        """carriers -> 1/0 on mutation-profiled lines, NaN elsewhere."""
        out = pd.Series(np.nan, index=self.ids)
        prof = self.ids.intersection(self.mut_profiled)
        out[prof] = 0.0
        hit = prof.intersection(s[s].index) if s.dtype == bool else prof.intersection(s)
        out[hit] = 1.0
        return out

    def sil(self, sym):
        x = self.expr(sym)
        lin = self.pop["OncotreeLineage"]
        dev = self.pop["role"] == "dev"
        res = pd.Series(np.nan, index=self.ids)
        ok = x.notna()
        res[ok] = 0.0
        cnt = lin[dev & ok].value_counts()
        big = set(cnt[cnt >= self.P["sil_min_lineage_lines"]].index)
        grp = lin.where(lin.isin(big), "_pan")
        for g in grp[ok].unique():
            members = ok & (grp == g)
            ref = members & dev if g != "_pan" else ok & dev
            vals = x[ref]
            if vals.empty:
                continue
            expressed = (vals >= 1.0).mean()
            med = vals.median()
            if expressed < self.P["sil_lineage_expressed_frac"]:
                continue
            hit = members & (x < self.P["sil_log2tpm_max"]) & (x < med - self.P["sil_below_lineage_median_by"])
            res[hit] = 1.0
        return res

    def lofbi(self, sym):
        d, c = self.dmg(sym), self.cn(sym)
        s = self.sil(sym)
        m = self.mut[(self.mut["HugoSymbol"] == sym) & (self.mut["LikelyLoF"].astype(str) == "True") &
                     (self.mut["AF"] >= self.P["lofbi_af_min"])]["ModelID"].unique()
        out = pd.Series(np.nan, index=self.ids)
        known = d.notna()
        out[known] = 0.0
        hit = known & ((d >= 2) | ((d >= 1) & ((c < self.P["cn_loss"]) | (s == 1))) | self.ids.isin(m))
        hit |= c < self.P["cn_del"]
        out[hit] = 1.0
        return out

    def fusions(self):
        if self._fus is None:
            f = p1.raw_matrix("fusion", columns=["ModelID", "Gene1", "Gene2", "FFPM", "IsDefaultEntryForModel"])
            f = f[(f["IsDefaultEntryForModel"] == "Yes") & (f["FFPM"] >= self.P["ffpm_min"])].copy()
            f["g5"] = f["Gene1"].str.split(" ").str[0]
            f["g3"] = f["Gene2"].str.split(" ").str[0]
            self._fus = f[f["g5"] != f["g3"]]
            self._fus_prof = pd.Index(p1.raw_matrix("fusion", columns=["ModelID"])["ModelID"].unique())
        return self._fus

    def msi(self):
        if self._msi is None:
            self._msi = molecular.msi_status(self.ids)["MSI"].reindex(self.ids)
        return self._msi

    def hypermutated(self):
        dev = self.pop.index[self.pop["role"] == "dev"]
        c = self.pa_count.reindex(self.ids)
        thr = c.reindex(dev).dropna().quantile(1 - self.P["hypermutated_top_frac"])
        return (self.msi() == 1) | (c >= thr)

    def subtype_tokens(self):
        st = self.mt["ModelSubtypeFeatures"].reindex(self.ids).fillna("")
        return st.map(lambda s: {t.strip() for t in re.split(r"\s*(?:,|;|\band\b)\s*", s) if t.strip()})

    # -- feature evaluation --
    def feature(self, fid: str) -> pd.Series:
        if fid in self.cache:
            return self.cache[fid]
        fam, name = fid.split(":", 1)
        R = self.rules["features"]
        if fid in R:
            v = self._rule(fid, R[fid])
        elif fam == "SIG" and name == "MSI":
            v = self.msi()
        elif fam == "AMP":
            c = self.cn(name); v = (c > self.P["cn_amp"]).astype(float).where(c.notna())
        elif fam == "DEL":
            c = self.cn(name); v = (c < self.P["cn_del"]).astype(float).where(c.notna())
        elif fam == "LOSS":
            c = self.cn(name); v = (c < self.P["cn_loss"]).astype(float).where(c.notna())
        elif fam == "LOFBI":
            v = self.lofbi(name)
        elif fam == "HOT":
            h = self.hot(name); v = (h > 0).astype(float).where(h.notna())
        else:
            raise KeyError(f"no rule for {fid}")
        self.cache[fid] = v.reindex(self.ids).astype(float)
        return self.cache[fid]

    def _rule(self, fid, r) -> pd.Series:
        t = r["type"]
        if t == "mutation":
            m = self.mut[self.mut["HugoSymbol"] == r["gene"]]
            keep = pd.Series(True, index=m.index)
            if "classes" in r:
                keep &= m["cls"].isin(r["classes"])
            if "protein_regex" in r:
                keep &= m["ProteinChange"].fillna("").str.contains(r["protein_regex"], regex=True)
            if "codon_min" in r:
                keep &= m["codon"].between(r["codon_min"], r["codon_max"])
            if "exon_regex" in r:
                keep &= m["Exon"].fillna("").astype(str).str.contains(r["exon_regex"], regex=True)
            if "dna_regex" in r:
                keep &= m["DNAChange"].fillna("").str.contains(r["dna_regex"], regex=True)
            if "hotspot_flag_or_regex" in r:
                hs = m["Hotspot"].astype(str).eq("True")
                keep &= hs | m["ProteinChange"].fillna("").str.contains(r["hotspot_flag_or_regex"], regex=True)
            ids = set(m.loc[keep, "ModelID"])
            if "merge_adjacent_frameshifts" in r:
                ids |= self._merged_inframe(m, r)
            return self._mprof(pd.Index(sorted(ids)))
        if t == "union":
            parts = [self.feature(f) for f in r["of"]]
            df = pd.concat(parts, axis=1)
            out = df.max(axis=1, skipna=True)
            return out.where(df.notna().any(axis=1))
        if t == "difference":
            a, b = self.feature(r["of"][0]), self.feature(r["of"][1])
            return (a * (1 - b.fillna(0))).where(a.notna())
        if t == "fusion":
            f = self.fusions()
            drv = r["driver"] if isinstance(r["driver"], list) else [r["driver"]]
            sel = fusion_select(f, r)
            out = pd.Series(np.nan, index=self.ids)
            prof = self.ids.intersection(self._fus_prof)
            out[prof] = 0.0
            out[prof.intersection(pd.Index(f.loc[sel, "ModelID"].unique()))] = 1.0
            return out
        if t == "her2_pos":
            c = self.cn("ERBB2")
            tok = self.subtype_tokens().map(lambda s: bool(s & set(self.P["her2_subtype_pos_tokens"])))
            out = (c > self.P["cn_amp"]).astype(float).where(c.notna())
            out[tok] = 1.0
            return out
        if t == "er_pos":
            return self._er_pos()
        if t == "er_pos_v3":
            return self._er_pos_v3()
        if t == "del_or_sil":
            c, s = self.cn(r["gene"]), self.sil(r["gene"])
            out = pd.concat([(c < self.P["cn_del"]).astype(float).where(c.notna()), s], axis=1)
            return out.max(axis=1, skipna=True).where(out.notna().any(axis=1))
        if t == "lofbi":
            return self.lofbi(r["gene"])
        if t == "brca_biallelic":
            parts = []
            for g in ("BRCA1", "BRCA2"):
                d, c = self.dmg(g), self.cn(g)
                af = self.mut[(self.mut["HugoSymbol"] == g) & (self.mut["LikelyLoF"].astype(str) == "True") &
                              (self.mut["AF"] >= self.P["lofbi_af_min"])]["ModelID"].unique()
                hit = (d >= 2) | ((d >= 1) & ((c < self.P["cn_loss"]) | self.ids.isin(af)))
                parts.append(hit.astype(float).where(d.notna()))
            parts.append(self.sil("BRCA1"))
            df = pd.concat(parts, axis=1)
            out = df.max(axis=1, skipna=True).where(df.notna().any(axis=1))
            out[self.hypermutated() & out.notna()] = 0.0
            return out
        if t == "cdk46_sel":
            c, d, cc = self.cn("CDKN2A"), self.dmg("CDKN2A"), self.cn("CCND1")
            df = pd.concat([(c < self.P["cn_del"]).astype(float).where(c.notna()), (d >= 1).astype(float).where(d.notna()),
                            (cc > self.P["cn_amp"]).astype(float).where(cc.notna())], axis=1)
            return df.max(axis=1, skipna=True).where(df.notna().any(axis=1))
        if t == "cn_below":
            c = self.cn(r["cn_gene"])
            return (c < self.P[r["threshold_key"]]).astype(float).where(c.notna())
        if t == "not_assayable":
            return pd.Series(np.nan, index=self.ids)
        raise KeyError(f"unknown rule type {t} for {fid}")

    def _merged_inframe(self, m: pd.DataFrame, r: dict) -> set:
        """Lines with >= 2 frameshift deletions in the rule's exon within max_gap_bp whose summed deleted length
        is a multiple of 3 (a complex in-frame indel split into two calls)."""
        gap = r["merge_adjacent_frameshifts"]["max_gap_bp"]
        x = m[(m["cls"] == "frameshift") & m["Exon"].fillna("").astype(str).str.contains(r.get("exon_regex", "."), regex=True)]
        pos = x["DNAChange"].fillna("").str.extract(r"c\.(\d+)(?:_(\d+))?del$")
        x = x.assign(a=pd.to_numeric(pos[0], errors="coerce"), b=pd.to_numeric(pos[1], errors="coerce")).dropna(subset=["a"])
        x["b"] = x["b"].fillna(x["a"])
        out = set()
        for mid, g in x.groupby("ModelID"):
            g = g.sort_values("a")
            if len(g) < 2:
                continue
            span_ok = (g["a"].iloc[1:].to_numpy() - g["b"].iloc[:-1].to_numpy() <= gap).all()
            if span_ok and int((g["b"] - g["a"] + 1).sum()) % 3 == 0:
                out.add(mid)
        self.merged_inframe = sorted(set(getattr(self, "merged_inframe", [])) | out)
        return out

    def _er_pos_v3(self) -> pd.Series:
        """AUD-4: Breast AND PAM50 luminal AND ESR1 log2(TPM+1) >= fixed cut; non-breast NA."""
        breast = self.pop["OncotreeLineage"] == "Breast"
        x = self.expr("ESR1")
        gs = molecular.global_signatures()
        cols = self.P["er_pam50_columns"]
        lum = pd.Series(np.nan, index=self.ids)
        if gs is not None and all(c in gs.columns for c in cols):
            sc = gs.reindex(self.ids)[cols].dropna(how="any")
            arg = sc.idxmax(axis=1).str.replace("Breast_PAM50_", "").str.replace("_score", "")
            lum[sc.index] = arg.isin(self.P["er_pam50_luminal"]).astype(float)
        hi = (x >= self.P["er_esr1_min_log2tpm"]).astype(float).where(x.notna())
        out = (hi * lum).where(breast & hi.notna() & lum.notna())
        self.er_expression_rule = hi.where(breast)
        self.er_pam50_luminal = lum
        self.er_d9 = {"breast_lines": int(breast.sum()), "breast_without_pam50_or_expr": int((breast & out.isna()).sum()),
                      "non_breast": "NA (AUD-4)"}
        return out

    def _er_pos(self) -> pd.Series:
        from sklearn.mixture import GaussianMixture
        tok = self.subtype_tokens()
        pos = tok.map(lambda s: bool(s & set(self.P["er_subtype_pos_tokens"])))
        neg = tok.map(lambda s: bool(s & set(self.P["er_subtype_neg_tokens"])))
        lin = self.pop["OncotreeLineage"]
        inl = lin.isin(self.P["er_mixture_lineages"])
        x = self.expr("ESR1")
        out = pd.Series(0.0, index=self.ids)
        fit = inl & x.notna() & (self.pop["role"] == "dev")
        seed = json.loads((CFG / "seeds.json").read_text(encoding="utf-8"))
        seed = int(seed["seeds"][0])
        gm = GaussianMixture(n_components=2, random_state=seed).fit(x[fit].to_numpy().reshape(-1, 1))
        hi = int(np.argmax(gm.means_.ravel()))
        ap = inl & x.notna()
        post = pd.Series(gm.predict_proba(x[ap].to_numpy().reshape(-1, 1))[:, hi], index=x[ap].index)
        out[post.index] = (post > self.P["er_mixture_posterior_cut"]).astype(float)
        self.er_gmm_only = out.copy()
        out[inl & x.isna()] = np.nan
        out[pos] = 1.0
        out[neg & ~pos] = 0.0
        self.er_fit = {"n_fit": int(fit.sum()), "means": [float(v) for v in gm.means_.ravel()], "high_component": hi}
        self.er_expression_rule = out.copy()
        # PI ruling D9: breast lines additionally require PAM50 luminal (argmax in {LumA, LumB})
        gate = lin == self.P["er_pam50_gate_lineage"]
        gs = molecular.global_signatures()
        cols = self.P["er_pam50_columns"]
        lum = pd.Series(np.nan, index=self.ids)
        if gs is not None and all(c in gs.columns for c in cols):
            sc = gs.reindex(self.ids)[cols].dropna(how="any")
            arg = sc.idxmax(axis=1).str.replace("Breast_PAM50_", "").str.replace("_score", "")
            lum[sc.index] = arg.isin(self.P["er_pam50_luminal"]).astype(float)
        self.er_pam50_luminal = lum
        out[gate] = (out[gate] * lum[gate]).where(lum[gate].notna() & out[gate].notna())
        self.er_d9 = {"breast_lines": int(gate.sum()), "breast_without_pam50": int((gate & lum.isna()).sum()),
                      "non_breast_er_pos_expression_rule_only": int(((~gate) & (out == 1)).sum()),
                      "flag_non_breast": "ER+ outside breast uses the expression rule only; no CGT test uses ER+ outside breast"}
        return out


def _two_by_two(a: pd.Series, b: pd.Series, la: str, lb: str) -> dict:
    j = pd.DataFrame({"a": a, "b": b}).dropna()
    t = {f"{la}={x}|{lb}={y}": int(((j["a"] == x) & (j["b"] == y)).sum()) for x in (1, 0) for y in (1, 0)}
    n = len(j)
    return {"n": int(n), "table": t, "agreement": (float((j["a"] == j["b"]).mean()) if n else None)}


def er_concordance(mol: "Molecular", pop: pd.DataFrame) -> dict:
    """PI ruling 3: concordance of the frozen ER_POS rule with PAM50 luminal (argmax of OmicsGlobalSignatures
    Breast_PAM50_* in {LumA, LumB}) and with ModelSubtypeFeatures 'ER+', on breast lines (D and R). Report only."""
    br = pop.index[pop["OncotreeLineage"] == "Breast"]
    er = mol.feature("CLIN:ER_POS").reindex(br)
    gmm = getattr(mol, "er_gmm_only", pd.Series(np.nan, index=mol.ids)).reindex(br)
    out = {"n_breast_lines": int(len(br)), "n_er_pos": int((er == 1).sum()),
           "n_er_pos_dev": int((er.reindex(pop.index[(pop["OncotreeLineage"] == "Breast") & (pop["role"] == "dev")]) == 1).sum())}
    gs = molecular.global_signatures()
    cols = ["Breast_PAM50_Basal_score", "Breast_PAM50_Her2_score", "Breast_PAM50_LumA_score", "Breast_PAM50_LumB_score",
            "Breast_PAM50_Normal_score"]
    if gs is not None and all(c in gs.columns for c in cols):
        sc = gs.reindex(br)[cols].dropna(how="any")
        arg = sc.idxmax(axis=1).str.replace("Breast_PAM50_", "").str.replace("_score", "")
        lum = arg.isin(["LumA", "LumB"]).astype(float).reindex(br)
        out["n_pam50_scored"] = int(len(sc))
        out["pam50_call_counts"] = {k: int(v) for k, v in arg.value_counts().items()}
        out["er_pos_vs_pam50_luminal"] = _two_by_two(er, lum, "ER_POS", "PAM50_luminal")
        out["gmm_only_vs_pam50_luminal"] = _two_by_two(gmm, lum, "GMM", "PAM50_luminal")
    else:
        out["pam50"] = "OmicsGlobalSignatures PAM50 columns absent"
    tok = mol.subtype_tokens().reindex(br)
    erplus = tok.map(lambda t: "ER+" in t).astype(float)
    out["n_subtype_ERplus"] = int(erplus.sum())
    out["er_pos_vs_subtype_ERplus"] = _two_by_two(er, erplus, "ER_POS", "subtype_ER+")
    out["gmm_only_vs_subtype_ERplus"] = _two_by_two(gmm, erplus, "GMM", "subtype_ER+")
    out["note"] = ("subtype 'ER+' absent is treated as 0 (unannotated lines included). ER_POS uses ModelSubtypeFeatures "
                   "first, so ER_POS=0|subtype=1 is 0 by construction; GMM-only rows show the mixture alone. "
                   "Flag only: the frozen rule is unchanged.")
    return out


# ------------------------------------------------------------------------------------------------- comparator
def indication_mask(codes: str, mt: pd.DataFrame) -> pd.Series:
    m = pd.Series(False, index=mt.index)
    for c in [c for c in codes.split(";") if c]:
        if c not in IND:
            raise KeyError(f"indication code {c!r} not in {INDICATION_MAP.name}")
        neg = [v for f, v in IND[c] if f == "not_OncotreeLineage"]
        if neg:
            m |= ~mt["OncotreeLineage"].isin(neg)
        for f, v in IND[c]:
            if f != "not_OncotreeLineage":
                m |= mt[f] == v
    return m


def comparator_mask(row, mt, mol: Molecular) -> pd.Series:
    comp = str(row["comparator"])
    if comp == "listed-only":
        return pd.Series(False, index=mt.index)
    if comp.startswith("pan-cancer"):
        m = pd.Series(True, index=mt.index)
        codes = [c for c in str(row["indication_codes"]).split(";") if c]
        if "SOLID" in codes:
            m &= indication_mask("SOLID", mt)
    else:
        m = indication_mask(str(row["indication_codes"]), mt)
    for tok in [t for t in str(row.get("comparator_restrict") or "").split(";") if t]:
        m &= mol.feature(tok).reindex(mt.index) == 1
    for tok in [t for t in str(row.get("comparator_exclude") or "").split(";") if t]:
        if tok.startswith("IND:"):
            m &= ~indication_mask(tok[4:], mt)
        else:
            m &= ~(mol.feature(tok).reindex(mt.index) == 1)
    return m


# ------------------------------------------------------------------------------------------------- main
def feature_d(feature: str, rules: dict) -> str:
    fam, name = feature.split(":", 1)
    if feature in rules["features"]:
        return "" if rules["features"][feature]["type"] == "not_assayable" else feature
    if fam == "SIG":
        return feature
    if fam in STANDARD_FAMILIES:
        ez = p1.entrez_of([name]).iloc[0]
        return "" if pd.isna(ez) else f"{fam}:{ez}"
    return ""


def norm_drugs(s: str) -> str:
    seen = []
    for d in str(s).split(";"):
        d = d.strip().lower()
        if d and d not in seen:
            seen.append(d)
    return ";".join(sorted(seen))


def main() -> None:
    rules = json.loads(RULES.read_text(encoding="utf-8"))
    mp = pd.read_csv(MAPPING, dtype=str, keep_default_na=False)
    corr = pd.read_csv(CORRECTIONS, dtype=str, keep_default_na=False)
    src = load_sources()
    log(f"{len(mp)} curated rows")

    # (1) resolve ids
    kept, dropped, id_rows = [], [], []
    for _, r in mp.iterrows():
        toks = [t.strip() for t in r["source_record_ids"].split(";") if t.strip()]
        genes = gene_tokens(r["feature"], rules, r["target_set"])
        ok_ids, reasons = [], []
        ctg_ok = False
        for t in toks:
            ok, why, text = resolve_id(t, src, r.get("ctg_biomarker_regex", ""), r["partition"] == "clinically-refuted")
            gm = ok and (t.startswith(("CTG:", "PMID:")) or any(g.lower() in text.lower() for g in genes))
            id_rows.append({"test_id": r["test_id"], "record_id": t, "resolved": ok, "reason": why, "gene_consistent": bool(gm),
                            "fda_kind": fda_kind(t.split("|", 1)[1]) if t.startswith("FDA:") else ""})
            if ok:
                ok_ids.append(t)
                ctg_ok |= t.startswith("CTG:")
            else:
                reasons.append(f"{t}: {why}")
        drop = None
        primary = [t for t in ok_ids if not (t.startswith("PMID:") and r["source_type"] != "pmid_only")]
        if not primary:
            drop = "no resolvable source record id"
        elif r["tier"] == "T3" and not ctg_ok:
            drop = "T3 row without a regex-matching cached CTG record (DESIGN 10.3)"
        elif r["source_type"] == "pmid_only" and r["partition"] != "clinically-refuted":
            drop = "pmid_only outside clinically-refuted (PI D5)"
        rr = dict(r)
        rr["source_ids"] = ";".join(ok_ids)
        rr["unresolved_ids"] = " | ".join(reasons)
        if drop:
            rr["drop_reason"] = drop
            dropped.append(rr)
        else:
            kept.append(rr)
    cg = pd.DataFrame(kept)
    log(f"kept {len(cg)}, dropped {len(dropped)}")

    # (2) collapse drugs; uniqueness
    cg["drugs"] = cg["drugs"].map(norm_drugs)
    cg["n_drugs"] = cg["drugs"].str.count(";") + 1
    dup = cg.duplicated(["feature", "target_set", "comparator"], keep=False)
    if dup.any():
        raise SystemExit(f"duplicate (feature, target_set, comparator): {cg.loc[dup, 'test_id'].tolist()}")

    # (3) carriers on dev + test lines (molecular data only)
    sp = p1.splits()
    mt = p1.model_table(["OncotreeLineage", "OncotreePrimaryDisease", "OncotreeCode", "ModelSubtypeFeatures",
                         "CellLineName", "StrippedCellLineName"])
    pop = sp[sp["role"].isin(["dev", "test"])].copy()
    pop = pop.join(mt.drop(columns=["OncotreeLineage"]), how="left")
    mtp = pop
    mol = Molecular(rules, pop, mt)
    gm_lib = p1.gene_map()
    in_lib = gm_lib["in_matrices"].fillna("").str.split(",").map(lambda x: "crispr" in x)
    lib_syms = set(gm_lib.loc[in_lib, "symbol"])

    feats = sorted(set(cg["feature"]) | {f for s in cg["comparator_restrict"] for f in s.split(";") if f} |
                   {f for s in cg["comparator_exclude"] for f in s.split(";") if f and not f.startswith("IND:")} |
                   {f for s in cg["feature_secondary"] for f in s.split(";") if f} |
                   {f for s in cg["secondary_comparator_exclude"] for f in s.split(";") if f})
    for f in feats:
        mol.feature(f)
    dev = pop["role"] == "dev"
    tst = pop["role"] == "test"
    carriers_long = []
    for f in feats:
        v = mol.feature(f)
        for mid in v.index[v == 1]:
            carriers_long.append({"feature": f, "ModelID": mid, "role": pop.at[mid, "role"]})

    tmin = 5
    cmin = 3
    s10 = json.loads((CFG / "design_params.json").read_text(encoding="utf-8")).get("clinical_s10", {})
    tmin = int(s10.get("testable_min_carriers", tmin))
    cmin = int(s10.get("ceiling_min_carriers", cmin))
    out = []
    for _, r in cg.iterrows():
        v = mol.feature(r["feature"])
        cm = comparator_mask(r, mtp, mol)
        fd = feature_d(r["feature"], rules)
        T = r["target_set"].split(";")
        tl = all(t in lib_syms for t in T)
        tl_any = any(t in lib_syms for t in T)
        rec = dict(r)
        rec.update({
            "feature_D": fd, "target_in_library": int(tl_any), "targets_missing_library": ";".join(t for t in T if t not in lib_syms),
            "n_carriers_D": int((v[dev] == 1).sum()), "n_carriers_R": int((v[tst] == 1).sum()),
            "n_profiled_D": int(v[dev].notna().sum()), "n_profiled_R": int(v[tst].notna().sum()),
            "n_comparator_D": int((cm & dev).sum()),
            "n_carriers_comparator_D": int(((v == 1) & cm & dev).sum()),
            "n_noncarriers_comparator_D": int(((v == 0) & cm & dev).sum()),
            "n_carriers_comparator_R": int(((v == 1) & cm & tst).sum()),
            "n_noncarriers_comparator_R": int(((v == 0) & cm & tst).sum()),
        })
        if not rec.get("scoring_rule"):
            rec["scoring_rule"] = ("any member, Bonferroni x |T| (family rule)" if len(T) >= 2 else "single target")
        # (4) testability (first failing reason)
        reason = ""
        if r["partition"] == "acquired-resistance" or r["comparator"] == "listed-only":
            reason = "acquired-resistance: listed only"
        elif not fd:
            reason = "feature not assayable (no feature_D)"
        elif not tl_any:
            reason = "target not in CRISPR library"
        elif rec["n_carriers_D"] < tmin:
            reason = f"< {tmin} D carriers (feature not in D)"
        elif rec["n_carriers_comparator_D"] < cmin:
            reason = f"< {cmin} D carriers in comparator"
        elif rec["n_noncarriers_comparator_D"] < cmin:
            reason = f"< {cmin} D non-carriers in comparator"
        rec["testable"] = int(reason == "")
        rec["untestable_reason"] = reason
        rec["ceiling_testable"] = int(reason == "" and r["sign"] in ("S", "W"))
        rec["p2_scored"] = int(rec["testable"] == 1 and r["partition"] in ("sensitivity", "bypass-resistance")
                               and r["ko_class"] in P2_KO_CLASSES and not r.get("headline_exclude")
                               and r["tier"] in P2_HEADLINE_TIERS)
        # AUD-7: merged-indication tests -> indication fixed effect; per-indication strata (beta per stratum downstream)
        strat = "stratified" in str(r["comparator"])
        rec["indication_fixed_effect"] = int(strat)
        if strat:
            st = {}
            for code in [c for c in str(r["indication_codes"]).split(";") if c]:
                im = indication_mask(code, mtp) & cm
                st[code] = [int(((v == 1) & im & dev).sum()), int(((v == 0) & im & dev).sum())]
            rec["indication_strata_D"] = json.dumps(st)
        else:
            rec["indication_strata_D"] = ""
        xl = re.findall(r"cross-lineage=([A-Za-z0-9_,]+)", str(r["comparator"]))
        rec["cross_lineage_codes"] = ";".join(xl)
        if xl:
            codes = xl[0].split(",")
            cc = cm & mtp["OncotreeCode"].isin(codes)
            rec["cross_lineage_codes_all_carriers"] = int(bool(((v[cc & dev] == 1).all()) and (cc & dev).any()))
        else:
            rec["cross_lineage_codes_all_carriers"] = ""
        sec = [t for t in str(r.get("secondary_comparator_exclude") or "").split(";") if t]
        if sec:
            cm2 = cm.copy()
            for tkn in sec:
                cm2 &= ~(mol.feature(tkn).reindex(mtp.index) == 1)
            rec["secondary_n_carriers_comparator_D"] = int(((v == 1) & cm2 & dev).sum())
            rec["secondary_n_noncarriers_comparator_D"] = int(((v == 0) & cm2 & dev).sum())
        else:
            rec["secondary_n_carriers_comparator_D"] = rec["secondary_n_noncarriers_comparator_D"] = ""
        out.append(rec)
    res = pd.DataFrame(out)

    # checks: corrections, known positives, ALK cap
    checks = {}
    kp = pd.read_csv(KNOWN_POS)
    name_map = {}
    for col in ("CellLineName", "StrippedCellLineName"):
        for mid, nm in mt[col].dropna().items():
            name_map.setdefault(str(nm).upper().replace("-", "").replace(" ", ""), mid)
    fus_feat = {"ALK": "CLIN:ALK_FUS", "ABL1": "CLIN:BCR_ABL1", "KMT2A": "CLIN:KMT2A_R", "NRG1": "CLIN:NRG1_FUS"}
    kp_res = []
    for _, k in kp.iterrows():
        mid = name_map.get(str(k["cell_line_name"]).upper().replace("-", "").replace(" ", ""))
        fv = mol.feature(fus_feat[k["driver_symbol"]]) if mid is not None else None
        val = None if (mid is None or mid not in fv.index) else fv.get(mid)
        if mid is not None and fv is not None and mid not in fv.index:  # line outside D/R: evaluate directly
            f = mol.fusions()
            rr = rules["features"][fus_feat[k["driver_symbol"]]]
            drv = rr["driver"] if isinstance(rr["driver"], list) else [rr["driver"]]
            x = f[f["ModelID"] == mid]
            sel = fusion_select(x, rr)
            val = float(sel.any())
        kp_res.append({"line": k["cell_line_name"], "ModelID": mid, "feature": fus_feat[k["driver_symbol"]],
                       "role": None if mid is None else str(p1.roles_of([mid]).iloc[0]), "carrier": None if val is None else bool(val == 1)})
    checks["fusion_known_positives"] = kp_res
    checks["fusion_known_positives_pass"] = all(x["carrier"] for x in kp_res)
    alk = mol.feature("CLIN:ALK_FUS")
    checks["alk_carriers_D"] = int((alk[dev] == 1).sum())
    fz = mol.fusions()
    alk_rule = rules["features"]["CLIN:ALK_FUS"]
    alk_rows = fz[fusion_select(fz, alk_rule) & fz["ModelID"].isin(alk.index[alk == 1])]
    checks["alk_carriers"] = [{"ModelID": mid, "role": pop.at[mid, "role"], "line": mt.at[mid, "StrippedCellLineName"],
                               "lineage": pop.at[mid, "OncotreeLineage"], "partners": sorted(set(g["g5"])),
                               "max_ffpm": float(g["FFPM"].max())} for mid, g in alk_rows.groupby("ModelID")]
    checks["alk_partner_rule"] = alk_rule["partners"]
    deviations = []
    cap = int(s10.get("alk_carriers_max", 10))
    if checks["alk_carriers_D"] > cap:
        deviations.append({"id": "DEV-CGT-1", "item": "ALK carriers <= 10 (R1-CL5, clinical_s10.alk_carriers_max)",
                           "observed": checks["alk_carriers_D"],
                           "ruling": "PI: the cap was a pre-data estimate, not a rule; all carriers with canonical partners kept, no tuning",
                           "check_instead": "every ALK carrier's partner is in the frozen CLIN:ALK_FUS partner rule"})
    kmt = mol.feature("CLIN:KMT2A_R")
    km_rows = fz[fusion_select(fz, rules["features"]["CLIN:KMT2A_R"]) & fz["ModelID"].isin(kmt.index[kmt == 1])]
    checks["kmt2a_carriers_D"] = int((kmt[dev] == 1).sum())
    checks["kmt2a_carriers_R"] = int((kmt[tst] == 1).sum())
    checks["kmt2a_carriers"] = [{"ModelID": mid, "role": pop.at[mid, "role"], "line": mt.at[mid, "StrippedCellLineName"],
                                 "lineage": pop.at[mid, "OncotreeLineage"], "partners": sorted(set(g["g3"]))}
                                for mid, g in km_rows.groupby("ModelID")]
    checks["er_concordance"] = er_concordance(mol, pop)
    checks["er_d9"] = dict(getattr(mol, "er_d9", {}))
    dv = pop["role"] == "dev"
    brd = (pop["OncotreeLineage"] == "Breast")
    for nm, ser in (("ESR1_ge_cut_breast", mol.er_expression_rule), ("ER_POS", mol.feature("CLIN:ER_POS")),
                    ("HR_POS_HER2_NEG", mol.feature("CLIN:HR_POS_HER2_NEG"))):
        checks["er_d9"][nm] = {"breast_D": int(((ser == 1) & brd & dv).sum()), "breast_R": int(((ser == 1) & brd & ~dv).sum())}
    esr1 = mol.expr("ESR1")
    san = []
    for want, names in ((1, rules["params"]["er_sanity_positive"]), (0, rules["params"]["er_sanity_negative"])):
        for nm in names:
            mid = name_map.get(nm.upper().replace("-", "").replace(" ", ""))
            v = mol.feature("CLIN:ER_POS").get(mid) if mid in mol.ids else None
            san.append({"line": nm, "ModelID": mid, "in_D_or_R": mid in mol.ids, "expected": want,
                        "ER_POS": None if v is None or pd.isna(v) else int(v),
                        "ESR1_log2tpm": None if mid not in mol.ids or pd.isna(esr1.get(mid)) else round(float(esr1[mid]), 2),
                        "pam50_luminal": None if mid not in mol.ids or pd.isna(mol.er_pam50_luminal.get(mid)) else int(mol.er_pam50_luminal[mid])})
    checks["er_sanity_lines"] = san
    checks["er_sanity_pass"] = all(x["ER_POS"] == x["expected"] for x in san if x["ER_POS"] is not None)
    checks["egfr_ex19_merged_frameshift_lines"] = getattr(mol, "merged_inframe", [])
    checks["er_mixture"] = getattr(mol, "er_fit", None)
    checks["corrections"] = corr["correction_id"].tolist()
    checks["unresolved_ids_in_kept_rows"] = int(sum(bool(x) for x in res["unresolved_ids"]))
    idd = pd.DataFrame(id_rows)
    checks["record_ids_total"] = int(len(idd))
    checks["record_ids_resolved"] = int(idd["resolved"].sum())
    checks["record_ids_gene_inconsistent"] = idd[idd["resolved"] & ~idd["gene_consistent"]][["test_id", "record_id"]].to_dict("records")

    # (5) outputs
    cols = ["test_id", "feature", "feature_D", "feature_secondary", "target_set", "comparator", "comparator_restrict",
            "comparator_exclude", "indication_codes", "sign", "partition", "ko_class", "scoring_rule", "tier", "cdx", "drugs",
            "n_drugs", "source_ids", "source_type", "in_gate", "tumour_agnostic", "excluded_indications",
            "allele_level_required", "proxy", "proxy_of", "correction", "target_in_library", "targets_missing_library",
            "n_carriers_D", "n_carriers_R", "n_profiled_D", "n_profiled_R", "n_comparator_D", "n_carriers_comparator_D",
            "n_noncarriers_comparator_D", "n_carriers_comparator_R", "n_noncarriers_comparator_R", "testable",
            "ceiling_testable", "p2_scored", "headline_exclude", "indication_fixed_effect", "indication_strata_D",
            "cross_lineage_codes", "cross_lineage_codes_all_carriers", "secondary_comparator_exclude",
            "secondary_n_carriers_comparator_D", "secondary_n_noncarriers_comparator_D", "untestable_reason",
            "unresolved_ids", "citation", "notes"]
    res[cols].to_csv(OUT_CGT, index=False)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    DROPPED.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(dropped, columns=list(mp.columns) + ["source_ids", "unresolved_ids", "drop_reason"]).to_csv(DROPPED, index=False)
    pd.DataFrame(carriers_long).to_csv(OUT_DIR / "cgt_feature_carriers.csv", index=False)
    idd.to_csv(OUT_DIR / "cgt_record_resolution.csv", index=False)

    head = res[res["proxy"].astype(str) == "0"]
    vc = lambda s: {str(k): int(v) for k, v in s.value_counts().sort_index().items()}  # noqa: E731
    p2_min = 15
    try:
        pe = json.loads((CFG / "primary_endpoints.json").read_text(encoding="utf-8"))
        p2_min = int(pe["P2_primary_clinical"]["min_testable_tests"])
    except Exception:
        pass
    summary = {
        "n_curated_rows": int(len(mp)), "n_tests": int(len(res)), "n_headline_tests (proxy=0)": int(len(head)),
        "n_proxy_tests": int(len(res) - len(head)), "n_dropped": len(dropped),
        "dropped": [{"test_id": d["test_id"], "reason": d["drop_reason"]} for d in dropped],
        "by_partition": vc(head["partition"]), "by_tier": vc(head["tier"]), "by_ko_class": vc(head["ko_class"]),
        "by_sign": vc(head["sign"]),
        "testable": int(head["testable"].sum()), "untestable": int((head["testable"] == 0).sum()),
        "untestable_by_reason": vc(head.loc[head["testable"] == 0, "untestable_reason"]),
        "testable_by_partition": vc(head.loc[head["testable"] == 1, "partition"]),
        "testable_by_tier": vc(head.loc[head["testable"] == 1, "tier"]),
        "testable_by_ko_class": vc(head.loc[head["testable"] == 1, "ko_class"]),
        "p2_scored_headline": int(head["p2_scored"].sum()),
        "p2_scored_headline_gate_disjoint": int(head.loc[head["in_gate"].astype(str) == "0", "p2_scored"].sum()),
        "p2_min_testable_tests": p2_min,
        "p2_status_expected": "primary" if int(head["p2_scored"].sum()) >= p2_min else f"descriptive (underpowered: {int(head['p2_scored'].sum())})",
        "ceiling_testable_headline": int(head["ceiling_testable"].sum()),
        "refuted_rows": int((head["partition"] == "clinically-refuted").sum()),
        "refuted_testable": int(((head["partition"] == "clinically-refuted") & (head["testable"] == 1)).sum()),
        "testable_SW_tests": head.loc[(head["testable"] == 1) & head["sign"].isin(["S", "W"]),
                                      ["test_id", "feature_D", "target_set", "sign", "partition", "tier", "n_carriers_D",
                                       "n_carriers_comparator_D", "n_carriers_R", "in_gate"]].to_dict("records"),
        "deviations": deviations,
        "headline_excluded": res.loc[res["headline_exclude"].astype(str) != "", ["test_id", "headline_exclude"]].to_dict("records"),
        "fda_id_kinds": vc(idd.loc[idd["fda_kind"] != "", "fda_kind"]),
        "recorded_gaps": [
            "ER+ -> endocrine therapy: no endocrine-monotherapy ER+ record in frozen CIViC/MOA/FDA (ER_POS__ESR1__BRCA rests on CIViC EID1541)",
            "avapritinib for KIT D816V systemic mastocytosis: no CIViC/MOA/FDA record (drug not listed)",
            "ribociclib/abemaciclib ER+/HER2- breast: no approval record in frozen sources (ER_POS__CDK46__BRCA lists palbociclib, T2)"],
        "checks": checks,
        "inputs": {"mapping": str(MAPPING.relative_to(ROOT)), "rules": str(RULES.relative_to(ROOT)),
                   "rules_version": rules.get("version")},
    }
    (OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    log(json.dumps({k: summary[k] for k in ("n_tests", "n_dropped", "testable", "untestable", "p2_scored_headline")}))
    if not checks["fusion_known_positives_pass"]:
        log("WARNING: fusion known-positive check failed; see summary.json checks")


if __name__ == "__main__":
    main()
