"""ALLELE and CLIN dictionary families (DESIGN v3 3, 10.3; PI ruling: the dictionary track builds them) from the
CGT track's read-only rule files ``configs/cgt_allele_rules.json`` (+ ``configs/cgt_extra_features.csv``, same
feature set). Molecular data only: no effects or predictions are read.

Rule types implemented (schema of cgt_allele_rules.json): mutation (gene, classes, protein_regex, codon_min/max,
exon_regex, dna_regex), union, difference, fusion (driver, driver_end, partners or "*"), lofbi, her2_pos,
er_pos (incl. PI ruling D9 PAM50 gate for breast), del_or_sil, brca_biallelic, cdk46_sel, cn_below;
not_assayable features are not built (recorded). Values: 1 / 0 / NaN (not profiled), as each rule states.
"""
from __future__ import annotations

import json
import re
import warnings

import numpy as np
import pandas as pd

from . import genesets, molecular, p1
from . import paths as P


def rules() -> dict:
    return json.loads((P.CONFIGS / "cgt_allele_rules.json").read_text(encoding="utf-8"))


def extra_features() -> pd.DataFrame:
    p = P.CONFIGS / "cgt_extra_features.csv"
    return pd.read_csv(p, dtype=str, keep_default_na=False) if p.exists() else pd.DataFrame()


def mut_class(m: pd.DataFrame) -> pd.Series:
    """missense / inframe / frameshift / splice / other, per the rules file's class definitions."""
    pc = m["ProteinChange"].fillna("").astype(str)
    vi = m["VariantInfo"].fillna("").astype(str)
    cls = pd.Series("other", index=m.index)
    cls[pc.str.match(r"^p\.[A-Z]\d+[A-Z]$")] = "missense"
    cls[pc.str.contains(r"(?:del|ins|dup)", regex=True) & ~pc.str.contains("fs")] = "inframe"
    cls[pc.str.contains("fs")] = "frameshift"
    cls[(cls == "other") & vi.str.contains(r"splice_(?:donor|acceptor)", regex=True)] = "splice"
    return cls


class Evaluator:
    def __init__(self, inp):
        """``inp``: dictionary.Inputs (population rows, raw caches)."""
        self.R = rules()
        self.Pm = self.R["params"]
        self.inp = inp
        self.ids = inp.ids
        m = molecular.mutations()
        self.mut_profiled = pd.Index(m["ModelID"].unique())
        m = m[m["ModelID"].isin(self.ids)].copy()
        m["cls"] = mut_class(m)
        m["codon"] = pd.to_numeric(m["ProteinChange"].astype(str).str.extract(r"^p\.[A-Z\*](\d+)")[0], errors="coerce")
        self.mut = m
        self.cache: dict[str, pd.Series] = {}
        self._fus = None

    # ---- matrices
    def _ez(self, sym):
        e = genesets.symbol_to_entrez([sym]).iloc[0]
        return e if isinstance(e, str) else None

    def _col(self, mat: pd.DataFrame, sym) -> pd.Series:
        ez = self._ez(sym)
        if ez is None or ez not in mat.columns:
            return pd.Series(np.nan, index=self.ids)
        return mat[ez].astype(float)

    def cn(self, sym):
        return self._col(self.inp.cn, sym)

    def dmg(self, sym):
        return self._col(self.inp.dmg, sym)

    def sil(self, sym):
        from . import dictionary as Dct
        ez = self._ez(sym)
        if ez is None or ez not in self.inp.expr.columns:
            return pd.Series(np.nan, index=self.ids)
        x, _ = Dct.sil_matrix(self.inp, [ez])
        return pd.Series(x[:, 0], index=self.ids).astype(float)

    def _mprof(self, lines) -> pd.Series:
        out = pd.Series(np.nan, index=self.ids)
        prof = self.ids.intersection(self.mut_profiled)
        out[prof] = 0.0
        out[prof.intersection(pd.Index(lines))] = 1.0
        return out

    def _lof_af_lines(self, sym):
        m = self.mut[(self.mut["HugoSymbol"] == sym) & (self.mut["LikelyLoF"].astype(str) == "True")
                     & (self.mut["AF"] >= self.Pm["lofbi_af_min"])]
        return set(m["ModelID"])

    def lofbi(self, sym):
        d, c, s = self.dmg(sym), self.cn(sym), self.sil(sym)
        out = pd.Series(np.nan, index=self.ids)
        known = d.notna()
        out[known] = 0.0
        hit = known & ((d >= 2) | ((d >= 1) & ((c < self.Pm["cn_loss"]) | (s == 1))) | self.ids.isin(self._lof_af_lines(sym)))
        hit |= c < self.Pm["cn_del"]
        out[hit] = 1.0
        return out

    def fusions(self):
        if self._fus is None:
            f = p1.raw_matrix("fusion", columns=["ModelID", "Gene1", "Gene2", "FFPM", "IsDefaultEntryForModel"])
            self._fus_prof = pd.Index(f["ModelID"].unique())
            f = f[(f["IsDefaultEntryForModel"] == "Yes") & (f["FFPM"] >= self.Pm["ffpm_min"])].copy()
            f["g5"] = f["Gene1"].astype(str).str.split(" ").str[0]
            f["g3"] = f["Gene2"].astype(str).str.split(" ").str[0]
            self._fus = f[f["g5"] != f["g3"]]
        return self._fus

    def msi(self):
        return molecular.msi_status(self.ids)["MSI"]

    def hypermutated(self):
        cnt = self.mut.groupby("ModelID").size().reindex(self.ids)
        cnt[self.ids.isin(self.mut_profiled) & cnt.isna()] = 0
        dev = self.inp.pop["role"] == "dev"
        thr = cnt[dev].dropna().quantile(1 - self.Pm["hypermutated_top_frac"])
        return (self.msi() == 1) | (cnt >= thr)

    def tokens(self):
        st = p1.model_table(["ModelSubtypeFeatures"])["ModelSubtypeFeatures"].reindex(self.ids).fillna("")
        return st.map(lambda s: {t.strip() for t in re.split(r"\s*(?:,|;|\band\b)\s*", s) if t.strip()})

    # ---- rules
    def feature(self, fid: str) -> pd.Series:
        if fid not in self.cache:
            self.cache[fid] = self._rule(self.R["features"][fid]).reindex(self.ids).astype(float)
        return self.cache[fid]

    def _rule(self, r) -> pd.Series:
        warnings.filterwarnings("ignore", message="This pattern is interpreted as a regular expression")
        t = r["type"]
        if t == "mutation":
            m = self.mut[self.mut["HugoSymbol"] == r["gene"]]
            keep = pd.Series(True, index=m.index)
            if "classes" in r:
                keep &= m["cls"].isin(r["classes"])
            if "protein_regex" in r:
                keep &= m["ProteinChange"].fillna("").astype(str).str.contains(r["protein_regex"], regex=True)
            if "codon_min" in r:
                keep &= m["codon"].between(r["codon_min"], r["codon_max"])
            if "exon_regex" in r:
                keep &= m["Exon"].fillna("").astype(str).str.contains(r["exon_regex"], regex=True)
            if "dna_regex" in r:
                keep &= m["DNAChange"].fillna("").astype(str).str.contains(r["dna_regex"], regex=True)
            if "hotspot_flag_or_regex" in r:
                hs = m["Hotspot"].astype(str).eq("True")
                keep &= hs | m["ProteinChange"].fillna("").astype(str).str.contains(r["hotspot_flag_or_regex"], regex=True)
            lines = set(m.loc[keep, "ModelID"])
            if "merge_adjacent_frameshifts" in r:
                lines |= self._merged_inframe(m, r)
            return self._mprof(sorted(lines))
        if t == "union":
            df = pd.concat([self.feature(f) for f in r["of"]], axis=1)
            return df.max(axis=1, skipna=True).where(df.notna().any(axis=1))
        if t == "difference":
            a, b = self.feature(r["of"][0]), self.feature(r["of"][1])
            return (a * (1 - b.fillna(0))).where(a.notna())
        if t == "fusion":
            f = self.fusions()
            drv = r["driver"] if isinstance(r["driver"], list) else [r["driver"]]
            own, other = ("g3", "g5") if r["driver_end"] == "3p" else ("g5", "g3")
            sel = f[own].isin(drv) & ~f[other].isin(drv)
            if r["partners"] != "*":
                sel &= f[other].isin(r["partners"])
            if r.get("exclude_partner_regex"):
                sel &= ~f[other].str.contains(r["exclude_partner_regex"], regex=True)
            if r.get("exclude_partner_without_ensembl"):
                col = "Gene2" if other == "g3" else "Gene1"
                sel &= ~f[col].astype(str).str.contains(r"\(\.\)", regex=True)
            if r.get("exclude_partner_within_bp"):
                co = genesets.gene_coords().drop_duplicates("symbol").set_index("symbol")[["chrom", "start", "end"]]
                w = int(r["exclude_partner_within_bp"])
                for dsym in drv:
                    if dsym not in co.index:
                        continue
                    ch, st, en = co.loc[dsym, ["chrom", "start", "end"]]
                    pc = co.reindex(f[other])
                    near = (pc["chrom"].to_numpy() == ch) & (pc["end"].to_numpy() >= st - w) & (pc["start"].to_numpy() <= en + w)
                    sel &= ~pd.Series(near, index=f.index)
            out = pd.Series(np.nan, index=self.ids)
            prof = self.ids.intersection(self._fus_prof)
            out[prof] = 0.0
            out[prof.intersection(pd.Index(f.loc[sel, "ModelID"].unique()))] = 1.0
            return out
        if t == "her2_pos":
            c = self.cn(r["gene"])
            tok = self.tokens().map(lambda s: bool(s & set(self.Pm["her2_subtype_pos_tokens"])))
            out = (c > self.Pm["cn_amp"]).astype(float).where(c.notna())
            out[tok] = 1.0
            return out
        if t == "er_pos":
            return self._er_pos(r["gene"])
        if t == "er_pos_v3":
            return self._er_pos_v3(r["gene"])
        if t == "del_or_sil":
            c, s = self.cn(r["gene"]), self.sil(r["gene"])
            df = pd.concat([(c < self.Pm["cn_del"]).astype(float).where(c.notna()), s], axis=1)
            return df.max(axis=1, skipna=True).where(df.notna().any(axis=1))
        if t == "lofbi":
            return self.lofbi(r["gene"])
        if t == "brca_biallelic":
            genes = [g for g in r["gene"].split(";") if g]
            parts = []
            for g in genes:
                d, c = self.dmg(g), self.cn(g)
                hit = (d >= 2) | ((d >= 1) & ((c < self.Pm["cn_loss"]) | self.ids.isin(self._lof_af_lines(g))))
                parts.append(hit.astype(float).where(d.notna()))
            parts.append(self.sil(genes[0]))
            df = pd.concat(parts, axis=1)
            out = df.max(axis=1, skipna=True).where(df.notna().any(axis=1))
            out[self.hypermutated() & out.notna()] = 0.0
            return out
        if t == "cdk46_sel":
            g1, g2 = [g for g in r["gene"].split(";") if g][:2]
            c, d, cc = self.cn(g1), self.dmg(g1), self.cn(g2)
            df = pd.concat([(c < self.Pm["cn_del"]).astype(float).where(c.notna()), (d >= 1).astype(float).where(d.notna()),
                            (cc > self.Pm["cn_amp"]).astype(float).where(cc.notna())], axis=1)
            return df.max(axis=1, skipna=True).where(df.notna().any(axis=1))
        if t == "cn_below":
            c = self.cn(r["cn_gene"])
            return (c < self.Pm[r["threshold_key"]]).astype(float).where(c.notna())
        raise KeyError(f"rule type {t} not built")

    def _merged_inframe(self, m: pd.DataFrame, r: dict) -> set:
        """Lines with >= 2 frameshift deletions in the rule's exon within max_gap_bp whose summed deleted length
        is a multiple of 3 (one complex in-frame indel split into two calls)."""
        gap = r["merge_adjacent_frameshifts"]["max_gap_bp"]
        x = m[(m["cls"] == "frameshift") & m["Exon"].fillna("").astype(str).str.contains(r.get("exon_regex", "."), regex=True)]
        pos = x["DNAChange"].fillna("").astype(str).str.extract(r"c\.(\d+)(?:_(\d+))?del$")
        x = x.assign(a=pd.to_numeric(pos[0], errors="coerce"), b=pd.to_numeric(pos[1], errors="coerce")).dropna(subset=["a"])
        x["b"] = x["b"].fillna(x["a"])
        out = set()
        for mid, g in x.groupby("ModelID"):
            g = g.sort_values("a")
            if len(g) < 2:
                continue
            ok = (g["a"].iloc[1:].to_numpy() - g["b"].iloc[:-1].to_numpy() <= gap).all()
            if ok and int((g["b"] - g["a"] + 1).sum()) % 3 == 0:  # literal-ok: codon length (in-frame = multiple of 3)
                out.add(mid)
        return out

    def _pam50_luminal(self) -> pd.Series:
        gs = molecular.global_signatures()
        cols = self.Pm["er_pam50_columns"]
        lum = pd.Series(np.nan, index=self.ids)
        if gs is not None and all(c in gs.columns for c in cols):
            sc = gs.reindex(self.ids)[cols].dropna(how="any")
            arg = sc.idxmax(axis=1).str.replace("Breast_PAM50_", "").str.replace("_score", "")
            lum[sc.index] = arg.isin(self.Pm["er_pam50_luminal"]).astype(float)
        return lum

    def _er_pos_v3(self, gene) -> pd.Series:
        """Gate lineage AND PAM50 luminal AND expression >= er_esr1_min_log2tpm; other lineages NA."""
        gate = self.inp.pop["OncotreeLineage"] == self.Pm["er_pam50_gate_lineage"]
        x = self._col(self.inp.expr, gene)
        lum = self._pam50_luminal()
        hi = (x >= self.Pm["er_esr1_min_log2tpm"]).astype(float).where(x.notna())
        return (hi * lum).where(gate & hi.notna() & lum.notna())

    def _er_pos(self, gene) -> pd.Series:
        from sklearn.mixture import GaussianMixture
        tok = self.tokens()
        pos = tok.map(lambda s: bool(s & set(self.Pm["er_subtype_pos_tokens"])))
        neg = tok.map(lambda s: bool(s & set(self.Pm["er_subtype_neg_tokens"])))
        lin = self.inp.pop["OncotreeLineage"]
        inl = lin.isin(self.Pm["er_mixture_lineages"])
        x = self._col(self.inp.expr, gene)
        out = pd.Series(0.0, index=self.ids)
        fit = inl & x.notna() & (self.inp.pop["role"] == "dev")
        gm = GaussianMixture(n_components=2, random_state=int(P.seeds()["seeds"][0])).fit(x[fit].to_numpy().reshape(-1, 1))
        hi = int(np.argmax(gm.means_.ravel()))
        ap = inl & x.notna()
        post = pd.Series(gm.predict_proba(x[ap].to_numpy().reshape(-1, 1))[:, hi], index=x[ap].index)
        out[post.index] = (post > self.Pm["er_mixture_posterior_cut"]).astype(float)
        out[inl & x.isna()] = np.nan
        out[pos] = 1.0
        out[neg & ~pos] = 0.0
        gate = lin == self.Pm["er_pam50_gate_lineage"]
        gs = molecular.global_signatures()
        cols = self.Pm["er_pam50_columns"]
        lum = pd.Series(np.nan, index=self.ids)
        if gs is not None and all(c in gs.columns for c in cols):
            sc = gs.reindex(self.ids)[cols].dropna(how="any")
            arg = sc.idxmax(axis=1).str.replace("Breast_PAM50_", "").str.replace("_score", "")
            lum[sc.index] = arg.isin(self.Pm["er_pam50_luminal"]).astype(float)
        out[gate] = (out[gate] * lum[gate]).where(lum[gate].notna() & out[gate].notna())
        return out


def build(inp) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """All buildable ALLELE:/CLIN: features -> (values population x features, list table, info)."""
    R = rules()
    ev = Evaluator(inp)
    ex = extra_features()
    cols, rows, skipped = {}, [], []
    for fid, r in R["features"].items():
        fam = fid.split(":", 1)[0]
        if fam not in ("ALLELE", "CLIN"):
            continue
        if r["type"] == "not_assayable":
            skipped.append({"feature_id": fid, "reason": "not_assayable: " + r.get("definition", "")})
            continue
        cols[fid] = ev.feature(fid)
        genes = [g for g in str(r.get("gene", "")).split(";") if g]
        rows.append({"feature_id": fid, "family": fam, "gene": ";".join(genes), "type": r["type"],
                     "definition": r.get("definition", ""), "source": r.get("source", "")})
    X = pd.DataFrame(cols, index=inp.ids).astype(np.float32)
    lst = pd.DataFrame(rows)
    info = {"rules_file": "configs/cgt_allele_rules.json", "rules_version": R.get("version"),
            "extra_file_rows": int(len(ex)),
            "extra_file_matches_rules": bool(len(ex) == 0 or set(ex["feature_id"]) == set(R["features"])),
            "not_built": skipped}
    return X, lst, info
