"""Gene annotation and gene lists for the dictionary (DESIGN 3): HGNC, Ensembl 116 coordinates, OncoKB cancer
genes, Part 1 paralog priors, fusion drivers, restricted gene set Gamma, partner priors.

All Part 1 files are read in place (read-only). Gene lists are written to ``configs/*.csv`` by
``scripts/build_dictionary.py`` with a ``source`` column; no gene symbol or Entrez literal lives in code.
"""
from __future__ import annotations

import functools
import gzip
import json
import re

import numpy as np
import pandas as pd

from . import p1
from . import paths as P

CACHE_DIR = P.CACHE / "dictionary"


def bp() -> dict:
    return P.load_config("build_params")


def resolve(rel: str):
    """Paths in build_params: ``data/...`` are Part 1 (read-only); everything else is project_b."""
    rel = str(rel)
    return (P.P1 / rel) if rel.startswith("data/") else (P.PROJECT / rel)


# ---------------------------------------------------------------------------------------------
# HGNC and symbol mapping
# ---------------------------------------------------------------------------------------------


@functools.lru_cache(maxsize=1)
def hgnc() -> pd.DataFrame:
    p1.install_p1_write_guard()
    h = pd.read_csv(resolve(bp()["dictionary"]["hgnc"]), sep="\t", low_memory=False,
                    usecols=["symbol", "entrez_id", "ensembl_gene_id", "location", "gene_group", "alias_symbol",
                             "prev_symbol", "locus_group"])
    h = h.dropna(subset=["entrez_id"]).copy()
    h["entrez"] = h["entrez_id"].astype(np.int64).astype(str)
    return h.drop(columns="entrez_id")


def symbol_to_entrez(symbols) -> pd.Series:
    """Symbol -> Entrez string: Part 1 gene_map first, then HGNC approved, previous and alias symbols."""
    syms = [str(s) for s in symbols]
    out = p1.entrez_of(syms).copy()
    h = hgnc()
    appr = h.drop_duplicates("symbol").set_index("symbol")["entrez"]
    miss = out.isna()
    out[miss] = appr.reindex(out.index[miss]).to_numpy()
    for col in ("prev_symbol", "alias_symbol"):
        miss = out.isna()
        if not miss.any():
            break
        e = h[["entrez", col]].dropna()
        e = e.assign(s=e[col].str.split("|")).explode("s")
        e = e.drop_duplicates("s", keep=False).set_index("s")["entrez"]  # unambiguous only
        out[miss] = e.reindex(out.index[miss]).to_numpy()
    return out


def ensembl_to_entrez() -> pd.Series:
    h = hgnc().dropna(subset=["ensembl_gene_id"])
    return h.drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id")["entrez"]


@functools.lru_cache(maxsize=1)
def symbol_of() -> pd.Series:
    gm = p1.gene_map()["symbol"]
    h = hgnc().drop_duplicates("entrez").set_index("entrez")["symbol"]
    return gm.combine_first(h)


# ---------------------------------------------------------------------------------------------
# Coordinates (Ensembl 116 GTF, GRCh38)
# ---------------------------------------------------------------------------------------------

_ATTR = re.compile(r'(gene_id|gene_name) "([^"]*)"')


def gene_coords(rebuild: bool = False) -> pd.DataFrame:
    """Entrez-indexed [chrom, start, end, symbol, ensembl] from the Ensembl 116 GTF gene records (cached)."""
    cache = CACHE_DIR / "gene_coords.parquet"
    if cache.exists() and not rebuild:
        return pd.read_parquet(cache)
    p1.install_p1_write_guard()
    rows = []
    with gzip.open(resolve(bp()["dictionary"]["gtf"]), "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)  # literal-ok: GTF has 9 tab fields
            if len(f) < 9 or f[2] != "gene":  # literal-ok: GTF has 9 tab fields
                continue
            a = dict(_ATTR.findall(f[8]))  # literal-ok: GTF attribute field
            rows.append((a.get("gene_id"), a.get("gene_name"), f[0], int(f[3]), int(f[4])))  # literal-ok: GTF start/end fields
    g = pd.DataFrame(rows, columns=["ensembl", "gene_name", "chrom", "start", "end"])
    e2z = ensembl_to_entrez()
    g["entrez"] = g["ensembl"].map(e2z)
    # genes on primary assembly only; one record per Entrez (first by chromosome order)
    g = g[g["chrom"].isin([str(i) for i in range(1, 23)] + ["X", "Y"])]  # literal-ok: autosomes 1-22
    g = g.dropna(subset=["entrez"]).drop_duplicates("entrez").set_index("entrez")
    g = g.rename(columns={"gene_name": "symbol"})[["chrom", "start", "end", "symbol", "ensembl"]]
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    g.to_parquet(cache)
    return g


# ---------------------------------------------------------------------------------------------
# Gene lists
# ---------------------------------------------------------------------------------------------


def oncokb_genes() -> pd.DataFrame:
    o = pd.read_csv(resolve(bp()["dictionary"]["oncokb_gene_list"]), sep="\t", dtype=str)
    o = o.rename(columns={"Hugo Symbol": "symbol", "Entrez Gene ID": "entrez"})
    o = o.dropna(subset=["entrez"])
    return o[["entrez", "symbol"]].drop_duplicates("entrez")


def paralog_top() -> pd.DataFrame:
    """Part 1 top-3 Ensembl Compara paralogs (P1/cache/priors/paralog_top.parquet; read-only)."""
    p1.install_p1_write_guard()
    t = pd.read_parquet(P.P1_CACHE / "priors" / "paralog_top.parquet")
    t.index = t.index.astype(str)
    return t


def top_paralogs_long(genes, k: int) -> pd.DataFrame:
    t = paralog_top()
    rows = []
    for g in genes:
        if g not in t.index:
            continue
        for r in range(1, k + 1):
            par = t.at[g, f"par{r}"]
            if par is not None and not (isinstance(par, float) and np.isnan(par)):
                rows.append((g, str(par), r, float(t.at[g, f"par{r}_score"])))
    return pd.DataFrame(rows, columns=["entrez", "partner_entrez", "rank", "score"])


def control_spec() -> pd.DataFrame:
    return pd.read_csv(resolve(bp()["controls"]["spec_file"]), dtype=str, keep_default_na=False)


def parse_feature_ref(ref: str) -> tuple[str, str]:
    """'LOSS:VPS4B' -> ('LOSS', 'VPS4B'). Named features (SIG:MSI) keep their name."""
    fam, _, name = ref.partition(":")
    return fam, name


GENE_FAMILIES = {"HOT", "LOFBI", "AMP", "DEL", "LOSS", "FUS", "SIL", "SIL-DEL", "SIL-MUT", "SIL-METH", "SIL-other",  # literal-ok: family ids (DESIGN 3)
                 "METH", "YABS", "EXPDOSE", "EXPGATE", "EXPPAR"}


def control_partner_genes() -> pd.DataFrame:
    """(dependency_symbol, partner_symbol, family) for every gene-level control feature (primary, secondary,
    fallback) of the spec."""
    s = control_spec()
    rows = []
    for _, r in s.iterrows():
        for col in ("feature_primary", "feature_secondary", "fallback_feature"):
            ref = r[col].strip()
            if not ref:
                continue
            fam, name = parse_feature_ref(ref)
            if fam in GENE_FAMILIES:
                rows.append((r["dependency_symbol"], name, fam, col))
    return pd.DataFrame(rows, columns=["dependency_symbol", "partner_symbol", "family", "column"])


def lineage_check_rows() -> pd.DataFrame:
    bpairs = pd.read_csv(P.P1_BIOMARKER_PAIRS, dtype=str)
    rel = bp()["controls"]["lineage_check_relation"]
    return bpairs[bpairs["relation"] == rel].copy()


def control_dependency_genes() -> pd.DataFrame:
    """Entrez of every positive-control and lineage-check dependency gene."""
    s = control_spec()
    lc = lineage_check_rows()
    syms = list(dict.fromkeys(list(s["dependency_symbol"]) + list(lc["dependency"])))
    ez = symbol_to_entrez(syms)
    return pd.DataFrame({"symbol": syms, "entrez": ez.to_numpy()})


def known_silenced() -> list[str]:
    return list(P.load_config("design_params")["dictionary_s3"]["known_silenced"])


def build_gamma(eligible_entrez) -> pd.DataFrame:
    """Gamma = OncoKB cancer genes U top-k paralogs of eligible genes U positive-control partners U known
    silenced (DESIGN 3). One row per gene with ';'-joined sources."""
    k = int(P.load_config("design_params")["dictionary_s3"]["n_paralogs"])
    src: dict[str, set] = {}

    def add(genes, s):
        for g in genes:
            if isinstance(g, str) and g:
                src.setdefault(g, set()).add(s)

    add(oncokb_genes()["entrez"], "oncokb_cancerGeneList")
    add(top_paralogs_long(list(eligible_entrez), k)["partner_entrez"], f"p1_paralog_top{k}_of_eligible")
    cp = control_partner_genes()
    add(symbol_to_entrez(cp["partner_symbol"]).dropna(), "positive_control_partners")
    add(symbol_to_entrez(known_silenced()).dropna(), "known_silenced")
    sym = symbol_of()
    out = pd.DataFrame({"entrez": sorted(src, key=int)})
    out["symbol"] = sym.reindex(out["entrez"]).to_numpy()
    out["source"] = [";".join(sorted(src[g])) for g in out["entrez"]]
    return out


def build_partner_priors(target_entrez) -> pd.DataFrame:
    """Per target gene: top-k Ensembl paralogs (Part 1 priors) and positive-control partners (DESIGN 9.2)."""
    k = int(P.load_config("design_params")["dictionary_s3"]["n_paralogs"])
    par = top_paralogs_long(list(target_entrez), k)
    par["relation"] = "paralog_rank" + par["rank"].astype(str)
    par["source"] = "P1 cache/priors/paralog_top.parquet (Ensembl Compara 116 protein_default)"
    cp = control_partner_genes()
    cp = cp[cp["dependency_symbol"] != cp["partner_symbol"]]
    cp = cp.assign(entrez=symbol_to_entrez(cp["dependency_symbol"]).to_numpy(),
                   partner_entrez=symbol_to_entrez(cp["partner_symbol"]).to_numpy())
    cp = cp.dropna(subset=["entrez", "partner_entrez"])
    cpr = pd.DataFrame({"entrez": cp["entrez"], "partner_entrez": cp["partner_entrez"], "rank": np.nan,
                        "score": np.nan, "relation": "control_partner:" + cp["family"],
                        "source": "configs/positive_controls_spec.csv (DESIGN 9.2)"})
    out = pd.concat([par, cpr], ignore_index=True)
    out = out.drop_duplicates(["entrez", "partner_entrez", "relation"])
    sym = symbol_of()
    out.insert(1, "symbol", sym.reindex(out["entrez"]).to_numpy())
    out.insert(3, "partner_symbol", sym.reindex(out["partner_entrez"]).to_numpy())  # literal-ok: column position
    return out.sort_values(["entrez", "relation", "partner_entrez"], key=lambda s: s.map(str)).reset_index(drop=True)


# ---------------------------------------------------------------------------------------------
# Fusion drivers
# ---------------------------------------------------------------------------------------------


def _civic_fusions() -> pd.DataFrame:
    v = pd.read_csv(P.EXTERNAL / "clinical" / "civic" / "nightly-VariantSummaries.tsv", sep="\t", low_memory=False)
    f = v[v["feature_type"] == "Fusion"]
    rows = []
    for _, r in f.iterrows():
        for side, col, st in (("5prime", "5_prime_partner", "5_prime_partner_status"),
                              ("3prime", "3_prime_partner", "3_prime_partner_status")):
            if r.get(st) == "known" and isinstance(r.get(col), str):
                rows.append((r[col], side, f"CIViC variant {r['variant_id']} ({r['feature_name']})"))
    return pd.DataFrame(rows, columns=["symbol", "side", "record"])


def _moa_fusions() -> pd.DataFrame:
    base = P.EXTERNAL / "clinical" / "moalmanac" / "extracted"
    hits = sorted(base.glob("*/molecular-oncology-almanac.json"))
    if not hits:
        return pd.DataFrame(columns=["symbol", "side", "record"])
    d = json.loads(hits[0].read_text(encoding="utf-8"))
    rows = []
    for r in d:
        if r.get("feature_type") != "Rearrangement":
            continue
        g1, g2 = (r.get("gene1") or "").strip(), (r.get("gene2") or "").strip()
        rec = f"MOAlmanac assertion {r.get('assertion_id')}"
        if g1 and g2:
            rows += [(g1, "5prime", rec), (g2, "3prime", rec)]
        elif g1:
            rows.append((g1, "any", rec))
    return pd.DataFrame(rows, columns=["symbol", "side", "record"])


def build_fusion_drivers() -> pd.DataFrame:
    """configs/fusion_drivers.csv: curated fusion partners (CIViC, MOAlmanac), OncoKB cancer genes only, that
    appear as the 3' gene of a curated fusion (or as a single-gene fusion record); orientation 3' except the
    design's 5' list (FGFR2/3). Genes curated only as 5' partners are not drivers."""
    cfg = bp()["dictionary"]
    ds = P.load_config("design_params")["dictionary_s3"]["fusion"]
    allrec = pd.concat([_civic_fusions().assign(db="CIViC"), _moa_fusions().assign(db="MOAlmanac")])
    allrec["entrez"] = symbol_to_entrez(allrec["symbol"]).to_numpy()
    allrec = allrec.dropna(subset=["entrez"])
    if cfg["fusion_driver_restrict_to_oncokb"]:
        allrec = allrec[allrec["entrez"].isin(set(oncokb_genes()["entrez"]))]
    h = hgnc().drop_duplicates("entrez").set_index("entrez")["gene_group"].fillna("")
    kin = h.str.contains(cfg["fusion_kinase_regex"], case=False, regex=True) & \
        ~h.str.contains(cfg["fusion_kinase_exclude_regex"], case=False, regex=True)
    five = set(symbol_to_entrez(list(ds["kinase_3prime_except_5prime"]) + list(cfg.get("fusion_5prime_drivers_extra", []))).dropna())
    sym = symbol_of()
    rows = []
    for ez, grp in allrec.groupby("entrez"):
        sides = set(grp["side"])
        is_k = bool(kin.get(ez, False))
        # DESIGN v3 3: the driver is the 3' gene (5' for the kinase_3prime_except_5prime list); a curated
        # single-gene record ("any") names the driver itself
        if ez in five:
            orient = "5prime"
        elif "3prime" in sides or sides == {"any"}:
            orient = "3prime"
        else:
            continue
        rows.append({"entrez": ez, "symbol": sym.get(ez, grp["symbol"].iloc[0]), "orientation": orient,
                     "is_kinase_hgnc": is_k, "curated_sides": ";".join(sorted(sides)),
                     "n_records": len(grp), "source": ";".join(sorted(set(grp["db"]))),
                     "example_record": grp["record"].iloc[0]})
    return pd.DataFrame(rows).sort_values("entrez", key=lambda s: s.astype(int)).reset_index(drop=True)
