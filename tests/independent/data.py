"""Data access for the independent checker (DESIGN 13.3).

Only data loaders are used: the cached covariate table (cache/covariates.parquet, written by
scripts/build_screen_covariates.py), the cached dictionary family matrices (cache/dictionary/*.parquet) and their
metadata (configs/dictionary_v1.json, configs/dictionary_gene_features.csv, results/dictionary/locus_groups.csv),
the effect NaN mask, configs/subsamples.json, and D effects through ``guard.load_effects(scope='dev')`` only.
No engine module (stats/en/arms/blocked/null/pipeline/endpoints/audit) is imported here.
"""
from __future__ import annotations

import json
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

OUT = Path(__file__).resolve().parent / "out"


def cfg(name: str) -> dict:
    with open(ROOT / "configs" / name, encoding="utf-8") as fh:
        return json.load(fh)


@lru_cache(None)
def design_params() -> dict:
    return cfg("design_params.json")


@lru_cache(None)
def pipeline_params() -> dict:
    return cfg("pipeline_params.json")


@lru_cache(None)
def seeds() -> dict:
    return cfg("seeds.json")


@lru_cache(None)
def covariates() -> pd.DataFrame:
    """D rows of the screen-covariate cache, plus METH_absent, which the cache file does not store. It is
    derived as the loader defines it (biomark.covariates.build: line not in the RRBS gene matrix,
    cache/dictionary/rrbs_gene.parquet). Fix 4: I first derived it from 'all METH features NaN', which gives 453
    D lines instead of 456."""
    c = pd.read_parquet(ROOT / "cache" / "covariates.parquet")
    c = c[c["role"] == "dev"].copy()
    rr = pd.read_parquet(ROOT / "cache" / "dictionary" / "rrbs_gene.parquet", columns=[]).index
    c["METH_absent"] = (~c.index.isin(rr)).astype(float)
    return c


@lru_cache(None)
def primary_disease() -> pd.Series:
    """OncotreePrimaryDisease (Model.csv via biomark.p1.model_table) for D lines; levels with < 5 D lines pooled
    (confounders_s5.primary_disease_sensitivity_min_lines)."""
    from biomark import p1
    k = int(design_params()["confounders_s5"]["primary_disease_sensitivity_min_lines"])
    pdz = p1.model_table(columns=["OncotreePrimaryDisease"])["OncotreePrimaryDisease"].reindex(dev_ids())
    pdz = pdz.fillna("_unknown").astype(str)
    n = pdz.value_counts()
    return pdz.where(pdz.map(n) >= k, "_rare_pd")


def dev_ids() -> pd.Index:
    return covariates().index


@lru_cache(None)
def meta() -> dict:
    d = cfg("dictionary_v1.json")
    return {f["id"]: f for f in d["features"]}


@lru_cache(None)
def gene_features() -> pd.DataFrame:
    return pd.read_csv(ROOT / "configs" / "dictionary_gene_features.csv", dtype={"entrez": str})


@lru_cache(None)
def eligible() -> list[str]:
    e = pd.read_csv(ROOT / "configs" / "eligible_genes.csv", dtype={"entrez": str})
    return sorted(e.loc[e["eligible"].astype(bool), "entrez"].tolist())


@lru_cache(None)
def shared_ids() -> tuple:
    return tuple(f for f, m in meta().items() if m["scope"] == "shared")


def universe(g: str) -> list[str]:
    """D_g = shared features + g's gene-specific features (dictionary_gene_features.csv), minus g's self-null rows
    (DESIGN 3 / 5.1). Gene-specific ids not in the dictionary are ignored."""
    gf = gene_features()
    gf = gf[gf["entrez"] == str(g)]
    sn = set(gf.loc[gf["relation"] == "self_null", "feature_id"])
    M = meta()
    out = [f for f in shared_ids() if f not in sn]
    seen = set(out)
    for f in gf.loc[gf["relation"] != "self_null", "feature_id"]:
        if f in M and f not in seen and f not in sn:
            out.append(f)
            seen.add(f)
    return out


def forced(g: str) -> set:
    """EN forced (still penalised) terms: g's gene-specific relations (own EXPDOSE, EXPPAR, own LOSS/AMP,
    partner/paralog LOSS, YABS, control partner feature); self-null rows excluded (DESIGN 4.1.3)."""
    gf = gene_features()
    gf = gf[(gf["entrez"] == str(g)) & (gf["relation"] != "self_null")]
    return set(gf["feature_id"])


@lru_cache(None)
def family_matrix(family: str) -> pd.DataFrame:
    X = pd.read_parquet(ROOT / "cache" / "dictionary" / f"{family}.parquet")
    return X.reindex(dev_ids()).astype(np.float64)


def feature_values(features) -> pd.DataFrame:
    """D x features (float64, NaN kept)."""
    M = meta()
    by = {}
    for f in features:
        by.setdefault(M[f]["family"], []).append(f)
    parts = [family_matrix(fam)[cols] for fam, cols in by.items()]
    return pd.concat(parts, axis=1)[list(features)]


@lru_cache(None)
def own_cn() -> pd.DataFrame:
    return pd.read_parquet(ROOT / "cache" / "dictionary" / "own_cn_log2.parquet").reindex(dev_ids()).astype(np.float64)


@lru_cache(None)
def expgate() -> pd.DataFrame:
    return pd.read_parquet(ROOT / "cache" / "dictionary" / "expgate.parquet").reindex(dev_ids()).astype(np.float64)


@lru_cache(None)
def locus_members() -> dict:
    """(family, locus_group) -> set of Entrez ids (results/dictionary/locus_groups.csv)."""
    lg = pd.read_csv(ROOT / "results" / "dictionary" / "locus_groups.csv", dtype={"entrez": str})
    return {k: set(v["entrez"]) for k, v in lg.groupby(["family", "locus_group"])}


def effects(genes) -> pd.DataFrame:
    """Observed D effects (Chronos) via the guard, scope 'dev' only."""
    from biomark import guard
    E = guard.load_effects(ids=list(dev_ids()), scope="dev", purpose="independent check 13.3 (dev only)",
                           genes=[str(g) for g in genes])
    return E.astype(np.float64)


@lru_cache(None)
def subsamples() -> list:
    s = cfg("subsamples.json")
    return s["subsamples"]
