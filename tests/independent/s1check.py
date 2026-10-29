"""Task 1: S1-obs p for 20 seeded genes vs the engine (|delta log10 p| <= indep_s13.s1_p_tol_abs_log10)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import data as DA
from .confounder import Options, s1_gene

ENG_COLS = ["gene", "feature", "p", "beta", "se", "df", "nu", "n_obs", "n_carriers", "n_carrier_clusters",
            "dropped_cov", "r2_design", "collinear_with_design", "n_lev1_excluded", "sd_dtilde", "df_rule", "t"]


def draw_genes(hit_genes, n=20, seed_key="seeds", seed_index=0):
    """Seeded stratified draw: n/2 eligible genes with >= 1 S1-obs F1 list hit (my task-2 BH), n/2 without."""
    s = DA.seeds()[seed_key]
    seed = s[seed_index] if isinstance(s, list) else s
    rng = np.random.default_rng(seed)
    elig = DA.eligible()
    hit = sorted(set(hit_genes) & set(elig))
    non = sorted(set(elig) - set(hit))
    a = rng.choice(hit, n // 2, replace=False).tolist()
    b = rng.choice(non, n - n // 2, replace=False).tolist()
    return sorted(a), sorted(b), {"seed_key": f"seeds.json:{seed_key}[{seed_index}]", "seed": int(seed)}


def compute(genes, opt=None):
    E = DA.effects(genes)
    rows = []
    for g in genes:
        U = DA.universe(g)
        rows += s1_gene(g, E[g].to_numpy(float), U, opt)
    return pd.DataFrame(rows)


def engine(genes, arm="S1-obs"):
    from biomark import endpoint_io as IO
    e = IO.read_arms(arms=[arm], genes=genes, columns=ENG_COLS)
    e["gene"] = e["gene"].astype(str)
    return e


def compare(mine: pd.DataFrame, eng: pd.DataFrame, tol: float) -> tuple[dict, pd.DataFrame]:
    m = mine.merge(eng, on=["gene", "feature"], how="outer", suffixes=("", "_eng"), indicator=True)
    both = m[m["_merge"] == "both"].copy()
    lp = np.log10(both["p"].astype(float))
    le = np.log10(both["p_eng"].astype(float))
    both["dlog10p"] = (lp - le).abs()
    both.loc[lp.isna() & le.isna(), "dlog10p"] = 0.0
    both["ok"] = both["dlog10p"] <= tol
    only_mine = m[m["_merge"] == "left_only"]
    only_eng = m[m["_merge"] == "right_only"]
    for c in ("beta", "se", "df"):
        both[f"rel_{c}"] = ((both[c] - both[f"{c}_eng"]) / both[f"{c}_eng"].abs()).abs()
    dc = (both["dropped_cov"].fillna("") != both["dropped_cov_eng"].fillna("").astype(str).map(
        lambda s: ";".join(sorted(x for x in s.replace(",", ";").replace("|", ";").split(";") if x))))
    both["drop_mismatch"] = dc
    res = {"n_compared": int(len(both)), "n_match": int(both["ok"].sum()),
           "match_rate": float(both["ok"].mean()) if len(both) else None,
           "max_dlog10p": float(both["dlog10p"].max()),
           "median_dlog10p": float(both["dlog10p"].median()),
           "n_only_mine": int(len(only_mine)), "n_only_mine_untested": int((~only_mine["tested"].astype(bool)).sum()),
           "n_only_engine": int(len(only_eng)),
           "n_df_mismatch_gt_1e-9": int((both["rel_df"] > 1e-9).sum()),
           "n_beta_rel_gt_1e-6": int((both["rel_beta"] > 1e-6).sum()),
           "n_se_rel_gt_1e-6": int((both["rel_se"] > 1e-6).sum()),
           "n_dropped_cov_mismatch": int(dc.sum()),
           "n_obs_mismatch": int((both["n_obs"] != both["n_obs_eng"]).sum())}
    return res, both
