"""Task 3 comparison: my EN-obs vs the engine (frequency within indep_s13.en_freq_tol; selected set identical
excluding features within the tolerance of 0.6; beta within en_beta_rel_tol relative)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import data as DA

ENG_COLS = ["gene", "feature", "selection_freq", "sign_consistency", "maj_sign", "beta", "beta_full", "selected",
            "forced", "stage"]


def draw_genes(n=20, seed_key="seeds", seed_index=1):
    from biomark import endpoint_io as IO
    s = DA.seeds()[seed_key]
    seed = s[seed_index]
    sel = IO.read_arms(arms=["EN-obs"], columns=["gene"], selected_only=True)
    hit = sorted(set(sel["gene"].astype(str)) & set(DA.eligible()))
    non = sorted(set(DA.eligible()) - set(hit))
    rng = np.random.default_rng(seed)
    a = sorted(rng.choice(hit, n // 2, replace=False).tolist())
    b = sorted(rng.choice(non, n - n // 2, replace=False).tolist())
    return a, b, {"seed_key": f"seeds.json:{seed_key}[{seed_index}]", "seed": int(seed),
                  "n_genes_with_EN_obs_selection": len(hit)}


def engine(genes):
    from biomark import endpoint_io as IO
    e = IO.read_arms(arms=["EN-obs"], genes=genes, columns=ENG_COLS)
    e["gene"] = e["gene"].astype(str)
    lam = pd.read_csv(DA.ROOT / "results" / "benchmark" / "lambda.csv", dtype={"gene": str})
    lam = lam[(lam["arm"] == "EN-obs") & lam["gene"].isin(genes)]
    return e, lam


def tidy(res: dict) -> pd.DataFrame:
    d = pd.DataFrame({k: res[k] for k in ("features", "freq", "sign_consistency", "maj_sign", "beta_mean",
                                          "beta_full", "selected")}).rename(columns={"features": "feature"})
    d.insert(0, "gene", res["gene"])
    return d


def compare(results: list, eng: pd.DataFrame, lam: pd.DataFrame) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    p = DA.design_params()["indep_s13"]
    ftol, btol = float(p["en_freq_tol"]), float(p["en_beta_rel_tol"])
    fmin = float(DA.design_params()["surrogate_s4"]["selection_freq_min"])
    mine = pd.concat([tidy(r) for r in results], ignore_index=True)
    m = mine.merge(eng, on=["gene", "feature"], how="outer", indicator=True)
    both = m[m["_merge"] == "both"].copy()
    both["dfreq"] = (both["freq"] - both["selection_freq"]).abs()
    near = ((both["freq"] - fmin).abs() <= ftol) | ((both["selection_freq"] - fmin).abs() <= ftol)
    both["near_threshold"] = near
    # 'selected' collides: mine -> selected_x, engine -> selected_y
    both["sel_mine"] = both["selected_x"].astype(bool)
    both["sel_eng"] = both["selected_y"].fillna(False).astype(bool)
    both["sel_diff"] = both["sel_mine"] != both["sel_eng"]
    sb = both[both["sel_mine"] & both["sel_eng"]].copy()
    sb["rel_beta_mean"] = (sb["beta_mean"] - sb["beta"]).abs() / sb["beta"].abs()
    sb["rel_beta_full"] = (sb["beta_full_x"] - sb["beta_full_y"]).abs() / sb["beta_full_y"].abs()
    lam_rows = []
    for r in results:
        L = lam[lam["gene"] == r["gene"]]
        if len(L):
            lam_rows.append({"gene": r["gene"], "lam_G_mine": r["lam_G"], "lam_G_eng": float(L["lam_G"].iloc[0]),
                             "lam_X_mine": r["lam_X"], "lam_X_eng": float(L["lam_X"].iloc[0]),
                             "n_rows_mine": r["n_rows"], "n_rows_eng": int(L["n_rows"].iloc[0])})
    lamdf = pd.DataFrame(lam_rows)
    lamdf["relG"] = (lamdf["lam_G_mine"] - lamdf["lam_G_eng"]).abs() / lamdf["lam_G_eng"]
    lamdf["relX"] = (lamdf["lam_X_mine"] - lamdf["lam_X_eng"]).abs() / lamdf["lam_X_eng"]
    sel_diff_far = both["sel_diff"] & ~both["near_threshold"]
    res = {"n_pairs_compared": int(len(both)), "n_only_mine": int((m["_merge"] == "left_only").sum()),
           "n_only_engine": int((m["_merge"] == "right_only").sum()),
           "freq": {"n_within_tol": int((both["dfreq"] <= ftol + 1e-12).sum()),
                    "n_outside_tol": int((both["dfreq"] > ftol + 1e-12).sum()),
                    "max_abs_diff": float(both["dfreq"].max()),
                    "n_nonzero_either": int(((both["freq"] > 0) | (both["selection_freq"] > 0)).sum()),
                    "n_exact_equal": int((both["dfreq"] < 1e-12).sum())},
           "selected": {"n_mine": int(both["sel_mine"].sum()), "n_engine": int(both["sel_eng"].sum()),
                        "n_both": int((both["sel_mine"] & both["sel_eng"]).sum()),
                        "n_diff_total": int(both["sel_diff"].sum()),
                        "n_diff_excluding_near_threshold": int(sel_diff_far.sum()),
                        "n_near_threshold_excluded": int(near.sum())},
           "beta": {"n_compared_selected_both": int(len(sb)),
                    "mean_beta_n_within_rel_tol": int((sb["rel_beta_mean"] <= btol).sum()),
                    "mean_beta_max_rel": float(sb["rel_beta_mean"].max()) if len(sb) else None,
                    "full_beta_n_within_rel_tol": int((sb["rel_beta_full"] <= btol).sum()),
                    "full_beta_max_rel": float(sb["rel_beta_full"].max()) if len(sb) else None},
           "lambda": {"n_genes": int(len(lamdf)), "lam_G_max_rel": float(lamdf["relG"].max()),
                      "lam_X_max_rel": float(lamdf["relX"].max()),
                      "n_genes_lam_G_rel_gt_1e-6": int((lamdf["relG"] > 1e-6).sum()),
                      "n_genes_lam_X_rel_gt_1e-6": int((lamdf["relX"] > 1e-6).sum()),
                      "n_rows_mismatch": int((lamdf["n_rows_mine"] != lamdf["n_rows_eng"]).sum())}}
    res["pass_freq"] = res["freq"]["n_outside_tol"] == 0
    res["pass_selected"] = res["selected"]["n_diff_excluding_near_threshold"] == 0 and res["n_only_mine"] == 0 \
        and res["n_only_engine"] == 0
    res["pass_beta"] = len(sb) > 0 and res["beta"]["mean_beta_n_within_rel_tol"] == len(sb)
    res["pass"] = res["pass_freq"] and res["pass_selected"] and res["pass_beta"]
    return res, both, lamdf
