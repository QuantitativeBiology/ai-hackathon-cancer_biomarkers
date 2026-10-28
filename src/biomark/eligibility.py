"""Eligibility (DESIGN 2) and the "dependent" / "accurate-dependent" rules.

Eligible gene = (variable OR selective) AND m5 tier in {high, medium} AND r_within >= threshold, where
r_within is the Pearson r on D between fold-centred OOF prediction and observed effect after centring BOTH
within lineage x library-group strata (strata >= min lines; rare lineages pooled). Observed effects are
read for D lines only (``guard.load_effects(scope="dev")``).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import guard, p1
from . import paths as P


def _e() -> dict:
    return P.load_config("design_params")["eligibility_s2"]


def _dep() -> dict:
    return P.load_config("design_params")["dependency_s2"]


def r_within_strata(pred: np.ndarray, obs: np.ndarray, strata: np.ndarray, min_lines: int,
                    min_total: int = 3) -> tuple[np.ndarray, np.ndarray]:  # literal-ok: default minimum lines for a defined r (overridden from build_params)
    """Per column: Pearson r of pred vs obs after centring both within strata on the column's non-NaN rows;
    strata with < ``min_lines`` non-NaN rows are excluded for that column. Returns (r, n_used)."""
    obs = obs.astype(np.float64)
    pred = pred.astype(np.float64)
    M = ~np.isnan(obs)
    _, code = np.unique(strata, return_inverse=True)
    S = np.zeros((len(code), code.max() + 1))
    S[np.arange(len(code)), code] = 1.0
    cnt = S.T @ M                                   # K x G
    keep = cnt >= min_lines
    Mk = M & (S @ keep > 0.5)                       # row kept for column
    y = np.where(Mk, obs, 0.0)
    p = np.where(Mk, pred, 0.0)
    ck = S.T @ Mk
    with np.errstate(invalid="ignore", divide="ignore"):
        my = (S.T @ y) / ck
        mp = (S.T @ p) / ck
    my = np.nan_to_num(my)
    mp = np.nan_to_num(mp)
    yc = np.where(Mk, y - S @ my, 0.0)
    pc = np.where(Mk, p - S @ mp, 0.0)
    num = (yc * pc).sum(0)
    den = np.sqrt((yc ** 2).sum(0) * (pc ** 2).sum(0))
    n = Mk.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where((den > 0) & (n >= min_total), num / den, np.nan)
    return r, n


def compute(cov: pd.DataFrame, chunk: int = 2000) -> pd.DataFrame:  # literal-ok: memory chunk size (genes)
    """Per gene (all genes of the m5 tier table): flags, tier, r_within, n_used, eligible, reason."""
    e = _e()
    min_r = float(e["r_within_lin_lib_min"])
    min_lines = int(e["r_within_stratum_min_lines"])
    min_total = int(P.load_config("build_params")["eligibility"]["min_lines_for_r"])
    gc = p1.gene_classes()
    tiers = p1.gene_tiers(e["tier_model"])
    dev = p1.dev_ids()
    strata = (cov.loc[dev, "lineage"].astype(str) + "|" + cov.loc[dev, "library_group"].astype(str)).to_numpy()
    oof = p1.oof_predictions("m5", ids=dev)
    genes = [g for g in tiers.index if g in oof.columns]
    r = np.full(len(genes), np.nan)
    n = np.zeros(len(genes), dtype=int)
    for i in range(0, len(genes), chunk):
        gs = genes[i:i + chunk]
        pc = p1.fold_centre(oof[gs]).to_numpy()
        y = guard.load_effects(list(dev), scope="dev", genes=gs).to_numpy()
        r[i:i + chunk], n[i:i + chunk] = r_within_strata(pc, y, strata, min_lines, min_total)
    out = pd.DataFrame(index=pd.Index(genes, name="entrez"))
    out["symbol"] = tiers["symbol"].reindex(out.index) if "symbol" in tiers else p1.symbols(out.index).to_numpy()
    out["tier"] = tiers["tier"].reindex(out.index)
    cls = gc.reindex(out.index)
    flag = np.zeros(len(out), dtype=bool)
    for k in e["gene_classes_any"]:
        out[k] = cls[k].fillna(False).astype(bool)
        flag |= out[k].to_numpy()
    out["class_ok"] = flag
    out["tier_ok"] = out["tier"].isin(e["tiers"])
    out["r_beyond"] = r
    out["n_lines_r"] = n
    out["r_ok"] = out["r_beyond"] >= min_r
    out["eligible"] = out["class_ok"] & out["tier_ok"] & out["r_ok"]
    reason = np.where(~out["class_ok"], "not_variable_or_selective",
                      np.where(~out["tier_ok"], "tier_not_" + "_or_".join(e["tiers"]),
                               np.where(out["r_beyond"].isna(), "r_beyond_not_computable",
                                        np.where(~out["r_ok"], f"r_beyond<{min_r}", "eligible"))))
    out["reason"] = reason
    return out


# ---------------------------------------------------------------------------------------------
# Dependent / accurate-dependent (DESIGN 2)
# ---------------------------------------------------------------------------------------------


def dependent_rule(effects_D: pd.DataFrame, sensitivity: bool = False) -> pd.DataFrame:
    """Per gene: rule (absolute / relative), cut, frac_below. ``effects_D``: D lines x genes (NaN allowed).

    absolute: effect < threshold (-0.5; sensitivity -1.0); relative (frac_below(threshold) > 0.5): effect <
    the gene's D quantile (25th percentile)."""
    d = _dep()
    thr = float(d["threshold_sensitivity"] if sensitivity else d["threshold"])
    X = effects_D.to_numpy(dtype=np.float64)
    n = np.sum(~np.isnan(X), 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = np.nansum(X < thr, 0) / n
    q = np.nanquantile(X, float(d["relative_quantile"]), axis=0) if X.shape[0] else np.full(X.shape[1], np.nan)
    rel = frac > float(d["relative_rule_if_frac_below_gt"])
    cut = np.where(rel, q, thr)
    return pd.DataFrame({"rule": np.where(rel, "relative", "absolute"), "cut": cut, "frac_below": frac,
                         "n_obs": n}, index=effects_D.columns)


def dependent(effects: pd.DataFrame, rule: pd.DataFrame) -> pd.DataFrame:
    """Boolean (NaN where effect is NaN) dependent matrix for any lines, with cuts fitted on D (``rule``)."""
    cut = rule["cut"].reindex(effects.columns).to_numpy()
    X = effects.to_numpy(dtype=np.float64)
    out = np.where(np.isnan(X), np.nan, (X < cut).astype(float))
    return pd.DataFrame(out, index=effects.index, columns=effects.columns)


def accurate_dependent_D(dep: pd.DataFrame, pred_centred: pd.DataFrame) -> pd.DataFrame:
    """D: dependent AND fold-centred prediction ranks within the n_g most negative (n_g = # dependent lines,
    ranked among lines with non-NaN effect)."""
    out = pd.DataFrame(0.0, index=dep.index, columns=dep.columns)
    for g in dep.columns:
        d = dep[g]
        ok = d.notna()
        ng = int(d[ok].sum())
        if ng == 0:
            continue
        rk = pred_centred.loc[ok[ok].index, g].rank(method="first")
        top = rk.index[rk <= ng]
        out.loc[top, g] = d[top]
        out.loc[~ok, g] = np.nan
    return out


def accurate_threshold_R(pred_R: pd.DataFrame, n_dep_D: pd.Series, n_D: pd.Series) -> pd.Series:
    """R: prediction cut = R quantile at n_g / n_D (target-free)."""
    q = (n_dep_D / n_D).reindex(pred_R.columns).clip(0, 1)
    return pd.Series({g: float(np.quantile(pred_R[g].to_numpy(), q[g])) for g in pred_R.columns})
