"""Screen covariates of the confounder model (DESIGN 5.3), for D, R and set A, from the *screened* condition
(CRISPRScreenMap -> ScreenSequenceMap -> ModelCondition). Reads no gene effects.

Categorical columns are clean strings with explicit levels (lineage, library_group, cellformat, serumfree,
medium_family, MSI, TP53); continuous columns are float (screen_days, screen_days_dev_humagne, doubling_time,
MSIScore, DMG_BURDEN, TMB, PLOIDY). Category maps and NA rules: configs/build_params.json:covariates.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

from . import molecular, p1
from . import paths as P

CATEGORICAL = ["lineage", "library_group", "cellformat", "serumfree", "medium_family", "MSI", "TP53"]  # literal-ok: covariate name (DESIGN 5.3)
CONTINUOUS = ["screen_days", "screen_days_dev_humagne", "doubling_time", "MSIScore", "DMG_BURDEN", "TMB", "PLOIDY"]
FLAGS = ["CN_absent", "MSI_unknown", "METH_absent", "DT_unknown"]


def _cfg() -> dict:
    return P.load_config("build_params")["covariates"]


def _dp() -> dict:
    return P.load_config("design_params")


def pooled_lineage(lineage: pd.Series, roles: pd.Series) -> pd.Series:
    """Lineages with < lineage_min_lines D lines -> rare label (DESIGN 2, 5.3)."""
    e = _dp()["eligibility_s2"]
    n = lineage[roles == "dev"].value_counts()
    keep = set(n.index[n >= int(e["lineage_min_lines"])])
    lin = lineage.fillna(e["rare_lineage_label"])
    return lin.where(lin.isin(keep), e["rare_lineage_label"])


def library_group(lib: pd.Series) -> pd.Series:
    c = _cfg()
    return lib.map(c["library_map"]).fillna(c["library_other"])


def _detect_formulation_col(df: pd.DataFrame, cols) -> str | None:
    """The column whose non-empty values match the formulation regex in >= min_frac of cells (DESIGN 5.3:
    OnboardedMedia/FormulationID appear swapped)."""
    c = _cfg()
    rx = re.compile(c["formulation_detect_regex"])
    best = None
    for col in cols:
        v = df[col].dropna().astype(str)
        if len(v) == 0:
            continue
        frac = v.map(lambda s: bool(rx.search(s))).mean()
        if frac >= float(c["formulation_detect_min_frac"]) and (best is None or frac > best[1]):
            best = (col, frac)
    return None if best is None else best[0]


def medium_family_of(text) -> str | None:
    c = _cfg()
    if text is None or (isinstance(text, float) and np.isnan(text)) or str(text).strip() == "":
        return None
    base = re.split(c["medium_base_split_regex"], str(text).strip())[0].strip()
    for fam, rx in c["medium_family_regex"]:
        if re.search(rx, base, flags=re.IGNORECASE):
            return fam
    return c["medium_other"]


def _combine(fams) -> str | None:
    f = {x for x in fams if x is not None}
    if not f:
        return None
    return f.pop() if len(f) == 1 else _cfg()["medium_multi_family"]


def culture(sc: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """cellformat, serumfree, medium_family (+ medium_source) per line of ``sc`` (p1.screen_covariates)."""
    c = _cfg()
    mc = p1.model_condition(["FormulationID", "GrowthMedia", "SerumFreeMedia", "CellFormat"])
    mcol = _detect_formulation_col(mc, ["FormulationID", "GrowthMedia"])
    model = p1.model_table(["OnboardedMedia", "FormulationID", "SerumFreeMedia"])
    modcol = _detect_formulation_col(model, ["OnboardedMedia", "FormulationID"])
    out = pd.DataFrame(index=sc.index)
    # cell format of the screened condition(s)
    cf = sc["CellFormat"]
    def _cf(v):
        if not isinstance(v, str):
            return c["cellformat_missing"]
        parts = {c["cellformat_map"].get(x, c["cellformat_missing"]) for x in v.split("+")}
        return parts.pop() if len(parts) == 1 else c["cellformat_multi"]
    out["cellformat"] = cf.map(_cf)
    # serum-free: screened condition, else Model.csv, else unknown
    def _sf(v):
        if not isinstance(v, str):
            return None
        parts = set(v.split("+"))
        return c["serumfree_map"].get(parts.pop()) if len(parts) == 1 else None
    sf = sc["SerumFreeMedia"].map(_sf)
    msf = model["SerumFreeMedia"].reindex(sc.index).map(lambda v: c["serumfree_map"].get(str(v)) if pd.notna(v) else None)
    out["serumfree"] = sf.where(sf.notna(), msf).fillna(c["serumfree_missing"])
    out["serumfree_source"] = np.where(sf.notna(), "condition", np.where(msf.notna(), "model", "none"))
    # base medium: per screened condition, then Model.csv
    cond_fam = mc[mcol].map(medium_family_of) if mcol else pd.Series(dtype=object)
    fam, src = [], []
    for line, conds in sc["ModelConditionID"].items():
        f = None
        if isinstance(conds, str) and mcol:
            f = _combine([cond_fam.get(x) for x in conds.split("+")])
        if f is not None:
            fam.append(f); src.append("condition"); continue
        f = medium_family_of(model[modcol].get(line)) if modcol else None
        if f is not None:
            fam.append(f); src.append("model"); continue
        fam.append(c["medium_unknown"]); src.append("none")
    out["medium_family"] = fam
    out["medium_source"] = src
    return out, {"condition_formulation_column": mcol, "model_formulation_column": modcol}


def vif(X: pd.DataFrame) -> pd.Series:
    """VIF_j = diag(inv(corr)) for non-constant columns."""
    X = X.loc[:, X.std(ddof=0) > 0].astype(float)
    C = np.corrcoef(X.to_numpy(), rowvar=False)
    inv = np.linalg.pinv(C)
    return pd.Series(np.diag(inv), index=X.columns)


def shared_design(cov: pd.DataFrame, include_days: bool = True) -> pd.DataFrame:
    """Dummy-coded shared design (first level of each categorical dropped) + continuous columns."""
    cats = [c for c in CATEGORICAL if c in cov]
    X = pd.get_dummies(cov[cats].astype(str), drop_first=True, dtype=float)
    cont = ["DMG_BURDEN"] + (["screen_days_dev_humagne"] if include_days else [])
    if cov["doubling_time"].notna().any():
        cont += ["doubling_time", "DT_unknown"]
    return pd.concat([X, cov[cont].astype(float)], axis=1)


def build() -> tuple[pd.DataFrame, dict]:
    """Covariate table for the population (D, R, set A) and a checks dict."""
    c = _cfg()
    pop = molecular.population()
    ids = pop.index
    roles = pop["role"]
    sc = p1.screen_covariates().reindex(ids)
    cov = pd.DataFrame(index=ids)
    cov["role"] = roles
    cov["PatientID"] = pop["PatientID"]
    cov["cv_fold"] = pop["cv_fold"]
    cov["lineage_raw"] = pop["OncotreeLineage"]
    cov["lineage"] = pooled_lineage(pop["OncotreeLineage"], roles)
    cov["Library_raw"] = sc["Library"]
    cov["library_group"] = library_group(sc["Library"])
    cov["screen_days"] = sc["ScreenDays"].astype(float)
    hum = cov["library_group"] == c["humagne_group"]
    med = cov.loc[hum & (roles == "dev"), "screen_days"].median()
    cov["screen_days_dev_humagne"] = np.where(hum, cov["screen_days"] - med, 0.0)
    cul, cul_info = culture(sc)
    cov = cov.join(cul)
    # burden
    msi = molecular.msi_status(ids)
    lv = c["msi_levels"]
    cov["MSI"] = np.where(msi["MSI"] == 1, lv["MSI"], np.where(msi["MSI"] == 0, lv["MSS"], lv["unknown"]))
    cov["MSI_unknown"] = msi["MSI"].isna().astype(float)
    cov["msi_source"] = msi["msi_source"].fillna("none")
    cov["MSIScore"] = msi["MSIScore"]
    mb = p1.feature_block("M_burden")["M_burden"].reindex(ids)
    cov["DMG_BURDEN"] = mb.astype(float)
    prof = p1.raw_matrix("damaging", columns=[]).index
    summ = molecular.mutation_summaries(ids, prof)
    cov["TMB"] = molecular.tmb(summ) if len(summ.columns) else np.nan
    gs = molecular.global_signatures()
    cov["PLOIDY"] = gs["Ploidy"].reindex(ids).astype(float) if gs is not None else np.nan
    comp = p1.feature_block("composite")
    tcol = P.load_config("build_params")["dictionary"]["composite_features"]["COMP:TP53MUT"]
    t = comp[tcol].reindex(ids)
    tl = c["tp53_levels"]
    cov["TP53"] = np.where(t == 1, tl["mut"], np.where(t == 0, tl["wt"], tl["unknown"]))  # literal-ok: covariate name (DESIGN 5.3)
    # CN availability (gene-level own CN lives in cache/dictionary/own_cn_log2.parquet)
    cn_ids = p1.raw_matrix("cn", columns=[]).index
    cov["CN_absent"] = (~ids.isin(cn_ids)).astype(float)
    meth, _ = molecular.rrbs_gene_matrix()
    cov["METH_absent"] = (~ids.isin(meth.index if meth is not None else [])).astype(float)
    # doubling time and screen QC (CRISPRConfounders)
    cf = molecular.crispr_confounders()
    if cf is not None:
        dt = cf[c["doubling_time_column"]].reindex(ids).astype(float)
        cov["doubling_time_raw"] = dt
        cov["DT_unknown"] = dt.isna().astype(float)
        cov["doubling_time"] = dt.fillna(dt[roles == "dev"].median())
        for q in c["screen_qc_aux_columns"]:
            cov[f"aux_{q}"] = cf[q].reindex(ids).astype(float)
    else:
        cov["doubling_time_raw"] = np.nan
        cov["DT_unknown"] = np.nan
        cov["doubling_time"] = np.nan
    # hypermutated (DESIGN v3 dictionary_s3.hypermutated_rule): MSI or top fraction of D TMB
    frac = float(_dp()["dictionary_s3"]["hypermutated_top_frac"])
    thr = cov.loc[roles == "dev", "TMB"].quantile(1 - frac)
    cov["hypermutated"] = ((cov["MSI"] == lv["MSI"]) | (cov["TMB"] > thr)).astype(float)
    # VIF of the shared design on D; screen-days deviation dropped if VIF > max
    D = cov[roles == "dev"]
    v = vif(shared_design(D))
    dropped = []
    vmax = float(c["vif_drop_screen_days_if_gt"])
    if v.get("screen_days_dev_humagne", 0) > vmax:
        dropped.append("screen_days_dev_humagne")
        v = vif(shared_design(D, include_days=False))
    cov.attrs["dropped"] = dropped
    for col in CATEGORICAL:
        cov[col] = cov[col].astype(str)
    checks = {"culture": cul_info, "screen_days_humagne_D_median": float(med), "dropped_by_vif": dropped,
              "vif_max": float(v.max()), "vif_max_column": str(v.idxmax()),
              "hypermutated_tmb_threshold_D": float(thr)}
    # library concordance with Part 1 NaN-pattern groups (ARI check; dev lines only)
    from sklearn.metrics import adjusted_rand_score
    lg = p1.library_groups()
    j = lg.join(cov["library_group"], how="inner")
    checks["library_ari_vs_p1_lib_group"] = float(adjusted_rand_score(j["lib_group"], j["library_group"]))
    checks["library_ari_min"] = float(c["lib_concordance_ari_min"])
    checks["library_ari_pass"] = checks["library_ari_vs_p1_lib_group"] >= checks["library_ari_min"]
    return cov, {"checks": checks, "vif": v.sort_values(ascending=False)}


def summary(cov: pd.DataFrame) -> pd.DataFrame:
    """Levels per categorical and missingness per continuous column, by role."""
    rows = []
    roles = ["dev", "test", "noncancerous"]
    for col in CATEGORICAL + ["msi_source", "medium_source", "serumfree_source"]:
        if col not in cov:
            continue
        ct = pd.crosstab(cov[col], cov["role"]).reindex(columns=roles, fill_value=0)
        for lev, r in ct.iterrows():
            rows.append({"covariate": col, "kind": "categorical", "level": lev,
                         **{f"n_{k}": int(r[k]) for k in roles}})
    for col in CONTINUOUS + FLAGS + ["hypermutated"] + [c for c in cov if c.startswith("aux_")]:
        if col not in cov:
            continue
        rec = {"covariate": col, "kind": "continuous" if col in CONTINUOUS or col.startswith("aux_") else "flag",
               "level": ""}
        for k in roles:
            x = cov.loc[cov["role"] == k, col].astype(float)
            rec[f"n_{k}"] = int(x.notna().sum())
            rec[f"missing_{k}"] = int(x.isna().sum())
            rec[f"mean_{k}"] = float(x.mean()) if x.notna().any() else np.nan
        rows.append(rec)
    return pd.DataFrame(rows)
