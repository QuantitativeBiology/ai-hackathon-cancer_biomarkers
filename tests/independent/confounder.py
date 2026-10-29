"""Independent implementation of the DESIGN 5.3 confounder model and the S1 test (5.4), written from the design
text and configs only.

Reading (R0, fixed before any comparison with the engine):
* fit rows: D lines with a non-NaN effect for g and a non-NaN value of d (1.6: "in S1 ... NA lines are excluded").
* Z = intercept + design blocks (never dropped): lineage (reference ``_rare``), library_group, CellFormat,
  SerumFree, medium family (reference = most frequent level in the fit rows), ScreenDoublingTime (NA -> D median,
  as cached) + DT_unknown; plus the droppable covariates MSI (reference MSS), DMG_BURDEN, TP53 (reference WT),
  screen days (the cached 'screen days within Humagne' deviation, screen_days_dev_humagne), own CN of g (log2,
  NaN -> 0, with CN_absent as its indicator) and EXPGATE of g.  Reference levels follow
  configs/pipeline_params.json:covariates.
* leverage rule: a categorical level with < 2 lines in the fit rows is merged into the block's reference level.
* drops per (g, d) (5.3): covariate is d (MSI<->SIG:MSI, DMG_BURDEN<->SIG:DMG_BURDEN, TP53<->COMP:TP53MUT, per
  pipeline_params covariate_is_feature); association R^2 > 0.49 on the fit rows (categorical: R^2 of d on the
  covariate's level dummies; numeric: r^2); own CN when d is CN-derived (AMP/DEL/LOSS/SIL-DEL) and d's group
  (its genes or its locus group) contains g; TP53 when d has TP53 among its genes or is an AMP feature of MDM2.
* beta: OLS partial coefficient by FWL; SE: CR2 with PatientID clusters (HC2 for singletons), hat matrix of the
  full model H = H_Z + dtilde dtilde'/sum dtilde^2; rows with 1 - h_ii < 1e-8 excluded from the meat.
* df: min(n - rank(Z) - 2, nu) for continuous d (rank(Z) counts the intercept), nu = max(1, (sum W)^2/sum W^2 - 1);
  binary d (fallback triggered, A7(iv)): min(n - rank(Z) - 2, carrier clusters - 1). Two-sided t p-value.
* tests are dropped if binary carriers < 5 in the fit rows or var(dtilde) < 1e-8.
"""
from __future__ import annotations

from collections import Counter

import numpy as np
import scipy.linalg as sla
from scipy import stats

from . import data as DA

RANK_TOL = 1e-10            # pipeline_params numerics.rank_tol_rel
LEV_EPS = 1e-8              # design 5.3 / confounders_s5.leverage.exclude_from_meat_if_1_minus_h_lt
MERGE_LT = 2                # confounders_s5.leverage.merge_levels_with_lt_lines_in_fit
EIG_FLOOR_REL = 1e-12       # pipeline_params numerics.eig_floor_rel

DESIGN_CATS = [("lineage", "_rare"), ("library_group", None), ("cellformat", None), ("serumfree", None),
               ("medium_family", None)]
DESIGN_NUMS = ["doubling_time", "DT_unknown"]
DROPPABLE = ["MSI", "DMG_BURDEN", "TP53", "screen_days", "own_CN", "EXPGATE"]
COV_IS_FEATURE = {"MSI": {"SIG:MSI"}, "DMG_BURDEN": {"SIG:DMG_BURDEN"}, "TP53": {"COMP:TP53MUT"}}
CN_FAMS = {"AMP", "DEL", "LOSS", "SIL-DEL"}
TP53_EZ, MDM2_EZ = "7157", "4193"


class Options(dict):
    """Reading switches (defaults = R0). Alternatives are used only for labelled diagnostics."""
    DEFAULT = {"screen_days_col": "screen_days_dev_humagne", "na_flags_always": False, "lev_shared_only": False}

    def __init__(self, **kw):
        super().__init__(self.DEFAULT)
        self.update(kw)


# ------------------------------------------------------------------------------------------------ design coding

def cat_block(vals: np.ndarray, ref) -> tuple[np.ndarray, int]:
    """Dummy columns (levels != reference) after merging < 2-line levels into the reference."""
    cnt = Counter(vals.tolist())
    if ref is None:
        ref = sorted(cnt.items(), key=lambda kv: (-kv[1], str(kv[0])))[0][0]
    v = vals.copy()
    merged = 0
    for lev, n in cnt.items():
        if n < MERGE_LT and lev != ref:
            v[v == lev] = ref
            merged += 1
    levels = sorted(set(v.tolist()) - {ref})
    if not levels:
        return np.zeros((len(v), 0)), merged
    return np.column_stack([(v == L).astype(float) for L in levels]), merged


def cov_columns(rows: np.ndarray, g: str, name: str, opt: Options) -> np.ndarray:
    c = DA.covariates()
    if name == "MSI":
        return cat_block(c["MSI"].to_numpy()[rows].astype(str), "MSS")[0]
    if name == "TP53":
        return cat_block(c["TP53"].to_numpy()[rows].astype(str), "WT")[0]
    if name == "DMG_BURDEN":
        return c["DMG_BURDEN"].to_numpy(float)[rows, None]
    if name == "screen_days":
        return c[opt["screen_days_col"]].to_numpy(float)[rows, None]
    if name == "own_CN":
        oc = DA.own_cn()
        if g not in oc.columns:
            return np.zeros((rows.sum(), 0))
        x = np.nan_to_num(oc[g].to_numpy(float), nan=0.0)[rows, None]
        if opt["na_flags_always"]:
            return x
        return np.column_stack([x, c["CN_absent"].to_numpy(float)[rows]])
    if name == "EXPGATE":
        eg = DA.expgate()
        if g not in eg.columns:
            return np.zeros((rows.sum(), 0))
        return eg[g].to_numpy(float)[rows, None]
    raise KeyError(name)


def design_columns(rows: np.ndarray, opt: Options) -> tuple[np.ndarray, int]:
    c = DA.covariates()
    cols = [np.ones((rows.sum(), 1))]
    merged = 0
    for name, ref in DESIGN_CATS:
        X, m = cat_block(c[name].to_numpy()[rows].astype(str), ref)
        cols.append(X)
        merged += m
    for name in DESIGN_NUMS:
        cols.append(c[name].to_numpy(float)[rows, None])
    if opt["na_flags_always"]:
        for name in ("CN_absent", "MSI_unknown", "METH_absent"):
            cols.append(c[name].to_numpy(float)[rows, None])
    return np.column_stack(cols), merged


def orth(Z: np.ndarray) -> np.ndarray:
    """Orthonormal basis of col(Z) (pivoted QR, relative rank tolerance)."""
    Q, R, _ = sla.qr(Z, mode="economic", pivoting=True)
    d = np.abs(np.diag(R))
    r = int((d > RANK_TOL * d[0]).sum()) if len(d) else 0
    return Q[:, :r]


def resid(Q: np.ndarray, X: np.ndarray) -> np.ndarray:
    return X - Q @ (Q.T @ X)


# ------------------------------------------------------------------------------------------------ drops

def r2_assoc(dcol: np.ndarray, C: np.ndarray) -> float:
    """R^2 of d on [1, C] (C = covariate dummies or a numeric column)."""
    if C.shape[1] == 0:
        return 0.0
    Q = orth(np.column_stack([np.ones(len(dcol)), C]))
    dc = dcol - dcol.mean()
    sst = float(dc @ dc)
    if sst <= 0:
        return 0.0
    e = resid(Q, dcol[:, None])[:, 0]
    return 1.0 - float(e @ e) / sst


def structural_drops(g: str, f: str) -> set:
    m = DA.meta()[f]
    fam = m["family"]
    genes = set(str(x) for x in (m.get("genes") or []))
    out = set()
    for cov, feats in COV_IS_FEATURE.items():
        if f in feats:
            out.add(cov)
    if fam in CN_FAMS:
        grp = set(genes)
        lg = m.get("locus_group")
        if lg:
            grp |= DA.locus_members().get((fam, lg), set())
        if g in grp:
            out.add("own_CN")
    if TP53_EZ in genes or (fam == "AMP" and MDM2_EZ in genes):
        out.add("TP53")
    return out


def assoc_cov_column(rows, g, name, opt):
    """The covariate as used for the association rule: own CN by its numeric value only."""
    C = cov_columns(rows, g, name, opt)
    if name == "own_CN" and C.shape[1] > 1:
        C = C[:, :1]
    return C


# ------------------------------------------------------------------------------------------------ CR2 test

def cr2_test(Q: np.ndarray, y_t: np.ndarray, Dt: np.ndarray, clusters: np.ndarray, Qlev: np.ndarray | None = None):
    """Vectorised over the columns of Dt (FWL residuals of the features on Z). y_t = FWL residual of y.
    Returns beta, se, n_lev1_excluded (arrays over features). Qlev: basis used for h_Z (default Q)."""
    Ql = Q if Qlev is None else Qlev
    S = (Dt ** 2).sum(0)
    beta = (Dt * y_t[:, None]).sum(0) / S
    E = y_t[:, None] - Dt * beta[None, :]
    hZ = (Ql ** 2).sum(1)
    H = hZ[:, None] + Dt ** 2 / S[None, :]
    one_m = 1.0 - H
    ok = one_m >= LEV_EPS
    n_lev = (~ok).sum(0)
    # clusters
    uniq, inv, cnt = np.unique(clusters, return_inverse=True, return_counts=True)
    single = cnt[inv] == 1
    with np.errstate(divide="ignore", invalid="ignore"):
        u = np.where(ok & single[:, None], Dt * E / np.sqrt(np.where(ok, one_m, 1.0)), 0.0)
    meat = (u ** 2).sum(0)
    F = Dt.shape[1]
    for ci in np.where(cnt > 1)[0]:
        idx = np.where(inv == ci)[0]
        okc = ok[idx]                                   # m x F
        allok = okc.all(0)
        Qc = Ql[idx]
        Hz = Qc @ Qc.T
        dc = Dt[idx]                                    # m x F
        ec = E[idx]
        m = len(idx)
        if allok.any():
            j = np.where(allok)[0]
            Hf = Hz[None] + np.einsum("if,jf->fij", dc[:, j], dc[:, j]) / S[j][:, None, None]
            Mm = np.eye(m)[None] - Hf
            w, V = np.linalg.eigh(Mm)
            wmax = w.max(axis=1, keepdims=True)
            winv = np.where(w > EIG_FLOOR_REL * wmax, 1.0 / np.sqrt(np.clip(w, 1e-300, None)), 0.0)
            A = np.einsum("fik,fk,fjk->fij", V, winv, V)
            Ae = np.einsum("fij,jf->fi", A, ec[:, j])
            meat[j] += np.einsum("if,fi->f", dc[:, j], Ae) ** 2
        for jj in np.where(~allok)[0]:
            keep = okc[:, jj]
            if not keep.any():
                continue
            ii = np.where(keep)[0]
            Hs = Hz[np.ix_(ii, ii)] + np.outer(dc[ii, jj], dc[ii, jj]) / S[jj]
            w, V = np.linalg.eigh(np.eye(len(ii)) - Hs)
            winv = np.where(w > EIG_FLOOR_REL * w.max(), 1.0 / np.sqrt(np.clip(w, 1e-300, None)), 0.0)
            A = (V * winv) @ V.T
            meat[jj] += float(dc[ii, jj] @ (A @ ec[ii, jj])) ** 2
    se = np.sqrt(meat) / S
    return beta, se, n_lev


def nu_df(Dt: np.ndarray, clusters: np.ndarray) -> np.ndarray:
    _, inv = np.unique(clusters, return_inverse=True)
    W = np.zeros((inv.max() + 1, Dt.shape[1]))
    np.add.at(W, inv, Dt ** 2)
    return np.maximum(1.0, W.sum(0) ** 2 / (W ** 2).sum(0) - 1.0)


# ------------------------------------------------------------------------------------------------ S1 for a gene

def s1_gene(g: str, y_full: np.ndarray, features: list, opt: Options | None = None, X=None,
            drop_extra: dict | None = None, return_rows=False):
    """S1 tests of y (D vector, NaN allowed) on every feature of ``features`` for gene g.
    ``drop_extra``: feature -> set of extra covariate/block names to remove (used by the delta refits).
    Returns list of dict rows."""
    opt = opt or Options()
    c = DA.covariates()
    M = DA.meta()
    X = DA.feature_values(features) if X is None else X
    Xv = X.to_numpy(float)
    pid = c["PatientID"].to_numpy().astype(str)
    ry = ~np.isnan(y_full)
    out = []
    # group features by row pattern
    patt = {}
    for j, f in enumerate(features):
        rows = ry & ~np.isnan(Xv[:, j])
        patt.setdefault(rows.tobytes(), (rows, []))[1].append(j)
    for rows, js in patt.values():
        n = int(rows.sum())
        Zd, merged = design_columns(rows, opt)
        Qd = orth(Zd)
        covC = {k: cov_columns(rows, g, k, opt) for k in DROPPABLE}
        acol = {k: assoc_cov_column(rows, g, k, opt) for k in DROPPABLE}
        cl = pid[rows]
        y = y_full[rows]
        # drops per feature
        groups = {}
        for j in js:
            f = features[j]
            d = Xv[rows, j]
            drops = structural_drops(g, f)
            for k in DROPPABLE:
                if k in drops or acol[k].shape[1] == 0:
                    continue
                if r2_assoc(d, acol[k]) > 0.49:
                    drops.add(k)
            extra = set(drop_extra.get(f, set())) if drop_extra else set()
            groups.setdefault((frozenset(drops), frozenset(extra)), []).append(j)
        for (drops, extra), jj in groups.items():
            cols = []
            if not ({"lineage", "library", "culture", "growth"} & extra):
                cols.append(Zd)
            else:
                cols.append(design_columns_minus(rows, opt, extra))
            kept = [k for k in DROPPABLE if k not in drops and k not in extra]
            cols += [covC[k] for k in kept]
            if "primary_disease" in extra:
                cols.append(cat_block(DA.primary_disease().to_numpy()[rows].astype(str), None)[0])
            Z = np.column_stack(cols)
            Q = orth(Z)
            rZ = Q.shape[1]
            Qlev = None
            if opt["lev_shared_only"]:
                Qlev = orth(np.column_stack([Zd] + [covC[k] for k in kept if k not in ("own_CN", "EXPGATE")]))
            D = Xv[np.ix_(rows, jj)]
            Dt = resid(Q, D)
            yt = resid(Q, y[:, None])[:, 0]
            var_d = (Dt ** 2).mean(0)
            # design-collinearity R^2
            Ddt = resid(Qd, D)
            Dc = D - D.mean(0)
            r2d = 1.0 - (Ddt ** 2).sum(0) / np.maximum((Dc ** 2).sum(0), 1e-300)
            binary = np.array([M[features[j]]["d_type"] == "binary" for j in jj])
            ncar = np.where(binary, (D == 1).sum(0), np.nan)
            keep = (var_d >= 1e-8) & (~binary | (ncar >= 5))
            beta = np.full(len(jj), np.nan); se = beta.copy(); nlev = beta.copy()
            if keep.any():
                b, s, nl = cr2_test(Q, yt, Dt[:, keep], cl, Qlev)
                beta[keep], se[keep], nlev[keep] = b, s, nl
            nu = nu_df(Dt, cl)
            dfres = n - rZ - 2
            ccl = np.array([len(set(cl[D[:, k] == 1])) if binary[k] else np.nan for k in range(len(jj))])
            df = np.where(binary, np.minimum(dfres, np.maximum(1.0, ccl - 1)), np.minimum(dfres, nu))
            t = beta / se
            p = 2 * stats.t.sf(np.abs(t), df)
            for k, j in enumerate(jj):
                out.append({"gene": g, "feature": features[j], "n_obs": n, "n_carriers": ncar[k],
                            "n_carrier_clusters": ccl[k], "dropped_cov": ";".join(sorted(drops)),
                            "r2_design": r2d[k], "beta": beta[k], "se": se[k], "t": t[k], "df": df[k],
                            "nu": nu[k], "p": p[k], "rank_Z": rZ, "n_lev1_excluded": nlev[k],
                            "tested": bool(keep[k]), "merged_levels": merged,
                            "sd_dtilde": float(np.sqrt(var_d[k]))})
    return out


BLOCK_CATS = {"lineage": ["lineage"], "library": ["library_group"], "culture": ["cellformat", "serumfree", "medium_family"]}


def design_columns_minus(rows, opt, extra):
    """Design block matrix with the named blocks removed (delta refits, 5.5)."""
    c = DA.covariates()
    cols = [np.ones((rows.sum(), 1))]
    for name, ref in DESIGN_CATS:
        blk = [b for b, cs in BLOCK_CATS.items() if name in cs][0]
        if blk in extra:
            continue
        cols.append(cat_block(c[name].to_numpy()[rows].astype(str), ref)[0])
    if "growth" not in extra:
        for name in DESIGN_NUMS:
            cols.append(c[name].to_numpy(float)[rows, None])
    if opt["na_flags_always"]:
        for name in ("CN_absent", "MSI_unknown", "METH_absent"):
            cols.append(c[name].to_numpy(float)[rows, None])
    return np.column_stack(cols)
