"""Confounder-model statistics (DESIGN v3.2 §5.3-§5.4, §6.1): vectorised FWL residualisation, the univariate
partial-regression test of many (g, d) at once with CR2 (HC2 for singleton patients) SEs and the carrier-based
small-sample df, the covariate-drop rule, BH, and the empirical-Bayes shrinkage power of §6.1.

Model for one (g, d) on rows r = (g non-NaN) ∩ (d non-NaN):

    y_g = [1, design blocks, droppable shared covariates \\ drops(d), gene covariates of g \\ drops(g, d)] γ + β d + ε

- design blocks (never dropped): lineage, library group, culture, growth (configs/pipeline_params.json:covariates);
  a level with < merge_min lines in the fit rows merges into the block's reference level.
- droppable shared covariates: MSI, DMG_BURDEN, TP53, screen days; dropped for d if d is (or its variant group
  contains) the covariate's feature, if R²(d ~ covariate dummies) > drop_assoc_r2_gt, or (TP53) if d derives from
  TP53 / MDM2 amplification.
- gene covariates: own CN of g (NaN -> 0 + CN_absent indicator) and EXPGATE of g; own CN dropped when d is
  CN-derived and d's gene or variant group contains g, or by the R² rule; EXPGATE by the R² rule.
- collinear_with_design: R²(d ~ design blocks) > collinear_with_design_r2_gt (test kept, flagged).
- SE: CR2 with PatientID clusters; singletons h_ii = h_Z,ii + d̃_i²/Σd̃²; multi-line clusters (I - H_cc)^(-1/2)
  (closed form for 2 x 2, batched eigh otherwise); rows with 1 - h_ii < tol are excluded from the meat.
- df = min(n - rank(Z) - 2, ν), ν = max(1, (Σ_c W_c)² / Σ_c W_c² - 1), W_c = Σ_{i in c} d̃_i².
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import linalg as sla
from scipy import stats as sps

from . import params as PR

# ------------------------------------------------------------------------------------------------
# Linear-algebra helpers
# ------------------------------------------------------------------------------------------------


def orth(A: np.ndarray, tol_rel: float | None = None) -> np.ndarray:
    """Orthonormal basis of span(A) (rank-revealing QR with column pivoting); n x 0 if A has no non-null column."""
    A = np.asarray(A, dtype=np.float64)
    if A.ndim == 1:
        A = A[:, None]
    if A.shape[1] == 0 or A.shape[0] == 0:
        return np.zeros((A.shape[0], 0))
    tol = PR.stats_params()["rank_tol"] if tol_rel is None else tol_rel
    Q, R, _ = sla.qr(A, mode="economic", pivoting=True, check_finite=False)
    d = np.abs(np.diag(R))
    if d.size == 0 or d[0] <= 0:
        return np.zeros((A.shape[0], 0))
    return Q[:, : int((d > d[0] * tol).sum())]


def resid(Q: np.ndarray, A: np.ndarray) -> np.ndarray:
    """A minus its projection on span(Q) (Q orthonormal)."""
    if Q.shape[1] == 0:
        return np.array(A, dtype=np.float64, copy=True)
    return A - Q @ (Q.T @ A)


def r2_on(Q: np.ndarray, X: np.ndarray, has_intercept: bool = False) -> np.ndarray:
    """R² of each column of X regressed on [1, span(Q)]; ``has_intercept`` = the constant is in span(Q)."""
    X = np.asarray(X, dtype=np.float64)
    if X.ndim == 1:
        X = X[:, None]
    Xc = X - X.mean(0)
    tot = (Xc ** 2).sum(0)
    if Q.shape[1] == 0:
        return np.zeros(X.shape[1])
    if has_intercept:
        Qc = Q
    else:
        one = np.ones((X.shape[0], 1)) / np.sqrt(X.shape[0])
        Qc = orth(resid(one, Q))
        if Qc.shape[1] == 0:
            return np.zeros(X.shape[1])
    ss = ((Qc.T @ Xc) ** 2).sum(0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(tot > 0, ss / tot, 0.0)


def bh(p) -> np.ndarray:
    """Benjamini-Hochberg q-values; NaN p stay NaN and do not count in m."""
    p = np.asarray(p, dtype=np.float64)
    q = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    m = pv.size
    if m == 0:
        return q
    o = np.argsort(pv, kind="stable")
    r = pv[o] * m / np.arange(1, m + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(m)
    out[o] = np.minimum(r, 1.0)
    q[ok] = out
    return q


def bh_cutoff(p, alpha: float) -> float:
    """Largest p with BH q <= alpha (the fixed threshold p* of DESIGN §7.1); alpha/m if nothing passes."""
    p = np.asarray(p, dtype=np.float64)
    p = p[~np.isnan(p)]
    if p.size == 0:
        return 0.0
    q = bh(p)
    ok = q <= alpha
    return float(p[ok].max()) if ok.any() else float(alpha / p.size)


# ------------------------------------------------------------------------------------------------
# Covariate design
# ------------------------------------------------------------------------------------------------


class CovDesign:
    """Dummy/numeric coding of the canonical covariate table with the per-fit level-merge rule."""

    def __init__(self, cov: pd.DataFrame, spec: dict | None = None):
        self.spec = spec or PR.pp()["covariates"]
        self.n = len(cov)
        self.index = cov.index
        self.items: dict[str, list] = {}
        self._cols: dict[str, tuple] = {}
        groups = dict(self.spec["design_blocks"])
        groups.update(self.spec["droppable_shared"])
        groups["primary_disease"] = [self.spec["primary_disease_block"]]
        for name, items in groups.items():
            keep = []
            for it in items:
                col = it["col"]
                if col not in cov.columns:
                    continue
                if it["kind"] == "cat":
                    codes, levels = pd.factorize(cov[col].astype(str), sort=True)
                    ref = it.get("reference")
                    ref_code = int(np.flatnonzero(levels == ref)[0]) if ref is not None and ref in set(levels) else -1
                    self._cols[col] = ("cat", codes.astype(np.int64), np.asarray(levels), ref_code)
                else:
                    v = pd.to_numeric(cov[col], errors="coerce").to_numpy(dtype=np.float64)
                    if np.isnan(v).any():
                        v = np.where(np.isnan(v), np.nanmedian(v) if np.isfinite(v).any() else 0.0, v)
                    self._cols[col] = ("num", v)
                keep.append(col)
            self.items[name] = keep
        self.design_names = [k for k in self.spec["design_blocks"] if self.items.get(k)]
        self.droppable_names = [k for k in self.spec["droppable_shared"] if self.items.get(k)]
        self.n_merged = 0

    def block(self, name: str, rows: np.ndarray, merge_min: int | None = None, apply_rows=None) -> np.ndarray:
        """Columns of covariate group ``name`` coded on ``rows`` (boolean or index array; levels with < merge_min
        lines there merge into the reference), evaluated on ``apply_rows`` (default ``rows``); no intercept."""
        mm = PR.stats_params()["merge_min"] if merge_min is None else merge_min
        ap = rows if apply_rows is None else apply_rows
        nap = int(np.sum(ap)) if np.asarray(ap).dtype == bool else len(ap)
        parts = []
        for col in self.items.get(name, []):
            info = self._cols[col]
            if info[0] == "num":
                v = info[1][rows]
                if v.size and np.ptp(v) > 0:
                    parts.append((info[1][ap] - v.mean())[:, None])
                continue
            _, codes, levels, ref_code = info
            c = codes[rows]
            cnt = np.bincount(c, minlength=len(levels))
            ref = ref_code if ref_code >= 0 and cnt[ref_code] > 0 else int(np.argmax(cnt))
            keep = np.flatnonzero((cnt >= mm) & (np.arange(len(levels)) != ref))
            self.n_merged += int(((cnt > 0) & (cnt < mm) & (np.arange(len(levels)) != ref)).sum())
            if keep.size:
                parts.append((codes[ap][:, None] == keep[None, :]).astype(np.float64))
        if not parts:
            return np.zeros((nap, 0))
        return np.hstack(parts)

    def Z(self, rows: np.ndarray, names, apply_rows=None) -> np.ndarray:
        ap = rows if apply_rows is None else apply_rows
        nr = int(np.sum(ap)) if np.asarray(ap).dtype == bool else len(ap)
        return np.hstack([np.ones((nr, 1))] + [self.block(k, rows, apply_rows=ap) for k in names])


# ------------------------------------------------------------------------------------------------
# Feature drop information (identity / derivation part of the §5.3 drop rule)
# ------------------------------------------------------------------------------------------------


@dataclass
class DropInfo:
    """Per-feature static drop facts. ``ident[j]`` = droppable shared names dropped by identity/derivation;
    ``cn_genes[j]`` = genes whose own-CN covariate is dropped for d_j (CN-derived d: own gene(s) + group members)."""
    ident: list = field(default_factory=list)
    cn_genes: list = field(default_factory=list)

    @staticmethod
    def build(fmeta: pd.DataFrame, tp53_genes=(), mdm2_genes=()) -> "DropInfo":
        spec = PR.pp()["covariates"]
        cov_feat = {k: set(v) for k, v in spec["covariate_is_feature"].items()}
        cn_fam = set(spec["cn_derived_families"])
        amp_fam = set(spec["tp53_derivation"]["amp_families"])
        tp53_genes, mdm2_genes = set(map(str, tp53_genes)), set(map(str, mdm2_genes))
        vg = fmeta["variant_group"].where(fmeta["variant_group"].notna(), None)
        members_f = {}
        for f, g in vg.items():
            if g:
                members_f.setdefault(g, set()).add(f)
        members_g = {}
        for f, g in vg.items():
            if g:
                members_g.setdefault(g, set()).update(fmeta.at[f, "genes"])
        ident, cn_genes = [], []
        for f, r in fmeta.iterrows():
            grp = vg.get(f)
            same = {f} | (members_f.get(grp, set()) if grp else set())
            s = {k for k, fs in cov_feat.items() if fs & same}
            genes = set(r["genes"])
            fam = r["family"]
            if genes & tp53_genes or (fam in amp_fam and genes & mdm2_genes):
                s.add(spec["tp53_derivation"]["covariate"])
            ident.append(frozenset(s))
            if fam in cn_fam:
                cn_genes.append(frozenset(genes | (members_g.get(grp, set()) if grp else set())))
            else:
                cn_genes.append(frozenset())
        return DropInfo(ident=ident, cn_genes=cn_genes)


# ------------------------------------------------------------------------------------------------
# Clusters and the CR2 core
# ------------------------------------------------------------------------------------------------


@dataclass
class Clusters:
    single: np.ndarray            # positions (within the fit rows) of singleton patients
    multi: list                   # list of position arrays (size >= 2)

    @staticmethod
    def of(patients_rows: np.ndarray) -> "Clusters":
        codes, uniq = pd.factorize(pd.Series(patients_rows).astype(str))
        cnt = np.bincount(codes)
        single = np.flatnonzero(cnt[codes] == 1)
        order = np.argsort(codes, kind="stable")
        multi = []
        for c in np.flatnonzero(cnt >= 2):
            multi.append(order[np.searchsorted(codes[order], c):np.searchsorted(codes[order], c, side="right")])
        return Clusters(single=single, multi=multi)


def _inv_sqrt_batch(M: np.ndarray, floor: float):
    """Symmetric inverse square root of a batch (F, s, s) of PSD matrices; eigenvalues <= floor are treated as
    zero (pseudo-inverse) and counted. Closed form for s == 2."""
    F, s, _ = M.shape
    if s == 2 and F:
        a, b, c = M[:, 0, 0], M[:, 0, 1], M[:, 1, 1]
        det = a * c - b * b
        good = det > floor
        out = np.empty_like(M)
        nbad = np.zeros(F, dtype=np.int64)
        if good.any():
            sd = np.sqrt(det[good])
            ag, bg, cg = a[good], b[good], c[good]
            t = np.sqrt(ag + cg + 2 * sd)
            den = (ag + sd) * (cg + sd) - bg * bg
            k = t / den
            out[good, 0, 0] = k * (cg + sd)
            out[good, 1, 1] = k * (ag + sd)
            out[good, 0, 1] = out[good, 1, 0] = -k * bg
        if (~good).any():
            o2, n2 = _inv_sqrt_eigh(M[~good], floor)
            out[~good], nbad[~good] = o2, n2
        return out, nbad
    return _inv_sqrt_eigh(M, floor)


def _inv_sqrt_eigh(M: np.ndarray, floor: float):
    w, V = np.linalg.eigh(M)
    ok = w > floor
    winv = np.where(ok, 1.0 / np.sqrt(np.where(ok, w, 1.0)), 0.0)
    return np.einsum("fij,fj,fkj->fik", V, winv, V), (~ok).sum(1)


def carrier_clusters(Xraw: np.ndarray, binary, cl: Clusters) -> np.ndarray:
    """Distinct patient clusters among the carriers (x == 1) of each column of Xraw (fit rows, same order as ``cl``);
    NaN for non-binary columns."""
    Xraw = np.asarray(Xraw, dtype=np.float64)
    if Xraw.ndim == 1:
        Xraw = Xraw[:, None]
    C = (Xraw == 1).astype(np.float64)
    out = C.sum(0)
    for idx in cl.multi:
        out = out - np.maximum(C[idx].sum(0) - 1.0, 0.0)
    return np.where(np.broadcast_to(np.asarray(binary, dtype=bool), out.shape), out, np.nan)


def shared_df(Xt: np.ndarray, cl: Clusters, n: int, rankZ: int, car_clusters=None, sp: dict | None = None) -> dict:
    """The df rule shared by every test and its power (DESIGN 5.3, 6.1). cap = n - rank(Z) - df_subtract;
    nu = max(1, (sum_c W_c)^2 / sum_c W_c^2 - 1) with W_c = sum over lines i of cluster c of dtilde_i^2;
    df_nu = min(cap, nu). Binary d with ``car_clusters`` given and df_rule.binary = fallback_carrier_clusters:
    df = min(cap, carrier clusters - 1) (the pre-declared calibration fallback); otherwise df = df_nu.
    Returns df, nu, df_nu, df_rule (per column), cap."""
    sp = sp or PR.stats_params()
    Xt = np.asarray(Xt, dtype=np.float64)
    if Xt.ndim == 1:
        Xt = Xt[:, None]
    X2 = Xt ** 2
    sxx = X2.sum(0)
    W2 = (X2[cl.single] ** 2).sum(0)
    for idx in cl.multi:
        W2 = W2 + X2[idx].sum(0) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        nu = np.maximum(1.0, sxx ** 2 / np.where(W2 > 0, W2, np.nan) - 1.0)
    cap = float(n - rankZ - sp["df_subtract"])
    df_nu = np.minimum(cap, nu)
    rule = np.full(sxx.shape, "nu", dtype=object)
    df = df_nu.copy()
    if car_clusters is not None and sp["df_rule_binary"] == "fallback_carrier_clusters":
        cc = np.broadcast_to(np.asarray(car_clusters, dtype=np.float64), sxx.shape)
        fb = np.isfinite(cc)
        df = np.where(fb, np.minimum(cap, cc - 1.0), df_nu)
        rule[fb] = "fallback_carrier_clusters"
    with np.errstate(invalid="ignore"):
        df = np.where(df > 0, df, np.nan)
    return {"df": df, "nu": nu, "df_nu": df_nu, "df_rule": rule, "cap": cap}


def cr2_core(Yt: np.ndarray, Xt: np.ndarray, hZ: np.ndarray, Qz: np.ndarray, cl: Clusters, rankZ: int,
             sp: dict | None = None, car_clusters=None) -> dict:
    """Partial-regression CR2 test of every column of Xt (n x F, residualised on Z) against every column of Yt
    (n x T, residualised on Z). ``Qz`` (n x k) orthonormal basis of Z (for the cluster blocks of H), ``hZ`` its row
    norms. ``car_clusters`` (F; NaN for continuous columns) switches binary columns to the fallback df (shared_df);
    without it every column uses the nu rule. Returns arrays (F, T): beta, se, t, df, nu, p, df_nu and p_nu (the
    nu-rule sensitivity); plus sxx, n_lev1, df_rule (F)."""
    sp = sp or PR.stats_params()
    tol, floor = sp["lev1_tol"], sp["eig_floor"]
    n, F = Xt.shape
    T = Yt.shape[1]
    sxx = (Xt ** 2).sum(0)
    safe = np.where(sxx > 0, sxx, 1.0)
    B = (Xt.T @ Yt) / safe[:, None]                          # F x T
    X2 = Xt ** 2
    dd = shared_df(Xt, cl, n, rankZ, car_clusters, sp)
    df, nu, df_nu = dd["df"], dd["nu"], dd["df_nu"]
    # leverages of singleton rows
    hs = hZ[cl.single][:, None] + X2[cl.single] / safe[None, :]
    oneh = 1.0 - hs
    okh = oneh >= tol
    inv1h = np.where(okh, 1.0 / np.where(okh, oneh, 1.0), 0.0)
    nlev = (~okh).sum(0)
    # multi-line cluster factors (target independent)
    Afac = []
    for idx in cl.multi:
        s = len(idx)
        Hz = Qz[idx] @ Qz[idx].T
        xc = Xt[idx].T                                       # F x s
        M = np.eye(s)[None] - Hz[None] - xc[:, :, None] * xc[:, None, :] / safe[:, None, None]
        A, nb = _inv_sqrt_batch(M, floor)
        nlev = nlev + nb
        Afac.append((idx, xc, A))
    se = np.empty((F, T))
    Xs = Xt[cl.single]
    for k in range(T):
        e_s = Yt[cl.single, k][:, None] - Xs * B[:, k][None, :]
        meat = (X2[cl.single] * e_s ** 2 * inv1h).sum(0)
        for idx, xc, A in Afac:
            ec = Yt[idx, k][None, :] - xc * B[:, k][:, None]          # F x s
            u = np.einsum("fij,fj->fi", A, ec)
            meat = meat + (xc * u).sum(1) ** 2
        se[:, k] = np.sqrt(meat) / safe
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(se > 0, B / se, np.nan)
    p = 2 * sps.t.sf(np.abs(t), df[:, None])
    p_nu = 2 * sps.t.sf(np.abs(t), df_nu[:, None])
    bad = sxx <= 0
    B[bad], se[bad], t[bad], p[bad], p_nu[bad] = np.nan, np.nan, np.nan, np.nan, np.nan

    def full(v):
        return np.broadcast_to(v[:, None], B.shape).copy()

    return {"beta": B, "se": se, "t": t, "df": full(df), "nu": full(nu), "p": p, "df_nu": full(df_nu), "p_nu": p_nu,
            "sxx": sxx, "n_lev1": nlev, "df_rule": dd["df_rule"]}


# ------------------------------------------------------------------------------------------------
# Single-pair test (reference implementation; reused by D-stage audits and replicate.py)
# ------------------------------------------------------------------------------------------------


def pair_test(y: np.ndarray, d: np.ndarray, Z: np.ndarray, patients: np.ndarray, one_sided: str | None = None,
              sp: dict | None = None, binary: bool | None = None) -> dict:
    """CR2 test of d in y ~ Z + d on the rows where y and d are finite (Z must include the intercept), df by
    ``shared_df``. ``binary`` True applies the binary df rule (fallback: carrier clusters - 1), False the nu rule
    (continuous d); None detects binary d from its finite values (all in {0, 1}). ``one_sided`` in {None, 'neg',
    'pos'}: p for H1 beta < 0 / beta > 0. Returns beta, se, t, df, nu, p, df_nu, p_nu, df_rule, n, n_carriers,
    n_carrier_clusters, n_lev1_excluded, se_type."""
    y = np.asarray(y, np.float64)
    d = np.asarray(d, np.float64)
    ok = np.isfinite(y) & np.isfinite(d)
    Zr = np.asarray(Z, np.float64)[ok]
    Q = orth(Zr)
    yt, dt = resid(Q, y[ok][:, None]), resid(Q, d[ok][:, None])
    cl = Clusters.of(np.asarray(patients)[ok])
    if binary is None:
        binary = bool(ok.any() and np.isin(d[ok], (0.0, 1.0)).all())
    cc = carrier_clusters(d[ok], True, cl) if binary else None
    r = cr2_core(yt, dt, (Q ** 2).sum(1), Q, cl, Q.shape[1], sp, car_clusters=cc)
    out = {k: float(r[k][0, 0]) for k in ("beta", "se", "t", "df", "nu", "p", "df_nu", "p_nu")}
    out["df_rule"] = str(r["df_rule"][0])
    out["n_carrier_clusters"] = float(cc[0]) if cc is not None else np.nan
    if one_sided in ("neg", "pos") and np.isfinite(out["t"]):
        tt = out["t"] if one_sided == "pos" else -out["t"]
        out["p"] = float(sps.t.sf(tt, out["df"]))
    e = yt[:, 0] - (out["beta"] if np.isfinite(out["beta"]) else 0.0) * dt[:, 0]
    dfres = max(int(ok.sum()) - Q.shape[1] - 1, 1)
    out.update(n=int(ok.sum()), n_carriers=int((d[ok] == 1).sum()), n_lev1_excluded=int(r["n_lev1"][0]),
               se_type=PR.pp()["se_type"], rankZ=int(Q.shape[1]), sigma_resid=float(np.sqrt((e ** 2).sum() / dfres)),
               sd_dtilde=float(np.sqrt((dt[:, 0] ** 2).mean())) if ok.any() else np.nan)
    return out


# ------------------------------------------------------------------------------------------------
# Vectorised S1 scan
# ------------------------------------------------------------------------------------------------


@dataclass
class ScanInputs:
    """Everything the scan needs (rows = D lines, n)."""
    X: np.ndarray                        # n x F raw feature values (NaN = not assayed)
    fmeta: pd.DataFrame                  # F rows (same order): family, d_type, scope, genes, variant_group
    Dg: dict                             # gene -> int array of feature positions (𝒟_g minus self-null)
    rows_ok: np.ndarray                  # n x G bool (effect non-NaN), column order = genes
    genes: list
    cov: CovDesign
    own_cn: np.ndarray | None            # n x G log2 own CN (NaN = not profiled) or None
    expgate: np.ndarray | None           # n x G 0/1 or None
    patients: np.ndarray                 # n
    drop: DropInfo


def _pattern_groups(mask: np.ndarray) -> dict:
    """columns grouped by identical boolean column pattern -> {key: (rows_bool, col_positions)}."""
    out = {}
    for j in range(mask.shape[1]):
        k = mask[:, j].tobytes()
        if k not in out:
            out[k] = (mask[:, j].copy(), [])
        out[k][1].append(j)
    return {k: (v[0], np.asarray(v[1], dtype=np.int64)) for k, v in out.items()}


def _counts(Xr: np.ndarray, binary: np.ndarray):
    car = np.where(binary, (Xr == 1).sum(0), np.nan)
    sd = Xr.std(0) if Xr.shape[0] else np.zeros(Xr.shape[1])
    return car, sd


def in_universe(n_obs, n_car, sd, binary, sp=None) -> np.ndarray:
    sp = sp or PR.stats_params()
    return np.where(binary, n_car >= sp["min_carriers"], (n_obs >= sp["cont_min_nonnan"]) & (sd > 0))


def s1_scan(inp: ScanInputs, targets: dict, compute: bool = True, genes_subset=None) -> pd.DataFrame:
    """Univariate confounder-model test of every (g, d) in U for every target matrix in ``targets``
    ({arm: n x G array, columns in inp.genes order}). ``compute=False`` returns only the universe rows (counts,
    drops, collinearity). Output: one row per (arm, gene, feature) [arm omitted when compute=False]."""
    sp = PR.stats_params()
    X, fm, cov = inp.X, inp.fmeta, inp.cov
    binary = (fm["d_type"].to_numpy() == "binary")
    shared = (fm["scope"].to_numpy() == "shared")
    fid = fm.index.to_numpy()
    gpos = {g: i for i, g in enumerate(inp.genes)}
    genes = list(inp.genes if genes_subset is None else genes_subset)
    arms = list(targets) if compute else []
    feat_nan = np.isnan(X)
    fgroups = _pattern_groups(~feat_nan)
    fkeys = list(fgroups)
    fgroup_of = np.empty(X.shape[1], dtype=np.int64)
    for i, k in enumerate(fkeys):
        fgroup_of[fgroups[k][1]] = i
    gmask = inp.rows_ok[:, [gpos[g] for g in genes]]
    ggroups = _pattern_groups(gmask)
    dnames, xnames = cov.design_names, cov.droppable_names
    out: dict = {}
    for gkey, (grows, gcols) in ggroups.items():
        gl = [genes[j] for j in gcols]
        need = np.unique(np.concatenate([inp.Dg[g] for g in gl])) if gl else np.zeros(0, int)
        memb = np.zeros((len(gl), X.shape[1]), dtype=bool)
        for jg, g in enumerate(gl):
            memb[jg, inp.Dg[g]] = True
        need_grp = fgroup_of[need]
        order = np.argsort(need_grp, kind="stable")
        bounds = np.flatnonzero(np.diff(need_grp[order])) + 1
        for chunk in np.split(order, bounds):
            if chunk.size == 0:
                continue
            fcols = need[chunk]
            frows = fgroups[fkeys[need_grp[chunk[0]]]][0]
            rows = grows & frows
            nr = int(rows.sum())
            if nr < 2:
                continue
            ri = np.flatnonzero(rows)
            Xr = X[np.ix_(ri, fcols)].astype(np.float64)
            car, sd = _counts(Xr, binary[fcols])
            inU = in_universe(nr, car, sd, binary[fcols], sp)
            if not inU.any():
                continue
            fcols, Xr, car, sd = fcols[inU], Xr[:, inU], car[inU], sd[inU]
            # collinear with design; R² drop rule for droppable shared covariates
            Qd = orth(cov.Z(rows, dnames))
            r2d = r2_on(Qd, Xr, has_intercept=True)
            drops = [set(inp.drop.ident[j]) for j in fcols]
            for nm in xnames:
                Bk = cov.block(nm, rows)
                if Bk.shape[1] == 0:
                    continue
                r2 = r2_on(orth(Bk), Xr)
                for i in np.flatnonzero(r2 > sp["drop_r2"]):
                    drops[i].add(nm)
            dkey = np.array(["|".join(sorted(s)) for s in drops], dtype=object)
            Qs, Xts = {}, {}
            for dk in np.unique(dkey):
                S = set(dk.split("|")) if dk else set()
                Q = orth(cov.Z(rows, dnames + [nm for nm in xnames if nm not in S]))
                sel = np.flatnonzero(dkey == dk)
                Qs[dk] = Q
                Xts[dk] = (sel, resid(Q, Xr[:, sel]))
            gidx = np.array([gpos[g] for g in gl], dtype=np.int64)
            Yr = {a: targets[a][np.ix_(ri, gidx)].astype(np.float64) for a in arms}
            Ytil = {dk: {a: resid(Qs[dk], Yr[a]) for a in arms} for dk in Qs}
            prow = inp.patients[rows]
            cl = Clusters.of(prow) if compute else None
            carcl = carrier_clusters(Xr, binary[fcols], cl) if compute else None
            nNA_g = int(grows.sum()) - nr
            Mg = memb[:, fcols]
            for jg in np.flatnonzero(Mg.any(1)):
                g = gl[jg]
                inDg = Mg[jg]
                gc = _gene_covs(inp, gpos[g], rows)
                for dk, (sel, Xt_all) in Xts.items():
                    m = inDg[sel]
                    if not m.any():
                        continue
                    loc = sel[m]                       # positions within fcols
                    Xt = Xt_all[:, m]
                    codes, gdrop = _gene_cov_codes(gc, Xr[:, loc], fcols[loc], g, inp.drop, sp)
                    for code in np.unique(codes):
                        cm = codes == code
                        cols_c = [c for c, keep in zip(gc["cols"], _code_keep(code, gc)) if keep]
                        Q = Qs[dk]
                        if cols_c:
                            Qc = orth(resid(Q, np.column_stack(cols_c)))
                        else:
                            Qc = np.zeros((nr, 0))
                        Xc = resid(Qc, Xt[:, cm])
                        rankZ = Q.shape[1] + Qc.shape[1]
                        Qfull = np.hstack([Q, Qc])
                        lc = loc[cm]
                        nf = lc.size
                        dco = ["|".join(sorted(set(dk.split("|") if dk else []) | gdrop[i])) for i in np.flatnonzero(cm)]
                        base = {"gene": np.full(nf, g, dtype=object), "feature": fid[fcols[lc]],
                                "n_obs": np.full(nf, nr), "n_carriers": car[lc], "n_NA": np.full(nf, nNA_g),
                                "collinear_with_design": r2d[lc] > sp["collinear_r2"], "r2_design": r2d[lc],
                                "dropped_cov": np.asarray(dco, dtype=object)}
                        if not compute:
                            _acc(out, base)
                            continue
                        Yt = np.column_stack([resid(Qc, Ytil[dk][a][:, [jg]]) for a in arms])
                        r = cr2_core(Yt, Xc, (Qfull ** 2).sum(1), Qfull, cl, rankZ, sp, car_clusters=carcl[lc])
                        lowvar = r["sxx"] / nr < sp["min_resid_var"]
                        for k, a in enumerate(arms):
                            rec = dict(base)
                            rec["arm"] = np.full(nf, a, dtype=object)
                            for key in ("beta", "se", "t", "df", "nu", "p", "df_nu", "p_nu"):
                                v = r[key][:, k].copy()
                                v[lowvar] = np.nan
                                rec[key] = v
                            rec["n_lev1_excluded"] = r["n_lev1"]
                            rec["sd_dtilde"] = np.sqrt(r["sxx"] / nr)
                            rec["df_rule"] = r["df_rule"]
                            rec["n_carrier_clusters"] = carcl[lc]
                            _acc(out, rec)
    cols = (["arm"] if compute else []) + ["gene", "feature", "n_obs", "n_carriers", "n_NA",
                                           "collinear_with_design", "r2_design", "dropped_cov"]
    if not out:
        return pd.DataFrame(columns=cols)
    res = pd.DataFrame({k: np.concatenate(v) for k, v in out.items()})
    res = res[cols + [c for c in res.columns if c not in cols]]
    res["se_type"] = PR.pp()["se_type"]
    return res


def _acc(store: dict, rec: dict) -> None:
    for k, v in rec.items():
        store.setdefault(k, []).append(np.asarray(v))


def _gene_covs(inp: ScanInputs, gj: int, rows: np.ndarray) -> dict:
    """Gene covariate columns of g on rows: own CN (NaN -> 0) [+ CN_absent indicator], EXPGATE; constant columns
    omitted. 'kind' tags which droppable each column belongs to."""
    cols, kind = [], []
    if inp.own_cn is not None:
        c = inp.own_cn[rows, gj].astype(np.float64)
        na = np.isnan(c)
        c0 = np.where(na, 0.0, c)
        if np.ptp(c0) > 0:
            cols.append(c0)
            kind.append("own_CN")
            if na.any() and not na.all():
                cols.append(na.astype(np.float64))
                kind.append("own_CN")
    if inp.expgate is not None:
        e = inp.expgate[rows, gj].astype(np.float64)
        e = np.where(np.isnan(e), 0.0, e)
        if np.ptp(e) > 0:
            cols.append(e)
            kind.append("EXPGATE")
    return {"cols": cols, "kind": kind}


def _code_keep(code: int, gc: dict) -> list:
    drop_cn, drop_eg = bool(code & 2), bool(code & 1)
    return [not ((k == "own_CN" and drop_cn) or (k == "EXPGATE" and drop_eg)) for k in gc["kind"]]


def _gene_cov_codes(gc: dict, Xr: np.ndarray, fglob: np.ndarray, g: str, drop: DropInfo, sp: dict):
    """Per feature: code = 2*drop(own CN) + drop(EXPGATE) and the set of dropped gene-covariate names."""
    F = Xr.shape[1]
    dcn = np.array([g in drop.cn_genes[j] for j in fglob], dtype=bool)
    deg = np.zeros(F, dtype=bool)
    kinds = np.asarray(gc["kind"], dtype=object)
    if (kinds == "own_CN").any():
        Q = orth(np.column_stack([c for c, k in zip(gc["cols"], kinds) if k == "own_CN"]))
        dcn |= r2_on(Q, Xr) > sp["drop_r2"]
    else:
        dcn[:] = False
    if (kinds == "EXPGATE").any():
        Q = orth(np.column_stack([c for c, k in zip(gc["cols"], kinds) if k == "EXPGATE"]))
        deg |= r2_on(Q, Xr) > sp["drop_r2"]
    codes = dcn.astype(int) * 2 + deg.astype(int)
    names = [set((["own_CN"] if a else []) + (["EXPGATE"] if b else [])) for a, b in zip(dcn, deg)]
    return codes, names


# ------------------------------------------------------------------------------------------------
# BH families, ranking
# ------------------------------------------------------------------------------------------------


def add_bh(df: pd.DataFrame, own_expdose: np.ndarray, alpha: float | None = None) -> pd.DataFrame:
    """test_family F1 / F2 (own EXPDOSE), BH q within (arm, family) over all U rows with a p-value; selected =
    F1 and q < alpha and not collinear_with_design (the S1 list); selected_F2 = the same rule within F2 (reported
    separately; DESIGN §5.6, amendment A6-iii)."""
    alpha = PR.stats_params()["bh_q"] if alpha is None else alpha
    df = df.copy()
    df["test_family"] = np.where(own_expdose, "F2", "F1")
    df["q"] = np.nan
    for (a, tf), idx in df.groupby(["arm", "test_family"]).groups.items():
        df.loc[idx, "q"] = bh(df.loc[idx, "p"].to_numpy())
    ok = (df["q"] < alpha) & ~df["collinear_with_design"].astype(bool)
    df["selected"] = ok & (df["test_family"] == "F1")
    df["selected_F2"] = ok & (df["test_family"] == "F2")
    return df


# ------------------------------------------------------------------------------------------------
# Empirical-Bayes shrinkage and power (DESIGN §6.1)
# ------------------------------------------------------------------------------------------------


def eb_scale_mixture(b: np.ndarray, s: np.ndarray, rp: dict | None = None) -> dict:
    """Normal scale-mixture prior g = Σ π_k N(0, σ_k²) fitted by EM on b ~ N(θ, s²); returns posterior means,
    the grid and the weights. Grid: n_sigma log-spaced σ from frac*min(s) to mult*sqrt(max(b² - s²))."""
    rp = rp or PR.replication_params()
    b, s = np.asarray(b, float), np.asarray(s, float)
    ok = np.isfinite(b) & np.isfinite(s) & (s > 0)
    bo, so = b[ok], s[ok]
    lo = rp["eb_sigma_min_frac"] * so.min()
    hi = rp["eb_sigma_max_mult"] * np.sqrt(max(np.max(bo ** 2 - so ** 2), (2 * lo) ** 2))
    sig = np.geomspace(lo, max(hi, lo * 2), rp["eb_n_sigma"])
    var = sig[None, :] ** 2 + so[:, None] ** 2
    L = np.exp(-0.5 * bo[:, None] ** 2 / var) / np.sqrt(var)  # literal-ok: normal log-density factor
    pi = np.full(sig.size, 1.0 / sig.size)
    for _ in range(rp["eb_max_iter"]):
        w = L * pi
        w /= np.maximum(w.sum(1, keepdims=True), PR.stats_params()["tiny"])
        new = w.mean(0)
        done = np.max(np.abs(new - pi)) < rp["eb_tol"]
        pi = new
        if done:
            break
    w = L * pi
    w /= np.maximum(w.sum(1, keepdims=True), PR.stats_params()["tiny"])
    shrink = (w * (sig[None, :] ** 2 / var)).sum(1)
    post = np.full(b.shape, np.nan)
    post[ok] = shrink * bo
    return {"posterior_mean": post, "sigma": sig, "pi": pi}


def shrunk_beta(beta: np.ndarray, se: np.ndarray, sd_dtilde: np.ndarray) -> np.ndarray:
    """§6.1: standardise (β* = β·SD_D(d̃)), shrink with one EB fit over all rows, back-transform."""
    sd = np.asarray(sd_dtilde, float)
    r = eb_scale_mixture(np.asarray(beta) * sd, np.asarray(se) * sd)
    with np.errstate(divide="ignore", invalid="ignore"):
        return r["posterior_mean"] / sd


def se_r_pred(sigma_resid_D: np.ndarray, binary: np.ndarray, n1_R=None, n0_R=None, sd_R=None, n_R=None):
    """SE_R_pred of §6.1 (binary: σ·sqrt(1/n1 + 1/n0); continuous: σ / (SD_R(d̃_R)·sqrt(n_R)))."""
    sig = np.asarray(sigma_resid_D, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        sb = sig * np.sqrt(1.0 / np.asarray(n1_R, float) + 1.0 / np.asarray(n0_R, float)) if n1_R is not None else np.nan
        sc = sig / (np.asarray(sd_R, float) * np.sqrt(np.asarray(n_R, float))) if sd_R is not None else np.nan
    return np.where(np.asarray(binary, bool), sb, sc)


def power_t(beta_shrunk, se_pred, nu, alpha: float | None = None) -> np.ndarray:
    """1 - F_t(t_{1-α, ν}; ncp = |β_shrunk| / SE_R_pred) (noncentral t; same df rule as the test)."""
    alpha = PR.replication_params()["alpha"] if alpha is None else alpha
    nu = np.asarray(nu, float)
    ncp = np.abs(np.asarray(beta_shrunk, float)) / np.asarray(se_pred, float)
    crit = sps.t.ppf(1 - alpha, nu)
    return sps.nct.sf(crit, nu, ncp)


def nu_from(dtilde: np.ndarray, clusters) -> float:
    """ν of §5.3 from residualised d̃ and cluster labels (for R-side power before any R effect is read)."""
    d2 = np.asarray(dtilde, float) ** 2
    W = pd.Series(d2).groupby(np.asarray(clusters)).sum().to_numpy()
    return float(max(1.0, W.sum() ** 2 / (W ** 2).sum() - 1.0))
