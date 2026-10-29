"""Task 3: independent two-stage fixed-lambda stability elastic net (DESIGN 4.1) for EN-obs.

Reading (fixed before comparing with the engine):
* rows: D lines with a non-NaN effect for g (identical for both EN arms); target = observed effect, not scaled.
* features: my D_g restricted to pairs meeting the 3 counts in g's rows (binary >= 5 carriers among non-NaN
  lines, continuous >= 20 non-NaN lines and SD > 0); standardised once on all D lines (mean, SD with ddof 0,
  non-NaN values), then NA -> 0 (i.e. NA -> standardised 0) - the CN_absent / MSI_unknown / METH_absent
  indicators are in the unpenalised block.
* FWL block per fit (estimated on that fit's rows): intercept, lineage (ref _rare), library group, CellFormat,
  SerumFree, medium family, ScreenDoublingTime, DT_unknown, CN_absent, MSI_unknown, METH_absent; < 2-line levels
  merged into the reference (5.3 leverage rule, referenced by 4.1.1). Held-out CV rows are residualised with
  the training-fold coefficients.
* column drop per fit: residual SD < 1e-6 or (binary) < 2 carriers in the fit's rows.
* SIS per stage and per fit: top K = 300 by |corr(residualised column, stage target)| (ties: feature order),
  union the forced terms (gene-specific relations except self-null). Stage G target = residualised y;
  stage X target = the same fit's stage-G residual.
* lambda: sklearn ElasticNet objective (1/2n)||r||^2 + a*l1*|b|_1 + a*(1-l1)/2*|b|^2, l1_ratio 0.9, no
  intercept, tol 1e-3, max_iter 10000, cyclic. 25-alpha log grid from alpha_max = max|X't|/(n*l1_ratio) of
  the full-D screened set down to eps*alpha_max. 5-fold CV over P1 cv_fold (SIS inside each training fold),
  1-SE rule: largest alpha with mean fold MSE <= min mean + SD(fold MSEs, ddof 1)/sqrt(5) at the minimum.
  Stage G first; stage X CV uses each training fold's stage-G fit at lambda_G.
* stability: the 50 subsamples of configs/subsamples.json intersected with g's rows; fixed alphas; frequency =
  share of the 50 fits with beta != 0; sign consistency = majority share among non-zero fits; selected =
  frequency >= 0.6 and sign consistency >= 0.9; mean beta over all 50 (zeros included).
"""
from __future__ import annotations

import os
import warnings

import numpy as np

from . import data as DA
from .confounder import cat_block, orth

L1, K_SIS, NLAM, EPS, TOL, MAXIT = 0.9, 300, 25, 1e-2, 1e-3, 10000
SD_MIN, CAR_MIN = 1e-6, 2
FREQ_MIN, SIGN_MIN = 0.6, 0.9
CATS = [("lineage", "_rare"), ("library_group", None), ("cellformat", None), ("serumfree", None),
        ("medium_family", None)]
NUMS = ["doubling_time", "DT_unknown", "CN_absent", "MSI_unknown", "METH_absent"]


def load_params():
    s = DA.design_params()["surrogate_s4"]
    pp = DA.pipeline_params()["en"]
    assert s["l1_ratio"] == L1 and s["sis_K_per_stage"] == K_SIS and s["n_lambda"] == NLAM
    assert s["eps"] == EPS and s["tol"] == TOL and s["selection_freq_min"] == FREQ_MIN
    assert s["sign_consistency_min"] == SIGN_MIN and pp["max_iter"] == MAXIT
    assert s["drop_column_in_fit"] == {"resid_sd_lt": SD_MIN, "carriers_lt": CAR_MIN}


class ZCoder:
    """Unpenalised block coded on a fit's rows; applies the same coding to other rows."""

    def __init__(self, fit_rows: np.ndarray):
        c = DA.covariates()
        self.spec = []
        for name, ref in CATS:
            v = c[name].to_numpy().astype(str)[fit_rows]
            from collections import Counter
            cnt = Counter(v.tolist())
            r = ref if ref is not None else sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
            levels = sorted(l for l, n in cnt.items() if n >= 2 and l != r)
            self.spec.append((name, levels))

    def __call__(self, rows: np.ndarray) -> np.ndarray:
        c = DA.covariates()
        cols = [np.ones((rows.sum(), 1))]
        for name, levels in self.spec:
            v = c[name].to_numpy().astype(str)[rows]
            if levels:
                cols.append(np.column_stack([(v == L).astype(float) for L in levels]))
        for name in NUMS:
            cols.append(c[name].to_numpy(float)[rows, None])
        return np.column_stack(cols)


class Residualiser:
    def __init__(self, fit_rows: np.ndarray):
        self.zc = ZCoder(fit_rows)
        Z = self.zc(fit_rows)
        self.Z = Z
        self.Q = orth(Z)
        self.pinv = np.linalg.pinv(Z)

    def fit(self, X: np.ndarray) -> np.ndarray:
        return X - self.Q @ (self.Q.T @ X)

    def apply(self, rows: np.ndarray, Xfit: np.ndarray, Xother: np.ndarray) -> np.ndarray:
        G = self.pinv @ Xfit
        return Xother - self.zc(rows) @ G


def enet(X, y, alphas):
    from sklearn.linear_model import enet_path
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, coefs, _ = enet_path(X, y, l1_ratio=L1, alphas=np.asarray(alphas, float), tol=TOL, max_iter=MAXIT,
                                precompute=False)
    return coefs  # p x n_alpha


def enet_one(X, y, alpha):
    from sklearn.linear_model import ElasticNet
    if X.shape[1] == 0:
        return np.zeros(0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = ElasticNet(alpha=alpha, l1_ratio=L1, fit_intercept=False, tol=TOL, max_iter=MAXIT,
                       selection="cyclic", precompute=False).fit(X, y)
    return m.coef_


class Gene:
    def __init__(self, g: str, y_full: np.ndarray):
        self.g = g
        M = DA.meta()
        self.rows = ~np.isnan(y_full)
        self.y = y_full
        U = DA.universe(g)
        X = DA.feature_values(U)
        Xv = X.to_numpy(float)
        rr = self.rows
        keep = []
        for j, f in enumerate(U):
            col = Xv[rr, j]
            nn = ~np.isnan(col)
            if M[f]["d_type"] == "binary":
                ok = (col[nn] == 1).sum() >= 5
            else:
                ok = nn.sum() >= 20 and np.nanstd(col) > 0
            keep.append(ok)
        keep = np.array(keep)
        self.features = [f for f, k in zip(U, keep) if k]
        Xv = Xv[:, keep]
        self.binary = np.array([M[f]["d_type"] == "binary" for f in self.features])
        self.carrier = (Xv == 1)
        mu = np.nanmean(Xv, axis=0)
        sd = np.nanstd(Xv, axis=0)
        self.Xs = np.nan_to_num((Xv - mu) / sd, nan=0.0)
        self.stage = np.array([M[f]["stage"] for f in self.features])
        self.forced = np.array([f in DA.forced(g) for f in self.features])
        self.folds = DA.covariates()["cv_fold"].to_numpy()
        self.n_universe_mine = len(U)

    # -------------------------------------------------------------- one fit
    def prep(self, rows):
        R = Residualiser(rows)
        Xt = R.fit(self.Xs[rows])
        yt = R.fit(self.y[rows][:, None])[:, 0]
        sd = Xt.std(0)
        car = self.carrier[rows].sum(0)
        usable = (sd >= SD_MIN) & (~self.binary | (car >= CAR_MIN))
        return R, Xt, yt, usable

    def screen(self, Xt, t, usable, stage):
        cand = np.where(usable & (self.stage == stage))[0]
        if not len(cand):
            return cand
        Xc = Xt[:, cand]
        nr = np.linalg.norm(Xc, axis=0) * max(np.linalg.norm(t), 1e-300)
        cor = np.abs(Xc.T @ t) / nr
        order = np.argsort(-cor, kind="stable")
        top = set(cand[order[:K_SIS]].tolist())
        top |= set(cand[self.forced[cand]].tolist())
        return np.array(sorted(top), dtype=int)

    def fit_fixed(self, rows, aG, aX):
        _, Xt, yt, usable = self.prep(rows)
        beta = np.zeros(len(self.features))
        sG = self.screen(Xt, yt, usable, "G")
        bG = enet_one(Xt[:, sG], yt, aG)
        beta[sG] = bG
        r = yt - Xt[:, sG] @ bG
        sX = self.screen(Xt, r, usable, "X")
        bX = enet_one(Xt[:, sX], r, aX)
        beta[sX] = bX
        return beta

    # -------------------------------------------------------------- lambda by CV + 1-SE
    @staticmethod
    def grid(X, t):
        n = len(t)
        amax = np.abs(X.T @ t).max() / (n * L1) if X.shape[1] else 1.0
        return np.logspace(np.log10(amax), np.log10(amax * EPS), NLAM)

    @staticmethod
    def one_se(mse):
        mse = np.asarray(mse)                      # folds x alphas
        m = mse.mean(0)
        se = mse.std(0, ddof=1) / np.sqrt(mse.shape[0])
        i = int(np.argmin(m))
        ok = np.where(m <= m[i] + se[i])[0]
        return int(ok.min())                       # alphas descending: smallest index = largest alpha

    def select_lambda(self):
        rows = self.rows
        _, Xt, yt, usable = self.prep(rows)
        sG = self.screen(Xt, yt, usable, "G")
        gridG = self.grid(Xt[:, sG], yt)
        folds = sorted(set(self.folds[rows].tolist()))
        cvG, trains = [], []
        for k in folds:
            tr = rows & (self.folds != k)
            ho = rows & (self.folds == k)
            R, Xtr, ytr, us = self.prep(tr)
            Xho = R.apply(ho, self.Xs[tr], self.Xs[ho])
            yho = R.apply(ho, self.y[tr][:, None], self.y[ho][:, None])[:, 0]
            s = self.screen(Xtr, ytr, us, "G")
            C = enet(Xtr[:, s], ytr, gridG)
            pred = Xho[:, s] @ C
            cvG.append(((yho[:, None] - pred) ** 2).mean(0))
            trains.append((tr, ho, Xtr, ytr, us, Xho, yho, s))
        iG = self.one_se(cvG)
        aG = gridG[iG]
        bG = enet_one(Xt[:, sG], yt, aG)
        r = yt - Xt[:, sG] @ bG
        sX = self.screen(Xt, r, usable, "X")
        gridX = self.grid(Xt[:, sX], r)
        cvX = []
        for tr, ho, Xtr, ytr, us, Xho, yho, s in trains:
            b = enet_one(Xtr[:, s], ytr, aG)
            rtr = ytr - Xtr[:, s] @ b
            rho = yho - Xho[:, s] @ b
            sx = self.screen(Xtr, rtr, us, "X")
            C = enet(Xtr[:, sx], rtr, gridX)
            cvX.append(((rho[:, None] - Xho[:, sx] @ C) ** 2).mean(0))
        iX = self.one_se(cvX)
        return aG, gridX[iX], {"iG": iG, "iX": iX, "gridG0": gridG[0], "gridX0": gridX[0]}

    def stability(self, aG, aX):
        ids = DA.dev_ids()
        B = []
        for sub in DA.subsamples():
            rows = self.rows & ids.isin(sub)
            B.append(self.fit_fixed(rows, aG, aX))
        return np.array(B)


def run_gene(g: str, lam: tuple | None = None) -> dict:
    """lam = (alpha_G, alpha_X) skips the CV (diagnostic only: stability conditional on given alphas)."""
    load_params()
    y = DA.effects([g])[g].to_numpy(float)
    G = Gene(g, y)
    if lam is None:
        aG, aX, info = G.select_lambda()
    else:
        (aG, aX), info = lam, {"lambda_source": "given (diagnostic)"}
    full = G.fit_fixed(G.rows, aG, aX)
    B = G.stability(aG, aX)
    nz = B != 0
    freq = nz.mean(0)
    pos = (B > 0).sum(0)
    neg = (B < 0).sum(0)
    nnz = nz.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        sc = np.where(nnz > 0, np.maximum(pos, neg) / np.maximum(nnz, 1), np.nan)
    maj = np.where(pos >= neg, 1.0, -1.0)
    sel = (freq >= FREQ_MIN) & (sc >= SIGN_MIN)
    return {"gene": g, "lam_G": float(aG), "lam_X": float(aX), "info": info, "n_rows": int(G.rows.sum()),
            "features": G.features, "freq": freq.tolist(), "sign_consistency": sc.tolist(),
            "maj_sign": maj.tolist(), "beta_mean": B.mean(0).tolist(), "beta_full": full.tolist(),
            "selected": sel.tolist(), "n_universe_mine": G.n_universe_mine}


def _init_worker():
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[v] = "1"


def _run_star(args):
    return run_gene(*args)


def run_many(genes, workers=3, lams=None):
    from multiprocessing import get_context
    args = [(g, None if lams is None else lams[g]) for g in genes]
    with get_context("spawn").Pool(workers, initializer=_init_worker) as pool:
        return pool.map(_run_star, args, chunksize=1)
