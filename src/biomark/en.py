"""Hierarchical two-stage (G then X) fixed-λ stability elastic net (DESIGN v3.2 §4.1). Identical code for EN-pred
(target = fold-centred prediction) and EN-obs (target = observed effect): the caller passes the target vector, the
gene's row mask (effect non-NaN; for EN-pred from the NaN mask only) and the gene's feature positions in U.

Per fit (full D, each CV training fold, each of the 50 subsamples), on that fit's rows:
1. FWL: target and 𝒟 (standardised once on D, NaN -> 0) are residualised on the unpenalised block [1, lineage,
   library group, culture, growth, CN_absent, MSI_unknown, METH_absent] estimated on the fit rows; no intercept.
   A column is unavailable in a fit if residual SD < resid_sd_lt or (binary) < carriers_lt carriers.
2. Stage G on the residualised target, stage X on the fit's own stage-G residual.
3. Per stage SIS to K by |corr| plus forced terms of that stage (still penalised).
4. λ: l1_ratio, n_lambda path (grid from full D, eps), CV over Part 1 cv_fold, 1-SE on full D; fixed per gene and
   stage in every subsample (and in the fixed-threshold null).
5. Stability: selected = frequency >= freq_min over the 50 subsamples and >= sign_min of selections with the
   majority sign; rank by frequency, |mean β over all 50 (zeros included)|, feature id; zero frequency unrankable.
6. Fidelity r2_oof_beyond_design: pooled out-of-fold R² of the surrogate on the held-out residualised target.

The penalised problem is solved by coordinate descent on the Gram matrix (scikit-learn's cd_fast, the solver behind
``enet_path``), objective (1/2n)||y - Xb||² + α(ρ||b||₁ + (1-ρ)/2 ||b||²).
"""
from __future__ import annotations

import hashlib
import json
import warnings
from collections import OrderedDict

import numpy as np
import pandas as pd
from sklearn.linear_model import _cd_fast as cd_fast

from . import params as PR
from . import stats as ST

# ------------------------------------------------------------------------------------------------
# Subsamples (configs/subsamples.json)
# ------------------------------------------------------------------------------------------------


def make_subsamples(lines, lineage, patients, seed: int, n_pairs: int) -> list:
    """``n_pairs`` complementary half-sample pairs, patient-grouped and lineage-stratified: each patient is
    stratified by its modal lineage; within a stratum the patients are shuffled and each is put in the half with
    fewer lines of that stratum so far (ties at random). Returns 2 * n_pairs lists of line ids (pair i -> 2i, 2i+1)."""
    lines = pd.Index(lines)
    lin = pd.Series(np.asarray(lineage), index=lines).astype(str)
    pat = pd.Series(np.asarray(patients), index=lines).astype(str)
    size = pat.value_counts()
    strat = lin.groupby(pat).agg(lambda v: v.value_counts().sort_index().idxmax())
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n_pairs):
        a = set()
        for L in sorted(strat.unique()):
            pts = np.array(sorted(strat.index[strat == L]), dtype=object)
            pts = pts[rng.permutation(len(pts))]
            na = nb = 0
            for p in pts:
                k = int(size[p])
                if na < nb or (na == nb and rng.random() < 0.5):
                    a.add(p)
                    na += k
                else:
                    nb += k
        ina = pat.isin(a).to_numpy()
        out.append(list(lines[ina]))
        out.append(list(lines[~ina]))
    return out


def subsample_manifest(subs: list, seed: int) -> dict:
    return {"format": PR.pp()["subsamples"]["format"], "seed_key": "subsample_seed", "seed": int(seed),
            "rule": PR.pp()["subsamples"]["half_rule"], "n_pairs": len(subs) // 2, "subsamples": subs}


def read_subsamples(path) -> list:
    """Accepts this module's format ({'subsamples': [[ids], ...]}) or {'pairs': [{'a': [...], 'b': [...]}]}."""
    with open(path, encoding="utf-8") as fh:
        m = json.load(fh)
    if "subsamples" in m:
        return [list(map(str, s)) for s in m["subsamples"]]
    if "pairs" in m:
        out = []
        for pr in m["pairs"]:
            out += [list(map(str, pr["a"])), list(map(str, pr["b"]))]
        return out
    raise ValueError(f"unrecognised subsample manifest {path}")


def max_attainable(binary_vals: np.ndarray, rows_ok: np.ndarray, sub_masks: list, carriers_min: int) -> np.ndarray:
    """Fraction of subsamples in which each binary column has >= carriers_min carriers among rows_ok rows."""
    B = np.nan_to_num(binary_vals) == 1
    cnt = np.array([(B[rows_ok & s]).sum(0) for s in sub_masks])
    return (cnt >= carriers_min).mean(0)


# ------------------------------------------------------------------------------------------------
# Coordinate descent on the Gram matrix
# ------------------------------------------------------------------------------------------------


def enet_gram(G: np.ndarray, Xy: np.ndarray, y: np.ndarray, n: int, alphas, l1_ratio: float, tol: float,
              max_iter: int, coef_init=None) -> np.ndarray:
    """Elastic-net path (alphas in the given order, warm-started) from the Gram matrix; returns p x len(alphas)."""
    p = G.shape[0]
    alphas = np.atleast_1d(np.asarray(alphas, dtype=np.float64))
    out = np.zeros((p, alphas.size))
    if p == 0:
        return out
    G = np.ascontiguousarray(G, dtype=np.float64)
    Xy = np.ascontiguousarray(Xy, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)
    w = np.zeros(p) if coef_init is None else np.array(coef_init, dtype=np.float64)
    rng = np.random.RandomState(0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for i, a in enumerate(alphas):
            w, _, _, _ = cd_fast.enet_coordinate_descent_gram(w, a * l1_ratio * n, a * (1.0 - l1_ratio) * n, G, Xy, y,
                                                               max_iter, tol, rng, 0, 0)
            out[:, i] = w
    return out


def alpha_grid(Xy: np.ndarray, n: int, l1_ratio: float, n_lambda: int, eps: float) -> np.ndarray:
    amax = np.max(np.abs(Xy)) / (n * l1_ratio) if Xy.size else 0.0
    if not np.isfinite(amax) or amax <= 0:
        return np.zeros(0)
    return np.geomspace(amax, amax * eps, n_lambda)


def one_se(mse: np.ndarray) -> int:
    """mse: folds x alphas (NaN folds ignored); index of the largest α within one SE of the minimum."""
    ok = ~np.isnan(mse).all(1)
    m = mse[ok]
    mean = m.mean(0)
    k = m.shape[0]
    se = m.std(0, ddof=1) / np.sqrt(k) if k > 1 else np.zeros(m.shape[1])
    j = int(np.argmin(mean))
    return int(np.flatnonzero(mean <= mean[j] + se[j])[0])


# ------------------------------------------------------------------------------------------------
# Engine
# ------------------------------------------------------------------------------------------------


class ENEngine:
    """Per-dataset EN engine. ``X`` raw feature values (n x F; NaN = not assayed), ``fmeta`` aligned (family, stage,
    d_type, scope), ``cov`` CovDesign over the same rows, ``na_flags`` n x k NA indicators (move with 𝒟 rows under the
    null), ``folds`` cv_fold per row, ``sub_masks`` list of 50 boolean row masks. ``Xstd`` (optional) = 𝒟
    standardised on D, NaN -> 0 (computed here when omitted)."""

    def __init__(self, X: np.ndarray, fmeta: pd.DataFrame, cov: ST.CovDesign, na_flags: np.ndarray, folds,
                 sub_masks: list, Xstd: np.ndarray | None = None, ep: dict | None = None):
        self.ep = ep or PR.en_params()
        self.n, self.F = X.shape
        self.binary = (fmeta["d_type"].to_numpy() == "binary")
        self.stage = fmeta["stage"].to_numpy()
        self.shared = np.flatnonzero(fmeta["scope"].to_numpy() == "shared")
        self.shared_pos = -np.ones(self.F, dtype=np.int64)
        self.shared_pos[self.shared] = np.arange(self.shared.size)
        self.fid = fmeta.index.to_numpy()
        if Xstd is None:
            Xstd = standardise(X)
        self.Xs = Xstd
        self.car01 = ((X == 1) & self.binary[None, :]).astype(np.float32)
        self.cov = cov
        self.na = np.asarray(na_flags, dtype=np.float64)
        self.folds = np.asarray(folds)
        self.fold_ids = sorted(np.unique(self.folds).tolist())
        self.sub = [np.asarray(s, dtype=bool) for s in sub_masks]
        self.cache: OrderedDict = OrderedDict()
        self.hits = self.miss = 0
        self.gram_shared = self.shared.size <= self.ep["gram_cache_max"]

    # ---- FWL ----------------------------------------------------------------------------------------------
    def U(self, rows_fit, rows_apply=None) -> np.ndarray:
        ap = rows_fit if rows_apply is None else rows_apply
        Z = self.cov.Z(rows_fit, self.cov.design_names, apply_rows=ap)
        return np.hstack([Z, self.na[ap]])

    def _entry(self, gkey: bytes, grows: np.ndarray, fit: str) -> dict:
        key = (gkey, fit)
        if key in self.cache:
            self.hits += 1
            self.cache.move_to_end(key)
            return self.cache[key]
        self.miss += 1
        if fit == "full":
            rows = grows
        elif fit.startswith("cv"):
            rows = grows & (self.folds != int(fit[2:]))
        else:
            rows = grows & self.sub[int(fit[2:])]
        Uf = self.U(rows)
        Q = ST.orth(Uf)
        e = {"rows": rows, "Q": Q, "n": int(rows.sum())}
        Xsh = ST.resid(Q, self.Xs[rows][:, self.shared])
        e["car"] = self.car01[rows][:, self.shared].sum(0)
        e["sd"] = np.sqrt((Xsh ** 2).mean(0)) if e["n"] else np.zeros(self.shared.size)
        if self.gram_shared:
            e["G"] = Xsh.T @ Xsh
        e["Xt"] = Xsh.astype(np.float32)
        if fit.startswith("cv"):
            te = grows & (self.folds == int(fit[2:]))
            e["te"] = te
            Ute = self.U(rows, te)
            coefU, *_ = np.linalg.lstsq(Uf, self.Xs[rows][:, self.shared], rcond=None)
            e["Xte"] = (self.Xs[te][:, self.shared] - Ute @ coefU).astype(np.float32)
            e["Ute"], e["Uf"] = Ute, Uf
        self.cache[key] = e
        while len(self.cache) > self.ep["cache_entries"]:
            self.cache.popitem(last=False)
        return e

    def _gene_block(self, e: dict, feats: np.ndarray) -> "_Block":
        return _Block(self, e, feats)

    def _target(self, e, y):
        yt = ST.resid(e["Q"], y[e["rows"]][:, None])[:, 0]
        yte = None
        if "te" in e:
            cy, *_ = np.linalg.lstsq(e["Uf"], y[e["rows"]], rcond=None)
            yte = y[e["te"]] - e["Ute"] @ cy
        return yt, yte

    @staticmethod
    def _sis(B, r, cand, forced_mask, K):
        """cand: positions available in this stage; returns sorted positions: top-K |corr| plus forced."""
        if cand.size == 0:
            return cand
        f = cand[forced_mask[cand]]
        rest = cand[~forced_mask[cand]]
        if rest.size > K:
            nrm = B.norm[rest]
            c = np.abs(B.dot(r)[rest]) / np.where(nrm > 0, nrm, np.inf)
            top = np.argsort(-c, kind="stable")[:K]
            rest = rest[top]
        return np.sort(np.concatenate([f, rest]))

    # ---- one fit ------------------------------------------------------------------------------------------
    def _stage(self, B, r, cols, alphas):
        n = B.n
        if cols.size == 0 or alphas.size == 0:
            return np.zeros((cols.size, max(alphas.size, 1)))
        ep = self.ep
        if alphas.size == 1 and not B.gram_cached(cols):
            return enet_x(B.cols(cols), r, n, alphas[0], ep["l1_ratio"], ep["tol"], ep["max_iter"])[:, None]
        return enet_gram(B.gram(cols), B.dot(r)[cols], r, n, alphas, ep["l1_ratio"], ep["tol"], ep["max_iter"])

    def _fixed(self, B, r, cols, alpha):
        if cols.size == 0 or not np.isfinite(alpha):
            return np.zeros(cols.size)
        return self._stage(B, r, cols, np.array([alpha]))[:, 0]

    def fit_gene(self, y: np.ndarray, grows: np.ndarray, feats: np.ndarray, forced: np.ndarray,
                 lam: dict | None = None, K: int | None = None) -> dict:
        """Fit one gene. ``y`` n-vector (NaN outside ``grows`` allowed), ``feats`` positions of U_g features,
        ``forced`` boolean over ``feats``. ``lam`` = {'G': α_G, 'X': α_X} to skip CV (fixed-threshold null).
        Returns per-feature arrays (freq, sign_cons, mean_beta, beta_full, n_avail, selected) and lam/fidelity."""
        ep = self.ep
        K = ep["K"] if K is None else K
        feats = np.asarray(feats, dtype=np.int64)
        forced = np.asarray(forced, dtype=bool)
        y = np.where(np.isfinite(y), y, 0.0).astype(np.float64)
        gkey = hashlib.blake2b(np.packbits(grows).tobytes(), digest_size=16).digest()  # literal-ok: hash width
        stG = self.stage[feats] == "G"
        stX = ~stG
        # ---- full-D fit and λ by CV ----------------------------------------------------------------------
        e0 = self._entry(gkey, grows, "full")
        B0 = self._gene_block(e0, feats)
        y0, _ = self._target(e0, y)
        n0 = e0["n"]
        out = {"r2_oof": np.nan, "K": K}
        selG0 = self._sis(B0, y0, np.flatnonzero(B0.avail & stG), forced, K)
        gridG = alpha_grid(B0.dot(y0)[selG0], n0, ep["l1_ratio"], ep["n_lambda"], ep["eps"])
        cv = []
        if lam is None:
            for k in self.fold_ids:
                ek = self._entry(gkey, grows, f"cv{k}")
                if ek["n"] == 0 or int(ek["te"].sum()) == 0:
                    continue
                Bk = self._gene_block(ek, feats)
                yt, yte = self._target(ek, y)
                sel = self._sis(Bk, yt, np.flatnonzero(Bk.avail & stG), forced, K)
                path = self._stage(Bk, yt, sel, gridG)
                cv.append({"B": Bk, "yt": yt, "yte": yte, "selG": sel, "pathG": path})
            if gridG.size and cv:
                mse = np.array([((c["yte"][:, None] - c["B"].te_cols(c["selG"]) @ c["pathG"]) ** 2).mean(0)
                                for c in cv])
                iG = one_se(mse)
                aG = float(gridG[iG])
            else:
                iG, aG = -1, np.inf
        else:
            aG = float(lam["G"])
        bG0 = self._fixed(B0, y0, selG0, aG)
        r0 = y0 - B0.matmul(selG0, bG0)
        selX0 = self._sis(B0, r0, np.flatnonzero(B0.avail & stX), forced, K)
        gridX = alpha_grid(B0.dot(r0)[selX0], n0, ep["l1_ratio"], ep["n_lambda"], ep["eps"])
        if lam is None:
            mseX = []
            for c in cv:
                Bk = c["B"]
                bG = c["pathG"][:, iG] if iG >= 0 else np.zeros(c["selG"].size)
                rt = c["yt"] - Bk.matmul(c["selG"], bG)
                rte = c["yte"] - Bk.te_matmul(c["selG"], bG)
                selX = self._sis(Bk, rt, np.flatnonzero(Bk.avail & stX), forced, K)
                pathX = self._stage(Bk, rt, selX, gridX) if gridX.size else None
                c.update(selX=selX, pathX=pathX, rte=rte)
                if pathX is not None:
                    mseX.append(((rte[:, None] - Bk.te_cols(selX) @ pathX) ** 2).mean(0))
            if gridX.size and mseX:
                iX = one_se(np.array(mseX))
                aX = float(gridX[iX])
            else:
                iX, aX = -1, np.inf
            sse = sst = 0.0
            for c in cv:
                pred = c["yte"] - c["rte"]
                if iX >= 0 and c.get("pathX") is not None:
                    pred = pred + c["B"].te_matmul(c["selX"], c["pathX"][:, iX])
                sse += float(((c["yte"] - pred) ** 2).sum())
                sst += float((c["yte"] ** 2).sum())
            out["r2_oof"] = 1.0 - sse / sst if sst > 0 else np.nan
        else:
            aX = float(lam["X"])
        bX0 = self._fixed(B0, r0, selX0, aX)
        beta_full = np.zeros(feats.size)
        beta_full[selG0] = bG0
        beta_full[selX0] = bX0
        out.update(lam_G=aG, lam_X=aX, beta_full=beta_full)
        # ---- stability --------------------------------------------------------------------------------------
        nsub = len(self.sub)
        Bs = np.zeros((nsub, feats.size))
        avail = np.zeros(feats.size)
        for s in range(nsub):
            es = self._entry(gkey, grows, f"ss{s}")
            if es["n"] == 0:
                continue
            Bb = self._gene_block(es, feats)
            avail += Bb.avail
            yt, _ = self._target(es, y)
            sel = self._sis(Bb, yt, np.flatnonzero(Bb.avail & stG), forced, K)
            bG = self._fixed(Bb, yt, sel, aG)
            Bs[s, sel] = bG
            rt = yt - Bb.matmul(sel, bG)
            selX = self._sis(Bb, rt, np.flatnonzero(Bb.avail & stX), forced, K)
            Bs[s, selX] = self._fixed(Bb, rt, selX, aX)
        nz = Bs != 0
        nsel = nz.sum(0)
        freq = nsel / nsub
        pos = (Bs > 0).sum(0)
        with np.errstate(invalid="ignore", divide="ignore"):
            maj = np.where(nsel > 0, np.maximum(pos, nsel - pos) / np.maximum(nsel, 1), np.nan)
        sign = np.where(pos * 2 >= nsel, 1.0, -1.0)
        sign[nsel == 0] = 0.0
        out.update(freq=freq, sign_cons=maj, maj_sign=sign, mean_beta=Bs.mean(0), n_avail=avail,
                   max_attainable=avail / nsub,
                   selected=(freq >= ep["freq_min"]) & (np.nan_to_num(maj) >= ep["sign_min"]), stage_G=stG,
                   betas=Bs)
        return out


class _Block:
    """Lazy residualised columns of one gene's U_g features in one fit: X'r for all columns by one mat-vec on the
    cached shared block; float64 columns materialised only for screened / non-zero sets."""

    def __init__(self, eng: "ENEngine", e: dict, feats: np.ndarray):
        self.e, self.n = e, e["n"]
        sp = eng.shared_pos[feats]
        self.is_sh = sp >= 0
        self.sp = sp
        self.spec_pos = np.flatnonzero(~self.is_sh)
        spec = feats[~self.is_sh]
        self.spec_idx = -np.ones(feats.size, dtype=np.int64)
        self.spec_idx[self.spec_pos] = np.arange(spec.size)
        car = np.empty(feats.size)
        sd = np.empty(feats.size)
        car[self.is_sh], sd[self.is_sh] = e["car"][sp[self.is_sh]], e["sd"][sp[self.is_sh]]
        self.Xsp = self.Xsp_te = None
        if spec.size:
            ri = e.setdefault("row_idx", np.flatnonzero(e["rows"]))
            Xr = eng.Xs[np.ix_(ri, spec)]
            self.Xsp = ST.resid(e["Q"], Xr)
            car[~self.is_sh] = eng.car01[np.ix_(ri, spec)].sum(0)
            sd[~self.is_sh] = np.sqrt((self.Xsp ** 2).mean(0)) if self.n else 0.0
            if "te" in e:
                coefU, *_ = np.linalg.lstsq(e["Uf"], Xr, rcond=None)
                ti = e.setdefault("te_idx", np.flatnonzero(e["te"]))
                self.Xsp_te = eng.Xs[np.ix_(ti, spec)] - e["Ute"] @ coefU
        self.avail = (sd >= eng.ep["resid_sd_min"]) & (~eng.binary[feats] | (car >= eng.ep["carriers_min"]))
        self.norm = np.sqrt(self.n) * sd
        self._dot_key, self._dot_val = None, None

    def dot(self, r: np.ndarray) -> np.ndarray:
        if self._dot_key is r:
            return self._dot_val
        out = np.empty(self.sp.size)
        if self.is_sh.any():
            full = self.e["Xt"].T @ r.astype(np.float32)
            out[self.is_sh] = full[self.sp[self.is_sh]]
        if self.Xsp is not None:
            out[self.spec_pos] = self.Xsp.T @ r
        self._dot_key, self._dot_val = r, out
        return out

    def _get(self, cols, test: bool) -> np.ndarray:
        cols = np.asarray(cols, dtype=np.int64)
        src = self.e["Xte"] if test else self.e["Xt"]
        out = np.empty((src.shape[0], cols.size))
        sh = self.is_sh[cols]
        if sh.any():
            out[:, sh] = src[:, self.sp[cols[sh]]]
        if (~sh).any():
            S = self.Xsp_te if test else self.Xsp
            out[:, ~sh] = S[:, self.spec_idx[cols[~sh]]]
        return out

    def cols(self, cols) -> np.ndarray:
        return self._get(cols, False)

    def te_cols(self, cols) -> np.ndarray:
        return self._get(cols, True)

    def matmul(self, cols, b) -> np.ndarray:
        nz = np.flatnonzero(b)
        return self.cols(cols[nz]) @ b[nz] if nz.size else np.zeros(self.n)

    def te_matmul(self, cols, b) -> np.ndarray:
        nz = np.flatnonzero(b)
        return self.te_cols(cols[nz]) @ b[nz] if nz.size else np.zeros(self.e["Xte"].shape[0])

    def gram_cached(self, cols) -> bool:
        """True when the shared block's Gram is cached (gene-specific columns are added by cross products)."""
        return "G" in self.e

    def gram(self, cols) -> np.ndarray:
        cols = np.asarray(cols, dtype=np.int64)
        if not self.gram_cached(cols):
            A = self.cols(cols)
            return A.T @ A
        sh = self.is_sh[cols]
        G = np.empty((cols.size, cols.size))
        s = self.sp[cols[sh]]
        G[np.ix_(sh, sh)] = self.e["G"][np.ix_(s, s)]
        if (~sh).any():
            P_ = self.Xsp[:, self.spec_idx[cols[~sh]]]
            A = self.cols(cols)
            C = A.T @ P_
            G[:, ~sh] = C
            G[~sh, :] = C.T
        return G


def enet_x(X: np.ndarray, y: np.ndarray, n: int, alpha: float, l1_ratio: float, tol: float, max_iter: int) -> np.ndarray:
    """Single-α elastic net by coordinate descent on X (same objective and convergence rule as ``enet_gram``)."""
    Xf = np.asfortranarray(X, dtype=np.float64)
    w = np.zeros(Xf.shape[1])
    rng = np.random.RandomState(0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        w, _, _, _ = cd_fast.enet_coordinate_descent(w, alpha * l1_ratio * n, alpha * (1.0 - l1_ratio) * n, Xf,
                                                     np.ascontiguousarray(y, dtype=np.float64), max_iter, tol, rng,
                                                     0, 0)
    return w


def standardise(X: np.ndarray) -> np.ndarray:
    """𝒟 standardised once on D (column mean/SD over non-NaN rows), NaN -> 0."""
    X = np.asarray(X, dtype=np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        mu = np.nanmean(X, 0)
        sd = np.nanstd(X, 0)
    sd = np.where(np.isfinite(sd) & (sd > 0), sd, 1.0)
    Z = (X - np.nan_to_num(mu)) / sd
    return np.where(np.isnan(Z), 0.0, Z)


def rank_en(freq: np.ndarray, mean_beta: np.ndarray, fid: np.ndarray) -> np.ndarray:
    """Rank within a gene-agnostic list: frequency desc, |mean β| desc, feature id asc; zero frequency -> NaN."""
    df = pd.DataFrame({"f": -np.asarray(freq, float), "b": -np.abs(np.asarray(mean_beta, float)),
                       "id": np.asarray(fid, dtype=object)})
    order = df.sort_values(["f", "b", "id"], kind="mergesort").index.to_numpy()
    rk = np.empty(len(df))
    rk[order] = np.arange(1, len(df) + 1)
    rk[np.asarray(freq) <= 0] = np.nan
    return rk
