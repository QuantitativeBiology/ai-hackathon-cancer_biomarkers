"""Fixed-threshold joint-row permutation null (DESIGN v3.2 §7.1).

One within-(lineage × library) permutation of line labels per draw (seed permutation_seed + b), applied jointly to
every 𝒟 row: shared and gene-specific feature columns and the NA indicators (CN_absent, MSI_unknown, METH_absent)
move together; targets, predictions, covariates (incl. own CN and EXPGATE) and the effect NaN mask stay fixed.
Strata with < strata_merge_min lines merge into the lineage's largest stratum. Each draw re-runs, on the null gene
subset and at the OBSERVED operating points:
- S1-obs and S1-pred with the observed F1 BH p-cutoffs (p*_F1 per arm; the S1 lists are F1, amendment A6-iii),
- EN-obs and EN-pred stability with the observed per-gene, per-stage λ_1SE, K and frequency/sign rules,
- MD-strict = permuted EN-pred selected ∧ D-stage p <= observed p*_sc ∧ same sign ∧ not collinear.
FDR_emp_fixed(list) = mean_b |hits_b| / max(1, |hits_obs|) on the subset.
"""
from __future__ import annotations

import os
import time

import numpy as np
import pandas as pd

from . import arms as A
from . import params as PR
from . import stats as ST

LISTS = ("S1-obs", "S1-pred", "EN-obs selected", "EN-pred selected", "MD-strict", "MD-global")


def null_genes(data, n: int | None = None, seed: int | None = None) -> list:
    """configs/null_genes.txt when present (real data), else a seeded subset of data.genes (all if fewer)."""
    np_ = PR.null_params()
    n = np_["gene_subset_n"] if n is None else n
    if not data.extras.get("toy"):
        from . import paths as P
        f = P.CONFIGS / "null_genes.txt"
        if f.exists():
            g = [x.strip() for x in f.read_text(encoding="utf-8").splitlines() if x.strip()]
            return [x for x in g if x in set(data.genes)]
    if len(data.genes) <= n:
        return list(data.genes)
    rng = np.random.default_rng(np_["gene_subset_seed"] if seed is None else seed)
    return sorted(rng.choice(np.asarray(data.genes, dtype=object), size=n, replace=False).tolist())


def strata(lineage, library) -> np.ndarray:
    """lineage x library strata; strata with < merge_min lines merge into the lineage's largest stratum."""
    mm = PR.null_params()["strata_merge_min"]
    lin = pd.Series(np.asarray(lineage)).astype(str)
    lib = pd.Series(np.asarray(library)).astype(str)
    s = lin + "|" + lib
    cnt = s.map(s.value_counts())
    big = s[cnt >= mm].groupby(lin[cnt >= mm]).agg(lambda v: v.value_counts().idxmax())
    small = cnt < mm
    s[small] = [big.get(l, f"{l}|*") for l in lin[small]]
    return s.to_numpy()


def permutation(strat: np.ndarray, seed: int) -> np.ndarray:
    """perm such that new_row[i] = old_row[perm[i]], permuting only within strata."""
    rng = np.random.default_rng(seed)
    perm = np.arange(len(strat))
    for k in pd.unique(strat):
        idx = np.flatnonzero(strat == k)
        perm[idx] = idx[rng.permutation(idx.size)]
    return perm


def _draw(data, res_obs: dict, genes: list, b: int, s1_only: bool = False) -> dict:
    """One fixed-threshold draw: returns hit counts per list (on ``genes``) and timing."""
    t0 = time.perf_counter()
    st = strata(data.cov_df["lineage"], data.cov_df["library_group"])
    perm = permutation(st, PR.null_params()["permutation_seed"] + b)
    Xp, NAp, Xsp = data.X[perm], data.na_flags[perm], data.Xstd[perm]
    cut = res_obs["p_cut"]
    tg = {a: data.target(A.TARGET_OF[a]) for a in A.S1_ARMS}
    s1 = ST.s1_scan(ST.ScanInputs(X=Xp, fmeta=data.fmeta, Dg=data.Dg, rows_ok=data.rows_ok, genes=data.genes,
                                  cov=data.cov, own_cn=data.own_cn, expgate=data.expgate, patients=data.patients,
                                  drop=data.drop), tg, compute=True, genes_subset=genes)
    s1["test_family"] = np.where(A.own_expdose_mask(data, s1), "F2", "F1")
    hits = {}
    for a in A.S1_ARMS:
        m = (s1["arm"] == a) & ~s1["collinear_with_design"].astype(bool)
        h = np.zeros(len(s1), dtype=bool)
        for tf in ("F1",):
            c = cut.get(f"{a}|{tf}")
            if c is None:
                continue
            mm = m & (s1["test_family"] == tf)
            h |= (mm & (s1["p"] <= c)).to_numpy()
        hits[a] = int(h.sum())
    out = {"b": b, "s1_rows": int(len(s1)), **{k: hits[k] for k in hits}}
    if not s1_only:
        lam_tab = res_obs["lam"]
        U = s1[s1["arm"] == "S1-obs"][["gene", "feature"]]
        eng = A.en_engine(data, X=Xp, Xstd=Xsp, na_flags=NAp)
        eng.max_patterns = int(PR.runtime()["blocked"]["en_max_patterns"])     # bounds worker memory; cache only
        gf = A.gene_feats(data, U)
        gi = {g: i for i, g in enumerate(data.genes)}
        genes = sorted(genes, key=lambda g: data.rows_ok[:, gi[g]].tobytes())
        en_parts = []
        for arm in A.EN_ARMS:
            lt = lam_tab[lam_tab["arm"] == arm].set_index("gene")
            lam = {g: {"G": lt.at[g, "lam_G"], "X": lt.at[g, "lam_X"]} for g in genes if g in lt.index}
            t, _, _, _ = A.run_en_genes(data, arm, [g for g in genes if g in lam], gf, engine=eng, lam=lam)
            en_parts.append(t)
            out[f"{arm} selected"] = int(t["selected"].sum()) if len(t) else 0
        en_tab = pd.concat([t for t in en_parts if len(t)], ignore_index=True) if any(len(t) for t in en_parts) \
            else pd.DataFrame(columns=["arm", "gene", "feature", "selected"])
        p_sc = res_obs.get("p_cut_sc", 0.0)
        if len(en_tab):
            md = A.strict(en_tab, s1, "EN-pred", p_cut=p_sc)
            out["MD-strict"] = int(md["pass"].sum())
            mg = A.md_global(en_tab, s1, "EN-pred", p_cut=cut.get("S1-obs|F1", 0.0))   # A7 at the observed p*_F1
            out["MD-global"] = int(mg["pass"].sum())
        else:
            out["MD-strict"] = out["MD-global"] = 0
    out["wall_s"] = time.perf_counter() - t0
    return out


_WN = {}


def _winit(payload):
    if isinstance(payload, str):                      # workerstore directory (file-backed, mmap)
        from . import workerstore as WS
        payload = WS.load(payload)
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = str(payload["blas_threads"])
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(payload["blas_threads"])
    except Exception:  # noqa: BLE001
        pass
    _WN.update(payload)


def _wdraw(args):
    b, s1_only = args
    return _draw(_WN["data"], _WN["res"], _WN["genes"], b, s1_only)


def observed_hits(res_obs: dict, genes: list) -> dict:
    gs = set(genes)
    a = res_obs["arms"]
    out = {}
    for arm in A.S1_ARMS:
        s = a[(a["arm"] == arm) & a["gene"].isin(gs)]
        out[arm] = int(s["selected"].sum())
    for arm in A.EN_ARMS:
        s = a[(a["arm"] == arm) & a["gene"].isin(gs)]
        out[f"{arm} selected"] = int(s["selected"].sum())
    md = res_obs.get("md_strict")
    out["MD-strict"] = int(md["gene"].isin(gs).sum()) if md is not None else 0
    mg = res_obs.get("md_global")
    out["MD-global"] = int(mg["gene"].isin(gs).sum()) if mg is not None else 0
    return out


def summarise(draws: pd.DataFrame, obs: dict) -> pd.DataFrame:
    npar = PR.null_params()
    rng = np.random.default_rng(npar["permutation_seed"])
    rows = []
    for L in LISTS:
        if L not in draws:
            continue
        h = draws[L].to_numpy(dtype=float)
        ho = obs.get(L, 0)
        den = max(1, ho)
        fdr = float(h.mean() / den)
        bs = [rng.choice(h, size=h.size, replace=True).mean() / den for _ in range(npar["boot_B"])]
        lo, hi = np.quantile(bs, [(1 - npar["ci_level"]) / 2, 1 - (1 - npar["ci_level"]) / 2])
        rows.append({"list": L, "fdr_emp_fixed": fdr, "ci_lo": float(lo), "ci_hi": float(hi), "hits_obs": ho,
                     "mean_hits_perm": float(h.mean()), "p_any_hit": float((h >= 1).mean()), "B": int(h.size),
                     "pass": fdr <= npar["fdr_max"]})
    return pd.DataFrame(rows)


def run_null(data, res_obs: dict, B: int | None = None, workers: int = 1, genes=None, decide_B: bool = False,
             s1_only: bool = False, log=None) -> dict:
    """Fixed-threshold null. ``decide_B``: B = B_raised iff B_raised x wall(draw 0) <= wall_cap (timing only).
    Returns {'draws': per-draw counts, 'summary': FDR table, 'B': B, 'wall_draw0_s': ...}."""
    npar = PR.null_params()
    genes = null_genes(data) if genes is None else list(genes)
    B = npar["B_base"] if B is None else int(B)
    draws = []
    data.target("obs")
    payload = {"data": data.for_workers(include_obs=True), "res": {k: res_obs[k] for k in ("p_cut", "lam", "p_cut_sc")
                                                                    if k in res_obs},
               "genes": genes, "blas_threads": npar["blas_threads"]}
    if workers <= 1:
        _WN.update(payload)
        d0 = _draw(data, payload["res"], genes, 0, s1_only)
        draws.append(d0)
        if decide_B and npar["B_raised"] * d0["wall_s"] / npar["sec_per_min"] <= npar["wall_cap_min"]:
            B = npar["B_raised"]
        for b in range(1, B):
            draws.append(_draw(data, payload["res"], genes, b, s1_only))
            if log:
                log(draws[-1])
    else:
        from concurrent.futures import ProcessPoolExecutor

        from . import workerstore as WS
        with WS.staged(payload, "null") as path, \
                ProcessPoolExecutor(max_workers=workers, initializer=_winit, initargs=(path,)) as ex:
            d0 = ex.submit(_wdraw, (0, s1_only)).result()
            draws.append(d0)
            if decide_B and npar["B_raised"] * d0["wall_s"] / npar["sec_per_min"] <= npar["wall_cap_min"]:
                B = npar["B_raised"]
            for d in ex.map(_wdraw, [(b, s1_only) for b in range(1, B)]):
                draws.append(d)
                if log:
                    log(d)
    dr = pd.DataFrame(draws)
    obs = observed_hits(res_obs, genes)
    return {"draws": dr, "summary": summarise(dr, obs), "B": B, "wall_draw0_s": draws[0]["wall_s"],
            "genes": genes, "hits_obs": obs}
