"""The four benchmark arms on one (gene, 𝒟) universe and MD-strict (DESIGN v3.2 §5.1-§5.2, §6 ranking rules).

- U = {(g, d): g in data.genes, d in 𝒟_g minus self-null, §3 counts on g's non-NaN rows (binary >= min carriers,
  continuous >= continuous_min_nonnan lines and SD > 0)}; identical for every arm.
- S1-obs / S1-pred: stats.s1_scan on observed effects / fold-centred predictions (same rows, drops, SEs);
  rank_stat = p (rank p asc, |β| desc, feature id asc); selected = BH q < s1_s5.bh_q within (arm, F1|F2), not
  collinear_with_design.
- EN-obs / EN-pred: en.ENEngine.fit_gene (same rows, subsamples, configs); rank_stat = selection frequency (rank
  frequency desc, |mean β over 50| desc, feature id asc; frequency 0 unrankable); selected = frequency and sign rule.
- MD-strict = EN-pred selected ∧ D-stage (the S1-obs confounder-model test of the same (g, d)) same sign ∧
  BH q_sc < md_strict_s5.q_sc_bh within EN-pred selected candidates ∧ not collinear_with_design. EN-obs-strict alike.

EN-pred never reads observed effects: it consumes ``data.target('pred')`` and the row mask only.
"""
from __future__ import annotations

import os
import time

import numpy as np
import pandas as pd

from . import en as EN
from . import params as PR
from . import stats as ST

ARMS = ("S1-obs", "S1-pred", "EN-obs", "EN-pred")
S1_ARMS = ("S1-obs", "S1-pred")
EN_ARMS = ("EN-obs", "EN-pred")
TARGET_OF = {"S1-obs": "obs", "S1-pred": "pred", "EN-obs": "obs", "EN-pred": "pred"}


# ------------------------------------------------------------------------------------------------
# Universe and S1
# ------------------------------------------------------------------------------------------------


def scan_inputs(data, X=None, na_flags=None) -> ST.ScanInputs:
    return ST.ScanInputs(X=data.X if X is None else X, fmeta=data.fmeta, Dg=data.Dg, rows_ok=data.rows_ok,
                         genes=data.genes, cov=data.cov, own_cn=data.own_cn, expgate=data.expgate,
                         patients=data.patients, drop=data.drop)


def run_s1(data, arms=S1_ARMS, genes=None, X=None, workers: int = 1, chunk: int | None = None) -> pd.DataFrame:
    """S1 scan for the requested S1 arms (or the D-stage 'S1-obs' only); one row per (arm, g, d) in U.
    ``workers`` > 1: gene chunks (sorted by NaN pattern) over a process pool; results identical to serial."""
    genes = list(data.genes if genes is None else genes)
    chunk = chunk or PR.null_params()["chunk_genes"]
    if workers <= 1 or len(genes) <= chunk or X is not None:
        tg = {a: data.target(TARGET_OF[a]) for a in arms}
        df = ST.s1_scan(scan_inputs(data, X), tg, compute=True, genes_subset=genes)
        return finish_s1(data, df)
    from concurrent.futures import ProcessPoolExecutor
    gi = {g: i for i, g in enumerate(data.genes)}
    order = sorted(genes, key=lambda g: data.rows_ok[:, gi[g]].tobytes())
    jobs = [(tuple(arms), order[i:i + chunk]) for i in range(0, len(order), chunk)]
    payload = {"data": data.for_workers(include_obs=any(TARGET_OF[a] == "obs" for a in arms)), "gf": None,
               "blas_threads": PR.null_params()["blas_threads"]}
    with ProcessPoolExecutor(max_workers=workers, initializer=_worker_init, initargs=(payload,)) as ex:
        parts = list(ex.map(_worker_s1, jobs))
    df = pd.concat([p for p in parts if len(p)], ignore_index=True)
    key = pd.Categorical(df["gene"], categories=genes, ordered=True)
    df = df.assign(_k=key).sort_values(["arm", "_k"], kind="mergesort").drop(columns="_k").reset_index(drop=True)
    return finish_s1(data, df)


def _worker_s1(args):
    arms, genes = args
    data = _W["data"]
    tg = {a: data.target(TARGET_OF[a]) for a in arms}
    return ST.s1_scan(scan_inputs(data), tg, compute=True, genes_subset=genes)


def finish_s1(data, df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    own = own_expdose_mask(data, df)
    df = ST.add_bh(df, own)
    df["rank_stat"] = df["p"]
    df["stage"] = ""
    df["rank"] = np.nan
    for a, idx in df.groupby("arm").groups.items():
        sub = df.loc[idx]
        ok = sub["p"].notna()
        o = sub[ok].assign(_b=-sub.loc[ok, "beta"].abs()).sort_values(["p", "_b", "feature"], kind="mergesort").index
        df.loc[o, "rank"] = np.arange(1, len(o) + 1)
    return df


def universe(data, genes=None) -> pd.DataFrame:
    return ST.s1_scan(scan_inputs(data), {}, compute=False, genes_subset=genes)


def own_expdose_mask(data, df: pd.DataFrame) -> np.ndarray:
    fam = data.fmeta["family"].reindex(df["feature"]).to_numpy()
    fg = data.fmeta["genes"].reindex(df["feature"]).to_numpy()
    expd = PR.pp().get("own_expdose_family", "EXPDOSE")
    return np.array([(f == expd) and (g in set(x)) for f, g, x in zip(fam, df["gene"], fg)], dtype=bool)


# ------------------------------------------------------------------------------------------------
# EN arms
# ------------------------------------------------------------------------------------------------


def gene_feats(data, U: pd.DataFrame) -> dict:
    """gene -> (feature positions in U_g, forced mask)."""
    pos = pd.Series(np.arange(len(data.fmeta)), index=data.fmeta.index)
    out = {}
    for g, sub in U.groupby("gene", sort=False):
        p = np.sort(pos.reindex(sub["feature"].unique()).to_numpy().astype(np.int64))
        fset = data.forced.get(g, set())
        out[g] = (p, np.array([data.fmeta.index[i] in fset for i in p], dtype=bool))
    return out


def en_engine(data, X=None, Xstd=None, na_flags=None) -> EN.ENEngine:
    return EN.ENEngine(data.X if X is None else X, data.fmeta, data.cov,
                       data.na_flags if na_flags is None else na_flags, data.folds, data.sub_masks,
                       Xstd=data.Xstd if Xstd is None else Xstd)


def run_en_genes(data, arm: str, genes, gf: dict, engine: EN.ENEngine | None = None, lam: dict | None = None,
                 K: int | None = None, keep_betas: bool = False):
    """Fit EN for ``genes``; returns (long table of U rows, per-gene lambda/fidelity table, timing dict)."""
    eng = engine or en_engine(data)
    Y = data.target(TARGET_OF[arm])
    gpos = {g: i for i, g in enumerate(data.genes)}
    rows, lt, betas = [], [], {}
    t0 = time.perf_counter()
    for g in genes:
        if g not in gf:
            continue
        feats, forced = gf[g]
        j = gpos[g]
        lg = None if lam is None else lam.get(g)
        if lam is not None and lg is None:
            continue
        r = eng.fit_gene(Y[:, j], data.rows_ok[:, j], feats, forced, lam=lg, K=K)
        fid = data.fmeta.index.to_numpy()[feats]
        rows.append(pd.DataFrame({
            "arm": arm, "gene": g, "feature": fid, "rank_stat": r["freq"], "selection_freq": r["freq"],
            "sign_consistency": r["sign_cons"], "maj_sign": r["maj_sign"], "beta": r["mean_beta"],
            "beta_full": r["beta_full"], "selected": r["selected"],
            "stage": np.where(r["stage_G"], "G", "X"), "en_max_attainable_freq": r["max_attainable"],
            "forced": forced}))
        lt.append({"arm": arm, "gene": g, "lam_G": r["lam_G"], "lam_X": r["lam_X"], "K": r["K"],
                   "r2_oof_beyond_design": r["r2_oof"], "n_rows": int(data.rows_ok[:, j].sum())})
        if keep_betas:
            betas[g] = r["betas"]
    tim = {"arm": arm, "n_genes": len(lt), "wall_s": time.perf_counter() - t0, "cache_hits": eng.hits,
           "cache_misses": eng.miss}
    tab = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    return tab, pd.DataFrame(lt), tim, betas


def finish_en(tab: pd.DataFrame) -> pd.DataFrame:
    if tab.empty:
        return tab
    tab = tab.copy()
    unreach = PR.en_params()["unreachable_lt"]
    tab["en_unreachable"] = tab["en_max_attainable_freq"] < unreach
    tab["rank"] = np.nan
    for a, idx in tab.groupby("arm").groups.items():
        s = tab.loc[idx]
        tab.loc[idx, "rank"] = EN.rank_en(s["selection_freq"].to_numpy(), s["beta"].to_numpy(),
                                          (s["gene"] + "|" + s["feature"]).to_numpy())
    return tab


# ------------------------------------------------------------------------------------------------
# Parallel execution (gene chunks over processes; arrays shared read-only)
# ------------------------------------------------------------------------------------------------

_W = {}


def _worker_init(payload):
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = str(payload.get("blas_threads", 1))
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(payload.get("blas_threads", 1))
    except Exception:  # noqa: BLE001 - optional
        pass
    _W["data"] = payload["data"]
    _W["gf"] = payload["gf"]
    _W["engine"] = None


def _worker_en(args):
    arm, genes, lam, K = args
    data = _W["data"]
    if _W["engine"] is None:
        _W["engine"] = en_engine(data)
    return run_en_genes(data, arm, genes, _W["gf"], engine=_W["engine"], lam=lam, K=K)[:3]  # literal-ok: first three return values


def run_en(data, arms=EN_ARMS, genes=None, U=None, workers: int = 1, lam: dict | None = None, K=None,
           chunk: int | None = None):
    """EN arms over ``genes`` (default data.genes). ``lam`` = {arm: {gene: {'G': α, 'X': α}}} for fixed-λ runs."""
    U = universe(data, genes) if U is None else U
    gf = gene_feats(data, U)
    genes = [g for g in (data.genes if genes is None else genes) if g in gf]
    chunk = chunk or PR.null_params()["chunk_genes"]
    tabs, lts, tims = [], [], []
    if workers <= 1 or len(genes) <= chunk:
        eng = en_engine(data)
        for a in arms:
            t_arm = time.time()
            t, l, tm, _ = run_en_genes(data, a, genes, gf, engine=eng, lam=None if lam is None else lam.get(a), K=K)
            tabs.append(t), lts.append(l), tims.append(tm)
            tims.append({"arm": a, "kind": "arm_wall", "start_epoch": t_arm, "end_epoch": time.time(), "workers": 1})
    else:
        from concurrent.futures import ProcessPoolExecutor
        gi = {g: i for i, g in enumerate(data.genes)}
        order = sorted(genes, key=lambda g: data.rows_ok[:, gi[g]].tobytes())
        for a in arms:
            t_arm = time.time()
            jobs = [(a, order[i:i + chunk], None if lam is None else lam.get(a), K) for i in range(0, len(order), chunk)]
            payload = {"data": data.for_workers(include_obs=TARGET_OF[a] == "obs"), "gf": gf,
                       "blas_threads": PR.null_params()["blas_threads"]}
            with ProcessPoolExecutor(max_workers=workers, initializer=_worker_init, initargs=(payload,)) as ex:
                for t, l, tm in ex.map(_worker_en, jobs):
                    tabs.append(t), lts.append(l), tims.append(tm)
            tims.append({"arm": a, "kind": "arm_wall", "start_epoch": t_arm, "end_epoch": time.time(), "workers": workers})
    tab = finish_en(pd.concat([t for t in tabs if len(t)], ignore_index=True)) if any(len(t) for t in tabs) \
        else pd.DataFrame()
    lt = pd.concat(lts, ignore_index=True) if lts else pd.DataFrame()
    return tab, lt, tims


# ------------------------------------------------------------------------------------------------
# MD-strict and assembly
# ------------------------------------------------------------------------------------------------


def strict(en_tab: pd.DataFrame, s1_obs: pd.DataFrame, arm: str = "EN-pred", p_cut: float | None = None) -> pd.DataFrame:
    """MD-strict (arm EN-pred) or EN-obs-strict: selected candidates of ``arm`` passing the D-stage screen.
    ``p_cut`` given (fixed-threshold null): pass = D-stage p <= p_cut instead of q_sc < q (same-sign rule kept)."""
    sp = PR.stats_params()
    sel = en_tab[(en_tab["arm"] == arm) & en_tab["selected"].astype(bool)]
    cols = ["gene", "feature", "beta", "se", "df", "nu", "p", "n_obs", "n_carriers", "dropped_cov",
            "collinear_with_design", "se_type", "df_rule", "df_nu", "p_nu"]
    d = s1_obs[s1_obs["arm"] == "S1-obs"][[c for c in cols if c in s1_obs]].rename(columns={"p": "p_sc", "p_nu": "p_sc_nu"})
    m = sel[["gene", "feature", "selection_freq", "beta", "maj_sign", "stage"]].rename(
        columns={"beta": "beta_en"}).merge(d, on=["gene", "feature"], how="left")
    m["d_stage_testable"] = m["p_sc"].notna()
    m["q_sc"] = ST.bh(m["p_sc"].to_numpy())
    m["same_sign"] = np.sign(m["beta"]) == m["maj_sign"]
    passed = m["same_sign"] & ~m["collinear_with_design"].fillna(False).astype(bool)
    passed &= (m["q_sc"] < sp["q_sc"]) if p_cut is None else (m["p_sc"] <= p_cut)
    m["pass"] = passed.fillna(False).astype(bool)
    m["list"] = "MD-strict" if arm == "EN-pred" else "EN-obs-strict"
    return m


def md_global(en_tab: pd.DataFrame, s1_obs: pd.DataFrame, arm: str = "EN-pred", p_cut: float | None = None) -> pd.DataFrame:
    """MD-global (amendment A7, the catalogue): EN-pred selected ∩ S1-obs F1 BH hits (full-family BH at s1_s5.bh_q,
    not collinear_with_design), same sign (S1-obs β vs the EN majority sign). ``p_cut`` given (fixed-threshold
    null): the S1-obs condition is F1 and p <= p_cut (the observed p*_F1). One row per EN-pred selected candidate
    with ``pass`` = in MD-global; ``list`` = 'MD-global'."""
    sel = en_tab[(en_tab["arm"] == arm) & en_tab["selected"].astype(bool)]
    cols = ["gene", "feature", "beta", "se", "df", "p", "q", "test_family", "selected", "n_obs", "n_carriers",
            "dropped_cov", "collinear_with_design", "se_type", "df_rule", "p_nu"]
    s = s1_obs[s1_obs["arm"] == "S1-obs"]
    d = s[[c for c in cols if c in s]].rename(columns={"selected": "s1_obs_selected"})
    m = sel[["gene", "feature", "selection_freq", "beta", "maj_sign", "stage"]].rename(
        columns={"beta": "beta_en"}).merge(d, on=["gene", "feature"], how="left")
    m["same_sign"] = np.sign(m["beta"]) == m["maj_sign"]
    f1 = (m["test_family"] == "F1") & ~m["collinear_with_design"].fillna(False).astype(bool)
    hit = m["s1_obs_selected"].fillna(False).astype(bool) if p_cut is None else (m["p"] <= p_cut)
    m["s1_f1_hit"] = (f1 & hit).fillna(False).astype(bool)
    m["pass"] = (m["same_sign"] & m["s1_f1_hit"]).fillna(False).astype(bool)
    m["list"] = "MD-global" if arm == "EN-pred" else "EN-obs-global"
    return m


def p_cut_sc(m: pd.DataFrame) -> float:
    return ST.bh_cutoff(m["p_sc"].to_numpy(), PR.stats_params()["q_sc"])


def assemble(data, s1: pd.DataFrame, en_tab: pd.DataFrame) -> pd.DataFrame:
    """Long arms table on the common universe with design columns (DESIGN §5.1)."""
    fm = data.fmeta
    parts = []
    if len(s1):
        parts.append(s1)
    if len(en_tab):
        u = s1[s1["arm"] == s1["arm"].iloc[0]][["gene", "feature", "n_obs", "n_carriers", "n_NA",
                                                "collinear_with_design", "dropped_cov"]] if len(s1) else None
        e = en_tab if u is None else en_tab.merge(u, on=["gene", "feature"], how="left")
        parts.append(e)
    a = pd.concat(parts, ignore_index=True, sort=False)
    a["family"] = fm["family"].reindex(a["feature"]).to_numpy()
    a["d_type"] = fm["d_type"].reindex(a["feature"]).to_numpy()
    a["variant_group"] = fm["variant_group"].reindex(a["feature"]).to_numpy()
    a["feature_stage"] = fm["stage"].reindex(a["feature"]).to_numpy()
    a["access_class"] = [data.access(g, f) for g, f in zip(a["gene"], a["feature"])]
    if len(en_tab):
        mx = en_tab.drop_duplicates(["gene", "feature"]).set_index(["gene", "feature"])["en_max_attainable_freq"]
        a["en_max_attainable_freq"] = mx.reindex(pd.MultiIndex.from_arrays([a["gene"], a["feature"]])).to_numpy()
    a["selected"] = a["selected"].fillna(False).astype(bool)
    a["se_type"] = PR.pp()["se_type"]
    return add_dedup(a)


def add_dedup(a: pd.DataFrame) -> pd.DataFrame:
    """Within (arm, gene, variant group) only the best-ranked member keeps a de-duplicated rank (DESIGN §5.1);
    rank_dedup re-ranks the kept rows 1..m; ``lost_to_group_member`` names the member that kept the slot."""
    a = a.copy()
    a["rank_dedup"] = np.nan
    a["lost_to_group_member"] = ""
    for arm, idx in a.groupby("arm").groups.items():
        s = a.loc[idx]
        s = s[s["rank"].notna()].sort_values("rank")
        key = np.where(s["variant_group"].notna() & (s["variant_group"].astype(str) != ""),
                       s["gene"] + "|" + s["variant_group"].astype(str), s["gene"] + "|" + s["feature"])
        first = ~pd.Series(key, index=s.index).duplicated()
        best = pd.Series(s["feature"].to_numpy(), index=key)[~pd.Index(key).duplicated()]
        kept = s.index[first.to_numpy()]
        a.loc[kept, "rank_dedup"] = np.arange(1, len(kept) + 1)
        lost = s.index[~first.to_numpy()]
        a.loc[lost, "lost_to_group_member"] = best.reindex(np.asarray(key)[~first.to_numpy()]).to_numpy()
    return a


def rank_of_feature(a: pd.DataFrame, arm: str, gene: str, feature: str) -> float:
    """§10.5: 1 + number of de-duplicated rows of ``arm`` outside (g, variant group of d) with a better rank;
    NaN when (g, d) is unrankable in the arm."""
    s = a[a["arm"] == arm]
    me = s[(s["gene"] == gene) & (s["feature"] == feature)]
    if me.empty or pd.isna(me["rank"].iloc[0]):
        return np.nan
    r = me["rank"].iloc[0]
    vg = me["variant_group"].iloc[0]
    same = (s["gene"] == gene) & ((s["feature"] == feature) | ((s["variant_group"] == vg) if pd.notna(vg) and vg
                                                              else False))
    better = s[(~same) & s["rank_dedup"].notna() & (s["rank"] < r)]
    return float(1 + len(better))


def run_arms(data, arms=ARMS, workers: int = 1, genes=None, timings: dict | None = None) -> dict:
    """All requested arms + MD-strict / EN-obs-strict. Returns dict(arms, s1, en, lam, md_strict, en_obs_strict,
    universe_hash, p_cut)."""
    timings = {} if timings is None else timings
    arms = tuple(arms)
    s1_arms = tuple(a for a in S1_ARMS if a in arms)
    if "EN-pred" in arms and "S1-obs" not in s1_arms:
        s1_arms = s1_arms + ("S1-obs",)              # D-stage screen of MD-strict
    t0 = time.perf_counter()
    if s1_arms:
        s1 = run_s1(data, s1_arms, genes=genes, workers=workers)
    else:
        s1 = pd.DataFrame()
    timings["s1_obs_pred"] = time.perf_counter() - t0
    timings["s1_window"] = (time.time() - timings["s1_obs_pred"], time.time())
    U = (s1[s1["arm"] == s1["arm"].iloc[0]][["gene", "feature"]] if len(s1) else universe(data, genes))
    en_arms = tuple(a for a in EN_ARMS if a in arms)
    t0 = time.perf_counter()
    en_tab, lam, tims = (run_en(data, en_arms, genes=genes, U=U, workers=workers) if en_arms
                         else (pd.DataFrame(), pd.DataFrame(), []))
    timings["en_obs_pred"] = time.perf_counter() - t0
    timings["en_detail"] = tims
    s1_out = s1[s1["arm"].isin([a for a in arms if a in S1_ARMS])] if len(s1) else s1
    a = assemble(data, s1_out, en_tab)
    out = {"arms": a, "s1": s1, "en": en_tab, "lam": lam, "universe": U}
    if "EN-pred" in arms:
        md = strict(en_tab, s1, "EN-pred")
        out["md_strict_all"] = md
        out["md_strict"] = md[md["pass"]].reset_index(drop=True)
        out["p_cut_sc"] = p_cut_sc(md)
        mg = md_global(en_tab, s1, "EN-pred")                # A7: the catalogue
        out["md_global_all"] = mg
        out["md_global"] = mg[mg["pass"]].reset_index(drop=True)
        mds = set(zip(out["md_strict"]["gene"], out["md_strict"]["feature"]))
        out["md_global"]["in_md_strict"] = [(g, f) in mds for g, f in zip(out["md_global"]["gene"], out["md_global"]["feature"])]
    if "EN-obs" in arms and len(s1) and (s1["arm"] == "S1-obs").any():
        mo = strict(en_tab, s1, "EN-obs")
        out["en_obs_strict_all"] = mo
        out["en_obs_strict"] = mo[mo["pass"]].reset_index(drop=True)
    out["p_cut"] = observed_cutoffs(s1)
    return out


def observed_cutoffs(s1: pd.DataFrame) -> dict:
    """Observed BH p-cutoffs per (arm, test family): p*_F1 (S1-obs), the S1-pred cutoff, and F2 alike."""
    cut = {}
    if not len(s1):
        return cut
    for (a, tf), sub in s1.groupby(["arm", "test_family"]):
        cut[f"{a}|{tf}"] = ST.bh_cutoff(sub["p"].to_numpy(), PR.stats_params()["bh_q"])
    return cut
