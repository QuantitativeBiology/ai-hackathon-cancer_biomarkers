"""Real-data I/O for the endpoint scripts (run_gate.py, build_p1_pool.py, replicate.py): frozen D-side tables,
dictionary values / covariates / own CN / EXPGATE for any line set (D, R, set A, pseudo-R), the gene-level D
residual SD of §6.1, and line-subsetting of ``pipeline.Data`` for the pseudo-R dry run. Effects are read only by
the scripts, through ``guard.load_effects``."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from . import params as PR
from . import paths as P
from . import pipeline as PL
from . import stats as ST


ARMS_DIR = P.CACHE / "benchmark" / "arms"
XLIB_ARMS_DIR = P.CACHE / "benchmark" / "xlib_arms"


def read_arms(columns=None, arms=None, genes=None, features=None, selected_only: bool = False, dataset_dir=None,
              categorical: bool = False) -> pd.DataFrame:
    """Reader of the streamed arms dataset (``<dataset_dir>/part-*.parquet``, written by blocked.run_streamed).

    Reads part by part with pyarrow row filters, so only the requested ``columns`` (None = all) and rows are
    materialised: ``arms`` (arm names), ``genes`` (Entrez strings), ``features`` (feature ids), ``selected_only``
    (rows with selected == True). String columns are stored as categoricals; they are returned as plain strings
    (shared objects, 8 bytes per row) unless ``categorical`` is True. The full table is ~26 M rows: always pass
    columns (and a row filter where possible)."""
    d = Path(dataset_dir or ARMS_DIR)
    parts = sorted(d.glob("part-*.parquet"))
    if not parts:
        raise FileNotFoundError(f"no arms parts under {d} (run scripts/run_arms.py)")
    filt = []
    if arms is not None:
        filt.append(("arm", "in", [str(a) for a in arms]))
    if genes is not None:
        filt.append(("gene", "in", [str(g) for g in genes]))
    if features is not None:
        filt.append(("feature", "in", [str(f) for f in features]))
    if selected_only:
        filt.append(("selected", "==", True))
    cols = None if columns is None else list(dict.fromkeys(columns))
    out = []
    for p in parts:
        df = pd.read_parquet(p, columns=cols, filters=filt or None)
        if len(df):
            out.append(df)
    if not out:
        return pd.read_parquet(parts[0], columns=cols).iloc[:0]
    df = pd.concat(out, ignore_index=True)
    if not categorical:
        for c in df.columns:
            if isinstance(df[c].dtype, pd.CategoricalDtype):
                df[c] = df[c].astype(object)
    return df


RANK_COLS = ["arm", "gene", "feature", "rank", "rank_stat", "variant_group", "beta", "maj_sign"]
S1O_COLS = ["arm", "gene", "feature", "beta", "se", "sd_dtilde", "dropped_cov", "test_family", "access_class",
            "n_carriers", "selected"]


def pool_inputs(dataset_dir=None, genes=None) -> tuple:
    """The two slim tables the P1 pool / lists need (instead of the whole arms table): ranks of every arm
    (RANK_COLS) and the S1-obs rows (S1O_COLS) used for EB shrinkage, drops, access class and testability."""
    a = read_arms(columns=RANK_COLS, genes=genes, dataset_dir=dataset_dir)
    s1o = read_arms(columns=S1O_COLS, arms=["S1-obs"], genes=genes, dataset_dir=dataset_dir)
    return a, s1o


def arms_table(columns=None, **kw) -> pd.DataFrame:
    """Back-compatible name for ``read_arms`` (pass ``columns``; the full table does not fit in memory)."""
    return read_arms(columns=columns, **kw)


def s1_table() -> pd.DataFrame:
    s = pd.read_parquet(P.CACHE / "benchmark" / "s1.parquet")
    s["gene"] = s["gene"].astype(str)
    return s


UNIVERSES = ("eligible", "ctrl")


def _sfx(universe: str) -> str:
    if universe not in UNIVERSES:
        raise ValueError(f"universe must be one of {UNIVERSES}")
    return "" if universe == "eligible" else f"_{universe}"


def arms_dir(universe: str = "eligible") -> Path:
    """Arms dataset of a universe (resolved at call time): cache/benchmark/arms (eligible) or .../arms_ctrl (U_ctrl,
    scripts/run_arms.py --universe ctrl)."""
    return P.CACHE / "benchmark" / f"arms{_sfx(universe)}"


def md_dir(universe: str = "eligible") -> Path:
    """results/md (eligible) or results/md_ctrl (U_ctrl; BH families within U_ctrl)."""
    return P.RESULTS / f"md{_sfx(universe)}"


def benchmark_dir(universe: str = "eligible") -> Path:
    return P.RESULTS / f"benchmark{_sfx(universe)}"


def universe_genes() -> dict:
    """{'eligible': set, 'ctrl': set} of Entrez strings from configs/eligible_genes.csv and configs/u_ctrl_genes.txt
    (empty if that file is absent). The two sets are disjoint (asserted)."""
    el = pd.read_csv(P.CONFIGS / "eligible_genes.csv", dtype={"entrez": str})
    elig = set(el.loc[el["eligible"].astype(str).str.lower().isin(["true", "1"]), "entrez"])
    f = P.CONFIGS / "u_ctrl_genes.txt"
    ctrl = {x.strip() for x in f.read_text(encoding="utf-8").splitlines() if x.strip()} if f.exists() else set()
    both = elig & ctrl
    if both:
        raise ValueError(f"u_ctrl_genes.txt overlaps the eligible set: {sorted(both)}")
    return {"eligible": elig, "ctrl": ctrl}


def universe_of(genes, ug: dict | None = None) -> dict:
    """gene -> 'eligible' | 'ctrl' | 'none'."""
    ug = universe_genes() if ug is None else ug
    return {str(g): ("eligible" if str(g) in ug["eligible"] else "ctrl" if str(g) in ug["ctrl"] else "none")
            for g in genes}


def md_strict(universe: str = "eligible") -> pd.DataFrame:
    return pd.read_csv(md_dir(universe) / "md_strict.csv", dtype={"gene": str})


def md_global(universe: str = "eligible") -> pd.DataFrame:
    """The catalogue (amendment A7): results/md/md_global.csv (U_ctrl: results/md_ctrl/md_global.csv)."""
    return pd.read_csv(md_dir(universe) / "md_global.csv", dtype={"gene": str})


def catalogue() -> pd.DataFrame:
    """Default catalogue of every downstream analysis (A7) = MD-global."""
    return md_global()


def en_obs_strict() -> pd.DataFrame:
    f = P.RESULTS / "md" / "en_obs_strict.csv"
    return pd.read_csv(f, dtype={"gene": str}) if f.exists() else pd.DataFrame(columns=["gene", "feature"])


def feature_meta() -> pd.DataFrame:
    from . import dictionary as DCT
    fl = pd.DataFrame(DCT.load_meta()["features"]).set_index("id")
    return fl


def line_values(ids, features) -> pd.DataFrame:
    """Dictionary feature values (ids x features) from cache/dictionary (population rows: D, R, set A)."""
    from . import dictionary as DCT
    fl = feature_meta()
    feats = [f for f in dict.fromkeys(features) if f in fl.index]
    parts = []
    for fam, sub in fl.loc[feats].groupby("family"):
        parts.append(DCT.load_family(fam, columns=list(sub.index)).reindex(list(ids)))
    X = pd.concat(parts, axis=1) if parts else pd.DataFrame(index=list(ids))
    return X.reindex(columns=feats)


def line_gene_covs(ids, genes):
    from . import dictionary as DCT
    oc = pd.read_parquet(DCT.CACHE_DIR / "own_cn_log2.parquet").reindex(index=list(ids), columns=list(genes))
    eg = pd.read_parquet(DCT.CACHE_DIR / "expgate.parquet").reindex(index=list(ids), columns=list(genes))
    return oc, eg


def line_covariates(ids) -> pd.DataFrame:
    return PL.real_covariates(pd.Index(list(ids), name="ModelID"))


def library_strings(ids) -> pd.Series:
    c = pd.read_parquet(P.CACHE / "covariates.parquet")
    return c["Library_raw"].reindex(list(ids)).fillna("")


def sigma_resid_genes(Y: np.ndarray, rows_ok: np.ndarray, cov: ST.CovDesign, genes) -> pd.Series:
    """σ̂_resid,D(g) (§6.1): residual SD of y_g on [1, design blocks, droppable covariates] over g's non-NaN rows."""
    out = {}
    names = cov.design_names + cov.droppable_names
    for pat in {rows_ok[:, j].tobytes() for j in range(rows_ok.shape[1])}:
        rows = np.frombuffer(pat, dtype=bool)
        cols = [j for j in range(rows_ok.shape[1]) if rows_ok[:, j].tobytes() == pat]
        Q = ST.orth(cov.Z(rows, names))
        R = ST.resid(Q, Y[rows][:, cols])
        dfr = max(int(rows.sum()) - Q.shape[1], 1)
        for j, c in zip(cols, np.sqrt((R ** 2).sum(0) / dfr)):
            out[str(genes[j])] = float(c)
    return pd.Series(out)


def subset_data(d: PL.Data, keep: np.ndarray) -> PL.Data:
    """``Data`` restricted to the lines in boolean ``keep`` (pseudo-R dry run: D' = folds != 0)."""
    keep = np.asarray(keep, dtype=bool)
    obs_full = d.obs_loader

    def obs():
        return np.asarray(obs_full(), dtype=np.float64)[keep]

    return PL.Data(lines=d.lines[keep], genes=list(d.genes), X=d.X[keep], fmeta=d.fmeta, Dg=d.Dg, forced=d.forced,
                   cov_df=d.cov_df.loc[keep], na_flags=d.na_flags[keep],
                   own_cn=None if d.own_cn is None else d.own_cn[keep],
                   expgate=None if d.expgate is None else d.expgate[keep], rows_ok=d.rows_ok[keep],
                   pred=d.pred[keep], folds=d.folds[keep], sub_masks=[m[keep] for m in d.sub_masks],
                   patients=d.patients[keep], access_map=d.access_map, access_fn=d.access_fn,
                   tp53_genes=d.tp53_genes, mdm2_genes=d.mdm2_genes, obs_loader=obs, extras=dict(d.extras))


def read_json(path) -> dict:
    return json.loads(open(path, encoding="utf-8").read())


def write_json(path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1, default=_js), encoding="utf-8")


def _js(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if np.isnan(o) else float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


# ------------------------------------------------------------------------------------------------
# P1 pool and lists (pre-tag; R molecular data only)
# ------------------------------------------------------------------------------------------------


EB_CHUNK_ROWS = 262144                     # literal-ok: row chunk of the low-memory EB (memory only, never results)


def shrunk_beta_lowmem(beta, se, sd_dtilde, chunk: int = EB_CHUNK_ROWS) -> np.ndarray:
    """``stats.shrunk_beta`` (§6.1 EB normal scale mixture; unchanged grid, start, EM, tolerance and cap) with the
    n x K temporaries built in row chunks: only the likelihood matrix L (n x K) is held. Every row-wise operation is
    the same elementwise arithmetic; the EM column sum accumulates the rows in the same order as the one-shot
    ``w.mean(0)`` (previous partial sum prepended to each chunk), so the result is bit-identical
    (tests/test_pool_fast.py). Peak memory ~ n x K x 8 bytes instead of ~5 such arrays."""
    rp = PR.replication_params()
    tiny = PR.stats_params()["tiny"]
    sd = np.asarray(sd_dtilde, float)
    b, s = np.asarray(beta, float) * sd, np.asarray(se, float) * sd
    ok = np.isfinite(b) & np.isfinite(s) & (s > 0)
    bo, so = b[ok], s[ok]
    lo = rp["eb_sigma_min_frac"] * so.min()
    hi = rp["eb_sigma_max_mult"] * np.sqrt(max(np.max(bo ** 2 - so ** 2), (2 * lo) ** 2))
    sig = np.geomspace(lo, max(hi, lo * 2), rp["eb_n_sigma"])
    n, K = bo.size, sig.size
    L = np.empty((n, K))
    for i in range(0, n, chunk):
        var = sig[None, :] ** 2 + so[i:i + chunk, None] ** 2
        L[i:i + chunk] = np.exp(-0.5 * bo[i:i + chunk, None] ** 2 / var) / np.sqrt(var)  # literal-ok: as stats
    pi = np.full(K, 1.0 / K)
    buf = np.empty((min(chunk, n) + 1, K))

    def weights(i, p):
        m = min(chunk, n - i)
        w = np.multiply(L[i:i + m], p, out=buf[1:m + 1])
        w /= np.maximum(w.sum(1, keepdims=True), tiny)
        return w, m

    for _ in range(rp["eb_max_iter"]):
        acc = None
        for i in range(0, n, chunk):
            w, m = weights(i, pi)
            if acc is None:
                acc = np.add.reduce(w, axis=0)
            else:
                buf[0] = acc                                   # same left-to-right row order as w.sum(0)
                acc = np.add.reduce(buf[:m + 1], axis=0)
        new = acc / n
        done = np.max(np.abs(new - pi)) < rp["eb_tol"]
        pi = new
        if done:
            break
    shrink = np.empty(n)
    for i in range(0, n, chunk):
        w, m = weights(i, pi)
        var = sig[None, :] ** 2 + so[i:i + m, None] ** 2
        shrink[i:i + m] = (w * (sig[None, :] ** 2 / var)).sum(1)
    post = np.full(b.shape, np.nan)
    post[ok] = shrink * bo
    with np.errstate(divide="ignore", invalid="ignore"):
        return post / sd


def make_pool(arms_tab: pd.DataFrame, s1_obs: pd.DataFrame, sigma_g: pd.Series, D_ids, R_ids,
              nan_D: pd.DataFrame, arms=None, fast: bool = True, workers: int = 1) -> tuple:
    """Pool_G and frozen top-k lists (DESIGN §6.1). ``s1_obs``: S1-obs rows (all U incl. F2; beta, se, sd_dtilde,
    dropped_cov, test_family) used for the EB shrinkage and the drops; ``arms_tab``: ranks of every arm on the same U.
    ``arms``: the arms to list (default the four benchmark arms; XLIB §6.3: S1-obs and EN-obs only).
    ``fast`` (default): low-memory EB (shrunk_beta_lowmem), cached design blocks / clusters in pool_power over
    ``workers`` processes, and the arms rows restricted to P_G once before the lists; results are identical to
    ``fast=False`` (tests/test_pool_fast.py).
    Returns (pool DataFrame over U_G, lists dict for p1_pool.json)."""
    from . import arms as A
    from . import endpoints as EP
    from . import params as PR
    e = PR.endpoints()["P1_primary"]
    s = s1_obs.copy()
    s["gene"] = s["gene"].astype(str)
    ok = s["beta"].notna() & s["se"].notna() & (s["sd_dtilde"] > 0)
    shrink = shrunk_beta_lowmem if fast else ST.shrunk_beta          # bit-identical; fast = low memory
    bsh = pd.Series(shrink(s.loc[ok, "beta"].to_numpy(), s.loc[ok, "se"].to_numpy(), s.loc[ok, "sd_dtilde"].to_numpy()),
                    index=pd.MultiIndex.from_arrays([s.loc[ok, "gene"], s.loc[ok, "feature"]]))
    fl = feature_meta()
    stage = fl["stage"].reindex(s["feature"]).to_numpy()
    UG = s[(stage == "G") & (s["test_family"] == "F1")][["gene", "feature"]].drop_duplicates()
    from . import endpoints as EP
    feats = list(dict.fromkeys(UG["feature"]))
    XR = line_values(R_ids, feats)
    fpos = pd.Series(np.arange(len(feats)), index=feats)
    covR_df = line_covariates(R_ids)
    ocR, egR = line_gene_covs(R_ids, sorted(UG["gene"].unique()))
    covR = EP.RDesign(covR_df, ocR, egR)
    patR = covR_df["PatientID"].astype(str).to_numpy()
    genes = sorted(UG["gene"].unique())
    expR = EP.expected_nonnan(nan_D.reindex(index=list(D_ids), columns=genes), library_strings(D_ids),
                              library_strings(R_ids))
    drops = s.set_index(["gene", "feature"])["dropped_cov"]
    binary = fl["d_type"].reindex(feats) == "binary"
    pool = EP.pool_power(UG, bsh, sigma_g, XR.to_numpy(dtype=np.float64), fpos, binary, covR, patR, expR, drops,
                         workers=workers, fast=fast)
    src = s1_obs if "access_class" in s1_obs else arms_tab[arms_tab["arm"] == "S1-obs"]
    acc = src.drop_duplicates(["gene", "feature"]).set_index(["gene", "feature"])["access_class"]
    pool["access_class"] = acc.reindex(pd.MultiIndex.from_arrays([pool["gene"], pool["feature"]])).to_numpy()
    P_G = pool[pool["powered"]]
    out = {"k_rule": e["k_rule"], "pool_size": int(len(P_G)), "k": EP.p1_k(len(P_G)), "lists": {}, "strata": {},
           "secondary_k": {}}
    arms = tuple(A.ARMS if arms is None else arms)
    if fast:     # every list pool is a subset of P_G: rows outside it never enter a list (identical lists, less work)
        pk = set(zip(P_G["gene"].astype(str), P_G["feature"]))
        keep = np.array([(str(g), f) in pk for g, f in zip(arms_tab["gene"], arms_tab["feature"])], dtype=bool)
        arms_tab = arms_tab[keep]
    out["lists"]["G"] = EP.build_lists(arms_tab, P_G, arms)
    for name, m in EP.strata_of(P_G["access_class"]).items():
        sub = P_G[m.to_numpy()]
        out["strata"][name] = {"pool_size": int(len(sub)), "k": EP.p1_k(len(sub))}
        out["lists"][name] = EP.build_lists(arms_tab, sub, arms)
    for k2 in e["secondary_k"]:
        out["secondary_k"][str(k2)] = EP.build_lists(arms_tab, P_G, arms, k=int(k2))
    return pool, out


REQUEST_COLS = ["gene", "feature", "reasons", "beta_D", "se_D", "dropped_cov", "sign_D", "expected_sign", "power",
                "family", "access_class"]


def build_requests(pool: pd.DataFrame, lists: dict, s1_obs: pd.DataFrame, md: pd.DataFrame, eos: pd.DataFrame,
                   s1_hits: pd.DataFrame, controls: pd.DataFrame | None, cgt_tests: pd.DataFrame | None,
                   extra: pd.DataFrame | None = None, mg: pd.DataFrame | None = None) -> pd.DataFrame:
    """configs/catalogue/r_requests.csv rows: every (g, d) needing R statistics, with '|'-joined reasons."""
    rows = []

    def add(g, f, why, exp=""):
        rows.append((str(g), f, why, exp))

    for stratum, L in lists["lists"].items():
        for arm, v in L.items():
            for g, f, _ in v["rows"]:
                add(g, f, f"P1:{stratum}:{arm}")
    for kk, L in lists.get("secondary_k", {}).items():
        for arm, v in L.items():
            for g, f, _ in v["rows"]:
                add(g, f, f"P1:k{kk}:{arm}")
    if mg is not None:
        for g, f in zip(mg.get("gene", []), mg.get("feature", [])):
            add(g, f, "MD-global")
    for g, f in zip(md.get("gene", []), md.get("feature", [])):
        add(g, f, "MD-strict")
    for g, f in zip(eos.get("gene", []), eos.get("feature", [])):
        add(g, f, "EN-obs-strict")
    for g, f in zip(s1_hits.get("gene", []), s1_hits.get("feature", [])):
        add(g, f, "S1-obs BH")
    if controls is not None:
        from . import endpoints as EP
        for _, c in controls.iterrows():
            for f in EP.control_features(c):
                add(c["dependency_entrez"], f, "control", c.get("expected_sign", ""))
    if cgt_tests is not None:
        for _, t in cgt_tests.iterrows():
            for g in t["targets"]:
                if isinstance(t["feature_D"], str):
                    add(g, t["feature_D"], f"CGT:{t['test_id']}", t["sign"])
    if extra is not None:
        for g, f, w in zip(extra["gene"], extra["feature"], extra["reason"]):
            add(g, f, w)
    R = pd.DataFrame(rows, columns=["gene", "feature", "reason", "expected_sign"])
    U = s1_obs.assign(gene=s1_obs["gene"].astype(str)).drop_duplicates(["gene", "feature"]).set_index(["gene", "feature"])
    R = R[np.array([k in U.index for k in zip(R["gene"], R["feature"])], dtype=bool)]
    agg = R.groupby(["gene", "feature"], sort=False).agg(
        reasons=("reason", lambda v: "|".join(dict.fromkeys(v))),
        expected_sign=("expected_sign", lambda v: next((x for x in v if x in ("S", "W")), ""))).reset_index()
    key = pd.MultiIndex.from_arrays([agg["gene"], agg["feature"]])
    agg["beta_D"] = U["beta"].reindex(key).to_numpy()
    agg["se_D"] = U["se"].reindex(key).to_numpy()
    agg["dropped_cov"] = U["dropped_cov"].reindex(key).fillna("").to_numpy()
    agg["sign_D"] = np.sign(agg["beta_D"])
    pw = pool.set_index(["gene", "feature"])["power"] if "power" in pool else pd.Series(dtype=float)
    agg["power"] = pw.reindex(key).to_numpy()
    fl = feature_meta()
    agg["family"] = fl["family"].reindex(agg["feature"]).to_numpy()
    acc = U["access_class"] if "access_class" in U else pd.Series(dtype=object)
    agg["access_class"] = acc.reindex(key).to_numpy() if len(acc) else ""
    return agg[REQUEST_COLS]


def fitter_for(ids, genes, features, Y: pd.DataFrame):
    """endpoints.RFitter on any line set (values, covariates, own CN, EXPGATE from the dictionary caches)."""
    from . import endpoints as EP
    X = line_values(ids, features)
    oc, eg = line_gene_covs(ids, genes)
    cov = line_covariates(ids)
    binary = feature_meta()["d_type"].reindex(X.columns) == "binary"
    return EP.RFitter(X.to_numpy(dtype=np.float64), pd.Series(np.arange(X.shape[1]), index=X.columns), cov, oc, eg,
                      Y.reindex(index=list(ids), columns=list(genes)), cov["PatientID"].astype(str).to_numpy(),
                      binary=binary)
