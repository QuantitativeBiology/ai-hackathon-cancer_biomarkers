"""Real-data I/O for the endpoint scripts (run_gate.py, build_p1_pool.py, replicate.py): frozen D-side tables,
dictionary values / covariates / own CN / EXPGATE for any line set (D, R, set A, pseudo-R), the gene-level D
residual SD of §6.1, and line-subsetting of ``pipeline.Data`` for the pseudo-R dry run. Effects are read only by
the scripts, through ``guard.load_effects``."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from . import paths as P
from . import pipeline as PL
from . import stats as ST


def arms_table() -> pd.DataFrame:
    a = pd.read_parquet(P.CACHE / "benchmark" / "arms.parquet")
    a["gene"] = a["gene"].astype(str)
    return a


def s1_table() -> pd.DataFrame:
    s = pd.read_parquet(P.CACHE / "benchmark" / "s1.parquet")
    s["gene"] = s["gene"].astype(str)
    return s


def md_strict() -> pd.DataFrame:
    return pd.read_csv(P.RESULTS / "md" / "md_strict.csv", dtype={"gene": str})


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


def make_pool(arms_tab: pd.DataFrame, s1_obs: pd.DataFrame, sigma_g: pd.Series, D_ids, R_ids,
              nan_D: pd.DataFrame) -> tuple:
    """Pool_G and frozen top-k lists (DESIGN §6.1). ``s1_obs``: S1-obs rows (all U incl. F2; beta, se, sd_dtilde,
    dropped_cov, test_family) used for the EB shrinkage and the drops; ``arms_tab``: ranks of every arm on the same U.
    Returns (pool DataFrame over U_G, lists dict for p1_pool.json)."""
    from . import arms as A
    from . import endpoints as EP
    from . import params as PR
    e = PR.endpoints()["P1_primary"]
    s = s1_obs.copy()
    s["gene"] = s["gene"].astype(str)
    ok = s["beta"].notna() & s["se"].notna() & (s["sd_dtilde"] > 0)
    bsh = pd.Series(ST.shrunk_beta(s.loc[ok, "beta"].to_numpy(), s.loc[ok, "se"].to_numpy(),
                                   s.loc[ok, "sd_dtilde"].to_numpy()),
                    index=pd.MultiIndex.from_arrays([s.loc[ok, "gene"], s.loc[ok, "feature"]]))
    fl = feature_meta()
    stage = fl["stage"].reindex(s["feature"]).to_numpy()
    UG = s[(stage == "G") & (s["test_family"] == "F1")][["gene", "feature"]].drop_duplicates()
    feats = list(dict.fromkeys(UG["feature"]))
    XR = line_values(R_ids, feats)
    fpos = pd.Series(np.arange(len(feats)), index=feats)
    covR = ST.CovDesign(line_covariates(R_ids))
    patR = line_covariates(R_ids)["PatientID"].astype(str).to_numpy()
    genes = sorted(UG["gene"].unique())
    expR = EP.expected_nonnan(nan_D.reindex(index=list(D_ids), columns=genes), library_strings(D_ids),
                              library_strings(R_ids))
    drops = s.set_index(["gene", "feature"])["dropped_cov"]
    binary = fl["d_type"].reindex(feats) == "binary"
    pool = EP.pool_power(UG, bsh, sigma_g, XR.to_numpy(dtype=np.float64), fpos, binary, covR, patR, expR, drops)
    acc = arms_tab[arms_tab["arm"] == "S1-obs"].drop_duplicates(["gene", "feature"]).set_index(["gene", "feature"])["access_class"]
    pool["access_class"] = acc.reindex(pd.MultiIndex.from_arrays([pool["gene"], pool["feature"]])).to_numpy()
    P_G = pool[pool["powered"]]
    out = {"k_rule": e["k_rule"], "pool_size": int(len(P_G)), "k": EP.p1_k(len(P_G)), "lists": {}, "strata": {},
           "secondary_k": {}}
    out["lists"]["G"] = EP.build_lists(arms_tab, P_G, A.ARMS)
    for name, m in EP.strata_of(P_G["access_class"]).items():
        sub = P_G[m.to_numpy()]
        out["strata"][name] = {"pool_size": int(len(sub)), "k": EP.p1_k(len(sub))}
        out["lists"][name] = EP.build_lists(arms_tab, sub, A.ARMS)
    for k2 in e["secondary_k"]:
        out["secondary_k"][str(k2)] = EP.build_lists(arms_tab, P_G, A.ARMS, k=int(k2))
    return pool, out


REQUEST_COLS = ["gene", "feature", "reasons", "beta_D", "se_D", "dropped_cov", "sign_D", "expected_sign", "power",
                "family", "access_class"]


def build_requests(pool: pd.DataFrame, lists: dict, s1_obs: pd.DataFrame, md: pd.DataFrame, eos: pd.DataFrame,
                   s1_hits: pd.DataFrame, controls: pd.DataFrame | None, cgt_tests: pd.DataFrame | None,
                   extra: pd.DataFrame | None = None) -> pd.DataFrame:
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
    R = R[[k in U.index for k in zip(R["gene"], R["feature"])]]
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
