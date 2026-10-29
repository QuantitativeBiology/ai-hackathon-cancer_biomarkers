"""Endpoint code: positive-control gate (DESIGN v3.2 §9.3), P1 replication (§6.1, primary_endpoints.P1_primary),
P2 CGT recall (§10.5, primary_endpoints.P2_primary_clinical, amendment A4), direction / actionability / population
tables (§10.5 Q4a-c) and the verdict branches. Pure functions on tables; the scripts (run_gate.py, build_p1_pool.py,
replicate.py) do the I/O and the guard checks. Nothing here reads effects: callers pass them in.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as sps

from . import gitutil
from . import params as PR
from . import paths as P
from . import stats as ST


class EndpointRefusal(PermissionError):
    """An endpoint script was asked to run without its preconditions (tag ancestry, single R read, hashes)."""


# ------------------------------------------------------------------------------------------------
# Preconditions
# ------------------------------------------------------------------------------------------------


def descends_from(tag: str, repo=None) -> bool:
    """True iff tag exists and is an ancestor of (or equal to) HEAD."""
    repo = repo or P.PROJECT
    if not gitutil.tag_exists(repo, tag):
        return False
    return gitutil.git(repo, "merge-base", "--is-ancestor", f"refs/tags/{tag}", "HEAD").returncode == 0


def require_descends(tag: str, repo=None, check=None) -> None:
    check = descends_from if check is None else check       # looked up at call time (patchable in tests)
    if not check(tag, repo):
        raise EndpointRefusal(f"HEAD does not descend from tag {tag!r}")


def prior_r_reads(log: pd.DataFrame) -> pd.DataFrame:
    """Granted replication (R) reads in the guard access log."""
    if log is None or log.empty or "scope" not in log:
        return pd.DataFrame(columns=getattr(log, "columns", []))
    return log[log["scope"].astype(str) == "replication"]


def require_single_r_read(log: pd.DataFrame, reproduce_of: dict | None = None) -> None:
    """Refuse a second R read. DESIGN §1.3 (h): afterwards only reproductions (same commit and config hash), which
    the caller must request explicitly with ``reproduce_of`` = {'commit': ..., 'config_sha256': ...}."""
    prev = prior_r_reads(log)
    if prev.empty:
        return
    if reproduce_of is None:
        raise EndpointRefusal(f"R was already read ({len(prev)} granted access(es)); a second R open is refused")
    if not (prev["commit"].astype(str) == str(reproduce_of.get("commit"))).all():
        raise EndpointRefusal("reproduction refused: commit differs from the first R access")
    if "config_sha256" in reproduce_of and "purpose" in prev:
        tag = f"cfg={reproduce_of['config_sha256']}"
        if not prev["purpose"].astype(str).str.contains(tag, regex=False).all():
            raise EndpointRefusal("reproduction refused: config hash differs from the first R access")


def catalogue_mismatches(manifest: dict, root=None) -> list:
    """Files of the frozen catalogue manifest ({'files': {relpath: sha256}}) whose sha256 differs or is missing."""
    import hashlib
    root = root or P.PROJECT
    bad = []
    for rel, sha in manifest.get("files", {}).items():
        f = root / rel
        if not f.exists() or hashlib.sha256(f.read_bytes()).hexdigest() != sha:
            bad.append(rel)
    return bad


def catalogue_manifest(relpaths, root=None) -> dict:
    import hashlib
    root = root or P.PROJECT
    return {"files": {r: hashlib.sha256((root / r).read_bytes()).hexdigest() for r in relpaths}}


# ------------------------------------------------------------------------------------------------
# Exact statistics and verdict branches
# ------------------------------------------------------------------------------------------------


def clopper_pearson(k: int, n: int, level: float | None = None) -> tuple:
    level = float(PR.endpoints()["P1_primary"]["ci"]["level"]) if level is None else level
    if n == 0:
        return (np.nan, np.nan)
    a = (1 - level) / 2
    lo = 0.0 if k == 0 else float(sps.beta.ppf(a, k, n - k + 1))
    hi = 1.0 if k == n else float(sps.beta.ppf(1 - a, k + 1, n - k))
    return (lo, hi)


def mcnemar_exact(b: int, c: int) -> float:
    """Exact two-sided McNemar p from the discordant counts b (arm 1 hit only) and c (arm 2 hit only)."""
    n = int(b + c)
    if n == 0:
        return 1.0
    return float(min(1.0, 2 * sps.binom.cdf(min(b, c), n, 0.5)))


def descriptive(n: int) -> str:
    return f"descriptive (underpowered: {n})"


def p1_branch(lo: float, hi: float, contrast: str) -> str:
    """Branch of a P1 contrast from its 95% CI (primary_endpoints.P1_primary.contrasts.<contrast>.branches)."""
    br = PR.endpoints()["P1_primary"]["contrasts"][contrast]["branches"]
    if lo > 0:
        return br["lower > 0"]
    if hi < 0:
        return br["upper < 0"]
    return br["else"]


def p1_verdict(k: int, lo: float, hi: float, contrast: str) -> str:
    kmin = int(PR.endpoints()["P1_primary"]["min_k"])
    return descriptive(k) if k < kmin else p1_branch(lo, hi, contrast)


def p2_verdict(n_testable: int, p: float, recall_en: float, recall_s1: float) -> str:
    pc = PR.endpoints()["P2_primary_clinical"]
    if n_testable < int(pc["min_testable_tests"]):
        return descriptive(n_testable)
    br = pc["primary_contrast"]["branches"]
    keys = list(br)
    alpha = float(pc["primary_contrast"]["alpha"])
    if p < alpha and recall_en > recall_s1:
        return br[keys[0]]
    if p < alpha and recall_s1 > recall_en:
        return br[keys[1]]
    return br[keys[2]]


# ------------------------------------------------------------------------------------------------
# Gate (§9.3)
# ------------------------------------------------------------------------------------------------


def _ids(v) -> list:
    if not isinstance(v, str) or not v.strip():
        return []
    return [x for x in v.replace("|", ";").split(";") if x]


def control_features(row) -> list:
    """Recoverable feature ids of a control: primary + secondary ids, and the fallback ids when the control's gate
    feature is the fallback (SIL:<x> already expanded by build_positive_controls.py)."""
    f = _ids(row.get("feature_primary_ids")) + _ids(row.get("feature_secondary_ids"))
    if str(row.get("gate_feature", "")) == "fallback":
        f += _ids(row.get("fallback_feature_ids"))
    return list(dict.fromkeys(f))


def _recovered(lst: pd.DataFrame, gene: str, feats: list, sign: str) -> bool:
    if lst is None or len(lst) == 0:
        return False
    m = lst[(lst["gene"].astype(str) == str(gene)) & lst["feature"].isin(feats)]
    if sign in ("S", "W"):
        m = m[(m["beta"] < 0) if sign == "S" else (m["beta"] > 0)]
    return bool(len(m))


def gate_attempt(pc_elig: pd.DataFrame, controls: pd.DataFrame, lists: dict, universe_of: dict | None = None) -> dict:
    """Gate on ``lists`` ({name: DataFrame(gene, feature, beta)}, must contain 'MD-strict'). Denominator = gate-eligible
    rows of pc_eligibility.csv. ``universe_of`` (gene -> 'eligible' / 'ctrl'), when given, is recorded per row: the
    caller reads each control's rows from the universe of its gene (U_ctrl for non-eligible control genes, DESIGN §2).
    Returns {'rows': per-control table, 'summary': dict}."""
    g9 = PR.dp()["gate_s9"]
    unreach = PR.en_params()["unreachable_lt"]
    e = pc_elig[pc_elig["gate_eligible"].astype(str).str.lower().isin(["true", "1"])].copy()
    c = controls.set_index("control_id")
    rows = []
    for _, r in e.iterrows():
        cr = c.loc[r["control_id"]]
        feats = control_features(cr)
        rec = {"control_id": r["control_id"], "family": r["family"], "gene": str(cr["dependency_entrez"]),
               "features": ";".join(feats), "expected_sign": cr["expected_sign"],
               "access_class": r.get("access_class", ""), "en_unreachable": float(r.get("en_max_attainable_freq", 1.0)) < unreach}
        if universe_of is not None:
            rec["universe"] = universe_of.get(rec["gene"], "none")
        for name, L in lists.items():
            rec[f"recovered_{name}"] = _recovered(L, cr["dependency_entrez"], feats, cr["expected_sign"])
        rec["recovered"] = rec["recovered_MD-strict"]
        rows.append(rec)
    T = pd.DataFrame(rows)

    def rate(sub, col="recovered"):
        k, n = int(sub[col].sum()), int(len(sub))
        lo, hi = clopper_pearson(k, n)
        return {"k": k, "n": n, "recall": k / n if n else np.nan, "ci95": [lo, hi]}

    fam = {f: rate(s) for f, s in T.groupby("family")}
    fam_ok = {f: v["recall"] >= g9["family_recall_min"] for f, v in fam.items() if v["n"] >= g9["family_min_eligible_rows"]}
    overall = rate(T)
    summ = {"decided_on": "MD-strict", "overall": overall, "by_family": fam, "family_pass": fam_ok,
            "pass": bool(overall["recall"] >= g9["overall_recall_min"] and all(fam_ok.values())),
            "by_list": {n: rate(T, f"recovered_{n}") for n in lists},
            "by_access_class": {a: rate(s) for a, s in T.groupby("access_class")},
            "PAR_non_direct": rate(T[(T["family"] == "PAR") & (T["access_class"] != "direct")]),
            "en_unreachable": rate(T[T["en_unreachable"]]), "n_denominator": int(len(T)),
            **({"by_universe": {u: rate(s) for u, s in T.groupby("universe")}} if "universe" in T else {}),
            "thresholds": {"overall_recall_min": g9["overall_recall_min"], "family_recall_min": g9["family_recall_min"],
                           "family_min_eligible_rows": g9["family_min_eligible_rows"]}}
    return {"rows": T, "summary": summ}


# ------------------------------------------------------------------------------------------------
# Ranks (fast rank_of_feature, same semantics as arms.rank_of_feature)
# ------------------------------------------------------------------------------------------------


class RankIndex:
    """Per-arm de-duplicated rank arrays for many rank_of_feature_D queries (DESIGN §10.5)."""

    def __init__(self, arms_tab: pd.DataFrame):
        self.by = {}
        for arm, s in arms_tab.groupby("arm"):
            s = s[["gene", "feature", "rank", "rank_dedup", "variant_group", "beta"]]
            key = pd.MultiIndex.from_arrays([s["gene"].astype(str), s["feature"]])
            d = s[s["rank_dedup"].notna()]
            vg = d["variant_group"].where(d["variant_group"].notna() & (d["variant_group"].astype(str) != ""), None)
            self.by[arm] = {"tab": s.set_index(key), "dr": np.sort(d["rank"].to_numpy(dtype=float)),
                            "dedup": d.assign(vg=vg.to_numpy())}

    @classmethod
    def slim(cls, dedup_ranks: dict, rows: pd.DataFrame) -> "RankIndex":
        """Memory-light index for the full real table: ``dedup_ranks`` = {arm: rank values of every de-duplicated
        row} and ``rows`` = all rows (gene, feature, rank, rank_dedup, variant_group, beta, arm) of the genes that
        will be queried. ``rank_of`` is identical to the full constructor for those genes."""
        self = cls.__new__(cls)
        self.by = {}
        for arm, dr in dedup_ranks.items():
            s = rows[rows["arm"] == arm][["gene", "feature", "rank", "rank_dedup", "variant_group", "beta"]]
            key = pd.MultiIndex.from_arrays([s["gene"].astype(str), s["feature"]])
            d = s[s["rank_dedup"].notna()]
            vg = d["variant_group"].where(d["variant_group"].notna() & (d["variant_group"].astype(str) != ""), None)
            self.by[arm] = {"tab": s.set_index(key), "dr": np.sort(np.asarray(dr, dtype=float)),
                            "dedup": d.assign(vg=vg.to_numpy())}
        return self

    def rank_of(self, arm: str, gene: str, feature: str) -> tuple:
        """(rank_of_feature_D, beta) or (nan, nan) when (g, d) is not rankable in the arm."""
        A = self.by.get(arm)
        if A is None or (str(gene), feature) not in A["tab"].index:
            return np.nan, np.nan
        me = A["tab"].loc[(str(gene), feature)]
        if isinstance(me, pd.DataFrame):
            me = me.iloc[0]
        r = me["rank"]
        if pd.isna(r):
            return np.nan, float(me["beta"])
        better = int(np.searchsorted(A["dr"], r, side="left"))
        d = A["dedup"]
        vg = me["variant_group"]
        same = (d["gene"].astype(str) == str(gene)) & ((d["feature"] == feature) |
                                                       ((d["vg"] == vg) if (isinstance(vg, str) and vg) else False))
        better -= int((same & (d["rank"] < r)).sum())
        return float(1 + better), float(me["beta"])


# ------------------------------------------------------------------------------------------------
# P2 (§10.5)
# ------------------------------------------------------------------------------------------------


def _targets(v) -> list:
    return [x.strip() for x in str(v).replace("|", ";").split(";") if x.strip() and x.strip().lower() != "nan"]


def p2_test_set(cgt: pd.DataFrame, U: pd.DataFrame, set_: str = "non_proxy", gate_disjoint: bool = False,
                entrez_of=None) -> pd.DataFrame:
    """Headline-eligible CGT tests (tiers T1-T3, A4; S/W; sensitivity / bypass-resistance; concordant KO classes;
    proxy = 0 for non_proxy; headline_exclude empty) with ``testable`` = feature_D in U for some g in T with >= 5
    carriers. ``entrez_of``: symbol -> entrez map for target sets given as symbols."""
    cs = PR.dp()["clinical_s10"]
    tiers = [t.split()[0] for t in cs["tiers"] if not t.startswith("refuted")]
    conc = [k for k in cs["ko_classes"] if k.startswith("concordant")]
    t = cgt.copy()
    m = t["tier"].isin(tiers) & t["sign"].isin(["S", "W"]) & t["partition"].isin(["sensitivity", "bypass-resistance"])
    m &= t["ko_class"].isin(conc)
    if "headline_exclude" in t:
        m &= t["headline_exclude"].isna() | (t["headline_exclude"].astype(str).str.strip() == "")
    if set_ == "non_proxy":
        m &= pd.to_numeric(t.get("proxy", 0), errors="coerce").fillna(0).astype(int) == 0
    if gate_disjoint:
        m &= pd.to_numeric(t.get("in_gate", 0), errors="coerce").fillna(0).astype(int) == 0
    t = t[m].copy()
    nmin = int(cs["testable_min_carriers"])
    Uc = U.set_index(["gene", "feature"])["n_carriers"] if "n_carriers" in U else None
    tg = []
    for v in t["target_set"]:
        ts = _targets(v)
        if entrez_of is not None:
            ts = [str(entrez_of.get(x, x)) for x in ts]
        tg.append(ts)
    t["targets"] = tg
    ok = []
    for ts, f in zip(t["targets"], t["feature_D"]):
        good = False
        if isinstance(f, str) and f and Uc is not None:
            for g in ts:
                nc = Uc.get((g, f))
                if nc is not None and (pd.isna(nc) or nc >= nmin):
                    good = True
        ok.append(good)
    t["testable"] = ok
    return t


def p2_score(arms_tab: pd.DataFrame, tests: pd.DataFrame, k: int | None = None, arms=None, ridx=None) -> pd.DataFrame:
    """Per (arm, test): rank_of_feature_D (best over g in T), hit (rank <= k, any sign), sign_correct."""
    k = int(PR.endpoints()["P2_primary_clinical"]["k"]) if k is None else int(k)
    ridx = ridx or RankIndex(arms_tab)
    arms = list(ridx.by) if arms is None else list(arms)
    rows = []
    for arm in arms:
        for _, t in tests.iterrows():
            best, beta, gbest = np.nan, np.nan, ""
            for g in t["targets"]:
                r, b = ridx.rank_of(arm, g, t["feature_D"]) if isinstance(t["feature_D"], str) else (np.nan, np.nan)
                if np.isfinite(r) and (not np.isfinite(best) or r < best):
                    best, beta, gbest = r, b, g
            hit = bool(np.isfinite(best) and best <= k)
            sc = (bool(beta < 0) if t["sign"] == "S" else bool(beta > 0)) if hit else np.nan
            rows.append({"arm": arm, "test_id": t["test_id"], "tier": t["tier"], "gene": gbest,
                         "feature_D": t["feature_D"], "rank_of_feature_D": best, "hit": hit, "beta": beta,
                         "sign": t["sign"], "sign_correct": sc, "testable": bool(t["testable"])})
    return pd.DataFrame(rows)


def p2_summary(scores: pd.DataFrame, arm1: str | None = None, arm2: str | None = None) -> dict:
    """Recall per arm (Clopper-Pearson) over testable tests, exact McNemar arm1 vs arm2, verdict."""
    pc = PR.endpoints()["P2_primary_clinical"]["primary_contrast"]
    a1, a2 = (pc["pair"].split(" vs ") if arm1 is None else (arm1, arm2))
    s = scores[scores["testable"]]
    n = int(s["test_id"].nunique())
    per = {}
    for arm, sub in s.groupby("arm"):
        kk = int(sub["hit"].sum())
        hits = sub[sub["hit"]]
        per[arm] = {"hits": kk, "n": n, "recall": kk / n if n else np.nan, "ci95": list(clopper_pearson(kk, n)),
                    "direction_accuracy": float(hits["sign_correct"].astype(float).mean()) if len(hits) else np.nan,
                    "sign_errors": int((hits["sign_correct"] == False).sum())}  # noqa: E712
    h1 = s[s["arm"] == a1].set_index("test_id")["hit"]
    h2 = s[s["arm"] == a2].set_index("test_id")["hit"].reindex(h1.index).fillna(False)
    b, c = int((h1 & ~h2).sum()), int((~h1 & h2).sum())
    p = mcnemar_exact(b, c)
    r1 = per.get(a1, {}).get("recall", np.nan)
    r2 = per.get(a2, {}).get("recall", np.nan)
    status = "descriptive" if n < int(PR.endpoints()["P2_primary_clinical"]["min_testable_tests"]) else "evaluated"
    return {"n_testable": n, "per_arm": per, "pair": [a1, a2], "mcnemar": {"b": b, "c": c, "p": p},
            "p2_status": status, "verdict": p2_verdict(n, p, r1, r2)}


def actionability(md: pd.DataFrame, s1_hits: pd.DataFrame, cgt: pd.DataFrame, fmeta: pd.DataFrame | None = None,
                  n_perm: int | None = None, seed: int | None = None) -> pd.DataFrame:
    """Q4b: pair-level (g in T and d = feature_D) and feature-level (d is any CGT feature_D) actionable fractions of
    MD-strict and S1-obs, each vs n_perm within-list gene-label permutations (PROG/PSCORE excluded)."""
    cs = PR.dp()["clinical_s10"]
    n_perm = int(cs["actionability_null_permutations"]) if n_perm is None else n_perm
    seed = int(PR.seeds()["actionability_null_seed"]) if seed is None else seed
    pairs = {(g, f) for ts, f in zip(cgt["target_set"].map(_targets), cgt["feature_D"]) if isinstance(f, str)
             for g in ts}
    feats = {f for f in cgt["feature_D"] if isinstance(f, str)}
    excl = ("PROG", "PSCORE")
    out = []
    for name, L in (("MD-strict", md), ("S1-obs", s1_hits)):
        L = L[~L["feature"].str.split(":").str[0].isin(excl)]
        g, f = L["gene"].astype(str).to_numpy(), L["feature"].to_numpy()
        n = len(L)
        obs_pair = float(np.mean([(a, b) in pairs for a, b in zip(g, f)])) if n else np.nan
        obs_feat = float(np.mean([b in feats for b in f])) if n else np.nan
        rng = np.random.default_rng(seed)
        perm = [np.mean([(a, b) in pairs for a, b in zip(rng.permutation(g), f)]) for _ in range(n_perm)] if n else []
        out.append({"list": name, "n": n, "pair_level": obs_pair, "feature_level": obs_feat,
                    "pair_level_perm_mean": float(np.mean(perm)) if perm else np.nan,
                    "pair_level_perm_p": float((1 + np.sum(np.asarray(perm) >= obs_pair)) / (1 + len(perm))) if perm else np.nan})
    return pd.DataFrame(out)


# ------------------------------------------------------------------------------------------------
# P1 (§6.1): R-evaluability, power, pool, lists, replication, bootstrap
# ------------------------------------------------------------------------------------------------


def expected_nonnan(nan_D: pd.DataFrame, lib_D: pd.Series, lib_R: pd.Series) -> pd.DataFrame:
    """R lines x genes: gene expected non-NaN for an R line if >= R_nan_frac_min of D lines with the same library
    string are non-NaN; unseen strings use the union NaN set of their component libraries ('+'-joined)."""
    fmin = float(PR.dp()["replication_s6"]["R_nan_frac_min"])
    ok_D = ~nan_D.astype(bool)
    frac = ok_D.groupby(lib_D.reindex(nan_D.index).astype(str)).mean()
    ok_lib = frac >= fmin
    comp = {}
    for s in ok_lib.index:
        for c in s.split("+"):
            comp.setdefault(c, []).append(s)
    rows = []
    for line, s in lib_R.astype(str).items():
        if s in ok_lib.index:
            rows.append(ok_lib.loc[s].to_numpy())
        else:
            parts = [c for c in s.split("+")]
            v = np.ones(nan_D.shape[1], dtype=bool)
            for c in parts:
                if c in ok_lib.index:
                    v &= ok_lib.loc[c].to_numpy()
                elif comp.get(c):
                    v &= np.all([ok_lib.loc[x].to_numpy() for x in comp[c]], axis=0)
            rows.append(v)
    return pd.DataFrame(np.array(rows), index=lib_R.index, columns=nan_D.columns)


POWER_COLS = ("n_R", "n1_R", "n0_R", "sd_R", "nu_R", "df_nu_R")
POWER_CHUNKS_PER_WORKER = 4                 # literal-ok: load balancing of the parallel power groups (no effect on values)
POWER_CACHE_MAX = 1024                      # literal-ok: LRU size of cached shared design blocks (memory only)
_PW = {}


def _power_init(path):
    from . import workerstore as WS
    _PW["ctx"] = WS.load(path)


def _power_chunk(groups):
    return _power_groups(_PW["ctx"], groups)


def _power_groups(ctx: dict, groups) -> list:
    """The per-(gene, drops) fits of pool_power for ``groups`` [(g, dkey, row indices)]: [(indices, values, rule)]."""
    from collections import OrderedDict
    XR, fp, bin_f, covR, pats, expR = (ctx[k] for k in ("XR", "fp", "bin_f", "covR", "pats", "expR"))
    fast = ctx.get("fast", True)
    finite = np.isfinite(XR)
    model_Z = getattr(covR, "model_Z", None)
    zc, cc = OrderedDict(), {}
    out = []
    for g, dk, idx in groups:
        rows_g = expR[g].to_numpy().astype(bool)
        S = set(dk.split("|")) if dk else set()
        pat = {}
        for i in idx:
            pat.setdefault((rows_g & finite[:, fp[i]]).tobytes(), []).append(i)
        for pk, ii in pat.items():
            rows = np.frombuffer(pk, dtype=bool)
            nr = int(rows.sum())
            ii = np.asarray(ii)
            Xr = XR[rows][:, fp[ii]]
            n1 = (Xr == 1).sum(0)
            vals = {k: np.full(len(ii), np.nan) for k in POWER_COLS}
            vals["n_R"][:], vals["n1_R"][:], vals["n0_R"][:] = nr, n1, nr - n1
            rl = np.full(len(ii), "", dtype=object)
            if nr < 2:
                out.append((ii, vals, rl))
                continue
            names = covR.design_names + [nm for nm in covR.droppable_names if nm not in S]
            if fast and (model_Z is None or isinstance(covR, RDesign)):
                zk = (pk, dk)
                Zs = zc.get(zk)
                if Zs is None:
                    Zs = covR.Z(rows, names)
                    zc[zk] = Zs
                    if len(zc) > POWER_CACHE_MAX:
                        zc.popitem(last=False)
                else:
                    zc.move_to_end(zk)
                extra = covR.gene_cols(str(g), rows, S) if isinstance(covR, RDesign) else []
                Z = np.hstack([Zs, np.column_stack(extra)]) if extra else Zs        # == RDesign.model_Z
                cl = cc.get(pk)
                if cl is None:
                    cl = cc[pk] = ST.Clusters.of(pats[rows])
            else:
                Z = model_Z(g, rows, S) if model_Z is not None else covR.Z(rows, names)
                cl = ST.Clusters.of(pats[rows])
            Q = ST.orth(Z)
            Dt = ST.resid(Q, Xr)
            ccl = ST.carrier_clusters(Xr, bin_f[ii], cl)
            dd = ST.shared_df(Dt, cl, nr, Q.shape[1], car_clusters=ccl)
            vals["sd_R"][:] = np.sqrt((Dt ** 2).mean(0))
            vals["nu_R"][:] = dd["df"]
            vals["df_nu_R"][:] = dd["df_nu"]
            rl[:] = dd["df_rule"]
            out.append((ii, vals, rl))
    return out


def pool_power(U_rows: pd.DataFrame, beta_shrunk: pd.Series, sigma_g: pd.Series, XR: np.ndarray, fpos: pd.Series,
               binary: pd.Series, covR, patientsR: np.ndarray, expR: pd.DataFrame, drops: pd.Series,
               workers: int = 1, fast: bool = True) -> pd.DataFrame:
    """R-evaluability and power of every (g, d) of ``U_rows`` from R molecular data only (no R effect).

    The df is the R test's own df (DESIGN §6.1 "the test and its power share one df rule"): the same §5.3 design as
    ``RFitter.fit`` (``covR.model_Z``: design blocks, droppable covariates and, when ``covR`` is an ``RDesign``, the
    gene covariates of g, all with the D drops), the same cap and the same binary fallback (``stats.shared_df``), on
    the R lines expected non-NaN for g and non-NaN for d. SD_R(d̃_R) for SE_R_pred uses that residual too.

    ``fast`` (default): the shared design block (rows, drops) and the patient clusters (rows) are cached across genes
    (the same arrays, so identical arithmetic), and the (gene, drops) groups are split over ``workers`` processes
    (workerstore payload); every group is computed independently and written to its own rows, so the result is
    identical to ``fast=False`` (tests/test_pool_fast.py)."""
    r6 = PR.dp()["replication_s6"]
    nmin1, nmin0 = int(r6["R_min_carriers"]), int(r6["R_min_noncarriers"])
    cmin = int(PR.dp()["dictionary_s3"]["continuous_min_nonnan"])
    pm = float(r6["power_min"])
    U = U_rows[["gene", "feature"]].copy().reset_index(drop=True)
    U["gene"] = U["gene"].astype(str)
    key = pd.MultiIndex.from_arrays([U["gene"], U["feature"]])
    dr = drops.reindex(key).to_numpy()
    U["dkey"] = ["|".join(sorted(set(str(x).split("|")) - {""})) if isinstance(x, str) else "" for x in dr]
    fp = fpos.reindex(U["feature"]).to_numpy().astype(np.int64)
    bin_f = binary.reindex(U["feature"]).fillna(True).to_numpy().astype(bool)
    groups = [(g, dk, np.asarray(idx)) for (g, dk), idx in U.groupby(["gene", "dkey"], sort=False).groups.items()]
    ctx = {"XR": np.asarray(XR, dtype=np.float64), "fp": fp, "bin_f": bin_f, "covR": covR,
           "pats": np.asarray(patientsR).astype(str), "expR": expR, "fast": fast}
    w = max(1, int(workers or 1)) if fast else 1
    if w <= 1 or len(groups) < 2:
        parts = [_power_groups(ctx, groups)]
    else:
        from concurrent.futures import ProcessPoolExecutor

        from . import workerstore as WS
        cuts = np.array_split(np.arange(len(groups)), w * POWER_CHUNKS_PER_WORKER)
        with WS.staged(ctx, "pool_power") as path, \
                ProcessPoolExecutor(max_workers=w, initializer=_power_init, initargs=(path,)) as ex:
            parts = list(ex.map(_power_chunk, [[groups[i] for i in c] for c in cuts if len(c)]))
    res = {k: np.full(len(U), np.nan) for k in POWER_COLS}
    rule = np.full(len(U), "", dtype=object)
    for p in parts:
        for ii, vals, rl in p:
            for k in POWER_COLS:
                res[k][ii] = vals[k]
            rule[ii] = rl
    for k, v in res.items():
        U[k] = v
    U["df_rule"] = rule
    U["binary"] = bin_f
    U["R_evaluable"] = np.where(bin_f, (U["n1_R"] >= nmin1) & (U["n0_R"] >= nmin0), (U["n_R"] >= cmin) & (U["sd_R"] > 0))
    sig = sigma_g.reindex(U["gene"]).to_numpy()
    U["se_R_pred"] = ST.se_r_pred(sig, bin_f, U["n1_R"].to_numpy(), U["n0_R"].to_numpy(), U["sd_R"].to_numpy(),
                                  U["n_R"].to_numpy())
    U["beta_shrunk"] = beta_shrunk.reindex(key).to_numpy()
    with np.errstate(invalid="ignore", divide="ignore"):
        U["power"] = ST.power_t(U["beta_shrunk"].to_numpy(), U["se_R_pred"].to_numpy(), U["nu_R"].to_numpy())
    U["powered"] = U["R_evaluable"] & (U["power"] >= pm)
    return U.drop(columns=["dkey"])


def p1_k(pool_size: int) -> int:
    e = PR.endpoints()["P1_primary"]
    return int(min(int(e["k_max"]), pool_size // int(e["k_pool_divisor"])))


def build_lists(arms_tab: pd.DataFrame, pool: pd.DataFrame, arms, k: int | None = None) -> dict:
    """Top-k list per arm within ``pool`` (unique (g, d)), after de-duplicating (g, variant group) to the arm's
    best-ranked member inside the pool; unrankable rows (S1 p NaN, EN frequency 0) never enter; n_tied_at_rank_k.
    An arm with no ranked rows in ``arms_tab`` gets an empty list."""
    key = set(zip(pool["gene"].astype(str), pool["feature"]))
    k = p1_k(len(key)) if k is None else k
    out = {}
    for arm in arms:
        s = arms_tab[(arms_tab["arm"] == arm) & arms_tab["rank"].notna()].copy()
        s["gene"] = s["gene"].astype(str)
        inpool = np.array([(g, f) in key for g, f in zip(s["gene"], s["feature"])], dtype=bool)
        s = s[inpool].sort_values("rank", kind="mergesort")      # boolean mask (an empty list would select columns)
        vg = s["variant_group"].where(s["variant_group"].notna() & (s["variant_group"].astype(str) != ""), None)
        gk = np.where(vg.notna(), s["gene"] + "|" + vg.astype(str), s["gene"] + "|" + s["feature"])
        s = s[~pd.Index(gk).duplicated()]
        top = s.head(k)
        tied = int((s["rank_stat"] == top["rank_stat"].iloc[-1]).sum()) if len(top) else 0
        sign = np.sign(top["maj_sign"].fillna(np.sign(top["beta"]))) if "maj_sign" in top else np.sign(top["beta"])
        out[arm] = {"k": int(k), "n_list": int(len(top)), "n_tied_at_rank_k": tied,
                    "rows": [[g, f, float(sg)] for g, f, sg in zip(top["gene"], top["feature"], sign)]}
    return out


def strata_of(access_class: pd.Series) -> dict:
    """Stratum masks of the P1 access strata (primary_endpoints.P1_primary.strata_reported)."""
    a = access_class.astype(str)
    return {"direct": a == "direct", "non_direct": a != "direct", "proxy": a == "proxy"}


class RDesign(ST.CovDesign):
    """CovDesign of one line set that also carries the gene covariates (own CN, EXPGATE), so the R test and its power
    build the identical §5.3 design (``model_Z``)."""

    def __init__(self, cov: pd.DataFrame, own_cn: pd.DataFrame | None = None, expgate: pd.DataFrame | None = None):
        super().__init__(cov)
        self.own_cn, self.expgate = own_cn, expgate

    def gene_cols(self, g: str, rows, drops: set) -> list:
        extra = []
        if self.own_cn is not None and "own_CN" not in drops and g in self.own_cn.columns:
            c = self.own_cn[g].to_numpy(dtype=np.float64)[rows]
            na = np.isnan(c)
            extra.append(np.where(na, 0.0, c))
            if na.any() and not na.all():
                extra.append(na.astype(float))
        if self.expgate is not None and "EXPGATE" not in drops and g in self.expgate.columns:
            extra.append(np.nan_to_num(self.expgate[g].to_numpy(dtype=np.float64)[rows]))
        return extra

    def model_Z(self, g: str, rows, drops: set) -> np.ndarray:
        """[1, design blocks, droppable shared covariates not in drops, gene covariates of g not in drops]."""
        Z = self.Z(rows, self.design_names + [nm for nm in self.droppable_names if nm not in drops])
        extra = self.gene_cols(str(g), rows, drops)
        return np.hstack([Z, np.column_stack(extra)]) if extra else Z


class RFitter:
    """§5.3 confounder model on any line set with FIXED drops (from D), CR2 + ν; rows may be an index array with
    repeats (bootstrap) with matching cluster ids."""

    def __init__(self, X: np.ndarray, fpos: pd.Series, cov: pd.DataFrame, own_cn: pd.DataFrame | None,
                 expgate: pd.DataFrame | None, Y: pd.DataFrame, patients: np.ndarray, binary: pd.Series | None = None):
        self.X, self.fpos, self.Y = X, fpos, Y
        self.cov = RDesign(cov, own_cn, expgate)
        self.own_cn, self.expgate = own_cn, expgate
        if binary is None:
            vals = [np.unique(X[np.isfinite(X[:, j]), j]) for j in range(X.shape[1])]
            binary = pd.Series([bool(len(v) and np.isin(v, [0.0, 1.0]).all()) for v in vals], index=fpos.index)
        self.binary = binary.reindex(fpos.index).fillna(True).astype(bool)
        self.patients = np.asarray(patients).astype(str)
        self.gpos = {g: i for i, g in enumerate(map(str, Y.columns))}
        self.Yv = Y.to_numpy(dtype=np.float64)

    def fit(self, g: str, f: str, drops: str, idx: np.ndarray | None = None, clusters: np.ndarray | None = None,
            sign: float = 0.0) -> dict:
        g = str(g)
        n = self.X.shape[0]
        idx = np.arange(n) if idx is None else np.asarray(idx)
        clusters = self.patients[idx] if clusters is None else np.asarray(clusters)
        y = self.Yv[idx, self.gpos[g]]
        x = self.X[idx, int(self.fpos[f])]
        ok = np.isfinite(y) & np.isfinite(x)
        S = set(str(drops).split("|")) if isinstance(drops, str) and drops else set()
        rows = idx[ok]
        Z = self.cov.model_Z(g, rows, S)
        xo = x[ok]
        car = xo == 1
        pat = self.patients[rows]
        out = {"n": int(ok.sum()), "n_carriers": int(car.sum()), "n_carrier_patients": int(len(set(pat[car]))),
               "n_noncarrier_patients": int(len(set(pat[~car])))}
        if ok.sum() < 2:
            out.update(beta=np.nan, se=np.nan, nu=np.nan, df=np.nan, p=np.nan, p_one=np.nan)
            return out
        r = ST.pair_test(y[ok], xo, Z, clusters[ok], binary=bool(self.binary[f]))
        out.update({k: r[k] for k in ("beta", "se", "nu", "df", "p", "df_nu", "p_nu", "df_rule", "n_lev1_excluded")})
        out["p_one"] = one_sided(r["t"], r["df"], sign)
        return out


def one_sided(t: float, df: float, sign: float) -> float:
    if not np.isfinite(t) or sign == 0:
        return np.nan
    return float(sps.t.sf(t, df) if sign > 0 else sps.t.cdf(t, df))


def replicated(fit: dict, sign: float, alpha: float | None = None, binary: bool = True) -> bool:
    alpha = PR.dp()["replication_s6"]["alpha_one_sided"] if alpha is None else alpha
    if binary and (fit.get("n_carrier_patients", 0) < 2 or fit.get("n_noncarrier_patients", 0) < 2):
        return False
    b, p = fit.get("beta", np.nan), fit.get("p_one", np.nan)
    return bool(np.isfinite(b) and np.sign(b) == np.sign(sign) and np.isfinite(p) and p < alpha)


def bootstrap_rates(fitter: RFitter, lists: dict, drops: dict, binary: dict, B: int | None = None,
                    seed: int | None = None, log=None) -> dict:
    """Paired patient bootstrap over R (DESIGN §6.1): lists and drops fixed; per draw every list row refitted; copies
    get new cluster ids; degenerate pairs (< 2 distinct carrier / non-carrier patients) = not replicated.
    ``lists``: {list_name: [[g, f, sign], ...]}. Returns {'rates': draws x lists, 'frac_degenerate': ...}."""
    r6 = PR.dp()["replication_s6"]
    B = int(r6["bootstrap_B"]) if B is None else int(B)
    seed = int(PR.seeds()[r6["bootstrap_seed_key"]]) if seed is None else int(seed)
    pats = fitter.patients
    up = np.unique(pats)
    rows_of = {p: np.flatnonzero(pats == p) for p in up}
    pairs = sorted({(g, f, s) for L in lists.values() for g, f, s in L})
    rng = np.random.default_rng(seed)
    rates = np.zeros((B, len(lists)))
    degen = np.zeros((B, len(lists)))
    names = list(lists)
    for b in range(B):
        draw = rng.choice(up, size=up.size, replace=True)
        idx = np.concatenate([rows_of[p] for p in draw])
        cl = np.concatenate([np.full(rows_of[p].size, f"{p}#{j}", dtype=object) for j, p in enumerate(draw)])
        rep, deg = {}, {}
        for g, f, s in pairs:
            fit = fitter.fit(g, f, drops.get((g, f), ""), idx=idx, clusters=cl, sign=s)
            bi = binary.get(f, True)
            deg[(g, f, s)] = bool(bi and (fit["n_carrier_patients"] < 2 or fit["n_noncarrier_patients"] < 2))
            rep[(g, f, s)] = replicated(fit, s, binary=bi)
        for j, nm in enumerate(names):
            L = lists[nm]
            if L:
                rates[b, j] = np.mean([rep[tuple(x)] for x in L])
                degen[b, j] = np.mean([deg[tuple(x)] for x in L])
        if log and (b + 1) % max(1, B // 10) == 0:  # literal-ok: progress log cadence
            log(b + 1)
    return {"rates": pd.DataFrame(rates, columns=names), "frac_degenerate": pd.DataFrame(degen, columns=names)}


def p1_contrasts(point: dict, boot: pd.DataFrame, k_of: dict, prefix: str = "") -> list:
    """Contrasts of primary_endpoints.P1_primary for one stratum: point diff, percentile CI, branch/descriptive."""
    e = PR.endpoints()["P1_primary"]
    level = float(e["ci"]["level"])
    out = []
    for name, c in e["contrasts"].items():
        a1, a2 = [x.strip() for x in c["diff"].split(" - ")]
        k1, k2 = prefix + a1, prefix + a2
        if k1 not in boot or k2 not in boot:
            continue
        d = (boot[k1] - boot[k2]).to_numpy()
        lo, hi = np.quantile(d, [(1 - level) / 2, 1 - (1 - level) / 2])
        k = min(k_of.get(k1, 0), k_of.get(k2, 0))
        out.append({"contrast": name, "diff": c["diff"], "estimate": point[k1] - point[k2], "ci_lo": float(lo),
                    "ci_hi": float(hi), "k": int(k), "verdict": p1_verdict(int(k), float(lo), float(hi), name)})
    return out


# ------------------------------------------------------------------------------------------------
# Q4c population by lineage
# ------------------------------------------------------------------------------------------------


def population_by_lineage(fitter: RFitter, tests: pd.DataFrame, drops: dict, lineage: np.ndarray,
                          oncotree_code: np.ndarray, indicated_of=None) -> tuple:
    """Per testable test and target g: within-indication β (lines whose OncotreeCode is in the indication codes, all
    lines when tumour-agnostic) and β in each non-indicated lineage with >= population_min_carriers carriers and
    non-carriers; test label on-label / beyond-label / off-label-only / indication-untestable / not supported
    (support = expected sign with one-sided p < alpha_one_sided). ``indicated_of(test_row)`` overrides the exact
    OncotreeCode match (bool per fitter line)."""
    cs = PR.dp()["clinical_s10"]
    nmin = int(cs["population_min_carriers"])
    alpha = float(PR.dp()["replication_s6"]["alpha_one_sided"])
    lineage = np.asarray(lineage).astype(str)
    code = np.asarray(oncotree_code).astype(str)
    per, lab = [], []
    for _, t in tests[tests["testable"]].iterrows():
        sgn = -1.0 if t["sign"] == "S" else 1.0
        codes = set(_targets(str(t.get("indication_codes", "")).replace(";", "|")))
        agnostic = str(t.get("tumour_agnostic", "0")) in ("1", "True", "true")
        if agnostic:
            ind = np.ones(len(code), dtype=bool)
        elif indicated_of is not None:
            ind = np.asarray(indicated_of(t), dtype=bool)          # e.g. analyses.indicated_lines (parent codes)
        else:
            ind = np.isin(code, list(codes))
        for g in t["targets"]:
            if g not in fitter.gpos or t["feature_D"] not in fitter.fpos.index:
                continue
            f = t["feature_D"]
            dr = drops.get((g, f), "")
            x = fitter.X[:, int(fitter.fpos[f])]
            fin = fitter.fit(g, f, dr, idx=np.flatnonzero(ind), sign=sgn)
            n1, n0 = int(np.sum(x[ind] == 1)), int(np.sum(x[ind] == 0))
            ind_ok = n1 >= nmin and n0 >= nmin
            ind_sup = ind_ok and replicated(fin, sgn, alpha)
            per.append({"test_id": t["test_id"], "gene": g, "scope": "within-indication", "lineage": "",
                        "beta": fin["beta"], "p_one": fin["p_one"], "n_carriers": n1, "testable": ind_ok})
            n_sup_out = 0
            for L in sorted(set(lineage[~ind])):
                m = (lineage == L) & ~ind
                c1, c0 = int(np.sum(x[m] == 1)), int(np.sum(x[m] == 0))
                if c1 < nmin or c0 < nmin:
                    continue
                r = fitter.fit(g, f, dr, idx=np.flatnonzero(m), sign=sgn)
                sup = replicated(r, sgn, alpha)
                n_sup_out += int(sup)
                per.append({"test_id": t["test_id"], "gene": g, "scope": "non-indicated lineage", "lineage": L,
                            "beta": r["beta"], "p_one": r["p_one"], "n_carriers": c1, "testable": True})
            if not ind_ok:
                label = "indication-untestable" if n_sup_out == 0 else "off-label-only"
            elif ind_sup:
                label = "beyond-label" if n_sup_out else "on-label"
            else:
                label = "off-label-only" if n_sup_out else "not supported"
            lab.append({"test_id": t["test_id"], "gene": g, "label": label, "n_nonindicated_supported": n_sup_out,
                        "tumour_agnostic": agnostic})
    return pd.DataFrame(per), pd.DataFrame(lab)
