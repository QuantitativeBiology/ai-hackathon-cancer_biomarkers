"""Toy diagnostics for the 5 failures of the flagged toy suite at 8cffc91 (reported, not pass rules).

Reads the replicate pickles of scripts/toy_nc_bounds.py (--cache) and runs extra engine/oracle work:

    python -X utf8 scripts/toy_diagnostics.py run-toys --kind pure_null --idx 0 1 2 3 --cache <dir>   # worker
    python -X utf8 scripts/toy_diagnostics.py run-toys --kind ridge --idx 0 --cache <dir>             # worker
    python -X utf8 scripts/toy_diagnostics.py launch --kind pure_null --n 20 --procs 5 --cache <dir>
    python -X utf8 scripts/toy_diagnostics.py a-null --r 0 --B 100 --workers 2 --cache <dir>         # S1 null draws
    python -X utf8 scripts/toy_diagnostics.py analyse --cache <dir>      # -> results/toy/diagnostics.json + CSVs

(a) S1-obs FDP_true vs fixed-threshold FDR_emp (vectorised S1 null), decomposition of false hits into planted /
    base (lineage + latent) / noise components of the observed effect; oracle FDP and oracle permutation null.
(b) pure-null hits per list (engine) and oracle S1-pred hits on the same 20 toys.
(c) fold-artefact hits without centring per list; EN-pred refit with folds independent of cv_fold.
(d) MD-global = EN-pred selected ∩ S1-obs F1 BH: FDP_true and recall, smooth and ridge toys.
(e) ridge MD-strict false hits: correlated with planted determinants or background.
"""
from __future__ import annotations

import argparse
import copy
import json
import os
import pickle
import subprocess
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from biomark import arms as A  # noqa: E402
from biomark import gitutil  # noqa: E402
from biomark import null as NL  # noqa: E402
from biomark import params as PR  # noqa: E402
from biomark import paths as P  # noqa: E402
from biomark import pipeline as PL  # noqa: E402
from biomark import toysim as ts  # noqa: E402

OUTD = ROOT / "results" / "toy"
LISTS = ts.LISTS
DET_KINDS = ("pan_lineage", "lineage_restricted", "own_locus_cn", "msi_like", "co_located", "co_located_tight",
             "genomic_with_expression_proxy", "prog_driven", "p53_w", "predictor_blind", "wrong_sign_cgt",
             "six_carrier_strong")


def toy_seed(r):
    return int(PR.seeds()["toy_seed"]) + r


def null_seed(b):
    return int(PR.seeds()["toy_null_seed"]) + b


def load_rep(cache: Path, r: int) -> dict:
    with open(cache / f"rep{r}.pkl", "rb") as fh:
        return pickle.load(fh)


def list_hits(out, name):
    return ts.list_hits(out, name)


def md_global(out):
    a = out["arms"]
    en = a[(a["arm"] == "EN-pred") & a["selected"].astype(bool)][["gene", "feature", "beta", "stage"]]
    s1 = a[(a["arm"] == "S1-obs") & a["selected"].astype(bool)][["gene", "feature", "beta", "p"]]
    m = en.merge(s1, on=["gene", "feature"], suffixes=("", "_s1"))
    return m[np.sign(m["beta"]) == np.sign(m["beta_s1"])].reset_index(drop=True), m


def recall(hits, toy, eligible=None):
    rec, n = 0, 0
    per = {}
    for k in DET_KINDS:
        pl = toy.plants(k)
        if not len(pl):
            continue
        r_ = ts.recovered(hits, toy, k)
        if eligible is not None:
            keep = [g in eligible for g in pl["gene"]]
            r_ = r_[np.asarray(keep)]
        per[k] = [int(r_.sum()), int(len(r_))]
        rec += int(r_.sum())
        n += int(len(r_))
    return {"recall": rec / max(1, n), "n_recovered": rec, "n_plants": n, "per_kind": per}


def roles(toy):
    return toy.answer_key["pairs"].set_index(["gene", "feature"])["role"]


def false_rows(hits, toy):
    rl = roles(toy)
    h = hits.copy()
    h["gene"] = h["gene"].astype(str)
    h["role"] = [rl.get((g, f), "none") for g, f in zip(h["gene"], h["feature"])]
    return h[~h["role"].isin(ts.TRUE_ROLES | ts.NEUTRAL_ROLES)]


# ------------------------------------------------------------------------------------------------
# Decomposition of an association into planted / base / noise parts (oracle full design, drop rules)
# ------------------------------------------------------------------------------------------------


class Decomposer:
    def __init__(self, toy):
        self.toy = toy
        self.nw = ts.make_toy(seed=toy.seed, scenario="null", predictor=toy.predictor)
        self.Z, self.drop = ts.oracle_design(toy, "full")
        self.c5 = PR.dp()["confounders_s5"]
        rl = toy.answer_key["pairs"]
        self.true_of = rl[rl["role"].isin(ts.TRUE_ROLES)].groupby("gene")["feature"].apply(list).to_dict()
        self.planted_feats = set(rl.loc[rl["role"].isin(ts.TRUE_ROLES), "feature"])

    def _Z(self, x):
        Z = self.Z
        for grp, cols in self.drop.items():
            Q = ts._orth(np.column_stack([np.ones(len(Z)), Z[cols].to_numpy(float)]))
            if ts._r2_on(Q, x[:, None])[0] > self.c5["drop_assoc_r2_gt"]:
                Z = Z.drop(columns=cols)
        return Z

    def __call__(self, g, f, target="obs"):
        toy, nw = self.toy, self.nw
        xr = toy.features[f]
        x = xr.fillna(0.0).to_numpy(float)
        y = (toy.effects if target == "obs" else ts.fold_centred(toy))[g]
        ok = (y.notna() & xr.notna()).to_numpy()
        Z = self._Z(x)
        Q = ts._orth(Z.to_numpy(float)[ok])
        xt = x[ok] - Q @ (Q.T @ x[ok])
        sxx = float(xt @ xt)
        T, T0 = toy.truth[g].to_numpy(float)[ok], nw.truth[g].to_numpy(float)[ok]
        yv = y.to_numpy(float)[ok]
        comp = {"planted": T - T0, "base": T0, "noise": yv - T}
        res = {}
        for k, v in comp.items():
            vt = v - Q @ (Q.T @ v)
            res[f"beta_{k}"] = float(xt @ vt / sxx) if sxx > 0 else np.nan
        yt = yv - Q @ (Q.T @ yv)
        b = float(xt @ yt / sxx)
        df = int(ok.sum()) - Q.shape[1] - 1
        se = float(np.sqrt(((yt - b * xt) ** 2).sum() / df / sxx))
        res.update(beta=b, se_ols=se, t_ols=b / se)
        for k in ("planted", "base", "noise"):
            res[f"t_{k}"] = res[f"beta_{k}"] / se
        # correlation with the gene's planted determinants and oracle t with them added
        tf = [t for t in self.true_of.get(g, []) if t in toy.features and t != f]
        res["gene_has_planted"] = bool(tf)
        if tf:
            Xd = np.column_stack([toy.features[t].fillna(0.0).to_numpy(float)[ok] for t in tf])
            Xdt = Xd - Q @ (Q.T @ Xd)
            nr = np.sqrt((Xdt ** 2).sum(0))
            r = (xt @ Xdt) / (np.sqrt(sxx) * np.where(nr > 0, nr, np.nan))
            j = int(np.nanargmax(np.abs(r)))
            res["max_abs_r_with_planted"] = float(np.abs(r[j]))
            res["planted_feature_max_r"] = tf[j]
            Q2 = ts._orth(np.column_stack([Z.to_numpy(float)[ok], Xd]))
            xt2 = x[ok] - Q2 @ (Q2.T @ x[ok])
            yt2 = yv - Q2 @ (Q2.T @ yv)
            s2 = float(xt2 @ xt2)
            if s2 > 1e-8:
                b2 = float(xt2 @ yt2 / s2)
                df2 = df - Q2.shape[1] + Q.shape[1]
                se2 = float(np.sqrt(((yt2 - b2 * xt2) ** 2).sum() / df2 / s2))
                res["t_adjusted_for_planted"] = b2 / se2
            else:
                res["t_adjusted_for_planted"] = np.nan
        else:
            res["max_abs_r_with_planted"] = np.nan
            res["planted_feature_max_r"] = ""
            res["t_adjusted_for_planted"] = np.nan
        # dominant component (largest |t| contribution)
        comps = {k: abs(res[f"t_{k}"]) for k in ("planted", "base", "noise")}
        res["dominant"] = max(comps, key=comps.get)
        res["feature_is_planted_elsewhere"] = f in self.planted_feats
        return res


def classify(d: dict) -> str:
    """'proxy_of_planted' if the gene has planted determinants and the planted component carries >= half of t;
    'base_structure' if the lineage/latent (null-world truth) component does; else 'noise'."""
    tt = abs(d["t_planted"]) + abs(d["t_base"]) + abs(d["t_noise"])
    if tt == 0:
        return "noise"
    if abs(d["t_planted"]) >= 0.5 * abs(d["t_ols"]) and np.sign(d["t_planted"]) == np.sign(d["t_ols"]):
        return "proxy_of_planted"
    if abs(d["t_base"]) >= 0.5 * abs(d["t_ols"]) and np.sign(d["t_base"]) == np.sign(d["t_ols"]):
        return "base_structure"
    return "noise"


# ------------------------------------------------------------------------------------------------
# Workers: extra engine runs
# ------------------------------------------------------------------------------------------------


def _reduce(out):
    a = out["arms"]
    keep = {"arms_sel": a[a["selected"].astype(bool)].copy(), "md_strict": out["md_strict"],
            "md_strict_all": out.get("md_strict_all"), "lam": out["lam"], "p_cut": out["p_cut"],
            "eligibility": out["eligibility"], "timings": out["timings"]}
    s1 = out["s1"]
    keep["s1_f2_sel"] = s1[s1["selected_F2"].astype(bool)].copy() if "selected_F2" in s1 else None
    keep["s1_pred_all"] = s1[s1["arm"] == "S1-pred"][["gene", "feature", "p", "q", "df", "nu", "beta", "se",
                                                       "n_carriers", "test_family", "collinear_with_design"]].copy()
    return keep


def run_toys(kind: str, idx: list, cache: Path, workers: int) -> None:
    for i in idx:
        t0 = time.perf_counter()
        if kind == "pure_null":
            toy = ts.make_toy(seed=null_seed(i), scenario="null")
        elif kind == "ridge":
            toy = ts.make_toy(seed=toy_seed(i), predictor="ridge")
        else:
            raise ValueError(kind)
        out = PL.run_pipeline(toy, null_B=0, workers=workers)
        red = _reduce(out)
        red["wall_s"] = time.perf_counter() - t0
        red["seed"] = toy.seed
        with open(cache / f"{kind}{i}.pkl", "wb") as fh:
            pickle.dump(red, fh)
        print(kind, i, f"{red['wall_s']:.0f}s", flush=True)


def launch(kind: str, n: int, procs: int, cache: Path, start: int = 0) -> float:
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               PYTHONPATH=str(ROOT / "src"))
    ids = list(range(start, start + n))
    parts = [ids[k::procs] for k in range(procs)]
    t0 = time.perf_counter()
    ps = []
    for k, part in enumerate(parts):
        if not part:
            continue
        log = open(cache / f"{kind}_launch{k}.log", "w", encoding="utf-8")
        ps.append((subprocess.Popen([sys.executable, "-X", "utf8", __file__, "run-toys", "--kind", kind, "--idx",
                                     *map(str, part), "--cache", str(cache)], env=env, stdout=log,
                                    stderr=subprocess.STDOUT, cwd=str(ROOT)), log))
    codes = [p.wait() for p, _ in ps]
    for _, lg in ps:
        lg.close()
    if any(codes):
        raise SystemExit(f"{kind} exit codes {codes}")
    return time.perf_counter() - t0


def a_null(r: int, B: int, workers: int, cache: Path) -> None:
    """Vectorised S1 fixed-threshold null (null.run_null s1_only) on replicate r; draws b = 0..B-1."""
    pk = load_rep(cache, r)
    toy = ts.make_toy(seed=pk["seed"])
    data = PL.Data.from_toy(toy, fold_centre=True)
    out = pk["R1"]
    res = {"arms": out["arms"], "lam": out["lam"], "md_strict": out["md_strict"],
           "p_cut": {k: v for k, v in out["p_cut"].items() if k != "MD-strict"}}
    t0 = time.perf_counter()
    nul = NL.run_null(data, res, B=B, workers=workers, s1_only=True)
    rec = {"r": r, "B": B, "wall_s": time.perf_counter() - t0, "genes": len(nul["genes"]),
           "draws": nul["draws"].to_dict("list"), "summary": nul["summary"].to_dict("records"),
           "hits_obs": nul["hits_obs"], "p_cut": res["p_cut"]}
    (cache / f"a_null{r}.json").write_text(json.dumps(rec, default=float), encoding="utf-8")
    print(json.dumps({k: rec[k] for k in ("r", "B", "wall_s", "hits_obs")}, default=float))


# ------------------------------------------------------------------------------------------------
# Oracle permutation null (same strata and seeds as the engine null)
# ------------------------------------------------------------------------------------------------


def oracle_perm_null(toy, B: int, p_cut: float) -> list:
    data = PL.Data.from_toy(toy, fold_centre=True)
    st = NL.strata(data.cov_df["lineage"], data.cov_df["library_group"])
    hits = []
    for b in range(B):
        perm = NL.permutation(st, PR.null_params()["permutation_seed"] + b)
        tp = copy.copy(toy)
        tp.features = pd.DataFrame(toy.features.to_numpy()[perm], index=toy.features.index,
                                   columns=toy.features.columns)
        cv = toy.covariates.copy()
        for c in ("CN_absent", "MSI_unknown"):
            cv[c] = cv[c].to_numpy()[perm]
        tp.covariates = cv
        so = ts.oracle_s1(tp, "obs", "full")
        m = (so["test_family"] == "F1") & ~so["collinear"]
        hits.append(int((m & (so["p"] <= p_cut)).sum()))
    return hits


# ------------------------------------------------------------------------------------------------
# Analyses
# ------------------------------------------------------------------------------------------------


def analyse_a(cache, reps, n_oracle_B):
    rows, per = [], []
    for r in reps:
        pk = load_rep(cache, r)
        toy = ts.make_toy(seed=pk["seed"])
        out = pk["R1"]
        dec = Decomposer(toy)
        rec = {"r": r, "seed": pk["seed"]}
        for name in ("S1-obs", "S1-pred"):
            h = list_hits(out, name)
            f = ts.fdp_true(h, toy)
            rec[name] = {"n_hits_nonneutral": f["n_hits"], "n_false": f["n_false"], "fdp_true": f["fdp"],
                         "n_hits_raw": int(len(h))}
        nf = cache / f"a_null{r}.json"
        if nf.exists():
            nj = json.loads(nf.read_text(encoding="utf-8"))
            for name in ("S1-obs", "S1-pred"):
                d = np.asarray(nj["draws"][name], float)
                ho = nj["hits_obs"][name]
                rec[name].update(null_B=int(d.size), null_mean_hits=float(d.mean()), null_sd_hits=float(d.std(ddof=1)),
                                 fdr_emp_fixed=float(d.mean() / max(1, ho)), hits_obs_null_subset=ho,
                                 p_null_ge_realised_false=float((d >= rec[name]["n_false"]).mean()),
                                 null_draw_wall_s=nj["wall_s"] / max(1, d.size))
                rec[name]["gap_fdr_minus_fdp"] = rec[name]["fdr_emp_fixed"] - rec[name]["fdp_true"]
        # decomposition of S1-obs false hits
        fr = false_rows(list_hits(out, "S1-obs"), toy)
        cls = []
        for _, x in fr.iterrows():
            d = dec(x["gene"], x["feature"], "obs")
            c = classify(d)
            cls.append(c)
            rows.append({"r": r, "list": "S1-obs", "gene": x["gene"], "feature": x["feature"],
                         "family": x.get("family", ""), "p": x["p"], "df": x.get("df", np.nan),
                         "n_carriers": x.get("n_carriers", np.nan), "beta": x["beta"], "role": x["role"],
                         "class": c, **{k: d[k] for k in ("t_ols", "t_planted", "t_base", "t_noise",
                                                          "gene_has_planted", "max_abs_r_with_planted",
                                                          "planted_feature_max_r", "t_adjusted_for_planted",
                                                          "feature_is_planted_elsewhere")}})
        rec["S1-obs"]["false_hit_classes"] = pd.Series(cls, dtype=object).value_counts().to_dict()
        # oracle on the same toy
        so = ts.oracle_s1(toy, "obs", "full")
        oh = so[so["hit"]]
        fo = ts.fdp_true(oh, toy)
        pc = float(oh["p"].max()) if len(oh) else 0.0
        rec["oracle_S1-obs"] = {"n_hits_nonneutral": fo["n_hits"], "n_false": fo["n_false"], "fdp_true": fo["fdp"],
                                "p_cut": pc}
        ofr = false_rows(oh, toy)
        rec["oracle_S1-obs"]["false_hit_classes"] = pd.Series([classify(dec(g, f)) for g, f in
                                                               zip(ofr["gene"], ofr["feature"])],
                                                              dtype=object).value_counts().to_dict()
        eng_false = set(zip(fr["gene"], fr["feature"]))
        rec["oracle_S1-obs"]["false_overlap_with_engine"] = len(eng_false & set(zip(ofr["gene"], ofr["feature"])))
        if n_oracle_B:
            oh_null = oracle_perm_null(toy, n_oracle_B, pc)
            rec["oracle_S1-obs"].update(null_B=n_oracle_B, null_mean_hits=float(np.mean(oh_null)),
                                        fdr_emp_fixed=float(np.mean(oh_null) / max(1, len(oh))),
                                        hits_obs_raw=int(len(oh)))
        per.append(rec)
        print("a", r, json.dumps({k: rec[k] for k in ("S1-obs", "oracle_S1-obs")}, default=float)[:600], flush=True)
    return per, pd.DataFrame(rows)


def analyse_b(cache, n):
    rows, ora_rows = [], []
    any_hit = {k: [] for k in list(LISTS) + ["MD-global", "S1-obs F2", "S1-pred F2"]}
    mean_raw = {k: [] for k in any_hit}
    oracle_any = {"S1-obs": [], "S1-pred": [], "MD-strict (oracle analogue)": []}
    for b in range(n):
        with open(cache / f"pure_null{b}.pkl", "rb") as fh:
            red = pickle.load(fh)
        toy = ts.make_toy(seed=red["seed"], scenario="null")
        a = red["arms_sel"]
        out = {"arms": a, "md_strict": red["md_strict"]}
        lists = {k: ts.list_hits(out, k) for k in LISTS}
        lists["MD-global"] = md_global(out)[0]
        f2 = red["s1_f2_sel"]
        if f2 is not None:
            lists["S1-obs F2"] = f2[f2["arm"] == "S1-obs"]
            lists["S1-pred F2"] = f2[f2["arm"] == "S1-pred"]
        for k, h in lists.items():
            f = ts.fdp_true(h, toy)
            any_hit[k].append(f["n_hits"])
            mean_raw[k].append(int(len(h)))
            for _, x in false_rows(h, toy).iterrows():
                rows.append({"toy": b, "seed": red["seed"], "list": k, "gene": x["gene"], "feature": x["feature"],
                             "family": x.get("family", ""), "p": x.get("p", np.nan), "q": x.get("q", np.nan),
                             "df": x.get("df", np.nan), "nu": x.get("nu", np.nan),
                             "n_carriers": x.get("n_carriers", np.nan), "beta": x.get("beta", np.nan),
                             "rank_stat": x.get("rank_stat", np.nan), "role": x["role"]})
        so = ts.oracle_s1(toy, "obs", "full")
        sp = ts.oracle_s1(toy, "pred", "full")
        omd = ts.oracle_md_strict(toy, sp, so)
        for k, h in (("S1-obs", so[so["hit"]]), ("S1-pred", sp[sp["hit"]]), ("MD-strict (oracle analogue)", omd)):
            oracle_any[k].append(ts.fdp_true(h, toy)["n_hits"])
        spi = sp.set_index(["gene", "feature"])
        eng = red["s1_pred_all"].set_index(["gene", "feature"])
        for _, x in false_rows(sp[sp["hit"]], toy).iterrows():
            e = eng.loc[(x["gene"], x["feature"])] if (x["gene"], x["feature"]) in eng.index else None
            ora_rows.append({"toy": b, "seed": red["seed"], "gene": x["gene"], "feature": x["feature"],
                             "family": x["family"], "p_oracle": x["p"], "df_oracle": x["df"],
                             "n_carriers": x["n_carriers"], "q_oracle": x["q"],
                             "p_engine": None if e is None else float(e["p"]),
                             "q_engine": None if e is None else float(e["q"]),
                             "df_engine": None if e is None else float(e["df"])})
        # oracle p for engine S1-pred false hits
        for rr in rows:
            if rr["toy"] == b and rr["list"] == "S1-pred" and "p_oracle" not in rr:
                k = (rr["gene"], rr["feature"])
                rr["p_oracle"] = float(spi.loc[k, "p"]) if k in spi.index else None
                rr["q_oracle"] = float(spi.loc[k, "q"]) if k in spi.index else None
                rr["df_oracle"] = float(spi.loc[k, "df"]) if k in spi.index else None
        print("b", b, {k: v[-1] for k, v in any_hit.items()}, flush=True)
    summ = {k: {"p_any_nonneutral_hit": float(np.mean([v > 0 for v in vals])), "n_toys_with_hit": int(sum(v > 0 for v in vals)),
                "mean_nonneutral_hits": float(np.mean(vals)), "mean_raw_hits": float(np.mean(mean_raw[k])),
                "per_toy": vals} for k, vals in any_hit.items()}
    osum = {k: {"p_any_nonneutral_hit": float(np.mean([v > 0 for v in vals])), "per_toy": vals}
            for k, vals in oracle_any.items()}
    return summ, osum, pd.DataFrame(rows), pd.DataFrame(ora_rows)


def analyse_c(cache, reps, workers):
    per = []
    for r in reps:
        pk = load_rep(cache, r)
        toy = ts.make_toy(seed=pk["seed"])
        fa = toy.answer_key["fold_artefact"]
        key = list(zip(fa["gene"], fa["feature"]))
        out_c, out_n = pk["R1"], pk["R3"]
        el = out_c["eligibility"]
        elig = set(el.loc[el["eligible"].astype(bool), "gene"].astype(str))
        dn = PL.Data.from_toy(toy, fold_centre=False)
        s1n = A.run_s1(dn, ("S1-pred",), workers=workers)
        rec = {"r": r, "seed": pk["seed"], "artefacts": [f"{g}|{f}" for g, f in key],
               "eligible": [g in elig for g, _ in key]}

        def cnt(df):
            if df is None or not len(df):
                return 0
            s = set(zip(df["gene"].astype(str), df["feature"]))
            return sum(k in s for k in key)

        a_c, a_n = out_c["arms"], out_n["arms"]
        rec["S1-pred BH"] = {"centred": cnt(a_c[(a_c["arm"] == "S1-pred") & a_c["selected"]]),
                             "without": cnt(s1n[s1n["selected"].astype(bool)])}
        rec["EN-pred selected"] = {"centred": cnt(a_c[(a_c["arm"] == "EN-pred") & a_c["selected"]]),
                                   "without": cnt(a_n[(a_n["arm"] == "EN-pred") & a_n["selected"]])}
        rec["MD-strict"] = {"centred": cnt(out_c["md_strict"]), "without": cnt(out_n["md_strict"])}
        rec["S1-obs BH (centring irrelevant)"] = cnt(a_c[(a_c["arm"] == "S1-obs") & a_c["selected"]])

        def look(a, arm, col):
            out = []
            for g, f in key:
                m = a[(a["arm"] == arm) & (a["gene"].astype(str) == g) & (a["feature"] == f)]
                out.append(None if m.empty else float(m[col].iloc[0]))
            return out

        rec["en_freq_centred"] = look(a_c, "EN-pred", "rank_stat")
        rec["en_freq_without"] = look(a_n, "EN-pred", "rank_stat")
        rec["en_sign_cons_without"] = look(a_n, "EN-pred", "sign_consistency")
        s1nl = s1n.assign(arm="S1-pred")
        rec["s1_pred_p_without"] = look(s1nl, "S1-pred", "p")
        rec["s1_pred_q_without"] = look(s1nl, "S1-pred", "q")
        rec["s1_obs_p"] = look(a_c, "S1-obs", "p")
        mda = out_n.get("md_strict_all")
        rec["md_q_sc_without"] = look(mda.assign(arm="EN-pred"), "EN-pred", "q_sc") if mda is not None else None
        lamn = out_n["lam"].set_index("gene")
        lamc = out_c["lam"][out_c["lam"]["arm"] == "EN-pred"].set_index("gene")
        gs = [g for g, _ in key if g in lamn.index]
        rec["lam_G_without"] = [float(lamn.at[g, "lam_G"]) for g in gs]
        rec["lam_G_centred"] = [float(lamc.at[g, "lam_G"]) if g in lamc.index else None for g in gs]
        # EN-pred without centring, lambda CV'd over folds independent of cv_fold (diagnostic only)
        if gs:
            dr = copy.copy(dn)
            rng = np.random.default_rng(int(PR.seeds()["toy_seed"]) + 1000 + r)
            dr.folds = rng.permutation(dn.folds)
            U = A.universe(dr, gs)
            tab, lt, _ = A.run_en(dr, ("EN-pred",), genes=gs, U=U, workers=1)
            rec["en_freq_without_random_cv_folds"] = look(tab, "EN-pred", "rank_stat")
            rec["en_selected_without_random_cv_folds"] = cnt(tab[tab["selected"].astype(bool)])
            rec["lam_G_without_random_cv_folds"] = [float(x) for x in lt.set_index("gene").reindex(gs)["lam_G"]]
            # and with the centred run's lambda held fixed
            lam = {"EN-pred": {g: {"G": float(lamc.at[g, "lam_G"]), "X": float(lamc.at[g, "lam_X"])} for g in gs
                               if g in lamc.index}}
            tab2, _, _ = A.run_en(dn, ("EN-pred",), genes=list(lam["EN-pred"]), U=A.universe(dn, list(lam["EN-pred"])),
                                  workers=1, lam=lam)
            rec["en_freq_without_centred_lambda"] = look(tab2, "EN-pred", "rank_stat")
            rec["en_selected_without_centred_lambda"] = cnt(tab2[tab2["selected"].astype(bool)])
        per.append(rec)
        print("c", r, json.dumps(rec, default=float)[:800], flush=True)
    return per


def fdp_rec(h, toy, elig):
    f = ts.fdp_true(h, toy)
    rc = recall(h, toy)
    rce = recall(h, toy, elig)
    return {"n_hits_nonneutral": f["n_hits"], "n_false": f["n_false"], "fdp_true": f["fdp"],
            "recall_all_plants": rc["recall"], "recall_eligible_plants": rce["recall"],
            "n_plants": rc["n_plants"], "n_plants_eligible": rce["n_plants"], "per_kind": rce["per_kind"]}


def analyse_d(cache, reps, ridge_ids):
    sm, rg = [], []
    for r in reps:
        pk = load_rep(cache, r)
        toy = ts.make_toy(seed=pk["seed"])
        out = pk["R1"]
        el = out["eligibility"]
        elig = set(el.loc[el["eligible"].astype(bool), "gene"].astype(str))
        mg, mg_any = md_global(out)
        rec = {"r": r, "seed": pk["seed"], "MD-global": fdp_rec(mg, toy, elig),
               "MD-global_ignoring_sign": fdp_rec(mg_any, toy, elig)}
        for k in LISTS:
            rec[k] = fdp_rec(ts.list_hits(out, k), toy, elig)
        sm.append(rec)
    for i in ridge_ids:
        p = cache / f"ridge{i}.pkl"
        if not p.exists():
            continue
        with open(p, "rb") as fh:
            red = pickle.load(fh)
        toy = ts.make_toy(seed=red["seed"], predictor="ridge")
        out = {"arms": red["arms_sel"], "md_strict": red["md_strict"]}
        el = red["eligibility"]
        elig = set(el.loc[el["eligible"].astype(bool), "gene"].astype(str))
        mg, mg_any = md_global(out)
        rec = {"i": i, "seed": red["seed"], "MD-global": fdp_rec(mg, toy, elig),
               "MD-global_ignoring_sign": fdp_rec(mg_any, toy, elig)}
        for k in LISTS:
            rec[k] = fdp_rec(ts.list_hits(out, k), toy, elig)
        rg.append(rec)
    return sm, rg


def analyse_e(cache, ridge_ids):
    per, rows = [], []
    for i in ridge_ids:
        p = cache / f"ridge{i}.pkl"
        if not p.exists():
            continue
        with open(p, "rb") as fh:
            red = pickle.load(fh)
        toy = ts.make_toy(seed=red["seed"], predictor="ridge")
        dec = Decomposer(toy)
        out = {"arms": red["arms_sel"], "md_strict": red["md_strict"]}
        rec = {"i": i, "seed": red["seed"]}
        for name in ("MD-strict", "S1-obs", "EN-pred selected"):
            h = ts.list_hits(out, name)
            f = ts.fdp_true(h, toy)
            rec[name] = {"n_hits_nonneutral": f["n_hits"], "n_false": f["n_false"], "fdp_true": f["fdp"]}
            if name == "MD-strict":
                fr = false_rows(h, toy)
                cl = []
                for _, x in fr.iterrows():
                    d = dec(x["gene"], x["feature"], "obs")
                    c = classify(d)
                    cl.append(c)
                    rows.append({"ridge_toy": i, "seed": red["seed"], "gene": x["gene"], "feature": x["feature"],
                                 "family": toy.feature_meta["family"].get(x["feature"], ""), "p_sc": x.get("p_sc"),
                                 "q_sc": x.get("q_sc"), "selection_freq": x.get("selection_freq"), "beta": x["beta"],
                                 "role": x["role"], "class": c,
                                 **{k: d[k] for k in ("t_ols", "t_planted", "t_base", "t_noise", "gene_has_planted",
                                                      "max_abs_r_with_planted", "planted_feature_max_r",
                                                      "t_adjusted_for_planted", "feature_is_planted_elsewhere")}})
                rec[name]["false_hit_classes"] = pd.Series(cl, dtype=object).value_counts().to_dict()
                rec[name]["false_in_planted_genes"] = int(sum(r_["gene_has_planted"] for r_ in rows
                                                              if r_["ridge_toy"] == i))
        per.append(rec)
        print("e", i, json.dumps(rec, default=float)[:600], flush=True)
    return per, pd.DataFrame(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run-toys", "launch", "a-null", "analyse"])
    ap.add_argument("--kind", default="pure_null")
    ap.add_argument("--idx", type=int, nargs="*", default=[])
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--procs", type=int, default=5)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--r", type=int, default=0)
    ap.add_argument("--B", type=int, default=100)
    ap.add_argument("--oracle-B", type=int, default=20)
    ap.add_argument("--ridge-ids", type=int, nargs="*", default=[0, 1, 2, 3, 4])
    ap.add_argument("--parts", default="abcde")
    ap.add_argument("--cache", type=Path, required=True)
    a = ap.parse_args(argv)
    if a.cmd == "run-toys":
        run_toys(a.kind, a.idx, a.cache, a.workers)
        return 0
    if a.cmd == "launch":
        w = launch(a.kind, a.n, a.procs, a.cache, a.start)
        (a.cache / f"{a.kind}_wall.json").write_text(json.dumps({"wall_s": w, "n": a.n}), encoding="utf-8")
        print(f"{a.kind}: {a.n} toys in {w / 60:.1f} min")
        return 0
    if a.cmd == "a-null":
        a_null(a.r, a.B, a.workers, a.cache)
        return 0
    reps = list(range(int(PR.dp()["toy_s13"]["nc_bound_replicates"])))
    dj = OUTD / "diagnostics.json"
    res = json.loads(dj.read_text(encoding="utf-8")) if dj.exists() else {}
    res.update({"commit": gitutil.head(P.PROJECT), "config_sha256": PR.config_hashes(),
                "toy_answer_key_sha256": ts.toy_config_sha256(), "status": "reported diagnostics, not pass rules"})
    if "a" in a.parts:
        per, fr = analyse_a(a.cache, reps, a.oracle_B)
        res["a_s1_obs_null_calibration"] = per
        fr.to_csv(OUTD / "diag_a_s1_false_hits.csv", index=False)
    if "b" in a.parts:
        summ, osum, hr, orr = analyse_b(a.cache, a.n)
        res["b_pure_null"] = {"engine": summ, "oracle": osum, "n_toys": a.n}
        hr.to_csv(OUTD / "diag_b_pure_null_hits.csv", index=False)
        orr.to_csv(OUTD / "diag_b_pure_null_oracle_s1_pred_hits.csv", index=False)
    if "c" in a.parts:
        res["c_fold_artefact_without_centring"] = analyse_c(a.cache, reps, a.workers)
    if "d" in a.parts:
        sm, rg = analyse_d(a.cache, reps, a.ridge_ids)
        res["d_md_global"] = {"definition": "EN-pred selected ∩ S1-obs F1 BH hits, same sign (MD-global_ignoring_sign "
                                            "alongside)", "smooth": sm, "ridge": rg}
    if "e" in a.parts:
        per, er = analyse_e(a.cache, a.ridge_ids)
        res["e_ridge_md_strict"] = per
        er.to_csv(OUTD / "diag_e_ridge_md_false_hits.csv", index=False)
    dj.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
