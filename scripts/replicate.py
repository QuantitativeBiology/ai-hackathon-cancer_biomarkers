"""The single replication run (DESIGN v3.2 §1.3, §6.1, primary_endpoints.P1_primary) and its dry run.

    PART2_UNBLIND=1 python -B scripts/replicate.py                 # the one R read (after tag freeze-catalogue)
    python -B scripts/replicate.py --pseudo-r fold0                # dry run: D fold 0 as pseudo-R -> results/dryrun/

Real run preconditions (refusal otherwise): HEAD descends from freeze-catalogue; the guard access log shows no earlier
granted R read (a reproduction needs --reproduce with the same commit and config hash, §1.3 h); every file of
configs/freeze_catalogue.json matches its sha256; guard.load_effects conditions (tag, PART2_UNBLIND=1, clean tree,
purpose). Reads configs/catalogue/r_requests.csv and results/benchmark/p1_pool.json (lists, drops and pool fixed),
R effects once (scope replication) and set-A effects once (scope limits).

Writes (results/replication/ or results/dryrun/): r_tests.csv (β, SE, ν, df, p two-sided and one-sided in the
expected direction, n, carriers, power, replicated), p1_summary.json (point rates, bootstrap CIs, contrasts and
verdicts per stratum), p1_by_access.csv, bootstrap_p1.csv (draws x lists), r_lines.csv (R accurate-dependent
membership of MD-strict genes), setA_lines.csv (set-A dependent fraction per requested gene).

Dry run: pseudo-R = D lines with cv_fold 0, D' = the other folds. Lists stay the frozen full-D lists (in-sample, as
§1.3 states; registry status dryrun, uncitable). The D-side quantities are refitted on D': S1-obs (β, SE, drops,
d̃ SD) for the requested genes with the unchanged engine, σ̂_resid,D'(g), the expected-NaN rule, EB shrinkage, the
pool and power; fold-0 effects are read with scope dev. Set A is not read (setA_lines.csv carries status dryrun).
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from biomark import endpoint_io as IO  # noqa: E402
from biomark import endpoints as EP  # noqa: E402
from biomark import eligibility as ELG  # noqa: E402
from biomark import gitutil  # noqa: E402
from biomark import guard  # noqa: E402
from biomark import p1  # noqa: E402
from biomark import params as PR  # noqa: E402
from biomark import paths as P  # noqa: E402
from biomark import stats as ST  # noqa: E402

TAG = "freeze-catalogue"
SIGN = {"S": -1.0, "W": 1.0}


def config_digest() -> str:
    h = PR.config_hashes()
    return hashlib.sha256("".join(f"{k}={h[k]}" for k in sorted(h)).encode()).hexdigest()[:16]  # literal-ok: id length


def fitter_for(ids, genes, features, Y: pd.DataFrame) -> EP.RFitter:
    X = IO.line_values(ids, features)
    oc, eg = IO.line_gene_covs(ids, genes)
    cov = IO.line_covariates(ids)
    return EP.RFitter(X.to_numpy(dtype=np.float64), pd.Series(np.arange(X.shape[1]), index=X.columns), cov, oc, eg,
                      Y.reindex(index=list(ids), columns=list(genes)), cov["PatientID"].astype(str).to_numpy())


def list_sets(pool_json: dict) -> tuple:
    lists, k_of = {}, {}
    for stratum, L in pool_json["lists"].items():
        for arm, v in L.items():
            nm = f"{stratum}|{arm}"
            lists[nm] = [(str(g), f, float(s)) for g, f, s in v["rows"]]
            k_of[nm] = int(v["k"])
    return lists, k_of


def core(fitter: EP.RFitter, req: pd.DataFrame, pool_json: dict, binary: dict, B: int, out: Path, mode: str,
         log=print) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for _, r in req.iterrows():
        s = SIGN.get(str(r.get("expected_sign", "")), float(np.sign(r["beta_D"])) if np.isfinite(r["beta_D"]) else 0.0)
        fit = fitter.fit(r["gene"], r["feature"], r["dropped_cov"] if isinstance(r["dropped_cov"], str) else "", sign=s)
        rows.append({"gene": r["gene"], "feature": r["feature"], "reasons": r["reasons"], "sign_tested": s,
                     "power": r.get("power", np.nan), **fit,
                     "replicated": EP.replicated(fit, s, binary=binary.get(r["feature"], True)),
                     "se_type": PR.pp()["se_type"]})
    rt = pd.DataFrame(rows)
    rt.to_csv(out / "r_tests.csv", index=False)
    lists, k_of = list_sets(pool_json)
    drops = {(str(g), f): (d if isinstance(d, str) else "") for g, f, d in zip(req["gene"], req["feature"], req["dropped_cov"])}
    point = {}
    for nm, L in lists.items():
        reps = [EP.replicated(fitter.fit(g, f, drops.get((g, f), ""), sign=s), s, binary=binary.get(f, True)) for g, f, s in L]
        point[nm] = float(np.mean(reps)) if reps else np.nan
    bt = EP.bootstrap_rates(fitter, lists, drops, binary, B=B, log=lambda b: log(f"bootstrap {b}/{B}"))
    bt["rates"].to_csv(out / "bootstrap_p1.csv", index=False)
    summ = {"mode": mode, "B": B, "strata": {}}
    by = []
    for stratum in pool_json["lists"]:
        k = pool_json["k"] if stratum == "G" else pool_json["strata"][stratum]["k"]
        kof = {nm: k for nm in lists if nm.startswith(stratum + "|")}
        cs = EP.p1_contrasts(point, bt["rates"], kof, prefix=f"{stratum}|")
        summ["strata"][stratum] = {"k": int(k), "rates": {nm.split("|", 1)[1]: point[nm] for nm in kof},
                                   "frac_degenerate": {nm.split("|", 1)[1]: float(bt["frac_degenerate"][nm].mean())
                                                       for nm in kof}, "contrasts": cs}
        for c in cs:
            by.append({"stratum": stratum, **c})
    pd.DataFrame(by).to_csv(out / "p1_by_access.csv", index=False)
    IO.write_json(out / "p1_summary.json", summ)
    return summ


def r_lines(md_genes, D_ids, R_ids, Y_D: pd.DataFrame, Y_R: pd.DataFrame, pred_R: pd.DataFrame) -> pd.DataFrame:
    """R accurate-dependent membership (DESIGN §2): dependent by the D-fitted cut and prediction <= the R quantile at
    n_g / n_D (target-free threshold)."""
    if not len(md_genes):
        return pd.DataFrame(columns=["gene", "ModelID", "dependent", "accurate_dependent"])
    rule = ELG.dependent_rule(Y_D[md_genes])
    dep = ELG.dependent(Y_R[md_genes], rule)
    thr = ELG.accurate_threshold_R(pred_R[md_genes], (ELG.dependent(Y_D[md_genes], rule) == 1).sum(), rule["n_obs"])
    out = []
    for g in md_genes:
        acc = (dep[g] == 1) & (pred_R[g] <= thr[g])
        out.append(pd.DataFrame({"gene": g, "ModelID": dep.index, "dependent": dep[g].to_numpy(),
                                 "accurate_dependent": acc.to_numpy()}))
    return pd.concat(out, ignore_index=True)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pseudo-r", choices=["fold0"], default=None)
    ap.add_argument("--reproduce", action="store_true")
    ap.add_argument("--boot-B", type=int, default=None, help="dry run only (the real run uses the design B)")
    ap.add_argument("--no-refit", action="store_true", help="dry run only: use full-D S1-obs instead of refitting on D'")
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args(argv)
    if not a.pseudo_r:
        if a.boot_B is not None or a.no_refit:
            raise SystemExit("--boot-B / --no-refit are dry-run options")
        EP.require_descends(TAG)
        cfg = config_digest()
        commit = gitutil.head(P.PROJECT)
        EP.require_single_r_read(guard.read_log(), {"commit": commit, "config_sha256": cfg} if a.reproduce else None)
        bad = EP.catalogue_mismatches(IO.read_json(P.CONFIGS / "freeze_catalogue.json"))
        if bad:
            raise EP.EndpointRefusal(f"frozen catalogue files changed since freeze_catalogue.json: {bad}")
    req = pd.read_csv(P.CONFIGS / "catalogue" / "r_requests.csv", dtype={"gene": str})
    pool_json = IO.read_json(P.RESULTS / "benchmark" / "p1_pool.json")
    fl = IO.feature_meta()
    binary = (fl["d_type"] == "binary").to_dict()
    lists, _ = list_sets(pool_json)
    genes = sorted(set(req["gene"]) | {g for L in lists.values() for g, _, _ in L})
    feats = sorted(set(req["feature"]) | {f for L in lists.values() for _, f, _ in L})
    md_genes = sorted(set(req.loc[req["reasons"].str.contains("MD-strict"), "gene"]))
    dev = list(p1.dev_ids())
    if a.pseudo_r:
        if a.reproduce:
            raise SystemExit("--reproduce applies to the real run only")
        out = P.RESULTS / "dryrun"
        folds = p1.folds().reindex(dev)
        R_ids = [i for i in dev if folds[i] == 0]
        D2 = [i for i in dev if folds[i] != 0]
        Y_all = guard.load_effects(dev, scope="dev", genes=genes)
        req2 = req.copy()
        if not a.no_refit:
            from biomark import arms as A
            from biomark import pipeline as PL
            d = PL.Data.from_real(genes=genes, save_subsamples=False)
            dD = IO.subset_data(d, np.isin(np.asarray(d.lines), D2))
            s1 = A.run_s1(dD, ("S1-obs",), workers=a.workers)
            u = s1.set_index(["gene", "feature"])
            key = pd.MultiIndex.from_arrays([req2["gene"], req2["feature"]])
            req2["beta_D"] = u["beta"].reindex(key).to_numpy()
            req2["dropped_cov"] = u["dropped_cov"].reindex(key).fillna("").to_numpy()
            nan_D = pd.read_parquet(P.CACHE / "effect_nan_mask.parquet").reindex(index=D2, columns=genes)
            covD = ST.CovDesign(IO.line_covariates(D2))
            sig = IO.sigma_resid_genes(Y_all.loc[D2].to_numpy(), ~nan_D.to_numpy(dtype=bool), covD, genes)
            arms_tab = IO.arms_table()
            s1o = s1.copy()                                  # finish_s1 already set test_family F1 / F2
            acc = arms_tab[arms_tab["arm"] == "S1-obs"].drop_duplicates(["gene", "feature"]).set_index(["gene", "feature"])["access_class"]
            s1o["access_class"] = acc.reindex(pd.MultiIndex.from_arrays([s1o["gene"], s1o["feature"]])).to_numpy()
            pool, pj = IO.make_pool(arms_tab[arms_tab["gene"].isin(genes)], s1o, sig, D2, R_ids, nan_D)
            out.mkdir(parents=True, exist_ok=True)
            pool.to_parquet(out / "p1_pool.parquet")
            pj["note"] = "dry run: D' = folds 1-4, pseudo-R = fold 0; lists frozen from full-D ranks (in-sample)"
            IO.write_json(out / "p1_pool.json", pj)
            pool_json = pj
        fitter = fitter_for(R_ids, genes, feats, Y_all.loc[R_ids])
        summ = core(fitter, req2, pool_json, binary, a.boot_B or PR.dp()["replication_s6"]["bootstrap_B"], out,
                    "dryrun")
        pr = p1.oof_predictions("m5", columns=md_genes, ids=R_ids) if md_genes else pd.DataFrame(index=R_ids)
        r_lines(md_genes, D2, R_ids, Y_all.loc[D2], Y_all.loc[R_ids], pr).to_csv(out / "r_lines.csv", index=False)
        pd.DataFrame({"gene": genes, "status": "dryrun: set A not read", "n_setA": np.nan, "n_dependent": np.nan,
                      "frac_dependent": np.nan}).to_csv(out / "setA_lines.csv", index=False)
        IO.write_json(out / "registry_status.json", {"status": "dryrun", "citable": False,
                                                     "utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")})
        print(summ["strata"].get("G", {}))
        return 0
    # ---- the real, single R run (preconditions checked above) -------------------------------------------------
    test, nc = list(p1.test_ids()), list(p1.noncancerous_ids())
    Y_R = guard.load_effects(test, scope="replication", purpose=f"replicate.py single R run cfg={cfg}", genes=genes)
    Y_A = guard.load_effects(nc, scope="limits", purpose=f"limits: set-A normal-line flag cfg={cfg}", genes=genes)
    Y_D = guard.load_effects(dev, scope="dev", genes=genes)
    out = P.RESULTS / "replication"
    fitter = fitter_for(test, genes, feats, Y_R)
    core(fitter, req, pool_json, binary, int(PR.dp()["replication_s6"]["bootstrap_B"]), out, "replication")
    pr = p1.test_predictions("m5_avg", columns=md_genes, ids=test) if md_genes else pd.DataFrame(index=test)
    r_lines(md_genes, dev, test, Y_D, Y_R, pr).to_csv(out / "r_lines.csv", index=False)
    rule = ELG.dependent_rule(Y_D)
    depA = ELG.dependent(Y_A, rule)
    pd.DataFrame({"gene": genes, "n_setA": depA.notna().sum().to_numpy(), "n_dependent": (depA == 1).sum().to_numpy(),
                  "frac_dependent": ((depA == 1).sum() / depA.notna().sum()).to_numpy()}).to_csv(out / "setA_lines.csv",
                                                                                                index=False)
    IO.write_json(out / "run.json", {"commit": commit, "config_digest": cfg, "status": "ok",
                                     "utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")})
    return 0


if __name__ == "__main__":
    sys.exit(main())
