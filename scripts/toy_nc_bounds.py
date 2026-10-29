"""NC2-NC6 bounds from 5 toy replicates (DESIGN v3.2 §7.3; definitions fixed in docs/toy_nc_bounds.md first).

    python -X utf8 scripts/toy_nc_bounds.py --cache <dir>            # launch 5 replicates, timebox, write JSON
    python -X utf8 scripts/toy_nc_bounds.py --replicate r --cache <dir>   # one replicate (used by the launcher)
    python -X utf8 scripts/toy_nc_bounds.py --aggregate-only --cache <dir> --wall-s <s>   # re-aggregate

Per replicate r (seed toy_seed + r, smooth truth toy): R1 run_pipeline(fold_centre=True, null_B=0); R2 low-tier
universe run_arms(("S1-obs", "EN-pred")); R3 run_pipeline(fold_centre=False, null_B=0, arms=("EN-pred",)).
Engine outputs are pickled to --cache (scratch; reused by scripts/toy_diagnostics.py). Writes
results/toy/nc_bounds.json and appends to results/deviations.csv when a fallback triggers.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
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
from biomark import params as PR  # noqa: E402
from biomark import paths as P  # noqa: E402
from biomark import pipeline as PL  # noqa: E402
from biomark import toysim as ts  # noqa: E402

N_REP = int(PR.dp()["toy_s13"]["nc_bound_replicates"])
NC = PR.dp()["negative_controls_s7"]
PASS = ts.load_toy_config()["pass"]
TIMEBOX_MIN = float(PR.tiers()["tiers"]["A"]["stages"]["toy_run_bounds"])
WORKERS_PER_REP = 2
DOC = ROOT / "docs" / "toy_nc_bounds.md"
OUT = ROOT / "results" / "toy" / "nc_bounds.json"
DEV = ROOT / "results" / "deviations.csv"
DEV_COLS = ["id", "utc", "stage", "item", "provisional", "toy_value", "final_value", "reason", "source", "agent"]


def seed_of(r: int) -> int:
    return int(PR.seeds()[PR.dp()["toy_s13"]["seed_key"]]) + r


# ------------------------------------------------------------------------------------------------
# One replicate
# ------------------------------------------------------------------------------------------------


def run_replicate(r: int, cache: Path, workers: int) -> dict:
    t0 = time.perf_counter()
    toy = ts.make_toy(seed=seed_of(r))
    tim = {}
    t = time.perf_counter()
    out1 = PL.run_pipeline(toy, fold_centre=True, null_B=0, workers=workers)
    tim["R1_s"] = time.perf_counter() - t
    t = time.perf_counter()
    low = list(toy.plants("low_tier")["gene"])
    d2 = PL.Data.from_toy(toy, fold_centre=True, genes=low)
    res2 = A.run_arms(d2, arms=("S1-obs", "EN-pred"), workers=workers)
    tim["R2_s"] = time.perf_counter() - t
    t = time.perf_counter()
    out3 = PL.run_pipeline(toy, fold_centre=False, null_B=0, arms=("EN-pred",), workers=workers)
    tim["R3_s"] = time.perf_counter() - t
    t = time.perf_counter()
    run_ridge_nc2(r, cache, workers)
    tim["R2_ridge_s"] = time.perf_counter() - t
    tim["total_s"] = time.perf_counter() - t0
    keep2 = {k: res2[k] for k in ("arms", "s1", "en", "lam", "md_strict", "md_strict_all", "p_cut", "p_cut_sc")
             if k in res2}
    with open(cache / f"rep{r}.pkl", "wb") as fh:
        pickle.dump({"r": r, "seed": seed_of(r), "R1": out1, "R2": keep2, "R3": out3, "timing": tim}, fh)
    return stats_from_cache(r, cache, toy)


def nc2_stats(toy, res2: dict) -> dict:
    low = list(toy.plants("low_tier")["gene"])
    lam = res2["lam"]
    fid = lam[lam["arm"] == "EN-pred"].set_index("gene")["r2_oof_beyond_design"].reindex(low)
    md2, en2 = res2["md_strict"], res2["en"]
    nsel = int(en2["selected"].astype(bool).sum()) if len(en2) else 0
    return {"NC2a": {"value": float(np.nanmedian(fid.to_numpy())) if fid.notna().any() else None,
                     "n": int(fid.notna().sum()), "count": nsel,
                     "per_gene_fidelity": {g: (None if pd.isna(v) else float(v)) for g, v in fid.items()}},
            "NC2b": {"value": float(len(md2)) / len(low), "n": len(low), "count": int(len(md2)),
                     "md_strict_rows": int(len(md2)), "en_pred_selected": nsel,
                     "md_strict_pairs": [f"{g}|{f}" for g, f in zip(md2["gene"], md2["feature"])]}}


def run_ridge_nc2(r: int, cache: Path, workers: int) -> dict:
    """Ruling R1: NC2 on the ridge toy (same seed), low-tier universe, S1-obs + EN-pred."""
    toy = ts.make_toy(seed=seed_of(r), predictor="ridge")
    low = list(toy.plants("low_tier")["gene"])
    d2 = PL.Data.from_toy(toy, fold_centre=True, genes=low)
    res2 = A.run_arms(d2, arms=("S1-obs", "EN-pred"), workers=workers)
    st = nc2_stats(toy, res2)
    st["seed"] = int(toy.seed)
    st["commit"] = gitutil.head(P.PROJECT)
    (cache / f"rep{r}_ridge_nc2.json").write_text(json.dumps(st, indent=1, default=float), encoding="utf-8")
    return st


def stats_from_cache(r: int, cache: Path, toy=None) -> dict:
    with open(cache / f"rep{r}.pkl", "rb") as fh:
        pk = pickle.load(fh)
    toy = ts.make_toy(seed=pk["seed"]) if toy is None else toy
    stats = replicate_stats(toy, pk["R1"], pk["R2"], pk["R3"])
    stats["timing"] = pk["timing"]
    (cache / f"rep{r}.json").write_text(json.dumps(stats, indent=1, default=float), encoding="utf-8")
    return stats


def _sel(out, arm):
    a = out["arms"]
    return a[(a["arm"] == arm) & a["selected"].astype(bool)]


def replicate_stats(toy, out1: dict, res2: dict, out3: dict) -> dict:
    ak = toy.answer_key
    s = {"seed": int(toy.seed)}
    # NC2 ------------------------------------------------------------------------------------------
    low = list(toy.plants("low_tier")["gene"])
    lam = res2["lam"]
    fid = lam[lam["arm"] == "EN-pred"].set_index("gene")["r2_oof_beyond_design"].reindex(low)
    md2 = res2["md_strict"]
    en2 = res2["en"]
    p_sc1 = out1["p_cut"].get("MD-strict")
    md2_fixed = A.strict(en2, res2["s1"], "EN-pred", p_cut=p_sc1) if len(en2) else pd.DataFrame(columns=["pass"])
    s["NC2a"] = {"value": float(np.nanmedian(fid.to_numpy())) if fid.notna().any() else None,
                 "n": int(fid.notna().sum()), "count": int(en2["selected"].astype(bool).sum()) if len(en2) else 0,
                 "per_gene_fidelity": {g: (None if pd.isna(v) else float(v))
                                                                   for g, v in fid.items()}}
    s["NC2b"] = {"value": float(len(md2)) / len(low), "n": len(low), "count": int(len(md2)),
                 "md_strict_rows": int(len(md2)),
                 "md_strict_pairs": [f"{g}|{f}" for g, f in zip(md2["gene"], md2["feature"])],
                 "en_pred_selected": int(en2["selected"].astype(bool).sum()) if len(en2) else 0,
                 "md_strict_rows_with_R1_p_cut_sc": int(md2_fixed["pass"].sum()) if len(md2_fixed) else 0,
                 "R1_p_cut_sc": p_sc1}
    # NC3 ------------------------------------------------------------------------------------------
    sn = out1["self_null"].copy()
    sn["gene"] = sn["gene"].astype(str)
    key = set(zip(ak["self_null"]["gene"], ak["self_null"]["feature"]))
    planted = sn[[(g, f) in key for g, f in zip(sn["gene"], sn["feature"])]]
    pt = planted[planted["testable"].astype(bool)]
    el = out1["eligibility"]
    elig = set(el.loc[el["eligible"].astype(bool), "gene"].astype(str))
    allt = sn[sn["testable"].astype(bool)]
    elt = allt[allt["gene"].isin(elig)]
    s["NC3"] = {"value": float(pt["frac_dependent"].mean()) if len(pt) else None, "n": int(len(pt)),
                "count": int(np.nansum(pt["frac_dependent"].to_numpy())) if len(pt) else 0,
                "planted_rows": planted[["gene", "feature", "n_carriers", "testable", "carrier_median",
                                         "frac_dependent"]].to_dict("records"),
                "secondary_all_genes": {"value": float(allt["frac_dependent"].mean()) if len(allt) else None,
                                        "n": int(len(allt))},
                "secondary_eligible_genes": {"value": float(elt["frac_dependent"].mean()) if len(elt) else None,
                                             "n": int(len(elt))}}
    # NC4 ------------------------------------------------------------------------------------------
    fa = ak["fold_artefact"]
    a1 = out1["arms"]
    a3 = out3["arms"]

    def freq(a, g, f):
        m = a[(a["arm"] == "EN-pred") & (a["gene"].astype(str) == g) & (a["feature"] == f)]
        return None if m.empty else float(m["rank_stat"].iloc[0])

    sc, sn_ = _sel(out1, "EN-pred"), _sel(out3, "EN-pred")
    after = sum(((sc["gene"].astype(str) == g) & (sc["feature"] == f)).any() for g, f in zip(fa["gene"], fa["feature"]))
    before = sum(((sn_["gene"].astype(str) == g) & (sn_["feature"] == f)).any() for g, f in zip(fa["gene"], fa["feature"]))
    s["NC4"] = {"after_centring": int(after), "without_centring": int(before), "n": int(len(fa)),
                "freq_after": [freq(a1, g, f) for g, f in zip(fa["gene"], fa["feature"])],
                "freq_without": [freq(a3, g, f) for g, f in zip(fa["gene"], fa["feature"])],
                "artefact_eligible": [g in elig for g in fa["gene"]]}
    # NC5 ------------------------------------------------------------------------------------------
    cgt = toy.cgt
    tests = cgt[(cgt["partition"] == "bypass-resistance") | (cgt["sign"] == "N")]
    cs = out1["cgt_scores"]
    per_arm = {}
    for arm in A.ARMS:
        aa = a1[a1["arm"] == arm]
        n_rank = n_fs = 0
        rows = []
        for _, t in tests.iterrows():
            c = cs[(cs["arm"] == arm) & (cs["test_id"] == t["test_id"])]
            rk = float(c["rank_of_feature_D"].iloc[0]) if len(c) else np.nan
            hit = bool(c["hit"].iloc[0]) if len(c) else False
            beta = np.nan
            for g in str(t["target_set"]).split("|"):
                m = aa[(aa["gene"].astype(str) == g) & (aa["feature"] == t["feature_D"])]
                if len(m):
                    beta = float(m["beta"].iloc[0])
            if np.isfinite(rk):
                n_rank += 1
                n_fs += int(hit and beta < 0)
            rows.append({"test_id": t["test_id"], "sign": t["sign"], "partition": t["partition"],
                         "rank_of_feature_D": None if not np.isfinite(rk) else rk, "hit": hit,
                         "beta": None if not np.isfinite(beta) else beta})
        per_arm[arm] = {"rate": (n_fs / n_rank) if n_rank else None, "n": n_rank, "false_sens": n_fs, "tests": rows}
    rates = [v["rate"] for v in per_arm.values() if v["rate"] is not None]
    s["NC5"] = {"value": max(rates) if rates else None, "n": min((v["n"] for v in per_arm.values()), default=0),
                "count": int(sum(v["false_sens"] for v in per_arm.values())),
                "per_arm": per_arm, "n_tests": int(len(tests))}
    # NC6 ------------------------------------------------------------------------------------------
    gc = out1["gene_classes"].astype({"gene": str}).set_index("gene")["failure_class"].astype(str)
    lo = sorted(toy.plants("lineage_only")["gene"])
    cls = gc.reindex(lo)
    s["NC6"] = {"value": float(cls.eq("3").mean()), "n": len(lo), "count": int(cls.eq("3").sum()),
                "classes": cls.fillna("NA").to_dict()}
    return s


# ------------------------------------------------------------------------------------------------
# Aggregation
# ------------------------------------------------------------------------------------------------


def _append_deviation(rows: list) -> None:
    new = not DEV.exists()
    with open(DEV, "a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=DEV_COLS)
        if new:
            w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in DEV_COLS})


def aggregate(reps: list, wall_s: float, write_dev: bool = True, ridge: list | None = None,
              commit: str = "") -> dict:
    over = wall_s / 60.0 > TIMEBOX_MIN
    smooth_nc2 = {k: [{"r": i, "value": rp[k]["value"], "count": rp[k].get("count")} for i, rp in enumerate(reps)]
                  for k in ("NC2a", "NC2b")}
    if ridge is not None:                       # ruling R1: NC2 from the ridge toy
        reps = [dict(rp, NC2a=rd["NC2a"], NC2b=rd["NC2b"]) for rp, rd in zip(reps, ridge)]
    spec = {
        "NC2a": dict(key="NC2_lowtier_median_fidelity_max", kind="max", continuous=True,
                     definition="median EN-pred r2_oof_beyond_design over the toy low-tier plant genes (R2)"),
        "NC2b": dict(key="NC2_lowtier_md_strict_per_gene_max", kind="max", continuous=False,
                     definition="MD-strict rows (R2, BH within the low-tier universe) / number of low-tier genes"),
        "NC3": dict(key="NC3_self_null_dependent_frac_max", kind="max", continuous=False,
                    definition="mean frac_dependent (carrier median < -0.5) over testable pairs of the answer-key "
                               "self-null set (planted self-null gene: DEL_g and DEL locus mates), engine self_null (R1)"),
        "NC5": dict(key="NC5_false_sensitivity_rate_max", kind="max", continuous=False,
                    definition="max over the four arms of (#bypass-resistance/N toy CGT tests hit at rank_of_feature_D "
                               "<= 200 with beta < 0) / (#such tests ranked in the arm), engine cgt_scores (R1)"),
        "NC6": dict(key="NC6_lineage_tf_class3_frac_min", kind="min", continuous=False,
                    definition="fraction of toy lineage_only plant genes with engine gene_classes failure_class 3 (R1)"),
    }
    utc = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    out_nc, dev_rows = {}, []
    for nc, sp in spec.items():
        prov = float(NC[sp["key"]])
        vals = [rp[nc]["value"] for rp in reps]
        ns = [rp[nc]["n"] for rp in reps]
        fin = [v for v in vals if v is not None]
        tb = (max(fin) if sp["kind"] == "max" else min(fin)) if len(fin) == len(vals) and fin else None
        lim = prov if sp["kind"] == "max" else 1.0 - prov
        f2_fail = [] if sp["continuous"] else [i for i, n in enumerate(ns) if n == 0 or 1.0 / n > lim + 1e-12]
        cnts = [rp[nc].get("count") for rp in reps]
        degen = False                           # ruling R2
        if fin and len(fin) == len(vals):
            degen = (float(np.ptp(fin)) == 0.0 and fin[0] in (0.0, 1.0)) or \
                (None not in cnts and all(c == 0 for c in cnts))
        faithful = tb is not None and not f2_fail and not degen
        if over:
            final, src = prov, "provisional"
            reason = f"timebox exceeded ({wall_s / 60:.1f} min > {TIMEBOX_MIN:.0f} min)"
        elif not faithful:
            final, src = prov, "provisional"
            reason = ("F2 resolution: per-replicate denominator n = " + str(ns) + f", 1/n > {lim:g}" if f2_fail
                      else (f"R2 degeneracy: values {fin}, underlying counts {cnts}" if degen
                            else "toy statistic undefined in a replicate"))
        else:
            final, src, reason = tb, "toy", ""
        out_nc[nc] = {"design_key": sp["key"], "definition": sp["definition"], "bound_kind": sp["kind"],
                      "per_replicate": [{"r": i, "seed": rp["seed"], "value": rp[nc]["value"], "n": rp[nc]["n"]}
                                        for i, rp in enumerate(reps)],
                      "underlying_counts": cnts,
                      "toy": "ridge (ruling R1)" if (nc in ("NC2a", "NC2b") and ridge is not None) else "smooth",
                      "toy_bound": tb, "provisional": prov,
                      "faithful": {"F1": True, "F2": not f2_fail, "F2_failing_replicates": f2_fail,
                                   "R2_degenerate": degen},
                      "final_value": final, "final_source": src, "reason": reason}
        if src == "provisional" and not over:
            dev_rows.append({"id": f"DEV-TOY-{nc}-{commit}", "utc": utc, "stage": "toy_run_bounds",
                             "item": f"negative_controls_s7.{sp['key']}", "provisional": prov,
                             "toy_value": "" if tb is None else tb, "final_value": prov,
                             "reason": f"no faithful toy analogue: {reason} (docs/toy_nc_bounds.md §2 and addendum R2; engine {commit})",
                             "source": "scripts/toy_nc_bounds.py", "agent": "toy-verifier"})
    # NC4: the toy pass bounds (not re-derived)
    a_max, w_min = int(PASS["fold_offset_selected_after_centring_max"]), int(PASS["fold_offset_selected_without_centring_min"])
    pr4 = [{"r": i, "seed": rp["seed"], "after_centring": rp["NC4"]["after_centring"],
            "without_centring": rp["NC4"]["without_centring"], "n": rp["NC4"]["n"],
            "pass_after": rp["NC4"]["after_centring"] <= a_max, "pass_without": rp["NC4"]["without_centring"] >= w_min,
            "freq_after": rp["NC4"]["freq_after"], "freq_without": rp["NC4"]["freq_without"]}
           for i, rp in enumerate(reps)]
    out_nc["NC4"] = {"design_key": "NC4_rule", "definition": "fold-artefact features EN-pred-selected after centring "
                     "(R1) / without centring (R3); bounds = toy_s13.pass fold_offset_*",
                     "per_replicate": pr4, "toy_bound": {"after_centring_max_observed": max(p["after_centring"] for p in pr4),
                                                         "without_centring_min_observed": min(p["without_centring"] for p in pr4)},
                     "provisional": {"after_centring_max": a_max, "without_centring_min": w_min},
                     "final_value": {"after_centring_max": a_max, "without_centring_min": w_min},
                     "final_source": "toy pass bounds (not re-derived)",
                     "replicates_passing": {"after": sum(p["pass_after"] for p in pr4),
                                            "without": sum(p["pass_without"] for p in pr4)}}
    if over:
        dev_rows.append({"id": f"DEV-TOY-TIMEBOX-{commit}", "utc": utc, "stage": "toy_run_bounds",
                         "item": "negative_controls_s7 NC2-NC6", "provisional": "all", "toy_value": "reported",
                         "final_value": "provisional",
                         "reason": f"5 replicates took {wall_s / 60:.1f} min > toy_run_bounds {TIMEBOX_MIN:.0f} min",
                         "source": "scripts/toy_nc_bounds.py", "agent": "toy-verifier"})
    if write_dev and dev_rows:
        _append_deviation(dev_rows)
    extras = {f"r{i}": {"NC2b": {k: v for k, v in rp["NC2b"].items() if k != "value"},
                        "NC2a_per_gene": rp["NC2a"]["per_gene_fidelity"],
                        "NC3": {k: v for k, v in rp["NC3"].items() if k != "value"},
                        "NC5": rp["NC5"]["per_arm"], "NC6_classes": rp["NC6"]["classes"],
                        "timing": rp.get("timing")} for i, rp in enumerate(reps)}
    if ridge is not None:
        for k in ("NC2a", "NC2b"):
            out_nc[k]["smooth_values_reported"] = smooth_nc2[k]
    return {"nc": out_nc, "replicate_details": extras, "deviation_rows_written": dev_rows if write_dev else [],
            "timebox_min": TIMEBOX_MIN, "wall_clock_min": wall_s / 60.0, "timebox_exceeded": over}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--replicate", type=int, default=None)
    ap.add_argument("--cache", type=Path, required=True)
    ap.add_argument("--aggregate-only", action="store_true")
    ap.add_argument("--wall-s", type=float, default=None)
    ap.add_argument("--no-deviations", action="store_true")
    ap.add_argument("--stats-only", action="store_true", help="recompute rep*.json from rep*.pkl")
    ap.add_argument("--ridge-nc2-only", action="store_true", help="only the ruling-R1 ridge NC2 run of --replicate")
    a = ap.parse_args(argv)
    a.cache.mkdir(parents=True, exist_ok=True)
    if a.replicate is not None:
        if a.ridge_nc2_only:
            run_ridge_nc2(a.replicate, a.cache, WORKERS_PER_REP)
        elif a.stats_only:
            stats_from_cache(a.replicate, a.cache)
        else:
            run_replicate(a.replicate, a.cache, WORKERS_PER_REP)
        return 0
    t_launch = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    if not a.aggregate_only:
        env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
                   PYTHONPATH=str(ROOT / "src"))
        t0 = time.perf_counter()
        procs = []
        for r in range(N_REP):
            log = open(a.cache / f"rep{r}.log", "w", encoding="utf-8")
            procs.append((subprocess.Popen([sys.executable, "-X", "utf8", __file__, "--replicate", str(r),
                                            "--cache", str(a.cache)], env=env, stdout=log, stderr=subprocess.STDOUT,
                                           cwd=str(ROOT)), log))
        codes = []
        for p, log in procs:
            codes.append(p.wait())
            log.close()
        wall = time.perf_counter() - t0
        (a.cache / "wall.json").write_text(json.dumps({"wall_s": wall, "codes": codes, "launched_utc": t_launch}),
                                           encoding="utf-8")
        if any(codes):
            raise SystemExit(f"replicate exit codes {codes}; see {a.cache}/rep*.log")
    else:
        wall = float(a.wall_s)
    reps = [json.loads((a.cache / f"rep{r}.json").read_text(encoding="utf-8")) for r in range(N_REP)]
    commit = gitutil.head(P.PROJECT)[:7]
    rf = [a.cache / f"rep{r}_ridge_nc2.json" for r in range(N_REP)]
    ridge = [json.loads(f.read_text(encoding="utf-8")) for f in rf] if all(f.exists() for f in rf) else None
    agg = aggregate(reps, wall, write_dev=not a.no_deviations, ridge=ridge, commit=commit)
    res = {"description": "NC2-NC6 toy bounds (DESIGN v3.2 §7.3); definitions: docs/toy_nc_bounds.md",
           "definitions_doc": "docs/toy_nc_bounds.md",
           "definitions_doc_sha256": hashlib.sha256(DOC.read_bytes()).hexdigest(),
           "seeds": [seed_of(r) for r in range(N_REP)], "scenario": "truth", "predictor": "smooth",
           "launched_utc": t_launch, "workers": {"replicates_concurrent": N_REP, "engine_workers_each": WORKERS_PER_REP},
           "commit": gitutil.head(P.PROJECT), "config_sha256": PR.config_hashes(),
           "toy_answer_key_sha256": ts.toy_config_sha256(), **agg}
    out_path = OUT.with_name(f"nc_bounds_{commit}.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    print(json.dumps({k: {kk: v.get(kk) for kk in ("toy_bound", "provisional", "final_value", "final_source", "reason")}
                      for k, v in res["nc"].items()}, indent=1, default=float))
    print(f"wall {wall / 60:.1f} min (timebox {TIMEBOX_MIN:.0f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
