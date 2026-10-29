"""XLIB secondary (DESIGN v3.2 §6.3, primary_endpoints.XLIB_secondary): S1-obs and EN-obs only.

    python -B scripts/run_xlib.py [--workers-s1 N --workers-en N]      # background, after the fixed-threshold null
    python -B scripts/run_xlib.py --reuse-arms [--expect-parts P --expect-rows R]   # reuse a finished arms dataset

Discovery = D lines whose library string is exactly the discovery library; replication = D lines with no such
component; others in neither. The two arms are refitted on discovery lines with their own λ CV (same params,
Part 1 folds and subsamples restricted to discovery lines); pool, power (discovery-set shrunk effects, replication
molecular data), k and lists follow the P1 rules; replication tests use the same confounder model and drops; paired
patient bootstrap over replication lines; contrast EN-obs − S1-obs with the P1 branches. u_rule all_D uses U and
eligibility as built on all D; discovery_only recomputes eligibility on discovery lines (no refit) and restricts the
pool. Lists (pool, k, de-duplication, ranks) are built for the XLIB arms only. ``--reuse-arms`` skips the streamed
arms when cache/benchmark/xlib_arms/ is complete (blocked.verify_dataset: its _manifest.json, or for a dataset written
before manifests existed the part / row counts of that run's log via --expect-parts / --expect-rows); the summary then
records reused_arms = true. Writes results/benchmark/replication_xlib.csv (+ xlib_summary.json),
cache/benchmark/xlib_arms/ and results/timing/{xlib_s1_obs, xlib_en_obs}.json. Refuses unless HEAD descends from
freeze-bounds.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from biomark import analyses as AN  # noqa: E402
from biomark import arms as A  # noqa: E402
from biomark import endpoint_io as IO  # noqa: E402
from biomark import endpoints as EP  # noqa: E402
from biomark import params as PR  # noqa: E402
from biomark import paths as P  # noqa: E402
from biomark import pipeline as PL  # noqa: E402
from biomark import stats as ST  # noqa: E402

XARMS = tuple(PR.dp()["replication_s6"]["xlib_arms"])


def reuse_dataset(dD, stage: Path, expect_parts=None, expect_rows=None) -> dict:
    """--reuse-arms: verify cache/benchmark/xlib_arms is a complete streamed dataset of these XLIB arms / genes /
    discovery lines (blocked.verify_dataset: its manifest, or for a dataset written before manifests the run's
    logged part and row counts) and write the manifest if it was missing."""
    from biomark import blocked as BL
    try:
        man = BL.verify_dataset(IO.XLIB_ARMS_DIR, XARMS, list(dD.genes), list(dD.lines), stage_dir=stage,
                                expect_parts=expect_parts, expect_rows=expect_rows)
    except ValueError as e:
        raise SystemExit(f"--reuse-arms refused: {e}")
    f = IO.XLIB_ARMS_DIR / BL.MANIFEST
    if not f.exists():
        IO.write_json(f, dict(man, written_by="scripts/run_xlib.py --reuse-arms", utc=AN.utc_now().isoformat()))
    print(f"reusing {IO.XLIB_ARMS_DIR}: {len(man['parts'])} parts, {man['n_rows']} rows "
          f"({man.get('verified_from', 'manifest')})")
    return man


def part_window(d: Path) -> tuple:
    mt = [p.stat().st_mtime for p in Path(d).glob("part-*.parquet")]
    return min(mt), max(mt)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers-s1", type=int, default=None)
    ap.add_argument("--workers-en", type=int, default=None)
    ap.add_argument("--reuse-arms", action="store_true",
                    help="skip run_streamed and use the complete dataset in cache/benchmark/xlib_arms")
    ap.add_argument("--expect-parts", type=int, default=None,
                    help="--reuse-arms on a dataset without _manifest.json: part count from the run's log")
    ap.add_argument("--expect-rows", type=int, default=None,
                    help="--reuse-arms on a dataset without _manifest.json: row count from the run's log")
    a = ap.parse_args(argv)
    EP.require_descends("freeze-bounds")
    from biomark import blocked as BL
    ws1, g1 = AN.guard_workers(a.workers_s1, "s1", log=print)
    wen, g2 = AN.guard_workers(a.workers_en, "en", log=print)
    data = PL.Data.from_real()
    split = AN.xlib_split(IO.library_strings(list(data.lines)))
    disc = (split == "discovery").to_numpy()
    repl_ids = list(data.lines[(split == "replication").to_numpy()])
    dD = IO.subset_data(data, disc)
    # streamed arms on discovery lines (same orchestration as run_arms; XLIB arms S1-obs and EN-obs only)
    stage = P.CACHE / "benchmark" / "_stage_xlib"
    if a.reuse_arms:
        man = reuse_dataset(dD, stage, a.expect_parts, a.expect_rows)
        win = man.get("timings", {}).get("s1_window") or part_window(IO.XLIB_ARMS_DIR)
        timing_note = ("reused arms dataset: " + ("S1 window from its manifest" if man.get("timings") else
                       "window = part write times of the original streamed run (S1/EN compute not re-timed; "
                       "see that run's log)"))
    else:
        res = BL.run_streamed(dD, XARMS, out_dir=IO.XLIB_ARMS_DIR, stage_dir=stage, workers_s1=ws1, workers_en=wen,
                              log=print)
        win, timing_note, man = res["timings"]["s1_window"], "", None
    s1_start, s1_end = [AN.dt.datetime.fromtimestamp(x, tz=AN.dt.timezone.utc) for x in win]
    en_start = AN.utc_now() if a.reuse_arms else s1_end
    arms_tab, s1o = IO.pool_inputs(dataset_dir=IO.XLIB_ARMS_DIR)
    genes = list(dD.genes)
    disc_ids = list(dD.lines)
    nan_D = pd.read_parquet(P.CACHE / "effect_nan_mask.parquet").reindex(index=list(data.lines), columns=genes)
    Yd = dD.target("obs")
    sig = IO.sigma_resid_genes(Yd, dD.rows_ok, ST.CovDesign(dD.cov_df), genes)
    pool, pj = IO.make_pool(arms_tab, s1o, sig, disc_ids, repl_ids, nan_D, arms=XARMS)   # §6.3: XLIB arms only
    # u_rule discovery_only
    el = pd.read_csv(P.CONFIGS / "eligible_genes.csv", dtype={"entrez": str})
    strata = (dD.cov_df["lineage"].astype(str) + "|" + dD.cov_df["library_group"].astype(str)).to_numpy()
    ok_disc = AN.discovery_eligibility(el, dD.pred, np.where(dD.rows_ok, Yd, np.nan), strata, genes)
    keep = set(ok_disc.index[ok_disc.to_numpy()])
    Pd = pool[pool["powered"] & pool["gene"].isin(keep)]
    pj_disc = {"pool_size": int(len(Pd)), "k": EP.p1_k(len(Pd)), "lists": {"G": EP.build_lists(arms_tab, Pd, XARMS)},
               "strata": {}}
    for name, m in EP.strata_of(Pd["access_class"]).items():
        sub = Pd[m.to_numpy()]
        pj_disc["strata"][name] = {"pool_size": int(len(sub)), "k": EP.p1_k(len(sub))}
        pj_disc["lists"][name] = EP.build_lists(arms_tab, sub, XARMS)
    # replication on non-Avana D lines (observed D effects, scope dev)
    Y_all = pd.DataFrame(data.target("obs"), index=data.lines, columns=genes)
    feats = sorted({f for pjx in (pj, pj_disc) for L in pjx["lists"].values() for arm in XARMS
                    for _, f, _ in L.get(arm, {}).get("rows", [])})
    fl = IO.feature_meta()
    binary = (fl["d_type"] == "binary").to_dict()
    drops = {(str(g), f): (d if isinstance(d, str) else "") for g, f, d in zip(s1o["gene"], s1o["feature"], s1o["dropped_cov"])}
    fitter = IO.fitter_for(repl_ids, genes, feats, Y_all.loc[repl_ids])
    rows, summ = [], {"discovery_lines": int(disc.sum()), "replication_lines": len(repl_ids), "u_rules": {}}
    for u_rule, pjx in (("all_D", pj), ("discovery_only", pj_disc)):
        lists, k_of = {}, {}
        for stratum, L in pjx["lists"].items():
            for arm in XARMS:
                v = L.get(arm, {"rows": [], "k": 0})
                lists[f"{stratum}|{arm}"] = [(str(g), f, float(s)) for g, f, s in v["rows"]]
                k_of[f"{stratum}|{arm}"] = int(v["k"])
        point = {nm: float(np.mean([EP.replicated(fitter.fit(g, f, drops.get((g, f), ""), sign=s), s,
                                                  binary=binary.get(f, True)) for g, f, s in L])) if L else np.nan
                 for nm, L in lists.items()}
        bt = EP.bootstrap_rates(fitter, lists, drops, binary)
        summ["u_rules"][u_rule] = {"pool_size": pjx["pool_size"], "k": pjx["k"]}
        for stratum in pjx["lists"]:
            for c in EP.p1_contrasts(point, bt["rates"], {nm: k_of[nm] for nm in lists if nm.startswith(stratum + "|")},
                                     prefix=f"{stratum}|"):
                if c["contrast"] != "multivariate_adds_value":
                    continue
                rows.append({"u_rule": u_rule, "stratum": stratum, **{f"rate_{arm}": point[f"{stratum}|{arm}"] for arm in XARMS},
                             **c, "frac_degenerate_EN-obs": float(bt["frac_degenerate"][f"{stratum}|EN-obs"].mean())})
    out = P.RESULTS / "benchmark"
    pd.DataFrame(rows).to_csv(out / "replication_xlib.csv", index=False)
    summ["caveat"] = PR.endpoints()["XLIB_secondary"]["caveat"]
    summ["reused_arms"] = bool(a.reuse_arms)
    if man is not None:
        summ["arms_dataset"] = {"n_rows": man.get("n_rows"), "n_parts": len(man.get("parts", {})),
                                "verified_from": man.get("verified_from", "manifest")}
    IO.write_json(out / "xlib_summary.json", summ)
    AN.write_timing("xlib_s1_obs", {"n_genes": len(genes), "n_discovery": int(disc.sum()), "memory_guard": g1,
                                    "reused_arms": bool(a.reuse_arms), **({"note": timing_note} if timing_note else {})},
                    workers=ws1, start_utc=s1_start, end_utc=s1_end)
    AN.write_timing("xlib_en_obs", {"note": "EN-obs refit + pool + replication bootstrap", "memory_guard": g2},
                    workers=wen, start_utc=en_start)
    AN.append_registry("scripts/run_xlib.py", "ok", ["results/benchmark/replication_xlib.csv", "results/benchmark/xlib_summary.json"])
    print(pd.DataFrame(rows).to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
