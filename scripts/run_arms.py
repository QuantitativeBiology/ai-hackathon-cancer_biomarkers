"""Run the four benchmark arms + MD-strict on real D data (DESIGN v3.2 §4-§5).

    python -B scripts/run_arms.py --pilot            # 50-gene timing pilot -> results/timing/surrogate_pilot.json
    python -B scripts/run_arms.py                    # full run (refused before tag freeze-bounds)
    python -B scripts/run_arms.py --universe ctrl    # U_ctrl genes (configs/u_ctrl_genes.txt; BH within U_ctrl)
    python -B scripts/run_arms.py --universe lowtier # NC2 low-tier genes (configs/lowtier_genes.txt; DESIGN §7.3)

Outputs (full run, eligible): cache/benchmark/arms/part-*.parquet (bulk; sha256 -> results/bulk_manifest.csv),
results/benchmark/{selected.csv, lambda.csv, p_cut.json, summary.json}, results/md/{md_global, md_global_all,
md_strict, md_strict_all, en_obs_strict}.csv, results/timing/{s1_obs,s1_pred,en_obs,en_pred,md_strict}.json.
``--universe ctrl`` writes the same files under ctrl-only roots (output_paths): cache/benchmark/arms_ctrl/,
results/benchmark_ctrl/ (timing in results/benchmark_ctrl/timing/, never budgeted) and results/md_ctrl/.
EN-pred never reads observed effects (it receives predictions and the NaN mask only).
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from biomark import arms as A  # noqa: E402
from biomark import gitutil  # noqa: E402
from biomark import params as PR  # noqa: E402
from biomark import paths as P  # noqa: E402
from biomark import pipeline as PL  # noqa: E402

FREEZE_BOUNDS = "freeze-bounds"


def sha(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


UNIVERSES = ("eligible", "ctrl", "lowtier")


def output_paths(universe: str = "eligible") -> dict:
    """Every location a full run writes, per universe. The ctrl universe (U_ctrl, configs/u_ctrl_genes.txt) writes
    only under ctrl-specific roots, so it never overwrites an eligible output: cache/benchmark/arms_ctrl/ (+ its stage
    dir), results/benchmark_ctrl/ (selected, lambda, p_cut, summary, timing/) and results/md_ctrl/. The one shared
    file is results/bulk_manifest.csv, keyed by path (ctrl rows have ctrl-specific paths)."""
    if universe not in UNIVERSES:
        raise ValueError(f"universe must be one of {UNIVERSES}")
    rt = PR.runtime()["blocked"]
    sfx = "" if universe == "eligible" else f"_{universe}"
    rb = P.RESULTS / f"benchmark{sfx}"
    return {"dataset": P.PROJECT / (rt["dataset_dir"] + sfx), "stage": P.PROJECT / (rt["stage_dir"] + sfx),
            "benchmark": rb, "md": P.RESULTS / f"md{sfx}",
            "timing": None if universe == "eligible" else rb / "timing",
            "bulk_manifest": P.RESULTS / "bulk_manifest.csv"}


def universe_genes(universe: str) -> list | None:
    """None = the eligible genes (Data.from_real default); ctrl = configs/u_ctrl_genes.txt, asserted disjoint from the
    eligible set (DESIGN §2: arms ⊆ eligible; U_ctrl runs separately)."""
    if universe == "eligible":
        return None
    name, builder = {"ctrl": ("u_ctrl_genes.txt", "build_u_ctrl_genes.py"),
                     "lowtier": ("lowtier_genes.txt", "build_lowtier_genes.py")}[universe]
    f = P.CONFIGS / name
    if not f.exists():
        raise SystemExit(f"{f} missing: run scripts/{builder}")
    genes = [x.strip() for x in f.read_text(encoding="utf-8").splitlines() if x.strip()]
    el = pd.read_csv(P.CONFIGS / "eligible_genes.csv", dtype={"entrez": str})
    both = set(genes) & set(el.loc[el["eligible"].astype(bool), "entrez"])
    if both:
        raise SystemExit(f"{name} overlaps the eligible set: {sorted(both)}")
    return genes


def write_timing(stem: str, rec: dict, timing_dir=None, **kw) -> None:
    """results/timing/<stem>.json (schema of analyses.write_timing) for the eligible universe (``timing_dir`` None);
    the ctrl universe writes to results/benchmark_ctrl/timing/ so budget_check.py never budgets it."""
    from biomark import analyses as AN
    AN.write_timing(stem, rec, timing_dir=timing_dir, **kw)


def _epoch(t: float) -> dt.datetime:
    return dt.datetime.fromtimestamp(t, tz=dt.timezone.utc)


def ensure_nan_mask(genes, path: Path) -> Path:
    """Write the effect NaN mask (D lines x genes, True = NaN) if absent. Only the mask leaves this function."""
    if path.exists():
        have = pd.read_parquet(path, columns=[])
        cols = pd.read_parquet(path).columns
        if set(genes) <= set(cols) and len(have):
            return path
    from biomark import guard, p1
    dev = list(p1.dev_ids())
    m = guard.load_effects(dev, scope="dev", genes=list(genes)).isna()
    path.parent.mkdir(parents=True, exist_ok=True)
    m.to_parquet(path)
    return path


def bulk_manifest(paths) -> None:
    f = P.RESULTS / "bulk_manifest.csv"
    rows = pd.read_csv(f) if f.exists() else pd.DataFrame(columns=["path", "sha256", "bytes", "utc", "script"])
    for p in paths:
        rel = Path(p).resolve().relative_to(P.PROJECT).as_posix()
        rows = rows[rows["path"] != rel]
        rows.loc[len(rows)] = {"path": rel, "sha256": sha(p), "bytes": Path(p).stat().st_size,
                               "utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                               "script": "scripts/run_arms.py"}
    rows.to_csv(f, index=False)


def pilot(args) -> int:
    ep, npar = PR.en_params(), PR.null_params()
    el = pd.read_csv(P.CONFIGS / "eligible_genes.csv", dtype={"entrez": str})
    elig = list(el.loc[el["eligible"].astype(bool), "entrez"])
    rng = np.random.default_rng(int(PR.seeds()["subsample_seed"]))
    genes = sorted(rng.choice(np.asarray(elig, dtype=object), size=min(ep["pilot_genes"], len(elig)),
                              replace=False).tolist())
    main_mask = P.CACHE / "effect_nan_mask.parquet"
    mpath = main_mask if main_mask.exists() else ensure_nan_mask(genes, P.CACHE / "pilot" / "effect_nan_mask.parquet")
    t0 = time.perf_counter()
    data = PL.Data.from_real(genes=genes, nan_mask_path=mpath, save_subsamples=False)
    t_load = time.perf_counter() - t0
    tim = {}
    t0 = time.perf_counter()
    s1 = A.run_s1(data, A.S1_ARMS)                      # serial: per-gene CPU incl. pattern-group overhead
    tim["s1_both"] = time.perf_counter() - t0
    U = s1[s1["arm"] == "S1-obs"][["gene", "feature"]]
    w = args.workers or 1
    chunk = max(1, -(-len(genes) // w))                  # spread the pilot genes over all workers
    res = {}
    for arm in A.EN_ARMS:
        t0 = time.perf_counter()
        tab, lam, det = A.run_en(data, (arm,), U=U, workers=w, chunk=chunk)
        tim[arm] = time.perf_counter() - t0
        res[arm] = (tab, lam, det)
    hits = sum(d.get("cache_hits", 0) for a in res for d in res[a][2])
    miss = sum(d.get("cache_misses", 0) for a in res for d in res[a][2])
    n_all = len(elig)
    full_workers = npar["all_workers"]
    cpu_en = sum(d["wall_s"] for a in res for d in res[a][2] if "wall_s" in d)   # summed worker seconds, both arms
    per_gene_en = cpu_en / len(genes)
    proj_en_min = per_gene_en * n_all / full_workers / npar["sec_per_min"]
    proj_s1_min = tim["s1_both"] / len(genes) * n_all / full_workers / npar["sec_per_min"]
    # fixed-threshold null cost: one permuted draw = S1 both arms + fixed-lambda EN both arms per null gene
    from biomark import null as NL
    lamd = {a: {r.gene: {"G": r.lam_G, "X": r.lam_X} for r in res[a][1].itertuples()} for a in res}
    obs = {"p_cut": A.observed_cutoffs(A.finish_s1(data, s1)), "lam": pd.concat([res[a][1] for a in res]),
           "arms": None, "p_cut_sc": 0.0}
    t0 = time.perf_counter()
    NL._draw(data, obs, genes, 0)
    t_draw = time.perf_counter() - t0
    per_gene_draw = t_draw / len(genes)
    draw_min = per_gene_draw * npar["gene_subset_n"] / npar["sec_per_min"]
    del lamd
    rec = {"n_genes": len(genes), "n_features": int(data.X.shape[1]), "n_lines": int(data.X.shape[0]),
           "workers": w, "chunk": chunk, "worker_cpu_s_en": cpu_en, "load_s": t_load, "stage_wall_s": tim,
           "wall_s": float(t_load + sum(tim.values()) + t_draw), "cache_hit_rate": hits / max(1, hits + miss),
           "cache_hits": hits, "cache_misses": miss, "n_eligible": n_all, "projection_workers": full_workers,
           "projected_wall_min": {"s1_obs_pred": proj_s1_min, "en_obs_pred": proj_en_min,
                                  "null_one_draw_single_worker": draw_min,
                                  "null_B_base": draw_min * npar["B_base"] / npar["workers"],
                                  "null_B_raised": draw_min * npar["B_raised"] / npar["workers"]},
           "null_B_rule": "B = B_raised iff B_raised x wall(draw 1) <= wall_cap_min (single draw wall)",
           "K": ep["K"], "pilot_max_projected_wall_min": ep["pilot_max_min"],
           "K_rule": "projected EN wall > pilot_max_projected_wall_min => K = sis_K_fallback",
           "K_decision": ep["K_fallback"] if proj_en_min > ep["pilot_max_min"] else ep["K"],
           "universe_rows": int(len(U)), "selected": {a: int(res[a][0]["selected"].sum()) for a in res},
           "fidelity_median": float(np.nanmedian(res["EN-pred"][1]["r2_oof_beyond_design"]))}
    write_timing("surrogate_pilot", rec)
    print(json.dumps(rec, indent=1, default=float))
    return 0


def full(args) -> int:
    """Streamed (blocked) arms on the full universe: S1 + EN per pattern chunk to a parquet dataset, global steps on
    slim columns (biomark.blocked). Peak memory ~ parent (slim columns) + workers x per-worker RSS, guarded."""
    from biomark import analyses as AN
    from biomark import blocked as BL
    if not gitutil.tag_exists(P.PROJECT, FREEZE_BOUNDS) and not args.force_before_bounds:
        raise SystemExit(f"tag {FREEZE_BOUNDS} missing: full arms are launched by the PI after freeze-bounds")
    genes = universe_genes(args.universe)
    out = output_paths(args.universe)
    data = PL.Data.from_real(genes=genes)
    ws1, g1 = AN.guard_workers(args.workers_s1, "s1", log=print)
    wen, g2 = AN.guard_workers(args.workers_en, "en", log=print)
    res = BL.run_streamed(data, A.ARMS, out_dir=out["dataset"], stage_dir=out["stage"],
                          workers_s1=ws1, workers_en=wen, log=print)
    summ = write_outputs(res, data, args.universe, {"s1": g1, "en": g2}, ws1)
    print(json.dumps(summ, indent=1, default=str))
    return 0


def write_outputs(res: dict, data, universe: str, guards: dict, ws1: int) -> dict:
    """Write the results of blocked.run_streamed for ``universe`` to output_paths(universe) only."""
    from biomark import endpoint_io as IO
    out = output_paths(universe)
    out_dir, rb, rm, td = out["dataset"], out["benchmark"], out["md"], out["timing"]
    tim = res["timings"]
    rb.mkdir(parents=True, exist_ok=True)
    rm.mkdir(parents=True, exist_ok=True)
    bulk_manifest([out_dir / p for p in res["parts"]])
    IO.read_arms(selected_only=True, dataset_dir=out_dir).to_csv(rb / "selected.csv", index=False)
    res["lam"].to_csv(rb / "lambda.csv", index=False)
    pc = dict(res["p_cut"], **({"MD-strict": res["p_cut_sc"]} if "p_cut_sc" in res else {}))
    (rb / "p_cut.json").write_text(json.dumps(pc, indent=1), encoding="utf-8")
    res["md_global"].to_csv(rm / "md_global.csv", index=False)          # A7: the catalogue
    res["md_global_all"].to_csv(rm / "md_global_all.csv", index=False)
    res["md_strict"].to_csv(rm / "md_strict.csv", index=False)          # frozen secondary
    res["md_strict_all"].to_csv(rm / "md_strict_all.csv", index=False)
    if "en_obs_strict" in res:
        res["en_obs_strict"].to_csv(rm / "en_obs_strict.csv", index=False)
    summ = {"universe": universe, "n_genes": len(data.genes), "n_rows_arms": int(res["n_rows"]),
            "n_parts": len(res["parts"]), "dataset_dir": str(out_dir.relative_to(P.PROJECT).as_posix()),
            "bh_families": "within this universe (all rows of the dataset)",
            "selected": res["selected"], "n_md_global": int(len(res["md_global"])), "n_md_strict": int(len(res["md_strict"])),
            "workers": guards, "config_sha256": PR.config_hashes(),
            "subsamples_sha256": sha(P.CONFIGS / "subsamples.json"),
            "dictionary_built_utc": data.extras.get("dictionary_built_utc"), "commit": gitutil.head(P.PROJECT)}
    if universe == "ctrl":
        summ["u_ctrl_genes_sha256"] = sha(P.CONFIGS / "u_ctrl_genes.txt")
    elif universe == "lowtier":
        summ["lowtier_genes_sha256"] = sha(P.CONFIGS / "lowtier_genes.txt")
    (rb / "summary.json").write_text(json.dumps(summ, indent=1, default=str), encoding="utf-8")
    n = len(data.genes)
    g1, g2 = guards.get("s1"), guards.get("en")
    s0, s1e = tim["s1_window"]
    for stem in ("s1_obs", "s1_pred"):          # one scan, two stems: same window and run_id (counted once)
        write_timing(stem, {"note": "S1-obs and S1-pred computed in one streamed scan", "n_genes": n, "memory_guard": g1},
                     td, workers=ws1, start_utc=_epoch(s0), end_utc=_epoch(s1e))
    for arm, stem in (("EN-obs", "en_obs"), ("EN-pred", "en_pred")):
        win = [d for d in tim["en_walls"] if d["arm"] == arm][0]
        write_timing(stem, {"n_genes": n, "memory_guard": g2}, td, workers=int(win["workers"]),
                     start_utc=_epoch(win["start_epoch"]), end_utc=_epoch(win["end_epoch"]))
    g0, g1e = tim["global_window"]
    write_timing("md_strict", {"note": "global BH / ranks / dedup / MD lists and final parts (blocked.run_streamed)"},
                 td, workers=1, start_utc=_epoch(g0), end_utc=_epoch(g1e))
    register_outputs(universe)
    return summ


def registry_outputs(universe: str, existing_only: bool = True) -> list:
    """The results files of a full run of ``universe`` (repo-relative; the bulk parts are in bulk_manifest.csv)."""
    sfx = "" if universe == "eligible" else f"_{universe}"
    files = [f"results/benchmark{sfx}/{f}" for f in ("selected.csv", "lambda.csv", "p_cut.json", "summary.json")]
    files += [f"results/md{sfx}/{f}.csv" for f in ("md_global", "md_global_all", "md_strict", "md_strict_all",
                                                    "en_obs_strict")]
    assert all((P.PROJECT / f).parent in (output_paths(universe)["benchmark"], output_paths(universe)["md"])
               for f in files)
    return [f for f in files if not existing_only or (P.PROJECT / f).exists()]


def register_outputs(universe: str) -> None:
    """results/registry.csv row (status ok) for the files this run wrote (DESIGN §13.5: citable outputs)."""
    from biomark import analyses as AN
    script = "scripts/run_arms.py" + ("" if universe == "eligible" else f" --universe {universe}")
    AN.append_registry(script, "ok", registry_outputs(universe), registry=P.RESULTS / "registry.csv")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", action="store_true")
    ap.add_argument("--workers", type=int, default=None, help="pilot only")
    ap.add_argument("--workers-s1", type=int, default=None, help="S1 pool (default runtime_params workers_default.s1)")
    ap.add_argument("--workers-en", type=int, default=None, help="EN pools (default runtime_params workers_default.en)")
    ap.add_argument("--universe", choices=list(UNIVERSES), default="eligible")
    ap.add_argument("--force-before-bounds", action="store_true", help="PI override of the freeze-bounds check")
    a = ap.parse_args(argv)
    return pilot(a) if a.pilot else full(a)


if __name__ == "__main__":
    sys.exit(main())
