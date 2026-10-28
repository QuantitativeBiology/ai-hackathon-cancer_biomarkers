"""Run the four benchmark arms + MD-strict on real D data (DESIGN v3.2 §4-§5).

    python -B scripts/run_arms.py --pilot            # 50-gene timing pilot -> results/timing/surrogate_pilot.json
    python -B scripts/run_arms.py                    # full run (refused before tag freeze-bounds)
    python -B scripts/run_arms.py --universe ctrl    # U_ctrl genes (non-eligible controls / CGT targets)

Outputs (full run): cache/benchmark/{arms,s1,en}.parquet (bulk; sha256 -> results/bulk_manifest.csv),
results/benchmark/{selected.csv, lambda.csv, p_cut.json, summary.json}, results/md/{md_strict.csv,
md_strict_all.csv, en_obs_strict.csv}, results/timing/{s1_obs,s1_pred,en_obs,en_pred,md_strict}.json.
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


def write_timing(stem: str, rec: dict, tag: str = "") -> None:
    """results/timing/<stem>.json for the eligible universe (stems of tiers.json stage_job_map); other universes
    write to results/benchmark/timing<tag>/ so budget_check.py never sees an unmapped stem."""
    d = P.PROJECT / PR.pp()["timing_dir"] if not tag else P.RESULTS / "benchmark" / f"timing{tag}"
    d.mkdir(parents=True, exist_ok=True)
    rec = {"job": stem, "utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"), **rec}
    (d / f"{stem}.json").write_text(json.dumps(rec, indent=1, default=float), encoding="utf-8")


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
    hits = sum(d["cache_hits"] for a in res for d in res[a][2])
    miss = sum(d["cache_misses"] for a in res for d in res[a][2])
    n_all = len(elig)
    full_workers = npar["all_workers"]
    cpu_en = sum(d["wall_s"] for a in res for d in res[a][2])            # summed worker seconds, both arms
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
           "workers_used": w, "chunk": chunk, "worker_cpu_s_en": cpu_en, "load_s": t_load, "wall_s": tim, "cache_hit_rate": hits / max(1, hits + miss),
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
    if not gitutil.tag_exists(P.PROJECT, FREEZE_BOUNDS) and not args.force_before_bounds:
        raise SystemExit(f"tag {FREEZE_BOUNDS} missing: full arms are launched by the PI after freeze-bounds")
    genes = None
    if args.universe == "ctrl":
        f = P.CONFIGS / "u_ctrl_genes.txt"
        genes = [x.strip() for x in f.read_text(encoding="utf-8").splitlines() if x.strip()]
    data = PL.Data.from_real(genes=genes)
    npar = PR.null_params()
    w = args.workers or npar["all_workers"]
    tim = {}
    res = A.run_arms(data, arms=A.ARMS, workers=w, timings=tim)
    tag = "" if args.universe == "eligible" else f"_{args.universe}"
    cb = P.CACHE / "benchmark"
    cb.mkdir(parents=True, exist_ok=True)
    rb, rm = P.RESULTS / "benchmark", P.RESULTS / "md"
    rb.mkdir(parents=True, exist_ok=True)
    rm.mkdir(parents=True, exist_ok=True)
    bulk = []
    for k in ("arms", "s1", "en"):
        f = cb / f"{k}{tag}.parquet"
        res[k].to_parquet(f)
        bulk.append(f)
    bulk_manifest(bulk)
    a = res["arms"]
    a[a["selected"]].to_csv(rb / f"selected{tag}.csv", index=False)
    res["lam"].to_csv(rb / f"lambda{tag}.csv", index=False)
    pc = dict(res["p_cut"], **({"MD-strict": res["p_cut_sc"]} if "p_cut_sc" in res else {}))
    (rb / f"p_cut{tag}.json").write_text(json.dumps(pc, indent=1), encoding="utf-8")
    res["md_strict"].to_csv(rm / f"md_strict{tag}.csv", index=False)
    res["md_strict_all"].to_csv(rm / f"md_strict_all{tag}.csv", index=False)
    if "en_obs_strict" in res:
        res["en_obs_strict"].to_csv(rm / f"en_obs_strict{tag}.csv", index=False)
    summ = {"universe": args.universe, "n_genes": len(data.genes), "n_rows_U": int(len(res["universe"])),
            "selected": a.groupby("arm")["selected"].sum().astype(int).to_dict(), "n_md_strict": int(len(res["md_strict"])),
            "config_sha256": PR.config_hashes(), "subsamples_sha256": sha(P.CONFIGS / "subsamples.json"),
            "dictionary_built_utc": data.extras.get("dictionary_built_utc"), "commit": gitutil.head(P.PROJECT)}
    (rb / f"summary{tag}.json").write_text(json.dumps(summ, indent=1, default=str), encoding="utf-8")
    n = len(data.genes)
    for stem, key in (("s1_obs", "s1_obs_pred"), ("s1_pred", "s1_obs_pred")):
        write_timing(stem, {"wall_s": tim[key], "note": "S1-obs and S1-pred computed in one scan", "n_genes": n}, tag)
    for arm, stem in (("EN-obs", "en_obs"), ("EN-pred", "en_pred")):
        cpu = sum(d["wall_s"] for d in tim["en_detail"] if d["arm"] == arm)
        write_timing(stem, {"worker_cpu_s": cpu, "wall_s_both_arms": tim["en_obs_pred"], "workers": w, "n_genes": n}, tag)
    write_timing("md_strict", {"wall_s": 0.0, "note": "computed from S1-obs and EN-pred tables inside run_arms"}, tag)
    print(json.dumps(summ, indent=1, default=str))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", action="store_true")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--universe", choices=["eligible", "ctrl"], default="eligible")
    ap.add_argument("--force-before-bounds", action="store_true", help="PI override of the freeze-bounds check")
    a = ap.parse_args(argv)
    return pilot(a) if a.pilot else full(a)


if __name__ == "__main__":
    sys.exit(main())
