"""Fixed-threshold joint-row permutation null on real D data (DESIGN v3.2 §7.1), after scripts/run_arms.py.

    python -B scripts/run_null.py                 # B = 20, raised to 50 iff 50 x wall(draw 1) <= 180 min; 6 workers
    python -B scripts/run_null.py --s1-vectorised # the additional 200 S1-only fixed-threshold permutations

Reads the observed operating points written by run_arms.py (cache/benchmark/arms.parquet, results/benchmark/
{lambda.csv, p_cut.json}, results/md/md_strict.csv) and configs/null_genes.txt (300-gene subset ⊆ eligible).
Writes results/null/{summary.json, fdr_fixed.csv, draws.csv} (or *_s1_vectorised*), results/timing/null_fixed.json.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd  # noqa: E402

from biomark import gitutil  # noqa: E402
from biomark import null as NL  # noqa: E402
from biomark import params as PR  # noqa: E402
from biomark import paths as P  # noqa: E402
from biomark import pipeline as PL  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--s1-vectorised", action="store_true")
    ap.add_argument("--force-before-bounds", action="store_true")
    a = ap.parse_args(argv)
    if not gitutil.tag_exists(P.PROJECT, "freeze-bounds") and not a.force_before_bounds:
        raise SystemExit("tag freeze-bounds missing: the null is launched by the PI after freeze-bounds")
    npar = PR.null_params()
    genes_file = P.CONFIGS / "null_genes.txt"
    if not genes_file.exists():
        raise SystemExit("configs/null_genes.txt missing (frozen at freeze-bounds)")
    genes = [x.strip() for x in genes_file.read_text(encoding="utf-8").splitlines() if x.strip()]
    el = pd.read_csv(P.CONFIGS / "eligible_genes.csv", dtype={"entrez": str})
    elig = set(el.loc[el["eligible"].astype(bool), "entrez"])
    assert set(genes) <= elig, "null_genes.txt must be a subset of the eligible genes (DESIGN §2)"
    data = PL.Data.from_real()
    rb = P.RESULTS / "benchmark"
    pc = json.loads((rb / "p_cut.json").read_text(encoding="utf-8"))
    from biomark import endpoint_io as IO
    res = {"arms": IO.read_arms(columns=["arm", "gene", "feature", "selected"], selected_only=True),
           "lam": pd.read_csv(rb / "lambda.csv", dtype={"gene": str}),
           "md_strict": pd.read_csv(P.RESULTS / "md" / "md_strict.csv", dtype={"gene": str}),
           "p_cut": {k: v for k, v in pc.items() if k != "MD-strict"}, "p_cut_sc": pc.get("MD-strict", 0.0)}
    from biomark import analyses as AN
    w, gnull = AN.guard_workers(a.workers, "null", log=print)
    out_dir = P.RESULTS / "null"
    out_dir.mkdir(parents=True, exist_ok=True)
    from biomark import analyses as AN
    t_start = AN.utc_now()
    t0 = time.perf_counter()
    if a.s1_vectorised:
        nul = NL.run_null(data, res, B=npar["s1_B"], workers=w, genes=genes, s1_only=True)
        sfx = "_s1_vectorised"
    else:
        nul = NL.run_null(data, res, B=npar["B_base"], workers=w, genes=genes, decide_B=True)
        sfx = ""
    wall = time.perf_counter() - t0
    nul["draws"].to_csv(out_dir / f"draws{sfx}.csv", index=False)
    nul["summary"].to_csv(out_dir / f"fdr_fixed{sfx}.csv", index=False)
    summ = {"fdr_emp_fixed": dict(zip(nul["summary"]["list"], nul["summary"]["fdr_emp_fixed"])),
            "pass": dict(zip(nul["summary"]["list"], nul["summary"]["pass"].astype(bool))),
            "fdr_emp_rerun": None, "p_cut_F1": {k: v for k, v in res["p_cut"].items() if k.endswith("|F1")},
            "p_cut_sc": res["p_cut_sc"], "lambda_source": "observed", "B": nul["B"],
            "wall_draw0_s": nul["wall_draw0_s"], "hits_obs": nul["hits_obs"], "n_genes": len(genes),
            "config_sha256": PR.config_hashes(), "commit": gitutil.head(P.PROJECT)}
    (out_dir / f"summary{sfx}.json").write_text(json.dumps(summ, indent=1, default=float), encoding="utf-8")
    from biomark import analyses as AN
    td = P.PROJECT / PR.pp()["timing_dir"]
    f_t = td / "null_fixed.json"
    parked = out_dir / "timing_s1_vectorised_subrecord.json"     # not budgeted until the main null record exists
    if a.s1_vectorised:
        srec = AN.timing_subrecord("s1_vectorised", t_start, AN.utc_now(), w, B=nul["B"],
                                   wall_draw0_s=nul["wall_draw0_s"], memory_guard=gnull)
        if f_t.exists():
            prev = json.loads(f_t.read_text(encoding="utf-8"))
            prev.setdefault("subrecords", {})["s1_vectorised"] = srec
            f_t.write_text(json.dumps(prev, indent=1, default=float), encoding="utf-8")
        else:
            parked.write_text(json.dumps(srec, indent=1, default=float), encoding="utf-8")
    else:
        subs = {}
        if f_t.exists():
            subs.update(json.loads(f_t.read_text(encoding="utf-8")).get("subrecords", {}))
        if parked.exists():
            subs["s1_vectorised"] = json.loads(parked.read_text(encoding="utf-8"))
        extra = {"subrecords": subs} if subs else {}
        AN.write_timing("null_fixed", {"B": nul["B"], "wall_draw0_s": nul["wall_draw0_s"], "memory_guard": gnull, **extra},
                        start_utc=t_start,
                        workers=w)
        if parked.exists():
            parked.unlink()
    AN.append_registry("scripts/run_null.py" + (" --s1-vectorised" if a.s1_vectorised else ""), "ok",
                       [f"results/null/summary{sfx}.json", f"results/null/fdr_fixed{sfx}.csv",
                        f"results/null/draws{sfx}.csv"])
    print(json.dumps(summ, indent=1, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
