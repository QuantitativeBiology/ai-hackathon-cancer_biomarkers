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
    res = {"arms": pd.read_parquet(P.CACHE / "benchmark" / "arms.parquet"),
           "lam": pd.read_csv(rb / "lambda.csv", dtype={"gene": str}),
           "md_strict": pd.read_csv(P.RESULTS / "md" / "md_strict.csv", dtype={"gene": str}),
           "p_cut": {k: v for k, v in pc.items() if k != "MD-strict"}, "p_cut_sc": pc.get("MD-strict", 0.0)}
    w = a.workers or npar["workers"]
    out_dir = P.RESULTS / "null"
    out_dir.mkdir(parents=True, exist_ok=True)
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
    td = P.PROJECT / PR.pp()["timing_dir"]
    td.mkdir(parents=True, exist_ok=True)
    (td / f"null_fixed{sfx}.json").write_text(json.dumps({
        "job": f"null_fixed{sfx}", "utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "wall_s": wall, "B": nul["B"], "workers": w, "wall_draw0_s": nul["wall_draw0_s"]}, indent=1), encoding="utf-8")
    print(json.dumps(summ, indent=1, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
