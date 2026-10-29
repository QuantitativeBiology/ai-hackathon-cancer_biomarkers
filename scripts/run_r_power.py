"""R power for every replication test (DESIGN v3.2 §6.1 power rule; §12.2 classes 5 / 7) -> results/replication/r_power.csv.

    python -B scripts/run_r_power.py [--workers N]

The frozen pool (build_p1_pool.py) computed R power only for its own universe U_G (stage-G F1 pairs), so the
replication tests of other pairs (stage-X families such as PROG, EXPPAR, PSCORE, EXPDOSE, SIL-other; F2) carry
power = NaN in results/replication/r_tests.csv. This script applies the same pre-declared formula
(design_params power_formula: 1 - F_t(t_0.95, nu_R; ncp = |beta_shrunk| / SE_R_pred), with the binary or continuous
SE_R_pred of stats.se_r_pred) through the same code (endpoint_io.pair_power -> endpoints.pool_power) to every
(gene, feature) of r_tests.csv. Inputs need no R effect: R molecular data (feature values, covariates, own CN, EXPGATE,
library strings), the D-side EB shrinkage over every S1-obs row and σ̂_resid,D(g) from D effects (guard scope dev).
Pairs whose frozen pool power exists are recomputed as a check and must agree (relative 1e-9; last-bit floating
differences only); their frozen value is kept
(``power_source`` = frozen pool). ``power_effective`` = power if the pair is R-evaluable, else 0 (the pool's rule:
powered = R-evaluable and power >= power_min). Registry row; no timing stem (not a budgeted stage).
Refuses unless HEAD descends from freeze-catalogue and r_tests.csv exists.
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
from biomark import endpoint_io as IO  # noqa: E402
from biomark import endpoints as EP  # noqa: E402
from biomark import paths as P  # noqa: E402
from biomark import stats as ST  # noqa: E402

OUT = P.RESULTS / "replication" / "r_power.csv"
RTOL = 1e-9          # recomputation vs frozen pool: last-bit floating differences only (BLAS / chunking)
COLS = ["gene", "feature", "power", "R_evaluable", "power_effective", "power_source", "n_R", "n1_R", "n0_R", "sd_R",
        "nu_R", "se_R_pred", "beta_shrunk", "binary"]


def merge_power(rt: pd.DataFrame, pw: pd.DataFrame) -> pd.DataFrame:
    """Final table: frozen pool power where r_tests has it (asserted equal to the recomputation), else recomputed."""
    rt = rt.assign(gene=rt["gene"].astype(str))
    m = rt[["gene", "feature", "power"]].rename(columns={"power": "power_frozen"}).merge(
        pw.assign(gene=pw["gene"].astype(str)), on=["gene", "feature"], how="left")
    fz = m["power_frozen"].notna()
    bad = fz & ~np.isclose(m["power_frozen"], m["power"], rtol=RTOL, atol=0, equal_nan=False)
    if bad.any():
        raise AssertionError(f"{int(bad.sum())} recomputed powers differ from the frozen pool, e.g. "
                             f"{m.loc[bad, ['gene', 'feature', 'power_frozen', 'power']].head(3).to_dict('records')}")
    m["power_source"] = np.where(fz, "frozen pool (recomputation identical)",
                                 "computed after the R read (R molecular data, frozen D side; no R effects)")
    m["power"] = np.where(fz, m["power_frozen"], m["power"])
    ev = m["R_evaluable"].astype("boolean").fillna(False).astype(bool)
    m["power_effective"] = np.where(ev, m["power"], 0.0)
    return m[COLS]


def run(workers: int = 1, check=None) -> pd.DataFrame:
    EP.require_descends("freeze-catalogue", check=check)
    from biomark import guard, p1
    t0 = time.perf_counter()
    rt = pd.read_csv(P.RESULTS / "replication" / "r_tests.csv", dtype={"gene": str})
    s1o = IO.read_arms(arms=["S1-obs"], columns=IO.S1O_COLS)
    dev, test = list(p1.dev_ids()), list(p1.test_ids())
    genes = sorted(rt["gene"].astype(str).unique())
    nan_D = pd.read_parquet(P.CACHE / "effect_nan_mask.parquet").reindex(index=dev, columns=genes)
    Y = guard.load_effects(dev, scope="dev", genes=genes).reindex(columns=genes).to_numpy()
    sig = IO.sigma_resid_genes(Y, ~nan_D.to_numpy(dtype=bool), ST.CovDesign(IO.line_covariates(dev)), genes)
    pw = IO.pair_power(rt[["gene", "feature"]], s1o, sig, dev, test, nan_D, workers=workers)
    (P.CACHE / "replication").mkdir(parents=True, exist_ok=True)
    pw.to_parquet(P.CACHE / "replication" / "r_power_raw.parquet")      # before the identity check (re-check without a rerun)
    out = merge_power(rt, pw)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    AN.append_registry("scripts/run_r_power.py", "ok", ["results/replication/r_power.csv"])
    print(f"r_power: {len(out)} tests; {int(out['power'].notna().sum())} with power "
          f"({int((out['power_source'].str.startswith('computed')).sum())} computed after R); "
          f"wall {time.perf_counter() - t0:.0f} s")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--workers", type=int, default=1)
    a = ap.parse_args(argv)
    run(workers=a.workers)
    return 0


if __name__ == "__main__":
    sys.exit(main())
