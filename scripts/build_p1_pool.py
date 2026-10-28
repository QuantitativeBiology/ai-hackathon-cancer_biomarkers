"""Freeze the P1 pool, lists and the R request list BEFORE tag freeze-catalogue (DESIGN v3.2 §1.3, §6.1).

    python -B scripts/build_p1_pool.py

Reads the observed arms (cache/benchmark/arms.parquet), MD-strict / EN-obs-strict, the positive controls and the CGT
table; R lines contribute MOLECULAR data only (features, covariates, library strings for the expected-NaN rule); no
R effect is read. D effects (dev scope) give σ̂_resid,D(g). Writes results/benchmark/p1_pool.json,
cache/benchmark/p1_pool.parquet (bulk), configs/catalogue/r_requests.csv and configs/freeze_catalogue.json (sha256 of
the catalogue files that replicate.py verifies). Refuses unless HEAD descends from freeze-bounds and before
freeze-catalogue exists.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd  # noqa: E402

from biomark import endpoint_io as IO  # noqa: E402
from biomark import endpoints as EP  # noqa: E402
from biomark import genesets as GS  # noqa: E402
from biomark import gitutil  # noqa: E402
from biomark import guard  # noqa: E402
from biomark import p1  # noqa: E402
from biomark import paths as P  # noqa: E402
from biomark import stats as ST  # noqa: E402

CATALOGUE_FILES = ["configs/catalogue/r_requests.csv", "results/benchmark/p1_pool.json", "results/md/md_global.csv",
                   "results/md/md_strict.csv",
                   "results/md/en_obs_strict.csv", "results/benchmark/p_cut.json", "results/benchmark/lambda.csv"]


def main() -> int:
    EP.require_descends("freeze-bounds")
    if gitutil.tag_exists(P.PROJECT, "freeze-catalogue"):
        raise EP.EndpointRefusal("freeze-catalogue exists: the pool and r_requests are frozen")
    a = IO.arms_table()
    s1o = a[a["arm"] == "S1-obs"].copy()
    dev, test = list(p1.dev_ids()), list(p1.test_ids())
    genes = sorted(s1o["gene"].unique())
    nan_D = pd.read_parquet(P.CACHE / "effect_nan_mask.parquet").reindex(index=dev, columns=genes)
    Y = guard.load_effects(dev, scope="dev", genes=genes).to_numpy()
    covD = ST.CovDesign(IO.line_covariates(dev))
    sig = IO.sigma_resid_genes(Y, ~nan_D.to_numpy(dtype=bool), covD, genes)
    pool, lists = IO.make_pool(a, s1o, sig, dev, test, nan_D)
    (P.CACHE / "benchmark").mkdir(parents=True, exist_ok=True)
    pool.to_parquet(P.CACHE / "benchmark" / "p1_pool.parquet")
    lists["pool_file_sha256"] = __import__("hashlib").sha256((P.CACHE / "benchmark" / "p1_pool.parquet").read_bytes()).hexdigest()
    IO.write_json(P.RESULTS / "benchmark" / "p1_pool.json", lists)
    ctl = pd.read_csv(P.CONFIGS / "positive_controls.csv", dtype=str)
    ctl = ctl[ctl["scope"] == "gate"]
    cgt = pd.read_csv(P.CONFIGS / "clinical_ground_truth.csv")
    ez = GS.symbol_to_entrez(sorted({x for v in cgt["target_set"] for x in EP._targets(v)}))
    tests = EP.p2_test_set(cgt, s1o, "incl_proxy", entrez_of=ez.dropna().astype(str).to_dict())
    s1h = a[(a["arm"] == "S1-obs") & a["selected"].astype(bool)]
    req = IO.build_requests(pool, lists, s1o, IO.md_strict(), IO.en_obs_strict(), s1h, ctl, tests, mg=IO.md_global())
    out = P.CONFIGS / "catalogue" / "r_requests.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    req.to_csv(out, index=False)
    present = [f for f in CATALOGUE_FILES if (P.PROJECT / f).exists()]
    IO.write_json(P.CONFIGS / "freeze_catalogue.json", EP.catalogue_manifest(present))
    print(f"pool {lists['pool_size']} (k_G = {lists['k']}); r_requests {len(req)} rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
