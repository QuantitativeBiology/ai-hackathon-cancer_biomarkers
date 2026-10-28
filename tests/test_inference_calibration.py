"""DESIGN §13.2 test_inference_calibration: size of the §5.3 CR2 + ν test under the null on the real covariate designs.

Per n ∈ {255 (R lines), 1022 (D lines)} × carriers ∈ {3, 5, 10} × {homoskedastic, carrier SD × 2}: 20,000 null
simulations (1,000 random carrier placements × 20 null outcome draws per placement). Each simulation is the
binary-d test of y ~ [shared design] + d with PatientID clusters, computed by the engine's own core
(``stats.orth`` / ``stats.resid`` / ``stats.cr2_core`` with ``stats.Clusters``, the path used by ``pair_test`` and
``s1_scan``). Pass: size at α 0.05 in [0.035, 0.065] and at α 1e-3 ≤ 2e-3.

Real inputs: covariates only (cache/covariates.parquet, built by scripts/build_screen_covariates.py from screen
metadata, no effects); the outcome is simulated noise. The shared design is the S1 model's design blocks plus the
droppable shared covariates (no drops; gene covariates are gene-specific and omitted).

Set BIOMARK_CALIB_SIMS to change the per-cell count for exploration only; the pass rule is evaluated only at the
design count (other counts skip).
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd
import pytest

from biomark import paths as P
from biomark import stats as ST

N_SIMS = 20_000                       # DESIGN §13.2
N_PLACEMENTS, N_Y = 1_000, 20         # 1,000 × 20 = 20,000
SIZE_005 = (0.035, 0.065)
SIZE_0001_MAX = 2e-3
CARRIERS = (3, 5, 10)
ROLES = {255: "test", 1022: "dev"}
HETERO = {"homoskedastic": 1.0, "carrier_sd_x2": 2.0}
RESULTS = {}


def _sims() -> int:
    return int(os.environ.get("BIOMARK_CALIB_SIMS", N_SIMS))


def _covariates(role: str) -> pd.DataFrame:
    """Canonical S1 covariate columns (as pipeline.real_covariates) for the lines of one role."""
    c = pd.read_parquet(P.CACHE / "covariates.parquet")
    c = c[c["role"] == role]
    out = pd.DataFrame(index=c.index)
    for k in ("lineage", "library_group", "cellformat", "serumfree", "medium_family", "doubling_time",
              "DT_unknown", "MSI", "DMG_BURDEN", "PatientID"):
        out[k] = c[k]
    out["tp53_status"] = c["TP53"]
    dropped = set(c.attrs.get("dropped", []))
    out["screen_days"] = 0.0 if "screen_days_dev_humagne" in dropped else c["screen_days_dev_humagne"]
    return out


@pytest.fixture(scope="module")
def designs():
    out = {}
    for n, role in ROLES.items():
        cov = _covariates(role)
        assert len(cov) == n, (role, len(cov))
        cd = ST.CovDesign(cov)
        rows = np.ones(n, dtype=bool)
        Z = cd.Z(rows, cd.design_names + cd.droppable_names)
        Q = ST.orth(Z)
        out[n] = {"Z": Z, "Q": Q, "hZ": (Q ** 2).sum(1), "cl": ST.Clusters.of(cov["PatientID"].astype(str).to_numpy()),
                  "patients": cov["PatientID"].astype(str).to_numpy()}
    return out


def _simulate(des: dict, k: int, sd_carrier: float, n_sims: int, seed: int) -> np.ndarray:
    Q, hZ, cl = des["Q"], des["hZ"], des["cl"]
    n = Q.shape[0]
    rng = np.random.default_rng(seed)
    n_pl = max(1, n_sims // N_Y)
    out = np.empty((n_pl, N_Y))
    for i in range(n_pl):
        d = np.zeros(n)
        d[rng.choice(n, size=k, replace=False)] = 1.0
        sd = np.where(d == 1, sd_carrier, 1.0)
        Y = rng.standard_normal((n, N_Y)) * sd[:, None]
        Xt = ST.resid(Q, d[:, None])
        Yt = ST.resid(Q, Y)
        r = ST.cr2_core(Yt, Xt, hZ, Q, cl, Q.shape[1], car_clusters=ST.carrier_clusters(d, True, cl))
        out[i] = r["p"][0]
    return out.ravel()


def test_vectorised_core_equals_pair_test(designs):
    """The simulation path (orth/resid/cr2_core) equals stats.pair_test on the same inputs."""
    des = designs[255]
    rng = np.random.default_rng(11)
    n = des["Z"].shape[0]
    for k in CARRIERS:
        d = np.zeros(n)
        d[rng.choice(n, size=k, replace=False)] = 1.0
        y = rng.standard_normal(n) * np.where(d == 1, 2.0, 1.0)
        r = ST.pair_test(y, d, des["Z"], des["patients"])
        Xt, Yt = ST.resid(des["Q"], d[:, None]), ST.resid(des["Q"], y[:, None])
        c = ST.cr2_core(Yt, Xt, des["hZ"], des["Q"], des["cl"], des["Q"].shape[1],
                        car_clusters=ST.carrier_clusters(d, True, des["cl"]))
        assert c["p"][0, 0] == pytest.approx(r["p"], rel=1e-10)


@pytest.mark.parametrize("hetero", list(HETERO))
@pytest.mark.parametrize("k", CARRIERS)
@pytest.mark.parametrize("n", list(ROLES))
def test_null_size(designs, n, k, hetero):
    n_sims = _sims()
    seed = 1000 * n + 10 * k + list(HETERO).index(hetero)
    p = _simulate(designs[n], k, HETERO[hetero], n_sims, seed)
    ok = np.isfinite(p)
    s05, s001 = float(np.mean(p[ok] < 0.05)), float(np.mean(p[ok] < 1e-3))
    RESULTS[(n, k, hetero)] = {"n_sims": int(p.size), "n_finite": int(ok.sum()), "size_0.05": s05, "size_0.001": s001}
    _dump()
    assert ok.all(), f"{(~ok).sum()} non-finite p"
    if n_sims != N_SIMS:
        pytest.skip(f"exploratory count {n_sims} != design {N_SIMS}: size_0.05={s05:.4f}, size_0.001={s001:.5f}")
    assert SIZE_005[0] <= s05 <= SIZE_005[1], f"size at 0.05 = {s05:.4f} (n={n}, carriers={k}, {hetero})"
    assert s001 <= SIZE_0001_MAX, f"size at 1e-3 = {s001:.5f} (n={n}, carriers={k}, {hetero})"


def _dump():
    path = os.environ.get("BIOMARK_CALIB_OUT")
    if path:
        rows = [{"n": n, "carriers": k, "variance": h, **v} for (n, k, h), v in sorted(RESULTS.items())]
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(rows, fh, indent=1)
