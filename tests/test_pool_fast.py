"""The fast make_pool path (endpoint_io.shrunk_beta_lowmem + arms rows restricted to P_G before the lists) gives output
identical to the slow path (stats.shrunk_beta, unfiltered lists): bit-identical EB on synthetic data across chunk
sizes, and make_pool on the toy (four-arm in-memory table and the two-arm streamed XLIB dataset) with identical pool
and lists. Toy / synthetic data only."""
import numpy as np
import pandas as pd
import pytest

from biomark import arms as A
from biomark import endpoint_io as IO
from biomark import stats as ST
from test_xlib_pool import _pool_args, toy_io, world  # noqa: F401  (module fixtures reused)


@pytest.mark.parametrize("chunk", [1, 257, 4096, 10 ** 6])
def test_lowmem_eb_bit_identical(chunk):
    rng = np.random.default_rng(3)
    n = 6000
    s = rng.uniform(0.05, 1.0, n)
    b = rng.normal(0, 1, n) * (rng.random(n) < 0.2) + rng.normal(0, s)
    sd = rng.uniform(0.2, 2.0, n)
    b[::97] = np.nan
    sd[::131] = 0.0
    ref = ST.shrunk_beta(b, s, sd)
    got = IO.shrunk_beta_lowmem(b, s, sd, chunk=chunk)
    assert np.array_equal(ref, got, equal_nan=True)


def _same(a, b):
    pa, la = a
    pb, lb = b
    pd.testing.assert_frame_equal(pa, pb, check_exact=True)
    assert la == lb


def test_make_pool_fast_equals_slow_four_arms(toy_io):
    w = toy_io
    full = A.run_arms(w["d"], arms=A.ARMS, workers=1)["arms"]
    tab = full[IO.RANK_COLS].copy()
    s1o = full.loc[full["arm"] == "S1-obs", [c for c in IO.S1O_COLS if c in full]].reset_index(drop=True)
    args = _pool_args(w)
    slow = IO.make_pool(tab, s1o, *args, fast=False)
    fast = IO.make_pool(tab, s1o, *args)
    _same(slow, fast)
    assert all(fast[1]["lists"]["G"][a]["n_list"] > 0 for a in A.ARMS)


def test_make_pool_fast_equals_slow_xlib_dataset(toy_io):
    w = toy_io
    tab, s1o = IO.pool_inputs(dataset_dir=w["root"] / "xlib_arms")
    xarms = ("S1-obs", "EN-obs")
    args = _pool_args(w)
    _same(IO.make_pool(tab, s1o, *args, arms=xarms, fast=False), IO.make_pool(tab, s1o, *args, arms=xarms))
