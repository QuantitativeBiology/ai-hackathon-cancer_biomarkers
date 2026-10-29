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


@pytest.fixture(scope="module")
def power_inputs():
    """pool_power inputs on the toy: 40 eligible genes, every second line as the 'R' molecular set, RDesign with own CN
    and EXPGATE, expected-NaN rows from the toy effects' NaN pattern, the S1-obs drops."""
    from biomark import endpoints as EP
    from biomark import params as PR
    from biomark import pipeline as PL
    from biomark import toysim as TS
    toy = TS.make_toy()
    el = PL.toy_eligibility(toy)
    genes = [g for g in toy.genes.index if el.loc[g, "eligible"]][:40]
    d = PL.Data.from_toy(toy, genes=genes)
    s1 = A.finish_s1(d, A.run_s1(d, ("S1-obs",)))
    fm = toy.feature_meta
    UG = s1[(fm["stage"].reindex(s1["feature"]).to_numpy() == "G") & (s1["test_family"] == "F1")][
        ["gene", "feature"]].drop_duplicates()
    R = list(d.lines)[::2]
    feats = list(dict.fromkeys(UG["feature"]))
    thr = PR.dp()["dictionary_s3"]["expressed_log2tpm"]
    covR = EP.RDesign(d.cov_df.loc[R], toy.cn_log2.reindex(index=R, columns=genes),
                      (toy.expression.reindex(index=R, columns=genes) >= thr).astype(float))
    key = pd.MultiIndex.from_arrays([UG["gene"], UG["feature"]])
    rng = np.random.default_rng(0)
    return dict(U_rows=UG, beta_shrunk=pd.Series(rng.normal(0, 0.3, len(UG)), index=key),
                sigma_g=pd.Series(0.3, index=genes), XR=toy.features.reindex(index=R, columns=feats).to_numpy(float),
                fpos=pd.Series(np.arange(len(feats)), index=feats), binary=fm["d_type"].reindex(feats) == "binary",
                covR=covR, patientsR=d.cov_df.loc[R, "PatientID"].astype(str).to_numpy(),
                expR=toy.effects.reindex(index=R, columns=genes).notna(),
                drops=s1.set_index(["gene", "feature"])["dropped_cov"])


@pytest.mark.parametrize("workers", [1, 2])
def test_pool_power_fast_equals_slow(power_inputs, workers):
    from biomark import endpoints as EP
    slow = EP.pool_power(**power_inputs, fast=False)
    fast = EP.pool_power(**power_inputs, workers=workers)
    pd.testing.assert_frame_equal(slow, fast, check_exact=True)
    assert slow["powered"].any() and (slow["df_rule"] != "").any() and slow["n_R"].nunique() > 1
