"""File-backed worker payloads (src/biomark/workerstore.py): the parallel path equals the serial path on toy data, pred
pools carry no observed-effect file, and the worker cache is removed after each pool."""
import numpy as np
import pandas as pd
import pytest

from biomark import arms as A
from biomark import guard
from biomark import pipeline as PL
from biomark import toysim as TS
from biomark import workerstore as WS

from synth_data import make_data


@pytest.fixture(scope="module")
def toy_data():
    return PL.Data.from_toy(TS.make_toy())


def _sorted(df, key):
    return df.sort_values(key, kind="mergesort").reset_index(drop=True)


def test_file_backed_parallel_s1_equals_serial(toy_data):
    d = toy_data
    genes = d.genes[:8]
    ser = A.run_s1(d, A.S1_ARMS, genes=genes)
    par = A.run_s1(d, A.S1_ARMS, genes=genes, workers=2, chunk=3)
    key = ["arm", "gene", "feature"]
    pd.testing.assert_frame_equal(_sorted(ser, key), _sorted(par, key))


def test_file_backed_parallel_en_equals_serial(toy_data):
    d = toy_data
    genes = d.genes[:6]
    U = A.universe(d, genes)
    ser, lam_s, _ = A.run_en(d, ("EN-pred", "EN-obs"), genes=genes, U=U, workers=1)
    par, lam_p, _ = A.run_en(d, ("EN-pred", "EN-obs"), genes=genes, U=U, workers=2, chunk=2)
    key = ["arm", "gene", "feature"]
    a, b = _sorted(ser, key), _sorted(par, key)
    pd.testing.assert_frame_equal(a.drop(columns=["rank"]), b.drop(columns=["rank"]))
    assert np.allclose(a["rank"].fillna(-1), b["rank"].fillna(-1))
    pd.testing.assert_frame_equal(_sorted(lam_s, ["arm", "gene"]), _sorted(lam_p, ["arm", "gene"]))


def test_large_arrays_are_memory_mapped_and_cleaned_up(toy_data):
    d = toy_data
    path = WS.stage({"data": d.for_workers(include_obs=True), "x": 1}, "t", min_bytes=1 << 16)
    fl = WS.files(path)
    assert "payload.pkl" in fl and any(f.endswith("_X.npy") for f in fl) and any(f.endswith("_obs.npy") for f in fl)
    back = WS.load(path)
    assert isinstance(back["data"].X, np.memmap) and not back["data"].X.flags.writeable
    np.testing.assert_array_equal(np.asarray(back["data"].X), d.X)
    np.testing.assert_array_equal(np.asarray(back["data"].target("obs")), d.target("obs"))
    assert back["x"] == 1 and back["data"].genes == d.genes
    del back
    WS.cleanup(path)
    import os
    assert not os.path.exists(path)


def test_pred_pool_store_has_no_observed_effect_file(toy_data):
    d = toy_data
    d.target("obs")                                          # observed effects loaded in the parent ...
    path = WS.stage({"data": d.for_workers(include_obs=False)}, "pred", min_bytes=1 << 16)   # small arrays to npy too
    try:
        fl = WS.files(path)
        back = WS.load(path)
        assert not any("obs" in f for f in fl), fl          # ... but never written for a pred pool
        assert back["data"]._obs is None and back["data"].obs_loader is None
        assert any(f.endswith("_pred.npy") for f in fl)     # the prediction target does travel
        del back
    finally:
        WS.cleanup(path)


def test_parallel_en_pred_is_blind_to_observed_effects(monkeypatch):
    def refuse(*a, **k):
        raise AssertionError("EN-pred read observed effects")

    monkeypatch.setattr(guard, "load_effects", refuse)
    data = make_data(seed=4, obs_loader=lambda: guard.load_effects(scope="dev"), G=4)
    U = A.universe(data)
    tab, lam, _ = A.run_en(data, ("EN-pred",), U=U, workers=2, chunk=1)
    assert len(tab) and data._obs is None
    ser, _, _ = A.run_en(data, ("EN-pred",), U=U, workers=1)
    key = ["gene", "feature"]
    assert np.array_equal(_sorted(tab, key)["selection_freq"].to_numpy(), _sorted(ser, key)["selection_freq"].to_numpy())
