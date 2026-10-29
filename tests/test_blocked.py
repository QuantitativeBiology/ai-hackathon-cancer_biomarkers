"""Streamed (blocked) arms orchestration == the in-memory run_arms on the toy: every arms-table column (BH sets, ranks,
de-duplicated ranks), MD-strict, MD-global, EN-obs-strict, the observed p-cutoffs and the lambda table; serial and with
2 workers. Also: the reader API, the memory guard and the FWL-cache eviction policy leave results unchanged."""
import numpy as np
import pandas as pd
import pytest

from biomark import analyses as AN
from biomark import arms as A
from biomark import blocked as BL
from biomark import pipeline as PL
from biomark import toysim as TS

KEY = ["arm", "gene", "feature"]
# discrete outputs compared exactly (BH sets, selections, ranks, de-duplicated ranks, counts)
EXACT = {"selected", "selected_F2", "rank", "rank_dedup", "n_obs", "n_carriers", "n_NA", "n_carrier_clusters",
         "n_lev1_excluded", "pass", "same_sign", "s1_f1_hit", "s1_obs_selected", "in_md_strict", "K", "n_rows",
         "selection_freq", "maj_sign", "sign_consistency", "en_max_attainable_freq", "d_stage_testable"}


@pytest.fixture(scope="module")
def toy():
    return PL.Data.from_toy(TS.make_toy())


@pytest.fixture(scope="module")
def ref(toy):
    return A.run_arms(toy, arms=A.ARMS, workers=1)


def _norm(df, key=KEY):
    df = df.copy()
    for c in df.columns:
        if isinstance(df[c].dtype, pd.CategoricalDtype):
            df[c] = df[c].astype(object)
    return df.sort_values(key, kind="mergesort").reset_index(drop=True)[sorted(df.columns)]


def _cmp(a, b, key=KEY):
    a, b = _norm(a, key), _norm(b, key)
    assert list(a.columns) == list(b.columns), set(a.columns) ^ set(b.columns)
    assert len(a) == len(b)
    for c in a.columns:
        x, y = a[c], b[c]
        if c in EXACT or pd.api.types.is_bool_dtype(x) or not pd.api.types.is_numeric_dtype(x):
            pass
        elif pd.api.types.is_numeric_dtype(y):
            # continuous statistics: last-bit BLAS differences between chunkings (column blocks) are allowed
            np.testing.assert_allclose(x.to_numpy(dtype=float), y.to_numpy(dtype=float), rtol=1e-9, atol=1e-12,
                                       err_msg=c)
            continue
        if pd.api.types.is_numeric_dtype(x) and not pd.api.types.is_bool_dtype(x):
            np.testing.assert_array_equal(x.to_numpy(dtype=float), y.to_numpy(dtype=float), err_msg=c)
        else:
            assert (x.astype(object).where(x.notna(), None).tolist() == y.astype(object).where(y.notna(), None).tolist()), c


@pytest.mark.parametrize("workers", [1, 2])
def test_blocked_equals_unblocked(toy, ref, tmp_path, workers):
    out = BL.run_streamed(toy, A.ARMS, out_dir=tmp_path / "arms", stage_dir=tmp_path / "stage", workers_s1=workers,
                          workers_en=workers, chunk=40, log=lambda *a: None)
    got = BL.read_dataset(tmp_path / "arms")
    assert out["n_rows"] == len(ref["arms"]) and len(out["parts"]) > 3
    _cmp(got, ref["arms"])
    for k in ("md_strict", "md_global", "en_obs_strict", "md_strict_all", "md_global_all"):
        _cmp(out[k], ref[k], ["gene", "feature"])
    _cmp(out["lam"], ref["lam"], ["arm", "gene"])
    assert out["p_cut"] == ref["p_cut"] and out["p_cut_sc"] == ref["p_cut_sc"]
    assert not (tmp_path / "stage").exists()
    # string columns are stored as categoricals
    raw = pd.read_parquet(sorted((tmp_path / "arms").glob("part-*.parquet"))[0])
    assert isinstance(raw["gene"].dtype, pd.CategoricalDtype) and isinstance(raw["feature"].dtype, pd.CategoricalDtype)


def test_reader_api_filters_and_columns(toy, ref, tmp_path, monkeypatch):
    from biomark import endpoint_io as IO
    BL.run_streamed(toy, A.ARMS, out_dir=tmp_path / "arms", stage_dir=tmp_path / "stage", chunk=60, log=lambda *a: None)
    monkeypatch.setattr(IO, "ARMS_DIR", tmp_path / "arms")
    sel = IO.read_arms(columns=["arm", "gene", "feature", "beta"], arms=["S1-obs"], selected_only=True)
    r = ref["arms"]
    exp = r[(r["arm"] == "S1-obs") & r["selected"]]
    assert list(sel.columns) == ["arm", "gene", "feature", "beta"] and len(sel) == len(exp)
    assert set(zip(sel["gene"], sel["feature"])) == set(zip(exp["gene"], exp["feature"]))
    g = toy.genes[:3]
    sub = IO.read_arms(columns=["arm", "gene", "rank"], genes=g)
    assert set(sub["gene"]) <= set(g) and len(sub) == int(r["gene"].isin(g).sum())
    assert sub["gene"].dtype == object                                  # plain strings for downstream code


def test_eviction_policy_does_not_change_results(toy):
    U = A.universe(toy, toy.genes[:5])
    a, _, _ = A.run_en(toy, ("EN-pred",), genes=toy.genes[:5], U=U)
    chunks = [[g] for g in toy.genes[:5]]
    parts = []
    A.run_en(toy, ("EN-pred",), U=U, chunks=chunks, on_chunk=lambda i, arm, t, l, tm: parts.append(t), finish=False,
             max_patterns=1)
    b = A.finish_en(pd.concat(parts, ignore_index=True))
    _cmp(a, b)


def test_memory_guard_caps_workers(monkeypatch):
    import psutil

    class VM:
        available = 3 * (1 << 30)

    monkeypatch.setattr(psutil, "virtual_memory", lambda: VM())
    eff, info = AN.guard_workers(8, "en")
    assert 1 <= eff < 8 and info["effective"] == eff and info["requested"] == 8
    VM.available = 64 * (1 << 30)
    eff2, _ = AN.guard_workers(6, "en")
    assert eff2 <= 6 and eff2 >= eff
    assert AN.guard_workers(1, "s1")[0] == 1
