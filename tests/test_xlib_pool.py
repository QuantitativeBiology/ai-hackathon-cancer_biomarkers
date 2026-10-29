"""XLIB (DESIGN v3.2 §6.3) pool and lists on a two-arm streamed dataset: the toy pipeline through
blocked.run_streamed(S1-obs, EN-obs) -> endpoint_io.pool_inputs -> endpoint_io.make_pool. Before the fix make_pool
listed all four benchmark arms and build_lists raised KeyError 'rank' for an arm absent from the dataset (an empty
row filter selected no columns); the lists of the XLIB arms must equal those built from the in-memory run. Also the
completed-dataset manifest used by run_xlib.py --reuse-arms. Toy data only (R molecular data = a toy line subset)."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from biomark import arms as A
from biomark import blocked as BL
from biomark import endpoint_io as IO
from biomark import endpoints as EP
from biomark import params as PR
from biomark import pipeline as PL
from biomark import stats as ST
from biomark import toysim as TS

ROOT = Path(__file__).resolve().parents[1]
XARMS = tuple(PR.dp()["replication_s6"]["xlib_arms"])
N_GENES = 12


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    toy = TS.make_toy()
    el = PL.toy_eligibility(toy)
    genes = [g for g in toy.genes.index if el.loc[g, "eligible"]][:N_GENES]
    d = PL.Data.from_toy(toy, genes=genes)
    root = tmp_path_factory.mktemp("xlib_pool")
    res = BL.run_streamed(d, XARMS, out_dir=root / "xlib_arms", stage_dir=root / "stage", workers_s1=1, workers_en=1,
                          log=lambda *a: None)
    ref = A.run_arms(d, arms=XARMS, workers=1)
    return {"toy": toy, "d": d, "root": root, "res": res, "ref": ref, "genes": list(d.genes)}


@pytest.fixture
def toy_io(world, monkeypatch):
    """endpoint_io's real-data readers (dictionary caches, covariates, own CN / EXPGATE, library strings) on the toy."""
    toy, d = world["toy"], world["d"]
    thr = PR.dp()["dictionary_s3"]["expressed_log2tpm"]
    cov = d.cov_df.copy()
    monkeypatch.setattr(IO, "feature_meta", lambda: toy.feature_meta[["family", "stage", "d_type"]])
    monkeypatch.setattr(IO, "line_values", lambda ids, feats: toy.features.reindex(index=list(ids), columns=list(feats)))
    monkeypatch.setattr(IO, "line_covariates", lambda ids: cov.loc[list(ids)])
    monkeypatch.setattr(IO, "line_gene_covs", lambda ids, genes: (
        toy.cn_log2.reindex(index=list(ids), columns=list(genes)),
        (toy.expression.reindex(index=list(ids), columns=list(genes)) >= thr).astype(float)))
    monkeypatch.setattr(IO, "library_strings", lambda ids: cov["library_group"].astype(str).reindex(list(ids)))
    return world


def _pool_args(world):
    d, genes = world["d"], world["genes"]
    Y = d.target("obs")
    sig = IO.sigma_resid_genes(Y, d.rows_ok, ST.CovDesign(d.cov_df), genes)
    D_ids = list(d.lines)
    R_ids = D_ids[::2]                                   # toy stand-in for the replication lines' molecular data
    nan_D = pd.DataFrame(~d.rows_ok, index=d.lines, columns=genes)
    return sig, D_ids, R_ids, nan_D


def test_make_pool_on_two_arm_streamed_dataset(toy_io):
    w = toy_io
    arms_tab, s1o = IO.pool_inputs(dataset_dir=w["root"] / "xlib_arms")
    assert set(arms_tab["arm"]) == set(XARMS) and "rank" in arms_tab
    sig, D_ids, R_ids, nan_D = _pool_args(w)
    pool, pj = IO.make_pool(arms_tab, s1o, sig, D_ids, R_ids, nan_D)              # the call that crashed
    for arm in set(A.ARMS) - set(XARMS):                                         # absent arms: empty, never an error
        assert pj["lists"]["G"][arm]["n_list"] == 0
    pool2, pj2 = IO.make_pool(arms_tab, s1o, sig, D_ids, R_ids, nan_D, arms=XARMS)   # run_xlib: XLIB arms only
    assert set(pj2["lists"]["G"]) == set(XARMS)
    pd.testing.assert_frame_equal(pool, pool2)
    P_G = pool2[pool2["powered"]]
    assert len(P_G) and pj2["k"] == EP.p1_k(len(P_G))
    # §6.3 lists follow the P1 rules: identical to the lists of the in-memory two-arm run on the same pool
    ref = w["ref"]["arms"]
    want = EP.build_lists(ref, P_G, XARMS)
    for arm in XARMS:
        assert pj2["lists"]["G"][arm] == want[arm], arm
        assert pj2["lists"]["G"][arm]["n_list"] > 0, arm
        assert pj2["lists"]["G"][arm] == pj["lists"]["G"][arm]


def test_build_lists_absent_arm_is_empty():
    tab = pd.DataFrame({"arm": "S1-obs", "gene": ["g1"], "feature": ["HOT:1"], "rank": [1.0], "rank_stat": [0.1],
                        "variant_group": [None], "beta": [-1.0]})
    pool = pd.DataFrame({"gene": ["g1"], "feature": ["HOT:1"]})
    L = EP.build_lists(tab, pool, ["S1-obs", "EN-pred"], k=1)
    assert L["S1-obs"]["n_list"] == 1 and L["EN-pred"] == {"k": 1, "n_list": 0, "n_tied_at_rank_k": 0, "rows": []}


def test_manifest_verifies_complete_dataset_and_reuse(world, tmp_path, monkeypatch):
    out, d = world["root"] / "xlib_arms", world["d"]
    man = BL.verify_dataset(out, XARMS, d.genes, d.lines, stage_dir=world["root"] / "stage")
    assert man["n_rows"] == world["res"]["n_rows"] and sorted(man["parts"]) == world["res"]["parts"]
    with pytest.raises(ValueError, match="arms"):
        BL.verify_dataset(out, A.ARMS, d.genes, d.lines)
    with pytest.raises(ValueError, match="genes / lines"):
        BL.verify_dataset(out, XARMS, d.genes[:-1], d.lines)
    # a copy without the manifest (dataset written before manifests): logged counts required and checked
    import shutil
    legacy = tmp_path / "legacy"
    shutil.copytree(out, legacy)
    (legacy / BL.MANIFEST).unlink()
    n, k = world["res"]["n_rows"], len(world["res"]["parts"])
    with pytest.raises(ValueError, match="no _manifest.json"):
        BL.verify_dataset(legacy, XARMS, d.genes)
    with pytest.raises(ValueError, match="expected"):
        BL.verify_dataset(legacy, XARMS, d.genes, expect_parts=k, expect_rows=n + 1)
    assert BL.verify_dataset(legacy, XARMS, d.genes, expect_parts=k, expect_rows=n)["n_rows"] == n
    with pytest.raises(ValueError, match="stage dir"):
        BL.verify_dataset(legacy, XARMS, d.genes, stage_dir=out, expect_parts=k, expect_rows=n)
    # a changed part fails the manifest check
    bad = tmp_path / "bad"
    shutil.copytree(out, bad)
    p = sorted(bad.glob("part-*.parquet"))[0]
    p.write_bytes(p.read_bytes() + b"x")
    with pytest.raises(ValueError, match="sha256"):
        BL.verify_dataset(bad, XARMS, d.genes, d.lines)
    # run_xlib --reuse-arms: verifies, writes the manifest of a legacy dataset, refuses an incomplete one
    spec = importlib.util.spec_from_file_location("run_xlib", ROOT / "scripts" / "run_xlib.py")
    rx = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rx)
    monkeypatch.setattr(IO, "XLIB_ARMS_DIR", legacy)
    m = rx.reuse_dataset(d, tmp_path / "no_stage", expect_parts=k, expect_rows=n)
    assert (legacy / BL.MANIFEST).exists() and m["n_rows"] == n
    assert rx.reuse_dataset(d, tmp_path / "no_stage")["n_rows"] == n                # manifest now present
    monkeypatch.setattr(IO, "XLIB_ARMS_DIR", bad)
    with pytest.raises(SystemExit, match="--reuse-arms refused"):
        rx.reuse_dataset(d, tmp_path / "no_stage")
    assert np.isfinite(n)
