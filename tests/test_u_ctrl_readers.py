"""U_ctrl readers (DESIGN v3.2 §2): controls of non-eligible genes are read from the U_ctrl run where the design uses
controls (§8.3 cross-lineage, §12.3 known_recovery), and U_ctrl never reaches a top-k list, the null subset, P1 or P2.
Toy world of tests/u_ctrl_world.py (biomark.toysim; no real data)."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import u_ctrl_world
from biomark import arms as A
from biomark import endpoint_io as IO
from biomark import m5 as M5

ROOT = Path(__file__).resolve().parents[1]


def _script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    w = u_ctrl_world.build(tmp_path_factory.mktemp("u_ctrl_readers"))
    ws = M5.Workspace(w["root"])
    return dict(w, ws=ws, src=M5.source_for(ws))


def _key(w):
    return (w["moved"]["dependency_entrez"], w["moved"]["feature_primary"])


def test_workspace_reads_u_ctrl_only_on_request(world):
    ws, g = world["ws"], world["moved"]["dependency_entrez"]
    assert ws.has_ctrl() and g in ws.ctrl_genes()
    assert not (ws.arms(arms=["S1-obs"], columns=["gene"])["gene"] == g).any()          # default: eligible only
    assert (ws.arms(arms=["S1-obs"], columns=["gene"], universe="ctrl")["gene"] == g).any()
    assert _key(world) in set(zip(*[M5.catalogue(ws, "ctrl")[c] for c in ("gene", "feature")]))
    assert _key(world) not in set(zip(*[M5.catalogue(ws)[c] for c in ("gene", "feature")]))


def test_analysis_pairs_add_u_ctrl_controls_only_when_asked(world):
    ws, src = world["ws"], world["src"]
    _, _, p0 = M5.analysis_pairs(ws, src)
    assert "universe" not in p0 and _key(world) not in set(zip(p0["gene"], p0["feature"]))
    _, _, p1 = M5.analysis_pairs(ws, src, u_ctrl_controls=True)
    k = p1.set_index(["gene", "feature"])
    assert k.loc[_key(world), "universe"] == "ctrl" and k.loc[_key(world), "source"] == "control"
    e = p1[p1["universe"] == "eligible"].drop(columns="universe").reset_index(drop=True)
    pd.testing.assert_frame_equal(e, p0.reset_index(drop=True))                      # eligible pairs unchanged


def test_cross_lineage_includes_u_ctrl_control(world):
    ws, src = world["ws"], world["src"]
    r = _script("run_cross_lineage").run(ws, src, 1)
    lab = r["labels"].set_index(["gene", "feature"])
    assert lab.loc[_key(world), "universe"] == "ctrl"
    assert r["summary"]["u_ctrl_controls"] == {"included": True, "n_pairs": 1}
    bsh = M5.shrunk_beta_ws(ws, {_key(world)})
    assert np.isfinite(bsh[_key(world)])                                              # U_ctrl S1-obs row used


def test_known_recovery_reads_u_ctrl_ranks_and_lists(world):
    ws, src = world["ws"], world["src"]
    pc = pd.DataFrame(columns=["gene", "feature", "pair_class", "in_md_global", "in_md_strict"])
    ws.write_csv(pc, "results/limits/pair_classes.csv")
    r = _script("build_candidate_cards").run(ws, src)
    kr = ws.read_csv("results/candidates/known_recovery.csv").set_index("id")
    mv = world["moved"]["control_id"]
    assert kr.at[mv, "universe"] == "ctrl"
    assert str(kr.at[mv, "recovered_MD-strict"]).lower() == "true"
    for arm in A.ARMS:
        assert np.isfinite(float(kr.at[mv, f"rank_{arm}"])), arm                      # rank within U_ctrl
    kept = world["kept"]["control_id"].tolist()
    assert (kr.loc[kept, "universe"] == "eligible").all()
    summ = ws.read_json("results/candidates/cards_summary.json")
    assert summ["u_ctrl_run_included"] and summ["by_universe"]["ctrl"]["recovered_md_strict"] >= 1


@pytest.mark.parametrize("script", ["run_null", "build_p1_pool", "run_q4", "replicate", "build_null_genes"])
def test_top_k_null_p1_p2_never_read_u_ctrl(script):
    """The null subset, P1 (pool / lists) and P2 (Q4) read the eligible dataset only (DESIGN §2)."""
    txt = (ROOT / "scripts" / f"{script}.py").read_text(encoding="utf-8")
    for bad in ("arms_ctrl", "md_ctrl", "benchmark_ctrl", "u_ctrl_genes", 'universe="ctrl"', "universe='ctrl'"):
        assert bad not in txt, (script, bad)
    assert IO.ARMS_DIR == IO.arms_dir("eligible") and IO.ARMS_DIR.name == "arms"     # default reader = eligible
