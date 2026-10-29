"""Gate bug fix (DESIGN v3.2 §2, §9.3): a gate-eligible positive control whose gene fails eligibility is analysed in
the separate universe U_ctrl and must be read from there by scripts/run_gate.py.

Toy world (biomark.toysim; no real data): a planted ADD-HOT control gene is declared non-eligible, so the eligible
arms run never sees it; run on U_ctrl (that gene + the toy's non-eligible LINEAGE-CHECK genes, BH within U_ctrl) it is
recovered on MD-strict. Inputs are written in the real on-disk layout (results/md/ vs results/md_ctrl/,
cache/benchmark/arms/ vs cache/benchmark/arms_ctrl/) and run_gate.run() is called on them. Before the fix run_gate read
the eligible lists only, so the control counted as not recovered."""
import importlib.util
from pathlib import Path

import pandas as pd
import pytest

import u_ctrl_world
from biomark import analyses as AN
from biomark import endpoint_io as IO
from biomark import endpoints as EP
from biomark import params as PR
from biomark import paths as P

ROOT = Path(__file__).resolve().parents[1]


def _script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    return u_ctrl_world.build(tmp_path_factory.mktemp("gate_u_ctrl"))


@pytest.fixture
def gate(world, monkeypatch):
    for f in (PR.dp, PR.pp, PR.en_params, PR.endpoints, PR.stats_params):
        f()                                                               # frozen params cached before CONFIGS moves
    root = world["root"]
    monkeypatch.setattr(P, "CONFIGS", root / "configs")
    monkeypatch.setattr(P, "RESULTS", root / "results")
    monkeypatch.setattr(P, "CACHE", root / "cache")
    monkeypatch.setattr(IO, "ARMS_DIR", root / "cache" / "benchmark" / "arms")
    monkeypatch.setattr(AN, "write_timing", lambda *a, **k: None)
    monkeypatch.setattr(AN, "append_registry", lambda *a, **k: None)
    for f in (root / "results" / "gates").glob("*.json") if (root / "results" / "gates").exists() else []:
        f.unlink()
    (root / "results" / "gate_iterations.csv").unlink(missing_ok=True)
    return _script("run_gate")


def test_gate_recovers_non_eligible_control_via_u_ctrl(world, gate):
    out = gate.run(None, None, check=lambda tag, repo=None: True)
    rows = pd.DataFrame(out["rows"]).set_index("control_id")
    mv = world["moved"]["control_id"]
    assert bool(rows.at[mv, "recovered"]), "non-eligible control must be read from U_ctrl (MD-strict)"
    assert bool(rows.at[mv, "recovered_MD-strict"]) and bool(rows.at[mv, "recovered_MD-global"])
    assert rows.at[mv, "universe"] == "ctrl"
    kept = world["kept"]["control_id"].tolist()
    assert (rows.loc[kept, "universe"] == "eligible").all()
    md = world["runs"]["eligible"]["md_strict"]                         # eligible rows still from the eligible run
    for cid, r in world["kept"].set_index("control_id").iterrows():
        exp = ((md["gene"].astype(str) == r["dependency_entrez"]) & (md["feature"] == r["feature_primary"])
               & ((md["beta"] < 0) if r["expected_sign"] == "S" else (md["beta"] > 0))).any()
        assert bool(rows.at[cid, "recovered"]) == bool(exp), cid
    h = out["inputs_sha256"]
    assert {"arms_dataset", "arms_dataset_ctrl", "md_strict.csv", "md_ctrl/md_strict.csv", "u_ctrl_genes.txt"} <= set(h)
    assert h["arms_dataset"] != h["arms_dataset_ctrl"]
    assert out["summary"]["by_universe"]["ctrl"]["k"] == 1


def test_gate_refuses_without_u_ctrl_outputs(world, gate, tmp_path):
    md = world["root"] / "results" / "md_ctrl" / "md_strict.csv"
    keep = tmp_path / "md_strict.csv"
    md.rename(keep)
    try:
        with pytest.raises(EP.EndpointRefusal, match="U_ctrl outputs missing"):
            gate.run(None, None, check=lambda tag, repo=None: True)
        assert not (world["root"] / "results" / "gates" / "G_pc_attempt1.json").exists()
    finally:
        keep.rename(md)


def test_gate_refuses_control_gene_in_no_universe(world, gate):
    f = world["root"] / "configs" / "u_ctrl_genes.txt"
    old = f.read_text(encoding="utf-8")
    f.write_text("".join(f"{g}\n" for g in old.split() if g != world["moved"]["dependency_entrez"]), encoding="utf-8")
    try:
        with pytest.raises(EP.EndpointRefusal, match="neither the eligible set nor U_ctrl"):
            gate.run(None, None, check=lambda tag, repo=None: True)
    finally:
        f.write_text(old, encoding="utf-8")
