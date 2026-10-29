"""Registry rows for citable outputs (DESIGN §13.5): run_arms.py registers its results files per universe, and
scripts/backfill_registry.py writes rows for earlier arms runs from their own metadata only (status ok only with a
recorded run commit that exists and MD row counts that agree). Synthetic files in a temporary git repo; no real data."""
import importlib.util
import json
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from biomark import paths as P

ROOT = Path(__file__).resolve().parents[1]


def _script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _git(root, *args):
    return subprocess.run(["git", "-C", str(root), *args], check=True, capture_output=True, text=True).stdout.strip()


def _arms_outputs(root, universe, commit, n_global=2, n_strict=1, utc="2026-10-29T03:33:56+00:00"):
    sfx = "" if universe == "eligible" else f"_{universe}"
    rb, md = root / "results" / f"benchmark{sfx}", root / "results" / f"md{sfx}"
    td = root / "results" / ("timing" if universe == "eligible" else f"benchmark{sfx}/timing")
    for d in (rb, md, td):
        d.mkdir(parents=True, exist_ok=True)
    (rb / "summary.json").write_text(json.dumps({
        "universe": universe, "commit": commit, "n_md_global": n_global, "n_md_strict": n_strict,
        "config_sha256": {"design_params.json": "d" * 64, "seeds.json": "s" * 64}}), encoding="utf-8")
    for f in ("selected.csv", "lambda.csv"):
        (rb / f).write_text("gene,feature\n1,HOT:1\n", encoding="utf-8")
    (rb / "p_cut.json").write_text("{}", encoding="utf-8")
    (md / "md_global.csv").write_text("gene,feature\n" + "".join(f"{i},HOT:{i}\n" for i in range(2)), encoding="utf-8")
    (md / "md_strict.csv").write_text("gene,feature\n1,HOT:1\n", encoding="utf-8")
    (td / "md_strict.json").write_text(json.dumps({"job": "md_strict", "utc": utc, "run_id": "run_arms-1-x"}),
                                       encoding="utf-8")


@pytest.fixture
def repo(tmp_path):
    _git(tmp_path, "init", "-q")
    (tmp_path / "README").write_text("x", encoding="utf-8")
    _git(tmp_path, "add", "-A")
    _git(tmp_path, "-c", "user.name=t", "-c", "user.email=t@localhost", "commit", "-qm", "c")
    return tmp_path, _git(tmp_path, "rev-parse", "HEAD")


def test_backfill_ok_only_with_recorded_commit_and_matching_counts(repo, capsys):
    root, head = repo
    bf = _script("backfill_registry")
    _arms_outputs(root, "eligible", head)
    _arms_outputs(root, "ctrl", "", utc="2026-10-29T04:36:41+00:00")                 # no recorded commit
    assert bf.main([], root=root) == 0
    reg = pd.read_csv(root / "results" / "registry.csv", dtype=str).set_index("run_id")
    e = reg.loc["20261029T033356_run_arms_backfill"]
    assert e["status"] == "ok" and e["commit"] == head and "[backfilled" in e["script"]
    assert set(e["outputs"].split(";")) == {"results/benchmark/selected.csv", "results/benchmark/lambda.csv",
                                            "results/benchmark/p_cut.json", "results/benchmark/summary.json",
                                            "results/md/md_global.csv", "results/md/md_strict.csv"}
    assert pd.isna(e["data_status_sha256"]) and pd.isna(e["dirty"])                   # unknown fields stay empty
    c = reg.loc["20261029T043641_run_arms_ctrl_backfill"]
    assert c["status"] == "backfilled_no_commit"
    assert all(o.startswith(("results/benchmark_ctrl/", "results/md_ctrl/")) for o in c["outputs"].split(";"))
    n = len(reg)
    assert bf.main([], root=root) == 0 and len(pd.read_csv(root / "results" / "registry.csv")) == n   # idempotent


@pytest.mark.parametrize("commit, n_global, status", [("deadbeef" * 5, 2, "backfilled_no_commit"),
                                                      (None, 3, "backfilled_mismatch")])
def test_backfill_never_ok_without_valid_metadata(repo, commit, n_global, status):
    root, head = repo
    _arms_outputs(root, "eligible", commit or head, n_global=n_global)
    row = _script("backfill_registry").backfill_row(
        "eligible", root, _script("run_arms").registry_outputs("eligible", existing_only=False))
    assert row["status"] == status


def test_backfill_dry_run_writes_nothing(repo):
    root, head = repo
    _arms_outputs(root, "eligible", head)
    assert _script("backfill_registry").main(["--dry-run"], root=root) == 0
    assert not (root / "results" / "registry.csv").exists()


def test_run_arms_registers_its_outputs(tmp_path, monkeypatch):
    ra = _script("run_arms")
    monkeypatch.setattr(P, "PROJECT", tmp_path)
    monkeypatch.setattr(P, "RESULTS", tmp_path / "results")
    for u in ("eligible", "ctrl"):
        _arms_outputs(tmp_path, u, "abc")
        ra.register_outputs(u)
    reg = pd.read_csv(tmp_path / "results" / "registry.csv", dtype=str)
    assert list(reg["status"]) == ["ok", "ok"]
    assert reg["script"].tolist() == ["scripts/run_arms.py", "scripts/run_arms.py --universe ctrl"]
    assert "results/md/md_global.csv" in reg["outputs"].iloc[0].split(";")
    assert all(o.startswith(("results/benchmark_ctrl/", "results/md_ctrl/")) for o in reg["outputs"].iloc[1].split(";"))


@pytest.mark.parametrize("script, needle", [
    ("run_null", "append_registry"), ("run_xlib", "append_registry"), ("run_confounder_audit", "append_registry"),
    ("run_component_readout", "append_registry"), ("run_q4", "append_registry"), ("build_p1_pool", "append_registry"),
    ("replicate", 'append_registry("scripts/replicate.py", "ok"'), ("run_arms", "register_outputs(universe)"),
    ("run_cross_lineage", "ws.register("), ("run_m5_audits", "ws.register("), ("run_limits", "ws.register("),
    ("build_candidate_cards", "ws.register(")])
def test_results_writing_scripts_register(script, needle):
    assert needle in (ROOT / "scripts" / f"{script}.py").read_text(encoding="utf-8")
