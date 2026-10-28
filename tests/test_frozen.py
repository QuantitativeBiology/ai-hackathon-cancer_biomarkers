"""scripts/check_frozen.py on throwaway git repositories (DESIGN §13.4)."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("check_frozen", ROOT / "scripts" / "check_frozen.py")
CF = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(CF)

DESIGN = "# DESIGN\n\nrule one\nrule two\n"


def _git(repo: Path, *args):
    r = subprocess.run(["git", "-c", "user.name=test", "-c", "user.email=test@example.invalid", *args], cwd=str(repo),
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r.stdout


def _h(text: str) -> str:
    return hashlib.sha256(text.replace("\r\n", "\n").encode()).hexdigest()


def _w(repo: Path, rel: str, text: str, crlf: bool = False):
    p = repo / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes((text.replace("\n", "\r\n") if crlf else text).encode())


def _commit(repo: Path, msg: str):
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", msg)


@pytest.fixture()
def repo(tmp_path):
    r = tmp_path / "repo"
    r.mkdir()
    _git(r, "init", "-q")
    _git(r, "config", "core.autocrlf", "false")
    _w(r, "DESIGN.md", DESIGN, crlf=True)                      # CRLF on disk: hashes are LF-normalised
    _w(r, "configs/design_params.json", '{"a": 1}\n')
    _w(r, "configs/frozen_hashes.json", json.dumps({"design": {
        "DESIGN.md": _h(DESIGN), "configs/design_params.json": _h('{"a": 1}\n')}, "_comment": "x"}))
    _commit(r, "design")
    _git(r, "tag", "freeze-design")
    return r


def _res(repo):
    return CF.evaluate(repo)


def test_clean_design_passes(repo):
    res = _res(repo)
    assert res["fails"] == []
    assert any("design_params.json: OK" in n for n in res["notes"])
    assert CF.main(["--root", str(repo)]) == 0


def test_changed_frozen_config_fails(repo):
    _w(repo, "configs/design_params.json", '{"a": 2}\n')
    res = _res(repo)
    assert any("design_params.json: current sha256" in f for f in res["fails"])
    assert CF.main(["--root", str(repo)]) == 1


def test_decision_is_named_but_does_not_pass(repo):
    _w(repo, "configs/design_params.json", '{"a": 2}\n')
    _w(repo, "RESEARCH_TREE.md", "- D3 (M4): changed design_params.json after review; outputs not seen\n")
    f = [x for x in _res(repo)["fails"] if "design_params.json" in x][0]
    assert "D3 (M4)" in f


def test_recorded_hash_must_match_tag(repo):
    rec = json.loads((repo / "configs/frozen_hashes.json").read_text())
    rec["design"]["configs/design_params.json"] = _h('{"a": 3}\n')
    _w(repo, "configs/frozen_hashes.json", json.dumps(rec))
    _w(repo, "configs/design_params.json", '{"a": 3}\n')
    res = _res(repo)
    assert any("recorded hash != file at tag freeze-design" in f for f in res["fails"])


def test_amendments_append_only(repo):
    amended = DESIGN + "\n## Amendments after freeze (logged)\n- A1 something\n"
    _w(repo, "DESIGN.md", amended)
    assert any("DESIGN.md" in f for f in _res(repo)["fails"])           # no design_amendments record yet
    rec = json.loads((repo / "configs/frozen_hashes.json").read_text())
    rec["design_amendments"] = {"DESIGN.md": _h(amended), "amendments": ["A1"]}
    _w(repo, "configs/frozen_hashes.json", json.dumps(rec))
    assert _res(repo)["fails"] == []
    # editing text above the heading, even with a matching amendments hash, fails
    edited = DESIGN.replace("rule two", "rule 2") + "\n## Amendments after freeze (logged)\n- A1 something\n"
    _w(repo, "DESIGN.md", edited)
    rec["design_amendments"]["DESIGN.md"] = _h(edited)
    _w(repo, "configs/frozen_hashes.json", json.dumps(rec))
    assert any("text above the amendments heading differs" in f for f in _res(repo)["fails"])


def test_protected_results_before_freeze_bounds_fail(repo):
    _w(repo, "results/null/summary.json", "{}")
    _commit(repo, "null results")
    res = _res(repo)
    assert any("tag freeze-bounds does not exist" in f for f in res["fails"])


def test_freeze_bounds_ancestor_rule(repo):
    _w(repo, "results/toy/nc_bounds.json", "{}")
    rec = json.loads((repo / "configs/frozen_hashes.json").read_text())
    rec["bounds"] = {"results/toy/nc_bounds.json": _h("{}")}
    _w(repo, "configs/frozen_hashes.json", json.dumps(rec))
    _commit(repo, "bounds")
    _git(repo, "tag", "freeze-bounds")
    assert _res(repo)["fails"] == []
    _w(repo, "results/benchmark/arms.csv", "a\n")
    _commit(repo, "arms")
    res = _res(repo)
    assert res["fails"] == [] and any("is an ancestor" in n for n in res["notes"])


def test_freeze_bounds_after_protected_commit_fails(repo):
    _w(repo, "results/gates/G_pc_attempt1.json", "{}")
    _commit(repo, "gate before bounds")
    _w(repo, "results/toy/nc_bounds.json", "{}")
    rec = json.loads((repo / "configs/frozen_hashes.json").read_text())
    rec["bounds"] = {"results/toy/nc_bounds.json": _h("{}")}
    _w(repo, "configs/frozen_hashes.json", json.dumps(rec))
    _commit(repo, "bounds")
    _git(repo, "tag", "freeze-bounds")
    res = _res(repo)
    assert any("is NOT an ancestor" in f for f in res["fails"])


def test_bounds_section_requires_its_tag_and_vice_versa(repo):
    rec = json.loads((repo / "configs/frozen_hashes.json").read_text())
    rec["bounds"] = {"configs/design_params.json": _h('{"a": 1}\n')}
    _w(repo, "configs/frozen_hashes.json", json.dumps(rec))
    assert any("tag freeze-bounds missing" in f for f in _res(repo)["fails"])
    del rec["bounds"]
    _w(repo, "configs/frozen_hashes.json", json.dumps(rec))
    _git(repo, "commit", "-q", "--allow-empty", "-m", "x")
    _git(repo, "tag", "freeze-bounds")
    assert any("no 'bounds' section" in f for f in _res(repo)["fails"])


def test_coverage_reported_and_strict(repo):
    _w(repo, "configs/build_params.json", "{}")
    res = _res(repo)
    assert res["fails"] == []
    assert any("configs/build_params.json" in c and "absent at tag freeze-design" in c for c in res["coverage"])
    assert any("configs/build_params.json" in f for f in CF.evaluate(repo, strict_coverage=True)["fails"])


def test_real_repository_record_is_consistent():
    """The project's own frozen_hashes.json passes (coverage gaps are reported, not failed)."""
    res = CF.evaluate(ROOT)
    assert res["fails"] == [], res["fails"]
