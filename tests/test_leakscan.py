"""scripts/leakscan.py (DESIGN §13.4, §13.2 "leakscan (planted "KRAS" and 0.7 give two hits)")."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("leakscan", ROOT / "scripts" / "leakscan.py")
LS = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(LS)

SYMS = {"KRAS", "TP53", "ERBB2"}


def _tree(tmp_path: Path, files: dict) -> Path:
    for rel, text in files.items():
        p = tmp_path / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding="utf-8")
    return tmp_path


def _kinds(hits):
    return sorted((h["kind"], h["value"]) for h in hits)


def test_planted_kras_and_0_7_give_two_hits(tmp_path):
    root = _tree(tmp_path, {"src/biomark/planted.py": 'GENE = "KRAS"\nTHR = 0.7\n'})
    hits = LS.scan(root, symbols=SYMS)
    assert _kinds(hits) == [("gene", "KRAS"), ("num", "0.7")]


def test_planted_with_real_gene_map(tmp_path):
    try:
        syms = LS.gene_symbols()
    except Exception as e:  # noqa: BLE001  (Part 1 gene_map not reachable)
        pytest.skip(f"gene_map unavailable: {e}")
    assert "KRAS" in syms and all(len(s) >= 3 for s in syms)
    root = _tree(tmp_path, {"src/biomark/planted.py": 'GENE = "KRAS"\nTHR = 0.7\n'})
    assert _kinds(LS.scan(root, symbols=syms)) == [("gene", "KRAS"), ("num", "0.7")]


def test_allowed_numbers_bools_and_negation(tmp_path):
    src = ("A = (0, 1, -1, 2, 0.5, 100, 1.0, 100.0, True, False)\n"
           "B = -2\n"
           "C = -0.5\n"
           "D = 1e-3\n"
           "E = +3\n")
    hits = LS.scan(_tree(tmp_path, {"src/biomark/m.py": src}), symbols=SYMS)
    assert _kinds(hits) == [("num", "-0.5"), ("num", "-2"), ("num", "0.001"), ("num", "3")]
    assert [h["line"] for h in sorted(hits, key=lambda h: h["line"])] == [2, 3, 4, 5]


def test_literal_ok_exempts_numbers_only(tmp_path):
    src = ('X = 0.7  # literal-ok: plotting margin\n'
           'Y = 0.7  # literal-ok:\n'                       # empty reason: not an exemption
           'Z = "TP53"  # literal-ok: covariate name\n')
    hits = LS.scan(_tree(tmp_path, {"src/biomark/m.py": src}), symbols=SYMS)
    assert _kinds(hits) == [("gene", "TP53"), ("num", "0.7")]
    assert {h["line"] for h in hits} == {2, 3}
    assert [h for h in hits if h["kind"] == "gene"][0]["literal_ok_line"]


def test_params_loader_and_scripts_exempt_from_numbers_not_strings(tmp_path):
    root = _tree(tmp_path, {
        "src/biomark/params.py": 'A = 0.49\nG = "ERBB2"\n',
        "scripts/run.py": 'B = 0.7\nG = "KRAS"\nL = "ACH-000123"\n',
        "src/other/x.py": 'C = 0.7\n',                     # outside src/biomark: numbers not scanned
    })
    hits = LS.scan(root, symbols=SYMS)
    got = sorted((h["file"], h["kind"], h["value"]) for h in hits)
    assert got == [("scripts/run.py", "ach", "ACH-000123"), ("scripts/run.py", "gene", "KRAS"),
                   ("src/biomark/params.py", "gene", "ERBB2")]


def test_string_rules_exact_symbol_and_ach_search(tmp_path):
    src = ('a = "kras"\n'               # case differs: not equal
           'b = "AMP:KRAS"\n'           # not equal to a symbol
           'c = "TP5"\n'
           'd = f"x {a} ACH-000001 y"\n'  # f-string constant part
           'e = "ACH-00001"\n'          # 5 digits: no hit
           'g = ["ERBB2"]\n')
    hits = LS.scan(_tree(tmp_path, {"scripts/s.py": src}), symbols=SYMS | {"TP5"})
    assert _kinds(hits) == [("ach", "ACH-000001"), ("gene", "ERBB2"), ("gene", "TP5")]


def test_short_symbols_never_loaded():
    """gene_symbols() keeps only symbols of >= 3 characters (DESIGN §13.4)."""
    assert LS.MIN_SYMBOL_LEN == 3


def test_cli_exit_codes(tmp_path, monkeypatch):
    monkeypatch.setattr(LS, "gene_symbols", lambda: SYMS)
    clean = _tree(tmp_path / "clean", {"src/biomark/m.py": "A = 1\n"})
    assert LS.main(["--root", str(clean)]) == 0
    dirty = _tree(tmp_path / "dirty", {"src/biomark/m.py": 'A = "KRAS"\n'})
    assert LS.main(["--root", str(dirty)]) == 1


def test_every_literal_ok_marker_is_listed(tmp_path):
    src = ('X = 0.7  # literal-ok: margin\n'
           'Y = "a"  # literal-ok: suppresses nothing\n'
           'Z = 0.3  # literal-ok\n'                      # malformed: no reason, exempts nothing
           'S = "# literal-ok: inside a string is not a marker"\n')
    ex = []
    hits = LS.scan(_tree(tmp_path, {"src/biomark/m.py": src}), symbols=SYMS, exemptions=ex)
    assert _kinds(hits) == [("num", "0.3")]
    by = {m["line"]: m for m in ex}
    assert set(by) == {1, 2, 3}
    assert by[1]["suppressed"] == ["0.7"] and by[1]["reason"] == "margin"
    assert by[2]["suppressed"] == [] and by[2]["well_formed"]
    assert not by[3]["well_formed"]
    out = LS.render(hits, ex)
    assert "EXEMPT src/biomark/m.py:1: suppressed ['0.7']" in out and "MALFORMED" in out


def _allow(root, entries):
    import json
    p = root / "configs" / "leakscan_allow.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps({"allow": entries}), encoding="utf-8")


def test_allowlisted_string_is_exempt_only_in_its_file(tmp_path):
    root = _tree(tmp_path, {"src/biomark/a.py": 'FAM = "TP53"\n', "src/biomark/b.py": 'G = "TP53"\n',
                            "scripts/c.py": 'H = "TP53"\n'})
    _allow(root, [{"file": "src/biomark/a.py", "literal": "TP53", "reason": "covariate column name"}])
    allowed, errors = [], []
    hits = LS.scan(root, symbols=SYMS, allowed=allowed, errors=errors)
    assert errors == []
    assert sorted(h["file"] for h in hits) == ["scripts/c.py", "src/biomark/b.py"]   # elsewhere still fails
    assert [(a["file"], a["value"], a["reason"]) for a in allowed] == [("src/biomark/a.py", "TP53", "covariate column name")]
    out = LS.render(hits, [], allowed, errors)
    assert "EXEMPT src/biomark/a.py:1:6: [gene] TP53 allowlisted  reason: covariate column name" in out
    assert out.endswith("RESULT: FAIL")


def test_allowlist_is_exact_literal_and_never_numbers_or_ach(tmp_path):
    root = _tree(tmp_path, {"src/biomark/a.py": 'A = "TP53"\nB = "KRAS"\nC = 0.7\nD = "ACH-000001"\n'})
    _allow(root, [{"file": "src/biomark/a.py", "literal": "TP53", "reason": "x"},
                  {"file": "src/biomark/a.py", "literal": "0.7", "reason": "numbers are not strings"}])
    allowed, errors = [], []
    hits = LS.scan(root, symbols=SYMS, allowed=allowed, errors=errors)
    assert _kinds(hits) == [("ach", "ACH-000001"), ("gene", "KRAS"), ("num", "0.7")]
    assert [a["literal"] for a in allowed if a.get("unused")] == ["0.7"]
    assert "matched nothing" in LS.render(hits, [], allowed, errors)


def test_invalid_allowlist_entries_fail(tmp_path, monkeypatch):
    root = _tree(tmp_path, {"src/biomark/a.py": 'A = "TP53"\n'})
    _allow(root, [{"file": "src/biomark/a.py", "literal": "TP53", "reason": "  "},
                  {"file": "src/biomark/a.py", "literal": "ACH-000001", "reason": "no"},
                  {"file": "src/biomark/gone.py", "literal": "TP53", "reason": "file absent"},
                  {"file": "src/biomark/a.py", "literal": "TP53"}])
    errors = []
    hits = LS.scan(root, symbols=SYMS, errors=errors)
    assert len(hits) == 1                                  # nothing valid exempts the TP53 hit
    joined = " ".join(errors)
    assert "empty reason" in joined and "ACH- ids cannot be allowlisted" in joined
    assert "does not exist" in joined and "needs string fields" in joined
    clean = _tree(tmp_path / "clean", {"src/biomark/a.py": "A = 1\n"})
    _allow(clean, [{"file": "src/biomark/a.py", "literal": "TP53", "reason": ""}])
    monkeypatch.setattr(LS, "gene_symbols", lambda: SYMS)
    assert LS.main(["--root", str(clean)]) == 1            # no hits, but an invalid entry fails


def test_project_allowlist_is_valid_and_has_no_gene_references():
    entries, errors = LS.load_allow(ROOT / "configs" / "leakscan_allow.json")
    assert errors == [] and entries
    assert all((ROOT / e["file"]).exists() for e in entries)
    assert not any(e["file"] == "scripts/build_cgt.py" for e in entries)   # real gene uses move to config
