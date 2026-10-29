"""scripts/check_report_numbers.py (DESIGN §13.5) on a synthetic report and synthetic results in a temporary git repo:
inline citations (sentence / table-row scope, rounding to the printed precision, %, signs, ranges, fractions),
uncited numbers under --strict, structural exemptions, source tracking and registry status, Part 1 pinned tags,
run ids and exit codes. No project results are read."""
import importlib.util
import json
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("check_report_numbers", ROOT / "scripts" / "check_report_numbers.py")
C = importlib.util.module_from_spec(spec)
spec.loader.exec_module(C)

REG = ("run_id,utc,script,status,outputs\n"
       "20261029T010000_run_gate,2026-10-29T01:00:00+00:00,scripts/run_gate.py,ok,results/gates/g.json\n"
       "20261029T020000_run_q4,2026-10-29T02:00:00+00:00,scripts/run_q4.py,ok,results/clinical\n"
       "20261029T030000_replicate,2026-10-29T03:00:00+00:00,scripts/replicate.py,dryrun,results/dry/d.json\n"
       "20261029T040000_old,2026-10-29T04:00:00+00:00,scripts/x.py,failed,results/old/o.json\n")


def _git(root, *args):
    subprocess.run(["git", "-C", str(root), *args], check=True, capture_output=True)


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    r = tmp_path_factory.mktemp("report_world")
    files = {
        "results/gates/g.json": {"summary": {"overall": {"k": 21, "n": 39, "recall": 0.5384615, "ci95": [0.3718, 0.7036]},
                                             "by_universe": {"ctrl": {"k": 0, "n": 7}}},
                                 "gate_after": "False (21/39)", "beta": -0.0196, "big": 11662, "frac": 0.9437},
        "results/dry/d.json": {"x": 0.25},
        "results/old/o.json": {"x": 0.75},
        "results/unreg/u.json": {"x": 0.125},
        "configs/design_params.json": {"q": 0.05},
    }
    for rel, obj in files.items():
        (r / rel).parent.mkdir(parents=True, exist_ok=True)
        (r / rel).write_text(json.dumps(obj), encoding="utf-8")
    (r / "results" / "clinical").mkdir(parents=True)
    (r / "results" / "clinical" / "p2.csv").write_text("arm,recall,n\nEN-pred,0.4167,24\nS1-obs,0.3333,24\n",
                                                       encoding="utf-8")
    (r / "results" / "registry.csv").write_text(REG, encoding="utf-8")
    (r / "DESIGN.md").write_text("Eligible: r_within >= 0.10 in strata of >= 5 lines; k = 200; B 20 (-> 50); "
                                 "beta -0.5 shift; 18,434 genes.\n", encoding="utf-8")
    (r / "RESEARCH_TREE.md").write_text("- D15 MD-global replaces MD-strict: FDP 0.228 vs 0.013.\n", encoding="utf-8")
    (r / "results" / "deviations.csv").write_text("id,item,toy_value,final_value\nDEV-1,budget,4.5,6.25\n",
                                                  encoding="utf-8")
    (r / "configs" / "report_number_exemptions.json").write_text(
        (ROOT / "configs" / "report_number_exemptions.json").read_text(encoding="utf-8"), encoding="utf-8")
    _git(r, "init", "-q")
    _git(r, "add", "-A")
    _git(r, "-c", "user.name=t", "-c", "user.email=t@localhost", "commit", "-qm", "w")
    (r / "results" / "untracked.json").write_text(json.dumps({"x": 0.5}), encoding="utf-8")   # never added
    return r


def _run(world, body, strict=False, name="REPORT.md"):
    rep = world / name
    rep.write_text(body, encoding="utf-8")
    return C.check(rep, world, strict=strict)


def _bad(rows):
    return [r for r in rows if not r["ok"]]


def test_inline_citations_match_at_printed_precision(world):
    body = ("The gate failed at attempt 2 with 21/39 controls (54%, 95% CI 0.37-0.70) (results/gates/g.json). "
            "Recall 0.54 and 0.538 round correctly, β = −0.020 and 11,662 pairs, frac 94.4% (results/gates/g.json).\n\n"
            "- U_ctrl recovered 0 of 7 controls (results/gates/g.json).\n\n"
            "P2 recall by arm (results/clinical/p2.csv):\n\n"
            "| arm | recall | n |\n|---|---|---|\n| EN-pred | 0.42 | 24 |\n| S1-obs | 33.3% | 24 |\n")
    rows = _run(world, body, strict=True)
    assert not _bad(rows), _bad(rows)
    cited = [r for r in rows if r["kind"] == "cited"]
    assert {"21", "39", "54%", "0.37", "0.70", "−0.020", "11,662", "94.4%", "0", "7", "0.42", "33.3%"} <= {
        r["printed"] for r in cited}
    assert any(r["kind"] == "exempt" and r["printed"] == "2" for r in rows)            # "attempt 2"


@pytest.mark.parametrize("body, printed", [
    ("Recall was 0.539 (results/gates/g.json).", "0.539"),                             # 0.5385 -> 0.538 at 3 dp
    ("The effect was 0.020 (results/gates/g.json).", "0.020"),                         # sign checked
    ("Recall was 61% (results/gates/g.json).", "61%"),
    ("We saw 22/39 (results/gates/g.json).", "22"),
    ("| arm | recall |\n|---|---|\n| EN-pred | 0.43 | (results/clinical/p2.csv) |", "0.43"),
])
def test_unmatched_cited_numbers_fail(world, body, printed):
    bad = _bad(_run(world, body))
    assert [r["printed"] for r in bad] == [printed] and "not found" in bad[0]["msg"]


def test_citation_scope_is_the_sentence(world):
    body = "Recall 0.54 (results/gates/g.json). The other value 0.125 is from elsewhere."
    rows = _run(world, body, strict=True)
    assert [(r["printed"], r["kind"]) for r in _bad(rows)] == [("0.125", "uncited")]
    assert not _bad(_run(world, body, strict=False))                                  # uncited only fails --strict
    # abbreviations do not end a sentence: the citation still covers 0.54
    assert not _bad(_run(world, "EN-pred vs. S1-obs, e.g. Recall 0.54 (results/gates/g.json).", strict=True))
    # a multi-line paragraph is one sentence until its full stop
    assert not _bad(_run(world, "Recall was\n0.54 over 39\ncontrols (results/gates/g.json).", strict=True))


def test_strict_exemptions_code_and_comments(world):
    body = ("## 3. Results\n\nAs in DESIGN §5.5 and Tsherniak 2017 (run on 2026-10-29), P2 uses k = 200 and B = 20; "
            "rules: q < 0.05, >= 5 carriers, 95% CI. Code `0.123` <!--c:claim_7 0.9--> [link](https://x.org/v1.2) "
            "and S1-obs, Q4a, class 3b, commit 415c593 are labels.\n\n```\n0.555\n```\n\n1. first item\n\n"
            "| # | rule |\n|---|---|\n| 12 | §8.3 labels |\n")
    rows = _run(world, body, strict=True)
    assert not _bad(rows), _bad(rows)
    assert {r["printed"] for r in rows if r["kind"] == "exempt"} >= {"2017", "200", "20", "0.05", "5", "95%"}


def test_source_tracking_and_registry(world):
    rows = _run(world, "A 0.5 (results/untracked.json). B 0.25 (results/dry/d.json). C 0.75 (results/old/o.json). "
                       "D 0.125 (results/unreg/u.json). E 9 (results/missing.json). F 0.05 (configs/design_params.json).")
    src = {r["sources"]: r for r in rows if r["kind"] == "source"}
    assert "not tracked" in src["results/untracked.json"]["msg"] and src["results/untracked.json"]["ok"]
    assert "uncitable" in src["results/dry/d.json"]["msg"]
    assert "'failed'" in src["results/old/o.json"]["msg"]
    assert "no results/registry.csv row" in src["results/unreg/u.json"]["msg"]
    assert not src["results/missing.json"]["ok"] and "does not exist" in src["results/missing.json"]["msg"]
    assert "configs/design_params.json" not in src                                    # configs: tracked is enough
    nonstrict = _bad(rows)
    assert {r["sources"] for r in nonstrict} == {"results/missing.json"}               # warnings only, plus missing
    strict = _bad(_run(world, (world / "REPORT.md").read_text(encoding="utf-8"), strict=True))
    assert {r["sources"] for r in strict if r["kind"] == "source"} == set(src)


def test_registry_directory_outputs_and_run_ids(world):
    rows = _run(world, "Recall 0.4167 (results/clinical/p2.csv) in run 20261029T020000_run_q4; "
                       "not run 20261029T030000_replicate nor 20261029T090000_ghost.", strict=False)
    bad = _bad(rows)
    assert [r["printed"] for r in bad] == ["20261029T030000_replicate", "20261029T090000_ghost"]
    assert "dryrun" in bad[0]["msg"] and "not in results/registry.csv" in bad[1]["msg"]


def test_pinned_tags_part1_form(world):
    ok = ("<!--n-src:g=results/gates/g.json-->Recall <!--n:g#summary/overall/recall-->0.538, "
          "<!--n:g#summary/overall/recall|pct-->53.8%, <!--n:results/gates/g.json#beta|abs-->0.02.")
    assert not _bad(_run(world, ok, strict=True))
    for body, msg in [("<!--n:results/gates/g.json#summary/overall/recall-->0.54%", "pct"),
                      ("<!--n:results/gates/g.json#summary/overall/recall-->0.537", "mismatch"),
                      ("<!--n:results/gates/g.json#summary/nope-->0.5", "unresolved"),
                      ("<!--n:results/gates/g.json#gate_after-->0.5", "not a number")]:
        bad = _bad(_run(world, body))
        assert bad and msg in bad[0]["msg"], (body, bad)


def test_exemption_needs_reason(world, tmp_path):
    f = tmp_path / "ex.json"
    f.write_text(json.dumps({"exemptions": [{"id": "x", "regex": "1", "reason": " "}]}), encoding="utf-8")
    with pytest.raises(ValueError, match="no reason"):
        C.load_exemptions(world, str(f))
    for e in C.load_exemptions(ROOT):                                                  # the project config
        assert e["reason"].strip() and e["id"]


def test_main_defaults_and_exit_codes(world, capsys):
    (world / "docs").mkdir(exist_ok=True)
    (world / "REPORT.md").write_text("Gate 21/39 (results/gates/g.json).", encoding="utf-8")
    (world / "docs" / "effective_rules.md").write_text("Threshold q = 0.05 (configs/design_params.json).",
                                                       encoding="utf-8")
    assert C.main(["--root", str(world)]) == 0
    out = capsys.readouterr().out
    assert "REPORT.md" in out and "effective_rules.md" in out
    (world / "docs" / "effective_rules.md").write_text("Threshold q = 0.06 (configs/design_params.json).",
                                                       encoding="utf-8")
    assert C.main(["--root", str(world)]) == 1
    (world / "docs" / "effective_rules.md").write_text("Uncited 0.07.", encoding="utf-8")
    assert C.main(["--root", str(world)]) == 0 and C.main(["--root", str(world), "--strict"]) == 1
    js = world / "out.json"
    C.main(["--root", str(world), "--json", str(js), "--quiet"])
    assert set(json.loads(js.read_text(encoding="utf-8"))) == {str(world / "REPORT.md"),
                                                               str(world / "docs" / "effective_rules.md")}


def test_text_sources_need_verbatim_numbers(world):
    body = ("The rule is r_within 0.10 over 18434 genes (DESIGN.md §2). The switch followed FDP 0.228 "
            "(RESEARCH_TREE.md D15). The budget moved from 4.5 to 6.25 (results/deviations.csv). "
            "Shift −0.5 and k 200 (DESIGN.md; results/gates/g.json).")
    rows = _run(world, body, strict=True)
    assert not _bad(rows), _bad(rows)
    assert not [r for r in rows if r["kind"] == "source"]                             # no registry row needed
    for bad_body, printed in [("The rule is r_within 0.1 (DESIGN.md).", "0.1"),        # verbatim, not rounded
                              ("FDP 0.23 (RESEARCH_TREE.md D15).", "0.23"),
                              ("Budget 6.2 (results/deviations.csv).", "6.2"),
                              ("Strata of 4 lines (DESIGN.md).", "4"),
                              ("We used 434 genes (DESIGN.md).", "434")]:               # whole token only
        bad = _bad(_run(world, bad_body))
        assert [r["printed"] for r in bad] == [printed], (bad_body, bad)
    # json / csv sources keep numeric matching at the printed precision
    assert not _bad(_run(world, "Recall 0.5 (results/gates/g.json).", strict=True))
