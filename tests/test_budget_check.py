"""scripts/budget_check.py on synthetic timing directories (real configs/tiers.json, read-only)."""
from __future__ import annotations

import csv
import datetime as dt
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("budget_check", ROOT / "scripts" / "budget_check.py")
BC = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(BC)

T0 = dt.datetime(2026, 11, 1, 8, 0, tzinfo=dt.timezone.utc)


def _write(d: Path, job: str, start_min: float, wall_min: float, workers: int = 1, explicit_start=True, **extra):
    end = T0 + dt.timedelta(minutes=start_min + wall_min)
    rec = {"job": job, "utc": end.isoformat(), "wall_s": wall_min * 60.0, "workers": workers, **extra}
    if explicit_start:
        rec["start_utc"] = (T0 + dt.timedelta(minutes=start_min)).isoformat()
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{job}.json").write_text(json.dumps(rec), encoding="utf-8")


def _m3b_all_done(d: Path, wall=1.0):
    """Every M3b Tier A job, sequential, each `wall` minutes."""
    t = 0.0
    for s, m in BC.STAGE_MILESTONE.items():
        if m != "M3b":
            continue
        for j in BC.load_tiers()["stage_job_map"][s]:
            _write(d, j, t, wall)
            t += wall
    return t


def test_stage_tables_cover_tiers_json():
    t = BC.load_tiers()
    assert BC.check_stage_map(t) == []


def test_empty_dir_passes_and_reports_not_done(tmp_path):
    res = BC.evaluate(tmp_path / "timing", milestone="M3b", deviations=tmp_path / "dev.csv")
    assert res["errors"] == []
    assert all(r["status"] == "not done: time" for r in res["stages"].values())
    assert "not done: time" in BC.render(res)


def test_unmapped_stem_fails(tmp_path):
    d = tmp_path / "timing"
    _write(d, "data_status", 0, 1)
    _write(d, "my_adhoc_job", 0, 1)
    res = BC.evaluate(d, milestone="M3b", deviations=tmp_path / "dev.csv")
    assert any("unmapped job stem 'my_adhoc_job'" in e for e in res["errors"])
    assert BC.main(["--timing-dir", str(d), "--deviations", str(tmp_path / "dev.csv"), "--dry-run"]) == 1


def test_schema_errors_fail(tmp_path):
    d = tmp_path / "timing"
    d.mkdir()
    (d / "eligibility.json").write_text(json.dumps({"job": "eligibility", "utc": T0.isoformat(), "wall_s": 5.0}),
                                        encoding="utf-8")                         # no workers
    (d / "build_dictionary.json").write_text(json.dumps({"job": "build_dictionary", "utc": T0.isoformat(), "workers": 1}),
                                       encoding="utf-8")                          # no wall_s, no start
    (d / "data_status.json").write_text(json.dumps({"job": "other", "utc": T0.isoformat(), "wall_s": 1, "workers": 1}),
                                        encoding="utf-8")                         # job != stem
    res = BC.evaluate(d, milestone="M3b", deviations=tmp_path / "dev.csv")
    errs = " ".join(res["errors"])
    assert "eligibility.json: missing/invalid 'workers'" in errs
    assert "build_dictionary.json: missing numeric 'wall_s'" in errs
    assert "data_status.json: field job='other'" in errs


def test_legacy_fields_warn_not_fail(tmp_path):
    d = tmp_path / "timing"
    d.mkdir()
    (d / "surrogate_pilot.json").write_text(json.dumps(
        {"job": "surrogate_pilot", "utc": T0.isoformat(), "workers_used": 5, "load_s": 10.0,
         "wall_s": {"a": 20.0, "b": 30.0}}), encoding="utf-8")
    res = BC.evaluate(d, milestone="M4", deviations=tmp_path / "dev.csv")
    assert res["errors"] == []
    assert len(res["warnings"]) == 2
    assert res["peak_workers"] == 5


def test_concurrency_over_12_fails_and_run_id_counts_once(tmp_path):
    d = tmp_path / "timing"
    _write(d, "null_fixed", 0, 100, workers=6)
    _write(d, "s1_obs", 10, 20, workers=6)
    res = BC.evaluate(d, milestone="M4", deviations=tmp_path / "dev.csv")
    assert res["peak_workers"] == 12 and res["errors"] == []                   # 12 is allowed
    _write(d, "confounder_delta", 15, 5, workers=1)
    res = BC.evaluate(d, milestone="M4", deviations=tmp_path / "dev.csv")
    assert res["peak_workers"] == 13
    assert any("13 concurrent workers > max_total_workers 12" in e for e in res["errors"])
    # one process writing two stems (same run_id) is counted once
    d2 = tmp_path / "t2"
    _write(d2, "s1_obs", 0, 20, workers=12, run_id="run_arms-1")
    _write(d2, "s1_pred", 0, 20, workers=12, run_id="run_arms-1")
    res = BC.evaluate(d2, milestone="M4", deviations=tmp_path / "dev.csv")
    assert res["peak_workers"] == 12 and res["errors"] == []
    # without run_id the same two files are two processes -> 24
    d3 = tmp_path / "t3"
    _write(d3, "s1_obs", 0, 20, workers=12)
    _write(d3, "s1_pred", 0, 20, workers=12)
    assert BC.evaluate(d3, milestone="M4", deviations=tmp_path / "dev.csv")["peak_workers"] == 24


def test_back_to_back_jobs_do_not_overlap(tmp_path):
    d = tmp_path / "timing"
    _write(d, "null_fixed", 0, 30, workers=12)
    _write(d, "s1_obs", 30, 10, workers=12)
    assert BC.evaluate(d, milestone="M4", deviations=tmp_path / "dev.csv")["peak_workers"] == 12


def test_start_derived_from_end_minus_wall(tmp_path):
    d = tmp_path / "timing"
    _write(d, "null_fixed", 0, 60, workers=8, explicit_start=False)
    _write(d, "en_obs", 30, 60, workers=8, explicit_start=False)
    res = BC.evaluate(d, milestone="M4", deviations=tmp_path / "dev.csv")
    assert res["peak_workers"] == 16 and res["errors"]


def test_stage_wall_is_union_and_lateness(tmp_path):
    d = tmp_path / "timing"
    _write(d, "s1_obs", 0, 25, workers=12, run_id="r")
    _write(d, "s1_pred", 0, 25, workers=12, run_id="r")                         # same scan: union = 25 min
    res = BC.evaluate(d, milestone="M4", deviations=tmp_path / "dev.csv")
    st = res["stages"]["s1_obs_pred"]
    assert st["done"] and st["wall_min"] == pytest.approx(25.0)
    assert st["late_min"] == pytest.approx(5.0) and st["status"] == "late"     # box 20
    assert res["stages"]["en_obs_pred"]["status"] == "not done: time"


def test_tier_b_started_before_tier_a_done_is_listed(tmp_path):
    d = tmp_path / "timing"
    _write(d, "component_readout", 0, 10)
    _write(d, "m1_proxy_attribution", 5, 10)                                   # Tier B of M5, starts at 5
    res = BC.evaluate(d, milestone="M5", deviations=tmp_path / "dev.csv")
    early = {x["stage"]: x for x in res["tier_b_early"]}
    assert "m1_proxy_attribution_panel" in early
    pend = set(early["m1_proxy_attribution_panel"]["tier_a_not_done_at_start"])
    assert {"component_readout", "subtypes_gene_module", "cross_lineage"} <= pend
    # all M5 Tier A done before the Tier B start -> not listed
    d2 = tmp_path / "t2"
    t = 0.0
    for s, m in BC.STAGE_MILESTONE.items():
        if m == "M5" and BC.load_tiers()["_stages"][s]["tier"] == "A":
            for j in BC.load_tiers()["stage_job_map"][s]:
                _write(d2, j, t, 1)
                t += 1
    _write(d2, "m1_proxy_attribution", t + 1, 10)
    assert BC.evaluate(d2, milestone="M5", deviations=tmp_path / "dev.csv")["tier_b_early"] == []


def test_cut_order_applied_once_per_milestone_and_logged(tmp_path):
    d, dev = tmp_path / "timing", tmp_path / "deviations.csv"
    dev.write_text(",".join(BC.DEV_COLS) + "\nDEV-TOY-NC3,x,toy_run_bounds,item,0.1,0,0.1,r,scripts/toy_nc_bounds.py,toy\n",
                   encoding="utf-8")
    t = _m3b_all_done(d)
    _write(d, "null_fixed", t, 180 + 70, workers=6)                              # 70 min late (> 1 h)
    res = BC.evaluate(d, milestone="M4", deviations=dev)
    assert res["cumulative_lateness_h"] == pytest.approx(70 / 60)
    assert res["cut"]["action"] == "cut" and res["cut"]["index"] == 0
    assert res["cut"]["entry"] == "tcga_prevalence_cbioportal"
    assert BC.main(["--timing-dir", str(d), "--deviations", str(dev), "--milestone", "M4"]) == 0
    rows = list(csv.DictReader(open(dev, encoding="utf-8")))
    assert [r["id"] for r in rows] == ["DEV-TOY-NC3", "DEV-BUDGET-1"]
    assert rows[1]["item"] == "tcga_prevalence_cbioportal" and rows[1]["stage"] == "M4"
    # same milestone again: no second cut
    BC.main(["--timing-dir", str(d), "--deviations", str(dev), "--milestone", "M4"])
    assert len(list(csv.DictReader(open(dev, encoding="utf-8")))) == 2
    # next milestone still > 1 h late: next unused entry
    res = BC.evaluate(d, milestone="M5", deviations=dev)
    assert res["cut"]["index"] == 1 and res["cut"]["entry"] == "null_rerun_selection_adaptive_B10"


def test_cut_skips_entries_whose_stage_is_done_and_no_cut_under_trigger(tmp_path):
    d, dev = tmp_path / "timing", tmp_path / "dev.csv"
    _write(d, "null_fixed", 0, 180 + 50, workers=6)                              # 50 min late: no cut
    assert BC.evaluate(d, milestone="M4", deviations=dev)["cut"] is None
    _write(d, "null_fixed", 0, 180 + 65, workers=6)
    _write(d, "tcga_prevalence", 300, 10)                                        # cut_order[0] already done
    res = BC.evaluate(d, milestone="M6", deviations=dev)
    assert res["cut"]["entry"] == "null_rerun_selection_adaptive_B10"
    assert any("cut_order[0]" in w and "already done" in w for w in res["warnings"])


def test_dry_run_never_writes(tmp_path):
    d, dev = tmp_path / "timing", tmp_path / "dev.csv"
    _write(d, "null_fixed", 0, 300, workers=6)
    assert BC.main(["--timing-dir", str(d), "--deviations", str(dev), "--milestone", "M4", "--dry-run"]) == 0
    assert not dev.exists()


def test_only_top_level_timing_is_budgeted_ctrl_and_dryrun_listed(tmp_path):
    """ctrl-universe (results/benchmark/timing_ctrl/) and dry-run (results/dryrun/) timings are never budgeted,
    even with mapped stems; sub-directories of results/timing are ignored; both are listed for visibility."""
    d = tmp_path / "results" / "timing"
    _write(d, "s1_obs", 0, 10, workers=12)
    _write(d / "old", "s1_obs", 0, 500, workers=99)                             # sub-directory: ignored
    _write(tmp_path / "results" / "benchmark" / "timing_ctrl", "s1_obs", 0, 500, workers=12)
    _write(tmp_path / "results" / "dryrun", "replicate", 0, 500, workers=12)
    (tmp_path / "results" / "benchmark" / "summary.json").write_text("{}", encoding="utf-8")   # not timing
    excl = (tmp_path / "results" / "benchmark", tmp_path / "results" / "dryrun")
    res = BC.evaluate(d, milestone="M4", deviations=tmp_path / "dev.csv", excluded_dirs=excl)
    assert res["errors"] == [] and res["n_records"] == 1 and res["peak_workers"] == 12
    nb = {Path(k).name: v for k, v in res["not_budgeted"].items()}
    assert nb == {"benchmark": [str(Path("timing_ctrl") / "s1_obs.json")], "dryrun": ["replicate.json"]}
    assert "not budgeted" in BC.render(res)


def test_new_m4_stems_are_mapped():
    """Stems written by the M4 analysis scripts (commit 82c3240) are in stage_job_map."""
    job_stage = {j for js in BC.load_tiers()["stage_job_map"].values() for j in js}
    assert {"xlib_s1_obs", "xlib_en_obs", "component_readout", "confounder_delta", "vif", "replicate",
            "q4_summary", "s1_obs", "s1_pred", "en_obs", "en_pred", "md_strict", "null_fixed",
            "surrogate_pilot"} <= job_stage
