"""budget_check counts timing sub-records (e.g. null_fixed.subrecords.s1_vectorised) toward the stage wall and
concurrency; every scripted write_timing call passes workers."""
import ast
import datetime as dt
import importlib.util
import json
from pathlib import Path

import pytest

from biomark import analyses as AN

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("budget_check", ROOT / "scripts" / "budget_check.py")
BC = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(BC)
T0 = dt.datetime(2026, 11, 1, 8, 0, tzinfo=dt.timezone.utc)


def _iso(m):
    return (T0 + dt.timedelta(minutes=m)).isoformat()


def _null(d, subs=None, workers=6):
    rec = {"job": "null_fixed", "start_utc": _iso(0), "utc": _iso(60), "wall_s": 3600.0, "workers": workers, "run_id": "main"}
    if subs is not None:
        rec["subrecords"] = subs
    d.mkdir(parents=True, exist_ok=True)
    (d / "null_fixed.json").write_text(json.dumps(rec), encoding="utf-8")


def test_subrecord_interval_counts_toward_stage_wall(tmp_path):
    d = tmp_path / "timing"
    _null(d)
    base = BC.evaluate(d, milestone="M4", deviations=tmp_path / "dev.csv")["stages"]["null_fixed_B"]
    assert base["wall_min"] == pytest.approx(60.0)
    _null(d, {"s1_vectorised": {"start_utc": _iso(90), "utc": _iso(150), "wall_s": 3600.0, "workers": 6, "run_id": "s1v"}})
    res = BC.evaluate(d, milestone="M4", deviations=tmp_path / "dev.csv")
    st = res["stages"]["null_fixed_B"]
    assert st["wall_min"] == pytest.approx(120.0)                      # union of [0,60] and [90,150]
    assert st["end"] == T0 + dt.timedelta(minutes=150)
    assert not res["errors"]


def test_subrecord_counts_toward_concurrency(tmp_path):
    d = tmp_path / "timing"
    _null(d, {"s1_vectorised": {"start_utc": _iso(30), "utc": _iso(50), "workers": 7, "run_id": "s1v"}})
    res = BC.evaluate(d, milestone="M4", deviations=tmp_path / "dev.csv")
    assert res["peak_workers"] == 13 and any("max_total_workers" in e for e in res["errors"])
    _null(d, {"s1_vectorised": {"start_utc": _iso(30), "utc": _iso(50), "workers": 7}})   # default own run id
    assert BC.evaluate(d, milestone="M4", deviations=tmp_path / "dev.csv")["peak_workers"] == 13
    _null(d, {"s1_vectorised": {"start_utc": _iso(30), "utc": _iso(50), "workers": 6, "run_id": "main"}})
    assert BC.evaluate(d, milestone="M4", deviations=tmp_path / "dev.csv")["peak_workers"] == 6   # same process


def test_subrecord_schema_errors(tmp_path):
    d = tmp_path / "timing"
    _null(d, {"s1_vectorised": {"utc": _iso(50), "wall_s": 60.0}})              # no workers
    assert any("sub-record" in e and "workers" in e for e in BC.evaluate(d, milestone="M4", deviations=tmp_path / "dev.csv")["errors"])
    _null(d, {"s1_vectorised": {"utc": _iso(50), "workers": 2}})                # no wall_s / start
    assert any("sub-record" in e for e in BC.evaluate(d, milestone="M4", deviations=tmp_path / "dev.csv")["errors"])


def test_timing_subrecord_helper_matches_budget_schema(tmp_path):
    s = AN.timing_subrecord("s1_vectorised", T0, T0 + dt.timedelta(minutes=20), 4, B=200)
    d = tmp_path / "timing"
    _null(d, {"s1_vectorised": s})
    res = BC.evaluate(d, milestone="M4", deviations=tmp_path / "dev.csv")
    assert not res["errors"] and res["stages"]["null_fixed_B"]["wall_min"] == pytest.approx(60.0)
    with pytest.raises(ValueError):
        AN.timing_subrecord("x", T0, T0, 0)


def test_every_scripted_write_timing_passes_workers():
    missing = []
    for f in (ROOT / "scripts").glob("*.py"):
        for n in ast.walk(ast.parse(f.read_text(encoding="utf-8"))):
            if isinstance(n, ast.Call) and getattr(n.func, "attr", getattr(n.func, "id", "")) == "write_timing":
                kws = {k.arg for k in n.keywords}
                rec = n.args[1] if len(n.args) > 1 else None
                in_rec = isinstance(rec, ast.Dict) and any(isinstance(k, ast.Constant) and k.value == "workers"
                                                          for k in rec.keys)
                passthrough = None in kws                              # **kw forwarding wrapper (run_arms)
                if not (kws & {"workers"} or in_rec or passthrough or isinstance(rec, ast.Name)):
                    missing.append(f"{f.name}:{n.lineno}")
    assert not missing, missing
