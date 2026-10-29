"""Milestone budget check (DESIGN §14, configs/tiers.json ``_comment``).

    python -X utf8 scripts/budget_check.py [--milestone M4] [--dry-run] [--json out.json]

Reads every ``results/timing/<job>.json``, maps job stems to stages through ``tiers.json.stage_job_map`` and
compares each stage's wall time with its Tier A / Tier B timebox. It

- FAILS (exit 1) on a timing file whose stem is not in ``stage_job_map``, on a timing record that does not meet the
  minimal schema (docs/process_checks.md §1), and when more than ``concurrency.max_total_workers`` workers ran at
  the same time;
- lists every Tier B stage that started before all Tier A stages of its milestone were done;
- reports stages without (all of) their timing files as ``not done: time``;
- when cumulative lateness (sum over done stages of milestones up to ``--milestone`` of max(0, wall − timebox))
  exceeds ``lateness_trigger_h``, applies the next unused ``cut_order`` entry (one per milestone run; entries whose
  stage is already done are skipped) and appends it to ``results/deviations.csv``.

Timing record (minimal schema): ``job`` (= file stem), ``utc`` (ISO end time, written when the job ends),
``wall_s`` (number, whole-job wall seconds) and ``workers`` (int, worker processes the job ran, 1 if serial).
Optional: ``start_utc`` (ISO start; otherwise start = utc − wall_s) and ``run_id`` (stems written by one process share
it, so concurrency counts that process once). Legacy forms accepted with a warning: ``workers_used`` for ``workers``;
a dict ``wall_s`` (summed, plus ``load_s``).

The stage -> milestone map is not in tiers.json; it is derived from the milestone ``content`` strings and DESIGN §14
and kept here (STAGE_MILESTONE), with a check that it covers every stage of tiers.json exactly once.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TIERS = ROOT / "configs" / "tiers.json"
TIMING = ROOT / "results" / "timing"
DEVIATIONS = ROOT / "results" / "deviations.csv"
# Timing written outside results/timing on purpose (never budgeted, only counted for visibility): ctrl-universe arms
# (scripts/run_arms.py --universe ctrl -> results/benchmark_ctrl/timing/; older layout results/benchmark/timing_ctrl/)
# and dry runs (replicate.py -> results/dryrun/).
NOT_BUDGETED = (ROOT / "results" / "benchmark", ROOT / "results" / "benchmark_ctrl",
                ROOT / "results" / "benchmark_lowtier", ROOT / "results" / "dryrun")
DEV_COLS = ["id", "utc", "stage", "item", "provisional", "toy_value", "final_value", "reason", "source", "agent"]
SOURCE = "scripts/budget_check.py"
AGENT = "budget_check"
SEC_PER_MIN, MIN_PER_H = 60.0, 60.0

# Derived from tiers.json milestones[].content and DESIGN §14 "Schedule" (docs/process_checks.md §1.3).
STAGE_MILESTONE = {
    "data_status": "M3b", "cgt_features": "M3b", "dictionary": "M3b", "screen_covariates": "M3b",
    "eligibility": "M3b", "nan_mask_subsamples": "M3b", "toy_run_bounds": "M3b",
    "s1_obs_pred": "M4", "en_obs_pred": "M4", "md_strict": "M4", "confounder_audit": "M4", "gate": "M4",
    "null_fixed_B": "M4", "xlib_arms": "M4", "independent_check": "M4", "null_rerun_selection_adaptive_B10": "M4",
    "component_readout": "M5", "subtypes_gene_module": "M5", "cross_lineage": "M5",
    "colocation_cutting_culture_hypermut": "M5", "m1_proxy_attribution_panel": "M5",
    "cgt_build_freeze": "M6", "cgt_scoring_q4": "M6", "combinatorial_lookup": "M6", "tcga_prevalence_cbioportal": "M6",
    "replication_R": "M7",
    "limits_classes": "M8", "candidate_cards": "M8",
}


def _ts(s) -> dt.datetime:
    t = dt.datetime.fromisoformat(str(s).replace("Z", "+00:00"))
    return t if t.tzinfo else t.replace(tzinfo=dt.timezone.utc)


def load_tiers(path: Path = TIERS) -> dict:
    t = json.loads(Path(path).read_text(encoding="utf-8"))
    stages = {}
    for tier in ("A", "B"):
        for s, m in t["tiers"][tier]["stages"].items():
            stages[s] = {"tier": tier, "box_min": float(m)}
    t["_stages"] = stages
    t["_order"] = [m["id"] for m in t["milestones"]]
    return t


def check_stage_map(t: dict) -> list:
    """Errors in the stage tables themselves (stage_job_map / STAGE_MILESTONE vs the tier stages)."""
    errs = []
    st = set(t["_stages"])
    if set(t["stage_job_map"]) != st:
        errs.append(f"stage_job_map stages != tier stages: missing {sorted(st - set(t['stage_job_map']))}, "
                    f"extra {sorted(set(t['stage_job_map']) - st)}")
    if set(STAGE_MILESTONE) != st:
        errs.append(f"STAGE_MILESTONE != tier stages: missing {sorted(st - set(STAGE_MILESTONE))}, "
                    f"extra {sorted(set(STAGE_MILESTONE) - st)}")
    bad = {m for m in STAGE_MILESTONE.values() if m not in t["_order"]}
    if bad:
        errs.append(f"STAGE_MILESTONE names unknown milestones {sorted(bad)}")
    seen = {}
    for s, jobs in t["stage_job_map"].items():
        for j in jobs:
            if j in seen:
                errs.append(f"job stem {j!r} mapped to two stages ({seen[j]}, {s})")
            seen[j] = s
    return errs


def read_record(path: Path) -> tuple[dict | None, list, list]:
    """Normalised record {job, start, end, wall_s, workers, run_id} plus (errors, warnings)."""
    errs, warns = [], []
    stem = path.stem
    try:
        r = json.loads(path.read_text(encoding="utf-8"))
    except Exception as e:  # noqa: BLE001
        return None, [f"{path.name}: unreadable JSON ({e})"], warns
    if not isinstance(r, dict):
        return None, [f"{path.name}: not a JSON object"], warns
    if r.get("job") != stem:
        errs.append(f"{path.name}: field job={r.get('job')!r} != file stem {stem!r}")
    wall = r.get("wall_s")
    if isinstance(wall, dict):
        vals = [v for v in wall.values() if isinstance(v, (int, float))]
        if len(vals) != len(wall):
            errs.append(f"{path.name}: wall_s dict has non-numeric values")
        wall = float(sum(vals)) + float(r.get("load_s", 0.0) or 0.0)
        warns.append(f"{path.name}: legacy wall_s dict summed (+load_s) = {wall:.1f} s")
    elif isinstance(wall, (int, float)) and not isinstance(wall, bool):
        wall = float(wall)
    else:
        wall = None
    workers = r.get("workers")
    if workers is None and "workers_used" in r:
        workers = r["workers_used"]
        warns.append(f"{path.name}: legacy field workers_used used as workers")
    if not (isinstance(workers, int) and not isinstance(workers, bool) and workers >= 1):
        errs.append(f"{path.name}: missing/invalid 'workers' (int >= 1)")
        workers = None
    end = start = None
    try:
        end = _ts(r["utc"]) if "utc" in r else None
    except ValueError:
        errs.append(f"{path.name}: unparsable utc {r.get('utc')!r}")
    if end is None:
        errs.append(f"{path.name}: missing 'utc' (end time)")
    if "start_utc" in r:
        try:
            start = _ts(r["start_utc"])
        except ValueError:
            errs.append(f"{path.name}: unparsable start_utc {r.get('start_utc')!r}")
    if wall is None:
        if start is not None and end is not None:
            wall = (end - start).total_seconds()
        else:
            errs.append(f"{path.name}: missing numeric 'wall_s' (and no start_utc)")
    if start is None and end is not None and wall is not None:
        start = end - dt.timedelta(seconds=wall)
    if start is not None and end is not None and start > end:
        errs.append(f"{path.name}: start after end")
    rec = {"job": stem, "start": start, "end": end, "wall_s": wall, "workers": workers,
           "run_id": str(r.get("run_id", stem)), "subs": []}
    subs = r.get("subrecords", {})
    if not isinstance(subs, dict):
        errs.append(f"{path.name}: 'subrecords' must be an object")
        subs = {}
    for name, s in subs.items():
        sr, se = _read_sub(path.name, stem, name, s)
        errs += se
        if sr is not None:
            rec["subs"].append(sr)
    return rec, errs, warns


def _read_sub(fname: str, stem: str, name: str, s) -> tuple[dict | None, list]:
    """A sub-record {start_utc?, utc, wall_s?, workers, run_id?}: same minimal schema as a record (no job); its
    interval counts toward the parent stage's wall and toward concurrency (own run_id, default <stem>:<name>)."""
    errs = []
    if not isinstance(s, dict):
        return None, [f"{fname}: sub-record {name!r} is not an object"]
    w = s.get("workers")
    if not (isinstance(w, int) and not isinstance(w, bool) and w >= 1):
        errs.append(f"{fname}: sub-record {name!r} missing/invalid 'workers' (int >= 1)")
        w = None
    try:
        end = _ts(s["utc"])
    except (KeyError, ValueError):
        return None, errs + [f"{fname}: sub-record {name!r} missing/unparsable 'utc' (end time)"]
    start = None
    if "start_utc" in s:
        try:
            start = _ts(s["start_utc"])
        except ValueError:
            errs.append(f"{fname}: sub-record {name!r} unparsable start_utc")
    wall = s.get("wall_s")
    if start is None:
        if isinstance(wall, (int, float)) and not isinstance(wall, bool):
            start = end - dt.timedelta(seconds=float(wall))
        else:
            return None, errs + [f"{fname}: sub-record {name!r} needs numeric 'wall_s' or 'start_utc'"]
    if start > end:
        errs.append(f"{fname}: sub-record {name!r} start after end")
    return {"job": f"{stem}:{name}", "start": start, "end": end, "workers": w,
            "run_id": str(s.get("run_id", f"{stem}:{name}"))}, errs


def _intervals(rec: dict) -> list:
    """(start, end) of a record and of its sub-records."""
    iv = [(rec["start"], rec["end"])] if rec["start"] is not None and rec["end"] is not None else []
    return iv + [(s["start"], s["end"]) for s in rec.get("subs", [])]


def _union_s(intervals) -> float:
    iv = sorted(intervals)
    tot, cur_s, cur_e = 0.0, None, None
    for s, e in iv:
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                tot += (cur_e - cur_s).total_seconds()
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        tot += (cur_e - cur_s).total_seconds()
    return tot


def peak_concurrency(recs: list) -> tuple[int, str]:
    """Max over time of Σ workers of running processes (a run_id counts once, with its max workers)."""
    runs = {}
    for r in recs:
        if r["start"] is None or r["end"] is None or r["workers"] is None:
            continue
        runs.setdefault(r["run_id"], []).append(r)
    events = []
    for rid, rs in runs.items():
        s, e, w = min(x["start"] for x in rs), max(x["end"] for x in rs), max(x["workers"] for x in rs)
        if e > s:
            events.append((s, 1, w, rid))
            events.append((e, 0, w, rid))      # ends sort before starts at the same instant
    events.sort(key=lambda x: (x[0], x[1]))
    cur, peak, when, live = 0, 0, "", set()
    for t, kind, w, rid in events:
        if kind == 1:
            cur += w
            live.add(rid)
            if cur > peak:
                peak, when = cur, f"{t.isoformat()} ({'+'.join(sorted(live))})"
        else:
            cur -= w
            live.discard(rid)
    return peak, when


def not_budgeted(dirs=NOT_BUDGETED) -> dict:
    """{dir: [json files]} under the excluded timing locations (results/benchmark*/timing*/, results/dryrun/)."""
    out = {}
    for d in dirs:
        d = Path(d)
        if not d.exists():
            continue
        files = [p for p in d.rglob("*.json") if not d.name.startswith("benchmark") or p.parent.name.startswith("timing")]
        if files:
            out[str(d)] = sorted(str(p.relative_to(d)) for p in files)
    return out


def evaluate(timing_dir: Path = TIMING, tiers_path: Path = TIERS, milestone: str | None = None,
             deviations: Path = DEVIATIONS, excluded_dirs=NOT_BUDGETED) -> dict:
    """Only ``timing_dir/*.json`` (top level) is budgeted; sub-directories and ``excluded_dirs`` are not."""
    t = load_tiers(tiers_path)
    errors = check_stage_map(t)
    warnings = []
    job_stage = {j: s for s, js in t["stage_job_map"].items() for j in js}
    recs = {}
    for p in sorted(Path(timing_dir).glob("*.json")) if Path(timing_dir).exists() else []:
        if p.stem not in job_stage:
            errors.append(f"unmapped job stem {p.stem!r} ({p.name}) not in tiers.json stage_job_map")
            continue
        rec, e, w = read_record(p)
        errors += e
        warnings += w
        if rec is not None:
            recs[p.stem] = rec
    order = t["_order"]
    if milestone is None:
        hit = [STAGE_MILESTONE[job_stage[j]] for j in recs]
        milestone = max(hit, key=order.index) if hit else order[0]
    if milestone not in order:
        raise SystemExit(f"unknown milestone {milestone!r}; known {order}")
    upto = set(order[: order.index(milestone) + 1])
    # stages
    stages = {}
    for s, info in t["_stages"].items():
        jobs = t["stage_job_map"][s]
        have = [recs[j] for j in jobs if j in recs]
        done = len(have) == len(jobs) and all(r["start"] is not None and r["end"] is not None for r in have)
        row = {"stage": s, "tier": info["tier"], "milestone": STAGE_MILESTONE.get(s), "box_min": info["box_min"],
               "jobs": jobs, "jobs_present": [r["job"] for r in have], "done": done}
        if have and all(r["start"] is not None for r in have):
            row["start"] = min(s for r in have for s, _ in _intervals(r))
        if done:
            row["end"] = max(e for r in have for _, e in _intervals(r))
            row["wall_min"] = _union_s([iv for r in have for iv in _intervals(r)]) / SEC_PER_MIN
            row["late_min"] = max(0.0, row["wall_min"] - info["box_min"])
            row["status"] = "late" if row["late_min"] > 0 else "ok"
        else:
            row["status"] = "not done: time"
        stages[s] = row
    # concurrency
    peak, when = peak_concurrency(list(recs.values()) + [s for r in recs.values() for s in r.get("subs", [])])
    cmax = int(t["concurrency"]["max_total_workers"])
    if peak > cmax:
        errors.append(f"{peak} concurrent workers > max_total_workers {cmax} at {when}")
    # Tier B started before all Tier A stages of its milestone were done
    tierb_early = []
    for s, row in stages.items():
        if row["tier"] != "B" or "start" not in row:
            continue
        a_rows = [r for r in stages.values() if r["tier"] == "A" and r["milestone"] == row["milestone"]]
        pending = [r["stage"] for r in a_rows if not r["done"] or r["end"] > row["start"]]
        if pending:
            tierb_early.append({"stage": s, "milestone": row["milestone"], "start": row["start"].isoformat(),
                                "tier_a_not_done_at_start": pending})
    # lateness and cuts
    late_h = sum(r.get("late_min", 0.0) for r in stages.values() if r["milestone"] in upto) / MIN_PER_H
    trig = float(t["lateness_trigger_h"])
    cut = None
    prior = _prior_cuts(deviations)
    if late_h > trig:
        if any(p["milestone"] == milestone for p in prior):
            cut = {"action": "none", "reason": f"cut already applied at {milestone}"}
        else:
            used = {p["entry"] for p in prior}
            for i, entry in enumerate(t["cut_order"]):
                if entry in used:
                    continue
                st = entry.split(" ")[0]
                if st in stages and stages[st]["done"]:
                    warnings.append(f"cut_order[{i}] {entry!r} skipped: stage already done")
                    continue
                cut = {"action": "cut", "index": i, "entry": entry, "milestone": milestone,
                       "cumulative_lateness_h": late_h}
                break
            else:
                cut = {"action": "none", "reason": "cut_order exhausted"}
                errors.append(f"cumulative lateness {late_h:.2f} h > {trig} h and cut_order exhausted")
    return {"milestone": milestone, "upto": [m for m in order if m in upto], "errors": errors, "warnings": warnings, "stages": stages,
            "peak_workers": peak, "peak_at": when, "max_total_workers": cmax, "tier_b_early": tierb_early,
            "cumulative_lateness_h": late_h, "lateness_trigger_h": trig, "cut": cut, "n_records": len(recs),
            "not_budgeted": not_budgeted(excluded_dirs)}


def _prior_cuts(deviations: Path) -> list:
    out = []
    if not Path(deviations).exists():
        return out
    with open(deviations, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r.get("source") == SOURCE and str(r.get("id", "")).startswith("DEV-BUDGET-"):
                out.append({"milestone": r.get("stage", ""), "entry": r.get("item", "")})
    return out


def log_cut(res: dict, deviations: Path = DEVIATIONS) -> dict | None:
    cut = res.get("cut")
    if not cut or cut.get("action") != "cut":
        return None
    deviations = Path(deviations)
    deviations.parent.mkdir(parents=True, exist_ok=True)
    cols = DEV_COLS
    if deviations.exists() and deviations.stat().st_size:
        with open(deviations, newline="", encoding="utf-8") as fh:
            cols = next(csv.reader(fh))
    n = len(_prior_cuts(deviations)) + 1
    row = {"id": f"DEV-BUDGET-{n}", "utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
           "stage": cut["milestone"], "item": cut["entry"], "provisional": "", "toy_value": "",
           "final_value": f"cut_order[{cut['index']}] applied",
           "reason": f"cumulative lateness {cut['cumulative_lateness_h']:.2f} h > lateness_trigger_h "
                     f"{res['lateness_trigger_h']} at end of {cut['milestone']}",
           "source": SOURCE, "agent": AGENT}
    new = not deviations.exists() or not deviations.stat().st_size
    with open(deviations, "a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        if new:
            w.writeheader()
        w.writerow(row)
    return row


def render(res: dict) -> str:
    L = [f"budget_check: milestone {res['milestone']}, {res['n_records']} timing record(s)"]
    for s, r in res["stages"].items():
        if r["milestone"] not in res["upto"] and not r["jobs_present"]:
            continue
        w = f"{r['wall_min']:.1f}" if "wall_min" in r else "-"
        L.append(f"  [{r['milestone']}/{r['tier']}] {s:<38} box {r['box_min']:>5.0f} min  wall {w:>6}  "
                 f"{r['status']}  jobs {len(r['jobs_present'])}/{len(r['jobs'])}")
    L.append(f"peak concurrent workers: {res['peak_workers']} (max {res['max_total_workers']}) {res['peak_at']}")
    L.append(f"cumulative lateness: {res['cumulative_lateness_h']:.2f} h (trigger {res['lateness_trigger_h']} h)")
    if res["cut"]:
        L.append(f"cut: {res['cut']}")
    for x in res["tier_b_early"]:
        L.append(f"TIER-B-EARLY: {x['stage']} ({x['milestone']}) started {x['start']} before Tier A "
                 f"{x['tier_a_not_done_at_start']} were done")
    for d, fs in res.get("not_budgeted", {}).items():
        L.append(f"not budgeted ({d}): {len(fs)} timing file(s): {', '.join(fs[:8])}{' ...' if len(fs) > 8 else ''}")
    for w in res["warnings"]:
        L.append(f"WARN: {w}")
    for e in res["errors"]:
        L.append(f"FAIL: {e}")
    L.append("RESULT: " + ("FAIL" if res["errors"] else "PASS"))
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--milestone", default=None, help="milestone ending now (default: latest with a timing file)")
    ap.add_argument("--timing-dir", type=Path, default=TIMING)
    ap.add_argument("--tiers", type=Path, default=TIERS)
    ap.add_argument("--deviations", type=Path, default=DEVIATIONS)
    ap.add_argument("--dry-run", action="store_true", help="never write deviations.csv")
    ap.add_argument("--json", type=Path, default=None)
    a = ap.parse_args(argv)
    res = evaluate(a.timing_dir, a.tiers, a.milestone, a.deviations)
    if not a.dry_run:
        row = log_cut(res, a.deviations)
        if row:
            res["deviation_row"] = row
    print(render(res))
    if a.json:
        a.json.write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
    return 1 if res["errors"] else 0


if __name__ == "__main__":
    sys.exit(main())
