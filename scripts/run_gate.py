"""Positive-control gate attempt (DESIGN v3.2 §9.3).

    python -B scripts/run_gate.py                                   # attempt 1
    python -B scripts/run_gate.py --reason "<bug fix>" --toy-test tests/test_x.py::test_y   # attempt n > 1

Denominator: gate-eligible rows of results/controls/pc_eligibility.csv (frozen at freeze-bounds). Recovered = (g,
primary or secondary feature) in the list with the expected sign; decided on MD-strict (results/md/md_strict.csv),
every arm reported (cache/benchmark/arms.parquet, selected rows). Writes results/gates/G_pc_attempt{n}.json (never
overwritten) and appends a row to results/gate_iterations.csv. Refuses unless HEAD descends from tag freeze-bounds.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd  # noqa: E402

from biomark import arms as A  # noqa: E402
from biomark import endpoint_io as IO  # noqa: E402
from biomark import endpoints as EP  # noqa: E402
from biomark import gitutil  # noqa: E402
from biomark import paths as P  # noqa: E402

TAG = "freeze-bounds"
GATES = P.RESULTS / "gates"
ITER = P.RESULTS / "gate_iterations.csv"


def sha(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def next_attempt() -> int:
    GATES.mkdir(parents=True, exist_ok=True)
    have = [int(f.stem.replace("G_pc_attempt", "")) for f in GATES.glob("G_pc_attempt*.json")]
    return max(have, default=0) + 1


def run(reason: str | None, toy_test: str | None, check=None) -> dict:
    from biomark import analyses as AN
    EP.require_descends(TAG, check=check)
    t_start = AN.utc_now()
    n = next_attempt()
    if n > 1 and not (reason and toy_test):
        raise EP.EndpointRefusal("attempt > 1 needs --reason and --toy-test (a failing-then-passing toy test)")
    pc = pd.read_csv(P.RESULTS / "controls" / "pc_eligibility.csv", dtype=str)
    ctl = pd.read_csv(P.CONFIGS / "positive_controls.csv", dtype=str)
    md = IO.md_strict()
    a = IO.arms_table()
    lists = {"MD-strict": md, "MD-global": IO.md_global()}     # decided on MD-strict (A7); MD-global reported
    for arm in A.ARMS:
        lists[arm] = a[(a["arm"] == arm) & a["selected"].astype(bool)][["gene", "feature", "beta"]]
    res = EP.gate_attempt(pc, ctl, lists)
    commit = gitutil.head(P.PROJECT)
    out = {"attempt": n, "utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"), "commit": commit,
           "reason": reason or "attempt 1", "summary": res["summary"],
           "rows": res["rows"].to_dict("records"),
           "inputs_sha256": {"pc_eligibility.csv": sha(P.RESULTS / "controls" / "pc_eligibility.csv"),
                             "positive_controls.csv": sha(P.CONFIGS / "positive_controls.csv"),
                             "md_strict.csv": sha(P.RESULTS / "md" / "md_strict.csv"),
                             "arms.parquet": sha(P.CACHE / "benchmark" / "arms.parquet")}}
    f = GATES / f"G_pc_attempt{n}.json"
    if f.exists():
        raise EP.EndpointRefusal(f"{f.name} exists and is never overwritten")
    IO.write_json(f, out)
    prev = None
    if n > 1:
        prev = IO.read_json(GATES / f"G_pc_attempt{n - 1}.json")["summary"]
    row = {"attempt": n, "utc": out["utc"], "commit": commit, "reason": out["reason"], "toy_test": toy_test or "",
           "gate_before": "" if prev is None else f"{prev['pass']} ({prev['overall']['k']}/{prev['overall']['n']})",
           "gate_after": f"{res['summary']['pass']} ({res['summary']['overall']['k']}/{res['summary']['overall']['n']})",
           "attempt_file_sha256": sha(f)}
    it = pd.read_csv(ITER) if ITER.exists() else pd.DataFrame(columns=list(row))
    it = pd.concat([it, pd.DataFrame([row])], ignore_index=True)
    it.to_csv(ITER, index=False)
    AN.write_timing("gate_pc", {"attempt": n}, start_utc=t_start, workers=1)
    AN.append_registry("scripts/run_gate.py", "ok", [f"results/gates/{f.name}", "results/gate_iterations.csv"])
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reason", default=None)
    ap.add_argument("--toy-test", default=None)
    a = ap.parse_args(argv)
    out = run(a.reason, a.toy_test)
    s = out["summary"]
    print(f"gate attempt {out['attempt']}: pass={s['pass']} overall {s['overall']['k']}/{s['overall']['n']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
