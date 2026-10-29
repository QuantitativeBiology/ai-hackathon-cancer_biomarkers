"""One-off: results/registry.csv rows for full arms runs made before run_arms.py registered its outputs (DESIGN §13.5).

    python -B scripts/backfill_registry.py --dry-run     # print the rows only
    python -B scripts/backfill_registry.py               # append them (idempotent: an existing run_id is skipped)

Covers the eligible universe (results/benchmark/*, results/md/*) and U_ctrl (results/benchmark_ctrl/*,
results/md_ctrl/*): the outputs are the files of run_arms.registry_outputs(universe). Every field comes from the run's
own metadata, never from the current state: commit and config hashes from results/benchmark{,_ctrl}/summary.json, the
end time and process run id from its md_strict timing JSON (the global step, the run's last stage). Unknown fields
(dirty flag, data_status hash) stay empty. The script column is marked "[backfilled]".

Status ok only if the summary records a run commit that exists in this repository AND the MD list row counts equal the
summary's n_md_global / n_md_strict; otherwise status ``backfilled_no_commit`` or ``backfilled_mismatch`` (not
citable: check_report_numbers.py accepts status ok only).
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import re
import subprocess
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd  # noqa: E402

from biomark import paths as P  # noqa: E402

COLS = ["run_id", "utc", "script", "commit", "dirty", "config_sha256", "design_params_sha256", "data_status_sha256",
        "status", "outputs"]
UNIVERSES = ("eligible", "ctrl")


def _run_arms():
    spec = importlib.util.spec_from_file_location("run_arms", Path(__file__).resolve().parent / "run_arms.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def commit_exists(commit: str, root: Path) -> bool:
    if not re.fullmatch(r"[0-9a-f]{7,40}", str(commit or "")):
        return False
    r = subprocess.run(["git", "-C", str(root), "cat-file", "-e", f"{commit}^{{commit}}"], capture_output=True)
    return r.returncode == 0


def _rows_csv(p: Path) -> int:
    with open(p, encoding="utf-8") as fh:
        return max(sum(1 for _ in fh) - 1, 0)


def backfill_row(universe: str, root: Path, outputs: list) -> dict | None:
    """The registry row for the existing full run of ``universe`` (None if it never ran)."""
    sfx = "" if universe == "eligible" else f"_{universe}"
    summ_f = root / "results" / f"benchmark{sfx}" / "summary.json"
    if not summ_f.exists():
        return None
    s = json.loads(summ_f.read_text(encoding="utf-8"))
    tdir = root / "results" / ("timing" if universe == "eligible" else f"benchmark{sfx}/timing")
    t = json.loads((tdir / "md_strict.json").read_text(encoding="utf-8")) if (tdir / "md_strict.json").exists() else {}
    utc = str(t.get("utc") or "")
    stamp = utc.replace(":", "").replace("-", "").split("+")[0] if utc else "00000000T000000"
    h = s.get("config_sha256") if isinstance(s.get("config_sha256"), dict) else {}
    commit = str(s.get("commit") or "")
    status = "ok"
    if not commit_exists(commit, root):
        status = "backfilled_no_commit"
    else:
        md = root / "results" / f"md{sfx}"
        for f, key in (("md_global.csv", "n_md_global"), ("md_strict.csv", "n_md_strict")):
            if key in s and (not (md / f).exists() or _rows_csv(md / f) != int(s[key])):
                status = "backfilled_mismatch"
    script = "scripts/run_arms.py" + ("" if universe == "eligible" else f" --universe {universe}")
    return {"run_id": f"{stamp}_run_arms{sfx}_backfill", "utc": utc,
            "script": f"{script} [backfilled from {summ_f.relative_to(root).as_posix()}; timing run_id "
                      f"{t.get('run_id', '')}]",
            "commit": commit, "dirty": "",
            "config_sha256": hashlib.sha256("".join(f"{k}={h[k]}" for k in sorted(h)).encode()).hexdigest() if h else "",
            "design_params_sha256": h.get("design_params.json", ""), "data_status_sha256": "",
            "status": status, "outputs": ";".join(outputs)}


def main(argv=None, root: Path | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)
    root = Path(root or P.PROJECT)
    ra = _run_arms()
    reg = root / "results" / "registry.csv"
    old = pd.read_csv(reg, dtype=str) if reg.exists() and reg.stat().st_size else pd.DataFrame(columns=COLS)
    new = []
    for u in UNIVERSES:
        outs = [o for o in ra.registry_outputs(u, existing_only=False) if (root / o).exists()]
        row = backfill_row(u, root, outs) if outs else None
        if row is None:
            print(f"{u}: no full run found")
            continue
        if row["run_id"] in set(old["run_id"].astype(str)):
            print(f"{u}: {row['run_id']} already in the registry")
            continue
        print(f"{u}: {row['run_id']} status {row['status']} commit {row['commit'][:10]} ({len(outs)} outputs)")
        new.append(row)
    if new and not a.dry_run:
        from filelock import FileLock
        with FileLock(str(reg) + ".lock", timeout=120):
            cur = pd.read_csv(reg, dtype=str) if reg.exists() and reg.stat().st_size else pd.DataFrame(columns=COLS)
            pd.concat([cur, pd.DataFrame(new, columns=COLS)], ignore_index=True).to_csv(reg, index=False)
        print(f"appended {len(new)} row(s) to {reg.relative_to(root).as_posix()} "
              f"({dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
