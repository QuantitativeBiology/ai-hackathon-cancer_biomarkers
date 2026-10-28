"""Pending-download status (DESIGN 1.5): which files are in external/depmap_extra and which fallbacks are active.

File names, tiers and fallback texts come from ``configs/build_params.json: data_status``; recorded checksums
from the file named there (``configs/depmap_extra_sha256.txt``).
"""
from __future__ import annotations

import datetime as _dt
import hashlib
from pathlib import Path

from . import paths as P


def sha256(path, chunk: int = 1 << 22) -> str:  # literal-ok: read buffer size
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def _cfg() -> dict:
    return P.load_config("build_params")["data_status"]


def extra_dir() -> Path:
    return P.PROJECT / _cfg()["dir"]


def recorded_checksums() -> dict[str, str]:
    p = P.PROJECT / _cfg()["checksums_file"]
    out = {}
    if p.exists():
        for ln in p.read_text(encoding="utf-8").splitlines():
            parts = ln.split()
            if len(parts) >= 2:
                out[Path(parts[1]).name] = parts[0].lower()
    return out


def find(name: str) -> Path | None:
    """Path of a pending download in external/depmap_extra (by file name or the entry's glob), else None."""
    d = extra_dir()
    for e in _cfg()["files"]:
        if e["file"] == name:
            if "glob" in e:
                hits = sorted(d.glob(e["glob"])) if d.exists() else []
                return hits[0] if hits else None
            p = d / e["file"]
            return p if p.exists() else None
    p = d / name
    return p if p.exists() else None


def present(name: str) -> bool:
    return find(name) is not None


def status(compute_sha: bool = True) -> dict:
    """Status dict written to results/data_status.json."""
    rec = recorded_checksums()
    rows, fallbacks = [], []
    for e in _cfg()["files"]:
        p = find(e["file"])
        row = {"file": e["file"], "tier": e["tier"], "present": p is not None, "path": None, "size": None,
               "sha256": None, "sha256_recorded": None, "sha256_matches": None, "used_for": e["used_for"],
               "fallback": None}
        if p is not None:
            row["path"] = p.relative_to(P.PROJECT).as_posix()
            row["size"] = p.stat().st_size
            if compute_sha:
                row["sha256"] = sha256(p)
            r = rec.get(p.name)
            row["sha256_recorded"] = r
            if r is not None and row["sha256"] is not None:
                row["sha256_matches"] = r == row["sha256"]
        else:
            row["fallback"] = e["fallback"]
            if e["tier"] not in ("dropped", "not_requested_in_design"):
                fallbacks.append({"file": e["file"], "fallback": e["fallback"]})
        rows.append(row)
    listed = {r["path"] for r in rows if r["path"]}
    others = []
    d = extra_dir()
    for p in sorted(d.iterdir()) if d.exists() else []:
        rel = p.relative_to(P.PROJECT).as_posix()
        if p.is_file() and rel not in listed:
            sh = sha256(p) if compute_sha else None
            others.append({"path": rel, "size": p.stat().st_size, "sha256": sh,
                           "sha256_recorded": rec.get(p.name),
                           "sha256_matches": (rec.get(p.name) == sh) if rec.get(p.name) and sh else None})
    return {"generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
            "other_files": others,
            "dir": _cfg()["dir"], "files": rows, "fallbacks_active": fallbacks,
            "all_present_match_recorded": all(r["sha256_matches"] is not False for r in rows + others)}
