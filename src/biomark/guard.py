"""Part 2 replication guard (DESIGN.md section 1).

Chronos gene effects (DepMap 26Q3) are read directly from ``P1/cache/crispr.parquet`` with a pyarrow
row filter, so rows that are refused never enter memory. Part 1's own guarded loader is NOT used (its
per-finalist test budget is spent, and it writes its log into Part 1).

Rules (roles from ``P1/configs/splits.csv``):

- **dev** (discovery, 1022 lines): always allowed, never logged.
- **test** (replication R, 255 lines; scope ``"replication"``): allowed only if ALL of
  (a) git tag ``freeze-catalogue`` exists in project_b;
  (b) env ``PART2_UNBLIND == "1"``;
  (c) the project_b working tree is clean (no modified/staged tracked file other than the access
      log; no untracked file outside ``results/``);
  (d) the caller passes a non-empty ``purpose`` string.
- **noncancerous** (15 lines; scope ``"limits"``): refused unless ``purpose`` starts with
  ``"limits:"`` and the tag exists.
- **unassigned** screened lines (not in splits.csv): always refused.

Every granted non-dev access appends one row (UTC time, commit, purpose, scope, n ids, roles, n genes,
caller) to ``results/replication_access_log.csv`` under a file lock. Refusals raise
``ReplicationGuardError`` and write nothing.

``p1_test_metrics(relpath, purpose)`` is the only reader for ``P1/results/test/`` and applies the
test rule (a)-(d).
"""
from __future__ import annotations

import csv
import datetime as _dt
import inspect
import json
import os
from pathlib import Path

import pandas as pd

from . import gitutil
from . import p1
from . import paths as P

CRISPR_PARQUET = P.P1_CRISPR_PARQUET
ACCESS_LOG = P.REPLICATION_ACCESS_LOG
LOCK_FILE = P.CACHE / "replication_access_log.lock"
REPO = P.PROJECT
FREEZE_TAG = P.FREEZE_TAG
UNBLIND_ENV = "PART2_UNBLIND"
LOG_COLS = ["utc_time", "commit", "purpose", "scope", "n_ids", "roles", "n_genes", "caller"]
SCOPES = {"dev": {"dev"}, "replication": {"dev", "test"}, "limits": {"dev", "noncancerous"}}
SCOPE_DEFAULT_ROLE = {"dev": "dev", "replication": "test", "limits": "noncancerous"}


class ReplicationGuardError(PermissionError):
    """Raised when replication / non-cancerous effects or P1 test metrics are requested without
    meeting the guard conditions."""


# ---------------------------------------------------------------------------------------------
# Conditions
# ---------------------------------------------------------------------------------------------


def tag_exists() -> bool:
    return gitutil.tag_exists(REPO, FREEZE_TAG)


def _rel(path) -> str:
    try:
        return Path(path).resolve().relative_to(Path(REPO).resolve()).as_posix()
    except ValueError:
        return Path(path).as_posix()


def dirty_entries() -> list[str]:
    """Porcelain entries that make the tree unclean for unblinding."""
    log_rel = _rel(ACCESS_LOG)
    bad = []
    for ln in gitutil.status_porcelain(REPO, untracked=True):
        code, path = ln[:2], ln[3:].strip().strip('"')  # literal-ok: git porcelain status columns
        if " -> " in path:
            path = path.split(" -> ", 1)[1]
        if path == log_rel:
            continue
        if code == "??" and (path == "results/" or path.startswith("results/")):
            continue
        bad.append(ln)
    return bad


def tree_clean() -> bool:
    return not dirty_entries()


def unblind_check(purpose=None) -> list[str]:
    """Reasons why replication access would be refused now (empty list = allowed)."""
    reasons = []
    if not tag_exists():
        reasons.append(f"git tag '{FREEZE_TAG}' does not exist in project_b")
    if os.environ.get(UNBLIND_ENV) != "1":
        reasons.append(f"env {UNBLIND_ENV} != '1'")
    if not tree_clean():
        reasons.append(f"working tree not clean: {dirty_entries()[:5]}")  # literal-ok: error-message preview length
    if not (isinstance(purpose, str) and purpose.strip()):
        reasons.append("no purpose string given")
    return reasons


def limits_check(purpose=None) -> list[str]:
    reasons = []
    if not (isinstance(purpose, str) and purpose.startswith("limits:") and purpose[len("limits:"):].strip()):
        reasons.append("non-cancerous lines need a purpose starting with 'limits:'")
    if not tag_exists():
        reasons.append(f"git tag '{FREEZE_TAG}' does not exist in project_b")
    return reasons


# ---------------------------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------------------------


def _caller() -> str:
    here = Path(__file__).resolve()
    for fr in inspect.stack()[1:]:
        if Path(fr.filename).resolve() != here:
            return f"{Path(fr.filename).name}:{fr.function}:{fr.lineno}"
    return "?"


def _commit() -> str:
    try:
        return gitutil.head(REPO)
    except Exception:  # noqa: BLE001 - repo without commits
        return ""


def _log(purpose, scope, roles: pd.Series, n_genes) -> None:
    from filelock import FileLock

    log = Path(ACCESS_LOG)
    log.parent.mkdir(parents=True, exist_ok=True)
    Path(LOCK_FILE).parent.mkdir(parents=True, exist_ok=True)
    vc = roles.value_counts()
    row = {"utc_time": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
           "commit": _commit(), "purpose": purpose, "scope": scope, "n_ids": int(len(roles)),
           "roles": ";".join(f"{k}:{int(v)}" for k, v in sorted(vc.items())), "n_genes": n_genes,
           "caller": _caller()}
    with FileLock(str(LOCK_FILE), timeout=120):  # literal-ok: file-lock wait seconds
        new = not log.exists() or log.stat().st_size == 0
        with open(log, "a", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=LOG_COLS)
            if new:
                w.writeheader()
            w.writerow(row)


def read_log() -> pd.DataFrame:
    p = Path(ACCESS_LOG)
    return pd.read_csv(p) if p.exists() else pd.DataFrame(columns=LOG_COLS)


# ---------------------------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------------------------


def load_effects(ids=None, scope: str = "dev", purpose: str | None = None, genes=None) -> pd.DataFrame:
    """Chronos gene effects (float32 DataFrame; rows = ``ids`` in the given order, index ModelID;
    columns = Entrez strings, all 18,434 or ``genes``).

    ``scope``: "dev" (dev IDs only), "replication" (test, optionally with dev IDs) or "limits"
    (non-cancerous, optionally with dev IDs). ``ids=None`` means all lines of the scope's role.
    """
    if scope not in SCOPES:
        raise ValueError(f"scope must be one of {sorted(SCOPES)}")
    if ids is None:
        ids = p1.ids_by_role(SCOPE_DEFAULT_ROLE[scope])
    ids = [str(i) for i in ids]
    if not ids:
        raise ValueError("no ids requested")
    roles = p1.roles_of(ids)
    present = set(roles)
    if "unassigned" in present:
        raise ReplicationGuardError(f"{int((roles == 'unassigned').sum())} IDs are not in the Part 1 population "
                                    "(splits.csv); refused")
    extra = present - SCOPES[scope]
    if extra:
        raise ReplicationGuardError(f"scope {scope!r} does not cover roles {sorted(extra)}")
    if "test" in present:
        why = unblind_check(purpose)
        if why:
            raise ReplicationGuardError("replication effects refused: " + "; ".join(why))
    if "noncancerous" in present:
        why = limits_check(purpose)
        if why:
            raise ReplicationGuardError("non-cancerous effects refused: " + "; ".join(why))
    cols = None if genes is None else [str(g) for g in genes]
    if present != {"dev"}:
        _log(purpose, scope, roles, "all" if cols is None else len(cols))
    df = pd.read_parquet(CRISPR_PARQUET, columns=cols, filters=[("ModelID", "in", sorted(set(ids)))])
    miss = sorted(set(ids) - set(df.index))
    if miss:
        raise KeyError(f"{len(miss)} ModelIDs have no CRISPR screen, e.g. {miss[:3]}")  # literal-ok: error-message preview length
    df = df.loc[ids]
    df.columns = df.columns.astype(str)
    return df


def screened_ids() -> pd.Index:
    """ModelIDs of the CRISPR matrix (index only; no values)."""
    return pd.Index(pd.read_parquet(CRISPR_PARQUET, columns=[]).index.astype(str), name="ModelID")


def p1_test_metrics(relpath: str = "", purpose: str | None = None):
    """Read a file (or list a directory) under ``P1/results/test/`` - refused unless the replication
    conditions hold (tag, PART2_UNBLIND=1, clean tree, purpose). Granted reads are logged.

    Returns a DataFrame (.csv/.csv.gz/.parquet), a parsed object (.json), text (other files) or a
    list of relative paths (directory)."""
    path = (P.P1_TEST_RESULTS / relpath).resolve()
    if not P.is_under(path, P.P1_TEST_RESULTS):
        raise ValueError("relpath must stay inside P1/results/test")
    why = unblind_check(purpose)
    if why:
        raise ReplicationGuardError("Part 1 test metrics refused: " + "; ".join(why))
    _log(purpose, f"p1_test_metrics:{Path(relpath).as_posix()}", pd.Series([], dtype=str), "")
    if path.is_dir():
        return sorted(p.relative_to(P.P1_TEST_RESULTS).as_posix() for p in path.rglob("*") if p.is_file())
    name = path.name.lower()
    if name.endswith(".csv") or name.endswith(".csv.gz"):
        return pd.read_csv(path)
    if name.endswith(".parquet"):
        return pd.read_parquet(path)
    if name.endswith(".json"):
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    return path.read_text(encoding="utf-8")
