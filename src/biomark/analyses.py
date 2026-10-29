"""Remaining Tier A analysis pieces: XLIB split and discovery-only eligibility (DESIGN v3.2 §6.3), component readout
(§4.2), design VIF / condition number (§5.3), Q4c population precision and carriers-only OR (§10.5), the OncoTree
parent-code mapping built from local files, and the registry / timing writers (§13.6, §14). Engine semantics are
unchanged: these functions call stats / en / arms / audit / endpoints as they are."""
from __future__ import annotations

import copy
import datetime as dt
import json
import re

import numpy as np
import pandas as pd
from scipy import stats as sps

from . import params as PR
from . import paths as P

# ------------------------------------------------------------------------------------------------
# Registry and timing
# ------------------------------------------------------------------------------------------------

REGISTRY_COLS = ["run_id", "utc", "script", "commit", "dirty", "config_sha256", "design_params_sha256",
                 "data_status_sha256", "status", "outputs"]


def _sha(path) -> str:
    import hashlib
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else ""


def append_registry(script: str, status: str, outputs, registry=None) -> dict:
    """Append one row to results/registry.csv under a file lock (DESIGN §13.6)."""
    from filelock import FileLock

    from . import gitutil
    reg = registry or (P.RESULTS / "registry.csv")
    reg.parent.mkdir(parents=True, exist_ok=True)
    try:
        commit = gitutil.head(P.PROJECT)
        dirty = bool(gitutil.status_porcelain(P.PROJECT, untracked=False))
    except Exception:  # noqa: BLE001
        commit, dirty = "", True
    h = PR.config_hashes()
    import hashlib
    utc = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    stamp = utc.replace(":", "").replace("-", "").split("+")[0]
    row = {"run_id": f"{stamp}_{script.rsplit('/', 1)[-1].replace('.py', '').replace(' ', '')}",
           "utc": utc, "script": script, "commit": commit, "dirty": dirty,
           "config_sha256": hashlib.sha256("".join(f"{k}={h[k]}" for k in sorted(h)).encode()).hexdigest(),
           "design_params_sha256": h["design_params.json"],
           "data_status_sha256": _sha(P.RESULTS / "data_status.json"), "status": status,
           "outputs": ";".join(str(o) for o in outputs)}
    with FileLock(str(reg) + ".lock", timeout=120):  # literal-ok: lock wait seconds
        old = pd.read_csv(reg) if reg.exists() and reg.stat().st_size else pd.DataFrame(columns=REGISTRY_COLS)
        pd.concat([old, pd.DataFrame([row])], ignore_index=True).to_csv(reg, index=False)
    return row


_RUN_ID = None


def run_id() -> str:
    """One id per process, shared by every timing stem the process writes (budget_check counts the process once)."""
    global _RUN_ID
    if _RUN_ID is None:
        import os
        import sys
        script = os.path.basename(sys.argv[0] or "python").replace(".py", "")
        _RUN_ID = f"{script}-{os.getpid()}-{dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%S')}"
    return _RUN_ID


def utc_now() -> dt.datetime:
    return dt.datetime.now(dt.timezone.utc)


def _check_workers(stem: str, wk) -> int:
    """workers must be an int >= 1 (no silent default; docs/process_checks.md §1.1)."""
    if isinstance(wk, bool) or wk is None:
        raise ValueError(f"timing {stem}: 'workers' (int >= 1) required")
    if isinstance(wk, float) and wk.is_integer():
        wk = int(wk)
    if not isinstance(wk, (int, np.integer)) or int(wk) < 1:
        raise ValueError(f"timing {stem}: invalid workers {wk!r} (int >= 1 required)")
    return int(wk)


def timing_subrecord(name: str, start_utc: dt.datetime, end_utc: dt.datetime, workers: int, **extra) -> dict:
    """A sub-record of a timing file (schema of a record minus job): budget_check counts its interval toward the
    parent stage's wall and toward concurrency."""
    return {"start_utc": start_utc.isoformat(timespec="seconds"), "utc": end_utc.isoformat(timespec="seconds"),
            "wall_s": float((end_utc - start_utc).total_seconds()), "workers": _check_workers(name, workers),
            "run_id": run_id(), **extra}


def guard_workers(requested: int | None, kind: str, log=None) -> tuple:
    """Pre-spawn memory guard: effective workers = min(requested, floor(room / per-worker RSS)), room =
    min(total budget - current parent RSS, available memory - safety margin), at least 1. Settings and the per-worker
    RSS estimates are in configs/runtime_params.json. Returns (effective, info) for the timing record."""
    import psutil
    rt = PR.runtime()
    m = rt["memory"]
    req = int(rt["workers_default"].get(kind, 1) if requested is None else requested)
    per = float(m["rss_per_worker_gb"].get(kind, m["rss_per_worker_gb"]["default"]))
    gib = float(1 << 30)  # literal-ok: bytes per GiB
    avail = psutil.virtual_memory().available / gib
    parent = psutil.Process().memory_info().rss / gib
    room = min(float(m["total_budget_gb"]) - parent, avail - float(m["safety_margin_gb"]))
    eff = max(1, min(req, int(room // per))) if req > 1 else 1
    info = {"kind": kind, "requested": req, "effective": eff, "per_worker_gb": per, "available_gb": round(avail, 2),
            "parent_rss_gb": round(parent, 2), "budget_gb": float(m["total_budget_gb"])}
    if log and eff < req:
        log(f"memory guard: {kind} workers {req} -> {eff} ({info})")
    return eff, info


def write_timing(stem: str, rec: dict, timing_dir=None, *, wall_s: float | None = None, workers: int | None = None,
                 start_utc: dt.datetime | None = None, end_utc: dt.datetime | None = None) -> dict:
    """results/timing/<stem>.json with the minimal schema of docs/process_checks.md §1.1: job (= stem), utc (end time),
    numeric wall_s, integer workers >= 1, start_utc and run_id. ``wall_s`` / ``workers`` may also come in ``rec``."""
    d = timing_dir or (P.PROJECT / PR.pp()["timing_dir"])
    d.mkdir(parents=True, exist_ok=True)
    rec = dict(rec)
    wall = float(rec.pop("wall_s")) if wall_s is None and "wall_s" in rec else wall_s
    wk = rec.pop("workers", None) if workers is None else workers
    wk = _check_workers(stem, wk)
    end = end_utc or utc_now()
    if wall is None and start_utc is not None:
        wall = (end - start_utc).total_seconds()
    if wall is None or not isinstance(wall, (int, float)):
        raise ValueError(f"timing {stem}: numeric wall_s required")
    start = start_utc or (end - dt.timedelta(seconds=float(wall)))
    out = {"job": stem, "utc": end.isoformat(timespec="seconds"), "start_utc": start.isoformat(timespec="seconds"),
           "wall_s": float(wall), "workers": wk, "run_id": rec.pop("run_id", run_id()), **rec}
    (d / f"{stem}.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    return out


# ------------------------------------------------------------------------------------------------
# XLIB (§6.3)
# ------------------------------------------------------------------------------------------------


def _quoted(text: str) -> str:
    m = re.search(r"'([^']+)'", text)
    if not m:
        raise ValueError(f"no quoted library string in {text!r}")
    return m.group(1)


def xlib_split(library: pd.Series) -> pd.Series:
    """'discovery' (library string exactly the discovery library), 'replication' (string without it), '' otherwise;
    strings from design_params.replication_s6.xlib_discovery / xlib_replication."""
    r6 = PR.dp()["replication_s6"]
    disc = _quoted(r6["xlib_discovery"])
    rep_excl = disc                     # xlib_replication: "D lines with no <discovery library> in the library string"
    if disc not in r6["xlib_replication"]:
        raise ValueError("xlib_replication does not name the discovery library")
    lib = library.fillna("").astype(str)
    has = lib.str.split("+").map(lambda parts: rep_excl in parts)
    return pd.Series(np.where(lib == disc, "discovery", np.where(~has & (lib != ""), "replication", "")),
                     index=library.index)


def discovery_eligibility(elig: pd.DataFrame, pred: np.ndarray, obs: np.ndarray, strata: np.ndarray, genes) -> pd.Series:
    """u_rule = discovery_only: r_within recomputed on discovery lines (no refit); class and tier criteria kept."""
    from . import eligibility as ELG
    e = PR.dp()["eligibility_s2"]
    min_total = int(P.load_config("build_params")["eligibility"]["min_lines_for_r"])
    r, _ = ELG.r_within_strata(pred, obs, strata, int(e["r_within_stratum_min_lines"]), min_total)
    el = elig.set_index("entrez").reindex(list(genes))
    ok = el["class_ok"].astype(bool) & el["tier_ok"].astype(bool) & (pd.Series(r, index=el.index) >= e["r_within_lin_lib_min"])
    return ok.fillna(False)


# ------------------------------------------------------------------------------------------------
# Component readout (§4.2)
# ------------------------------------------------------------------------------------------------

COMPONENTS = {"M1": "M1 expression ridge (E)",
              "M2n-M1": "M2 local inputs (own / paralog expression, CN, deep deletion, mutation)",
              "M3n-M1": "M3 pool (E, C, M, H, meta)"}


def component_targets(oof: dict, fold_of: pd.Series) -> dict:
    """Fold-centred ŷ_M1, ŷ_M2n − ŷ_M1, ŷ_M3n − ŷ_M1 (each frame lines x genes)."""
    from . import p1 as P1
    m1 = P1.fold_centre(oof["m1"], fold_of)
    return {"M1": m1, "M2n-M1": P1.fold_centre(oof["m2n"] - oof["m1"], fold_of),
            "M3n-M1": P1.fold_centre(oof["m3n"] - oof["m1"], fold_of)}


def data_with_pred(data, pred: np.ndarray):
    """Shallow copy of a pipeline.Data whose prediction target is ``pred`` (n x G, data.genes order)."""
    d = copy.copy(data)
    d.pred = np.asarray(pred, dtype=np.float64)
    d._obs = None
    return d


def component_selection(md: pd.DataFrame, comp_tabs: dict, access: pd.Series) -> pd.DataFrame:
    """Per MD-strict (g, d): full-D β and selection frequency per component; label expression-mediated (non-zero
    full-D β in M1 only), genomics-layer (non-zero in M2n−M1 or M3n−M1 only), both, or none; carrier_input names
    the component inputs carrying a proxy-class feature."""
    key = pd.MultiIndex.from_arrays([md["gene"].astype(str), md["feature"]])
    out = md[["gene", "feature"]].copy()
    out["gene"] = out["gene"].astype(str)
    nz = {}
    for c, t in comp_tabs.items():
        t = t.assign(gene=t["gene"].astype(str)).set_index(["gene", "feature"])
        out[f"beta_full_{c}"] = t["beta_full"].reindex(key).to_numpy()
        out[f"freq_{c}"] = t["selection_freq"].reindex(key).to_numpy()
        out[f"selected_{c}"] = t["selected"].reindex(key).fillna(False).astype(bool).to_numpy()
        nz[c] = np.nan_to_num(out[f"beta_full_{c}"].to_numpy()) != 0
    ex = nz.get("M1", np.zeros(len(out), bool))
    ge = nz.get("M2n-M1", np.zeros(len(out), bool)) | nz.get("M3n-M1", np.zeros(len(out), bool))
    out["layer"] = np.where(ex & ge, "both", np.where(ex, "expression-mediated", np.where(ge, "genomics-layer", "none")))
    acc = access.reindex(key).to_numpy()
    out["access_class"] = acc
    out["carrier_input"] = ["; ".join(COMPONENTS[c] for c in comp_tabs if nz[c][i]) if acc[i] == "proxy" else ""
                            for i in range(len(out))]
    return out


# ------------------------------------------------------------------------------------------------
# VIF / condition number of the D design (§5.3)
# ------------------------------------------------------------------------------------------------


def design_vif(cov) -> tuple:
    """VIF of every column of the shared D design (design blocks + droppable covariates, dummy-coded with the per-fit
    merge rule on all D rows) and the condition number of the standardised design, after removing exactly aliased
    columns (identical dummies; the estimating code is unaffected because it uses a rank-revealing basis).
    Returns (vif Series, condition number, list of (kept, aliased) column pairs)."""
    from . import covariates as CV
    rows = np.ones(cov.n, dtype=bool)
    cols, names = [], []
    for nm in cov.design_names + cov.droppable_names:
        B = cov.block(nm, rows)
        for j in range(B.shape[1]):
            cols.append(B[:, j])
            names.append(f"{nm}_{j}")
    X = pd.DataFrame(np.column_stack(cols), columns=names)
    aliases, keep = [], []
    for c in X.columns:
        twin = next((k for k in keep if np.array_equal(X[k].to_numpy(), X[c].to_numpy())), None)
        if twin is None:
            keep.append(c)
        else:
            aliases.append((twin, c))
    X = X[keep]
    v = CV.vif(X)
    Z = (X - X.mean()) / X.std(ddof=0).replace(0, np.nan)
    s = np.linalg.svd(Z.dropna(axis=1).to_numpy(), compute_uv=False)
    return v.sort_values(ascending=False), float(s[0] / s[-1]), aliases


# ------------------------------------------------------------------------------------------------
# OncoTree parent codes from local files; Q4c population precision and carriers-only OR
# ------------------------------------------------------------------------------------------------


def oncotree_terms(moa_json=None) -> dict:
    """OncoTree code -> set of terms, from the frozen MOAlmanac release (oncotree_code / oncotree_term pairs)."""
    p = moa_json or (P.EXTERNAL / "clinical" / "moalmanac" / "extracted" / "vanallenlab-moalmanac-db-6d3bd09"
                     / "molecular-oncology-almanac.json")
    pairs: dict = {}
    if not p.exists():
        return pairs

    def walk(o):
        if isinstance(o, dict):
            c, t = o.get("oncotree_code"), o.get("oncotree_term")
            if c and t:
                pairs.setdefault(str(c), set()).add(str(t))
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)

    walk(json.loads(p.read_text(encoding="utf-8")))
    return pairs


def indicated_lines(codes, model: pd.DataFrame, terms: dict) -> tuple:
    """Boolean per line (model rows): OncotreeCode in ``codes`` or Subtype / PrimaryDisease / Lineage name equal to an
    OncoTree term of one of the codes (partial parent mapping from Model.csv's own chain). Returns (mask, unmapped)."""
    codes = [c for c in codes if c]
    names = set().union(*[terms.get(c, set()) for c in codes]) if codes else set()
    m = model["OncotreeCode"].isin(codes)
    for col in ("OncotreeSubtype", "OncotreePrimaryDisease", "OncotreeLineage"):
        if col in model:
            m |= model[col].isin(names)
    hit = {c for c in codes if model["OncotreeCode"].eq(c).any() or
           any(model[col].isin(terms.get(c, set())).any() for col in ("OncotreeSubtype", "OncotreePrimaryDisease",
                                                                      "OncotreeLineage") if col in model)}
    return m.to_numpy(), sorted(set(codes) - hit)


def population_precision(x: np.ndarray, acc_dep: np.ndarray, indicated: np.ndarray) -> dict:
    """Among carriers: population precision = accurate-dependent fraction within the indication (and outside), and the
    carriers-only odds ratio of accurate-dependent membership, indicated vs not (conditional MLE, exact 95% CI)."""
    car = (np.asarray(x) == 1) & np.isfinite(np.asarray(acc_dep, dtype=float))
    a = np.asarray(acc_dep, dtype=float) == 1
    ind = np.asarray(indicated, dtype=bool)
    t = np.array([[int((car & ind & a).sum()), int((car & ind & ~a).sum())],
                  [int((car & ~ind & a).sum()), int((car & ~ind & ~a).sum())]])
    out = {"carriers_indicated": int(t[0].sum()), "carriers_not_indicated": int(t[1].sum()),
           "precision_indicated": t[0, 0] / t[0].sum() if t[0].sum() else np.nan,
           "precision_not_indicated": t[1, 0] / t[1].sum() if t[1].sum() else np.nan}
    if t[0].sum() and t[1].sum():
        r = sps.contingency.odds_ratio(t)
        ci = r.confidence_interval(confidence_level=float(PR.endpoints()["P1_primary"]["ci"]["level"]))
        out.update(or_carriers=float(r.statistic), or_ci_lo=float(ci.low), or_ci_hi=float(ci.high),
                   or_fisher_p=float(sps.fisher_exact(t)[1]))
    else:
        out.update(or_carriers=np.nan, or_ci_lo=np.nan, or_ci_hi=np.nan, or_fisher_p=np.nan)
    return out
