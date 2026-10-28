"""Read-only loaders for Part 1 inputs and outputs (docs/part1_brief.md).

Every loader returns pandas objects indexed by ``ModelID`` (``ACH-xxxxxx``) and/or Entrez ID strings
(``"7105"``). Nothing here writes into Part 1, and nothing here reads Chronos gene-effect values:
targets go through ``biomark.guard.load_effects`` and Part 1 test metrics through
``biomark.guard.p1_test_metrics``.

``install_p1_write_guard()`` (called by every loader and by ``import_essml``) adds a process-wide
audit hook that refuses any write-mode ``open``/mkdir/rename/remove under P1 and wraps
``pyarrow.parquet.write_table`` (pandas' parquet writer bypasses Python ``open``).
"""
from __future__ import annotations

import functools
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from . import paths as P

# ---------------------------------------------------------------------------------------------
# Write protection for Part 1
# ---------------------------------------------------------------------------------------------
_GUARD_INSTALLED = False
_P1_PREFIX = os.path.normcase(str(P.P1)).rstrip("\\/") + os.sep
_P1_ROOT_NC = os.path.normcase(str(P.P1)).rstrip("\\/")
_WRITE_FLAGS = os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_APPEND | os.O_TRUNC
_PATH_EVENTS = {  # event -> indices of path arguments that are write targets
    "os.mkdir": (0,), "os.rename": (0, 1), "os.remove": (0,), "os.rmdir": (0,), "os.truncate": (0,),
    "os.chmod": (0,), "os.utime": (0,), "os.symlink": (1,), "os.link": (1,), "shutil.rmtree": (0,),
    "shutil.copyfile": (1,), "shutil.move": (0, 1), "shutil.copytree": (1,),
}


class Part1WriteError(PermissionError):
    """Raised on any attempt to write under the read-only Part 1 tree."""


def _in_p1(path) -> bool:
    try:
        if path is None or isinstance(path, int):
            return False
        s = os.path.normcase(os.path.abspath(os.fsdecode(os.fspath(path))))
    except TypeError:
        return False
    return s == _P1_ROOT_NC or s.startswith(_P1_PREFIX)


def _audit(event, args):
    if event == "open":
        path, mode, flags = (tuple(args) + (None, None, None))[:3]
        writing = (isinstance(mode, str) and any(c in mode for c in "wax+")) or \
                  (isinstance(flags, int) and bool(flags & _WRITE_FLAGS))
        if writing and _in_p1(path):
            raise Part1WriteError(f"refused write into Part 1: {path}")
    elif event in _PATH_EVENTS:
        for i in _PATH_EVENTS[event]:
            if i < len(args) and _in_p1(args[i]):
                raise Part1WriteError(f"refused {event} under Part 1: {args[i]}")


def install_p1_write_guard() -> None:
    """Idempotently install the Part 1 write blocker for this process."""
    global _GUARD_INSTALLED
    if _GUARD_INSTALLED:
        return
    sys.addaudithook(_audit)
    orig = pq.write_table

    @functools.wraps(orig)
    def write_table(table, where, *a, **k):
        if _in_p1(where):
            raise Part1WriteError(f"refused parquet write into Part 1: {where}")
        return orig(table, where, *a, **k)

    pq.write_table = write_table
    _GUARD_INSTALLED = True


install_p1_write_guard()

# ---------------------------------------------------------------------------------------------
# Safe importer for Part 1's essml package
# ---------------------------------------------------------------------------------------------


def _refuse(name):
    def f(*a, **k):
        raise Part1WriteError(f"essml.{name} writes into Part 1 and is disabled in Part 2")
    f.__name__ = name.rsplit(".", 1)[-1]
    return f


def _rebind_defaults(module, old, new, only=("*",)) -> int:
    """Replace default argument values equal to ``old`` by ``new`` in the module's functions/methods
    (defaults such as ``model_dir=MODEL_DIR`` are bound at definition time)."""
    n = 0

    def fix(fn):
        nonlocal n
        if getattr(fn, "__module__", None) != module.__name__:
            return
        if "*" not in only and fn.__name__ not in only:
            return
        if fn.__defaults__:
            d = tuple(new if (isinstance(v, Path) and v == old) else v for v in fn.__defaults__)
            n += sum(a is not b for a, b in zip(d, fn.__defaults__))
            fn.__defaults__ = d
        if fn.__kwdefaults__:
            for k, v in fn.__kwdefaults__.items():
                if isinstance(v, Path) and v == old:
                    fn.__kwdefaults__[k] = new
                    n += 1

    for obj in list(vars(module).values()):
        if isinstance(obj, type):
            for m in vars(obj).values():
                f = getattr(m, "__func__", m)
                if callable(f) and hasattr(f, "__defaults__"):
                    fix(f)
        elif callable(obj) and hasattr(obj, "__defaults__"):
            fix(obj)
    return n


def import_essml(models: bool = True):
    """Import Part 1's ``essml`` read-only and return the package.

    - sets ``sys.dont_write_bytecode = True`` BEFORE adding ``P1/src`` to ``sys.path``;
    - installs the Part 1 write blocker;
    - redirects the write hazards of brief section 6 (``m2_local.INNER_CACHE``,
      ``m5_stack.INNER_CACHE/OUTER_CACHE``, ``modelstore.MODEL_DIR`` incl. the
      ``save_fold_model`` default, ``guard.ACCESS_LOG`` incl. defaults) to
      ``project_b/cache/p1_redirect/``; readers such as ``modelstore.load_fold_models`` keep
      reading P1;
    - disables ``io.build_cache``, ``io.write_gene_map`` and ``features.save_raw``.

    Reusable helpers then include ``essml.io.parse_gene_column``, ``essml.metrics``,
    ``essml.mutdeps``, ``essml.predictability`` and ``essml.features.FoldTransformer``.
    Use ``biomark.guard.load_effects`` for targets, not ``essml.guard``.
    """
    sys.dont_write_bytecode = True
    install_p1_write_guard()
    src = str(P.P1_SRC)
    if src not in sys.path:
        sys.path.insert(0, src)
    import importlib

    essml = importlib.import_module("essml")
    if not P.is_under(essml.__file__, P.P1_SRC):
        raise ImportError(f"'essml' resolved to {essml.__file__}, not Part 1's {P.P1_SRC}")
    mods = ["essml.io", "essml.guard", "essml.modelstore", "essml.features"]
    if models:
        mods += ["essml.models.m2_local", "essml.models.m5_stack"]
    for m in mods:
        importlib.import_module(m)
    P.ensure_dirs()
    for dotted, (new, rebind) in P.P1_WRITE_HAZARDS.items():
        modname, attr = dotted.rsplit(".", 1)
        mod = sys.modules.get(modname)
        if mod is None:
            continue
        old = getattr(mod, attr)
        if P.is_under(old, P.P1):  # first call: redirect; later calls are no-ops
            setattr(mod, attr, Path(new))
            if rebind:
                _rebind_defaults(mod, Path(old), Path(new), rebind)
    io_ = sys.modules["essml.io"]
    io_.build_cache = _refuse("io.build_cache")
    io_.write_gene_map = _refuse("io.write_gene_map")
    sys.modules["essml.features"].save_raw = _refuse("features.save_raw")
    return essml


# ---------------------------------------------------------------------------------------------
# Generic parquet reader with column / row selection
# ---------------------------------------------------------------------------------------------


def _schema_columns(path) -> list[str]:
    s = pq.read_schema(path)
    idx = set()
    if s.metadata and b"pandas" in s.metadata:
        import json
        meta = json.loads(s.metadata[b"pandas"])
        idx = {c for c in meta.get("index_columns", []) if isinstance(c, str)}
    return [c for c in s.names if c not in idx]


def read_parquet(path, columns=None, ids=None, strict: bool = True) -> pd.DataFrame:
    """Read a ModelID-indexed P1 parquet with optional column (pyarrow) and row (filter) selection.

    ``strict=True`` raises KeyError for requested columns/IDs that are absent; otherwise they are
    dropped silently. Rows come back in the order of ``ids`` when given.
    """
    if P.is_under(path, P.P1_TEST_RESULTS):
        raise PermissionError("Part 1 test results are read only via biomark.guard.p1_test_metrics")
    if P.is_under(path, P.P1_CRISPR_PARQUET):
        raise PermissionError("Chronos effects are read only via biomark.guard.load_effects")
    if columns is not None:
        columns = [str(c) for c in columns]
        have = set(_schema_columns(path))
        miss = [c for c in columns if c not in have]
        if miss and strict:
            raise KeyError(f"{len(miss)} columns not in {Path(path).name}, e.g. {miss[:5]}")
        columns = [c for c in columns if c in have]
    filters = None
    if ids is not None:
        ids = [str(i) for i in ids]
        filters = [("ModelID", "in", ids)]
    df = pd.read_parquet(path, columns=columns, filters=filters)
    if ids is not None:
        miss = [i for i in ids if i not in df.index]
        if miss and strict:
            raise KeyError(f"{len(miss)} ModelIDs not in {Path(path).name}, e.g. {miss[:5]}")
        df = df.loc[[i for i in ids if i in df.index]]
    df.columns = df.columns.astype(str)
    return df


# ---------------------------------------------------------------------------------------------
# Splits and gene annotations
# ---------------------------------------------------------------------------------------------


@functools.lru_cache(maxsize=1)
def _splits_cached() -> pd.DataFrame:
    sp = pd.read_csv(P.P1_SPLITS, dtype={"ModelID": str, "PatientID": str, "role": str})
    return sp.set_index("ModelID")


def splits() -> pd.DataFrame:
    """1292 x [PatientID, OncotreeLineage, role, cv_fold]; role in dev (1022) / test (255) /
    noncancerous (15); cv_fold 0-4 for dev, -1 otherwise."""
    return _splits_cached().copy()


def ids_by_role(role: str) -> pd.Index:
    sp = _splits_cached()
    return sp.index[sp["role"] == role]


def dev_ids() -> pd.Index:
    return ids_by_role("dev")


def test_ids() -> pd.Index:
    """Replication (R) line IDs. IDs only: no effect values."""
    return ids_by_role("test")


def noncancerous_ids() -> pd.Index:
    return ids_by_role("noncancerous")


def roles_of(ids) -> pd.Series:
    """Role per ModelID; 'unassigned' for lines not in splits.csv."""
    ids = pd.Index([str(i) for i in ids], name="ModelID")
    return _splits_cached()["role"].reindex(ids).fillna("unassigned").rename("role")


def folds() -> pd.Series:
    """Dev CV fold (0-4) per dev ModelID."""
    sp = _splits_cached()
    return sp.loc[sp["role"] == "dev", "cv_fold"].astype(int).rename("cv_fold")


def gene_map() -> pd.DataFrame:
    """21,269 rows indexed by Entrez string: [symbol, in_matrices]."""
    install_p1_write_guard()
    gm = pd.read_parquet(P.P1_GENE_MAP)
    gm["entrez"] = gm["entrez"].astype(str)
    return gm.set_index("entrez")


def symbols(entrez) -> pd.Series:
    """Entrez -> symbol (NaN where unknown)."""
    gm = gene_map()["symbol"]
    return gm.reindex(pd.Index([str(e) for e in entrez], name="entrez"))


def entrez_of(syms) -> pd.Series:
    """Symbol -> Entrez (first match in gene_map; NaN where unknown)."""
    gm = gene_map().reset_index()
    first = gm.drop_duplicates("symbol").set_index("symbol")["entrez"]
    return first.reindex(pd.Index(list(syms), name="symbol"))


def gene_classes() -> pd.DataFrame:
    """P1 configs/gene_classes.csv indexed by Entrez: symbol, all, variable, selective,
    rare_selective, common_essential, nonessential, n_obs, sd, n_below, frac_below."""
    return pd.read_csv(P.P1_GENE_CLASSES, dtype={"entrez": str}).set_index("entrez")


def gene_tiers(model: str = "m5") -> pd.DataFrame:
    """P1 predictable_genes{,_m1}.csv indexed by Entrez (tier in high/medium/low/not evaluable)."""
    return pd.read_csv(P.P1_PREDICTABLE_GENES[model], dtype={"entrez": str}).set_index("entrez")


def eval_genes() -> pd.DataFrame:
    return pd.read_csv(P.P1_EVAL_GENES, dtype={"entrez": str})


def library_groups() -> pd.DataFrame:
    """P1 NaN-pattern library groups (dev lines only): [lib_group, n_nan_genes]."""
    return pd.read_csv(P.P1_LIBRARY_GROUPS, dtype={"ModelID": str}).set_index("ModelID")


# ---------------------------------------------------------------------------------------------
# Feature blocks and raw caches
# ---------------------------------------------------------------------------------------------


def feature_block(name: str, columns=None, ids=None, strict: bool = True) -> pd.DataFrame:
    """P1 raw feature block ``cache/features/<name>.parquet`` (E, C, C_instab, M, M_burden, H,
    F_gene, F_pair, meta, meta_prov, composite, P), ModelID index, string columns."""
    if name not in P.P1_FEATURE_BLOCKS:
        raise KeyError(f"unknown feature block {name!r}; expected one of {P.P1_FEATURE_BLOCKS}")
    return read_parquet(P.P1_FEATURES / f"{name}.parquet", columns=columns, ids=ids, strict=strict)


def raw_matrix(name: str, columns=None, ids=None, strict: bool = True) -> pd.DataFrame:
    """P1 raw cache: expression (2130x19218), cn (1763x18778), damaging (2395x19657),
    hotspot (2395x165) float32; fusion = long table (94,718 rows). CRISPR is refused here."""
    if name == "crispr":
        raise PermissionError("Chronos effects are read only via biomark.guard.load_effects")
    if name not in P.P1_RAW_CACHES:
        raise KeyError(f"unknown raw cache {name!r}")
    if name == "fusion":
        df = pd.read_parquet(P.P1_RAW_CACHES[name], columns=columns)
        return df if ids is None else df[df["ModelID"].isin([str(i) for i in ids])]
    return read_parquet(P.P1_RAW_CACHES[name], columns=columns, ids=ids, strict=strict)


def model_table(columns=None) -> pd.DataFrame:
    """Model.csv (2554 lines) indexed by ModelID, optionally a column subset. Only empty cells are NaN
    (literal strings such as "None" are kept)."""
    install_p1_write_guard()
    use = None if columns is None else ["ModelID"] + [c for c in columns if c != "ModelID"]
    return pd.read_csv(P.P1_MODEL_CSV, usecols=use, low_memory=False, keep_default_na=False,
                       na_values=[""]).set_index("ModelID")


def model_condition(columns=None) -> pd.DataFrame:
    """ModelCondition.csv (3385 conditions) indexed by ModelConditionID, optionally a column subset.
    Only empty cells are NaN (PlateCoating "None" is a real level)."""
    install_p1_write_guard()
    use = None if columns is None else ["ModelConditionID"] + [c for c in columns if c != "ModelConditionID"]
    return pd.read_csv(P.P1_MODEL_CONDITION_CSV, usecols=use, low_memory=False, keep_default_na=False,
                       na_values=[""]).set_index("ModelConditionID")


# ---------------------------------------------------------------------------------------------
# Predictions
# ---------------------------------------------------------------------------------------------


def run_id(model: str) -> str:
    return P.P1_RUN_IDS.get(model, model)


def oof_predictions(model: str = "m5", columns=None, ids=None) -> pd.DataFrame:
    """Dev out-of-fold predictions (1022 x 18434 float32, no NaN). ``model`` is a key of
    ``paths.P1_RUN_IDS`` (m5 = headline, m1, m2n, m3n, t) or a full run id."""
    return read_parquet(P.P1_RUNS / run_id(model) / "pred.parquet", columns=columns, ids=ids)


def component_oof(columns=None) -> dict[str, pd.DataFrame]:
    """OOF predictions of the components / variants: {m1, m2n, m3n, t}."""
    return {k: oof_predictions(k, columns=columns) for k in ("m1", "m2n", "m3n", "t")}


def per_gene_fold(model: str = "m5") -> pd.DataFrame:
    """Per (fold, gene) dev-CV metrics of a run (``cache/runs/<run>/per_gene_fold.parquet``)."""
    df = pd.read_parquet(P.P1_RUNS / run_id(model) / "per_gene_fold.parquet")
    if "gene" in df:
        df["gene"] = df["gene"].astype(str)
    return df


def test_predictions(model: str = "m5_avg", columns=None, ids=None) -> pd.DataFrame:
    """Replication-line predictions from all-dev refits (255 x 18434), ``cache/test/preds/<model>``.
    Predictions only: no targets, no metrics."""
    return read_parquet(P.P1_TEST_PREDS / f"{model}.parquet", columns=columns, ids=ids)


def fold_centre(pred: pd.DataFrame, fold_of: pd.Series | None = None) -> pd.DataFrame:
    """Centre each gene within each CV fold (Part 1 D3 pooled-OOF artefact fix).

    ``fold_of``: fold label per ModelID (default: dev ``folds()``); every row of ``pred`` must have
    one. NaN-aware (nanmean). Returns a frame of the same shape, index, columns and dtype.
    """
    f = (folds() if fold_of is None else fold_of).reindex(pred.index)
    if f.isna().any():
        raise KeyError(f"{int(f.isna().sum())} rows of pred have no fold, e.g. {list(f.index[f.isna()][:3])}")
    X = pred.to_numpy(copy=True)
    fv = f.to_numpy()
    for k in np.unique(fv):
        m = fv == k
        with np.errstate(all="ignore"):
            mu = np.nanmean(X[m], axis=0)
        X[m] -= mu.astype(X.dtype, copy=False)
    return pd.DataFrame(X, index=pred.index, columns=pred.columns)


# ---------------------------------------------------------------------------------------------
# Screen and culture covariates (no gene-effect values are read)
# ---------------------------------------------------------------------------------------------

_COND_COLS = ["CellFormat", "GrowthMedia", "SerumFreeMedia", "PlateCoating"]
_MODEL_COLS = ["GrowthPattern", "ModelType", "OnboardedMedia", "FormulationID", "SerumFreeMedia"]


def _join(values) -> object:
    v = sorted({str(x) for x in values if pd.notna(x) and str(x) != ""})
    return "+".join(v) if v else np.nan


def screen_covariates(extra_condition_cols=("FormulationID",)) -> pd.DataFrame:
    """Screen/culture covariates for every CRISPR-screened line (dev, test, noncancerous, unassigned).

    ``extra_condition_cols``: further ModelCondition.csv columns of the screened condition(s), '+'-joined
    like the default ones (default adds FormulationID, needed by the base-medium parser of DESIGN 5.3).

    Sources: CRISPRScreenMap (ScreenID -> ModelID) -> ScreenSequenceMap rows of those screens
    (Library, Days, ScreenType, ModelConditionID) -> ModelCondition.csv; plus Model.csv.
    Columns:
      role, n_screens, Library (sorted '+'-joined if several), ScreenDays (median Days over the
      sequences of the line's screens), ScreenType, ModelConditionID, CellFormat, GrowthMedia,
      SerumFreeMedia, PlateCoating (condition-level, '+'-joined over the screened conditions),
      model_GrowthPattern, model_ModelType, model_OnboardedMedia, model_FormulationID,
      model_SerumFreeMedia (Model.csv).
    """
    install_p1_write_guard()
    sm = pd.read_csv(P.P1_SCREEN_MAP_CSV, dtype=str)
    ss = pd.read_csv(P.P1_SCREEN_SEQUENCE_MAP_CSV,
                     usecols=["ScreenID", "Days", "Library", "ScreenType", "ModelConditionID"],
                     dtype={"ScreenID": str, "Library": str, "ScreenType": str, "ModelConditionID": str})
    ss = ss.merge(sm, on="ScreenID", how="inner")
    cond_cols = list(dict.fromkeys(list(_COND_COLS) + list(extra_condition_cols)))
    mc = model_condition(cond_cols)
    ss = ss.join(mc, on="ModelConditionID")
    g = ss.groupby("ModelID")
    out = pd.DataFrame(index=pd.Index(sorted(sm["ModelID"].unique()), name="ModelID"))
    out["role"] = roles_of(out.index).to_numpy()
    out["n_screens"] = sm.groupby("ModelID")["ScreenID"].nunique().reindex(out.index).astype(int)
    out["Library"] = g["Library"].agg(_join).reindex(out.index)
    out["ScreenDays"] = g["Days"].median().reindex(out.index)
    for c in ["ScreenType", "ModelConditionID"] + cond_cols:
        out[c] = g[c].agg(_join).reindex(out.index)
    mt = model_table(_MODEL_COLS).reindex(out.index)
    for c in _MODEL_COLS:
        out[f"model_{c}"] = mt[c]
    return out
