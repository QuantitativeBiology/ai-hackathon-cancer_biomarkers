"""Path constants for Part 2 (``PROJECT``) and Part 1 (``P1``, strictly read-only), plus config loading.

Part 1 paths follow docs/part1_brief.md. Nothing in this module writes into Part 1; ``ensure_dirs()``
creates Part 2 cache/result directories only.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

# ---------------------------------------------------------------------------------------------
# Part 2 (this repository)
# ---------------------------------------------------------------------------------------------
PROJECT = Path(__file__).resolve().parents[2]
SRC = PROJECT / "src"
CONFIGS = PROJECT / "configs"
SCRIPTS = PROJECT / "scripts"
DOCS = PROJECT / "docs"
EXTERNAL = PROJECT / "external"          # new downloads (git-ignored)

CACHE = PROJECT / "cache"                # git-ignored
CACHE_P1_REDIRECT = CACHE / "p1_redirect"  # destinations for Part 1 code's hard-coded write paths
M2_INNER_CACHE = CACHE_P1_REDIRECT / "m2_inner"
M5_INNER_CACHE = CACHE_P1_REDIRECT / "m5_inner"
M5_OUTER_CACHE = CACHE_P1_REDIRECT / "m5_outer"
P1_MODEL_DIR_REDIRECT = CACHE_P1_REDIRECT / "models"
P1_GUARD_LOG_REDIRECT = CACHE_P1_REDIRECT / "p1_guard_access_log.csv"

RESULTS = PROJECT / "results"
REPLICATION_ACCESS_LOG = RESULTS / "replication_access_log.csv"
INPUT_VERIFICATION = RESULTS / "input_verification.json"

PART1_COMMIT_FILE = CONFIGS / "part1_commit.txt"
SEEDS_FILE = CONFIGS / "seeds.json"
FREEZE_TAG = "freeze-catalogue"

# ---------------------------------------------------------------------------------------------
# Part 1 (READ-ONLY). Override the root with env BIOMARK_P1_ROOT only for relocation.
# ---------------------------------------------------------------------------------------------
P1 = Path(os.environ.get("BIOMARK_P1_ROOT", str(PROJECT.parent / "project"))).resolve()
P1_SRC = P1 / "src"
P1_DATA = P1 / "data"
P1_CACHE = P1 / "cache"
P1_FEATURES = P1_CACHE / "features"
P1_RUNS = P1_CACHE / "runs"
P1_MODELS = P1_CACHE / "models"
P1_TEST_PREDS = P1_CACHE / "test" / "preds"
P1_CONFIGS = P1 / "configs"
P1_RESULTS = P1 / "results"
P1_RESULTS_RUNS = P1_RESULTS / "runs"
P1_TEST_RESULTS = P1_RESULTS / "test"     # guarded: biomark.guard.p1_test_metrics only

# raw DepMap 26Q3 files
P1_RAW_FILES = {
    "crispr": P1_DATA / "CRISPRGeneEffect.csv",
    "expression": P1_DATA / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
    "cn": P1_DATA / "OmicsCNGeneWGS.csv",
    "damaging": P1_DATA / "OmicsSomaticMutationsMatrixDamaging.csv",
    "hotspot": P1_DATA / "OmicsSomaticMutationsMatrixHotspot.csv",
    "fusion": P1_DATA / "OmicsFusionFiltered.csv",
}
P1_MODEL_CSV = P1_DATA / "Model.csv"
P1_MODEL_CONDITION_CSV = P1_DATA / "ModelCondition.csv"
P1_SCREEN_MAP_CSV = P1_DATA / "CRISPRScreenMap.csv"
P1_SCREEN_SEQUENCE_MAP_CSV = P1_DATA / "ScreenSequenceMap.csv"
P1_SHA256 = P1_CONFIGS / "depmap_local_sha256.txt"

# parquet caches (Entrez-string columns, ModelID index)
P1_CRISPR_PARQUET = P1_CACHE / "crispr.parquet"   # read ONLY via biomark.guard.load_effects
P1_RAW_CACHES = {k: P1_CACHE / f"{k}.parquet" for k in ("expression", "cn", "damaging", "hotspot", "fusion")}
P1_GENE_MAP = P1_CACHE / "gene_map.parquet"
P1_FEATURE_BLOCKS = ("E", "C", "C_instab", "M", "M_burden", "H", "F_gene", "F_pair", "meta", "meta_prov",
                     "composite", "P")

# configs
P1_SPLITS = P1_CONFIGS / "splits.csv"
P1_GENE_CLASSES = P1_CONFIGS / "gene_classes.csv"
P1_EVAL_GENES = P1_CONFIGS / "eval_genes.csv"
P1_PREDICTABLE_GENES = {"m5": P1_CONFIGS / "predictable_genes.csv", "m1": P1_CONFIGS / "predictable_genes_m1.csv"}
P1_PREDICTABILITY_RULES = P1_CONFIGS / "predictability.json"
P1_LIBRARY_GROUPS = P1_CONFIGS / "library_groups.csv"
P1_ONCOGENE_PAIRS = P1_CONFIGS / "oncogene_pairs.csv"
P1_BIOMARKER_PAIRS = P1_CONFIGS / "biomarker_pairs.csv"

# run IDs (dev out-of-fold predictions; cache/runs/<run_id>/pred.parquet, 1022 x 18434)
P1_RUN_IDS = {
    "m5": "20261027-161041_m5_avg_6dbebcf4",      # headline
    "m1": "20261026-222117_m1_E_resB0_91e52f9d",
    "m2n": "20261027-044450_m2n_E_paralogs_62ce1575",
    "m3n": "20261027-014824_g3_m3n_601e3f5e",
    "t": "20261026-224125_t_E_status_24b46e2a",
}
HEADLINE_RUN = P1_RUN_IDS["m5"]

# Part 1 write hazards (brief section 6): "module.GLOBAL" -> (Part 2 destination, functions whose
# definition-time defaults are rebound; "*" = every function/method of the module).
# biomark.p1.import_essml applies them. Readers keep their P1 defaults (e.g. modelstore.load_fold_models
# still reads P1/cache/models); features.save_raw / io.build_cache / io.write_gene_map are disabled.
P1_WRITE_HAZARDS = {
    "essml.models.m2_local.INNER_CACHE": (M2_INNER_CACHE, ()),
    "essml.models.m5_stack.INNER_CACHE": (M5_INNER_CACHE, ()),
    "essml.models.m5_stack.OUTER_CACHE": (M5_OUTER_CACHE, ()),
    "essml.modelstore.MODEL_DIR": (P1_MODEL_DIR_REDIRECT, ("save_fold_model",)),
    "essml.guard.ACCESS_LOG": (P1_GUARD_LOG_REDIRECT, ("*",)),
}


def ensure_dirs() -> None:
    """Create Part 2 cache/result directories (never anything under P1)."""
    for d in (CACHE, CACHE_P1_REDIRECT, M2_INNER_CACHE, M5_INNER_CACHE, M5_OUTER_CACHE,
              P1_MODEL_DIR_REDIRECT, RESULTS):
        d.mkdir(parents=True, exist_ok=True)


def is_under(path, root) -> bool:
    """True if ``path`` is ``root`` or inside it (resolved, case-insensitive on Windows)."""
    p = os.path.normcase(str(Path(path).resolve()))
    r = os.path.normcase(str(Path(root).resolve()))
    return p == r or p.startswith(r.rstrip("\\/") + os.sep)


def load_config(name: str):
    """Load ``configs/<name>`` (``.json`` appended if no suffix) from Part 2."""
    p = CONFIGS / name
    if not p.suffix:
        p = p.with_suffix(".json")
    with open(p, encoding="utf-8") as fh:
        return json.load(fh)


def seeds() -> dict:
    """Fixed seeds (configs/seeds.json): keys ``seeds``, ``subsample_seed``, ``permutation_seed``."""
    return load_config("seeds.json")


def part1_commit() -> str:
    """The Part 1 commit Part 2 is pinned to (configs/part1_commit.txt)."""
    return PART1_COMMIT_FILE.read_text(encoding="utf-8").strip()
