"""Single parameter loader for the Tier A engine (stats, en, arms, null, audit, pipeline).

Design numbers come only from ``configs/design_params.json`` (DESIGN §F), ``primary_endpoints.json``, ``seeds.json``
and ``tiers.json``; implementation tolerances from ``configs/pipeline_params.json``. This module is the leakscan's
"params loader" exemption: everything else in ``src/biomark`` reads numbers through the accessors below.
"""
from __future__ import annotations

import functools
import hashlib
import json

from . import paths as P


@functools.lru_cache(maxsize=None)
def _load(name: str) -> dict:
    return P.load_config(name)


def dp() -> dict:
    """configs/design_params.json."""
    return _load("design_params.json")


def pp() -> dict:
    """configs/pipeline_params.json (implementation constants)."""
    return _load("pipeline_params.json")


def endpoints() -> dict:
    return _load("primary_endpoints.json")


def seeds() -> dict:
    return _load("seeds.json")


def tiers() -> dict:
    return _load("tiers.json")


def sha256(name: str) -> str:
    return hashlib.sha256((P.CONFIGS / name).read_bytes()).hexdigest()


def config_hashes() -> dict:
    return {n: sha256(n) for n in ("design_params.json", "pipeline_params.json", "seeds.json",
                                    "primary_endpoints.json")}


# ------------------------------------------------------------------------------------------------
# Typed accessors (one place where every number used by the engine is named)
# ------------------------------------------------------------------------------------------------


@functools.lru_cache(maxsize=1)
def stats_params() -> dict:
    d, p = dp(), pp()
    c5, s1 = d["confounders_s5"], d["s1_s5"]
    return {
        "drop_r2": float(c5["drop_assoc_r2_gt"]),
        "collinear_r2": float(c5["collinear_with_design_r2_gt"]),
        "merge_min": int(c5["leverage"]["merge_levels_with_lt_lines_in_fit"]),
        "lev1_tol": float(c5["leverage"]["exclude_from_meat_if_1_minus_h_lt"]),
        "min_carriers": int(s1["min_carriers"]),
        "min_resid_var": float(s1["min_resid_var"]),
        "bh_q": float(s1["bh_q"]),
        "cont_min_nonnan": int(d["dictionary_s3"]["continuous_min_nonnan"]),
        "q_sc": float(d["md_strict_s5"]["q_sc_bh"]),
        "rank_tol": float(p["numerics"]["rank_tol_rel"]),
        "eig_floor": float(p["numerics"]["eig_floor_rel"]),
        "tiny": float(p["numerics"]["tiny"]),
        "eigh_min_size": int(p["numerics"]["cluster_eigh_min_size"]),
        "df_subtract": 2,  # DESIGN 5.3: df = min(n - rank(Z) - 2, nu)
    }


@functools.lru_cache(maxsize=1)
def en_params() -> dict:
    s4, p = dp()["surrogate_s4"], pp()["en"]
    return {
        "l1_ratio": float(s4["l1_ratio"]),
        "n_lambda": int(s4["n_lambda"]),
        "eps": float(s4["eps"]),
        "tol": float(s4["tol"]),
        "K": int(s4["sis_K_per_stage"]),
        "K_fallback": int(s4["sis_K_fallback"]),
        "n_pairs": int(s4["n_complementary_pairs"]),
        "freq_min": float(s4["selection_freq_min"]),
        "sign_min": float(s4["sign_consistency_min"]),
        "resid_sd_min": float(s4["drop_column_in_fit"]["resid_sd_lt"]),
        "carriers_min": int(s4["drop_column_in_fit"]["carriers_lt"]),
        "fidelity_class3": float(s4["fidelity_class3_max"]),
        "unreachable_lt": float(s4["en_unreachable_if_max_attainable_freq_lt"]),
        "pilot_genes": int(s4["pilot_genes"]),
        "pilot_max_min": float(s4["pilot_max_projected_wall_min"]),
        "stages_G": list(s4["stages"]["G"]),
        "stages_X": list(s4["stages"]["X"]),
        "max_iter": int(p["max_iter"]),
        "cache_entries": int(p["fwl_cache_max_entries"]),
        "gram_cache_max": int(p["gram_cache_max_features"]),
        "forced_exclude": list(p["forced_relations_exclude"]),
    }


@functools.lru_cache(maxsize=1)
def null_params() -> dict:
    n7, c14, p = dp()["null_s7"], dp()["compute_s14"], pp()["null"]
    return {
        "B_base": int(n7["B_base"]), "B_raised": int(n7["B_raised"]), "wall_cap_min": float(n7["wall_cap_min"]),
        "fdr_max": float(n7["fdr_emp_max"]), "s1_B": int(n7["s1_vectorised_B"]), "workers": int(n7["workers"]),
        "gene_subset_n": int(n7["gene_subset_n"]), "strata_merge_min": int(p["strata_merge_min"]),
        "ci_level": float(p["ci_level"]), "boot_B": int(p["bootstrap_B"]), "s1_chunk": int(p["s1_vectorised_chunk"]),
        "fg_workers": int(c14["foreground_workers_during_null"]), "all_workers": int(c14["workers"]),
        "chunk_genes": int(c14["chunk_genes"]), "blas_threads": int(c14["blas_threads_in_worker"]),
        "permutation_seed": int(seeds()["permutation_seed"]),
        "gene_subset_seed": int(seeds()[n7["gene_subset_seed_key"]]),
        "sec_per_min": 60.0,
    }


@functools.lru_cache(maxsize=1)
def replication_params() -> dict:
    r6, eb = dp()["replication_s6"], pp()["eb"]
    return {"alpha": float(r6["alpha_one_sided"]), "power_min": float(r6["power_min"]),
            "eb_n_sigma": int(eb["n_sigma"]), "eb_max_iter": int(eb["max_iter"]), "eb_tol": float(eb["tol"]),
            "eb_sigma_min_frac": float(eb["sigma_min_frac_of_min_se"]), "eb_sigma_max_mult": float(eb["sigma_max_mult"])}


def lists_json(obj) -> str:
    return json.dumps(obj, sort_keys=True, default=str)
