"""Run the independent reimplementation check (DESIGN 13.3) end to end and write results/gates/G_indep.json.

    python -X utf8 -m tests.independent.run_all [--workers 3]

Detailed per-test comparisons go to tests/independent/out/ (CSV)."""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import datetime as dt
import json
import subprocess
import time

import numpy as np
import pandas as pd

from . import counts, data as DA, delta, encheck, enet, s1check, selftest


def _js(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, float) and not np.isfinite(o):
        return None
    if isinstance(o, (set, tuple)):
        return list(o)
    raise TypeError(type(o))


def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, list):
        return [_clean(v) for v in o]
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def annotations() -> dict:
    """Static notes: readings where the design text is ambiguous, fixes to my own code, diagnoses."""
    return {
        "readings_R0": [
            "5.3 screen days: one droppable covariate = cached screen_days_dev_humagne ('screen days within Humagne'); "
            "pipeline_params lists screen_days as droppable, DESIGN lists it both in the library block and as droppable",
            "5.3 S1 design has no METH_absent / MSI_unknown (1.6: S1 drops NA lines per feature); own CN NaN -> 0 with CN_absent",
            "5.3 'Leverages depend on (NaN group, d) only': read as a computational remark; my h_ii uses the full Z "
            "(incl. own CN, EXPGATE). Engine p agree to 1e-13, so the engine also uses the exact full-model hat",
            "5.3 df 'n - rank(Z) - 2': rank(Z) counts the intercept; binary d min(., carrier clusters - 1) (A7 iv)",
            "5.6 F1: BH over all U tests except own EXPDOSE incl. collinear_with_design tests; list hits exclude "
            "collinear_with_design (5.3 'never counted as determinants'); same rule for F2, S1-pred, MD-strict, MD-global",
            "4.1 'NA -> 0' after standardisation once on D (mean/SD ddof 0 over non-NaN D values)",
            "4.1 FWL block = the 4.1 list only (no screen days, no droppable covariates); 5.3 < 2-line level merge applied",
            "4.1 SIS: top-300 |corr| union forced (gene-specific relations except self-null); 25-alpha grid from the "
            "full-D screened set; 1-SE with SD(fold MSE, ddof 1)/sqrt(5); stage X CV on each fold's own stage-G residual",
            "13.3 'beta within 1e-3 relative': which beta is unstated; pass/fail uses mean beta over 50 subsamples "
            "(the EN rank statistic, engine column 'beta'); full-D beta reported alongside",
            "5.5 library block also removes screen days (5.3 puts them in the library block); burden = MSI + DMG_BURDEN; "
            "primary-disease sensitivity = block added to the full model, same flag rule (5.5 does not restate it)",
        ],
        "fixes_to_my_code": [
            "1 counts.py: list hits had omitted the collinear_with_design exclusion (BH cut unchanged and identical to engine)",
            "2 s1check.py: dropped_cov comparison accepts the engine's '|' separator (reporting only)",
            "3 encheck.py: removed two dead lines referencing a non-existent column (crash)",
            "4 data.py: METH_absent derived as 'all METH features NaN' (453 D lines); the loader biomark.covariates.build "
            "defines it as 'not in the RRBS gene matrix' (456). Changed to the loader definition. Found because with the "
            "engine's lambdas plugged in (diagnostic) full-D betas still differed by 0.1-5%; three other EN readings "
            "(screen days in FWL, SD ddof 1, raw NA->0) were tried as diagnostics and did not explain it; not adopted",
            "5 enet.py: optional fixed-lambda argument for that diagnostic only",
        ],
        "task3_diagnosis": "after fix 4: lambda_G/X agree to 4e-7, full-D beta to 4e-5 for all 30 selected features, "
            "frequency max |diff| 0.02 (one subsample fit in 3 features of genes 5045/6927), selected sets identical. "
            "The only out-of-tolerance item is mean beta of 5045 x SIL-METH:604 (rel 1.25e-3 vs 1e-3): its full-D beta "
            "agrees to 6e-10, so the gap comes from the same 5045 subsample fit(s) whose active set differs "
            "(AMP:6427 0.08 vs 0.06, AMP:9123 0.10 vs 0.12) - a solver-tolerance (tol 1e-3) / SIS-boundary effect in "
            "one half-sample, not a rule difference. Strictly a FAIL of the beta tolerance; PASS if 'beta' means full-D beta.",
    }


def main(workers=3):
    DA.OUT.mkdir(exist_ok=True)
    tol = DA.design_params()["indep_s13"]
    t0 = time.time()
    rep = {"design_ref": "DESIGN 13.3; tolerances design_params.indep_s13", "tolerances": tol,
           "utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
           "commit": subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                                    cwd=DA.ROOT).stdout.strip(),
           "agent": "independent-checker", "reads": "D effects via guard.load_effects(scope='dev') only; no R",
           "self_test_cr2_max_rel_diff_vs_bruteforce": selftest.run()}

    # ---------------- task 2 first (its BH hits stratify the task-1 gene draw)
    c = counts.run()
    hit_genes = c.pop("_s1obs_hit_genes")
    s1o = c["S1-obs"]
    task2 = {"family_definition": counts.__doc__.strip(),
             "S1-obs": s1o, "S1-pred_extra": c["S1-pred"],
             "MD_lists_recomputed": c["md"], "gate": c["gate"]}
    task2["pass"] = bool(s1o["F1_match_exact"] and c["gate"]["exact_match"])
    task2["key_numbers"] = {"F1_tests": s1o["F1"]["n_tests"], "F1_bh_rejections": s1o["F1"]["bh_rejections"],
                            "F1_p_cut": s1o["F1"]["p_cut"], "F1_list_hits_mine": s1o["F1"]["list_hits"],
                            "F1_selected_engine": s1o["engine"]["F1_selected"],
                            "gate_recall_mine": f"{c['gate']['MD-strict_engine_file']['k']}/{c['gate']['MD-strict_engine_file']['n']}",
                            "gate_recall_attempt1": f"{c['gate']['attempt1']['k']}/{c['gate']['attempt1']['n']}"}

    # ---------------- task 1
    a, b, info = s1check.draw_genes(hit_genes)
    genes = a + b
    mine = s1check.compute(genes)
    eng = s1check.engine(genes)
    r1, both = s1check.compare(mine, eng, float(tol["s1_p_tol_abs_log10"]))
    worst = both.sort_values("dlog10p", ascending=False).head(10)
    r1.update(info)
    r1["genes_with_F1_hits"], r1["genes_without_F1_hits"] = a, b
    r1["worst"] = worst[["gene", "feature", "p", "p_eng", "dlog10p", "beta", "beta_eng", "se", "se_eng", "df",
                         "df_eng", "dropped_cov"]].to_dict("records")
    r1["drop_paths_exercised"] = both["dropped_cov"].value_counts().to_dict()
    r1["n_tests_with_merged_levels"] = int((both["merged_levels"] > 0).sum())
    r1["pass"] = bool(r1["n_match"] == r1["n_compared"] and r1["n_only_engine"] == 0
                      and r1["n_only_mine"] == r1["n_only_mine_untested"])
    both.drop(columns=["_merge"]).to_csv(DA.OUT / "task1_s1_obs_comparison.csv.gz", index=False)

    # supplement: 2 seeded genes whose engine universe is smaller than my D_g (genes with NaN effect rows);
    # checks that the extra pairs are exactly my untested ones (< 5 carriers in the gene's rows)
    gd = s1o["universe"]["genes_differing"]
    rng = np.random.default_rng(DA.seeds()["seeds"][3])
    sup = sorted(rng.choice(gd, 2, replace=False).tolist())
    m2 = s1check.compute(sup)
    rs, b2 = s1check.compare(m2, s1check.engine(sup), float(tol["s1_p_tol_abs_log10"]))
    rs.update({"genes": sup, "seed_key": "seeds.json:seeds[3]"})
    rs["pass"] = bool(rs["n_match"] == rs["n_compared"] and rs["n_only_engine"] == 0
                      and rs["n_only_mine"] == rs["n_only_mine_untested"])
    r1["supplement_nan_row_genes"] = rs
    s1o["universe"]["genes_differing"] = len(gd)

    # ---------------- task 3
    ea, eb, einfo = encheck.draw_genes()
    egenes = ea + eb
    R = enet.run_many(egenes, workers=workers)
    eeng, lam = encheck.engine(egenes)
    r3, eboth, lamdf = encheck.compare(R, eeng, lam)
    r3.update(einfo)
    r3["genes_with_EN_obs_selection"], r3["genes_without"] = ea, eb
    nz = eboth[(eboth["freq"] > 0) | (eboth["selection_freq"] > 0)]
    r3["worst_freq"] = nz.sort_values("dfreq", ascending=False).head(10)[
        ["gene", "feature", "freq", "selection_freq", "sel_mine", "sel_eng"]].to_dict("records")
    sb = eboth[eboth["sel_mine"] & eboth["sel_eng"]].copy()
    sb["rel_beta_mean"] = (sb["beta_mean"] - sb["beta"]).abs() / sb["beta"].abs()
    r3["worst_beta"] = sb.sort_values("rel_beta_mean", ascending=False).head(10)[
        ["gene", "feature", "beta_mean", "beta", "rel_beta_mean", "beta_full_x", "beta_full_y"]].rename(
        columns={"beta": "beta_mean_eng", "beta_full_x": "beta_full_mine", "beta_full_y": "beta_full_eng"}).to_dict("records")
    r3["selected_differences"] = eboth[eboth["sel_diff"]][["gene", "feature", "freq", "selection_freq",
                                                            "sign_consistency_x", "sign_consistency_y",
                                                            "near_threshold"]].to_dict("records")
    r3["lambda_table"] = lamdf.to_dict("records")
    eboth.drop(columns=["_merge"]).to_csv(DA.OUT / "task3_en_obs_comparison.csv.gz", index=False)

    # ---------------- task 4
    pairs, dinfo = delta.draw_pairs()
    dres = delta.run(pairs)
    dres.to_csv(DA.OUT / "task4_delta.csv", index=False)
    dpath = DA.ROOT / "results" / "confounders" / "delta.csv"
    r4 = dict(dinfo)
    r4["engine_file_exists"] = dpath.exists()
    r4["pairs"] = dres.to_dict("records")
    r4["n_pairs_flagged_any_block"] = int((dres["flagged_blocks"] != "").sum())
    r4["full_model_vs_md_strict_beta_max_rel"] = float(
        ((dres["beta_full"] - dres["beta_md_strict_engine"]).abs() / dres["beta_md_strict_engine"].abs()).max())
    r4["full_model_vs_md_strict_p_max_dlog10"] = float(
        (np.log10(dres["p_full"]) - np.log10(dres["p_md_strict_engine"])).abs().max())
    if dpath.exists():
        r4["comparison"] = "engine delta.csv exists - compare manually"
        r4["pass"] = None
    else:
        r4["pass"] = None
        r4["status"] = "reported only (results/confounders/delta.csv does not exist)"

    rep["task1_s1_obs"] = r1
    rep["task2_bh_and_gate"] = task2
    rep["task3_en_obs"] = r3
    rep["task4_confounding_delta"] = r4
    rep["overall"] = {"task1": r1["pass"], "task2": task2["pass"], "task3": r3["pass"],
                      "task4": "reported only (no engine file)"}
    rep["notes"] = annotations()
    rep["wall_s"] = time.time() - t0
    out = DA.ROOT / "results" / "gates" / "G_indep.json"
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(_clean(json.loads(json.dumps(rep, default=_js))), fh, indent=1)
    print(json.dumps(rep["overall"]), rep["wall_s"])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=3)
    main(ap.parse_args().workers)
