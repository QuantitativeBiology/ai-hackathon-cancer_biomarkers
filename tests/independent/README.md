# Independent reimplementation check (DESIGN §13.3)

This directory holds a second implementation of the core statistics. It was written only from DESIGN.md (§1–§6, §9.3, §13.3,
amendments A1–A7), `configs/design_params.json`, `configs/pipeline_params.json` (implementation constants only),
`results/deviations.csv` (DEV-DF-FALLBACK, DEV-CAL-CONSERVATIVE) and the data loaders. None of the engine
modules (`stats.py`, `en.py`, `arms.py`, `blocked.py`, `null.py`, `pipeline.py`, `endpoints.py`, `audit.py`),
`docs/impl_notes_pipeline.md` or their tests were opened. Engine outputs are read only for comparison, through
`biomark.endpoint_io.read_arms` and `results/{md,gates,benchmark}`. Effects come only from
`guard.load_effects(scope='dev')`. R is never read.

Run (about 15 min, at most 3 worker processes):

    python -X utf8 -m tests.independent.run_all --workers 3

This writes `results/gates/G_indep.json` and the detailed comparison tables in `tests/independent/out/` (the large ones as `.csv.gz`, which git ignores; only `task4_delta.csv` is small enough to commit).

| module | content |
|---|---|
| `data.py` | data access: covariate cache, dictionary caches and metadata, D effects through the guard, subsamples |
| `confounder.py` | §5.3 confounder design, per-test drops, CR2 (PatientID clusters, HC2 singletons), df rule with the A7(iv) fallback, S1 p |
| `selftest.py` | synthetic check of my CR2 against a brute-force full-matrix CR2 and statsmodels HC2 |
| `counts.py` | task 2: F1/F2 BH from the engine's p values with my family definition; MD-strict and MD-global recomputed; gate recall |
| `s1check.py` | task 1: seeded gene draw, my S1-obs against the engine's |
| `enet.py`, `encheck.py` | task 3: two-stage fixed-λ stability elastic net (§4.1) for EN-obs, and the comparison |
| `delta.py` | task 4: §5.5 drop-one-block confounding deltas |
| `run_all.py` | runs everything and writes G_indep.json |

## Readings fixed before any comparison (R0)

The full statements are in the module docstrings. The points where the design text leaves a choice:

1. **S1 design and droppable covariates.** Z is the intercept, lineage (reference `_rare`), library group,
   CellFormat, SerumFree, medium family (reference: the most frequent level in the fit rows), doubling time and
   DT_unknown. To that I add MSI (reference MSS), DMG_BURDEN, TP53 (reference WT), screen days, own CN (NaN → 0,
   with CN_absent) and EXPGATE. Reference levels and the droppable grouping come from `pipeline_params.covariates`.
   - "Screen days" is the cached within-Humagne deviation (`screen_days_dev_humagne`). §5.3 names only that
     covariate, and it lists screen days as droppable.
   - METH_absent and MSI_unknown are not in the S1 design. §1.6 says S1 excludes NA lines per feature.
2. **Leverage for CR2.** I use the exact hat matrix of the full model, own CN and EXPGATE included. This means I
   read "leverages depend on (NaN group, d) only" (§5.3) as a statement about computation, not as a change to
   the estimator.
3. **df.** `n − rank(Z) − 2`, where rank(Z) counts the intercept. Binary d uses `min(that, carrier clusters − 1)`
   (A7(iv)). Continuous d uses `min(that, ν)`.
4. **F1.** Every (g, d) in U with d ≠ own EXPDOSE and a finite p. BH runs over the whole family, so tests that are
   `collinear_with_design` take part in the BH step-up. They are then removed from the *list*, because §5.3 says
   they are "never counted as determinants".
5. **EN.** Features are standardised once on all D lines (mean, and SD with ddof 0, over non-NaN values), and NA
   then becomes standardised 0.
   - The FWL block is the §4.1 list. It has no droppable covariates and no screen days.
   - SIS takes the top 300 by |corr| and adds the forced terms (their union).
   - λ uses the sklearn objective, a 25-α grid starting from α_max of the full-D screened set, and 5-fold CV
     over `cv_fold` with SIS inside each fold. The 1-SE rule uses SD(fold MSE, ddof 1)/√5.
   - Stage X is cross-validated on each training fold's own stage-G residual at λ_G.
6. **Delta.** The library block also removes screen days, which §5.3 places in that block. Burden is MSI plus
   DMG_BURDEN. The OncotreePrimaryDisease sensitivity block (levels with ≥ 5 D lines) is added to the full model,
   and I apply the same flag rule to it.

## Fixes to my own code (each a clear error, made without tuning toward the engine's numbers)

1. `counts.py`: I first counted every BH rejection as a list hit. This left out the §5.3 rule that
   `collinear_with_design` pairs are never determinants. The BH cut was unchanged: it was identical to the
   engine's before and after the fix.
2. `s1check.py`: the `dropped_cov` comparison now accepts the engine's `|` separator. This affects reporting only.
3. `encheck.py`: removed two leftover lines that referenced a non-existent column, which had crashed the script.
4. `data.py`: the covariate cache lacks METH_absent, so I derived it myself. I first defined it as "every METH
   feature is NaN", which gives 453 D lines with RRBS. The loader I am allowed to read, `biomark.covariates.build`,
   defines it as "line not in the RRBS gene matrix" (`cache/dictionary/rrbs_gene.parquet`), which gives 456.
   That is a data-input error on my side, so I now use the loader's definition. It touches only the EN FWL block:
   S1 under R0 does not use METH_absent.
   - What exposed it: with the engine's λ plugged in, as a diagnostic, my full-D β still differed by 0.1–5%.
   - I tried three other EN readings on the same diagnostic (screen days in the FWL block, SD with ddof 1, raw
     NA→0 before standardising). None explained the gap, so I did not adopt any of them.
   - After the fix, gene 8880's full-D β agrees to 3e-8.
5. `enet.py`: added an optional fixed-λ argument, used for that diagnostic only. The pass/fail numbers always use
   my own CV λ.

## Results (run of 2026-10-29; details in `results/gates/G_indep.json`)

- **Task 1, S1-obs: PASS.**
  - Sample: 20 genes drawn with `seeds.json:seeds[0]`. Ten have F1 hits in my task-2 BH and ten do not.
  - All 59,955 (g, d) tests match. The largest |Δlog10 p| is 8.5e-14.
  - β, SE, df, ν, n, carriers, carrier clusters, `dropped_cov`, `r2_design` and leverage exclusions are identical.
  - Supplement: 2 genes with NaN effect rows (`seeds[3]`). 5,733 tests match. The 264 pairs that only I have
    are exactly my untested ones (< 5 carriers).
- **Task 2, BH counts and gate: PASS (exact).**
  - F1 has 6,559,561 tests and a p cut of 2.556e-5, identical to the engine.
  - There are 3,358 BH rejections. 216 of them are `collinear_with_design`, which leaves 3,142 list hits: the same
    set as the engine's.
  - F2 gives 333. S1-pred gives 276,557.
  - MD-strict (11,662) and MD-global (1,380) recompute as identical sets.
  - Gate recall is 21/39, matching attempt 1 per family and per row.
- **Task 3, EN-obs: FAIL on the strict β tolerance only.**
  - Sample: 20 genes drawn with `seeds[1]`; 10 have EN-obs selections.
  - λ agrees within 4e-7.
  - Frequencies differ by at most 0.02.
  - The selected sets are identical: 30 features in both.
  - Full-D β agrees within 4e-5.
  - Mean β over the 50 subsamples is within 1e-3 for 29 of the 30 selected features. The exception is gene 5045 ×
    SIL-METH:604, at 1.25e-3. This traces to one half-sample fit with a different active set.
- **Task 4, delta: reported only.** `results/confounders/delta.csv` does not exist.
  - Sample: 10 MD-strict pairs drawn with `seeds[2]`.
  - My full-model β equals the engine's D-stage β within 2e-15.
  - One pair is flagged: 79897 × SIL-other:1289 on the lineage block (β 0.200 full vs 0.056 without lineage).

### Task 4 follow-up: engine comparison (`python -X utf8 -m tests.independent.delta_compare`)

- **Engine file.** `results/confounders/delta.csv` has the same sha256 as its row in `bulk_manifest.csv`.
  `delta_flags.csv` has no manifest row.
- **Pairs.** All 10 seeded pairs and all 80 (pair, block) rows are present in the engine's file.
- **Flags: PASS.** Flags agree for all 80 block rows and all 10 pairs. Both sides flag 79897 × SIL-other:1289 on
  lineage, and neither marks any pair subtype-confounded.
- **Numbers, 7 of 8 blocks: PASS.** For lineage, culture, growth, burden, TP53, own CN and primary disease,
  β_full, β_minus, SE_minus and Δβ agree within 1.7e-12 relative (70/70 rows).
- **Library block: FAIL (10/10 rows).** Relative β_minus differs by up to 0.029 and Δβ by up to 0.91 (where Δβ is
  small).
  - Cause: the engine drops only library_group, while I also drop screen days, following the §5.3 block
    definition (reading 6 in the readings list above).
  - Re-running the engine's reading as a diagnostic (not adopted) matches it to 8e-15.
  - This is a choice of reading, not a numerical error on either side.
  - The design text supports my reading. pipeline_params supports the engine's: it lists screen days as droppable
    and keeps the library block to library_group.
  - It changes no flag in this sample.
