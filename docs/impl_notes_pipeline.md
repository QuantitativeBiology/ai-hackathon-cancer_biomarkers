# Implementation notes: Tier A engine (stats, en, arms, null, audit, pipeline)

Author: pipeline-author agent. Design: DESIGN.md v3.2 (tag `freeze-design`). Numbers: `configs/design_params.json`,
`primary_endpoints.json`, `seeds.json`; implementation constants only in the new `configs/pipeline_params.json`
(tolerances, covariate coding, solver limits, name maps), loaded through `src/biomark/params.py` (the leakscan's
params-loader exemption). No gene symbol or ACH id is a string literal in the new modules (checked with an AST
scan against `gene_map` symbols; note `CR2` and `TP53` are HGNC symbols, so `se_type` and the TP53 covariate name
come from config).

## Modules and API

| Module | Entry points |
|---|---|
| `stats.py` | `orth`, `resid`, `r2_on`, `bh`, `bh_cutoff`, `CovDesign` (coding + per-fit level merge), `DropInfo` (identity part of the drop rule), `Clusters`, `cr2_core` (vectorised CR2/HC2 + ν), `pair_test` (one (g, d); reusable by `replicate.py`; `one_sided='neg'/'pos'`), `s1_scan` (all (g, d) of U for one or more targets), `add_bh`, `eb_scale_mixture`, `shrunk_beta`, `se_r_pred`, `power_t`, `nu_from` |
| `en.py` | `make_subsamples`, `read_subsamples`, `enet_gram`, `alpha_grid`, `one_se`, `ENEngine.fit_gene`, `standardise`, `rank_en` |
| `arms.py` | `run_s1`, `universe`, `run_en` (process pool over gene chunks), `strict` (MD-strict / EN-obs-strict), `assemble`, `add_dedup`, `rank_of_feature`, `run_arms`, `observed_cutoffs` |
| `null.py` | `null_genes`, `strata`, `permutation`, `run_null` (fixed thresholds, B rule, workers), `summarise` |
| `audit.py` | `pair_fit`, `confounding_delta`, `locus_resolution`, `amplicon_test`, `cross_lineage`, `self_null`, `pair_classes`, `gene_classes`, `cgt_scores`, `gate`, `run_audits` |
| `pipeline.py` | `Data` (one container: `from_toy`, `from_real`), `toy_eligibility`, `run_pipeline(data_or_toy, *, fold_centre, null_B, arms, workers, null_workers, audits)` |
| scripts | `run_arms.py [--pilot] [--universe ctrl]`, `run_null.py [--s1-vectorised]` (both refuse to run before tag `freeze-bounds` except the pilot) |

## Design ambiguities resolved (all reversible; none changes a frozen number)

1. **S1 design vs NA indicators.** S1 excludes NA lines per feature family (§1.6), so CN_absent / MSI_unknown /
   METH_absent are *not* in the S1 design; they are in the EN unpenalised block only (§4.1). Under the null they
   move with 𝒟 rows in both arms (S1 through the NaN pattern of the permuted feature).
2. **Own CN as a gene covariate** is coded NaN → 0 plus a CN_absent indicator column (part of the own-CN covariate,
   dropped with it). The R² drop rule for own CN uses both columns; for EXPGATE its single column.
3. **Covariate-drop R²** is computed on each test's own fit rows (gene non-NaN ∩ feature non-NaN). "Variant group
   contains the covariate's feature" uses `variant_group` of `dictionary_v1.json`; TP53/MDM2 derivation resolves the
   symbols in `pipeline_params.json` to Entrez at run time. "CN-derived" = AMP, DEL, LOSS, SIL-DEL.
4. **Level merge** (< 2 lines in the fit rows → reference). Reference = configured level (lineage `_rare`, MSI `MSS`,
   TP53 `WT`) when present in the rows, else the most frequent level in the rows.
5. **df = min(n − rank(Z) − 2, ν)** exactly as written, with rank(Z) including the gene covariates.
6. **CR2 clusters**: symmetric (I − H_cc)^(−1/2); 2 × 2 closed form, eigh for larger blocks; eigen-directions with
   eigenvalue ≤ `eig_floor_rel` are excluded from the meat and counted in `n_lev1_excluded` (with singleton rows
   1 − h < 1e-8).
7. **U** = §3 counts on the gene's non-NaN rows (binary ≥ `s1_s5.min_carriers`, continuous ≥ 20 non-NaN and SD > 0),
   identical for all arms. S1 rows with residual variance < `min_resid_var` stay in U with p = NaN (unrankable).
8. **BH families**: F1 and F2 (own EXPDOSE) within each S1 arm over all U rows with a p (collinear rows count in m,
   but are never selected). MD-strict q_sc = BH over all EN-pred selected candidates with a D-stage p (collinear
   count in m, never pass; candidates without a D-stage test, e.g. < 5 carriers on the gene's rows, are
   `d_stage_testable = False` and fail).
9. **Null cutoffs**: p* = largest p with q ≤ α in the observed family; if nothing passed, α/m (a single test at the
   first BH step). F2 gets its own cutoff (design names only p*_F1). The S1-only vectorised draws (B = 200) reuse
   seeds `permutation_seed + b`, so draws 0..B−1 are shared with the EN draws. The 50-draw rule uses the single-draw
   wall of draw 0, literally `50 × wall ≤ 180 min`.
10. **EN details**: features standardised once on all D rows (non-NaN mean/SD) then NaN → 0, never rescaled after
   FWL (needed for the exact FWL ≡ huge-scaled-unpenalised equivalence, tested). λ grid from the full-D fit of each
   stage (after SIS) and reused for the CV folds; held-out fold residualised with the training-fold FWL
   coefficients; stage X is CV'd on each fold's stage-G residual at λ_G. Fidelity = 1 − SSE/Σ ỹ_test² pooled over
   folds. `en_max_attainable_freq` = share of the 50 subsample fits in which the column passed the in-fit
   availability rule (exact, gene-specific). Rank ties broken by `gene|feature` id.
11. **Solver**: coordinate descent on the Gram matrix via scikit-learn's `_cd_fast.enet_coordinate_descent_gram`
   (the routine behind `enet_path`; private API, pinned sklearn 1.5.2; equality with `enet_path` is unit-tested).
12. **FWL cache**: in-memory LRU per worker (key = gene NaN pattern × fit), storing the residualised shared block
   (float32) and, when the shared block has ≤ 1,500 columns, its Gram. Memory-mapped files were not needed at the
   measured sizes; same semantics (read-only, computed on that fit's rows).
13. **Subsamples**: `configs/subsamples.json` is read if present (the dictionary agent's file has the same
   `subsamples` key). Otherwise each patient is stratified by its modal lineage and patients are assigned greedily
   to the half with fewer lines of that stratum, in a seeded random order.
14. **EN-pred blindness**: EN-pred consumes `Data.target('pred')` and `rows_ok` (NaN mask) only. Worker payloads for
   pred arms carry no observed effects (`Data.for_workers(include_obs=False)`). Test: `test_en.py::
   test_en_pred_is_blind_to_observed_effects` (monkeypatched `guard.load_effects` raises; EN-obs then fails, EN-pred
   runs, row sets equal).
15. **Audits (first implementation, for the toy interface and M5)**: co-location refit = model + neighbour indicator
   on rows without concordant carriers (d = nb = 1), d tested (a discordant-only refit cannot carry the §5.3
   design); cross-lineage uses raw lineages with ≥ 10 lines, per-lineage power from the §6.1 rule with σ̂_resid of
   the full model, `lineage-restricted` needs ≥ 1 lineage with F5 BH q < 0.05 and an OLS F-test interaction p < 0.05;
   NC3 `frac_dependent` = 1 if the carrier median effect < −0.5 (its mean over pairs is the NC3 statistic),
   `frac_carriers_dependent` alongside; pair classes 7 (R) and 8 (KO–drug) and the culture/hypermutation flags of
   §8.6–8.7 are not computed yet (pre-R classes only).
16. **EB shrinkage grid**: 20 log-spaced σ from 0.1·min(s) to 2·sqrt(max(b² − s²)), no point mass at 0 (design:
   "normal scale mixture, 20 log-spaced sigma").

17. **S1 list = F1** (`selected`); own-EXPDOSE BH hits are `selected_F2` (toysim interface "BH q < bh_q (F1)",
   amendment A6-iii). The null therefore uses p*_F1 only.
18. **Feature NaN patterns.** The real 𝒟 has 630 distinct NaN patterns (gene-level CN coverage differs by gene), so
   the S1 scan works per (gene NaN pattern × feature NaN pattern) group with integer group ids and a gene ×
   feature membership matrix. An early draft compared byte-string pattern keys in an object array and silently
   dropped features from U on real data (13,108 instead of 17,900 rows for 3 genes × 2 arms); the toy (few
   patterns) did not show it. Fixed before any real result was written; only the timing pilot touched real effects.
19. **Performance choices** (no effect on results): pivoted-QR `orth`; no full-matrix boolean row copies
   (`np.ix_`); EN screening from one mat-vec on the cached residualised block, columns materialised only for the
   screened/non-zero sets; Gram composed from the cached shared Gram plus gene-specific cross products; single-α
   fits without a cached Gram use the X-based CD of the same solver; S1 and EN parallel over gene chunks sorted by
   NaN pattern (`compute_s14.chunk_genes`), equality with serial tested.

## Timing (real D: 1,022 lines, 12,813 feature columns = 2,993 shared + gene-specific; 2,193 eligible genes)

`results/timing/surrogate_pilot.json` (50 seeded eligible genes, 5 workers, run while the 6-worker toy suite shared
the 12-core machine): S1 both arms 203 s serial (4.1 s/gene incl. ~15 s per NaN-pattern group) -> 12.4 min projected
at 12 workers; EN both arms 683 worker-s (13.7 s/gene incl. ~30 s cold FWL cache per worker per pattern; warm
≈ 2 s/gene/arm) -> 41.6 min projected (< 45 min, so K stays 300; warm estimate ≈ 25 min); one fixed-threshold null
draw on 300 genes ≈ 29 min single worker -> B = 20 ≈ 96 min on 6 workers; B = 50 would need 50 × 29 min > 180 min,
so the rule keeps B = 20. A first pilot run put all 50 genes in one worker chunk and over-projected EN by ~4×
(fixed: chunk = genes / workers, projection from summed worker CPU). Selections were identical across both pilot
runs (deterministic). Real EN-pred selects ~20 features per gene (1,025 / 50 genes) versus 6 for EN-obs; the median
EN-pred fidelity is 0.29 (consistent with amendment A6: prediction arms are dense; MD-strict is the catalogue).

## Toy answer-key run (`BIOMARK_TOY_PIPELINE=1 pytest tests/test_toy_answers.py -k pipeline`, 74 min)

14 passed, 5 failed (bounds and answer key untouched):
- `fold_offset_artefact` (positive control, without centring >= 3/5): 0/5. The artefact features reach selection
  frequency 0.42-0.50 (far above background) but < 0.6: λ is CV'd over Part 1 cv_fold, and a fold-specific offset
  is by construction unpredictable out of fold, so λ_1SE is high. Design property, not a code fault.
- `gate`: KeyError from a column collision (my gate table repeated `gate_eligible`/`family`); fixed (gate output =
  control_id + recovered*); recall 1.0 on the same toy outputs (re-run of the flagged suite pending).
- `error_rates_smooth_toy`: S1-obs FDR_emp_fixed 0.036 < FDP_true 0.090 − 0.05. The 10 false S1-obs hits are
  background pairs (not proxies of planted determinants); the null's ~4.8 hits/draw match nominal BH. No code cause
  found; looks like realised-FDP variability of the univariate list.
- `error_rates_ridge_toy`: MD-strict FDP_true 0.228 > 0.15 (23 background pairs; EN-pred FDP 0.51, S1-pred 0.79).
  Mechanism as in A6: ridge OOF predictions inherit training-fold noise of y, and the D-stage screen uses the same D
  effects, so the screen does not remove these candidates.
- `pure_null`: S1-pred P(any non-neutral hit) 0.25 > 0.10 (5/20 toys, 1-2 background hits each); S1-obs 0/20. The
  assertion stops at S1-pred, so MD-strict was not evaluated by the test.
