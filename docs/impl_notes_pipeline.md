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

# Endpoint code (gate, P1 replication, P2 CGT recall)

| File | Content |
|---|---|
| `src/biomark/endpoints.py` | pure functions: `descends_from` / `require_descends` (tag is an ancestor of HEAD), `require_single_r_read`, `catalogue_manifest` / `catalogue_mismatches`, `clopper_pearson`, `mcnemar_exact`, `p1_verdict`, `p2_verdict`, `gate_attempt`, `RankIndex` (fast `rank_of_feature_D`, equality with `arms.rank_of_feature` tested), `p2_test_set`, `p2_score`, `p2_summary` (recall + CI, direction accuracy, sign errors, McNemar, verdict), `actionability`, `population_by_lineage`, `expected_nonnan`, `pool_power`, `p1_k`, `build_lists`, `strata_of`, `RFitter` (§5.3 model with fixed drops on any line set, CR2 + ν), `replicated`, `bootstrap_rates`, `p1_contrasts` |
| `src/biomark/endpoint_io.py` | real-data loaders for any line set (values, covariates, own CN, EXPGATE, library strings), σ̂_resid,D(g), `subset_data`, `make_pool`, `build_requests` |
| `scripts/run_gate.py` | §9.3 attempt n -> `results/gates/G_pc_attempt{n}.json` (never overwritten) + row in `results/gate_iterations.csv` (commit, reason, toy test, gate before/after, file sha256); attempt > 1 needs `--reason` and `--toy-test`; refuses unless HEAD descends from `freeze-bounds` |
| `scripts/build_p1_pool.py` | pre-tag: Pool_G, k_G, k_s, frozen lists (k_G, strata, secondary k) -> `results/benchmark/p1_pool.json` (+ bulk pool parquet), `configs/catalogue/r_requests.csv`, `configs/freeze_catalogue.json`; R contributes molecular data only; refuses before `freeze-bounds` or after `freeze-catalogue` |
| `scripts/replicate.py` | the single R run and `--pseudo-r fold0` (writes only `results/dryrun/`) |
| `tests/test_endpoint_branches.py` | 16 tests: every P1 branch (both contrasts, three branches + descriptive), end-to-end bootstrap reaching "predictor adds value", every P2 branch + descriptive, headline set (A4), gate pass/fail branches, refusals (pre-tag replicate / gate / pool build, second R open after the tag before any effect read, reproduction rules, catalogue hash mismatch), expected-NaN rule, pool/power/list builder |

Ambiguities resolved:
1. **Who freezes the pool.** §6.1 fixes k_G, the pool and power "before the tag", so they are built by a separate pre-tag
   script (`build_p1_pool.py`) that also writes `r_requests.csv` (reasons: P1 lists at k_G, every stratum and
   secondary k; MD-strict; EN-obs-strict; S1-obs BH; gate controls; CGT tests incl. proxy copies; an `extra` hook
   for cross-lineage / shortlist rows from M5) and the `freeze_catalogue.json` manifest that `replicate.py` checks.
2. **Pseudo-R dry run.** §1.3 calls it in-sample; the coordinator asked for refits on the other folds. Lists stay the
   frozen full-D ranks (re-running the arms on D′ is out of scope); everything else on the D side is refitted on
   D′ = folds 1–4: S1-obs (β, SE, drops, SD(d̃)) for the requested genes with the unchanged engine, σ̂_resid, the
   expected-NaN rule, EB shrinkage, pool and power. Fold-0 effects are read with scope dev; set A is not read
   (`setA_lines.csv` has status dryrun); `registry_status.json` = dryrun, uncitable. `--no-refit` and `--boot-B`
   exist for the dry run only.
3. **σ̂_resid,D(g)** = gene-level residual SD of y_g on [1, design blocks, droppable covariates]. **d̃_R for power**
   omits the gene covariates (own CN, EXPGATE); the R test itself uses the full model with the D drops.
4. **Replication sign.** A list row replicates with its arm's sign (EN majority sign; S1 sign of β); `r_tests.csv`
   tests CGT/control rows in their pre-registered direction, others in the direction of β_D. The degenerate rule
   (< 2 distinct carrier or non-carrier patients) applies to binary features, counting original patients.
5. **Single R read.** Any granted `replication` row in the guard log refuses a new run; `--reproduce` allows only the
   same commit and config digest (stored in the guard purpose string), per §1.3 (h). Set-A (`limits`) reads do not
   count as R reads.
6. **P2 headline set** is recomputed (tiers T1–T3, A4; S/W; sensitivity + bypass-resistance; concordant KO classes;
   proxy = 0; `headline_exclude` empty; `gate_disjoint` drops `in_gate`), with testable = `feature_D` in U with ≥ 5
   carriers for some g in T (CGT symbols mapped to Entrez). The CGT file's `p2_scored` column is not used, so any
   disagreement shows up.
7. **Gate features** = primary + secondary ids, plus the fallback ids when `gate_feature = fallback`.
8. **Q4c population**: indicated lines = OncotreeCode exactly in `indication_codes` (no OncoTree hierarchy, so parent
   codes such as SOLID match nothing unless tumour-agnostic); population precision and the carriers-only OR are
   not implemented yet. P1 secondary variants (U_all, gene-matched, Deming, XLIB) are not implemented; the
   secondary-k lists are frozen and requested.
