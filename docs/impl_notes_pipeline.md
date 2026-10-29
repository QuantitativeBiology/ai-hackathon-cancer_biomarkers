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
| `scripts/run_gate.py` | §9.3 attempt n -> `results/gates/G_pc_attempt{n}.json` (never overwritten) + row in `results/gate_iterations.csv` (commit, reason, toy test, gate before/after, file sha256); attempt > 1 needs `--reason` and `--toy-test`; refuses unless HEAD descends from `freeze-bounds`. Each control is read from the universe of its gene (eligible lists, or U_ctrl: `results/md_ctrl/`, `cache/benchmark/arms_ctrl/`) for every arm and both MD lists; `universe` per row, `summary.by_universe`, input hashes of both datasets; refuses if a control gene is in neither universe or U_ctrl outputs are missing / stale |
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

# Remaining Tier A pieces (XLIB, component readout, confounder audit, Q4, registry)

| File | Content | Timing stem(s) |
|---|---|---|
| `src/biomark/analyses.py` | `append_registry` (file-locked `results/registry.csv`, §13.6 columns), `write_timing`, `xlib_split`, `discovery_eligibility`, `component_targets`, `data_with_pred`, `component_selection`, `design_vif` (alias-aware), `oncotree_terms`, `indicated_lines`, `population_precision` | — |
| `scripts/run_xlib.py` | §6.3 S1-obs + EN-obs refitted on discovery lines (own λ CV, folds / subsamples restricted), pool / power / k / lists by the P1 rules (replication molecular data), replication fits on non-Avana D lines, paired patient bootstrap, contrast EN-obs − S1-obs, `u_rule` all_D and discovery_only -> `results/benchmark/replication_xlib.csv` (+ `xlib_summary.json`, `cache/benchmark/xlib_arms.parquet`) | xlib_s1_obs, xlib_en_obs |
| `scripts/run_component_readout.py` | §4.2 EN (unchanged engine) on fold-centred ŷ_M1, ŷ_M2n − ŷ_M1, ŷ_M3n − ŷ_M1 for MD-strict genes -> `results/components/component_selection.csv` | component_readout |
| `scripts/run_confounder_audit.py` | §5.5 deltas for MD-strict ∪ S1-obs BH ∪ gate controls (parallel chunks of `audit.confounding_delta`) -> `results/confounders/delta.csv`, `delta_flags.csv`; §5.3 VIF / condition number -> `results/covariates/vif_check.json` | confounder_delta, vif |
| `scripts/run_q4.py` | §10.5 Q4a/P2 (all sets / subsets / k), direction, Q4b actionability, Q4c population labels + precision + carriers-only OR -> `results/clinical/{q4_summary.json, p2_scores.csv, actionability.csv, q4c_population.csv, q4c_labels.csv}` | q4_summary |
| `scripts/replicate.py` | now appends a registry row (status ok, or dryrun for `--pseudo-r`) and writes stem `replicate` (the dry run writes its timing into `results/dryrun/`, never an unmapped stem) | replicate |
| `scripts/run_arms.py` | `--universe ctrl` (genes = `configs/u_ctrl_genes.txt`, asserted disjoint from the eligible set; BH within U_ctrl) writes only ctrl-specific paths (`output_paths`): `cache/benchmark/arms_ctrl/` (+ `_stage_arms_ctrl/`), `results/benchmark_ctrl/` (selected, lambda, p_cut, summary incl. `u_ctrl_genes_sha256`; timings in `results/benchmark_ctrl/timing/`, never budgeted) and `results/md_ctrl/`; `results/bulk_manifest.csv` rows are keyed by the ctrl part paths | unchanged |
| `scripts/build_u_ctrl_genes.py` | U_ctrl = (positive-control dependency genes ∪ CGT target genes) minus eligible, from frozen configs only (CGT symbols mapped through `eligible_genes.csv`) -> `configs/u_ctrl_genes.txt`; `--check` compares | — |
| `tests/test_tier_a_analyses.py` | 13 tests (XLIB split, discovery eligibility, component labels, component target with the unchanged engine, VIF + alias detection, precision / OR vs scipy, parent-code mapping, registry / timing writers, every scripted stem is in `stage_job_map`, all four scripts refuse before freeze-bounds) | — |

Ambiguities resolved:
1. **XLIB refit.** The two arms are refitted on discovery lines with the engine unchanged; the 50 subsamples and the
   Part 1 folds are the full-D ones restricted to discovery lines (still complementary and patient-grouped);
   standardisation is on discovery lines. The replication string rule ("no Avana in the library string") is read
   as no `+`-component equal to the discovery library. `discovery_only` recomputes r_within on discovery lines
   with the class and tier criteria of `eligible_genes.csv` and restricts the pool; the lists are re-cut in that
   pool from the same discovery ranks (no refit).
2. **Component readout.** "Full-D fits": a component pair counts as present when its full-D β at λ_1SE ≠ 0
   (frequency and the stability selection are reported alongside). layer = expression-mediated (M1 only),
   genomics-layer (M2n − M1 and/or M3n − M1 only), both, none; `carrier_input` names the component inputs for
   proxy-class features.
3. **VIF.** Computed on the design as the scan codes it (merge rule on all D rows). Real D: serumfree = unknown and
   medium_family = unknown are the same 14 lines, so two culture dummies are exactly aliased (condition number
   ~1e16). The estimating code is unaffected (rank-revealing QR); `design_vif` reports the alias and drops it,
   giving max VIF 6.17 (≤ 10) and condition number 10.8. `results/covariates/vif.csv` belongs to the covariate builder
   and is written only if absent.
4. **OncoTree parents.** There is no OncoTree hierarchy file locally. A partial parent mapping is built from local
   files: a line is indicated if its OncotreeCode equals the code, or its Model.csv Subtype / PrimaryDisease /
   Lineage name equals an OncoTree term paired with that code in the frozen MOAlmanac release. This raises e.g.
   NSCLC 5 → 80, COADREAD 15 → 66, BRCA 20 → 51 and OVT 0 → 40 D lines. Codes without a term or a name match
   (e.g. SOLID, BTC, CESC, GB, TNBC, DLBCL, ALCL) fall back to exact matching and are listed as
   `unmapped_indication_codes` in `q4_summary.json`.
5. **Population precision** = share of carriers that are accurate-dependent (D rule), within and outside the
   indication. **Carriers-only OR** = conditional-MLE OR of accurate-dependent membership, indicated vs not, with an
   exact 95% CI and Fisher p (NaN when one side has no carriers).
6. **Not yet done in Q4**: the tier trend test (reported as not run when a tier has < 5 tests), the miss
   decomposition, and the tumour-agnostic ⇒ cross-lineage check (needs the §8.3 labels from M5).
7. Comment-only change in `stats.py` (`# literal-ok` on the EB Gaussian density) so the project leakscan passes on
   the engine; no semantic change.

# Fixes before freeze-bounds (F1–F6) and amendment A7

- **F1 df fallback (DESIGN §5.3, A7 iv; deviations.csv DEV-DF-FALLBACK).** `stats.shared_df` is the single df rule of
  every test and of power. cap = n − rank(Z) − 2. Binary d: df = min(cap, carrier clusters − 1), where carrier
  clusters = distinct PatientID clusters among the carriers in the fit rows (`stats.carrier_clusters`). Continuous d:
  min(cap, ν). df ≤ 0 gives p = NaN. The ν-rule df and p stay as the logged sensitivity (`df_nu`, `p_nu`), and
  `df_rule` is a column of every S1 / MD / R table. Selected by `pipeline_params.json:df_rule.binary`.
  - Where it applies: `s1_scan` (D, XLIB, the null); `pair_test` (`binary=None` detects 0/1 d from the raw
    values); `audit.pair_fit`; `RFitter` (R tests, bootstrap; feature types from the dictionary).
  - `cr2_core` only sees the residualised d̃, so it applies the fallback only when given `car_clusters`. Without
    them it is the ν rule, which is what `tests/test_inference_calibration.py` measures (see the report: that test
    must pass carrier counts to measure the active rule).
  - The same simulation with the fallback (seeds and counts unchanged) has max size 0.0587 at α 0.05 and 0.0011 at
    1e-3. Four cells are below the 0.035 floor: (255, 3, SD × 2) 0.0318; (255, 3, homoskedastic) 0.0104;
    (255, 5, homoskedastic) 0.0293; (1022, 3, homoskedastic) 0.0276.
  - 50-gene pilot, F1 BH hits: S1-obs 16 → 12 and S1-pred 5,996 → 5,264 (strictly nested). EN is unaffected.
- **F2 power df = test df.** `endpoints.pool_power` builds the R test's own design: `RDesign.model_Z`, shared with
  `RFitter.fit`, with the gene covariates and the D drops. It then calls `stats.shared_df` with the cap and the
  fallback. SD_R(d̃_R) uses the same residual. `RFitter.cov` is now an `RDesign`, a `CovDesign` that carries own
  CN / EXPGATE.
- **F3 NMF D-only invariance.** `dictionary.fam_prog` fits on D (unchanged) and runs `transform` separately per role
  (D, R, set A), in float64 with one BLAS thread (`threadpoolctl`). `tests/test_nmf.py` passes 3/3.
  - Dictionary rebuild: only `PROG.parquet` and `nmf_components.parquet` changed (50 PROG feature hashes; `sd_D` of
    11 PROG features changed at the 6th decimal). PSCORE and every other family are unchanged.
  - `eligible_genes.csv`, `eligibility.csv`, `positive_controls.csv` and `pc_eligibility.csv` are byte-identical
    after re-running `build_eligibility.py` and `build_positive_controls.py`.
- **F4 timing schema.** `analyses.write_timing` writes job, utc (end), start_utc, numeric wall_s, workers and a
  per-process run_id; it raises if wall_s cannot be determined.
  - `run_arms`: s1_obs and s1_pred share one window and run_id; en_obs and en_pred get per-arm windows recorded in
    `arms.run_en` (timing only); md_strict has workers 1. `pilot` records numeric wall_s and workers.
  - XLIB: per-stem windows. `run_gate`: stem gate_pc plus a registry row.
  - `run_null --s1-vectorised`: a sub-record `s1_vectorised` of `null_fixed.json`, never a new stem.
- **F5** gene constants of `scripts/build_cgt.py` -> `configs/cgt_build_rules.json` (`params.cgt_build_rules`).
  The unmodified script first reproduced every CGT output byte for byte, and so did the refactored one: the
  sha256 of `clinical_ground_truth.csv`, `cgt_curated_mapping.csv`, `cgt_features.csv`, `results/cgt/*`,
  `cgt_dropped.csv` and `cache/cgt/*` is unchanged.
- **F6** `literal-ok` reasons on guard.py 79/102/148/203 and p1.py 55/209/219/399. Project leakscan: 0 hits.
- **A7 MD-global** = EN-pred selected ∩ S1-obs F1 BH hits (not collinear), same sign (`arms.md_global`).
  - Outputs: `res/out["md_global"]` and `["md_global_all"]`, with flags `pass`, `s1_f1_hit`, `same_sign` and
    `in_md_strict`; `results/md/md_global{,_all}.csv`. In the fixed-threshold null it is the list `MD-global`
    (F1 p ≤ observed p*_F1, same sign).
  - Downstream defaults (`IO.catalogue()` = MD-global): audits' pair and gene classes; `r_requests` (reason
    MD-global; MD-strict kept); Q4b actionability; component readout and confounder audit (both lists, flagged);
    replicate `r_lines`.
  - The gate is decided on MD-strict, with MD-global reported. The P1/P2 contrasts are unchanged.
- **Re-verify fixes.** `analyses.write_timing` raises on a missing or invalid `workers`; there is no default, and every
  scripted call passes it (checked by AST in `tests/test_budget_subrecords.py`). Timing files may carry `subrecords`
  (`{name: {start_utc, utc, wall_s, workers, run_id}}`, built by `analyses.timing_subrecord`).
  - `budget_check` counts their intervals in the stage wall (union) and in concurrency (own run_id, default
    `<stem>:<name>`), and reports schema errors.
  - `run_null --s1-vectorised` writes `null_fixed.subrecords.s1_vectorised`. Before the main null exists, it is
    parked in `results/null/timing_s1_vectorised_subrecord.json` and merged by the main run, so a placeholder never
    marks the stage done.

# Worker payloads on Windows spawn (infrastructure fix after the first full arms launch)

The full real-D arms run crashed at worker spawn: the `Data` payload passed as `ProcessPoolExecutor(initargs=...)` is
pickled through the spawn pipe (OSError 22 / truncated pickle). `src/biomark/workerstore.py` now stages every pool
payload to `cache/workers/<run_id>/<pool>-<id>/`.
- Arrays ≥ 1 MiB (found recursively, including `Data` fields: X, Xstd, pred, own CN, EXPGATE, rows_ok, and `_obs` only
  for obs pools) go to `.npy`; everything else goes to `payload.pkl`. Only the directory path crosses the pipe.
- Workers rebuild the payload with `np.load(mmap_mode="r")` (read-only; `Data` is restored without `__post_init__`).
- `WS.staged(...)` deletes the directory when the pool exits.
- Used by `arms.run_s1`, `arms.run_en`, `null.run_null` and `run_confounder_audit`. XLIB and the component readout
  go through `arms`.
- The null payload no longer carries the observed arms table, which `_draw` never used.
- Pred-arm pools are staged from `Data.for_workers(include_obs=False)`, so no observed-effect file exists.
  `tests/test_workerstore.py` checks this, checks parallel (2 workers, mmap) == serial for S1 and EN on the toy,
  memmaps and cleanup, and runs parallel EN-pred under a raising `guard.load_effects`.

# M5 / M6 / M8 analysis scripts (subtypes, audits, combinatorial lookup, limits, candidate cards)

Author: M5-M8 analysis agent, after `freeze-bounds`, before any M5-M8 output existed. Catalogue = MD-global
(amendment A7) everywhere; MD-strict is computed alongside (flags `in_md_global` / `in_md_strict`, `lists`).
Implementation constants (not design numbers) are in the new `configs/m5_m8_params.json` (`m5.cfg()`); design
numbers come from `design_params.json` through `params`.

| File | Content | Timing stem(s) |
|---|---|---|
| `src/biomark/m5.py` | `Workspace` (project layout, or a toy root with the same layout; arms read through `endpoint_io.read_arms` with column / row pushdown, never whole), `RealSource` / `ToySource` (data loaders: `pipeline.Data`, D effects via `guard.load_effects(scope='dev')`, expression, controls, gene info, CGT, paralogs), `prepare_toy`, `catalogue`, `analysis_pairs`, `subtypes`, `modules`, `culture_audit`, `culture_overlap`, `hypermut_tp53`, `cross_lineage` (chunked), `eb_fit_streamed` / `shrunk_beta_ws`, `RankLookup`, `map_pairs` (workerstore pools) | — |
| `src/biomark/combinatorial.py` | readers of the Esmaeili Anvar 2024 harmonised table and the Flister 2024 supplement (to Entrez pair keys only), `label_of`, `candidates`, `background`, `validation`, `lookup` | — |
| `src/biomark/limits.py` | `pair_classes` (on `audit.pair_classes`), `gene_classes` (on `audit.gene_classes`), `pair_labels`, `q1b`, `known_recovery`, `shortlist` | — |
| `src/biomark/audit.py` | `cross_lineage` split into `cross_lineage_stats` (per pair) + `cross_lineage_labels` (BH F4/F5 + labels); `cross_lineage()` = labels(stats(...)), behaviour unchanged | — |
| `scripts/run_subtypes.py` | §8.1 / §8.2 -> `results/subtypes/{subtypes,subtype_lines,subtype_genes,modules,module_lines}.csv`, `summary.json` | subtypes, modules |
| `scripts/run_cross_lineage.py` | §8.3 -> `results/lineage/{cross_lineage,cross_lineage_per_lineage}.csv`, `summary.json` | cross_lineage |
| `scripts/run_m5_audits.py` | §8.4-§8.7 -> `results/locus/{resolution,amplicon}.csv`, `results/subtypes/{culture_audit,culture_overlap,hypermut_tp53}.csv`, `results/m5_audits_summary.json` | locus_resolution, amplicon_test, culture_audit, hypermut_tp53 |
| `scripts/run_combinatorial_lookup.py` | §11 -> `results/combinatorial/{lookup.csv, validation.json}` | combinatorial_lookup |
| `scripts/run_limits.py` | §12.1 / §12.2, Q1b, NC6 -> `results/limits/{pair_classes,gene_classes,pair_labels}.csv`, `limits_summary.json` | limits |
| `scripts/build_candidate_cards.py` | §12.3 -> `results/candidates/{known_recovery,shortlist}.csv`, `cards_summary.json` | candidate_cards |
| `tests/test_m5_m8_toy.py` | the six scripts end to end on a toy workspace + synthetic unit checks | — |

Every script refuses unless HEAD descends from `freeze-bounds` (`EP.require_descends`, as in `run_gate.py`), before
anything is read or written. `--toy DIR` runs on a toy workspace (prepared on first use by the unchanged engine;
never the project's `results/`). Each script writes timing through `analyses.write_timing` and a registry row (status
`toy` in a toy workspace). No R: no scope `replication`, no `p1.test_predictions`; R-dependent outputs read
`results/replication/*.csv` only if it exists.

Ambiguities resolved:
1. **Line share (§8.1).** β = the full-D EN-pred coefficient (`beta_full`), Σ over the EN-pred *stability-selected*
   features of g (both stages); x̃ = 𝒟 standardised on D and residualised on the EN unpenalised block (design blocks
   + NA indicators) on g's rows, x̄ = its mean. A catalogue feature with `beta_full = 0` (selected by stability only)
   gets no members and a note (counted in `summary.json:beta_full_zero`). Membership is computed for MD-global ∪
   MD-strict pairs. Per gene: `n_unexplained` = dependent but not accurate-dependent; `n_accurate_unassigned`.
2. **Module carriers (§8.2).** Binary d: x = 1. Continuous d: the dependency side of the D median
   (sign(β)·(x − median) < 0). "Accurate-dependent on ≥ 50% of module genes" uses the genes with a non-NaN status for
   that line. Effective clusters: average linkage on 1 − |r| (pairwise-complete D effects), cut at distance 1 − 0.5.
   Modules per list (MD-global; MD-strict alongside) and per sign.
3. **Cross-lineage pairs (§8.3).** MD-global ∪ MD-strict ∪ all positive-control features in U ∪ the S1-obs BH top 200
   by rank (`pipeline_params.audit.cross_lineage_top_s1`). One BH family (F4, F5) per run over all those pairs: pairs
   are fitted in parallel chunks (`cross_lineage_stats`) and BH + labels are applied once to the concatenated tables
   (`cross_lineage_labels`; tested equal to a single call). Shrunk β for the power rule: one EB fit over every S1-obs
   U row, streamed in row chunks (`eb_fit_streamed`: same grid, start, tolerance and cap as `stats.eb_scale_mixture`;
   tested equal); the posterior is read for the needed pairs only.
   **A7-iii limitation:** `single_lineage_driven` = cross-lineage label with ≤ 1 powered lineage that is nominally
   significant (p < 0.05) with the pooled sign. It is reported in `summary.json`; the label itself is not changed.
4. **Co-location / cutting pairs (§8.4-§8.5).** MD-global ∪ MD-strict ∪ S1-obs BH ∪ positive controls (as in the toy
   adapter), so passengers that appear only as S1-obs hits are resolved too; `audit.locus_resolution` and
   `audit.amplicon_test` are unchanged.
5. **Culture (§8.6).** Membership vector = subtype members vs the other lines with an effect for g. Variables:
   CellFormat, 3D (CellFormat = 3D), SerumFree, medium family, library group; flag on max V ≥ 0.3. Rescue: the
   confounder model restricted to the level of the flagging variable that holds most members keeps the catalogue
   sign at p < 0.05 (`culture_confounded` = flagged and not rescued -> class 6). Module overlap: one-sided
   hypergeometric over the eligible-gene universe, BH over all (module, set) tests.
6. **Hypermutation / TP53 (§8.7).** Hypermutated = MSI covariate level MSI ∪ SIG:TMB ≥ the D 95th percentile (toy:
   the TMB covariate). HOT/ALLELE/LOFBI must keep sign with |β| ≥ 0.5 |β_full| in non-hypermutated lines; SIG must
   keep sign in TP53-mutant lines (`tp53_status = mut`). A refit that cannot be estimated gets `passed = NaN` (not a
   failure).
7. **Pair classes (§12.2).** Base = `audit.pair_classes` (any flagged delta block counts as confounding, as in the toy
   suite). Added class-6 reasons: culture flag, hypermutation, TP53, subtype-confounded (primary-disease delta).
   Order 4 → 5 → 6 → 8 → 7. Class 5 also applies when R power < 0.5 (pre-tag power from
   `configs/catalogue/r_requests.csv`, or from `r_tests.csv`). Class 8 = a catalogue pair matching a CGT test of KO
   class `discordant-mechanism`, with ceiling label data-null / data-contradicted, or with the opposite sign. Class 7
   (R powered, not replicated) overrides 9/8 once R exists. Before R, classes 9/8 (and 6 without power) carry
   `class_status = pending R`. Gene-level class-8 rows (`source = cgt_ceiling`) are added for concordant CGT tests
   whose ceiling is data-null / contradicted (`results/clinical/cgt_ceiling.csv`, from the M6 ceiling script; until
   it exists, class 8 uses KO class and sign only).
8. **Gene classes (§12.1).** Universe = variable ∪ selective (`eligible_genes.csv`) ∪ control dependency genes. D SD
   comes from Part 1 `gene_classes.csv` (same D Chronos effects; avoids reading every gene's effects); r_within =
   `r_beyond`; fidelity = EN-pred `r2_oof_beyond_design`. The M1 named share is applied only when the Tier B file
   `results/attribution/m1_named_share.csv` exists. `failure_class_md_strict` is reported alongside. Pan-lethal =
   Part 1 DepMap `common_essential` flag (Krill-Burger 2023 method).
9. **Extra labels.** processing-limited = a LOSS/AMP/DEL feature of g's own gene (own-arm CN; CYCLOPS).
   predictor-blind = S1-obs BH ∧ not EN-pred selected (∧ replicated in R: pending; `run_xlib` produces no per-pair XLIB
   replication, so XLIB does not enter this label). predictor-only not replicated needs R.
10. **Combinatorial (§11).** Esmaeili Anvar 2024 is the harmonised source for Dede, Parrish, Thompson, CHyMErA and
   Ito; the raw Dede / Thompson supplements are not re-read, so nothing is counted twice. Parrish enters through the
   Esmaeili re-analysis (DESIGN excluded only the unprocessed counts). Screen-level hit = the authors' per-study hit
   list (`Data_for_Figure_01:Figure1_E`); tested = `Figure1_B`; line-level hit = a `Supp_table1` row with CohenD ≥ 0.8
   and deltaMean ≤ −0.5, counted only in studies whose hit list contains the pair (`m5_m8_params`). Flister 2024: only
   the preprint's own screens (NCIH1299 and MDAMB231 columns of media-3), counted as one screen with two lines. Pairs
   are keyed by Entrez (symbols mapped with `genesets.symbol_to_entrez`, aliases included). Candidates: list pairs
   whose feature is a loss family of another gene A (LOFBI, LOSS, DEL, SIL*, METH, YABS, and EXPPAR with β > 0);
   relation paralog (Part 1 top-3 Compara paralogs, either direction) / spec partner / other; validation uses the
   SL-direction paralog candidates. Background = every screened pair with a known identity, minus the list's own
   pairs, restricted to the strata occupied by the list (identity decile of the background × family-size bin
   1/2/3-4/5-9/10+, family size = max of the two genes). Fisher exact OR with exact CI, plus the
   stratum-standardised expected rate. `lookup.csv` holds identifiers and our counts only.
11. **Candidate cards (§12.3).** known_recovery = every positive control (all scopes) and every testable CGT test of
   the headline-eligible set incl. proxy copies (`endpoints.p2_test_set(..., 'incl_proxy')`). Rank = §10.5
   rank_of_feature_D via `m5.RankLookup` (same semantics as `endpoints.RankIndex`, tested equal, without loading the
   arms table). The shortlist is empty with status `pending R` before the single R run; `check_shortlist` asserts
   shortlist ⊆ class 9 ∧ R replicated. not_in_curated_lists: neither (g, d) nor (g, gene of d) is in controls, CGT,
   P1 biomarker_pairs (first word of the biomarker text read as a symbol) or oncogene_pairs.
12. **Memory.** The arms table (~26 M rows) is never loaded whole. Scripts read S1-obs selected rows, the rows of the
   needed genes, and, for the EB fit and the rank lookup, only 2-3 float columns of one arm at a time. `map_pairs`
   applies `analyses.guard_workers` when it is available.

# Streamed (blocked) arms orchestration and the reader API (memory fix before the full launch)

A single in-memory arms table for 2,193 genes is ~26 M rows: ~25 GiB as one parquet-read copy, and > 30 GB peak in
the old `run_arms` parent. `src/biomark/blocked.py:run_streamed` replaces the in-memory assembly.

**Orchestration.** It is used by `scripts/run_arms.py`, the `ctrl` universe and `scripts/run_xlib.py`. Results equal
`arms.run_arms` (`tests/test_blocked.py`, serial and 2 workers).
- Genes are sorted by NaN pattern and cut into chunks at every pattern change and every `chunk_genes` genes.
- One S1 pool scans all chunks. Each chunk's raw scan is written immediately to a stage part
  (`arms.finish_s1_local`), and the parent keeps only slim categorical columns.
- One pool per EN arm does the same. Workers evict FWL caches to one row pattern
  (`ENEngine.max_patterns = runtime_params.blocked.en_max_patterns`; cache only — the frozen `fwl_cache_max_entries`
  stays the upper bound).
- Global steps run on the slim columns with the same rules as `finish_s1` / `finish_en` / `add_dedup` /
  `run_arms`: BH per (arm, family), S1 ranks (p, −|β|, feature; stable), EN ranks (`en.rank_en`), de-duplicated
  ranks, observed cut-offs.
- MD-strict / MD-global / EN-obs-strict use `arms.strict` / `arms.md_global` on the EN-selected candidates and
  their S1-obs rows, read back from the stage parts.
- Each final part is `arms.annotate` plus the global columns, written as `<dataset>/part-NNNNN.parquet` with
  categorical string columns. Each part's sha256 goes to `results/bulk_manifest.csv`.
- Equality: discrete outputs are exactly equal (BH sets, selections, ranks, de-duplicated ranks, MD lists,
  strings). Continuous statistics agree to rtol 1e-9: BLAS/QR rounding differs at the last bit between column blocks,
  so there is no bitwise identity across chunkings.

**Reader API (use this for every arms read; never load the whole table).**
```python
from biomark import endpoint_io as IO
IO.read_arms(columns=[...], arms=None | [arm, ...], genes=None | [entrez, ...], features=None | [id, ...],
             selected_only=False, dataset_dir=None, categorical=False)       # -> pandas DataFrame
IO.pool_inputs(dataset_dir=None, genes=None)    # (ranks of every arm: RANK_COLS, S1-obs rows: S1O_COLS)
```
- Datasets: `IO.ARMS_DIR` = `cache/benchmark/arms` (eligible universe; the default), `IO.arms_dir("ctrl")` =
  `cache/benchmark/arms_ctrl` (U_ctrl; read only by `run_gate.py`, `run_cross_lineage.py` and
  `build_candidate_cards.py`, for controls / CGT targets — never by top-k lists, the null, P1 or P2), and
  `IO.XLIB_ARMS_DIR` = `cache/benchmark/xlib_arms`. `IO.md_dir(u)`, `IO.benchmark_dir(u)`, `IO.universe_genes()`,
  `IO.universe_of(genes)` resolve the per-universe locations at call time.
- Reads go part by part with pyarrow row filters. Categoricals come back as plain strings (shared objects) unless
  `categorical=True`.
- Columns of a part: the `assemble` columns of the former table (arm, gene, feature, S1 statistics incl. df_rule /
  df_nu / p_nu, q, selected, selected_F2, test_family, EN columns, family, d_type, variant_group, feature_stage,
  access_class, en_max_attainable_freq, rank, rank_stat, rank_dedup, lost_to_group_member, se_type).
- P2 ranks: `endpoints.RankIndex.slim(dedup_ranks_by_arm, rows_of_queried_genes)` gives the same `rank_of` as the
  full index (tested).

**Memory guard and worker counts** (`configs/runtime_params.json`, not frozen).
- Defaults: S1 8, EN 6, null 6 (flags `--workers-s1`, `--workers-en`; `run_null --workers`).
- `analyses.guard_workers` lowers a pool's workers to floor(room / per-worker RSS), where room = min(10 GB budget
  − parent RSS, available − 1 GB). Per-worker estimates come from the 100-gene smoke: S1 0.5, EN 0.9, null 1.5
  (unmeasured, conservative).
- The effective count and the guard inputs are written into the timing records (`memory_guard`).
- Null draws process genes by NaN pattern with the same one-pattern cache.
