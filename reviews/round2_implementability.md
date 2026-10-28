# Round 2 review: implementability and agent failure modes (DESIGN v2)

Reviewer: senior research engineer and AI-agent-failure-mode reviewer (independent; did not write DESIGN v2).
What I checked: DESIGN.md v2, configs/{design_params,primary_endpoints,tiers,seeds}.json, src/biomark/*.py, tests/
(31 tests pass in 33 s), docs/data_availability.md, docs/external_data_inventory.md, configs/external_manifest.json,
external/, and small reads of Part 1 configs and schemas.

Facts behind the items below (all checked on 2026-10-28):
- `external/depmap_extra/` is empty. OmicsGlobalSignatures, OmicsSomaticMutations, CRISPRConfounders and OmicsProfiles
  are all absent.
- MSI fallback coverage: Cell Model Passports `msi_status` is known for 441 of 1022 D lines (36 MSI) and 111 of 255 R
  lines. `ModelSubtypeFeatures` contains "MSI" for only 18 lines.
- statsmodels 0.14.4 has no CR2 sandwich. `combine_effects` supports only Paule–Mandel ("iterated") and
  DerSimonian–Laird, not REML. R is not installed, so no clubSandwich or metafor reference exists.
- Candidate genes before the `r_within` filter: 3,062 (variable ∨ selective, tier high or medium).
- Measured: one gene × one EN arm (2 stages × (5 CV paths + 1 full path + 50 fixed-λ fits), K = 310, n = 1000) takes
  2.2 s single-threaded. Residualising the 3,000-column 𝒟 against the unpenalised block anew in each of the 112 fits
  adds about 11 s per gene-arm.
- ModelCondition `CellFormat` levels: adherent, suspension, dome, mixed, spheroid, unknown, NaN.
- Screen-map library strings in D ∪ R: 13 distinct values, including Avana+Humagne-CD (20), Brunello+KY (1) and
  TKOv3 combinations.
- All commits are authored by `ritabaiao`, so git author cannot separate agents.

---

### R2-IM1 (critical): the pre-freeze steps depend on code and data that do not exist until M3b and M4
**Problem.** Three steps are required before `freeze-design`, which tiers.json places in M2 (1.5 h):
- §7.2 needs `simulate_control_bounds.py`, i.e. "the toy null's 95th percentile at real D dimensions". That requires
  the complete pipeline (arms, stability, D-stage, BH) and the toy generator, both scheduled for M3b and M4.
- §9.3 needs `pc_eligibility.csv`, which requires the built dictionary (M3b).
- §7.1 and §13.4 put `configs/null_genes.txt` among the frozen hashes, but it requires the eligible gene list (M3b).

An agent at M2 must either freeze provisional bounds and quietly "replace" them later, or skip the freeze. Both
outcomes are unlogged deviations of exactly the kind §13.4 forbids.

**Change.** Use two tags:
- `freeze-design`, now: design_params (with NC bounds marked `provisional`), primary_endpoints, tiers, seeds, the
  controls *specification* (R2-IM5), culture and lineage lists, and the toy plants with `toy_seed`.
- `freeze-bounds`, at the end of M3b: NC1–NC6 bounds from the toy simulation, `pc_eligibility.csv`,
  `null_genes.txt`, `dictionary_v1.json` and `data_status.json`.

`freeze-bounds` must precede the first commit that writes any of `results/benchmark/`, `results/md/`,
`results/controls/G_*` or `results/null/`. Update tiers.json: the M2 content line drops "sim bounds, pc_eligibility",
and M3b ends with `freeze-bounds`.

**Check.** `gate.py frozen`:
- for each of those result paths, `git merge-base --is-ancestor freeze-bounds <first commit touching path>`;
- `configs/frozen_hashes.json` has two sections (`design`, `bounds`), and no key under `negative_controls_s7` still
  carries `provisional` at `freeze-bounds`.

### R2-IM2 (critical): R can be opened only once, but the R tests needed downstream are not listed
**Problem.** Guard condition (h) allows only reproduction after the first R access (same config hash and code
commit). Many deliverables need R statistics, and v2 does not list them as `replicate.py` outputs:
- §8.3 Q1b: R replication of cross-lineage vs lineage-restricted MD-strict determinants, which are not necessarily in
  any top-100.
- §12: class 5 (R-underpowered) and class 7 (powered, not replicated) for every MD-strict pair, and pair class
  "predictor-blind", which needs R tests on all S1-obs BH hits.
- §6.2: the R non-Avana subset, gene-matched lists and k = 500.
- Tier B LMM: §6.2 does not say whether it runs on D or R.
- §2: "R accurate-dependent" is defined but nothing consumes it.

An agent at M8 that finds a missing R number can only drop the analysis silently or break the guard (new commit plus
new access). A Tier B item finished after the tag is locked out permanently.

**Change.**
- At `freeze-catalogue`, commit `configs/catalogue/r_requests.csv` (g, d, reason ∈ {P1 list k ∈ {50, 100, 200, 500}
  per arm and variant, MD-strict, EN-obs-strict, S1-obs BH, cross-lineage or restricted label, control, CGT
  feature}).
- `replicate.py` computes, in one run, R β, SE, p, n, carriers and power for every row, plus R accurate-dependent
  membership for MD-strict subtypes. It writes `results/replication/r_tests.csv` and `r_lines.csv`.
- Every later R-dependent script reads only these files.
- Declare LMM a D-only robustness check, so it needs no R access.
- Declare that set A (scope `limits`) is read by `replicate.py` in the same run.
- The `--pseudo-R fold0` dry run must produce every downstream table that consumes R (Q1b, pair_classes,
  gene_classes classes 5 and 7, the R non-Avana subset) from `results/replication_dryrun/`.

**Check.**
- `gate.py G_rep`: every (g, d) in `results/limits/pair_classes.csv` with class 5 or 7, every Q1b row and every P1
  list row is present in `r_requests.csv`.
- No file under `scripts/` other than `replicate.py` calls `guard.load_effects(scope="replication"|"limits")`
  (AST scan).
- The dry run's output directory contains the same file set as the planned real outputs.

### R2-IM3 (critical): essential inputs are missing, and the data deadline falls after the step that uses them
**Problem.** `external/depmap_extra/` is empty. §1.5 sets the data deadline at the start of M4, but `build_dictionary.py`
(M3b) and `pc_eligibility` (before freeze) already branch on these files, through WGD, LoH, ploidy, SPEC_*, TMB, MSI
and alleles. A file that arrives mid-M3b would force an unlogged rebuild or an ad hoc mixed state.

Without OmicsSomaticMutations, every allele-level CGT test is `untestable: allele`: BRAF V600, KRAS G12C, EGFR
activating, PIK3CA, ERBB2 TKD, AKT1 E17K, IDH1 R132, KIT and FGFR3. That is most of the S tests behind P2. No rule
says what happens to P2 when only a handful of tests remain testable. §10.3 also leaves it undefined whether
"RAS-mutant" (the EGFR W test) and "damaging + LoH" (BRCA) are allele-level or can be computed from local files.
MSI from fallbacks is NA for about 56% of D lines, which affects the WRN←SIG:MSI control and the MSI covariate.

**Change.**
- Make `data_status.json` an input to `build_dictionary.py`: the deadline becomes the first dictionary build, and its
  sha256 is stored in `dictionary_v1.json`. Late files feed only labelled secondaries, never a rebuild.
- Ask the user now for the four Priority-1 files, which are the project's critical path.
- In `build_cgt_features.py`, pre-declare:
  - RAS-mutant = HOT:KRAS ∨ HOT:NRAS, labelled `gene-level proxy`;
  - LoH = own CN < 0.75 when GlobalSignatures is absent, labelled proxy.
- Add `P2_primary_clinical.min_testable_tests = 15`. Below it, P2 is reported as "descriptive (underpowered: n
  testable tests)", with a secondary recall in which allele tests use gene-level HOT, labelled `HOT proxy` and never
  pooled with the headline.

**Check.**
- `gate.py data`: `dictionary_v1.json.data_status_sha256` equals the sha256 of the committed `data_status.json`, and
  `data_status.json` was committed before `dictionary_v1.json`.
- `q4_summary.json` has `n_testable` and `p2_status ∈ {primary, descriptive}`, recomputed by script from
  `min_testable_tests`.

### R2-IM4 (major): EN-pred's blindness, its row sets, subsamples and cost are underspecified
**Problems.**
- (a) §1.6 says "each (g, d) uses lines with non-NaN effect", but EN-pred runs with `load_effects` patched to raise.
  It cannot learn the NaN mask, so it must either fit all 1022 lines (rows differ from EN-obs, and the arms are no
  longer paired) or break blindness.
- (b) Seeds conflict. `subsample_design` uses `subsample_seed`, while `compute_s14.gene_seed` is "seed_list[i] XOR
  entrez" with i undefined. It is unclear whether subsamples are shared across genes and arms.
- (c) The FWL residualisation of 𝒟 is specified per fit. Recomputed per gene, that costs about 13 s instead of 2.2 s
  per gene-arm. That is roughly 1.8 h for the main arms, 1.5 h for XLIB and 3.6 h for null B = 20, against a 180-min
  cap. Shared across genes, those figures fall to about 15 min, 10 min and 40 min. The pilot rule changes only K, so
  it would cut K when the real cost is residualisation.
- (d) Unspecified: the stage-X target inside CV folds and subsamples (the subsample's own stage-G residual or the
  full-D one?).
- (e) Unspecified: under §7.1 permutation, whether `CN_absent` and `MSI_unknown` move with the 𝒟 rows. If they stay
  fixed, NA → 0 cells land on lines that have data.

**Change.**
- (a) A separate process writes `cache/effect_nan_mask.parquet` (boolean) and its sha256 to `results/bulk_manifest.csv`.
  This is the only effect-derived input EN-pred may read, and both EN arms use identical row sets.
- (b) Draw the 25 complementary pairs once from `subsample_seed` over D lines, patient-grouped and lineage-stratified,
  into `configs/subsamples.json`, shared by all genes and arms. Each fit uses subsample ∩ non-NaN rows. Delete
  `gene_seed`; coordinate descent is deterministic.
- (c) Cache the FWL-residualised 𝒟 per (NaN-pattern group × fold or subsample), memmapped and read-only. The pilot
  times 50 genes with this cache.
- (d) Every fit uses its own stage-G residual.
- (e) NA indicators are 𝒟 attributes and are permuted with the 𝒟 rows; covariates are not.

**Check.**
- `tests/test_en_pred_blind.py` also asserts that the row sets of EN-pred and EN-obs are equal per gene.
- `tests/test_subsamples.py` checks the manifest hash equality across arms.
- A unit test confirms the permuted `CN_absent` stays aligned with the permuted CN rows.
- `results/timing/surrogate_pilot.json` reports `fwl_cache_hit_rate`.

### R2-IM5 (major): the leakscan rule forbids the scripts the design requires, and "decision literal" is undefined
**Problem.** §13.4 bans gene symbols and `ACH-` ids as literals in `scripts/`, and bans "decision literals" in `src/`,
but several required scripts would contain them:
- `build_positive_controls.py` is "built from this section", which means hard-coding 60+ symbol pairs;
- `build_cgt_features.py` hard-codes allele lists (BRAF V600, EGFR set, KIT D816, …);
- the known-positive fusion lines are written as line names;
- `build_cgt.py` hard-codes corrections (KIT D816V, del(17p)).

An agent will either hit leakscan failures or weaken leakscan. "Decision literal" has no operational definition, so
no script can check it.

**Change.**
- Move the lists into frozen configs: `configs/positive_controls_spec.csv` (the §9.2 rows), `configs/cgt_allele_rules.json`,
  `configs/fusion_known_positives.csv` (cell-line names resolved to ModelID by Model.csv), and
  `configs/cgt_corrections.csv`. The build scripts only read and join them.
- Define leakscan as an AST scan of `src/` and `scripts/`:
  - (i) string constants matching an HGNC symbol of 3 or more characters from `gene_map`, or matching
    `ACH-\d{6}`;
  - (ii) numeric constants outside {0, 1, -1, 2, 0.5, 100} in `src/biomark/` other than `params.py`, unless the line
    carries `# literal-ok: <reason>`.

**Check.**
- `gate.py leakscan` exits 0 on HEAD.
- `tests/test_leakscan.py` plants a file with `"KRAS"` and `0.7` and expects two hits.
- `frozen_hashes.json` includes the four new config files.

### R2-IM6 (major): undefined terms and thresholds that force unlogged choices
**Problem.** Each of the following leaves an agent to choose:
- §8.1: "surrogate deviation" is undefined. Which stages count, residualised or raw features, and selected or all
  non-zero β?
- §8.5: the α for "(ii) not significantly negative" is not given, and neither is the meaning of "g is depleted".
- §8.6: the Cramér's V threshold that flags a subtype as culture-confounded is missing, so the class-6 culture rule
  has no input.
- §12: it is unclear how pair-level classes 4–8 become a gene class. Eligible genes with fidelity ≥ 0.10 and no
  MD-strict pair have no class. The §12 gene universe and the class-5 "too rare" cut are not given.
- §5.3: the CellFormat → {adherent, suspension, 3D, mixed/unknown} map is missing. Data have dome, spheroid, mixed,
  unknown and NaN. The "formulation regex" and the medium-family map are also missing.
- §3: the "partner priors" list for LOSS scope and the kinase list for FUS orientation are not given.
- §6.1: the SE_R formula is missing, and so is the expected non-NaN rule for R library strings absent from D or mixed
  (13 strings exist).
- §6.3: Brunello+KY is neither Avana-only nor KY-only.
- §13.1: "different agent by git author" cannot be checked, because every commit has the same author.
- §1.3: "catalogue agents do not read P1 REPORT.md test sections" cannot be scripted.

**Change.** Add keys, with the proposed values below, to design_params before `freeze-design`:
- `subtypes_s8.line_share_def` = "β_gd·(x̃_ℓd − x̄_d) / Σ_{selected d'} β_gd'·(x̃_ℓd' − x̄_d') over both stages,
  FWL-residualised features; only lines with negative sum";
- `cutting.alpha` = 0.05; `cutting.depleted` = "carrier median effect of g < −0.5";
- `culture_flag_cramers_v` = 0.3;
- `limits_s12`:
  - universe = variable ∪ selective genes ∪ controls;
  - gene class = minimum over its pairs' classes under the precedence 4 < 5 < 6 < 7 < 8 < 9 (class 9 if any pair is
    valid);
  - new class `3b` "predictable, nameable, no MD-strict determinant";
  - class 5 = < 5 D carriers or R power < 0.5;
- `culture_maps` = {dome, spheroid → 3D; mixed, unknown, NaN → mixed/unknown}, with an explicit regex table for the
  medium family;
- `fusion.orientation` = "driver from fusion_drivers.csv is the 3′ gene (5′ for FGFR2/3)", which removes the need for
  a kinase list;
- `partner_priors` = the partner column of `positive_controls_spec.csv`;
- `replication_s6.se_r` = "σ̂_resid,D(g) · sqrt(1/n1_R + 1/n0_R)";
- `R_nan_rule` = "gene expected non-NaN if ≥ 95% of D lines with the same library string are non-NaN; strings absent
  in D use the union NaN set of their component libraries";
- `xlib_replication` = "no Avana in library string".

Declare as **procedural, not script-checked**: the P1 REPORT.md/RESEARCH_TREE.md reading rule, and toy authorship
(the `Agent: <name>` trailer is required; checked by trailer, not author).

**Check.** `gate.py params` asserts every listed key exists and is non-null. `tests/test_culture_map.py` asserts every
observed CellFormat level maps to exactly one class.

### R2-IM7 (major): hand-written estimators have no reference implementation
**Problem.** The following estimators must be written by hand and cannot be tested against R:
- CR2, with its df convention (Bell–McCaffrey or n − p) unstated;
- REML τ² with HKSJ;
- the EB normal scale mixture;
- Deming regression;
- the FaST-LMM.

`test_s1_vectorised` compares only CR1 against statsmodels. The patient bootstrap of R has further gaps:
- it does not say how resampled duplicate patients are clustered;
- it does not say whether covariate drops (|ρ| > 0.7) are re-derived per draw;
- it does not say how empty covariate levels are handled.

§13.3 asks for 1e-6 agreement on EN selection from an independent reimplementation, which is impossible with
coordinate descent at tol 1e-3, unless the reimplementation reuses sklearn, in which case it is not independent.

**Change.**
- CR2 with residual df = n_clusters − 1 (t). Test: with singleton clusters it equals statsmodels HC2 to 1e-10.
- Replace REML with Paule–Mandel (`combine_effects(method_re="iterated")`) plus the HKSJ variance from the same
  function, so no hand-written τ² estimator is needed.
- EB mixture: toy known answer (simulate from a 3-component scale mixture, recover the posterior means; r > 0.99 vs
  the oracle).
- Deming: use `scipy.odr` with sx = SE_D and sy = SE_R.
- Bootstrap:
  - each resampled patient copy gets a new cluster id;
  - drops and pool membership are fixed from the full-R fit;
  - zero-variance columns are dropped by pivoted QR and logged.
- §13.3 tolerances: S1 p within 1e-6, counts exact, EN selection frequency within 2/50 with an identical selected set
  excluding features within 0.04 of the 0.6 cut, and β within 1e-3 relative.

**Check.**
- `tests/test_cr2.py`, `tests/test_eb_mixture.py` and `tests/test_bootstrap_clusters.py` exist and pass.
- `gate.py G_indep` reads the tolerances from design_params, not from code.

### R2-IM8 (major): the schedule is over-committed and late slippage has no cut trigger
**Problem.**
- About 18 h remain, not the 19 in tiers.json, so the reserve is about 1.5 h.
- XLIB re-runs all four arms on about 620 Avana-only lines, but only `replication_xlib` (20 min) is budgeted, inside
  the 1 h M7.
- The null (12 workers, background) overlaps the M4 foreground arms (12 workers). The concurrency rule covers only the
  Tier B refit, so both jobs slow by about 2× without a log entry.
- The cut order triggers only "if M4 ends > 1 h late", so slippage in M5–M8 has no rule.
- The independent reimplementation sits in M7, on the critical path, although its inputs exist after gate attempt 1.
- `budget_check.py` stage names have no mapping to job names in `results/timing/`.

**Critical path:**
1. User downloads.
2. Dictionary, covariates and eligibility (M3b, about 1.5 h).
3. Toy and bounds (`freeze-bounds`).
4. S1 and EN arms.
5. MD-strict, then gate attempt 1.
6. Subtypes and cross-lineage.
7. CGT scoring, which needs `freeze-cgt`.
8. `freeze-catalogue`, then the single R run.
9. Limits, claims and REPORT.

The toy, the CGT curation, the null and the independent check are parallel branches that must not take cores from
this path.

**Change.**
- Set `hours_remaining_at_v2: 18.0` and `reserve_h: 1.5`.
- Add the Tier A stage `xlib_arms: 30` in M4 (background, after the main arms) and move it out of M7.
- Allocate cores: null 4 workers, foreground 8 (`compute_s14.null_workers`).
- Move the independent check (13.3) to start right after `G_pc_attempt1.json`, in parallel with M5 and M6.
- Run `budget_check.py` at every milestone end. If cumulative lateness exceeds 1 h, apply the next unused entry of
  `cut_order` and log it to `results/deviations.csv`.
- Add `tiers.json.stage_job_map` (stage → timing file stems).
- Prioritise `gate.py` subcommands: leakscan, frozen, data, controls and G_pc in M3b; the rest by M7.

**Check.**
- `budget_check.py` exits non-zero if any `results/timing/*.json` stem is unmapped, if the sum of concurrent workers
  in the timing logs exceeds 12, or if a milestone end exists without a matching `deviations.csv` decision when
  lateness > 1 h.
