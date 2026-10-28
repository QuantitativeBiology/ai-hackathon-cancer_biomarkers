# Review round 1 — statistics (reviewer: statistician)

Scope: selective inference, stability selection, surrogate specification, confounder model, collinearity,
non-independence, permutation null, BH families, cross-lineage rule, replication power, and whether MD vs S1 is a fair
comparison. Reviewed DESIGN.md v1 against PROMPT.md, docs/part1_brief.md and docs/lit/*.md. I made small read-only
reads in Part 1 (configs only). Facts used below:
- **Eligible set.** 3,062 eligible genes (1,331 high, 1,731 medium), from `P1/configs/gene_classes.csv` ×
  `predictable_genes.csv` with the §2 rule.
- **Lineages.** Dev has 6 lineages with < 10 lines (Thyroid 9, Ampulla 8, Eye 6, Testis 3, Vulva 1, Adrenal 1).
- **Library × lineage.** Library is partly nested in lineage: Bladder, Kidney, Skin and Lymphoid are about 100%
  Avana-only, and Humagne-CD is concentrated in Bowel, CNS, Esophagus/Stomach and Pancreas.
- **R lineages.** The R lineage mix is roughly proportional to D (CNS 24, Lung 21, Ovary 18, ...).

Item format: ID · severity · problem · change (text or rule for DESIGN.md) · script check.

---

### R1-ST1 · critical · MD's D-stage test is double dipping; its BH q-values are not valid

**Problem.** MD step (ii) (§5.1) tests d against the observed effect y_g in D after d was selected on ŷ_g (§4.1).
Each OOF prediction of a line in fold j is a function of X and y_{−j}, so the predictions of all D lines together use
every y in D. Noise in y that happens to correlate with d is learnt into ŷ, which makes d more likely to be selected.
That same noise is then scored again in the observed-effect test. The test is conditional on selection, so its
p-values are anti-conservative. BH "over all MD candidates" is applied after selection, so the q-values carry no FDR
guarantee. This also biases the benchmark in MD's favour (see ST4).

**Change (§5.1, new paragraph "Inferential status").**
- "The D-stage association in MD is a screening filter, not a confirmatory test. Its p-values are reported as
  *selection-conditional (nominal)*. The FDR attached to MD in D is the **empirical end-to-end FDR** from the null in
  ST2. Confirmatory claims about MD candidates come only from R (§6) and the external sources (§7)."
- Optional cross-fit sensitivity, if time allows: for each fold k, select on the fold-k lines' ŷ only and test on the
  fold-k lines' y only. Those predictions do not use y_k, and y_k ⟂ y_{−k} given X if lines are independent. Combine
  the five fold-wise tests by Stouffer's method for the (g, d) pairs selected in ≥ 3 folds.
- Wording rule: REPORT.md must not use "significant" for MD D-stage results without the qualifier.

**Check.** A grep test fails if `results/md/*.csv` has a column named `q_bh` without the companion columns
`p_status == "selection_conditional"` and `fdr_empirical`. A REPORT lint flags "significant" near "MD" unless a
q-value from R or an empirical FDR is cited.

### R1-ST2 · critical · The permutation null as specified does not reproduce the pipeline and is not affordable

**Problem 1: the null is too easy.** §5.4 permutes each d within lineage. That breaks the correlations among 𝒟
features, which the elastic net and the locus collapsing depend on. It also breaks d's correlation with library and
culture, and §7.2 says something different ("permutation of 𝒟").

**Problem 2: the selection step is not re-run.** If the null permutes d but does not re-run stability selection, it
cannot measure the selection-induced inflation (ST1).

**Problem 3: compute.** 200 permutations × 3,062 genes × 50 subsamples of the full elastic net does not fit in the
budget. A 1022 × several-thousand-feature EN fit costs about 1 s, so one pass is about 150 k fits ≈ 2 h wall on
12 cores, and one pass per permutation is impossible.

**Change (§5.4 multiple testing and §7.2, unified).**
- **Permutation unit.** Permute **rows of the whole 𝒟 block jointly** within lineage × library-group strata. Strata
  with fewer than 2 lines are merged into the lineage's largest stratum. y, ŷ and the covariates stay fixed. This
  keeps the feature–feature correlations and locus structure. The null is "no molecular determinant beyond
  covariates".
- **What is re-run.** Every permutation re-runs the full MD pipeline: SIS screen (ST6), stability selection, the
  D-stage test and the BH threshold. It also re-runs S1 and EN-obs (ST4) with their own thresholds.
- **Size.** Use a pre-declared random subset of 300 eligible genes (seed from `configs/seeds.json`, drawn before any
  result) and B = 50 permutations. If the pilot timing (ST6) allows, raise B to 200 on the same subset. Do not adapt B
  after seeing results.
- **Estimator.**
  - FDR_emp(list) = mean_b |hits_b| / |hits_obs| on the 300-gene subset, with a bootstrap CI over permutations.
  - Also report P(|hits_b| ≥ 1).
  - Apply the subset's rate to the full list as an *estimate* and say that it is one.
- **Pass rule (pre-declared).** For each list, FDR_emp ≤ 0.10. If FDR_emp > 0.10, the list's D-stage threshold is
  reported as miscalibrated, not retuned. Retuning would change a test after seeing results, which PROMPT.md forbids.

**Check.**
- `tests/test_permutation_null.py` asserts three things: the permuted 𝒟 has the same feature–feature correlation
  matrix as the original (up to row order); within each stratum, the multiset of rows is unchanged; and the pipeline
  function called under permutation is the same function object as the one used for the observed run.
- `results/null/summary.json` must contain `B`, `gene_subset_sha`, `fdr_emp` and `ci` for MD, S1 and EN-obs.

### R1-ST3 · critical · The confounder model adjusts away the determinant for several positive controls

**Problem.** §5.3 always includes "own-gene CN of g" and "TMB/MSI".
- When d is AMP or DEL at g's own locus (ERBB2, MYCN, EGFR amp; in effect also MTAP/9p21 for co-located genes), the
  own-CN covariate is close to collinear with d, so β_d collapses. Three of the §7.1 gate controls are therefore
  adjusted away by design.
- When d is MSI (WRN–MSI), the MSI covariate is d itself. TMB is also strongly correlated with MSI and with LOF calls
  in hypermutated lines.
- The same applies to EXP-own, whose correlation with own CN can be high.

**Change (§5.3, rule).**
- "For each (g, d), drop from the covariate set any covariate that (a) is d or a component of d's locus feature, or
  (b) has |Spearman r| > 0.7 with d in the test lines. Record the dropped covariates per test in a `dropped_cov`
  column."
- "For SIG determinants, TMB and MSI are not both covariates. MSI as d is adjusted for TMB only if |r| ≤ 0.7, and vice
  versa."
- "For AMP/DEL determinants that contain g, the own-CN-adjusted estimate is a *sensitivity* result. The primary
  estimate omits own CN. The Cas9-cutting rule (§5.4) handles the artefact instead."
- Apply the same rule to MD, S1 and EN-obs.

**Check.** A unit test fits the confounder model for ERBB2 × ERBB2-AMP and WRN × MSI on D. It asserts that
`dropped_cov` contains own-CN and MSI respectively, and that neither design matrix has a column with |r| > 0.7 to d.
A script scans all fitted tests and asserts max |r(d, covariate)| ≤ 0.7.

### R1-ST4 · major · MD vs S1 is not a like-for-like comparison

**Problem.** MD gets two advantages:
- an extra filter (stability selection) before the observed-effect test;
- a BH family of a few hundred selected candidates, compared with S1's millions of tests (3,062 genes × |𝒟|, and
  |𝒟| is likely in the thousands, since SIL alone can contribute thousands).

**What "MD > S1" could mean.** It could reflect any of four things: (1) the predictor's denoising, (2) multivariate
sparse selection, (3) the stability filter, or (4) a lighter multiplicity correction. Q2 asks about (1). The S1-pred
ablation separates (1) only for univariate tests.

**Other asymmetries.**
- The test statistic differs: S1 uses an LRT, and MD's test is unspecified.
- The LMM sensitivity is run for S1 only.
- Matched-size lists need a ranking for MD that is not defined.

**Change (§6).**
- Add the competitor **EN-obs**: the identical stability-selection pipeline (same 𝒟, SIS screen, λ rule, subsample
  count and seeds) fitted on the **observed** y_g, followed by the same D-stage test. The primary Q2 contrast is
  **MD vs EN-obs**, which isolates the predictor. MD vs S1 is secondary, and S1 vs S1-pred vs EN-obs vs MD is
  reported as an ablation ladder.
- **One test for all lists:** OLS LRT for d in the §5.3 model (after the ST3 drops), with the same covariates. Run the
  LMM sensitivity on the top 1,000 hits of **each** list.
- **Pre-declared rankings for matched-size lists:**
  - MD and EN-obs: selection frequency, ties broken by D-stage p;
  - S1 and S1-pred: p-value.
  - Use k ∈ {50, 100, 200, 500} ∩ [1, min list size]. The matched-k comparison is primary for replication and CGT
    recall; the threshold-based lists are secondary.
- **Gene-matched variant.** Strong genes (lineage TFs) can dominate a top-k by p. Also compare at most m = 2 features
  per gene, the top m by each list's own ranking.
- **Overlap strata.** Report replication separately for MD∩S1, MD-only and S1-only hits. The difference between lists
  is identified only by the non-shared strata.
- **CGT recall** is reported as recall-at-k curves (k as above) with an exact McNemar test on paired pair-level
  outcomes, not as raw list recall. Larger lists recall more trivially.

**Check.**
- `results/benchmark/lists.parquet` must hold one row per (list, g, d) with `rank_stat` and `rank`. A script asserts
  that all four lists use the same covariate-set hash, the same test function, and the same subsample count and seed
  list (read from run logs).
- `results/benchmark/replication_by_k.csv` has all lists at every feasible k.
- `results/benchmark/overlap_strata.csv` exists.

### R1-ST5 · major · Replication in 255 lines is mostly underpowered, and power from D effects is inflated

**Problem: few carriers.** A feature with 5 dev carriers has about 1.25 expected R carriers, so P(≥ 3 carriers in R)
≈ 0.13. A one-sided α = 0.05 test with 80% power needs Δ ≈ 1.43 σ at 3 carriers and ≈ 0.79 σ at 10 carriers.

**Problem: winner's curse.** D effects of selected hits are biased upward. The bias is larger for lists with stricter
or double-dipped selection (MD; S1 at BH over millions). Power computed from D effects is therefore optimistic, and
optimistic to different degrees across lists.

**Problem: hit-level counts.** Hit-level "same sign and p < 0.05" fractions over mostly underpowered hits will be
dominated by noise.

**Change (§6 Power and comparison 1).**
- **Evaluable in R:** ≥ 3 carriers and ≥ 3 non-carriers in R after the §5.3 covariate handling. Non-evaluable hits are
  counted and reported, not dropped silently.
- **Power from a shrunk effect.** Use β_D shrunk by an empirical-Bayes normal-mixture prior fitted to all tests of that
  list's family (e.g. ashr-style; or simply β_D·max(0, 1 − SE²/β_D²)). The same method is used for all lists.
- **List-level primary statistics,** each with a bootstrap CI over R lines (resample R lines, refit, 1,000 reps, fixed
  seed):
  - (a) sign concordance among evaluable hits, binomial test against 0.5;
  - (b) replication slope, the Deming or OLS slope of β_R on β_D;
  - (c) observed/expected replications, with the expected count from the shrunk power.
- **MD vs EN-obs vs S1 differences:** a paired bootstrap over R lines, where both lists are re-evaluated on the same
  resample. Report the difference with its CI; no claim without a CI that excludes 0.

**Check.** `results/replication/summary.json` must contain the fields `n_evaluable`, `n_nonevaluable`,
`sign_concordance`, `slope`, `obs_over_exp` and `ci_*` per list, plus `delta_*` with a CI for each pairwise
comparison. A script asserts that the power was computed with the declared shrinkage function for every list (one
function hash) and before the `freeze-catalogue` tag (file mtime and git log).

### R1-ST6 · major · Stability selection as specified gives no error control; the surrogate is under-specified

**Problem 1: no error bound.** Choosing λ by CV/1-SE *within* each half-sample voids the Meinshausen–Bühlmann /
Shah–Samworth bound, so the 60% threshold has no error interpretation.

**Problem 2: unspecified details.**
- Standardisation of rare binary features is not stated.
- Fidelity "R² beyond lineage" is in-sample. With thousands of features it is inflated, and it drives failure class 3
  in §8.

**Problem 3: compute.** Feasibility is unknown (see ST2).

**Problem 4: covariate asymmetry.** The surrogate has lineage unpenalised but not library or culture. PROG or SIL
features that proxy organoid or library status can be selected as "determinants" and only partly caught later.

**Change (§4.1).**
- **Complementary-pairs stability selection** (Shah & Samworth 2013): 25 complementary half-sample pairs (= 50
  subsamples), patient-grouped and lineage-stratified.
- **Selection size fixed per gene by a q-rule:** in each subsample take the first q = 20 features to enter the EN
  path. Selected = frequency ≥ π = 0.6 with a consistent sign. Report the bound E[V] ≤ q² / ((2π − 1)·p) per gene with
  p = |𝒟_g|. The CV/1-SE λ is used only for the fidelity fit.
- **SIS pre-screen inside each subsample:** the top 300 𝒟 features by |partial correlation| with the target given the
  unpenalised covariates. This makes 50 × 3,062 × 3 pipelines (MD, EN-obs, permutations) feasible. Run a timed pilot
  on 100 genes before M4 and record the result. The fallback to 20 subsamples, if needed, applies to **all** EN lists
  equally.
- **Unpenalised covariates:** lineage + library group + growth pattern (the same block as §5.3, minus gene-specific
  terms).
- **Standardisation:** all features z-scored on D. Binary features with < 1% prevalence are flagged in the output, and
  rare features have low selection stability by construction (≈ 2–3 carriers per half). A selection-frequency vs
  carrier-count plot is reported.
- **Fidelity = cross-validated partial R²** (patient-grouped 5-fold) of the surrogate beyond the unpenalised
  covariates. The class-3 threshold of 0.1 applies to this OOS value.

**Check.**
- `results/surrogate/selection.parquet` contains `q`, `pi`, `ev_bound` and `n_subsamples`.
- A test asserts that every subsample's selected set has size ≤ q.
- `fidelity_oos` exists, and §8 class assignment reads `fidelity_oos`, not the in-sample R² (grep the limits script).
- `results/timing/pilot_surrogate.json` exists with a commit earlier than the first M4 result.

### R1-ST7 · major · Collinear covariates and rare lineages make the confounder model and the "confounding delta" unstable

**Problem 1: screen length.** Screen length is nearly a deterministic function of library group (median 21 / 24 / 14
days for Avana / Humagne / KY).

**Problem 2: growth pattern.** Suspension ≈ haematopoietic lineages and Dome ≡ organoid, which is concentrated in a
few lineages and library groups. Library is fully nested in lineage for Bladder, Kidney, Skin and Lymphoid.

**Problem 3: rare lineages.** Lineage fixed effects for the 1–9-line lineages consume degrees of freedom. A lineage
with one line has its line fitted perfectly, so that line is uninformative and silently dropped.

**Problem 4: the delta flag.** "|β| drops > 50%" is unstable when β is small. "With and without each block" depends on
the order in which blocks are added.

**Problem 5: within-lineage subtypes.** These are not adjusted for (e.g. luminal vs basal breast, B vs T lymphoid).

**Change (§5.3).**
- **Screen length:** enter as the within-library deviation in days (days − library median). If its VIF is still > 10,
  drop it and state that library group absorbs it.
- **Growth pattern:** collapse to {adherent, suspension, 3D (dome/spheroid/organoid), mixed/unknown}.
- **Lineage:** lineages with < 10 D lines are pooled into `_rare` for the covariate (the Part 1 convention).
- **Gate on the design:** before any test, compute the condition number and VIFs of the covariate design. Write
  `results/covariates/vif.csv`. Any covariate with VIF > 10 is resolved by the declared rule above, not ad hoc.
- **Confounding delta = drop-one-block from the full model.** Flag a pair if (|β_full − β_minus| > 0.5·|β_minus|) and
  (|β_full − β_minus| > 1·SE_minus). Blocks: lineage, library + length, culture, burden (TMB/MSI), own-CN.
- **Sensitivity:** add OncotreePrimaryDisease (diseases with ≥ 5 D lines; others take the lineage level) as an extra
  block. A hit whose β flags on this block is labelled "subtype-confounded" in §8 class 6.

**Check.**
- `vif.csv` exists and its max VIF is ≤ 10 in the final design.
- A test asserts that no lineage level in the design matrix has fewer than 10 lines.
- `results/confounders/delta.csv` has one row per (g, d, block) with `beta_full, beta_minus, se_minus, flag`. A script
  recomputes the flag rule from those columns and asserts that it matches.

### R1-ST8 · major · The cross-lineage rule is undefined at key points and is not calibrated

**Problem: undefined terms.**
- "Meta effect significant" has no α or family.
- DerSimonian–Laird τ² is unstable with strata of 3 carriers.

**Problem: the "≥ 2 lineages nominal p < 0.05" rule is not calibrated.**
- It favours features present in many lineages. Under the null, with 10 testable lineages, P(≥ 2 nominal hits) ≈
  0.09.
- Absence of significance elsewhere is labelled "lineage-restricted" even when it is only lack of power or lack of
  carriers (MYCN AMP is almost all neuroblastoma).
- Leave-one-lineage-out drops only the largest lineage.

**Change (§5.2 within vs across lineage).**
- **Meta-analysis method:** REML τ² with a Hartung–Knapp–Sidik–Jonkman CI. Report I², τ and the 95% prediction
  interval.
- **Labels:**
  - **cross-lineage:** (a) meta p passes BH at 0.05 over all (g, d) entering the cross-lineage stage (its own family);
    (b) same sign with one-sided p < 0.05 in ≥ 2 lineages; (c) the sign is kept under leave-one-lineage-out for
    **every** contributing lineage.
  - **lineage-restricted:** the per-lineage effect passes BH within its lineage family **and** the d × lineage
    interaction LRT passes p < 0.05.
  - **single-lineage carriers (untestable):** < 2 lineages meet the ≥ 3/≥ 3 carrier rule.
  - **undetermined:** everything else (replaces "pooled-only").
- **Calibration:** apply the same labelling to the ST2 permutations and report the false cross-lineage rate per list.

**Check.** `results/subtypes/cross_lineage.csv` has the columns `tau2_method == "REML"`, `ci_method == "HKSJ"`,
`n_lineages_testable`, `pi_lo`, `pi_hi`, `label`, `lolo_all_sign_kept`. A script recomputes each label from the stored
per-lineage statistics. `results/null/summary.json` includes `false_cross_lineage_rate`.

### R1-ST9 · major · BH families and test counts are ambiguous

**Problem.**
- "Each at its full test count" is ambiguous for MD (post-selection count vs |G|·|𝒟|).
- EXP-own tests (own expression vs own dependency) are numerous, near-tautological and dominate S1's discoveries,
  which distorts both burden and recall.
- The cross-lineage stage, the per-lineage tests, the CGT and the positive controls have no stated families.

**Change (§5.4 multiple testing): a family table**, written to `configs/test_families.json` before any test is run.

| Family | Rule |
|---|---|
| F1 S1-main | All (g, d) with d ∉ {EXP-own}. BH 0.05 at the count of tests actually run. |
| F2 S1-EXP-own | Its own family, reported separately. Never counted toward recall of non-expression biomarkers. |
| F3 MD/EN-obs | BH over the selected candidates. Labelled screening-only (ST1); empirical FDR from ST2 is the reported error rate. |
| F4 cross-lineage meta | Its own family (ST8). |
| F5 per-lineage tests | Within each lineage family. |
| F6 R replication | No BH. List-level statistics (ST5); hit-level one-sided p is reported as descriptive. |

- BH is valid under positive dependence (PRDS), which is plausible here. Storey q is reported as a sensitivity only.
- The test count n_F used for each family is logged and must equal the number of rows in that family's result table.

**Check.** A script loads `test_families.json`. For each family it recomputes the BH q from the stored p-values and
asserts that it equals the stored q, and that n_F equals the row count. It also asserts that the commit of
`test_families.json` precedes the first result commit.

### R1-ST10 · major · Cell lines are not independent; module counts double-count correlated genes

**Problem.**
- 104 screened lines share a patient, and OLS treats them as independent.
- Residual line similarity beyond lineage (Wainberg; Gonçalves) inflates OLS p-values.
- The LMM sensitivity is S1-only (see ST4).
- Cross-gene modules M(d) and hit counts treat the genes of one complex (e.g. mito-ribosome, CCT, SAGA) as independent
  evidence. A "recurs across 12 genes" module can be one biological signal.

**Change.**
- **§5.3:** cluster-robust (CR2) SEs with PatientID as the cluster for every OLS test in every list. In R, also
  cluster by patient if any R lines share one.
- **§6:** the LMM (line-similarity kernel from **D-only** dependency profiles of non-eligible genes, to avoid using
  the tested gene) is run on the top 1,000 of each list. Report the fraction that stays at nominal p < 0.001 under
  the LMM.
- **§5.2:** for M(d), cluster genes by dependency correlation in D (|r| > 0.5, average linkage). Report a module's size
  both as genes and as the **effective number of gene clusters**. Benchmark counts (§6.3–6.4) are also reported at the
  gene-cluster level.

**Check.**
- Every result table has a `se_type == "CR2_patient"` column.
- `results/benchmark/lmm_sensitivity.csv` has rows for all lists.
- The module table has `n_genes` and `n_gene_clusters`. A test checks n_gene_clusters ≤ n_genes and that the
  clustering threshold is read from config.

### R1-ST11 · major · The positive-control gate can tune the pipeline on the clinical ground truth

**Problem.**
- About half of the §7.1 controls (ERBB2, EGFR, BRAF, KRAS, ESR1, PRMT5/MAT2A–MTAP, WRN–MSI) are also CGT pairs.
- "If the gate fails, stop and diagnose" invites pipeline changes that are scored on those same pairs, so CGT recall
  becomes partly in-sample.
- The ≥ 60% threshold on about 15 evaluable controls has a wide binomial CI (e.g. 9/15 = 0.60, 95% CI 0.32–0.84).

**Change (§7.1).**
- "Any pipeline change made after the first gate evaluation is logged in `results/gate_log.csv` (commit, change,
  reason) and applied to **all** lists."
- "If any change is logged, CGT recall (§7.3 Q4a) is reported on two sets: all pairs, and CGT pairs **not** among the
  §7.1 controls. The latter is the primary set for the MD vs S1/EN-obs contrast."
- "The gate result is reported with an exact binomial CI. Gate pass/fail is not used to choose between MD variants."

**Check.**
- `gate_log.csv` exists (it may be empty).
- If it is non-empty, `results/cgt/recall.csv` must contain the column `set ∈ {all, excl_controls}` and the
  benchmark script must read `excl_controls` as primary. A script asserts the overlap list between §7.1 and the CGT
  table is computed from both config files, not hard-coded.

### R1-ST12 · minor · The toy simulation cannot test the risks above; add a fold-artefact check

**Problem.**
- §9 builds the toy "predictor" as truth + noise, independent of the observed noise. That cannot reproduce the
  selection-test dependence (ST1) that most needs testing.
- The toy has no collinear 𝒟 blocks and no library–lineage nesting.
- The Leiden check does not test association with the CV fold (Part 1 pitfall D3).

**Change (§9 and §5.2).**
- **Toy predictor:** the toy predictor is a **ridge model trained by 5-fold CV on the noisy toy observations** from
  features that overlap 𝒟. Its fold-centred OOF output is the input to MD.
- **Toy design additions:** add correlated feature blocks (r ≈ 0.9, a co-location analogue) and a library factor
  nested in 3 of the 12 lineages with its own effect.
- **Full-null pass rules:** under a full null, the ST2 estimator must return FDR_emp ≤ 0.10 for MD, EN-obs and S1,
  and MD's nominal D-stage q must be shown to be anti-conservative. The latter is a positive control for ST1: if it is
  not anti-conservative, the toy is too weak.
- **§5.2 fold check:** report AMI(Leiden clusters, cv_fold). It must be < 0.02, or the clustering is flagged as an
  OOF artefact.

**Check.** `tests/test_toy_pipeline.py` asserts four things: the toy predictor is fitted (not truth + noise), by
checking the presence of a fitted model object and the CV fold vector; null FDR_emp ≤ 0.10 for each list; the nominal
MD q under the null yields more hits than expected (the anti-conservativeness flag is True); and
`results/subtypes/leiden_ami.csv` has a `cv_fold` row.

---

**Out of my area / no change needed:**
- The §1 R-access guard and the freeze tags are statistically sound as specified.
- The −0.5 dependency threshold with a −1.0 sensitivity is fine.
- The §7.4 matched-background comparison for the combinatorial screens is adequate. Add a Fisher exact test with a CI
  for the hit-rate contrast.
