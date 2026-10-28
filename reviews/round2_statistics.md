# Review round 2: statistics (targeted review of the new v2 mechanisms)

Reviewer: statistician (round 2). Scope: the 2×2 arms, the hierarchical stability EN, shrinkage power and the top-k
P1 endpoint with its paired bootstrap, the joint-row permutation null and the FDR_emp pass rule, the covariate-drop
rule, the cross-lineage rule, CR1/CR2 SEs and the cross-library leak note. I read DESIGN.md v2, configs/design_params.json,
primary_endpoints.json, tiers.json, docs/part1_brief.md, reviews/round1_statistics.md and round1_response.md. The only
Part 1 access was one read-only carrier count on dev lines from P1 caches (pandas, no essml import):
STAG2 ~17, CREBBP ~14, MLH1 dosage-2 11, HOT:EGFR 13, ERBB2 AMP 41, MYCN AMP 23. Most gate controls are above
10 carriers. FUS, ALK, rare SIL and allele features are not.

For each verdict branch I checked by dry-run reasoning whether it can be reached. Items 1, 2, 4 and 7 contain
conditions that are either unreachable or vacuous as written.

---

### R2-ST1 · critical · Rare-carrier inference: robust SEs with df = n − p (and z in power) are anti-conservative; CR2 undefined at leverage 1

**Problem.**
- *Effective df.* For a binary d, the robust SE of β_d is set almost entirely by the carriers' residuals. Its
  sampling distribution therefore has about n_c − 1 effective df (the Welch case), not n − rank − 2.
- *Size at nominal 0.05 (one-sided).*
  - R tests with 3 carriers: actual size ≈ P(t_2 > 1.645) ≈ 0.12.
  - D tests with 5 carriers: actual size ≈ 0.09.
- *Size at BH-scale thresholds.* The inflation grows by orders of magnitude: P(t_4 > 4.9) ≈ 4e-3, against a nominal
  1e-6.
- *CR1 = HC1 has no leverage correction.* For singleton patients CR1 reduces to HC1, and a carrier's leverage is about
  1/n_c.
- *Consequences.*
  - F1 BH, the S1 top-k by p and the D-stage q_sc fill up with 5–9-carrier features carried by a single outlying
    carrier.
  - The power formula Φ(|β|/SE − z_0.95) overstates power for exactly these pairs.
  - The R "replicated" rule then has a false-positive rate near 0.12.
- *Leverage 1 in R.* Lineage, culture and medium levels are defined on D, and a level can have only 1 of the 255 R
  lines. That line has h_ii = 1, so (I − H_gg)^−1/2 in CR2 is undefined.
- *Bootstrap gaps.*
  - Duplicated patients in a resample are not given new cluster ids.
  - Resamples with < 2 distinct carrier patients have no stated rule.
- *Conflict 9's premise is wrong.* It says CR2 is too slow for S1. For singleton patients CR2 = HC2, and
  h_ii = h_Z,ii + d̃_i²/Σd̃². Here h_Z,ii is shared within a NaN group and d̃ is the FWL residual already computed, so
  HC2 is vectorisable. The ~50 multi-line patient blocks are 2×2 to 4×4 matrices.

**Change (§5.3, §5.4, §6.1; `confounders_s5.se_type`, `replication_s6.power_formula`).**
1. **SE.** CR2 for every test (HC2 for singleton patients), vectorised as above, with `se_type = "CR2"` everywhere.
   CR1 remains only as a logged sensitivity.
2. **Reference distribution.** Use t with df = min(n − rank(Z) − 2, ν), where ν = max(1, (Σ_c W_c)² / Σ_c W_c² − 1)
   and W_c = Σ_{i∈c} d̃_i², summed over the lines i of patient cluster c (d̃ = the FWL-residualised feature).
   - This is a Satterthwaite approximation for the HC variance. For a rare binary d it gives ν ≈ n_c − 1; for a
     continuous d it gives ν ≈ n/3 or more.
   - Pre-declared fallback, if the calibration test fails: for binary d, df = (number of carrier clusters) − 1.
3. **Power.** power = 1 − F_t(t_{0.95, ν_R}; ncp = |β_shrunk| / SE_R_pred), where ν_R is computed from the R carrier
   counts (molecular data only). This is the same df rule as the test it predicts.
4. **Leverage.** In each fit, merge any design level with < 2 lines in the fit's rows into its block's reference level
   (lineage → `_rare`), and log the merge. Exclude rows with 1 − h_ii < 1e-8 from the CR2 meat, and record
   `n_lev1_excluded`.
5. **Bootstrap.**
   - Give resampled patient copies distinct cluster ids.
   - A pair with < 2 distinct carrier or < 2 distinct non-carrier patients in a resample counts as not replicated.
   - Report `frac_degenerate` per arm.

**Check.**
- `tests/test_inference_calibration.py`: 20,000 null simulations for each n ∈ {255, 1022} × n_c ∈ {3, 5, 10} ×
  {homoskedastic, carrier SD × 2}, using the real D/R covariate design. Empirical size at α = 0.05 must lie within
  [0.035, 0.065], and size at α = 1e-3 must be ≤ 2e-3.
- `tests/test_s1_vectorised.py` compares against a brute-force numpy CR2 + ν on 200 tests (statsmodels has no CR2):
  |Δlog10 p| < 1e-6.
- A script asserts that the power function's df equals the R test's df for every pool row.

### R2-ST2 · critical · The null re-runs adaptive thresholds under the complete null, so FDR_emp is biased low; the toy and NC1 bounds are ill-posed

**Problem.**
- *Adaptive thresholds go strict under the complete null.* §7.1 re-derives λ_1SE by CV and re-runs BH inside every
  permutation draw.
  - Every 𝒟 row is permuted, so CV picks λ near λ_max, and BH under the complete null makes a discovery with
    probability ≤ q.
  - So |hits_b| ≈ 0 by construction, while the observed run operates at a looser λ and a looser BH cutoff.
  - FDR_emp = mean|hits_b| / |hits_obs| therefore measures false positives at a stricter threshold than the one that
    produced the observed hits. The bias runs in MD-strict's favour, and MD-strict is the list whose only error rate
    is this number.
- *The toy bound is vacuous or ill-posed.* The toy check "full null FDR_emp ≤ 0.10" fails in either reading:
  - On a pure-null toy, every observed hit is false. FDR_emp is then 0/1 (vacuous pass) or ≈ 1 (fail) by chance.
  - On the planted toy, the true FDP is known and is the quantity that should be checked.
- *NC1's replacement bound is degenerate.* `simulate_control_bounds.py` would replace the NC1 bound (0.10) with the
  95th percentile of the toy null's FDR_emp, which is degenerate (0 or ≈ 1). The calibration target would become
  vacuous.
- *NA masks.* The design does not say whether `CN_absent`, `MSI_unknown` and the S1 NA masks move with the permuted
  𝒟 rows. If they stay fixed, imputed zeros land on lines whose flag says "profiled".

**Change (§7.1, §13.1, `null_s7`, `toy_s13.pass`, `negative_controls_s7`).**
1. **Primary null = fixed-threshold null.** In each draw, apply the observed run's operating points:
   - the per-gene, per-stage λ_1SE;
   - SIS K, freq ≥ 0.6 and sign ≥ 0.9;
   - the observed F1 BH p-cutoff p*_F1, taken from the full observed family, not the subset;
   - the observed MD-strict cutoff p*_sc, with the same-sign rule.

   FDR_emp_fixed = mean_b |hits_b| / max(1, |hits_obs|), with a bootstrap CI over draws, as in SAM. The pass rule
   (≤ 0.10) applies to FDR_emp_fixed. This version skips CV, so it costs less than the current one.
2. **Secondary null.** The current full re-run becomes `FDR_emp_rerun`, with B = 10 in the same 180-min cap after the
   primary draws. It is reported as "selection-adaptive" and is not the pass rule.
3. **Toy pass rules.**
   - On the planted toy: FDP_true(list) ≤ 0.10 + 0.05 (MC tolerance), using plant labels, and
     FDR_emp_fixed ≥ FDP_true − 0.05.
   - On a pure-null toy: P(|hits| ≥ 1) ≤ 0.10 for the BH lists, and mean |hits| reported for the EN lists.
4. **NC1** stays at 0.10 and is excluded from simulated replacement. It is a pre-declared target, not a data bound.
5. **NA handling.** NA indicators and NA masks permute with their 𝒟 rows. The covariates own CN and EXPGATE stay
   fixed, and the design says so explicitly.

**Check.**
- `results/null/summary.json` has keys `fdr_emp_fixed`, `fdr_emp_rerun`, `p_cut_F1`, `p_cut_sc` and `lambda_source`
  ("observed").
- A script recomputes p_cut from the observed F1 table.
- `simulate_control_bounds.py` refuses to touch the `NC1_fdr_emp_max` key.
- `tests/test_permutation_null.py` asserts that `CN_absent` equals the permuted CN row's NA pattern.
- The toy test asserts both planted-toy inequalities.

### R2-ST3 · major · The covariate-drop rule is undefined for categorical blocks and would remove lineage/library/culture adjustment from proxy features

**Problem.**
- *Undefined for categorical blocks.* "|Spearman(cov, d)| > 0.7" has no definition for lineage, library group,
  CellFormat, medium family, MSI {3 levels} or TP53 {3 levels}. An implementer will compute it on integer codes, which
  are arbitrary.
- *Wrong when it can be defined.* For a lineage proxy (SIL with Cramér's V > 0.5, a PROG lineage factor, an
  Ewing-type FUS), the rule deletes the confounder and keeps the proxy. Lineage would become a "determinant". This
  is the opposite of the intent, which was to stop alternative names for d (MSI, DMG_BURDEN, TP53, own CN) from
  adjusting d away.

**Change (§5.3, `confounders_s5.drop_rules`).**
- Droppable covariates are only {MSI, DMG_BURDEN, TP53 status, own CN, EXPGATE own, screen-days}.
- For a categorical covariate, the association is R²(d ~ the covariate's dummies), with drop at R² > 0.49 (≡ |r| >
  0.7).
- Lineage, library and culture are never dropped. If R²(d ~ lineage + library + culture) > 0.49, the test runs with
  them and is labelled `collinear_with_design`. These pairs go to §12 class 6 (or to "lineage-unresolvable" in §8.3)
  and are not counted as determinants.

**Check.**
- A scan of all test tables asserts `dropped_cov ⊆ {MSI, DMG_BURDEN, TP53, own_CN, EXPGATE, screen_days}`.
- It also asserts that every row with R²(design) > 0.49 carries `collinear_with_design = True`.
- Unit test: a synthetic d equal to a lineage dummy keeps lineage and is flagged.

### R2-ST4 · major · The in-fit carrier floor makes 5–9-carrier features unselectable by EN (an unreachable branch)

**Problem.** §4.1 step 1 drops a column with < 5 carriers in the fit's rows.
- *Complementary halves cap the frequency.* The two halves of a pair split the m D carriers.
  - For m ≤ 9, at most one half of each pair can hold ≥ 5 carriers, so selection frequency ≤ 0.5 < 0.6 always.
  - For m = 10, the expected maximum frequency is ≈ 0.62.
- *Effect on EN-obs, EN-pred and MD-strict.* These lists can never contain features with 5–9 carriers, although they
  are in U and gate-eligible (≥ 5).
  - S1 can rank such features, so the arms differ structurally on rare features: FUS, ALK (≤ 10 by design), rare SIL
    and allele-level CGT features. This affects P2 (CGT recall over U) and the gate's denominators.
  - The gate is decided on MD-strict, so any 5–9-carrier control is a guaranteed miss.
- *Even 3 is too high.* A floor of 3 still caps m = 5 at 0.5.

**Change (§4.1 step 1; `surrogate_s4`, new key `min_carriers_in_fit: 2`).**
- Drop a column in a fit only if it has < 2 carriers in that fit's rows. Two is the minimum for a non-degenerate
  column after FWL. With a floor of 2, the expected maximum frequency is ≈ 0.81 at m = 5.
- `pc_eligibility.csv` and `arms.csv` gain an `en_max_attainable_freq` column, computed from the committed subsample
  manifest. Rows with a value < 0.6 are flagged `en_unreachable`, and the gate and P2 report them separately for every
  arm.

**Check.**
- From the manifest, a script asserts that `en_max_attainable_freq ≥ 0.6` for every feature with ≥ 5 D carriers
  (expected to hold at a floor of 2), or lists the exceptions.
- The toy gains one 6-carrier planted feature with a strong effect. It must be selectable, i.e. frequency > 0 in
  > 50% of subsamples.

### R2-ST5 · major · Continuous features have no definition in U, R-evaluability, power or the shrinkage prior

**Problem.**
- *No carrier definition.* "≥ 5 carriers", "≥ 3 carriers and ≥ 3 non-carriers with expected non-NaN effect" and
  "SE_R_pred from … R carrier counts" are defined only for binary d. Continuous d have no carriers: EXPDOSE, EXPPAR,
  PROG, PSCORE, DMG_BURDEN, ploidy, CIN and SPEC_*.
- *Consequence for the pool.* These families dominate EN selections, so pool membership and power for a large share of
  every arm's top-k are undefined.
- *Mixed scales in the prior.* The EB prior is "fitted once on all F1 (β, SE)". That pools mean differences (binary)
  with per-unit slopes (continuous, standardised or not), so the scale mixture is fitted on non-exchangeable effects.
- *Prior coverage.* F2 rows (EXPDOSE own) are in U but not in the prior fit.

**Change (§1.6, §6.1; `replication_s6`, `primary_endpoints.universe`).**
- *Continuous d.*
  - In U if ≥ 20 non-NaN D lines and SD > 0 among them.
  - EXPDOSE counts as "carriers" its expressed lines (log2TPM+1 ≥ 1), with ≥ 5 required.
  - R-evaluable if ≥ 20 non-NaN R lines and SD_R(d̃_R) > 0.
  - SE_R_pred = σ_D,resid / (SD_R(d̃_R) · √n_R), where d̃_R is d residualised on the R covariate design (molecular
    data only).
- *Prior.* Fit the EB prior on standardised effects β* = β · SD_D(d̃) (with matching SE*). Apply it to every row of U,
  including F2, and back-transform.

**Check.**
- A script asserts that every U row has non-null `power`, `evaluable` and `d_type ∈ {binary, continuous}`.
- `power` uses one function hash.
- A unit test: rescaling a continuous feature by 10 leaves its power unchanged.

### R2-ST6 · major · P1 fairness: covariate sets differ between EN and S1, F2 rows can saturate top-k, and EN ranks on ties

**Problem.**
- *(a) Covariate sets differ.* EN arms are residualised on lineage + library only. S1 arms, the D-stage screen and
  every R test also adjust for culture. An EN-only selection that proxies culture (a PROG organoid/suspension factor)
  is then penalised in R for a covariate-set difference, not for being multivariate. The contrast "EN-obs − S1-obs =
  multivariate adds value" is confounded by this. MSI, burden and TP53 cannot enter the EN FWL block, because they are
  𝒟 features (WRN←MSI would be FWL'd to zero). Culture is not in 𝒟, so nothing forces its omission.
- *(b) F2 rows in the top-k.* U includes F2 (own EXPDOSE), and EN forces it. These near-tautological pairs replicate
  almost surely, so all four top-100 lists converge, the rates saturate, and the two contrasts lose power. F2 was
  separated in §5.6 precisely because these rows dominate S1.
- *(c) Ties at frequency 1.0.* With a fixed λ, many (g, d) reach freq = 1.0. The EN top-100 is then decided by the
  tie-break "|mean β|", whose averaging set (all 50 subsamples, or selections only) is not defined.

**Change (§4.1, §6.1; `surrogate_s4.fwl_unpenalised`, `P1_primary.candidate_pool`, `rank_rule_EN`).**
- (a) Add the culture block (CellFormat, SerumFree, medium family) to `fwl_unpenalised` for both EN arms. If PI
  ruling 3 is kept instead, rename the contrast to "EN pipeline vs S1 (full model)" and report the culture-block delta
  (§5.5) for the EN-only strata.
- (b) The P1 pool excludes F2 rows for all arms. Top-k including F2 becomes a secondary variant.
- (c) Define mean β as the mean over all 50 subsamples with zeros included. Report `n_tied_at_rank_k` per arm.

**Check.**
- A config test checks that `fwl_unpenalised ⊇ culture` (or that the renamed contrast key is present).
- A script asserts that no P1 primary list row has `family == F2`.
- `arms.csv` has `mean_beta_def = "all50"`.
- `replication_by_k.csv` has `n_tied_at_rank_k`.

### R2-ST7 · major · The cross-lineage label is unreachable at k = 2–3 lineages, and label precedence is undefined

**Problem.**
- *Unreachable at small k.* The meta test uses HKSJ, a t test with k − 1 df over the contributing lineages. With
  k = 2 (the minimum the rule allows), a two-sided p of 1e-2 needs t ≈ 64; with k = 3 it needs t ≈ 10. Under F4 BH
  thresholds, "cross-lineage" is effectively unreachable unless ≥ 4–5 lineages are testable. That contradicts the
  stated minimum of "≥ 2 powered lineages", and it biases Q1b.
- *Precedence.* The labels are not mutually exclusive. A feature can have 1 powered lineage that is significant
  (q < 0.05) with interaction p < 0.05, and so qualify as both "lineage-restricted" and "lineage-unresolvable".
- *Interaction not estimable.* The interaction LRT cannot be estimated when only one lineage has carriers (MYCN).
- *"Powered" has no α.* The α behind "powered" is not stated.

**Change (§8.3; `subtypes_s8.cross_lineage`).**
- *Pooled significance.* Take it from the pooled confounder model restricted to the testable lineages (lineage fixed
  effects, CR2 + ν from R2-ST1), with BH over F4.
- *Descriptive meta.* REML + HKSJ CI, I², τ and the prediction interval are reported, not gated.
- *Powered.* Per-lineage one-sided α = 0.05 under the R2-ST1 power rule.
- *Precedence.* cross-lineage > lineage-restricted > lineage-unresolvable > undetermined.
- *Lineage-restricted* additionally requires ≥ 2 testable lineages, so that the interaction is estimable. A feature
  with carriers in a single lineage gets `lineage-unresolvable (single-lineage carriers)`.

**Check.**
- `tests/test_cross_lineage_labels.py` feeds synthetic statistics and asserts that each of the 4 labels is reachable.
  In particular, a 2-lineage, 2-organ-group consistent effect (β = −1, 15 carriers each) is labelled cross-lineage, and
  each case gets exactly one label.
- The label script asserts precedence.

### R2-ST8 · minor · Leakage notes: XLIB universe selected on replication-line effects; pseudo-R dry run is in-sample

**Problem.**
- *XLIB.* The XLIB caveat covers only the pred arms. Gene eligibility (`r_within`), the variable/selective classes,
  `frac_below` (relative rule) and the m5 tiers all use effects of all D lines, including the no-Avana replication
  lines. So U, which all arms share, is pre-selected on replication-set effects. This inflates the absolute XLIB rates
  for every arm, including EN-obs and S1-obs. Paired contrasts survive, but absolute rates do not.
- *Pseudo-R dry run.* With D fold 0 as stand-in R, fold-0 effects were used to build every arm's list. The fold-0
  lines' OOF predictions are blind to them, but the other folds' predictions were trained on them. So the dry run is
  in-sample and would show artefactual arm differences.

**Change (§6.3, §1.3).**
- Add to the XLIB caveat: "U is selected on all-D effects; absolute XLIB rates are optimistic."
- Add the sensitivity `u_rule = discovery_only`: eligibility recomputed on Avana-only lines (OOF predictions plus
  discovery effects, no refit).
- Dry-run outputs get registry status `dryrun`. They are uncitable and are not inspected for arm differences.
  Branch reachability is shown only by the synthetic unit test.

**Check.**
- `replication_xlib.csv` has rows with `u_rule ∈ {all_D, discovery_only}`.
- `check_report_numbers.py --strict` rejects any tag resolving to status `dryrun`.

---

**Adequate as written.** The following need no change beyond the items above:
- the P1 verdict branches in `primary_endpoints.json`, which are mutually exclusive and exhaustive, with all three
  reachable;
- the paired-bootstrap design with fixed lists (a conditional inference on the frozen lists, which REPORT.md should
  state as such);
- the EB-shrunk power as a shared, pre-tag property of each (g, d);
- the joint-row permutation unit and its strata;
- the MD-strict wording rule;
- the drop-one-block delta rule.
