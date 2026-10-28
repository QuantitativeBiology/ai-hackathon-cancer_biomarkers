# DESIGN — Predictor-derived dependency subtypes and their molecular biomarkers

Version: v3.2 (after review round 3; includes PI amendment D7: methylation). Status: frozen at tag `freeze-design`;
bounds frozen at tag `freeze-bounds` at the end of M3b (§13.4).
Inputs: PROMPT.md, docs/part1_brief.md, LITERATURE.md, reviews/round{1,2,3}_*.md and their response files. Decision
numbers live in §F and its config files; prose numbers are descriptive.

## 0. Questions, claims, endpoints

- **Q1 subtypes.** For dependencies the Part 1 predictor predicts accurately *beyond lineage and library*, can its
  reasons be read out as assayable molecular determinants that group lines dependent for the same reason, within and
  across lineages? **Q1b**: do cross-lineage determinants replicate, or match clinical biomarkers, more often than
  lineage-restricted ones (§8.3)?
- **Q2 benchmark.** Does the predictor add information beyond the observed data (EN-pred vs EN-obs), and does sparse
  multivariate selection add beyond univariate tests (EN-obs vs S1-obs)?
- **Q3 confounders.** Which candidates survive lineage, library, culture/growth, passengers, TP53, co-location,
  cutting artefacts and multiplicity?
- **Q4 clinical.** Recovery, direction, actionability and patient population against a frozen clinical ground truth.
- **Q5 limits.** Every eligible or control dependency gets a failure class (§12).

**Falsifiable novelty claims** (vs Tsherniak 2017, Dempster 2020, Behan 2019, Pacini 2024): (a) predictor-derived
lists replicate better than, or name determinants missed by, EN-obs and S1-obs; (b) line-level subtypes with a tested
cross-lineage label; (c) scoring a dependency catalogue against tiered clinical ground truth with a decomposition of
misses.

**Endpoints** (`configs/primary_endpoints.json`): P1 (Q2) = R replication of the top k_G per arm in the *genomic*
universe U_G, reported by access stratum (§6.1); P2 (Q4a) = CGT recall at k = 200 (§10.5). Everything else is
secondary. **The first deliverables a reader sees** are `results/candidates/known_recovery.csv` and `shortlist.csv`
(§12.3); REPORT.md opens with them. Non-goals: a new predictor; outcome or drug-response modelling.

## 1. Data, roles, protection

| Role | Lines | Use |
|---|---|---|
| D (discovery) | 1022 Part 1 dev lines, `cv_fold` 0–4 | everything before `freeze-catalogue` |
| R (replication) | 255 Part 1 test lines | read once, after `freeze-catalogue`, by one `scripts/replicate.py` run |
| set A | 15 non-cancerous | same `replicate.py` run (scope `limits`): normal-line dependency flag (§12.3) |

1.1 **Predictions.** D: Part 1 OOF headline `m5_avg` (`p1.oof_predictions("m5")`), fold-centred per gene
(`p1.fold_centre`); component OOF runs m1/m2n/m3n (mean = M5 to 1e-6). R: all-dev refit `p1.test_predictions("m5_avg")`.
Targets: Chronos 26Q3 via `guard.load_effects` only.

1.2 **What R is.** Part 1 opened R's effects once to score its predictor; no biomarker analysis has used R. R shares
the Chronos fit with D, so replication tests generalisation across lines, not across processing; it is also reported
for R lines screened without Avana.

1.3 **Guard (extends `src/biomark/guard.py`).**
- *Single entry point.* Only `guard.py` may reference CRISPR effect files (`crispr.parquet`, `CRISPRGene*`,
  `ScreenGeneEffect`), `results/test`, `cache/test`, `essml.guard` or `io.load_matrix`; only `p1.py` imports `essml`; no `subprocess` targets ../project; only `scripts/replicate.py`
  calls `load_effects` with scope `replication` or `limits` (AST scan). The deny-list adds `P1/cache/test/**` except
  `preds/`.
- *Unblind conditions* (added to (a) tag, (b) `PART2_UNBLIND=1`, (c) clean tree, (d) purpose): (e) caller is
  `replicate.py`; (f) sha256 of every frozen catalogue file, including `configs/catalogue/r_requests.csv`, matches
  `configs/freeze_catalogue.json`; (g) the tag target equals the one recorded at the first access; (h) afterwards only
  reproductions (same config hash and commit).
- *One R run.* `r_requests.csv` lists every (g, d) needing R statistics, with reason ∈ {P1/P2 list row at any k, arm or
  variant; MD-strict; EN-obs-strict; S1-obs BH; cross-lineage/restricted label; control; CGT test; shortlist}.
  `replicate.py` writes, in one run, `results/replication/r_tests.csv` (β, SE, ν, p, n, carriers, power, bootstrap
  draws for the endpoint lists), `r_lines.csv` (R accurate-dependent membership for MD-strict subtypes) and
  `setA_lines.csv`. Every later R-dependent script reads only these files.
- *Dry run.* Before the tag, `replicate.py --pseudo-R fold0` (D fold 0 as stand-in R) must produce the same file set
  and every downstream R-consuming table (Q1b, pair and gene classes, cards). Its registry status is `dryrun`:
  uncitable and not inspected for arm differences (it is in-sample). Verdict-branch reachability is shown by a
  synthetic unit test.
- *Immutability.* After the first access, the diff of `configs/catalogue/`, `replicate.py` and replication modules since
  the tag must be empty.

1.4 **Part 1 untouched.** `scripts/snapshot_part1.py` writes `configs/part1_manifest.csv` (path, size, mtime_ns,
excluding `data/`, which `verify_inputs.py` sha256-checks); `tests/test_part1_untouched.py` compares it, HEAD and
cleanliness at every gate.

1.5 **Downloads.** Obtained (`external/depmap_extra/`, sha256 in `configs/depmap_extra_sha256.txt`):
OmicsSomaticMutations (ProteinChange, Hotspot, LikelyLoF, AF), OmicsGlobalSignatures (MSIScore, LoHFraction, WGD, CIN,
Ploidy, Aneuploidy), CRISPRConfounders (ScreenDoublingTime, screen QC, CasActivity), OmicsCNSegmentsWGS (unused),
OmicsProfiles (keys only), and CCLE RRBS promoter methylation (CCLE_RRBS_TSS1kb_20181022; amendment D7). Not obtained:
PRISM (unused). **Data deadline = the first dictionary build**: `build_dictionary.py` reads `results/data_status.json`
(committed first) and stores its sha256 in `dictionary_v1.json`; a later file feeds only labelled secondaries, never a
rebuild. MSI = MSIScore ≥ 20; Cell Model Passports `msi_status` is a concordance check only.

1.6 **Missing data and row sets.** 𝒟 is built from raw caches. CN-derived features are NA for lines without CN; in EN
fits NA → 0 plus `CN_absent`/`MSI_unknown`/`METH_absent` indicators, which are attributes of 𝒟 rows and move with them under
permutation; in S1, D-stage, subtype and carrier counts NA lines are excluded for that family. A separate process writes
`cache/effect_nan_mask.parquet` (sha256 in `results/bulk_manifest.csv`), the only effect-derived input EN-pred may read;
both EN arms fit identical row sets per gene. `n_obs, n_carriers, n_NA` are written per test.

1.7 **D-only fitting.** NMF, scalers, SIL medians, locus groups, programme scores and thresholds are fit on D and
frozen; other lines are transformed. A test rebuilds them with R rows deleted and asserts identical hashes.

## 2. Eligibility and "dependent"

- **Eligible gene** = (variable ∨ selective) ∧ m5 tier ∈ {high, medium} ∧ `r_within ≥ 0.10`, where `r_within` is the
  D Pearson r between fold-centred OOF prediction and observed effect after centring *both* within lineage ×
  library-group strata (≥ 5 lines; rare lineages pooled). Failing genes go to §12 class 3 (`results/eligibility.csv`).
- **Non-eligible control and CGT-target genes** run in a separate universe `U_ctrl`, excluded from every top-k list,
  the null subset and P1/P2 (scripts assert arms ⊆ eligible and `null_genes.txt` ⊆ eligible).
- **Associations use continuous effects.** "Dependent" serves subtype membership only: effect < −0.5; for genes with
  frac_below(−0.5) > 0.5, effect < the gene's D 25th percentile (`dependent_rule = relative`). Sensitivity −1.0.
- **Accurate-dependent** (ℓ, g): dependent ∧ fold-centred prediction ranks within the n_g most negative (n_g = dependent
  D lines). In R: prediction ≤ the R quantile at n_g/n_D (target-free).

## 3. Interpretable feature dictionary 𝒟

`scripts/build_dictionary.py` → `configs/dictionary_v1.json` (id `FAMILY:<entrez or name>`, family, gene(s),
definition, d_type, carriers or non-NaN counts D/R, NA count, Cramér's V with lineage, locus group, data_status
sha256). Γ = OncoKB cancer genes ∪ Part 1 top-3 Ensembl paralogs of eligible genes ∪ partners in
`configs/positive_controls_spec.csv` ∪ known-silenced list. Binary features need ≥ 5 D carriers; continuous ones ≥ 20
non-NaN D lines and SD > 0 (EXPDOSE: ≥ 5 expressed lines).

| Family | Definition | Scope |
|---|---|---|
| HOT | hotspot > 0, gene-level | all |
| ALLELE | allele classes from OmicsSomaticMutations per `configs/cgt_allele_rules.json` (e.g. BRAF_V600, BRAF_nonV600, KRAS_G12C, KRAS_G12D/V, KRAS_G13, KRAS_Q61, NRAS_Q61, EGFR_activating, PIK3CA_helical, PIK3CA_kinase, AKT1_E17K, IDH1_R132, ERBB2_TKD, KIT_ex11, KIT_D816, FGFR3_hot) | listed genes |
| LOFBI | dosage 2; or dosage ≥ 1 with own CN < 0.75, SIL, or a LikelyLoF call with AF ≥ 0.7; or X-linked dosage ≥ 1 in males | Γ |
| AMP / DEL | relative CN > 2 / < 0.25 | AMP: Γ ∪ own gene; DEL: all |
| LOSS | relative CN < 0.75 | gene-specific: own gene, paralogs, spec partners |
| FUS | DepMap filtered fusion, FFPM ≥ 0.1, driver from `configs/fusion_drivers.csv` as 3′ gene (5′ for FGFR2/3) | drivers |
| SIL | log2(TPM+1) < 1 ∧ expressed in ≥ 50% of own lineage (≥ 10 D lines, else pan-D) ∧ < lineage median − 2; mechanism SIL-DEL (CN < 0.25), SIL-MUT (dosage ≥ 1), SIL-METH (promoter methylation ≥ 0.5, D7), SIL-other | Γ; Y genes within males |
| METH | promoter hypermethylation (max TSS±1kb RRBS fraction ≥ 0.5 and ≥ 0.3 above the D lineage median), restricted gene set, NA without RRBS (`dictionary_s3.methylation`, D7) | Γ |
| YABS | Y paralog absent (female, or male with SIL of the Y gene) | gene-specific |
| CLIN | clinical composites from `build_cgt_features.py` (molecular data only): ER_POS, HER2_POS, HRPOS_HER2NEG, BRCA_BIALLELIC, RAS_MUT | all |
| EXPDOSE | own expression among expressed lines (0 elsewhere, centred) | gene-specific; the only own-expression determinant |
| EXPGATE | own expressed vs not | covariate only |
| EXPPAR | expression of top-3 paralogs | gene-specific |
| PROG | 50 NMF factors (5,000 most variable genes, nndsvda, cd, 500 iter, `nmf_seed`) | all; "expression programme (no clinical assay)" |
| PSCORE | ISG (Hallmark IFN-α), NRF2 target score | all |
| SIG | MSI (MSIScore ≥ 20), MSIScore, Ploidy, WGD, CIN, Aneuploidy, LoHFraction, TMB (log1p protein-altering calls), DMG_BURDEN; SPEC_MMR/POLE/APOBEC/TOBACCO/UV (coding calls, ≥ 10 SNVs; sanity checks §F) | all |
| COMP | TP53MUT (absence biomarkers coded as presence with flipped sign), RB_ALT | all |

- **Stage and universe families.** U_G (genomic universe) = HOT, ALLELE, LOFBI, AMP, DEL, LOSS, FUS, SIL-DEL, SIL-MUT, SIL-METH, METH,
  YABS, CLIN, SIG, COMP (= EN stage G). Expression families EXPDOSE, EXPPAR, PROG, PSCORE, SIL-other (= stage X).
- **Self-null features** are removed from 𝒟_g and kept as checks (NC3): DEL_g, SIL_g (all mechanisms incl. SIL-METH_g),
  METH_g, LOFBI_g, DEL_h for h in g's DEL locus group. LOSS_g (CYCLOPS) and AMP_g stay.
- **Variant groups** for de-duplication (§5.1): locus groups (genes within 10 Mb whose binary calls have Jaccard ≥ 0.5
  among CN-profiled D lines, per AMP/DEL/LOSS family, Ensembl 116 coordinates) and allele groups (HOT:x with its
  ALLELE:x_* and CLIN features built from them) and SIL groups (the SIL mechanisms of one gene, which are mutually exclusive
  features with precedence DEL > MUT > METH > other). Features stay separate; CDKN2A/CDKN2B/MTAP form one group.
- **Cap.** Shared 𝒟 ≤ 3,000; if exceeded, the largest family's minimum frequency rises to 10, then 20 (logged). SIL with
  Cramér's V(lineage) > 0.5 is flagged.
- **Not assayed:** methylation for lines without RRBS (`METH_absent`; record in `dictionary_v1.json`); coverage is
  reported (§12.2) and REPORT.md lists it as a limit.

## 4. Reading determinants out of the predictor

### 4.1 Hierarchical stability elastic net (EN-pred, EN-obs)
Target: fold-centred ŷ_g (EN-pred) or observed y_g (EN-obs); identical code, configs, subsamples and row sets.
1. **FWL.** 𝒟 standardised once on D. Every fit residualises target and features on the unpenalised block (lineage,
   library group, culture block, growth block, `CN_absent`, `MSI_unknown`, `METH_absent`; §5.3) estimated on that fit's rows; no
   intercept. Residualised 𝒟 is cached per (NaN-pattern group × fold or subsample), memmapped read-only. A column is
   dropped in a fit if residual SD < 1e-6 or < 2 carriers in the fit's rows.
2. **Stages.** Stage G on the residualised target; stage X on that fit's own stage-G residual.
3. **Screen.** Per stage, SIS to K = 300 by |corr| inside each fit; forced (still penalised): own EXPDOSE, EXPPAR, own
   LOSS/AMP, paralog/partner LOSS, YABS, the control's partner feature.
4. **λ.** l1_ratio 0.9; 25-λ path, eps 1e-2, tol 1e-3; CV over Part 1 `cv_fold`; 1-SE on full D, held fixed per gene
   and stage in all subsamples.
5. **Stability.** 25 complementary half-sample pairs drawn once (`subsample_seed`; patient-grouped, lineage-stratified)
   into `configs/subsamples.json`, shared by all genes and arms; each fit uses subsample ∩ non-NaN rows. Selected =
   frequency ≥ 0.6 and ≥ 90% of selections with the majority sign. Rank: frequency, then |mean β| over all 50
   subsamples (zeros included), then feature id; `n_tied_at_rank_k` reported. `en_max_attainable_freq` per feature
   (from the manifest) is written; rows < 0.6 are `en_unreachable` and reported separately.
6. **Fidelity** `r2_oof_beyond_design` = out-of-fold R² of the surrogate on the residualised target (folds above);
   class-3 threshold 0.10.
7. **Pilot.** 50 genes with the FWL cache (`results/timing/surrogate_pilot.json`, incl. cache hit rate); projected wall
   > 45 min ⇒ K = 150. Subsamples are never reduced.
8. **Blindness.** EN-pred runs with `load_effects` patched to raise and an audit hook on `crispr.parquet`
   (`tests/test_en_pred_blind.py`, which also asserts equal row sets with EN-obs).

### 4.2 Component readout (Tier A)
Full-D fits on fold-centred ŷ_M1, ŷ_M2n − ŷ_M1, ŷ_M3n − ŷ_M1 → `component_selection.csv`: "expression-mediated" vs
"genomics-layer". It names `carrier_input` for proxy features and is the fallback when §4.3 is cut.

### 4.3 M1-proxy direct attribution (Tier B, ≤ 300-gene panel: controls, CGT targets, top MD-strict genes)
Exact linear attribution from saved `t_E_status_24b46e2a` folds (T ≈ M1): per fold, W in gene blocks ≤ 500, aggregated
on the fly into EXP-own, EXP-paralog, PROG_k (least-squares coefficients of w_g on the NMF loadings in z-space) and
EXP-other; each line attributed by the fold model that did not train on it; contributions fold-centred; additivity
within 1e-5. A named share < 0.5 means "M1 unnameable" (class 3). Agreement = fraction of panel MD-strict (g, d) whose
family is among the top-5 terms of the accurate-dependent lines. No M3 refit.

### 4.4 Access class (every (g, d) in U)
- `direct`: d's source matrix and gene equal one of M2's 13 inputs for g (own expression/not-expressed/hinge, own CN,
  own deep deletion, own damaging mutation; paralog expression/minimum, deep deletion, damaging mutation).
- `potentially-direct`: d derives from a matrix in M3's pool (E, C, M, H) but is not an M2 input (selection
  unverifiable without the cut M3 refit).
- `expression-derived`: EXP/PROG/PSCORE/SIL-other/SIL-METH. `proxy`: not derived from any predictor input (MSI and
  other SIG, FUS, ALLELE, CLIN, COMP, METH).
Wording: "learnt d" only for direct and expression-derived; potentially-direct carries a qualifier; proxy reads "learnt
a correlate of d (via X)". **Claims that the predictor finds determinants beyond its hand-coded inputs come only from
the non-direct stratum; "beyond any input" only from proxy.**

## 5. Benchmark arms, MD-strict, confounder model

5.1 **Arms** on the identical universe U = {(g, d): g eligible, d ∈ 𝒟_g minus self-null, meeting §3 counts}:
**S1-obs** (univariate confounder-model test on observed effects; the standard), **S1-pred** (same on fold-centred
predictions), **EN-obs**, **EN-pred** (= MD-pure). `results/benchmark/arms.csv` holds (arm, g, d, family, access_class,
rank_stat, rank, d_type, en_max_attainable_freq); scripts assert equal universes, covariate hash, subsamples and seeds.
Before any top-k, features of one g in one variant group are de-duplicated to the best-ranked member (all arms). A CGT
`feature_D` is scored by its own `rank_of_feature_D` (§10.5): the group's best member never lends its rank.

5.2 **MD-strict** = EN-pred selected candidates that pass the D-stage screen on observed effects: confounder model,
same sign, BH q_sc < 0.05 within EN-pred selected candidates. The p and q are *selection-conditional (nominal)*; the
error rate is the empirical FDR (§7); confirmatory claims come only from R and external sources; REPORT.md never calls
it "significant". EN-obs gets the same screen (EN-obs-strict). MD-strict is the catalogue for subtypes and clinical
claims.

5.3 **Confounder model** (S1 arms, D-stage, R tests, ceiling); built by `scripts/build_screen_covariates.py` from the
*screened* condition (CRISPRScreenMap → ScreenSequenceMap without pDNA → ModelCondition), D ∪ R ∪ set A, no effects.
Category maps, regexes and NA levels are in `configs/build_params.json`.
- *Design blocks (never dropped):* lineage (< 10 D lines → `_rare`); library group (Avana-only, KY-only,
  Humagne-only, Avana+KY, other) + screen days within Humagne (dropped if VIF > 10); culture (screened CellFormat
  {adherent, suspension, 3D, mixed/unknown}, SerumFree {T, F, unknown}, base-medium family; the parser picks the
  formulation column by content because Model.csv OnboardedMedia/FormulationID appear swapped; lines without media →
  unknown); growth (ScreenDoublingTime, NA → D median + `DT_unknown`).
- *Droppable covariates:* MSI {MSI, MSS, unknown}, DMG_BURDEN, TP53 {mut, WT, unknown}, own CN of g, EXPGATE of g,
  screen days. Per (g, d), drop one if it is d or in d's variant group, or if its association with d exceeds 0.49 (R²
  of d on the covariate's dummies; ≡ |r| > 0.7 for numeric), or own CN when d is CN-derived and its group contains g
  (own-CN-adjusted = sensitivity), or TP53 when d derives from TP53 or MDM2 AMP. Logged in `dropped_cov`; a scan asserts
  `dropped_cov` ⊆ the droppable set.
- *Collinear with design:* if R²(d ~ lineage + library + culture + growth) > 0.49 the test keeps the design and is
  labelled `collinear_with_design`; such pairs go to §12 class 6 (or lineage-unresolvable) and are never counted as
  determinants.
- *Leverage:* in each fit, a design level with < 2 lines in the fit's rows merges into its block's reference level
  (logged); rows with 1 − h_ii < 1e-8 are excluded from the meat (`n_lev1_excluded`). VIF/condition number →
  `results/covariates/vif.csv` (max ≤ 10).
- *SE and df (every test):* CR2 with PatientID clusters (HC2 for singletons). For singletons
  h_ii = h_Z,ii + d̃_i²/Σd̃², with d̃ the FWL residual of d; the ~50 multi-line patient blocks (2–4 lines) use batched
  closed-form (I − H_cc)^−1/2 with the rank-1 update. Leverages depend on (NaN group, d) only, so they are computed once
  per feature and reused across genes. Reference t with df = min(n − rank(Z) − 2, ν), ν = max(1, (Σ_c W_c)²/Σ_c W_c² −
  1), W_c = Σ_{i∈c} d̃_i². Fallback if the calibration test fails: binary d uses df = carrier clusters − 1. CR1 is a
  logged sensitivity. `se_type = CR2` is a column of every table.

5.4 **S1 implementation.** Genes grouped by NaN pattern; QR-residualise Y and 𝒟 on the shared design; own CN and
EXPGATE added by FWL update; partial-regression t with §5.3 SE and df; F = t². Drop tests with < 5 carriers among non-NaN
lines or residual variance < 1e-8. `tests/test_s1_vectorised.py`: 200 random tests vs a brute-force numpy CR2 + ν
(|Δlog10 p| < 1e-6) and, with singleton clusters, vs statsmodels HC2 (1e-10).

5.5 **Confounding delta** (MD-strict ∪ S1-obs BH hits ∪ controls): drop-one-block refits (lineage, library, culture,
growth, burden, TP53, own CN); flag if |β_full − β_minus| > 0.5·|β_minus| and > SE_minus. Sensitivity block
OncotreePrimaryDisease (≥ 5 D lines) → "subtype-confounded" (`results/confounders/delta.csv`).

5.6 **Test families** (`test_families_s5`; n_F = rows, q recomputed by script): F1 S1 main; F2 S1 own-EXPDOSE
(excluded from P1 pools and from recall of non-expression biomarkers); F3 EN D-stage (screening only); F4
cross-lineage pooled tests; F5 per-lineage; F6 R replication (no BH); F7 CGT ceiling. BH under PRDS.

## 6. Replication endpoints

6.1 **P1 (primary).** Pool_G = (g, d) ∈ U_G, not F2, R-evaluable and R power ≥ 0.5. R-evaluable: binary ≥ 3 carriers and
≥ 3 non-carriers, continuous ≥ 20 non-NaN R lines with SD_R(d̃_R) > 0, using expected non-NaN effects (a gene is
expected non-NaN for an R line if ≥ 95% of D lines with the same library string are non-NaN; unseen strings use the
union NaN set of their components). k_G = min(100, ⌊|Pool_G|/3⌋), fixed in `results/benchmark/p1_pool.json` before the
tag. Top-k_G of each arm within Pool_G; replicated = same sign and one-sided p < 0.05 in R (§5.3 model, same drops,
CR2 + ν). Contrasts EN-pred − EN-obs ("predictor adds value") and EN-obs − S1-obs ("multivariate adds value"), each by a
paired patient bootstrap over R (1,000, `bootstrap_seed`; lists, drops and pool fixed from the full-R fit; resampled
patient copies get new cluster ids; zero-variance columns dropped by pivoted QR and logged; a pair with < 2 distinct
carrier or non-carrier patients in a draw counts as not replicated; `frac_degenerate` reported). Claim only if the 95%
CI excludes 0 in the claimed direction; this is inference conditional on the frozen lists. Reported by access stratum
(`direct`, `non_direct`, and `proxy` for "beyond any input") in `results/benchmark/p1_by_access.csv`: each arm is
re-ranked within Pool_G ∩ stratum with k_s = min(100, ⌊|Pool_s|/3⌋) (in `p1_pool.json`), same bootstrap and branches.
If k_G or any k_s < 20 the result is "descriptive (underpowered: k)" with no claim. EN rows with selection frequency 0
are not rankable; if fewer than k remain the list is shorter and `n_list` is reported.
- **Power** (a property of (g, d), identical across arms, committed before the tag): β_D from S1-obs, standardised
  (β* = β·SD_D(d̃)), shrunk by an empirical-Bayes normal scale mixture fitted once on all U rows (incl. F2) and
  back-transformed. SE_R_pred = σ̂_resid,D(g)·√(1/n1_R + 1/n0_R) (binary) or σ̂_resid,D(g)/(SD_R(d̃_R)·√n_R)
  (continuous), from R molecular data only. Power = 1 − F_t(t_0.95,ν_R; ncp = |β_shrunk|/SE_R_pred) with ν_R from §5.3
  on R carriers, so the test and its power share one df rule.

6.2 **Secondary.** U_all with the same rules (F2 still excluded); top-k incl. F2; k = 50/200/500; top-k over all of U
(non-evaluable = not replicated); gene-matched (≤ 2 features per gene); overlap strata EN-pred∩S1-obs, EN-pred-only,
S1-obs-only; set differences; sign concordance, Deming slope (`scipy.odr`, sx = SE_D, sy = SE_R), observed/expected;
S1-pred contrasts; R non-Avana subset.

6.3 **XLIB (secondary).** S1-obs and EN-obs only (the valid contrast): discovery = D lines whose library string is
exactly Avana; replication = D lines with no Avana in the string; others in neither. Same pool, k and bootstrap rules.
Caveat: U and eligibility were selected on all-D effects, so absolute XLIB rates are optimistic; sensitivity
`u_rule = discovery_only` recomputes eligibility on discovery lines (no refit). `replication_xlib.csv` has `u_rule`.

## 7. Null calibration and negative controls

7.1 **Fixed-threshold null (primary).** One within-(lineage × library) permutation of line labels per draw, applied
jointly to every 𝒟 row (shared, gene-specific and NA indicators); targets, predictions and covariates (incl. own CN,
EXPGATE) stay fixed; strata < 2 lines merge into the lineage's largest stratum. Each draw re-runs screen and stability
for both EN arms and S1 for both S1 arms on the 300-gene subset (`null_gene_subset_seed` → `configs/null_genes.txt`),
at the observed operating points: per-gene, per-stage λ_1SE, K, frequency/sign rules, the observed F1 BH p-cutoff
p*_F1 (from the full observed family), the S1-pred cutoff, and the MD-strict cutoff p*_sc with the same-sign rule.
FDR_emp_fixed(list) = mean_b |hits_b| / max(1, |hits_obs|) on the subset, bootstrap CI over draws, plus P(|hits_b| ≥ 1)
and false cross-lineage rate. **Pass: FDR_emp_fixed ≤ 0.10 per thresholded list**; failure is reported as
miscalibration, never retuned. B = 20, raised to 50 iff 50 × wall(draw 1) ≤ 180 min (timing only); per-draw seed
`permutation_seed + b`; 6 workers. S1 arms additionally get 200 vectorised fixed-threshold permutations.
7.2 **Selection-adaptive null (secondary).** Full re-run including λ CV and BH (`FDR_emp_rerun`), B = 10 within the same
cap after the primary draws; reported, not a pass rule. `results/null/summary.json` holds `fdr_emp_fixed`,
`fdr_emp_rerun`, `p_cut_F1`, `p_cut_sc`, `lambda_source = observed`; a script recomputes p*_F1.
7.3 **Negative controls** (`negative_controls_s7`): NC1 FDR_emp_fixed ≤ 0.10 (pre-declared target, never replaced);
NC2 100 matched-SD low-tier genes (`lowtier_seed`): median fidelity, MD-strict hits per gene; NC3 self-null features:
fraction of testable pairs whose carriers are dependent; NC4 fold artefact (toy plant); NC5 false-sensitivity rate on
resistance/N CGT tests; NC6 lineage-TF genes in class 3. NC2–NC6 provisional bounds are replaced by the max (min for NC6)
of the toy statistic over 5 toy replicates (`toy_seed + r`, r = 0..4) at toy dimensions, written to
`results/toy/nc_bounds.json` and frozen at `freeze-bounds` with the provisional values alongside; if the replicates
exceed the `toy_run_bounds` timebox, the provisional values become final (logged in `deviations.csv`). Outside its
bound = failed control.

## 8. Subtypes and confounder audit

8.1 **Per-gene subtype** S(g, d), d ∈ MD-strict: accurate-dependent lines of g whose line share
s_ℓd = β_gd(x̃_ℓd − x̄_d) / Σ_{selected d′} β_gd′(x̃_ℓd′ − x̄_d′) (both stages, FWL-residualised features, full-D EN-pred
fit; lines with negative sum only) is ≥ 0.25; binary d also requires carrier status. Multi-label allowed; dependent but
mispredicted lines are "unexplained".

8.2 **Biomarker module** M(d): genes with (g, d, same sign) ∈ MD-strict; member lines = carriers accurate-dependent on ≥
50% of them; size as genes and effective gene clusters (D dependency |r| > 0.5, average linkage).

8.3 **Within vs across lineage** (MD-strict ∪ controls ∪ S1-obs BH top 200). Testable lineage: ≥ 10 D lines, ≥ 3
carriers, ≥ 3 non-carriers; powered: per-lineage power ≥ 0.5 at one-sided α 0.05 under the §6.1 power rule. Pooled
significance comes from the confounder model restricted to testable lineages (CR2 + ν), BH over F4. Paule–Mandel τ²
with HKSJ CI (`statsmodels combine_effects`, method_re="iterated"), I², τ and prediction interval are descriptive.
Labels, in precedence order: **cross-lineage** (pooled q < 0.05 ∧ ≥ 2 powered lineages in ≥ 2 organ-system groups ∧ same
sign in ≥ 80% of powered lineages ∧ sign kept under leave-one-lineage-out for every contributing lineage) >
**lineage-restricted** (≥ 2 testable lineages ∧ per-lineage BH q < 0.05 in F5 ∧ d × lineage interaction p < 0.05) >
**lineage-unresolvable** (< 2 powered lineages; "single-lineage carriers" when only one lineage has carriers) >
**undetermined**. `tests/test_cross_lineage_labels.py` shows each label reachable (a 2-lineage, 2-group effect β = −1 with
15 carriers each is cross-lineage) and exactly one label per case. **Q1b**: Fisher OR (95% CI) of R replication and of
CGT membership, cross-lineage vs lineage-restricted MD-strict determinants; a null result is a finding.

8.4 **Co-location.** For every CN/SIL determinant in a locus group, refit on lines discordant between d and each
neighbour; ≥ 5 discordant carriers and same sign at p < 0.05 ⇒ gene-level name, else locus-level
(`results/locus/resolution.csv`).

8.5 **CN-cutting artefact** (AMP determinants whose group contains g): (i) carriers-minus-non-carriers median effect of
unexpressed amplicon genes (log2TPM+1 < 1); (ii) within carriers, r(g effect, g expression), negative-significant at
α 0.05. "Cutting artefact" if (i) < −0.2 and not (ii); "amplified-driver addiction" if (i) ≥ −0.2 and carrier median
effect of g < −0.5; else "unresolved". No literature exceptions.

8.6 **Culture.** Cramér's V of subtype membership vs CellFormat, 3D, serum-free, medium family and library; V ≥ 0.3 flags
the subtype (→ class 6 unless a within-stratum test keeps the sign at p < 0.05). Module genes are tested for overlap
with the culture gene sets (hypergeometric BH < 0.05).

8.7 **Hypermutation and TP53.** Hypermutated = MSI ∪ top 5% TMB. HOT/ALLELE/LOFBI determinants must keep sign with
|β| ≥ 50% of full in non-hypermutated lines; SIG determinants must keep sign within TP53-mutant lines. Failures → class 6.

## 9. Positive controls and the gate

9.1 **Files.** `configs/positive_controls_spec.csv` (family, dependency_symbol, feature_primary, feature_secondary,
fallback_feature, expected_sign, citation, note; frozen at `freeze-design`) → `scripts/build_positive_controls.py` →
`configs/positive_controls.csv` (adds control_id, dependency_entrez, scope ∈ {gate, class_check, untestable}, in_cgt,
n_carriers_D/R, access_class, en_max_attainable_freq, gate_eligible), frozen at `freeze-bounds` and an ancestor of the
first commit writing `results/md/`. A spec id `SIL:<x>` expands to every SIL-DEL/SIL-MUT/SIL-METH/SIL-other:<x> feature in 𝒟;
the control is recovered if any expanded feature is recovered with the expected sign. `build_positive_controls.py` fails on
any spec id (after expansion) that is not a `dictionary_v1.json` id.

9.2 **Controls** (S = β < 0 for feature presence; W = β > 0):
- ADD-HOT (S): KRAS←HOT:KRAS; NRAS←HOT:NRAS; BRAF←ALLELE:BRAF_V600 (HOT:BRAF as labelled secondary proxy);
  PIK3CA←HOT:PIK3CA; CTNNB1←HOT:CTNNB1; EGFR←ALLELE:EGFR_activating (sec. HOT:EGFR; AMP:EGFR "non-clinical").
- ADD-CN (S): ERBB2←AMP:ERBB2; MYCN←AMP:MYCN; CCND1←AMP:CCND1; PKMYT1←AMP:CCNE1; PSMC2←LOSS:PSMC2 (CYCLOPS,
  "processing-limited").
- ADD-EXP (S): ESR1←EXPDOSE:ESR1 (clinical expression biomarker).
- PAR (S): SMARCA2←LOFBI:SMARCA4; ARID1B←LOFBI:ARID1A; STAG1←LOFBI:STAG2; VPS4A←LOSS:VPS4B; VPS4B←LOSS:VPS4A;
  ENO2←LOSS:ENO1; PAPSS1←LOSS:PAPSS2; PRMT5←DEL:MTAP; MAT2A←DEL:MTAP; ME3←LOSS:ME2; CREBBP←LOFBI:EP300;
  EP300←LOFBI:CREBBP; DDX3X←YABS:DDX3Y (secondaries in the spec file).
- LOF-SL (S): PIK3CB←LOFBI:PTEN; NFE2L2←LOFBI:KEAP1. SIL (S): WRN←SIL:MLH1; UBC←SIL:UBB.
- SIG (S): WRN←SIG:MSI; KIF18A←SIG:WGD. PROG (S): ADAR←PSCORE:ISG; NFE2L2←PSCORE:NRF2.
- P53 (W on COMP:TP53MUT): MDM2, MDM4, PPM1D, USP7. W: CDK4, CDK6, CCND1 ← LOFBI:RB1.
- LINEAGE-CHECK (class_check): the 18 lineage-TF rows of P1 `biomarker_pairs.csv`, expected in class 3.

9.3 **Gate.** `results/controls/pc_eligibility.csv` (frozen at `freeze-bounds`) fixes the denominator. Recovered = (g,
primary or secondary feature) in the list with the expected sign (EN: selected; S1: BH; MD-strict: catalogue). Pass
(decided on MD-strict; all arms reported, exact binomial CIs): ≥ 60% of gate-eligible controls and ≥ 50% within each
family with ≥ 2 eligible rows. Also reported: recall by access class, PAR recall on rows with a non-direct encoding, and
`en_unreachable` rows separately. Run once: `results/gates/G_pc_attempt1.json` (sha256 logged, never overwritten); the
verdict never chooses between variants. Afterwards only bug fixes, each a row of `results/gate_iterations.csv` (commit,
reason, failing-then-passing toy test, gate before/after), applied to all arms; commits touching `src/` must touch
`tests/`. With > 1 row the CGT headline is gate-disjoint and REPORT.md states "gate failed at attempt 1 (x/y); passed at
attempt n after Z".

## 10. Clinical ground truth (CGT)

10.1 **Build.** `scripts/build_cgt.py` reads frozen CIViC, MOAlmanac, FDA CDx and cached CTG records plus
`configs/cgt_mapping.csv` (curated; every row ≥ 1 resolvable CIViC EID/AID, MOAlmanac id, FDA PMA/BLA/510(k) or NCT;
`clinically-refuted` rows may instead cite a PubMed-verified PMID, `source_type = pmid_only`) and
`configs/cgt_corrections.csv`. Unresolvable rows → `results/clinical/cgt_dropped.csv`. `scripts/build_cgt_features.py`
(molecular data only) applies `configs/cgt_allele_rules.json` and `configs/fusion_known_positives.csv` (cell-line names
resolved via Model.csv) and writes `configs/cgt_features.csv`, which also supplies the ALLELE and CLIN 𝒟 families and
therefore runs before the first dictionary build. Neither script reads effects or predictions (static import scan).
Both outputs are tagged `freeze-cgt` before any dependency data are joined.

10.2 **Scoring unit, sign, partitions.** A CGT *test* = (feature, target set T, comparator, indication set); drugs
mapping to one test are collapsed. Comparator `within-indication` or `pan-cancer + lineage covariate`. Sign S ⇔ β < 0 on
the target's effect for feature presence, W ⇔ β > 0; absence biomarkers coded as presence of the mutation with flipped
sign. Partitions: `sensitivity`, `bypass-resistance` (scored W), `acquired-resistance` (listed, untestable),
`expected-null`, `clinically-refuted`. An S hit on a W test is a sign error. Columns: `test_id, feature, feature_D,
target_set, comparator, indication_codes, sign, partition, ko_class, tier, drugs, source_ids, source_type, in_gate,
tumour_agnostic, excluded_indications, allele_level_required, proxy`; (feature, target_set, comparator) unique.
`feature_D` is the 𝒟 id the arms can rank; tests without one are a P2 miss category.

10.3 **Content rules (before `freeze-cgt`).**
- Corrections (`cgt_corrections.csv`): KIT D816V avapritinib/midostaurin S, imatinib acquired-resistance; del(17p)→BCL2
  population-only N; BRAF sensitivity requires V600; EGFR amplification is not a CGT test; EGFR feature = activating
  alleles; PTEN→AKT target {AKT1, AKT2, AKT3}; RAS/NRAS WT/mutant→EGFR merged into one W test (CLIN:RAS_MUT, comparator
  RAS/BRAF-WT CRC); ESR1 LBD N, untestable; CDKN2A±MTAP merged into the MTAP test with the §8.4 resolution test.
- **Proxy rule (shared with §9.2):** allele-defined tests use ALLELE features (`proxy = 0`). A labelled secondary copy
  with `HOT:<gene>` and `proxy = 1` exists only where the gene-level call approximates the allele; it never enters the
  headline. Never proxied: KRAS G12C vs other KRAS, EGFR T790M/C797S, KIT D816V vs exon 11, ERBB2 TKD.
- Features: fusions with canonical partners, known positives must be carriers (NCI-H3122, NCI-H2228 ALK; K562
  BCR::ABL1; MV4-11 KMT2A; MDA-MB-175-VII NRG1), ALK carriers ≤ 10; BRCA_BIALLELIC = dosage 2, or damaging with AF ≥ 0.7
  or own CN < 0.75, or SIL BRCA1, excluding hypermutated lines; ER_POS = `ModelSubtypeFeatures`, else a 2-component
  Gaussian mixture on ESR1 within Breast ∪ Uterus ∪ Ovary (posterior 0.5); HER2_POS = AMP:ERBB2 ∨ "HER2+"; MTAP-loss = DEL
  ∨ SIL; MSI per §1.5.
- KO classes: `concordant` (incl. allele-specific inhibitors of addicted oncogenes vs WT); `concordant-family` (|T| ≥ 2,
  pre-registered members; hit if any member; ceiling Bonferroni × |T|); `concordant-KO/combination-required`;
  `discordant-mechanism`; `population-only`; `non-cell-autonomous`; `untestable`. Combinations scored on the anchor
  drug's target.
- Tiers (exclusive): T1 approval with label biomarker; T2 CIViC A/B Supports or NCCN, not approved; T3 trial-stage
  biomarker-selected NCT only; refuted. Forward-looking T3 rows (CCNE1→PKMYT1, ATM→ATR, SMARCA4→SMARCA2, BRCA→POLQ,
  BRCA→USP1) only with a matching cached CTG record.

10.4 **Data ceiling** = S1-obs on that test: comparator-restricted §5.3 model (Wilcoxon robustness), one-sided in the
expected direction, BH q < 0.1 (F7), for S/W tests with ≥ 3 D carriers; labels `data-supported`, `data-contradicted`,
`data-null`, `untestable` (`cgt_ceiling.csv`; recomputed from numbers). Concordant data-null/contradicted tests feed §12
class 8.

10.5 **Q4 metrics** (`results/clinical/q4_summary.json`).
- **Q4a / P2**: recall of testable S/W tests (concordant classes; testable = `feature_D` in U with ≥ 5 carriers) in each
  arm's top-200 over U_all (secondary 100/500; hit = g ∈ T and d = `feature_D`, any sign); exact McNemar between arms;
  sets `non_proxy` (headline) and `incl_proxy` (secondary), each on `all` and `gate_disjoint`; per tier with exact CIs;
  trend test only if every tier has ≥ 5 tests. If fewer than 15 testable headline tests, P2 is "descriptive
  (underpowered: n)" (`p2_status`). Also recall over data-supported tests and `miss_decomposition` per arm (ineligible
  gene class 1/2/3; no `feature_D` or feature not in 𝒟_g; eligible and ranked below k; `lost_to_group_member`; ceiling
  data-null), summing to the misses. A test is hit iff g ∈ T and `rank_of_feature_D` ≤ k, where `rank_of_feature_D` = 1 +
  the number of de-duplicated rows outside (g, variant group of `feature_D`) with a better rank_stat; the same rank is
  used in `known_recovery.csv`. **Primary contrast:** EN-pred vs S1-obs on `non_proxy` (`gate_disjoint` if > 1 gate
  attempt), exact two-sided McNemar at α 0.05, three branches (`primary_endpoints.json`), not evaluated when
  descriptive; per-arm recall with Clopper–Pearson 95% CI; other pairs secondary and unadjusted.
- Direction accuracy (hits with the expected sign / all hits) and sign-error rate for every arm and the ceiling;
  N-partition hit ratio vs 1,000 matched random pairs (flag > 2); false-sensitivity rate.
- **Q4b actionability**: pair- and feature-level fractions of MD-strict vs S1-obs and vs 200 within-list g-label
  permutations; PROG/PSCORE excluded; line coverage by lineage.
- **Q4c population** (`q4c_population.csv`): within-indication β; per non-indicated lineage (≥ 3/≥ 3) β; on-label /
  beyond-label / off-label-only / indication-untestable; population precision; carriers-only OR of accurate-dependent
  membership, indicated vs not; tumour-agnostic tests must be cross-lineage; exclusions reported.
- **Q4d**: recovering discordant or non-cell-autonomous tests is a flag. **Prevalence (Tier B)**: D vs TCGA (cBioPortal,
  cached), Wilson CIs, log2-ratio flag > 1.

10.6 **Clinically refuted set** (expected KO sign S; ≥ 3 rows): KRAS-mutant NSCLC → MEK1/2; CCND1 amp/CDKN2A loss →
CDK4/6 as selection biomarker; BRAF V600E CRC → BRAF monotherapy; EGFR-amplified GBM → EGFR TKI. `recall_refuted` vs
`recall_approved` (exact CIs); similar rates are a Q5 finding.

## 11. Combinatorial screens (Tier A, 45 min, never cut)

For MD-strict and S1-obs candidates "loss of A → dependency on paralog B", look up A–B in the downloaded Dede 2020,
Thompson 2021 and Esmaeili Anvar 2024 tables and the Flister et al. bioRxiv supplement (analysis only). Parrish 2021
(counts only) is not reprocessed; predicted-SL resources are excluded. Labels from file rows only: replicated (hit in ≥ 2
screens or ≥ 2 lines), tested-not-replicated, untested. Validation rate of MD-strict vs S1-obs paralog candidates vs a
screened background matched on sequence-identity decile and family size (Fisher, CI); supporting evidence
(`results/combinatorial/lookup.csv`).

## 12. Limits, candidate cards

12.1 **Gene classes** (`results/limits/gene_classes.csv`; universe = variable ∪ selective genes ∪ controls): gene-level
first — 1 no between-line contrast (D SD < 0.15); 2 not predictable (low tier); 3 accurate only via lineage/library
(`r_within` < 0.10) or unnameable (fidelity < 0.10 or M1 named share < 0.5); 3b predictable and nameable but no
MD-strict determinant. Otherwise 9 if any pair is class 9, else its furthest pair class (max); all are listed.

12.2 **Pair classes** (`pair_classes.csv`, MD-strict pairs): 4 expression gating only, or own EXPDOSE without AMP/SIL
backing unless the gene is in the clinical expression list (ESR1, ERBB2, AR, PGR, FOLR1); 5 too rare (< 5 D carriers)
or R power < 0.5; 6 confounded (culture flag, library, cutting artefact, unresolved co-location, hypermutation, TP53,
subtype-confounded, collinear_with_design); 7 powered but not replicated; 8 KO–drug discordance; 9 valid candidate.
Extra labels: predictor-blind (S1-obs hit, replicated in R or XLIB, not EN-pred-selected); predictor-only not replicated;
processing-limited (own-arm CN, CYCLOPS); pan-lethal by CRISPR (Krill-Burger 2023). Methylation is assayed only for the ~843 RRBS-profiled CCLE lines (D7); its coverage is reported.

12.3 **Candidate cards** (`scripts/build_candidate_cards.py`, frozen outputs only):
`results/candidates/known_recovery.csv` (one row per control and testable CGT test: rank and sign per arm, R result,
ceiling label, access class) and `results/candidates/shortlist.csv` (class-9 pairs replicated in R, ranked by
cross-lineage label, then R one-sided p, then number of powered lineages; columns g, d, access_class, lineages and
per-lineage β, carriers D/R, delta flags, locus resolution, combinatorial label, CGT/actionability, set-A dependent
fraction as a normal-line flag, and `not_in_curated_lists` = not in controls ∪ CGT ∪ P1 biomarker_pairs ∪
oncogene_pairs). REPORT.md opens with both tables (top 15 shortlist rows) and three case studies = the top three
shortlist rows; the word "novel" is never applied to a gene. A test asserts shortlist ⊆ (class 9 ∧ R replicated).

## 13. Verification and anti-drift

13.1 **Toy simulation** (`src/biomark/toy.py`, `tests/test_toy_pipeline.py`), written by a different agent than the
pipeline author (commits carry an `Agent: <name>` trailer; checked by trailer). Same entry points and
`design_params.json` (hash asserted); plants and `toy_seed` frozen at `freeze-design`. 600 lines × 12 lineages, 300
genes, 400 features with correlated blocks (r ≈ 0.9), library nested in 3 lineages; predictor = ridge trained by 5-fold
CV on noisy toy observations plus per-fold offsets. Pass bounds (`toy_s13`): pan-lineage recall ≥ 0.8; own-locus CN
(own-CN covariate dropped by rule) ≥ 0.8; MSI-like signature ≥ 4/5 with passenger LOF ≤ 1/5; co-located passengers at
locus level ≥ 4/5; culture-confounded in MD-strict ≤ 1/5; lineage-only genes 0 in MD-strict, ≥ 80% class 3; fold-offset
feature ≤ 1/5 after centring and ≥ 3/5 without; genomic feature with expression proxy recovered in stage G ≥ 4/5; a
6-carrier strong plant selected in > 50% of subsamples; wrong-sign CGT pair = sign error. Error rates: on the planted toy
FDP_true ≤ 0.15 per list and FDR_emp_fixed ≥ FDP_true − 0.05; on a pure-null toy (`toy_null_seed`) P(|hits| ≥ 1) ≤ 0.10
for BH lists, EN mean hits reported. Five replicates of the toy run give the NC2–NC6 bounds (§7.3).

13.2 **Unit tests**: `test_inference_calibration.py` (20,000 null simulations per n ∈ {255, 1022} × carriers ∈ {3, 5,
10} × {homoskedastic, carrier SD × 2} on the real covariate designs: size at α 0.05 in [0.035, 0.065], at α 1e-3 ≤
2e-3); `test_cr2.py`; `test_eb_mixture.py` (3-component known answer, r > 0.99 vs oracle); `test_bootstrap_clusters.py`;
power df = test df; power invariant to feature rescaling; cross-lineage labels; culture map; leakscan (planted `"KRAS"`
and `0.7` give two hits); subsamples; permuted NA flags aligned; dictionary/locus counts; self-null removal; covariate
drops (ERBB2 × AMP:ERBB2 drops own CN; WRN × SIG:MSI drops MSI; a lineage-dummy d is flagged); FWL; determinism;
additivity; NMF; CGT schema; EN-pred blind; guard refusals.

13.3 **Independent reimplementation**, started right after `G_pc_attempt1.json` by a separate agent in
`tests/independent/`: S1-obs p for 20 genes (1e-6), BH counts and gate recall (exact), EN selection for 20 genes
(frequency within 2/50; identical selected set excluding features within 0.04 of 0.6; β within 1e-3 relative), deltas for
10 pairs → `results/gates/G_indep.json`; tolerances read from `design_params.indep_s13`.

13.4 **Two freeze tags and deviations.** `freeze-design` (now): design_params (NC2–NC6 `provisional`), primary_endpoints,
tiers, seeds, build_params, positive_controls_spec, cgt_allele_rules, fusion_known_positives, cgt_corrections, culture and
lineage lists, toy plants. `freeze-bounds` (end of M3b): NC bounds from the toy run, `pc_eligibility.csv`,
`positive_controls.csv`, `null_genes.txt`, `subsamples.json`, `dictionary_v1.json`, `data_status.json`.
`configs/frozen_hashes.json` has sections `design` and `bounds`. `freeze-bounds` must be an ancestor of the first commit
touching `results/benchmark/`, `results/md/`, `results/gates/` or `results/null/`. A later change to a frozen file needs a
decision `D<n>` in RESEARCH_TREE.md naming the file and whether dependent outputs were seen; triggered fallbacks write
`results/deviations.csv`. **Leakscan** = AST scan of `src/` and `scripts/` for (i) string constants equal to an HGNC
symbol of ≥ 3 characters from `gene_map` or matching `ACH-\d{6}`, (ii) numeric constants outside {0, 1, −1, 2, 0.5,
100} in `src/biomark/` other than the params loader, unless the line carries `# literal-ok: <reason>`.

13.5 **Report numbers and claims.** `scripts/check_report_numbers.py` (ported from Part 1) `--strict` at the final gate;
tagged sources must be tracked and resolve to registry rows with status ok (`exploratory`, `dryrun` uncitable). `results/claims.csv` (claim id, class incl. `beyond_inputs`,
evidence, stratum, rule) is generated from §6, §8.3 and §12; a sentence applying a controlled word (validated,
replicated, cross-lineage, valid biomarker, recovers, significant) to a named gene or pair carries a matching
`<!--c:claim_id-->`; a `beyond_inputs` claim needs `stratum = non_direct`. `.gitignore` ignores `*.gz`/`*.parquet`, so
committed results are plain `.csv`/`.json` ≤ 5 MB; bulk tables stay in `cache/` with sha256 in
`results/bulk_manifest.csv`; `gate.py tracked` asserts every registry output is tracked.

13.6 **Ledger, registry, gates.** `configs/design_ledger.csv` (one row per deliverable and accepted review item) is
generated from both response files; no row `planned` at the final gate. `results/registry.csv` (file-locked: run_id,
script, commit, dirty, config/design_params/data_status hashes, status, outputs). Tracks use own output folders and
D-id ranges. `gate.py` subcommands: `leakscan, frozen, data, controls, G_pc, params` in M3b; `tracked, cgt, ledger,
G_rep, G_indep, part1` by M7. `gate.py params` asserts every §F key exists and is non-null. **Procedural (not
script-checked):** catalogue agents do not read P1 REPORT.md test sections or P1 RESEARCH_TREE.md.

## 14. Tiers, compute, schedule (≈ 16 h left; `configs/tiers.json`)

- **Tier A (must):** 𝒟, four arms, MD-strict, confounder audit, gate, fixed-threshold null, subtypes and modules,
  cross-lineage, co-location/cutting/culture/hypermutation, component readout, CGT recovery/direction/actionability/
  population, combinatorial lookup, replication in R, XLIB (2 arms), limits, candidate cards. **Tier B:** M1-proxy
  attribution panel; TCGA prevalence; selection-adaptive null. **Cut:** LMM, EN-pred-clean, M3 refit, clustering,
  CTG harvest, cross-fit, separate bound-simulation script.
- **Compute:** foreground 6 workers while the null runs on 6 (else 12), 1 BLAS thread each, peak RSS ≤ 12 GB. Each job writes `results/timing/<job>.json`; `tiers.json.stage_job_map` maps stages to job
  stems; `scripts/budget_check.py` runs at every milestone end, fails on unmapped stems or > 12 concurrent workers, and
  when cumulative lateness exceeds 1 h applies the next unused `cut_order` entry, logged in `results/deviations.csv`.
- **Schedule (h):** M2 v3 + `freeze-design` 0.5; M3b dictionary, covariates, eligibility, controls, CGT features, guard/
  gate, toy → `freeze-bounds` 3.0; M4 arms, MD-strict, gate, null and XLIB in background 3.0; M5 subtypes, audit,
  attribution 2.0; M6 CGT scoring, combinatorial 2.0; M7 `freeze-catalogue`, single R run 0.75; M8 limits, cards,
  report, review 2.75. Total 14.0 + 2.0 reserve.

## §F Frozen parameters

Each row maps 1:1 to a key in the named file; code reads only the file.

| § | Parameter | Value | File : key |
|---|---|---|---|
| 1–2 | data deadline, MSI; eligibility; dependent | first dictionary build, MSIScore ≥ 20; tier high/medium ∧ r_within ≥ 0.10, U_ctrl; −0.5 (−1.0), relative 25th pct, SD < 0.15 | design_params : data_s1, eligibility_s2, dependency_s2 |
| 3 | dictionary | binary ≥ 5, continuous ≥ 20; AMP > 2, DEL < 0.25, LOSS < 0.75; SIL; groups (Jaccard 0.5, 10 Mb); cap 3,000; NMF 50; spectra | design_params : dictionary_s3 |
| 4 | EN; ranking | l1 0.9; FWL design; stages G/X; K 300 (150); 25 λ; 1-SE fixed; 25 pairs; freq 0.6; sign 0.9; in-fit carriers 2; mean β over 50 | design_params : surrogate_s4 |
| 4 | attribution | M1-proxy ≤ 300 genes; additivity 1e-5; named share 0.5 | design_params : attribution_s4 |
| 5 | covariates, SE; S1; MD-strict; families | droppable set; R² > 0.49; CR2 + ν; leverage; VIF 10; BH 0.05; q_sc 0.05; F1–F7 | design_params : confounders_s5, s1_s5, md_strict_s5, test_families_s5 |
| 6 | P1; XLIB | U_G, k_G = min(100, floor(pool/3)), power 0.5, one-sided 0.05, bootstrap 1,000, strata; 2 arms, u_rule | primary_endpoints : P1_primary, XLIB_secondary |
| 6 | power; R rules | standardised EB; noncentral t; SE_R; NaN rule | design_params : replication_s6 |
| 7 | null; controls | fixed thresholds, 300 genes, B 20 (→ 50), 6 workers, FDR 0.10, rerun B 10; NC1 fixed, NC2–6 from toy | design_params : null_s7, negative_controls_s7 |
| 8 | subtypes, lineage, cutting, culture | share 0.25; module 50%; 3/3; power 0.5; q 0.05; 2 groups; 80%; −0.2; V 0.3 | design_params : subtypes_s8 |
| 9–11 | gate; CGT; P2; combinatorial | 60% / 50%; ceiling ≥ 3, q 0.1, testable ≥ 5; k 200, ≥ 15 tests; 45 min | design_params : gate_s9, clinical_s10, combinatorial_s11; primary_endpoints : P2_primary_clinical |
| 12–14 | limits; toy; independent; compute, tiers; seeds | aggregation, 3b; plants, bounds; tolerances; 6 + 6 workers, budgets | design_params : limits_s12, toy_s13, indep_s13, compute_s14; tiers.json; seeds.json |
