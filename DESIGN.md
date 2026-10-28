# DESIGN — Predictor-derived dependency subtypes and their molecular biomarkers

Version: v2 (after review round 1). Status: not frozen (round 2 pending).
Inputs: PROMPT.md, docs/part1_brief.md, LITERATURE.md, reviews/round1_*.md, reviews/round1_response.md. Decision
numbers live in §F and its config files; prose numbers are descriptive.

## 0. Questions, claims, endpoints

- **Q1 subtypes.** For dependencies the Part 1 predictor predicts accurately *beyond lineage and library*, can its
  reasons be read out as assayable molecular determinants that group lines dependent for the same reason, within and
  across lineages? **Q1b**: do cross-lineage determinants replicate, or match clinical biomarkers, more often than
  lineage-restricted ones (§8.3)?
- **Q2 benchmark.** Does the predictor add information beyond the observed data (EN-pred vs EN-obs), and does sparse
  multivariate selection add beyond univariate tests (EN-obs vs S1-obs)?
- **Q3 confounders.** Which candidates survive lineage, library/screen, culture, passenger/hypermutation, TP53/Cas9
  response, co-location, CN/cutting artefacts and multiplicity?
- **Q4 clinical.** Recovery, direction, actionability and patient population against a frozen clinical ground truth.
- **Q5 limits.** Every eligible or control dependency gets a failure class (§12).

**Falsifiable novelty claims** (vs Tsherniak 2017, Dempster 2020, Behan 2019, Pacini 2024): (a) predictor-derived
lists replicate better than, or name determinants missed by, EN-obs and S1-obs; (b) line-level subtype assignment with
a tested cross-lineage label; (c) scoring a dependency catalogue against tiered CIViC/FDA ground truth with a
data-ceiling decomposition of misses.

**Endpoints** (`configs/primary_endpoints.json`): P1 (Q2) = R replication of top-100 per arm (§6.1); P2 (Q4a) =
CGT recall at matched k (§10.5). XLIB within-D cross-library replication is a secondary co-endpoint (§6.3).
Everything else is secondary. Non-goals: a new predictor; outcome or drug-response modelling.

## 1. Data, roles, protection

| Role | Lines | Use |
|---|---|---|
| D (discovery) | 1022 Part 1 dev lines, `cv_fold` 0–4 | everything before `freeze-catalogue` |
| R (replication) | 255 Part 1 test lines | opened once, after `freeze-catalogue`, by `scripts/replicate.py` |
| set A | 15 non-cancerous | §12 limits only (guard scope `limits`) |

1.1 **Predictions.** D: Part 1 OOF headline `m5_avg` (`p1.oof_predictions("m5")`), fold-centred per gene
(`p1.fold_centre`) before any cross-line use; component OOF runs m1/m2n/m3n (their mean equals M5 to 1e-6). R: all-dev
refit `p1.test_predictions("m5_avg")`. Targets: Chronos 26Q3 via `guard.load_effects` only.

1.2 **What R is.** Part 1 opened R's effects once to score its predictor; no biomarker analysis has used R. R shares
the Chronos fit with D (guide efficacies, CN/proximity correction, RCPC), so replication tests generalisation across
lines, not across processing. R replication is also reported for the R lines screened without Avana.

1.3 **Guard (extends `src/biomark/guard.py`).**
- *Single entry point.* Only `guard.py` may contain `crispr.parquet`, `CRISPRGeneEffect`, `results/test`,
  `cache/test`, `CRISPRGeneDependency`, `CRISPRGeneEffectUncorrected`, `ScreenGeneEffect` or call `essml.guard` /
  `io.load_matrix`; only `p1.py` may `import essml`; no `subprocess` targets ../project. `p1.read_parquet` already
  refuses the CRISPR cache and P1 test results; the deny-list adds `P1/cache/test/**` except `preds/` and any
  effect-derived DepMap download. Catalogue agents do not read P1 REPORT.md test sections or RESEARCH_TREE.md.
- *New unblind conditions* (added to the current (a) tag, (b) `PART2_UNBLIND=1`, (c) clean tree, (d) purpose):
  (e) caller is `scripts/replicate.py`; (f) sha256 of every frozen catalogue file matches
  `configs/freeze_catalogue.json` (committed in the tagged commit); (g) the tag target equals the one recorded at the
  first granted access (moved tag ⇒ refusal); (h) after the first access, a request is granted only with the same
  config hash and code commit (reproduction only).
- *Dry run.* Before the tag, `replicate.py --pseudo-R fold0` runs end to end with D fold 0 as stand-in R (catalogue
  reused) into `results/replication_dryrun/`; every verdict branch of `primary_endpoints.json` must be shown reachable
  (unit test with synthetic inputs).
- *Immutability.* After the first real access, `git diff freeze-catalogue..HEAD -- configs/catalogue/
  scripts/replicate.py src/biomark/replication*.py` must be empty.

1.4 **Part 1 untouched.** `scripts/snapshot_part1.py` writes `configs/part1_manifest.csv` (path, size, mtime_ns; all of
../project except `data/`, which is sha256-checked by `verify_inputs.py`) at M3b start; `tests/test_part1_untouched.py`
adds a manifest comparison to its HEAD/clean/pycache checks. It runs at every gate, not only in pytest.

1.5 **Downloads** (deadline: start of M4). `scripts/data_status.py` writes `results/data_status.json` (file, present,
sha256, fallback) once, committed before the first M4 result; every analysis branches on it, and late files feed only
labelled secondary analyses.
- Essential: OmicsGlobalSignatures (MSI, ploidy, WGD, CIN, LoH); OmicsSomaticMutations (alleles, LikelyLoF/AF, TMB,
  spectra); CRISPRConfounders (doubling time); OmicsProfiles only if needed for keys.
- Nice-to-have: OmicsInferredMolecularSubtypes, RRBS. Dropped: CNSegments, CRISPRGeneDependency, Gene.csv.
- Fallbacks: MSI = Cell Model Passports `msi_status` (local `external/probe/cmp_models.csv.gz`, 1,044 models, BROAD_ID
  → ModelID; D coverage reported) ∪ `ModelSubtypeFeatures` "MSI"; unknown is NA, never MSS. No doubling-time proxy is
  invented. Without the mutation file, allele-defined CGT tests are `untestable: allele` and spectra are "not tested".

1.6 **Missing data.** 𝒟 is built from raw P1 caches, never fold-transformed blocks. CN-derived features are NA for the
55 D / 16 R lines without CN. In the EN arms NA → 0 plus unpenalised `CN_absent` / `MSI_unknown` indicators; in S1,
D-stage, subtype and carrier counts NA lines are excluded for that family. Each (g, d) uses lines with non-NaN effect;
`n_obs, n_carriers, n_NA` are written per test.

1.7 **D-only fitting.** NMF basis, z-scalers, SIL lineage medians, locus groups, programme scores and thresholds are fit
on D and frozen; R/prediction-only lines are transformed, never refit. A test rebuilds them with R rows deleted from all
inputs and asserts identical hashes.

## 2. Eligibility and "dependent"

- **Eligible gene** = (variable ∨ selective) ∧ m5 tier ∈ {high, medium} ∧ `r_within ≥ 0.10`, where `r_within` is the
  Pearson r on D between fold-centred OOF prediction and observed effect after centring *both* within lineage ×
  library-group strata (strata ≥ 5 lines; rare lineages pooled). Failing genes go to §12 class 3. Positive-control
  genes are analysed regardless and reported separately. Output `results/eligibility.csv`.
- **Associations use continuous effects.** "Dependent" is used only for subtype membership: effect < −0.5; for genes
  with frac_below(−0.5) > 0.5 (PRMT5, MAT2A, PKMYT1, …), effect < the gene's D 25th percentile
  (`dependent_rule = relative`). Sensitivity: −1.0.
- **Accurate-dependent** (ℓ, g): dependent ∧ fold-centred prediction ranks within the n_g most negative, n_g = number of
  dependent D lines. In R: prediction ≤ the R quantile at n_g/n_D (target-free, computed before the tag).

## 3. Interpretable feature dictionary 𝒟

Built by `scripts/build_dictionary.py` → `configs/dictionary_v1.json` (feature id, family, gene(s), definition,
n_carriers_D, n_carriers_R, NA count, Cramér's V with lineage, locus group). Feature ids are `FAMILY:<entrez or name>`.
Restricted gene set Γ = OncoKB cancer genes ∪ Part 1 top-3 Ensembl paralogs of eligible genes ∪ positive-control
partners ∪ known-silenced list (§F). Binary features need ≥ 5 D carriers.

| Family | Definition | Scope |
|---|---|---|
| HOT | hotspot > 0 (gene-level) | all |
| LOFBI | biallelic LOF: dosage 2, or ≥ 1 with own CN < 0.75 or SIL, or X-linked dosage ≥ 1 in males | Γ |
| AMP / DEL | relative CN > 2 / < 0.25, gene-level | AMP: Γ ∪ own gene; DEL: all |
| LOSS | relative CN < 0.75 (hemizygous) | gene-specific: own gene and paralog/partner priors |
| FUS | DepMap filtered fusion, FFPM ≥ 0.1, gene in `configs/fusion_drivers.csv`, kinase as 3′ partner (5′ for FGFR2/3) | drivers |
| SIL | log2(TPM+1) < 1 ∧ gene expressed in ≥ 50% of the line's own lineage (lineages ≥ 10 D lines; else pan-D) ∧ value < lineage median − 2; mechanism in order SIL-DEL (CN < 0.25), SIL-MUT (dosage ≥ 1), SIL-METH (RRBS only), SIL-other | Γ; Y-linked genes defined within males |
| YABS | Y-paralog absent: female, or male with SIL of the Y gene (e.g. DDX3Y for DDX3X) | gene-specific |
| EXPDOSE | own expression among expressed lines (log2TPM+1 ≥ 1; 0 elsewhere, centred) | gene-specific; the only own-expression determinant |
| EXPGATE | own expressed vs not | covariate only (§5.3), never a determinant |
| EXPPAR | expression of top-3 paralogs | gene-specific |
| PROG | 50 NMF factors (5,000 most variable genes, `nndsvda`, cd, 500 iter, `nmf_seed`), annotated by Hallmark/lineage enrichment | all; "expression programme (no clinical assay)" |
| PSCORE | ISG (Hallmark IFN-α, local P1 copy) and NRF2 target score | all |
| SIG | MSI; ploidy, WGD, CIN, LoH (GlobalSignatures; else `C_instab` as CIN proxy labelled proxy, others not done); DMG_BURDEN = log1p damaged-gene count (renamed from "TMB"); TMB only from the mutation file | all |
| SIG-SPEC | mutation file only: SPEC_MMR, SPEC_POLE, SPEC_APOBEC, SPEC_TOBACCO, SPEC_UV (reference-free definitions, §F) with sanity checks (MMR–MSI ρ > 0.5; UV highest in Skin; tobacco Lung > rest) | all |
| COMP | TP53MUT (absence biomarkers are coded as presence of the mutation with the sign flipped), RB_ALT | all |

- **Self-null features** are removed from 𝒟_g and kept as sanity checks: DEL_g, SIL_g, LOFBI_g, and DEL_h for h in g's
  DEL locus group. Their carriers must not be dependent (§7.3 NC3). LOSS_g (hemizygous, CYCLOPS) and AMP_g stay.
- **Locus groups** (no merging): connected components of genes within 10 Mb whose binary calls have Jaccard ≥ 0.5 among
  CN-profiled D lines, per family (AMP, DEL, LOSS), GRCh38 coordinates from Ensembl 116. Groups serve reporting,
  covariate drops, top-k dedup and the resolution test (§8.4). CDKN2A/CDKN2B/MTAP fall in one group and stay separate
  features (r(CDKN2A, MTAP) = 0.79; 128 D lines CDKN2A-del/MTAP-intact).
- **Cap.** Shared 𝒟 ≤ 3,000 features; if exceeded, the largest family's minimum frequency rises to 10, then 20 (logged
  in `results/dictionary/family_counts.csv`). SIL features with Cramér's V(lineage) > 0.5 are flagged as lineage proxies.
- **Not assayed** families get a record with reason in `dictionary_v1.json`: methylation (unless RRBS; then SIL-METH =
  promoter β ≥ 0.7 ∧ SIL), mutational signatures (unless the mutation file arrives), allele-level variants (mutation
  file; HOT stays gene-level in 𝒟). MLH1 SIL-other vs MSI status is the known-positive check of the methylation proxy.
  REPORT.md lists the unassayed families as limits.

## 4. Reading determinants out of the predictor

### 4.1 Hierarchical stability elastic net (EN-pred, EN-obs)
Target: fold-centred OOF prediction ŷ_g (EN-pred) or observed effect y_g (EN-obs); identical code, configs and seeds
otherwise (test asserts config equality apart from the target).
1. **FWL.** Standardise 𝒟 once on D. In every fit (full D, each CV fold, each subsample), residualise target and
   features on the unpenalised block (lineage, library group, `CN_absent`, `MSI_unknown`) estimated on that fit's rows;
   fit with no intercept. Drop columns with residual SD < 1e-6 or < 5 carriers in the fit's rows.
2. **Stages.** Stage G: genomic/signature families (HOT, LOFBI, AMP, DEL, LOSS, FUS, SIL-DEL/MUT/METH, YABS, SIG,
   COMP). Stage X: expression families (EXPDOSE, EXPPAR, PROG, PSCORE, SIL-other) on the stage-G residual. Genomic
   determinants are therefore not absorbed by their expression proxies. Report `fidelity_G` and `dfidelity_X`.
3. **Screen.** Per stage, SIS to the K = 300 features with largest |corr| with the residualised target, inside each fit.
   Forced (still penalised): own EXPDOSE, EXPPAR, own LOSS/AMP, paralog/partner LOSS, YABS, the positive-control
   partner feature.
4. **λ.** l1_ratio 0.9; 25-λ path, eps 1e-2, tol 1e-3; CV over Part 1 `cv_fold`; 1-SE rule on full D; that λ is held
   fixed in all subsamples.
5. **Stability.** 25 complementary pairs of half-samples (50 subsamples), patient-grouped, lineage-stratified
   (`subsample_seed`). Selected = frequency ≥ 0.6 and ≥ 90% of selections with the majority sign. Rank (all of U):
   frequency, then |mean β|, then feature id.
6. **Fidelity** `r2_oof_beyond_lineage_library` = out-of-fold R² of the surrogate (CV folds above) on the residualised
   target. Class-3 threshold 0.10 applies to this value.
7. **Pilot.** Time 50 random eligible genes first (`results/timing/surrogate_pilot.json`); if projected wall > 45 min,
   K = 150. The number of subsamples is never reduced.
8. **EN-pred never reads observed effects.** The EN-pred job runs with `guard.load_effects` patched to raise and an audit
   hook on opening `crispr.parquet` (`tests/test_en_pred_blind.py`). Its eligibility input is a precomputed gene list.

### 4.2 Component readout (Tier A, refit-free)
Full-D fits (no subsampling) on fold-centred ŷ_M1, ŷ_M2n − ŷ_M1, ŷ_M3n − ŷ_M1 → `component_selection.csv`:
"expression-mediated" vs "genomics-layer"; the pre-declared fallback for §4.3.

### 4.3 Direct attribution (Tier B, ≤ 300-gene panel: controls, CGT targets, top MD-strict genes)
- *M1-proxy*: exact linear attribution from saved `t_E_status_24b46e2a` folds (T ≈ M1). Per fold, W in blocks of ≤ 500
  genes, aggregated on the fly into EXP-own, EXP-paralog, PROG_k (least-squares coefficients of w_g on the NMF loadings
  mapped to z-space) and EXP-other; never materialise the tensor. A named share < 0.5 means "M1 unnameable".
- *M3*: refit panel genes on Part 1 dev folds (inner cache redirected to project_b; copying P1's `m2_inner` read-only
  is allowed); `pred_contrib=True`; features mapped E→EXP-own/paralog/other, C→locus, M→LOF, H→HOT, meta→"lineage/meta"
  (not a determinant), `m1`→"M1-mediated". Reproduction: per-gene r ≥ 0.999 vs `g3_m3n` OOF, else report the shortfall
  and continue (one ≤ 30-min diagnosis).
- OOF rule: each line attributed by the fold model that did not train on it; contributions fold-centred. Additivity
  (sum + bias = prediction) within 1e-5 is tested.
- *Agreement*: fraction of panel MD-strict (g, d) whose family is among the top-5 attribution terms of the
  accurate-dependent lines.

### 4.4 Access class (every selected (g, d))
`direct` (d's raw feature is one of M2's 13 inputs for g), `potentially-direct` (in M3's pool: E/C/M/H; verified via
`sel_full_` for panel genes), `expression-derived` (EXP/PROG/PSCORE/SIL), `proxy` (not an input: MSI, SIG, FUS, COMP,
DMG_BURDEN; `carrier_input` names the carrying component from §4.2/§4.3). "The predictor learnt d" is allowed only for
direct and expression-derived; potentially-direct carries a qualifier; proxy reads "learnt a correlate of d (via X)".

## 5. Benchmark arms, MD-strict, confounder model

5.1 **Arms** on the identical universe U = {(g, d): g eligible, d ∈ 𝒟_g (shared + gene-specific − self-null), ≥ 5
carriers among non-NaN lines}: **S1-obs** (univariate confounder-model test on observed effects; the standard),
**S1-pred** (same on fold-centred predictions), **EN-obs**, **EN-pred** (= MD-pure). `results/benchmark/arms.csv` holds
one row per (arm, g, d) with `rank_stat, rank`; a script asserts equal (g, d) universes and identical covariate hash,
subsample count and seeds across arms. Before any top-k, CN features of the same g in one locus group are deduplicated
to the best-ranked member (all arms).

5.2 **MD-strict** = EN-pred selected candidates that pass the D-stage observed-effect screen: confounder model (§5.3),
same sign, BH q_sc < 0.05 within EN-pred selected candidates. Its p and q are *selection-conditional (nominal)*: OOF
predictions are functions of D's effects, so the screen double-dips. Its error rate is the empirical FDR (§7);
confirmatory claims come only from R and external sources. REPORT.md never calls it "significant". EN-obs gets the same
screen (EN-obs-strict) for symmetric reporting. MD-strict is the catalogue for subtypes and clinical claims.

5.3 **Confounder model** (S1 arms, D-stage screen, R tests, ceiling). Built by `scripts/build_screen_covariates.py`
from CRISPRScreenMap → ScreenSequenceMap (no pDNA) → ModelCondition of the *screened* condition (extends
`p1.screen_covariates`), for D ∪ R ∪ set A, reading no effects. Blocks:
- lineage (< 10 D lines → `_rare`);
- library group from the screen map (Avana-only, KY-only, Humagne-only, Avana+KY, other; P1 `lib_group` only for an
  ARI ≥ 0.9 concordance check) + screen-days deviation within Humagne-CD (dropped if VIF > 10; Avana is always 21 d and
  KY 14 d);
- culture: screened CellFormat {adherent, suspension, 3D, mixed/unknown}, SerumFree {T, F, unknown}, base-medium family
  {RPMI, DMEM, IMDM, DMEM-F12, other, unknown}. The parser picks the column by content because Model.csv
  `OnboardedMedia`/`FormulationID` appear swapped; the 69 D+R lines without media (Sanger condition) are `unknown`;
- burden: MSI {MSI, MSS, unknown}, DMG_BURDEN; TP53 status {mut, WT, unknown};
- gene-specific: own CN of g (log2, NA → 0 + `CN_absent`), EXPGATE of g;
- doubling time only if CRISPRConfounders arrives.
Per-(g, d) drop rule, logged in `dropped_cov`: drop any covariate that is d, is in d's locus group, or has |Spearman|
> 0.7 with d on the test lines; omit own CN when d is CN-derived and its locus group contains g (own-CN-adjusted =
sensitivity, feeds §8.5); omit TP53 when d derives from TP53 or MDM2 amplification; MSI and DMG_BURDEN are not both
kept when |ρ| > 0.7. VIF/condition number of the shared design → `results/covariates/vif.csv` (max ≤ 10 after the
declared rules). SEs are patient-clustered (CR1 vectorised for all tests; CR2 for D-stage, ceiling and R tests);
`se_type` is a column of every table.

5.4 **S1 implementation.** Genes grouped by NaN pattern; within a group, QR-residualise Y and 𝒟 on the shared design,
add own CN and EXPGATE by FWL update; t from the partial regression, df = n − rank(Z) − 2, F = t² (= Gaussian LRT);
cluster sandwich from per-line and per-cluster cross-products. Drop tests with < 5 carriers among non-NaN lines or
residual variance < 1e-8. `tests/test_s1_vectorised.py`: 200 random (g, d), incl. NaN-pattern genes and dropped
covariates, vs statsmodels OLS (cluster cov): |Δlog10 p| < 1e-6, identical df.

5.5 **Confounding delta** (MD-strict ∪ S1-obs BH hits ∪ controls only): drop-one-block refits (lineage, library,
culture, burden, TP53, own CN); flag if |β_full − β_minus| > 0.5·|β_minus| and > 1·SE_minus. Sensitivity block
OncotreePrimaryDisease (≥ 5 D lines) → "subtype-confounded". `results/confounders/delta.csv`, flag recomputed by script.

5.6 **Test families** (`design_params.json: test_families_s5`, logged n_F = table rows, q recomputed by script): F1
S1 main; F2 S1 own-EXPDOSE (never counted toward recall of non-expression biomarkers); F3 EN D-stage (screening only);
F4 cross-lineage meta; F5 per-lineage; F6 R replication (no BH; list-level); F7 CGT ceiling. BH under PRDS; Storey q
as sensitivity only.

## 6. Replication endpoints

6.1 **P1 (primary).** Pool = (g, d) ∈ U that are R-evaluable (≥ 3 carriers and ≥ 3 non-carriers with expected
non-NaN effect from their library) and have R power ≥ 0.5. Top-100 of each arm's ranking within the pool. Replicated =
same sign and one-sided p < 0.05 in R (confounder model with the same drops, CR2). Rate differences EN-pred − EN-obs
("predictor adds value") and EN-obs − S1-obs ("multivariate adds value"), each with a paired patient-grouped bootstrap
over R lines (1,000, `bootstrap_seed`; lists fixed, tests refit); claim only if the 95% CI excludes 0 in the claimed
direction, branches as in `primary_endpoints.json`.
- **Power is a property of (g, d), identical across arms**: β_D from the S1-obs model, shrunk by an empirical-Bayes
  normal scale mixture fitted once on all F1 tests (limits winner's curse); SE_R predicted from D residual SD and R
  carrier counts (molecular data only). Computed and committed before the tag; one function hash for all arms.

6.2 **Secondary.** k = 50/200/500; top-k over all of U (non-evaluable = not replicated, counts reported); gene-matched
lists (≤ 2 features per gene); overlap strata EN-pred∩S1-obs / EN-pred-only / S1-obs-only; set differences
EN-pred\S1-obs, S1-obs\EN-pred, EN-pred\EN-obs (`rep_rate, n_pos_ctrl, n_cgt`); list-level sign concordance (binomial
vs 0.5), Deming slope β_R on β_D, observed/expected replications; S1-pred contrasts; R non-Avana subset. Tier B: LMM (FaST-LMM-style: one eigendecomposition of a D-only kernel from non-eligible
genes' dependency profiles, δ grid per gene, GLS) on the top 200 of each list; fraction keeping p < 0.001.

6.3 **XLIB co-endpoint.** Re-run all four arms on Avana-only D lines (own λ CV, same params); replicate on D lines with
no Avana screen; Avana+KY and other mixed lines are in neither (`results/benchmark/replication_xlib.csv`; script asserts
disjoint IDs matching the screen map). Same top-100, power and bootstrap rules. **Caveat:** discovery-line OOF
predictions come from fold models trained partly on replication-line effects, so pred arms are not independent here;
EN-pred vs EN-obs on XLIB is descriptive, EN-obs vs S1-obs is a valid co-endpoint. Tier B adds **EN-pred-clean** from
an M1 refit cross-validated within the discovery lines only (~10 min).

## 7. Null calibration and negative controls

7.1 **Permutation null.** One within-(lineage × library) permutation of line labels per draw, applied jointly to every
row of 𝒟 (shared and gene-specific); targets, predictions and covariates stay fixed; strata < 2 lines merge into the
lineage's largest stratum. Each draw re-runs the full pipeline for all four arms (screen, stability, D-stage filter,
BH) on a pre-drawn 300-eligible-gene subset (`null_gene_subset_seed` → `configs/null_genes.txt`, committed before any
result); draw b uses `permutation_seed + b`. B = 20; raised to 50 only if 50 × wall(draw 1) ≤ 180 min, decided from
timing before any observed output is read. FDR_emp(list) = mean_b |hits_b| / max(1, |hits_obs|) on the subset, with
bootstrap CI over draws, plus P(|hits_b| ≥ 1) and the false cross-lineage rate. **Pass: FDR_emp ≤ 0.10 per thresholded
list** (S1-obs/S1-pred BH, EN selected, MD-strict); failure is reported as miscalibration, never retuned. S1 arms
additionally get 200 vectorised permutations on the same subset (QQ and per-test calibration). Background job, 3 h cap
(`results/null/perm_manifest.json`). Tests: the permuted 𝒟 has the same feature correlation matrix and per-stratum row
multiset; the pipeline function under permutation is the same object as in the observed run.

7.2 **Negative controls with numeric bounds** (`negative_controls_s7`; provisional values replaced by
`scripts/simulate_control_bounds.py` with the toy null's 95th percentile at real D dimensions before `freeze-design`):
NC1 permutation FDR_emp; NC2 100 matched-SD low-tier genes (median fidelity, MD-strict hits per gene); NC3 self-null
features (fraction of testable pairs whose carriers are dependent); NC4 fold artefact (AMI with `cv_fold`, and the toy
fold-offset plant); NC5 false-sensitivity rate on resistance/N CGT tests; NC6 lineage-TF genes (SOX10, MITF, PAX8, …
from P1 `biomarker_pairs.csv`) landing in class 3. Outside its bound = failed control.

## 8. Subtypes and confounder audit

8.1 **Per-gene subtype** S(g, d), d ∈ MD-strict: accurate-dependent lines of g for which d is a line-level determinant
(binary d: carriers; contribution β_gd(x_ℓd − x̄_d) of the full-D EN-pred fit is in the dependency direction and ≥ 25%
of the line's surrogate deviation). Lines may be multi-label; dependent but mispredicted lines are "unexplained" (§12).

8.2 **Biomarker module** M(d): genes with (g, d, same sign) ∈ MD-strict; member lines = carriers accurate-dependent on ≥
50% of them. Size reported as genes and as effective gene clusters (D dependency correlation |r| > 0.5, average
linkage).

8.3 **Within vs across lineage** (MD-strict ∪ controls ∪ S1-obs BH top 200). Testable lineage: ≥ 10 D lines, ≥ 3
carriers and ≥ 3 non-carriers; powered: power ≥ 0.5 for the pooled shrunk effect. Meta-analysis REML τ² + HKSJ CI; I²,
τ, 95% prediction interval.
- **cross-lineage**: meta BH q < 0.05 (F4) ∧ ≥ 2 powered lineages from ≥ 2 organ-system groups (§F `lineage_groups`) ∧
  same sign in ≥ 80% of powered lineages ∧ sign kept under leave-one-lineage-out for every contributing lineage;
  d × lineage interaction reported;
- **lineage-restricted**: per-lineage BH q < 0.05 (F5) ∧ interaction LRT p < 0.05;
- **lineage-unresolvable**: < 2 powered lineages; **undetermined**: otherwise.
Labels are recomputed by script from stored statistics; the same labelling on the permutations gives the false
cross-lineage rate. **Q1b**: Fisher OR (95% CI) of R replication and of CGT membership, cross-lineage vs
lineage-restricted MD-strict determinants (`results/subtypes/cross_lineage_test.csv`); a null result is a finding.

8.4 **Co-location.** For every CN/SIL determinant in a locus group, refit on lines discordant between d and each
neighbour; ≥ 5 discordant carriers and same sign at p < 0.05 ⇒ named at gene level, else at locus level
(`results/locus/resolution.csv`). MTAP-specific evidence comes from the CDKN2A-del/MTAP-intact discordant lines.

8.5 **CN-cutting artefact** (AMP determinants whose locus group contains g): (i) median effect of unexpressed amplicon
genes (log2TPM+1 < 1) in carriers minus non-carriers; (ii) within carriers, r(g effect, g expression). Label "cutting
artefact" if (i) < −0.2 and (ii) is not significantly negative; "amplified-driver addiction" if (i) ≥ −0.2 and g is
depleted; else "unresolved" (`results/confounders/amplicon_test.csv`). No literature exceptions; ERBB2, MYCN, CCND1 are
expected to pass (positive controls of the rule).

8.6 **Culture.** Cramér's V of subtype membership vs screened CellFormat, 3D, serum-free, medium family and library.
Module genes are tested for overlap with the pre-declared culture gene sets (§F; hypergeometric BH < 0.05); flagged
subtypes get a medium/serum-free stratified test. Where culture nests in lineage, tests are within stratum.

8.7 **Hypermutation and TP53.** Hypermutated = MSI ∪ top 5% DMG_BURDEN. LOF/HOT determinants must keep sign with
abs β ≥ 50% of full in non-hypermutated lines; SIG/CIN/WGD determinants must keep sign within TP53-mutant lines. Failures →
§12 class 6.

8.8 **Unsupervised check (Tier C).** Fold-centred M5 OOF (raw, lineage-residualised) → PCA → cosine kNN → networkx
Louvain (leidenalg/igraph not installed), 3 resolutions × 10 seeds; AMI with lineage, library, MD-strict subtypes and
`cv_fold` (≥ 0.02 ⇒ OOF artefact); clusters not explained by 𝒟 are "unnamed".

## 9. Positive controls and the gate

9.1 **File.** `configs/positive_controls.csv`, built by `scripts/build_positive_controls.py` from this section,
committed before any association (its commit must be an ancestor of the first commit writing `results/md/`). Schema:

| Column | Content |
|---|---|
| control_id | e.g. `PAR_SMARCA2_SMARCA4` |
| family | ADD-HOT, ADD-CN, ADD-EXP, PAR, LOF-SL, SIL, SIG, PROG, P53, W, LINEAGE-CHECK |
| dependency_symbol, dependency_entrez | dependency gene |
| feature_primary, feature_secondary | 𝒟 ids (`FAMILY:<entrez or name>`); secondary may be empty |
| expected_sign | S (β < 0 for feature presence) or W (β > 0) |
| scope | gate, class_check, untestable |
| citation | literature key (resolved in LITERATURE.md / docs/lit) |
| in_cgt | 1 if the pair is a CGT test (filled after `freeze-cgt`) |
| n_carriers_D, n_carriers_R | recomputed from raw caches (molecular data only) |
| gate_eligible | dependency ≥ 100 non-NaN D lines ∧ feature present with ≥ 5 carriers (binary) |
| note | e.g. "non-clinical", "gene-level proxy", "processing-limited" |

9.2 **Controls and encodings.**
- ADD-HOT (S): KRAS←HOT:KRAS; NRAS←HOT:NRAS; BRAF←HOT:BRAF (cell-line control; gene-level proxy for V600);
  PIK3CA←HOT:PIK3CA; CTNNB1←HOT:CTNNB1; EGFR←HOT:EGFR (secondary AMP:EGFR, "non-clinical").
- ADD-CN (S): ERBB2←AMP:ERBB2; MYCN←AMP:MYCN; CCND1←AMP:CCND1; PKMYT1←AMP:CCNE1; PSMC2←LOSS:PSMC2 (CYCLOPS;
  "processing-limited": Chronos arm-level correction biases it to null).
- ADD-EXP (S): ESR1←EXPDOSE:ESR1 (clinical expression biomarker; exempt from class 4).
- PAR (S): SMARCA2←LOFBI:SMARCA4; ARID1B←LOFBI:ARID1A; STAG1←LOFBI:STAG2; VPS4A←LOSS:VPS4B (sec. EXPPAR:VPS4B);
  VPS4B←LOSS:VPS4A; ENO2←LOSS:ENO1 (sec. EXPPAR:ENO1); PAPSS1←LOSS:PAPSS2 (sec. DEL:PAPSS2); PRMT5←DEL:MTAP
  (sec. SIL:MTAP); MAT2A←DEL:MTAP (sec. SIL:MTAP); ME3←LOSS:ME2; CREBBP←LOFBI:EP300; EP300←LOFBI:CREBBP;
  DDX3X←YABS:DDX3Y.
- LOF-SL (S): PIK3CB←LOFBI:PTEN (sec. DEL:PTEN); NFE2L2←LOFBI:KEAP1.
- SIL (S): WRN←SIL:MLH1; UBC←SIL:UBB.
- SIG (S): WRN←SIG:MSI; KIF18A←SIG:WGD (fallback SIG:CIN proxy, labelled).
- PROG (S): ADAR←PSCORE:ISG; NFE2L2←PSCORE:NRF2.
- P53 (W, coded on COMP:TP53MUT): MDM2, MDM4, PPM1D, USP7.
- W (W): CDK4←LOFBI:RB1; CDK6←LOFBI:RB1; CCND1←LOFBI:RB1.
- LINEAGE-CHECK (scope class_check): the 18 lineage-TF rows of P1 `biomarker_pairs.csv`, expected in class 3.
Controls with < 5 carriers are listed as `untestable` (scope) and in the report.

9.3 **Gate.** Carrier counts and `gate_eligible` are written to `results/controls/pc_eligibility.csv` before
`freeze-design`, which fixes the denominator. Recovered = (g, primary or secondary feature) in the list with the
expected sign (EN: selected; S1: BH q < 0.05; MD-strict: in catalogue). Pass (decided on MD-strict; all arms reported
with exact binomial CIs): ≥ 60% of gate-eligible controls overall and ≥ 50% within each family with ≥ 2 eligible rows.
Run once: `results/gates/G_pc_attempt1.json` (sha256 logged, never overwritten). The gate's verdict is never used to
choose between variants. After it, only bug fixes are allowed, each a row in `results/gate_iterations.csv` (row 0 =
first run; commit, reason, failing toy test that now passes, gate before/after), applied to all arms; every later
commit touching `src/` must also touch `tests/`. If more than one row exists, the CGT headline is the gate-disjoint set
and REPORT.md states "gate failed at attempt 1 (x/y); passed at attempt n after Z". The control list cannot change
after `freeze-design`.

## 10. Clinical ground truth (CGT)

10.1 **Build.** `scripts/build_cgt.py` reads CIViC nightly, MOAlmanac, the FDA CDx scrape (all frozen under
`external/`, sha256 in `configs/external_sha256.txt`) and the curated `configs/cgt_mapping.csv`, where every row carries
≥ 1 resolvable source record id (CIViC EID/AID, MOAlmanac id, FDA PMA/BLA/510(k) number in the frozen scrape, or NCT id
in a frozen CTG JSON; PMIDs supplementary, fetched without user identifiers). Unresolvable rows are dropped to
`results/clinical/cgt_dropped.csv`. A targeted CTG fetch of the NCT ids named in the mapping (≤ 30 records) is Tier A;
a broad harvest is Tier C. The script reads no effects or predictions (static import scan). Output
`configs/clinical_ground_truth.csv`, plus `configs/indication_map.csv` (indication → OncotreeCode/PrimaryDisease/CGT
feature, TCGA study ids) and `configs/cgt_features.csv` (§10.3), all committed and tagged `freeze-cgt` before any
dependency data are joined (`git merge-base --is-ancestor`). The docs/lit draft table is a starting list, not a source.

10.2 **Scoring unit, sign, partitions.** A CGT *test* = (𝒟/CGT feature, target set T, comparator, indication set);
drugs mapping to one test are collapsed (drug list kept). Comparator: `within-indication` (default for
lineage-restricted biomarkers) or `pan-cancer + lineage covariate` (tumour-agnostic). Sign: S ⇔ β < 0 on the target's
Chronos effect for feature presence; W ⇔ β > 0; absence biomarkers coded as presence of the mutation with the sign
flipped. Partitions: `sensitivity` (S), `bypass-resistance` (W; scored, e.g. RAS-mut → EGFR in CRC, RB1 LOF → CDK4/6),
`acquired-resistance` (T315I, C797S, KIT D816V for imatinib; listed, untestable), `expected-null` (N),
`clinically-refuted` (§10.6). An S hit on a W test is a sign error. Columns: `test_id, feature, target_set, comparator,
indication_codes, sign, partition, ko_class, tier, drugs, source_ids, in_gate, tumour_agnostic, excluded_indications,
allele_level_required`; (feature, target_set, comparator) unique.

10.3 **Corrections and operationalisation (before `freeze-cgt`).**
- KIT D816V: avapritinib/midostaurin S; imatinib → acquired-resistance. del(17p)/TP53 → BCL2: `population-only`, N.
  BRAF sensitivity tests require V600. EGFR amplification is not a CGT test (cell-line control only); EGFR CGT feature =
  activating alleles (L858R, exon-19 del, G719X, L861Q, exon-20 ins, T790M). PTEN → AKT: T = {AKT1, AKT2, AKT3}. KRAS/NRAS
  WT/mutant → EGFR merge into one W test (comparator RAS/BRAF-WT CRC). ESR1 LBD: N vs ER+ WT, untestable. CDKN2A±MTAP
  merge into the 9p21/MTAP test with the MTAP-specific resolution test (§8.4).
- `scripts/build_cgt_features.py` (molecular data only): alleles from OmicsSomaticMutations (BRAF V600; KRAS G12C,
  G12/G13/Q61; EGFR activating set; PIK3CA E542K/E545X/H1047X/C420R; ERBB2 TKD; AKT1 E17K; IDH1 R132; KIT exon 11,
  D816; FGFR3 R248/S249/G370/Y373) — absent file ⇒ `untestable: allele`, never mapped to gene-level HOT; fusions with
  canonical partners (CIViC), orientation and FFPM rule of §3, known positives must be carriers (NCI-H3122, NCI-H2228
  ALK; K562 BCR::ABL1; MV4-11 KMT2A; MDA-MB-175-VII NRG1), ALK carriers ≤ 10; BRCA1/2 biallelic (dosage 2, or damaging +
  LoH, or SIL BRCA1) excluding hypermutated lines; ER+ = `ModelSubtypeFeatures`, else a 2-component Gaussian mixture on
  ESR1 within Breast ∪ Uterus ∪ Ovary (posterior 0.5); HER2+ = AMP:ERBB2 ∨ "HER2+"; MTAP-loss = DEL ∨ SIL; MSI per §1.5.
- KO-concordance classes: `concordant` (incl. allele-specific inhibitors of addicted oncogenes vs WT, e.g. KRAS G12C);
  `concordant-family` (|T| ≥ 2, pre-registered members, e.g. CDK4/6, MEK1/2, pan-AKT, FGFR1–4, PARP1/2, NTRK1–3; hit if
  any member, ceiling Bonferroni × |T|); `concordant-KO/combination-required` (BRAF V600E CRC); `discordant-mechanism`
  (PARP trapping, IDH neomorphs, degraders); `population-only`; `non-cell-autonomous` (ICI, ADCs incl. HER2-low T-DXd);
  `untestable`. Combinations are scored on the anchor drug's target.
- Tiers (mutually exclusive, highest wins): T1 regulatory approval with label biomarker; T2 CIViC A/B Supports or
  NCCN-cited, not approved; T3 trial-stage biomarker-selected NCT only; refuted. CDx presence is an attribute.
- Forward-looking trial-stage rows only if the cached CTG record matches the biomarker regex: CCNE1 AMP → PKMYT1,
  ATM LOF → ATR, SMARCA4 LOF → SMARCA2 (`in_gate`), BRCA1/2 biallelic → POLQ and → USP1.

10.4 **Q4-0 data ceiling** = S1-obs on that test: for each S/W test with ≥ 3 D carriers in the comparator, the §5.3
model restricted to the comparator (Wilcoxon as robustness), one-sided in the expected direction, BH q < 0.1 over CGT
tests (F7). Labels `data-supported`, `data-contradicted`, `data-null`, `untestable` (`results/clinical/cgt_ceiling.csv`;
labels recomputed from numbers). Recall is reported over testable tests (≥ 5 carriers) and over data-supported tests;
the gap is the method's miss rate; concordant data-null/contradicted tests feed §12 class 8.

10.5 **Q4 metrics** (`results/clinical/q4_summary.json`).
- **Q4a / P2**: recall of testable S/W tests (concordant classes) in each arm's top-200 (secondary 100/500; hit = g ∈ T
  and d = feature, any sign), exact McNemar between arms, on `all` and `gate_disjoint` sets (`recall_all`,
  `recall_gate_disjoint`); per tier with exact binomial CIs; Cochran–Armitage trend only if every tier has ≥ 5 testable
  tests, else "descriptive only".
- **Direction accuracy** = hits with the pre-registered sign / all hits; sign-error rate = 1 − that; for every arm and
  the ceiling. A toy planted wrong-sign pair must count as a sign error, not a miss.
- **N-partition hit ratio** vs 1,000 random (g, d) matched for carrier count and lineage spread (flag > 2);
  false-sensitivity rate on resistance/N tests.
- **Q4b actionability**: pair-level (g ∈ T and d = a sensitivity/bypass test feature) and feature-level fractions of
  MD-strict, vs S1-obs and vs 200 within-list permutations of g labels; PROG/PSCORE determinants excluded from the
  numerators; line coverage by lineage.
- **Q4c population** (`results/clinical/q4c_population.csv`): within-indication β_in; per non-indicated lineage (≥ 3/≥
  3) β_L; classes on-label / beyond-label / off-label-only / indication-untestable; population precision; carriers-only
  OR of accurate-dependent membership, indicated vs non-indicated; tumour-agnostic tests must be cross-lineage (§8.3);
  excluded indications reported, not scored.
- **Q4d**: recovery of discordant or non-cell-autonomous tests is a flag, not a success.
- **Prevalence (Tier B)**: D vs TCGA PanCancer Atlas (cBioPortal, ≈ 40 cached checksummed calls; NSCLC = LUAD + LUSC)
  per (biomarker, indication); Wilson CIs, log2-ratio flag > 1; acquired biomarkers `tcga_na`.

10.6 **Clinically refuted set** (expected KO sign S): KRAS-mutant NSCLC → MEK1/2 (SELECT-1); CCND1 amp/CDKN2A loss →
CDK4/6 as a selection biomarker (PALOMA-1 cohort 2); BRAF V600E CRC → BRAF monotherapy; EGFR-amplified GBM → EGFR TKI.
Each needs a resolvable source id before `freeze-cgt` or is dropped; ≥ 3 rows required. `recall_refuted` vs
`recall_approved` (exact CIs); similar rates are a Q5 finding.

## 11. Combinatorial screens (Tier B, 1 h)

For MD-strict and S1-obs candidates of the form "loss of A → dependency on paralog B", look up A–B in the downloaded
Dede 2020, Thompson 2021 and Esmaeili Anvar 2024 tables and the Klingbeil/Flister bioRxiv supplement (analysis only, no
redistribution). Parrish 2021 has counts only (processed GI tables blocked) and is not reprocessed; predicted-SL
resources are excluded. Labels from file rows only: replicated (hit in ≥ 2 screens or ≥ 2 lines), tested-not-replicated,
untested. Validation rate of MD-strict vs S1-obs paralog candidates vs matched screened background (Fisher, CI);
supporting evidence only (`results/combinatorial/lookup.csv`).

## 12. Limits: failure classes

Per gene, first applicable class (`results/limits/gene_classes.csv`): 1 no between-line contrast (D SD < 0.15; not the
common-essential flag); 2 not predictable (low tier); 3 accurate only via lineage/library (`r_within` < 0.10) or
unnameable in 𝒟 (fidelity < 0.10, or M1 named share < 0.5); 4 determinant is only expression gating, or own EXPDOSE
without AMP/SIL backing, unless the gene is in the clinical expression-biomarker list (ESR1, ERBB2, AR, PGR, FOLR1); 5
determinant too rare or R-underpowered; 6 confounded (culture, library, cutting artefact, unresolved co-location,
hypermutation, TP53, subtype-confounded); 7 powered but not replicated; 8 KO–drug discordance (CGT class or ceiling
data-null/contradicted); 9 valid candidate. A gene with `dependent_rule = relative` is never class 1 while SD ≥ 0.15.

Per (g, d) (`results/limits/pair_classes.csv`): **predictor-blind** (S1-obs hit, replicated in R or XLIB, not
EN-pred-selected); **predictor-only, not replicated** (EN-pred hit, powered, failed); **processing-limited** (own-arm CN
and CYCLOPS under Chronos correction); **pan-lethal by CRISPR** (Krill-Burger 2023 caveat). Unassayed determinant types
(§3) are listed. Each class gets counts, examples and evidence columns.

## 13. Verification and anti-drift

13.1 **Toy simulation** (`tests/test_toy_pipeline.py`; generator `src/biomark/toy.py`), **written by a different agent
than the pipeline author** (checked by git author or `Agent:` trailer). It calls the same entry points with the same
`design_params.json` (hash asserted); planted effects and `toy_seed` fixed at `freeze-design`. 600 lines × 12
lineages, 300 genes, 400 features with correlated blocks (r ≈ 0.9), a library factor nested in 3 lineages; the toy
predictor is a ridge model trained by 5-fold CV on the noisy toy observations, plus per-fold offsets. Plants and pass
bounds (§F `toy_s13`): pan-lineage recall ≥ 0.8; own-locus CN (own-CN covariate dropped by rule) recall ≥ 0.8; MSI-like signature
(with enriched passenger LOF) recovered in ≥ 4/5 and passenger LOF selected in ≤ 1/5; co-located passengers at locus
level in ≥ 4/5; culture-confounded in MD-strict ≤ 1/5; lineage-only genes 0 in MD-strict and ≥ 80% class 3; fold-offset
feature selected ≤ 1/5 after centring and ≥ 3/5 without (positive control for the artefact); genomic feature with an
expression proxy recovered in stage G ≥ 4/5; full null FDR_emp ≤ 0.10 per list (B = 20); a wrong-sign CGT pair counts as
a sign error. The D-stage nominal anti-conservativeness under the null is reported as a diagnostic.

13.2 **Unit tests**: dictionary and locus-group counts vs raw caches; self-null removal; covariate drops (ERBB2 ×
AMP:ERBB2 drops own CN, WRN × SIG:MSI drops MSI); FWL vs lineage scaled by 1e4 (|Δβ| < 1e-3); subsample manifest
(disjoint pairs, no split patient); determinism (1 vs 12 workers); additivity; NMF reproducibility; CGT schema; control
features exist; EN-pred blind; guard refusals (deny-list, moved tag, changed hash).

13.3 **Independent reimplementation.** Before `freeze-catalogue`, a separate agent reimplements, from raw inputs in
`tests/independent/`, five quantities: EN-pred selection for 20 genes, S1-obs p for those genes, gate recall, BH counts
and confounding deltas for 10 pairs; tolerance 1e-6 (exact for counts) → `results/gates/G_indep.json` with agent name.

13.4 **Frozen files and deviations.** At `freeze-design`, `configs/frozen_hashes.json` records sha256 of §F configs,
the dictionary build config, controls, culture/lineage lists, `null_genes.txt` and `cutting`/locus rules. A later
change needs a decision `D<n>` in RESEARCH_TREE.md naming the file and whether dependent outputs were seen (if so the
report gives both values); triggered fallbacks write `results/deviations.csv` automatically. No decision literals in
`src/biomark/*.py` outside the params loader, and no gene symbols, Entrez ids or `ACH-` ids as literals in `src/` or
`scripts/` (leakscan; `configs/` and `tests/` whitelisted).

13.5 **Report numbers and claims.** `scripts/check_report_numbers.py` is ported from Part 1 (copied, not imported; tag
`<!--n:results/...json#path-->`), run `--strict` at the final gate; every tagged source must be git-tracked and resolve
to a `results/registry.csv` row with `status = ok`, clean tree and the frozen design_params hash (dirty-tree runs are
`exploratory` and uncitable). `results/claims.csv` (claim id, class, evidence file, rule) is generated from §8.3, §6
and §12; a REPORT.md sentence applying a controlled word (validated, replicated, cross-lineage, valid biomarker,
recovers, significant) to a named gene or pair needs `<!--c:claim_id-->` of a matching class.
- **.gitignore** ignores `*.gz`/`*.parquet` everywhere (checked: `results/x.csv.gz` is ignored), so such results
  would silently stay uncommitted. Rule: committed results are plain `.csv`/`.json` ≤ 5 MB; bulk per-test tables live in
  `cache/` with sha256 in `results/bulk_manifest.csv`; `gate.py tracked` asserts every registry output is tracked.

13.6 **Ledger, registry, gates.** `configs/design_ledger.csv` (one row per deliverable and per accepted review item,
status ∈ planned/done/declined/not done) is generated from reviews/round1_response.md at `freeze-design`; no row may be
`planned` at the final gate. `results/registry.csv` (file-locked append: run_id, script, commit, dirty, config hash,
design_params hash, data_status hash, status, outputs). Parallel tracks use their own output folders and D-id ranges
(M4 D100–, M5 D200–, M6 D300–). `pip freeze` → `configs/pip_freeze.txt`. All checks are subcommands of
`scripts/gate.py` (modelled on Part 1's): `leakscan, frozen, tracked, controls, G_pc, cgt, data, ledger, G_rep, G_indep,
part1`.

## 14. Tiers, compute, schedule (≈ 19 h left; `configs/tiers.json`)

- **Tier A (must):** 𝒟; four arms; MD-strict; confounder audit; gate; null; subtypes per gene and per module;
  cross-lineage; co-location/cutting/culture/hypermutation; component readout; CGT recovery, direction, actionability,
  population by lineage; replication in R and XLIB; limits. **Tier B:** attribution panel; combinatorial lookup (1 h);
  TCGA prevalence; LMM top-200; EN-pred-clean. **Tier C:** clustering; CTG harvest; cross-fit sensitivity (select on
  fold-k predictions, test on fold-k effects, Stouffer over pairs selected in ≥ 3 folds).
- **Compute:** loky, 12 workers × 50-gene chunks, 1 BLAS thread each, read-only memmaps, per-gene seed
  `seed_list[i] ^ int(entrez)`, peak RSS ≤ 12 GB; refit never concurrent with 12-worker jobs. Each job writes
  `results/timing/<job>.json`; `scripts/budget_check.py` compares with the timeboxes.
- **Schedule (h):** M2 finish + round 2 + freeze 1.5; M3b infrastructure + toy 3.0; M4 discovery 3.0; M5 subtypes 2.5;
  M6 CGT + clinical 2.0; M7 independent check, freeze, R, XLIB 1.0; M8 limits + report 3.5; total 16.5 + 2.5 reserve.
  The cut order and never-cut list are in `tiers.json`.

## §F Frozen parameters

Each row maps 1:1 to a key in the named file; code reads only the file.

| § | Parameter | Value | File : key |
|---|---|---|---|
| 2 | eligibility | (variable ∨ selective) ∧ m5 tier ∈ {high, medium} ∧ r_within ≥ 0.10 (strata ≥ 5 lines; lineage ≥ 10 else `_rare`) | design_params : eligibility_s2 |
| 2 | dependent | < −0.5 (sens. −1.0); relative < D 25th pct if frac_below > 0.5; no contrast SD < 0.15 | design_params : dependency_s2 |
| 2 | accurate-dependent | rank ≤ n_g (D); R quantile n_g/n_D | design_params : dependency_s2 |
| 3 | carriers; CN cuts | ≥ 5; AMP > 2, DEL < 0.25, LOSS < 0.75 | design_params : dictionary_s3 |
| 3 | SIL | < 1 log2TPM; ≥ 50% expressed in lineage; < median − 2; V > 0.5 flagged | design_params : dictionary_s3 |
| 3 | locus groups; resolution | 10 Mb, Jaccard ≥ 0.5; ≥ 5 discordant, p < 0.05 | design_params : dictionary_s3.locus_* |
| 3 | cap | shared 𝒟 ≤ 3,000; escalate min carriers 10, 20 | design_params : dictionary_s3 |
| 3 | NMF; programmes | 50, 5,000 genes, nndsvda/cd/500, nmf_seed; ISG, NRF2 | design_params : dictionary_s3.nmf / programme_scores |
| 3 | MSI order; hypermutated | GlobalSignatures → CMP → ModelSubtypeFeatures; MSI ∪ top 5% DMG_BURDEN | design_params : dictionary_s3 |
| 3 | spectra; methylation | reference-free SPEC_* + sanity; β ≥ 0.7 | design_params : dictionary_s3 |
| 4 | EN | l1 0.9; FWL lineage+library(+NA flags); stages G→X; SIS K 300 (fallback 150, pilot 50 genes, 45 min); 25 λ, eps 1e-2, tol 1e-3; P1 cv_fold; 1-SE fixed; 25 pairs; freq ≥ 0.6; sign ≥ 0.9 | design_params : surrogate_s4 |
| 4 | fidelity | r2_oof_beyond_lineage_library; class 3 < 0.10 | design_params : surrogate_s4 |
| 4 | ranking | EN freq, abs β, id; S1 p, abs β, id; locus dedup | design_params : surrogate_s4; primary_endpoints : arms |
| 4 | attribution | ≤ 300 genes; r ≥ 0.999; additivity 1e-5; named share 0.5; top-5; line share 0.25 | design_params : attribution_s4 |
| 5 | covariates | blocks; drop abs ρ > 0.7 etc.; VIF ≤ 10; CR1/CR2 patient; delta rule; disease ≥ 5 | design_params : confounders_s5 |
| 5 | S1; MD-strict | partial t, ≥ 5 carriers, var ≥ 1e-8, BH 0.05; same sign ∧ q_sc < 0.05 | design_params : s1_s5, md_strict_s5 |
| 5 | test families | F1–F7 | design_params : test_families_s5 |
| 6 | P1 | k 100; power ≥ 0.5; R ≥ 3/3; one-sided 0.05; bootstrap 1,000 patient; contrasts and branches | primary_endpoints : P1_primary |
| 6 | power/shrinkage | EB normal scale mixture on F1, shared across arms | design_params : replication_s6 |
| 6 | secondary; XLIB | k 50/200/500, variants; Avana-only → no-Avana | primary_endpoints : P1_primary.secondary_*, XLIB_co_endpoint |
| 6 | LMM | top 200, p 0.001 | design_params : lmm_s6 |
| 7 | null | joint row permutation in lineage×library; 300 genes; B 20 (→ 50 by timing rule, cap 180 min); FDR_emp ≤ 0.10; S1 200 | design_params : null_s7 |
| 7 | negative controls | NC1–NC6 bounds (simulated before freeze) | design_params : negative_controls_s7 |
| 8 | subtypes; cross-lineage | module ≥ 50%; clusters abs r > 0.5 avg; testable ≥ 3/3; powered 0.5; q 0.05; ≥ 2 lineages, ≥ 2 organ groups; ≥ 80% sign; LOLO | design_params : subtypes_s8 |
| 8 | lineage groups; culture sets | 8 organ groups; 4 culture gene sets, BH 0.05 | design_params : subtypes_s8 |
| 8 | cutting; hypermutation; clustering | −0.2; ≥ 50% β; PCA 50, k 15, Louvain, AMI fold < 0.02 | design_params : subtypes_s8 |
| 9 | gate | ≥ 60% overall, ≥ 50% per family (≥ 2 rows); MD-strict decides; exact CI | design_params : gate_s9 |
| 9 | controls | schema §9.1, list §9.2 | positive_controls.csv (built by script) |
| 10 | CGT | ceiling ≥ 3 carriers, q < 0.1; testable ≥ 5; trend ≥ 5/tier; N ratio 2 (1,000 pairs); actionability null 200 | design_params : clinical_s10 |
| 10 | P2 | top-200 (100/500); McNemar; all + gate-disjoint | primary_endpoints : P2_primary_clinical |
| 11 | combinatorial | 4 sources; ≥ 2 screens/lines; matched background; Fisher; 60 min | design_params : combinatorial_s11 |
| 13 | toy | dimensions, plants, pass bounds | design_params : toy_s13 |
| 14 | compute; tiers | 12 workers, chunk 50, 12 GB; timeboxes, budgets, cut order | design_params : compute_s14; tiers.json |
| all | seeds | seeds, subsample, permutation, null_gene_subset, bootstrap, toy, nmf | seeds.json |
