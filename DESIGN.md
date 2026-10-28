# DESIGN — Predictor-derived dependency subtypes and their molecular biomarkers

Version: v1 (draft for review round 1). Status: NOT frozen.
Inputs: PROMPT.md (goal, fixed), docs/part1_brief.md (Part 1 interfaces), docs/lit/*.md and LITERATURE.md (prior work),
docs/data_availability.md (data access).

## 0. Question and claims we want to be able to make

The goal is fixed by PROMPT.md. Operationally, we ask:

- **Q1 (subtypes).** For gene dependencies that the Part 1 predictor predicts accurately, can the predictor's own
  reasons be read out as molecular determinants, and do they partition dependent cell lines into subtypes that share the
  dependency *for the same molecular reason*, within and across lineages?
- **Q2 (benchmark).** Are predictor-derived biomarkers better than standard single-feature association tests, judged by
  out-of-sample replication, recovery of positive controls and clinical ground truth, robustness to confounders, and
  testing burden?
- **Q3 (confounders).** Which candidate biomarkers survive control for lineage, culture conditions, screen library and
  screen length, passenger alterations and hypermutation, genomic co-location (including Cas9 cutting toxicity), and
  multiple testing?
- **Q4 (clinical).** Do the subtypes recover established predictive biomarkers with drugs that reached trials, with the
  correct direction of effect, the clinically actionable subtypes and the correct patient population?
- **Q5 (limits).** Which dependencies cannot be reliably used as cancer biomarkers, and why? Each eligible dependency is
  assigned a failure class.

Non-goals: building a new predictor (Part 1's is used as frozen); patient outcome modelling; drug-response prediction.

## 1. Data and roles

| Role | Lines | Use |
|---|---|---|
| Discovery (D) | 1022 Part 1 dev lines | Everything in §§2–6: attribution, surrogates, subtypes, single-feature tests, confounder audit |
| Replication (R) | 255 Part 1 test lines | Opened **once**, after the discovery catalogue and the benchmark lists are frozen (tag `freeze-catalogue`) |
| Prediction-only | 838 expressed, unscreened lines | Not used in primary analyses (optional: subtype prevalence) |

- **R is not pristine.** Part 1 opened R's CRISPR data once, to score its predictor; no biomarker analysis was ever run on
  R. Part 1's published per-gene test metrics (`P1/results/test/`) are **not read** by any Part 2 script before
  `freeze-catalogue` (enforced: `src/biomark/guard.py` refuses R targets and P1 test metrics before the tag exists and
  logs each access to `results/replication_access_log.csv`).
- **Predictions:** D uses Part 1's out-of-fold (OOF) headline predictions (`m5_avg`, run `20261027-161041_m5_avg_6dbebcf4`,
  1022 × 18434); R uses Part 1's all-dev refit predictions (`P1/cache/test/preds/m5_avg.parquet`). OOF predictions are
  centred per gene within fold before any cross-line analysis (Part 1 D3 pooled-OOF artefact).
- **Targets:** Chronos gene effect, DepMap 26Q3 (`P1/cache/crispr.parquet`, read directly in Part 2 under the Part 2
  guard). "Dependent" = effect < −0.5 (pre-set; sensitivity −1.0).
- **Molecular data:** Part 1 caches (expression, CN, damaging, hotspot, fusion, Model, ModelCondition, screen maps);
  pending user download (DATA.md): OmicsGlobalSignatures (MSI, ploidy, WGD, CIN, LoH), CRISPRConfounders (doubling time,
  screen QC), OmicsSomaticMutations (full calls; mutation-level filtering), OmicsCNSegmentsWGS (segments, arm-level
  events), legacy CCLE RRBS methylation (optional), PRISM (optional). Any analysis that needs a file that does not
  arrive by the deadline (§10) is recorded with status "not done: data unavailable".
- **Part 1 is read-only.** Part 1 code is imported only with bytecode writing off and every write-path redirected into
  this folder (part1_brief §6); Part 1 scripts are never run. A test (`tests/test_part1_untouched.py`) checks that
  `git status` in ../project is clean and that no file under ../project is newer than the start of this project.

## 2. Eligible dependencies

- **Eligible genes (primary):** (variable ∨ selective) ∧ Part 1 m5 tier ∈ {high, medium}, using Part 1's frozen dev-CV
  tables (`P1/configs/gene_classes.csv`, `P1/configs/predictable_genes.csv`). This operationalises "accurate prediction".
- **Positive-control genes** (§7) are analysed regardless of tier and reported separately.
- **Accurate line-level predictions:** for subtype assignment, a (line, gene) prediction is "accurate-dependent" if the
  OOF prediction is below the gene's dev 20th percentile of predictions **and** the observed effect < −0.5. Lines that
  are dependent but predicted as not dependent are "unexplained dependencies" and feed §8.

## 3. Interpretable feature dictionary 𝒟

The determinants must be assayable molecular features. 𝒟 is built from dev lines only, fixed before any association is
computed, and versioned (`configs/dictionary_v1.json`, built by `scripts/build_dictionary.py`).

| Family | Definition | Minimum frequency |
|---|---|---|
| HOT | hotspot mutation per gene | ≥ 5 dev lines |
| LOF | damaging mutation (dosage > 0) per gene, restricted to OncoKB cancer genes ∪ genes with a paralog/partner prior | ≥ 5 |
| AMP / DEL | relative CN > 2 / < 0.25 (Part 1 thresholds), then **locus-collapsed** (§5.4) | ≥ 5 |
| FUS | fusion partner gene (Part 1 F_gene) | ≥ 5 |
| SIL | gene not expressed (log2 TPM+1 < 1) where the gene is expressed in ≥ 50% of dev lines; candidate methylation/deletion silencing; split into SIL-DEL (deep deletion) vs SIL-other (methylation if RRBS arrives) | ≥ 5 |
| EXP-own, EXP-paralog | target gene's own expression; expression of its top-3 Ensembl paralogs (Part 1 priors) | gene-specific |
| PROG | expression programmes: 50 NMF factors of dev expression (non-negative, 5,000 most variable genes, seed list), each annotated by Hallmark/lineage enrichment | — |
| SIG | MSI score, TMB (log1p damaging count), ploidy, WGD, CIN, LoH fraction (OmicsGlobalSignatures; Part 1 C_instab/M_burden proxies if not downloaded) | — |
| COMP | TP53 status, RB pathway status (Part 1 composite flags) | — |

Lineage (OncotreeLineage) and culture/screen covariates (§5) are **not** determinants; they enter as covariates.

## 4. Reading determinants out of the predictor

Two complementary readouts. Both are run; agreement between them is a reported quality measure.

### 4.1 Primary: surrogate (distillation) of the predictor in 𝒟
For each eligible gene g, fit a sparse linear surrogate of the fold-centred OOF headline prediction ŷ_g on 𝒟 plus
unpenalised lineage indicators:

ŷ_g ~ lineage (unpenalised) + Σ_d β_gd · d   (elastic net, α = 0.9; λ by patient-grouped 5-fold CV, 1-SE rule)

- **Why predictions, not observations:** the predictor is a fixed function of molecular inputs, so its prediction is a
  denoised readout of what it learnt; the surrogate expresses that function in assayable terms. Its **fidelity**
  (R² of surrogate vs ŷ_g, beyond lineage) measures how much of the predictor's reasoning 𝒟 can name.
- **Selection stability:** 50 patient-grouped half-sample subsamples (seed list); a determinant is *selected* if it is
  chosen in ≥ 60% of subsamples with a consistent sign (stability selection). Compute guard: if the 50× scheme exceeds
  4 h on 12 cores for the eligible set, fall back to 20 subsamples (pre-declared).
- **Line-level determinant:** for an accurate-dependent line, its determinant(s) = selected d with β_gd · d_line
  contributing ≥ 25% of the line's surrogate deviation from the lineage mean in the dependency direction.

### 4.2 Secondary: direct attribution of the Part 1 model
- Refit M1, M2, M3 on Part 1's dev folds in this folder (caches redirected); verify the refit reproduces Part 1's OOF
  predictions (max |Δ| ≤ 1e-4 per component, else stop and investigate).
- M1 (linear kernel ridge on expression): exact per-feature contributions (z − z̄)·w; aggregated to 𝒟 terms EXP-own,
  EXP-paralog, PROG (projection of w_g onto the NMF basis; residual reported).
- M2 (13-feature ridge): exact contributions for own/paralog expression, CN, deletion, mutation.
- M3 (LightGBM): TreeSHAP (`pred_contrib=True`) over its 51 features, mapped to 𝒟 families by feature identity.
- M5 = mean of the three (additive).
- **Agreement metric:** for each selected (g, d) from 4.1, is d (or its 𝒟 family for the same gene) among the top-5
  direct-attribution terms for the accurate-dependent lines? Reported as a fraction; disagreements are listed.

## 5. Candidate biomarkers, subtypes, and confounder control

### 5.1 Candidate biomarker (model-derived list, MD)
A candidate is (g, d, sign). It enters MD if (i) selected by 4.1, and (ii) d is associated with the **observed** effect
of g in D, in the full confounder model (§5.3), at BH-FDR < 0.05 over all MD candidates.

### 5.2 Subtypes
- **Per-gene subtype** S(g, d): accurate-dependent lines of g whose line-level determinant includes d.
- **Cross-gene subtype (biomarker module)** M(d): the set of genes {g : (g, d, same sign) ∈ MD} and the lines carrying d
  that are accurate-dependent on ≥ 50% of them. Lines may carry several subtypes (multi-label).
- **Unsupervised check:** Leiden clustering (kNN graph, k = 15) of dev lines on their direct-attribution profiles over
  eligible genes (fold-centred); report adjusted mutual information with (a) lineage, (b) library group, (c) MD subtype
  labels; clusters not explained by 𝒟 are listed as "unnamed".
- **Within vs across lineage:** for each (g, d), within-lineage effects in every lineage with ≥ 3 carriers and ≥ 3
  non-carriers, combined by random-effects meta-analysis; a determinant is **cross-lineage** if the meta effect is
  significant and has the same sign with nominal p < 0.05 in ≥ 2 lineages; otherwise **lineage-restricted** (named) or
  **pooled-only**. Leave-one-lineage-out: the effect must keep its sign when the most-represented lineage is dropped.

### 5.3 Confounder model
Observed effect_g ~ d + lineage + library group + screen length (days) + growth pattern + serum-free medium + TMB/MSI +
own-gene CN of g (+ doubling time and screen QC if CRISPRConfounders arrives).
Report the coefficient of d with and without each covariate block ("confounding delta"); flag if |β| drops > 50%.

### 5.4 Specific confounders
- **Culture conditions:** cross-tabulate subtype membership vs GrowthPattern, organoid ModelType, medium formulation
  and serum-free status (Cramér's V); where culture nests within lineage (haematopoietic ≈ suspension), the test is run
  within stratum only. Subtypes enriched for integrin/actin or lipid/cholesterol programmes are flagged as possible culture
  artefacts.
- **Library / screen length:** library group as covariate; plus a library-stratified replication within D (discover
  in Avana-only lines, check sign in the rest).
- **Passenger alterations and hypermutation:** LOF restricted to cancer/prior genes; TMB/MSI covariate; sensitivity
  analysis excluding hypermutated lines (top 5% TMB or MSI-high); for LOF determinants, the effect must persist in MSS
  lines.
- **Genomic co-location:** AMP/DEL features whose CN calls correlate r > 0.8 across dev lines within 10 Mb are merged
  into a locus feature (e.g., 9p21: CDKN2A/CDKN2B/MTAP). The mechanistic gene within a locus is named only if it is the
  dependency's paralog, pathway partner or the dependency itself; otherwise the determinant is reported at locus level.
  Conditional tests (d | neighbour) for every CN or SIL determinant.
- **Cas9 cutting / CN artefact:** for AMP determinants, test whether the dependency gene lies in the amplified locus;
  if yes, label "possible cutting artefact" unless the dependency is the amplified oncogene itself with literature
  support (e.g., ERBB2).
- **Multiple testing:** BH within each list (MD and S1 separately, each at its full test count); a lineage-stratified
  permutation null (d permuted within lineage, 200 permutations) calibrates the end-to-end pipeline (number of MD
  candidates expected under the null); report empirical FDR.

## 6. Benchmark against single-feature association tests (S1)

- **S1:** for every eligible gene g and every d ∈ 𝒟 (with gene-specific EXP terms), the same confounder model (§5.3)
  fitted on **observed** effects; likelihood-ratio test; BH-FDR < 0.05 over all S1 tests. GDSC-ANOVA-style with
  lineage covariate. A Gonçalves-style linear mixed model with a cell-line similarity random effect is run as a
  sensitivity for the top 1,000 S1 hits.
- **S1-pred (ablation):** S1 on fold-centred predictions instead of observations, to separate denoising from sparse
  multivariate selection.
- **Comparison (same eligible genes; same 𝒟):**
  1. Replication in R (once, frozen lists): fraction of hits with the same sign and one-sided p < 0.05 in R, and the
     effect-size correlation D vs R; also at matched list sizes (top-k by p-value, k = 50, 100, 200, 500).
  2. Recall of positive controls (§7) and of clinical ground truth (§7).
  3. Fraction cross-lineage, fraction flagged by the confounder audit.
  4. Number of tests (burden).
- **Power:** before opening R, compute per candidate the power to detect the D effect size in R at one-sided α 0.05;
  replication rates are reported over candidates with power ≥ 0.5 and over all.

## 7. Ground truth and controls

### 7.1 Positive controls (pipeline gate, pre-declared)
WRN–MSI; PRMT5 and MAT2A–MTAP/9p21 deletion; SMARCA2–SMARCA4 LOF; ARID1B–ARID1A LOF; STAG1–STAG2 LOF; VPS4A–VPS4B
loss; ENO2–ENO1 loss; PAPSS1–PAPSS2 loss; KRAS/NRAS/BRAF/PIK3CA/CTNNB1 hotspot → own dependency; ERBB2 AMP → ERBB2;
EGFR mut/amp → EGFR; MYCN AMP → MYCN; ESR1 expression → ESR1; MDM2 dependency in TP53-WT. **Gate:** MD must recover
≥ 60% of those with ≥ 5 carriers in D (recovered = (g, d) in MD with the expected sign). If not, stop and diagnose before
any biological claim.

### 7.2 Negative controls
- Lineage-stratified permutation of 𝒟 (§5.4) → expected MD size near 0.
- Random genes of matched variance from the low-tier set → surrogate fidelity near 0, no MD candidates.
- Known artefact classes must be labelled as such: CN-amplicon cutting, ancestry-linked guide SNPs (if detectable),
  organoid/library proxies.

### 7.3 Clinical ground truth (CGT)
- Sources (scriptable, licences in DATA.md; docs/lit/clinical_ground_truth.md): CIViC nightly evidence + assertions
  (CC0; Predictive, levels A/B, "Supports"), MOAlmanac (CC BY 4.0), FDA CDx list (HTML scrape, cross-validated against
  CIViC `fda_companion_test` and MOAlmanac FDA indications), OncoKB free `cancerGeneList`/`allCuratedGenes` (gene-level
  only), ClinicalTrials.gov v2 for ≥ 1 biomarker-selected NCT per pair. Cancer Genome Interpreter (no scriptable
  download; last update 2022) and JAX-CKB (licence) are excluded. Every download is frozen with its version/date in
  `configs/external_sha256.txt`.
- **Pre-registration:** the expected sign (S = stronger dependency in biomarker+, W = weaker, N = KO not informative)
  of every pair is recorded before any dependency data are joined. Recall is scored on S/W rows; N rows are an
  expected-null set. Resistance pairs are negative controls (must not be assigned a sensitivity direction).
- **Tiers:** recall is reported per evidence tier (CIViC assertion Tier I-A; FDA-CDx-backed; CIViC A/B; trial-stage
  only). Recall is expected to rise with tier; if it does not, that is reported as evidence the subtypes track something
  other than actionability.
- **Forward-looking set:** MTAP-del → PRMT5/MAT2A and MSI → WRN were discovered in this data type and reached phase I;
  recovering them with tumour-agnostic breadth is the strongest available prospective check.
- **Mapping:** (biomarker alteration → 𝒟 feature, drug → target gene(s) = dependency, expected sign, indication(s),
  tumour-agnostic flag, KO-concordance class). KO-concordance classes: *concordant* (target inhibition ≈ loss:
  kinase addiction), *discordant-mechanism* (PARP trapping, neomorphic IDH1/2, degraders of mutant protein,
  allele-specific inhibitors where KO hits both alleles), *non-cell-autonomous* (immunotherapy, antibody–drug
  conjugates), *untestable* (target not in CRISPR library or biomarker absent from D). Table:
  `configs/clinical_ground_truth.csv`, built by script from sources plus a curated mapping file with a citation per row;
  frozen before any CGT comparison (tag `freeze-cgt`).
- **Questions:**
  - Q4a recovery: fraction of *concordant, testable* (S/W) CGT pairs recovered by MD vs S1; direction-of-effect
    accuracy is reported as a separate number (a sign error is worse than a miss).
  - Q4b actionability: fraction of MD subtypes (and of dependent lines) whose determinant is a CGT biomarker.
  - Q4c patient population: for each recovered pair, the lineage composition of the subtype vs the approved
    indication(s); enrichment of indicated lineages (Fisher); tumour-agnostic approvals should appear cross-lineage.
    Cell-line prevalence of the biomarker per lineage vs TCGA PanCancer Atlas prevalence (cBioPortal API).
  - Q4d the discordant and non-cell-autonomous classes are expected **not** to be recovered; recovery there is a
    flag, not a success.

### 7.4 Combinatorial CRISPR screens
For MD candidates of the form "loss/silencing of A → dependency on paralog B": look up the A–B pair in Parrish 2021
(GSE178179), Dede 2020, Esmaeili Anvar 2024 (in4mer), Klingbeil 2024 (bioRxiv) and Thompson 2021 supplementary pair lists;
label replicated (synthetic-lethal hit in ≥ 1 screen), tested-not-replicated, untested. Compare to the hit rate of the
other screened pairs (matched background; inter-screen precision ≈ 0.4 is the ceiling). Where screened lines are in
DepMap, check concordance per line.

## 8. Limits: failure classes for every eligible-or-control dependency

Each gene gets one primary class (first that applies, order pre-set):
1. not variable / common essential (no between-line contrast);
2. not predictable (low tier) → no readout;
3. predicted but surrogate fidelity beyond lineage < 0.1 → "lineage-only" or "unnameable in 𝒟";
4. determinant is only the gene's own expression (EXP-own) → tautological unless backed by AMP/SIL;
5. determinant too rare (< 5 carriers) or underpowered in R;
6. confounded (culture, library, CN-cutting, co-location unresolved, hypermutation);
7. not replicated in R;
8. KO–drug discordance (from CGT classes);
9. valid biomarker candidate.
Report counts, examples, and for classes 3–8 the evidence per gene (`results/limits/gene_classes.csv`).

## 9. Verification (before relying on any result)
- **Toy data with known answers** (`tests/test_toy_pipeline.py`): simulate 600 lines × 12 lineages, 300 genes, 𝒟 of
  400 features with planted (g, d) effects (pan-lineage, lineage-restricted, co-located passenger, culture-confounded,
  hypermutation-passenger), noisy "observations" and a "predictor" = truth + noise. The pipeline must recover planted
  pan-lineage pairs (recall ≥ 0.8), label co-located passengers at locus level, flag culture-confounded pairs, and keep
  empirical FDR ≤ 0.1 under a full null.
- Refit reproduces Part 1 OOF predictions (§4.2). Unit tests on 𝒟 construction (counts match raw matrices), on locus
  collapsing (9p21 merges), and on the guard.
- Every headline number is produced by a script into `results/` and checked by an agent that did not produce it.

## 10. Milestones, compute, deadlines (24 h budget)
| Milestone | Content | Budget |
|---|---|---|
| M2 | design + ≤ 3 review rounds; freeze `freeze-design` | 3 h |
| M3 | infrastructure: guard, dictionary, covariates, toy tests, refit check | 3 h |
| M4 | surrogates + S1 + confounder audit on D; positive-control gate | 4 h |
| M5 | subtypes, cross-lineage, direct attribution agreement, combinatorial lookup | 3 h |
| M6 | CGT table frozen; clinical comparison | 3 h |
| M7 | freeze catalogue → replication in R (once) | 1 h |
| M8 | limits; REPORT.md; review | 4 h |
- Data deadline: files requested from the user that have not arrived by the start of M4 are "not available" and their
  analyses use the declared proxies (C_instab, M_burden) or are marked not done.
- Seeds: `configs/seeds.json` (fixed list). Every run logs config hash, commit, input checksums.
