# Predictor-derived dependency subtypes and their molecular biomarkers

**Status: final.** Written on 2026-10-29. The single replication in R (255 lines, one logged read, after tag `freeze-catalogue`) is done (DESIGN.md §1.3). A red-team review's corrections are applied. The budget was extended by the user (DEV-BUDGET-OVERRUN).

This is Part 2 of a two-part project. Part 1 built a CRISPR gene-essentiality predictor on DepMap 26Q3 (`../project`, read-only).

- The rules in force are listed in `docs/effective_rules.md`, cited here as ER#. Each number cites its source file.
- Claims of replication refer to the single run in R (DESIGN.md §1.3: 255 held-out lines). They are not clinical-validity claims.

## 1. Summary

1. **Question.** Can the reasons behind an accurate essentiality prediction be read out as determinants that can be assayed? And do those determinants group cell lines that share a dependency for the same reason, within and across lineages? (PROMPT.md; DESIGN §0)

2. **Catalogue.** MD-global is the set of EN-pred selected pairs that are also S1-obs full-family BH hits with the same sign.
   - It holds 1,380 gene-determinant pairs on the 1,022 dev lines (`results/benchmark/summary.json`; DESIGN.md §1.3). MD-strict, the frozen secondary, holds 11,662 (`results/benchmark/summary.json`).
   - The switch to MD-global was made because MD-strict failed its toy error bounds (A7; Section 4.3). It came before any real association result, the one exception being aggregate selection counts from a 50-gene timing pilot (`configs/design_params.json`: pilot_genes); pilot record: `results/timing/surrogate_pilot.json`.
   - By construction, MD-global is a subset of the S1-obs hits: EN-pred selects among them. Its error control is therefore inherited from S1-obs.
   - The catalogue pairs form 558 modules, 112 of them with more than one gene (`results/subtypes/summary.json`).

3. **Benchmark against single-feature tests (Q2).**
   - Selected pairs: S1-obs 3,142; EN-obs 738; EN-pred 50,728; S1-pred 276,557 (`results/benchmark/summary.json`).
   - **P1, the primary endpoint: neither "adds value" claim is supported** (`results/replication/p1_summary.json`, k_G = 100).
   - Does the predictor add value? EN-pred − EN-obs = −0.07 (bootstrap 95% CI −0.12 to 0.02): no detectable difference (`results/replication/p1_summary.json`).
   - Does the multivariate model add value? EN-obs − S1-obs = −0.10 (−0.17 to −0.05): the univariate test replicates better (`results/replication/p1_summary.json`).
   - **The predictor verdict is fragile.**
     - The bootstrap rates sit below the observed ones, so the percentile CI is biased. A basic bootstrap CI for EN-pred − EN-obs is entirely negative, (−0.16, −0.02) (`results/replication/p1_secondary.json`, bootstrap_diagnostic).
     - At deeper lists the sign flips: EN-pred replicates better than EN-obs at k = 200 (0.79 vs 0.715) and at k = 500 (0.634 vs 0.60) (`results/replication/p1_secondary.json`, secondary_k).
     - S1-obs is never beaten at k = 50, 100 or 200 (0.98, 0.99, 0.91) (`results/replication/p1_summary.json`; `results/replication/p1_secondary.json`). At k = 500, S1-pred is level with it (`results/replication/p1_secondary.json`: 0.708 vs 0.702).
     - The pre-declared verdicts stand as stated.

4. **Combinatorial screens.**
   - 31 of 183 tested MD-global paralog candidates replicate: 16.9% (CI 11.8–23.2%), against an expected 2.5% for paralog pairs matched on sequence identity and family size, about 6.8× (`results/combinatorial/validation.json`; `results/report_counts.json`: combinatorial). The unmatched background rate is 1.0% (OR 19.6) (`results/combinatorial/validation.json`).
   - S1-obs candidates replicate at the same rate (16.7%) (`results/combinatorial/validation.json`). MD-global is a subset of S1-obs, so the two lists are nested and the between-list Fisher test is descriptive only.
   - The candidates are real paralog dependencies, but the predictor adds nothing over S1-obs here. The screens are not independent of DepMap: the pairs were nominated from DepMap, and the cell lines overlap.

5. **Clinical ground truth (Q4).**
   - The positive-control gate **failed twice: 21/39 (54%) at attempt 1, and 21/39 at attempt 2** (`results/gates/G_pc_attempt1.json`, `G_pc_attempt2.json`).
   - **P2 primary (non-proxy, gate-disjoint, k = 200):** 6 testable tests, and **all four arms recover 0 of 6** (95% CI 0–0.459) (`results/clinical/q4_summary.json`, non_proxy|gate_disjoint|k200). So P2 is "descriptive (underpowered: 6)", with no McNemar test (`results/clinical/q4_summary.json`, non_proxy|gate_disjoint|k200).
   - **Secondary (16 non-proxy tests, gate controls included):** S1-obs recovers 4, EN-obs 1, EN-pred 0 and S1-pred 0 (`results/report_counts.json`: p2 hits, from `results/clinical/p2_scores.csv`). Exact McNemar EN-pred vs S1-obs p = 0.125, "no detectable difference" (`results/clinical/q4_summary.json`, non_proxy|all|k200).
   - EN-pred recovers no test at k = 200 in any non-proxy subset.

6. **Confounders (Q3).**
   - The drop-one-block audit flags 1,856 of 13,447 tested pairs, 98 of them among the 1,380 MD-global pairs (`results/report_counts.json`: confounders, counted from `results/confounders/delta_flags.csv`).
   - The fixed-threshold null gives FDR_emp 0.0009 for MD-global and 0.022 for S1-obs (`results/null/summary.json`). **NC1 fails for EN-obs (0.200)**, labelled "FDR not controlled" (`results/null/summary.json`; `configs/outcome_actions.csv`). EN-pred (0.327) is labelled "not calibrated" and exceeds 0.10 even on a null that cannot see the ridge FDP (`results/null/summary.json`; `configs/design_params.json`: NC1).
   - NC6 also fails: 0 of 18 lineage-TF genes are in class 3 (`results/limits/limits_summary.json`), so class-3 (lineage) labels are flagged unreliable.
   - NC2a fails (0.104 against ≤ 0.0011; no bug found) and NC3 fails (0.424 against ≤ 0.1; its genes are labelled CN/cutting-artefact suspects) (`results/controls/nc_real.json`). NC2b (0.26, ≤ 1.2) and NC5 (0.0, ≤ 0.2) pass (`results/controls/nc_real.json`).

7. **Cross-lineage (Q1b).**
   - Of the 1,380 MD-global pairs, 430 are cross-lineage, 22 lineage-restricted, 894 lineage-unresolvable and 34 undetermined (`results/lineage/summary.json`).
   - 139 of the 430 cross-lineage labels are "single-lineage driven": at most one nominally significant powered lineage shares the pooled sign (`results/lineage/summary.json`; A7-iii).
   - The CGT-membership odds ratio is uninformative: 5 cross-lineage vs 0 lineage-restricted CGT pairs; OR infinite (CI 0.045 to infinity), p = 1.0 (`results/limits/limits_summary.json:Q1b`).
   - R replication, cross-lineage vs lineage-restricted: 387 vs 17 replicated pairs, OR 2.64 (CI 0.72–7.95), Fisher p = 0.072 (`results/limits/limits_summary.json:Q1b`). No detectable difference (`results/limits/limits_summary.json:Q1b`).

8. **Replication in R (P1, 255 lines).** Same-sign replication at one-sided p < 0.05 of each arm's top 100 genomic candidates (`results/replication/p1_summary.json`; DESIGN.md §1.3):
   - S1-obs 0.99, S1-pred 0.93, EN-obs 0.89, EN-pred 0.82 (`results/replication/p1_summary.json`).
   - In no access stratum does EN-pred replicate detectably better than EN-obs (`results/replication/p1_by_access.csv`). Direct pairs: 0.59 vs 0.77 (`results/replication/p1_summary.json`, k = 22). Non-direct: 0.88 vs 0.86 (`results/replication/p1_summary.json`). Proxy: 0.84 vs 0.79 (`results/replication/p1_summary.json`, k = 94).
   - XLIB (S1-obs and EN-obs only; discovery on 638 D lines, replication on 200 D lines from another library): S1-obs 0.82, EN-obs 0.62 (k = 100; pool 594) (`results/benchmark/replication_xlib.csv`, `results/benchmark/xlib_summary.json`).

9. **Limits (Q5).** Gene and pair failure classes after R (`results/limits/limits_summary.json`):
   - genes (3,501): 3b 1,415; 3 884; 2 417; 9 364; 5 282; 7 70; 6 46; 1 21; 8 2 (`results/limits/limits_summary.json`; codes in `results/limits/gene_classes.csv`);
   - MD-global pairs: 9 (valid, replicated) = 552; 5 (rare or under-powered) = 537; 6 (confounded) = 160; 7 (powered, not replicated) = 131 (`results/limits/limits_summary.json`; codes in `results/limits/pair_classes.csv`).
   - The earlier 1,031 class-9 pairs and 735-pair shortlist came from missing R power for continuous features, now fixed (`results/report_counts.json`: pre_power_fix; `tests/test_limits_power.py`). Every class-9 pair is now replicated in R with power ≥ 0.5 (`results/limits/pair_classes.csv`).

   Class codes (DESIGN §12.1–12.2):
   - Gene classes: 1 no between-line contrast; 2 not predictable; 3 accurate only via lineage or library, or unnameable; 3b predictable and nameable but with no MD-strict determinant (DESIGN.md §12.1). Otherwise a gene takes 9 if any of its pairs is class 9, else its furthest pair class (DESIGN.md §12.1).
   - Pair classes: 4 expression gating only; 5 too rare or under-powered; 6 confounded; 7 not replicated; 8 KO–drug discordant; 9 valid candidate (DESIGN.md §12.2). Classes 5, 7 and 8 and class 9 are final after the R run and the data ceiling (`results/limits/limits_summary.json`).

   Limits already visible:
   - 109,919 accurate-but-unassigned dependent lines in the MD-strict gene set (`results/subtypes/summary.json`) contradict the brief's premise that every accurate prediction implies a nameable determinant.
   - 136 of the 154 CN/SIL pairs in MD-global (88%) resolve only to a locus, not to a gene (`results/locus/resolution.csv`).
   - 771 of the 1,380 MD-global pairs (56%) are expression-derived and only 285 are direct (`results/report_counts.json`: components, counted from `results/components/component_selection.csv`).
   - the gate failed;
   - every non-eligible control was missed (0/7) (`results/gates/G_pc_attempt2.json`: by_universe);
   - the df fallback leaves few-carrier binary features under-powered (ER11);
   - the own-expression control ESR1 could not be recovered under the list definitions (Section 6.2).

10. **Candidates** (`results/candidates/shortlist.csv`, `cards_summary.json`).
    - The shortlist holds 552 class-9 MD-global pairs replicated in R (`results/candidates/cards_summary.json`).
    - Expression-level determinants make up most of it: 171 expression-programme and 162 paralog-expression features (`results/report_counts.json`: shortlist). Only 15 are hotspot and 37 loss-of-function features (`results/report_counts.json`: shortlist).
    - The top hits are paralog pairs: FAM50A–FAM50B, DDX3X–DDX3Y, DNAJC19–DNAJC15, EIF1AX, RPP25L, INTS6, RAB6A, ATP6V0E1 and EFR3A. Most are not in the curated lists the cards check (`results/candidates/shortlist.csv`: not_in_curated_lists). MDM2 dependency in TP53 wild-type lines also appears.
    - Several top hits are frequently dependent in the 15 non-cancerous lines: RAB6A 0.73, MDM2 0.67, FAM50A 0.47 (`results/candidates/shortlist.csv:setA_dependent_frac`; `results/report_counts.json`: set_a). That flags possible normal-cell dependency. Set A is an unconditional fraction over 15 immortalised lines (DESIGN.md §1.3), not split by determinant, so this is not a therapeutic-window measurement.

## 2. Question and approach

Part 1's predictor, M5, is the average of three models: a kernel ridge M1, a paralog stage M2 and a LightGBM M3 (RESEARCH_TREE.md, M0). It predicts DepMap 26Q3 CRISPR effects out of fold. The brief (PROMPT.md) asks that every accurate prediction should imply a determinant. To test this, the project:

- reads determinants out of the out-of-fold predictions with a stability elastic net over an interpretable feature dictionary (EN-pred);
- compares EN-pred in a 2×2 design against the same elastic net on observed effects (EN-obs), and against a univariate confounder-adjusted test on observed (S1-obs) and predicted (S1-pred) effects (DESIGN §5.1, D1);
- names a catalogue, MD-global, and groups lines into subtypes and modules by it (DESIGN §8, A7);
- audits confounders, scores the catalogue against clinical ground truth (CIViC, MOAlmanac, FDA companion diagnostics) and combinatorial CRISPR screens, and assigns each dependency a failure class (DESIGN §9–§12);
- decides the primary endpoints on 255 held-out lines (R), read once after tag `freeze-catalogue` (DESIGN.md §1.3, §6).

There are two primary endpoints (ER24, ER25); everything else is secondary.
- **P1:** replication in R of each arm's top k_G genomic candidates (Q2).
- **P2:** CGT recall at k = 200 (Q4a).

Claimed novelty is limited to three falsifiable claims (DESIGN §0):
1. predictor-derived lists replicate better than the observed-data arms, or name determinants those arms miss;
2. line-level subtypes carry a tested cross-lineage label;
3. a dependency catalogue is scored against tiered clinical ground truth, with a decomposition of the misses.

None of the three is supported yet.

## 3. Data

The full table of sources, versions, checksums and licences is in DATA.md.

- **DepMap Public 26Q3.** CRISPRGeneEffect, Model, ModelCondition, expression, CN, the damaging and hotspot mutation matrices, fusions, CRISPRScreenMap and ScreenSequenceMap. Read in place from Part 1; checksums in `../project/configs/depmap_local_sha256.txt`; CC BY 4.0 (DATA.md, "Inherited").
- **DepMap extras supplied by the user.** OmicsSomaticMutations, OmicsGlobalSignatures, OmicsProfiles, CRISPRConfounders, OmicsCNSegmentsWGS and CCLE RRBS promoter methylation (2018-10-22), all registered before the first dictionary build (D7; DATA.md, "New"). RRBS covers 456 of the 1,022 D lines: these are in the RRBS gene matrix, and 453 of them have at least one METH feature value (`results/limits/limits_summary.json`; `results/gates/G_indep.json` notes), so methylation features are not assayed for the rest.
- **Clinical ground truth.**
  - CIViC nightly 2026-10-28 (CC0).
  - MOAlmanac v.2026-08-02. Its repository states GPL-2.0; the reuse licence for the data is unconfirmed in DATA.md.
  - The FDA companion-diagnostic list (public domain).
  - OncoKB gene lists (not redistributed).
  - In total: 64 CGT tests, 27 testable (`results/cgt/summary.json`).
- **Combinatorial screens.**
  - Esmaeili Anvar 2024, a harmonised table covering Dede, Parrish, Thompson, CHyMErA and Ito. It is the single source used; the Dede 2020 and Thompson 2021 supplements were downloaded but are not counted again (`configs/m5_m8_params.json`).
  - The Flister et al. 2024 bioRxiv supplement. Its licence is "cc_no": analysis only, and no row or value of it is written out (`results/combinatorial/validation.json:_licence_note`).
  - The Parrish processed tables could not be downloaded; counts only (DATA.md).
- **Line roles.** D = the 1,022 Part 1 dev lines; R = 255 held-out lines, never analysed; set A = 15 non-cancerous lines (DESIGN.md §1.3; ER1).
- **Not committed.** Raw and restricted data. Bulk tables stay in `cache/`, with their sha256 in `results/bulk_manifest.csv`.

## 4. Methods

### 4.1 Predictor and dictionary
- **Predictions.** Part 1's out-of-fold M5 values, fold-centred per gene. R uses the all-dev refit (DESIGN §1.1).
- **Eligible genes.** A gene is eligible if it is variable or selective, has tier high or medium, and has `r_within` ≥ 0.10 beyond lineage × library. 2,193 of 18,434 genes qualify: 1,172 high, 1,021 medium (`results/report_counts.json`: eligibility, counted from `configs/eligible_genes.csv`).
- **The dictionary 𝒟 (DESIGN §3).** Shared features are capped at 3,000 (`configs/design_params.json`). It holds these feature families:
  - hotspot, allele and loss-of-function mutations;
  - amplification, deletion and shallow loss;
  - fusions, silencing and methylation;
  - clinical composites;
  - expression dose and paralog expression;
  - expression programmes (`configs/design_params.json`: 50 NMF factors);
  - signatures, and TP53/RB composites.
- **D only.** NMF, scalers and thresholds are fitted on D only (DESIGN §1.7).
- **Access class.** Every (gene, feature) pair has one: `direct`, `potentially-direct`, `expression-derived` or `proxy` (DESIGN §4.4). It limits what "the predictor learnt d" can mean.

### 4.2 Confounder model and the four arms
**Confounder model** (DESIGN §5.3, ER9):
- Never dropped: the design blocks lineage, library group, culture block and growth.
- Droppable per test: MSI, mutation burden, TP53, own copy number and the own-expression gate.
- CR2 standard errors, clustered by patient.

**df rule.** The pre-declared reference-t df is min(n − rank − 2, ν) (`results/deviations.csv`: DEV-DF-FALLBACK). It was shown to be anti-conservative at n = 255 (DESIGN.md §13.2; `results/deviations.csv`: DEV-DF-FALLBACK). The pre-declared fallback therefore applies to binary features: df = carrier clusters − 1 (`results/deviations.csv`: DEV-DF-FALLBACK). It is used in every S1 test, in power and in the R tests (ER11; Section 5.5).

**The four arms** (DESIGN §5.1):
- **S1-obs and S1-pred.** One test per gene-feature pair, on observed or predicted effects. BH runs over the main family F1, 6,559,561 tests (`results/gates/G_indep.json`).
- **EN-obs and EN-pred.** A two-stage hierarchical elastic net: stage G genomic, stage X expression, with Frisch–Waugh–Lovell residualisation on the design block. Selection uses 25 complementary half-sample pairs, frequency ≥ 0.6 and ≥ 90% sign agreement (`configs/design_params.json`; DESIGN §4.1).
- EN-pred never reads observed effects, only the NaN mask.

All four arms ran on identical universes: 2,193 genes, 26,246,924 rows in 54 parts (`results/benchmark/summary.json`).

### 4.3 MD-global and MD-strict
**MD-strict** is the set of EN-pred selected pairs that pass an observed-effect screen run inside the EN-pred selection (BH q < 0.05, same sign; DESIGN §5.2). It was the frozen catalogue until A7.
- On the engine's toy runs with a ridge-type predictor, its true FDP was 0.228, 0.255, 0.237, 0.200 and 0.160, against a limit of 0.15 (DESIGN.md A7).
- On pure-null toys, its P(any hit) was 0.30, against a limit of 0.10 (DESIGN.md A7).
- The cause: the observed screen re-uses the noise that the ridge predictions absorbed (A7-i, A7-v).

**MD-global** is the set of EN-pred selected pairs that are also S1-obs F1 BH hits with the same sign. On the same toys it had FDP 0, 0, 0.013, 0 and 0.025, and P(any hit) 0.00, at equal recall (0.815–0.926 on the ridge toys) (DESIGN.md A7-v).

How the two lists are used:
- Both are reported. MD-global is the catalogue for everything downstream; the gate and NC2b stay on MD-strict (ER2, ER3, ER12, ER17).
- MD-strict p-values are selection-conditional, and neither list is called "significant". Confirmation comes from R and from external sources.

### 4.4 Null, negative controls and positive controls
**Fixed-threshold null** (ER23, ER15):
- Line labels are permuted within lineage × library, jointly across all dictionary rows.
- 300 pre-drawn genes; B = 20, raised to 50 only if timing allows (`configs/design_params.json`: null_s7).
- A thresholded list passes if FDR_emp_fixed ≤ 0.10.

**Negative controls.** The NC2–NC6 bounds come from five toy replicates (ER16–ER22).

**Positive controls.** 40 were specified, giving 39 gate rows in 10 families (`results/report_counts.json`: controls_spec, counted from `configs/positive_controls.csv`; DESIGN §9; RESEARCH_TREE.md D10 for 39/40):
- hotspot additions and copy-number additions;
- paralog and synthetic-lethal pairs;
- MSI/WGD signatures;
- TP53 and RB1 opposite-sign controls;
- the ESR1 expression biomarker.

The gate passes at ≥ 60% overall and ≥ 50% in each family that has at least two rows. It is decided on MD-strict (ER12).

### 4.5 Clinical and combinatorial validation
**CGT** (DESIGN §10).
- A pre-registered table of feature, target set, comparator, indication, sign and partition, built from CIViC, MOAlmanac and FDA records.
- Frozen (tag `freeze-cgt`) before any dependency data were joined.
- Headline tiers are T1–T3. Allele-defined tests use allele features, with gene-level proxies as secondary. The ceiling is S1-obs on the same test (ER26, ER28).
- Q4 metrics: recall at k = 200 per arm with exact McNemar; direction accuracy and sign-error rate; actionability; population by lineage.

**Combinatorial** (DESIGN §11).
- A candidate of the form "loss of A → dependency on paralog B" is looked up in the harmonised screen tables.
- **The per-line hit rule is an agent implementation choice (D18), not a design value.** The rule is Cohen's d ≥ 0.8 and Δmean ≤ −0.5, counted only in studies whose own hit list contains the pair. It was set before any M5 output existed.
- `replicated` means a hit in ≥ 2 screens or ≥ 2 lines (ER29).
- The comparison background is screened paralog pairs, matched on sequence-identity decile and family size (ER30).

### 4.6 Limits
Each eligible or control gene gets a failure class (DESIGN §12):
- no contrast;
- not predictable;
- accurate only through lineage or library, or unnameable;
- predictable, but with no determinant.

Each catalogue pair gets a pair class:
- expression gating only;
- too rare or under-powered;
- confounded;
- not replicated;
- KO–drug discordant;
- valid candidate.

### 4.7 Compute
- The arms run in gene blocks on a parquet dataset, under a memory guard (D17; `results/benchmark/summary.json:workers`).
- Workers: S1 8; EN 6 requested and 5 effective (`results/benchmark/summary.json:workers`); null 6 requested (`configs/runtime_params.json`: workers_default).
- The machine has 15.7 GB of RAM (RESEARCH_TREE.md D17).
- Wall times exceed the tier timeboxes. This is logged as lateness, not treated as a design change.

## 5. Verification

### 5.1 Toy simulation with an answer key
- `src/biomark/toysim.py` and `tests/test_toy_answers.py` were written by a verification agent, separate from the pipeline author. Commits carry an `Agent:` trailer (DESIGN §13.1, A6-iv, LESSONS 9).
- The toy contains planted determinants, co-located passengers, culture-confounded and lineage-only genes, and a fold-offset artefact.
- The answer key was frozen before the engine ever ran on it.

### 5.2 Suite on engine commit 04aa982: 19 pass, 3 fail (`results/toy/suite_run_04aa982.json`)
Source (`results/toy/suite_run_04aa982.json`): 3,565 s. **Failed checks:**
1. **Fold-offset artefact positive control.** It needed ≥ 3 of 5 without centring and got 0 (`results/toy/suite_run_04aa982.json`; `configs/design_params.json`: toy replicates). The bound was derived from S1-pred BH, but EN selection frequency without centring is 0.36–0.50, below 0.6 (RESEARCH_TREE.md D14; `configs/design_params.json`: selection frequency; A7-iii). Fold centring stays.
2. **Ridge-toy MD-strict.** FDP_true 0.228 (101 hits, 23 false), against a limit of 0.15 (`results/toy/suite_run_04aa982.json`). This is the double dip (A7-i).
3. **Pure-null toy, S1-pred.** P(any hit) 0.2, against a limit of 0.1 (`results/toy/suite_run_04aa982.json`). The earlier engine commit gave 0.25 (RESEARCH_TREE.md D14); the file above is for 04aa982.

MD-global passes its three error checks on the same runs (same file). These checks are in-sample: they use the same toy seeds on which MD-global was chosen (`tests/test_toy_md_global.py`, written after A7).
- ridge toy: 78 hits, 0 false (`results/toy/suite_run_04aa982.json`);
- smooth toy: 60 hits, 0 false (`results/toy/suite_run_04aa982.json`);
- pure null: P(any) = 0.0 over 20 toys (`results/toy/suite_run_04aa982.json`; `configs/design_params.json`: pure-null toys).

The fixed-threshold null cannot see the ridge-predictor FDP of the prediction-based lists (D11). EN-pred and S1-pred FDR_emp will therefore be labelled "not calibrated for ridge-type predictors", and R replication arbitrates (A6-i).

The suite was re-run on the final engine (b125298) with the same outcome: 19 passed, 3 failed, the same three failures, and MD-global 3/3 passing (RESEARCH_TREE.md); run record: `results/toy/suite_run_final.json`.

### 5.3 The MD-global switch
- The switch was decided on the toy numbers in Section 4.3 (A7 header; D15). At that point no real association result existed, apart from aggregate selection counts from a 50-gene timing pilot (`configs/design_params.json`: pilot_genes).
- Its real consequence at the gate: MD-global recovers 17/39 controls, against 21/39 for MD-strict (`results/gates/G_pc_attempt2.json:by_list`). The switch did not improve control recovery.

### 5.4 Process checks
Source: `docs/process_checks.md`. The checks are:
- `budget_check`;
- `leakscan`: it caught the planted gene string and the planted constant, and found no hits at 04aa982 (docs/process_checks.md);
- `check_frozen`: PASS. Six files that DESIGN §13.4 lists under `freeze-design` were created later; they are hashed in the `bounds` section (D13-vii);
- the DESIGN §13.2 unit tests.

Before any real result, these checks found four problems, all fixed (D13):
- the calibration failure (Section 5.5);
- a mismatch between power df and test df;
- NMF scores for D rows that were not invariant;
- hard-coded gene symbols in the clinical builder.

Remaining gaps:
- The additivity test for the Tier B attribution has no code.
- The M3b builders' timing writers never recorded wall times (D13-viii).

The fast test suite passed at the re-verification commit, 83e575b (git log).

### 5.5 Calibration fallback
The calibration test runs 20,000 null simulations per cell (`configs/design_params.json`: n_sims).
- **Under the ν rule:** at n = 255 with doubled carrier SD, size at α 0.05 was 0.072, 0.073 and 0.069 for 3, 5 and 10 carriers (`results/deviations.csv`: DEV-DF-FALLBACK). This triggered the pre-declared fallback: for binary d, df = carrier clusters − 1 (`results/deviations.csv`: DEV-DF-FALLBACK).
- **Under the fallback:** the maximum size is 0.0587 at α 0.05 and 0.0011 at α 1e-3 (`results/deviations.csv`: DEV-CAL-CONSERVATIVE). But 4 of 12 cells are conservative, below the 0.035 floor: 0.03175, 0.01040, 0.02930 and 0.02760 (`results/deviations.csv`: DEV-DF-FALLBACK, DEV-CAL-CONSERVATIVE).

The test is left failing. The cost is lower power for binary features with 3–5 carriers (`results/deviations.csv`: DEV-CAL-CONSERVATIVE). One consequence is already visible: EGFR-activating has 8 carriers and df 7, so β = −0.63 gives p = 0.13, and the control is missed at the gate (RESEARCH_TREE.md D19).

### 5.6 Independent reimplementation (DESIGN §13.3)
A separate agent, which did not read the engine modules, reimplemented S1-obs, the BH families, the gate and EN-obs selection; EN-pred was not reimplemented (`results/gates/G_indep.json`, commit 4c902fb). The independence is partial. The checker adopted the engine loader's METH_absent definition and, for one diagnostic only, plugged in the engine's λ (`G_indep.json` notes).

| Check | Result | Tolerance |
|---|---|---|
| S1-obs p, 20 genes | 59,955 of 59,955 pairs match; max abs(Δlog10 p) 8.5e-14; df, β, SE all within 1e-6 | 1e-6 |
| F1 BH | 6,559,561 tests; 3,358 rejections; 3,142 list hits (216 collinear-with-design removed); identical hit sets | exact |
| MD lists | MD-strict 11,662 and MD-global 1,380, recomputed identically | exact |
| Gate | 21/39, same rows as attempt 1 | exact |
| EN-obs frequency, 20 genes | 59,303 pairs within tolerance; max difference 0.02 | 0.04 |
| EN-obs selected set | identical (30 of 30; 11 near-threshold features excluded) | identical |
| λ (G and X) | max relative difference 4e-7 | none stated |
| **EN-obs mean β over 50 subsamples** | **FAIL: 29 of 30 within tolerance; worst 1.25e-3 relative (gene 5045 × SIL-METH:604)** | 1e-3 relative |
| Full-D β | 30 of 30 within 4e-5 | n/a |
| Confounding delta, 10 pairs (80 block rows) | 0 flag disagreements with the engine; 7 of 8 blocks match to 1.7e-12; the library block differs (DEV-DELTA-LIBRARY-BLOCK) (`results/gates/G_indep.json`; `results/deviations.csv`) | flags exact |

- **The β-tolerance check failed.** The checker attributes it to the subsample fits for gene 5045, where the active set differed (solver tolerance 1e-3 and the screening boundary) (`results/gates/G_indep.json:notes`). The full-D β of that pair agrees to 6e-10 (`results/gates/G_indep.json:notes`). No tolerance was changed.
- The file records tasks 1 and 2 as pass and task 3 as fail (`results/gates/G_indep.json`: overall).
- **Pair universe.** The checker's universe is larger than the engine's: 6,574,471 vs 6,561,731 pairs, i.e. 12,740 pairs in 75 genes that only the checker has (`results/gates/G_indep.json`). All of them are untested, because they have fewer than 5 carriers on the gene's rows (`results/gates/G_indep.json`; `configs/design_params.json`: carrier minimum). The checker carries such pairs as untested rows; the engine leaves them out of its table. The S1 p-values and the BH counts match exactly.
- **Delta check against the engine** (`results/gates/G_indep.json`: task4_engine_comparison): 0 flag disagreements over 80 block rows. The library block differs because the design is ambiguous about screen days (`results/deviations.csv`: DEV-DELTA-LIBRARY-BLOCK).

## 6. Results on the dev set D

### 6.1 Benchmark list sizes
The run covered 2,193 genes and 26.2 million (g, d) rows, all four arms in 45 minutes (RESEARCH_TREE.md M4; D19). Selected pairs (`results/benchmark/summary.json`):

| List | Pairs |
|---|---|
| S1-obs | 3,142 |
| EN-obs | 738 |
| EN-pred | 50,728 |
| S1-pred | 276,557 |
| MD-global | 1,380 |
| MD-strict | 11,662 |

The prediction arms select far more than the observed arms. With a ridge-type predictor their error rates are not controlled and cannot be read from the null (A6-i), so list size is not a performance claim. Replication of these lists in R: Section 7.7.

### 6.2 Positive-control gate: failed at attempt 1 and at attempt 2

**Attempt 1** (`results/gates/G_pc_attempt1.json`): MD-strict recovered **21/39 = 0.538 (95% CI 0.372–0.699)**, below the 60% bar.
- Families below 50%: ADD-CN 1/5, SIG 0/2, W 1/3 (`results/gates/G_pc_attempt1.json`: by_family; `configs/design_params.json`: gate_s9).
- Families passing: ADD-HOT 5/6, LOF-SL 1/2, P53 4/4, PAR 7/13, PROG 2/2 (`results/gates/G_pc_attempt1.json`: by_family).
- ADD-EXP (0/1) and SIL (0/1) have fewer than two rows each (`results/gates/G_pc_attempt1.json`: by_family).

**Bug search (D19).** The U_ctrl universe, which holds the non-eligible control and CGT-target genes, had never been built or read. As a result SMARCA2, ARID1B, MYCN and ESR1 had no rows. This was fixed as a DESIGN §9.3 bug, with a toy test that fails before the fix and passes after it (`test_gate_u_ctrl.py`). The 41 U_ctrl genes were run through the arms, and the gate was rerun as attempt 2 (`results/benchmark_ctrl/summary.json`).

**Attempt 2** (`results/gates/G_pc_attempt2.json`): MD-strict **21/39 again**, with the same three failing families.

| Breakdown | Recall |
|---|---|
| Eligible universe | 21/32 (0.656, CI 0.468–0.814) |
| U_ctrl universe | **0/7** (CI 0–0.410) |
| S1-pred | 29/39 (list membership; rewards list size: 276,557 pairs) (`results/gates/G_pc_attempt2.json`; `results/benchmark/summary.json`) |
| MD-strict | 21/39 |
| S1-obs | 21/39 |
| EN-pred | 21/39 |
| MD-global | 17/39 |
| EN-obs | 15/39 |
| Access class: direct | 4/10 |
| Access class: potentially-direct | 9/17 |
| Access class: expression-derived | 3/4 |
| Access class: proxy | 5/8 |
| PAR rows with a non-direct encoding | 4/8 |

**The seven U_ctrl controls are missed by MD-strict, MD-global and S1-obs. There is no further bug (D20, corrected in D21).** Four are recovered by other arms (`results/gates/G_pc_attempt2.json`): S1-pred recovers MYCN, SMARCA2 and ARID1B, and EN-obs recovers ESR1 and SMARCA2.
- **SMARCA2 × SMARCA4 loss-of-function:** β = −0.48, p = 1e-5, q = 0.15 (RESEARCH_TREE.md D20).
- **ARID1B × ARID1A loss-of-function:** q = 0.099 in the U_ctrl BH family (RESEARCH_TREE.md D20). For MD-strict, both misses arise because EN-pred did not select the pair.
- **ENO2, PAPSS1 and ME3 × shallow loss:** β ≈ 0 (RESEARCH_TREE.md D20). Shallow loss does not capture the homozygous deletions these pairs need.
- **MYCN amplification:** p = 0.01, absorbed by lineage (RESEARCH_TREE.md D20).
- **ESR1 own-expression:** q = 0.0008 (RESEARCH_TREE.md D20), but own-expression belongs to test family F2, so it cannot enter the F1 lists (S1-obs, MD-strict, MD-global). EN-obs does recover it. It is reported as a specification gap and not re-scored.

Other misses outside U_ctrl: MAT2A × MTAP deletion (q = 0.15), and EGFR-activating as in Section 5.5 (RESEARCH_TREE.md D19).

Eligible-only recall (21/32) is above 60%, but it is not the gate's denominator, and it was not used to change the verdict (`results/gates/G_pc_attempt2.json`: by_universe; `configs/design_params.json`: gate_s9). The wording was fixed in advance (`configs/outcome_actions.csv`): "gate failed at attempt 1 (21/39) and at attempt 2 after the U_ctrl fix (21/39); CGT headline restricted to gate-disjoint tests." MD-global's 17/39 is reported for information only (`configs/outcome_actions.csv`: gate_md_global; `results/gates/G_pc_attempt2.json`: by_list).

### 6.3 Subtypes and modules (`results/subtypes/summary.json`, MD-global and MD-strict)

| | MD-global | MD-strict |
|---|---|---|
| Pairs | 1,380 | 11,662 |
| Pairs with member lines | 1,014 | 8,299 |
| Modules | 558 | 2,221 |
| Multi-gene modules | 112 | 694 |

- The file's gene-level block (2,134 genes; 279,074 unexplained lines; 109,919 accurate-but-unassigned lines) describes the MD-strict gene set (`results/subtypes/summary.json`): `results/subtypes/subtype_genes.csv` has 2,134 rows, while only 764 genes carry an MD-global pair (`results/report_counts.json`: md_lists). There is no per-list split of these line counts. The top-level `n_pairs` (11,662) and `n_pairs_with_members` (8,299) are MD-strict as well (`results/subtypes/summary.json`).
- Subtype membership uses the share rule (DESIGN §8.1). An "unexplained" line is dependent but mispredicted.
- Cross-lineage labels and the co-location, amplicon, culture and hypermutation audits are in Sections 7.2 and 7.3.

### 6.4 Combinatorial validation (`results/combinatorial/validation.json`)
The test covers paralog candidates (loss of A → dependency on B), screened in 23 cell lines across six sources (`results/combinatorial/validation.json`: _n_cell_lines; `results/report_counts.json`: combinatorial).

| List | Candidates | Tested | Replicated | Rate (95% CI) | Background rate | OR vs background (95% CI) |
|---|---|---|---|---|---|---|
| MD-global | 191 | 183 | 31 | 16.9% (11.8–23.2) | 1.03% (110/10,717) | 19.6 (12.3–30.6) |
| MD-strict | 363 | 317 | 37 | 11.7% (8.4–15.7) | 0.78% (114/14,539) | 16.7 (11.0–24.9) |
| S1-obs | 210 | 198 | 33 | 16.7% (11.8–22.6) | 1.01% (108/10,702) | 19.6 (12.5–30.2) |

- **Between lists** (descriptive only: MD-global and MD-strict are nested in S1-obs, so the lists are not independent samples). MD-global vs S1-obs: OR 1.02 (0.57–1.81), Fisher p = 1.0 (`results/combinatorial/validation.json`). MD-strict vs S1-obs: OR 0.66 (0.39–1.14), p = 0.11 (`results/combinatorial/validation.json`).
- **Matched expectation.** For MD-global the sequence-identity-matched expected rate is 2.5%, so the enrichment over that expectation (16.9% vs 2.5%, a derived ratio of 6.8) is smaller than the OR over the unmatched background suggests (`results/report_counts.json`: combinatorial; `results/combinatorial/validation.json`).
- The hit rule is an agent choice (Section 4.5).
- **Conclusion.** Catalogue paralog candidates are experimentally supported combinatorial dependencies, but the test shows no advantage over the univariate S1-obs list.

### 6.5 Worked subtype examples

The five MD-global pairs below show what a recovered subtype looks like. All five are class 9 and replicated in R, and each is a row of `results/candidates/shortlist.csv`.
- They were picked to cover the cases the brief asks about, not by effect size. The exception is FAM50A, which is shortlist rank 1 (`results/candidates/shortlist.csv`).
- **Members** are dependent lines that the predictor ranks as dependent and whose share of the explained signal is at least 0.25 (`results/subtypes/subtypes.csv`; DESIGN.md §8.1).
- **Lineage composition** gives the five largest lineages in that file.
- **Powered lineages** counts lineages with per-lineage power of at least 0.5 (`results/lineage/cross_lineage.csv`; DESIGN.md §8.3).
- **Set-A fractions** are over the 15 non-cancerous lines (DESIGN.md §1.3).
- **R effects** are one-sided, from `results/replication/r_tests.csv`.

**Example 1: FAM50A × expression of its paralog FAM50B (cross-lineage, expression determinant).**
- **Determinant.** FAM50A is less essential where its paralog FAM50B is highly expressed: EXPPAR:26240, pooled β +0.277 (`results/lineage/cross_lineage.csv`). The access class is direct, because paralog expression is a hand-coded input of the predictor's paralog stage (`results/candidates/shortlist.csv`; DESIGN §4.4).
- **Members.** 174 of 176 accurate-dependent lines, out of 256 dependent lines (`results/subtypes/subtypes.csv`). They span 21 lineages; the largest are Lymphoid 32, CNS/Brain 25, Bowel 24, Esophagus/Stomach 17 and Myeloid 13 (`results/subtypes/subtypes.csv`).
- **Lineage label.** Cross-lineage: 21 powered lineages, same-sign fraction 1.0 (`results/lineage/cross_lineage.csv`). The interaction p of 0.00099 shows that the effect size varies by lineage (`results/lineage/cross_lineage.csv`).
- **R.** β = 0.278 (SE 0.020), p = 1.6e-26 (`results/replication/r_tests.csv`).
- **Other support.** Combinatorial screens: replicated (`results/candidates/shortlist.csv`). The silencing and methylation variants of the same paralog are also on the shortlist (`results/candidates/shortlist.csv`).
- **Set A.** Dependent fraction 0.467 (`results/candidates/shortlist.csv`).
- **Reading.** Low FAM50B marks FAM50A-dependent lines in 21 lineages (`results/lineage/cross_lineage.csv`). But the dependency is also common in non-cancerous lines, so it is not tumour-selective. The "subtype" is nearly the whole dependent set. Because the access class is direct, the result shows that the predictor used an input it was given, not a determinant beyond its inputs.

**Example 2: MDM2 × TP53 mutation (cross-lineage; the determinant is shared by a 71-gene module) (`results/subtypes/modules.csv`).**
- **Determinant.** MDM2 is less essential in TP53-mutant lines: COMP:TP53MUT, pooled β +0.740 (`results/lineage/cross_lineage.csv`). The access class is proxy, because TP53 status is not a predictor input (`results/candidates/shortlist.csv`).
- **Members.** The subtype files list no member lines (`results/subtypes/subtypes.csv`). They define members among carriers of the feature, and here the dependent lines are the non-carriers (TP53 wild type).
- **Module.** COMP:TP53MUT with positive sign covers 71 genes (`results/subtypes/modules.csv`). These include MDM2, MDM4, PPM1D and USP7, but also ribosomal and exosome genes (`results/subtypes/modules.csv`).
- **Lineage label.** Cross-lineage across 19 powered lineages (`results/lineage/cross_lineage.csv`). The interaction p is 1.2e-12 (`results/lineage/cross_lineage.csv`).
- **R.** β = 0.834 (SE 0.089), p = 2.9e-17 (`results/replication/r_tests.csv`).
- **Other support.** MDM2 is a P53-family gate control (`results/gates/G_pc_attempt2.json`).
- **Set A.** Dependent fraction 0.667 (`results/candidates/shortlist.csv`).
- **Reading.** This recovers a known TP53-status dependency. It is strongly heterogeneous across lineages, and the module mixes p53-pathway genes with general ribosome and exosome genes, so it is not 71 independent dependencies (`results/subtypes/modules.csv`).

**Example 3: KRAS × KRAS hotspot mutation (cross-lineage, genomic determinant).**
- **Determinant.** Lines with a KRAS hotspot mutation depend more on KRAS: HOT:3845, pooled β −1.03 (`results/lineage/cross_lineage.csv`). The access class is potentially-direct (`results/candidates/shortlist.csv`).
- **Members.** 95 of 187 accurate-dependent lines (`results/subtypes/subtypes.csv`). The largest lineages are Bowel 30, Lung 21, Pancreas 14, Biliary Tract 7 and Esophagus/Stomach 6 (`results/subtypes/subtypes.csv`).
- **Lineage label.** Cross-lineage, with 10 powered lineages (`results/lineage/cross_lineage.csv`).
- **R.** β = −1.080 (SE 0.140), p = 3.0e-9 (`results/replication/r_tests.csv`).
- **Other support.** KRAS is a gate control (`results/candidates/known_recovery.csv`).
- **Set A.** Dependent fraction 0.4 (`results/candidates/shortlist.csv`).
- **Reading.** This is the cleanest genomic subtype, and it shows the pipeline returns a known biology. But the call is gene-level: the allele-specific KRAS G12C clinical test is not recovered at k = 200 (`results/candidates/known_recovery.csv`). The 92 accurate-dependent non-member lines (derived) are not explained by this feature (`results/report_counts.json`: subtype_examples).

**Example 4: PRMT5 × MTAP loss (cross-lineage, multi-gene module with a clinical test).**
- **Determinant.** MTAP loss (CLIN:MTAP_LOSS) makes lines more dependent on PRMT5, pooled β −0.329 (`results/lineage/cross_lineage.csv`). The access class is proxy (`results/candidates/shortlist.csv`).
- **Members.** 47 of 112 accurate-dependent lines (`results/subtypes/subtypes.csv`). The largest lineages are Pancreas 11 and Esophagus/Stomach 8 (`results/subtypes/subtypes.csv`).
- **Module.** The MTAP-loss module groups PRMT5, PRMT1, PELO, METTL16 and WDR77 (`results/subtypes/modules.csv`).
- **Lineage label.** Cross-lineage, with 10 powered lineages and interaction p = 0.386, so the effect is homogeneous across lineages (`results/lineage/cross_lineage.csv`).
- **R.** β = −0.160 (SE 0.061), p = 0.0062 (`results/replication/r_tests.csv`). That is about half the D estimate.
- **Clinical.** CGT test MTAP_LOSS__PRMT5__pan, tier T2. S1-obs ranks it 139 (a hit at k = 200); EN-pred ranks it 7,961 (`results/clinical/p2_scores.csv`).
- **Set A.** Dependent fraction 0.067 (`results/candidates/shortlist.csv`).
- **Reading.** This is the one example with a curated clinical test. Its MD-global entry comes from the S1-obs half of the intersection, not from a determinant that the predictor revealed.

**Example 5: CCND1 × expression of CCND3 (lineage-restricted).**
- **Determinant.** CCND1 is less essential where its paralog CCND3 is highly expressed: EXPPAR:896, pooled β +0.127 (`results/lineage/cross_lineage.csv`).
- **Members.** 23 of 147 accurate-dependent lines (`results/subtypes/subtypes.csv`).
- **Lineage label.** Lineage-restricted, with only 1 powered lineage (`results/lineage/cross_lineage.csv`). The significant lineage is Lung (β 0.307, q 0.0064), but its power of 0.42 is below the 0.5 cut (`results/lineage/cross_lineage.csv`).
- **R.** β = 0.157 (SE 0.044), p = 4.4e-4 (`results/replication/r_tests.csv`).
- **Set A.** Dependent fraction 0.133 (`results/candidates/shortlist.csv`).
- **Reading.** The label rests on a lineage that fails the power cut, so this is a lineage-concentrated candidate, not a confirmed restriction.

**Across vs within lineage** (`results/report_counts.json`: lineage_by_class; MD-global pairs by cross-lineage label, all pairs vs class 9):

| Label | MD-global pairs | Class 9 | Class-9 share |
|---|---|---|---|
| Cross-lineage | 430 | 314 | 73% |
| Lineage-restricted | 22 | 14 | 64% |
| Lineage-unresolvable | 894 | 196 | 22% |
| Undetermined | 34 | 28 | 82% |

- **Caveat.** The cross-lineage label needs ≥ 2 powered lineages, and class 9 needs R power ≥ 0.5, so the higher class-9 share of cross-lineage pairs is confounded with power (DESIGN.md §8.3, §12.2).
- DDX3X, at shortlist ranks 2, 4 and 14, is not one of the five examples; it carries the `cn_cutting_suspect` label (`results/candidates/shortlist.csv`; Section 7.8).

## 7. Further results, including the replication in R

### 7.1 Confounder audit
Drop-one-block deltas on MD-global ∪ S1-obs hits ∪ controls, including the U_ctrl gate controls.
- Blocks: lineage, library, culture, growth, burden, TP53 and own CN.
- A pair is flagged if abs(Δβ) > 0.5·abs(β_minus) and > SE_minus.
- Also VIF and subtype-confounding.

The audit covered 13,447 pairs: 13,438 eligible plus 9 U_ctrl gate-control pairs (`results/report_counts.json`: confounders). **1,856 are delta-flagged, 610 subtype-confounded and 2,319 flagged by either** (`results/report_counts.json`: confounders, counted from `results/confounders/delta_flags.csv`).

| Source list (`results/report_counts.json`: confounders.by_list) | Pairs | Delta-flagged | Subtype-confounded |
|---|---|---|---|
| MD-global | 1,380 | 98 | 23 |
| MD-strict | 11,662 | 1,571 | 596 |
| S1-obs BH hits | 3,142 | 382 | 37 |
| Controls | 51 | 1 | 0 |

(Same file. The lists overlap, so the rows do not sum.)

- Of the 1,856 flagged pairs, 1,797 are flagged on the lineage block alone and 32 on own CN alone (`results/report_counts.json`: confounders.n_delta_flagged_by_blocks, from `results/confounders/delta_flags.csv`: blocks_flagged).
- VIF: max 6.17 (library_1) against a limit of 10; condition number 10.84; 41 columns; culture_4 and culture_9 aliased (`results/covariates/vif_check.json`).

### 7.2 Cross-lineage labels (Q1b)
Labels follow DESIGN §8.3, over 11,757 pairs (`results/lineage/summary.json`):

| List (pairs) | Cross-lineage | Lineage-restricted | Lineage-unresolvable | Undetermined |
|---|---|---|---|---|
| MD-global (1,380) | 430 | 22 | 894 | 34 |
| MD-strict (11,662) | 486 | 60 | 11,075 | 41 |
| Controls (51) | 24 | — | 27 | — |
| Top 200 S1-obs BH hits | 200 | — | — | — |

- 894 of the 1,380 MD-global pairs (65%, derived) are lineage-unresolvable, because fewer than 2 lineages are powered (`results/lineage/summary.json`; `results/report_counts.json`: lineage; DESIGN.md §8.3).
- Single-lineage-driven cross-lineage labels: 139 (MD-global), 160 (MD-strict), 5 (controls), 1 (S1-obs BH) (`results/lineage/summary.json`).
- Q1b, CGT membership, cross-lineage vs lineage-restricted: 5 vs 0 pairs; OR infinite (CI 0.045 to infinity); Fisher p = 1.0 (`results/limits/limits_summary.json:Q1b`). Uninformative (`results/limits/limits_summary.json:Q1b`).
- R replication, cross-lineage vs lineage-restricted: 387 vs 17 replicated; OR 2.64 (0.72–7.95), p = 0.072 (`results/limits/limits_summary.json:Q1b`).
- Known weakness of the rule: a one-lineage effect can be labelled cross-lineage when the null lineages happen to share its sign (RESEARCH_TREE.md D13-vi: 1 of 5 toy seeds).

### 7.3 Co-location, cutting artefact, culture, hypermutation
DESIGN §8.4–8.7:
- **Co-location** (`results/m5_audits_summary.json`, `results/locus/resolution.csv`): of 1,114 CN/SIL determinants in locus groups, 64 resolve to a gene and 1,050 stay at locus level. For MD-global pairs, 18 of 154 resolve to a gene and 136 stay at locus level (`results/report_counts.json`: locus, counted from `results/locus/resolution.csv`).
- **CN-cutting artefact** (`results/m5_audits_summary.json`; `results/locus/amplicon.csv`): of 190 amplicon determinants, 57 are amplified-driver addiction, 133 unresolved and **0 cutting artefacts**. For MD-global: 12 addiction and 7 unresolved, of 19 (`results/report_counts.json`: locus).
- **Culture** (`results/m5_audits_summary.json`; run on the MD-strict subtype set, which contains MD-global): of 11,662 subtypes, 60 are flagged (Cramér's V ≥ 0.3), 51 rescued by the within-stratum test and 9 confounded; 14 module-overlap flags.
- **Hypermutation and TP53** (`results/m5_audits_summary.json`): 466 determinants tested, 5 failed (4 TP53-mutant, 1 hypermutation) (`results/m5_audits_summary.json`; `results/report_counts.json`: m5_audits).
- **Component readout** (`results/components/component_selection.csv`).
  - MD-global pairs by layer: expression-mediated 807, both 487, genomics layer 68, none 18 (of 1,380) (`results/report_counts.json`: components).
  - By access class: expression-derived 771, direct 285, proxy 168, potentially-direct 156 (`results/report_counts.json`: components).
  - Tier B direct attribution was not run.

### 7.4 Clinical ground truth: recovery, direction, actionability, population
The headline set is gate-disjoint (ER14, ER26).
**P2 primary (non-proxy, gate-disjoint, k = 200).** 6 testable tests (5 T1, 1 T3) (`results/clinical/q4_summary.json`; `results/report_counts.json`: p2).
- Every arm recovers 0 of 6: recall 0 (95% CI 0–0.459) (`results/clinical/q4_summary.json`, non_proxy|gate_disjoint|k200).
- Verdict "descriptive (underpowered: 6)": no McNemar, and the trend test was not run (`results/clinical/q4_summary.json`).
- At k = 100 every arm is also 0/6 (`results/report_counts.json`: p2). At k = 500, only S1-obs recovers one test (1/6 = 0.167, CI 0.004–0.641) (`results/clinical/q4_summary.json`, non_proxy|gate_disjoint|k500).

**Secondary sets** (`results/clinical/q4_summary.json`; `results/report_counts.json`: p2; `configs/primary_endpoints.json`: k; not gate-disjoint, unadjusted, gate controls included):

| Set, k | Tests | S1-obs | EN-obs | EN-pred | S1-pred | McNemar EN-pred vs S1-obs |
|---|---|---|---|---|---|---|
| non-proxy all, 200 | 16 | 4 (0.25) | 1 (0.063) | 0 | 0 | p = 0.125, no detectable difference |
| non-proxy all, 500 | 16 | 5 (0.313) | 2 (0.125) | 0 | 0 | p = 0.0625 |
| incl-proxy all, 200 | 22 | 7 (0.318) | 2 (0.091) | 0 | 0 | p = 0.0156, univariate recovers more |
| incl-proxy all, 500 | 22 | 9 (0.409) | 5 (0.227) | 3 (0.136) | 0 | p = 0.031, univariate recovers more |

(Same file.)
- The significant McNemar results are secondary and include gate controls, so they do not enter the P2 claim. They point away from the predictor-derived arm.
- EN-pred and S1-pred recover no test at k = 200 in any non-proxy subset.

**Direction.** Every recovered test has the expected sign: direction accuracy 1.0 and 0 sign errors for S1-obs, EN-obs, and EN-pred at k = 500 incl-proxy (`results/clinical/q4_summary.json`). Direction is not estimable for arms with no hits (`results/clinical/q4_summary.json`).

**Data ceiling** (`results/clinical/cgt_ceiling.csv`; `results/clinical/ceiling_summary.json`; DESIGN §10.4). S1-obs is run on each ceiling-testable test in its own comparator (§5.3 model, one-sided, BH q < 0.1 over the tests, F7).
- 28 tests are ceiling-testable: 19 data-supported, 9 data-null and 0 data-contradicted (`results/clinical/ceiling_summary.json`). The other 36 are untestable (`results/clinical/ceiling_summary.json`).
- The one-sided Wilcoxon robustness check agrees with the label in 21 of 28 (`results/clinical/ceiling_summary.json`).
- Data-null tests: EGFR-activating in NSCLC (and its HOT proxy), ALK fusion, BCR::ABL1, RAS-mutant → EGFR in CRC, ER+ → ESR1, CCNE1 amplification → PKMYT1, and the refuted KRAS-mutant → MEK and EGFR-amplified GBM tests (`results/clinical/cgt_ceiling.csv`).
- **Caveat.** For BCR::ABL1 and ALK fusion the carrier indicator is collinear with the design blocks inside the comparator, so the fit is degenerate (BCR::ABL1 β ≈ −2e14) (`results/clinical/cgt_ceiling.csv`). Both are flagged and, as in S1, cannot be data-supported; for BCR::ABL1 the Wilcoxon q is 0.002 (`results/clinical/cgt_ceiling.csv`).
- Stratified comparators get no separate indication fixed effect; the §5.3 lineage block is used (`results/clinical/ceiling_summary.json`: rules).
- **Class 8 from the ceiling:** concordant tests with a data-null ceiling add 7 gene-level class-8 rows on 4 genes; ESR1 and PKMYT1 end in gene class 8 (`results/report_counts.json`: ceiling_class8).

**Miss decomposition** (`results/clinical/miss_decomposition.csv`; `results/report_counts.json`: miss_testable). The category order (ineligible gene, no feature_D, ceiling data-null, lost to a group member, ranked below k) is an implementation choice; the design lists the categories only (`results/clinical/cgt_secondary.json`: rules). Misses among testable tests at k = 200:

| Set (testable tests) | Arm | Misses | Ceiling data-null | Lost to a group member | Ranked below k |
|---|---|---|---|---|---|
| headline, gate-disjoint (6) | S1-obs | 6 | 1 | 2 | 3 |
| headline, gate-disjoint (6) | EN-obs | 6 | 1 | 2 | 3 |
| headline, gate-disjoint (6) | EN-pred | 6 | 1 | 1 | 4 |
| headline, gate-disjoint (6) | S1-pred | 6 | 1 | 1 | 4 |
| non-proxy all (16) | S1-obs | 12 | 3 | 3 | 6 |
| non-proxy all (16) | EN-obs | 15 | 3 | 5 | 7 |
| non-proxy all (16) | EN-pred | 16 | 3 | 5 | 8 |
| non-proxy all (16) | S1-pred | 16 | 3 | 5 | 8 |

- The untestable headline tests add 14 ineligible-gene and 5 no-feature misses per arm (`results/clinical/cgt_secondary.json`: miss_decomposition).
- **Recall over data-supported tests** (non-proxy all, k = 200; 11 tests): S1-obs 3/11, every other arm 0/11 (`results/clinical/cgt_secondary.json`: recall_data_supported).
- **NC5** (false-sensitivity rate): 0.0, a pass against ≤ 0.2 (`results/clinical/cgt_secondary.json`: NC5). Only 2–3 tests per arm are ranked (`results/clinical/cgt_secondary.json`: NC5; `results/controls/nc_real.json`).
- **Refuted vs approved** (k = 200): refuted 0/3 in every arm; approved (T1 sensitivity, non-proxy, testable) S1-obs 3/10, EN-obs 1/10, EN-pred 0/10, S1-pred 0/10 (`results/clinical/cgt_secondary.json`: refuted_vs_approved).
- The TCGA prevalence (Tier B) was cut.

**Known-recovery rows** (`results/candidates/cards_summary.json`; `results/report_counts.json`: known_recovery): 88 rows, 58 controls and 30 CGT tests.
- Gate controls: MD-global 17/39, MD-strict 21/39 (`results/report_counts.json`: known_recovery).
- Lineage-TF class-check rows: 0 of 18 recovered (`results/report_counts.json`: known_recovery).
- CGT rows: MD-global 5/30, MD-strict 7/30 (`results/report_counts.json`: known_recovery).

**Actionability** (`results/clinical/actionability.csv`; `configs/design_params.json`: actionability_null_permutations; 200 within-list gene-label permutations):

| List | n | Pair-level | Permutation mean | Permutation p | Feature-level |
|---|---|---|---|---|---|
| MD-global | 754 | 0.0080 | 0.00030 | 0.005 | 0.069 |
| S1-obs | 2,464 | 0.0057 | 0.00017 | 0.005 | 0.051 |
| MD-strict | 5,044 | 0.0014 | 0.000048 | 0.005 | 0.054 |

A permutation p of 0.005 is the floor for 200 permutations (`results/clinical/actionability.csv`; `configs/design_params.json`: actionability_null_permutations). No test of MD-global vs S1-obs was run.

**Population** (`results/clinical/q4c_labels.csv`, `results/clinical/q4c_population.csv`). Of 16 tests: 6 on-label, 4 beyond-label, 3 not supported, 2 off-label only, 1 indication-untestable (`results/report_counts.json`: q4c).
- Six indication codes could not be mapped (BALL, BNHL, BTC, DLBCL, EGC, FL), because there is no full OncoTree file.
- BRAF V600, pan-cancer: β = −0.870, one-sided p = 2.9e-16, 61 carriers, on-label (`results/clinical/q4c_population.csv`; `results/clinical/q4c_labels.csv`).
- EGFR-activating in NSCLC is "not supported" (indicated-lineage precision 0.0) (`results/clinical/q4c_labels.csv`), as are CCNE1 amplification → PKMYT1 and KMT2A-r → MEN1.

Tier counts among the 64 tests: 41 T1, 7 T2, 5 T3, 1 preclinical, 3 refuted, plus 7 proxy tests (`results/cgt/summary.json`). The headline set shrank from 18 scored tests (`results/cgt/summary.json`) to 16 non-proxy testable, and then to 6 gate-disjoint non-proxy (`results/clinical/q4_summary.json`). `cgt/summary.json` expected 8 gate-disjoint tests; the difference of 2 is the proxy tests (`results/cgt/summary.json`; `results/report_counts.json`: p2).

### 7.5 Null calibration (FDR_emp)
Fixed-threshold null on the 300-gene subset, B = 20 draws, 3 workers (`results/null/summary.json`; `results/null/fdr_fixed.csv`). B stayed at 20: draw 0 took 706 s, so 50 draws would take about 589 min (derived), over the 180-min cap (`results/report_counts.json`: null; `configs/design_params.json`: null_s7).

| List | Observed hits (subset) | Mean hits per permutation | FDR_emp_fixed (95% CI) | P(any hit in a draw) | ≤ 0.10 | Label (`configs/outcome_actions.csv`) |
|---|---|---|---|---|---|---|
| S1-obs | 515 | 11.55 | 0.0224 (0.0185–0.0269) | 1.0 | yes | — |
| S1-pred | 38,482 | 1,570.05 | 0.0408 (0.0396–0.0420) | 1.0 | yes | not calibrated (A6/A7) |
| EN-obs | 107 | 21.35 | **0.1995 (0.1813–0.2196)** | 1.0 | **no** | **FDR not controlled** |
| EN-pred | 6,979 | 2,284.45 | **0.3273 (0.3229–0.3318)** | 1.0 | **no** | not calibrated (A6/A7) |
| MD-strict | 1,701 | 63.1 | 0.0371 (0.0350–0.0389) | 1.0 | yes | not calibrated (A6/A7) |
| MD-global | 213 | 0.2 | 0.0009 (0–0.0021) | 0.15 | yes | — |

**NC1.**
- **Failed for EN-obs** (`results/null/summary.json`; `configs/design_params.json`: NC1; FDR_emp 0.20 > 0.10). The list is labelled "FDR not controlled"; no threshold was changed.
- EN-pred also exceeds 0.10 (`results/null/summary.json`; `configs/design_params.json`: NC1; 0.327).
- For S1-pred, EN-pred and MD-strict, the null cannot be read as an FDR. On ridge toys their true FDP was 0.7–0.8 (RESEARCH_TREE.md D11: EN-pred and S1-pred) and 0.16–0.26 (RESEARCH_TREE.md A7: MD-strict), and the permutation null does not see it.
- The selection-adaptive re-run was cut (Tier B; `fdr_emp_rerun: null`).

**Observed-count bug, found and fixed.** The first summary had MD-global observed hits = 0, because `md_global` was missing from the observed inputs (`results/null/summary.json:resummarised`).
- The summary was recomputed from the same draws (draws committed at a7c6c8a, resummarise commit 732e50e; `summary.json:resummarised`).
- The only changed value is MD-global's observed count (0 → 213) (`results/null/summary.json:resummarised`). The permutation counts and the other lists are unchanged.

**Negative-control bounds** (`results/toy/nc_bounds.json`):
- NC2a ≤ 0.0011 and NC2b ≤ **1.2** (`configs/frozen_hashes.json`: bounds_values). NC2b is set by one ridge replicate; its per-replicate values were 0, 0.2, 0.2, 0.2 and 1.2, so the bound is very loose (`results/toy/nc_bounds.json`).
- NC3 0.1, NC5 0.2 and NC6 0.7 keep their provisional values (`configs/frozen_hashes.json`: bounds_values).

**Real-data values.**
- **NC2a fails**: the median EN-pred fidelity of the 100 low-tier genes is 0.104, against ≤ 0.0011 (`results/controls/nc_real.json`). Per `outcome_actions` a bug check was done and found none: the statistic runs through the same code as for the eligible genes, whose median fidelity is 0.339 (`results/report_counts.json`: nc_detail). The toy bound came from low-tier genes whose signal is latent only, which real low-tier predictions are not.
- **NC2b passes**: 0.26 MD-strict hits per low-tier gene, against ≤ 1.2; MD-global 0/100 for information (provisional 0.1) (`results/controls/nc_real.json`).
- The NC2 genes are "matched-SD" by an implementation choice: each of 100 eligible genes drawn with `lowtier_seed` is matched to the unused low-tier gene with the nearest D effect SD (`results/controls/lowtier_genes_summary.json`; list in `configs/lowtier_genes.txt`).
- **NC3 fails**: 0.424 of the 684 testable self-null pairs are dependent, against ≤ 0.1, driven by DEL self-null pairs (0.527 of 478) (`results/controls/nc_real.json`; `results/report_counts.json`: nc_detail). Per `outcome_actions` the affected genes are labelled CN/cutting-artefact suspects: 46 genes, 33 MD-global pairs and 13 of the 552 shortlist pairs (`results/limits/limits_summary.json`: cn_cutting_suspect).
- **NC5 passes**: 0.0, against ≤ 0.2 (`results/controls/nc_real.json`).
- NC4 is the toy fold-offset control (A7-iii).
- **NC6 fails**: 0 of 18 lineage-TF genes are in class 3, against a bound of ≥ 0.7 (`results/limits/limits_summary.json`). The 17 genes with r_within ≥ 0.10 are classes 3b or 9, so the class-3 rule does not fire (`results/report_counts.json`: lineage_tf). Per `outcome_actions`, class-3 (lineage) labels are flagged unreliable.

### 7.6 XLIB (library-independent replication, S1-obs and EN-obs only)
Rates, S1-obs vs EN-obs, with the EN-obs − S1-obs contrast (`results/benchmark/replication_xlib.csv`; `results/benchmark/xlib_summary.json`):
- G: 0.82 vs 0.62, EN-obs − S1-obs −0.20 (−0.23, −0.07), univariate better (`results/benchmark/replication_xlib.csv`).
- direct (k 16): 0.625 vs 0.4375, descriptive (`results/benchmark/replication_xlib.csv`).
- non-direct: 0.83 vs 0.63, −0.20 (−0.25, −0.07), univariate better (`results/benchmark/replication_xlib.csv`).
- proxy (k 50): 0.78 vs 0.66, −0.12 (−0.18, 0.02), no detectable difference (`results/benchmark/replication_xlib.csv`).

The arms were reused from a legacy dataset, whose lines are not verifiable from a manifest (`results/benchmark/xlib_summary.json:arms_dataset`). The discovery_only sensitivity (`u_rule`) was cut under the pre-declared cut order (DEV-CUT-XLIB-DISC, `results/deviations.csv`). It had also crashed: `eligible_genes.csv` lacks columns that the discovery step expects. The main all-D rule is unaffected. Caveat: eligibility was selected on all-D effects, so absolute rates are optimistic (DESIGN §6.3).

### 7.7 P1 and P2 replication in R
A single run after `freeze-catalogue`.
The run read R once: one replication read of 255 lines and one set-A read of 15 lines, both at commit ae6aa84 (DESIGN.md §1.3); the read is logged in `results/replication_access_log.csv`. The pool was frozen before the tag: 1,020 powered U_G candidates, so k_G = 100 (`results/benchmark/p1_pool.json`). B = 1,000 paired, patient-grouped bootstrap draws (`results/replication/p1_summary.json`).

Rates and contrasts (`results/replication/p1_summary.json`; `results/replication/p1_by_access.csv`):

| Stratum (k) | S1-obs | S1-pred | EN-obs | EN-pred | EN-pred − EN-obs (95% CI) | EN-obs − S1-obs (95% CI) |
|---|---|---|---|---|---|---|
| G, headline (100) | 0.99 | 0.93 | 0.89 | 0.82 | −0.07 (−0.12, 0.02): no detectable difference | −0.10 (−0.17, −0.05): univariate better |
| direct (22) | 0.955 | 0.727 | 0.773 | 0.591 | −0.18 (−0.32, 0.00): no detectable difference | −0.18 (−0.23, 0.09): no detectable difference |
| non-direct (100) | 0.99 | 0.97 | 0.86 | 0.88 | 0.02 (−0.08, 0.06): no detectable difference | −0.13 (−0.20, −0.07): univariate better |
| proxy (94) | 0.904 | 0.819 | 0.787 | 0.84 | 0.05 (−0.05, 0.10): no detectable difference | −0.12 (−0.14, 0.00): no detectable difference |

(`results/replication/p1_summary.json`, `p1_by_access.csv`.)

**Verdicts.**
- "Predictor adds value": not supported in any stratum.
- "Multivariate adds value": the reverse holds in the headline and non-direct strata, where the univariate test is better.
- The high absolute rates (`results/replication/p1_summary.json`: 0.82–0.99) reflect lists of powered, strong D associations re-tested in lines that share the Chronos processing with D (DESIGN §1.2). They measure generalisation across lines, not across processing.
- The rates also track list composition:
  - Lists ranked by p pick the strongest effects, so the contrasts partly compare ranking rules.
  - The S1-obs top 100 is dominated by TP53-status features and spans few genes. EN-pred covers more genes.
  - Part 1 had already scored its predictor on R (DESIGN §1.2), so R is held out from this project but not untouched.
- **Bootstrap bias.** The bootstrap mean sits below the observed rates in every list (`results/replication/p1_secondary.json`, bootstrap_diagnostic). The percentile CIs are therefore unreliable, and the basic CIs are reported alongside as a diagnostic.
- **Secondary k** (`results/replication/p1_secondary.json`):

| k | S1-obs | S1-pred | EN-obs | EN-pred |
|---|---|---|---|---|
| 50 | 0.98 | 0.98 | 0.92 | 0.90 |
| 200 | 0.91 | 0.845 | 0.715 | 0.79 |
| 500 | 0.702 | 0.708 | 0.60 | 0.634 |

- Not computed: the remaining P1 secondaries (gene-matched, observed/expected, sign concordance, Deming) were cut. The R non-Avana subset would have needed a second R read (DEV-P1-SECONDARIES, `results/deviations.csv`).

**Set A** (`results/report_counts.json`: set_a, counted from `results/replication/setA_lines.csv`; 15 non-cancerous lines). Across the 2,142 requested genes, the median dependent fraction is 0.20, and 144 genes are dependent in at least half of the set-A lines (`results/report_counts.json`: set_a).

### 7.8 Limits classes
Final state, after R (`results/limits/limits_summary.json`, R_status available), with R power for every R test (`results/replication/r_power.csv`) and the data ceiling for class 8 (Section 7.4). The Tier B attribution is missing.
- Genes (3,501): 3b 1,415; 3 884; 2 417; 9 364; 5 282; 7 70; 6 46; 1 21; 8 2 (`results/limits/limits_summary.json`; codes in `results/limits/gene_classes.csv`).
- MD-global pairs: class 9 = 552; class 5 = 537; class 6 = 160; class 7 = 131 (`results/limits/limits_summary.json`).
- MD-strict pairs: class 9 = 615; class 5 = 10,623; class 6 = 174; class 7 = 148; class 4 = 102 (`results/limits/limits_summary.json`).
- Class 8 from the ceiling: 7 gene-level rows on 4 genes; ESR1 and PKMYT1 end in gene class 8 (`results/report_counts.json`: ceiling_class8).
- **CN/cutting-artefact suspects (NC3 outcome action).** 46 genes with at least one dependent self-null pair carry the label `cn_cutting_suspect`; the label changes no class. 33 MD-global pairs and 13 of the 552 shortlist pairs sit on these genes, including DDX3X at shortlist ranks 2, 4 and 14 (`results/limits/limits_summary.json`: cn_cutting_suspect; `results/candidates/shortlist.csv`).
- Labels: 231 processing-limited pairs and 1,761 predictor-blind rows (`results/limits/limits_summary.json`).
- NC6 fails (Section 7.5).

### 7.9 Candidate cards
`known_recovery.csv` has 88 rows (`results/candidates/cards_summary.json`; Section 7.4). With R: MD-global recovers 22 and MD-strict 28 (`results/candidates/cards_summary.json`).

`shortlist.csv` (552 rows, top 15 below) holds class-9 MD-global pairs that replicated in R (`results/candidates/cards_summary.json`). Its features are PROG 171, EXPPAR 162, COMP 61, LOFBI 37, SIL-other 26, CLIN 23, HOT 15, PSCORE 14 (`results/report_counts.json`: shortlist, counted from `results/candidates/shortlist.csv`).

Top rows (`results/candidates/shortlist.csv`):

| Rank | Dependency | Determinant | Access | Powered lineages | R one-sided p | Combinatorial | Set-A dependent fraction |
|---|---|---|---|---|---|---|---|
| 1 | FAM50A | low FAM50B expression (EXPPAR:26240) | direct | 21 | 1.6e-26 | replicated | 0.47 |
| 2 | DDX3X (cn_cutting_suspect) | low DDX3Y expression (EXPPAR:8653) | direct | 16 | 1.7e-21 | tested, not replicated | 0.07 |
| 3 | DNAJC19 | low DNAJC15 expression (EXPPAR:29103) | direct | 21 | 1.2e-19 | tested, not replicated | 0.53 |
| 4 | DDX3X (cn_cutting_suspect) | expression programme NMF49 | expression-derived | 16 | 2.8e-19 | — | 0.07 |
| 5 | MYBL2 | EXPPAR:4603 | direct | 14 | 1.8e-17 | tested, not replicated | 0.33 |
| 6 | MDM2 | TP53 mutation composite | proxy | 19 | 2.9e-17 | — | 0.67 |
| 7 | EIF1AX | EXPPAR:9086 | direct | 17 | 1.3e-16 | tested, not replicated | 0.33 |
| 8 | TP53 | TP53 mutation composite | proxy | 19 | 4.1e-16 | — | 0.00 |
| 9 | MDM2 | TP53 loss of function (LOFBI:7157) | potentially-direct | 20 | 3.4e-15 | tested, not replicated | 0.67 |
| 10 | RPP25L | EXPPAR:54913 | direct | 20 | 7.5e-15 | tested, not replicated | 0.47 |
| 11 | INTS6 | EXPPAR:203522 | direct | 19 | 9.1e-15 | tested, not replicated | 0.20 |
| 12 | RAB6A | EXPPAR:51560 | direct | 17 | 6.0e-14 | tested, not replicated | 0.73 |
| 13 | ATP6V0E1 | EXPPAR:155066 | direct | 13 | 1.0e-13 | tested, not replicated | 0.07 |
| 14 | DDX3X (cn_cutting_suspect) | Y-chromosome loss (YABS:8653) | expression-derived | 16 | 3.1e-13 | tested, not replicated | 0.07 |
| 15 | EFR3A | EXPPAR:22979 | direct | 10 | 1.1e-12 | tested, not replicated | 0.47 |

(`results/candidates/shortlist.csv`.)

**Case studies** (chosen by `build_candidate_cards.py`):
1. **FAM50A, dependent when the paralog FAM50B is low.** Cross-lineage across 21 lineages, replicated in R and in the combinatorial screens (`results/candidates/shortlist.csv`). But set-A dependency is 0.47 (`results/candidates/shortlist.csv`).
2. **DDX3X, dependent when DDX3Y is low or the Y chromosome is lost.** This is the same biology captured three times: paralog expression, an NMF programme, and Y loss. The determinant is consistent with sex and Y-chromosome status of the line rather than a lineage-specific tumour alteration; sex was not tested directly. The pair was tested in the combinatorial screens and did not replicate there (`results/candidates/shortlist.csv`). Set-A dependency is low (`results/candidates/shortlist.csv`: 0.07). DDX3X carries the `cn_cutting_suspect` label from the NC3 outcome action (Section 7.8; `results/candidates/shortlist.csv`).
3. **DNAJC19, dependent when DNAJC15 is low.** Whether low DNAJC15 reflects promoter silencing was not tested here, because the determinant is the expression feature. It was tested in combinatorial screens but not replicated there, and set-A dependency is 0.53 (`results/candidates/shortlist.csv`).

None of these is named "novel". All three are paralog-buffering dependencies, the class the predictor's M2 stage encodes directly. The word "novel" is not applied to any gene.

## 8. Deviations and amendments

### 8.1 Deviations (`results/deviations.csv`)
| Id | Item | Provisional → final | Reason |
|---|---|---|---|
| DEV-TOY-NC3 and DEV-TOY-NC3-04aa982 | NC3 bound | 0.1 → 0.1 | no faithful toy analogue; per-replicate denominator 3, so 1/n > 0.1 |
| DEV-TOY-NC5 and DEV-TOY-NC5-04aa982 | NC5 bound | 0.2 → 0.2 | no toy analogue; denominators 1–2, so 1/n > 0.2 |
| DEV-TOY-NC6-04aa982 | NC6 bound | 0.7 → 0.7 | toy value 1.0 in 5/5 replicates: zero variance at the edge of the range (rule R2) |
| DEV-DF-FALLBACK | S1 df rule | ν rule → binary df = carrier clusters − 1 | pre-declared fallback, triggered because the calibration test failed at n = 255 |
| DEV-CAL-CONSERVATIVE | calibration size floor 0.035 | kept; 4/12 cells conservative | no further fallback is pre-declared; test left failing; lower power for features with 3–5 carriers |
| DEV-DICT-1 | fusion 5′-driver exception | — → KMT2A added | D10, A5 |
| DEV-CGT-1 | ALK carrier cap | 10 → 11 carriers kept | pre-data estimate (D9) |
| DEV-BUDGET-M3B-M4 | tier timeboxes | — | M3b builder walls not recorded (D13-viii); M4 walls over their timeboxes on a 15.7 GB machine (D17); the overnight interruption |

### 8.2 Amendments after freeze (DESIGN.md)
| Id | Content | Did result data exist? |
|---|---|---|
| A1 (D9) | CGT operationalisation: ER+ = breast ∧ luminal ∧ ESR1 ≥ 3.0; KMT2A-r; EGFR_activating excludes T790M; FGFR targets; indication fixed effect | no (before `freeze-cgt`) |
| A2 | FDA De Novo / 513(f)(2) grants allowed as source records | no |
| A3 | The refuted set has 3 rows by construction | no |
| A4 | Tier `preclinical` reported as secondary, outside the P2 headline | no |
| A5 (D10) | KMT2A added to the 5′ fusion exception; SPEC_MMR excluded after failing its sanity check (ρ 0.27 vs MSIScore, threshold 0.5) (RESEARCH_TREE.md D10) | no |
| A6 (D11) | Ridge predictor: prediction-arm FDP is not detectable by the null, so it is reported as uncalibrated; eligibility strata with < 5 lines excluded; pure-null count rule; toy file names | no (toy only) |
| A7 (D12–D15) | Catalogue = MD-global; NC bound rulings; fold-offset and denominator findings reported; df fallback triggered; ridge-replicate correction | no (toy runs of the real engine; no real association result) |

## 9. Limitations

- **The gate failed twice (21/39).** The pipeline recovers 54% of the specified positive controls on MD-strict, and 17/39 on MD-global (`results/gates/G_pc_attempt2.json`). Families ADD-CN, SIG and W fall below 50% (`results/gates/G_pc_attempt2.json`; `configs/design_params.json`: gate_s9). The gate is a pre-declared bar and was not re-scored after the fact.
- **The gate denominator mixes universes.** The 7 U_ctrl controls are 0/7, and ESR1 could not be recovered by construction (own-expression is family F2) (`results/gates/G_pc_attempt2.json`: by_universe). The eligible-only recall (21/32) does not replace the verdict (`results/gates/G_pc_attempt2.json`: by_universe).
- **Negative-control bounds.** NC2b = 1.2 is very loose (provisional value 0.1) (`configs/frozen_hashes.json`: bounds_values). NC3, NC5 and NC6 keep their provisional bounds, because the toy has no faithful analogue (ER17–ER21).
- **NC1 failed for EN-obs** (`results/null/summary.json`; FDR_emp 0.20). EN-pred (0.327), S1-pred and MD-strict are "not calibrated" (`results/null/summary.json`). Only S1-obs (0.022) and MD-global (0.0009) are within 0.10 with no label (`results/null/summary.json`; `configs/design_params.json`: NC1). The null covers 300 genes with B = 20 (`results/null/summary.json`).
- **NC6 failed.** Lineage-TF genes are not class 3 under the frozen rule (0/18) (`results/limits/limits_summary.json`), so class-3 (lineage) labels are flagged unreliable.
- **P2 is descriptive.** The gate-disjoint headline has 6 testable tests, and every arm recovers 0 (`results/clinical/q4_summary.json`, non_proxy|gate_disjoint|k200). EN-pred recovers no CGT test at k = 200 in any non-proxy subset.
- **65% of MD-global pairs are lineage-unresolvable** (`results/lineage/summary.json`; `results/report_counts.json`: lineage; 894/1,380), and 139 of the 430 cross-lineage labels are single-lineage driven (`results/lineage/summary.json`).
- **Completed after the R read** (DEV-POST-R-COMPLETION, `results/deviations.csv`), with definitions frozen beforehand and D data only: the data ceiling (19 of 28 ceiling-testable tests data-supported, 9 data-null) (`results/clinical/ceiling_summary.json`); NC2/NC3/NC5 on real data, of which NC2a and NC3 fail and NC5 passes (`results/controls/nc_real.json`); refuted-vs-approved recall, refuted 0/3 in every arm (`results/clinical/cgt_secondary.json`); and the independent delta check, 0 flag disagreements over 80 block rows (`results/gates/G_indep.json`: task4_engine_comparison).
- **NC2a failed** (0.104 against ≤ 0.0011; outcome action: "NC2a failed", then a bug check, which found none) (`results/controls/nc_real.json`).
- **NC3 failed** (0.424 against ≤ 0.1; outcome action: "NC3 failed", affected genes labelled CN/cutting-artefact suspects): 46 genes, 33 MD-global pairs and 13 of the 552 shortlist pairs, including DDX3X at ranks 2, 4 and 14 (`results/controls/nc_real.json`; `results/limits/limits_summary.json`: cn_cutting_suspect; `results/candidates/shortlist.csv`).
- **Confounder audit, library block:** the design is ambiguous on whether screen days belong to the library block (DEV-DELTA-LIBRARY-BLOCK, `results/deviations.csv`). The independent check agrees with the engine on every flag.
- **Three toy checks fail** (Section 5.2), one of them the positive control for the fold artefact. The ridge-predictor FDP of the prediction-based lists is invisible to the permutation null.
- **The β tolerance failed in the independent check** (Section 5.6), on one of 30 pairs, by 25% over the tolerance (`results/report_counts.json`: independent, from `results/gates/G_indep.json`). The cause is attributed, not removed.
- **Calibration.** CR2 + ν was anti-conservative at n = 255 (`results/deviations.csv`: DEV-DF-FALLBACK). The fallback is conservative for binary features with 3–5 carriers, so rare determinants are under-powered (`results/deviations.csv`: DEV-CAL-CONSERVATIVE; ER11). Example: EGFR-activating, 8 carriers (RESEARCH_TREE.md D19).
- **The catalogue is a D-only, selection-conditional list.** MD-global p-values are not confirmatory. With a BH family of about 6.6 million tests (`results/report_counts.json`: bh_family), borderline known pairs fail: SMARCA2, ARID1B, MAT2A.
- **Methylation** is assayed for only 456 of the 1,022 D lines (`results/limits/limits_summary.json`; `results/gates/G_indep.json`: notes; DESIGN.md §1.3; 456 in the RRBS matrix, 453 with a METH feature value), and methylation findings carry that coverage.
- **Proxy and access class.** Only non-direct pairs can support "the predictor learnt more than its hand-coded inputs". Most gate controls are not direct (`results/gates/G_pc_attempt2.json`: by_access_class; direct 4/10).
- **Agent-chosen rules.** The cross-lineage rule can mislabel a one-lineage effect (D13-vi). The combinatorial hit rule is an agent choice (D18).
- **Shared processing.** R shares the Chronos fit with D, so replication tests generalisation across lines, not across processing (DESIGN §1.2).
- **External data.** Flister 2024 is a preprint and analysis-only. The MOAlmanac reuse licence is unconfirmed (DATA.md).
- **Not done.** LMM, the M3 refit, clustering and the CTG harvest were cut (DESIGN §14). The Tier B M1-proxy attribution is not yet coded (`docs/process_checks.md` §4).

## 10. Reproducibility

- **Tags.** `freeze-design` (4a079e8, DESIGN v3.2), `freeze-cgt` (681bf39), `freeze-bounds` (4235153). `freeze-catalogue` (ae6aa84; the R read was logged at this commit).
- **Key commits.**
  - engine: 04aa982 (toy suite);
  - arms: f673762, with code at 8220ab1 (`results/benchmark/summary.json:commit`);
  - gate attempts: 0346d69 and 415c593 (the gate JSONs record f6737621 and b4874674);
  - independent check: 4c902fb;
  - subtypes and combinatorial: dcbacd4;
  - confounder-audit change: 19e4b31;
  - null draws: a7c6c8a (run started 08:49; resummarised at 732e50e);
  - M5 audits: f42ddf4;
  - pre-R limits: b125298;
  - XLIB: 2d3675a and 8989602;
  - toy suite on the final engine: dfcd474.
- **Frozen configs and hashes.** `configs/frozen_hashes.json` has sections `design`, `bounds`, `bounds_values` and `actions`; they are checked by `scripts/check_frozen.py`. `configs/outcome_actions.csv` was hashed before any real output.
- **Seeds and fixed inputs.** `configs/seeds.json`; subsamples in `configs/subsamples.json`; null genes in `configs/null_genes.txt`.
- **Pipeline scripts (`scripts/`):**
  - inputs: `verify_inputs`, `data_status`;
  - builders: `build_*` (dictionary, covariates, eligibility, controls, CGT, null genes, subsamples, U_ctrl genes);
  - arms and gate: `run_arms`, `run_gate`, `run_null`, `run_xlib`;
  - audits and summaries: `run_component_readout`, `run_confounder_audit`, `run_subtypes`, `run_cross_lineage`, `run_m5_audits`, `run_q4`, `run_combinatorial_lookup`, `run_limits`;
  - replication and cards: `build_p1_pool`, `replicate`, `build_candidate_cards`.

  Execution order follows RESEARCH_TREE.md. The commit and config hashes of each run are written to `results/registry.csv`. Implementation notes are in `docs/impl_notes_*.md`.
- **Tests.**
  - `tests/`: the fast suite passed at 83e575b.
  - Toy suite: `BIOMARK_TOY_PIPELINE=1 python -X utf8 -m pytest tests/test_toy_answers.py -k pipeline` (`results/toy/suite_run_04aa982.json`, wall time 3,565 s at 04aa982).
  - R is unlocked only by `replicate.py` with `PART2_UNBLIND=1`, after the tag (DESIGN §1.3).
- **Number check.** `scripts/check_report_numbers.py --strict`, final run:
  - Every cited number matches its cited source at the printed precision. Structural numbers are exempt only through `configs/report_number_exemptions.json`, each with a reason, and no number is uncited.
  - The strict gate still fails on **four source checks**: `results/gates/G_indep.json`, `results/cgt/summary.json`, `results/toy/suite_run_04aa982.json` and `results/toy/nc_bounds.json` have no `results/registry.csv` row. Their writers (the independent check, `build_cgt.py` at `freeze-cgt`, and the toy runs) predate the run registry, and rows were not back-filled after the fact.
  - The non-strict check passes, with only those four warnings.
