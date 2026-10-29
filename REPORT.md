# Predictor-derived dependency subtypes and their molecular biomarkers (DRAFT)

**Status: draft.** Written on 2026-10-29, after gate attempt 2.

This is Part 2 of a two-part project. Part 1 built a CRISPR gene-essentiality predictor on DepMap 26Q3 (`../project`, read-only).

- `{{...}}` placeholders mark numbers that need the replication run in R, or an analysis that is still running. No number behind a placeholder has been seen.
- The rules in force are listed in `docs/effective_rules.md`, cited here as ER#. Each number cites its source file.
- No claim in this draft is a replicated or "valid biomarker" claim: R has not been read.

## 1. Summary

1. **Question.** Can the reasons behind an accurate essentiality prediction be read out as determinants that can be assayed? And do those determinants group cell lines that share a dependency for the same reason, within and across lineages? (PROMPT.md; DESIGN §0)

2. **Catalogue.** MD-global is the set of EN-pred selected pairs that are also S1-obs full-family BH hits with the same sign.
   - It holds 1,380 gene-determinant pairs on the 1,022 dev lines (`results/benchmark/summary.json`). MD-strict, the frozen secondary, holds 11,662.
   - The switch to MD-global was made before any real result existed, because MD-strict failed its toy error bounds (A7; Section 4.3).
   - The catalogue pairs form 558 modules, 112 of them with more than one gene (`results/subtypes/summary.json`).

3. **Benchmark against single-feature tests (Q2).**
   - Selected pairs: S1-obs 3,142; EN-obs 738; EN-pred 50,728; S1-pred 276,557 (`results/benchmark/summary.json`).
   - Whether the predictor adds replicated candidates beyond observed data, and whether the elastic net adds beyond S1: P1 contrasts {{P1_EN-pred_minus_EN-obs}} and {{P1_EN-obs_minus_S1-obs}}, at k_G = {{P1_k_G}}.

4. **Combinatorial screens.**
   - 31 of 183 tested MD-global paralog candidates replicate: 17% (CI 11.8–23.2%), against 1.0% of the screened background. OR 19.6 (12.3–30.6) (`results/combinatorial/validation.json`).
   - S1-obs candidates do as well: MD-global vs S1-obs OR 1.02 (CI 0.57–1.81), same file.
   - So the enrichment shows the candidates are real paralog pairs. It does not show that the predictor adds anything to S1-obs.

5. **Clinical ground truth (Q4).**
   - The positive-control gate **failed twice: 21/39 (54%) at attempt 1, and 21/39 at attempt 2** (`results/gates/G_pc_attempt1.json`, `G_pc_attempt2.json`).
   - Still pending: CGT recall at k = 200, direction, actionability and population. P2 recall {{P2_recall_EN-pred}} vs {{P2_recall_S1-obs}}; McNemar {{P2_p}}.
   - The gate-disjoint headline set has 8 scored tests (`results/cgt/summary.json`), so P2 is expected to be "descriptive (underpowered: 8)" (ER27).

6. **Confounders (Q3).** Drop-one-block deltas and the co-location, cutting-artefact, culture and hypermutation audits: {{confounder_audit_summary}}. Fixed-threshold null FDR_emp: {{null_FDR_emp_per_list}}.

7. **Cross-lineage (Q1b).** Cross-lineage vs lineage-restricted MD-global determinants: {{n_cross_lineage}} / {{n_lineage_restricted}}. Replication odds ratio: {{Q1b_OR}}.

8. **Replication in R (P1, 255 lines).** Same-sign replication at one-sided p < 0.05 of each arm's top k_G candidates: {{P1_rates_by_arm}}. By access stratum: {{P1_by_access}}. XLIB: {{XLIB_rates}}.

9. **Limits (Q5).** Gene and pair failure classes: {{limits_class_counts}}. Limits already visible:
   - the gate failed;
   - every non-eligible control was missed (0/7);
   - the df fallback leaves few-carrier binary features under-powered (ER11);
   - the own-expression control ESR1 could not be recovered under the list definitions (Section 6.2).

10. **Candidates.** The shortlist (class 9, replicated in R) and three case studies: {{shortlist_top15}}, {{case_studies}}.

## 2. Question and approach

Part 1's predictor, M5, is the average of three models: a kernel ridge M1, a paralog stage M2 and a LightGBM M3 (RESEARCH_TREE.md, M0). It predicts DepMap 26Q3 CRISPR effects out of fold. The brief (PROMPT.md) asks that every accurate prediction should imply a determinant. To test this, the project:

- reads determinants out of the out-of-fold predictions with a stability elastic net over an interpretable feature dictionary (EN-pred);
- compares EN-pred in a 2×2 design against the same elastic net on observed effects (EN-obs), and against a univariate confounder-adjusted test on observed (S1-obs) and predicted (S1-pred) effects (DESIGN §5.1, D1);
- names a catalogue, MD-global, and groups lines into subtypes and modules by it (DESIGN §8, A7);
- audits confounders, scores the catalogue against clinical ground truth (CIViC, MOAlmanac, FDA companion diagnostics) and combinatorial CRISPR screens, and assigns each dependency a failure class (DESIGN §9–§12);
- decides the primary endpoints on 255 held-out lines (R), read once after tag `freeze-catalogue` (DESIGN §1.3, §6).

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
- **DepMap extras supplied by the user.** OmicsSomaticMutations, OmicsGlobalSignatures, OmicsProfiles, CRISPRConfounders, OmicsCNSegmentsWGS and CCLE RRBS promoter methylation (2018-10-22), all registered before the first dictionary build (D7; DATA.md, "New"). RRBS covers 453 of the 1,022 D lines (D10), so methylation features are not assayed for the rest.
- **Clinical ground truth.**
  - CIViC nightly 2026-10-28 (CC0).
  - MOAlmanac v.2026-08-02. Its repository states GPL-2.0; the reuse licence for the data is unconfirmed in DATA.md.
  - The FDA companion-diagnostic list (public domain).
  - OncoKB gene lists (not redistributed).
  - In total: 64 CGT tests, 27 testable (`results/cgt/summary.json`).
- **Combinatorial screens.**
  - Esmaeili Anvar 2024, a harmonised table covering Dede, Parrish, Thompson, CHyMErA and Ito; plus Dede 2020 and Thompson 2021.
  - The Flister et al. 2024 bioRxiv supplement. Its licence is "cc_no": analysis only, and no row or value of it is written out (`results/combinatorial/validation.json:_licence_note`).
  - The Parrish processed tables could not be downloaded; counts only (DATA.md).
- **Line roles.** D = the 1,022 Part 1 dev lines; R = 255 held-out lines, never analysed; set A = 15 non-cancerous lines (ER1).
- **Not committed.** Raw and restricted data. Bulk tables stay in `cache/`, with their sha256 in `results/bulk_manifest.csv`.

## 4. Methods

### 4.1 Predictor and dictionary
- **Predictions.** Part 1's out-of-fold M5 values, fold-centred per gene. R uses the all-dev refit (DESIGN §1.1).
- **Eligible genes.** A gene is eligible if it is variable or selective, has tier high or medium, and has `r_within` ≥ 0.10 beyond lineage × library. 2,193 of 18,434 genes qualify: 1,172 high, 1,021 medium (`results/dictionary/eligibility_summary.json`).
- **The dictionary 𝒟 (DESIGN §3).** Shared features are capped at 3,000. It holds these feature families:
  - hotspot, allele and loss-of-function mutations;
  - amplification, deletion and shallow loss;
  - fusions, silencing and methylation;
  - clinical composites;
  - expression dose and paralog expression;
  - expression programmes (50 NMF factors);
  - signatures, and TP53/RB composites.
- **D only.** NMF, scalers and thresholds are fitted on D only (DESIGN §1.7).
- **Access class.** Every (gene, feature) pair has one: `direct`, `potentially-direct`, `expression-derived` or `proxy` (DESIGN §4.4). It limits what "the predictor learnt d" can mean.

### 4.2 Confounder model and the four arms
**Confounder model** (DESIGN §5.3, ER9):
- Never dropped: the design blocks lineage, library group, culture block and growth.
- Droppable per test: MSI, mutation burden, TP53, own copy number and the own-expression gate.
- CR2 standard errors, clustered by patient.

**df rule.** The pre-declared reference-t df is min(n − rank − 2, ν). It was shown to be anti-conservative at n = 255. The pre-declared fallback therefore applies to binary features: df = carrier clusters − 1. It is used in every S1 test, in power and in the R tests (ER11; Section 5.5).

**The four arms** (DESIGN §5.1):
- **S1-obs and S1-pred.** One test per gene-feature pair, on observed or predicted effects. BH runs over the main family F1, 6,559,561 tests (`results/gates/G_indep.json`).
- **EN-obs and EN-pred.** A two-stage hierarchical elastic net: stage G genomic, stage X expression, with Frisch–Waugh–Lovell residualisation on the design block. Selection uses 25 complementary half-sample pairs, frequency ≥ 0.6 and ≥ 90% sign agreement (DESIGN §4.1).
- EN-pred never reads observed effects, only the NaN mask.

All four arms ran on identical universes: 2,193 genes, 26,246,924 rows in 54 parts (`results/benchmark/summary.json`).

### 4.3 MD-global and MD-strict
**MD-strict** is the set of EN-pred selected pairs that pass an observed-effect screen run inside the EN-pred selection (BH q < 0.05, same sign; DESIGN §5.2). It was the frozen catalogue until A7.
- On the engine's toy runs with a ridge-type predictor, its true FDP was 0.228, 0.255, 0.237, 0.200 and 0.160, against a limit of 0.15.
- On pure-null toys, its P(any hit) was 0.30, against a limit of 0.10.
- The cause: the observed screen re-uses the noise that the ridge predictions absorbed (A7-i, A7-v).

**MD-global** is the set of EN-pred selected pairs that are also S1-obs F1 BH hits with the same sign. On the same toys it had FDP 0, 0, 0.013, 0 and 0.025, and P(any hit) 0.00, at equal recall (0.815–0.926 on the ridge toys) (A7-v).

How the two lists are used:
- Both are reported. MD-global is the catalogue for everything downstream; the gate and NC2b stay on MD-strict (ER2, ER3, ER12, ER17).
- MD-strict p-values are selection-conditional, and neither list is called "significant". Confirmation comes from R and from external sources.

### 4.4 Null, negative controls and positive controls
**Fixed-threshold null** (ER23, ER15):
- Line labels are permuted within lineage × library, jointly across all dictionary rows.
- 300 pre-drawn genes; B = 20, raised to 50 only if timing allows.
- A thresholded list passes if FDR_emp_fixed ≤ 0.10.

**Negative controls.** The NC2–NC6 bounds come from five toy replicates (ER16–ER22).

**Positive controls.** 40 were specified, giving 39 gate rows in 10 families (DESIGN §9; D10 for 39/40):
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
- Workers: S1 8; EN 6 requested and 5 effective; null 6 requested.
- The machine has 15.7 GB of RAM.
- Wall times exceed the tier timeboxes. This is logged as lateness, not treated as a design change.

## 5. Verification

### 5.1 Toy simulation with an answer key
- `src/biomark/toysim.py` and `tests/test_toy_answers.py` were written by a verification agent, separate from the pipeline author. Commits carry an `Agent:` trailer (DESIGN §13.1, A6-iv, LESSONS 9).
- The toy contains planted determinants, co-located passengers, culture-confounded and lineage-only genes, and a fold-offset artefact.
- The answer key was frozen before the engine ever ran on it.

### 5.2 Suite on engine commit 04aa982: 19 pass, 3 fail
Source: `results/toy/suite_run_04aa982.json` (3,565 s). **Failed checks:**
1. **Fold-offset artefact positive control.** It needed ≥ 3 of 5 without centring and got 0. The bound was derived from S1-pred BH, but EN selection frequency without centring is 0.36–0.50, below 0.6 (A7-iii, D14). Fold centring stays.
2. **Ridge-toy MD-strict.** FDP_true 0.228 (101 hits, 23 false), against a limit of 0.15. This is the double dip (A7-i).
3. **Pure-null toy, S1-pred.** P(any hit) 0.2, against a limit of 0.1. The earlier engine commit gave 0.25 (D14); the file above is for 04aa982.

MD-global passes its three error checks on the same runs (same file):
- ridge toy: 78 hits, 0 false;
- smooth toy: 60 hits, 0 false;
- pure null: P(any) = 0.0 over 20 toys.

The fixed-threshold null cannot see the ridge-predictor FDP of the prediction-based lists (D11). EN-pred and S1-pred FDR_emp will therefore be labelled "not calibrated for ridge-type predictors", and R replication arbitrates (A6-i).

The suite ran on 04aa982 and `src/` has changed since, so a re-run is planned (RESEARCH_TREE.md, 08:40 note): {{toy_suite_rerun_result}}.

### 5.3 The MD-global switch
- The switch was decided before any real association result existed (A7 header; D15), on the toy numbers in Section 4.3.
- Its real consequence at the gate: MD-global recovers 17/39 controls, against 21/39 for MD-strict (`results/gates/G_pc_attempt2.json:by_list`). The switch did not improve control recovery.

### 5.4 Process checks
Source: `docs/process_checks.md`. The checks are:
- `budget_check`;
- `leakscan`: it caught the planted gene string and the planted constant, and found 0 hits at 04aa982;
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

The fast test suite passed 236 tests at the re-verification commit (git log, 83e575b).

### 5.5 Calibration fallback
The calibration test runs 20,000 null simulations per cell.
- **Under the ν rule:** at n = 255 with doubled carrier SD, size at α 0.05 was 0.072, 0.073 and 0.069 for 3, 5 and 10 carriers. This triggered the pre-declared fallback: for binary d, df = carrier clusters − 1.
- **Under the fallback:** the maximum size is 0.0587 at α 0.05 and 0.0011 at α 1e-3. But 4 of 12 cells are conservative, below the 0.035 floor: 0.03175, 0.01040, 0.02930 and 0.02760 (`results/deviations.csv`: DEV-DF-FALLBACK, DEV-CAL-CONSERVATIVE).

The test is left failing. The cost is lower power for binary features with 3–5 carriers. One consequence is already visible: EGFR-activating has 8 carriers and df 7, so β = −0.63 gives p = 0.13, and the control is missed at the gate (D19).

### 5.6 Independent reimplementation (DESIGN §13.3)
A separate agent that did not read the engine code reimplemented S1-obs, the BH families, the gate and EN selection (`results/gates/G_indep.json`, commit 4c902fb).

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
| Confounding delta, 10 pairs | not scored: `results/confounders/delta.csv` did not exist yet; the checker's own refit flagged 1 of 10 | reported only |

- **The β-tolerance check failed.** The checker attributes it to the subsample fits for gene 5045, where the active set differed (solver tolerance 1e-3 and the screening boundary). The full-D β of that pair agrees to 6e-10 (`G_indep.json:notes`). No tolerance was changed.
- The file records tasks 1 and 2 as pass and task 3 as fail (`overall`).
- **Pair universe.** The checker's universe is larger than the engine's: 6,574,471 vs 6,561,731 pairs, i.e. 12,740 pairs in 75 genes that only the checker has. All of them are untested, because they have fewer than 5 carriers on the gene's rows. The checker carries such pairs as untested rows; the engine leaves them out of its table. The S1 p-values and the BH counts match exactly.
- Independent confounder-delta check against the engine file: {{indep_delta_check}}.

## 6. Results so far (dev set D only; R has not been read)

### 6.1 Benchmark list sizes
The run covered 2,193 genes and 26.2 million (g, d) rows, all four arms in 45 minutes (D19). Selected pairs (`results/benchmark/summary.json`):

| List | Pairs |
|---|---|
| S1-obs | 3,142 |
| EN-obs | 738 |
| EN-pred | 50,728 |
| S1-pred | 276,557 |
| MD-global | 1,380 |
| MD-strict | 11,662 |

The prediction arms select far more than the observed arms. With a ridge-type predictor their error rates are not controlled and cannot be read from the null (A6-i), so list size is not a performance claim. Replication of these lists: P1 {{P1}}.

### 6.2 Positive-control gate: failed at attempt 1 and at attempt 2

**Attempt 1** (`results/gates/G_pc_attempt1.json`): MD-strict recovered **21/39 = 0.538 (95% CI 0.372–0.699)**, below the 60% bar.
- Families below 50%: ADD-CN 1/5, SIG 0/2, W 1/3.
- Families passing: ADD-HOT 5/6, LOF-SL 1/2, P53 4/4, PAR 7/13, PROG 2/2.
- ADD-EXP (0/1) and SIL (0/1) have fewer than two rows each.

**Bug search (D19).** The U_ctrl universe, which holds the non-eligible control and CGT-target genes, had never been built or read. As a result SMARCA2, ARID1B, MYCN and ESR1 had no rows. This was fixed as a DESIGN §9.3 bug, with a toy test that fails before the fix and passes after it (`test_gate_u_ctrl.py`). The 41 U_ctrl genes were run through the arms, and the gate was rerun as attempt 2.

**Attempt 2** (`results/gates/G_pc_attempt2.json`): MD-strict **21/39 again**, with the same three failing families.

| Breakdown | Recall |
|---|---|
| Eligible universe | 21/32 (0.656, CI 0.468–0.814) |
| U_ctrl universe | **0/7** (CI 0–0.410) |
| S1-pred | 29/39 |
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

**All seven U_ctrl misses are genuine; there is no further bug (D20).**
- **SMARCA2 × SMARCA4 loss-of-function:** β = −0.48, p = 1e-5, q = 0.15.
- **ARID1B × ARID1A loss-of-function:** q = 0.099. This and SMARCA2 are lost to the BH family of about 123,000 tests.
- **ENO2, PAPSS1 and ME3 × shallow loss:** β ≈ 0. Shallow loss does not capture the homozygous deletions these pairs need.
- **MYCN amplification:** p = 0.01, absorbed by lineage.
- **ESR1 own-expression:** q = 0.0008, but own-expression belongs to test family F2. The gate counts only F1 lists, so this control could never be recovered under the list definitions. It is reported as a specification gap and not re-scored.

Other misses outside U_ctrl: MAT2A × MTAP deletion (q = 0.15), and EGFR-activating as in Section 5.5 (D19).

Eligible-only recall (21/32) is above 60%, but it is not the gate's denominator, and it was not used to change the verdict. The wording was fixed in advance in `configs/outcome_actions.csv`: "gate failed at attempt 1 (21/39) and at attempt 2 after the U_ctrl fix (21/39); CGT headline restricted to gate-disjoint tests." MD-global's 17/39 is reported for information only (outcome_actions, `gate_md_global`).

### 6.3 Subtypes and modules (MD-global, `results/subtypes/summary.json`)

| | MD-global | MD-strict |
|---|---|---|
| Pairs | 1,380 | 11,662 |
| Pairs with member lines | 1,014 | 8,299 |
| Modules | 558 | 2,221 |
| Multi-gene modules | 112 | 694 |

- The same file also has gene-level totals: 2,134 genes, 279,074 unexplained lines and 109,919 accurate-but-unassigned lines. These totals are not labelled by list in the file.
- Subtype membership uses the share rule (DESIGN §8.1). An "unexplained" line is dependent but mispredicted.
- Cross-lineage labels, co-location and culture checks are pending (Section 7).

### 6.4 Combinatorial validation (`results/combinatorial/validation.json`)
The test covers paralog candidates (loss of A → dependency on B), screened in 23 cell lines across six sources.

| List | Candidates | Tested | Replicated | Rate (95% CI) | Background rate | OR vs background (95% CI) |
|---|---|---|---|---|---|---|
| MD-global | 191 | 183 | 31 | 16.9% (11.8–23.2) | 1.03% (110/10,717) | 19.6 (12.3–30.6) |
| MD-strict | 363 | 317 | 37 | 11.7% (8.4–15.7) | 0.78% (114/14,539) | 16.7 (11.0–24.9) |
| S1-obs | 210 | 198 | 33 | 16.7% (11.8–22.6) | 1.01% (108/10,702) | 19.6 (12.5–30.2) |

- **Between lists.** MD-global vs S1-obs: OR 1.02 (0.57–1.81), Fisher p = 1.0. MD-strict vs S1-obs: OR 0.66 (0.39–1.14), p = 0.11.
- **Matched expectation.** For MD-global the sequence-identity-matched expected rate is 2.5%, so the enrichment over that expectation (16.9% vs 2.5%, a derived ratio of 6.8) is smaller than the OR over the unmatched background suggests.
- The hit rule is an agent choice (Section 4.5).
- **Conclusion.** Catalogue paralog candidates are experimentally supported combinatorial dependencies, but the test shows no advantage over the univariate S1-obs list.

## 7. Pending sections

### 7.1 Confounder audit
Drop-one-block deltas on MD-global ∪ S1-obs hits ∪ controls, including the U_ctrl gate controls.
- Blocks: lineage, library, culture, growth, burden, TP53 and own CN.
- A pair is flagged if abs(Δβ) > 0.5·abs(β_minus) and > SE_minus.
- Also VIF and subtype-confounding.

Result: {{confounder_audit}}. The job is re-running after the overnight interruption (RESEARCH_TREE.md).

### 7.2 Cross-lineage labels (Q1b)
- Labels follow DESIGN §8.3: counts {{cross_lineage_counts}}.
- Fisher OR of R replication and of CGT membership, cross-lineage vs lineage-restricted: {{Q1b_OR}}.
- Known weakness of the rule: a one-lineage effect can be labelled cross-lineage when the null lineages happen to share its sign (1 of 5 toy seeds, D13-vi).

### 7.3 Co-location, cutting artefact, culture, hypermutation
DESIGN §8.4–8.7:
- co-location resolution: {{colocation_resolution}}
- cutting-artefact calls: {{cutting_artefact_calls}}
- culture flags: {{culture_flags}}
- hypermutation/TP53 failures: {{hypermutation_tp53_failures}}

### 7.4 Clinical ground truth: recovery, direction, actionability, population
The headline set is gate-disjoint (ER14, ER26).
- Recall at k = 200 per arm: {{P2_recall_by_arm}}, with exact McNemar {{P2_mcnemar}}, or "descriptive (underpowered: n)": {{p2_status}}.
- Miss decomposition: {{miss_decomposition}}.
- Ceiling labels: {{ceiling_labels}}.
- Direction accuracy and sign-error rate: {{direction_accuracy}}.
- NC5 false-sensitivity rate: {{nc5_value}}.
- Recall on refuted vs approved tests: {{recall_refuted_vs_approved}}.
- Actionability: {{actionability}}.
- Population (on-label, beyond-label, off-label-only, indication-untestable): {{q4c_population}}.
- TCGA prevalence (Tier B): {{tcga_prevalence}}.

Tier counts available now, among the 64 tests: 41 T1, 7 T2, 5 T3, 1 preclinical, 3 refuted (`results/cgt/summary.json`).

### 7.5 Null calibration (FDR_emp)
- FDR_emp_fixed per list: {{null_FDR_emp_per_list}}.
- P(any hit): {{null_p_any}}.
- Selection-adaptive re-run: {{null_FDR_emp_rerun}}.
- B = {{null_B}}.
- NC1 pass/fail: {{NC1_status}}.
- NC2a and NC2b: {{NC2_values}}. Their bounds are 0.0011 and 1.2; 1.2 is very loose (ER17).
- NC3: {{NC3_value}}. NC4: {{NC4_value}}. NC6: {{NC6_value}}.

### 7.6 XLIB (library-independent replication, S1-obs and EN-obs only)
Rates {{XLIB_rates}}; contrast {{XLIB_contrast}}. Caveat: eligibility was selected on all-D effects, so absolute rates are optimistic (DESIGN §6.3).

### 7.7 P1 and P2 replication in R
A single run after `freeze-catalogue`.
- P1 pool size: {{P1_pool}}; k_G: {{P1_k_G}}.
- Replication by arm: {{P1_rates_by_arm}}.
- Contrasts with bootstrap 95% CI: {{P1_contrasts}}.
- By access stratum (`direct`, `non_direct`, `proxy`): {{P1_by_access}}.
- Verdict branch: {{P1_verdict}}.
- Set A normal-line dependency flags: {{setA}}.

### 7.8 Limits classes
Gene classes {{gene_class_counts}}; pair classes {{pair_class_counts}}; predictor-blind and processing-limited labels {{limit_labels}}.

### 7.9 Candidate cards
`known_recovery.csv` {{known_recovery}}; `shortlist.csv` top 15 {{shortlist_top15}}; three case studies {{case_studies}}. The word "novel" is not applied to any gene.

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
| A5 (D10) | KMT2A added to the 5′ fusion exception; SPEC_MMR excluded after failing its sanity check (ρ 0.27 vs MSIScore, threshold 0.5) | no |
| A6 (D11) | Ridge predictor: prediction-arm FDP is not detectable by the null, so it is reported as uncalibrated; eligibility strata with < 5 lines excluded; pure-null count rule; toy file names | no (toy only) |
| A7 (D12–D15) | Catalogue = MD-global; NC bound rulings; fold-offset and denominator findings reported; df fallback triggered; ridge-replicate correction | no (toy runs of the real engine; no real association result) |

## 9. Limitations

- **The gate failed twice (21/39).** The pipeline recovers 54% of the specified positive controls on MD-strict, and 17/39 on MD-global. Families ADD-CN, SIG and W fall below 50%. The gate is a pre-declared bar and was not re-scored after the fact.
- **The gate denominator mixes universes.** The 7 U_ctrl controls are 0/7, and ESR1 could not be recovered by construction (own-expression is family F2). The eligible-only recall (21/32) does not replace the verdict.
- **Negative-control bounds.** NC2b = 1.2 is very loose (provisional value 0.1). NC3, NC5 and NC6 keep their provisional bounds, because the toy has no faithful analogue (ER17–ER21).
- **Three toy checks fail** (Section 5.2), one of them the positive control for the fold artefact. The ridge-predictor FDP of the prediction-based lists is invisible to the permutation null.
- **The β tolerance failed in the independent check** (Section 5.6), on one of 30 pairs, by 25% over the tolerance. The cause is attributed, not removed.
- **Calibration.** CR2 + ν was anti-conservative at n = 255. The fallback is conservative for binary features with 3–5 carriers, so rare determinants are under-powered (ER11). Example: EGFR-activating, 8 carriers.
- **The catalogue is a D-only, selection-conditional list.** MD-global p-values are not confirmatory. With a BH family of about 6.6 million tests, borderline known pairs fail: SMARCA2, ARID1B, MAT2A.
- **Methylation** is assayed for only 453 of the 1,022 D lines (D10), and methylation findings carry that coverage.
- **Proxy and access class.** Only non-direct pairs can support "the predictor learnt more than its hand-coded inputs". Most gate controls are not direct (direct 4/10).
- **Agent-chosen rules.** The cross-lineage rule can mislabel a one-lineage effect (D13-vi). The combinatorial hit rule is an agent choice (D18).
- **Shared processing.** R shares the Chronos fit with D, so replication tests generalisation across lines, not across processing (DESIGN §1.2).
- **External data.** Flister 2024 is a preprint and analysis-only. The MOAlmanac reuse licence is unconfirmed (DATA.md).
- **Not done.** LMM, the M3 refit, clustering and the CTG harvest were cut (DESIGN §14). The Tier B M1-proxy attribution is not yet coded (`docs/process_checks.md` §4).

## 10. Reproducibility

- **Tags.** `freeze-design` (4a079e8, DESIGN v3.2), `freeze-cgt` (681bf39), `freeze-bounds` (4235153). `freeze-catalogue`: {{freeze_catalogue_commit}}.
- **Key commits.**
  - engine: 04aa982 (toy suite);
  - arms: f673762, with code at 8220ab1 (`results/benchmark/summary.json:commit`);
  - gate attempts: 0346d69 and 415c593 (the gate JSONs record f6737621 and b4874674);
  - independent check: 4c902fb;
  - subtypes and combinatorial: dcbacd4;
  - confounder-audit change: 19e4b31.
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
  - `tests/`: the fast suite passed 236 tests at 83e575b.
  - Toy suite: `BIOMARK_TOY_PIPELINE=1 python -X utf8 -m pytest tests/test_toy_answers.py -k pipeline` (about 1 h).
  - R is unlocked only by `replicate.py` with `PART2_UNBLIND=1`, after the tag (DESIGN §1.3).
- **Number check.** `scripts/check_report_numbers.py --strict`: {{number_check_status}}.
