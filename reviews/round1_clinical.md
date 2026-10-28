# Review round 1: clinical ground truth and translational validity (R1-CL)

Reviewer role: translational oncologist / clinical biomarker expert. Scope: DESIGN.md v1 §7.1, §7.3, §8 and the CGT
draft in docs/lit/clinical_ground_truth.md §2 (the source for `configs/clinical_ground_truth.csv`).
Read-only checks on Part 1 dev lines (1022): gene-level hotspot carriers BRAF 72, KRAS 167, PIK3CA 103, EGFR 13,
ERBB2 11, ESR1 1, KIT 3, FLT3 4, IDH2 0. Raw fusion rows with ALK as 3' partner: 18 lines, FGFR2 as 5' partner: 37,
KMT2A as 5' partner: 20, BCR::ABL1: 11, RET 3': 1. `ModelSubtypeFeatures` gives HER2+ 13, MSI 9, and ER+ only ~5
lines. Breast has 54 dev lines and NSCLC 80. The hotspot matrix is **gene-level only**, so it cannot tell alleles
apart.
None of the changes below adds heavy compute. Together they cost under 1 h of CPU, mostly per-row tests on fewer than
100 CGT tests.

---

### R1-CL1 (critical): clinically wrong or mis-specified expectations in the CGT draft
**Problem.** Several rows would score the method against a wrong clinical truth:
- (a) **KIT D816V → imatinib, S (row 18).** This is wrong. Imatinib is labelled for ASM *without* D816V (or with
  unknown KIT status), and the D816V assay (H140006) selects patients by *excluding* D816V. D816V is a
  drug-resistance allele. On the KO side, D816V mast-cell lines are still KIT-dependent.
- (b) **del(17p)/TP53 loss → venetoclax/BCL2, S (row 51).** This is wrong. del(17p) defined an unmet-need
  population (the first accelerated approval). It does not predict greater BCL2 dependency, and venetoclax is now
  labelled for CLL regardless of 17p.
- (c) **"BRAF hotspot → BRAF" (§7.1, row 1).** Only V600 (class I) predicts response to BRAF inhibitors. Class III
  (kinase-impaired, e.g. D594/G466) mutants are RAS/CRAF-dependent, not BRAF-dependent. A gene-level hotspot flag
  (72 dev lines) pools all of these.
- (d) **"EGFR mut/amp → EGFR" (§7.1).** EGFR amplification is not a clinical predictive biomarker in any indication:
  EGFR TKIs failed in EGFR-amplified GBM, and EGFR amplification is not a TKI selection criterion in NSCLC. Also,
  several of the 13 EGFR "hotspot" dev lines may be passengers in hypermutated lines.
- (e) **PTEN loss → AKT1 (row 27).** Capivasertib is pan-AKT, so a single-gene AKT1 KO is not the clinical
  perturbation. The preclinical dependency of PTEN-null cells is on PIK3CB.
- (f) **Rows 23 and 24** (KRAS/NRAS WT → EGFR, S; mutant → EGFR, W) are the *same contrast* counted twice.
- (g) **ESR1 LBD missense (row 29).** Measured against ER+ ESR1-WT lines, the expected KO sign is N: both groups are
  ER-dependent, and the clinical claim is SERD over AI, not "more ESR1 dependency". There is also only 1 dev carrier.
- (h) **CDKN2A±MTAP (row 41)** duplicates row 39 once 9p21 is locus-collapsed (§5.4).

**Change (add to §7.3 "Mapping").** Add a rule: "Before `freeze-cgt`, the CGT table is corrected as follows:
- (a) Split row 18. KIT D816V → avapritinib/midostaurin, target KIT, S. KIT D816V → imatinib moves to the resistance
  partition.
- (b) Row 51 becomes class `population-only`, sign N.
- (c) BRAF rows require V600, both in the CGT and in the §7.1 positive control. Non-V600 BRAF is not scored as
  sensitivity.
- (d) EGFR amplification is removed from the CGT and stays only as a cell-line positive control, labelled
  'non-clinical'. The EGFR CGT feature is restricted to activating alleles (L858R, exon-19 del, G719X, L861Q,
  exon-20 ins, T790M).
- (e) Row 27 gets target set {AKT1, AKT2, AKT3} under the family rule (R1-CL6).
- (f) Rows 23 and 24 merge into one test: 'RAS-mut → EGFR, W, comparator RAS/BRAF-WT CRC'.
- (g) Row 29 gets sign N vs ER+ WT and is marked untestable (< 3 carriers).
- (h) Row 41 merges into row 39 as the 9p21 locus feature, and the MTAP-specific conditional test is reported (§5.4)."

**Check.** `tests/test_cgt_table.py` asserts that:
- no row has (KIT D816V, imatinib, S);
- no row has (del17p|TP53, BCL2, S);
- every BRAF S row has `feature_def` matching `V600`;
- no CGT row has feature `AMP:EGFR`;
- the AKT row has |target_set| = 3;
- no two rows share (feature, target_set, comparator) (see R1-CL3).

### R1-CL2 (critical): the positive-control gate and the CGT recall are not independent
**Problem.** The §7.1 gate includes ERBB2-AMP, EGFR, BRAF, KRAS, PIK3CA, ESR1, MSI→WRN and MTAP→PRMT5/MAT2A. These
are most of the *testable* concordant CGT pairs. The pipeline must pass the gate, and "stop and diagnose" invites
fixing it until it does. The Q4a recall on those pairs is then partly guaranteed by construction. It would overstate
clinical validation, and the forward-looking WRN/PRMT5 claim (§7.3) would be circular.

**Change (add to §7.3).** "Each CGT test carries `in_gate ∈ {0,1}`. Q4a is reported on (i) all testable tests and
(ii) **gate-disjoint** tests. (ii) is the headline whenever the pipeline changed after the first gate run. Every
post-gate change is logged in `results/gate_iterations.csv` (commit, reason, gate result before/after). Only bug
fixes are allowed there, and thresholds and 𝒟 are not tuned. The forward-looking claim (MSI→WRN, MTAP→PRMT5/MAT2A)
is stated as 'recovered, and also used as a pipeline gate', never as independent validation."

**Check.** The Q4a summary JSON has keys `recall_all` and `recall_gate_disjoint`. `results/gate_iterations.csv`
exists, with one row or more (first run = row 0). If it has more than one row, the report generator emits the
gate-disjoint number first. A grep test confirms REPORT.md does not call WRN/PRMT5 recovery "independent".

### R1-CL3 (critical): define the scoring unit, the comparator, the sign convention and the resistance partition
**Problem.** The CGT is drug-row based. BRAF V600E appears 3 times, PIK3CA 2, PARP/BRCA 3 and MSI 2, so recall is
inflated by duplicate drugs. "Expected sign" has no stated comparator: ER+ vs all lines is a lineage contrast, while
ER+ vs ER− breast is the clinical contrast. It also has no stated mapping to the sign of the Chronos coefficient.
Finally, "resistance pairs must not be assigned a sensitivity direction" (§7.3) contradicts row 24, which is a W
recall target. Bypass resistance is testable. Acquired on-target resistance alleles are not.

**Change (replace the §7.3 "Pre-registration" bullet).**
- The scoring unit is a **CGT test** = (𝒟 feature, target set T, comparator population, indication set). Drugs that
  map to the same test are collapsed, and the drug list is kept as an attribute.
- Every test pre-registers:
  - a comparator: `within-indication` (default for lineage-restricted biomarkers: ER/HER2 within breast and
    gastric, EGFR within NSCLC, RAS within CRC), or `pan-cancer + lineage covariate` (tumour-agnostic biomarkers);
  - a sign: S ≡ the coefficient of the biomarker on the target's Chronos effect is < 0 (more negative in
    biomarker+). W ≡ > 0. For absence biomarkers (WT), code the feature as presence of the mutation and flip the sign;
  - a partition:
    - `sensitivity` (S rows);
    - `bypass-resistance` (W rows, scored: RAS-mut → EGFR in CRC; RB1 LOF → CDK4/CDK6);
    - `acquired-resistance` (ABL1 T315I, EGFR C797S, KIT D816V for imatinib: untestable, listed only);
    - `expected-null` (N rows).
- The success criterion for a `bypass-resistance` row is the W sign. Any S-sign hit on a W row is a **sign error**.

**Check.** The CSV schema is validated in `tests/test_cgt_table.py`:
- required columns `test_id, feature, target_set, comparator, indication_codes, sign, partition, drugs, citation`;
- `sign` ∈ {S, W, N};
- `partition` ∈ the 4 values above (+ `clinically-refuted`, R1-CL7);
- unique (feature, target_set, comparator);
- every `sensitivity` row has sign S and every `bypass-resistance` row has sign W;
- the `freeze-cgt` tag precedes the first commit of any file under `results/clinical/`
  (`git merge-base --is-ancestor`).

### R1-CL4 (major): add a method-independent "data ceiling" test per CGT row, so misses are attributed correctly
**Problem.** A CGT miss by MD could mean three things: the method failed; CRISPR KO does not reproduce the drug
(discordance); or the data cannot show it (too few carriers). The design cannot tell these apart, but Q5 needs
exactly this distinction. Several of the strongest clinical addictions are rare in D (ROS1, RET, FLT3-ITD, NTRK,
KIT) and fall below the 𝒟 minimum of 5. MD and S1 then miss them by construction.

**Change (add as §7.3 "Q4-0 data ceiling").** "For every S/W test with ≥ 3 carriers in D, run one pre-specified test:
the §5.3 confounder model restricted to the comparator, plus a Wilcoxon rank-sum test as a robustness check.
Use one-sided p in the expected direction, with BH over CGT tests only. Each test is labelled:
- `data-supported` (q < 0.1 in the expected direction);
- `data-contradicted` (q < 0.1 in the opposite direction);
- `data-null`;
- `untestable` (< 3 carriers, target not in the library, or feature not resolvable; R1-CL5).

Q4a recall is reported over (i) testable tests (carriers ≥ 5, so representable in 𝒟) and (ii) data-supported tests.
The gap between them is the method's miss rate. `data-null` and `data-contradicted` concordant tests are assigned to
§8 class 8 as empirical KO–drug discordance. This gives the 'limits' answer with evidence."

**Check.** `results/clinical/cgt_ceiling.csv` has one row per S/W test, with `n_carriers`, `beta`, `p_one_sided`,
`q`, `label`. A script asserts that the labels follow from the numbers and that recall over data-supported tests is
≥ recall over testable tests. A toy test with planted rows checks the label logic.

### R1-CL5 (major): operationalise each clinical biomarker as an assayable cell-line feature, frozen before the join
**Problem.** Clinical biomarkers are allele-, fusion- or IHC-defined. 𝒟 is gene-level:
- HOT pools BRAF class I/III, KRAS G12C with other KRAS alleles, and activating EGFR alleles with passengers.
- FUS uses any partner call, so 18 dev lines carry ALK as 3' partner although only about 6–8 DepMap lines are known
  ALK-rearranged. The 37 FGFR2-5' calls are mostly not kinase-retaining oncogenic fusions.
- BRCA1/2 "damaging" includes monoallelic truncations in MSI lines.
- ER+ has no direct label: `ModelSubtypeFeatures` covers only ~5 ER+ dev lines.
- §8 class 4 would demote ER (EXP-own) as "tautological", although ER protein/mRNA *is* the clinical biomarker.

**Change (add `scripts/build_cgt_features.py` to §7.3; output `configs/cgt_features.parquet`, built from molecular
data only, committed before `freeze-cgt`).**
- **Alleles** come from OmicsSomaticMutations: BRAF V600; KRAS G12C, and G12/G13/Q61 separately; EGFR activating set
  (R1-CL1d); PIK3CA (E542K, E545X, H1047X, C420R); ERBB2 TKD; AKT1 E17K; IDH1 R132; KIT exon 11 and D816;
  FGFR3 R248/S249/G370/Y373.
  If that file has not arrived by M4, allele-defined tests are labelled `untestable: allele not resolvable`. They are
  not mapped to the gene-level hotspot flag.
- **Fusions:** canonical partner list from CIViC VariantSummaries 5'/3' fields. The kinase gene must be the 3' partner
  (5' for FGFR2/3), and the call must pass the DepMap filtered fusion file and FFPM ≥ 0.1 (pre-set).
  Known positives must be carriers: NCI-H3122 and NCI-H2228 for ALK, K562 for BCR::ABL1, MV4-11 for KMT2A, and
  MDA-MB-175-VII for NRG1.
- **BRCA1/2 and HR genes:** biallelic = damaging dosage 2, or damaging plus LoH (if signatures arrive), or SIL of
  BRCA1. Hypermutated lines (§5.4 definition) are excluded from the carrier set.
- **ER+/HER2+:** first `ModelSubtypeFeatures`. Otherwise, a 2-component Gaussian mixture on ESR1 log2 TPM within
  Breast ∪ Uterus ∪ Ovary dev lines (cut at posterior 0.5, fitted on expression only). HER2+ = ERBB2 AMP or "HER2+".
  HR+/HER2− = ER+ ∧ ¬HER2+.
- **MTAP-loss** = deep deletion ∨ SIL. **MSI** = OmicsGlobalSignatures, else `ModelSubtypeFeatures` "MSI" ∪ top-TMB
  proxy, labelled as proxy.
- **§8 class 4 exemption:** a gene with a CGT test whose clinical biomarker is expression/protein of the target
  (ESR1, ERBB2) is not demoted to "tautological".

**Check.** A unit test asserts that the named reference lines are carriers, and that the ALK-fusion carrier count is
≤ 10. The `cgt_features.parquet` commit precedes the first read of CRISPR data in any CGT script (guard log
timestamp later than the file's commit time). §8 output has no ESR1 in class 4.

### R1-CL6 (major): fix and extend the KO–drug concordance classes
**Problem.**
- "Allele-specific inhibitors where KO hits both alleles" is listed as *discordant*. That misclassifies KRAS G12C:
  KRAS-mutant lines are KRAS-addicted, so the contrast against KRAS-WT is concordant (S). Only "G12C vs other KRAS
  alleles" is N.
- Paralog-family drugs (pan-AKT, MEK1/2, CDK4/6, FGFR1–4, PARP1/2, pan-RAF) are scored on one gene, which
  under-calls them.
- Clinically discordant contexts with concordant KO have no class. The main example is BRAF V600E CRC: KO is
  expected S, but monotherapy fails because of EGFR feedback, and encorafenib needs cetuximab.
- The lit survey (§3.10) says ESR1 KO may be "weakly selective in culture". The expectation should still be S:
  ER+ breast lines in FBS media are ESR1-dependent. Do not pre-weaken it.

**Change (replace the §7.3 KO-concordance list).** Classes:
- `concordant` (includes allele-specific inhibitors of addicted oncogenes, compared against WT);
- `concordant-family`: the target set has ≥ 2 paralogs. The test is recovered if any member is significant with the
  expected sign after Bonferroni ×|T| within the test. The member list is pre-registered, e.g. CDK4/6 → {CDK4, CDK6},
  with CDK4 expected in breast and CDK6 in haematopoietic lines; MEK → {MAP2K1, MAP2K2};
- `concordant-KO/combination-required`: BRAF V600E in CRC. It is scored for KO direction, and its clinical context is
  noted;
- `discordant-mechanism` (PARP trapping, IDH neomorphs, degraders);
- `population-only` (row 51);
- `non-cell-autonomous` (ICI, ADCs incl. HER2-low T-DXd);
- `untestable`.

Combination regimens are scored on the anchor drug's target only (encorafenib → BRAF; inavolisib → PIK3CA).

**Check.** The `ko_class` vocabulary is validated. Every `concordant-family` row has |target_set| ≥ 2 and a
`scoring_rule` field. No KRAS G12C row has class `discordant-mechanism`.

### R1-CL7 (major): add a "clinically refuted" set as a translational specificity control
**Problem.** The CGT contains only successes. Q5 asks which dependencies *cannot* be used as biomarkers. The most
informative cases are pairs that are real in cell lines but failed in patients. Without them, the design cannot show
whether dependency subtypes discriminate translational success from in-vitro-only biology.

**Change (add §7.3 partition `clinically-refuted`, expected KO sign S, citation per row; verify each citation before
`freeze-cgt` (all from memory) and drop any that cannot be verified).**
- KRAS-mutant NSCLC → MEK1/2 (selumetinib + docetaxel, SELECT-1, Jänne 2017 JAMA, negative);
- CCND1 amp / CDKN2A loss → CDK4/6 as a *selection* biomarker (PALOMA-1 cohort 2, Finn 2015 Lancet Oncol: no added
  benefit over unselected ER+);
- BRAF V600E CRC → BRAF monotherapy (vemurafenib, Kopetz 2015 JCO, ORR ≈ 5%);
- EGFR-amplified GBM → EGFR TKI (EORTC 26034, van den Bent 2009 JCO).

Report the MD and S1 recovery rate in the refuted set next to the approved set. If the refuted rate is similar to
the approved rate, that is a Q5 finding: in-vitro dependency does not encode translational success
(feedback, heterogeneity, therapeutic index).

**Check.** `partition == clinically-refuted` rows ≥ 3, each with a PMID. `results/clinical/q4_summary.json` contains
`recall_refuted` and `recall_approved`, with exact binomial CIs.

### R1-CL8 (minor): extend the forward-looking (trial-stage, CRISPR-discovered) set with KO-concordant analogues
**Problem.** The forward-looking set has only 2 pairs, both in the gate (R1-CL2). PARP/BRCA is N, but several
synthetic-lethal targets in BRCA/ATM/CCNE1/SMARCA4 contexts were found in functional-genomic screens, are KO-concordant
and are in biomarker-selected trials.

**Change.** Add trial-stage rows only if the ClinicalTrials.gov v2 eligibility regex confirms them before
`freeze-cgt`. Each is from memory until then:
- CCNE1 amp → PKMYT1 (lunresertib, NCT04855656);
- ATM LOF → ATR (camonsertib, NCT04497116);
- SMARCA4 LOF → SMARCA2 (PRT3789, NCT05639751; this one is also in the gate, flag `in_gate`);
- BRCA1/2 biallelic → POLQ (ART4215, NCT04991480) and → USP1 (KSQ-4279, NCT05240898).

**Check.** Each such row has an NCT id whose cached CTG JSON (`data/clinical/ctg/<NCT>.json`, checksummed) matches
the biomarker regex. Rows without a matching cached record fail the table test.

### R1-CL9 (major): concrete definition of "correct patient population" (Q4c)
**Problem.**
- "Lineage composition of the subtype vs approved indication (Fisher)" is confounded by biomarker prevalence. HER2+
  lines are breast/gastric because ERBB2 amplification is, not because the subtype learnt the indication.
- OncotreeLineage is too coarse for most labels: Lung mixes NSCLC and SCLC/neuroendocrine, Myeloid mixes CML and
  AML, and Breast carries no HR/HER2 status.
- Indication lists are incomplete. HER2-amp also includes CRC (tucatinib + trastuzumab) and biliary (zanidatamab).
- Tumour-agnostic labels carry exclusions: the dabrafenib + trametinib agnostic label excludes CRC.

**Change (replace the Q4c bullet).** Indications are mapped in `configs/indication_map.csv` (CGT indication →
OncotreeCode / OncotreePrimaryDisease / CGT subtype feature, plus TCGA study IDs), frozen with `freeze-cgt`. For each
data-supported or recovered test:
- (i) Within-indication effect β_in: carriers vs non-carriers among indicated lines.
- (ii) Per non-indicated lineage with ≥ 3 carriers and ≥ 3 non-carriers: β_L.
- (iii) Class:
  - `on-label`: β_in expected sign, one-sided p < 0.05;
  - `beyond-label`: on-label, plus ≥ 1 non-indicated lineage significant. Each such lineage is checked against the
    refuted set and against agnostic labels and exclusions;
  - `off-label-only`;
  - `indication-untestable` (< 3 carriers in the indication).
- (iv) **Population precision** = on-label / (on-label + off-label-only).
- (v) Among carriers only: the odds ratio of accurate-dependent membership for indicated vs non-indicated lineage
  (Fisher). This replaces the subtype-vs-all-lines Fisher test.
- (vi) Rows flagged `tumour_agnostic` must be cross-lineage (§5.2). The `excluded_indications` field (e.g. CRC for
  BRAF V600E) is reported, not scored.

**Check.** `results/clinical/q4c_population.csv` has the columns above, one row per test. Every test with an
indication has ≥ 1 mapped code in `indication_map.csv`. A script recomputes population precision from the CSV.

### R1-CL10 (minor): define the TCGA prevalence comparison
**Problem.** "Cell-line prevalence vs TCGA prevalence" has no metric, denominator or handling of acquired biomarkers.
EGFR T790M, ESR1 LBD and ABL1 T315I are absent from primary TCGA by design.

**Change.**
- Per (biomarker, indication), compare prevalence in D (carriers / indicated D lines) with TCGA PanCancer Atlas
  prevalence from cBioPortal. Denominators are the `_sequenced` sample list for mutations and the `_cna` list
  (GISTIC −2/+2) for CN. ER/HER2 come from the BRCA `SUBTYPE` clinical attribute. NSCLC = LUAD + LUSC.
- Report the Wilson 95% CIs and log2(ratio). Flag |log2 ratio| > 1 as "cell-line over/under-representation", and
  report the Spearman correlation across cells with ≥ 20 D lines.
- Acquired biomarkers get `tcga_na: acquired`.
- Cache API responses under `data/clinical/cbioportal/`, checksummed. Query only the genes in the CGT, about 40
  calls.

**Check.** `results/clinical/prevalence.csv` has columns `dep_prev, dep_lo, dep_hi, tcga_prev, tcga_lo, tcga_hi,
log2_ratio, flag`. A test recomputes one cell (BRAF V600 in SKCM) from the cached JSON.

### R1-CL11 (minor): make the evidence tiers clinically meaningful and the trend test honest
**Problem.**
- "FDA-CDx-backed" is not a clinical evidence tier. ER and CDK4/6 (standard of care, highest evidence) have no CDx.
- The tiers as written overlap.
- With fewer than about 10 testable tests per tier, "recall rises with tier" cannot be tested.

**Change.** Tiers are hierarchical and mutually exclusive (each test takes its highest tier):
- T1: regulatory approval with a label-specified biomarker in that indication (FDA/EMA via MOAlmanac indications;
  CDx presence kept as an attribute);
- T2: CIViC A/B Supports, or an NCCN-cited assertion, without approval;
- T3: trial-stage, with a biomarker-selected NCT only;
- refuted (R1-CL7).

Recall is reported per tier with an exact binomial 95% CI. A Cochran–Armitage trend test is run only if every tier
has ≥ 5 testable tests; otherwise the result is "descriptive only".

**Check.** Each test has exactly one `tier`. `q4_summary.json` has per-tier n, k and CI. The trend p-value is null
when any tier has n < 5.

### R1-CL12 (minor): concrete actionability (Q4b) and direction-of-effect metrics with baselines
**Problem.** "Fraction of MD subtypes whose determinant is a CGT biomarker" is inflated by frequent CGT features
(KRAS, PIK3CA, TP53, BRAF), which can be linked to non-target genes. "Direction-of-effect accuracy" has no
denominator.

**Change (definitions to add to Q4a/Q4b).**
- **Direction accuracy** = (# S/W tests hit by the method with the expected sign) / (# S/W tests hit with any sign).
  **Sign-error rate** = 1 − direction accuracy. Both are reported for MD, S1 and the ceiling test.
- **N-partition hit rate** = fraction of N tests with an MD hit, compared with the MD hit rate over 1,000 random
  (g, d) pairs matched for carrier count and lineage spread. A ratio > 2 is flagged.
- **Pair-level actionability** = |MD candidates where g ∈ T and d = feature of a CGT sensitivity or bypass test| /
  |MD|. **Feature-level actionability** = |MD candidates whose d is any CGT feature| / |MD|. Both are compared with
  S1 and with a null that permutes g labels within MD (200×, preserving feature frequencies).
- **Line coverage** = fraction of D lines in ≥ 1 pair-level actionable subtype, by lineage.

**Check.** `q4_summary.json` contains `direction_accuracy_{MD,S1,ceiling}`, `n_partition_hit_ratio`,
`actionability_pair`, `actionability_feature`, `actionability_null_p95` and `line_coverage_by_lineage`.
A toy-data test with one planted sign-flipped CGT pair confirms that it counts as a sign error, not a miss.
