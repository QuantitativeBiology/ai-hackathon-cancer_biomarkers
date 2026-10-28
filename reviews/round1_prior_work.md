# Review round 1: prior work and biology (reviewer R1-PW)

Scope: novelty relative to Tsherniak 2017, Dempster 2020, Behan 2019, Pacini 2024, Gonçalves 2020, Lord 2020, Wainberg 2021 and Pan 2022; biological correctness of the positive controls, biomarker families, paralog and collateral logic, and mechanism attribution; missing biology and artefacts. I checked Part 1 tables read-only (`P1/configs/gene_classes.csv`, `predictable_genes.csv`, `oncogene_pairs.csv`, `biomarker_pairs.csv`). Facts from those tables are quoted below. Literature tagged [from-memory] was not retrieved in this session.

Overall: the design is sound and confounder-aware. Two problems would undermine its central claims, and several biomarker families are defined in ways that make well-known biology unrecoverable by construction.

---

### R1-PW1 (critical): no baseline separates "the predictor" from "a sparse multivariate model"; novelty claims are unstated
**Problem.** Prior work already reads biomarkers out of multivariate dependency models fitted on *observed* effects: Tsherniak 2017 (RF top features), Dempster 2020 (the DepMap "predictive features" ensemble, few features), and Cai 2024 (MOVE, SHAP). §6 compares MD only with univariate S1 and with S1-pred. Any MD advantage could therefore come from fitting a sparse multivariate model, not from the Part 1 predictor, and Q2 cannot then support the project's premise ("an accurate predictor implies a determinant"). §0 also never states what is new relative to Behan 2019, Pacini 2024 and Lord 2020.
**Change.**
- Add the benchmark arm **EN-obs**: §4.1 run unchanged (same 𝒟, unpenalised lineage, α = 0.9, the same stability selection, the same §5.1 observed-effect filter), with observed `y_g` in place of `ŷ_g`.
- In §6, report MD, EN-obs, S1 and S1-pred on the same genes, with matched top-k lists.
- Report the set differences MD\S1, S1\MD and MD\EN-obs with their replication rate in R and their positive-control and CGT content.
- Add to §0 the novelty claims, stated so they can be falsified:
  - (a) the predictor-derived lists replicate better than, or name determinants that are missed by, EN-obs and S1 (Behan/Pacini-style ANOVA);
  - (b) the line-level subtype assignment and its cross-lineage recurrence are new;
  - (c) the scoring of a dependency-derived catalogue against tiered CIViC/FDA ground truth is new (no published precision/recall was found; clinical_ground_truth.md §4).
**Check.** `results/benchmark/lists.csv` has method ∈ {MD, EN-obs, S1, S1-pred} for the same gene set. `results/benchmark/setdiff.csv` exists with columns `set, n, rep_rate, n_pos_ctrl, n_cgt`. A test asserts that EN-obs and MD use identical configs apart from the target.

### R1-PW2 (critical): own-gene tautologies and a collinear covariate block positive controls by construction
**Problem.**
1. §5.3 always includes "own-gene CN of g". For ERBB2 AMP → ERBB2, MYCN AMP → MYCN and EGFR amp, the determinant *is* own-gene CN (or a locus containing g), so the covariate absorbs the effect. Three §7.1 controls then fail by design.
2. 𝒟 includes features that make KO trivially uninformative for g itself: SIL_g, DEL/locus features containing g, and LOF_g (a gene that is not expressed or is deleted cannot be a dependency).
3. EXP-own mixes two things. The first is expression gating: an unexpressed gene is never a dependency, which is trivial and covers most of the EXP-own signal for lineage genes. The second is dosage addiction (more ESR1/MYCN gives a stronger dependency), which is a real biomarker class (Pacini 2024; ER+ is a clinical biomarker). Failure class 4 treats both as tautological unless AMP or SIL backs them, which wrongly demotes ESR1-type addiction.
**Change.**
- In §5.3, drop the own-gene CN covariate whenever d is a CN or locus feature whose span contains g. Use the R1-PW9 test instead.
- Remove SIL_g, DEL/LOSS loci containing g and LOF_g from 𝒟_g. Log them as "self-null" sanity checks: carriers must have effect ≈ 0, which also catches guide and CN artefacts.
- Split EXP-own into EXP-own-gate (expressed vs not, at TPM ≥ 1) and EXP-own-dose (EXP-own refitted within lines with log2(TPM+1) ≥ 1). Only EXP-own-dose can be a determinant. Failure class 4 then reads "determinant is only expression gating".
**Check.** A unit test asserts that for every (g, d) the covariate matrix contains no CN column from a locus containing g when d is CN-derived. It also asserts that 𝒟_g has no SIL_g, LOF_g or locus-containing-g features. `results/md.csv` has an `exp_own_kind` column.

### R1-PW3 (major): "determinant" can name a feature the predictor never saw
**Problem.** Part 1 inputs are restricted.
- M1 uses expression only.
- M2 uses own/paralog expression, CN, deep deletion and mutation.
- M3 uses the top-50 E/C/M/H/meta columns per gene, and meta includes lineage one-hot.
- MSI, ploidy/WGD/CIN, TMB, fusions and most LOF features are not inputs.
The surrogate in 𝒟 will still select them whenever the prediction correlates with them through proxies, for example MSI through MLH1 silencing or WGD through proliferation programmes. Calling such a d "what the predictor learnt" is then a mechanism-attribution error. Prior work on attributing models to biomarkers (Dempster 2020 curator caveat; Tsherniak 2017, where 82% of top features are expression) shows that this proxying is the common case, not the exception.
**Change.** Give every selected (g, d) an `access_class`:
- *direct*: d's raw feature is among M2's 13 features or M3's `sel_full_` for g;
- *expression-derived*: EXP, PROG or SIL;
- *proxy*: d is not an input. Name the top direct-attribution input that carries it (from §4.2).
Claims of the form "the predictor learnt d" are allowed only for direct and expression-derived d. A proxy d is reported as "the predictor learnt a correlate of d (via X)". Report the counts per class.
**Check.** `results/md.csv` has a non-null `access_class`, plus a `carrier_input` column for proxy rows. A test checks that the direct set equals M2 features ∪ M3 `sel_full_` for g in each fold.

### R1-PW4 (major): binary "dependent < −0.5" and failure class 1 discard the flagship clinical pairs
**Problem.** In Part 1 tables, PRMT5 has 98% of dev lines below −0.5, MAT2A 95%, PKMYT1 95% and KIF18A 86%. PRMT5, MAT2A and PKMYT1 are also flagged common-essential. The biomarker for these genes is *graded* (MTAP-null lines are more dependent: Kryukov/Mavrakis 2016; CCNE1-amplified lines are more PKMYT1-dependent: Gallo 2022 [from-memory]). That graded difference is exactly the therapeutic-window logic behind the trials. Under §2 and §8, "accurate-dependent" is nearly everyone, and failure class 1 ("common essential") fires first. The forward-looking MTAP → PRMT5/MAT2A pair would be labelled "no contrast" even if it were recovered.
**Change.**
- Failure class 1 = no between-line contrast (dev SD of effect < 0.15, pre-set), *not* the common-essential flag.
- For genes with frac_below(−0.5) > 0.5, define "dependent" relative to the gene: effect below the gene's dev 25th percentile (pre-declared). §5.1 tests already use continuous effects; state that explicitly.
- The Krill-Burger 2023 caveat (pan-lethal by CRISPR, selective by RNAi) goes into §8 as a named limit.
**Check.** In `results/limits/gene_classes.csv`, no gene with SD ≥ 0.15 is in class 1. PRMT5 and MAT2A carry `dependent_rule = relative`. The gate script reports MTAP → PRMT5 as recovered or not recovered, never as "ineligible".

### R1-PW5 (major): copy-loss biology is limited to deep deletion, which misses hemizygous paralog loss and CYCLOPS
**Problem.**
- VPS4A/VPS4B collateral vulnerability (Neggers 2020) is driven by *hemizygous* loss (18q with SMAD4, 16q with CDH1) and low expression, not by homozygous deletion. Part 1 counted only 4 deep-deleted VPS4B lines, and ENO1 has 1. Under "≥ 5 carriers" these controls silently drop out of the gate denominator.
- CYCLOPS genes (Nijhawan 2012 [from-memory]: PSMC2, SF3B1, POLR2A) are partial-loss vulnerabilities of essential genes. They are invisible with CN < 0.25.
- Chronos's arm-level proximity-bias alignment (part1_brief §7.1) partly removes own-arm CN effects by construction.
**Change.**
- Add a LOSS family: relative CN < 0.7 (hemizygous), locus-collapsed like AMP/DEL. Add continuous CN of the top-3 paralogs to EXP-paralog.
- Re-specify each paralog/collateral control by its biologically correct feature: VPS4A ← VPS4B LOSS or low EXP; VPS4B ← VPS4A LOSS (16q); PAPSS1 ← PAPSS2 DEL/LOSS (10q23 with PTEN); ENO2 ← ENO1 DEL, declared *untestable* (1 carrier).
- Add the CYCLOPS control PSMC2 ← own-locus LOSS (sign: more dependent), with a note that the arm correction biases it toward null. A miss is then reported as a "processing-limited" limit, not as model failure.
**Check.** `configs/dictionary_v1.json` contains LOSS features. `configs/positive_controls.csv` has `feature, sign, n_carriers_dev, testable` for each control. The gate denominator equals the rows with testable == True, and untestable rows are listed in the report.

### R1-PW6 (major): positive controls cover only oncogene addiction and paralog loss; they test none of the signature, programme, silencing or W-direction families
**Problem.** The §7.1 gate can pass on hotspot/amplification addiction alone, and these are the easiest cases (Lord 2020). The families PROMPT.md names (methylation/silencing, signatures, expression programmes) and the W direction have no controls. Part 1 already curated 43 + 19 pairs (D22/D24: EP300–CREBBP, ME2–ME3, PTEN–PIK3CB, TP53-WT → MDM4/PPM1D, lineage TFs), and the design ignores them. In Part 1 tables, all of the genes below are variable and high-tier, so they are testable: KIF18A, ADAR, NFE2L2, UBC, DDX3X, MDM2, PPM1D, USP7, CDK6 and CDK4.
**Change.** Load `P1/configs/biomarker_pairs.csv` and `oncogene_pairs.csv` as the base control set, then add pre-declared controls by family, each with feature, sign and citation:
- **SIG:** WGD/aneuploidy → KIF18A (S; Quinton 2021, Marquis 2021, Cohen-Sharir 2021 [from-memory]; KIF18A inhibitors are in phase I, so this is also a forward-looking CGT row).
- **PROG:** interferon/ISG programme → ADAR (S; Gannon 2018, Liu 2019 [from-memory]); NRF2 programme or KEAP1 LOF → NFE2L2 (S; non-paralog LOF).
- **SIL:** MLH1 SIL → WRN (S; the methylation route to MSI); UBB SIL → UBC (S; Tsherniak 2017).
- **Sex chromosomes:** DDX3Y absent (female, or loss of Y in males) → DDX3X (S) [from-memory]. The SIL rule "expressed in ≥ 50% of lines" excludes Y genes, so define SIL for Y-linked genes within males only and add female sex as a 𝒟 feature for Y-paralog tests. Sex otherwise stays a covariate.
- **AMP (clinical, forward-looking):** CCNE1 AMP → PKMYT1 (S; Gallo 2022; lunresertib trial).
- **p53 axis:** TP53-WT → MDM2, MDM4, PPM1D, USP7 (S).
- **W direction:** RB1 LOF → CDK4, CDK6 and CCND1 (W). RB1-null lines do not need CDK4/6. This is also a clinical resistance biomarker.
- **Lineage-only expected:** SOX10, MITF, PAX8 and similar lineage TFs are expected to land in failure class 3. They are a check on class assignment, not on recovery.
Score the gate per family (≥ 50% within each family with ≥ 2 testable controls, plus the overall 60%), so that addiction controls cannot carry it.
**Check.** `configs/positive_controls.csv` has a `family` column with ≥ 2 testable rows in each of HOT/AMP, paralog/collateral, SIL, SIG, PROG and W. `results/gate.json` reports recall per family and a pass/fail per family.

### R1-PW7 (major): the 𝒟 has no mutational signatures, although PROMPT.md names them
**Problem.** SIG holds only MSI, TMB, ploidy/WGD/CIN and LoH. The listed gaps (methods_prior_work.md, "Gaps") include MMR, APOBEC, POLE and HRD signatures as biomarkers beyond the MSI label. As written, the design cannot answer the "mutational signature" part of the goal.
**Change.** If `OmicsSomaticMutations.csv` arrives:
- compute per-line SBS spectra and fit exposures by NNLS to a pre-declared COSMIC SBS subset: APOBEC (SBS2 + 13), MMR (SBS6/15/26/44), POLE (SBS10a/b), HRD (SBS3), UV (SBS7a–d) and tobacco (SBS4);
- if the reference catalogue or trinucleotide context (hg38 ±1 bp) cannot be obtained within the budget, fall back to the 6-class substitution spectrum plus indel fraction from REF/ALT, which needs no reference;
- restrict the calls to somatic-filtered, non-gnomAD-common variants, because cell lines have no matched normal.
Sanity checks (pre-set): MMR exposure correlates with MSIScore (ρ > 0.5); UV exposure is highest in skin lines; POLE exposure is enriched in POLE-exonuclease-mutant lines. Test the signatures as 𝒟 features. Expected classes: HRD → PARP1 is N; MMR/indel → WRN is S. Otherwise mark the family "not done: data unavailable".
**Check.** `configs/dictionary_v1.json` has SIG-SBS features, or a `not_done` record with its reason. `results/sanity/signatures.json` holds the three sanity statistics and their pass/fail.

### R1-PW8 (major): TP53 status and the Cas9 DNA-damage response are missing from the confounder model
**Problem.** Cas9 cutting triggers a p53 response. TP53-WT lines show stronger depletion and selection against edited cells (Haapaniemi 2018, Ihry 2018 [from-memory]). TP53 status is also correlated with lineage, WGD/CIN and TMB. §5.3 includes TMB/MSI but not TP53, although COMP has it. The combinatorial_confounders survey lists it as a covariate.
**Change.** Add TP53 status (Part 1 composite `TP53_mut`, unknowns kept as their own level) to the §5.3 covariates and report its confounding delta. Omit it only when d is TP53-derived (TP53-WT controls, or MDM2 amp). Add a sensitivity check: any SIG/CIN/WGD determinant must keep its sign within TP53-mutant lines.
**Check.** The covariate spec in `configs/confounders.json` includes `TP53_mut`. A test asserts it is excluded exactly when d's source gene is TP53 or MDM2.

### R1-PW9 (major): the CN-cutting artefact rule relies on "literature support" rather than an empirical test
**Problem.** §5.4 exempts an amplified oncogene "with literature support", which is unscriptable and favours known biology. The established diagnostic (Aguirre 2016; Munoz 2016) is that *unexpressed* genes in the same amplicon also deplete. Chronos's CN correction can also over-correct true amplified-driver addiction (CRISPRcleanR was designed to avoid that).
**Change.** For each AMP/locus determinant whose span contains g, take carrier lines and compute:
- (i) the median effect of the unexpressed amplicon genes (log2(TPM+1) < 1) in carriers vs non-carriers;
- (ii) within carriers, the correlation of g's effect with g's expression.
Label the determinant "cutting artefact" if (i) is < −0.2 and (ii) is not negative-significant. Label it "amplified-driver addiction" if (i) is ≈ 0 and g is depleted. Otherwise label it "unresolved". Drop the literature exception. The rule is pre-set. ERBB2, MYCN and CCND1 are expected to pass, and are therefore positive controls for the test itself.
**Check.** `results/confounders/amplicon_test.csv` has one row per AMP determinant, with columns `n_unexpressed_neighbours, neighbour_effect_delta, r_expr_effect, label`, and every label belongs to the three-value set.

### R1-PW10 (major): some clinical-ground-truth mappings are biologically wrong for single-gene KO
**Problem.** The table in clinical_ground_truth.md §2, which feeds §7.3, has these problems:
- **Paralog-buffered drug targets.** Capivasertib inhibits AKT1/2/3, trametinib MEK1/2, larotrectinib NTRK1/2/3, CDK4/6 inhibitors both kinases, and erdafitinib FGFR1–4. A single-gene KO is buffered by the paralog, and the relevant kinase is lineage-specific (CDK6 in haematopoietic lines, CDK4 in breast). Scoring AKT1 or MAP2K1 alone would book a buffering miss as a model failure.
- **del(17p)/TP53 → venetoclax/BCL2 = S.** del(17p) defines a high-risk CLL population. It does not mark BCL2 dependence, and TP53 loss, if anything, lowers apoptotic priming. The row is not an S pair.
- **ER+ → CDK4.** This is a pathway-node pair with lineage-specific kinase usage.
**Change.**
- Add the KO-concordance class *paralog-buffered target*. Its dependency readout is the most negative effect across the drug's target family, pre-declared per row, and recall is reported separately from the single-target rows.
- Reclassify del(17p) → BCL2 as N, or exclude it.
- Add RB1 loss → CDK4/6 (W, resistance) and CCNE1 AMP → PKMYT1 (trial-stage S).
- Keep the `freeze-cgt` tag, and make these changes before it.
**Check.** In `configs/clinical_ground_truth.csv`, the rows for capivasertib, trametinib, larotrectinib/entrectinib, CDK4/6 inhibitors and erdafitinib have `ko_class = paralog_buffered` and a non-empty `target_family`. The venetoclax-del(17p) row is N or absent. A test asserts both.

### R1-PW11 (minor): LOF zygosity, allele classes and fusion features are too coarse for mechanism claims
**Problem.**
- Paralog synthetic lethality needs complete loss of A, but "dosage > 0" includes heterozygous passengers. STAG2 is X-linked, so one hit suffices in male lines.
- HOT collapses alleles. BRAF class II/III (non-V600) and KRAS G13D behave differently from V600 and G12/Q61, and the CGT rows are allele-specific (BRAF V600E, KRAS G12C).
- F_gene partner flags include read-throughs and passengers (Picco 2019: most fusions are non-functional). A fusion is only hit by a KO of its partner when the guides target exons retained in the fusion.
**Change.**
- Add LOF-bi: dosage 2, or dosage ≥ 1 with LOSS/LoH at the locus, or dosage ≥ 1 with SIL. Count X-linked dosage 1 in male lines as biallelic. Use LOF-bi for the paralog/TSG controls (SMARCA4, ARID1A, STAG2, KEAP1).
- If the full mutation file arrives, add allele features BRAF_V600, BRAF_nonV600, KRAS_G12C, KRAS_G12D/V, KRAS_G13D, KRAS_Q61, PIK3CA_helical and PIK3CA_kinase.
- Restrict FUS to in-frame fusions whose gene is in a curated fusion-driver list (CIViC/MOAlmanac fusion profiles ∪ OncoKB fusion genes). Mark fusion misses as "guide position unverified".
**Check.** The dictionary has LOF-bi features. The BRAF_V600 feature exists if the mutation file is present (otherwise there is a `not_done` record). The FUS features are a subset of `configs/fusion_drivers.csv`.

### R1-PW12 (minor): culture-artefact gene sets and the "cross-lineage" definition are under-specified
**Problem.**
- "Integrin/actin or lipid/cholesterol programmes" is a free-text judgement, and it omits the best-documented medium-sensitive class: ferroptosis/selenoprotein genes (GPX4, SEPHS2, SEPSECS, PSTK, EEFSEC, SECISBP2, TXNRD1, SLC7A11). The mesenchymal → GPX4 link (Viswanathan 2017 [from-memory]) is real biology, but its magnitude depends on the medium. Amino-acid genes (ASNS, PYCR1) are also medium-sensitive (Rossiter 2021).
- "Cross-lineage" can be satisfied by OncotreeLineages that are one tumour biology split in two (Bowel and Esophagus/Stomach adenocarcinoma; Myeloid and Lymphoid). Chan 2019 also reports an MSI × lineage interaction for WRN.
**Change.**
- Create `configs/culture_genesets.json` with pre-declared lists: integrin/FAK/actin, cholesterol/SREBP, ferroptosis/selenoprotein and amino-acid synthesis. Flag any subtype whose genes overlap a list (hypergeometric, BH < 0.05). For flagged subtypes, run the medium/serum-free stratified test.
- Define "cross-lineage" as ≥ 2 lineages from distinct pre-declared organ-system groups (`configs/lineage_groups.csv`: haematopoietic, GI, neural crest, and so on).
- Report the d × lineage interaction for every cross-lineage determinant.
**Check.** Both config files exist and are committed before `freeze-catalogue`. The `cross_lineage` flag in `results/subtypes.csv` is recomputable from those files. Every cross-lineage row has ≥ 2 distinct `organ_group` values.
