# LITERATURE.md: synthesis for Part 2 (dependency subtypes explained by molecular biomarkers)

Compiled 2026-10-28 from the four surveys and the Part 1 brief. Tags on references are those of the source surveys: [verified] retrieved, [verified-meta] bibliographic record only, [search] search snippet only, [from-memory]/[memory] not checked.

## 1. Where the field is

**Biomarker discovery from dependency screens.** The dominant design is a per-gene association test adjusted for context. GDSC-style ANOVA (Garnett 2012; Iorio 2016) was carried into CRISPR work by Behan 2019 (324 lines, WRN nominated for MSI) and Pacini 2024 (930 lines, 370 prioritised targets across 27 cancer types). Pacini 2024, the closest prior work, found many addictions are gain-of-function rather than synthetic lethal. Gonçalves 2020 used a linear mixed model (institute, growth condition, growth proxy fixed; cell-line similarity random; likelihood-ratio test, BH per drug), a reusable template. Predictive models agree that expression carries most signal: expression was the top feature in 82% of random-forest models (Tsherniak 2017), and expression beat genomics across five datasets (Dempster 2020), with a caveat that part of the gain is tissue identity. Lord 2020 showed that only 220 driver-associated interactions replicate across non-overlapping panels; oncogene addiction replicates better than oncogene-linked synthetic lethality, and power favours frequent drivers.

**Dependency-based subtyping.** Structure in dependency space has mostly been studied at gene level (co-essential modules: Pan 2018, Boyle 2018, Kim 2019, Rauscher 2018, Wainberg 2021, Pan 2022) rather than cell-line level. Cell-line subtyping is small or lineage-bound: Ras-dependent versus independent AML (Wang 2017), NMF clusters in 47 breast lines (Ding 2025), and the cross-lineage mesenchymal/undifferentiated state with distinct dependencies (Warren 2021), which is partly a culture-adaptation axis (Neiswender 2026). Pearson co-essentiality gives ubiquitous false positives unless non-independence between lines is modelled (GLS whitening; Wainberg 2021; Gheorghe 2022).

**Explainable dependency prediction.** DeepDEP (Chiu 2021), expression-only deep models (Itzhacky 2020), the pathway-structured BioVNN (Lin 2021, prospective test on 142 later lines), GATDep with GNNExplainer (Yu 2025), MOVE with SHAP (Cai 2024) and the TCGA-transfer elastic nets of Shi 2024 attribute single-gene predictions. Part 1 agrees: expression-only ridge misses RAS/RAF hotspot addiction, recovered only by the LightGBM stage (NRAS AUPRC 0.21 versus 0.66).

## 2. Statistical practice and confounders

| Confounder | Mechanism and evidence | Control adopted |
|---|---|---|
| Lineage | Biomarker and dependency co-enriched in a lineage (Simpson); lineage dominates drug response (Iorio 2016) | Pooled model with lineage term, within-lineage tests, leave-one-lineage-out; require direction in at least two lineages |
| Culture format, medium, serum | Format shifts integrin/actin dependencies, medium shifts lipid/cholesterol ones (Neiswender 2026); physiologic medium changes hundreds of essentials (Rossiter 2021); heme is suspension | Covariates from Model/ModelCondition; flag subtypes enriched for these modules; test within strata where nested |
| Library, screen quality, batch | Library tracks lineage and screen length (Part 1: P1 0.317 overall, 0.289 within library, 0.211 within lineage x library); institute/batch (Pacini 2021; McFarland 2018) | Library covariate, within-library checks; QC covariates from CRISPRConfounders (user download) |
| CN and cutting toxicity | Cuts at amplified loci cause gene-independent arrest (Aguirre 2016; Munoz 2016); corrected by CERES, CRISPRcleanR, Chronos (Meyers 2017; Iorio 2018; Dempster 2021) | Chronos scores plus residual check: same-locus CN gain must not predict dependency. Part 1: Chronos removes CN association partly by construction |
| Passengers, hypermutation, MSI | Neutral mutations co-occur with drivers; WRN-MSI adjusted for indel-driven paralog loss and shows lineage x MSI interaction (Chan 2019, secondhand); fusions mostly non-functional (Picco 2019) | TMB/MSI covariates; hotspot or damaging calls; recurrent fusions only; hypermutator-excluded sensitivity |
| Co-location and collateral deletion | MTAP is about 100 kb from CDKN2A (Muller 2012; Marjon 2016; Neggers 2020; Shi 2024) | Collapse co-deleted neighbours to locus features; conditional tests; separate "collateral" class |
| Non-independence | Lines share lineage and patients | LMM/GLS or within-lineage permutation; patient-grouped folds |
| Ancestry, growth, TP53 | Guide-SNP false negatives (Misek 2024); growth rate and TP53 [memory only] | PatientRace (sparse), ScreenDoublingTime, TP53 as covariates |
| Multiple testing | About 18k genes x many features, correlated | BH within biomarker family, hierarchical gene-then-feature control, lineage-stratified permutation as primary calibration, Storey q as sensitivity; simulate null bounds first (Part 1 permutation nulls missed their bands) |

Of 1,867 CRISPR pan-lethal genes about half are RNAi-selective (Krill-Burger 2023).

## 3. Combinatorial CRISPR screens as validation

| Dataset | Cells and scale | Access |
|---|---|---|
| Klingbeil 2024 (preprint) | NCI-H1299, MDA-MB-231 genome-wide, 36,648 pairs; 20-line hybrid; 461-pair meta-analysis | bioRxiv supplementary tables: scriptable |
| Parrish 2021 | PC9, HeLa; 1,030 pairs, 122 synthetic lethal in at least one line | GEO GSE178179: scriptable |
| Dede 2020 | A549, HT29, OVCAR8; 403 pairs, 24 synthetic lethal | Figshare, PMC supplements: scriptable |
| Esmaeili Anvar 2024 | K-562, A549, A375, MELJUSO; 4,435 pairs plus higher-order | Figshare: scriptable |
| Thompson 2021 | A375, MeWo, RPE-1; 1,191 pairs | Supplements scriptable; raw EGA controlled |
| Ito 2021, CHyMErA, Horlbeck, Najm | various | Locations unverified: manual |

Achievable concordance is modest. Per-screen hit rates are 6 to 15%; 68 to 86% of hits are line-specific; mean pairwise overlap between screens is about 31% (RPE-1: 2 of 26 hits shared); inter-screen reproducibility is precision 0.40 and recall 0.32; a feature-based predictor reached average precision about 49% (chance about 11%) and, in a 2026 preprint, ROC AUC 0.81 on unseen pairs. Screened pairs are biased toward close paralogs. So: replication is strong evidence, non-replication weak; use hits in two or more screens or lines as gold standard; compare with a matched background; only the paralog-buffering subset is testable.

## 4. Clinical ground truth

| Resource | Licence | Use |
|---|---|---|
| CIViC nightly TSV | CC0 | Core. 2,912 accepted predictive items; level A/B with "supports" gives 849 items and 630 profile x therapy pairs; NCT ids and DOID present |
| MOAlmanac | CC BY 4.0 | Core. GitHub JSON or REST; FDA indications and biomarker classes |
| FDA CDx list | US public domain | Core. HTML table only, needs parsing and cross-check against CIViC/MOAlmanac |
| OncoKB cancerGeneList / allCuratedGenes | Academic, redistribution restricted | Gene-level annotation only; variant levels need a licence token (HTTP 401) |
| ClinicalTrials.gov v2, cBioPortal | Public | Companion trials by eligibility regex; TCGA prevalence denominators |
| AACR GENIE | Registration; no redistribution | Prevalence, aggregate counts only |
| Cancer Genome Interpreter | CC0 per one survey, CC BY-NC per another; stale (Oct 2022) | Manual download, secondary cross-check |
| ESCAT; JAX-CKB | No file; paid | Manual tiers; CKB excluded |

**Mapping drug to target.** Benchmark rows are tuples (biomarker feature, type, drug, target gene, indication) with a pre-registered sign of the target-knockout effect in biomarker-positive versus negative lines: stronger (S), weaker (W, resistance), or not informative (N). Targets come from the label (via MOAlmanac), not drug names. Pathway drugs map to the drugged node; synthetic-lethal pairs invert roles (MTAP deletion to PRMT5/MAT2A; MSI to WRN). The survey's 52 rows are described as 45 distinct pairs, and its subcounts are inconsistent, so recount when freezing the table.

**KO-drug discordance classes** (to be pre-registered, so failure is not blamed on the model): PARP trapping rather than loss of function; neomorphic IDH1/2; immune and non-cell-autonomous biomarkers (MSI-H, TMB, PD-L1 for checkpoint inhibitors, while MSI is testable for WRN inhibitors); antibody-drug conjugates; partial inhibition and pan-lethality; off-target drug mechanism (Lin 2019); dominant-negative alleles; gatekeeper and resistance alleles; upstream biomarker with downstream drug node; loss of hormone dependence in culture (ER); and 9p21 co-deletion.

**Prior quantification of agreement.** Gonçalves 2020 reported enrichment of drug-sensitivity associations with knockout of the nominal target. SELECT scored a cell-line-derived model against 35 published trials (80% predicted response; Lee 2021). No published study was found reporting precision or recall of DepMap-derived biomarkers against OncoKB or CIViC levels; the metric must be defined here.

## 5. Available data for this project

| Status | Items |
|---|---|
| Present | DepMap 26Q3 CRISPRGeneEffect, screen maps, Model, ModelCondition, expression, WGS CN, fusions, damaging and hotspot matrices; Nusinow 2020 proteomics; BioGRID ORCS 2.0.18 and interactions; Ensembl 116 GTF, HGNC, paralogs, Hart core/non-essential sets, HuRI, MSigDB Hallmark, Reactome, BioPlex; Part 1 caches and out-of-fold predictions |
| Scriptable | CIViC, MOAlmanac, FDA CDx page, OncoKB cancerGeneList (annotation only; do not commit), cBioPortal, ClinicalTrials.gov, GDSC1/GDSC2 release 8.5, Cell Model Passports model list (header dated 2021, stale; has MSI status, ploidy and the ACH/COSMIC/Sanger id map), Parrish GEO, Dede and Esmaeili Anvar Figshare, Klingbeil supplements, 2026 paralog predictions (Figshare) |
| Needs user download (DepMap portal is blocked for scripts) | Priority 1: OmicsGlobalSignatures (MSI, ploidy, CIN, WGD), OmicsProfiles, CRISPRConfounders (only source of doubling time and screen QC), OmicsSomaticMutations, OmicsInferredMolecularSubtypes, SubtypeTree/Matrix. Priority 2: OmicsCNSegmentsWGS (arm-level CN), QC reports, CRISPRGeneDependency, control gene sets, Gene.csv, PortalCompounds. Optional: PRISM and legacy RRBS (not in 26Q3), CGI biomarkers, COSMIC CGC (account required, do not commit) |
| Excluded or unavailable | JAX-CKB (licence); OncoKB variant-level levels (token); GENIE rows (no redistribution; aggregates only); Thompson raw EGA data; methylation and PRISM within 26Q3; arm-level CNA file (derive from segments); WES CN |

Part 1 constraints: 838 lines have expression but no screen (prediction-only subtyping candidates); the test targets were opened once in Part 1 and the guard now refuses them, so using them is a declared protocol decision.

## 6. Gaps this project addresses

1. Subtypes derived from model-explained (attribution-based) dependency structure, with the defining biomarkers named; prior work clusters observed dependencies or attributes single genes.
2. A matched comparison of multivariate and single-feature discovery on shared positive and negative controls, with calibration under lineage-stratified permutation.
3. A joint confounder audit (lineage, culture, library, CN cutting, co-location) with subtypes re-derived after each correction, rather than one confounder at a time.
4. Subtype-level, indication-aware clinical concordance (direction, actionability, population) against CIViC, MOAlmanac and FDA CDx, with a pre-registered discordance partition and recall by evidence tier.
5. A quantified account of which dependencies cannot serve as biomarkers (the N classes, pan-lethal genes, low-predictability tiers).

Positive-control gate before biological claims: WRN-MSI, PRMT5/MAT2A-MTAP, SMARCA2-SMARCA4, ARID1B-ARID1A, STAG1-STAG2, VPS4A-VPS4B, ENO2-ENO1 (amplicon and ancestry artefacts as negative controls).

## 7. References

**Dependency maps, scoring, prediction**
- Tsherniak 2017, Cell, 10.1016/j.cell.2017.06.010 [verified]
- McFarland 2018, Nat Commun, 10.1038/s41467-018-06916-5 [verified]
- Meyers 2017, Nat Genet, 10.1038/ng.3984 [verified-meta]
- Dempster 2021 (Chronos), Genome Biol, 10.1186/s13059-021-02540-7 [verified]
- Dempster 2020, bioRxiv, 10.1101/2020.02.21.959627 [verified]
- Behan 2019, Nature, 10.1038/s41586-019-1103-9 [verified]
- Pacini 2021, Nat Commun, 10.1038/s41467-021-21898-7 [verified]
- Pacini 2024, Cancer Cell, 10.1016/j.ccell.2023.12.016 [verified]
- Ghandi 2019, Nature, 10.1038/s41586-019-1186-3 [verified]
- Gonçalves 2020, Mol Syst Biol, 10.15252/msb.20199405 [verified]
- Picco 2019, Nat Commun, 10.1038/s41467-019-09940-1 [verified]
- Krill-Burger 2023, Genome Biol, 10.1186/s13059-023-03020-w [verified]
- Misek 2024, Nat Commun, 10.1038/s41467-024-48957-z [verified]
- Shi 2024, Nat Cancer, 10.1038/s43018-024-00789-y [verified]
- Chiu 2021, Sci Adv, 10.1126/sciadv.abh1275 [verified]
- Itzhacky 2020, Mol Omics, 10.1039/d0mo00042f [verified]
- Lin 2021 (BioVNN), Bioinformatics, 10.1093/bioinformatics/btab137 [verified]
- Cai 2024 (MOVE), Nat Commun, 10.1038/s41467-024-54771-4 [verified]
- Yu 2025 (GATDep), J Transl Med, 10.1186/s12967-025-07501-3 [verified]
- Herranz-Ors 2026, Nature, 10.1038/s41586-026-10830-y [verified]
- Neiswender 2026, Nature, 10.1038/s41586-026-10843-7 [verified]
- Warren 2021 (Celligner), Nat Commun, 10.1038/s41467-020-20294-x [verified]

**Modules, subtypes, co-essentiality**
- Pan 2018, Cell Syst, 10.1016/j.cels.2018.04.011 [verified]
- Boyle 2018, Mol Syst Biol, 10.15252/msb.20188594 [verified]
- Kim 2019, Life Sci Alliance, 10.26508/lsa.201800278 [verified]
- Rauscher 2018, Mol Syst Biol, 10.15252/msb.20177656 [verified]
- Wainberg 2021, Nat Genet, 10.1038/s41588-021-00840-z [verified]
- Gheorghe 2022, BMC Bioinf, 10.1186/s12859-022-05078-y [verified]
- Pan 2022, Cell Syst, 10.1016/j.cels.2021.12.005 [verified]
- Wang 2017, Cell, 10.1016/j.cell.2017.01.013 [verified]
- Ding 2025, Cancer Biol Med, 10.20892/j.issn.2095-3941.2025.0290 [verified]

**Statistics and confounders**
- Garnett 2012, Nature, 10.1038/nature11005 [verified]
- Iorio 2016, Cell, 10.1016/j.cell.2016.06.017 [verified]
- Aguirre 2016, Cancer Discov, 10.1158/2159-8290.CD-16-0154 [verified]
- Munoz 2016, Cancer Discov, 10.1158/2159-8290.CD-16-0178 [verified]
- Iorio 2018 (CRISPRcleanR), BMC Genomics, 10.1186/s12864-018-4989-y [verified]
- Lord 2020, eLife, 10.7554/eLife.58925 [verified]
- Rossiter 2021, Cell Metab, 10.1016/j.cmet.2021.02.005 [verified-meta]
- Lin 2019, Sci Transl Med, 10.1126/scitranslmed.aaw8412 [verified]

**Positive controls and paralogs**
- Chan 2019 (WRN), Nature, 10.1038/s41586-019-1102-x [verified]
- Lieb 2019, eLife, 10.7554/eLife.43333 [verified]
- Kategaya 2019, iScience, 10.1016/j.isci.2019.02.006 [verified]
- Marjon 2016, Cell Rep, 10.1016/j.celrep.2016.03.043 [verified]
- Mavrakis 2016, Science, 10.1126/science.aad5944 [verified-meta]
- Kryukov 2016, Science [from-memory]
- Hoffman 2014, PNAS, 10.1073/pnas.1316793111 [verified-meta]
- Helming 2014, Nat Med, 10.1038/nm.3480 [verified]
- van der Lelij 2017, eLife, 10.7554/eLife.26980 [verified]
- Neggers 2020, Cell Rep, 10.1016/j.celrep.2020.108493 [verified]
- Muller 2012, Nature, 10.1038/nature11331 [verified]
- De Kegel 2019, PLoS Genet, 10.1371/journal.pgen.1008466 [verified]
- De Kegel 2021, Cell Syst, 10.1016/j.cels.2021.08.006 [verified-meta]

**Combinatorial screens**
- Dede 2020, Genome Biol, 10.1186/s13059-020-02173-2 [verified]
- Thompson 2021, Nat Commun, 10.1038/s41467-021-21478-9 [verified]
- Parrish 2021, Cell Rep, 10.1016/j.celrep.2021.109597 [verified]
- Esmaeili Anvar 2024, Nat Commun, PMC11055879, Figshare 10.6084/m9.figshare.24243832.v1 [verified]
- Klingbeil 2024, bioRxiv, 10.1101/2024.07.16.603642 [verified]
- 2026 predicted paralog dependency map, bioRxiv, data 10.6084/m9.figshare.31058182 [verified]
- Gonatopoulos-Pournatzis 2020 (CHyMErA), Nat Biotechnol [search]
- Ito 2021, Nat Genet [search]
- Horlbeck 2018, Cell [search]

**Clinical ground truth**
- Tamborero 2018 (CGI), Genome Med 10:25 [from-memory]
- Mateo 2018 (ESCAT), Ann Oncol, PMID 30137196 [verified]
- Mosele 2024, Ann Oncol, PMID 38834388 [verified]
- Reardon 2021 (MOAlmanac), Nat Cancer [verified]
- Ferretti 2024, PMID 38658754 [verified]
- Murai 2012, PMID 23118055 [verified]
- Pettitt 2018, PMID 29748565 [verified]
- Jaaks 2022, PMID 35197630 [verified]
- Lee 2018 (ISLE), PMID 29959327 [verified]
- Lee 2021 (SELECT), PMID 33857424 [verified]
- Sinha 2024 (PERCEPTION), PMID 38637658 [verified]
