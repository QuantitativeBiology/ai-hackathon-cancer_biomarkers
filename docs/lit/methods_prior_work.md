# Prior work: dependency biomarkers, subtyping, association statistics, positive controls

Compiled 2026-10-28.

**Tags.** [verified] = abstract or full text retrieved this session. [verified-meta] = bibliographic record retrieved, content from search summaries. [from-memory] = not retrieved.

Not found: specific "Sarthy", "Gerstung-type" or "Kevin" papers. The Pacini 2024 and Behan 2019 methods were not accessible, so their covariate details are [from-memory].

---

## 1. Biomarker discovery from dependency screens

**Tsherniak 2017.** *Defining a Cancer Dependency Map.* Cell. 10.1016/j.cell.2017.06.010 [verified]
501 RNAi screens analysed with DEMETER. It found 769 differentially required genes, and nonlinear (random-forest-type) models over 66,646 features predicted 426 of them (55%). **The top biomarker was expression-based in 82% of models**, including methylation-driven paralog loss (UBB hypermethylation leading to UBC dependency). Relevance: expression dominates predictive models, and that makes interpretation harder (see §3).

**McFarland 2018.** *Improved estimation of cancer dependencies from large-scale RNAi screens…* Nat Commun. 10.1038/s41467-018-06916-5 [verified]
DEMETER2 is a hierarchical Bayesian model with screen-quality parameters and batch normalisation across 3 RNAi datasets. Relevance: screen quality and batch are first-class confounders.

**Meyers 2017.** *Computational correction of copy-number effect improves specificity of CRISPR-Cas9 essentiality screens.* Nat Genet. 10.1038/ng.3984 [verified-meta]
CERES models the gene-independent Cas9 cutting toxicity that scales with copy number (342 lines). Relevance: raw CN-dependency associations are artefact-prone.

**Dempster 2019.** *Extracting Biological Insights from the Project Achilles Genome-Scale CRISPR Screens.* bioRxiv. 10.1101/720243 [verified]
Achilles pipeline: QC, dependency probabilities, screen-quality correction.

**Dempster 2021.** *Chronos: a cell population dynamics model of CRISPR experiments…* Genome Biol. 10.1186/s13059-021-02540-7 [verified]
Chronos is the current DepMap scoring method. It gives better control separation and **stronger biomarker associations**, and it has the lowest CN and screen-quality bias among the methods compared. Relevance: use Chronos scores, but still test residual CN bias.

**Dempster 2020.** *Gene expression has more power for predicting in vitro cancer cell vulnerabilities than genomics.* bioRxiv. 10.1101/2020.02.21.959627 [verified]
Systematic DNA-vs-RNA comparison across 5 viability datasets. Expression wins even for vulnerabilities with canonical DNA markers, and the best models use few features. A curator caveat says part of the gain may be tissue identity encoded in expression. Relevance: our kernel-ridge expression model will partly learn lineage.

**Behan 2019.** *Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens.* Nature. 10.1038/s41586-019-1103-9 [verified]
Project Score: 324 lines, 30 cancer types. Fitness effects were combined with genomic biomarkers and tractability, and WRN was nominated for MSI tumours. Methods [from-memory]: GDSC-style ANOVA on binary cancer functional events, with MSI status and tissue/growth covariates, pan-cancer plus per-cancer-type analyses, and CRISPRcleanR CN correction.

**Pacini 2021.** *Integrated cross-study datasets of genetic dependencies in cancer.* Nat Commun. 10.1038/s41467-021-21898-7 [verified]
Integrates the Broad and Sanger screens, correcting for sgRNA efficiency, CN-driven cutting and batch. The combined data give more power and extra biomarkers. Relevance: institute/library is a batch covariate.

**Pacini 2024.** *A comprehensive clinically informed map of dependencies in cancer cells and framework for target prioritization.* Cancer Cell. 10.1016/j.ccell.2023.12.016 [verified]
930 lines with multi-omics. It finds expression markers beyond driver genes, and **many addictions are gain-of-function rather than synthetic lethal**. PPI-network support was used to prioritise 370 targets across 27 cancer types, which were then mapped to tumour cohorts. Relevance: the closest prior work to our "clinically informed" validation. A secondary Cell Reports study using this resource fit `dependency ~ marker + cancer type`, with BH correction after excluding pan- and never-essential genes [search summary only].

**Ghandi 2019.** *Next-generation characterization of the CCLE.* Nature. 10.1038/s41586-019-1186-3 [verified]
Covers mutation, splicing, methylation, histone marks, miRNA and RPPA for 1,072 lines. These are our feature sources.

**Gonçalves 2022.** *Pan-cancer proteomic map of 949 human cell lines.* Cancer Cell. 10.1016/j.ccell.2022.06.010 [verified]
Proteomics adds protein-specific dependency biomarkers.

**Gonçalves 2020.** *Drug mechanism-of-action discovery through the integration of pharmacological and CRISPR screens.* Mol Syst Biol. 10.15252/msb.20199405 [verified, full-text methods]
**Linear mixed model (LIMIX) per drug–gene pair.** Fixed covariates were institute (Broad/Sanger, which also stands in for library), PC1 of drug response (a growth-rate proxy) and growth condition (adherent/suspension/semi-adherent). The **random effect was a cell-line similarity matrix from fitness profiles**. Significance used a likelihood-ratio test with BH per drug at FDR < 10%. CN was handled upstream with CRISPRcleanR. Relevance: a directly reusable template for our single-feature benchmark.

**Picco 2019.** *Functional linkage of gene fusions to cancer cell fitness…* Nat Commun. 10.1038/s41467-019-09940-1 [verified]
8,354 fusions in 1,011 lines, **mostly non-functional**. Use only recurrent or functional fusions as features.

**Shimada 2021.** *shinyDepMap…* eLife. 10.7554/eLife.57116 [verified]
Combined CRISPR+shRNA efficacy/selectivity, with gene clustering.

**Krill-Burger 2023.** *Partial gene suppression improves identification of cancer vulnerabilities when CRISPR-Cas9 knockout is pan-lethal.* Genome Biol. 10.1186/s13059-023-03020-w [verified]
1,867 genes are pan-lethal by CRISPR, but about 50% of them are selective by RNAi. Relevance: pan-essential filtering removes some selective vulnerabilities, so we should consider RNAi as a secondary readout.

**Misek 2024.** *Germline variation contributes to false negatives in CRISPR-based experiments with varying burden across ancestries.* Nat Commun. 10.1038/s41467-024-48957-z [verified]
Most ancestry–dependency associations are artefacts. Germline SNPs in 1.2–2.5% of sgRNA targets reduce cutting, disproportionately in African-ancestry lines. Relevance: an artefact class to screen for (guide-level SNP masking).

**Shi 2024.** *Building a translational cancer dependency map for The Cancer Genome Atlas.* Nat Cancer. 10.1038/s43018-024-00789-y [verified]
Elastic-net DepMap models were transferred to TCGA after expression alignment and stroma removal. Expression-only models and multi-omic models performed comparably for about 76% of genes. They recovered SMARCA2/4, STAG1/2 and EP300/CREBBP, and validated PAPSS1 lethality driven by PAPSS2 co-deletion with PTEN, which was linked to survival. Relevance: a template for clinical validation.

**Chiu 2021.** *Predicting and characterizing a cancer dependency map of tumors with deep learning (DeepDEP).* Sci Adv. 10.1126/sciadv.abh1275 [verified]
An autoencoder pretrained on tumours, then trained on DepMap and projected onto about 8,000 TCGA tumours. Interpretation was used to propose synthetic essentiality.

**Itzhacky 2020.** *Prediction of cancer dependencies from expression data using deep learning.* Mol Omics. 10.1039/d0mo00042f [verified]
Expression-only deep models.

**Lin 2021.** *Using interpretable deep learning to model cancer dependencies (BioVNN).* Bioinformatics. 10.1093/bioinformatics/btab137 [verified]
A pathway-structured "visible" neural network, evaluated prospectively on 142 lines from a later release (a good design).

**Cai 2024.** *Synthetic augmentation of cancer cell line multi-omic datasets using unsupervised deep learning (MOVE).* Nat Commun. 10.1038/s41467-024-54771-4 [verified]
A VAE imputes missing omics for 1,523 lines, and **SHAP** identifies features that drive clustering and dependency biomarkers. Relevance: the nearest precedent for SHAP-based clustering and biomarker attribution.

**Yu 2025.** *Context-aware gene dependency modeling via graph attention networks (GATDep).* J Transl Med. 10.1186/s12967-025-07501-3 [verified]
A GNN over GSVA features with GNNExplainer attributions. On TCGA, predicted dependencies associate with MSI and chemotherapy response.

**Yin 2025/26.** *Interpretable neural networks prioritize cancer driver genes from genome-wide dependency landscapes.* PLoS Comput Biol. 10.1371/journal.pcbi.1014869 [verified]
Tests whether a mutation is encoded in the dependency landscape, the inverse of our task.

**Gong 2025.** *Pan-Cancer Biomarker Analysis from the Cancer Dependency Map: A Blueprint for Precision Oncology.* bioRxiv. 10.1101/2025.02.07.637152 [verified]
1,150 profiles. Validates DepMap against known therapeutic windows, then tests single- and multi-gene biomarkers, including IHC-ready proteins.

**Herranz-Ors 2026.** *A tumour-derived organoid biobank maps cancer gene dependencies.* Nature. 10.1038/s41586-026-10830-y [verified]
256 clinically annotated organoids, with CRISPR screens in 162. It reports genomic and clinical dependency markers, organoid-specific essentials and KRAS allele effects. Relevance: an external validation set with linked patient context.


---

## 2. Dependency-based structure: co-essentiality, modules, subtypes

**Pan 2018.** *Interrogation of Mammalian Protein Complex Structure, Function, and Membership Using Genome-Scale Fitness Screens.* Cell Syst. 10.1016/j.cels.2018.04.011 [verified]
Co-essentiality recovers complexes and SWI/SNF sub-modules.

**Boyle 2018.** *High-resolution mapping of cancer cell networks using co-functional interactions.* Mol Syst Biol. 10.15252/msb.20188594 [verified]
**PCA-based removal of technical artefacts** before correlation (>300 screens). Community detection compresses co-functional genes into one score. Relevance: compress dependencies into module scores before clustering cell lines.

**Kim 2019.** *A network of human functional gene interactions from knockout fitness screens in cancer cells.* Life Sci Alliance. 10.26508/lsa.201800278 [verified]. **Rauscher 2018.** *Toward an integrated map of genetic interactions in cancer cells (MINGLE).* Mol Syst Biol. 10.15252/msb.20177656 [verified]
The first is a co-essentiality network from 276 screens. MINGLE integrates 85 screens across libraries.

**Wainberg 2021.** *A genome-wide atlas of co-essential modules assigns function to uncharacterized genes.* Nat Genet. 10.1038/s41588-021-00840-z [verified]
Pearson co-essentiality has "ubiquitous false positives". **GLS** whitening, which accounts for non-independence between cell lines (lineage relatedness), gives calibrated p-values. Relevance: cell lines are not exchangeable, and the same problem affects our association tests.

**Gheorghe 2022.** *Optimal construction of a functional interaction network from pooled library CRISPR fitness screens.* BMC Bioinf. 10.1186/s12859-022-05078-y [verified]
**Covariance normalisation is critical**, and CERES/BF scores outperformed Chronos for network construction. Partial correlation was also evaluated. Relevance: the choice of score matters for gene-gene structure.

**Pan 2022.** *Sparse dictionary learning recovers pleiotropy from human cell fitness screens (Webster).* Cell Syst. 10.1016/j.cels.2021.12.005 [verified]
Graph-regularised sparse dictionary learning lets a gene load on several functions. Relevance: build cell-line subtypes on the latent function loadings rather than on raw genes.

**Wang 2017.** *Gene Essentiality Profiling Reveals Gene Networks and Synthetic Lethal Interactions with Oncogenic Ras.* Cell. 10.1016/j.cell.2017.01.013 [verified]
AML within-lineage example: 14 lines split into Ras-dependent and Ras-independent groups.

**Ding 2025.** *Comprehensive investigation of the molecular basis of cancer dependencies… breast cancer.* Cancer Biol Med. 10.20892/j.issn.2095-3941.2025.0290 [verified]
47 breast lines. Per-omic linear regression with intrinsic subtype as a covariate, plus **NMF clustering of dependencies** into 2 clusters (mitochondrial vs signalling/replication), followed by cluster-specific association tests. Relevance: a small-scale version of our design.

**Warren 2021.** *Global computational alignment of tumor and cell line transcriptional profiles (Celligner).* Nat Commun. 10.1038/s41467-020-20294-x [verified]
Several hundred lines from diverse lineages share a **mesenchymal/undifferentiated state with distinct dependencies**, so a cross-lineage subtype exists and is partly a culture-adaptation axis.

**Neiswender 2026.** *A dependency map enhanced with next-generation 3D cancer models.* Nature. 10.1038/s41586-026-10843-7 [verified, full text]
147 organoid/spheroid screens. Minipool screens separate effects: **growth format altered integrin/actin dependencies, while culture medium altered lipid/cholesterol dependencies**. Under serum-free media, format barely mattered (r≈0.99). Organoid library effects were confounded with model type, which they handled by pretraining Chronos. GBM NextGen models clustered as "glial" and traditional lines as "mesenchymal". Biomarker classes reported: expression addiction, paralog-expression dependency, GOF/LOF.

**Dharia 2021.** *A first-generation pediatric cancer dependency map.* Nat Genet. 10.1038/s41588-021-00819-w [verified]
82 low-mutation paediatric lines with adult-like dependency complexity, a natural test of expression-defined biomarkers.

---

## 3. Statistics and confounders

**GDSC ANOVA.** Garnett 2012, *Systematic identification of genomic markers of drug sensitivity in cancer cells*, Nature, 10.1038/nature11005 [verified]. Iorio 2016, *A Landscape of Pharmacogenomic Interactions in Cancer*, Cell, 10.1016/j.cell.2016.06.017 [verified]
Iorio stresses that tissue lineage mediates response. The GDSCTools formulation `Y ~ C(tissue) + C(MSI) + [media/growth] + feature` is described in the GDSCTools preprint [search summary]. The FDR approach (BH, or permutation in some analyses) is [from-memory].

**CN-driven cutting toxicity.**
- Aguirre 2016, *Genomic copy number dictates a gene-independent cell response to CRISPR-Cas9 targeting*, Cancer Discov, 10.1158/2159-8290.CD-16-0154 [verified]. Cuts in amplified regions cause G2 arrest even for unexpressed and intergenic targets.
- Munoz 2016, *CRISPR Screens… Generate False-Positive Hits for Highly Amplified Genomic Regions*, Cancer Discov, 10.1158/2159-8290.CD-16-0178 [verified].
- Iorio 2018, *Unsupervised correction of gene-independent cell responses to CRISPR-Cas9 targeting (CRISPRcleanR)*, BMC Genomics, 10.1186/s12864-018-4989-y [verified]. Segment-based correction that keeps amplified-driver signal.

**Integration and reproducibility.** Lord 2020, *Integrative analysis of large-scale loss-of-function screens identifies robust cancer-associated genetic interactions*, eLife, 10.7554/eLife.58925 [verified]
Only 220 driver-associated interactions replicate across independent screens and **non-overlapping cell-line panels**. Oncogene addiction is more robust than oncogene-associated synthetic lethality, and robust pairs are enriched for PPI. The reviewer noted that power biases the results toward frequent drivers.

**Medium.** Rossiter 2021, *CRISPR screens in physiologic medium reveal conditionally essential genes in human cells*, Cell Metab, 10.1016/j.cmet.2021.02.005 [verified-meta]. HPLM vs RPMI changes hundreds of essential genes in K562 cells.

**Hypermutation and MSI.** Chan 2019 (below) adjusted the WRN–MSI association for indel-driven loss of other genes and paralogs. It also found that MSI lines from atypical lineages were less WRN-dependent. This is a lineage × MSI interaction, not a main effect [search summary of extended data].

**Drug-target sanity check.** Lin 2019, *Off-target toxicity is a common mechanism of action of cancer drugs undergoing clinical trials*, Sci Transl Med, 10.1126/scitranslmed.aaw8412 [verified]. Clinical "targets" are often dispensable in CRISPR screens.

---

## 4. Positive controls (mechanistic dependency–biomarker pairs)

| Pair | Mechanism | Citation |
|---|---|---|
| WRN ↔ MSI | Synthetic lethality with MMR deficiency, needs helicase activity | Chan 2019 Nature 10.1038/s41586-019-1102-x [verified]; Behan 2019; Lieb 2019 eLife 10.7554/eLife.43333 [verified]; Kategaya 2019 iScience 10.1016/j.isci.2019.02.006 [verified] |
| PRMT5 / MAT2A / RIOK1 ↔ MTAP del (CDKN2A co-deletion) | MTA accumulation inhibits PRMT5 | Marjon 2016 Cell Rep 10.1016/j.celrep.2016.03.043 [verified]; Mavrakis 2016 Science 10.1126/science.aad5944 [verified-meta]; Kryukov 2016 Science [from-memory]; Kalev 2021 Cancer Cell 10.1016/j.ccell.2020.12.010 [verified] |
| SMARCA2 ↔ SMARCA4 loss | Paralog SWI/SNF ATPase | Hoffman 2014 PNAS 10.1073/pnas.1316793111 [verified-meta] |
| ARID1B ↔ ARID1A loss | Paralog SWI/SNF subunit | Helming 2014 Nat Med 10.1038/nm.3480 [verified] |
| STAG1 ↔ STAG2 loss | Cohesin paralog | van der Lelij 2017 eLife 10.7554/eLife.26980 [verified] |
| VPS4A ↔ VPS4B loss (18q/SMAD4); VPS4B ↔ VPS4A loss (16q/CDH1) | Collateral paralog loss, >30% of cancers | Neggers 2020 Cell Rep 10.1016/j.celrep.2020.108493 [verified] |
| ENO2 ↔ ENO1 del (1p36) | Collateral lethality | Muller 2012 Nature 10.1038/nature11331 [verified] |
| PAPSS1 ↔ PAPSS2 del (with PTEN) | Collateral paralog loss | Shi 2024 [verified] |
| UBC ↔ UBB hypermethylation | Methylation-silenced paralog | Tsherniak 2017 [verified] |
| KRAS/BRAF/PIK3CA/ESR1-type addiction | Oncogene or expression addiction | Lord 2020; Pacini 2024; Neiswender 2026 [verified] |

Paralog resources: De Kegel 2019, PLoS Genet, 10.1371/journal.pgen.1008466 [verified] (paralog buffering explains 13–17% of variable essentiality). De Kegel 2021, Cell Syst, 10.1016/j.cels.2021.08.006 [verified-meta]. Dede 2020, Genome Biol, 10.1186/s13059-020-02173-2 [verified] (half of expressed genes are never essential and are enriched for paralogs). Thompson 2021, Nat Commun, 10.1038/s41467-021-21478-9 [verified]. Parrish 2021, Cell Rep, 10.1016/j.celrep.2021.109597 [verified]. Kebabci 2026, bioRxiv, 10.64898/2026.01.19.700065 [verified] (cell-line-specific paralog SL predictor).

**Known artefacts and negative controls:**
1. Amplicon genes, including unexpressed genes, scoring as dependent through CN cutting (Aguirre, Munoz).
2. Guide-SNP false negatives tied to ancestry (Misek).
3. Lineage-confounded associations, where a marker that is merely enriched in a lineage "predicts" a lineage dependency (Iorio 2016; Dempster 2020 caveat; Wainberg GLS).
4. Non-functional fusions and passenger mutations in hypermutators (Picco; Chan).
5. Medium- and format-driven dependencies: lipid/cholesterol and integrin/actin genes (Neiswender; Rossiter).

---

## Design implications

1. **Use Chronos gene effects and keep a CN-residual check.** After scoring, regress each gene's effect on local CN, and flag any biomarker where CN gain of a region predicts dependency on genes in that same region (Meyers; Aguirre; Dempster 2021).
2. **Single-feature benchmark = GDSC ANOVA plus an LMM variant.** Fit `dep ~ lineage + MSI + growth_type + institute/library + screen_quality + feature`. Also fit Gonçalves-style LIMIX with a random effect from a cell-line similarity kernel and a likelihood-ratio test. Report both, plus BH-FDR, effect size and minimum group size (n≥3–5 mutants).
3. **Treat cell lines as non-exchangeable.** Use GLS/LMM or within-lineage permutation (permute labels inside lineage strata) for empirical FDR, because Pearson/OLS p-values are inflated (Wainberg; Gheorghe).
4. **Separate lineage from biomarker signal explicitly.** For every subtype–biomarker link, report pan-cancer, lineage-adjusted and leave-one-lineage-out results. Expression-based biomarkers (which top Tsherniak and Dempster 2020) are partly lineage proxies.
5. **Cluster on denoised, compressed dependency space.** Remove pan-essential and never-essential genes, then apply PCA-artefact removal (Boyle) or Webster/NMF loadings (Pan 2022; Ding 2025) before clustering. Check stability with bootstrap/consensus clustering and across Broad/Sanger (Pacini 2021).
6. **Report the lineage composition of each subtype.** Require cross-lineage subtypes to persist after lineage-centering dependencies. Check the Celligner "mesenchymal/undifferentiated" axis as a likely culture-adaptation cluster (Warren; Neiswender).
7. **Model culture covariates.** Encode growth pattern (adherent/suspension/spheroid/organoid), medium (serum vs serum-free) and library. Flag subtypes driven by integrin/actin or lipid/cholesterol modules as possible culture artefacts (Neiswender 2026; Rossiter 2021).
8. **Handle co-deletions at the locus level.** Collapse co-deleted neighbours (CDKN2A/MTAP/9p21, PTEN/PAPSS2, SMAD4/VPS4B, CDH1/VPS4A) into locus features. Attribute mechanism only when the paralog or metabolic partner is the deleted gene (Marjon; Neggers; Shi).
9. **MSI and hypermutation.** Include MSI/TMB as covariates. For mutation features, use putative-damaging or hotspot filters and recurrence rather than any mutation. Test MSI × lineage interactions (Chan 2019).
10. **Attribution for the predictor.** Use SHAP (exact TreeSHAP for LightGBM). For kernel ridge, use either linear-kernel weights or SHAP/permutation importance on held-out folds. Aggregate attributions per subtype and require agreement with the univariate benchmark (Cai 2024; Lin 2021).
11. **Positive-control recall as a gate.** Before biological claims, the pipeline must recover WRN–MSI, PRMT5/MAT2A–MTAP, SMARCA2–SMARCA4, ARID1B–ARID1A, STAG1–STAG2, VPS4A–VPS4B and ENO2–ENO1. It must also not "discover" amplicon or ancestry artefacts.
12. **Reproducibility across panels.** Discover in one institute or library and replicate in the other, using non-overlapping lines (Lord 2020). Use NextGen/organoid screens as an external set (Neiswender; Herranz-Ors 2026).
13. **Clinical validation.** Align expression to TCGA (Celligner or the Shi 2024 pipeline), project the subtypes onto tumours, and test concordance with clinical biomarkers: MSI-H, MTAP IHC/9p21 loss, SMARCA4/ARID1A loss and ER/HER2. Then test association with outcome or response, adjusting for stage and lineage.
14. **Probe pan-lethal genes with RNAi.** Use DEMETER2 data as a secondary readout for selectivity among CRISPR pan-lethal genes (Krill-Burger 2023).

---

## Gaps our project could fill

- **Predictor-derived subtypes are rare.** Prior work clusters observed dependencies (Ding 2025; Neiswender 2026) or attributes single-gene predictions (DeepDEP, BioVNN, MOVE). Few works cluster cell lines on **model-explained (attribution-based) dependency structure** and then name the biomarkers that define each subtype.
- **No systematic comparison of multivariate vs single-feature biomarker discovery on shared positive and negative controls.** That comparison should include calibration under lineage-stratified permutation.
- **Confounder-aware cross-lineage subtypes.** Lineage, culture format and medium, and CN-cutting effects are usually handled one at a time. A joint audit, with subtypes re-derived after each correction, is missing.
- **Biclustering of dependencies together with omics is underexplored.** I found SynLeGG-type tools and NMF, but no established biclustering benchmark.
- **Mutational and expression signatures as biomarkers** (e.g., MMR/APOBEC signatures beyond the MSI label, methylation-silenced paralogs) are seldom tested systematically against dependencies.
- **Clinical concordance at subtype level.** Most translational work maps single dependencies to TCGA (Shi 2024; DeepDEP). Validating whole subtypes against clinical biomarkers and outcome is open, and organoid biobanks with patient data (Herranz-Ors 2026) now make it feasible.
