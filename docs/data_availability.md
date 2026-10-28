# Data availability (probed 2026-10-28)

Probes used python requests / curl. Small probe files are in `external/probe/` (git-ignored). Status codes: present, scriptable, needs-user-download, licensed-restricted, unavailable.

## 1. DepMap 26Q3 files not yet present

The 26Q3 `README.txt` (`C:\Users\ritab\Desktop\project\data\README.txt`) lists about 100 files. It contains no methylation, no proteomics, no PRISM or other drug-response file, and no arm-level CNA file. Part 1 found that depmap.org is Cloudflare-blocked for scripts. This session confirmed it: `depmap.org/portal/...` returns an HTML challenge page, `downloads.depmap.org` does not resolve, and the GCS bucket `depmap-external-downloads` returns 403. Every DepMap file below therefore needs a manual download.

Already present: CRISPRGeneEffect, CRISPRScreenMap, ScreenSequenceMap, Model, ModelCondition, OmicsCNGeneWGS, OmicsExpressionTPMLogp1HumanProteinCodingGenes, OmicsFusionFiltered, OmicsSomaticMutationsMatrixDamaging, OmicsSomaticMutationsMatrixHotspot.

| Wanted | 26Q3 file name | Status | Fields / notes |
|---|---|---|---|
| Ploidy, CIN, WGD, LOH fraction, MSI score, aneuploidy | `OmicsGlobalSignatures.csv` (README name; not "OmicsSignatures") | needs-user-download | MSIScore, Ploidy, CIN, WGD, LoHFraction, Aneuploidy. PureCN and MSIsensor2. Small. **Priority 1.** |
| Raw MSI repeats | `OmicsMicrosatelliteRepeats.csv` | needs-user-download (optional) | Upstream of MSIScore. Large. Skip unless needed. |
| Full somatic calls | `OmicsSomaticMutations.csv` (also `OmicsSomaticMutationsMAF.maf`) | needs-user-download | Mutect2, hg38. Includes VepImpact, hotspot flags, VAF, protein change, oncogene/TSG and LoF annotations. Large (hundreds of MB). **Priority 1.** |
| Copy-number segments | `OmicsCNSegmentsWGS.csv` | needs-user-download | Columns: SequencingID, ModelID, IsDefault*, CONTIG, START, END, state, NUM_POINTS_COPY_RATIO, SEGMENT_COPY_NUMBER. Arm-level CN can be derived from this. |
| Arm-level CNA | none in 26Q3 | unavailable | Derive from segments. `OmicsGlobalSignatures.Aneuploidy` is the summary. |
| Sequencing and model mapping | `OmicsProfiles.csv` | needs-user-download | ModelID, ProfileID, IsDefaultEntryForModel, IsDefaultEntryForMC. Needed to map profile IDs. Small. **Priority 1.** |
| Molecular subtypes | `OmicsInferredMolecularSubtypes.csv`, `SubtypeTree.csv`, `SubtypeMatrix.csv` | needs-user-download | One-hot subtype and context membership (MSI, PAM50, KRAS allele, etc.). Small. |
| CRISPR screen-level QC and confounders | `CRISPRConfounders.csv`, `AchillesScreenQCReport.csv`, `AchillesSequenceQCReport.csv` | needs-user-download | Includes ScreenNNMD, ScreenROCAUC, CasActivity, **ScreenDoublingTime**, library, screen-type and growth-pattern booleans. This is the only doubling-time source. Small. **Priority 1.** |
| Gene-level dependency probability and controls | `CRISPRGeneDependency.csv`, `CRISPRInferredCommonEssentials.csv`, `AchillesCommonEssentialControls.csv`, `AchillesNonessentialControls.csv` | needs-user-download | Medium or small. |
| Uncorrected effect and screen-level effect | `CRISPRGeneEffectUncorrected.csv`, `ScreenGeneEffect.csv` | needs-user-download (optional) | Large. |
| Gene annotation | `Gene.csv` | needs-user-download | Small. |
| PRISM or drug screens | not in 26Q3 | unavailable | Only `PortalCompounds.csv` is listed (compound metadata, targets, SMILES). PRISM Repurposing is a separate portal release (Repurposing Public 24Q2), and its files were not listed in the 26Q3 README. Use GDSC as the scriptable substitute. |
| Methylation, RRBS | not in 26Q3 | unavailable | Legacy CCLE RRBS is on the portal "Legacy" tab and was not verified. Cell Model Passports has methylation (see section 3). |
| Proteomics | not in 26Q3 | unavailable | Nusinow 2020 is already present at `data/external_omics/nusinow2020`. |
| WES copy number | `OmicsCNGeneMC_WES.csv` | not wanted | Not maintained. |

## 2. Model.csv and ModelCondition.csv columns (26Q3)

Model.csv has 2,528 rows. Fill rates in brackets.

- Culture: `GrowthPattern` [100%]. Values: Adherent 1472, Suspension 466, Dome 241, Unknown 157, Mixed 139, Spheroid 79. Also `OnboardedMedia` (media ID) [82%], `FormulationID` (text, e.g. "RPMI + 20% FBS + ...") [82%], `SerumFreeMedia` (bool) [82%], `PlateCoating` [3%], `ModelType`, `EngineeredModel` [1%].
- Lineage: `OncotreeLineage` [99%], `OncotreePrimaryDisease`, `OncotreeSubtype`, `OncotreeCode` [91%], `DepmapModelType` [100%].
- Patient: `Sex` (Male 1361, Female 1035, Unknown 129, NaN 29), `PatientRace` [32%, the ancestry proxy; sparse], `Age` [81%], `AgeCategory`, `PrimaryOrMetastasis` [84%], `SampleCollectionSite`, `Stage` [4%], treatment fields (sparse).
- IDs: `ModelID`, `PatientID`, `CCLEName`, `SangerModelID` [48%], `COSMICID` [38%], `RRID`.
- **Not in Model.csv:** MSI, ploidy, doubling time. These come from `OmicsGlobalSignatures.csv` (MSI, ploidy) and `CRISPRConfounders.csv` (ScreenDoublingTime). Cell Model Passports also has msi_status and ploidy.
- Genetic ancestry (inferred) is not in the file list. Only the self-reported `PatientRace` is available.

ModelCondition.csv columns: ModelConditionID, ModelID, ParentModelConditionID, DataSource, **CellFormat** (adherent, suspension and others; the 2D/3D indicator, e.g. "adherent"), Morphology, PassageNumber, **GrowthMedia**, **FormulationID**, **SerumFreeMedia**, PlateCoating, PrescreenTreatmentDrug(+Days), AnchorDrug(+Concentration, DaysWithDrug), `model_id`. Screen type (2DS, 2DO, 3DN...) is in `ScreenSequenceMap.ScreenType` and in CRISPRConfounders.

## 3. External datasets

| Dataset | Status | Exact URL | Version / date | Size | Licence | Fields / notes |
|---|---|---|---|---|---|---|
| GDSC1 fitted dose response | scriptable (verified: HTTP 200, range requests work) | https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC1_fitted_dose_response_27Oct23.xlsx | release 8.5, file dated 2023-10-27 | 29.4 MB | Free for academic use; the site says CC BY 4.0, which should be checked before publication | Cell line, COSMIC_ID, drug, LN_IC50, AUC, RMSE, Z_SCORE, max conc. 8.5 is the latest release found; the directory listing is 404, so newer releases were not checked. |
| GDSC2 fitted dose response | scriptable (verified) | https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC2_fitted_dose_response_27Oct23.xlsx | as above | 21.3 MB | as above | same fields |
| Cell Model Passports model list | scriptable (downloaded: 1,984 models, 51 cols) | https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz | file header dated 2021-04-27 (stale for the model list) | 145 KB | Cell Model Passports terms (academic, CC BY 4.0 for data) | Maps COSMIC_ID, BROAD_ID (ACH-), model_id (SIDM). Has msi_status, ploidy, growth_properties, gender, ethnicity, tissue, cancer_type. Needed to join GDSC to DepMap. Methylation and omics matrices are on cellmodelpassports.sanger.ac.uk/downloads (not probed in depth). |
| CIViC nightly | scriptable (downloaded) | https://civicdb.org/downloads/nightly/nightly-ClinicalEvidenceSummaries.tsv, `nightly-VariantSummaries.tsv`, `nightly-MolecularProfileSummaries.tsv` | nightly, 2026-10-28 | 4.2 MB, 0.6 MB, 0.7 MB | CC0 | Evidence columns: molecular_profile, disease, doid, therapies, evidence_type, evidence_direction, evidence_level, significance, rating, nct_ids, variant_origin. Variants: gene, entrez_id, chromosome, start, stop, hgvs, allele_registry_id. About 4,970 evidence lines. The API `/api/genes/1` returns 410; use the TSVs. |
| FDA companion diagnostics | scriptable but HTML only (downloaded page, 2 tables, 243 rows) | https://www.fda.gov/medical-devices/in-vitro-diagnostics/list-cleared-or-approved-companion-diagnostic-devices-in-vitro-and-imaging-tools | page current to 2026-10-28 | 179 KB | US government public domain | Columns: PMA/BLA number, device, indication, drug, biomarker/claim. Requires an HTML table parse and manual cleaning of the biomarker text. There is no CSV. The pharmacogenomic-biomarkers page URL guessed in the probe returned 404 (not re-searched). |
| OncoKB cancerGeneList | scriptable, no token (downloaded: 1,247 genes) | https://www.oncokb.org/api/v1/utils/cancerGeneList | API info shows oncotree_2025_10_03 | 467 KB | OncoKB terms: academic non-commercial; commercial use needs a licence. Redistribution is restricted. | Fields: hugoSymbol, entrezGeneId, oncokbAnnotated, geneType, occurrenceCount, mSKImpact, foundation, vogelstein, sangerCGC. Use as an annotation only; do not commit the file. |
| OncoKB allAnnotatedVariants and actionable levels | licensed-restricted | https://www.oncokb.org/api/v1/utils/allAnnotatedVariants | n/a | n/a | needs a free academic token (registration); HTTP 401 without one | Without a token, only cancerGeneList and `/info` work. |
| Cancer Genome Interpreter biomarkers | unverified, likely needs-user-download | https://www.cancergenomeinterpreter.org/biomarkers | n/a | n/a | CC BY-NC 4.0 (stated on the site, not re-verified) | The site is now a single-page app. Every guessed path, including cgi_biomarkers_latest.zip, returned the same 1.4 KB HTML shell with 200, so no file URL was confirmed. The user should download the biomarkers file from the Biomarkers page. |
| MOAlmanac database | scriptable (GitHub) | https://raw.githubusercontent.com/vanallenlab/moalmanac-db/main/molecular-oncology-almanac.json (1.49 MB); releases https://github.com/vanallenlab/moalmanac-db/releases (latest tag `v.2026-08-02`) | 2026-08 | 1.5 MB | CC BY 4.0 (data), per repo | The releases carry no assets, so use the repo file or the tarball. The old `moalmanac.org/api` returns 404. |
| cBioPortal public API | scriptable (verified) | https://www.cbioportal.org/api/... ; BRAF in skcm_tcga_pan_can_atlas_2018: `GET /molecular-profiles/skcm_tcga_pan_can_atlas_2018_mutations/mutations?sampleListId=skcm_tcga_pan_can_atlas_2018_all&entrezGeneId=673&projection=SUMMARY` returned 245 mutations, with 448 samples in the study | portal v7.2.0 | small per query | Open, no key; study data under the terms of TCGA (public) | Frequency = distinct mutated samples / 448 (the sample list includes non-primary samples; use the `_sequenced` list for the denominator). Alterations in `proteinChange`, `mutationType`. CNA via the `_gistic` profile. |
| ClinicalTrials.gov API v2 | scriptable (verified) | https://clinicaltrials.gov/api/v2/studies?query.term=MTAP&countTotal=true | API 2.0.5, data 2026-10-28 | JSON, small | Public domain; the terms ask for attribution | 71 studies match "MTAP". 36 are recruiting. Fields can be limited with `fields=NCTId,BriefTitle,InterventionName,Condition,OverallStatus`. Paging uses `nextPageToken`. |
| BioGRID ORCS | present | local: `data/external_screens/biogrid_orcs/` (tarball `BIOGRID-ORCS-ALL-homo_sapiens-2.0.18.screens.tar.gz`, 111 extracted screen files) | ORCS 2.0.18 | n/a | MIT | Also `data/raw_external/biogrid/BIOGRID-ALL-LATEST.tab3.zip` (183 MB, interaction data). |
| Ensembl gene coordinates | present | `data/raw_external/ensembl116_gtf/Homo_sapiens.GRCh38.116.gtf.gz` | Ensembl 116, GRCh38 | n/a | open (EMBL-EBI) | Also present: HGNC, ensembl_paralogs, hart_ceg2/hart_neg, huri, msigdb_hallmark, reactome_gmt, bioplex. |
| COSMIC Cancer Gene Census | licensed-restricted | https://cancer.sanger.ac.uk/census | n/a | n/a | requires a login (302 to /cosmic/login); academic free, commercial paid; redistribution not allowed | Substitutes: OncoKB cancerGeneList has a `sangerCGC` flag (a Boolean derived from CGC); Vogelstein flags are also in that file. Alternatively the user can download CGC manually with their own COSMIC account. |
| UniProt, ProteomicsDB | not probed in depth | | | | | ProteomicsDB site answers 200 (no download probed). |

## 4. Python packages

Python 3.13.1. Installed: pandas 2.2.3, numpy 2.1.3, scipy 1.15.2, statsmodels 0.14.4, lightgbm 4.7.0, torch 2.6.0, pyarrow 19.0.1, networkx 3.4.2, matplotlib 3.10.1, seaborn 0.13.2, openpyxl 3.1.5, requests 2.32.3. Not installed: scikit-learn, xgboost, shap.
`pip install --dry-run shap` succeeds and would install shap 0.54.1 plus slicer 0.0.8 (a cp312-abi3 wheel is available). It would pull in scikit-learn if that is missing; the dry-run output shown listed only shap and slicer, so scikit-learn may be resolved separately and is worth checking before use. Nothing was installed.

## USER DOWNLOAD REQUEST

Please download the following from the DepMap portal (https://depmap.org/portal/data_page/, release **DepMap Public 26Q3**, "All Downloads" tab; search the file name). Save each file unchanged to:

`C:\Users\ritab\Desktop\project_b\external\depmap_extra\`

Priority 1 (needed for covariates and analysis):
1. `OmicsGlobalSignatures.csv` (ploidy, CIN, WGD, LoH fraction, MSI score, aneuploidy)
2. `OmicsProfiles.csv` (ID mapping)
3. `CRISPRConfounders.csv` (includes ScreenDoublingTime and screen QC)
4. `OmicsSomaticMutations.csv` (full mutation calls; large)
5. `OmicsInferredMolecularSubtypes.csv`
6. `SubtypeTree.csv` and `SubtypeMatrix.csv`

Priority 2:
7. `OmicsCNSegmentsWGS.csv` (for arm-level CN)
8. `AchillesScreenQCReport.csv` and `AchillesSequenceQCReport.csv`
9. `CRISPRGeneDependency.csv`
10. `AchillesCommonEssentialControls.csv`, `AchillesNonessentialControls.csv`, `CRISPRInferredCommonEssentials.csv`
11. `Gene.csv`
12. `PortalCompounds.csv` (small; compound-target annotations)

Optional, only if the portal offers them in a different tab (they are not in the 26Q3 README):
13. PRISM Repurposing secondary screen dose-response and treatment-info files (latest release: "Repurposing Public 24Q2" or newer).
14. Legacy CCLE RRBS methylation (TSS 1kb upstream average), from the portal "Legacy" or "Downloads" archive.

Also manual: the Cancer Genome Interpreter biomarkers file (https://www.cancergenomeinterpreter.org/biomarkers, destination `external\cgi\`) and, if the user holds a COSMIC account, `cancer_gene_census.csv` (https://cancer.sanger.ac.uk/census, destination `external\cosmic\`; do not commit).
