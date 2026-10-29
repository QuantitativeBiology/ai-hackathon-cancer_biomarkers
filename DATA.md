# DATA — sources, versions, checksums, licences

No raw or restricted data are committed. Part 1 inputs are read in place from `../project/data/` (read-only);
their checksums are in `../project/configs/depmap_local_sha256.txt` and are re-verified by `scripts/verify_inputs.py`.
New downloads go to `external/` (git-ignored); their checksums are recorded in `configs/external_sha256.txt`.

## Inherited from Part 1 (read-only)
| Dataset | Source | Version | Checksum | Licence |
|---|---|---|---|---|
| DepMap CRISPRGeneEffect, Model, ModelCondition, OmicsCNGeneWGS, OmicsExpressionTPMLogp1HumanProteinCodingGenes, OmicsFusionFiltered, OmicsSomaticMutationsMatrixDamaging, OmicsSomaticMutationsMatrixHotspot, CRISPRScreenMap, ScreenSequenceMap | DepMap portal | Public 26Q3 | ../project/configs/depmap_local_sha256.txt | CC BY 4.0 |
| Part 1 frozen model outputs (CV predictions, splits, tiers) | ../project (git commit recorded in configs/part1_commit.txt) | Part 1 HEAD | git commit hash | as Part 1 |

## New in this project
| Dataset | Source URL | Version / date | Checksum ref | Licence | Status |
|---|---|---|---|---|---|
| CIViC nightly (ClinicalEvidence, Assertion, MolecularProfile, Variant summaries) | https://civicdb.org/downloads/nightly/ | nightly 2026-10-28 | configs/external_sha256.txt (external/clinical/civic/) | CC0 1.0 | downloaded |
| MOAlmanac database release | https://doi.org/10.5281/zenodo.17089309 (record 21992952) | v.2026-08-02, DOI 10.5281/zenodo.21992952 | configs/external_sha256.txt (external/clinical/moalmanac/) | Repo and Zenodo state GPL-2.0; CC BY 4.0 for the data is unconfirmed, verify before reuse | downloaded |
| FDA list of cleared/approved companion diagnostics (raw HTML + parsed CSVs) | https://www.fda.gov/medical-devices/in-vitro-diagnostics/list-cleared-or-approved-companion-diagnostic-devices-in-vitro-and-imaging-tools | content current as of 2026-09-28; fetched 2026-10-28 | configs/external_sha256.txt (external/clinical/fda_cdx/) | US government public domain | downloaded; 232 rows parsed, 1 flagged |
| OncoKB cancerGeneList, allCuratedGenes (public endpoints) | https://www.oncokb.org/api/v1/utils/ | data v7.6 (2026-09-30) | configs/external_sha256.txt (external/clinical/oncokb/) | Academic non-commercial terms unconfirmed; do not redistribute or commit | downloaded |
| Parrish 2021 pgPEN counts (GEO GSE178179) + analysis code (Zenodo 10.5281/zenodo.5081113) | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178179/suppl/ | GEO update 2021-10-27 | configs/external_sha256.txt (external/combinatorial/parrish2021/) | Article CC BY-NC-ND; code MIT | partial: counts only; processed GI tables (Cell Rep supplement) blocked |
| Dede 2020 Genome Biol supplement (PMC7558751) | https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7558751/supplementaryFiles | Genome Biol 21:262 | configs/external_sha256.txt (external/combinatorial/dede2020/) | CC BY 4.0; data CC0 waiver | downloaded |
| Esmaeili Anvar 2024 in4mer (figshare) | https://doi.org/10.6084/m9.figshare.24243832.v1 | v1, modified 2024-04-08 | configs/external_sha256.txt (external/combinatorial/esmaeili2024/) | CC BY 4.0 | downloaded (notebooks skipped) |
| bioRxiv 10.1101/2024.07.16.603642 supplement (brief: "Klingbeil"; record lists Flister et al., Pfizer) | https://www.biorxiv.org/content/10.1101/2024.07.16.603642v3.supplementary-material | v3, 2025-01-16 | configs/external_sha256.txt (external/combinatorial/klingbeil2024_biorxiv/) | bioRxiv licence "cc_no": no reuse licence; analysis only, do not redistribute | downloaded |
| Thompson 2021 Nat Commun supplement (PMC7910459) | https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7910459/supplementaryFiles | Nat Commun 12:1302 | configs/external_sha256.txt (external/combinatorial/thompson2021/) | CC BY 4.0 | downloaded |

| DepMap OmicsSomaticMutations, OmicsGlobalSignatures, OmicsProfiles, CRISPRConfounders, OmicsCNSegmentsWGS | DepMap portal (manual download by the user, 2026-10-28; portal blocks scripts) | Public 26Q3 (as requested; files not self-versioned) | configs/depmap_extra_sha256.txt | CC BY 4.0 | present in external/depmap_extra/ |
| CCLE RRBS promoter methylation (CCLE_RRBS_TSS1kb_20181022.txt.gz; 21,338 TSS±1kb loci × ~843 lines, CCLE names) | DepMap portal, CCLE 2019 legacy release (manual download by the user, 2026-10-28) | 2018-10-22 (Ghandi 2019) | configs/depmap_extra_sha256.txt | CC BY 4.0 | present; registered before the first dictionary build |
| PRISM Repurposing | DepMap portal | — | — | CC BY 4.0 | not obtained: not used |
| OncoTree hierarchy | oncotree.mskcc.org | — | — | — | not obtained: indication matching used MOAlmanac OncoTree term names plus DepMap Model.csv names (partial); 6 unmapped codes listed in results/clinical/q4_summary.json |
| TCGA prevalence (cBioPortal) | cbioportal.org | — | — | — | not used: Tier B item cut under tiers.json cut_order (DEV-BUDGET-M3B-M4) |

All of the above except the DepMap rows are fetched by `scripts/fetch_external.py` (manifest: configs/external_manifest.json).
