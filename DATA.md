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
(filled as data are acquired)
