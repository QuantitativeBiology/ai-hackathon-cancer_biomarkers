# Review round 1 — data science (DepMap/omics), DESIGN v1

Reviewer role: data, feature dictionary, missing data, IDs and versions, downloads, replication leakage, licensing.
Checked read-only against Part 1 (P1 HEAD 7f94b5a = configs/part1_commit.txt; `git status` in P1 clean):
Model.csv (2554 rows), ModelCondition.csv, ScreenSequenceMap.csv, CRISPRScreenMap.csv, configs/splits.csv,
gene_classes.csv, predictable_genes.csv and the cn/damaging/hotspot/fusion caches (few columns only).
`external/depmap_extra/` is empty: none of the requested downloads has arrived yet.

Facts found while checking, which the items below rely on:
- 9p21 in D (967 lines with CN): r(CDKN2A, MTAP) = 0.79 on linear CN and 0.76 on log CN. The deep-deletion phi
  coefficient is 0.68: 174 lines have both deleted, 128 have CDKN2A deleted with MTAP intact, 4 have MTAP deleted
  with CDKN2A intact. For CDKN2A vs CDKN2B, r = 0.97 and phi = 0.95.
- VPS4B, ENO1 and ENO2 have **0** dev lines with relative CN < 0.25. With CN < 0.75 the counts are 340, 116 and 75.
- In ScreenSequenceMap (non-pDNA screens used in CRISPRScreenMap), screen length is set by the library: Avana is
  always 21 d (1051 screens) and KY always 14 d (322). Only Humagne-CD varies (14–25 d). 227 of 1361 lines have
  2–4 screens, and 220 of them have more than one Days value.
- `CPD` is filled for 2 of 3623 sequences, so there is no local doubling-time proxy.
- The screened `CellFormat` (via ModelConditionID) disagrees with `Model.GrowthPattern` for at least 38 lines:
  24 lines annotated Spheroid and 14 annotated Dome were screened as adherent. `SerumFreeMedia` is missing for 21%
  of screens, and `FormulationID` has 139 levels.
- P1 `library_groups.csv` is inferred (ARI 0.949) and covers dev lines only. True library composition per line is
  available for D and R from the screen maps. Among R lines, 168 are Avana-only, 41 Humagne-CD-only and 9 KY-only.
- MSI fallback: `external/probe/cmp_models.csv.gz` (Cell Model Passports) gives msi_status for 441/1022 dev lines
  (36 MSI) and 111/255 R lines (7 MSI). `Model.ModelSubtypeFeatures` mentions MSI for 18 lines.
- Dictionary sizes: 45 of 165 hotspot genes have ≥ 5 dev carriers. There are 3,062 eligible genes
  ((variable ∨ selective) ∧ tier high/medium; 1,331 high).
- No dev/test PatientID overlap. Fusion calls cover all 1022 D and 255 R lines.
- data_availability.md has two errors: it says scikit-learn is not installed (1.5.2 is installed) and gives
  Model.csv as 2,528 rows (it has 2,554). igraph and leidenalg are **not** installed, which affects §5.2 Leiden.
- OmicsSomaticMutations.csv (per the 26Q3 README) carries ModelID, IsDefaultEntryForModel, AF, ProteinChange,
  LikelyLoF, Hotspot, OncogeneHighImpact/TumorSuppressorHighImpact and AMClass.

---

### R1-DA1 — critical — The locus-collapse rule does not merge 9p21, and merging would erase the MTAP control
**Problem.** §5.4 merges AMP/DEL features with r > 0.8 within 10 Mb, and §9 asserts a unit test that "9p21 merges".
In fact r(CDKN2A, MTAP) = 0.79, so MTAP is not merged and the test as written fails. Lowering the threshold would
merge MTAP into a CDKN2A locus feature. That would destroy the MTAP-specific PRMT5/MAT2A positive control, whose
evidence is the 128 D lines that have CDKN2A deleted but MTAP intact.
**Change.** Keep the gene-level AMP/DEL/LOSS features in 𝒟. Define **locus groups** separately:
- a group is a connected component of genes within 10 Mb with binary-call Jaccard ≥ 0.5 among CN-profiled D lines;
- groups are used for reporting, for the "mechanistic gene" rule, and for a mandatory **resolution test**: for
  determinant d in group L, refit the model on the lines that are discordant between d and each neighbour;
- if there are ≥ 5 discordant carriers and the effect keeps its sign at p < 0.05, d is named at gene level;
  otherwise the determinant is reported at locus level.

The positive-control scoring uses the gene-level MTAP-DEL feature (or MTAP SIL).
**Check.** `tests/test_locus.py` asserts: CDKN2A and CDKN2B are in one group; MTAP is in the same group but stays a
separate feature; discordant counts are written to `results/locus/resolution.csv`; and on toy data a planted
passenger in a group is reported at locus level.

### R1-DA2 — critical — Positive-control feature encodings are undefined; several cannot exist under current 𝒟
**Problem.** §7.1 names the controls by biology, not by 𝒟 feature. Three problems follow:
- VPS4B and ENO1 have 0 deep deletions in D, so "VPS4A–VPS4B loss" and "ENO2–ENO1 loss" cannot be DEL features.
  Their known mechanism is hemizygous loss or low expression (VPS4B CN < 0.75 in 340 lines).
- MSI is not in local data.
- The gate denominator ("controls with ≥ 5 carriers") therefore depends on choices that are made later, which makes
  the gate tunable.

𝒟 also has no family for hemizygous/shallow loss, so CYCLOPS-type and paralog-dosage dependencies cannot be named
at all.
**Change.**
1. Add family **LOSS**: relative CN < 0.75, built only for (a) the dependency gene's own locus and (b) its
   paralog/partner priors, with ≥ 5 carriers.
2. Commit `configs/positive_controls.csv` before any association is computed. Columns: `control_id,
   dependency_entrez, feature_id` (exactly one primary 𝒟 feature, plus optional secondary), `expected_sign,
   n_carriers_D, n_carriers_R` (molecular data only), `gate_eligible`.
3. Fix the encodings, e.g.:
   - VPS4A ← LOSS(VPS4B) primary, EXP-paralog(VPS4B) secondary;
   - ENO2 ← EXP-paralog(ENO1) or LOSS(ENO1);
   - PRMT5/MAT2A ← DEL(MTAP);
   - WRN ← MSI (source per R1-DA5).

**Check.** A script verifies that the commit introducing `positive_controls.csv` is an ancestor of the first commit
that writes `results/md/`. It also recomputes the carrier counts from the raw caches and asserts they equal the file.

### R1-DA3 — major — The own-gene CN covariate absorbs own-locus determinants
**Problem.** §5.3 always adjusts for the own-gene CN of g. For d = AMP(ERBB2) → ERBB2, AMP(MYCN) → MYCN, or
LOSS/DEL of g's own locus (CYCLOPS), the determinant and the covariate are nearly collinear. These positive controls
would then fail by construction, or show huge "confounding deltas". Chronos already CN-corrects, so own-CN
adjustment is a sensitivity analysis, not the primary model.
**Change.** When d is a CN-derived feature (AMP/DEL/LOSS/SIL-DEL) whose gene or locus group contains g, omit the
own-CN covariate from the primary model. Report the own-CN-adjusted coefficient as the cutting/CN-artefact
sensitivity, which feeds the §5.4 Cas9 label.
**Check.** `results/models/covariate_spec.csv` lists the covariates for each (g, d). A test asserts that own-CN is
absent iff d is in g's own locus group.

### R1-DA4 — major — Culture and screen covariates must come from the screened condition, encoded per line
**Problem.** "Library group + screen length (days)" are nearly collinear: Avana = 21 d and KY = 14 d always.
227 lines have multiple screens. P1's `lib_group` is inferred, covers dev lines only, and cannot be reused for R.
"Growth pattern" from Model.csv is not the screened format. SerumFreeMedia has 21% missing values, which would
silently become False or be dropped. Medium formulation has 139 levels. No doubling time exists locally.
**Change.** `scripts/build_screen_covariates.py` writes `configs/screen_covariates.csv` for D ∪ R ∪ set A from
CRISPRScreenMap → ScreenSequenceMap (excluding pDNA) → ModelCondition. It reads no effect values. Columns:
- multi-hot library (Avana, KY, Humagne-CD, other) and `n_screens`;
- `screen_days` = mean over the line's screens, entered **only** as the within-Humagne-CD deviation (0 elsewhere),
  and dropped if VIF > 10;
- `screen_3d` (ScreenType ∈ {3DO, 3DN});
- screened CellFormat (adherent / suspension / 3D / mixed / unknown);
- SerumFree (True / False / unknown, with unknown as an explicit level);
- base-medium family (regex on the FormulationID first token: RPMI / DMEM / IMDM / DMEM-F12 / other / unknown);
- `any_failed_qc` (PassesQC False).

Model.GrowthPattern and P1 `lib_group` are used only for a concordance check. Doubling time enters only if
CRISPRConfounders arrives; otherwise it is "not done: data unavailable" (no proxy is invented).
**Check.** The script asserts no NaN in the covariates (only explicit levels), ARI with P1 `lib_group` ≥ 0.9 on D,
and writes the VIF table to `results/confounders/vif.csv`.

### R1-DA5 — major — Classify the pending downloads and pre-declare the fallbacks
**Problem.** The DATA/§1 request lists 14 items without saying which ones the design needs. Two consequences:
- CNSegments is requested, but 𝒟 has no ARM family.
- CRISPRGeneDependency is target-derived, has no use here, and enlarges the leakage surface.

The MSI fallback "C_instab/M_burden" is not an MSI measure.
**Change.** State in §1 which downloads are which.

**Essential** (each blocks a pre-declared analysis):
- (1) OmicsGlobalSignatures, for MSI (WRN control, hypermutation sensitivity) and ploidy/WGD;
- (2) OmicsSomaticMutations.csv, for allele-level CGT mapping (KRAS G12C; BRAF V600 vs class 2/3; EGFR exon19del/
  L858R vs T790M/ex20ins; IDH1 R132), LikelyLoF/AF for LOF-strong, and real TMB;
- (3) CRISPRConfounders, the only doubling-time source;
- (4) OmicsProfiles, only if (1) is keyed by SequencingID.

**Nice-to-have:** OmicsInferredMolecularSubtypes/SubtypeMatrix (PAM50, HER2 labels for Q4c), PortalCompounds
(drug→target), Achilles QC reports (redundant with (3)), RRBS methylation.

**Drop:** CNSegments (unless an ARM family is added; not advised under the time budget), CRISPRGeneDependency,
Achilles control lists and CRISPRInferredCommonEssentials (P1 gene_classes covers them), and Gene.csv (HGNC is
local).

**Fallbacks:**
- MSI = CMP msi_status ∪ ModelSubtypeFeatures "MSI", with unknown left as NA (never MSS). Coverage is 441/1022 D
  lines and must be reported.
- Without (2), allele-specific CGT pairs are classed "untestable (allele level)" and are never mapped to gene-level
  HOT.

**Check.** `scripts/check_inputs.py` writes `results/data_status.json` (file, present, sha256, fallback used). Every
"not done: data unavailable" row in the report must match a row with present = false.

### R1-DA6 — major — Replication set R: close every leakage path and state what the joint Chronos fit implies
**Problem.** The guard covers only `crispr.parquet` and P1 `results/test/`. R gene effects also reach an agent
through four other paths:
- CRISPRGeneDependency, CRISPRGeneEffectUncorrected and ScreenGeneEffect (if downloaded);
- P1 `REPORT.md` (G4 test section, per-tier test numbers);
- P1 `RESEARCH_TREE.md`;
- `results/test/*/`.

Separately, R effects were co-estimated with D in the same Chronos run: guide efficacies, CN and proximity-bias
correction, and RCPC components were fitted over all lines. R is therefore "held-out lines, same processing", not an
independent cohort. The design does not say this (only P1 pitfall 5 does).
**Change.**
1. The guard deny-list covers every file whose values derive from R gene effects (the three DepMap files above,
   P1 `results/test/**`, `cache/test/**` except `preds/`). Agents building the catalogue are instructed not to read
   the P1 REPORT.md test sections or P1 RESEARCH_TREE before `freeze-catalogue`.
2. Add one sentence to §1: R shares the Chronos fit with D, so replication tests generalisation across lines, not
   across processing.
3. Report replication separately for the R lines screened without Avana (50 lines, Humagne-CD/KY only), as a
   partial library-independence check.
4. Optional, only if a ≤ 30-min probe finds a scriptable file: use Sanger Project Score fold changes (Cell Model
   Passports, CC BY 4.0, independent CRISPRcleanR processing) on the R lines that have Sanger screens.

**Check.** `tests/test_guard.py` asserts a refusal for every deny-listed path before the tag.
`results/replication_access_log.csv` must have no R-target rows with timestamps earlier than the `freeze-catalogue`
tag commit time.

### R1-DA7 — major — The SIL definition is lineage-confounded
**Problem.** "log2(TPM+1) < 1 where the gene is expressed in ≥ 50% of dev lines" marks lineage-restricted genes,
such as epithelial genes, as "silenced" in every line of lineages that never express them. SIL then becomes a
lineage proxy that is collinear with the unpenalised lineage term, or it produces spurious cross-lineage modules.
**Change.** SIL_g(line) = 1 iff all three hold:
- log2(TPM+1) < 1;
- g is expressed (≥ 1) in ≥ 50% of D lines **of the same lineage** (lineages with ≥ 10 D lines; otherwise the
  pan-D rule);
- the line's value is < the lineage median − 2.

Split SIL into exclusive mechanisms, in this order:
- SIL-DEL: own CN < 0.25;
- SIL-MUT: damaging dosage ≥ 1, i.e. possible NMD;
- SIL-METH: only if methylation arrives;
- SIL-other.

Lineage medians are taken from D only.
**Check.** A dictionary test recomputes the counts. For each SIL feature, Cramér's V with lineage is written out,
and features with V > 0.5 are flagged as lineage proxies in `dictionary_v1.json`.

### R1-DA8 — major — Missing-data rules for 𝒟 and targets
**Problem.** Several kinds of missing value are not handled:
- 55 D and 16 R lines have no CN profile, and P1 feature blocks fill these as diploid (0 plus a flag).
- MSI will be missing for most lines under the fallback.
- Chronos effects are NaN for some genes in KY/Humagne-only lines (library NaN patterns).
- The OOF predictions have no NaN, so the surrogate (fit on predictions) and MD/S1 (fit on observations) silently
  use different line sets.

**Change.**
1. Build 𝒟 from the raw P1 caches (`cn`, `damaging`, `hotspot`, `fusion`, `expression`), never from the
   fold-transformed `cache/features` blocks.
2. CN-derived features (AMP, DEL, LOSS, SIL-DEL, locus) are NA for lines without a CN profile. Within such a gene's
   CN, NaN follows the P1 `mutdeps` convention.
3. In the elastic-net surrogate, NA is set to 0 and one unpenalised `CN_absent` indicator is added (likewise
   `MSI_unknown`).
4. In S1, MD, subtype and carrier counts, NA lines are excluded for that family.
5. Every per-(g, d) test uses lines with a non-NaN observed effect, and the ≥ 5-carrier rule is applied on those
   lines. `n_obs`, `n_carriers` and `n_NA` are written per test.

**Check.** A test asserts that every CN-family column is NA for the 55 D lines. The S1/MD output tables must carry
`n_obs/n_carriers/n_NA`, and a script asserts that n_carriers ≥ 5 on every reported row.

### R1-DA9 — major — Freeze every R-side transform from D, and define "accurate-dependent" on R's prediction scale
**Problem.** Fitting the NMF programmes, z-scaling, SIL lineage medians, locus groups and the 20th-percentile
threshold on all lines would leak R into the dictionary. Second, the R predictions come from the all-dev refit,
whose offset differs from the fold-centred OOF scale, so reusing D's 20th-percentile cut on R misclassifies lines.
**Change.**
1. `build_dictionary.py` fits everything on D only and freezes it: the 5,000 most variable genes, the NMF basis
   W (50 factors, seed list), scalers, SIL medians and locus groups.
2. R programme scores are obtained by NNLS projection onto the fixed W.
3. The R accurate-dependent threshold is the per-gene 20th percentile of R's own refit predictions. This is
   target-free and may be computed before the freeze.
4. R carrier counts (molecular data only) are computed before the freeze. They feed the §6 power analysis and a
   pre-declared benchmark restriction to features with ≥ 3 R carriers.

**Check.** A test rebuilds the dictionary artefacts with R rows deleted from all inputs and asserts identical
hashes. `dictionary_v1.json` records the D-only fit hash.

### R1-DA10 — minor — Mutation encodings: allele, zygosity, burden name
**Problem.** HOT is gene-level, so KRAS G12C, G12D and Q61 are pooled, and BRAF class 1 is pooled with classes 2/3.
LOF = damaging dosage > 0 includes heterozygous passenger truncations, which are frequent in MSI lines. Paralog SL
usually needs biallelic loss: ARID1A has 114 D carriers with dosage > 0 but only 52 with dosage 2. "TMB = log1p
damaging count" is actually a damaged-gene count.
**Change.**
- Rename the burden feature `DMG_BURDEN`, and use true TMB only from OmicsSomaticMutations.
- For paralog/partner-prior genes, LOF-strong (dosage 2, or dosage ≥ 1 with own LOSS/DEL, or LikelyLoF with
  AF ≥ 0.7 when the MAF arrives) is the primary feature; LOF-any is a sensitivity.
- HOT stays gene-level in 𝒟. Allele-level variables exist only in the CGT mapping and require the MAF
  (R1-DA5).
- Define hypermutated lines = MSI (R1-DA5 source) ∪ top 5% DMG_BURDEN in D.

**Check.** A dictionary test checks the LOF-strong counts against the raw damaging and CN caches, and the CGT table
has an `allele_level_required` column whose rows are "untestable" when the MAF is absent.

### R1-DA11 — minor — ID, version and environment consistency
**Problem.** data_availability.md misstates the Model.csv row count and the scikit-learn status. Leiden needs
igraph/leidenalg, which are not installed. External sources key on symbols (FDA, MOAlmanac), on Entrez (OncoKB,
CIViC variants, cBioPortal) or on BROAD_ID (CMP), with no declared mapping rule. `crispr.parquet` is a cache whose
checksum is not recorded.
**Change.**
- Correct data_availability.md.
- Either install leidenalg/igraph and pin them in `requirements.txt`, or pre-declare networkx Louvain
  (seed list) as the clustering method.
- All joins go on ModelID and Entrez strings. Symbol → Entrez uses local HGNC (with previous and alias symbols),
  and unmapped rows go to `results/idmap_unmapped.csv`.
- CMP BROAD_ID → ModelID.
- Record the sha256 of every P1 cache read, and spot-check `crispr.parquet` against CRISPRGeneEffect.csv on D
  rows only (100 random dev lines × 100 genes, equal to 1e-6).

**Check.** `scripts/verify_inputs.py` covers the cache hashes and the spot check. The unmapped fraction per source is
reported and must be < 5% for curated CGT rows.

### R1-DA12 — minor — DATA.md and licence compliance
**Problem.**
- DATA.md has no rows yet for files already downloaded to `external/probe/`: CIViC (CC0), OncoKB cancerGeneList
  (non-commercial, no redistribution), CMP model list (CC BY 4.0) and the FDA CDx page (public domain).
- P1 DATA.md never filled its DepMap licence section.
- Parrish 2021 is CC BY-NC-ND (article) and the Klingbeil bioRxiv licence is unchecked.
- The 2026 "predicted paralog dependency map" preprint is model predictions on DepMap lines, so using it as
  §7.4 "screen replication" would be circular.

**Change.**
- `scripts/verify_data_md.py` requires every file under `external/` to have a DATA.md row with source URL,
  version/date, sha256 (matching `configs/external_sha256.txt`), licence and a redistribution flag.
- Record the DepMap terms URL and the CC BY 4.0 citation requirement.
- Committed derived tables from OncoKB, CGI or COSMIC may hold only IDs and booleans, no copied text.
- §7.4 counts only measured combinatorial screens. Predicted-SL resources are excluded or reported separately as
  "predicted".

**Check.** `verify_data_md.py` runs in CI/tests and fails on any unlisted or mismatched file. A grep test checks
that no committed file under `configs/` or `results/` contains OncoKB `oncokbAnnotated`/description fields.
