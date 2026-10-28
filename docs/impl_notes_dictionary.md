# Implementation notes — M3b dictionary, covariates, eligibility, controls (Tier A)

One line per choice the design (v3.1 + rulings R3-1/R3-2) leaves open: **choice — reason**. Numbers live in
`configs/build_params.json` (implementation) or `configs/design_params.json` (decisions).

## Data, rows, deadline
- Dictionary rows = Part 1 population (1022 D, 255 R, 15 set A), molecular data only; all fitting on D rows — DESIGN 1.7; R/A transformed.
- `results/data_status.json` written by `scripts/data_status.py` before the dictionary build; its sha256 is stored in `dictionary_v1.json` — data deadline (DESIGN 1.5, PI ruling 3).
- data_status lists sha256 of every file in external/depmap_extra (configured files + `other_files`) vs configs/depmap_extra_sha256.txt — R3-1(d).
- RRBS (CCLE_RRBS_TSS1kb_20181022) is USED (METH, SIL-METH) — amendment D7 / R3-1; D coverage 453–456/1022 lines (gene-level after mapping), 2 of 843 columns unmapped (NCIH1339_LUNG, NCIH684_LIVER).
- OmicsGlobalSignatures: one row per ModelID = first non-null value over its IsDefaultEntryForModel rows — RNA rows carry PAM50 only, DNA rows carry MSI/CIN; 102 models have both.
- MSI = MSIScore ≥ data_s1.msi_threshold (20); no CMP/ModelSubtypeFeatures fallback; CMP concordance (504 MSS/MSS, 43 MSI/MSI, 8 score-MSI vs CMP-MSS, 0 score-MSS vs CMP-MSI) in results/dictionary/checks.json — DESIGN v3 1.5.
- OmicsSomaticMutations: default entries only; the file holds MODERATE/HIGH (protein-altering) calls only, so TMB and spectra are coding-only — as DESIGN v3 3 states.
- OmicsCNSegmentsWGS, OmicsProfiles unused — DESIGN 1.5.

## Covariates (DESIGN 5.3; `cache/covariates.parquet`)
- Library group from the '+'-joined screened libraries; exact maps Avana / KY / Humagne-CD / Avana+KY, everything else (incl. Brunello, TKOv3, Avana+Humagne-CD) → `other` — design lists five levels.
- Library ARI vs P1 lib_group on D = 0.950 (≥ 0.9 check passes).
- Screen days within Humagne = days − D Humagne-only median (24), 0 for other libraries; VIF of the shared D design max 9.4 → not dropped.
- CellFormat of the screened condition(s): dome/spheroid/organoid → 3D; mixed/unknown/missing/multi-format → mixed/unknown; no Model.csv fallback — design says "screened CellFormat".
- SerumFree and base medium: screened condition, else Model.csv, else unknown — design_params media_column_rule.
- Base medium family parsed per screened condition from the first component of the formulation text (split on '+'); conditions disagreeing → `other`; mixtures (e.g. "IMDM:RPMI (4:1)") → `other` — no family is defined for mixtures.
- Formulation column chosen by content (≥ 90% regex match): FormulationID in both ModelCondition and Model.csv in this release.
- Doubling time NA → D median + DT_unknown — design 5.3; CasActivity/ScreenNNMD/ScreenROCAUC stored as `aux_*` only (not in the confounder model; design lists only doubling time).
- TP53 status from P1 composite TP53_mut; no unknown level occurs (all lines have mutation data).
- Continuous covariates stored: screen_days, screen_days_dev_humagne, doubling_time, MSIScore, DMG_BURDEN, TMB, PLOIDY — PI ruling 4; flags CN_absent, MSI_unknown, METH_absent (R3-1 e), DT_unknown.
- hypermutated = MSI ∨ TMB > D 95th percentile — design_params hypermutated_rule (v3: TMB).
- Own CN covariate (`cache/dictionary/own_cn_log2.parquet`) = clip(log2((CN+0.1)/1.1), ±3), NaN kept for consumers to set 0 + CN_absent — P1 C-block transform without the fill.

## Eligibility (DESIGN 2)
- Strata = pooled lineage × screen-map library group; stratum size counted on the gene's non-NaN D lines; lines in strata < 5 excluded from r_within — design says "strata ≥ 5 lines" only.
- r_within undefined (NaN, not eligible) if < 10 lines remain — build_params eligibility.min_lines_for_r.
- Output `configs/eligible_genes.csv` (task) and the identical `results/eligibility.csv` (DESIGN 2); `r_beyond` = r_within.
- dependent / accurate-dependent rules implemented as functions in `src/biomark/eligibility.py` (D-fitted cuts applied to any lines; R quantile rule target-free).

## Gene lists
- Γ = OncoKB cancerGeneList ∪ top-3 P1 paralogs of eligible genes ∪ spec partners ∪ known_silenced → configs/gamma_genes.csv (3464 genes, source column).
- Partner priors per target = top-3 Ensembl paralogs (P1 priors) ∪ spec feature genes other than the dependency (incl. SIL/YABS feature genes) — design_params partner_priors; PPI partners not used (design names paralogs and spec partners only).
- Targets of gene-specific features = eligible genes ∪ positive-control and lineage-check dependency genes — "controls analysed regardless".
- **DEV-DICT-1** (PI ruling 4): KMT2A added to the 5' driver list (build_params dictionary.fusion_5prime_drivers_extra), extending design_params dictionary_s3.fusion.kinase_3prime_except_5prime (key unchanged); consistent with CGT CLIN:KMT2A_R.
- Fusion drivers = CIViC fusion partners + MOAlmanac rearrangement genes, OncoKB genes only, that appear as the 3' gene of a curated fusion or as a single-gene record; FGFR2/3 and KMT2A 5'; other genes curated only as 5' partners (EWSR1, BCR, EML4…) are not drivers — DESIGN v3 "driver as 3' gene"; OncoKB has no fusion flag, so it serves only as the gene filter.
- DepMap fusion orientation: Gene1 = 5', Gene2 = 3' (verified on EML4::ALK, FGFR3::TACC3); Gene1 == Gene2 rows dropped.

## Families
- HOT over all 165 hotspot genes; AMP over Γ ∪ targets (AMP_g of a non-Γ target is gene-specific); DEL over all genes; LOSS gene-specific (own + partners).
- LOFBI NA rule: 1 if any clause true; NA if dosage ≥ 1 and a clause needing missing CN/sex could make it true; LikelyLoF clause = ≥ 1 LikelyLoF call with AF ≥ 0.7 in the gene.
- SIL lineage stats on D lines of the raw OncotreeLineage if ≥ 10 D lines, else pan-D; Y genes: male D lines only (≥ 10 males in lineage, else pan-male), females/unknown sex NA.
- SIL split (R3-2): mutually exclusive SIL-DEL / SIL-MUT / SIL-METH / SIL-other features, precedence DEL > MUT > METH > other; umbrella SIL kept internal (LOFBI clause, YABS); lines without CN cannot be SIL-DEL and lines without RRBS cannot be SIL-METH, so they fall to the next mechanism.
- METH: gene value = max methylation fraction over the gene's TSS±1kb loci (symbol prefix → Entrez via P1 gene_map/HGNC; 839 locus symbols unmapped); lineage median over RRBS-profiled D lines (lineage ≥ 10 such lines, else pan-D); `meth_above_lineage_median` 0.3 from the design rule text (build_params, as no design_params key exists).
- YABS: Y paralog = a chrY gene among the target's top-3 paralogs; female = 1, male = SIL of the Y gene (within males), unknown sex NA.
- EXPDOSE centre = mean over expressed D lines; EXPPAR raw log2(TPM+1) of top-3 paralogs; EXPGATE stored as covariate (cache/dictionary/expgate.parquet).
- PROG: sklearn NMF on D log2(TPM+1) of the 5,000 highest-D-variance genes (pinned settings, nmf_seed 0); every row (incl. D) scored with model.transform for identical D/R treatment; annotation = Hallmark hypergeometric on top-50 loading genes (BH within component) else top lineage (results/dictionary/prog_annotation.csv).
- PSCORE: member symbols → Entrez; genes with D SD 0 dropped.
- SIG: MSI, MSISCORE, PLOIDY, WGD (binary), CIN, ANEUPLOIDY, LOH from GlobalSignatures (NA for lines without WGS-derived values), TMB = log1p(# default calls), DMG_BURDEN = P1 M_burden, SPEC_* per design (lines < 10 SNVs NA for every spectrum); POLE threshold = D 99th percentile of SNV count; CC>TT = 2-bp `substitution` calls CC>TT or GG>AA.
- COMP: TP53MUT, RB_ALT = P1 composite TP53_mut, RB_altered.
- ALLELE / CLIN (PI ruling 1): built in the dictionary by `src/biomark/cgt_features.py` from the CGT track's read-only configs/cgt_allele_rules.json (+ cgt_extra_features.csv, same feature set), its rule semantics reproduced (mutation classes/regex/codons/exon/DNA, union, difference, fusion with driver end and partners, lofbi, her2_pos, er_pos with the D9 PAM50 gate, del_or_sil, brca_biallelic with hypermutated set to 0, cdk46_sel, cn_below); not_assayable rules are recorded, not built; same ≥ 5 D-carrier floor; `configs/cgt_features.csv` = built list with counts and in_D. HELD by PI until the CGT rule audit is final (one rebuild then).
- Continuous minima: ≥ 20 non-NaN D lines and SD > 0; EXPDOSE ≥ 5 expressed D lines — design v3.
- Cap: shared binary families only; escalate the currently largest shared family not yet at 20 to the next level (10, then 20), repeat until ≤ 3000; gene-specific features keep the floor 5 — "largest family's minimum frequency rises".
- Cramér's V on D with pooled lineage; lineage_proxy flag for SIL-* and METH with V > 0.5.
- Locus groups over kept AMP/DEL/LOSS features plus (DEL) every target with ≥ 1 D DEL carrier; distance = interval gap ≤ 10 Mb; edges need Jaccard ≥ 0.5 on CN-profiled D lines; connected components (groups may span > 10 Mb by chaining); genes without coordinates are singletons.
- Variant groups: CN locus group; VG_SIL:<g> for the SIL mechanisms of a gene; VG_ALLELE:<g> for HOT/ALLELE/CLIN of a gene.
- Self-null per target: DEL_g, LOFBI_g, METH_g, SIL-*_g and DEL_h for h in g's DEL group (R3-1 f).
- Access class (per (g, d), `dictionary.access_class`, column in dictionary_gene_features.csv and positive_controls.csv): direct = own EXPDOSE/AMP/DEL/LOSS/LOFBI or paralog EXPPAR/DEL/LOFBI; potentially-direct = other HOT/LOFBI/AMP/DEL/LOSS; expression-derived = EXP*/PROG/PSCORE/SIL-*/YABS; proxy = SIG/FUS/ALLELE/CLIN/COMP/METH (R3-1).
- dictionary_v1.json omits null/false/empty keys per feature and stores sources per family — keeps the committed file < 5 MB.
- Feature hash = sha256 (16 hex) of the float32 column over population rows (`hash`) and over D rows (`hash_D`).

## Controls (DESIGN 9)
- Spec SIL:<x> expands to SIL-DEL/MUT/METH/other:<x>; carriers = union; build fails (UnresolvableFeature) on unknown family, unresolvable symbol or undeclared named feature (R3-2); ALLELE/CLIN ids are valid if declared in design_params allele_classes / clin_features even if not yet built.
- n_obs_D from cache/effect_nan_mask.parquet (the only effect-derived input; scripts/build_effect_nan_mask.py, dev scope).
- gate_eligible = n_obs_D ≥ 100 ∧ (primary, else secondary, else fallback) in 𝒟_g with ≥ 5 carriers among non-NaN D lines (binary) or present (continuous); `gate_feature` names which one.
- en_max_attainable_freq = fraction of the 50 subsamples where the feature has ≥ min_carriers_in_fit (2) carriers among subsample ∩ non-NaN rows (1.0 for continuous).
- Subsamples (configs/subsamples.json): patients assigned to the half with fewer lines within each pooled-lineage stratum (random order, random tie-break), patient lineage = modal pooled lineage.
- Citation keys checked against LITERATURE.md + docs/lit (`citation_resolved`); unresolved keys are kept and noted.
- No negative-controls file: DESIGN §9 specifies none (NC1–NC6 are §7 bounds).
- Spectra sanity (PI ruling 5): SIG:SPEC_MMR failed its pre-declared check (Spearman with MSIScore 0.27 < 0.5; coding-only calls) and is excluded from 𝒟 as a determinant, recorded under `failed_sanity` in dictionary_v1.json and in summary.csv; the check is unchanged. UV (top lineage Skin) and tobacco (Lung > rest, p 4e-37) pass.
- Fusion known positives (configs/fusion_known_positives.csv): NCI-H3122 and K562 carriers; NCI-H2228 and MDA-MB-175-VII are not in the Part 1 population; MV4-11 becomes a FUS:KMT2A carrier only after DEV-DICT-1 (rebuild pending).
- RRBS coverage (PI ruling 6) reported in summary.csv (METH, SIL-METH rows) and dictionary_v1.json `rrbs_coverage`.

## Build log
- Dictionary builds so far: 3 (debugging of cap scope, access class, JSON size; all with the same results/data_status.json sha256, no association results existed). `build_count` is stored in dictionary_v1.json from the next build on (next = 4: ALLELE/CLIN, DEV-DICT-1, SPEC_MMR exclusion), to run once the PI declares the CGT rule files final.
