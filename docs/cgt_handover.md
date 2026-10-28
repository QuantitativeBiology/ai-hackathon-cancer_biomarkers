# CGT build — handover from the first (stopped) agent

The first CGT agent read the spec and the sources but wrote no deliverables (its run was halted mid-query). No CRISPR
effects or predictions were read. Its findings, kept for the resuming agent:

## Sources and record-ID fields
- CIViC `external/clinical/civic/nightly-ClinicalEvidenceSummaries.tsv`: evidence_id (EID), evidence_level A–E,
  evidence_direction, significance, nct_ids, doid, therapies, evidence_status; 2,912 Predictive.
- CIViC `nightly-AssertionSummaries.tsv`: assertion_id (AID), amp_category, nccn_guideline, regulatory_approval,
  fda_companion_test; 30 Predictive (AID 85 = KIT D816V / imatinib / Resistance → CL1a).
- CIViC `nightly-VariantSummaries.tsv` (fusion partners), `nightly-MolecularProfileSummaries.tsv`.
- MOAlmanac `external/clinical/moalmanac/.../molecular-oncology-almanac.json`: 1,018 assertions (assertion_id,
  feature_type, oncotree_code, therapy_name, therapy_sensitivity/resistance, predictive_implication, pmid, nct,
  _deprecated).
- FDA `external/clinical/fda_cdx/fda_cdx_parsed.csv` (232 rows; row_index + pma_ids) and
  `fda_cdx_group_indications_parsed.csv` (9 group rows).
- OncoKB gene lists (gene level only).

## Record IDs already located
- BRAF V600E: CIViC EID 1409/1411/7287/3017/9851/13309/12161/12162, AID 7/10/11/13; FDA rows 39/40/72/159/160/161,
  CRC 6/156/184/204, NSCLC 66/85/123, ATC 188, LGG 213; MOA 144–161, 811, 865, 911.
- EGFR activating: CIViC EID 2994–2997/11219/11227/11228/12718/1592/1867, AID 5/6/105/130; FDA rows 41–51, 59–62,
  75, 98/99, 125/127; MOA 272–275, 947–950. EGFR amp: only MOA 713 (clinical evidence, NSCLC) → consistent with
  removing EGFR amp from the sensitivity partition (CL1d).
- PTEN/AKT1/PIK3CA → capivasertib: CIViC EID 12020, 3039, 12034; MOA 937–942; FDA rows 87, 14 → target_set
  {AKT1, AKT2, AKT3} (CL1e).
- Alpelisib/inavolisib PIK3CA: FDA rows 74/97/158/202, 13/198; MOA 1000.
- Refuted: KRAS-mut NSCLC → MEK: CIViC EID 2998 (SELECT-1). CDKN2A loss → CDK4/6: CIViC EID 7444, 11653, 11654, 7299,
  7300.
- MTAP deletion: CIViC EID 12158 (anvumetostat). MAT2A/IDE397 and WRN/HRO761 have no CIViC/MOA/FDA record → cached
  ClinicalTrials.gov records NCT05245500, NCT05094336, NCT04794699, NCT05838768, NCT06004245.
- ESR1: CIViC EID 11443, 13394; MOA 864 (ESR1 LBD = N vs ER+ WT; 1 dev carrier).
- RB1 → CDK4/6 (bypass-resistance W): CIViC EID 1535/1595/1630; MOA 824/892.
- KIT D816V: AID 85 (imatinib resistance); the avapritinib/midostaurin S rows still need record IDs.

## Outstanding
ER+/CDK4-6 IDs (probably the FDA group table); RAS-mut → EGFR in CRC merge (CL1f); MSI → pembrolizumab (MOA 768/769)
vs MSI → WRN split; fusion tests (ALK/ROS1/RET/NTRK/FGFR2/FGFR3/BCR-ABL1/KMT2A/NRG1); BRCA1/2 → PARP1 (N,
discordant mechanism); ADC/ICI non-cell-autonomous rows; del(17p) → BCL2 (population-only N, CL1b); remaining refuted
rows (BRAF V600E CRC monotherapy; EGFR-amp GBM, which has no CIViC/MOA/FDA record).

## Practical notes
- Use `python -X utf8` (configs are UTF-8; cp1252 decoding fails). pandas fails under `python -I` (dateutil missing).
- Query one source at a time with narrow column output.
- PI decision D5: a `clinically-refuted` row may cite a PMID as its only source record if the PMID is verified to exist
  via the PubMed E-utilities API (no email or identifiers sent); such rows carry `source_type = pmid_only`.
