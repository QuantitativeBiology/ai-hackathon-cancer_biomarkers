# CGT implementation notes (M6, DESIGN v3 §10)

Files: `configs/cgt_curated_mapping.csv` (= DESIGN's `cgt_mapping.csv`; one row per test), `configs/cgt_allele_rules.json`
(machine rules for every ALLELE/CLIN feature), `configs/cgt_extra_features.csv` (table view of the same rules),
`configs/cgt_corrections.csv` (CL1a–h, each with the test ids and the assertion the build checks),
`scripts/fetch_cgt_records.py` (CTG v2 + PubMed esummary cache), `scripts/build_cgt.py`, `tests/test_cgt_table.py`.
No CRISPR effects, predictions or P1 metrics were read at any point.

## Step 1: curation decisions (one line each)

- Record id syntax: `CIViC:EID<n>` / `CIViC:AID<n>`, `MOA:<assertion_id>` (not list index; the handover's "MOA 144–161" were checked against assertion_id), `FDA:<row_index>|<PMA>` (both must match the parsed scrape), `CTG:<NCT>`, `PMID:<n>`.
- Tier is assigned from the frozen sources only (T1 needs an FDA CDx row or MOA FDA-Approved record naming the label biomarker); approvals absent from the sources lower the tier (ER+→ESR1, ER+→CDK4/6, KIT D816V→midostaurin, KIT exon 11→imatinib are T2).
- KIT D816V (CL1a): S test with midostaurin only (CIViC EID1725); avapritinib has no CIViC/MOA/FDA record → not listed; imatinib row is acquired-resistance (AID85, FDA H140006).
- Acquired-resistance rows (KIT D816V/imatinib, EGFR C797S, ABL1 T315I) use comparator `listed-only`, sign N, ko_class `untestable` (keeps (feature, T, comparator) unique vs the KIT S row).
- del(17p) (CL1b): `population-only`, N, operationalised as TP53 hemizygous loss in CLL (untestable in practice).
- BRAF (CL1c): V600 allele feature; the tumour-agnostic test excludes CRC (label exclusion, `comparator_exclude=IND:COADREAD`); BRAF V600E CRC is its own `concordant-KO/combination-required` test.
- BRAF V600E CRC monotherapy (§10.6 refuted item) has the same (feature, target, comparator) as the approved encorafenib+cetuximab test; kept as one T1 test carrying the refuting record EID1405 (Kopetz 2015, PMID 26460303 cached) → refuted set = 3 rows (≥ 3 required).
- EGFR (CL1d): one test on ALLELE:EGFR_activating (L858R, ex19 in-frame del, G719X, L861Q, S768I, ex20 ins, T790M); ALLELE:EGFR_L858R_EX19DEL kept as `feature_secondary`; T790M not a separate test (would double count).
- AMP:EGFR appears only in the clinically-refuted GBM row (DESIGN §10.6 lists it; CL1d's "no AMP:EGFR row" is enforced for all non-refuted partitions).
- EGFR-amp GBM refuted row: no refuting CIViC/MOA/FDA record (CIViC EID772/1128 support EGFRvIII→TKI) → PI D5 `pmid_only` with PMID 19204207 (EORTC 26034), verified by esummary.
- PTEN/capivasertib (CL1e): breast test on the label composite CLIN:CAPI_ALT (PIK3CA/AKT1/PTEN), T = AKT1–3; separate PTEN_LOF→AKT test for mHSPC (PTEN IHC label, genomic biallelic loss as approximation).
- RAS→EGFR (CL1f): one W test, CLIN:RAS_MUT (KRAS/NRAS codons 12/13/59/61/117/146) within CRC, BRAF V600E lines excluded so the reference is RAS/BRAF-WT CRC.
- ESR1 LBD (CL1g): N vs ER+ WT (comparator_restrict CLIN:ER_POS); testability computed, not hard-coded.
- 9p21 (CL1h): CDKN2A±MTAP merged into CLIN:MTAP_LOSS (DEL ∨ SIL) tests for PRMT5 (T2, CIViC EID12158 B) and MAT2A (T3, CTG only); DEL:MTAP as `feature_secondary`.
- KRAS G12C (CL6): S test vs KRAS-WT (other KRAS alleles excluded from the comparator), plus an N test G12C vs other KRAS-mutant lines (class `concordant`, never `discordant-mechanism`).
- Families pre-registered: MEK {MAP2K1, MAP2K2}; AKT {AKT1–3}; CDK4/6 {CDK4, CDK6}; NTRK {1–3}; FGFR {FGFR1–4} for both pemigatinib/futibatinib and erdafitinib (one family definition); PARP {PARP1, PARP2}; NRG1 fusion {ERBB2, ERBB3}; ICI {PDCD1} (MSI, TMB) and {PDCD1, CD274} (PD-L1 IHC).
- Indications sharing feature/target/comparator are merged into one test (ALK: NSCLC ∪ ALCL ∪ IMT; HER2+: breast ∪ GEA ∪ CRC ∪ biliary; BRCA: ovary ∪ breast ∪ pancreas ∪ prostate).
- ER+ → ESR1: the only ER+-selection record in the frozen sources is CIViC EID1541 (PALOMA-1, letrozole partner); used, tier T2, flagged unresolved (no endocrine-monotherapy record).
- ER+ → CDK4/6 uses CLIN:HR_POS_HER2_NEG within breast (EID1541); ribociclib/abemaciclib approvals are not in the sources.
- RB1 → CDK4/6 (W): pan-cancer + lineage comparator because the records span breast, GBM, prostate and osteosarcoma (EID12045 B; EID1535/1595/1630 D; MOA 824/892 preclinical); tier T2 from EID12045.
- KIT exon 11 GIST: CIViC A (EID654) + MOA imatinib GIST label (CD117+, not genotype) → T2.
- CDKN2A/CCND1 → CDK4/6 refuted: union feature CLIN:CDK46_SEL; records are basket/TAPUR "Does Not Support" items outside breast → pan-cancer comparator (PALOMA-1 cohort 2 is the rationale only, not a cited record).
- KRAS-mut NSCLC → MEK refuted: ALLELE:KRAS_MUT within NSCLC, CIViC EID2998 (SELECT-1, level A) + EID4867.
- Forward-looking T3 rows (CL8): CTG records cached; CCNE1 (NCT04855656), SMARCA4 (NCT05639751), BRCA→POLQ (NCT04991480), MSI→WRN (NCT05838768, NCT06004245), MTAP (NCT05245500, NCT05094336, NCT04794699, NCT03435250) match their regex.
- ATM→ATR (NCT04497116, eligibility says only "ATR inhibitor-sensitizing mutations") and BRCA→USP1 (NCT05240898, no BRCA/HRD text) fail the regex; kept in the curated mapping so the build drops them visibly to `results/clinical/cgt_dropped.csv` (DESIGN §10.3 requires a matching CTG record; MOA 123–125 elimusertib is not used to rescue ATM→ATR).
- Non-cell-autonomous / IHC rows (HER2-low T-DXd, MET IHC ADC, PD-L1, FOLR1, CLDN18, TMB-H) are N with features `not_assayable` → no `feature_D`, untestable; MSI→PD-1 is N with SIG:MSI.
- IDH1/IDH2 → own gene: N, `discordant-mechanism`; BRCA → PARP1/2: N, `discordant-mechanism`; HRR-panel and HRD-score labels not added as separate tests.
- Proxy rule (DESIGN v3 §10.3): labelled `__proxyHOT` copies (proxy = 1) for BRAF V600 (×3), EGFR_activating, PIK3CA_ACT, EZH2_GOF, IDH1_R132; none for KRAS G12C, EGFR T790M/C797S, KIT D816V/exon 11, ERBB2 TKD.
- in_gate = 1 when the test's target and feature overlap a §9.2 gate control (BRAF, EGFR, KRAS G12C, PIK3CA, ERBB2/HER2+, ER+/ESR1, RB1→CDK4/6, MSI→WRN, MTAP→PRMT5/MAT2A, CCNE1→PKMYT1, SMARCA4→SMARCA2).
- Fusions: P1 fusion cache, default entries, Gene1 = 5′ / Gene2 = 3′ (verified on NCI-H3122 EML4→ALK, K562 BCR→ABL1, RT112 FGFR3→TACC3), FFPM ≥ 0.1, partner ≠ driver, partner ∈ CIViC VariantSummaries canonical partners; curator additions: ROS1 + SLC34A2, EZR (HCC78); RET + NCOA4; NRG1 + TENM4 (MDA-MB-175-VII, `fusion_known_positives.csv`); KMT2A (PI ruling 2): mechanistic rule, any fusion call with KMT2A as 5′ gene, KMT2A::KMT2A excluded (`partners: "*"`; replaces the memory-based list).
- MET exon 14: splice-site variants at c.2888-1/-2 and c.3028+1/+2 or deletions inside c.2888–3028 (intronic −18/−20 variants excluded).
- FLT3 ITD: in-frame ins/dup codons 572–630 from short-read calls (under-detection noted).
- ER_POS: ModelSubtypeFeatures tokens first, else 2-component GMM on ESR1 log2(TPM+1) fitted on dev lines of Breast/Uterus/Ovary; HER2_POS: CN > 2 ∨ "HER2+"; MTAP_LOSS: DEL ∨ SIL; MSI via `biomark.molecular.msi_status` (MSIScore ≥ 20).
- Feature ids: dictionary families use `FAMILY:<symbol>` in the mapping and `FAMILY:<entrez>` in `feature_D`; ALLELE/CLIN ids are their own `feature_D`; `not_assayable` features have no `feature_D`.

## Step 2: build decisions (one line each)

- `scripts/build_cgt.py` resolves 526 record tokens (524 resolved; the 2 unresolved are the CTG records of the dropped ATM→ATR and BRCA→USP1 rows); a soft gene-consistency check (feature gene or target named in the record text) passes for every resolved CIViC/MOA/FDA id.
- Drop rules: no resolvable primary id; T3 without a regex-matching cached CTG record; `pmid_only` outside clinically-refuted. PMIDs on non-pmid_only rows are supplementary and never count as the resolving id. Dropped rows → `results/clinical/cgt_dropped.csv`.
- Drugs: lower-cased, de-duplicated, sorted per test (`n_drugs`).
- Population: P1 dev (D) and test (R) lines; carriers 1/0 on profiled lines, NaN on unprofiled; mutation calls = OmicsSomaticMutations default entries (cached subset `cache/cgt/`), CN/damaging/expression/hotspot/fusion = P1 raw caches, MSI via `biomark.molecular.msi_status`.
- Indication → Model.csv mapping lives in `configs/indication_map.csv` (PI ruling 4; one row per indication code × OncoTree match value, with match_field, source and citation); build_cgt.py reads it and fails on unknown codes; `SOLID` = `not_OncotreeLineage` Myeloid/Lymphoid; pan-cancer tests with `SOLID` are restricted to solid lineages.
- Testable (P2) = feature_D present ∧ ≥ 1 target in the CRISPR library (gene_map `in_matrices` metadata, no effect values) ∧ ≥ 5 D carriers ∧ ≥ 3 D carriers and ≥ 3 D non-carriers in the comparator; `p2_scored` additionally requires partition sensitivity/bypass and a concordant class; `ceiling_testable` = testable ∧ sign S/W.
- ER mixture: 126 dev Breast/Uterus/Ovary lines, component means 0.07 and 1.84 log2(TPM+1), seed = seeds.json `seeds[0]`; 23/54 dev breast lines called ER+ (higher than the ~13 subtype-annotated; GMM may over-call — check before freeze).
- Known fusion positives (NCI-H3122, NCI-H2228, K562, MV4-11, MDA-MB-175-VII) are all carriers (H2228 and MDA-MB-175-VII are outside D/R; evaluated directly).
- DEVIATION DEV-CGT-1 (PI ruling 1): 11 D ALK carriers > pre-data cap of 10; all canonical (7 NPM1::ALK ALCL + 4 EML4::ALK, FFPM ≥ 1.9); kept, not tuned; logged in `summary.json` `deviations`; the test now asserts every ALK carrier's partner is in the frozen CLIN:ALK_FUS partner list.
- Proxy copies get their own `feature_D` (HOT:<entrez>); the summary's headline counts use proxy = 0 only.
- KMT2A recount (PI ruling 2): 19 D / 2 R carriers (was 11 D); 12 D carriers in the acute-leukaemia comparator (adds EOL-1, whose call is KMT2A::RP11-770J1.8; EOL-1 is a known KMT2A-PTD line); 8 of the 19 are solid-tumour calls with non-canonical partners (outside the comparator; they do enter the pan-D feature count).
- ER concordance (PI ruling 3, flag only, rule unchanged): breast D+R lines n = 65, ER_POS 29. vs PAM50 luminal (argmax ∈ {LumA, LumB}; 21 luminal): 15/14/6/30 (ER+lum / ER+non-lum / ER−lum / ER−non-lum), agreement 0.69; GMM alone 0.54. vs ModelSubtypeFeatures exact 'ER+' (3 lines): all 3 ER_POS, but 26 ER_POS lines lack the token, agreement 0.60. Flag: the mixture over-calls ER+ relative to PAM50 (14 ER_POS non-luminal lines).
- Recorded gaps (PI ruling 5, `summary.json` `recorded_gaps`): ER+ endocrine-therapy record, avapritinib (KIT D816V SM), ribociclib/abemaciclib approvals.
- PI ruling D9 (operationalisation before any dependency join): CLIN:ER_POS for Breast lines = pre-declared expression rule AND PAM50 luminal (argmax of OmicsGlobalSignatures Breast_PAM50_* ∈ {LumA, LumB}; breast lines without PAM50 = NA, none occur); non-breast lines keep the expression rule only (flagged; 29 such D+R lines, no CGT test uses ER+ outside breast). CLIN:HR_POS_HER2_NEG inherits it. Breast ER+ falls from 23 D / 6 R (pre-D9) to 11 D / 4 R; HR+/HER2− breast is 6 D / 3 R. Tests: ER+→ESR1 11 comparator carriers (testable); HR+/HER2−→CDK4/6 6 (testable); CAPI_ALT in HR+/HER2− 3 vs 3 (testable, minimal); PIK3CA_ACT in HR+/HER2− 2 (untestable; HOT proxy 3, testable); ESR1 LBD in ER+ 1 (untestable).

## Pre-freeze audit rulings (reviews/cgt_audit.md; PI: all items accepted; molecular-only operationalisation, no dependency data)

- AUD-1 (critical) KMT2A_R: KMT2A 5′ fusion, partner ≠ KMT2A, excluding partners within ±2 Mb of KMT2A (Ensembl 116 coordinates via `biomark.genesets.gene_coords`) and lncRNA-type partners (RP11-/LINC/AC…/-AS, or no Ensembl id); no partner list. Carriers 19 D / 2 R → 12 D / 0 R; EOL-1 (KMT2A-PTD, lncRNA call) and the six 11q23.3 neighbour calls drop out; ONS76 (KMT2A::WNT3, CNS) remains, outside the leukaemia comparator (11 comparator carriers).
- AUD-2 (critical) CDK46_SEL refuted: EID7444 (Does Not Support Resistance) removed; build now rejects any CIViC EID on a clinically-refuted row unless it is Does Not Support + Sensitivity/Response or Supports + Resistance (`REFUTING`); comparator excludes Breast and LOFBI:RB1 lines (388 carriers / 408 non-carriers in D).
- AUD-3 (major) ER_POS→ESR1: CIViC/MOAlmanac/FDA searched for an ER+ fulvestrant/elacestrant/tamoxifen/SERM record: none (only PGR/TFF3/AGR2/TP53 modifiers and ESR1-mutant SERD records). Row flagged unsupported (`headline_exclude`), p2_scored = 0; letrozole dropped (drugs empty); EID1541 kept only as the ER+ population record.
- AUD-4 (major) CLIN:ER_POS := Breast ∧ PAM50 luminal ∧ ESR1 log2(TPM+1) ≥ 3.0 (fixed); non-breast NA; HR+/HER2− := ER_POS ∧ ¬HER2+. Sanity (reported, cut not tuned): MCF7 5.62, T47D 5.64, ZR-75-1 5.04 → ER+; MDA-MB-231 0.14, BT-549 0.02 → ER−. Breast ER+ 11 D / 4 R (D9) → 6 D / 4 R; HR+/HER2− 6 D / 3 R → 3 D / 3 R. Supersedes the D9 rule and the GMM (no longer used).
- AUD-5 (major) RB1→CDK4/6: no RB1-only A/B record (CIViC RB1 A–C profiles are phosphorylation/alpelisib, SCLC chemotherapy, RB1-WT∧CDKN2A composite); demoted from T2 to tier `preclinical` (outside the clinical_s10 tier list; PI to confirm the label); in_gate = 1 stated in notes.
- AUD-6 (major) FGFR: T = {FGFR3} bladder (FGFR3_ALT__FGFR3__BLCA), {FGFR2} cholangiocarcinoma (FGFR2_FUS__FGFR2__CHOL), class `concordant`; drug lists kept.
- AUD-7 (major) merged-indication tests (≥ 2 indication codes, within-indication): comparator suffix `:stratified` (indication fixed effect), `indication_fixed_effect = 1`, per-indication D carrier/non-carrier counts in `indication_strata_D`; per-indication β to be produced by the ceiling (§10.4). Carrier-defined OncoTree codes recorded in the comparator as `cross-lineage=ALCLALKP` (ALK) and `cross-lineage=CMLBCRABL1` (BCR::ABL1); the build verifies all D lines of those codes are carriers.
- AUD-8 (minor) T790M removed from ALLELE:EGFR_activating (CL1d text updated); the HOT:EGFR proxy no longer approximates T790M.
- AUD-9 (minor) EGFR ex19: frameshift deletions in exon 19 of one line within 10 bp whose summed length is a multiple of 3 are merged as in-frame: ACH-002529 (2+16 bp), ACH-000150 (2+16), ACH-002650 (2+2+8, dev); EGFR_activating D carriers 7 → 8 (NSCLC comparator 4 → 5).
- AUD-10 (minor) CAPI_ALT uses ALLELE:PIK3CA_CAPI = PIK3CA missense flagged Hotspot (OmicsSomaticMutations) or R88Q (CIViC EID12034); FDA row 87 gives no allele list. CAPI_ALT D carriers 180 → 197 (HR+/HER2− comparator 2 carriers / 1 non-carrier after AUD-4: untestable).
- AUD-11 (minor) BRAF_V600 adds V600_K601delinsE/D/K (no new carriers).
- AUD-12 (minor) ESR1_LBD restricted to SERD-label alleles (K303R, E380Q, V422del, S463P, L469X, L536X, Y537X, D538G, L549P): 4 → 1 D carrier.
- AUD-13 (minor) BRAF_V600→MEK secondary reference excluding NRAS-mutant and LOFBI:NF1 lines: 38 carriers vs 4 non-carriers (primary 38 vs 14).
- AUD-14 (minor) EID1405 / PMID 26460303 moved from source ids to notes of BRAF_V600E__BRAF__CRC; the refuted set is 3 rows by construction.
- AUD-15 (minor) De Novo (DEN) and 513(f)(2)/510(k)/HDE FDA ids accepted; resolution table and summary report `fda_kind` (177 PMA, 1 De Novo, 1 510(k), 1 HDE).
- AUD-16 (minor) ABL1 T315I listed drugs trimmed to imatinib, dasatinib (records EID234, 4365/4368, MOA 69).
- AUD-17 (minor) MTAP→PRMT5/MAT2A: pre-registered expectation of a data-null ceiling (drug-specific MTA-cooperative mechanism; PRMT5/MAT2A KO broadly essential); §8.4 MTAP-specific resolution test is part of DESIGN §8.4 (not built here).
- AUD-18 (minor) EGFR-amp GBM citation states that EORTC 26034 was not amplification-selected.
- DESIGN text items proposed by the audit (not edited here; for the DESIGN owner): §10.1 add De Novo/513(f)(2) to allowed FDA ids; §10.6 refuted set = 3 rows by construction; §10.3 avapritinib recorded gap; `preclinical` tier label.
- Concurrent edit adopted: ALLELE:EGFR_ACT was renamed ALLELE:EGFR_activating by another agent (DESIGN §3 / positive controls); mapping generator and tests follow the new name; `test_cgt_table` EGFR assertion fixed to match test id `EGFR_ACT__EGFR__NSCLC` (the renamed check had matched no rows).
- Totals after audit: 64 tests (57 headline), 27 testable, 19 P2-scored headline (8 gate-disjoint; ≥ 15 minimum still met).
- AUD-5 follow-up (PI): tier `preclinical` confirmed for RB1→CDK4/6; reported as a secondary and excluded from the P2 headline (build rule `P2_HEADLINE_TIERS = {T1, T2, T3}`; `test_p2_headline_tiers`).
