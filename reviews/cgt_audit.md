# CGT pre-freeze audit (independent clinical-oncology review, 2026-10-28)

## Scope

**Read:** DESIGN §10, `clinical_s10`, `configs/clinical_ground_truth.csv` (current on-disk version, 64 rows), the curated mapping, allele rules, corrections, `indication_map.csv`, the implementation notes, `summary.json` and the record-resolution table.

**Records opened:** CIViC evidence and assertion rows; MOAlmanac by `assertion_id`; FDA CDx by `row_index`; the CTG and PubMed caches.

**Allele calls:** grepped in `OmicsSomaticMutations.csv` (default entries only). No CRISPR effects or predictions were read.

## Coverage

- All 28 testable sensitivity/bypass rows (22 headline + 6 proxy).
- The 3 clinically-refuted rows.
- 10 of the remaining 33 rows, drawn with `random.Random(20261028).sample` over the sorted test_ids: BRCA_BI__PARP, IDH1_R132__proxyHOT, BRAF_FUS__LGG, KIT_D816V__SM, KIT_EX11__GIST, DEL17P__BCL2, PDL1_IHC, FGFR2_FUS__CHOL, ABL1_T315I__listed, IDH2_MUT.

## Verdict

Record support is good overall. For the BRAF, EGFR, HER2, ALK, BCR::ABL1, KRAS G12C, RAS→EGFR, capivasertib, FGFR3, KMT2A, EZH2, MTAP, CCNE1 and sampled rows, the cited records name the claimed biomarker, drug, indication and direction.

The allele rules behave correctly on real calls:
- BRAF class II calls are excluded (K601N, L485W, N486_P490del).
- EGFR exon-19 in-frame deletions are captured and frameshifts excluded; the Exon field reads `19/28`.
- The KRAS codon set is correct.
- The PIK3CA set matches the alpelisib label.
- MET c.3028+1 is captured, and intronic c.2888-18/-20 SNVs are excluded.
- NPM1 W288fs is captured.

CL1a–h are applied as their assertions state.

**Before `freeze-cgt`, fix:** the 2 critical findings and the 5 major ones.

## Findings

| row / test id | severity | problem | proposed fix |
|---|---|---|---|
| KMT2A_R__MEN1__acute_leuk (CLIN:KMT2A_R) | **critical** | "Any partner, KMT2A 5′" admits 8 solid-tumour calls. 6 have 11q23.3 neighbours as partner, which points to read-through or local artefacts: ATP5MG (HCC38), CD3G (TEN), UBE4A (EFE184), SCN2B (HCC95), PHLDB1 (CAOV3), TMPRSS4 (ESO26). The others are WNT3 (ONS76) and lncRNAs (C32, NCI-H1568). All 8 enter the pan-D carrier count that the arms rank. EOL-1 (KMT2A::RP11-770J1.8) is a KMT2A-PTD line. The revumenib label (FDA row 226, MOA 1001/1002) is "KMT2A translocation" by break-apart FISH, so EOL-1 is not a label carrier. | Exclude partners within ±2 Mb of KMT2A and RP11/lncRNA partners, or require a canonical partner list (AFF1, MLLT3, MLLT1, MLLT10, AFDN, ELL, EPS15, SEPTIN6/9…). Exclude EOL-1. Recount and log the change as a pre-freeze rule change. |
| CDK46_SEL__CDK46__refuted | **critical** | CIViC EID7444 (CDKN2A loss, NSCLC, palbociclib) is "Does Not Support **Resistance**". That is the wrong direction: it does not refute sensitivity. EID7299/7300/11653/11654 ("Does Not Support Sensitivity") and EID11673 (CCND1-3 amp, TAPUR) are valid. | Remove EID7444. Add a build check that refuted rows cite only "Does Not Support Sensitivity" or "Supports Resistance" records. |
| ER_POS__ESR1__BRCA | major | Wrong drug→target. The only record, EID1541 ("ESR1 Overexpression", palbociclib+letrozole, PALOMA-1), supports the CDK4/6 row. Letrozole inhibits aromatase (CYP19A1). It works by systemic oestrogen deprivation, which is non-cell-autonomous in culture, not by blocking ESR1. ER+→ESR1 KO is mechanistically right (SERD/SERM), but no cited record supports it. | Re-source to a fulvestrant or elacestrant ER+ record if MOAlmanac has one. Otherwise flag it as unsupported and exclude it from the headline P2 set. Drop letrozole from `drugs`. |
| RB1_LOF__CDK46__pan | major | T2 rests on EID12045, a **composite** profile ("PIK3CA OR RB1 OR ESR1 Mutation", PALOMA-3 ctDNA) that does not isolate RB1. The remaining records (EID1535/1595/1630, MOA 824/892) are level D or preclinical. The W sign and family class are correct. | Cite an RB1-only A/B record or demote the row out of T2. The row is in_gate = 1, so state this. |
| FGFR3_ALT__FGFR__BLCA; FGFR2_FUS__FGFR__CHOL | major | T = {FGFR1–4} copies the drugs' pan-FGFR selectivity. For a single-gene KO in oncogene addiction, the dependency is the altered gene. A 4-gene family credits spurious FGFR1/4 hits and costs Bonferroni ×4. This differs from CAPI_ALT, where an upstream biomarker names no AKT isoform. | Set T = {FGFR3} for the bladder row and {FGFR2} for the cholangiocarcinoma row, class `concordant`. Keep the drug lists. |
| ALK_FUS__ALK__NSCLC_ALCL_IMT; BCR_ABL1__ABL1__CML_ALL; HER2_POS (4 indications) | major | When indications are merged into one "within-indication" comparator, the contrast is confounded with lineage. 7 of the 11 ALK carriers are ALCL (Lymphoid), while the non-carriers are mostly NSCLC. The `ALCLALKP` and `CMLBCRABL1` OncoTree codes are defined by the biomarker itself, so they contain no non-carriers. | Use an indication fixed effect (stratified contrast) for merged tests and report per-indication β. Treat carrier-defined codes as cross-lineage contrasts. |
| CLIN:ER_POS / CLIN:HR_POS_HER2_NEG | major | The GMM component means are 0.07 and 1.84 log2(TPM+1). ER+ breast lines are typically above 4, so the mixture splits "off" from "low". The D9 PAM50 gate protects breast lines only. ER_POS__CDK46 has 26 D carriers but only 6 in breast, so most of the carriers the arms rank are uterus or ovary lines called by the mixture alone. | Verify the scale (ESR1 in MCF7 and T47D). Fit on breast only or use a pre-declared fixed cut. Set non-breast lines to NA. Re-run before freeze. This is a molecular-only rule, so changing it is not tuning against effects. |
| CDK46_SEL__CDK46__refuted (comparator) | major | The pan-cancer comparator includes breast, where CDK4/6 inhibition is approved, not refuted. CDKN2A loss and RB1 loss are mutually exclusive, so the non-carriers are enriched for RB1-null lines. That inflates S for a reason unrelated to CDKN2A as a selector. | Exclude Breast (the records are all outside breast). Add LOFBI:RB1 to `comparator_exclude`. |
| ALLELE:EGFR_ACT (+ proxyHOT) | minor | T790M is in the "activating" union. It adds no carriers (it co-occurs with a sensitiser, e.g. NCI-H1975 L858R+T790M). Its inclusion means the HOT proxy approximates an allele that §10.3 says is never proxied. | Drop T790M from EGFR_ACT and update the CL1d text. |
| ALLELE:EGFR_EX19DEL | minor | ACH-002529 has two adjacent frameshifts (E746VfsTer16 and L747FfsTer14, AF ≈ 0.98). This may be a complex in-frame indel split into two calls. | Inspect DNAChange. If the net change is in-frame, call it ex19del. |
| CLIN:CAPI_ALT (PIK3CA part) | minor | Uses the alpelisib allele set. The capivasertib CDx (FDA row 87) is broader, and real calls of Q546K ×3, G1049R ×4, N345K, E726K and M1043I are missed. | Give CAPI_ALT a capivasertib-specific PIK3CA set. |
| ALLELE:BRAF_V600 | minor | The rule misses p.V600_K601delinsE, which is class I. | Add `^p\.V600_K601delins[EDK]$`. |
| ALLELE:ESR1_LBD | minor | Any missense in codons 303–549 admits VUS. Impact is low because the row is N and untestable. | Restrict to the label alleles (E380Q, S463P, L536X, Y537X, D538G…). |
| BRAF_V600__MEK__MEL | minor | The BRAF-WT melanoma reference includes NRAS- and NF1-mutant lines, which are also MEK-dependent. This is still the correct clinical contrast. | Keep it. Add a secondary analysis with an NRAS/NF1-excluded reference. |
| BRAF_V600E__BRAF__CRC | minor | The §10.6 monotherapy-refuted claim (EID1405, PMID 26460303) is folded into the T1 combination test. The class is right, because the EGFR feedback that defeats monotherapy is drug-specific and the KO is still S. As a result the refuted set has 3 rows. | Move EID1405 to `notes`. State in §10.6 that the refuted set has 3 rows by construction. |
| KMT2A_R row source | minor | FDA row 226 is a De Novo grant (DEN240067), which is not in the §10.1 list of allowed source types. The row also resolves via MOA. | Add De Novo/513(f)(2) to §10.1. |
| ABL1_T315I__listed | minor | Bosutinib and nilotinib are listed, but the opened records cover only imatinib and dasatinib. | Cite a record for each drug or trim the list. |
| MTAP_LOSS__PRMT5/MAT2A | minor | MTA-cooperative PRMT5 inhibitors work through a drug-specific mechanism, and PRMT5 KO is near-pan-lethal, so a data-null ceiling is likely. | Pre-register that expectation. Confirm the CL1h §8.4 resolution test exists. |

## Rows confirmed correct (no change)

- **N and discordant rows:** BRCA→PARP1/2 N `discordant-mechanism` (PARP1 loss confers resistance to trapping). IDH1/IDH2 N. ESR1 LBD N within ER+. DEL17P N `population-only`. IHC, ADC and immune-checkpoint rows N, non-cell-autonomous.
- **KRAS and RAS rows:** KRAS G12C S vs KRAS-WT, plus an N test vs other KRAS-mutant lines (`concordant`). RAS→EGFR W vs a RAS/BRAF-WT CRC reference.
- **Families:** pan-AKT → AKT1–3, for both capivasertib and PTEN. CDK4/6 pair, any member (CDK4 in breast, CDK6 in haematological lines). MEK1/2.
- **Synthetic-lethal and menin rows:** MEN1 for KMT2A-r and NPM1c. POLQ for BRCA (NCT04991480). PKMYT1 for CCNE1-amp (NCT04855656, EID12761/2). SMARCA2 for SMARCA4. WRN for MSI.
- **Dropped rows:** ATM→ATR and BRCA→USP1 were dropped correctly.
- **Refuted rows:** KRAS-mut NSCLC→MEK (EID2998 SELECT-1 "Supports Resistance" A; EID4867) and EGFR-amp GBM (PMID 19204207, EORTC 26034) are adequate. Neither trial was amplification-selected; say so in the citation.

## Corrections CL1a–h

All are applied. Residuals:
- **CL1a:** avapritinib is absent (a recorded gap). Note the gap in §10.3.
- **CL1d:** fix the T790M wording (see the EGFR_ACT row above).
