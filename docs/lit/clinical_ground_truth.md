# Clinical ground truth for dependency-subtype validation

Survey date 2026-10-28. `[verified]` = resource/record opened or downloaded in this session; `[from-memory]` = not re-checked here, treat as a claim to confirm before use. Nothing below is fabricated; where a lookup failed I say so.

---

## 1. Resource table

| Resource | URL / download path | Licence | Scriptable? | Usable fields |
|---|---|---|---|---|
| **FDA list of authorized companion diagnostics (CDx)** | `fda.gov/medical-devices/in-vitro-diagnostics/list-fda-authorized-companion-diagnostic-devices-in-vitro-and-imaging-tools` (the older `list-cleared-or-approved-...` slug 301-redirects here) `[verified]` | US Gov public domain `[from-memory]` | **Partly.** HTML table only; **no CSV/XLSX/API link found** `[verified]`. Scrape + parse (~122 kB page, "Content current as of: 09/28/2026") `[verified]` | Diagnostic Name (Manufacturer); Indication and Sample Type; Drug Trade Name (Generic) + NDA/BLA; Biomarker(s); Biomarker details; PMA/510(k)/513(f)(2)/HDE number + approval date `[verified]` |
| **OncoKB** | `oncokb.org`; API `oncokb.org/api/v1` | Free for academic research, licence required for hospital/commercial; token from Account Settings after licence request `[from-memory via cBioPortal docs, verified page]` | **Mixed.** `GET /api/v1/utils/cancerGeneList.txt` → 200, 139 kB, no token `[verified]`; `/utils/allCuratedGenes.txt` → 200, 1.68 MB, no token `[verified]`; **`/utils/allActionableVariants.txt` → HTTP 401** `[verified]`. `oncokb.org/terms` and `/api-access` are JS-only and did not render for WebFetch `[verified]` | cancerGeneList: Hugo Symbol, Entrez, GRCh37/38 isoform+RefSeq, Gene Type (ONCOGENE/TSG), panel membership (MSK-IMPACT, FOUNDATION ONE, Vogelstein, COSMIC CGC v99), aliases `[verified]`. allCuratedGenes adds **Highest Level of Evidence (sensitivity)** and **(resistance)**, Setting, Summary, Background `[verified]`. Per-variant levels 1/2/3A/3B/4, R1/R2 need the licensed endpoint |
| **CIViC** | `civicdb.org`; `https://civicdb.org/downloads/nightly/nightly-<X>.tsv` | **CC0 1.0** (data); MIT (code) `[verified, civic.readthedocs.io FAQ]` | **Yes**, fully. Nightly + historical monthly, TSV and VCF `[verified]`. `/releases/main` and `/downloads/nightly` directory listings are JS/403 but the per-file paths work `[verified]` | `ClinicalEvidenceSummaries.tsv` (4.2 MB, 4,972 rows): molecular_profile(+id), disease, **doid**, phenotypes, therapies, therapy_interaction_type, evidence_type, evidence_direction, evidence_level (A–E), significance, evidence_statement, citation_id/source_type/citation, **nct_ids**, rating, evidence_status, variant_origin, URLs `[verified]`. `AssertionSummaries.tsv` adds **amp_category** (Tier I–IV / Level A–D), **nccn_guideline(+version)**, **regulatory_approval**, **fda_companion_test** `[verified]`. Also `VariantSummaries` (coordinates, HGVS, 5'/3' fusion partners), `FeatureSummaries`, `MolecularProfileSummaries` `[verified]` |
| **Cancer Genome Interpreter — Cancer Biomarkers DB** (Tamborero 2018, *Genome Med* 10:25) | `cancergenomeinterpreter.org/biomarkers` | UPF repository record (DOI 10.34810/data405) lists **CC0** `[from-memory, search result not opened]` | **No stable path found this session.** Site is a JS SPA; `/biomarkers`, `/api/biomarkers`, `/data/biomarkers/cgi_biomarkers_20221017.zip` and `_latest.zip` all return the same 1.4 kB HTML shell `[verified]`. Last dataset update Oct 2022 `[from-memory, bbglab wiki]` — **stale; use as a secondary cross-check only** | Documented classification: gene/alteration, drug, cancer type, drug response (sensitivity / resistance / toxicity), level of clinical evidence `[from-memory]` |
| **ESMO ESCAT** | Mateo et al., *Ann Oncol* 2018;29(9):1895–1902, PMID 30137196 `[verified]`; refreshed in Mosele et al., *Ann Oncol* 2024 ESMO PMWG NGS recommendations, PMID 38834388 `[verified]` | Journal copyright; **no machine-readable file** | **No** — manual tier assignment from the papers | Tiers I-A/I-B/I-C, II, III, IV, V, X; I-A = randomised trial, survival endpoint, in that tumour type `[from-memory for exact wording]`. Useful as the *indication-matched* axis |
| **JAX-CKB / CKB (now Genomenon)** | `ckb.jax.org`, `ckb.genomenon.com`; EULA at `ckbhome.jax.org/legal/eula` | **Restrictive, paid.** CKB BOOST = subscription; licence grants internal use only and excludes rights to download/redistribute; CKB CORE limited to ~50 genes; FLEX = flat-file/JSON for integrators. Became part of Genomenon in May 2024 `[from-memory, search results]` | Only under FLEX subscription | — **exclude from this project** (cannot redistribute a derived benchmark) |
| **Molecular Oncology Almanac (MOAlmanac)** (Reardon 2021, *Nat Cancer*) | DB repo `github.com/vanallenlab/moalmanac-db`; API `api.moalmanac.org` (OpenAPI at `/openapi.json`) | **Database CC BY 4.0** (commercial use allowed with credit), code GPL-2.0; Zenodo concept DOI 10.5281/zenodo.17089309 `[verified]` | **Yes.** Git clone, or REST: `/propositions`, `/statements`, `/biomarkers`, `/therapies`, `/indications`, `/diseases`, `/documents`, `/genes`, `/search`, `/strengths`, `/mappings`, `/contributions`; filters AND across keys, OR within (`?gene=BRAF&gene=KRAS`); last_updated 2026-10-08 `[verified]` | GA4GH-style: `predicate` (e.g. `predictsSensitivityTo`), disease (NCIt id), therapy, `Indication.FDAIndication`/`EMAIndication`, document type "Regulatory approval" with the label PDF URL, biomarker classes incl. `SomaticVariant`, `CopyNumber`, `GeneFusion`, `GermlineVariant`, `MicrosatelliteStability`, `MismatchRepair`, `HomologousRecombination`, `ProteinExpression`, `TumorMutationalBurden`, `WildType` `[verified]` |
| **cBioPortal public API** | `https://www.cbioportal.org/api` (Swagger `/api/swagger-ui/index.html`) | Site content © cBioPortal; per-study terms from the source cohort `[verified — no licence statement on the API docs page]` | **Yes, no auth.** `GET /api/studies` returned 551 studies incl. all 33 `*_tcga_pan_can_atlas_2018` `[verified]` | Study/sample/patient clinical attributes, mutations, CNA, expression — for lineage and biomarker-prevalence denominators |
| **AACR Project GENIE** | `synapse.org/genie`; mirror `genie.cbioportal.org` (Bearer token) | **Registration + click-through terms**; no re-identification, **no redistribution without written permission from the GENIE Coordinating Center**; cite consortium + version. v19.0-public ≈ 271,837 samples / 227,696 patients `[from-memory, AACR FAQ/data guide via search]` | Yes once registered (Synapse client / tokened cBioPortal) | Sample-level mutation/CNA/fusion + OncoTree code, age, sex, sample type — best real-world prevalence denominator |
| **ClinicalTrials.gov API v2** | `https://clinicaltrials.gov/api/v2/studies` | US Gov public domain `[from-memory]` | **Yes, no auth.** `/api/v2/version` → `{"apiVersion":"2.0.5","dataTimestamp":"2026-10-28T09:00:04"}` `[verified]`; `?query.term=...&fields=...` and `/studies/{NCTId}` both work `[verified]` | NCTId, briefTitle, overallStatus, phases, interventions, conditions, **eligibilityCriteria free text** — biomarker-selected enrolment is recoverable by regex (verified for MTAP, MSI/dMMR, KRAS G12C) `[verified]` |

**CIViC content (this session's nightly dump)** `[verified]`: 2,912 accepted Predictive items — significance Sensitivity/Response 1,835, Resistance 1,051, Reduced Sensitivity 16, Adverse Response 10; direction Supports 2,520 / Does Not Support 392; levels D 1,101, B 800, C 790, A 193, E 28; origin Somatic 2,492, Rare Germline 67. Restricting to **level A/B + Supports** gives 849 items, 630 distinct (molecular_profile × therapy) pairs and 174 distinct (profile × therapy × disease) triples at level A — the right size for a benchmark. 668 Predictive rows carry an NCT id; 2,874 carry a DOID. Fusion profiles use `v::ALK Fusion`, `BCR::ABL1 Fusion`, `FGFR2::v Fusion` conventions (30 ALK, 15 BCR::ABL1, 13 NTRK, 11 RET, 7 ROS1 at A/B) `[verified]`. The 30 Predictive **assertions** are near-pure Tier I-Level A with `fda_companion_test` populated `[verified]`.

---

## 2. Mapping biomarker→drug pairs onto (dependency gene, biomarker feature)

A CRISPR screen measures loss of fitness on **knockout of a gene**; a clinical biomarker is **an alteration that predicts response to a drug**. The bridge is the drug's *target gene*, so build the benchmark as a 5-tuple:

`(biomarker_feature, feature_type, drug, target_gene, indication)` → `expected_sign(CRISPR effect of target_gene KO | biomarker+ vs biomarker−)`

Rules used to construct the table:
1. **Target assignment.** Take the molecular target from the drug label / MOAlmanac `statements` (`reportedIn` → Regulatory approval document) rather than from drug-name parsing. Multi-target TKIs get the *indication-relevant* target (crizotinib → ALK in ALK+ NSCLC, ROS1 in ROS1+).
2. **Pathway drugs map to the node, not the biomarker.** BRAF V600E → trametinib tests `MAP2K1/MAP2K2` dependency; ER+ → palbociclib tests `CDK4`/`CCND1`; BCL2 for venetoclax.
3. **Synthetic-lethal pairs invert the roles**: the biomarker gene is *lost*, the target is a *different* gene (MTAP-del → PRMT5/MAT2A; MSI → WRN; BRCA1/2 → PARP1).
4. **Direction.** `stronger` = target KO more lethal in biomarker+ lines (oncogene addiction and synthetic lethality); `weaker` = biomarker predicts *resistance* to the drug (KRAS mut vs anti-EGFR), so target KO should be *less* consequential, or the dependency sits elsewhere; `KO-not-informative` = the drug's effect is not reproduced by gene loss (§3).
5. **Population.** Keep the label indication (and OncoTree/DOID code) so lineage distribution of a recovered subtype can be scored against the approved population, not just against the gene.
6. **Resistance alleles** (ABL1 T315I, EGFR C797S, PARP1 mutations) are **not** testable as dependency direction; keep them in a separate "resistance" partition used only for negative controls.

### Candidate benchmark pairs (45)

Expected direction: **S** = target dependency stronger in biomarker+ lines; **W** = weaker / biomarker predicts resistance; **N** = KO not informative. CDx = on the FDA CDx list as fetched 2026-10-28.

| # | Biomarker | Alt. type | Drug(s) | Target gene | Exp. | Indication(s) | Citation | CDx |
|---|---|---|---|---|---|---|---|---|
|1|BRAF V600E|missense GoF|vemurafenib, dabrafenib|BRAF|S|melanoma, NSCLC, ATC|cobas 4800 BRAF V600 P110020 `[verified]`|yes|
|2|BRAF V600E|missense|encorafenib + cetuximab|BRAF|S|CRC|FoundationOne CDx P200010/S026 `[verified]`|yes|
|3|BRAF V600E|missense|trametinib (± dabrafenib)|MAP2K1|S|melanoma|FDA CDx list `[verified]`|yes|
|4|BRAF fusion/alteration|fusion|tovorafenib (Ojemda)|BRAF|S|paediatric LGG|FDA CDx list `[verified]`|yes|
|5|EGFR ex19del / L858R|indel / missense|osimertinib|EGFR|S|NSCLC|cobas EGFR P150047 `[verified]`|yes|
|6|EGFR T790M|missense|osimertinib|EGFR|S|NSCLC (post-TKI)|cobas EGFR `[verified]`|yes|
|7|EGFR exon-20 ins|indel|sunvozertinib (Zegfrovy)|EGFR|S|NSCLC|FDA CDx list `[verified]`|yes|
|8|ERBB2 amplification|amp|trastuzumab, pertuzumab|ERBB2|S|breast, gastric|HER2 IHC/ISH P090015 `[verified]`|yes|
|9|ERBB2 TKD activating mut|missense/ins|zongertinib, sevabertinib|ERBB2|S|NSCLC|P160045/S049, P240026 `[verified]`|yes|
|10|HER2-low/ultralow (IHC 1+/0+)|protein expr|trastuzumab deruxtecan|ERBB2|N (ADC payload)|breast|PATHWAY 4B5 `[verified]`|yes|
|11|MET protein ≥50% 3+|protein expr|telisotuzumab vedotin|MET|N (ADC)|NSCLC|VENTANA MET SP44 P250002 `[verified]`|yes|
|12|MET exon-14 skipping|splice|capmatinib, tepotinib|MET|S|NSCLC|FoundationOne P170019/S011 `[verified]`|yes|
|13|ALK rearrangement|fusion|crizotinib, alectinib, lorlatinib|ALK|S|NSCLC|Ventana ALK D5F3, Vysis ALK `[verified]`|yes|
|14|ROS1 fusion|fusion|entrectinib, crizotinib|ROS1|S|NSCLC|FoundationOne CDx `[verified]`|yes|
|15|RET fusion|fusion|selpercatinib, pralsetinib|RET|S|NSCLC, thyroid|P230011/S006, P160045/S019 `[verified]`|yes|
|16|NTRK1/2/3 fusion|fusion|larotrectinib, entrectinib|NTRK1/2/3|S|solid tumours (agnostic)|FoundationOne P170019/S017 `[verified]`|yes|
|17|BCR::ABL1|fusion|imatinib, nilotinib|ABL1|S|CML, Ph+ ALL|MRDx BCR-ABL K173492 `[verified]`|yes|
|18|KIT D816V|missense|imatinib|KIT|S|aggressive systemic mastocytosis|H140006 `[verified]`|yes|
|19|KIT exon 11 mut|indel/missense|imatinib|KIT|S|GIST|CIViC A/B (10 items) `[verified]`|yes `[from-memory]`|
|20|PDGFRA D842V|missense|avapritinib|PDGFRA|S|GIST|therascreen PDGFRA P210002 `[verified]`|yes|
|21|FLT3 ITD / TKD D835,I836|indel/missense|midostaurin, gilteritinib, quizartinib|FLT3|S|AML|LeukoStrat CDx P160040(+S002,S011) `[verified]`|yes|
|22|KRAS G12C|missense|sotorasib, adagrasib|KRAS|S|NSCLC, CRC|Guardant360 CDx, therascreen KRAS `[verified]`|yes|
|23|KRAS/NRAS wild-type|absence of mut|cetuximab, panitumumab|EGFR|S (in WT)|CRC|FoundationOne CDx `[verified]`|yes|
|24|KRAS/NRAS mutant|missense|cetuximab, panitumumab|EGFR|W (resistance)|CRC|same label, exclusionary `[verified]`; CIViC 38 Resistance items `[verified]`|yes|
|25|PIK3CA mut (C420R,E542K,E545x,H1047x)|missense|alpelisib|PIK3CA|S|HR+/HER2− breast|therascreen PIK3CA, FoundationOne P170019/S006 `[verified]`|yes|
|26|PIK3CA mut|missense|inavolisib + palbociclib + fulvestrant|PIK3CA|S|HR+/HER2− breast|P170019/S062 `[verified]`|yes|
|27|PTEN protein loss|protein loss|capivasertib|AKT1|S|breast, mHSPC|VENTANA PTEN SP218 P250031 `[verified]`|yes|
|28|AKT1 E17K|missense GoF|capivasertib|AKT1|S|breast|CIViC A/B `[verified]`|yes `[from-memory]`|
|29|ESR1 missense (E380,S463,L536,Y537,D538…)|missense GoF|elacestrant, imlunestrant, camizestrant, vepdegestrant|ESR1|S|HR+ breast (post-AI)|Guardant360 CDx P250027 `[verified]`|yes|
|30|ER+ (ESR1/PGR expression)|protein/expr|tamoxifen, AI, fulvestrant|ESR1|S|breast|IHC standard of care; **not on CDx list** `[verified]`|no|
|31|ER+/HER2−|expr|palbociclib, ribociclib, abemaciclib|CDK4 (and CCND1)|S|breast|approved; **no CDx listed** `[verified absent]`|no|
|32|BRCA1 mut (germline/somatic)|LoF|olaparib, talazoparib, rucaparib|PARP1|N (trapping)|ovarian, breast, pancreas, mCRPC|BRACAnalysis CDx P140020, P140020/S015 `[verified]`|yes|
|33|BRCA2 mut|LoF|olaparib, niraparib|PARP1|N|as above|FoundationOne/BRACAnalysis `[verified]`|yes|
|34|HRD / genomic instability score|signature|niraparib|PARP1|N|ovarian|Myriad myChoice CDx `[verified]`|yes|
|35|MSI-H / dMMR|signature|pembrolizumab, dostarlimab, nivolumab(±ipi)|— (immune)|N (non-cell-autonomous)|solid, EC, CRC|Ventana MMR RxDx P200019; Idylla MSI P250004 `[verified]`|yes|
|36|TMB ≥10 mut/Mb|burden|pembrolizumab|— (immune)|N|solid tumours|FoundationOne P170019/S016 `[verified]`|yes|
|37|PD-L1 CPS/TPS|protein expr|pembrolizumab|CD274|N|multiple|PD-L1 IHC 22C3 pharmDx `[verified]`|yes|
|38|MSI-H / dMMR|signature|HRO761; VVD-133214|**WRN**|S|MSI solid tumours (ph I)|NCT05838768, NCT06004245 `[verified]`; Chan 2019 PMID 30971823, Ferretti 2024 PMID 38658754 `[verified]`|no|
|39|MTAP homozygous deletion|deep deletion|MRTX1719; anvumetostat (AMG 193); BMS-986504; AZD3470; BAY 3713372|**PRMT5**|S|solid tumours (ph I/II)|NCT05245500, NCT05094336, NCT07492680, NCT06130553, NCT06914128 `[verified]`; Kryukov/Mavrakis 2016 PMID 26912360/26912361 `[verified]`|no|
|40|MTAP deletion|deep deletion|IDE397; AG-270 (terminated)|**MAT2A**|S|solid tumours|NCT04794699, NCT03435250 `[verified]`|no|
|41|CDKN2A±MTAP homozygous loss|co-deleted locus|anvumetostat|PRMT5|S (co-location caveat)|solid, HNSCC, thoracic|NCT05094336, NCT06333951 eligibility text `[verified]`|no|
|42|IDH1 R132C/H/G/S/L|neomorphic|ivosidenib, olutasidenib|IDH1|N (neomorph)|AML, MDS, cholangiocarcinoma|Abbott RealTime IDH1 P170041(+S006) `[verified]`|yes|
|43|IDH1 R132 / IDH2 R172|neomorphic|vorasidenib|IDH1/IDH2|N|grade 2 astrocytoma/oligodendroglioma|Oncomine Dx, NDA 218784 `[verified]`|yes|
|44|IDH2 R140/R172|neomorphic|enasidenib|IDH2|N|AML|Abbott RealTime IDH2 P170005 `[verified]`|yes|
|45|FGFR2 fusion/rearrangement|fusion|pemigatinib|FGFR2|S|cholangiocarcinoma|FoundationOne P170019/S013 `[verified]`|yes|
|46|FGFR3 R248C/S249C/G370C/Y373C, fusions|missense/fusion|erdafitinib|FGFR3|S|urothelial|therascreen FGFR RGQ RT-PCR `[verified]`|yes|
|47|KMT2A rearrangement (11q23.3)|fusion|revumenib|MEN1 (menin)|S|acute leukaemia|Cytocell KMT2A FISH DEN240067 `[verified]`|yes|
|48|NPM1 mutation|indel|revumenib, ziftomenib|MEN1|S|R/R AML|FDA approvals Oct/Nov 2025 `[from-memory, search]`; **not on CDx list** `[verified absent]`|no|
|49|EZH2 Y646x / A682G / A692V|missense GoF|tazemetostat|EZH2|S|follicular lymphoma|cobas EZH2, NDA 213400 `[verified]`|yes|
|50|NRG1 fusion|fusion|zenocutuzumab|ERBB3 (ERBB2/3 heterodimer)|S|NSCLC, pancreatic (2024), cholangiocarcinoma (2026)|FDA approvals `[from-memory, search]`; **not on CDx list** `[verified absent]`|no|
|51|del(17p) / TP53 loss|deletion|venetoclax|BCL2|S|CLL|Vysis CLL FISH P150041 `[verified]`|yes|
|52|FOLR1 expression|protein expr|mirvetuximab soravtansine|FOLR1|N (ADC)|ovarian/FT/primary peritoneal|VENTANA FOLR1 P220006 `[verified]`|yes|

(Rows 1–52 with four numbering collisions removed in review ⇒ 45 distinct benchmark pairs; 31 carry an FDA CDx, 7 are trial-stage only.)

---

## 3. KO–drug discordance classes

These are the cases where "subtype did not recover the clinical biomarker" is a **property of CRISPR**, not a failure of the model. They must be pre-registered, not discovered post hoc.

1. **Trapping / poisoning, not loss of function.** PARP inhibitors kill by trapping PARP1–DNA complexes; `PARP1` knockout is not equivalent and can even *rescue* (Murai 2012 PMID 23118055; Pettitt 2018 PMID 29748565, PARP1 point mutants confer PARPi resistance) `[verified]`. So BRCA1/2+ lines need not show stronger PARP1 KO dependency. Same logic for TOP1/TOP2 poisons.
2. **Neomorphic enzymes.** Mutant IDH1/IDH2 gain a new activity (2-HG); the drug blocks the neoactivity whereas KO removes the wild-type enzyme — direction is undefined or opposite (rows 42–44).
3. **Non-cell-autonomous / immune biomarkers.** MSI-H, dMMR, TMB-high, PD-L1 predict checkpoint-inhibitor benefit through T cells; a monoculture CRISPR screen cannot test them (rows 35–37). Note the *same* biomarker **is** testable for WRN (row 38) — MSI is non-testable for pembrolizumab and testable for HRO761. Direction depends on the drug, not on the biomarker.
4. **Antibody–drug conjugates.** Target expression is a delivery requirement; cytotoxicity comes from the payload. ERBB2 KO ≠ T-DXd in HER2-low; likewise MET/telisotuzumab, FOLR1/mirvetuximab (rows 10, 11, 52).
5. **Partial vs complete inhibition / pan-lethality.** Many drugs partially inhibit targets. 1,867 genes are pan-lethal in CRISPR across the DepMap collection, masking selectivity that RNAi retains for ~50% of them (Krill-Burger 2023 PMID 37612728) `[verified]`. Pan-essential targets (e.g. RNA-Pol, proteasome, some spliceosome drugs) will show no biomarker-conditioned signal.
6. **Off-target drug mechanism.** For several clinical-stage agents the nominal target is non-essential and the drug kills via another protein (Lin 2019 PMID 31511426) `[verified]`. A missing dependency here invalidates the *drug annotation*, not the subtype.
7. **Dominant-negative / structural alleles** (e.g. mutant TP53, some ESR1 fusions): KO of the gene is not the clinical perturbation.
8. **Resistance biomarkers and gatekeeper mutations** (ABL1 T315I, EGFR C797S, KRAS mut vs anti-EGFR): CRISPR KO removes the whole protein, so acquired-resistance alleles are invisible by design. Use these only as negative controls.
9. **Downstream / pathway-node mismatch.** The biomarker is upstream of the drugged node (BRAF V600E → MEK inhibitor; ER+ → CDK4/6). Expect the dependency on the *node*, with weaker effect size than on the mutated gene itself.
10. **Hormonal/paracrine dependence.** ER+ lines in standard medium often lose oestrogen dependence; ESR1 KO may be weakly selective in culture even though ER+ is the strongest predictive biomarker in breast cancer.
11. **Co-located passenger confound.** MTAP sits beside CDKN2A at 9p21; MTAP-null is frequently a passenger of CDKN2A deletion (rows 39–41). Guide-level copy-number effects and co-deletion must be modelled or the biomarker is attributable to the wrong gene.
12. **Lineage confounding masquerading as a biomarker.** Iorio 2016 (PMID 27397505) emphasised tissue lineage as a dominant mediator of drug response `[verified]`; Pacini 2024 (PMID 38215750) found many dependency–marker links are **gain-of-function addiction rather than synthetic lethality** `[verified]`.

---

## 4. Prior work and how it quantified agreement

- **Behan 2019** (PMID 30971826) `[verified]`: 324 lines / 30 cancer types; integrated fitness + genomic biomarkers + tractability into a priority score; validated WRN–MSI prospectively. Precedent for *framework-level* rather than recall-level validation.
- **Pacini 2024** (PMID 38215750) `[verified]`: 930 lines; "clinically informed" dependency–marker associations + PPI networks → 370 priority targets for 27 cancer types; targets mapped back to sequenced tumour cohorts for prevalence. Closest template for our validation design.
- **Iorio 2016** (PMID 27397505) `[verified]`: 11,289 tumours mapped onto 1,001 lines × 265 drugs; showed lines recapitulate tumour oncogenic alterations; logic models for combination sensitising events.
- **Gonçalves 2020** (PMID 32627965) `[verified]`: 199,219 drug-response measurements × genome-wide CRISPR in 484 lines; **reported an enrichment of positive associations between drug-sensitivity profile and KO of the drug's nominal target** — the published quantitative handle on target-KO↔drug concordance, and the method to separate on- from off-target activity.
- **Lin 2019** (PMID 31511426) `[verified]`: the discordance counter-example set.
- **Jaaks 2022** (PMID 35197630) `[verified]`: 2,025 clinically relevant combinations across 125 breast/colon/pancreatic lines; synergy is rare and context-dependent — relevant to any combination-biomarker claim.
- **Lee 2018 ISLE** (PMID 29959327) and **Lee 2021 SELECT** (PMID 33857424) `[verified]`: SELECT was predictive of response in **80% of 35 published trials across 10 cancer types** plus WINTHER. The explicit precedent for scoring a cell-line-derived model against trial-level outcomes; Sinha 2024 PERCEPTION (PMID 38637658) `[verified]` extends it to single-cell.
- **Picco 2019** (PMID 31097696) `[verified]`: 8,354 fusion events in 1,011 lines functionally linked by CRISPR + drug screens — the method for fusion biomarkers (rows 13–17, 45, 47, 50).
- No published study was found this session that reports **precision/recall of DepMap-derived biomarkers against OncoKB/CIViC levels**; several targeted PubMed queries returned nothing. Treat that as a gap our analysis can fill (and as a reason to define the metric ourselves).

---

## 5. Design implications

1. **Build the ground truth from CIViC + MOAlmanac + the FDA CDx scrape only.** All three are redistributable (CC0, CC BY 4.0, public domain). Use OncoKB's free `cancerGeneList.txt`/`allCuratedGenes.txt` for *gene-level* highest-evidence annotation, and do not depend on `allActionableVariants` (401 without a licence). Exclude CKB.
2. **Freeze and version every download** (CIViC nightly date, MOAlmanac Zenodo DOI, FDA page "Content current as of" date, CTG `dataTimestamp`) into `data/clinical/` with a manifest; CIViC changes nightly and results will not reproduce otherwise.
3. **Define the benchmark as the 5-tuple** `(biomarker_feature, feature_type, drug, target_gene, indication)` plus a pre-registered expected sign. Store the expected sign *before* looking at any dependency data.
4. **Pre-register the discordance partition.** Score recall only on the `S`/`W` rows (38 of 45 here); report the `N` rows separately as an expected-null set. Recovering a biomarker in the `N` partition is a red flag for confounding, not a success.
5. **Two-level metric.** (a) *Gene-level recall*: fraction of benchmark target genes whose dependency is explained by the matched biomarker feature in any subtype. (b) *Indication-level precision*: of subtypes recovering a pair, the fraction whose lineage composition is enriched for the approved indication (DOID/OncoTree-matched). Report direction-of-effect accuracy as a third, separate number — a sign error is worse than a miss.
6. **Stratify by evidence tier** and report recall per tier: CIViC A/B + Supports (849 items, 630 pairs), CIViC assertions Tier I-Level A (30 Predictive), FDA-CDx-backed, ESCAT I-A/I-B. Recall should rise monotonically with tier; if it does not, the subtype model is picking up something other than clinical actionability.
7. **Keep resistance pairs as negative controls**, not as recall targets, and verify the model does *not* assign them a sensitivity-like direction (CIViC supplies 1,051 Resistance items with direction labels).
8. **Model the 9p21 co-deletion explicitly** for MTAP/CDKN2A and, generally, regress out copy-number and chromosomal-proximity effects before attributing a dependency to a deleted gene. This is the single most likely source of a spurious "clinically actionable" subtype in our setting.
9. **Treat pan-essential genes as out of scope** for biomarker claims (cf. 1,867 pan-lethals); pre-filter them and say so, rather than reporting them as unexplained misses.
10. **Use GENIE (registered) for prevalence and TCGA via the open cBioPortal API for lineage denominators.** Report, for each recovered subtype, biomarker prevalence in the matched indication alongside the cell-line prevalence — the Pacini 2024 "map targets to sequenced tumour cohorts" step. Do not redistribute GENIE rows; ship only aggregate counts.
11. **Harvest companion trials programmatically from CTG v2** (`query.term` + regex over `eligibilityCriteria`) so each benchmark row carries at least one biomarker-selected NCT id; this is how the trial-stage rows (MTAP/PRMT5, MSI/WRN) enter the ground truth with evidence, and it future-proofs the table against new approvals.
12. **Include the trial-stage synthetic-lethal set as the forward-looking test.** MTAP-del→PRMT5/MAT2A and MSI→WRN were *discovered in exactly our data type* and then reached phase I; if the subtypes recover them with correct lineage breadth (tumour-agnostic), that is the strongest available evidence that the pipeline would have found a real drug target prospectively.
13. **Report the limits section as a result, not a caveat**: the 7 `N` rows here (PARP/BRCA, IDH neomorphs, immune biomarkers, ADCs) are a concrete, citable answer to "which gene dependencies cannot be used as cancer biomarkers and why".

---

### Confirmed gaps / to-do
- OncoKB `oncokb.org/terms` and `/api-access` render only via JavaScript; the exact academic-licence text was **not** read this session — confirm before relying on redistribution of any OncoKB-derived field.
- No scriptable Cancer Genome Interpreter biomarker download was found (SPA; all candidate paths return the HTML shell) and the dataset is from Oct 2022. Either obtain the file manually or drop CGI.
- FDA CDx list has no bulk file; the scraper must tolerate a table whose rows do not always align in extracted text (observed mis-pairings in the tail of the page) — validate parsed rows against CIViC `fda_companion_test` assertions and MOAlmanac FDA indications.
