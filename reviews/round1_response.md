# Round 1 response — DESIGN v1 → v2

84 items (7 reviewers × 12). Decisions follow the binding PI decisions D1–D3 (RESEARCH_TREE.md) and the PI's 12
round-1 rulings. "§" refers to DESIGN.md v2. "Owner" is the milestone that implements the check (M2 = before
`freeze-design`; M3b = infrastructure; M4 = discovery; M5 = subtypes; M6 = clinical; M7 = freeze + replication; M8 =
report). Every accepted or modified row becomes a row of `configs/design_ledger.csv` (§13.6).

## Conflicts resolved beyond the PI rulings

1. **Self-null vs CYCLOPS** (PW2 vs PW5/DA2): DEL_g, SIL_g, LOFBI_g and DEL of g's locus group are removed from 𝒟_g;
   hemizygous LOSS_g stays a determinant (CYCLOPS), AMP_g stays (addiction).
2. **Hierarchical surrogate (PI-4) inside the AC1 algorithm:** two FWL stages (genomic G, then expression X on the
   G residual), SIS K = 300 per stage, same λ/stability rules.
3. **Unpenalised EN block:** ST6 asked for growth pattern too; PI ruling 3 fixes lineage + library (plus technical NA
   flags). Culture enters through the D-stage/R confounder model and the §8.6 audit.
4. **EN ranking tie-break:** ST4 suggested D-stage p; that would make EN-pred read observed effects, so ties break on
   |mean β| then feature id.
5. **XLIB independence:** discovery-line OOF predictions were trained partly on replication-line effects, so pred arms
   are not independent on XLIB; it is a valid co-endpoint for EN-obs vs S1-obs and descriptive for EN-pred vs EN-obs;
   Tier B adds EN-pred-clean (M1 refit within discovery lines).
6. **P1 pool:** "top-100 over powered candidates" is read as powered-first: the top 100 of each arm among R-evaluable
   pairs with power ≥ 0.5. Power is a property of (g, d) from one EB prior shared by all arms.
7. **BRAF V600** (CL1c) vs a gate runnable without downloads (PI-6): CGT tests require V600 (else untestable); the
   cell-line positive control uses gene-level HOT:BRAF, labelled "gene-level proxy".
8. **LOSS threshold:** 0.75 (DA2, PI ruling 6) over 0.7 (PW5).
9. **Cluster-robust SEs:** CR2 for every test (ST10) is too slow for ~9M S1 tests; CR1 vectorised for all, CR2 for the
   D-stage screen, ceiling and R tests.
10. **Output formats:** `.gitignore` ignores `*.gz`/`*.parquet`; CL5's `cgt_features.parquet` becomes `.csv`, and
    committed results are plain CSV/JSON (§13.5).
11. **P2 k:** fixed at 200 (100/500 secondary); PI-2 listed three k without a primary.
12. **EXPGATE** is a covariate in all four arms, including S1-obs, so the arms share one model (S1 is "standard" plus
    the shared covariates).
13. **Seeds:** four keys added to `configs/seeds.json` before any use (null_gene_subset, bootstrap, toy, nmf); the
    existing keys are unchanged.

## Declined (in whole or part)

- ST6 q-rule stability (fixed q = 20 per subsample, E[V] bound): replaced by PI ruling 3 (fixed λ_1SE); error rate comes
  from the permutation null.
- PW5 continuous paralog CN: redundant with LOSS-paralog and EXPPAR; adds untested continuous terms.
- PW7 COSMIC NNLS exposures: no reference genome/catalogue in budget; reference-free spectra proxies instead (PI ruling 6).
- DA6.4 Sanger Project Score probe: no time; library independence is covered by XLIB and the R non-Avana subset.
- AC4 full M1/M2 refit: replaced by PI ruling 11 (M1-proxy from T folds, M3 panel refit only).
- AC7 trying `pip install leidenalg`: Tier C uses networkx Louvain directly (PI ruling 11).
- ST12 anti-conservativeness as a pass bound: kept as a reported diagnostic (making it a bound would invite tuning the toy).

## Item table

| ID | Sev. | Decision | What changed in DESIGN (section) | Check | Owner |
|---|---|---|---|---|---|
| R1-PI-1 | critical | accept | 2×2 arms S1-obs, S1-pred, EN-obs, EN-pred (= MD-pure) on identical U; MD-strict = catalogue; claims rules (§5.1, §5.2, §0) | `results/benchmark/arms.csv` universe-equality assert; `tests/test_en_pred_blind.py` | M4 |
| R1-PI-2 | critical | modify | P1: top-100 among R-powered pairs, patient-grouped paired bootstrap, verdict branches; P2 at k = 200 (§6.1, §10.5; `primary_endpoints.json`) | report script reads only `primary_endpoints.json`, hash vs `frozen_hashes.json` | M2/M8 |
| R1-PI-3 | critical | accept | Tiers A/B/C, timeboxes, cut order, fixed-λ stability, vectorised S1 null, SIL restricted, cap 3,000 (§14, §4.1, §7.1, §3; `tiers.json`) | `scripts/budget_check.py` vs `results/timing/*.json`; dictionary cap assert | M3b–M8 |
| R1-PI-4 | major | accept | Hierarchical EN (stage G genomic, stage X expression on residual); `fidelity_G`, `dfidelity_X`; clinical expression exemption list; PROG excluded from actionability (§4.1, §12, §10.5) | toy plant "genomic with expression proxy" ≥ 4/5 in stage G; script asserts control families are valid biomarker families | M4 |
| R1-PI-5 | major | accept | Eligibility adds `r_within ≥ 0.10` after centring pred and obs within lineage × library; library in FWL block (§2, §4.1) | `results/eligibility.csv`; script asserts no gene < 0.10 in any list and all such genes class 3 | M3b/M8 |
| R1-PI-6 | major | accept | MSI fallback CMP `msi_status` ∪ ModelSubtypeFeatures; gate for all arms; bug-fix-only rule; gate-disjoint recall (§1.5, §9.3, §10.5) | `results/gate/…` rows per (control, arm); `gate_iterations.csv`; src-touch-requires-tests check | M4 |
| R1-PI-7 | major | modify | Recovered regardless of sign; direction accuracy and false-sensitivity rate separate; CGT built from sources + curated mapping with source ids (not the 45-row draft); Q4c replaced by CL9 metric (§10) | `tests/test_cgt_table.py`; toy wrong-sign pair = sign error | M6 |
| R1-PI-8 | major | accept | Q1b Fisher test; cross-lineage = powered lineages + meta + ≥ 80% sign + LOLO; lineage-unresolvable label (§8.3) | `results/subtypes/cross_lineage_test.csv`; label recompute script | M5/M7 |
| R1-PI-9 | major | modify | XLIB co-endpoint: Avana-only discovery → no-Avana replication, mixed excluded; pred-arm caveat; EN-pred-clean Tier B (§6.3) | `replication_xlib.csv`; disjoint-ID assert | M7 |
| R1-PI-10 | minor | accept | Accurate-dependent = rank ≤ n_g on fold-centred predictions; R quantile n_g/n_D (§2) | toy test: prevalence 40%, perfect predictor labels all dependent lines | M3b |
| R1-PI-11 | minor | modify | Sources already downloaded (Dede, Thompson, Esmaeili, Klingbeil); Parrish counts not reprocessed; 1 h Tier B; matched background, Fisher (§11) | `results/combinatorial/lookup.csv` columns; sources in `external_sha256.txt` | M6 |
| R1-PI-12 | minor | accept | Pair classes predictor-blind / predictor-only-not-replicated; not-assayed records for methylation and signatures; MLH1 SIL vs MSI check (§12, §3) | `pair_classes.csv`; script asserts every PROMPT family in `dictionary_v1.json` with count or `not_assayed` | M8 |
| R1-PW1 | critical | accept | EN-obs arm; set differences; falsifiable novelty claims (§5.1, §6.2, §0) | `setdiff.csv`; test: EN-obs and EN-pred configs identical except target | M4/M7 |
| R1-PW2 | critical | accept | Own-CN omitted for own-locus CN d (sensitivity); self-null features removed and used as NC3; EXPGATE covariate vs EXPDOSE determinant (§3, §5.3, §7.2) | unit test on covariate spec and 𝒟_g contents; `exp_own_kind` column | M3b |
| R1-PW3 | major | modify | `access_class` direct / potentially-direct / expression-derived / proxy with `carrier_input`; M3 `sel_full_` verified for panel genes only (no saved M3 folds) (§4.4) | non-null `access_class`; direct set = M2 inputs test | M5 |
| R1-PW4 | major | accept | Class 1 = SD < 0.15; relative dependent rule (25th pct) for frac_below > 0.5; continuous effects stated; Krill-Burger limit (§2, §12) | `gene_classes.csv` assert; PRMT5/MAT2A `dependent_rule = relative`; gate reports MTAP→PRMT5 | M3b/M8 |
| R1-PW5 | major | modify | LOSS (< 0.75) for own and paralog/partner; controls re-encoded (VPS4A/B, PAPSS1, ENO2); PSMC2 CYCLOPS control "processing-limited"; continuous paralog CN declined (§3, §9.2) | `positive_controls.csv` testable column; dictionary has LOSS | M3b |
| R1-PW6 | major | accept | Controls by family incl. SIG, PROG, SIL, P53, W, LINEAGE-CHECK; per-family gate rule (§9.2, §9.3) | ≥ 2 testable rows per family; `gate.json` per-family pass | M2/M4 |
| R1-PW7 | major | modify | Mutation-file-only reference-free spectra (MMR, POLE, APOBEC, tobacco, UV) with sanity checks; else not tested (§3) | `results/sanity/signatures.json`; `not_done` record | M3b |
| R1-PW8 | major | accept | TP53 status covariate with drop rule; SIG determinants keep sign in TP53-mut lines (§5.3, §8.7) | test: TP53 dropped exactly when d is TP53/MDM2-derived | M3b/M5 |
| R1-PW9 | major | accept | Empirical amplicon test, three labels, no literature exception (§8.5) | `amplicon_test.csv` label set | M5 |
| R1-PW10 | major | accept | `concordant-family` class with target sets; del(17p)→BCL2 population-only N; RB1→CDK4/6 W; CCNE1→PKMYT1 (§10.3) | `test_cgt_table.py` family rows, del17p row | M6 |
| R1-PW11 | minor | modify | LOFBI (dosage 2, LOSS/SIL, X-linked male); alleles only in CGT features (DA10); FUS curated drivers, orientation, FFPM (§3, §10.3) | LOFBI count test; FUS ⊆ `fusion_drivers.csv` | M3b |
| R1-PW12 | minor | accept | Culture gene sets incl. ferroptosis and amino acid; organ-system groups for cross-lineage; interaction reported (§8.3, §8.6; `design_params.subtypes_s8`) | labels recomputable from config; ≥ 2 organ groups per cross-lineage row | M5 |
| R1-ST1 | critical | accept | D-stage = selection-conditional screen; empirical FDR is the error rate; wording rule; cross-fit Tier C (§5.2, §14) | grep test: `q_sc` with `p_status`/`fdr_empirical`; REPORT lint | M4/M8 |
| R1-ST2 | critical | modify | Joint row permutation in lineage × library, full pipeline re-run, 300 genes, B = 20 (50 by timing rule), FDR ≤ 0.10 pass, not retuned (§7.1) | `tests/test_permutation_null.py`; `results/null/summary.json` | M4 |
| R1-ST3 | critical | accept | Per-test drop rule (abs ρ > 0.7, locus, own CN, MSI/TMB) with `dropped_cov`, all arms (§5.3) | unit test ERBB2/WRN; scan max abs r(d, cov) ≤ 0.7 | M3b |
| R1-ST4 | major | modify | Primary contrast EN-pred vs EN-obs; one test for all lists; rankings; gene-matched and overlap strata; McNemar for CGT; LMM top-200 Tier B (§5, §6.2, §10.5) | `arms.csv` hash/seed assert; `replication_by_k.csv`; `overlap_strata.csv` | M4/M7 |
| R1-ST5 | major | accept | R-evaluable rule; EB shrunk power; list-level statistics; paired bootstrap (§6.1, §6.2) | `results/replication/summary.json` fields; one shrinkage function hash, committed before tag | M7 |
| R1-ST6 | major | modify | CPSS 25 pairs; SIS in each fit; OOS fidelity; z-scoring; pilot. q-rule declined (fixed λ per PI ruling 3); growth pattern not unpenalised (§4.1) | `subsample_manifest.json`; limits script reads `r2_oof_beyond_lineage_library`; pilot timing commit before M4 | M4 |
| R1-ST7 | major | accept | Days as within-Humagne deviation; CellFormat collapse; `_rare` lineage; VIF gate; drop-one-block delta; disease sensitivity (§5.3, §5.5) | `vif.csv` max ≤ 10; no lineage level < 10; delta flag recompute | M3b/M4 |
| R1-ST8 | major | modify | REML + HKSJ, PI, labels merged with PI-8 rule (powered lineages, organ groups, LOLO all), calibrated on null (§8.3) | `cross_lineage.csv` columns; `false_cross_lineage_rate` | M5 |
| R1-ST9 | major | accept | Families F1–F7 in `design_params.test_families_s5` (§5.6) | script recomputes q and n_F; commit order | M2/M4 |
| R1-ST10 | major | modify | Patient-clustered CR1 (all) / CR2 (screen, R, ceiling); LMM Tier B top-200 with D-only kernel; module gene clusters (§5.3, §6.2, §8.2) | `se_type` column; `lmm_sensitivity.csv`; `n_gene_clusters ≤ n_genes` | M4/M5 |
| R1-ST11 | major | accept | Gate log, applied to all arms, gate-disjoint primary if changed, exact CI, gate not used to pick variants (§9.3, §10.5) | `gate_iterations.csv`; overlap computed from both configs | M4/M6 |
| R1-ST12 | minor | modify | Toy predictor = CV ridge on noisy obs + fold offsets; correlated blocks; nested library; fold AMI; anti-conservativeness reported, not a bound (§13.1, §8.8) | `tests/test_toy_pipeline.py` | M3b |
| R1-DA1 | critical | accept | Gene-level features kept; locus groups (Jaccard ≥ 0.5, 10 Mb); resolution test on discordant lines (§3, §8.4) | `tests/test_locus.py`; `results/locus/resolution.csv` | M3b/M5 |
| R1-DA2 | critical | accept | LOSS family; `positive_controls.csv` with exact feature ids committed before associations (§3, §9.1, §9.2) | ancestor-commit check; carrier counts recomputed | M2/M3b |
| R1-DA3 | major | accept | Own CN omitted iff d in g's own locus group; adjusted = sensitivity (§5.3) | `covariate_spec.csv` test | M3b |
| R1-DA4 | major | accept | Screened-condition covariates for D ∪ R ∪ A; explicit unknown levels; media column chosen by content (swap); 69 lines unknown; no doubling-time proxy (§5.3) | `build_screen_covariates.py` asserts no NaN, ARI ≥ 0.9 vs P1 lib_group | M3b |
| R1-DA5 | major | accept | Essential / nice / dropped downloads; fallbacks; data_status at M4 start (§1.5) | `results/data_status.json`; not-done rows match present = false | M3b/M4 |
| R1-DA6 | major | accept | Deny-list extended; same-Chronos-fit statement; R non-Avana subset reported; Sanger probe declined (time) (§1.2, §1.3) | guard tests per path; access log timestamps after tag | M3b/M7 |
| R1-DA7 | major | accept | Lineage-relative SIL with mechanism split; Cramér's V flag (§3) | dictionary test; V column | M3b |
| R1-DA8 | major | accept | Raw caches; CN NA; NA flags unpenalised; per-family exclusion; n_obs/n_carriers/n_NA (§1.6) | test: CN NA for 55 D lines; n_carriers ≥ 5 assert | M3b |
| R1-DA9 | major | accept | D-only fitting and frozen transforms; NMF transform for R; R accurate threshold by quantile; R carrier counts pre-freeze (§1.7, §2) | rebuild-without-R hash test | M3b |
| R1-DA10 | minor | accept | DMG_BURDEN rename; LOFBI primary, LOF-any sensitivity; HOT gene-level, alleles in CGT; hypermutated rule (§3, §8.7, §10.3) | LOFBI count test; `allele_level_required` column | M3b/M6 |
| R1-DA11 | minor | modify | Joins on ModelID/Entrez; HGNC symbol map; networkx Louvain declared; data_availability.md fixed; crispr.parquet D-row spot check (§3, §8.8, §13.6) | `verify_inputs.py` spot check; `idmap_unmapped.csv` < 5% | M3b |
| R1-DA12 | minor | accept | DATA.md rows for every external file; licence flags; measured screens only (§11, §13.6) | `scripts/verify_data_md.py`; OncoKB text grep | M3b |
| R1-AC1 | critical | accept | FWL, SIS K = 300 (per stage) with forced terms, 25-λ path, P1 cv_fold, 1-SE fixed λ, 50 subsamples never reduced, pilot → K = 150 (§4.1) | `surrogate_pilot.json`; FWL vs scaled-lineage test; config hash keys | M4 |
| R1-AC2 | critical | modify | B = 20 joint permutation; strata lineage × library (PI ruling 4); 300 genes; 3 h cap (§7.1) | `perm_manifest.json`; invariance test; toy null B = 20 | M4 |
| R1-AC3 | major | modify | Cap ≤ 3,000 with escalation; SIL restricted; locus collapse replaced by locus groups (DA1) (§3) | `family_counts.csv` cap assert; locus test | M3b |
| R1-AC4 | major | modify | No full refit (PI ruling 11): M1-proxy from T folds, M3 panel refit, r ≥ 0.999, one 30-min diagnosis, additivity (§4.3) | `results/refit/repro_check.csv`; additivity test | M5 |
| R1-AC5 | major | accept | Component readout Tier A, fallback for §4.3 (§4.2) | `component_selection.csv` | M4 |
| R1-AC6 | major | accept | OOF rule, fold-centred contributions, gene blocks ≤ 500, PROG coefficients, M3 mapping, M1-mediated bucket (§4.3) | family-sum test; RSS < 6 GB log | M5 |
| R1-AC7 | major | modify | Tier C; fold-centred OOF input; networkx Louvain; seed stability; fold test (§8.8) | `results/clusters/*.csv`, `backend.txt` | M5 (Tier C) |
| R1-AC8 | major | accept | Vectorised S1 by NaN pattern, partial t, FWL own-gene terms; extras restricted to hits ∪ controls; LMM Tier B or dropped (§5.4, §5.5, §6.2) | `tests/test_s1_vectorised.py` | M4 |
| R1-AC9 | major | accept | Compute and schedule table, cut order, never-cut list (§14; `tiers.json`) | `budget_check.py` | all |
| R1-AC10 | minor | accept | Complementary pairs; OOF fidelity beyond lineage + library (§4.1) | `subsample_manifest.json` | M4 |
| R1-AC11 | minor | accept | loky, 1 BLAS thread, memmaps, per-gene seeds, no concurrent refit, 12 GB (§14) | `tests/test_determinism.py` | M3b |
| R1-AC12 | minor | accept | NMF settings pinned, D-only, transform for R (§3) | NMF reproducibility test | M3b |
| R1-AF1 | critical | accept | Single entry point, extra unblind conditions (caller, catalogue hashes, tag target, reproduction-only), pseudo-R dry run, immutability diff (§1.3) | `gate.py leakscan`, `gate.py G_rep`; extended `tests/test_guard.py` | M3b/M7 |
| R1-AF2 | critical | accept | §F table and config files; frozen hashes; D-entry deviation rule; automatic `deviations.csv`; no decision literals (§13.4, §F) | `gate.py frozen`; literal-grep test | M2 |
| R1-AF3 | critical | accept | Every CGT row needs a resolvable source id; build script reads no effects; freeze-cgt ancestry; combinatorial labels from files only (§10.1, §11) | `gate.py cgt`; import scan | M6 |
| R1-AF4 | major | accept | Controls by feature id; pc_eligibility before freeze; immutable attempt 1; reporting rule (§9) | `tests/test_controls.py`; `gate.py G_pc` | M2/M4 |
| R1-AF5 | major | accept | Report-number check ported; claims table and `c:` tags; .gitignore checked (`*.gz` ignores results) → plain CSV/JSON rule + tracked check (§13.5) | `check_report_numbers.py --strict`; `gate.py tracked`; claim-word regex | M3b/M8 |
| R1-AF6 | major | accept | Part 1 manifest snapshot; HEAD check; essml only in p1.py; no subprocess into P1; run at every gate (§1.3, §1.4) | `tests/test_part1_untouched.py` (manifest); `gate.py part1` | M3b |
| R1-AF7 | major | accept | Compliance ledger from this table; pip freeze (§13.6) | `gate.py ledger` | M2/M8 |
| R1-AF8 | major | accept | Numeric bounds NC1–NC6 simulated before freeze (§7.2; `negative_controls_s7`) | `gate.py controls`; bounds commit precedes real runs | M2/M4 |
| R1-AF9 | major | accept | Toy uses same entry points/params, written by a different agent; independent reimplementation of 5 quantities (§13.1, §13.3) | `gate.py G_indep` (author check) | M3b/M7 |
| R1-AF10 | minor | accept | data_status at M4 start; late files secondary only (§1.5) | `gate.py data` | M4 |
| R1-AF11 | minor | accept | Registry with dirty flag and hashes; exploratory status; D-id ranges per track (§13.5, §13.6) | report-number check resolves registry rows | M3b |
| R1-AF12 | minor | accept | No gene/line literals in code; exception lists in configs, frozen (§13.4, §8.5 has no exceptions) | leakscan for `ACH-\d{6}` and gene symbols | M3b |
| R1-CL1 | critical | accept | KIT D816V split, del17p N, BRAF V600, EGFR amp removed, AKT family, RAS merge, ESR1 LBD N, 9p21 merge (§10.3) | `tests/test_cgt_table.py` assertions | M6 |
| R1-CL2 | critical | accept | `in_gate`; recall all and gate-disjoint; gate iterations; WRN/PRMT5 never "independent" (§9.3, §10.5; `primary_endpoints.wording_rules`) | `q4_summary.json` keys; REPORT grep | M6/M8 |
| R1-CL3 | critical | accept | Scoring unit = test; comparator; sign convention; five partitions; S on W = sign error (§10.2) | schema validation in `test_cgt_table.py` | M6 |
| R1-CL4 | major | accept | Data ceiling = S1-obs on that test, comparator-restricted, BH q < 0.1, four labels; recall over testable and data-supported (§10.4) | `cgt_ceiling.csv`; label recompute; toy rows | M6 |
| R1-CL5 | major | accept | `build_cgt_features.py` (alleles, fusions with reference lines, BRCA biallelic, ER mixture, HER2, MTAP, MSI); output CSV, committed before freeze-cgt; ESR1 exempt from class 4 (§10.3, §12) | reference-line carrier test; ALK ≤ 10; commit before any CGT CRISPR read | M6 |
| R1-CL6 | major | accept | KO-concordance classes incl. concordant-family and combination-required; KRAS G12C concordant (§10.3) | vocabulary test; family rows size(T) ≥ 2 | M6 |
| R1-CL7 | major | accept | Clinically refuted partition, ≥ 3 rows with resolvable ids, recall_refuted vs approved (§10.6) | `q4_summary.json` keys; row-count test | M6 |
| R1-CL8 | minor | accept | Forward-looking trial rows only with matching cached CTG record (§10.3) | cached CTG JSON regex test | M6 |
| R1-CL9 | major | accept | `indication_map.csv`; on/beyond/off-label classes; population precision; carriers-only OR; agnostic exclusions (§10.5) | `q4c_population.csv`; precision recompute | M6 |
| R1-CL10 | minor | modify | TCGA prevalence as Tier B with Wilson CIs and log2 flag (details in `clinical_s10.prevalence`) (§10.5) | `prevalence.csv`; one cell recomputed from cached JSON | M6 (Tier B) |
| R1-CL11 | minor | accept | Mutually exclusive tiers T1–T3 + refuted; trend test only if ≥ 5 per tier (§10.3, §10.5) | one tier per test; null trend p when n < 5 | M6 |
| R1-CL12 | minor | accept | Direction accuracy, sign-error rate, N-partition ratio, pair/feature actionability with null, line coverage (§10.5) | `q4_summary.json` keys; toy sign-flip test | M6 |
