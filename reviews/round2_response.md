# Round 2 response — DESIGN v2 → v3

23 items (fresh PI 7, statistics 8, implementability 8). These decisions follow the PI round-2 rulings (decision D6 in
RESEARCH_TREE.md) and the new data fact: OmicsSomaticMutations, OmicsGlobalSignatures, CRISPRConfounders,
OmicsCNSegmentsWGS and OmicsProfiles arrived; RRBS and PRISM were not obtained. "§" refers to DESIGN.md v3.

**Owner** is the milestone that implements the check: M2 = before `freeze-design`; M3b = infrastructure, ending in
`freeze-bounds`; M4 = discovery; M5 = subtypes; M6 = clinical; M7 = freeze + R; M8 = report.

**Impl impact** says whether the item changes outputs of the tracks already implementing v2:
- **DICT** = dictionary / eligibility / controls agent;
- **CGT** = CGT agent;
- **TOY** = toy-simulation agent.

## Conflicts resolved beyond the PI rulings

1. **ALLELE and CLIN become 𝒟 families** (FP-1d, FP-6). P2 hits require d ∈ U. Allele-defined and composite clinical
   tests (ER+, HER2+, BRCA biallelic, RAS-mutant) would otherwise be unhittable by every arm. `build_cgt_features.py`
   therefore runs before the first dictionary build and supplies these features. HOT stays gene-level. HOT and its
   alleles form a variant group, de-duplicated before top-k exactly like locus groups.
2. **BRAF and EGFR gate controls** now use ALLELE primaries, with HOT as a labelled secondary proxy. This is the same
   proxy rule as the CGT (FP-1d).
3. **SE rule:** IM7 proposed df = n_clusters − 1; the PI ruling makes ST1 override it. CR2 + ν (Satterthwaite) is
   used everywhere. IM7's statsmodels-HC2 equivalence test is kept as a unit test.
4. **Budget arithmetic:** "≤ 15 h + ≥ 3 h reserve" exceeds the ~17 h left. Milestones are set to 14.0 h plus a 3.0 h
   reserve.
5. **FP-4 strata.** `potentially-direct` (M3 pool) counts as non-direct for the PI-ruled "beyond hand-coded inputs"
   claim. A stricter "beyond any input" claim is allowed only from the proxy class.
6. **Doubling time** (now available) enters as a never-dropped *growth* design block. It is added to the EN FWL block
   together with culture (ST6).
7. **XLIB** runs only S1-obs and EN-obs (FP-2). IM8's xlib_arms stage is therefore 30 min for two arms.
8. **Methylation stays "not tested"** per the PI ruling. A file `external/depmap_extra/CCLE_RRBS_TSS1kb_20181022.txt.gz`
   is present but has no checksum entry, and DATA.md says "not obtained". It is not used unless the PI registers it
   before the first dictionary build.

## Declined

None. Parts of FP-2 are superseded by ST1 (CR1 → CR2) and by IM8 (null 4 workers → 6 workers, per the PI ruling).

## Item table

| ID | Sev. | Decision | What changed in DESIGN (section) | Check | Owner | Impl impact |
|---|---|---|---|---|---|---|
| R2-FP-1 | critical | accept (resolved by data) | Mutation file arrived. ALLELE family from `cgt_allele_rules.json`; allele tests use ALLELE (`proxy = 0`). HOT:gene copies are a labelled secondary (`proxy = 1`) with exclusions. Gate and CGT share one proxy rule. P2 sets `non_proxy` (headline) / `incl_proxy` (§3, §9.2, §10.3, §10.5; `clinical_s10.proxy_rule`, `P2_primary_clinical.sets`) | script: every `proxy = 1` row has a `HOT:` feature and is not in `allele_proxy_excluded`; `q4_summary.json` has `n_testable_proxy`, `n_testable_allele` | M3b/M6 | **DICT**: add ALLELE family; BRAF/EGFR spec rows → ALLELE primary. **CGT**: `proxy` and `feature_D` columns; allele rules config |
| R2-FP-2 | major | accept (CR1 part overridden by ST1) | Cut: LMM, EN-pred-clean, M3 refit, all Tier C, `simulate_control_bounds.py`. Tier B = M1-proxy attribution, TCGA, adaptive null. XLIB 2 arms. NC bounds from the toy run. Independent-check tolerances revised. Null 6 workers. Budget 14 + 3 h (§4.3, §6.3, §7.3, §13.3, §14; `tiers.json`) | `tiers.json` has no removed keys; `milestone_total_h + reserve_h ≤ hours_remaining`; workers sum ≤ 12; `budget_check.py` | M2 | **TOY**: toy run must write NC2–NC6 bounds. DICT, CGT: none |
| R2-FP-3 | major | accept | P1 headline on U_G (stage-G families); U_all secondary; own EXPDOSE (F2) excluded from all P1 pools; k_G = min(100, ⌊pool/3⌋) fixed before the tag (§6.1; `P1_primary.universe`, `k_rule`) | `p1_pool.json` committed before `freeze-catalogue`; script: no P1 row has an expression family or F2 | M7 | none |
| R2-FP-4 | major | accept | P1 and gate reported by access class; "beyond inputs" only from non-direct, "beyond any input" only from proxy; PAR recall on non-direct encodings (§4.4, §6.1, §9.3; `P1_primary.strata_reported`, `claim_rules`) | `p1_by_access.csv`; `claims.csv` has no `beyond_inputs` row without `stratum = non_direct` | M4/M7/M8 | **DICT**: `access_class` column in `positive_controls.csv` and `arms.csv` |
| R2-FP-5 | major | accept | Candidate cards `known_recovery.csv` and `shortlist.csv` (ranked; `not_in_curated_lists`; set-A flag); REPORT opens with them plus 3 rule-chosen case studies (§0, §12.3; `limits_s12.cards`, `report_opening`) | shortlist ⊆ class 9 ∧ R replicated; both files tracked and registered; no bare "novel" in REPORT.md | M8 | none |
| R2-FP-6 | minor | accept | Non-eligible controls and CGT targets go to `U_ctrl`, excluded from top-k, null and P1/P2; P2 miss decomposition; CGT `feature_D` (§2, §10.2, §10.5) | arms ⊆ eligible; `null_genes.txt` ⊆ eligible; `miss_decomposition` sums to misses | M3b/M6 | **DICT**: eligibility marks `U_ctrl`. **CGT**: `feature_D` column |
| R2-FP-7 | minor | accept | Combinatorial lookup moved to Tier A, 45 min, never cut (§11; `tiers.json`) | `combinatorial_lookup` in `tiers.A.stages` and `never_cut`; `lookup.csv` exists at final gate | M6 | none |
| R2-ST1 | critical | accept | CR2 everywhere (HC2 for singletons), vectorised per (NaN group, d). t with df = min(n − rank − 2, ν), ν Satterthwaite from W_c; carrier-cluster fallback. Leverage-1 rule (merge levels < 2 lines; exclude 1 − h < 1e-8 from the meat). Power by noncentral t with ν_R. Bootstrap: new cluster ids per copy; degenerate rule; `frac_degenerate` (§5.3, §6.1; `confounders_s5.df_rule`, `leverage`, `replication_s6`) | `test_inference_calibration.py` (size bounds); `test_s1_vectorised.py` vs brute-force CR2 and statsmodels HC2; power df = test df | M3b/M4 | **DICT**: `build_screen_covariates` level-merge logging only |
| R2-ST2 | critical | accept | Primary null = fixed-threshold (observed λ, p*_F1, p*_sc); selection-adaptive re-run secondary (B = 10, Tier B); NA flags permute with rows; toy checks true FDP on the planted toy and P(any hit) on a pure null; NC1 fixed at 0.10 (§7.1–7.3, §13.1; `null_s7`, `toy_s13.pass`) | `null/summary.json` keys; p*_F1 recomputed; NA-alignment test; toy inequalities | M3b/M4 | **TOY**: planted-toy FDP and pure-null tests. **DICT**: NA indicators stored as 𝒟 row attributes |
| R2-ST3 | major | accept | Droppable set {MSI, DMG_BURDEN, TP53, own CN, EXPGATE, screen days}; categorical association = R² on dummies > 0.49; design blocks never dropped; `collinear_with_design` label → class 6 (§5.3; `confounders_s5`) | scan: `dropped_cov` ⊆ droppable; flagged rows have R²(design) > 0.49; lineage-dummy unit test | M3b/M4 | **DICT**: covariate spec separates design and droppable blocks |
| R2-ST4 | major | accept | In-fit carrier floor 2; `en_max_attainable_freq` and `en_unreachable` reported for gate and P2; toy 6-carrier plant (§4.1, §9.3; `surrogate_s4.min_carriers_in_fit`) | manifest script lists features ≥ 5 carriers with attainable freq < 0.6; toy plant selected in > 50% of subsamples | M3b/M4 | **DICT**: `en_max_attainable_freq` in `positive_controls.csv`. **TOY**: add plant |
| R2-ST5 | major | accept | Continuous d: ≥ 20 non-NaN lines and SD > 0 (EXPDOSE ≥ 5 expressed); R-evaluable and SE_R formulas for continuous d; EB prior on standardised effects over all U incl. F2 (§3, §6.1; `replication_s6`) | every U row has `power`, `evaluable`, `d_type`; rescaling-invariance test | M3b/M7 | **DICT**: `d_type` and non-NaN counts in `dictionary_v1.json` |
| R2-ST6 | major | accept | Culture (and growth) block added to EN FWL; F2 excluded from P1 pools; tie-break = mean β over all 50 subsamples with zeros; `n_tied_at_rank_k` (§4.1, §6.1; `surrogate_s4`) | config test `fwl_unpenalised ⊇ culture`; no F2 row in P1 lists; `mean_beta_def = all50` | M4 | **DICT**: culture/growth covariates must be exported for D lines in FWL form |
| R2-ST7 | major | accept | Pooled significance from the restricted confounder model (CR2 + ν, BH F4). Meta descriptive with Paule–Mandel + HKSJ (statsmodels). Powered α stated. Precedence order. Lineage-restricted needs ≥ 2 testable lineages; single-lineage carriers → unresolvable (§8.3; `subtypes_s8.cross_lineage`) | `test_cross_lineage_labels.py`: all 4 labels reachable, one label per case | M5 | none |
| R2-ST8 | minor | accept | XLIB caveat: U selected on all-D effects; `u_rule = discovery_only` sensitivity. Dry run status `dryrun`, uncitable, not inspected; reachability via synthetic test (§1.3, §6.3) | `replication_xlib.csv` has `u_rule`; report checker rejects `dryrun` tags | M4/M7 | **DICT**: eligibility function must accept a line subset (discovery-only rerun) |
| R2-IM1 | critical | accept | Two tags: `freeze-design` now (NC2–NC6 provisional) and `freeze-bounds` at end of M3b (bounds, `pc_eligibility`, `positive_controls.csv`, `null_genes.txt`, `subsamples.json`, dictionary, data_status); ancestry rule for result paths (§13.4; `tiers.json` M2/M3b) | `gate.py frozen`: two hash sections; `freeze-bounds` ancestor of first commit under benchmark/md/gates/null; no `provisional` left at `freeze-bounds` | M2/M3b | **DICT**: outputs go into `freeze-bounds`, not `freeze-design` |
| R2-IM2 | critical | accept | `configs/catalogue/r_requests.csv`; one `replicate.py` run writes `r_tests.csv`, `r_lines.csv`, `setA_lines.csv` (incl. set A); downstream reads only these; dry run produces the full R-consuming file set (§1.3, §12.3) | `gate.py G_rep`: class 5/7 pairs, Q1b rows and P1/P2 rows ⊆ `r_requests.csv`; AST scan: only `replicate.py` uses replication/limits scopes; dry-run file set = planned set | M7 | none |
| R2-IM3 | critical | accept (resolved by data) | Data deadline = first dictionary build; data_status sha256 stored in `dictionary_v1.json`; late files secondary only. RAS_MUT from alleles (CLIN), not a proxy. LoH for BRCA via AF ≥ 0.7 or own CN < 0.75. P2 `min_testable_tests = 15` → descriptive below (§1.5, §3, §10.3, §10.5) | `gate.py data`: hash equality and commit order; `q4_summary.json` `p2_status` recomputed from the threshold | M3b/M6 | **DICT**: store data_status hash; `build_cgt_features` runs before the dictionary. **CGT**: `min_testable` reporting |
| R2-IM4 | major | accept | `cache/effect_nan_mask.parquet` is the only effect-derived EN-pred input (equal row sets). `configs/subsamples.json` shared by all genes/arms; per-gene seeds removed. FWL cache per (NaN group × fold/subsample). Stage X on each fit's own stage-G residual. NA indicators permute (§1.6, §4.1, §7.1) | `test_en_pred_blind.py` (row sets), `test_subsamples.py`, NA-alignment test; pilot reports cache hit rate | M3b/M4 | **DICT**: NaN mask and subsamples files (M3b) |
| R2-IM5 | major | accept | Lists in frozen configs: `positive_controls_spec.csv` (exists), `cgt_allele_rules.json`, `fusion_known_positives.csv` (exists), `cgt_corrections.csv`; AST leakscan definition with `# literal-ok` escape (§9.1, §10.1, §13.4) | `gate.py leakscan` exit 0; `test_leakscan.py` planted hits; new configs in `frozen_hashes.json` | M2/M3b | **DICT**: none beyond spec edits. **CGT**: write `cgt_allele_rules.json`, `cgt_corrections.csv`; build scripts read them only |
| R2-IM6 | major | accept | Keys added: line share definition, cutting α/depleted, culture V 0.3, `limits_s12` (universe, aggregation, class 3b, class 5), culture maps (in `build_params.json`), fusion orientation, partner priors, SE_R formulas, R NaN rule, XLIB library-string rules. Procedural rules declared; toy authorship by `Agent:` trailer (§8.1, §8.5, §8.6, §12.1, §13.1, §13.6) | `gate.py params` (keys present, non-null); `test_culture_map.py` | M2/M3b | **DICT**: culture map already in `build_params.json` (consistent). **TOY**: `Agent:` trailer |
| R2-IM7 | major | modify (df per ST1) | CR2 + ν (ST1 rule, not n_clusters − 1); Paule–Mandel replaces REML; EB toy known-answer test; Deming via `scipy.odr`; bootstrap rules; independent-check tolerances in `indep_s13` (§5.3, §6.1, §8.3, §13.2, §13.3) | `test_cr2.py`, `test_eb_mixture.py`, `test_bootstrap_clusters.py`; `G_indep` reads tolerances from config | M3b/M7 | none |
| R2-IM8 | major | accept (workers per PI) | Hours 17, milestones 14 + 3 reserve. Tier A `xlib_arms` (30 min, M4 background). Null 6 / foreground 6 workers. Independent check starts after gate attempt 1. Lateness trigger 1 h → next cut. `stage_job_map`. `gate.py` priorities (§14; `tiers.json`) | `budget_check.py` fails on unmapped stems, > 12 workers, or lateness without a deviations entry | all | none |
