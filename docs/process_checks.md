# Process checks (LESSONS rule 3)

Author: process-guard agent. Scripts: `scripts/budget_check.py`, `scripts/leakscan.py`, `scripts/check_frozen.py`;
tests: `tests/test_budget_check.py`, `tests/test_leakscan.py`, `tests/test_frozen.py`, plus the §13.2 unit tests
listed in §4. None of these modifies the engine, DESIGN.md or a frozen config.

## 1. `budget_check.py` (DESIGN §14, `tiers.json._comment`)

`python -X utf8 scripts/budget_check.py [--milestone M4] [--dry-run] [--json out.json]`; exit 1 on FAIL.

- Reads `results/timing/*.json` (top level only). Stems map to stages through `tiers.json.stage_job_map`.
- **FAIL**: a stem not in `stage_job_map`; a record that misses the minimal schema (§1.1); more than
  `concurrency.max_total_workers` (12) workers running at one time; `cut_order` exhausted while still late.
- **Lists**: every Tier B stage whose earliest job started before all Tier A stages of its milestone had ended
  (`TIER-B-EARLY`); stages without all their timing files as `not done: time`.
- **Lateness / cuts**: stage wall = union of its jobs' [start, end] intervals; lateness = max(0, wall − timebox),
  summed over done stages of all milestones up to `--milestone` (no credit for early stages). If it exceeds
  `lateness_trigger_h` (1 h), the next `cut_order` entry not yet logged by budget_check is applied, once per milestone
  run; an entry whose stage is already done is skipped (cutting it recovers nothing). The cut is appended to
  `results/deviations.csv` as `DEV-BUDGET-<n>` (existing header kept; `stage` = milestone, `item` = the entry).
- **Not budgeted** (listed for visibility only): `results/benchmark_ctrl/timing/` (ctrl-universe arms; older layout `results/benchmark/timing*/`) and
  `results/dryrun/` (pseudo-R dry run), as the pipeline author decided in 82c3240.

### 1.1 Minimal timing schema

| field | type | required | meaning |
|---|---|---|---|
| `job` | str | yes | equals the file stem |
| `utc` | ISO-8601 | yes | **end** time, written when the job ends |
| `wall_s` | number | yes (unless `start_utc`) | whole-job wall seconds |
| `workers` | int ≥ 1 | yes | worker processes the job ran (1 if serial) |
| `start_utc` | ISO-8601 | recommended | start time; otherwise start = `utc` − `wall_s` |
| `run_id` | str | recommended | shared by all stems one process writes, so concurrency counts the process once |

Accepted with a WARN (legacy): `workers_used` for `workers`; a dict `wall_s` (values summed, plus `load_s`).

### 1.2 Writers that do not meet the schema (as of 82c3240; not fixed here)

| writer | stem(s) | gap | effect in budget_check |
|---|---|---|---|
| `scripts/run_arms.py` pilot | `surrogate_pilot` | dict `wall_s`, `workers_used` | WARN only |
| `scripts/run_arms.py` full (l.179) | `s1_obs`, `s1_pred` | no `workers`; no `run_id` (one scan, two stems) | FAIL (schema); once fixed without `run_id` the shared scan counts twice (12 + 12 = 24 > 12) |
| `scripts/run_arms.py` full (l.182) | `en_obs`, `en_pred` | no numeric `wall_s` (`wall_s_both_arms`, `worker_cpu_s` instead), no `start_utc` | FAIL (schema) |
| `scripts/run_arms.py` full (l.183) | `md_strict` | no `workers` | FAIL (schema) |
| `scripts/run_null.py` (l.75) | `null_fixed_s1_vectorised` | stem not in `stage_job_map` (frozen tiers.json) | FAIL (unmapped) when `--s1-vectorised` runs |
| `scripts/run_xlib.py` (l.110-111) | `xlib_s1_obs`, `xlib_en_obs` | both written at the end with their own `wall_s`, so start = end − wall is wrong for `xlib_s1_obs` (it ran first); no `start_utc`/`run_id` | concurrency counts both at the end of the job (2 × w) |
| `scripts/run_confounder_audit.py` (l.100) | `vif` | no `workers` | FAIL (schema) |
| `scripts/run_q4.py` (l.115) | `q4_summary` | no `workers` | FAIL (schema) |
| `scripts/replicate.py` (l.216) | `replicate` | no `workers` | FAIL (schema) |
| M3b builders (`data_status.py`, `build_dictionary.py`, `build_screen_covariates.py`, `build_eligibility.py`, `build_effect_nan_mask.py`, `build_subsamples.py`, CGT features), `toy_nc_bounds.py`, `run_gate.py` | all M3b stems, `toy_pipeline`, `gate_pc` | write no timing file | stage always `not done: time` |

`src/biomark/analyses.write_timing` is the natural place to add `start_utc`, `workers` and `run_id` for every caller.

### 1.3 Stage → milestone map

tiers.json has no stage → milestone field; `STAGE_MILESTONE` in the script derives it from `milestones[].content`
and §14, and a test asserts it covers exactly the tier stages: M3b data_status … toy_run_bounds; M4 s1_obs_pred,
en_obs_pred, md_strict, confounder_audit, gate, null_fixed_B, xlib_arms, independent_check, Tier B
null_rerun_selection_adaptive_B10; M5 component_readout, subtypes_gene_module, cross_lineage,
colocation_cutting_culture_hypermut, Tier B m1_proxy_attribution_panel; M6 cgt_build_freeze, cgt_scoring_q4,
combinatorial_lookup, Tier B tcga_prevalence_cbioportal; M7 replication_R; M8 limits_classes, candidate_cards.

## 2. `leakscan.py` (DESIGN §13.4)

`python -X utf8 scripts/leakscan.py [--json out.json] [--allow FILE]`; exit 1 on any hit or invalid allowlist entry.

- (i) string constants equal to an HGNC symbol of ≥ 3 characters from Part 1 `gene_map`
  (`biomark.p1.gene_map()`, read-only), or containing `ACH-\d{6}`, in `src/` and `scripts/`;
- (ii) numeric constants (`-x` is one constant) outside {0, 1, −1, 2, 0.5, 100} in `src/biomark/` except
  `src/biomark/params.py` (the params loader).
- `# literal-ok: <reason>` (non-empty reason) exempts (ii) only, as §13.4 attaches it to (ii) (coordinator ruling).
- String exemptions come only from `configs/leakscan_allow.json` (`{"allow": [{file, literal, reason}]}`): an entry
  exempts that exact literal in that exact file and nowhere else, never a number and never an `ACH-` id. An entry
  with a missing field, an empty reason, an `ACH-` literal or a missing file fails the scan; an entry that matches
  nothing is a WARN (stale). Seeded with the identifiers that collide with HGNC symbols (family ids `FUS`, `COMP`,
  TP53-status covariate name `TP53`, spectrum suffix `POLE`) and the `BAD` status word in `verify_inputs.py`:
  9 entries covering 16 occurrences. Gene symbols used as genes (`scripts/build_cgt.py`) are not allowlisted;
  they move to config.
- Every exemption is printed as an `EXEMPT` line: each `literal-ok` comment (what it suppressed, reason,
  well-formed or not; 66 markers, 63 numeric constants suppressed at 82c3240) and each allowlisted occurrence with
  its reason, so no exemption is silent.
- Tests: planted `"KRAS"` and `0.7` give exactly two hits (with injected and with the real gene_map); an
  allowlisted literal in another file still fails; the allowlist never covers numbers or `ACH-` ids; invalid
  entries fail; the project allowlist is valid and has no `build_cgt.py` entry.

## 3. `check_frozen.py` (DESIGN §13.4)

`python -X utf8 scripts/check_frozen.py [--strict-coverage] [--json out.json]`; exit 1 on FAIL.

1. Each `frozen_hashes.json` section (`design`, `bounds` when present) is checked against the current files and the
   file at its tag (sha256 of LF-normalised bytes). DESIGN.md may differ from `design` only through
   `design_amendments`: its hash must match the current file and the text above `## Amendments after freeze` must
   equal DESIGN.md at `freeze-design`. A mismatch prints the RESEARCH_TREE decisions (`- D<n>`) naming the file but
   still fails.
2. `freeze-bounds` must be an ancestor of the first HEAD commit touching `results/{benchmark,md,gates,null}/`; with
   no tag, no such commit may exist. A `bounds` section without the tag (or the tag without a section) fails.
3. Coverage (reported; fails with `--strict-coverage`): §13.4 files missing from the record.

Current result: PASS. Coverage gaps: none of `build_params.json`, `positive_controls_spec.csv`,
`cgt_allele_rules.json`, `fusion_known_positives.csv`, `cgt_corrections.csv`, `toy_answer_key.json` (toy plants)
existed at tag `freeze-design`, so §13.4's freeze-design list is not actually frozen by that tag; the `bounds`
section and tag do not exist yet.

## 4. DESIGN §13.2 unit-test coverage

| §13.2 item | test | status |
|---|---|---|
| inference calibration (20,000 sims × n {255, 1022} × carriers {3, 5, 10} × {homosk., carrier SD × 2}) | `tests/test_inference_calibration.py` (new; 1-7 min) | ν rule (82c3240): 3/12 FAIL, anti-conservative at n = 255, SD × 2 (0.072 / 0.073 / 0.069; 1e-3: 0.0002 / 0.0021 / 0.0025). **Deployed fallback (04aa982, df = carrier clusters − 1 for binary d): 8/12 pass, 4 FAIL conservative** (size 0.05 < 0.035): n255/3/SD×2 0.03175, n255/3/homo 0.01040, n255/5/homo 0.02930, n1022/3/homo 0.02760; none anti-conservative (max 0.0587, 1e-3 max 0.0011). Accepted by the coordinator as the pre-declared fallback; left failing |
| `test_cr2.py` | `test_stats.py::test_cr2_and_nu_match_brute_force`, `::test_hc2_matches_statsmodels_singletons` | covered under another name |
| (§5.4 `test_s1_vectorised.py`: 200 random tests vs brute CR2) | `test_stats.py::test_vectorised_scan_equals_pair_test` (≥ 15 of 25 vs `pair_test`) | partial (fewer tests; reference is `pair_test`, which is itself brute-force-checked) |
| `test_eb_mixture.py` (3 components, r > 0.99 vs oracle) | `test_stats.py::test_eb_mixture_recovers_oracle` | covered under another name |
| `test_bootstrap_clusters.py` | `tests/test_bootstrap_clusters.py` (new) | 3/3 pass |
| power df = test df | `tests/test_power_df.py` (new) | 2/3 FAIL at 82c3240; **3/3 pass at 04aa982** (`pool_power` uses `stats.shared_df`) |
| power invariant to feature rescaling | `test_stats.py::test_power_invariant_to_feature_rescaling_and_df_shared` | covered |
| cross-lineage labels | `tests/test_cross_lineage_labels.py` (new; each label reachable, one label per case); toy: `test_toy_answers.py::test_pipeline_cross_lineage_labels` (env-gated) | 4/4 pass |
| culture map | `tests/test_culture_map.py` (new; §5.3 culture block maps, formulation-column detection, built levels) | 21/21 pass |
| leakscan | `tests/test_leakscan.py` (new) | 9/9 pass |
| subsamples | `test_en.py::test_subsamples_complementary_patient_grouped_lineage_stratified` | covered |
| permuted NA flags aligned | `test_arms.py::test_permutation_keeps_strata_column_correlations_and_na_alignment` | covered |
| dictionary / locus counts | `test_dictionary.py` (CN/HOT/LOSS counts vs raw, CDKN2A/CDKN2B/MTAP locus, Jaccard window, cap) | covered |
| self-null removal | `test_dictionary.py::test_self_null_exclusions`; toy `test_pipeline_self_null` | covered |
| covariate drops (ERBB2 × AMP:ERBB2, WRN × SIG:MSI, lineage dummy) | synthetic `test_stats.py::test_covariate_drop_rules`; real genes `tests/test_covariate_drops_real.py` (new) | 4/4 pass |
| FWL | `test_en.py::test_fwl_en_equals_en_with_huge_scaled_unpenalised_columns` | covered |
| determinism | `test_en.py::test_planted_feature_selected_and_deterministic`, `test_arms.py::test_parallel_s1_equals_serial`, `test_toy_answers.py::test_determinism_and_seed_dependence` | covered (no end-to-end byte-identical rerun of `run_arms`) |
| additivity (attribution 1e-5) | — | **missing**: M1-proxy attribution (§4.3, Tier B, M5) has no code yet |
| NMF | `tests/test_nmf.py` (new; D-only fit, D-row invariance to deleting R/set A, reproduces cached PROG) | 2/3 FAIL at 82c3240 (D rows moved ≤ 8.3e-4; cache ≠ refit ≤ 2.4e-4); **3/3 pass at 04aa982** (per-role transform, float64, 1 BLAS thread; dictionary rebuilt) |
| CGT schema | `test_cgt_table.py::test_schema`, `::test_vocabularies` | covered |
| EN-pred blind | `test_en.py::test_en_pred_is_blind_to_observed_effects` | covered |
| guard refusals | `test_guard.py`; `test_endpoint_branches.py` refusal tests | covered (note: `test_guard.py::test_dev_allowed_real_data_no_log` and `::test_dev_default_all_dev` read real dev effects, and `test_loaders.py` reads OOF/test predictions; agents barred from effects must deselect them) |

## 5. Findings for the engine owner (as found at 82c3240)

Re-verified at 04aa982: items 1, 2 and 6 are fixed (calibration now conservative-only, see §4; power df and NMF
tests pass), leakscan has 0 hits (item 4), the timing writers emit the §1.1 schema (§1.2 table is historical).
Remaining remarks on 04aa982: `analyses.write_timing` defaults a missing `workers` to 1 (`max(1, int(wk or 1))`),
so a writer that forgets it is recorded as serial rather than rejected; `run_null.py --s1-vectorised` stores its
run as a sub-record of `null_fixed`, so its interval is not counted in the stage wall or in concurrency.

1. **Calibration (§5.3 fallback trigger).** With R's 255-line design and carrier SD × 2, the CR2 + ν test is
   anti-conservative at α 0.05 (0.069–0.073) and at 1e-3 for 5 and 10 carriers. §5.3 pre-declares the fallback
   "binary d uses df = carrier clusters − 1" if this test fails; that is a PI decision (and affects `replicate.py`'s
   R tests, which use this test). D (1,022 lines) passes in all six cells.
2. **Power df ≠ test df** (`src/biomark/endpoints.py:pool_power`, l.396-447): ν_R is `stats.nu_from(d̃_R)` with d̃_R
   residualised on the shared design only (l.401 "gene covariates are not used for power (approximation)", l.428-431),
   while the R test (`RFitter.fit` → `stats.pair_test`) uses df = min(n − rank(Z) − 2, ν) with own CN / EXPGATE in Z.
   (a) With gene covariates the dfs differ (synthetic: 11.7 vs 17.7/24.9 for a feature correlated with own CN;
   44.8 vs 42.6 in the other direction); (b) the n − rank(Z) − 2 cap is never applied to power (4.24 vs 4.0 in a small
   design). §6.1 says "the test and its power share one df rule". (b) rarely binds at n = 255; (a) can move a pair
   across the power ≥ 0.5 pool threshold.
3. **Cross-lineage rule (design-level, not a code bug).** Per §8.3's precedence, a plant confined to one of four
   lineages is labelled cross-lineage when the three null lineages happen to share its sign (seed 21 of 5 tried;
   interaction p = 1e-14), because "powered" uses the pooled β and same-sign/LOLO use signs only.
   Conversely, opposite-sign lineage effects give a pooled β ≈ 0 and so no powered lineage.
4. **Leakscan hits at 82c3240** (70): `scripts/build_cgt.py` 46 gene symbols (CGT rule aliases and feature builders:
   pre-registered clinical definitions, not effect leaks, but hard-coded rather than in `cgt_allele_rules.json`);
   family / covariate / spectrum ids colliding with HGNC symbols (`FUS`, `COMP`, `TP53`, `POLE`, 15 hits, most with a
   `literal-ok` marker that does not apply to (i)); `BAD` in a `verify_inputs.py` print; 8 numeric slice / timeout
   constants in `guard.py` and `p1.py` without a marker. No `ACH-` id. After the allowlist (16 occurrences exempt):
   54 hits = the 46 `build_cgt.py` symbols and the 8 `guard.py` / `p1.py` numbers, which the pipeline author is
   moving to config and marking.
5. **freeze-design coverage**: six §13.4 freeze-design files were added after the tag (§3).
6. **NMF scores are not D-only invariant** (`src/biomark/dictionary.py:426`, `model.transform` on all population
   rows at once): the basis H is fit on D only and is identical with R/set A deleted, but sklearn's
   `NMF.transform` solves all rows jointly (shared initial scale from the mean of the whole matrix, one stopping
   rule), so D-row PROG scores depend on which other rows are transformed (max |Δ| 8.3e-4; 6 % of D cells > 1e-6).
   §1.7 asks for identical hashes. Transforming D rows on their own (and other rows separately) would make D rows
   invariant. Separately, a refit does not reproduce the cached PROG family bit-for-bit (max |Δ| 2.4e-4), so a
   dictionary rebuild would not give identical hashes either.
