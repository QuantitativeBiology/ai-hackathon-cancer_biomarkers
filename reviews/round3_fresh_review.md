# Round 3 (final) — fresh independent review of DESIGN v3.1

Reviewer: fresh agent (Opus 5.5), no access to earlier versions or reviews/. Read: PROMPT.md, DESIGN.md v3.1,
configs/{design_params,primary_endpoints,tiers,seeds}.json, docs/part1_brief.md (skim), plus spot checks of
configs/positive_controls_spec.csv, configs/build_params.json, results/data_status.json and Part 1 `configs/splits.csv`
(read-only, lineage names only).

**Verdict: not ready to freeze as is. Ready once R3-1 to R3-4 are applied.** All six are text or config edits; none
changes scope or adds compute. PROMPT.md coverage is complete: subtypes within and across lineages, a single-feature
benchmark, combinatorial screens, confounders (culture, passengers, co-location, multiplicity), clinical recovery with
direction, actionability and population, and limits. Methylation (D7) is in the design text but not yet in the configs
(R3-1). The compute is feasible. A back-of-envelope estimate for one null draw is about 61k fixed-λ EN fits on 300 genes,
roughly 1–3 min on 6 workers, so B = 20 fits in 180 min. The 200 vectorised S1 permutations are matrix products costing a
few Tflop in total, and nothing needs more than 12 GB. The schedule assumed 17 h; with about 16 h left, the reserve drops
to 2.0 h (`tiers.json: hours_remaining_at_v3` → 16.0, `reserve_h` → 2.0). That is acceptable.

Lineage map check: all 28 dev OncotreeLineage values map to one `subtypes_s8.lineage_groups` group. The cross-lineage
label precedence was dry-run and every branch is reachable (`undetermined` catches, for example, 3 powered same-sign
lineages in one organ group with no interaction). The gate (60% / 50%) is reachable both ways with plausible family
outcomes (ADD-CN 3/5, SIG 1/2, PSMC2 failing).

---

## R3-1 (critical) — Amendment D7 (methylation) is in DESIGN text only; the configs still encode "methylation not tested"

**Problem.**
- (a) `design_params.surrogate_s4.stages.G` lacks `SIL-METH` and `METH`. `primary_endpoints.P1_primary.universe_families`
  points to that key, so the code-defined U_G for the primary endpoint excludes methylation, while DESIGN §3 includes it.
- (b) `design_params.data_s1.not_obtained` still lists "RRBS methylation (SIL-METH absent; methylation not tested)", and
  `obtained` does not list RRBS.
- (c) `build_params.data_status.files[RRBS]` has `tier = not_obtained_per_design_v3` and `used_for` = "none … recorded but
  unused".
- (d) `results/data_status.json` (untracked) is stale: its RRBS row has `sha256 = null` and says "promoter beta ≥ 0.7",
  but `meth_high` is 0.5. This file is the data-deadline anchor (§1.5).
- (e) Missing-data handling has no RRBS indicator. About 1/3 of D lines lack RRBS (non-CCLE, newer and organoid lines,
  which correlate with library and culture). §1.6 sets NA → 0 with indicators only for `CN_absent`/`MSI_unknown`, so
  imputed METH zeros are confounded.
- (f) Own-gene methylation is not self-null, although SIL_g is. METH_g and SIL-METH_g silence g itself and would yield
  trivial W-sign hits.
- (g) §3 "Not assayed" bullet is a garbled sentence.

**Change.**
- `design_params.surrogate_s4.stages.G`: append `"SIL-METH", "METH"`.
- `data_s1.obtained`: append `"CCLE RRBS TSS1kb 20181022 (D7)"`.
- `data_s1.not_obtained`: `["PRISM (unused)"]`.
- `surrogate_s4.fwl_unpenalised`: append `"METH_absent"`. §1.6: "indicators `CN_absent`, `MSI_unknown`, `METH_absent`".
- `dictionary_s3.methylation.self_null`: `"METH_g and SIL-METH_g are self-null (removed from 𝒟_g, kept as NC3 checks)"`.
  Add both to the §3 self-null sentence.
- `dictionary_s3.methylation.access_class`: `{"METH": "proxy", "SIL-METH": "expression-derived"}`. Add the matching
  §4.4 wording.
- `build_params.data_status.files[RRBS]`: `tier = "essential"`,
  `used_for = "METH, SIL-METH (design_params dictionary_s3.methylation, meth_high 0.5)"`, `fallback = null`.
- Regenerate `results/data_status.json` with `compute_sha=True` and commit it before `build_dictionary.py`.
- §3 last bullet: "**Not assayed:** methylation for lines without RRBS (`METH_absent`); coverage reported (§12.2)."

**Script check.** `gate.py params` gains these asserts:
- `{"METH","SIL-METH"} <= set(dp["surrogate_s4"]["stages"]["G"])`
- `not any("RRBS" in x for x in dp["data_s1"]["not_obtained"])`
- `"METH_absent" in dp["surrogate_s4"]["fwl_unpenalised"]`
- in `data_status.json`, the row whose path matches `*RRBS*` has `sha256_matches is True` and its `used_for` contains no
  `0.7`

`tests/test_dictionary.py` asserts that METH_g and SIL-METH_g are absent from 𝒟_g.

## R3-2 (major) — positive_controls_spec.csv contradicts DESIGN §9.2 and the shared proxy rule

**Problem.** The spec (frozen at `freeze-design`) has four mismatches with DESIGN §9.2 and §10.3:
- BRAF has `feature_primary = HOT:BRAF`; §9.2 requires `ALLELE:BRAF_V600`, with HOT:BRAF as the labelled secondary proxy.
- EGFR has `HOT:EGFR` primary and `AMP:EGFR` secondary; §9.2 requires `ALLELE:EGFR_activating` primary, HOT:EGFR
  secondary, and AMP:EGFR "non-clinical".
- `SIL:MLH1`, `SIL:UBB` and `SIL:MTAP` are not resolvable ids, because §3 splits SIL into SIL-DEL / SIL-MUT / SIL-METH
  (stage G) and SIL-other (stage X).
- KIF18A's fallback `SIG:CIN_PROXY` is not a `dictionary_s3.sig_features` id.

The gate denominator and recovery would therefore be computed on features the design says are proxies or do not exist.

**Change** (spec rows):
- `ADD-HOT,BRAF,ALLELE:BRAF_V600,HOT:BRAF,,S,…,HOT:BRAF = labelled proxy`
- `ADD-HOT,EGFR,ALLELE:EGFR_activating,HOT:EGFR,AMP:EGFR,S,…,AMP:EGFR non-clinical`
- KIF18A fallback → `SIG:CIN`.
- Add to §9.1 and to `gate_s9`: `"sil_expansion": "a spec id SIL:<x> expands in build_positive_controls.py to every
  SIL-DEL/SIL-MUT/SIL-METH/SIL-other:<x> feature in 𝒟; recovered if any expanded feature is recovered with the expected
  sign"`.
- §3: SIL mechanisms are mutually exclusive features with precedence DEL > MUT > METH > other, and form one variant group
  per gene.

**Script check.** `build_positive_controls.py` fails if any `feature_primary`/`feature_secondary`/`fallback_feature`
(after SIL expansion) is not a `dictionary_v1.json` id or a declared family pattern. `tests/test_controls_spec.py`
asserts that the BRAF and EGFR primaries start with `ALLELE:`, and that every `HOT:` entry for a gene with an
`allele_classes` member sits in `feature_secondary`.

## R3-3 (major) — De-duplication before top-k silently removes allele-level and CLIN headline CGT features (P2)

**Problem.**
- §5.1 de-duplicates each (g, variant group) to its best-ranked member before any top-k. Allele groups contain
  HOT:x, ALLELE:x_* and the CLIN features built from them.
- The P2 hit rule requires `d = feature_D` exactly, and the non-proxy headline uses ALLELE/CLIN features.
- Whenever HOT:BRAF outranks ALLELE:BRAF_V600 (likely: almost the same carriers, more of them), or AMP:ERBB2 outranks
  CLIN:HER2_POS, or HOT:KRAS outranks ALLELE:KRAS_G12C, the headline feature is deleted before ranking. The test then
  scores as a miss for reasons unrelated to the arm. In the EN arms, stability selection also splits frequency across the
  collinear members.
- P2 can thus be driven by which group member wins, and `miss_decomposition` has no category for it.

**Change** (`primary_endpoints.P2_primary_clinical`, new key):
- `"rank_of_feature_D": "1 + number of rows of the arm's de-duplicated U_all ranking, excluding (g, variant group of
  feature_D), whose rank_stat is better than feature_D's own rank_stat; hit iff rank_of_feature_D <= k and g in T. The
  group's best member never lends its rank (no proxying)."`
- Add `"lost_to_group_member"` to `clinical_s10.miss_decomposition` for tests where a group member ranks ≤ k but
  feature_D does not.
- Apply the same rank to `known_recovery.csv`.

**Script check.** `tests/test_p2_rank.py`, synthetic case: HOT:BRAF ranked 3, ALLELE:BRAF_V600 ranked 150 in the raw list,
150 other groups ranked 4–149 with duplicates. Expected: rank_of_feature_D equals 1 + the de-duplicated count of better
rows, the hit is true for k = 200, and a KRAS_G12C test is never credited from HOT:KRAS's rank.

## R3-4 (major) — P2 has no pre-registered contrast or verdict branches

**Problem.** `P2_primary_clinical` defines a statistic (recall at k = 200) and "exact McNemar per pair of arms". It names
no primary pair, α or claim rule, so 6 arm pairs × 2 sets × 2 subsets leave the verdict open to choice after results are
seen. P1, by contrast, has explicit branches.

**Change** (`primary_endpoints.P2_primary_clinical`):
- `"primary_contrast": {"pair": "EN-pred vs S1-obs", "set": "non_proxy", "subset": "all (gate_disjoint if
  gate_iterations.csv has > 1 row)", "test": "exact McNemar, two-sided", "alpha": 0.05, "branches": {"p < 0.05 and
  EN-pred recall higher": "predictor recovers more clinical biomarkers", "p < 0.05 and S1-obs higher": "univariate
  recovers more", "else": "no detectable difference"}}`
- `"per_arm": "recall with exact (Clopper–Pearson) 95% CI"`
- `"other_pairs": "secondary, unadjusted, labelled"`
- If `p2_status = descriptive`, branches are not evaluated.
- Mirror this in DESIGN §10.5 in one sentence.

**Script check.** `tests/test_endpoint_branches.py` feeds synthetic discordant counts (b, c) = (8, 0), (0, 8) and (3, 3)
with 20 tests. It asserts the three P2 branches in that order, and that n = 12 yields `descriptive` with no branch.

## R3-5 (major) — P1 edge cases: the stratum contrast is not defined, and k_G can be degenerate

**Problem.**
- The beyond-inputs claim "uses the non-direct stratum with the same rule". It is unspecified whether that means (i)
  splitting each arm's frozen top-k_G by access class or (ii) re-ranking within Pool_G ∩ non_direct. Under (i), arms get
  different, possibly tiny or zero, denominators and the contrast can be undefined.
- `k_G = min(100, ⌊|Pool_G|/3⌋)` can be 0 or very small, with no fallback (P2 has one).
- In EN arms, if fewer than k_G pool rows have selection frequency > 0, the list is filled by zero-frequency rows ordered
  by feature_id. Their replication is then inherited from the power-filtered pool, not from the arm.

**Change** (`primary_endpoints.P1_primary`):
- `"strata_rule": "re-rank within Pool_G ∩ stratum; k_s = min(100, floor(|Pool_s|/3)), committed in p1_pool.json; same
  bootstrap and branches"`
- `"min_k": 20, "below_min": "descriptive (underpowered: k); no claim"`, applied to k_G and every k_s, including proxy
  for `beyond_any_input`
- `"en_zero_freq": "rows with selection_freq = 0 are not rankable; if fewer than k remain, the list is shorter and n_list
  is reported"`

**Script check.** `tests/test_endpoint_branches.py` adds three cases:
- a pool of 45 gives k = 15, so the result is descriptive with no claim;
- a non_direct pool of 300 gives k_s = 100;
- an EN arm with 60 nonzero-frequency rows and k = 100 gives n_list = 60.

The existing synthetic P1 test confirms all three contrast branches are reachable on the re-ranked strata.

## R3-6 (minor) — NC2–NC6 bound derivation is not executable as written

**Problem.** §7.3 and `toy_s13.writes_nc_bounds` replace the provisional bounds with "the 95th percentile (5th for NC6) of
the toy statistic at real D dimensions". This has three problems:
- The toy is 600 lines × 300 genes.
- The separate bound-simulation script was cut (`tiers.removed_in_v3`).
- One toy run gives one value, not a distribution.

The freeze-bounds step therefore has no defined computation. Nothing stops a replaced bound from being looser than the
provisional one.

**Change.**
- `toy_s13`: `"nc_bound_replicates": 5, "nc_bound_seeds": "toy_seed + r, r = 0..4", "nc_bound_rule": "max over
  replicates (min for NC6) at toy dimensions; if the 5 replicates exceed the toy_run_bounds timebox, the provisional
  values become final (deviations.csv)"`.
- Delete "at real D dimensions" in §7.3 and in `writes_nc_bounds`.
- Record both the provisional and the toy-derived bound in `frozen_hashes.json:bounds`.

**Script check.** `gate.py params` asserts that `negative_controls_s7.status_NC2_NC6 == "frozen"` after `freeze-bounds`,
and that each NC2–NC6 value equals the recomputed max/min over the 5 replicate outputs in `results/toy/nc_bounds.json`.

---

**Not raised (checked, fine):**
- Each §F row has a config key, apart from the D7 keys in R3-1.
- The gate family thresholds are reachable.
- The 𝒟 binary count rule (≥ 5) makes the §12.2 class-5 branch "< 5 D carriers" vacuous for MD-strict pairs. This is
  harmless, because class 5 stays reachable via R power < 0.5.
- The guard, freeze-tag ordering and single-R-run rules are internally consistent.
