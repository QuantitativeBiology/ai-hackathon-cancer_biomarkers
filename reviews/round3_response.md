# Round 3 response — DESIGN v3.1 → v3.2 (frozen at `freeze-design`)

Six items from `reviews/round3_fresh_review.md`. The PI ruled to accept all six. "§" refers to DESIGN.md v3.2.

**Owner** is the milestone that implements the check: M2 = before `freeze-design`; M3b = infrastructure, ending in
`freeze-bounds`; M6 = clinical; M7 = freeze + R.

**Impl impact** says whether the item changes outputs of the tracks already implementing v3:
- **DICT** = dictionary / eligibility / controls agent, which owns `configs/build_params.json`,
  `configs/positive_controls_spec.csv` and `results/data_status.json`. Those three files were not edited here; the PI
  has told the agent its parts of R3-1 and R3-2.
- **CGT** = CGT agent.
- **TOY** = toy-simulation agent.

| ID | Sev. | Decision | What changed in DESIGN (section) / configs | Check | Owner | Impl impact |
|---|---|---|---|---|---|---|
| R3-1 | critical | accept | §1.6 and §4.1: `METH_absent` indicator. §3: METH_g and SIL-METH_g are self-null; garbled "Not assayed" sentence replaced. §4.4: METH = proxy, SIL-METH = expression-derived. Config: `surrogate_s4.stages.G` += SIL-METH, METH; `fwl_unpenalised` += METH_absent; `data_s1.obtained` += RRBS (D7); `not_obtained` = [PRISM]; `dictionary_s3.methylation.self_null`, `access_class`, `not_assayed` | `gate.py params` asserts on stages.G, not_obtained, METH_absent and the data_status RRBS row (sha256 matches; no 0.7); `tests/test_dictionary.py`: METH_g and SIL-METH_g ∉ 𝒟_g | M2/M3b | **DICT**: RRBS row in `build_params.json` (tier essential, used_for, fallback null); regenerate `data_status.json` with sha256 before `build_dictionary.py`; METH_absent column; self-null removal |
| R3-2 | major | accept (DESIGN/config parts here) | §3: SIL mechanisms are mutually exclusive features (precedence DEL > MUT > METH > other) and form one variant group per gene. §9.1: `SIL:<x>` spec ids expand; the build fails on unresolvable ids. Config: `gate_s9.sil_expansion`, `dictionary_s3.sil_mechanism_precedence`, `variant_groups` | `build_positive_controls.py` fails on unresolved ids; `tests/test_controls_spec.py` (BRAF and EGFR primaries start with `ALLELE:`; HOT for allele genes only as secondary) | M3b | **DICT**: spec rows BRAF → `ALLELE:BRAF_V600` / HOT:BRAF; EGFR → `ALLELE:EGFR_activating` / HOT:EGFR / AMP:EGFR non-clinical; KIF18A fallback → `SIG:CIN` |
| R3-3 | major | accept | §5.1 and §10.5: P2 hit uses `rank_of_feature_D` (the group's best member never lends its rank); new miss category `lost_to_group_member`; same rank in `known_recovery.csv`. Config: `P2_primary_clinical.rank_of_feature_D`, `hit`; `clinical_s10.miss_decomposition` | `tests/test_p2_rank.py` (synthetic HOT:BRAF rank 3 / ALLELE:BRAF_V600 rank 150; no KRAS_G12C credit from HOT:KRAS) | M6/M8 | **CGT**: none (scoring-time). DICT: none |
| R3-4 | major | accept | §10.5: P2 primary contrast EN-pred vs S1-obs on `non_proxy` (`gate_disjoint` if > 1 gate attempt), exact two-sided McNemar α 0.05, three branches, not evaluated when descriptive; per-arm Clopper–Pearson. Config: `P2_primary_clinical.primary_contrast`, `per_arm`, `other_pairs` (keys exactly as proposed) | `tests/test_endpoint_branches.py` ((b, c) = (8, 0), (0, 8), (3, 3) with 20 tests → three branches; n = 12 → descriptive, no branch) | M2/M6 | none |
| R3-5 | major | accept | §6.1: strata re-ranked within Pool_G ∩ stratum with k_s = min(100, ⌊\|Pool_s\|/3⌋); `min_k` 20 → descriptive, no claim (k_G and every k_s, incl. proxy); zero-frequency EN rows not rankable; `n_list` reported. Config: `P1_primary.strata_reported` (+ proxy), `strata_rule`, `min_k`, `below_min`, `en_zero_freq` | `tests/test_endpoint_branches.py` cases: pool 45 → k 15 descriptive; non_direct pool 300 → k_s 100; EN with 60 nonzero rows, k 100 → n_list 60; branches reachable on re-ranked strata | M7 | none |
| R3-6 | minor | accept | §7.3 and §13.1: NC2–NC6 bounds = max (min for NC6) over 5 toy replicates (`toy_seed + r`) at toy dimensions, written to `results/toy/nc_bounds.json`; provisional values become final if the timebox is exceeded (logged); both values recorded at `freeze-bounds`; "at real D dimensions" removed. Config: `toy_s13.nc_bound_replicates`, `nc_bound_seeds`, `nc_bound_rule`, `writes_nc_bounds`; `negative_controls_s7._comment` | `gate.py params`: `status_NC2_NC6 == "frozen"` after `freeze-bounds`, and each value equals the recomputed max/min over the 5 replicate outputs | M3b | **TOY**: run 5 replicates and write `nc_bounds.json` |

`configs/tiers.json` was updated per the PI ruling: `hours_remaining_at_v3` 16.0 and `reserve_h` 2.0. The milestone
budgets already summed to 14.0 h, so none was rescaled.
