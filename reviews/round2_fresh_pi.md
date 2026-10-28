# Round 2 — fresh external PI review of DESIGN v2

Reviewer: fresh external PI (has not read v1 or round-1 reviews). Inputs read: PROMPT.md, DESIGN.md v2,
configs/{design_params,primary_endpoints,tiers}.json, docs/part1_brief.md, RESEARCH_TREE.md, DATA.md, and the
contents of `external/` and `P1/data/` (directory listings only).

## Overall verdict

v2 covers every part of PROMPT.md: subtypes within and across lineages (§8), a single-feature benchmark (§5–6),
combinatorial screens (§11), confounders including culture, passengers, co-location and multiplicity (§5.3, §8.4–8.7),
clinical recovery with direction, actionability and population (§10), and limits (§12). The protections against R
leakage and the pre-registration are sound. Three things stop it from giving a convincing answer in 18 h:
(1) the allele-level mutation file has not arrived (`external/depmap_extra/` is empty), so most canonical T1 clinical
tests become `untestable` under the current rules, and those rules contradict the positive-control encoding;
(2) the plan is over budget (16.5 h + 2.5 h reserve = 19 h against about 18 h left) and carries several analyses that do
not change the scientific answer; (3) the primary replication endpoint will most likely be decided by expression
features and by the predictor's own hand-coded inputs, which a cancer-biology audience will discount. Items follow,
most important first.

---

### R2-FP-1 (critical): clinical ground truth collapses without the allele file, and the rules contradict themselves

**Problem.** §10.3 says that without OmicsSomaticMutations every allele-defined CGT test is `untestable: allele` and is
"never mapped to gene-level HOT". The file is absent (`external/depmap_extra/` is empty; DATA.md has no row for it).
Under that rule the core T1 sensitivity tests become untestable: BRAF V600 → BRAF (melanoma, NSCLC, tumour-agnostic),
EGFR activating → EGFR, KRAS G12C, PIK3CA → PIK3CA, IDH1, KIT, FGFR3 and AKT1 E17K. So does the refuted row BRAF V600E
CRC → BRAF, which puts the "≥ 3 refuted rows" requirement at risk. P2 would then rest mostly on CN, fusion and
LOF tests. Meanwhile §9.2 accepts `HOT:BRAF` as a "gene-level proxy for V600" for the gate. The same proxy cannot be
valid for the gate and invalid for the CGT.

**Change.**
(a) Today, ask the user once more for one file only, OmicsSomaticMutations.csv (26Q3). It is the single download that
unlocks alleles, TMB and spectra.
(b) Before `freeze-cgt`, pre-declare a fallback branch. When `data_status.json` shows the file absent, each
allele-defined test whose gene is in the hotspot matrix becomes an `allele_proxy` test with feature `HOT:<gene>` and
`proxy=1`. Exceptions where the gene-level call is not the clinical allele stay `untestable: allele`: KRAS G12C vs
other KRAS, EGFR T790M/C797S, KIT D816V vs exon 11, ERBB2 TKD.
(c) Report P2 on all testable tests including proxies (headline) and on non-proxy tests only (secondary). Both are
declared in `primary_endpoints.json` before the freeze.
(d) Make §9.2 and §10.3 cite the same proxy rule.

**Compliance check.** `clinical_ground_truth.csv` has columns `proxy` and `allele_level_required`. A script asserts:
every row with `proxy=1` has a `HOT:` feature, and its gene is not in `design_params.clinical_s10.allele_proxy_excluded`.
`primary_endpoints.json: P2_primary_clinical.sets` contains both `all_incl_proxy` and `non_proxy`. The `freeze-cgt` tag
is an ancestor of the first M4 result commit. `q4_summary.json` reports `n_testable_proxy` and `n_testable_allele`.

### R2-FP-2 (major): over budget; cut what does not change the answer

**Problem.** `tiers.json` has a 16.5 h plan plus a 2.5 h reserve against about 18 h left (§14 still says about 19 h).
M3b, at 3 h, holds the dictionary, covariates, CGT features, guard extensions, the toy simulation and the bound
simulation. That is the highest-risk milestone and it has no slack. The null runs as a 3 h "background job", but the
compute rule gives every foreground job 12 workers, so the two compete for the same cores. These analyses add cost
without changing any answer to PROMPT.md:
- XLIB with all four arms: the pred arms are already declared non-independent there, and EN-pred-clean needs an M1
  refit.
- The M3 refit for the Tier B attribution panel: refits, a reproduction check to r ≥ 0.999, and a 30-min diagnosis.
- The FaST-LMM top-200.
- A separate `simulate_control_bounds.py` step, which duplicates the toy null.
- Tier C as a whole.
- CR2 alongside CR1, when only 104 of 1,292 lines share a patient.
- A 1e-6 tolerance for the independent reimplementation of EN selection, which is fitted at tol 1e-3.

**Change.**
- XLIB: secondary, S1-obs and EN-obs only (the valid contrast). Drop EN-pred-clean.
- Attribution: Tier B = M1-proxy only, which is exact from the saved folds. The M3 refit moves to Tier C.
- Delete the LMM and all of Tier C.
- Derive the NC bounds from the toy run that §13.1 already requires. No separate script.
- Use CR1 everywhere, plus a 200-test CR1-vs-CR2 check that is reported, not gating.
- Independent check: S1 p-values (1e-6), BH counts and gate recall (exact), and EN selection sets for 20 genes
  (Jaccard ≥ 0.95).
- Null: 6 workers while foreground jobs use ≤ 6.
- New budget: milestone total ≤ 15 h with a reserve ≥ 3 h.

**Compliance check.** `tiers.json`: `milestone_total_h + reserve_h ≤ hours_remaining`, and `hours_remaining ≤ 18`.
The keys `lmm_top200`, `xlib_en_pred_clean_m1_refit` and the Tier C stages are absent.
`design_params.compute_s14.null_workers + foreground_workers ≤ 12`. `scripts/simulate_control_bounds.py` does not
exist, and NC bounds are read from the toy output. `budget_check.py` passes at each milestone.

### R2-FP-3 (major): P1 will be decided by expression features, not by assayable biomarkers

**Problem.** P1 ranks the top 100 over all of U. U includes EXPPAR, PROG (50 NMF factors), PSCORE and, for the EN arms,
own EXPDOSE. M5 is two-thirds expression-driven (M1 is expression-only, and M2 is mostly expression plus paralogs).
EN-pred will therefore rank expression programmes and paralog expression highly, and expression-to-dependency
relations replicate easily in R because they share the same Chronos fit. A positive "predictor adds value" could then
mean only "the predictor's expression readout replicates". That says little about the biomarkers the prompt asks for
(mutation, CN, methylation, signature) and nothing about clinical translation. Own EXPDOSE is already demoted (F2,
class 4), but only in S1.

**Change.**
- Make the P1 headline the genomic universe U_G (stage-G families: HOT, LOFBI, AMP, DEL, LOSS, FUS, SIL-DEL/MUT/METH,
  YABS, SIG, COMP). U_all becomes a pre-declared secondary with identical rules.
- Exclude own EXPDOSE from every arm's top-k for P1, matching S1/F2.
- Because the U_G pool may be small, fix `k_G = min(100, floor(|pool_G| / 3))` from the pool computed before the tag.

**Compliance check.** `primary_endpoints.json: P1_primary.universe == "U_G"`, with `secondary_universes` containing
`"U_all"` and the k_G rule present. `results/benchmark/p1_pool.json` (committed before `freeze-catalogue`) records
`|pool_G|` and `k_G`. A script asserts that no P1 top-k row has family ∈ {EXPDOSE, EXPPAR, PROG, PSCORE, SIL-other},
and that no row is own EXPDOSE in any list.

### R2-FP-4 (major): circularity through Part 1's hand-coded inputs

**Problem.** M2 is fed own CN, deep deletion, mutation and the top-3 paralogs' expression, deletion and mutation. M3
screens E, C, M and H. So the PAR and ADD-CN controls, and much of the paralog-SL catalogue, are **direct inputs** of
the predictor. When EN-pred recovers them, that partly re-reads Part 1's design choices (themselves taken from the
literature). This favours EN-pred both in the gate (PAR is the largest control family) and in P1. §4.4 fixes the
wording but not the endpoint.

**Change.**
- Report P1 and gate recall stratified by access class (§4.4): `direct` vs non-direct (potentially-direct,
  expression-derived, proxy).
- Pre-declare that "the predictor discovers determinants beyond its inputs" may be claimed only from the
  non-direct stratum, with the same bootstrap rule.
- Give the gate's PAR family a second recall computed on non-direct encodings where they exist, e.g. `LOFBI:SMARCA4`
  is not an M2 input but `EXPPAR` is.

**Compliance check.** `primary_endpoints.json: P1_primary.strata_reported == ["direct", "non_direct"]`, with a
`claim_rules.beyond_inputs` entry. `results/benchmark/p1_by_access.csv` exists. `results/claims.csv` contains no claim of
class `beyond_inputs` unless its evidence row has `stratum == non_direct`.

### R2-FP-5 (major): no deliverable a cancer biologist would read first

**Problem.** The outputs are mostly rates, CIs and class counts. A cancer-biology reader will judge the work on two
tables: (i) does it recover the known biology (each control or CGT pair, with rank, sign and R replication per arm),
and (ii) which new, replicated, cross-lineage subtypes does it propose, with the evidence for each. Neither is defined.
"Class 9 valid candidate" is not a ranked, evidence-carrying shortlist, and novelty is not defined at all.

**Change.**
- Add `scripts/build_candidate_cards.py`. It reads frozen outputs only and writes two tables:
  - `results/candidates/known_recovery.csv`: one row per positive control and per testable CGT test, with the rank in
    each arm, sign, R result and ceiling label.
  - `results/candidates/shortlist.csv`: class-9 pairs that replicated in R, ranked by (cross-lineage label, R one-sided
    p, number of lineages). Columns: g, d, access_class, lineages and per-lineage β, n_carriers D/R, confounding-delta
    flags, locus resolution, combinatorial label, CGT/actionability, and `novel` (not in positive_controls ∪ CGT ∪
    P1 biomarker_pairs ∪ oncogene_pairs; worded "not in curated lists", never "novel").
- REPORT.md opens with these two tables (top 15 shortlist rows) and 3 case studies chosen by rule: the top three
  shortlist rows.

**Compliance check.** Both CSVs exist, are git-tracked and are listed in `results/registry.csv` with status ok. A test
asserts the shortlist rows ⊆ (pair_classes class 9 ∧ R replicated). REPORT.md contains `<!--n:results/candidates/` tags.
A grep finds no bare "novel" applied to a gene in REPORT.md.

### R2-FP-6 (minor): eligibility vs controls and CGT targets is inconsistent, and misses are not decomposed

**Problem.** §5.1 defines U over eligible genes, but §2 and `controls_bypass_eligibility: true` analyse control genes
regardless. It is unclear whether non-eligible controls enter the arms' rankings. If they do, P1 and P2 are inflated.
CGT targets that fail eligibility (e.g. rare ALK) silently count as P2 misses, although that is a predictor limit, not a
readout failure. Claim (c) promises this decomposition but no table defines it.

**Change.**
- Non-eligible control genes are run in a separate `U_ctrl` and excluded from all top-k lists, the null subset and
  P1/P2.
- `q4_summary.json` decomposes every P2 miss as {ineligible: class 1/2/3, eligible but feature absent from 𝒟_g,
  eligible and ranked below k, ceiling data-null}.

**Compliance check.** A script asserts the genes in `results/benchmark/arms.csv` ⊆ `eligibility.csv[eligible]` and
`configs/null_genes.txt` ⊆ eligible. `q4_summary.json.miss_decomposition` sums to the number of misses per arm.

### R2-FP-7 (minor): the combinatorial screen lookup is a named prompt deliverable but sits in Tier B

**Problem.** PROMPT.md asks to "complement with combinatorial CRISPR screens where possible". The data are already
downloaded (DATA.md), and the lookup is cheap. As Tier B it runs only after every Tier A stage of M6, so it is the
first casualty of a late M6.

**Change.** Move `combinatorial_lookup` to Tier A (45 min) and add it to `never_cut`.

**Compliance check.** In `tiers.json`, `combinatorial_lookup` is in `tiers.A.stages` and in `never_cut`.
`results/combinatorial/lookup.csv` exists at the final gate.

---

No further changes are needed in my area. The remaining circularity (MD-strict double-dipping, R sharing the Chronos
fit, positive controls doubling as gate and CGT) is correctly declared and handled by the R-only confirmatory rule, the
gate-disjoint set and the wording rules.
