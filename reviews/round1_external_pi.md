# Review round 1: external PI (scientific value and feasibility)

Reviewer role: external PI, senior computational cancer biologist. Inputs read: PROMPT.md, DESIGN.md v1,
docs/part1_brief.md, docs/data_availability.md, and the implications sections of docs/lit/*.md. I also made small
read-only checks of Part 1 configs and the local environment:
- Eligible set = (variable ∨ selective) ∧ m5 tier ∈ {high, medium} = **3,062 genes** (1,331 high, 1,731 medium).
- Dev lineages: 18 have ≥ 20 lines.
- scikit-learn 1.5.2 **is** installed (data_availability.md §4 is wrong on this). igraph, leidenalg and shap are not.
- `external/probe/cmp_models.csv.gz` (Cell Model Passports, already local) has `msi_status` for 1,044 models (63 MSI)
  and BROAD_ID for 1,715.

## Overall verdict

The design covers every part of PROMPT.md: subtypes within and across lineages, the single-feature benchmark,
combinatorial screens, confounders, clinical validation (direction, actionability, population) and limits. The
confounder and clinical ground-truth sections are strong, and the pre-registration of expected signs and KO-drug
discordance classes is exactly right.

Two problems would stop a reader from being convinced.

1. **Circular benchmark.** As written, a model-derived candidate (MD) must pass the same observed-effect test as S1
   (§5.1 ii uses the §5.3 model, which is S1's model). So MD is, almost by construction, S1 filtered by the surrogate
   with a smaller BH denominator. "MD beats S1" would then mostly reflect the smaller test count, not anything the
   predictor learnt.
2. **The premise is only partly true for this predictor.** Part 1's accuracy is largely lineage and library: P1 is
   0.3175 overall and 0.211 within lineage × library. The predictor is also nearly expression-only (M1 alone reaches
   0.308). Determinants read from it will be dominated by expression programmes unless the design forces the
   genomic, clinically assayable readout to the front.

The plan is also over budget: §10 sums to 21 h against about 20 h left, including design. It needs a cut list.

Items are ordered by importance.

---

### R1-PI-1 (critical): The MD arm is nested in S1, so the benchmark is circular

**Problem.** MD requires (i) surrogate selection on predictions and (ii) significance in the §5.3 observed-effect model,
which is the S1 test. Any MD hit is therefore an S1-type hit with a smaller BH family. The comparison cannot show that
the predictor adds information. The design also lacks the ablation that separates "the predictor's knowledge" from
"sparse multivariate selection", because S1-pred is univariate only.

**Change.** Insert into §5.1 and §6:
- Define **MD-pure** as the benchmark arm. It is the (g, d, sign) set selected by the stability surrogate on fold-centred
  OOF predictions, with lineage and library group as unpenalised covariates. It is built **without reading observed
  effects**.
- The current MD (surrogate plus observed-effect test) is kept as **MD-strict**. MD-strict is the catalogue used for
  subtypes and clinical claims, not the benchmark arm.
- The benchmark is a pre-declared 2 × 2 of {univariate, sparse multivariate (same elastic net and stability rule)} ×
  {observed effects, OOF predictions}: S1-obs, S1-pred, EN-obs, EN-pred (= MD-pure).
- "The predictor adds value" is claimed only if EN-pred beats EN-obs (see R1-PI-2). "Multivariate adds value" is
  claimed only if EN-obs beats S1-obs.

**Check.**
- The surrogate script imports no CRISPR loader. A test asserts that `crispr.parquet` is not opened while MD-pure is
  built, using the guard's access log or by monkeypatching `pd.read_parquet` to record paths.
- `results/benchmark/arms.csv` has exactly the four arms with identical gene and 𝒟 universes. The script asserts
  equal (g, d) universes.

### R1-PI-2 (critical): No primary endpoint or decision rule for Q2 and Q4

**Problem.** §6 lists about 8 comparison metrics and §7.3 about 6 more, but none is primary. A positive story can be
picked post hoc, which PROMPT.md forbids ("never change metrics ... to make results look better").

**Change.** Add §6.0 "Primary endpoints (frozen at `freeze-design`)":
- **P1 (Q2):** replication rate in R of the top-k = 100 candidates per arm, ranked by stability frequency then |β|
  for EN arms and by p for S1 arms. "Replicated" = same sign and one-sided p < 0.05 in R. Computed over candidates
  with R power ≥ 0.5. Use a paired line-bootstrap 95% CI of the difference (1,000 resamples, patient-grouped).
- **P2 (Q4a):** recall of concordant, testable S/W CGT pairs at matched k = 100, 200 and 500, with the same CI.
- **Decision rule:** an arm is "better" only if the CI of the difference excludes 0 on P1. P2 is reported but is
  underpowered, since roughly 20–30 testable pairs are expected. Everything else is secondary and labelled as such in
  REPORT.md.

**Check.**
- `configs/primary_endpoints.json` is committed before the `freeze-catalogue` tag. `git log` shows it unchanged
  afterwards.
- The report script reads k, the ranking rule and the CI settings only from that file and fails if the file's hash
  differs from the one recorded at the tag.

### R1-PI-3 (critical): Budget exceeds the time left; cut list and compute-feasible statistics

**Problem.**
- §10 totals 21 h against about 20 h left.
- Several steps are infeasible as specified:
  - Stability selection with a λ-CV inside each of 50 subsamples × 3,062 genes × about 2,000 features is roughly
    10⁶ elastic-net path fits, many hours on 12 cores.
  - The "200-permutation end-to-end null" in §5.4 would multiply that by 200.
  - §4.2 needs a full M1/M2/M3 CV refit (1–2 h, no saved folds for M1, M2 or M3) plus a 1e-4 reproduction check.
    LightGBM thread nondeterminism can fail that check and trigger an open-ended "stop and investigate".
  - Leiden is not installed.

**Change.** Add §10.1 "Tiers and timeboxes" (pre-declared; anything not finished is reported as "not done: time"):
- **Tier A (must):**
  - 𝒟, S1, and MD-pure/MD-strict with fixed-λ stability selection: λ chosen once per gene by grouped CV on all D, then
    50 half-samples at that λ (Meinshausen–Bühlmann).
  - The confounder audit, the gate, CGT recovery and direction, replication in R, and the limits classes.
- **Tier B:**
  - Direct attribution on a ≤ 300-gene panel (positive controls, CGT target genes, top MD-strict genes):
    - M1 component: exact linear attribution from the saved `t_E_status_24b46e2a` folds (no refit; labelled
      "M1-proxy", since T ≈ M1, P1 0.3082 vs 0.3081).
    - M3 component: refit only for the panel genes with LightGBM `pred_contrib`.
    - Reproduction check: per-gene r ≥ 0.999 against Part 1 OOF, not max |Δ| ≤ 1e-4. On failure, report the
      shortfall and continue.
  - Combinatorial lookup (see R1-PI-11).
  - The LMM sensitivity on the top 200 S1 hits, not 1,000.
- **Tier C (only if ahead of schedule):** unsupervised clustering (networkx Louvain, since Leiden is not installed),
  TCGA prevalence via cBioPortal, and ClinicalTrials.gov harvesting.
- **Null calibration:**
  - S1: 200 within-lineage permutations, vectorised by residualising on covariates (Frisch–Waugh).
  - MD-pure: 20 within-lineage permutations of 𝒟 on a random 200-gene subset (seed list), with the fixed-λ stability
    pipeline. Report empirical FDR with its binomial CI.
- **SIL family:** restrict to OncoKB cancer genes ∪ paralog/partner-prior genes, as for LOF. Otherwise SIL alone can
  reach thousands of features.

**Check.**
- `results/timing/stages.csv` logs wall time per stage. A script compares it against the timeboxes in
  `configs/tiers.json` and lists any Tier B/C stage started before every Tier A stage was marked done.
- `configs/dictionary_v1.json` feature count is ≤ a declared cap (e.g. 3,000).

### R1-PI-4 (major): Expression families will crowd out the genomic, clinically assayable determinants

**Problem.**
- The predictor is almost entirely expression-based. In a joint elastic net, EXP-own, EXP-paralog and 50 PROG factors
  will absorb the variance and suppress correlated genomic features. Examples: SMARCA4 LOF vs low SMARCA4 expression;
  ERBB2 AMP vs ERBB2 expression.
- Nearly all CGT biomarkers are genomic, so clinical recall would be depressed by construction.
- §7.1 also uses "ESR1 expression → ESR1" as a positive control, while §8 class 4 calls EXP-own determinants
  tautological. The two rules contradict each other.

**Change.**
- Make the surrogate hierarchical. Stage 1 uses genomic and signature families only (HOT, LOF, AMP/DEL loci, FUS,
  SIL-DEL, SIG, COMP). Stage 2 uses expression families (EXP-own, EXP-paralog, PROG, SIL-other) on the stage-1
  residual.
- Report fidelity_genomic and Δfidelity_expression separately.
- An EXP-own determinant is a biomarker (not class 4) only when the gene is a CGT biomarker measured clinically at
  protein or RNA level. Pre-list these: ESR1, ERBB2, AR, PGR, FOLR1. Otherwise it is class 4 unless backed by AMP
  or SIL.
- PROG determinants are labelled "expression programme (no clinical assay)" and are excluded from Q4b actionability
  numerators.

**Check.**
- The surrogate output has columns `fidelity_genomic` and `dfidelity_expr`.
- A test plants a genomic feature with a correlated expression proxy in the toy data (§9) and asserts that the genomic
  feature is recovered in stage 1.
- A script asserts that every §7.1 control's expected determinant family is one the §8 rules count as a valid
  biomarker.

### R1-PI-5 (major): "Accurate prediction" must mean accurate beyond lineage and library, and library must be a surrogate covariate

**Problem.**
- Eligibility uses Part 1's pooled tiers, which reward lineage and library structure.
- The surrogate adjusts for lineage but not for library group. Library is confounded with lineage and screen length,
  and expression partly encodes culture.
- A gene can be "accurately predicted" only because the model recognises lineage or library. Reading determinants from
  such a gene is the circularity PROMPT.md's premise risks.

**Change.**
- In §2, add a pre-set second eligibility criterion computed on D only: the within-(lineage × library) Pearson r of
  fold-centred OOF predictions vs observed effects, centred within strata with ≥ 5 lines, must be ≥ 0.10.
- Genes failing it go directly to §8 class 3 ("accurate only via lineage/library"), with their counts reported.
- Add library group as an unpenalised covariate in §4.1, beside lineage.

**Check.** `results/eligibility.csv` has `r_within_lin_lib`. A script asserts that no gene below 0.10 appears in any MD
list and that every such gene has class 3 in `results/limits/gene_classes.csv`.

### R1-PI-6 (major): The positive-control gate must be symmetric, tamper-proof and runnable without pending downloads

**Problem.**
- WRN–MSI, the strongest forward-looking control, depends on OmicsGlobalSignatures, which may never arrive.
  M_burden is not MSI.
- The gate is evaluated only for MD. If the gate fails, "stop and diagnose" invites tuning on the same controls that
  §6 and §7.3 then report as recall. That is double use of the controls.

**Change.**
- **MSI fallback, pre-declared:** Cell Model Passports `msi_status` (local probe file, CC BY 4.0; joined via BROAD_ID;
  add it to DATA.md). This is used if OmicsGlobalSignatures is absent by M4. If both exist, report their agreement.
- **Symmetric gate:** evaluate the gate for all four arms (R1-PI-1).
- **No tuning after the gate:** the only changes allowed after gate evaluation are bug fixes, each demonstrated by a
  failing toy test that passes afterwards, logged in RESEARCH_TREE.md and applied to all arms.
- Recall tables report gate controls in a separate row block. The headline CGT recall is computed with and without
  the gate pairs.

**Check.**
- `results/gate/gate_report.csv` has one row per (control, arm).
- Every commit touching `src/` after the `gate-evaluated` tag must also touch `tests/` (checked with `git diff --stat`).
- The CGT table has an `in_gate_set` column, and the report shows both recall numbers.

### R1-PI-7 (major): The CGT recovery definition makes direction accuracy tautological; the CGT scope needs a cap

**Problem.**
- §7.1 defines "recovered" as present with the expected sign. If Q4a uses the same definition, direction accuracy is
  100% by construction.
- The CGT build (CIViC + MOAlmanac + FDA scrape + CTG + curated mapping) is open-ended for a 3 h milestone.

**Change.**
- **Recovered** = (g, d) in the arm's list, regardless of sign.
- **Direction accuracy** = the fraction of recovered pairs with the pre-registered sign. Also report it over all
  testable S/W pairs with nominal p < 0.05 in D, which is method-independent.
- Resistance and N-class pairs: report the fraction "recovered with a sensitivity direction". This is the false
  positive rate.
- Scope: start from the 45-pair curated table in docs/lit/clinical_ground_truth.md. Extend it automatically only with
  CIViC Predictive A/B "Supports" rows and FDA CDx rows that map to a 𝒟 feature and a library gene, within a 1.5 h
  timebox. Use CTG only to attach an NCT ID to trial-stage pairs.
- **Q4c population**, pre-declared metric: for each recovered pair, the odds ratio of indicated lineages among
  subtype lines vs all carriers. Tumour-agnostic pairs must be cross-lineage under §5.2.

**Check.**
- `configs/clinical_ground_truth.csv` has `expected_sign` and is frozen at `freeze-cgt` before the join. The commit
  order is checked.
- The Q4 script outputs `recovered`, `sign_correct` and `false_sensitivity_rate` as separate columns.
- A unit test with a planted wrong-sign pair asserts recovered = 1 and sign_correct = 0.

### R1-PI-8 (major): Test PROMPT.md's cross-lineage claim instead of assuming it; fix the low-power bias of the cross-lineage rule

**Problem.**
- PROMPT.md asserts that cross-lineage determinants are stronger translational candidates. The design labels
  determinants cross-lineage but never tests whether the label predicts anything.
- The rule (meta significant + nominal p < 0.05 in ≥ 2 lineages, each with ≥ 3 carriers) mostly measures biomarker
  frequency, not consistency.

**Change.**
- Add Q1b: a pre-declared comparison of cross-lineage vs lineage-restricted MD-strict determinants on (a) replication
  in R and (b) CGT membership, by Fisher test with ORs. A null result is reported as a finding.
- Redefine **cross-lineage** as: ≥ 2 lineages each with power ≥ 0.5 for the pooled effect; random-effects meta
  p < 0.05; and same sign in ≥ 80% of powered lineages (I² reported).
- Determinants without ≥ 2 powered lineages are "lineage-unresolvable", not "pooled-only".

**Check.** `results/subtypes/cross_lineage_test.csv` exists with ORs and CIs. A script asserts that every "cross-lineage"
label has n_powered_lineages ≥ 2.

### R1-PI-9 (major): R is underpowered for rare biomarkers; add an independent-library replication inside D

**Problem.**
- With 255 R lines, a determinant with 5–15 D carriers has about 1–4 carriers in R. Most fusion, rare hotspot and
  locus determinants will be untestable.
- §5.4's "discover in Avana-only, check sign in the rest" is ambiguous. If discovery used all of D, the check is not
  independent.

**Change.**
- Pre-declare a second replication: re-run S1-obs and MD-pure (fixed λ, no extra stability runs) on Avana-only D lines
  (lib_group 0, n = 621). Replicate on non-Avana D lines (n = 401: KY/Sanger-derived, Humagne and mixed), using the
  same top-k and replication definition as R1-PI-2.
- Report it beside R as "cross-library replication". It is a different lab and library, which is the strongest
  available independence.
- Pool R with the cross-library set only in a clearly labelled secondary analysis.

**Check.** `results/benchmark/replication_xlib.csv` exists. A script asserts that the discovery-set ModelIDs and the
replication-set ModelIDs are disjoint and match lib_group.

### R1-PI-10 (minor): The "accurate-dependent" threshold ignores dependency prevalence

**Problem.** "Prediction below the dev 20th percentile and observed < −0.5" mislabels lines whenever more than 20% of
lines are dependent: correct predictions become "unexplained" and feed §8 spuriously. It also mixes fold offsets
unless percentiles are taken on fold-centred predictions.

**Change.**
- Line (ℓ, g) is accurate-dependent if observed effect < −0.5 and the fold-centred prediction ranks within the n_g most
  dependent predictions, where n_g = the number of observed dependent D lines for g.
- Use the −1.0 sensitivity threshold as already declared.

**Check.** On toy data with prevalence 40%, a test asserts that a perfect predictor labels all dependent lines
accurate-dependent.

### R1-PI-11 (minor): Scope the combinatorial lookup to what can validate, and timebox it

**Problem.**
- Five sources are named. Several need manual or bioRxiv supplementary downloads.
- Combinatorial screens can only validate paralog-buffering candidates. Inter-screen precision is only about 0.4.

**Change.**
- Use at most two scriptable sources (Parrish GSE178179 plus one consensus table, e.g. Esmaeili Anvar 2024 figshare),
  after first checking whether local BioGRID ORCS 2.0.18 already contains them. Timebox: 1.5 h.
- Gold standard: hit in ≥ 2 screens or ≥ 2 lines.
- Report the validation rate of MD-strict paralog candidates vs S1-obs paralog candidates vs a background of screened
  pairs matched on sequence identity and family size. Present it as supporting evidence with N stated, never as a
  primary endpoint.

**Check.** `results/combinatorial/lookup.csv` has columns `source, n_screens_hit, n_lines_hit, matched_background_rate`,
and `configs/external_sha256.txt` lists each source file.

### R1-PI-12 (minor): Make the limits answer quantitative about the predictor itself, and name unassayed determinant types

**Problem.**
- §8 classifies genes but never states where the predictor-based readout fails relative to direct association. That
  is the most informative limit for this project.
- PROMPT.md names methylation and mutational signatures. 𝒟 has neither, unless RRBS arrives, and that gap is silent.

**Change.** Add two classes to §8, assigned per (g, d), and report their counts:
- **"predictor-blind":** S1-obs hit, replicated (R or cross-library), not selected by MD-pure. Example: features absent
  from the predictor's inputs, such as MSI, which it sees only through expression.
- **"predictor-only, not replicated":** MD-pure hit that fails replication while powered.

Also add a §3 note:
- Methylation is represented only by the SIL-other expression proxy. It is validated by checking MLH1 SIL-other
  against CMP MSI status as a known-positive.
- Mutational signatures (APOBEC, POLE, HRD) are "not assayed". LoH fraction is the HRD proxy, if downloaded.
- REPORT.md lists both as limits.

**Check.** `results/limits/pair_classes.csv` contains both classes with counts. A script asserts that every
dictionary family named in PROMPT.md (mutation, CN, methylation, mutational signature, expression signature) appears
in `configs/dictionary_v1.json`, either with a feature count or with status `not_assayed` and a reason.

---

**Cut summary (to fit about 20 h):**
- Drop the full M1/M2/M3 refit, giving a ≤ 300-gene attribution panel.
- Leiden becomes optional Tier C Louvain, and TCGA prevalence moves to Tier C.
- LMM: top 1,000 → top 200.
- MD permutation null: 200 → 20 permutations on 200 genes.
- Stability selection uses fixed λ.
- Combinatorial: 5 sources → ≤ 2, with a timebox.
- CGT: start from the 45 curated pairs plus a capped automatic extension.

Suggested M2–M8 budget: 2 / 2.5 / 4 / 2.5 / 2.5 / 1 / 3.5 h, plus a 2 h reserve.
