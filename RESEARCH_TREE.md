# Research tree — Dependency subtypes and molecular biomarkers from the Part 1 essentiality predictor

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (latest first)

### 2026-10-29 ~09:50 — SUSPEND (user suspending the machine); resume point
Committed through 3833bd0 (REPORT.md draft with placeholders, docs/effective_rules.md, check_report_numbers.py), and
nothing is uncommitted. The single job running is the fixed-threshold null, started at 08:49
(`scripts/run_null.py`, 3 workers, log cache/run_null3.log). It has no checkpoint. If VS Code stays open, sleep pauses it
and it continues on wake. If its process is gone on resume, re-run it (about 3.2 h at 3 workers; free RAM first).
Interrupted earlier and not yet re-run: the confounder audit and cross-lineage, neither of which has output.

Unverified drafts: REPORT.md and docs/effective_rules.md, which still need the strict number check and the citation fixes.

Next on resume:
1. Check the null process and its log.
2. Run the confounder audit with 1–2 workers, then cross-lineage, run_m5_audits, run_q4 and run_component_readout.
3. Run run_limits and build_candidate_cards.
4. Re-run the toy suite (src changed after 04aa982).
5. After the null: run_xlib.
6. build_p1_pool, the pseudo-R dry run, then tag `freeze-catalogue`.
7. replicate.py, a single R run.
8. Re-run limits and cards, then fill REPORT.md, run the strict number check and the red-team review, prune LESSONS,
   and do the final commit.

Agents: the author (ab46581f) and the writer (ab7cb7af) are idle and can be resumed via SendMessage. Two open requests
for the author: registry rows from every results script, and making the number checker accept DESIGN.md and
RESEARCH_TREE.md as text sources for effective_rules.md.

### 2026-10-29 08:40 — INTERRUPTION (VS Code closed overnight); resumed
Everything up to 19e4b31 is committed and intact. Five running jobs were lost, and none of them had written any output:
- the fixed-threshold null, which had been running since 03:38 and has no per-draw checkpoint;
- cross-lineage;
- the confounder audit;
- the writer (effective_rules.md and the REPORT.md draft);
- check_report_numbers.py.
Re-run at 08:40: the null (background), then the confounder audit; the writer and author agents were resumed via
SendMessage. The machine has about 3.5 GB of RAM free before any job starts, so the memory guard cuts workers and the
null will be slower than projected. Remaining order: confounder audit → cross-lineage → M5 audits → Q4 → limits/cards →
toy suite re-run (src changed after 04aa982) → null → XLIB → freeze-catalogue (P1 pool, r_requests, dry run) →
single R replication → limits/cards re-run → REPORT.md, number check, red-team review, LESSONS prune, final commit.

### 2026-10-29 ~00:35 — M3b milestone: bounds frozen (tag `freeze-bounds`), engine verified on the toy
The Tier A engine (S1 with CR2 and a carrier df, a two-stage stability elastic net, the four arms, the MD lists, the
fixed-threshold null) was written by one agent. Two others checked it on the toy simulation and with process checks.
Neither found an engine bug in the final version, but the checks changed four things before any real association
result:
- The calibration test showed CR2 + ν is anti-conservative at n = 255, so the pre-declared df fallback now applies.
- Power and test df were made consistent.
- Expression programmes are now computed without the held-out lines.
- The catalogue became MD-global (A7). MD-strict, the frozen catalogue, had a true FDP of 0.16–0.26 on all five
  ridge-predictor toys and a false hit in 30% of pure-null toys. The double dip causes this: the observed screen re-uses
  the dev-set noise the ridge predictions absorbed. MD-global had FDP ≤ 0.025 and no pure-null hits, at the same recall.

The final toy suite gives 19 pass, 3 fail. The fails are the fold-offset positive control (the bound was derived on a
different list), ridge MD-strict, and pure-null S1-pred (4 marginal hits). All are reported, not changed. The NC bounds
come from 5 toy replicates under rulings fixed in advance: NC2a 0.0011, NC2b 1.2, and the provisional values for NC3,
NC5 and NC6. `configs/frozen_hashes.json:bounds` covers the dictionary, controls, null genes, subsamples, data status,
the CGT files and the pipeline constants. Next is M4 on real D: the four arms + MD lists, the confounder audit, gate
attempt 1, then the null in the background, XLIB and the independent check. Time used ≈ 7.5 h of 24.

### 2026-10-28 ~19:40 — M2 milestone: design frozen (tag `freeze-design`, DESIGN v3.2)
The design went through three review rounds: 7 reviewers / 84 items (R1), a fresh PI plus statistics and implementability
reviewers / 23 items (R2), and one fresh reviewer / 6 items (R3). Every item has a recorded decision
(reviews/round{1,2,3}_response.md). The core design: read determinants out of the Part 1 predictor with a stability
elastic net fitted to its out-of-fold predictions (EN-pred), and compare it in a 2×2 benchmark against the same model on
observed effects and univariate confounder-adjusted tests. The primary endpoint is replication of the top genomic
candidates in the 255 never-analysed Part 1 test lines (opened once), stratified by whether the determinant was a
hand-coded predictor input. Clinical validation uses a pre-registered CIViC/MOAlmanac/FDA test table with signs and
partitions fixed before any dependency data are joined. Frozen parameters are in configs/ and hashed in
configs/frozen_hashes.json.

All user-requested DepMap extras arrived, including RRBS methylation (amendment D7). The core library and the external
data are committed. Three build tracks are running (dictionary/eligibility/controls; clinical ground-truth table;
toy simulation by a separate verification agent). Next: finish M3b, run the toy at real settings, tag `freeze-bounds`,
then M4 (the four arms, MD-strict, gate, null). Time used ≈ 2.6 h of 24; the frozen schedule budgets 14 h + 2 h reserve.

### 2026-10-28 — SUSPEND (machine suspended by the user mid-M2/M3b) — resumed 19:02; all 4 agents resumed via SendMessage
Done and committed: M0 setup; M1 literature (LITERATURE.md, docs/lit/*, data availability, Part 1 brief); DESIGN v1 → v2
after round 1 (84 items, reviews/round1_response.md); round-2 reviews (reviews/round2_*.md, decisions D5–D6); core
library src/biomark/{paths,p1,guard,gitutil}.py with 31 passing tests; input checksums verified; external clinical and
combinatorial data fetched; user-supplied DepMap extras registered (configs/depmap_extra_sha256.txt). No analysis
result exists yet; the replication set R has never been read; no freeze tag exists.

Interrupted (uncommitted files left on disk; treat as unverified drafts, re-check before use):
(1) DESIGN v3 + reviews/round2_response.md + configs/{design_params,primary_endpoints,tiers,seeds}.json edits — the
editor was applying D6; check that round2_response.md has a row for all 23 items and configs match DESIGN §F, then commit
and tag `freeze-design`. (2) M3b dictionary/eligibility/controls track: src/biomark/{dictionary,covariates,eligibility,
molecular,genesets,datastatus}.py, scripts/data_status.py, configs/{build_params.json,positive_controls_spec.csv,
fusion_known_positives.csv}, results/data_status.json, p1.py edit — re-dispatch with "finish and test" against v3.
(3) CGT build (2nd attempt, two-step; docs/cgt_handover.md): scripts/fetch_cgt_records.py only — re-dispatch.
(4) Toy simulation (verification agent): configs/toy_answer_key.json only — re-dispatch; must stay a different agent from
the pipeline author. Next after these: freeze-bounds → M4 (four arms, MD-strict, gate, null).

### 2026-10-28 — M0: setup
Project started. Prompt saved verbatim (PROMPT.md), `.gitignore` set (no raw/derived bulk data), tracking files created.
Part 1 (../project, DepMap 26Q3, headline model M5 = average of kernel ridge M1, paralog stage M2, LightGBM M3) is read
in place and never modified.

Next: parallel literature surveys and a technical brief of Part 1, then DESIGN.md and up to three review rounds.

## Decision log
- D0 (M0): Part 1 is read-only; all caches/outputs live in this folder (`cache/`, `results/`). Inherited process rules
  from Part 1 LESSONS are seeded into LESSONS.md and marked as inherited.
- D1 (M2, R1): benchmark becomes a 2×2 {univariate S1, stability elastic net EN} × {observed, OOF-predicted};
  EN-pred (= MD-pure) is built without observed effects; "predictor adds value" only if EN-pred > EN-obs on the frozen
  primary endpoint (replication in R of top-100, power ≥ 0.5). Reason: v1's MD was S1 with a smaller BH family (R1-PI-1, ST1, PW1).
- D2 (M2, R1): covariates equal to / collinear with the determinant (|ρ| > 0.7) are dropped per test; self-null
  features removed; 9p21 handled by locus groups + resolution test, not merging; LOSS (CN < 0.75) family added (R1-ST3, PW2, DA1, DA2).
- D20 (M4, 04:40): [!] **Gate attempt 2 also failed, 21/39**, after the U_ctrl bug fix (41 genes; toy test
  test_gate_u_ctrl.py fails on 0346d69 and passes after). All 7 non-eligible gate controls are missed by every arm. These
  are genuine misses with no further bug:
  - SMARCA2 × SMARCA4-LOFBI: β −0.48, p 1e-5, q 0.15.
  - ARID1B × ARID1A-LOFBI: q 0.099. Both fail on the BH family size of about 123k tests.
  - ESR1 own-EXPDOSE: q 0.0008, but the own-EXPDOSE feature belongs to F2. The gate counts F1 lists, so this control
    could never be recovered under the list definitions, a specification gap reported as such; it is not re-scored
    after the fact.
  - ENO2, PAPSS1, ME3 × LOSS: β ≈ 0. LOSS (CN < 0.75) does not capture the homozygous deletions these pairs need.
  - MYCN AMP: p 0.01, absorbed by lineage.
  - Eligible controls alone: 21/32 (0.66).
  Per outcome_actions the project proceeds: the report states "gate failed at attempt 1 (21/39) and at attempt 2 after
  the U_ctrl fix (21/39)", and the CGT headline is restricted to gate-disjoint tests. S1-pred recovers the most controls
  (29/39), then MD-strict, S1-obs and EN-pred (21), MD-global (17) and EN-obs (15).
- D19 (M4, 03:40, first real results):
  - **Arms.** The four arms ran on D: 2,193 genes, 26.2 M rows, 45 min with S1 at 8 and EN at 5 workers. Selected:
    S1-obs 3,142, EN-obs 738, EN-pred 50,728, S1-pred 276,557. The catalogue MD-global has 1,380 pairs; MD-strict has 11,662.
  - [!] **Gate attempt 1 failed: 21/39** (0.54, CI 0.37–0.70). ADD-CN 1/5, SIG 0/2 and W 1/3 are below 0.5.
  - **Bug search, per outcome_actions.** U_ctrl (§2, the non-eligible control and CGT-target genes) was never built or
    read. As a result, SMARCA2, ARID1B, MYCN and ESR1 had no rows, which is a pipeline omission. This is fixed as a §9.3
    bug fix with a toy test, leading to attempt 2. The report will state "gate failed at attempt 1 (21/39)".
  - **Misses that are not bugs.**
    - EGFR-activating: 8 carriers, so the triggered df fallback gives df 7 and p = 0.13 for β −0.63. This is the cost
      of DEV-CAL-CONSERVATIVE.
    - MAT2A × MTAP-del: q = 0.15.
  - **Other steps.** The null runs in the background; the independent §13.3 check has started (separate agent, no
    engine code read).
- D18 (M5/M6 code, ~01:45, before any real catalogue exists): the M5/M6/M8 Tier A scripts are built and toy-tested
  (24/24). The implementation choices that DESIGN §8/§11/§12 leaves open are fixed in configs/m5_m8_params.json and hashed
  under frozen_hashes `actions` before any real M5 output. Among them, the combinatorial per-line hit rule (Cohen's
  d ≥ 0.8 and Δmean ≤ −0.5, only within studies that call the pair a hit) is an agent choice, not a design value, and
  will be stated as such in the report. Esmaeili 2024 is the single harmonised source for Dede, Parrish, Thompson,
  CHyMErA and Ito, so no screen is counted twice.
- D17 (M4, ~01:15, compute; no real result exists): the machine has 15.7 GB of RAM, but the full arms table is about
  26 M rows (≈ 25 GB as one object-dtype frame), and EN workers hold ≈ 0.8–1.7 GB each.
  - The arms run in gene blocks with a parquet dataset and slim global steps (blocked = unblocked, tested).
  - Workers are capped by a memory guard: S1 8, EN 6, null 6.
  - The null runs after the arms instead of beside them.
  - The frozen `fwl_cache_max_entries` is unchanged.
  - Walls will exceed the tier timeboxes; this is logged by budget_check as lateness, not as a design change.
- D16 (M3b end, 00:35): NC2b = 1.2 (MD-strict hits per low-tier gene) is applied as the mechanical result of ruling R1.
  R1 was fixed before any ridge value existed, so it is not revisited. MD-global's low-tier rate is reported next to the
  provisional 0.1 for information. NC values are recorded in frozen_hashes `bounds_values`; design_params.json is not
  edited. The calibration fallback is kept although 4 few-carrier cells are conservative (DEV-CAL-CONSERVATIVE).
  Tag `freeze-bounds`.
- D15 (M4, ~23:30): [!] The catalogue becomes MD-global (EN-pred ∩ S1-obs full-family BH, same sign); this is amendment
  A7, made before any real association result. On the engine's ridge toys MD-strict FDP was 0.228, 0.255 and 0.237 and
  MD-global's was 0, 0 and 0.013. On the pure nulls P(any hit) was 0.30 for MD-strict and 0.00 for MD-global. Recall was
  identical on all toys. MD-strict stays as a reported secondary; the gate, NC2b and the P1/P2 contrasts are unchanged.
- D14 (M4, ~23:15, NC bound rulings before any ridge NC2 value or real result):
  - R1: NC2 (low-tier false hits) is an error-type control, so per A6 its toy value is the max over 5 ridge-toy replicates.
    On the smooth toy EN-pred selects nothing on low-tier genes, so the value is degenerate (≈ 0).
  - R2: a toy NC statistic with zero between-replicate variance at the edge of its range, or with a zero underlying count,
    has no resolution. Its provisional value becomes final (DEV-TOY row), as for NC3/NC5. NC6 (1.0 in 5/5) → 0.7 stays.
  - Interim diagnostics on 8cffc91 found no engine bug.
    - Fold artefact: the EN frequency is 0.36–0.50 without centring and ≤ 0.04 with it. The bound was derived from S1-pred
      BH, not from EN.
    - S1-obs smooth: false hits 10/4/8/3/9 over 5 replicates; the null gives 5.5/draw. This is variability plus a
      denominator mismatch in the test (FDR_emp counts neutral hits, FDP_true does not).
    - [!] Pure null P(any hit): MD-strict 0.30, S1-pred 0.25. The candidate secondary list MD-global (EN-pred ∩ S1-obs BH)
      gives 0.00, with FDP 0 and recall 0.87–0.88 on the smooth toys (MD-strict: FDP ≤ 0.046, same recall). Whether it
      becomes the error-controlled catalogue (amendment A7) waits for the ridge-toy numbers.
- D13 (M4, ~23:00, process-guard findings before freeze-bounds, no association result yet):
  - (i) The CR2 + ν test is anti-conservative at n = 255 with heteroskedastic carriers (size 0.069–0.073 at α 0.05). This
    triggers the pre-declared §5.3 fallback, binary d → df = carrier clusters − 1, applied in every S1 test and in power
    and logged in deviations.csv. It is not a choice.
  - (ii) Power df ≠ test df, a bug; it is being fixed.
  - (iii) NMF scores were not D-only invariant (8e-4). D, R and set A will be transformed separately and the dictionary
    rebuilt before freeze-bounds.
  - (iv) build_cgt.py gene literals move to config, with a byte-identical CGT required.
  - (v) String exemptions from the leakscan go through an explicit allowlist (configs/leakscan_allow.json); `literal-ok`
    covers numbers only.
  - (vi) [!] The cross-lineage rule (§8.3) labelled a one-lineage plant cross-lineage in 1/5 toy seeds when the null
    lineages shared its sign by chance. This is kept as a reported limitation of the label; the rule is unchanged.
  - (vii) Six files that §13.4 lists under freeze-design were created later in M3b. They are hashed in the `bounds` section
    instead.
  - (viii) M3b builders wrote no timing records, so their wall time is unknown. This is logged as a deviation, not
    reconstructed.
- D12 (M4, 21:50, toy run of the real engine before freeze-bounds): 14/19 toy checks pass, 5 fail; no bound or answer
  key changed. [!] Ridge toy: MD-strict FDP_true 0.228 > 0.15, the double dip the answer key predicted (§4.1). The
  D-stage observed screen re-uses the D noise that the ridge OOF predictions absorbed, so the catalogue's D-based error
  rate is not a valid FDR claim with a ridge-type predictor (Part 1 M1/M2 are ridge); R replication (P1) is the arbiter.
  The other failures (fold-offset positive control, smooth S1-obs FDR_emp −0.004 beyond tolerance, pure-null S1-pred
  0.25, gate crash, fixed) go to an independent verifier for diagnosis; any resulting amendment (A7) will be written
  before freeze-bounds and before any real association result.
- D11 (M3b, toy): [!] with a ridge-type predictor, prediction-based lists have true FDP ≈ 0.7–0.8 that the
  permutation null cannot see (ridge spreads a determinant's signal onto correlated features). Response: MD-strict stays
  the catalogue, EN-pred FDR_emp labelled not calibrated, R replication arbitrates; amendment A6 (also eligibility strata
  < 5 lines excluded; pure-null count rule).
- D10 (M3b): dictionary track outcomes — shared 𝒟 2,965 features (cap 3,000; AMP/METH minimum carriers raised to 20,
  SIL-other to 10 by the frozen cap rule); 2,193 eligible genes; 39/40 controls gate-eligible. Rulings: ALLELE/CLIN features
  built into 𝒟 by the dictionary owner from the CGT rules; canonical name ALLELE:EGFR_activating; allele primaries only for
  BRAF and EGFR (R3-2), KRAS/NRAS/PIK3CA stay HOT (§9.2); KMT2A added to the 5′-driver fusion exception (DEV-DICT-1);
  SPEC_MMR excluded (failed its pre-declared ρ ≥ 0.5 sanity check vs MSIScore: 0.27); RRBS covers 453/1022 D lines.
- D9 (M3b): ER+ (and HR+/HER2−) operationalised as ESR1-expression rule AND PAM50 luminal call, decided before any
  dependency data were joined: the expression-only rule over-called (29/65 breast lines ER+, 14 non-luminal). ALK cap of 10
  was a pre-data estimate; 11 real carriers kept (DEV-CGT-1). KMT2A-r = any KMT2A-5′ fusion (no memory-based partner list).
- D8 (M2 milestone): design frozen at v3.2 after 3 rounds (113 items total, none open). Further changes only as logged
  amendments with reason; parameters cannot change after results are seen.
- D7 (M2): the user supplied CCLE RRBS promoter methylation before the first dictionary build, so methylation is in scope
  (families METH and SIL-METH, U_G, access class proxy); coverage is limited to RRBS-profiled CCLE lines.
- D6 (M2, R2): accept round 2 nearly in full: HC2/CR2 with carrier-based df (small-carrier size inflation, R2-ST1);
  fixed-threshold null (R2-ST2); EN in-fit carrier floor 2 (R2-ST4); P1 headline on genomic universe U_G with access-class
  strata (R2-FP-3/4); two freeze tags and a pre-listed R request file (R2-IM1/2); budget cuts (LMM, Tier C, M3 refit).
- D5 (M3b): CGT agent #1 halted before writing files; re-dispatched in two steps (curated mapping → mechanical build).
  Clinically-refuted rows may cite a PubMed-verified PMID only (flag pmid_only); all other rows need CIViC/MOA/FDA/NCT IDs.
- D4 (M2): M3b implementation of design-stable parts starts in parallel with round 2 to save wall-clock; round-2 changes are applied as diffs.
- D3 (M2, R1): compute-feasible statistics (fixed-λ stability selection, SIS K=300; joint-row permutation null B=20 on a
  pre-drawn 300-gene subset) and Tier A/B/C time boxes (R1-PI-3, AC1, AC2, ST2).

## Tree
- [x] M0 Setup — PROMPT.md, .gitignore, DATA.md, LESSONS.md, RESEARCH_TREE.md
- [~] M1 Literature and data availability → LITERATURE.md
  - [x] docs/part1_brief.md — OOF/test preds cached; M1/M2/T linear (exact attribution), M3 TreeSHAP; 5 write hazards into P1
  - [x] docs/lit/methods_prior_work.md, clinical_ground_truth.md, combinatorial_confounders.md; docs/data_availability.md
  - [x] LITERATURE.md synthesis (writer); clinical survey pair counts inconsistent → recount by script at CGT freeze
  - [~] User download request sent (DepMap extras: GlobalSignatures, CRISPRConfounders, full mutations, CN segments; optional RRBS, PRISM)
  - [~] External fetch (clinical: CIViC, MOAlmanac, FDA CDx, OncoKB gene lists; combinatorial screens)
- [x] M2 Design
  - [x] DESIGN v1; round 1: 7 reviewers, 84 items (18 critical) in reviews/round1_*.md
  - [x] DESIGN v2 + reviews/round1_response.md (PI decisions D1–D3); configs/{design_params,primary_endpoints,tiers}.json
  - [x] Round 2: fresh PI + statistics + implementability, 23 items (6 critical) in reviews/round2_*.md
  - [x] DESIGN v3 + reviews/round2_response.md (D6), all 23 accepted; v3.1 amendment D7 (RRBS methylation: METH, SIL-METH)
  - [x] Round 3: single fresh reviewer, 6 items accepted → v3.2; tag freeze-design
- [~] M3 Infrastructure
  - [x] src/biomark (paths, p1 loaders, guard, P1 write-block); 31 tests pass; input sha256 11/11 OK
  - [x] External fetch: CIViC, MOAlmanac (GPL-2.0), FDA CDx (232 rows), OncoKB lists; combinatorial: Esmaeili 2024 harmonised table (covers Dede/Parrish/Thompson/CHyMErA/Ito), Dede, Thompson, Flister 2024 bioRxiv (analysis-only licence); Parrish processed tables blocked (not needed)
  - [x] CGT table: 64 tests; independent audit (2 critical, 5 major, 11 minor, all fixed: AUD-1..18) → 27 testable,
    18 P2-headline, 8 gate-disjoint; tag `freeze-cgt` (no dependency data joined); DESIGN amendments A1–A5 (append-only)
  - [x] Dictionary track committed; [~] single rebuild #4 with ALLELE/CLIN
  - [x] M3b parallel tracks (started during round 2, D4): dictionary/eligibility/controls · CGT table (no dependency data) · toy simulation (separate verification agent); tag `freeze-bounds`
- [~] M4 Core pipeline engine (stats/HC2, two-stage stability EN, 4 arms, MD-strict, fixed-threshold null) — pipeline author agent; unit tests on synthetic data; real run after freeze-bounds
  - [x] Engine committed (8cffc91): 19/19 unit tests; real pilot on 50 genes → K = 300, null B = 20 (≈ 96 min on 6 workers)
  - [!] Toy run 14/19 on 8cffc91: ridge-toy MD-strict FDP 0.228 (D12) → [!] catalogue = MD-global (A7, D15)
  - [x] Process checks (budget_check, leakscan + allowlist, check_frozen, §13.2 tests) → df fallback triggered, power-df and NMF fixes (D13)
  - [x] Independent toy verifier: NC bounds (9.2 min), null_genes.txt, diagnostics; final suite on 04aa982 19/22, MD-global 3/3
  - [x] Endpoint code (author): run_gate.py, replicate.py (+ pseudo-R fold-0 dry run), P2 McNemar, branch tests; XLIB, component readout, confounder audit, Q4
  - [x] Tag `freeze-bounds` (frozen_hashes `bounds` + `bounds_values`)
  - [x] configs/outcome_actions.csv: pre-committed action per gate/NC/catalogue/replication outcome, hashed before any real output (LESSONS 12)
  - [-] Real run attempt 1 (00:40): crashed at worker spawn because the Windows pickle pipe limit was hit by the full-D initargs payload; no output written. Fix: memory-mapped worker files (author)
  - [ ] Real run: 4 arms + MD lists, confounder audit, gate attempt 1, null (background), XLIB, dry run
  - [ ] Independent reimplementation check §13.3 (separate agent, after G_pc_attempt1)
- [ ] M5 Subtypes: per gene and per module, cross-lineage labels, co-location/cutting/culture/hypermutation audits
- [ ] M6 Clinical and combinatorial: CGT recovery, direction, actionability, population by lineage; combinatorial-screen lookup (Esmaeili harmonised table)
- [ ] M7 `freeze-catalogue` + r_requests.csv, then a single replication run on R (+ set A)
- [ ] M8 Limits classes, candidate cards, REPORT.md (writer), number check, red-team review, LESSONS prune
- [ ] Final REPORT.md
