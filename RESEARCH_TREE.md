# Research tree — Dependency subtypes and molecular biomarkers from the Part 1 essentiality predictor

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (latest first)

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
  - [x] CGT table: 64 tests (29 testable, 22 P2-scored), 524/526 record IDs resolve; [~] independent audit → freeze-cgt
  - [~] M3b parallel tracks (started during round 2, D4): dictionary/eligibility/controls · CGT table (no dependency data) · toy simulation (separate verification agent)
- [~] M4 Core pipeline engine (stats/HC2, two-stage stability EN, 4 arms, MD-strict, fixed-threshold null) — pipeline author agent; unit tests on synthetic data; real run after freeze-bounds
- [ ] Final REPORT.md
