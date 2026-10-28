# Research tree — Dependency subtypes and molecular biomarkers from the Part 1 essentiality predictor

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (latest first)

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
- [~] M2 Design
  - [x] DESIGN v1; round 1: 7 reviewers, 84 items (18 critical) in reviews/round1_*.md
  - [x] DESIGN v2 + reviews/round1_response.md (PI decisions D1–D3); configs/{design_params,primary_endpoints,tiers}.json
  - [x] Round 2: fresh PI + statistics + implementability, 23 items (6 critical) in reviews/round2_*.md
  - [~] DESIGN v3 + reviews/round2_response.md (D6)
- [~] M3 Infrastructure
  - [x] src/biomark (paths, p1 loaders, guard, P1 write-block); 31 tests pass; input sha256 11/11 OK
  - [x] External fetch: CIViC, MOAlmanac (GPL-2.0), FDA CDx (232 rows), OncoKB lists; combinatorial: Esmaeili 2024 harmonised table (covers Dede/Parrish/Thompson/CHyMErA/Ito), Dede, Thompson, Flister 2024 bioRxiv (analysis-only licence); Parrish processed tables blocked (not needed)
  - [~] M3b parallel tracks (started during round 2, D4): dictionary/eligibility/controls · CGT table (no dependency data) · toy simulation (separate verification agent)
- [ ] M3+ Implementation (defined by DESIGN.md)
- [ ] Final REPORT.md
