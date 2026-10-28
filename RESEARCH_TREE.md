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

## Tree
- [x] M0 Setup — PROMPT.md, .gitignore, DATA.md, LESSONS.md, RESEARCH_TREE.md
- [~] M1 Literature and data availability → LITERATURE.md
  - [x] docs/part1_brief.md — OOF/test preds cached; M1/M2/T linear (exact attribution), M3 TreeSHAP; 5 write hazards into P1
  - [x] docs/lit/methods_prior_work.md, clinical_ground_truth.md, combinatorial_confounders.md; docs/data_availability.md
  - [x] LITERATURE.md synthesis (writer); clinical survey pair counts inconsistent → recount by script at CGT freeze
  - [~] User download request sent (DepMap extras: GlobalSignatures, CRISPRConfounders, full mutations, CN segments; optional RRBS, PRISM)
  - [~] External fetch (clinical: CIViC, MOAlmanac, FDA CDx, OncoKB gene lists; combinatorial screens)
- [~] M2 Design → DESIGN.md v1 drafted; review round 1 (7 reviewers: PI, prior work, stats, data, algorithms/compute, agent failure modes, clinical) running
- [~] M3 Infrastructure (design-independent): src/biomark loaders, replication guard, Part-1-untouched test
- [ ] M3+ Implementation (defined by DESIGN.md)
- [ ] Final REPORT.md
