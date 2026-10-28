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
- [~] M1 Literature and data availability → LITERATURE.md (5 parallel agents: part1 brief, methods, clinical, combinatorial+confounders, data access)
- [ ] M2 Design → DESIGN.md, review rounds (≤3) in `reviews/`
- [ ] M3+ Implementation (defined by DESIGN.md)
- [ ] Final REPORT.md
