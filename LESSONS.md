# LESSONS — process rules for running this project

One rule per line, with evidence in brackets. (P1) = inherited from Part 1 (../project/LESSONS.md). Pruned at M2 (2026-10-28)
after a fresh-agent "stop doing / worked well" review (reviews/lessons_review_m2.md).

1. (P1) Check data access (APIs, bot walls, licences) before planning around a source; give every "if obtainable" input a deadline and a recorded fate. [Part 1: DepMap bot wall; here: R2-IM3 data deadline = first dictionary build, RRBS file present without checksum until D7 registered it]
2. (P1) Lock held-out data before any modelling and commit the lock. [Part 1: split locked day 1; here: replication set R never read, opened once, R2-IM2 `r_requests.csv`]
3. (P1) A process promise is a check a script enforces, or it is cut. [Part 1; here: every round-response row has a Check column, e.g. `gate.py params`, `budget_check.py`]
4. (P1) Before freezing a decision rule, dry-run it on non-blinded inputs so every verdict branch is reachable. [Part 1: unreachable "fails"; here: R3-4/R3-5 added reachability tests `test_endpoint_branches.py`, R2-ST7 `test_cross_lineage_labels.py`]
5. (P1) Set a control's expected bound by simulating it under the null before the run. [Part 1: permutation bands missed; here: R3-6 NC2-NC6 bounds = max over 5 toy replicates]
6. (P1) Parallel tracks get their own files and ID ranges, and each track gets design changes as a written diff list. [Part 1: registry corruption; here: D4 parallel start needed per-agent "Impl impact" in round2/3 responses]
7. Cap review at the point of diminishing returns (stop when a round yields <= 1 critical item), and take new scope only as a logged amendment with a named re-check of the touched configs. [84/23/6 items, critical 18/6/1; D7 caused R3-1]
8. Commit only finished work; before a suspend write the resume point and the list of unverified drafts into the tree, then resume stalled agents with SendMessage. [SUSPEND note; all 4 agents resumed with no rework]
9. Keep the author of a pipeline and its verifier as different agents, recorded by an `Agent:` commit trailer. [toy simulation agent; R2-IM6]
10. Split curation work over many clinical sources into a curated-mapping step and a mechanical-build step, one narrow query per source, with a handover file for any halted agent. [D5; cgt_handover.md]
11. Have a technical-brief agent map the upstream codebase before design. [M1: part1_brief found exact attribution, cached OOF/test predictions and 5 write hazards into Part 1]
12. Tell every web-facing agent never to send the user's email or other identifiers to external services. [M1: methods-survey agent passed the email to NCBI once; D5 PubMed check explicitly "no email or identifiers sent"]
13. Audit curated ground truth with an independent agent before freezing it, opening the cited records. [M3b: the CGT audit found a refuted row citing a wrong-direction record, read-through KMT2A "fusions" and an ER+ cut at 1.84 log2TPM, none caught by the builder's 19 passing tests.]
