# Lessons review at M2 (fresh process reviewer)

Sources read: LESSONS.md, RESEARCH_TREE.md, full git log (14 commits), reviews/round{1,2,3}_response.md, docs/cgt_handover.md.

## What should we stop doing?

- Stop piling up review rounds with ever-larger frozen specs before any result exists: 113 review items across 3 rounds and 14 commits, and no analysis result yet. [round1 84 items (18 critical), round2 23 (6 critical), round3 6; RESEARCH_TREE "No analysis result exists yet"; 2.6 h of 24 used]
- Stop letting scope amendments land after reviewers have signed off: they re-open the frozen design. [D7 RRBS added mid-round 3 produced R3-1, a critical item, and forced configs edits to stages.G, fwl_unpenalised, not_obtained]
- Stop starting implementation tracks against a moving design without a diff protocol: every round creates per-agent rework. [D4 started M3b during round 2; round2/3 "Impl impact" columns assign DICT/CGT/TOY edits for 10+ items, e.g. R2-FP-1, R2-ST4, R3-1, R3-2, R3-6]
- Stop leaving unfinished agent output uncommitted on disk across a suspend, and stop relying on the next session to recognise it as a draft. [RESEARCH_TREE SUSPEND: four tracks left "unverified drafts"; git status still shows ~30 untracked files and a modified p1.py]
- Stop doing resource and budget arithmetic in prose; reviewers found it late. [round2 conflict 4: "15 h + 3 h reserve" exceeded the ~17 h left; round3 had to re-set hours_remaining_at_v3 16.0; clinical survey pair counts inconsistent, LITERATURE note]

## What worked well?

- Item tables with ID, severity, decision, DESIGN section, a machine-checkable Check, owner and implementation impact; every item got a recorded decision, so none were silently dropped. [round1/2/3_response.md; 113 items, 0 open, D8]
- Fresh reviewers and diminishing returns as a stopping signal: later rounds found fewer and milder problems, and round 3 stopped at one reviewer. [84 -> 23 -> 6 items; critical 18 -> 6 -> 1]
- Freezing by hash and tags with an ancestry rule (freeze-design, then freeze-bounds), plus a PI-ruled declined/accepted list. [4a079e8 "frozen hashes (design)"; R2-IM1; D8 "parameters cannot change after results are seen"]
- Splitting hard work into a curated step and a mechanical step, with a handover file that preserves what a halted agent learned, and resuming by SendMessage. [D5; docs/cgt_handover.md lists located record IDs; SUSPEND "all 4 agents resumed via SendMessage"]
- Keeping the pipeline author and the verifier separate, and mapping upstream before designing. [toy simulation by a separate verification agent, R2-IM6 `Agent:` trailer rule; part1_brief found 5 write hazards]

## Proposed LESSONS.md

(12 rules; (P1) marks inherited rules.)

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

Dropped: old P1 rule 7 (reviewers return reports as text) because all reviews here landed as files in reviews/ with no refusal; P1 rules 1-6 kept (rule 6 extended). Old rules 8-11 kept as 12, 8 (merged with 10 of old), 11, 10 respectively.
