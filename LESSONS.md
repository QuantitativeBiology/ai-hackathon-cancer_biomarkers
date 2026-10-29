# LESSONS — process rules for running this project

One rule per line, with evidence in brackets. (P1) = inherited from Part 1 (../project/LESSONS.md). Pruned at M2 and at
M3b (2026-10-29) after fresh-agent "stop doing / worked well" reviews (reviews/lessons_review_m2.md, lessons_review_m3b.md).
A rule that is one incident long stays in the review file until it recurs.

1. (P1) Check data access (APIs, bot walls, licences) before planning around a source; give every "if obtainable" input a deadline and a recorded fate. [Part 1: DepMap bot wall; here: RRBS registered by D7 before the first dictionary build]
2. (P1) Lock held-out data before any modelling and commit the lock. [Part 1: split locked day 1; here: R sealed by the tag guard; no access log exists before freeze-catalogue]
3. (P1) A process promise is a check a script enforces, and the checker lands in the same commit as the stage it checks. [here: budget_check retrofitted after the engine, 9 non-conforming timing writers, M3b wall time unknown (D13-viii)]
4. (P1) Before freezing a rule or a bound, run it on toy or non-blinded inputs: every verdict branch reachable, and every control bound set from its null simulation and able to fail. [Part 1: unreachable fails; here: test_endpoint_branches; NC2b = 1.2 and provisional NC3/5/6 (D14, D16)]
5. No design tag until the real engine has run on the toy and `check_frozen --strict-coverage` is clean; later changes are amendments that list their effect on every frozen bound. [MD-strict frozen as the catalogue, then FDP 0.16–0.26 on ridge toys → A7; six §13.4 files missing at freeze-design (D13-vii)]
6. (P1) Parallel tracks get their own files and ID ranges, and each track gets design changes as a written diff list. [Part 1: registry corruption; here: D4]
7. Cap review at the point of diminishing returns (stop when a round yields <= 1 critical item), and take new scope only as a logged amendment with a named re-check of the touched configs. [84/23/6 items, critical 18/6/1; D7 caused R3-1]
8. Commit only finished work; before a suspend write the resume point and the list of unverified drafts into the tree, then resume stalled agents with SendMessage. [SUSPEND note; all 4 agents resumed with no rework]
9. Author and verifier are different agents (`Agent:` commit trailer); verifiers of curated truth open the cited records; a brief barring a data scope names the tests that read it. [toy verifier found the double dip; CGT audit AUD-1..18; guard agent read dev effects via test_guard.py]
10. Test the deployed rule through the real entry point, not a function's defaults. [calibration test called cr2_core without carrier counts and kept measuring the old ν rule]
11. Have a technical-brief agent map the upstream codebase before design. [M1: part1_brief found exact attribution, cached OOF/test predictions and 5 write hazards into Part 1]
12. Pre-commit the action for every control and gate outcome before the first real result, and quote rules only from the current amendments. [M3b review risk; configs/outcome_actions.csv written before any real output was opened]

Standing brief clause for web-facing agents (formerly a rule): never send the user's email or other identifiers to external services. [M1: NCBI incident]
