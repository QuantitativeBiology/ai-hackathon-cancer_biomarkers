# Lessons review at M3b (fresh process reviewer)

Sources read: LESSONS.md (15 rules), RESEARCH_TREE.md (D0-D16), `git log -60`, DESIGN.md "Amendments after freeze" (A1-A7),
reviews/lessons_review_m2.md, docs/process_checks.md, docs/toy_nc_bounds.md, docs/impl_notes_pipeline.md (skim).
No data or results were opened. Not repeated from the M2 review: review-round caps, suspend/draft handling, parallel-track
diffs (the M2 "stop" items were in fact stopped: git status is clean, no review round after M2, amendments carry reasons).

## 1. Stop doing

1. **Freezing a decision rule before the engine has run on the toy.** `freeze-design` (v3.2) fixed MD-strict as the
   catalogue (A6-i repeated it at D11). Only afterwards did the toy show MD-strict FDP 0.16-0.26 on 5/5 ridge toys and
   P(any hit) 0.30 on pure nulls (D12, D14), forcing A7/D15 (catalogue = MD-global) plus a correction clause A7(v) that
   fixes A7's own numbers (0.85-0.93 -> 0.815-0.926). Rule 4 says dry-run before freezing, but the dry run was scheduled
   after the freeze. Cost: five decisions (D12-D16) and four DESIGN edits in about four hours; the rescue was only clean
   because no real association result existed yet. Next time the toy engine run is a precondition of the design tag.
2. **Tagging a freeze whose file list is not complete.** `check_frozen` found six DESIGN 13.4 files created after
   `freeze-design` (process_checks 3; D13-vii), so that tag did not freeze what DESIGN says it freezes. The patch was a
   second hash section (`bounds`). Cost: the "frozen" claim in DESIGN is partly false and relies on a reader knowing
   D13-vii. Run `check_frozen.py --strict-coverage` before creating a tag, and fail the tag if it is not clean.
3. **Writing the enforcement for a process promise after the work it governs.** `budget_check` was built after the
   engine; process_checks 1.2 lists 9+ writers that did not meet the timing schema, and D13-viii records that all M3b
   builders wrote no timing, so their wall time is "unknown, logged as deviation". The "7.5 h of 24" figure in the tree
   is therefore not machine-verified for M3b. Rule 3 says a promise is a script check or is cut; here the promise
   (budgeted tiers with cut order) existed for hours with no check. Cost: unrecoverable time data, a schema retrofit
   (f2abd93), and a still-open hole (`write_timing` defaults missing `workers` to 1; the vectorised-null sub-record is not
   counted in concurrency). Write the checker and the writer helper in the same commit that adds a budgeted stage.
4. **Setting bounds that cannot fail, then calling the control a bound.** After rulings R1/R2, NC2b is 1.2 MD-strict hits
   per low-tier gene against a provisional 0.1 (12x looser); NC3, NC5, NC6 keep provisional values because the toy was
   degenerate (D14-R2, toy_nc_bounds addendum); 3 of 22 toy checks fail and are "reported, not changed" (D16). Also
   NC2b is applied to MD-strict, which after A7 is no longer the catalogue (A7: gate and NC2b "stay on MD-strict"). Rule 5
   is satisfied in form (simulated, max of 5) and not in function: a control that passes whatever happens guards nothing.
   Cost: after M4 a failed negative control will have no pre-agreed meaning, which invites post-hoc reading.
5. **Accreting rulings as nested prose instead of a current-state table.** D13 has 8 sub-items, D14 has R1/R2 plus three
   interim findings, A7 has (i)-(v) with a self-correction, and DESIGN body text still describes MD-strict as the
   catalogue while the amendment overrides it. `impl_notes_pipeline.md` lists 13+ "ambiguities resolved" that also change
   behaviour. Cost: the report writer (M8) and the independent re-implementer (13.3) must reconcile body plus amendments
   plus decisions; that is where a number or rule will be quoted stale. Add one generated "effective parameters and rules"
   page (catalogue = MD-global, gate = MD-strict, df rule = fallback, NC values) checked against `frozen_hashes`.
6. **Letting lessons record single incidents as rules mid-milestone.** Rules 13-15 are each one event with an "M3b"/"M4"
   tag, written in the same commit as the fix (3b6a237, 94cbedf). They are useful, but 14 and 15 are narrow
   (a test file name; one function signature) and sit beside general rules; the file reached its 15-rule ceiling in one
   day. Cost is small today, but the file stops being read when it grows by incident. See section 3.

## 2. Worked well

1. **A toy simulation with a pre-written answer key, built by a different agent from the pipeline author.** It exposed
   the double-dip FDP (0.16-0.26 vs 0 with MD-global at equal recall), found before any real result (D12, D15, A7 header
   "none made after any association result existed"). This is the single most valuable control in the project.
2. **Process-guard agent on the engine.** Within one pass it found a size problem (CR2 + nu size 0.069-0.073 at
   n = 255), power df != test df, non-invariant NMF scores, and 54 leakscan hits (process_checks 5, D13). All were fixed
   and re-verified at 04aa982 (83e575b: suite 236 passed). Fixes were applied as the pre-declared 5.3 fallback, not as a
   choice, which keeps them non-post-hoc.
3. **Writing definitions before running them, with the hash recorded.** `toy_nc_bounds.md` states analogues, units,
   denominators and fallbacks before the script ran (sha256 is in the output JSON), and rulings R1/R2 are dated with
   "no ridge value existed yet". D16 then applied R1 mechanically instead of revisiting it. This is the form rule 5 should
   require.
4. **Reporting failures instead of tuning to them.** 19 pass / 3 fail on the final suite is stated in the tree, with the
   cause of each failure (different list for the bound; different denominators; 4 marginal hits) and no bound edited.
   The same applies to DEV-CAL-CONSERVATIVE (4 conservative cells kept).
5. **Independent audit of curated ground truth before its freeze.** The CGT audit (reviews/cgt_audit.md) found 2 critical
   and 5 major issues that 19 builder tests missed (a refuted row citing a wrong-direction record; ER+ cut). Result:
   27 testable / 18 headline tests, `freeze-cgt` tag, no dependency data joined.
6. **Mechanical freeze/ancestry and leakage checks with allowlists.** `check_frozen` (hash plus "freeze-bounds is an
   ancestor of the first results/ commit"), `leakscan` (symbols and constants, exemptions only via a reasoned allowlist
   file and printed as `EXEMPT`), and scripts that refuse to run before the tag (impl_notes: `run_arms.py`, `run_null.py`).
   They turn rules 2 and 3 into exit codes.

## 3. LESSONS.md pruning proposal (15 -> 13 rules)

- **Merge 4 + 5** into one: "Before freezing a rule or a bound, run it on the toy/non-blinded inputs: every verdict branch
  reachable, every control bound set from its null simulation and able to fail." Evidence: unreachable fails (P1);
  NC2b 1.2 and provisional NC3/5/6 cannot fail (D14, D16).
- **Merge 9 + 13 + 14** into: "Author and verifier are different agents (`Agent:` trailer); verifiers of curated truth open
  the cited records; a brief barring a data scope names the tests that read it." Evidence: toy verifier, CGT audit,
  test_guard.py reading dev effects.
- **Reword 15**: "Test the deployed rule through the real entry point, not the function default." Keep; evidence: calibration
  test ran cr2_core without carrier counts.
- **Reword 3**: add "and the checker lands in the same commit as the stage it checks." Evidence: budget_check retrofit,
  9 non-conforming timing writers, M3b time unknown.
- **Drop 10** (curated-mapping vs mechanical-build split): task-specific to the finished CGT build; the general part is
  already in rule 8 (handover file) and the new freeze-order rule.
- **Drop 12 as a rule**; move it into a standing brief template ("web agents: send no email or identifiers") that every
  web-facing dispatch reuses; as a lesson it is a one-line prompt clause, not a process rule.
- **New rule A (order of freezes):** "No design tag until the engine has run on the toy and `check_frozen --strict-coverage`
  is clean; later changes are amendments that list their effect on every frozen bound." Evidence: A6/A7, D13-vii.
- **New rule B (one effective-rules page):** "Every amendment updates a generated current-rules table; the report cites only
  that table." Evidence: MD-strict/MD-global split across DESIGN body, A7, D15, D16.
- Keep 1, 2, 6, 7, 8, 11 unchanged. Result: 1, 2, 3', 4+5, 6, 7, 8, 9+13+14, 11, 15', A, B = 12 rules; one slot remains
  free for the next incident. Add a standing sentence: "a rule that is one incident long stays in the review file until
  it recurs."

## 4. One process risk for the next milestones

**Post-hoc interpretation freedom at the first real result.** The controls and gate are weaker than they look: NC2b is
1.2 and lives on a list that is no longer the catalogue; NC3/NC5/NC6 are provisional; 3 toy checks fail; the gate stays
on MD-strict while all downstream analyses use MD-global; and the P1 replication in R is the one irreversible step
(255 lines, opened once). If the first real run misses any control or the gate, the project has no pre-agreed action, and
the D12-D16 pattern (a ruling every ~40 minutes) will recur on real data, now with outcomes visible. The same applies
inside the R replication and the report, where a quoted number could come from a pre-amendment rule.

**Check that would catch it:** before the first real `run_arms`, add and tag `freeze-actions`: a table (hashed in
`frozen_hashes`) mapping every gate/NC/calibration outcome (pass, fail, NaN, degenerate) to one pre-committed action
(proceed, report-only, stop and escalate), including "NC2b on MD-global reported next to MD-strict". `run_gate.py` and
`replicate.py` read the table and refuse to start R unless every real outcome is mapped and `check_frozen` passes. Add
a report-time script that diffs every number in REPORT.md against the effective-rules page, so a stale MD-strict or
pre-fallback value fails the build.
