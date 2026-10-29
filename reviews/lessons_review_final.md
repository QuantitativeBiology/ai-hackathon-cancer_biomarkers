# Lessons review at the final milestone (fresh process reviewer)

Sources read: LESSONS.md (12 rules), reviews/lessons_review_m2.md and lessons_review_m3b.md, RESEARCH_TREE.md (status
notes, D0-D22), results/deviations.csv (17 rows), `git log --oneline -80`. No data opened, nothing run. Not repeated from
M2/M3b: review-round caps, freeze order, enforcement-after-stage, accreting rulings. Those M3b items were acted on
(pre-committed outcome_actions, effective_rules.md, check_report_numbers.py), and the M3b risk (post-hoc rulings on first
real results) did not materialise: no ruling was written after a real result except bug fixes and corrections (D19-D21).

## 1. Stop doing

1. **Running multi-hour jobs with no checkpoint and no restart rule.** The fixed-threshold null (started 03:38) had no
   per-draw checkpoint; the overnight VS Code closure killed it and four other jobs with no output (INTERRUPTION note,
   08:40). It was re-run from zero at 08:49 (about 3.2 h at 3 workers). The SUSPEND note then states "it has no
   checkpoint" again, and the P1 pool build ran 2 h 06 of CPU with no output before a fast path replaced it. Cost: roughly
   5 h of null compute lost, the null finished late, DEV-BUDGET-M3B-M4 and DEV-CUT-XLIB-DISC (a cut forced by lateness),
   and the 24 h budget was passed (DEV-BUDGET-OVERRUN, user decision at 16:19). Fix: a job over ~20 min writes per-unit
   output (per gene block or per draw) and `--resume`; this is how the arms run (parquet blocks) and the null did not.
2. **Not budgeting memory against the real machine before writing the engine.** D17 found, after the engine existed, a
   26 M-row arms table of about 25 GB on a 15.7 GB machine; real run attempt 1 crashed at worker spawn (pickle pipe
   limit); the first null restart got 1 worker at 2.1 GB free; the pool build thrashed swap and was rewritten twice
   (4253173 low-memory EB shrinkage, 74d0536 fast path). Cost: dfcf424 plus the two rewrites, walls over every timebox,
   and about 3.5 GB free at restart cut the null to 3 workers. The time budget was done in hours; memory was never done
   in GB per worker times workers against free RAM. Fix: a peak-RSS line per stage in the pilot, checked by
   budget_check against `psutil` free memory at launch.
3. **Logging a cut before the failure is understood and before checking the design.** DEV-CUT-DRYRUN is stamped
   15:11:50, the same minute the rebuilt pool job started, citing the swap-thrashed pool. It was withdrawn 7 minutes later
   (DEV-CUT-DRYRUN-WITHDRAWN) because DESIGN 1.3 requires the dry run before the tag, it was not in `cut_order`, and
   `make_pool` had been rewritten after the last toy suite. Only the red-team review stopped the irreversible R read from
   being the first real execution of that path. DEV-CUT-XLIB-DISC likewise combined a lateness cut with a crash that was
   masked by it. Cost: one round of rework and a near miss on the one irreversible step. Fix: a cut is allowed only if it is
   listed in `cut_order`, a failed job is first rerun or diagnosed, and the cut row cites the DESIGN line it relaxes.
4. **Writing a universal claim before reading every arm.** D20 (04:40) said all 7 non-eligible controls were "missed by
   every arm" and that ESR1 "could never be recovered". D21 (15:20) found S1-pred recovers MYCN, SMARCA2, ARID1B and
   EN-obs recovers ESR1 and SMARCA2. The false sentence stood about 10 h and fed the REPORT draft, DEV-CUT-DRYRUN's
   context and df7b7c6. Same family: the null summary shipped with MD-global `hits_obs = 0` (a7c6c8a), caught by reading, then
   recomputed (DEV-NULL-RECOUNT). Cost: D21, DEV rows, red-team fixes. Fix: "every / none / never" in the tree or report
   is generated from a per-arm table (script output), not typed; `check_report_numbers` already exists for REPORT, so run
   it on RESEARCH_TREE decisions too.
5. **Back-filling bookkeeping in bulk.** DEV-DICT-1 and DEV-CGT-1 were written at 08:47 on 10-29 for D9/D10 events from
   10-28; registry rows were back-filled (c5e0e54, 073c9af `backfill_registry.py`); timing for M3b builders is permanently
   unknown (D13-viii). Cost: small but real: the deviation file is not a contemporaneous record for 5 of 17 rows. Fix: the
   decision that changes a frozen value and its deviation row land in one commit (a pre-commit check can enforce it).
   Note: the "scripts that ignore `--help` and run for real" candidate is not evidenced in any file I was allowed to read
   (no deviation, decision or commit message mentions it); I did not include it. If the author has an incident, it stays
   in the review file until it recurs (LESSONS preamble).

## 2. Worked well

1. **Pre-committing actions for outcomes (`outcome_actions.csv`, 56752ee), written before any real output.** The gate
   failed twice (21/39, D19, D20), NC1 failed for EN-obs (1cb695e), and P2 was underpowered (6 < 15); each time the
   action was a table lookup, not a new ruling. The D12-D16 pattern (a ruling every ~40 min) did not recur on real data.
2. **Reporting failures instead of tuning.** Gate failure stated as "21/39 at attempts 1 and 2"; EN-pred no better than
   EN-obs (-0.07 (-0.12, 0.02)) and univariate better (-0.10) reported as the primary endpoint outcome (20cd7f7);
   DEV-CAL-CONSERVATIVE left as a limit; toy 19/22 kept.
3. **The single irreversible read was protected by mechanism and by one fresh reviewer.** The R set was read once, on a
   clean tree after `freeze-catalogue` (d7d7835), the dry run was restored on red-team advice (cba4ba1, ae6aa84), the
   dry-run numbers were registry-marked uncitable, and the real read is logged once (20cd7f7). Fast paths were accepted
   only with `check_exact` equality to the slow path (74d0536, 4253173 "bit-identical").
4. **Independent verification of outputs, not just code.** Reimplementation check 4c902fb: 59,955/59,955 S1 p-values,
   gate 21/39 exact, one tolerance failure reported as a failure. The red-team review found D20's error, the RRBS
   456/453 mismatch and the pilot qualifier on A7 (df7b7c6). `check_report_numbers --strict` and the registry make a stale
   number a build error.
5. **Suspend and interruption notes that a cold session can execute.** Each gives the last commit hash, the running PID
   and log, the unverified drafts, and numbered next steps (SUSPEND #2 includes the exact fast-path command and a 20-min
   restart rule); agents were resumed by SendMessage and the tree was clean after each. The interruption cost five jobs
   but no committed work and no lost decisions.

## 3. LESSONS.md final edit (12 rules -> 13)

Merge, reword or add exactly as follows; rules not listed (1, 4, 5, 6, 7, 10, 11) stay unchanged except the bracketed
evidence where shown. The standing web-agent clause and the "one incident long" preamble line stay.

- **Rule 2 (reword evidence):** `2. (P1) Lock held-out data before any modelling and commit the lock. [Part 1: split locked
  day 1; here: R sealed by the tag guard and read once on a clean tree after freeze-catalogue, dry run restored before the
  tag (DEV-CUT-DRYRUN-WITHDRAWN, 20cd7f7)]`
- **Rule 3 (unchanged text, add evidence):** `... [here: budget_check retrofitted after the engine, 9 non-conforming
  timing writers, M3b wall time unknown (D13-viii); memory had no check at all until D17]`
- **Rule 8 (merge checkpointing in):** `8. Commit only finished work. Every job over 20 min writes per-unit output and has
  a resume flag; before a suspend write the resume point (last hash, PID, log, exact restart command, unverified drafts)
  into the tree, then resume stalled agents with SendMessage. [null had no checkpoint: lost 03:38 run, re-run 3.2 h
  (INTERRUPTION note); SUSPEND #2 note with a 20-min restart rule; 4 agents resumed with no rework]`
- **Rule 9 (add a clause):** `9. Author and verifier are different agents (`Agent:` trailer); verifiers of curated truth
  open the cited records; a brief barring a data scope names the tests that read it; a claim of "every / none / never" is
  generated from a per-arm table, not typed. [toy verifier found the double dip; CGT audit AUD-1..18; guard agent read dev
  effects via test_guard.py; D20 "missed by every arm" false, corrected in D21]`
- **Rule 12 (merge the cut rule):** `12. Pre-commit the action for every control, gate and budget outcome before the first
  real result, and quote rules only from the current amendments. A budget cut is allowed only if it is in `cut_order`,
  the failure behind it is first diagnosed or rerun, and the row names the DESIGN line it relaxes. [outcome_actions.csv
  written before any real output; gate failed twice with no new ruling; DEV-CUT-DRYRUN logged 15:11:50, withdrawn 7 min
  later]`
- **New rule 13 (memory and time budget in the unit that binds):** `13. Before writing a stage, measure its peak memory in
  the pilot and budget workers x GB against free RAM on this machine, in the same check as the time budget; a stage that
  does not fit is blocked or split before launch, not capped afterwards. [26 M rows ~25 GB on 15.7 GB (D17); attempt 1
  crashed at spawn; pool thrashed swap, rewritten twice (4253173, 74d0536); DEV-BUDGET-M3B-M4, DEV-BUDGET-OVERRUN]`
- **Add to rule 9 or the preamble, one line:** `A decision that changes a frozen value and its deviation row land in one
  commit.` (evidence: DEV-DICT-1, DEV-CGT-1 back-filled at 08:47). If a slot is needed, this is the first line to drop.
- **Preamble:** change "Pruned at M2 and at M3b" to "Pruned at M2, M3b and the final milestone (reviews/lessons_review_final.md)".

Result: rules 1-7 and 10-11 as now, 8, 9 and 12 widened, 13 new = 13 rules, no free slot. Rule 11 (technical brief) and
rule 1 are the weakest evidenced for future work; if a new incident needs a slot, drop 11 first (one incident, M1, already
paid for).
