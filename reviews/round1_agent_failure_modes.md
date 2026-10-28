# Review round 1 — AI-agent failure modes (reviewer: R1-AF)

Scope: places in DESIGN v1 where LLM agents could, without anyone noticing, move a threshold after seeing results, leak
replication data, invent ground-truth rows or citations, over-claim, lose decisions, write into Part 1, or put numbers in
the report that no results file backs. Evidence comes from Part 1 LESSONS (#5: accepted items were "done" without being
done; #6: prose safeguards lapsed and only gate scripts held; #7: unreachable verdict branches; #10: control bands never
simulated) and from inspecting this repo and ../project read-only.

What the design already does well: R is opened once behind a tag, the access log, sign pre-registration for the CGT,
a positive-control gate, toy data with planted answers, and seeds and config hashes. The problem is that most of these
are written as prose. Below, each one becomes a check that a script runs. All the checks fit in one
`scripts/gate.py` with subcommands, using Part 1's `gate.py` / `check_report_numbers.py` as the template (copy them; do
not import them). Estimated build cost for all items: about 1.5 to 2 h inside M3, plus about 0.5 h inside M8.

---

### R1-AF1 — critical — The replication guard is advisory; R can leak through direct reads, and the one-shot opening is not protected
**Problem.** `P1/cache/crispr.parquet` is a plain parquet file, and part1_brief §2 says Part 2 must read it directly. Any
script that calls `pd.read_parquet` on it gets R targets, and `src/biomark/guard.py` never sees the call. The same is true
of `P1/results/test/*` (per-gene test metrics for m5_avg, m1, b0–b2) and `P1/cache/test/`. Tags can be moved with
`git tag -f`. Nothing stops an agent from editing the catalogue or the replication script after R is opened and running
it again. Nothing tests the replication code before the single run, so a crash or a bug found after opening invites a
second, adapted opening.

**Change (add to §1).**
1. *Single entry point.* Only `src/biomark/guard.py` may contain the strings `crispr.parquet`, `CRISPRGeneEffect`,
   `results/test`, `results\test`, `cache/test` or `cache\test`, or call `essml.guard` / `io.load_matrix`.
   `guard.load_targets(ids)` returns rows only for dev IDs (from `P1/configs/splits.csv`, role == dev). If any requested
   ID is not dev, it raises unless all of the following hold: (a) tag `freeze-catalogue` exists; (b) its commit equals
   the hash recorded in `configs/freeze_catalogue.json` (committed in the tagged commit); (c) the working tree is clean;
   (d) the caller is `scripts/replicate.py`; (e) the sha256 of every frozen catalogue file matches that JSON.
2. *Re-runs are reproductions only.* A second R access is allowed only with the identical config hash and code commit
   (deterministic reproduction). Any other hash is refused. Every access appends a row to
   `results/replication_access_log.csv`.
3. *Dry run before freezing.* `scripts/replicate.py` runs end to end with D fold 0 standing in for R
   (`--pseudo-R fold0`, catalogue rebuilt on folds 1–4 or simply reused) before `freeze-catalogue`. Its output goes to
   `results/replication_dryrun/`, and every verdict branch is reached or shown to be reachable (P1 LESSONS #7).
4. *Immutable after opening.* After the first real access,
   `git diff <freeze-catalogue>..HEAD -- configs/catalogue/ scripts/replicate.py src/biomark/replication*.py` must be empty.

**Check.** `gate.py leakscan`: a grep over `src/ scripts/ tests/ notebooks/` for the forbidden strings outside guard.py,
plus a grep for `read_parquet(`/`read_csv(` on any path under `../project/cache` or `../project/results` that is not on
a whitelist (`configs/p1_read_whitelist.txt`). `tests/test_guard.py` must show that a test ID raises before the tag,
that a moved tag raises, and that a changed hash raises. `gate.py G_rep` must confirm the access log has at least one row,
every row has the same config hash, the frozen-path diff is empty, and the dry-run output exists in a commit that is an
ancestor of the tag.

### R1-AF2 — critical — Thresholds and the primary comparisons are scattered in prose and can be changed after results are seen
**Problem.** About 25 decision numbers are only in prose: −0.5/−1.0, the 20th percentile, α = 0.9, 1-SE, 50/20
subsamples, 60% stability, 25% line contribution, BH 0.05, |β| drop > 50%, ≥ 5 carriers, r > 0.8 within 10 Mb, the 60%
gate, power ≥ 0.5, k ∈ {50, 100, 200, 500}, fidelity < 0.1, k = 15 for Leiden, ≥ 3/≥ 3 lineage carriers, 200
permutations, the top 5% TMB, the top-5 agreement. Q2 lists four comparison families with four list sizes each and names
no primary one, so an agent can report whichever comparison favours MD. The "compute guard" fallback (50 → 20
subsamples) triggers on a measured runtime, which is fine only if it is logged.

**Change (add a §0.1 "Frozen parameters and endpoints").**
- All decision numbers live in `configs/design_params.json`, keyed by design section. Code reads them only from there.
  The file is hashed into `configs/frozen_hashes.json` at tag `freeze-design`.
- Declare **one primary endpoint per question**, with its verdict rule and every branch. Suggested for Q2: the R
  replication rate (same sign, one-sided p < 0.05) of MD vs S1 at matched k = 100 among candidates with power ≥ 0.5,
  with a paired bootstrap CI over genes. Verdict: "MD better" if the CI is above 0, "S1 better" if it is below, else
  "no difference". Everything else is labelled secondary in the report. The same is done for Q4a (recall on concordant,
  testable S/W rows at the CGT tier "CIViC A/B ∪ FDA-CDx").
- *Deviation rule.* A post-freeze change to `design_params.json`, the dictionary, the controls or the CGT needs a
  decision-log entry `D<n>` in RESEARCH_TREE.md. The entry says whether any output depending on the change had already
  been seen. If it had, the report gives the headline under both the old and the new value. Fallbacks such as 50 → 20
  subsamples are written automatically to `results/deviations.csv` by the code that triggers them.

**Check.** `gate.py frozen`: the sha256 of each frozen file must match `frozen_hashes.json`, or a `D<n>` must exist whose
text names that file and whose ID appears in the commit message that changed it (`git log --format=%B -- <file>`). The
report-number check (R1-AF5) must find each primary endpoint tagged in REPORT.md's first results section. A unit test
must confirm that grep finds no literal `-0.5`, `0.05`, `0.6`, `0.25` or `0.8` in `src/biomark/*.py` outside
`params.py` (whitelist for unrelated uses).

### R1-AF3 — critical — Clinical ground truth and combinatorial-screen labels can be fabricated or written from memory
**Problem.** The 45-row candidate CGT table in docs/lit/clinical_ground_truth.md was written by an agent. Its
`[verified]` tags are self-reported, and the survey itself says several claims are `[from-memory]`. DESIGN §7.3 builds
the CGT from "sources plus a curated mapping file with a citation per row", and that curated file is where an invented
PMID, NCT or sign could enter unnoticed. Pre-registering the sign "before any dependency data are joined" has no
mechanism behind it. §7.4 names five combinatorial screens (Parrish GSE178179, Dede 2020, in4mer, Klingbeil 2024,
Thompson 2021) that appear neither in data_availability.md nor in DATA.md, so an agent may fill in "replicated" labels
from memory.

**Change (add to §7.3 and §7.4).**
- Every row of `configs/cgt_mapping.csv` carries at least one **resolvable source key**: a CIViC EID/AID, a MOAlmanac
  statement id, an FDA PMA/BLA number present in the frozen CDx scrape, or an NCT id present in a frozen CTG JSON. All of
  these sit under `external/` with their sha256 in `configs/external_sha256.txt`. A PMID is supplementary and counts only
  if it resolves in a frozen PubMed esummary download, fetched **without** any email/tool parameter (LESSONS #8). Rows
  with no resolvable key are dropped and listed in `results/clinical/cgt_dropped.csv`. They are not kept "from
  knowledge".
- The expected sign and the KO-concordance class are columns of `cgt_mapping.csv`, frozen by tag `freeze-cgt`. The CGT
  build script must not import guard.py or read any effect or prediction file. The tag must be an ancestor of the first
  commit that adds anything under `results/clinical/`.
- §7.4: a pair label (replicated / tested-not-replicated / untested) may only come from a downloaded supplementary file
  with a row in DATA.md. If no file has arrived by the data deadline, §7.4 has status "not done: data unavailable". Give
  the user a download request now (LESSONS #1).

**Check.** `gate.py cgt`: every source key in cgt_mapping.csv must resolve in the frozen files (exact string match), and
the file hashes must match external_sha256.txt. A static import scan of `scripts/build_cgt.py` must find no reference to
the guard, crispr or preds. `git merge-base --is-ancestor freeze-cgt <first commit touching results/clinical>` must
succeed. §7.4 outputs must reference only files listed in DATA.md.

### R1-AF4 — major — The positive-control gate is ambiguous in feature terms, and "stop and diagnose" invites tuning until it passes
**Problem.** The controls are named biologically ("ENO1 loss", "VPS4B loss", "ESR1 expression", "MDM2 in TP53-WT", "EGFR
mut/amp"), not as 𝒟 feature IDs. After locus collapsing (§5.4), the ENO1 deletion becomes a 1p36 locus feature, and VPS4B
(18q) and 9p21 MTAP have the same problem. The denominator ("those with ≥ 5 carriers in D") is unknown until the data are
seen. "Stop and diagnose" has no bound, so the natural agent behaviour is to tweak α, stability or BH until 60% is
reached. Part 1 also found named entities missing at run time (LESSONS #7: POU5F1, the EC lines).

**Change (add to §7.1).** `configs/positive_controls.csv` lists (control_id, dependency Entrez, the 𝒟 feature ID(s)
accepted as a match, including the locus feature that contains the named gene, expected sign, and a citation). A script
resolves it against `dictionary_v1.json`. Carrier counts in D (no association computed) are written to
`results/controls/pc_eligibility.csv` **before** `freeze-design`, which fixes the denominator. The first gate evaluation
is written to `results/gates/G_pc_attempt1.json` and is never overwritten. If it fails, any fix is a D-entry under
R1-AF2, the report states "gate failed at attempt 1 (x/y); passed at attempt n after change Z", and the control list
cannot change after freeze.

**Check.** `tests/test_controls.py`: every control's feature IDs exist in the dictionary (or the control is listed as
ineligible with a reason). `gate.py G_pc` reads attempt 1 and later attempts, and fails if attempt 1 is missing or was
modified (its sha256 is recorded in the gate log at creation).

### R1-AF5 — major — No report-number or claim check, and the .gitignore would silently drop results files
**Problem.** DESIGN §9 says "every headline number is produced by a script into results/", but nothing ties REPORT.md
numbers to files. Part 1's `check_report_numbers.py` did exactly this and worked (LESSONS #6). There is also a concrete
gap: this repo's `.gitignore` ignores `*.gz` everywhere, so any `results/*.csv.gz` (Part 1 wrote `per_gene.csv.gz`) will
never be committed, and the report would then cite numbers that git cannot reproduce. Claim words ("validated",
"replicated", "cross-lineage", "valid biomarker", "recovers") are defined by rules in §§5.2, 6 and 8, but nothing stops
the writer agent from using them loosely.

**Change (add to §9).**
- Port Part 1's tag syntax (`<!--n:results/...json#path-->0.31`) into `scripts/check_report_numbers.py`. Run `--strict`
  at the final gate. Every tagged source file must be tracked by git. Add `!results/**` and `!configs/**` exceptions to
  .gitignore, or ban .gz under results/.
- Add a claims table, `results/claims.csv` (gene/pair, claim class, evidence file, rule), generated by script from the
  §8 failure classes, the §5.2 lineage labels and the R1 replication output. In REPORT.md, a sentence that applies a
  controlled claim word to a named gene or pair must carry `<!--c:claim_id-->`, and the claim class must match.

**Check.** `check_report_numbers.py --strict` exits 0. `git ls-files --error-unmatch` succeeds for each tagged source. A
regex over REPORT.md finds no controlled claim word plus a gene symbol on a line without a `c:` tag.

### R1-AF6 — major — The Part 1 read-only test misses the likely write paths
**Problem.** The proposed test runs `git status` in ../project and checks that no file is newer than project start.
`git status` cannot see the realistic writes: `__pycache__` in `P1/src/essml/`, `cache/m2_inner/*.npy`,
`cache/m5_*` and `cache/models/` are all gitignored in Part 1. Also, Part 1 itself was edited on 2026-10-28, the same day
this project started, so "project start" needs an exact baseline. `essml.guard._log` writes to
`P1/results/test_access_log.csv` on any non-dev request.

**Change (add to §1).**
1. At M3 start, `scripts/snapshot_part1.py` writes `configs/part1_manifest.csv` (path, size, mtime_ns) for every file
   under ../project except `data/` (whose sha256 is already in Part 1's configs). It also asserts that P1 `HEAD` equals
   `configs/part1_commit.txt` (currently 7f94b5a; the P1 tree is clean).
2. Only `src/biomark/p1.py` may import `essml`. It sets `sys.dont_write_bytecode = True` before the import, sets
   `m2_local.INNER_CACHE`, `m5_stack.INNER_CACHE/OUTER_CACHE` and the modelstore `MODEL_DIR` to project_b paths, and then
   asserts that every one of them resolves inside project_b. It must never call `essml.guard.load_targets` with non-dev
   IDs.
3. No `subprocess` call may target a path under ../project.

**Check.** `tests/test_part1_untouched.py` compares the live tree with the manifest: no new, removed or changed
(size/mtime) file outside `data/`, and HEAD unchanged. The leakscan (R1-AF1) also greps for `import essml` outside p1.py
and for `subprocess` plus `project\\scripts`. Run it at every gate, not only in pytest.

### R1-AF7 — major — Prose promises with no output, owner or fate (Part 1 lost 40 of 174 accepted items this way)
**Problem.** DESIGN v1 has about 40 deliverables ("flagged as possible culture artefacts", Leiden "unnamed" clusters, LMM
sensitivity on the top 1,000, TCGA prevalence via cBioPortal, per-line combinatorial concordance, the agreement metric,
the leave-one-lineage-out sign check, the MSS-only persistence check, and others). None of them has an output file or a
gate. Some depend on things that are absent:
- `leidenalg`/`igraph` are not installed (checked);
- data_availability.md says scikit-learn is absent, although it does import, so the environment record is unreliable;
- the combinatorial screen tables have not been requested;
- DESIGN cites `LITERATURE.md`, which does not exist, although PROMPT.md requires it.

**Change (add §9.1 "Compliance ledger").** `configs/design_ledger.csv` has one row per deliverable: (id = §-number +
letter, output path, milestone, status ∈ {planned, done, declined: reason, not done: data unavailable}). The ledger is
created at `freeze-design` and its status column is updated in place. Review items accepted in rounds 1–3 get rows too,
with the review ID (LESSONS #5). REPORT.md includes the ledger summary table. `pip freeze` is written to
`configs/pip_freeze.txt` at freeze, and packages missing for a deliverable are installed or the deliverable is declined.

**Check.** `gate.py ledger`: every `done` row's output path exists, is tracked by git and is non-empty. No row is still
`planned` at the final gate. Every file named in DESIGN.md as an input (regex over backticked paths) exists. Every R1–R3
"accepted" ID from reviews/ appears in the ledger.

### R1-AF8 — major — Negative controls and nulls have no numeric bounds ("near 0")
**Problem.** §7.2 says "expected MD size near 0" and "fidelity near 0", and the toy test says "empirical FDR ≤ 0.1" with
no tolerance for how many toy replicates. Part 1 had two permutation controls miss bands that were never simulated, and
it "accepted them as interpreted" (LESSONS #10). A "near 0" bound can be redefined after the fact.

**Change (add to §7.2).** Before `freeze-design`, simulate each negative control on the toy generator (R1-AF9) at real D
dimensions. Write the 95th percentile of the null statistic as the pass bound in `design_params.json`:
- MD size under within-lineage permutation;
- surrogate fidelity for matched-variance low-tier genes;
- empirical FDR over 20 toy nulls.

A real-data control outside its bound is reported as a failed control. The bound is not reinterpreted.

**Check.** `gate.py controls` compares `results/controls/*.json` against the bounds in `design_params.json` and records
pass or fail. The bounds' commit must be an ancestor of the first real-data control run (git log order).

### R1-AF9 — major — Verification can become self-confirming
**Problem.** The toy test (§9) and the "agent that did not produce it checks it" rule are both stated, but:
- an agent can write a toy test that exercises a toy-specific code path or config;
- an agent can tune the pipeline until the toy passes;
- an agent can "check" a result by rereading the producer's output instead of recomputing it.

**Change (add to §9).**
1. The toy test calls the same entry points with the same `design_params.json`; the test asserts the hash. The planted
   effects and the toy seed are fixed at `freeze-design`.
2. Before `freeze-catalogue`, a separate agent reimplements five headline quantities from raw inputs with its own code,
   in `tests/independent/`, and a script compares them with the pipeline (tolerance 1e-6, or exact for counts). The five
   are: the MD list for 20 genes, the S1 p-values for those genes, the positive-control gate recall, the BH counts, and
   the confounding delta for 10 pairs. Do this once (LESSONS #8). After that, the reviewer effort goes to red-team
   questions.
3. The outcome is recorded in `results/gates/G_indep.json` with the agent name.

**Check.** `gate.py G_indep`: the comparison file exists, the maximum deviation is within tolerance, and the author of
`tests/independent/` (git log `--format=%an` or a trailer `Agent:`) differs from the author of `src/biomark/`.

### R1-AF10 — minor — Late or missing user downloads can silently change which analysis is primary
**Problem.** Several covariates depend on files the user may not provide: MSI, doubling time, full mutation calls,
segments and RRBS. The deadline is "start of M4". If a file arrives during M6, an agent may swap the primary analysis
after results have already been seen, or fall back to a proxy without saying so.

**Change (add to §10).** At the start of M4, `scripts/data_status.py` writes `results/data_status.json`, recording each
requested file as present (with sha256) or absent. Every analysis chooses between the real file and the proxy by
reading this JSON. A file that arrives after the deadline may only feed secondary analyses, which are labelled as such.

**Check.** `gate.py data`: the JSON is committed before the first M4 results commit, and the primary result files record
the `data_status` hash they used, which must equal the frozen one.

### R1-AF11 — minor — Run provenance and decision memory
**Problem.** "Every run logs config hash, commit, input checksums" has no registry or format. Agents running parallel
tracks (M4/M5) will collide on shared logs and decision IDs (P1 LESSONS #12: registry corruption, D21 renumbered to
D23). After a context reset, an agent cannot tell which result was produced by which rule version.

**Change (add to §10).**
- `results/registry.csv` (file-locked append) has columns run_id, script, commit, dirty flag, config hash,
  design_params hash, data_status hash, status and outputs.
- Results produced from a dirty tree get `status = exploratory` and cannot be cited in REPORT.md.
- Each parallel track gets its own range of D-IDs (M4: D100–, M5: D200–, and so on) and its own output folder.

**Check.** The report-number check (R1-AF5) resolves each tagged file to a registry row with `status = ok`, dirty
flag = 0 and the frozen design_params hash.

### R1-AF12 — minor — Hard-coded gene, line and locus exceptions
**Problem.**
- §5.4 lets the agent name a mechanistic gene in a locus "with literature support (e.g., ERBB2)".
- 9p21 is named in prose.
- The positive and negative control lists live in the text.

Ad-hoc exceptions like these are where agents rescue favoured hits, and literals in the code break the rule against
hard-coding data.

**Change.** The cutting-artefact exceptions (`configs/cutting_exceptions.csv`: gene, citation key that resolves under
R1-AF3) and any named loci are config files frozen at `freeze-design`. Code must contain no gene symbols, Entrez IDs or
`ACH-` IDs.

**Check.** A leakscan extension greps `src/` and `scripts/` for `ACH-\d{6}` and for any symbol from
`P1/cache/gene_map.parquet` used as a string literal (whitelist: `tests/`, `configs/`). The `gate.py frozen` check
covers the exception files.

---

**Summary for the PI.** Items AF1–AF3 should be in place before any association is computed (M3). AF4, AF7 and AF8 must
be finished before `freeze-design`, because each needs a dry run on non-blinded inputs. AF5 and AF9 are needed before
the report. None of these items changes the scientific scope. They turn §§1, 7 and 9 of the design into checks that a
script can enforce.
