*(Human-written text starts)*

This is an exercise in assessing the bioinformatics capabilities of agentic AI. It is part of a hackathon run by the Beltrao Lab at ETH Zurich, to which we (INESC-ID & IST – Universidade de Lisboa) contribute, in which Claude was used to run bioinformatics projects (https://github.com/evocellnet/ai-hackathon/). All results in this project were generated entirely by AI, with minimal human intervention. This repository is not intended to be read as a traditional publication or used as a validated bioinformatics tool.

The goal of this project was to use a CRISPR gene essentiality predictor (built in a previous project: https://github.com/QuantitativeBiology/ai-hackathon-cancer_essentiality) to identify molecular biomarkers that explain why groups of cancer cell lines depend on specific genes, and to validate them against established clinical biomarkers. It was executed by Claude Opus 5.5, with delegation to Opus 5.5 and Sonnet 5.5 subagents.

None of the results have yet been verified by a human. We plan to verify the agent's methodology and benchmarking, and will update this section accordingly.

The rest of this repository is entirely AI-generated. It includes the code implementing the analysis and an AI-written report in the style of a journal article, describing the findings, the path taken to reach them, and benchmarking against other methods.

*(Human-written text ends)*

---

# Predictor-derived dependency subtypes and their molecular biomarkers

This is Part 2 of a two-part project. Part 1 built a CRISPR gene-essentiality predictor on DepMap Public 26Q3. Here we try to read molecular determinants out of that predictor and test whether they work as biomarkers. A determinant is a mutation, copy-number or methylation event, or an expression feature.

We use the determinants to group cancer cell lines into gene-dependency subtypes, within and across lineages. We then test the subtypes in three ways:
- against single-feature association tests;
- against combinatorial CRISPR screens;
- against clinical biomarkers and their companion trials.

## Where to start

| File | What it is |
|---|---|
| [REPORT.md](REPORT.md) | The findings in full: methods, verification, results, deviations and limitations. Every number cites its source file. |
| [RESEARCH_TREE.md](RESEARCH_TREE.md) | The research log: status summaries at each milestone and decisions D0–D22. |
| [DESIGN.md](DESIGN.md) | The pre-registered design (v3.2, three review rounds), plus the amendments A1–A7 made after the freeze. |
| [docs/effective_rules.md](docs/effective_rules.md) | One table listing every rule in force at the end, with its source. |
| [results/deviations.csv](results/deviations.csv) | Every departure from the frozen design, with its reason. |
| [LITERATURE.md](LITERATURE.md), [DATA.md](DATA.md) | The literature survey, and the data sources with versions, checksums and licences. |
| [LESSONS.md](LESSONS.md) | Process rules learnt while running the project. |
| [PROMPT.md](PROMPT.md) | The original task prompt, verbatim. |

## Main findings

These are taken from REPORT.md, which cites the source of each number. They come from the dev set D (1,022 lines) and from one replication run in 255 held-out lines (R).

- **The predictor-derived catalogue is not better than simpler methods.** We counted how many of each method's top 100 genomic candidates replicate in R:

  | Method | Replicated in R |
  |---|---|
  | Univariate confounder-adjusted test on observed effects (S1-obs) | 0.99 |
  | Elastic net on observed effects (EN-obs) | 0.89 |
  | Elastic net on predicted effects (EN-pred) | 0.82 |

  - **Predicted vs observed elastic net:** the difference is −0.07 (95% CI −0.12 to 0.02), so it is not detectable. This verdict is fragile: a basic bootstrap interval is entirely negative, and EN-pred replicates better than EN-obs at 200 and 500 candidates.
  - **Elastic net vs univariate test:** the univariate test replicates better, by 0.10 (CI 0.05 to 0.17). Replication across libraries points the same way: 0.82 vs 0.62.
- **The positive-control gate failed.** The pipeline recovered 21 of 39 known controls, against a pre-declared bar of 60%, at both attempts.
- **The predictor-derived lists recover almost no clinical biomarkers.**
  - **Primary clinical endpoint:** only 6 tests were eligible, below the minimum of 15, so the result is descriptive.
  - **Detectable biomarkers:** of the 11 clinical biomarkers detectable in cell-line data, EN-pred recovers 0 and S1-obs recovers 3.
  - **Refuted biomarkers:** none is recovered by any method.
- **The combinatorial CRISPR screens confirm the paralog candidates, but don't favour the predictor.** Catalogue candidates replicate at 16.9%, against 2.5% expected for matched paralog pairs, which is the same rate as the univariate list.
- **Shortlist: 552 candidates that replicated in R with adequate power.** Most are expression-level paralog dependencies (for example FAM50A–FAM50B, DDX3X–DDX3Y and DNAJC19–DNAJC15), and several are also frequent dependencies in non-cancerous lines. REPORT.md §6.5 works through five subtypes, giving their member lines and lineage make-up.
- **Several pre-declared controls failed.** These are reported as failures:
  - the empirical FDR for EN-obs;
  - predictor fidelity on low-tier genes;
  - the self-null control, mostly from deletion features;
  - the lineage transcription-factor class check.

## Repository layout

```
src/biomark/    analysis library: loaders and the R-access guard, dictionary, confounder model and S1 tests,
                stability elastic net, four arms, null, endpoints, audits, limits classes
scripts/        one script per pipeline stage (build_*, run_*, replicate.py, report_counts.py) and the
                process checks (check_frozen.py, leakscan.py, budget_check.py, check_report_numbers.py)
tests/          unit tests; toy simulation with an answer key (test_toy_answers.py, written by a separate
                verification agent); tests/independent/ (independent reimplementation of the core statistics)
configs/        frozen parameters (design_params, primary_endpoints, tiers, seeds), the clinical ground-truth
                table, pre-committed outcome actions, frozen hashes (configs/frozen_hashes.json)
results/        committed result tables (plain CSV/JSON), with results/registry.csv logging every run
reviews/        design reviews (3 rounds), red-team and LESSONS reviews
docs/           implementation notes, literature notes, data availability, effective rules
cache/, external/   not committed: bulk intermediates and downloaded external data (git-ignored)
```

## Data

- **Part 1 must be present.** The code expects the Part 1 repository at `../project/`, with DepMap 26Q3 under `../project/data/`. Part 1 is read in place and never modified; `tests/test_part1_untouched.py` checks this.
- **Some DepMap extras were downloaded by hand.** These are OmicsSomaticMutations, OmicsGlobalSignatures, OmicsProfiles, CRISPRConfounders, OmicsCNSegmentsWGS and CCLE RRBS methylation. They go in `external/depmap_extra/`.
- **External clinical and combinatorial data are fetched by script.** `scripts/fetch_external.py` downloads CIViC, MOAlmanac, FDA companion diagnostics, OncoKB gene lists and the combinatorial screen tables into `external/`.
- **Checksums and licences** are in [DATA.md](DATA.md).
- **Nothing raw or restricted is committed.** No raw or restricted data, and no bulk intermediate, is in git. The sha256 of each bulk file is in `results/bulk_manifest.csv`.

## Reproducing the run

Python 3.13.1, CPU only. Install the pinned packages with `pip install -r requirements.txt`. Run all commands from the repository root with `python -X utf8 -B`. The machine used had 15.7 GB of RAM, and `configs/runtime_params.json` caps workers to fit memory.

The run is organised around four git tags. Scripts refuse to run before the tag they depend on.

| Tag | What it freezes |
|---|---|
| `freeze-design` | design parameters |
| `freeze-cgt` | the clinical ground-truth table, frozen before any dependency data were joined |
| `freeze-bounds` | negative-control bounds, dictionary, controls, null genes |
| `freeze-catalogue` | catalogue and replication requests, frozen before the single read of R |

1. **Inputs and builds:**
   - `verify_inputs.py`, `fetch_external.py`, `data_status.py`
   - `build_screen_covariates.py`, `build_eligibility.py`, `build_gene_lists.py`, `build_effect_nan_mask.py`, `build_subsamples.py`
   - `build_dictionary.py`, `build_positive_controls.py`, `fetch_cgt_records.py`, `build_cgt.py`
   - `build_null_genes.py`, `build_u_ctrl_genes.py`, `build_lowtier_genes.py`
2. **Arms and null on D:**
   - `run_arms.py`, `run_arms.py --universe ctrl`, `run_arms.py --universe lowtier`
   - `run_gate.py`, `run_null.py`, `run_confounder_audit.py`, `run_component_readout.py`, `run_xlib.py --cut-discovery-only`
3. **Subtypes, audits, clinical and combinatorial:**
   - `run_subtypes.py`, `run_cross_lineage.py`, `run_m5_audits.py`
   - `run_q4.py`, `run_ceiling.py`, `run_cgt_secondary.py`, `run_combinatorial_lookup.py`
4. **Freeze and single replication:**
   - `build_p1_pool.py`, then the dry run `replicate.py --pseudo-r fold0`, then tag `freeze-catalogue`.
   - Then `PART2_UNBLIND=1 replicate.py` on a clean tree. The guard allows exactly one logged read of R; a reproduction needs `--reproduce`, with the same commit and config.
5. **After R:**
   - `p1_secondary.py`, `run_r_power.py`, `run_nc_real.py`, `run_limits.py`, `build_candidate_cards.py`
   - `report_counts.py`, `check_report_numbers.py REPORT.md`

The exact commands, worker counts and the order that was actually run (including interruptions and re-runs) are in RESEARCH_TREE.md and `results/registry.csv`.

## Verification

- **Unit and process tests:** `python -X utf8 -m pytest tests`. The tests `test_guard.py` and `test_loaders.py` read dev data.
- **Toy simulation with an answer key:** `BIOMARK_TOY_PIPELINE=1 python -X utf8 -m pytest tests/test_toy_answers.py tests/test_toy_md_global.py -k pipeline` (about 1 h). On the final engine, 19 checks pass and 3 fail; REPORT.md §5.2 describes the failures.
- **Independent reimplementation:** `python -X utf8 -m tests.independent.run_all`. The results are in `results/gates/G_indep.json`.
- **Integrity checks:** `scripts/check_frozen.py` (frozen hashes and tag order), `scripts/leakscan.py` (no hard-coded genes, line IDs or thresholds) and `scripts/check_report_numbers.py` (every number in REPORT.md against its cited file).

## Known limitations of this repository

- **Not run or not computed:**
  - some pre-declared analyses were cut under the design's cut order or not computed: TCGA prevalence, the predictor-attribution panel, the selection-adaptive null, and several P1 secondary variants;
  - the non-Avana subset of R, which would have needed a second read of R.
- **Strict number check:** `check_report_numbers.py --strict` fails on four result files that were written before the run registry existed.
- **Post-R items:** some D-only items were completed after the R read, using definitions frozen beforehand (DEV-POST-R-COMPLETION).
- **The full list** is in REPORT.md §9 and `results/deviations.csv`.
