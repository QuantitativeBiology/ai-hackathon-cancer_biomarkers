# Toy NC2–NC6 bounds: definitions (written before any replicate was run)

Author: toy verification agent (not the pipeline author). DESIGN v3.2 §7.3, §13.1, amendment A6;
`design_params.json: negative_controls_s7`, `toy_s13.nc_bound_rule`; timebox `tiers.json: tiers.A.stages.toy_run_bounds`
= 45 min. Script: `scripts/toy_nc_bounds.py`; output: `results/toy/nc_bounds.json`.

This file was written before `scripts/toy_nc_bounds.py` was first run (its sha256 is recorded in the output JSON). Every analogue, unit,
denominator rule and fallback below is fixed here; none was chosen after seeing a value.

## 1. Replicates and runs

- Seeds `toy_seed + r`, r = 0..4 (`seeds.json`), scenario `truth`, toy dimensions (`toy_s13`: 600 lines × 12 lineages,
  300 genes, 400 shared features). Predictor = the toy default `smooth` (amendment A6-i: recovery-type bounds use the
  smooth toy; the ridge toy is the unbounded error-rate check). NC bounds are not error-rate bounds.
- Engine only (`biomark.pipeline.run_pipeline`, `biomark.arms.run_arms`); no oracle numbers enter a bound. Per replicate:
  - **R1** `run_pipeline(toy_r, fold_centre=True, null_B=0)` (all four arms, MD-strict, audits). No NC2–NC6 statistic
    needs the null.
  - **R2** low-tier universe: `Data.from_toy(toy_r, fold_centre=True, genes=<low-tier genes>)`, then
    `arms.run_arms(data, arms=("S1-obs", "EN-pred"))` (EN-pred fidelity and MD-strict with the D-stage screen). This is a
    separate universe like `run_arms.py --universe ctrl`: MD-strict's BH q_sc is computed within the low-tier universe's
    own EN-pred selections, as a separate real NC2 run would.
  - **R3** `run_pipeline(toy_r, fold_centre=False, null_B=0, arms=("EN-pred",))` (exactly the answer-key test's
    no-centring run; NC4 only).
- Workers: the 5 replicates run as 5 concurrent processes with 2 engine workers each (10 in total).
- **Timebox**: wall-clock from launching the 5 replicates (R1–R3 of all r) to the last one finishing. If it exceeds
  45 min, every NC2–NC6 provisional value becomes final, the toy values are still reported, and one row is appended to
  `results/deviations.csv`.

## 2. Faithfulness rule (applied mechanically by the script)

An NC receives a toy-derived bound only if both hold; otherwise it keeps its provisional value, the toy values are
reported alongside, and a row naming the reason is appended to `results/deviations.csv`.

- **F1 (same unit, same code path).** The toy contains units of the real statistic's kind (genes, pairs or CGT tests of
  that class) and the statistic is computed from the same engine outputs the real run will use.
- **F2 (resolution).** For a rate or fraction with per-replicate denominator n_r, 1/n_r ≤ the provisional value (for
  the min-bound NC6: 1/n_r ≤ 1 − provisional) in every replicate. Otherwise the max over replicates can only be 0 or
  ≥ 1/n_r: it is quantised, not an estimate of the statistic's null spread. Continuous statistics satisfy F2. A
  per-gene count (hits / genes) uses n_r = number of genes.

Bound = max over the 5 replicates (min for NC6), as `nc_bound_rule` states. No headroom is added.

## 3. Analogues

| NC | Real statistic (DESIGN §7.3) | Toy analogue | Unit and denominator | Provisional |
|---|---|---|---|---|
| NC2a | 100 matched-SD low-tier genes (`lowtier_seed`): median fidelity | median over the toy's `low_tier` plant genes of EN-pred `r2_oof_beyond_design` (R2 `lam`) | genes; continuous | ≤ 0.02 |
| NC2b | same genes: MD-strict hits per gene | number of MD-strict rows (R2) / number of low-tier genes | genes, n_r = 10 | ≤ 0.1 |
| NC3 | self-null features: fraction of testable pairs whose carriers are dependent (carrier median effect < −0.5) | mean `frac_dependent` over `testable` rows of the engine's `self_null` table (R1) restricted to the answer-key self-null set (`answer_key["self_null"]`: the planted self-null gene's DEL_g and its DEL locus mates) | pairs; n_r = testable planted self-null pairs | ≤ 0.1 |
| NC4 | fold artefact (toy plant) | fold-artefact features EN-pred-selected after centring (R1) and without centring (R3) | 5 plants | existing toy pass bounds: after ≤ 1, without ≥ 3 |
| NC5 | false-sensitivity rate on resistance/N CGT tests | per arm: toy CGT tests with partition `bypass-resistance` or sign `N` whose `feature_D` is ranked in the arm; false sensitivity = hit (`rank_of_feature_D` ≤ k = 200, engine `cgt_scores`) with β < 0; statistic = max over the four arms | tests; n_r = number of such tests | ≤ 0.2 |
| NC6 | lineage-TF genes: fraction in class 3 | fraction of the toy's `lineage_only` plant genes with engine `gene_classes.failure_class == "3"` (R1) | genes, n_r = 10 | ≥ 0.7 |

Notes fixed in advance:

- **NC2.** The toy has exactly 10 low-tier genes; all are used, so SD matching and `lowtier_seed` have no toy analogue
  (there is nothing to sample from). Toy low-tier genes carry hidden-latent signal that the predictor sees but no 𝒟
  feature represents. They model "no nameable determinant", which is what both NC2 statistics test. They do not model
  "unpredictable". Reported alongside, not bounded: per-gene fidelity, the number of EN-pred selections, and the
  MD-strict count obtained with R1's observed p*_sc instead of a within-universe BH.
- **NC3.** Only the planted self-null gene follows the self-null biology in the generator (its DEL carriers have true
  effect N(0, 0.05)). For every other toy gene, DEL_g/LOFBI_g/SIL_g carriers keep the gene's ordinary effect. An
  all-gene toy statistic would measure the generator's gene-mean distribution rather than the pipeline, so it is
  reported as secondary only (all genes and eligible genes), not bounded.
- **NC4** is a toy-only control. Its bounds are the existing `toy_s13.pass` / answer-key values; per-replicate
  pass/fail is reported. It is not re-derived as a max over replicates.
- **NC5.** The toy CGT has one `bypass-resistance` W test (the p53_w plant) and at most one `expected-null` N test (a
  background HOT feature). Per arm values are reported.
- **NC6.** Lineage-only toy genes have strong lineage effects in 1–3 lineages (one with a lineage-proxy SIL feature) and
  nothing else. They are the analogue of the P1 lineage-TF rows.

## 4. Output (`results/toy/nc_bounds.json`)

Per NC: `definition`, `unit`, `per_replicate` (value, denominator, extras), `toy_bound` (max, or min for NC6),
`provisional`, `faithful` (F1, F2 with the failing replicates), `final_value` (toy bound if faithful and inside the
timebox, else provisional), `final_source` (`toy` / `provisional`) and `reason`. Top level: seeds, commit,
config sha256s, timebox and measured wall-clock, workers, and a pointer to this file.

## Addendum 2026-10-28 (≈22:45 UTC): coordinator rulings R1 and R2, decided before ridge NC2 values were computed

Context: the first run (engine 8cffc91, `results/toy/nc_bounds_8cffc91.json`) gave degenerate smooth-toy NC2 values.
EN-pred selected 0 features on every low-tier gene, so median fidelity ≈ 0⁻ and MD-strict per gene = 0. NC6 = 1.0 in
every replicate. When this addendum was written, no ridge-toy NC2 value had been computed. §1–§4 above stay unchanged
except where these rulings override them.

- **Ruling R1 (NC2 on the ridge toy).** NC2a/NC2b are error-type controls. Per amendment A6, error rates are judged on
  the ridge toy, which is also the analogue of Part 1's ridge-type predictor. The NC2 toy statistic is therefore the
  max over 5 ridge-toy replicates (`make_toy(toy_seed + r, predictor="ridge")`, r = 0..4) at toy dimensions, with the
  same R2 run (low-tier universe, `run_arms(("S1-obs", "EN-pred"))`, BH q_sc within the low-tier universe) and the same
  statistics as §3. The smooth-toy NC2 values are reported alongside, not bounded. NC3–NC6 stay on the smooth toy.
- **Ruling R2 (general degeneracy rule).** This applies to every NC2–NC6 statistic, including NC6 and any NC on the
  post-fix run. A statistic has no resolution, and its provisional value becomes final (with a `DEV-TOY` row in
  `results/deviations.csv` reporting the toy values, as for NC3/NC5), if either holds:
  - it has zero between-replicate variance and sits at the boundary of its range (0 or 1); or
  - its underlying count is 0 in every replicate.
  Underlying counts: NC2a = EN-pred selections on the low-tier genes; NC2b = MD-strict rows; NC3 = dependent testable
  pairs; NC5 = false-sensitivity hits; NC6 = class-3 genes. R2 is applied mechanically by `scripts/toy_nc_bounds.py`
  after the F1/F2 rule of §2. If ridge NC2 is also degenerate, R2 applies to it.
- **Labelling.** Every NC output names the engine commit it ran on (`nc_bounds_<commit>.json`). Deviation-row ids carry
  the commit.
