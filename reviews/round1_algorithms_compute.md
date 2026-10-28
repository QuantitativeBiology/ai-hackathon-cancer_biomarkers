# Review round 1: algorithms and compute (R1-AC)

Reviewer scope: whether the surrogate and stability selection, the S1 tests, the permutation null, the Part 1 refit with direct attribution, Leiden clustering and memory are feasible on 12 cores, 16 GB, CPU only, with about 20 h left for the whole project.

## Measurements behind this review

All of these are read-only on Part 1 or use synthetic data in the scratchpad.

- **Eligible genes.** (variable or selective) and m5 tier in {high, medium} gives 3,062 genes (1,331 high, 1,731 medium), counted from `P1/configs/gene_classes.csv` and `predictable_genes.csv`.
- **Raw 𝒟 family sizes before collapse or restriction** (dev lines, at least 5 carriers):
  - HOT: 45;
  - damaging-mutated genes: 4,375 (before the OncoKB/prior restriction);
  - AMP (CN > 2): 11,776 genes;
  - DEL (CN < 0.25): 261 genes;
  - SIL candidates (expressed in at least 50% of lines, not expressed in at least 5): 5,207.
  - Unless these are capped, p will be 3k to 10k.
- **Component OOF predictions.** M5 OOF equals the mean of the three component OOF runs (m1_E_resB0, m2n_E_paralogs, g3_m3n). The largest difference is 9.5e-7. So all component OOF predictions already exist, and a refit is needed only to get model internals.
- **Elastic net benchmark** (sklearn 1.5.2, one core, n = 1022, p = 3000, l1_ratio 0.9):

  | Setting | Time per gene | Total for 3,062 genes |
  |---|---|---|
  | Default `ElasticNetCV` (5 folds, 50 λ) | 51 s | about 43 core-h, or 3.6 h wall on 12 cores |
  | Same, eps 1e-2 and 25 λ | 24 s | |
  | Shared Gram matrix | 8 s | |
  | Marginal screen to the top 300 features, then CV | 0.15 s | |
  | 50 half-sample fits at a fixed λ, with screening | 0.6 s | |

  - A single half-sample fit at a fixed λ takes 0.07 s with all 3,000 features.
  - Screening recovered 8 of 10 planted features at the 60% stability threshold.
  - **The cost is in the per-gene λ cross-validation, not in the 50 subsamples.**
- **Environment.**
  - leidenalg and igraph are **not installed**. networkx 3.4.2 and lightgbm 4.7.0 (native `pred_contrib`) are present. shap is not installed and is not needed.
  - `P1/cache/m2_inner` is 3.8 GB.
  - Part 1 all-dev refit of M2 without the inner cache took 1,384 s, and M5 peak RSS was 3.3 GB (`P1/results/timing/refit_full_dev.json`).

## Items (most important first)

### R1-AC1 (critical): Make the surrogate algorithm explicit and fast, and fix the wrong fallback
- **Problem.** With default sklearn `ElasticNetCV` on p ≈ 3,000, §4.1 costs about 3.6 h wall, and more if p is larger.
  - The declared fallback ("if more than 4 h, use 20 subsamples") cuts the cheap part. The subsamples cost about 4 s per gene; the per-gene λ CV costs about 50 s.
  - sklearn cannot leave lineage unpenalised directly.
  - The rule for λ inside the subsamples is unspecified. Re-running CV in each subsample would multiply the cost by 50.
- **Change.** Add to §4.1:
  1. **Lineage.** Handle the unpenalised lineage term by Frisch–Waugh–Lovell (FWL). Standardise the 𝒟 columns once on D. Then, inside each fit (full D, each CV fold, each subsample), residualise ŷ_g and 𝒟 on lineage indicators estimated on that fit's rows, and fit the elastic net with no intercept. This is exactly equivalent to an unpenalised lineage term with the stated penalty. Drop a column whose residual SD is below 1e-6 or which has fewer than 5 carriers in the fit's rows.
  2. **Screening.** Inside every fit, screen to the K = 300 features with the largest |corr| with the lineage-residualised ŷ_g. The gene-specific EXP-own and EXP-paralog terms and the pre-declared partner feature of each positive-control pair always pass the screen. They are still penalised.
  3. **λ path and CV.** Use 25 λ values with eps = 1e-2 and tol = 1e-3. The CV folds are Part 1's `cv_fold`, which are already patient-grouped. Use the 1-SE rule.
  4. **λ in subsamples.** The full-D λ_1SE is held fixed in all subsamples (sklearn's alpha is per-sample, so it transfers).
  5. **Fallback.** Replace the fallback with this rule: time a pilot of 50 random eligible genes first. If the projected wall time is above 45 min, set K = 150. The number of subsamples is never reduced.
- **Expected cost.** About 1 s per gene, so 3,062 genes take about 5 min wall on 12 cores.
- **Check.**
  - `results/timing/surrogate_pilot.json` holds the projected and actual wall time.
  - A unit test on toy data compares FWL elastic-net coefficients with a fit in which the lineage columns are scaled by 1e4, so they are effectively unpenalised. The maximum |Δβ| should be below 1e-3.
  - The config hash includes `K`, `n_lambda`, `eps` and `lambda_rule = full_D_1se_fixed`.

### R1-AC2 (critical): 200 end-to-end permutations are infeasible; use B = 20 with a joint row permutation
- **Problem.** One end-to-end pass takes about 10 to 15 min wall, even after AC1 and AC8: surrogate plus stability selection about 5 min, vectorised S1 about 5 min, and the MD confounder filter. So 200 passes take 33 to 50 h. With the design as written (default ElasticNetCV) they would take more than 700 h.
  - "d permuted within lineage" done feature by feature would destroy the correlation between features. The locus collapse and the co-location tests would then behave differently under the null than on real data.
- **Change.** Rewrite §5.4 "Multiple testing" and §7.2 as follows:
  - **Permutation.** Each permutation draws **one** within-lineage permutation of the line labels and applies it to all rows of 𝒟 jointly, including the gene-specific EXP terms. Lineage, the covariates, ŷ and the observed effects stay attached to the original lines.
  - **Scope.** Run the full pipeline: dictionary features already frozen, surrogate, stability selection, the MD filter, and S1 with BH.
  - **B = 20.** The 20 seeds are fixed in `configs/seeds.json`.
  - **Output.** Report the mean null count of MD and S1 hits ± SE, and the empirical FDR (mean null count divided by the observed count). The goal is to estimate a mean count, so 20 draws are enough. The null does not give per-candidate permutation p-values.
  - **Scheduling.** Run it as a background job with a 3 h wall cap.
  - **Fallback.** If the pilot shows more than 9 min per permutation, run B = 20 on a seeded stratified subset of 600 eligible genes (tier × class) and report the counts per gene tested.
- **Check.**
  - `results/null/perm_manifest.json` lists the seeds, B, the gene set and the wall time per permutation.
  - A test asserts that a permuted 𝒟 keeps the column correlation matrix (permutation invariance) and keeps the lineage counts.
  - The toy full-null test uses the same code with B = 20.

### R1-AC3 (major): Cap and specify the dictionary size
- **Problem.** The raw counts are SIL 5,207, AMP 11,776 genes before collapse, and 4,375 mutated genes. Without explicit caps, |𝒟| can reach 6k to 10k. That raises surrogate and S1 cost and the test burden.
  - If a Gram-based solver is used, each worker needs 0.3 to 0.8 GB, or 4 to 10 GB across 12 workers. That is a risk of running out of memory.
  - The locus-collapse algorithm is not specified, and gene coordinates are not named. They are available from `data/raw_external/ensembl116_gtf`.
- **Change.** Add to §3:
  - **SIL.** Restrict SIL to the same prior gene set as LOF: OncoKB, plus paralogs and partners of eligible genes, plus known silenced genes such as MGMT, MLH1, CDKN2A, SMARCA4 and VHL.
  - **Locus collapse** (deterministic):
    1. Sort genes by GRCh38 start within each chromosome (Ensembl 116).
    2. Chain adjacent AMP (or DEL) features while the Pearson r of the 0/1 calls with the current locus seed is above 0.8 and the span is at most 10 Mb.
    3. The locus feature is the majority call.
  - **Ceiling.** The shared part of 𝒟 must have at most 3,000 features. If it is larger, raise the minimum frequency of the largest family to 10, then 20, and log it.
  - **Logging.** Log the per-family counts before and after each step.
- **Check.**
  - `scripts/build_dictionary.py` writes `results/dictionary/family_counts.csv` and asserts the 3,000 ceiling.
  - A unit test asserts that CDKN2A, CDKN2B and MTAP map to one locus, and that the ENO1 locus at 1p36 is collapsed.

### R1-AC4 (major): Reproduce the refit against the component runs; copy the inner cache; time-box the refit
- **Problem.**
  - §4.2 says "refit reproduces Part 1 OOF predictions … else stop and investigate". That has no time box, and the reference is left implicit.
  - LightGBM run with a different thread count may not reproduce to 1e-4.
  - Without `m2_inner`, M2 needs about 1,400 s per fit, about 2 h for 5 folds.
- **Change.** Edit §4.2:
  1. **Reference.** The reproduction target is each component's OOF `pred.parquet`:
     - m1: `20261026-222117_m1_E_resB0_91e52f9d`;
     - m2: `20261027-044450_m2n_E_paralogs_62ce1575`;
     - m3: `20261027-014824_g3_m3n_601e3f5e`.
     Their mean equals the M5 OOF to within 1e-6, which was verified.
  2. **Inner cache.** Copy `P1/cache/m2_inner` (read-only source) into `cache/m2_inner` here, which is gitignored. Record a SHA-256 manifest and set `m2_local.INNER_CACHE` to it.
  3. **Configuration.** Use the same n_jobs and thread counts as the Part 1 configs.
  4. **Tolerances.**
     - M1 and M2: max |Δ| ≤ 1e-4.
     - M3: max |Δ| ≤ 1e-3 and per-gene r ≥ 0.9999.
  5. **Expected cost.** M1 is 5 × about 80 s. M2 is 5 × about 3 min with the cache. M3 is 5 × 3 to 10 min. The total is at most 75 min wall, run in the background during M4. It must not run at the same time as 12-worker jobs.
  6. **Time box.** If the refit is not reproduced after one diagnosis attempt of at most 30 min, §4.2 is reported as "not done: refit not reproduced", and the AC5 fallback is used. The project does not stop, because §4.2 is secondary.
  7. **Additivity check.** Every attribution must satisfy: sum of contributions + bias = component prediction, within 1e-5, for each (line, gene).
- **Check.**
  - `results/refit/repro_check.csv` holds max |Δ| and r for each component and fold.
  - `tests/test_attribution_additivity.py`.
  - `cache/m2_inner/MANIFEST.sha256`.

### R1-AC5 (major): A refit-free component-level readout, used as complement and fallback
- **Problem.** The question "which part of the predictor uses d?" can be answered largely without a refit, because the component OOF predictions exist.
- **Change.** Add §4.3:
  - Run the §4.1 surrogate, without stability subsampling (full-D fit only), on ŷ_M1, on ŷ_M2 − ŷ_M1 and on ŷ_M3 − ŷ_M1 (the genomic and paralog layers). Each is fold-centred. This costs about 3 × 1 min.
  - A determinant is "expression-mediated" if it is selected for M1, and "genomics-layer" if it is selected only for the M2 or M3 increments.
  - This is reported always. If AC4 fails, it replaces the direct-attribution agreement metric. The fallback is pre-declared.
- **Check.** `results/surrogate/component_selection.csv` exists, with columns `gene, feature, M1_sel, M2inc_sel, M3inc_sel`.

### R1-AC6 (major): Specify the attribution mapping so the attribution tensor is never materialised
- **Problem.** M1 contributions per (line, gene, feature) would be 1022 × 3062 × 18,071 ≈ 5.7e10 values, or 226 GB in float32.
  - The "projection of w_g onto the NMF basis" mixes spaces. The NMF is on 5,000 raw log-TPM genes, while M1 works on about 18k z-scored features.
  - The OOF rule for attributions and the fold-centring are not stated.
  - Pushing the M3 `m1` SHAP value down through M1 is an approximation.
- **Change.** Add to §4.2:
  1. **OOF rule.** Each line's attribution comes from the fold model that did not train on it. Contributions are fold-centred: the fold mean contribution per term is subtracted, so the terms sum to the fold-centred prediction.
  2. **Block computation for M1.** Per fold, compute W in gene blocks of at most 500 genes (18,071 × 500 float32 = 36 MB). Aggregate on the fly into these terms:
     - EXP-own;
     - EXP-paralog (3);
     - PROG_k, defined as the least-squares coefficients of w_g on the NMF loadings mapped to z-space (H_k / sd_train), over the NMF genes;
     - EXP-other (residual).
  3. **Unnameable genes.** Report, per gene, the share of M1 contribution variance captured by own + paralog + PROG. A share below 0.5 means "M1 unnameable", which feeds §8 class 3.
  4. **M3.** Use native LightGBM `pred_contrib=True`. Map the features as follows:
     - E: EXP-own, EXP-paralog or EXP-other;
     - C: the locus of that gene;
     - M: LOF, or LOF-other;
     - H: HOT;
     - meta one-hot: "lineage/meta", which is not a determinant;
     - the `m1` feature: a separate "M1-mediated" bucket, with no push-down.
- **Check.** A test asserts that the families plus the residual sum to the total M1 contribution (within 1e-4) and that SHAP rows sum to the booster prediction. Peak RSS is logged below 6 GB.

### R1-AC7 (major): Leiden is not installed and the attribution profile is undefined; cluster prediction profiles first
- **Problem.**
  - leidenalg and igraph are missing.
  - "Direct-attribution profiles over eligible genes" has no defined dimension.
  - The check depends on the AC4 refit succeeding.
  - Raw profiles will cluster trivially by lineage, and pooled OOF predictions carry fold offsets (the Part 1 D3 artefact).
- **Change.** Edit §5.2 "Unsupervised check":
  - **Primary input.** Fold-centred M5 OOF predictions over the eligible genes (1022 × 3062). No refit is needed. Run it twice: raw, and lineage-residualised.
  - **Method.** PCA to 50 components, then a cosine kNN graph with k = 15, then Leiden (`leidenalg` if `pip install leidenalg igraph` succeeds from wheels within 10 min; otherwise networkx `louvain_communities`, as a pre-declared fallback).
  - **Resolution and stability.** Resolutions {0.5, 1, 2} × 10 seeds. Report the ARI across seeds.
  - **Fold artefact test.** Test cluster membership against `cv_fold` with a chi-square test.
  - **Secondary input.** Attribution-family profiles (lines × genes × about 8 families), only if AC4 succeeds.
  - **Cost.** Seconds.
- **Check.** `results/clusters/stability.csv` (seed ARI), `ami.csv` (lineage, library, MD, fold) and `backend.txt`.

### R1-AC8 (major): S1 must be vectorised; define "LRT" as a partial F/t test; implement the LMM FaST-LMM style or drop it
- **Problem.**
  - S1 is about 3,062 genes × up to 3,000 features, or about 9M tests. Fitted one at a time with statsmodels (2 to 5 ms per fit, and an LRT needs 2 fits), that is 10 to 25 core-h per pass. Inside the null that is impossible.
  - Chronos NaN patterns differ by library, so the lines used vary by gene.
  - statsmodels MixedLM cannot take an arbitrary similarity matrix and would take hours for 1,000 tests.
- **Change.** Add to §6:
  - **Vectorised test.**
    1. Group genes by their NaN pattern.
    2. Within each group, residualise Y and 𝒟 on the shared covariate design Z by QR.
    3. Add the gene-specific own-CN covariate by a rank-1 FWL update.
    4. Compute t from the partial correlation, with df = n − rank(Z) − 2.
    5. Report F = t², which equals the Gaussian-OLS LRT.
  - **Drops.** Drop a test if fewer than 5 carriers remain among the non-NaN lines, or if the residual variance is below 1e-8.
  - **LMM.**
    1. Eigendecompose the line-similarity matrix once.
    2. Rotate.
    3. Estimate δ per gene on the null model by a grid search.
    4. Test by GLS.
    If it is not implemented by the end of M4, drop it (pre-declared).
  - **Restrict the per-candidate extras** (random-effects meta-analysis, LOLO, confounding delta per covariate block, conditional co-location tests) to MD ∪ S1 FDR hits ∪ controls. Do not run them for every (g, d).
  - **Cost.** Each S1 pass takes at most 5 min wall.
- **Check.** `tests/test_s1_vectorised.py` compares 200 random (g, d) pairs, including NaN-pattern genes, against statsmodels OLS. It requires |Δ log10 p| < 1e-6 and identical df.

### R1-AC9 (major): The time plan overruns; add a compute schedule and a pre-declared cut order
- **Problem.** Milestones M2 to M8 sum to 21 h, but about 20 h are left, including the review that is under way. Nothing says what is cut first.
- **Change.** Add to §10:
  - **Compute table.**

    | Job | Wall time |
    |---|---|
    | Dictionary + NMF | 5 min |
    | Surrogate (AC1) | 5 to 10 min |
    | S1 (AC8) | 5 min |
    | Refit + attribution | up to 75 min, background in M4 |
    | Null (AC2) | up to 3 h, background in M5–M6 while CGT curation runs |
    | Clustering | 1 min |
    | Replication | 10 min |

    Total heavy CPU is about 5 h of wall time, overlapped with work that is not CPU-bound.
  - **Cut order** if M4 ends more than 1 h late:
    1. LMM sensitivity;
    2. direct-attribution refit (use AC5);
    3. attribution-profile clustering;
    4. B = 20 → 10;
    5. TCGA prevalence and per-line combinatorial concordance.
  - **Never cut:** the positive-control gate, the confounder model, CGT pre-registration, the single R replication, and §8.
- **Check.** Every job writes `results/timing/<job>.json` with wall time and peak RSS. `scripts/budget_check.py` compares the times with the table and prints any cuts it triggers.

### R1-AC10 (minor): Pin the stability-selection mechanics and out-of-sample fidelity
- **Problem.**
  - Half-samples drawn independently give no error control.
  - Fidelity R² measured in-sample is optimistic, especially with p > n.
- **Change.**
  - **Subsamples.** Use 25 complementary pairs (Shah & Samworth), which gives 50 half-samples. They are patient-grouped and lineage-stratified, with seeds from `seeds.json`.
  - **Fidelity.** Fidelity is the out-of-fold R² from the λ CV, minus the lineage-only R², computed on the same folds.
- **Check.** `results/surrogate/subsample_manifest.json`: each pair is disjoint, the two halves of a pair cover D, and no patient is split. Fidelity in the per-gene table is labelled `r2_oof_beyond_lineage`.

### R1-AC11 (minor): Rules for parallelism, memory and determinism
- **Problem.** On Windows, 12 loky workers each using multithreaded BLAS oversubscribe the CPU. If seeds depend on worker scheduling, results will not reproduce.
- **Change.**
  - Use joblib loky with 12 workers over chunks of 50 genes.
  - Set `threadpool_limits(1)` inside the workers.
  - Pass 𝒟 and Y as read-only memmaps.
  - Derive per-gene seeds as `seed_list[i] ^ int(entrez)`.
  - Do not run the refit (4 to 6 processes, about 3.3 GB) at the same time as any 12-worker job.
  - The peak RSS ceiling is 12 GB, logged with psutil.
- **Check.** Re-run the surrogate on 20 genes with 1 and with 12 workers. The selection frequencies must be identical (`tests/test_determinism.py`).

### R1-AC12 (minor): Pin the NMF settings (PROG)
- **Problem.** A 50-factor NMF on 1022 × 5,000 takes about 1 to 3 min. With random initialisation the factors are not reproducible, and R lines need the frozen basis.
- **Change.**
  - Fit `NMF(n_components=50, init="nndsvda", solver="cd", max_iter=500, random_state=seeds[0])` on D only.
  - Freeze H in `cache/nmf_H.npy` with its checksum.
  - R and prediction-only lines get factors by `transform` with the frozen H. They are never refit.
- **Check.** A unit test refits NMF on D and checks that the factors are identical to the stored H, and that the reconstruction error is logged.

## Areas where no change is needed

- **TreeSHAP cost.** About 1,557 boosters per fold × about 204 rows, with native `pred_contrib`, takes under 5 min.
- **Memory.** Apart from the M1 tensor (AC6) and Gram matrices (AC3), memory is not binding: 𝒟 and Y are under 50 MB each.
- **Replication in R.** It is only scoring of the frozen lists and takes minutes.
- **Toy pipeline.** At 600 × 300 × 400 it should run the same code in under 3 min. Keep it as the gate.
