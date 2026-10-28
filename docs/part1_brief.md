# Part 1 technical brief (for Part 2 design and code)

Source: `C:\Users\ritab\Desktop\project` (abbreviated `P1/` below). It is **read-only**. Everything here was checked
against the code and caches on 2026-10-28. Python 3.13.1, CPU only. `P1/requirements.txt` pins numpy 2.1.3,
scipy 1.15.2, pandas 2.2.3, pyarrow 19.0.1, scikit-learn 1.5.2, lightgbm 4.7.0 and statsmodels 0.14.4.

## 1. Data

**IDs.**
- Lines are DepMap `ModelID` (`ACH-xxxxxx`).
- Genes are **Entrez ID strings** (`"7105"`) in every cached matrix. Raw DepMap headers `"SYMBOL (Entrez)"` are
  parsed by `essml.io.parse_gene_column`.
- Symbols live in `P1/cache/gene_map.parquet` (columns `entrez, symbol, in_matrices`; 21,269 rows).
- Omics use the row with `IsDefaultEntryForModel == "Yes"`. Values are float32.

**Raw DepMap files** (`P1/data/`, 26Q3, checksums in `P1/configs/depmap_local_sha256.txt`) and their caches:

| Name | File | Cache (`P1/cache/`) | Lines × genes |
|---|---|---|---|
| crispr (target) | CRISPRGeneEffect.csv | `crispr.parquet` (guarded) | 1361 × 18434 |
| expression | OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv | `expression.parquet` | 2130 × 19218 |
| cn | OmicsCNGeneWGS.csv (relative CN, linear) | `cn.parquet` | 1763 × 18778 |
| damaging | OmicsSomaticMutationsMatrixDamaging.csv (0/1/2) | `damaging.parquet` | 2395 × 19657 |
| hotspot | OmicsSomaticMutationsMatrixHotspot.csv | `hotspot.parquet` | 2395 × 165 |
| fusion | OmicsFusionFiltered.csv (long table) | `fusion.parquet` | 94,718 rows |

Also present: Model.csv, ModelCondition.csv, and `CRISPRScreenMap.csv` / `ScreenSequenceMap.csv` (used post hoc
only). **Not local:** OmicsGlobalSignatures (MSI, ploidy, WGD, CIN), SubtypeMatrix/SubtypeTree, the full mutation
MAF and OmicsProfiles. If Part 2 needs them for clinical biomarkers they must be downloaded manually; the portal
blocks scripts.

**Raw feature blocks.** `essml.features.build_all` builds them, `scripts/build_features.py` runs it, and the output
is `P1/cache/features/<block>.parquet`, indexed by ModelID with string columns.

| Block | Content | Rows |
|---|---|---|
| E | log2(TPM+1), as supplied | 2130 |
| C | clip(log2((CN+0.1)/1.1), ±3); NaN → 0 (diploid); all-NaN genes dropped | 1763 |
| C_instab | fraction of autosomal genes with \|CN−1\| > 0.3 | 1763 |
| M | damaging dosage 0/1/2 | 2395 |
| M_burden | log1p(number of damaged genes) | 2395 |
| H | hotspot > 0 | 2395 (165 genes) |
| F_pair / F_gene | fusion pair `"a--b"` (Entrez) / partner indicators, present in ≥ 5 lines | 2130 |
| meta | OncotreeLineage, OncotreePrimaryDisease, Sex, Age | 2554 |
| meta_prov | GrowthPattern, ModelType, PrimaryOrMetastasis (**ablation only**: library proxy) | 2554 |
| composite | p53/Rb flags (TP53_mut, MDM2_amp, CDKN2A_deepdel, RB_altered, …) | 2395 |
| P | Nusinow 2020 proteomics | 375 |

**Fold-fitted transforms.** `features.FoldTransformer(specs).fit(raw, train_ids).transform(raw, ids)` produces
columns `"<block>:<feature>"`:
- E drops columns with training SD < 0.1, leaving about 18,071;
- count blocks keep columns present in ≥ 5 training lines;
- every block is z-scored on training lines, then clipped to ±8;
- a line absent from a genomic block gets 0 plus an `:absent` flag; absent expression gets the training mean.

**Coverage.**
- 1,292 screened lines have expression; 69 screened lines lack it and are outside the population.
- **838 lines have expression but no CRISPR screen**: 710 cell lines and 128 organoids, of which 85 are
  non-cancerous. They are candidates for prediction-only subtyping.
- 55 dev and 16 test lines have no CN profile, so their CN is treated as diploid. All lines have mutation data.

**Model.csv** (2554 × 48), counts over the 1,292 population lines:
- lineage and subtype: `OncotreeLineage` (29), `OncotreePrimaryDisease` (86), `OncotreeSubtype` (175), `Sex`,
  `Age`, `PrimaryOrMetastasis`;
- culture: `GrowthPattern` (Adherent 910, Suspension 167, Dome 108, Mixed 63, Spheroid 40), `ModelType` (Organoid
  108), `OnboardedMedia` (146 levels), `SerumFreeMedia`;
- `ModelSubtypeFeatures` (221 lines, e.g. "HER2+", "ER+", "EWS::FLI"): candidate clinical labels.

**ModelCondition.csv** (3385 × 17; about 1.5 conditions per model) has `CellFormat`, `GrowthMedia`,
`SerumFreeMedia`, `PlateCoating`, `DataSource` and `AnchorDrug`. Part 1 **never used it**, and DESIGN excluded
media and DataSource as screen-status proxies. To get the screened condition, join
`ScreenSequenceMap.ModelConditionID`; that file also gives `Library`, `Days` and `ScreenType`.

**Library groups.**
- `P1/configs/library_groups.csv` has columns `ModelID, lib_group, n_nan_genes` and covers **dev lines only**.
- `lib_group` is the rank of a line's gene-effect NaN pattern; 4 = other. Group sizes: 0 = Avana-only (621),
  1 = Avana+KY (154), 2 = Humagne-CD (143), 3 = KY (47), 4 = mixed (57). Adjusted Rand index against the true
  library map is 0.949 (`docs/library_map_check.md`).
- Per-line groups for test lines were not saved. Derive them for any line from CRISPRScreenMap →
  ScreenSequenceMap.Library; this reads no effect values.
- Library also tracks screen length: median 21 days for Avana, 24 for Humagne-CD, 14 for KY.

**Patient grouping.** `PatientID` is in `P1/configs/splits.csv`. 104 screened lines share a patient. Every split
and inner CV is a StratifiedGroupKFold grouped by patient and stratified by lineage.

## 2. Splits and the guard

**`P1/configs/splits.csv`** has columns `ModelID, PatientID, OncotreeLineage, role, cv_fold` and 1,292 rows:
- role dev 1022, `cv_fold` 0–4 with sizes 208/202/201/202/209;
- role test 255, `cv_fold` −1;
- role noncancerous 15, `cv_fold` −1 (set A).

The split seed is 20261026 (`configs/seeds.json`) and the generator is `scripts/make_splits.py`. Lineages with
fewer than 10 lines are pooled as `_rare` for stratification only.

**Guard** (`src/essml/guard.py`):
- `io.load_matrix('crispr')` raises without `unguarded=True`.
- `guard.load_targets(ids)` returns dev rows freely and writes no log.
- Any test or set-A ID needs `ESSML_UNBLIND=1`, a committed `configs/finalists.json`, a clean tree, a listed config
  hash, and **no earlier granted access for that (hash, scope)**.
- Every non-dev request is appended to `P1/results/test_access_log.csv`, which is a write into P1.

**Test status.** Yes, the test CRISPR data were opened once in Part 1, once per finalist (log rows 917–923,
2026-10-28T00:54–01:26). Those were: m5_avg, m1_E_resB0, t_E_status, t_E_status_noEC, b1, b2 and b0. Every
finalist's (hash, scope="test") is now used up, so the guard **will refuse test targets to everyone**. Part 2
would have to read `crispr.parquet` directly. That is a protocol decision: declare it, and lock the Part 2 design
before reading.

**Cached test predictions** (255 × 18434, models refit on all 1022 dev lines) are in `P1/cache/test/preds/`:
- `{m5_avg, m1_E_resB0, t_E_status, t_E_status_noEC, b0, b1, b2}.parquet` and `t_E_status_fold{0-4}.parquet`;
- `*_with_confidence.parquet` in long format.

The manifest is `P1/results/test/predictions_manifest.json` and test metrics are in `P1/results/test/<model>/`.

## 3. Models

All models are in `P1/src/essml/models/`. Registry keys are in `scripts/run.py:MODELS`.

| Model (run name) | Class | What it is |
|---|---|---|
| **M1** `m1_E_resB0` | `kernel_ridge.KernelRidgeM1(modalities=["E"], residual="B0")` | Multi-output ridge with a **linear** kernel (K = Z Zᵀ; no RBF in M1/M2/M5) on standardised E. Target = y − training gene mean (B0), so **lineage is not used**. One `eigh` per NaN-pattern gene cluster; leave-patient-out λ from 1e2–1e7, shared within variable vs other genes. Genes in patterns with < 50 rows get the gene mean. |
| **M2** `m2n_E_paralogs` | `m2_local.LocalStageM2(priors="paralogs", local_genes="variable", nested_lambda=True)` | Per **training-fold variable gene** (SD > 0.25), a ridge on 13 features: `m1_oof` (inner 5-fold OOF M1) plus own expression, not-expressed, hinge, CN, deep deletion and mutation, and for the top-3 Ensembl paralogs their expression, minimum, deep deletion and mutation. Other genes pass M1 through. Stored as `beta_` (G × 13) and `icpt_`. |
| **M3** `g3_m3n` | `lgbm.LGBMResidualM3(pool_blocks=[E,C,M,H,meta], k_screen=50, nested_es=True)` | **One LightGBM booster per training-fold variable gene** (about 2,186 per fold; about 1,557 non-empty) on y − M1_oof. Features: the own M1 prediction plus the top 50 pool columns by \|r\| with the residual. The pool is E, C, M, H and meta one-hot (lineage, disease, sex, age). `LGB_DEFAULT`: lr 0.05, 15 leaves, min_data_in_leaf 20, lambda_l2 1, max_bin 63; ≤ 500 rounds, early stopping 30 on one patient-grouped 80/20 split. Stored as `boosters_` (strings), `model_genes_` and `sel_full_` (indices into `cf_.columns_`). |
| **M5 headline** `m5_avg` | `m5_stack.StackM5(mode="average")` | One-third weights on M1, M2 and M3, each refitted on the training fold. |
| **T** `t_E_status` | `KernelRidgeT(modalities=["E","status"])` | M1 plus `TP53_mut` and `RB_altered`. Fold models are saved **with `Ztr_`**. |

**Attribution feasibility:**
- **M1 and T: exact and linear.** For pattern cluster c,
  `pred_g(x) = mu_g + rbar_cg + (z_x − zbar_c)ᵀ w_cg`, where `w_c = Z_cᵀ alpha_c` (features × genes in c) and
  `zbar_c` is the mean of the training rows of c.
  - Per-feature contribution: `(z_xf − zbar_cf) · w_cfg`.
  - `Z` = `model._train_design()`, i.e. `FoldTransformer` output, clipped ±8, times `col_w_` (all 1 for E-only).
  - Rebuilding Z needs `cache/features`. Each cluster holds `rows`, `genes`, `alpha` (float32) and `kbar`.
  - The largest cluster is about 814 rows × 17,087 genes. Full W is about 18k × 18k float32, roughly 1.3 GB per
    fold, so compute it per gene subset (variable genes: about 150 MB).
  - Ridge weights on correlated expression are dense. Aggregate them over gene sets or modules.
- **M2: exact.** Linear in 13 features, and the `m1` term decomposes through M1.
- **M3: feasible.** `lgb.Booster(model_str=s).predict(F, pred_contrib=True)` gives TreeSHAP over the 51 features
  plus a bias term.
  - Build `F` with `LGBMResidualM3._feats(True, cf_.transform(ids), rows, m1_pred[:, g], sel_full_[j])`.
  - The `m1` feature's SHAP value can be pushed down through M1's linear weights; this is an approximation.
- **M5:** one third of each component's attribution; the result is additive.

**Fitted model objects** are stored at `P1/cache/models/<name>_<hash8>/fold{k}.pkl` plus `meta.json`, as pickled
dicts with keys `cls, state, dropped, fold, train_ids, features, run_id` (`essml.modelstore`).

| Directory | Contents | Usable for attribution? |
|---|---|---|
| `t_E_status_24b46e2a` (5 folds, 60 MB each) and `t_E_status_noEC_fa0181c1` | Complete, `Ztr_` included | Yes, **no refit needed** |
| `m5_avg_6dbebcf4` (5 × 0.4 MB) | Weights and metadata only; `models_ = None` | No: prediction would refit all three inputs |
| `b2_289f5170` | kNN baseline | Not relevant |
| `dry_*` | Dry-run models; folds 1–4 trained on about 612 lines (fold 0 excluded) | No: not the CV models |

There are **no saved CV fold models for m1_E_resB0, M2 or M3**, and no saved all-dev refit of m5_avg (only its
test predictions).

**Refit cost** (measured; 12-core, 16 GB machine):

| Model | Cost | Notes |
|---|---|---|
| M1 | 50–80 s per fold | 12 threads |
| M2 | about 3 min per fold | With the inner-OOF cache; uncached all-dev refit 1384 s |
| M3 | 3–10 min per fold | 6 processes |
| m5_avg, all 1022 dev lines | 27–29 min | 4 threads, fresh caches, peak RSS about 3.3 GB (`P1/results/timing/refit_full_dev.json`) |
| Full headline CV | about 1–2 h | 5 folds |

## 4. Out-of-fold (dev CV) predictions and per-gene metrics

**Per run outputs:**
- `P1/cache/runs/<run_id>/pred.parquet`: 1022 × 18434 float32, index `ModelID`, **no NaN**. Each line is predicted
  by the one fold model that did not train on it.
- Same folder: `per_gene_fold.parquet` (`fold, gene, n_obs, r, rho, ap, nap, auroc, prev, r_B0, r_B1, is_<class>`)
  and `per_line.parquet`.
- `P1/results/runs/<run_id>/`: `summary.json`, `per_fold.csv` and `per_gene.csv.gz` (`gene, r_mean, r_sd, nap_mean,
  …`). The registry is `P1/results/registry.csv`.

**Key run IDs:**

| Model | run_id | dev P1 (variable genes) |
|---|---|---|
| **m5_avg (headline)** | `20261027-161041_m5_avg_6dbebcf4` | 0.3175 |
| m1_E_resB0 | `20261026-222117_m1_E_resB0_91e52f9d` | 0.3081 |
| m2n_E_paralogs | `20261027-044450_m2n_E_paralogs_62ce1575` | 0.3115 |
| g3_m3n | `20261027-014824_g3_m3n_601e3f5e` | 0.310 |
| t_E_status | `20261026-224125_t_E_status_24b46e2a` | 0.3082 |
| m1bt_full (E+genomics) | `20261027-002223_m1bt_full_48508032` | 0.3096 |
| b1 / b2 | `20261026-220921_b1_0e2a7a27` / `20261026-221456_b2_289f5170` | 0.184 / 0.253 |

Also available: the permuted-expression null `ctl_resB0_permuted` (`20261026-231349_…`), and variants run on
leave-one-lineage-out (`*_lolo`) and with the EC lines excluded (`*_noEC`).

```python
import pandas as pd
P = r"C:\Users\ritab\Desktop\project"
oof = pd.read_parquet(fr"{P}\cache\runs\20261027-161041_m5_avg_6dbebcf4\pred.parquet")  # 1022 x 18434
sp = pd.read_csv(fr"{P}\configs\splits.csv", index_col=0)
fold = sp.loc[oof.index, "cv_fold"]
```

**Per-gene tables:**
- `configs/predictable_genes.csv` (primary, m5_avg) and `predictable_genes_m1.csv`: `entrez, symbol,
  r_fold0..4, mean_r, min_r, sd_r, n_folds_evaluable, confidence, tier, …`.
- Tiers: high = mean_r ≥ 0.30 and min_r ≥ 0.20; medium = mean_r ≥ 0.15 and min_r > 0; confidence = mean_r − sd_r.
  Rules are in `configs/predictability.json`.
- m5_avg: 1,632 high, 6,198 medium, 10,085 low, 519 not evaluable (variable genes 1,054 / 921 / 143 / 41). The
  permuted null gives 0 high genes.
- `configs/gene_classes.csv` (dev-wide booleans: variable 2159, selective 2610, rare_selective, common_essential
  848, nonessential; plus `sd, frac_below`) and `configs/eval_genes.csv` (17,915 genes).

## 5. Existing biomarker work

**D24 mutation-driven dependencies** (`docs/mutation_dependencies.md`, `src/essml/mutdeps.py`,
`scripts/mutation_dependencies.py`):
- Pairs: `configs/oncogene_pairs.csv` (19 pairs, 17 eligible with ≥ 5 positive dev lines; columns `pair_id, group,
  dependency_gene/entrez, biomarker_gene/entrez, biomarker_type, threshold, n_mutant_dev, n_wt_dev,
  n_unknown_dev, eligible`). Operational choices: `configs/mutation_dependencies.json`.
- Results: `results/mutation_dependencies.json` and `_folds.csv`.
- Status thresholds: hotspot > 0, damaging > 0, amp CN > 2, deep deletion CN < 0.25, TP53-WT.
- A built-in "biomarker-only" reference uses the training-fold mean effect by status.

Findings (AUPRC for effect < −0.5):
- **E-only misses RAS/RAF hotspot addiction.** NRAS: E-only 0.21, biomarker-only 0.60, m5_avg 0.63, M3 0.66.
  BRAF goes 0.56 → 0.66 and KRAS 0.76 → 0.81. NRAS within-lineage r rises from 0.19 to 0.73 with M3.
- Expression already captures amplification, p53-axis and paralog pairs: genomics changes AUPRC by −0.02 to +0.05.
- Biomarker-only is worse than E-only on 13 of 15 pairs. Pooled ΔAUPRC: m5_avg +0.068 [0.051, 0.086], M3 +0.062,
  E+genomics +0.029.
- The pre-registered "genomics clearly helps" rule was not met; only NRAS qualifies.

**D22 literature pairs.**
- `configs/biomarker_pairs.csv` has 43 pairs: oncogene addiction, paralog loss, MSI→WRN and lineage TFs.
- Results: `results/predictability/biomarker_pairs.csv`. m5 tiers are 33 high, 7 medium and 3 low (ALK, ME3, AR,
  all rare).

## 6. Reusable code (read-only import)

```python
import sys; sys.dont_write_bytecode = True        # or PYTHONDONTWRITEBYTECODE=1 / python -B
sys.path.insert(0, r"C:\Users\ritab\Desktop\project\src")
```

Bytecode matters: importing would write `__pycache__/*.pyc` into `P1/src/essml/` whenever a .pyc is stale or
missing.

**Safe to use as is:** `metrics` (`per_gene_pearson`, `per_gene_ap_norm`, …), `evaluate.evaluate`,
`genesets.gene_classes`, `baselines` (GeneMean, LineageMean, KNNExpression), `heldout` (patient bootstrap),
`mutdeps` (`biomarker_status`, `fold_metrics`, `biomarker_only`, `paired_ci`), `predictability`
(`fold_stats`, `assign_tiers`, `load_table`), `priors.LocalFeatureBuilder.from_cache()` (reads only),
`features.load_raw`, `FoldTransformer`, and `kernel_ridge.KernelRidgeM1` with `fit` and `predict`.

**Hard-coded paths that WRITE into P1. Redirect or avoid them:**

| Path | Writer | When | Fix |
|---|---|---|---|
| `cache/m2_inner/*.npy` | `m2_local.INNER_CACHE` | M2/M3 fit on a cache miss | Set `m2_local.INNER_CACHE = <project_b>/cache/m2_inner` (module global, read at call time). Optionally copy P1's 3.8 GB for exact-key hits on the same dev folds. |
| `cache/m5_inner`, `cache/m5_outer` | `m5_stack.INNER_CACHE/OUTER_CACHE` | M5 fit/predict | Redirect, or pass `cache=False` (as `evaluate_test.py` does) |
| `cache/models/<run>/` | `modelstore.save_fold_model` (MODEL_DIR) | When saving | Pass `model_dir=` |
| `results/test_access_log.csv` | `guard._log` | Any non-dev target/count request | Request dev IDs only |
| `cache/*.parquet`, `cache/features/` | `io.build_cache`, `features.save_raw` | Only if missing, or if called | Do not call |
| `results/registry.csv`, `cache/runs/` | `scripts/run.py` and the other scripts | Always | **Never run P1 scripts.** Reimplement the CV loop. |

**Other gotchas:**
- `features.check_fit_ids` raises if any training ID has role test or noncancerous. Fitting on dev+test needs a
  monkeypatch of `features.check_fit_ids`.
- Models read `configs/splits.csv` for PatientID groups when `groups=None`.
- M5 fold pickles restore `Y_` through `guard.load_targets(train_ids)`, which works for dev IDs.

## 7. Pitfalls documented in Part 1

1. **Chronos is CN-corrected and post-processed.** It applies CN correction, proximity-bias (chromosome-arm)
   alignment and RCPC removal of PCs correlated with common-essential dropout. CN→effect associations are partly
   removed by construction.
2. **Library confound.**
   - Library is confounded with lineage and screen length.
   - P1 is 0.3175 overall, 0.289 within library and 0.211 within lineage × library.
   - GrowthPattern/ModelType ("Dome" ≡ organoid) proxy library; they added +0.014 to P1 and were excluded.
   - Check Part 2 subtypes against library and organoid status.
3. **Missing genomics = diploid wild type.** Missing CN (55 dev lines) is filled with 0 plus an `:absent` flag.
   `mutdeps` instead treats it as unknown; use that convention for biomarker strata.
4. **Pooled OOF artefact (D3).** Fold models differ in offsets, so pooling OOF predictions gives spurious r. Centre
   per gene within fold before clustering, or test whether clusters associate with fold.
5. **Joint Chronos fit.** Test-line effects come from the same Chronos run as the dev lines.
6. **Expression-only blind spot.** Only M3 recovers NRAS/KRAS/BRAF/PIK3CA hotspot dependencies.
7. **Low tiers are noise-limited.** Restrict subtype features to high/medium-tier variable genes.
8. **No probabilities.** The models are point regressors.
9. **Permutation nulls missed their bands** (D20b). Simulate control bounds first (P1 LESSONS #10).
