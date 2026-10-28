# Toy answer key (DESIGN v3 §13.1)

Written by the verification agent (not the pipeline author). Code: `src/biomark/toysim.py`; tests:
`tests/test_toy_answers.py`; parameters, bounds and derivation: `configs/toy_answer_key.json`
(`generator`, `pass`, `derivation`). No real data are read: dimensions come from `design_params.json:toy_s13`,
seeds from `seeds.json` (`toy_seed` truth, `toy_null_seed` pure null).

**Naming.** DESIGN v3 §13.1 names `src/biomark/toy.py` / `tests/test_toy_pipeline.py`; the brief for this track
named `toysim.py` / `test_toy_answers.py`. They are the same deliverable; DESIGN or the files should be renamed
before `freeze-design`.

## 1. What the generator produces

`make_toy(seed, scenario)` returns a `Toy` in the pipeline's conventions: 600 lines (`TOYL####`) × 12 lineages
(round-robin over the 8 organ-system groups of `subtypes_s8.lineage_groups`; one lineage has 8 lines → `_rare`),
300 genes (Entrez-like strings `9000001…`, symbols `TG001…`, 10 chromosomes, 4.5 Mb spacing), 400 shared
dictionary features + 849 gene-specific columns (EXPDOSE own, LOSS own/paralog, EXPPAR paralog), feature
metadata (family, gene, d_type, stage G/X, scope, locus group, carriers, NA), `gene_features` (𝒟_g membership),
covariates (lineage, library group with KY-only nested in 3 lineages, CellFormat, SerumFree, medium,
ScreenDoublingTime/DT_unknown, MSI/DMG_BURDEN/TMB, TP53_status, CN_absent, MSI_unknown, PatientID with ~15 shared-
patient pairs, OncotreePrimaryDisease), expression, log2 CN, observed effects (NaN for 8 genes in KY-only lines),
OOF-like predictions **with per-fold offsets** (the pipeline must fold-centre), patient-grouped lineage-stratified
folds, a positive-control table, a toy CGT table, and `answer_key`.

Correlated blocks: AMP/DEL loci of 3 genes (neighbour co-event probability 0.8, or 0.985 for "tight" loci), arm-level
LOSS (whole-arm locus groups), PROG pairs with r = 0.9, PLOIDY~WGD, NRF2 score ~ a LOFBI gene, MSI passengers.
Locus groups are computed with the DESIGN rule (≤ 10 Mb, Jaccard ≥ 0.5, per family) from the realised calls.

**Scenarios.** `truth` plants everything below. `null` is the identical world (same draws) with every molecular,
culture and cutting determinant set to zero; lineage-only genes, the lineage-dummy gene, the self-null gene, library
offsets, latent programmes and fold offsets remain. The pure-null toys use `toy_null_seed + b`.

**Predictors.** `smooth` (default): prediction = μ + 0.9·1.2·tanh((T_visible − μ)/1.2) + 0.6·library offset +
N(0, 0.12) + fold offset, where T_visible is the truth without the predictor-blind plants. `ridge` (DESIGN §13.1):
5-fold CV ridge (per-gene GCV α) on the noisy observations using all dictionary features (minus predictor-blind
features), lineage and library one-hots and the hidden latent programmes, plus fold offsets. Why both: see §4.

## 2. Plants (`answer_key["plants"]`; one row each, with sign, β, proxies and neutral features)

| Kind | n | Truth | What the pipeline must do |
|---|---|---|---|
| pan_lineage | 20 | HOT own ×6, HOT pathway ×3, LOFBI paralog ×4, LOSS paralog (CYCLOPS) ×4, SIL partner ×2 (stage X), AMP partner ×1; β −0.45…−0.9 (LOSS −0.4…−0.6), uniform across lineages | recall, cross-lineage label |
| lineage_restricted | 5 | HOT present pan-lineage, effect only in 1–2 lineages of one organ group | not cross-lineage |
| own_locus_cn | 5 | AMP of the gene's own locus (neighbours co-amplified, one neighbour unexpressed) | recovered with own CN dropped by rule; "amplified-driver addiction" |
| msi_like | 5 | dependency on true MSI; DMG_BURDEN = 5 + 2.5·MSI (r² ≈ 0.6 > 0.49 → dropped too); 12 LOFBI passengers enriched in MSI lines | SIG:MSI recovered, MSI in `dropped_cov`; passengers not selected |
| culture_confounded | 5 | dependency on suspension (3) or serum-free (2); a trap HOT/SIL feature with P(carrier) 0.5 vs 0.03 | trap not in MD-strict (culture is a design block) |
| co_located | 5 | DEL of core → dependency of another gene; passenger neighbour Jaccard 0.54–0.71 with ≥ 6 discordant carriers each side | passenger never named at gene level; true gene named |
| co_located_tight | 3 | as above, Jaccard ≥ 0.96, < 5 discordant carriers | locus level only |
| fold_offset_artefact | 5 | prediction offset −0.8 in one fold; HOT feature with 70% of 30 carriers in that fold; no effect | not selected after centring; selected without |
| genomic_with_expression_proxy | 5 | LOFBI of A → B; EXPPAR:A (or PSCORE:NRF2) strongly correlated | selected in stage G |
| prog_driven | 5 | −0.2…−0.3 per SD of a non-lineage PROG factor (partner r 0.9 accepted) | stage X recovery |
| p53_w | 2 | TP53-mutant lines less dependent (W, +0.5…+0.7) | recovered with TP53 dropped by rule |
| predictor_blind | 3 | HOT → gene in truth/observations, absent from every predictor input | S1-obs finds; EN-pred must not (blindness) |
| wrong_sign_cgt | 1 | truth W (+0.7); toy CGT declares S | counted as a sign error, not a miss |
| six_carrier_strong | 1 | exactly 6 carriers in 6 lineages, β −1.5 | EN-pred selection frequency > 0.5 |
| lineage_dummy_feature | 1 | SIL feature identical to a lineage indicator; gene dependent in that lineage | `collinear_with_design`, never a determinant |
| lineage_only | 10 | strong lineage effects (1–3 lineages, one with a lineage-proxy SIL), nothing else | 0 in MD-strict, class 3 |
| self_null | 1 | gene deleted (DEL, locus group of 3) → carriers' effect ≈ 0 | DEL_g and locus mates outside U; NC3 |
| cutting_artefact | 2 loci | every amplicon gene (core expressed, neighbours unexpressed) depleted −0.5 in carriers | core labelled "cutting artefact" |
| low_tier | 10 | latent-only genes with tier `low` | not eligible |

`answer_key["pairs"]` gives every (gene, feature) a role: `true`, `proxy_ok` (locus-group mates, PROG partner),
`neutral` (expression proxies of a true determinant, DMG_BURDEN for MSI genes, CIN for TP53 genes, LOSS_g of the
self-null gene, cutting-artefact pairs), or a trap role (`culture_trap`, `passenger`, `fold_artefact`, `self_null`,
`lineage_collinear`). **FDP_true** = hits that are neither true/proxy_ok nor neutral, over non-neutral hits.

## 3. Adapter interface the pipeline must expose

`biomark.pipeline.run_pipeline(toy, *, fold_centre=True, null_B=None, arms=("S1-obs","S1-pred","EN-obs","EN-pred"))`
returns a dict with (columns are the minimum; extra columns allowed):

| key | columns |
|---|---|
| meta | dict with `design_params_sha256` |
| eligibility | gene, eligible, r_within |
| arms | arm, gene, feature, rank, rank_stat (EN: selection frequency; S1: p), selected (EN rule / S1 BH F1), beta, stage |
| md_strict | gene, feature, beta, q_sc, dropped_cov |
| collinear_with_design | feature (gene optional) |
| locus_resolution | gene, feature, locus_group, level (gene/locus), n_discordant_carriers |
| confounding_delta | gene, feature, block, flagged |
| amplicon_test | gene, feature, label |
| cross_lineage | gene, feature, label |
| gene_classes | gene, failure_class |
| self_null | gene, feature, frac_dependent |
| null_calibration | list (S1-obs, S1-pred, EN-obs selected, EN-pred selected, MD-strict), fdr_emp_fixed, hits_obs, mean_hits_perm |
| cgt_scores | arm, test_id, hit, sign_correct |
| gate | control_id, recovered |

Part-B tests are skipped while `src/biomark/pipeline.py` lacks `run_pipeline`. Because they are slow, they also run
only with `BIOMARK_TOY_PIPELINE=1`: full arms, null B = 20, plus 20 pure-null toys. Run them before `freeze-bounds`
and after every `src/` fix (DESIGN §9.3 gate iterations). Runs used: smooth truth toy
(default seed), the same toy with `fold_centre=False, null_B=0, arms=("EN-pred",)`, the ridge truth toy, and
`pure_null_B` = 20 pure-null toys with `null_B=0` (override with `BIOMARK_TOY_PURE_NULL_B`).

## 4. Findings while building the key (for the PI / pipeline author)

1. **Dense predictors break the prediction lists' FDP.** With the DESIGN ridge predictor, OOF predictions are a
   linear function of every dictionary feature, so S1-pred flags whatever the predictor uses, including features that
   correlate by chance with hidden latent biology in these 600 lines. Oracle S1-pred FDP_true is ≈ 0.8 and a pure-null
   toy gives hundreds of S1-pred hits. Raising the ridge penalty makes it worse (residual variance shrinks faster than
   the noise weights). The fixed-threshold permutation null cannot see this, because permuting 𝒟 breaks the
   prediction–feature link, so `FDR_emp_fixed ≥ FDP_true − 0.05` cannot hold for S1-pred. An MD-strict analogue
   (S1-pred hits screened on observations, BH within the selection) still has FDP ≈ 0.7. The chance associations are
   real in-sample, so the observed-effect screen passes them, and the lenient within-selection BH threshold lets them
   through. Only R replication can catch these. Consequence for the key: prediction lists are hard-bounded on the
   smooth toy; on the ridge toy the bound is hard only for S1-obs, EN-obs and MD-strict, with S1-pred and EN-pred
   reported. Expect MD-strict to be at risk on the ridge toy. That would be a real result about the double dip, not a
   toy bug.
2. **The ridge toy cannot test most plant classes.** Ridge on toy features cannot learn lineage-restricted
   (interaction) or culture effects. With it, 4/5 fold-artefact, 5/5 lineage-restricted, 3/5 culture, 2/3
   predictor-blind genes and both cutting cores have r_within < 0.10 and drop out of U. So recovery bounds use the
   smooth predictor (the brief's definition), and the ridge toy is used only for error rates. This departs from DESIGN
   §13.1 and needs a D-entry.
3. **r_within with small strata.** DESIGN §2 does not say what happens to lineage × library strata with < 5 lines.
   Pooling them across libraries leaks library offsets into r_within and made 3/10 lineage-only genes eligible in the
   toy. The key (`oracle_r_within`) excludes such lines. The pipeline should do the same, or lineage-only bounds may
   fail.
4. **Covariate-drop rule (v3).** Satisfied by construction. DMG_BURDEN is dropped with SIG:MSI because r² > 0.49;
   MSI passengers keep the MSI covariate (R² ≈ 0.2). Lineage, library and culture are never dropped. Of the 10
   lineage-proxy SIL features, some exceed R²(d ~ design) = 0.49 (`collinear_with_design`) and some do not, so both
   branches are exercised. The lineage-dummy SIL feature has R² = 1.
5. **Pure null is not hit-free by design.** The self-null gene keeps LOSS_g (CYCLOPS rule) and EXPDOSE_g in U.
   Both contain its DEL carriers, whose effect is ≈ 0, so they are real associations in every scenario. Pure-null
   "P(|hits| ≥ 1)" therefore counts only non-neutral hits (`fdp_true(...)["n_hits"]`). Counting raw hits would give 1.0
   for any correct pipeline. Likewise S1 lists are F1 only. Own-EXPDOSE tests form the separate F2 family (≈ 150
   tests, BH within), which alone produces a chance hit in ~5% of pure nulls. Merging F1 and F2 doubles P(any hit).
6. **Fold-artefact plant strengthened.** At an offset of −0.5 the oracle 5th percentile of live artefacts (raw
   predictions, BH) was exactly 3/5. The bound rule would then have lowered DESIGN's ≥ 3/5 positive control to 2/5,
   so the plant was strengthened to −0.8 (live ≥ 4/5 in 8/8 test seeds) rather than weakening the bound.
7. **Permuted NA masks.** CN_absent (5%) and MSI_unknown (5%) lines carry NaN in their families. A permutation that
   moves 𝒟 rows but leaves the indicator columns behind changes which rows are imputed. The pure-null P(any hit) and
   the FDR_emp_fixed check are the toy's sensitivity to that bug.

## 5. Pass bounds and their derivation

Rule (`toysim.propose_bounds`, applied to `derivation.smooth`, 20 truth seeds `toy_seed + i` and 20 pure-null
seeds `toy_null_seed + i`). The oracle is a homoskedastic OLS partial-regression test on the DESIGN v3 design blocks,
with the v3 drop rules and BH. Selection is approximated by S1-pred BH hits; MD-strict is approximated by those hits
screened on observations. The oracle never runs the pipeline.

- fraction "min" bound = min(DESIGN value, floor₀.₀₅(oracle 5th percentile − 0.10));
- count "min" bound over n plants = min(DESIGN value, ⌊n·(p05/n − 0.10)⌋);
- "max" bound = max(DESIGN value, ⌈oracle 95th percentile⌉);
- error-rate bounds = DESIGN values (FDP_true ≤ 0.15, FDR_emp_fixed ≥ FDP_true − 0.05, pure-null P(any hit) ≤ 0.10);
  the oracle only checks that they are attainable.

The 0.10 headroom allows for EN stability selection and CR2 + ν tests being less powerful than the oracle.

<!-- BOUNDS-TABLE -->
| bound | value | DESIGN toy_s13 | oracle (smooth) mean / p05 / p95 | checked on |
|---|---|---|---|---|
| `pan_lineage_recall_min` | 0.8 | 0.8 | `pan_recall_md_strict_oracle` 1.0 / 1.0 / 1.0 | EN-pred selected and MD-strict |
| `pan_lineage_recall_s1_obs_min` | 0.9 | — | `pan_recall_s1_obs` 1.0 / 1.0 / 1.0 | S1-obs BH |
| `own_locus_cn_recall_min` | 0.8 | 0.8 | `own_locus_recall_s1_obs` 1.0 / 1.0 / 1.0 | MD-strict |
| `msi_like_signature_recovered_min` | 4 | 4 | `msi_recovered_s1_obs` 5.0 / 5.0 / 5.0 | MD-strict, + `dropped_cov` contains MSI |
| `msi_like_passenger_lof_selected_max` | 1 | 1 | `msi_passenger_genes_s1_obs` 0.05 / 0.0 / 0.05 | EN-pred selected (MSI genes with a passenger) |
| `co_located_locus_level_min` | 4 | 4 | `co_located_true_gene` 4.95 / 4.95 / 5.0 | reported in `locus_resolution`, passenger not gene-level |
| `co_located_true_gene_level_min` | 4 | — | `co_located_true_gene` 4.95 / 4.95 / 5.0 | true DEL named at gene level |
| `co_located_passenger_gene_level_max` | 1 | — | `co_located_passenger_gene` 0.0 / 0.0 / 0.0 | passenger named at gene level |
| `co_located_tight_gene_level_max` | 0 | — | `tight_gene` 0.0 / 0.0 / 0.0 | any gene-level call in tight loci |
| `culture_confounded_in_md_strict_max` | 1 | 1 | `culture_trap_md_strict_oracle` 0.0 / 0.0 / 0.0 | MD-strict and S1-obs BH |
| `culture_flagged_if_reported_frac_min` | 0.8 | — | — | trap rows in `confounding_delta` (culture block) flagged |
| `lineage_only_in_md_strict_max` | 0 | 0 | `lineage_only_hits_s1_obs` 0.0 / 0.0 / 0.0 | MD-strict rows of lineage-only genes |
| `lineage_only_class3_frac_min` | 0.8 | 0.8 | `lineage_only_ineligible_frac` 0.975 / 0.9 / 1.0 | `gene_classes` == 3 |
| `lineage_dummy_in_md_strict_max` | 0 | — | `lineage_dummy_collinear` 1.0 / 1.0 / 1.0 | MD-strict |
| `lineage_dummy_flagged_collinear` | True | — | `lineage_dummy_collinear` 1.0 / 1.0 / 1.0 | `collinear_with_design` |
| `fold_offset_selected_after_centring_max` | 1 | 1 | `fold_hits_centred` 0.0 / 0.0 / 0.0 | EN-pred selected (centred) |
| `fold_offset_selected_without_centring_min` | 3 | 3 | `fold_live_raw_bh` 4.95 / 4.95 / 5.0 | EN-pred selected with `fold_centre=False` |
| `genomic_stage_G_recovered_min` | 4 | 4 | `proxy_stage_G_s1_pred` 5.0 / 5.0 / 5.0 | EN-pred selected, stage G, primary feature |
| `prog_driven_en_pred_min` | 4 | — | `prog_recovered_s1_pred` 5.0 / 5.0 / 5.0 | EN-pred selected |
| `p53_md_strict_min` | 1 | — | `p53_recovered_s1_obs` 2.0 / 2.0 / 2.0 | MD-strict, W sign |
| `six_carrier_selected_subsample_frac_min` | 0.5 | 0.5 | `six_carrier_s1_obs` 1.0 / 1.0 / 1.0 | EN-pred `rank_stat` (frequency) > |
| `predictor_blind_en_pred_selected_max` | 0 | — | `blind_hits_s1_pred` 0.0 / 0.0 / 0.0 | EN-pred selected (blindness) |
| `predictor_blind_s1_obs_recall_min` | 0.9 | — | `blind_recall_s1_obs` 1.0 / 1.0 / 1.0 | S1-obs BH |
| `wrong_sign_counts_as_sign_error` | True | True | `wrong_sign_hit_w_s1_obs` 1.0 / 1.0 / 1.0 | `cgt_scores`: hit and not sign_correct |
| `cross_lineage_pan_frac_min` | 0.9 | — | `cross_lineage_pan_frac` 1.0 / 1.0 / 1.0 | `cross_lineage` label of labelled pan pairs |
| `cross_lineage_restricted_max` | 2 | — | `cross_lineage_restricted_n` 0.7 / 0.0 / 2.0 | restricted pairs labelled cross-lineage |
| `cutting_core_label_min` | 1 | — | `cutting_core_labelled` 2.0 / 2.0 / 2.0 | `amplicon_test` = cutting artefact |
| `own_locus_addiction_label_min` | 3 | — | `own_locus_addiction` 4.75 / 4.0 / 5.0 | `amplicon_test` = amplified-driver addiction |
| `self_null_in_universe_max` | 0 | — | — | self-null rows in `arms` |
| `nc3_frac_dependent_max` | 0.1 | — | `nc3_frac_dependent` 0.0 / 0.0 / 0.0 | `self_null.frac_dependent` |
| `gate_overall_recall_min` | 0.6 | — | `gate_recall_s1_obs` 1.0 / 1.0 / 1.0 | `gate` on MD-strict, gate-eligible |
| `gate_family_recall_min` | 0.5 | — | `gate_recall_s1_obs` 1.0 / 1.0 / 1.0 | per family with >= 2 rows |
| `planted_eligible_frac_min` | 0.85 | — | `planted_eligible_frac` 0.9974 / 0.9825 / 1.0 | oracle eligibility (self-test) |
| `planted_fdp_true_max` | 0.15 | 0.15 | `fdp_md_strict_oracle` 0.0543 / 0.0196 / 0.0921 | FDP_true per list (smooth: all 5; ridge: S1-obs, EN-obs, MD-strict) |
| `fdr_emp_fixed_minus_fdp_true_min` | -0.05 | -0.05 | — | `null_calibration.fdr_emp_fixed` >= FDP_true + bound (same lists) |
| `pure_null_bh_lists` | ['S1-obs', 'S1-pred', 'MD-strict'] | — | — | lists for the pure-null check |
| `pure_null_bh_lists_p_any_hit_max` | 0.1 | 0.1 | `null_hits_s1_obs` 0.05 / 0.0 / 0.05 | P(non-neutral hits >= 1) over pure-null toys |
| `pure_null_B` | 20 | 20 | — | number of pure-null toys |

Feasibility of the error-rate bounds (oracle, 20 seeds): smooth toy FDP_true p95 = {'fdp_s1_obs': 0.1135, 'fdp_s1_pred': 0.1242, 'fdp_md_strict_oracle': 0.0921}; pure-null
P(non-neutral hit) = {'s1_obs': 0.05, 's1_pred': 0.1, 'md_strict_oracle': 0.05}. Ridge toy FDP_true p05 = {'fdp_s1_obs': 0.0468, 'fdp_s1_pred': 0.8117, 'fdp_md_strict_oracle': 0.6717}; pure-null
P(non-neutral hit) = {'s1_obs': 0.1, 's1_pred': 1.0, 'md_strict_oracle': 1.0}. S1-pred and the MD-strict analogue are therefore structurally unable
to meet the error-rate bounds on the ridge toy (§4.1).

Lists hard-bounded on the smooth toy: ['S1-obs', 'S1-pred', 'EN-obs selected', 'EN-pred selected', 'MD-strict']; on the ridge toy: ['S1-obs', 'EN-obs selected', 'MD-strict']
(reported only: ['S1-pred', 'EN-pred selected']).

Reproduce: `PYTHONPATH=src python -m biomark.toysim derive --n-seeds 20 --write` (≈ 12 min; rewrites `derivation` and
`pass`). Any change after `freeze-design` needs a D-entry (DESIGN §13.4).
