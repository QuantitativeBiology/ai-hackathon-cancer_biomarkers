"""Toy answer key (DESIGN v3 §13.1), written by the verification agent, not the pipeline author.

Part A (runs now): generator self-tests -- shapes and conventions, determinism, planted signal recoverable by a
simple oracle, pure null without signal, fold offsets present and removed by centring, traps live.
Part B (skipped until ``biomark.pipeline.run_pipeline`` exists): the pipeline must meet the pre-declared bounds in
configs/toy_answer_key.json:pass (derivation: docs/toy_answer_key.md).

Env: BIOMARK_TOY_PURE_NULL_B overrides the number of pure-null toys in part B (default toy_s13.pass.pure_null_B).
"""
import os
import re

import numpy as np
import pandas as pd
import pytest

from biomark import paths as P
from biomark import toysim as ts

DP = P.load_config("design_params")
CFG = ts.load_toy_config()
BOUNDS = CFG["pass"]
FAM = CFG["families"]
NF = CFG["named_features"]


# =========================================================================================================
# Part A: generator self-tests
# =========================================================================================================


@pytest.fixture(scope="module")
def toy():
    return ts.make_toy()


@pytest.fixture(scope="module")
def null_world(toy):
    """Same seed, determinants removed (identical world)."""
    return ts.make_toy(seed=toy.seed, scenario="null")


@pytest.fixture(scope="module")
def pure_null():
    return ts.make_toy(scenario="null")


@pytest.fixture(scope="module")
def s1_obs(toy):
    return ts.oracle_s1(toy, "obs", "full")


@pytest.fixture(scope="module")
def s1_pred(toy):
    return ts.oracle_s1(toy, "pred", "full")


def _hits(df):
    return df[df["hit"]]


def test_config_dimensions_and_seeds():
    t13 = DP["toy_s13"]
    seeds = P.seeds()
    assert t13["seed_key"] in seeds and t13["null_seed_key"] in seeds
    assert CFG["generator"]["predictor"]["default"] in ("smooth", "ridge")
    assert set(ts.REQUIRED_OUTPUTS) >= {"arms", "md_strict", "locus_resolution", "confounding_delta", "null_calibration"}
    assert BOUNDS, "pass bounds not yet derived (python -m biomark.toysim derive --n-seeds 20 --write)"
    for k, v in BOUNDS.items():
        if not k.startswith("_"):
            assert isinstance(v, (int, float, bool, list)), k


def test_shapes_and_conventions(toy):
    t13 = DP["toy_s13"]
    n, G = t13["n_lines"], t13["n_genes"]
    assert toy.features.shape[0] == n and toy.effects.shape == (n, G) and toy.predictions.shape == (n, G)
    assert toy.truth.shape == (n, G) and toy.expression.shape == (n, G) and toy.cn_log2.shape == (n, G)
    assert list(toy.effects.columns) == list(toy.genes.index) == list(toy.predictions.columns)
    assert toy.covariates["lineage"].nunique() == t13["n_lineages"]
    assert all(re.fullmatch(r"TOYL\d{4}", i) for i in toy.lines)
    fm = toy.feature_meta
    assert list(fm.index) == list(toy.features.columns)
    assert (fm["scope"] == "shared").sum() == t13["n_features"]
    assert all(re.fullmatch(r"[A-Z]+:[A-Za-z0-9_]+", f) for f in fm.index)
    assert set(fm["family"]) <= set(FAM.values())
    assert set(fm["stage"]) == {"G", "X"}
    stages = DP["surrogate_s4"]["stages"]
    assert all(sf in stages[st] for sf, st in zip(fm["stage_family"], fm["stage"]))
    assert toy.gene_features["feature"].isin(fm.index).all()
    assert toy.gene_features["gene"].isin(toy.genes.index).all()
    need = {"lineage", "library_group", "CellFormat", "SerumFree", "medium_family", "ScreenDoublingTime", "DT_unknown",
            "MSI", "DMG_BURDEN", "TMB", "TP53_status", "CN_absent", "MSI_unknown", "PatientID",
            "OncotreePrimaryDisease"}
    assert need <= set(toy.covariates.columns)
    assert not toy.predictions.isna().any().any()
    assert toy.folds.isin(range(t13["n_folds"])).all() and toy.folds.nunique() == t13["n_folds"]
    # NaN conventions
    cn_feats = fm.index[fm["family"].isin([FAM["amp"], FAM["del"], FAM["loss"]])]
    absent = toy.covariates["CN_absent"].to_numpy()
    assert toy.features.loc[absent, cn_feats].isna().all().all()
    assert not toy.features.loc[~absent, cn_feats].isna().any().any()
    assert toy.cn_log2.loc[absent].isna().all().all()
    msi = toy.features[f"{FAM['sig']}:{NF['msi']}"]
    assert (msi.isna() == toy.covariates["MSI_unknown"]).all()
    nan_lines = toy.effects.isna().any(axis=1)
    assert (toy.covariates.loc[nan_lines, "library_group"] == CFG["generator"]["noise"]["library_nan_group"]).all()
    assert set(toy.covariates["library_group"]) <= set(CFG["generator"]["library_groups"])


def test_patients_within_one_fold_and_lineage(toy):
    cov = toy.covariates.assign(fold=toy.folds)
    g = cov.groupby("PatientID")
    assert (g["fold"].nunique() == 1).all() and (g["lineage"].nunique() == 1).all()
    assert (g.size() > 1).sum() >= 5


def test_library_partly_confounded_with_lineage(toy):
    tab = pd.crosstab(toy.covariates["lineage"], toy.covariates["library_group"]).to_numpy(float)
    n = tab.sum()
    exp = tab.sum(1, keepdims=True) * tab.sum(0, keepdims=True) / n
    chi2 = ((tab - exp) ** 2 / exp).sum()
    v = np.sqrt(chi2 / (n * (min(tab.shape) - 1)))
    assert 0.2 < v < 0.9, v
    g = CFG["generator"]
    nested_share = (toy.covariates["library_group"] == g["library_nested_group"]).groupby(toy.covariates["lineage"]).mean()
    assert (nested_share.sort_values().iloc[-3:] > 0.7).all() and nested_share.median() < 0.3


def test_determinism_and_seed_dependence(toy):
    again = ts.make_toy(seed=toy.seed)
    pd.testing.assert_frame_equal(again.features, toy.features)
    pd.testing.assert_frame_equal(again.effects, toy.effects)
    pd.testing.assert_frame_equal(again.predictions, toy.predictions)
    other = ts.make_toy(seed=toy.seed + 1)
    assert not np.allclose(other.effects.fillna(0).to_numpy(), toy.effects.fillna(0).to_numpy())


def test_null_world_differs_only_in_planted_effects(toy, null_world):
    pd.testing.assert_frame_equal(null_world.features, toy.features)
    pd.testing.assert_frame_equal(null_world.covariates, toy.covariates)
    planted = set(toy.plants()["gene"])
    d = (toy.truth - null_world.truth).abs().max()
    assert set(d.index[d > 1e-12]) <= planted
    assert d[[g for g in toy.genes.index if g not in planted]].max() == 0


def test_locus_groups_and_colocation_plants(toy):
    co = toy.answer_key["colocated"]
    md = DP["dictionary_s3"]["locus_resolution"]["min_discordant_carriers"]
    jmin = DP["dictionary_s3"]["locus_groups"]["binary_jaccard_min"]
    assert co["same_group"].all()
    assert (co["jaccard_true_passenger"] >= jmin).all()
    res = co[co["kind"] == "co_located"]
    tight = co[co["kind"] == "co_located_tight"]
    assert (res[["n_true_only_vs_left", "n_true_only_vs_right", "n_passenger_only"]] >= md).all().all()
    assert (tight[["n_true_only_vs_left", "n_true_only_vs_right"]].max(axis=1) < md).all()
    # self-null set of the self-null gene includes its DEL locus mates
    sn = toy.answer_key["self_null"]
    g = sn["gene"].iloc[0]
    assert f"{FAM['del']}:{g}" in set(sn["feature"]) and len(sn) >= 2


def test_oracle_recovers_planted_truth(toy, s1_obs, s1_pred):
    """The simple oracle meets every recovery bound that has an oracle analogue (bounds are achievable)."""
    ho, hp = _hits(s1_obs), _hits(s1_pred)
    assert ts.recovered(ho, toy, "pan_lineage").mean() >= BOUNDS["pan_lineage_recall_s1_obs_min"]
    assert ts.recovered(hp, toy, "pan_lineage").mean() >= BOUNDS["pan_lineage_recall_min"]
    assert ts.recovered(ho, toy, "own_locus_cn").mean() >= BOUNDS["own_locus_cn_recall_min"]
    assert ts.recovered(ho, toy, "msi_like", accept_proxies=False).sum() >= BOUNDS["msi_like_signature_recovered_min"]
    assert ts.recovered(hp, toy, "genomic_with_expression_proxy", stage="G").sum() >= BOUNDS["genomic_stage_G_recovered_min"]
    assert ts.recovered(hp, toy, "prog_driven").sum() >= BOUNDS["prog_driven_en_pred_min"]
    assert ts.recovered(ho, toy, "p53_w").sum() >= BOUNDS["p53_md_strict_min"]
    assert ts.recovered(ho, toy, "six_carrier_strong").all()
    assert ts.recovered(ho, toy, "predictor_blind").mean() >= BOUNDS["predictor_blind_s1_obs_recall_min"]
    assert ts.recovered(hp, toy, "predictor_blind").sum() <= BOUNDS["predictor_blind_en_pred_selected_max"]
    md = ts.oracle_md_strict(toy, s1_pred, s1_obs)
    for lst in (ho, hp, md):
        assert ts.fdp_true(lst, toy)["fdp"] <= BOUNDS["planted_fdp_true_max"]
    # wrong-sign CGT pair is a hit with the W sign (a sign error against the declared S), not a miss
    ws = toy.plants("wrong_sign_cgt").iloc[0]
    row = ho[(ho["gene"] == ws["gene"]) & (ho["feature"] == ws["feature"])]
    assert len(row) == 1 and row["beta"].iloc[0] > 0
    assert toy.cgt.set_index("feature").loc[ws["feature"], "sign"] == "S"


def test_eligibility_oracle(toy):
    el = ts.oracle_eligible(toy)
    det = toy.plants().query("kind not in ['lineage_only', 'low_tier', 'self_null', 'lineage_dummy_feature', "
                             "'cutting_artefact']")
    assert el.reindex(det["gene"]).mean() >= BOUNDS["planted_eligible_frac_min"]
    lo = toy.plants("lineage_only")["gene"]
    assert (~el.reindex(lo)).mean() >= BOUNDS["lineage_only_class3_frac_min"]
    assert not el.reindex(toy.plants("low_tier")["gene"]).any()


def test_pure_null_has_no_signal(pure_null):
    so = ts.oracle_s1(pure_null, "obs", "full")
    sp = ts.oracle_s1(pure_null, "pred", "full")
    f1 = so[(so["test_family"] == "F1") & ~so["collinear"]]
    assert 0.035 <= (f1["p"] < 0.05).mean() <= 0.065
    n_obs = ts.fdp_true(so[so["hit"]], pure_null)["n_hits"]
    n_pred = ts.fdp_true(sp[sp["hit"]], pure_null)["n_hits"]
    n_md = ts.fdp_true(ts.oracle_md_strict(pure_null, sp, so), pure_null)["n_hits"]
    assert n_obs <= 1 and n_pred <= 1 and n_md <= 1, (n_obs, n_pred, n_md)
    # the only real associations left in a pure null are neutral (self-null gene's own LOSS / EXPDOSE)
    assert ts.fdp_true(so[so["hit"]], pure_null)["n_neutral"] >= 1


def test_fold_offsets_present_and_removed_by_centring(toy):
    fa = toy.answer_key["fold_artefact"]
    for _, r in fa.iterrows():
        means = toy.predictions[r["gene"]].groupby(toy.folds).mean()
        assert means[r["fold"]] - means.drop(r["fold"]).mean() < r["offset"] / 2
    fc = ts.fold_centred(toy)
    assert fc.groupby(toy.folds).mean().abs().max().max() < 1e-12
    raw = ts.oracle_s1(toy, "pred_raw", "full")
    cen = ts.oracle_s1(toy, "pred", "full")
    key = set(zip(fa["gene"], fa["feature"]))
    live = sum((g, f) in key for g, f in zip(raw.loc[raw["hit"], "gene"], raw.loc[raw["hit"], "feature"]))
    gone = sum((g, f) in key for g, f in zip(cen.loc[cen["hit"], "gene"], cen.loc[cen["hit"], "feature"]))
    assert live >= BOUNDS["fold_offset_selected_without_centring_min"]
    assert gone <= BOUNDS["fold_offset_selected_after_centring_max"]
    # the fold-artefact features carry no observed-effect signal
    obs = ts.oracle_s1(toy, "obs", "full")
    assert sum((g, f) in key for g, f in zip(obs.loc[obs["hit"], "gene"], obs.loc[obs["hit"], "feature"])) == 0


def test_culture_trap_live_without_culture_and_removed_with_it(toy, s1_obs):
    traps = set(map(tuple, toy.answer_key["culture"][["gene", "feature"]].to_numpy()))
    ll = ts.oracle_s1(toy, "pred", "lineage_library")
    live = sum((g, f) in traps for g, f in zip(ll.loc[ll["hit"], "gene"], ll.loc[ll["hit"], "feature"]))
    assert live >= 3, "culture traps must be live when culture is not adjusted"
    adj = sum((g, f) in traps for g, f in zip(s1_obs.loc[s1_obs["hit"], "gene"], s1_obs.loc[s1_obs["hit"], "feature"]))
    assert adj <= BOUNDS["culture_confounded_in_md_strict_max"]


def test_msi_burden_and_drop_rule(toy, s1_obs):
    """DMG_BURDEN is collinear enough with MSI (R² > drop threshold) to be dropped when d = SIG:MSI; MSI passengers
    are not, so the MSI covariate adjusts them away."""
    thr = DP["confounders_s5"]["drop_assoc_r2_gt"]
    msi = toy.features[f"{FAM['sig']}:{NF['msi']}"]
    ok = msi.notna()
    r2 = np.corrcoef(msi[ok], toy.covariates.loc[ok, "DMG_BURDEN"])[0, 1] ** 2
    assert r2 > thr
    rows = s1_obs[s1_obs["feature"] == f"{FAM['sig']}:{NF['msi']}"]
    assert rows["dropped"].str.contains("MSI").all() and rows["dropped"].str.contains("DMG_BURDEN").all()
    pas = s1_obs[s1_obs["feature"].isin(toy.answer_key["msi_passengers"])]
    assert not pas["dropped"].str.contains("MSI").any()
    msi_genes = set(toy.plants("msi_like")["gene"])
    ph = pas[pas["hit"] & pas["gene"].isin(msi_genes)]
    assert ph["gene"].nunique() <= BOUNDS["msi_like_passenger_lof_selected_max"]


def test_lineage_dummy_feature_is_collinear_with_design(toy, s1_obs):
    ld = toy.plants("lineage_dummy_feature").iloc[0]
    rows = s1_obs[s1_obs["feature"] == ld["feature"]]
    assert len(rows) and rows["collinear"].all() and not rows["hit"].any()
    # lineage-proxy SIL features straddle the collinearity threshold (some flagged, some tested)
    lp = s1_obs[s1_obs["feature"].isin(toy.answer_key["lineage_proxy_features"])].drop_duplicates("feature")
    assert lp["collinear"].any() and (~lp["collinear"]).any()
    lo = set(toy.plants("lineage_only")["gene"]) | {ld["gene"]}
    assert not s1_obs.loc[s1_obs["hit"], "gene"].isin(lo).any()


def test_self_null_carriers_not_dependent(toy):
    thr = DP["dependency_s2"]["threshold"]
    sn = toy.answer_key["self_null"]
    for g, f in zip(sn["gene"], sn["feature"]):
        car = toy.features[f] == 1
        assert car.sum() >= DP["dictionary_s3"]["min_carriers"]
        assert toy.effects.loc[car, g].median() > thr
        assert toy.effects.loc[~car & toy.features[f].notna(), g].median() < thr
        assert f not in toy.features_of(g)


def test_six_carrier_plant(toy):
    r = toy.plants("six_carrier_strong").iloc[0]
    car = toy.features[r["feature"]] == 1
    assert car.sum() == CFG["generator"]["effects"]["six_carrier_n"]
    assert toy.covariates.loc[car, "lineage"].nunique() == car.sum()


def test_cutting_and_addiction_oracle(toy):
    cut = toy.plants("cutting_artefact").query("subkind == 'core'")
    got = sum(ts.oracle_amplicon(toy, g, f) == "cutting artefact" for g, f in zip(cut["gene"], cut["feature"]))
    assert got >= BOUNDS["cutting_core_label_min"]
    own = toy.plants("own_locus_cn")
    got = sum(ts.oracle_amplicon(toy, g, f) == "amplified-driver addiction" for g, f in zip(own["gene"], own["feature"]))
    assert got >= BOUNDS["own_locus_addiction_label_min"]


def test_cross_lineage_oracle(toy):
    pan = toy.plants("pan_lineage")
    lab = [ts.oracle_cross_lineage(toy, g, f) for g, f in zip(pan["gene"], pan["feature"])]
    assert np.mean([x == "cross-lineage" for x in lab]) >= BOUNDS["cross_lineage_pan_frac_min"]
    rs = toy.plants("lineage_restricted")
    lab = [ts.oracle_cross_lineage(toy, g, f) for g, f in zip(rs["gene"], rs["feature"])]
    assert sum(x == "cross-lineage" for x in lab) <= BOUNDS["cross_lineage_restricted_max"]


def test_ridge_variant_shows_dense_predictor_inflation():
    """DESIGN §13.1 ridge predictor: S1-pred on a dense linear predictor flags features it uses on in-sample chance
    structure. This is why prediction lists are hard-bounded only on the smooth toy (docs/toy_answer_key.md)."""
    t = ts.make_toy(scenario="null", predictor="ridge")
    sp = ts.oracle_s1(t, "pred", "full")
    so = ts.oracle_s1(t, "obs", "full")
    assert sp["hit"].sum() > 10 * max(1, so["hit"].sum())


def test_scoring_helpers(toy):
    pl = toy.plants("pan_lineage").head(2)
    hits = pd.DataFrame({"gene": pl["gene"], "feature": pl["feature"], "beta": [-1.0, 1.0]})
    rec = ts.recovered(hits, toy, "pan_lineage")
    assert rec.sum() == 1  # the second has the wrong sign
    bg = toy.genes.index.difference(toy.plants()["gene"])[0]
    junk = pd.DataFrame({"gene": [bg], "feature": [toy.feature_meta.index[0]], "beta": [-1.0]})
    f = ts.fdp_true(pd.concat([hits, junk]), toy)
    assert f["n_hits"] == 3 and f["n_false"] == 1
    assert ts.validate_outputs({}) and len(ts.validate_outputs({})) == len(ts.REQUIRED_OUTPUTS)


def test_docs_list_every_bound():
    doc = (P.DOCS / "toy_answer_key.md").read_text(encoding="utf-8")
    missing = [k for k in BOUNDS if not k.startswith("_") and k not in doc]
    assert not missing, missing


# =========================================================================================================
# Part B: pipeline answer key (skipped until biomark.pipeline exists)
# =========================================================================================================

PIPE = ts.get_pipeline()
RUN_PIPE = os.environ.get("BIOMARK_TOY_PIPELINE") == "1"
needs_pipeline = pytest.mark.skipif(
    PIPE is None or not RUN_PIPE,
    reason=(f"{ts.PIPELINE_MODULE}.{ts.PIPELINE_ENTRY} not implemented yet (toy answer key waits for it)"
            if PIPE is None else
            "toy pipeline answer key is slow (full arms + fixed-threshold null + 20 pure-null toys); "
            "run with BIOMARK_TOY_PIPELINE=1 (gate runs: before freeze-bounds and after every src/ fix)"))


@pytest.fixture(scope="module")
def out_truth(toy):
    return PIPE(toy)


@pytest.fixture(scope="module")
def out_nocentre(toy):
    return PIPE(toy, fold_centre=False, null_B=0, arms=("EN-pred",))


@pytest.fixture(scope="module")
def ridge_toy():
    return ts.make_toy(predictor="ridge")


@pytest.fixture(scope="module")
def out_ridge(ridge_toy):
    return PIPE(ridge_toy)


@pytest.fixture(scope="module")
def pure_null_outs():
    B = int(os.environ.get("BIOMARK_TOY_PURE_NULL_B", BOUNDS["pure_null_B"]))
    base = P.seeds()[DP["toy_s13"]["null_seed_key"]]
    outs = []
    for b in range(B):
        t = ts.make_toy(seed=base + b, scenario="null")
        outs.append((t, PIPE(t, null_B=0)))
    return outs


def _sel(out, arm):
    return ts.list_hits(out, arm)


def _has(df, gene, feats):
    return ((df["gene"].astype(str) == gene) & df["feature"].astype(str).isin(feats)).any()


@needs_pipeline
def test_pipeline_interface(out_truth, toy):
    assert not ts.validate_outputs(out_truth), ts.validate_outputs(out_truth)
    assert out_truth["meta"]["design_params_sha256"] == ts.design_params_sha256()
    el = out_truth["eligibility"]
    elig = set(el.loc[el["eligible"].astype(bool), "gene"].astype(str))
    assert set(out_truth["arms"]["gene"].astype(str)) <= elig
    assert set(out_truth["arms"]["arm"]) == set(ts.ARMS)


@needs_pipeline
def test_pipeline_pan_lineage_recall(out_truth, toy):
    for name in ("EN-pred selected", "MD-strict"):
        assert ts.recovered(_sel(out_truth, name), toy, "pan_lineage").mean() >= BOUNDS["pan_lineage_recall_min"], name
    assert ts.recovered(_sel(out_truth, "S1-obs"), toy, "pan_lineage").mean() >= BOUNDS["pan_lineage_recall_s1_obs_min"]


@needs_pipeline
def test_pipeline_own_locus_cn(out_truth, toy):
    assert ts.recovered(out_truth["md_strict"], toy, "own_locus_cn").mean() >= BOUNDS["own_locus_cn_recall_min"]


@needs_pipeline
def test_pipeline_msi_signature_and_drop_rule(out_truth, toy):
    md = out_truth["md_strict"]
    assert ts.recovered(md, toy, "msi_like", accept_proxies=False).sum() >= BOUNDS["msi_like_signature_recovered_min"]
    rows = md[md["feature"] == f"{FAM['sig']}:{NF['msi']}"]
    assert rows["dropped_cov"].astype(str).str.contains("MSI").all(), "SIG:MSI must drop the MSI covariate"
    en = _sel(out_truth, "EN-pred selected")
    msi_genes = set(toy.plants("msi_like")["gene"])
    pas = en[en["feature"].isin(toy.answer_key["msi_passengers"]) & en["gene"].astype(str).isin(msi_genes)]
    assert pas["gene"].nunique() <= BOUNDS["msi_like_passenger_lof_selected_max"]


@needs_pipeline
def test_pipeline_co_location(out_truth, toy):
    lr = out_truth["locus_resolution"].astype({"gene": str, "feature": str})
    co = toy.answer_key["colocated"]
    ok = pass_gene = true_gene = tight_gene = 0
    for _, r in co.iterrows():
        sub = lr[(lr["gene"] == r["gene"]) & (lr["locus_group"] == r["locus_group"])]
        pas_gene = ((sub["feature"] == r["passenger"]) & (sub["level"] == "gene")).any()
        tr_gene = ((sub["feature"] == r["true_feature"]) & (sub["level"] == "gene")).any()
        if r["kind"] == "co_located":
            ok += len(sub) > 0 and not pas_gene
            pass_gene += pas_gene
            true_gene += tr_gene
        else:
            tight_gene += int((sub["level"] == "gene").sum())
    assert ok >= BOUNDS["co_located_locus_level_min"]
    assert pass_gene <= BOUNDS["co_located_passenger_gene_level_max"]
    assert true_gene >= BOUNDS["co_located_true_gene_level_min"]
    assert tight_gene <= BOUNDS["co_located_tight_gene_level_max"]


@needs_pipeline
def test_pipeline_culture_confounded(out_truth, toy):
    traps = toy.answer_key["culture"]
    key = set(zip(traps["gene"], traps["feature"]))
    md = out_truth["md_strict"]
    assert sum((str(g), f) in key for g, f in zip(md["gene"], md["feature"])) <= BOUNDS["culture_confounded_in_md_strict_max"]
    s1 = _sel(out_truth, "S1-obs")
    assert sum((str(g), f) in key for g, f in zip(s1["gene"], s1["feature"])) <= BOUNDS["culture_confounded_in_md_strict_max"]
    de = out_truth["confounding_delta"]
    rep = de[[(str(g), f) in key for g, f in zip(de["gene"], de["feature"])] & (de["block"] == "culture")]
    if len(rep):
        assert rep["flagged"].astype(bool).mean() >= BOUNDS["culture_flagged_if_reported_frac_min"]


@needs_pipeline
def test_pipeline_lineage_only_and_lineage_dummy(out_truth, toy):
    md = out_truth["md_strict"].astype({"gene": str})
    lo = set(toy.plants("lineage_only")["gene"])
    assert md["gene"].isin(lo).sum() <= BOUNDS["lineage_only_in_md_strict_max"]
    gc = out_truth["gene_classes"].astype({"gene": str}).set_index("gene")["failure_class"].astype(str)
    assert gc.reindex(sorted(lo)).eq("3").mean() >= BOUNDS["lineage_only_class3_frac_min"]
    ld = toy.plants("lineage_dummy_feature").iloc[0]
    assert ld["feature"] in set(out_truth["collinear_with_design"]["feature"].astype(str))
    assert (md["feature"] == ld["feature"]).sum() <= BOUNDS["lineage_dummy_in_md_strict_max"]


@needs_pipeline
def test_pipeline_fold_offset_artefact(out_truth, out_nocentre, toy):
    fa = toy.answer_key["fold_artefact"]
    sel_c = _sel(out_truth, "EN-pred selected")
    sel_n = _sel(out_nocentre, "EN-pred selected")
    after = sum(_has(sel_c, g, [f]) for g, f in zip(fa["gene"], fa["feature"]))
    before = sum(_has(sel_n, g, [f]) for g, f in zip(fa["gene"], fa["feature"]))
    assert after <= BOUNDS["fold_offset_selected_after_centring_max"]
    assert before >= BOUNDS["fold_offset_selected_without_centring_min"], "positive control of the artefact failed"


@needs_pipeline
def test_pipeline_genomic_with_expression_proxy_stage_G(out_truth, toy):
    en = _sel(out_truth, "EN-pred selected")
    got = ts.recovered(en, toy, "genomic_with_expression_proxy", accept_proxies=False, stage="G").sum()
    assert got >= BOUNDS["genomic_stage_G_recovered_min"]
    assert ts.recovered(en, toy, "prog_driven").sum() >= BOUNDS["prog_driven_en_pred_min"]


@needs_pipeline
def test_pipeline_six_carrier_strong(out_truth, toy):
    r = toy.plants("six_carrier_strong").iloc[0]
    a = out_truth["arms"]
    row = a[(a["arm"] == "EN-pred") & (a["gene"].astype(str) == r["gene"]) & (a["feature"] == r["feature"])]
    assert len(row) == 1, "six-carrier plant missing from U (in-fit carrier floor 2 must keep it)"
    assert float(row["rank_stat"].iloc[0]) > BOUNDS["six_carrier_selected_subsample_frac_min"]


@needs_pipeline
def test_pipeline_self_null(out_truth, toy):
    sn = toy.answer_key["self_null"]
    a = out_truth["arms"]
    inU = a[[(str(g), f) in set(zip(sn["gene"], sn["feature"])) for g, f in zip(a["gene"], a["feature"])]]
    assert len(inU) <= BOUNDS["self_null_in_universe_max"]
    s = out_truth["self_null"].astype({"gene": str})
    s = s[[(g, f) in set(zip(sn["gene"], sn["feature"])) for g, f in zip(s["gene"], s["feature"])]]
    assert len(s) >= 1 and (s["frac_dependent"] <= BOUNDS["nc3_frac_dependent_max"]).all()


@needs_pipeline
def test_pipeline_wrong_sign_cgt_is_sign_error(out_truth, toy):
    truth = toy.answer_key["cgt_truth"]
    wrong = truth[(truth["sign"] == "S") & (truth["true_sign"] == "W")]["test_id"]
    cs = out_truth["cgt_scores"]
    rows = cs[cs["test_id"].isin(wrong) & cs["hit"].astype(bool)]
    assert len(rows) >= 1, "the wrong-sign pair must be a hit (any sign) in at least one arm, not a miss"
    assert not rows["sign_correct"].astype(bool).any()


@needs_pipeline
def test_pipeline_cross_lineage_labels(out_truth, toy):
    cl = out_truth["cross_lineage"].astype({"gene": str}).set_index(["gene", "feature"])["label"]
    pan = toy.plants("pan_lineage")
    lab = [cl.get((g, f)) for g, f in zip(pan["gene"], pan["feature"])]
    lab = [x for x in lab if x is not None]
    assert lab and np.mean([x == "cross-lineage" for x in lab]) >= BOUNDS["cross_lineage_pan_frac_min"]
    rs = toy.plants("lineage_restricted")
    assert sum(cl.get((g, f)) == "cross-lineage" for g, f in zip(rs["gene"], rs["feature"])) <= BOUNDS["cross_lineage_restricted_max"]


@needs_pipeline
def test_pipeline_cutting_and_addiction(out_truth, toy):
    am = out_truth["amplicon_test"].astype({"gene": str}).set_index(["gene", "feature"])["label"]
    cut = toy.plants("cutting_artefact").query("subkind == 'core'")
    assert sum(am.get((g, f)) == "cutting artefact" for g, f in zip(cut["gene"], cut["feature"])) >= BOUNDS["cutting_core_label_min"]
    own = toy.plants("own_locus_cn")
    assert sum(am.get((g, f)) == "amplified-driver addiction" for g, f in zip(own["gene"], own["feature"])) >= BOUNDS["own_locus_addiction_label_min"]


@needs_pipeline
def test_pipeline_predictor_blind_and_p53(out_truth, toy):
    assert ts.recovered(_sel(out_truth, "EN-pred selected"), toy, "predictor_blind").sum() <= BOUNDS["predictor_blind_en_pred_selected_max"]
    assert ts.recovered(_sel(out_truth, "S1-obs"), toy, "predictor_blind").mean() >= BOUNDS["predictor_blind_s1_obs_recall_min"]
    assert ts.recovered(out_truth["md_strict"], toy, "p53_w").sum() >= BOUNDS["p53_md_strict_min"]


@needs_pipeline
def test_pipeline_gate(out_truth, toy):
    gate = out_truth["gate"].merge(toy.controls, on="control_id")
    gate = gate[gate["gate_eligible"]]
    assert gate["recovered"].astype(bool).mean() >= BOUNDS["gate_overall_recall_min"]
    for fam_, sub in gate.groupby("family"):
        if len(sub) >= DP["gate_s9"]["family_min_eligible_rows"]:
            assert sub["recovered"].astype(bool).mean() >= BOUNDS["gate_family_recall_min"], fam_


def _error_rate_checks(out, toy, hard_lists):
    nc = out["null_calibration"].set_index("list")
    report = {}
    for name in ts.LISTS:
        f = ts.fdp_true(ts.list_hits(out, name), toy)
        report[name] = f["fdp"]
        if name in hard_lists:
            assert f["fdp"] <= BOUNDS["planted_fdp_true_max"], (name, f)
            if name in nc.index:
                assert nc.at[name, "fdr_emp_fixed"] >= f["fdp"] + BOUNDS["fdr_emp_fixed_minus_fdp_true_min"], name
    print("FDP_true per list:", report)


@needs_pipeline
def test_pipeline_error_rates_smooth_toy(out_truth, toy):
    _error_rate_checks(out_truth, toy, BOUNDS["_lists_hard_on_smooth_toy"])


@needs_pipeline
def test_pipeline_error_rates_ridge_toy(out_ridge, ridge_toy):
    _error_rate_checks(out_ridge, ridge_toy, BOUNDS["_lists_hard_on_ridge_toy"])


@needs_pipeline
def test_pipeline_pure_null(pure_null_outs):
    any_hit = {name: [] for name in ts.LISTS}
    for t, out in pure_null_outs:
        for name in ts.LISTS:
            any_hit[name].append(ts.fdp_true(ts.list_hits(out, name), t)["n_hits"])  # neutral pairs excluded
    for name in BOUNDS["pure_null_bh_lists"]:
        p_any = np.mean([h > 0 for h in any_hit[name]])
        assert p_any <= BOUNDS["pure_null_bh_lists_p_any_hit_max"], (name, p_any)
    print("pure-null mean hits (EN lists reported):", {k: float(np.mean(v)) for k, v in any_hit.items()})
