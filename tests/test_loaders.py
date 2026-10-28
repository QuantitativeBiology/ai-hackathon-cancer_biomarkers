"""Shapes and invariants of the Part 1 read-only loaders."""
import numpy as np
import pytest

from biomark import p1
from biomark import paths as P


@pytest.fixture(scope="module")
def oof():
    return p1.oof_predictions("m5")


def test_splits_counts():
    sp = p1.splits()
    assert sp.shape[0] == 1292
    assert sp["role"].value_counts().to_dict() == {"dev": 1022, "test": 255, "noncancerous": 15}
    assert sorted(p1.folds().value_counts().to_dict().items()) == [(0, 208), (1, 202), (2, 201), (3, 202),
                                                                     (4, 209)]
    assert (sp.loc[sp.role != "dev", "cv_fold"] == -1).all()
    assert {"PatientID", "OncotreeLineage"} <= set(sp.columns)


def test_headline_oof_shape_no_nan(oof):
    assert oof.shape == (1022, 18434)
    assert not oof.isna().any().any()
    assert set(oof.index) == set(p1.dev_ids())
    assert oof.columns.map(str.isdigit).all()


def test_fold_centre(oof):
    cols = list(oof.columns[:500])
    c = p1.fold_centre(oof[cols])
    assert c.shape == (1022, 500) and (c.dtypes == oof.dtypes.iloc[0]).all()
    means = c.groupby(p1.folds().reindex(c.index)).mean()
    assert np.abs(means.to_numpy()).max() < 1e-5
    # within-fold differences are preserved
    f = p1.folds().reindex(oof.index)
    a = oof[cols][f == 0]
    b = c[f == 0]
    np.testing.assert_allclose((a - a.iloc[0]).to_numpy(), (b - b.iloc[0]).to_numpy(), atol=1e-5)
    with pytest.raises(KeyError):
        p1.fold_centre(p1.test_predictions(columns=cols[:2]))


def test_component_and_test_predictions():
    cols = ["7157", "4609"]
    comp = p1.component_oof(columns=cols)
    assert set(comp) == {"m1", "m2n", "m3n", "t"}
    assert all(v.shape == (1022, 2) for v in comp.values())
    tp = p1.test_predictions()
    assert tp.shape == (255, 18434) and set(tp.index) == set(p1.test_ids())
    assert not tp.isna().any().any()


def test_gene_tables():
    gm = p1.gene_map()
    assert gm.shape[0] == 21269 and gm.loc["7157", "symbol"] == "TP53"
    assert p1.entrez_of(["TP53"]).iloc[0] == "7157"
    gc = p1.gene_classes()
    assert gc["variable"].sum() == 2159 and gc.index.dtype == object
    t = p1.gene_tiers("m5")
    assert t["tier"].value_counts()["high"] == 1632
    assert t.index.isin(gc.index).all()


def test_feature_blocks_and_raw():
    ids = list(p1.dev_ids()[:3])
    for name in ["E", "C", "M", "H", "F_gene", "meta", "meta_prov", "composite", "C_instab", "M_burden"]:
        b = p1.feature_block(name, ids=ids, columns=["7157"] if name in ("E", "C", "M") else None)
        assert b.shape[0] == 3 and list(b.index) == ids
    assert p1.feature_block("composite").shape[0] == 2395
    assert p1.raw_matrix("expression", columns=["7157"]).shape == (2130, 1)
    assert p1.raw_matrix("cn", columns=["7157"]).shape == (1763, 1)
    assert p1.raw_matrix("damaging", columns=["7157"]).shape == (2395, 1)
    assert p1.raw_matrix("hotspot").shape == (2395, 165)
    assert p1.raw_matrix("fusion", columns=["ModelID"]).shape[0] == 94718
    with pytest.raises(KeyError):
        p1.feature_block("E", columns=["not_a_gene"])


def test_model_tables():
    m = p1.model_table(["OncotreeLineage", "GrowthPattern"])
    assert m.shape == (2554, 2)
    mc = p1.model_condition(["CellFormat", "PlateCoating"])
    assert mc.shape[0] == 3385 and "None" in set(mc["PlateCoating"])


def test_screen_covariates_cover_population():
    sc = p1.screen_covariates()
    dt = sc[sc["role"].isin(["dev", "test"])]
    assert len(dt) == 1277
    need = ["Library", "ScreenDays", "ScreenType", "ModelConditionID", "CellFormat", "GrowthMedia",
            "SerumFreeMedia", "PlateCoating", "model_GrowthPattern", "model_ModelType",
            "model_OnboardedMedia", "model_FormulationID", "model_SerumFreeMedia"]
    assert set(need) <= set(sc.columns)
    n_no_lib = int(dt["Library"].isna().sum())
    print(f"\n[screen_covariates] dev+test lines lacking library: {n_no_lib}; "
          f"lacking GrowthMedia (screened condition): {int(dt['GrowthMedia'].isna().sum())}; "
          f"multi-screen lines: {int((dt['n_screens'] > 1).sum())}")
    assert n_no_lib == 0
    # the screen-derived library agrees with Part 1's NaN-pattern groups (ARI 0.949 in the brief)
    from sklearn.metrics import adjusted_rand_score
    lg = p1.library_groups()
    j = lg.join(sc["Library"], how="inner")
    assert len(j) == 1022
    ari = adjusted_rand_score(j["lib_group"], j["Library"])
    print(f"[screen_covariates] ARI(lib_group, Library) on dev = {ari:.3f}")
    assert ari > 0.8
    # library-linked screen length (brief: median 21 Avana, 24 Humagne-CD, 14 KY)
    med = dt.groupby("Library")["ScreenDays"].median()
    assert med["Avana"] == 21 and med["KY"] == 14


def test_paths_and_seeds():
    s = P.seeds()
    assert s["seeds"] == [0, 1, 2, 3, 4] and s["subsample_seed"] == 20261028 and s["permutation_seed"] == 20261029
    assert P.P1.exists() and P.P1_CRISPR_PARQUET.exists()
    assert all((P.P1_RUNS / r / "pred.parquet").exists() for r in P.P1_RUN_IDS.values())
