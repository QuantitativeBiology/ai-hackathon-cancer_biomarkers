"""Hierarchical stability EN (src/biomark/en.py): FWL equivalence, subsamples, selection, ranking, EN-pred blindness."""
import numpy as np
import pandas as pd
import pytest
from sklearn.linear_model import ElasticNet

from biomark import arms as A
from biomark import en as EN
from biomark import guard
from biomark import stats as ST

from synth_data import make_data


def test_fwl_en_equals_en_with_huge_scaled_unpenalised_columns():
    rng = np.random.default_rng(0)
    n, p = 160, 8
    lin = rng.integers(0, 3, n)
    U = np.column_stack([np.ones(n), lin == 1, lin == 2, rng.normal(size=n)]).astype(float)
    X = rng.normal(size=(n, p)) + 0.5 * U[:, [3]]
    y = X[:, 0] * 0.8 - X[:, 3] * 0.5 + U @ np.array([0.3, 1.0, -1.0, 0.7]) + rng.normal(size=n)
    alpha, l1 = 0.05, 0.9
    Q = ST.orth(U)
    Xt, yt = ST.resid(Q, X), ST.resid(Q, y[:, None])[:, 0]
    b_fwl = EN.enet_gram(Xt.T @ Xt, Xt.T @ yt, yt, n, [alpha], l1, 1e-14, 100000)[:, 0]
    c = 1e4
    big = ElasticNet(alpha=alpha, l1_ratio=l1, fit_intercept=False, tol=1e-14, max_iter=1000000)
    big.fit(np.column_stack([X, c * U]), y)
    assert np.allclose(big.coef_[:p], b_fwl, atol=1e-5)
    assert (np.abs(b_fwl) > 0).sum() >= 2


def test_enet_gram_matches_sklearn_path():
    rng = np.random.default_rng(1)
    n, p = 100, 12
    X = rng.normal(size=(n, p))
    y = X[:, 1] - 0.5 * X[:, 4] + rng.normal(size=n)
    grid = EN.alpha_grid(X.T @ y, n, 0.9, 25, 1e-2)
    ours = EN.enet_gram(X.T @ X, X.T @ y, y, n, grid, 0.9, 1e-10, 100000)
    from sklearn.linear_model import enet_path
    _, ref, _ = enet_path(X, y, l1_ratio=0.9, alphas=grid, tol=1e-10, max_iter=100000)
    assert np.allclose(ours, ref, atol=1e-6)


def test_subsamples_complementary_patient_grouped_lineage_stratified():
    data = make_data(seed=2)
    subs = EN.make_subsamples(data.lines, data.cov_df["lineage"], data.cov_df["PatientID"], 11, 25)
    assert len(subs) == 50
    pat = data.cov_df["PatientID"]
    lin = data.cov_df["lineage"]
    for i in range(25):
        a, b = set(subs[2 * i]), set(subs[2 * i + 1])
        assert not a & b and a | b == set(data.lines)
        assert not set(pat[list(a)]) & set(pat[list(b)])            # patients never split
        modal = lin.groupby(pat).agg(lambda v: v.value_counts().sort_index().idxmax())
        for L in lin.unique():
            pa = [p for p in modal.index[modal == L] if p in set(pat[list(a)])]
            na = pat.isin(pa).sum()
            nb = pat.isin(modal.index[modal == L]).sum() - na
            assert abs(na - nb) <= pat.value_counts().max()          # balanced per stratum (<= one patient)
    assert subs == EN.make_subsamples(data.lines, lin, pat, 11, 25)  # deterministic


def test_planted_feature_selected_and_deterministic():
    data = make_data(seed=3, plant_beta=-1.2)
    U = A.universe(data)
    tab1, lam1, _ = A.run_en(data, ("EN-obs",), U=U)
    tab2, lam2, _ = A.run_en(data, ("EN-obs",), U=U)
    pd.testing.assert_frame_equal(tab1, tab2)
    g0 = tab1[(tab1["gene"] == data.genes[0])].set_index("feature")
    assert bool(g0.at["HOT:9200001", "selected"]) and g0.at["HOT:9200001", "beta"] < 0
    assert g0.at["HOT:9200001", "selection_freq"] >= 0.6
    # stage of each row = stage of its feature; frequency within [0, 1]; ranks unique among rankable rows
    assert (tab1["stage"] == data.fmeta["stage"].reindex(tab1["feature"]).to_numpy()).all()
    assert tab1["selection_freq"].between(0, 1).all()
    rk = tab1["rank"].dropna()
    assert rk.is_unique and (tab1.loc[tab1["selection_freq"] == 0, "rank"].isna()).all()
    # fixed-lambda refit reproduces the selection exactly
    lam = {"EN-obs": {r.gene: {"G": r.lam_G, "X": r.lam_X} for r in lam1.itertuples()}}
    tab3, _, _ = A.run_en(data, ("EN-obs",), U=U, lam=lam)
    assert (tab3["selection_freq"].to_numpy() == tab1["selection_freq"].to_numpy()).all()


def test_rank_en_rules():
    rk = EN.rank_en(np.array([0.8, 0.8, 0.0, 0.5, 0.8]), np.array([0.1, -0.3, 0.9, 0.2, 0.1]),
                    np.array(["b", "a", "c", "d", "a2"]))
    assert np.isnan(rk[2])
    assert list(rk[[1, 4, 0, 3]]) == [1, 2, 3, 4]       # freq desc, |beta| desc, id asc


def test_en_pred_is_blind_to_observed_effects(monkeypatch):
    def refuse(*a, **k):
        raise AssertionError("EN-pred read observed effects")

    monkeypatch.setattr(guard, "load_effects", refuse)
    data = make_data(seed=4, obs_loader=lambda: guard.load_effects(scope="dev"))
    U = A.universe(data)
    tab, lam, _ = A.run_en(data, ("EN-pred",), U=U)
    assert len(tab) and data._obs is None
    with pytest.raises(AssertionError):                    # the patch is live: an obs arm does read effects
        A.run_en(data, ("EN-obs",), U=U)
    # identical row sets for both EN arms: the rows come from the NaN mask only
    assert (lam["n_rows"].to_numpy() == data.rows_ok.sum(0)[[data.genes.index(g) for g in lam["gene"]]]).all()
