"""Tier A statistics engine (src/biomark/stats.py) on synthetic data only."""
import numpy as np
import pandas as pd
import pytest
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

from biomark import arms as A
from biomark import stats as ST

from synth_data import make_data


def _design(n, rng):
    lin = rng.integers(0, 4, n)
    Z = np.column_stack([np.ones(n), *(lin == k for k in (1, 2, 3)), rng.normal(size=n)]).astype(float)
    return Z


def test_hc2_matches_statsmodels_singletons():
    rng = np.random.default_rng(1)
    n = 150
    Z = _design(n, rng)
    for trial in range(5):
        d = (rng.random(n) < 0.1 + 0.05 * trial).astype(float)
        y = Z @ rng.normal(size=Z.shape[1]) - 0.7 * d + rng.normal(size=n) * (1 + d)
        r = ST.pair_test(y, d, Z, np.arange(n).astype(str))
        fit = sm.OLS(y, np.column_stack([Z, d])).fit(cov_type="HC2")
        assert r["beta"] == pytest.approx(fit.params[-1], rel=1e-10, abs=1e-12)
        assert r["se"] == pytest.approx(fit.bse[-1], rel=1e-10)


def _cr2_brute(y, X, cl):
    """CR2 SE of the last coefficient: (X'X)^-1 Σ_c X_c' A_c e_c e_c' A_c X_c (X'X)^-1, A_c = (I - H_cc)^(-1/2)."""
    XtXi = np.linalg.inv(X.T @ X)
    H = X @ XtXi @ X.T
    e = y - X @ (XtXi @ X.T @ y)
    meat = np.zeros((X.shape[1], X.shape[1]))
    for c in np.unique(cl):
        idx = np.flatnonzero(cl == c)
        M = np.eye(idx.size) - H[np.ix_(idx, idx)]
        w, V = np.linalg.eigh(M)
        Ac = V @ np.diag(w ** -0.5) @ V.T
        u = X[idx].T @ Ac @ e[idx]
        meat += np.outer(u, u)
    V = XtXi @ meat @ XtXi
    return np.sqrt(V[-1, -1])


def test_cr2_and_nu_match_brute_force():
    rng = np.random.default_rng(2)
    n = 120
    Z = _design(n, rng)
    cl = np.arange(n).astype(object)
    cl[1] = cl[0]
    cl[3] = cl[2]
    cl[5] = cl[6] = cl[4]
    cl[10:14] = "four"
    for _ in range(4):
        d = (rng.random(n) < 0.15).astype(float)
        d[[0, 4]] = 1.0
        y = Z @ rng.normal(size=Z.shape[1]) + 0.5 * d + rng.normal(size=n)
        r = ST.pair_test(y, d, Z, cl.astype(str))
        se = _cr2_brute(y, np.column_stack([Z, d]), cl.astype(str))
        assert r["se"] == pytest.approx(se, rel=1e-8)
        Q = ST.orth(Z)
        dt = ST.resid(Q, d[:, None])[:, 0]
        W = pd.Series(dt ** 2).groupby(cl.astype(str)).sum().to_numpy()
        nu = max(1.0, W.sum() ** 2 / (W ** 2).sum() - 1)
        assert r["nu"] == pytest.approx(nu, rel=1e-10)
        assert r["df_nu"] == pytest.approx(min(n - Z.shape[1] - 2, nu))          # nu rule (sensitivity)
        car_cl = len(set(cl.astype(str)[d == 1]))                                 # distinct carrier patients
        assert r["n_carrier_clusters"] == car_cl
        assert r["df"] == pytest.approx(min(n - Z.shape[1] - 2, car_cl - 1))     # fallback (DESIGN 5.3)
        assert r["df_rule"] == "fallback_carrier_clusters"
        rc = ST.pair_test(y, d, Z, cl.astype(str), binary=False)                  # continuous rule unchanged
        assert rc["df"] == pytest.approx(r["df_nu"]) and rc["df_rule"] == "nu"


def test_vectorised_scan_equals_pair_test():
    data = make_data(seed=3)
    res = A.run_s1(data, ("S1-obs",))
    Y = data.target("obs")
    checked = 0
    for _, r in res.sample(25, random_state=0).iterrows():
        gj = data.genes.index(r["gene"])
        fj = data.fmeta.index.get_loc(r["feature"])
        from biomark import audit
        p = audit.pair_fit(data, r["gene"], r["feature"])
        if np.isnan(r["p"]):
            continue
        assert r["beta"] == pytest.approx(p["beta"], rel=1e-8, abs=1e-10)
        assert r["se"] == pytest.approx(p["se"], rel=1e-8)
        assert abs(np.log10(r["p"]) - np.log10(p["p"])) < 1e-6
        checked += 1
    assert checked >= 15
    assert Y.shape[1] == len(data.genes)


def test_covariate_drop_rules():
    data = make_data(seed=4)
    res = A.run_s1(data, ("S1-obs",))
    msi = res[res["feature"] == "SIG:MSI"]
    assert len(msi) and msi["dropped_cov"].str.contains("MSI").all()          # determinant == covariate
    dmg = res[res["feature"] == "SIG:DMG_BURDEN"]
    assert dmg["dropped_cov"].str.contains("DMG_BURDEN").all()
    lind = res[res["feature"] == "SIL:LINDUMMY"]
    assert len(lind) and lind["collinear_with_design"].all()                 # lineage dummy flagged ...
    assert not lind["dropped_cov"].str.contains("lineage").any()             # ... never dropped
    droppable = {"MSI", "DMG_BURDEN", "TP53", "screen_days", "own_CN", "EXPGATE"}
    for s in res["dropped_cov"]:
        assert set(filter(None, s.split("|"))) <= droppable                   # dropped_cov ⊆ droppable set
    own = res[(res["feature"] == f"AMP:{data.genes[1]}") & (res["gene"] == data.genes[1])]
    assert len(own) and own["dropped_cov"].str.contains("own_CN").all()      # own-locus CN drops own CN
    hot = res[res["feature"] == "HOT:9200001"]
    assert not hot["dropped_cov"].str.contains("MSI").any()


def test_bh_matches_statsmodels():
    rng = np.random.default_rng(5)
    p = rng.random(300) ** 3
    p[::17] = np.nan
    q = ST.bh(p)
    ok = ~np.isnan(p)
    assert np.allclose(q[ok], multipletests(p[ok], method="fdr_bh")[1])
    assert np.isnan(q[~ok]).all()
    c = ST.bh_cutoff(p, 0.05)
    assert ((q <= 0.05) == (p <= c))[ok].all()


def test_eb_mixture_recovers_oracle():
    rng = np.random.default_rng(6)
    m = 3000
    comp = rng.choice(3, size=m, p=[0.6, 0.3, 0.1])
    sig = np.array([0.05, 0.5, 2.0])[comp]
    theta = rng.normal(0, sig)
    s = rng.uniform(0.2, 0.6, m)
    b = theta + rng.normal(0, s)
    post = ST.eb_scale_mixture(b, s)["posterior_mean"]
    # oracle posterior mean with the true prior
    pis, sg = np.array([0.6, 0.3, 0.1]), np.array([0.05, 0.5, 2.0])
    var = sg[None] ** 2 + s[:, None] ** 2
    w = pis * np.exp(-0.5 * b[:, None] ** 2 / var) / np.sqrt(var)
    w /= w.sum(1, keepdims=True)
    oracle = (w * sg[None] ** 2 / var).sum(1) * b
    assert np.corrcoef(post, oracle)[0, 1] > 0.99


def test_power_invariant_to_feature_rescaling_and_df_shared():
    sig, n1, n0 = 0.4, np.array([5.0, 10.0]), np.array([200.0, 300.0])
    se = ST.se_r_pred(sig, np.array([True, True]), n1, n0)
    p1 = ST.power_t(np.array([-0.6, -0.3]), se, np.array([4.0, 9.0]))
    # continuous feature: rescaling d by c rescales beta by 1/c and SD_R(d~) by c -> same power
    c = 7.0
    se_c1 = ST.se_r_pred(sig, np.array([False]), sd_R=np.array([0.8]), n_R=np.array([200.0]))
    se_c2 = ST.se_r_pred(sig, np.array([False]), sd_R=np.array([0.8 * c]), n_R=np.array([200.0]))
    pa = ST.power_t(np.array([0.2]), se_c1, np.array([150.0]))
    pb = ST.power_t(np.array([0.2 / c]), se_c2, np.array([150.0]))
    assert np.allclose(pa, pb)
    assert np.all((p1 > 0) & (p1 < 1))
    rng = np.random.default_rng(7)
    dt = rng.normal(size=50)
    assert ST.nu_from(dt, np.arange(50)) == pytest.approx(max(1, (dt ** 2).sum() ** 2 / (dt ** 4).sum() - 1))


def test_scan_uses_the_binary_fallback_df():
    data = make_data(seed=9)
    res = A.run_s1(data, ("S1-obs",))
    b = res[res["feature"].str.startswith(("HOT:", "LOFBI:"))]
    c = res[res["feature"].str.startswith(("PROG:", "EXPDOSE:", "SIG:DMG"))]
    assert (b["df_rule"] == "fallback_carrier_clusters").all() and (c["df_rule"] == "nu").all()
    cap = b["n_obs"] - 0                                                      # df <= carrier clusters - 1
    assert (b["df"] <= b["n_carrier_clusters"] - 1 + 1e-9).all() and (b["df"] <= cap).all()
    assert np.allclose(c["df"], c["df_nu"])
    from biomark import audit
    r = b.iloc[0]
    p = audit.pair_fit(data, r["gene"], r["feature"])
    assert p["df"] == pytest.approx(r["df"]) and p["df_rule"] == "fallback_carrier_clusters"
