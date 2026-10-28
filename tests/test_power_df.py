"""DESIGN §13.2 "power df = test df" (§6.1: "Power = 1 − F_t(t_0.95,ν_R; ncp) with ν_R from §5.3 on R carriers, so
the test and its power share one df rule"; §5.3: df = min(n − rank(Z) − 2, ν) for every test).

For each (g, d), the df used by ``endpoints.pool_power`` (column ``nu_R``) must equal the df of the R test that the
replication run will compute (``endpoints.RFitter.fit(...)["df"]``) on the same R lines with the same drops.
Synthetic R-like data only.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from biomark import endpoints as EP


def _setup(n_pat=90, seed=7, gene_covs=False, n_cont=1):
    rng = np.random.default_rng(seed)
    pats = np.repeat(np.arange(n_pat), rng.choice([1, 2], size=n_pat, p=[0.8, 0.2]))
    n = len(pats)
    ids = [f"R{i:04d}" for i in range(n)]
    cov = pd.DataFrame({"lineage": rng.choice(["L1", "L2", "L3", "L4"], n),
                        "library_group": rng.choice(["Avana-only", "KY-only"], n), "cellformat": "adherent",
                        "serumfree": "F", "medium_family": rng.choice(["RPMI", "DMEM"], n),
                        "doubling_time": rng.normal(40, 5, n), "DT_unknown": 0.0, "MSI": "MSS",
                        "DMG_BURDEN": rng.normal(5, 0.5, n), "tp53_status": rng.choice(["mut", "WT"], n),
                        "PatientID": [f"P{p}" for p in pats]}, index=ids)
    Xb = (rng.random((n, 3)) < np.array([0.04, 0.1, 0.3])).astype(float)
    Xc = rng.uniform(-1, 1, (n, n_cont))
    X = np.column_stack([Xb, Xc])
    feats = ["HOT:1", "HOT:2", "HOT:3"] + [f"EXPDOSE:{j}" for j in range(n_cont)]
    genes = ["g1", "g2"]
    Y = pd.DataFrame(rng.normal(0, 0.4, (n, 2)), index=ids, columns=genes)
    own = egate = None
    if gene_covs:
        own = pd.DataFrame(rng.normal(0, 0.3, (n, 2)) + 0.8 * X[:, [1, 1]], index=ids, columns=genes)
        egate = pd.DataFrame((rng.random((n, 2)) < 0.85).astype(float), index=ids, columns=genes)
    fit = EP.RFitter(X, pd.Series(np.arange(len(feats)), index=feats), cov, own, egate, Y, cov["PatientID"].to_numpy())
    UG = pd.DataFrame([(g, f) for g in genes for f in feats], columns=["gene", "feature"])
    key = pd.MultiIndex.from_arrays([UG["gene"], UG["feature"]])
    pool = EP.pool_power(UG, pd.Series(-0.3, index=key), pd.Series(0.4, index=genes), X, fit.fpos,
                         pd.Series([True] * 3 + [False] * n_cont, index=feats), fit.cov, fit.patients,
                         pd.DataFrame(True, index=range(n), columns=genes), pd.Series("", index=key))
    return fit, pool


def _pairs(fit, pool):
    out = []
    for r in pool.itertuples():
        t = fit.fit(r.gene, r.feature, "")
        out.append((r.gene, r.feature, float(r.nu_R), float(t["df"])))
    return out


def test_power_df_equals_test_df_shared_design_only():
    """No gene covariates: the R test's design equals the power design."""
    fit, pool = _setup(gene_covs=False)
    bad = [(g, f, a, b) for g, f, a, b in _pairs(fit, pool) if not np.isclose(a, b, rtol=1e-10)]
    assert not bad, f"power df != test df: {bad}"


def test_power_df_equals_test_df_with_gene_covariates():
    """With own CN / EXPGATE in the R test (§5.3 model), the power df must still be the test's df."""
    fit, pool = _setup(gene_covs=True)
    bad = [(g, f, round(a, 4), round(b, 4)) for g, f, a, b in _pairs(fit, pool) if not np.isclose(a, b, rtol=1e-6)]
    assert not bad, f"power df != test df ({len(bad)} pairs): {bad}"


def test_power_df_applies_the_n_minus_rank_cap():
    """df = min(n − rank(Z) − 2, ν): with few lines and a continuous d, ν exceeds n − rank − 2."""
    fit, pool = _setup(n_pat=14, seed=3, gene_covs=False)
    pairs = _pairs(fit, pool)
    capped = [(g, f, a, b) for g, f, a, b in pairs if f.startswith("EXPDOSE")]
    assert capped and all(b < a + 1e-12 for _, _, a, b in capped)       # the test df is the cap here
    bad = [(g, f, round(a, 4), round(b, 4)) for g, f, a, b in capped if not np.isclose(a, b, rtol=1e-10)]
    assert not bad, f"power uses nu without the n - rank(Z) - 2 cap: {bad}"
