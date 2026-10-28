"""DESIGN §13.2 test_bootstrap_clusters: the paired patient bootstrap of §6.1 (endpoints.bootstrap_rates).

Checks, on synthetic R-like data (no real effects): patients are the resampling unit (all lines of a drawn patient
come together), each drawn copy gets its own cluster id, the draws are seeded and paired across lists, a fit on a
bootstrap index equals the CR2 + ν test on the physically duplicated data with copy-level clusters, and a pair with
< 2 distinct carrier patients in a draw is degenerate and not replicated.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from biomark import endpoints as EP
from biomark import stats as ST


def _R(n_pat=120, seed=5, extra_carrier_patients=None):
    rng = np.random.default_rng(seed)
    sizes = rng.choice([1, 2, 3], size=n_pat, p=[0.7, 0.2, 0.1])
    pats = np.repeat(np.arange(n_pat), sizes)
    n = len(pats)
    ids = [f"R{i:04d}" for i in range(n)]
    cov = pd.DataFrame({"lineage": rng.choice(["L1", "L2", "L3"], n),
                        "library_group": rng.choice(["Avana-only", "KY-only"], n), "cellformat": "adherent",
                        "serumfree": "F", "medium_family": rng.choice(["RPMI", "DMEM"], n),
                        "doubling_time": rng.normal(40, 5, n), "DT_unknown": 0.0, "MSI": "MSS",
                        "DMG_BURDEN": rng.normal(5, 0.5, n), "tp53_status": rng.choice(["mut", "WT"], n),
                        "PatientID": [f"P{p}" for p in pats]}, index=ids)
    X = (rng.random((n, 3)) < 0.2).astype(float)
    X[:, 2] = 0.0
    X[pats == 0, 2] = 1.0                       # HOT:3 carried by one patient only (all its lines)
    Y = rng.normal(0, 0.4, (n, 2))
    Y[:, 0] += -1.0 * X[:, 0]
    feats = ["HOT:1", "HOT:2", "HOT:3"]
    Yd = pd.DataFrame(Y, index=ids, columns=["g1", "g2"])
    return EP.RFitter(X, pd.Series(np.arange(3), index=feats), cov, None, None, Yd, cov["PatientID"].to_numpy())


class _Spy:
    def __init__(self, fitter):
        self.f, self.calls = fitter, []

    def __getattr__(self, k):
        return getattr(self.f, k)

    def fit(self, g, f, drops, idx=None, clusters=None, sign=0.0):
        self.calls.append((g, f, np.asarray(idx).copy(), np.asarray(clusters).copy()))
        return self.f.fit(g, f, drops, idx=idx, clusters=clusters, sign=sign)


def test_patients_are_the_unit_and_copies_get_new_cluster_ids():
    fit = _R()
    spy = _Spy(fit)
    lists = {"A": [("g1", "HOT:1", -1.0)], "B": [("g1", "HOT:2", -1.0)]}
    EP.bootstrap_rates(spy, lists, {}, {"HOT:1": True, "HOT:2": True}, B=15, seed=11)
    pats = fit.patients
    rows_of = {p: set(np.flatnonzero(pats == p)) for p in np.unique(pats)}
    n_pat = len(rows_of)
    by_draw = {}
    for g, f, idx, cl in spy.calls:
        by_draw.setdefault(f, []).append((idx, cl))
    assert len(by_draw["HOT:1"]) == len(by_draw["HOT:2"]) == 15
    for (i1, c1), (i2, c2) in zip(by_draw["HOT:1"], by_draw["HOT:2"]):
        np.testing.assert_array_equal(i1, i2)                     # paired: same draw for every list row
        np.testing.assert_array_equal(c1, c2)
    for idx, cl in by_draw["HOT:1"]:
        assert len(np.unique(cl)) == n_pat                        # one copy per draw slot = n_pat clusters
        for c in np.unique(cl):
            rows = idx[cl == c]
            p = pats[rows]
            assert len(set(p)) == 1                               # a cluster id never mixes patients
            assert set(rows) == rows_of[p[0]] and len(rows) == len(rows_of[p[0]])   # whole patient, once
        # a patient drawn k times appears as k distinct clusters
        n_copies = pd.Series([pats[idx[cl == c][0]] for c in np.unique(cl)]).value_counts()
        assert n_copies.max() >= 2                                # with replacement (n_pat = 120)
    # seeded: same seed -> same draws; different seed -> different draws
    spy2 = _Spy(fit)
    EP.bootstrap_rates(spy2, lists, {}, {"HOT:1": True, "HOT:2": True}, B=15, seed=11)
    for (_, _, i1, _), (_, _, i2, _) in zip(spy.calls, spy2.calls):
        np.testing.assert_array_equal(i1, i2)
    spy3 = _Spy(fit)
    EP.bootstrap_rates(spy3, lists, {}, {"HOT:1": True, "HOT:2": True}, B=3, seed=12)
    assert not np.array_equal(spy.calls[0][2], spy3.calls[0][2])


def test_bootstrap_fit_equals_duplicated_data_with_copy_clusters():
    fit = _R()
    pats = fit.patients
    up = np.unique(pats)
    rng = np.random.default_rng(3)
    draw = rng.choice(up, size=up.size, replace=True)
    idx = np.concatenate([np.flatnonzero(pats == p) for p in draw])
    cl = np.concatenate([np.full((pats == p).sum(), f"{p}#{j}", dtype=object) for j, p in enumerate(draw)])
    r = fit.fit("g1", "HOT:1", "", idx=idx, clusters=cl, sign=-1.0)
    # physically duplicated data, explicit design, brute CR2 via pair_test
    y = fit.Yv[idx, 0]
    d = fit.X[idx, 0]
    Z = fit.cov.Z(idx, fit.cov.design_names + fit.cov.droppable_names)
    ref = ST.pair_test(y, d, Z, cl.astype(str))
    for k in ("beta", "se", "nu", "df", "p"):
        assert r[k] == pytest.approx(ref[k], rel=1e-10), k
    # treating copies as one cluster (the wrong rule) gives a different SE whenever a patient was drawn twice
    wrong = ST.pair_test(y, d, Z, pats[idx])
    assert not np.isclose(wrong["se"], ref["se"], rtol=1e-6)


def test_single_carrier_patient_is_degenerate_and_never_replicated():
    fit = _R()
    lists = {"deg": [("g1", "HOT:3", -1.0)], "ok": [("g1", "HOT:1", -1.0)]}
    bt = EP.bootstrap_rates(fit, lists, {}, {"HOT:1": True, "HOT:3": True}, B=20, seed=2)
    assert (bt["rates"]["deg"] == 0).all()
    assert (bt["frac_degenerate"]["deg"] == 1).all()             # 0 or 1 carrier patients in every draw
    assert (bt["frac_degenerate"]["ok"] == 0).all()
    assert bt["rates"]["ok"].mean() > 0.9                         # planted determinant replicates
