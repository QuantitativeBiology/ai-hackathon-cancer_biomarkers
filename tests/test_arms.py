"""Arms assembly, MD-strict, ranking/dedup and the joint-row permutation null (synthetic data only)."""
import numpy as np
import pandas as pd

from biomark import arms as A
from biomark import null as NL

from synth_data import make_data


def test_four_arms_share_one_universe_and_md_strict_rule():
    data = make_data(seed=5, plant_beta=-1.5)
    res = A.run_arms(data, arms=A.ARMS, workers=1)
    a = res["arms"]
    keys = {arm: set(zip(s["gene"], s["feature"])) for arm, s in a.groupby("arm")}
    assert set(keys) == set(A.ARMS)
    first = keys["S1-obs"]
    assert all(k == first for k in keys.values())                          # identical (g, d) universes
    assert not any(f.startswith("EXPDOSE:") and f != "EXPDOSE:" + g for g, f in first)  # own EXPDOSE only
    for g, f in first:                                                     # U ⊆ 𝒟_g
        assert data.fmeta.index.get_loc(f) in set(data.Dg[g])
    md_all = res["md_strict_all"]
    sel = a[(a["arm"] == "EN-pred") & a["selected"]]
    assert len(md_all) == len(sel)
    exp = md_all["same_sign"] & (md_all["q_sc"] < 0.05) & ~md_all["collinear_with_design"].astype(bool)
    assert (md_all["pass"] == exp.fillna(False)).all()
    assert set(zip(res["md_strict"]["gene"], res["md_strict"]["feature"])) <= set(zip(sel["gene"], sel["feature"]))
    assert ("9100001", "HOT:9200001") in set(zip(res["md_strict"]["gene"], res["md_strict"]["feature"]))
    for col in ("arm", "gene", "feature", "family", "access_class", "rank_stat", "rank", "d_type",
                "en_max_attainable_freq", "selected", "beta", "stage", "se_type", "rank_dedup"):
        assert col in a.columns
    s1 = a[a["arm"] == "S1-obs"]
    assert (s1.loc[s1["collinear_with_design"].astype(bool), "selected"] == False).all()  # noqa: E712


def test_strict_on_constructed_tables():
    en = pd.DataFrame({"arm": "EN-pred", "gene": ["g1", "g1", "g2", "g3"], "feature": ["a", "b", "a", "c"],
                       "selected": True, "selection_freq": 0.9, "beta": [-0.2, 0.3, -0.1, -0.4],
                       "maj_sign": [-1.0, 1.0, -1.0, -1.0], "stage": "G"})
    s1 = pd.DataFrame({"arm": "S1-obs", "gene": ["g1", "g1", "g2", "g3"], "feature": ["a", "b", "a", "c"],
                       "beta": [-0.5, -0.5, -0.3, -0.6], "se": 0.1, "df": 50.0, "nu": 20.0,
                       "p": [1e-6, 1e-6, 0.5, 1e-5], "n_obs": 100, "n_carriers": 10, "dropped_cov": "",
                       "collinear_with_design": [False, False, False, True], "se_type": "CR2"})
    m = A.strict(en, s1)
    assert list(m["pass"]) == [True, False, False, False]   # sign flip; not significant; collinear
    m2 = A.strict(en, s1, p_cut=0.6)
    assert list(m2["pass"]) == [True, False, True, False]


def test_dedup_and_rank_of_feature():
    a = pd.DataFrame({"arm": "S1-obs", "gene": ["g", "g", "g", "h"], "feature": ["AMP:1", "AMP:2", "HOT:3", "HOT:3"],
                      "variant_group": ["LG", "LG", None, None], "rank": [2.0, 1.0, 3.0, 4.0]})
    d = A.add_dedup(a)
    r = d.set_index(["gene", "feature"])
    assert np.isnan(r.at[("g", "AMP:1"), "rank_dedup"]) and r.at[("g", "AMP:1"), "lost_to_group_member"] == "AMP:2"
    assert r.at[("g", "AMP:2"), "rank_dedup"] == 1 and r.at[("h", "HOT:3"), "rank_dedup"] == 3
    # the group's best member never lends its rank: AMP:1 is scored by its own rank_stat
    assert A.rank_of_feature(d, "S1-obs", "g", "AMP:1") == 1.0      # only rows outside (g, LG) count
    assert A.rank_of_feature(d, "S1-obs", "h", "HOT:3") == 3.0


def test_permutation_keeps_strata_column_correlations_and_na_alignment():
    data = make_data(seed=6)
    st = NL.strata(data.cov_df["lineage"], data.cov_df["library_group"])
    perm = NL.permutation(st, 99)
    assert (st[perm] == st).all()                                        # strata counts and membership kept
    assert sorted(perm) == list(range(len(perm)))
    Xp, NAp = data.X[perm], data.na_flags[perm]
    C0 = pd.DataFrame(data.X).corr().to_numpy()
    C1 = pd.DataFrame(Xp).corr().to_numpy()
    assert np.allclose(np.nan_to_num(C0), np.nan_to_num(C1))           # joint row move keeps correlations
    amp = data.fmeta.index.get_loc(f"AMP:{data.genes[1]}")
    assert (np.isnan(Xp[:, amp]) == (NAp[:, 0] == 1)).all()              # NA flags move with 𝒟 rows
    assert not (perm == np.arange(len(perm))).all()
    # a stratum with a single line is merged into the lineage's largest stratum
    lin = np.array(["A"] * 5 + ["B"] * 4)
    lib = np.array(["x"] * 4 + ["y"] + ["x"] * 4)
    s = NL.strata(lin, lib)
    assert s[4] == s[0]


def test_null_draw_runs_with_observed_thresholds():
    data = make_data(seed=7, plant_beta=-1.5)
    res = A.run_arms(data, arms=A.ARMS, workers=1)
    out = NL.run_null(data, res, B=2, workers=1, genes=data.genes)
    dr = out["draws"]
    assert len(dr) == 2 and set(NL.LISTS) <= set(dr.columns)
    sm = out["summary"]
    assert set(sm["list"]) == set(NL.LISTS)
    assert (sm["fdr_emp_fixed"] >= 0).all()


def test_parallel_s1_equals_serial():
    data = make_data(seed=8)
    b = A.run_s1(data, A.S1_ARMS, workers=2, chunk=1)        # first: observed effects not yet loaded
    a = A.run_s1(data, A.S1_ARMS)
    key = ["arm", "gene", "feature"]
    a, b = a.sort_values(key).reset_index(drop=True), b.sort_values(key).reset_index(drop=True)
    pd.testing.assert_frame_equal(a, b)
    assert not a.loc[a["test_family"] == "F2", "selected"].any()          # S1 list = F1 (A6-iii)


def test_md_global_rule_and_outputs():
    en = pd.DataFrame({"arm": "EN-pred", "gene": ["g1", "g1", "g2", "g3", "g4"], "feature": ["a", "b", "a", "c", "d"],
                       "selected": True, "selection_freq": 0.9, "beta": -0.2, "maj_sign": [-1.0, 1.0, -1.0, -1.0, -1.0],
                       "stage": "G"})
    s1 = pd.DataFrame({"arm": "S1-obs", "gene": ["g1", "g1", "g2", "g3", "g4"], "feature": ["a", "b", "a", "c", "d"],
                       "beta": [-0.5, -0.5, -0.3, -0.6, -0.4], "se": 0.1, "df": 50.0, "p": [1e-6, 1e-6, 0.5, 1e-5, 1e-4],
                       "q": [1e-4, 1e-4, 0.6, 1e-3, 0.01], "test_family": ["F1", "F1", "F1", "F1", "F2"],
                       "selected": [True, True, False, False, False], "n_obs": 100, "n_carriers": 10, "dropped_cov": "",
                       "collinear_with_design": [False, False, False, True, False], "se_type": "CR2"})
    m = A.md_global(en, s1)
    assert list(m["pass"]) == [True, False, False, False, False]   # sign flip; not a BH hit; collinear; F2 excluded
    m2 = A.md_global(en, s1, p_cut=0.6)
    assert list(m2["pass"]) == [True, False, True, False, False]   # fixed-threshold null rule: F1 and p <= p*_F1
    data = make_data(seed=5, plant_beta=-1.5)
    res = A.run_arms(data, arms=A.ARMS, workers=1)
    mg = res["md_global"]
    a = res["arms"]
    en_sel = set(zip(*[a.loc[(a.arm == "EN-pred") & a.selected, c] for c in ("gene", "feature")]))
    s1_sel = set(zip(*[a.loc[(a.arm == "S1-obs") & a.selected, c] for c in ("gene", "feature")]))
    got = set(zip(mg["gene"], mg["feature"]))
    assert got <= (en_sel & s1_sel) and ("9100001", "HOT:9200001") in got
    assert "in_md_strict" in mg.columns
