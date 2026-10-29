"""Every verdict branch of P1 and P2 is reachable, the gate decides as specified, and the endpoint scripts refuse
without their preconditions (LESSONS 4; R3-4/R3-5). Synthetic data only; no effects are read."""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats as sps
from statsmodels.stats.proportion import proportion_confint

from biomark import arms as A
from biomark import endpoints as EP
from biomark import params as PR

ROOT = Path(__file__).resolve().parents[1]


def _script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


BR1 = PR.endpoints()["P1_primary"]["contrasts"]
BR2 = PR.endpoints()["P2_primary_clinical"]["primary_contrast"]["branches"]
KMIN = int(PR.endpoints()["P1_primary"]["min_k"])
NMIN = int(PR.endpoints()["P2_primary_clinical"]["min_testable_tests"])


# ------------------------------------------------------------------------------------------------
# exact statistics
# ------------------------------------------------------------------------------------------------


def test_exact_statistics_match_references():
    for b, c in ((0, 0), (3, 9), (12, 2), (7, 7)):
        ref = 1.0 if b + c == 0 else sps.binomtest(b, b + c, 0.5).pvalue
        assert EP.mcnemar_exact(b, c) == pytest.approx(ref, rel=1e-12)
    for k, n in ((0, 20), (7, 20), (20, 20)):
        lo, hi = EP.clopper_pearson(k, n)
        rl, rh = proportion_confint(k, n, alpha=0.05, method="beta")
        assert (lo, hi) == pytest.approx((rl, rh), abs=1e-12)


# ------------------------------------------------------------------------------------------------
# P1 branches
# ------------------------------------------------------------------------------------------------


@pytest.mark.parametrize("contrast", list(BR1))
def test_p1_every_branch_reachable(contrast):
    br = BR1[contrast]["branches"]
    a1, a2 = [x.strip() for x in BR1[contrast]["diff"].split(" - ")]
    rng = np.random.default_rng(0)
    got = set()
    for shift, want in ((0.3, br["lower > 0"]), (-0.3, br["upper < 0"]), (0.0, br["else"])):
        boot = pd.DataFrame({f"G|{a1}": 0.5 + shift + rng.normal(0, 0.02, 400), f"G|{a2}": 0.5 + rng.normal(0, 0.02, 400)})
        point = {k: float(v.mean()) for k, v in boot.items()}
        c = [x for x in EP.p1_contrasts(point, boot, {k: KMIN + 5 for k in boot}, prefix="G|") if x["contrast"] == contrast][0]
        assert c["verdict"] == want
        got.add(c["verdict"])
        small = [x for x in EP.p1_contrasts(point, boot, {k: KMIN - 1 for k in boot}, prefix="G|")
                 if x["contrast"] == contrast][0]
        assert small["verdict"] == EP.descriptive(KMIN - 1)
    assert len(got) == 3


def _synthetic_R(n_pat=150, seed=1):
    """R-like line set: two-line patients, lineages, a planted determinant (feature 'HOT:1' on gene 'g1')."""
    rng = np.random.default_rng(seed)
    pats = np.repeat(np.arange(n_pat), rng.choice([1, 2], size=n_pat, p=[0.8, 0.2]))
    n = len(pats)
    ids = [f"R{i:04d}" for i in range(n)]
    cov = pd.DataFrame({"lineage": rng.choice(["L1", "L2", "L3"], n), "library_group": rng.choice(["Avana-only", "KY-only"], n),
                        "cellformat": "adherent", "serumfree": "F", "medium_family": rng.choice(["RPMI", "DMEM"], n),
                        "doubling_time": rng.normal(40, 5, n), "DT_unknown": 0.0, "MSI": "MSS",
                        "DMG_BURDEN": rng.normal(5, 0.5, n), "tp53_status": rng.choice(["mut", "WT"], n),
                        "PatientID": [f"P{p}" for p in pats]}, index=ids)
    X = (rng.random((n, 6)) < 0.15).astype(float)
    feats = [f"HOT:{j}" for j in range(1, 7)]
    Y = rng.normal(0, 0.4, (n, 3))
    Y[:, 0] += -1.0 * X[:, 0]
    Y[:, 1] += -0.9 * X[:, 1]
    Yd = pd.DataFrame(Y, index=ids, columns=["g1", "g2", "g3"])
    return EP.RFitter(X, pd.Series(np.arange(6), index=feats), cov, None, None, Yd, cov["PatientID"].to_numpy()), feats


def test_p1_end_to_end_bootstrap_reaches_predictor_adds_value():
    fit, feats = _synthetic_R()
    lists = {"G|EN-pred": [("g1", "HOT:1", -1.0), ("g2", "HOT:2", -1.0)],
             "G|EN-obs": [("g1", "HOT:5", -1.0), ("g3", "HOT:6", -1.0)],
             "G|S1-obs": [("g1", "HOT:1", -1.0), ("g3", "HOT:4", -1.0)]}
    binary = {f: True for f in feats}
    point = {nm: float(np.mean([EP.replicated(fit.fit(g, f, "", sign=s), s) for g, f, s in L])) for nm, L in lists.items()}
    assert point["G|EN-pred"] == 1.0 and point["G|EN-obs"] < 1.0
    bt = EP.bootstrap_rates(fit, lists, {}, binary, B=60, seed=3)
    cs = {c["contrast"]: c for c in EP.p1_contrasts(point, bt["rates"], {k: KMIN for k in lists}, prefix="G|")}
    assert cs["predictor_adds_value"]["verdict"] == BR1["predictor_adds_value"]["branches"]["lower > 0"]
    assert set(bt["frac_degenerate"].columns) == set(lists)


def test_replicated_rules_and_degenerate_pairs():
    fit, _ = _synthetic_R()
    r = fit.fit("g1", "HOT:1", "", sign=-1.0)
    assert EP.replicated(r, -1.0) and not EP.replicated(r, 1.0)          # same sign required
    one = r.copy()
    one["n_carrier_patients"] = 1
    assert not EP.replicated(one, -1.0)                                   # < 2 distinct carrier patients
    assert r["p_one"] == pytest.approx(sps.t.cdf(r["beta"] / r["se"], r["df"]))


def test_expected_nonnan_library_rule():
    D = [f"D{i}" for i in range(40)]
    lib_D = pd.Series(["Avana"] * 20 + ["KY"] * 20, index=D)
    nan = pd.DataFrame(False, index=D, columns=["a", "b"])
    nan.loc[D[20:], "b"] = True                       # gene b missing in KY
    nan.loc[D[0], "a"] = True                         # 1/20 Avana NaN -> 95 % non-NaN: still expected
    lib_R = pd.Series(["Avana", "KY", "Avana+KY"], index=["r1", "r2", "r3"])
    E = EP.expected_nonnan(nan, lib_D, lib_R)
    assert E.loc["r1"].tolist() == [True, True]
    assert E.loc["r2"].tolist() == [True, False]
    assert E.loc["r3"].tolist() == [True, False]      # unseen string: union NaN set of its components


# ------------------------------------------------------------------------------------------------
# P2 branches
# ------------------------------------------------------------------------------------------------


def _scores(n, en_hits, s1_hits):
    rows = []
    for i in range(n):
        for arm, hits in (("EN-pred", en_hits), ("S1-obs", s1_hits)):
            h = i in hits
            rows.append({"arm": arm, "test_id": f"t{i}", "hit": h, "testable": True, "sign": "S",
                         "sign_correct": True if h else np.nan})
    return pd.DataFrame(rows)


def test_p2_every_branch_reachable():
    keys = list(BR2)
    n = NMIN + 5
    s = EP.p2_summary(_scores(n, set(range(12)), set()))
    assert s["verdict"] == BR2[keys[0]] and s["mcnemar"]["p"] < 0.05
    s = EP.p2_summary(_scores(n, set(), set(range(12))))
    assert s["verdict"] == BR2[keys[1]]
    s = EP.p2_summary(_scores(n, {0, 1, 2}, {1, 2, 3}))
    assert s["verdict"] == BR2[keys[2]]
    s = EP.p2_summary(_scores(NMIN - 1, set(range(10)), set()))
    assert s["verdict"] == EP.descriptive(NMIN - 1) and s["p2_status"] == "descriptive"


def test_p2_headline_set_and_hit_rule():
    cgt = pd.DataFrame({
        "test_id": ["a", "b", "c", "d", "e", "f"], "feature_D": ["HOT:1", "HOT:2", "HOT:3", "HOT:1", "HOT:4", "HOT:5"],
        "target_set": ["g1", "g2;g3", "g1", "g1", "g1", "g1"], "sign": ["S", "S", "W", "S", "S", "N"],
        "partition": ["sensitivity", "sensitivity", "bypass-resistance", "sensitivity", "sensitivity", "expected-null"],
        "ko_class": ["concordant", "concordant-family", "concordant", "concordant", "concordant", "concordant"],
        "tier": ["T1", "T3", "T2", "preclinical", "T1", "T1"], "proxy": [0, 0, 0, 0, 1, 0],
        "headline_exclude": [np.nan] * 6, "in_gate": [1, 0, 0, 0, 0, 0]})
    U = pd.DataFrame({"gene": ["g1", "g3", "g1", "g1"], "feature": ["HOT:1", "HOT:2", "HOT:3", "HOT:4"],
                      "n_carriers": [10, 6, 4, 9]})
    T = EP.p2_test_set(cgt, U)
    assert list(T["test_id"]) == ["a", "b", "c"]                       # preclinical (A4), proxy, N partition out
    assert T.set_index("test_id")["testable"].to_dict() == {"a": True, "b": True, "c": False}
    assert list(EP.p2_test_set(cgt, U, gate_disjoint=True)["test_id"]) == ["b", "c"]
    arms_tab = pd.DataFrame({"arm": "EN-pred", "gene": ["g1", "g3", "g9"], "feature": ["HOT:1", "HOT:2", "HOT:9"],
                             "rank": [3.0, 1.0, 2.0], "rank_dedup": [3.0, 1.0, 2.0], "variant_group": [None] * 3,
                             "beta": [-0.5, 0.4, -0.1]})
    sc = EP.p2_score(arms_tab, T, k=2).set_index("test_id")
    assert sc.at["a", "rank_of_feature_D"] == 3 and not sc.at["a", "hit"]
    assert sc.at["b", "hit"] and sc.at["b", "gene"] == "g3" and sc.at["b", "sign_correct"] == False  # noqa: E712
    assert np.isnan(sc.at["c", "rank_of_feature_D"])


def test_rank_index_equals_arms_rank_of_feature():
    rng = np.random.default_rng(4)
    n = 60
    a = pd.DataFrame({"arm": "S1-obs", "gene": rng.choice(["g1", "g2", "g3"], n),
                      "feature": [f"F:{i}" for i in range(n)], "rank_stat": rng.random(n),
                      "variant_group": rng.choice([None, "LG1", "LG2"], n), "beta": rng.normal(size=n)})
    a["rank"] = a["rank_stat"].rank(method="first")
    a.loc[a.index[::11], "rank"] = np.nan
    a = A.add_dedup(a)
    idx = EP.RankIndex(a)
    for g, f in zip(a["gene"], a["feature"]):
        r1 = A.rank_of_feature(a, "S1-obs", g, f)
        r2, _ = idx.rank_of("S1-obs", g, f)
        assert (np.isnan(r1) and np.isnan(r2)) or r1 == r2


def test_actionability_and_population_tables():
    cgt = pd.DataFrame({"test_id": ["a"], "feature_D": ["HOT:1"], "target_set": ["g1"]})
    md = pd.DataFrame({"gene": ["g1", "g2", "g1"], "feature": ["HOT:1", "HOT:1", "PROG:3"]})
    act = EP.actionability(md, md.iloc[:1], cgt, n_perm=50, seed=1).set_index("list")
    assert act.at["MD-strict", "pair_level"] == 0.5 and act.at["MD-strict", "feature_level"] == 1.0
    fit, _ = _synthetic_R()
    tests = pd.DataFrame({"test_id": ["a"], "feature_D": ["HOT:1"], "targets": [["g1"]], "sign": ["S"],
                          "testable": [True], "indication_codes": ["C1"], "tumour_agnostic": [0]})
    lineage = np.array(["L1", "L2", "L3"] * (len(fit.patients) // 3 + 1))[: len(fit.patients)]
    code = np.where(lineage == "L1", "C1", "C9")
    per, lab = EP.population_by_lineage(fit, tests, {}, lineage, code)
    assert lab["label"].iloc[0] in {"on-label", "beyond-label", "off-label-only", "indication-untestable", "not supported"}
    assert lab["label"].iloc[0] == "beyond-label"                       # planted pan-lineage effect
    assert (per["scope"] == "within-indication").sum() == 1


# ------------------------------------------------------------------------------------------------
# Gate
# ------------------------------------------------------------------------------------------------


def _gate_inputs(n_fam=4, rec=None):
    ids = [f"C{i}" for i in range(2 * n_fam)]
    pc = pd.DataFrame({"control_id": ids, "family": [f"F{i // 2}" for i in range(2 * n_fam)], "gate_eligible": "True",
                       "access_class": "proxy", "en_max_attainable_freq": "1.0"})
    ctl = pd.DataFrame({"control_id": ids, "dependency_entrez": [f"g{i}" for i in range(2 * n_fam)],
                        "feature_primary_ids": [f"HOT:{i}" for i in range(2 * n_fam)], "feature_secondary_ids": "",
                        "fallback_feature_ids": "", "gate_feature": "primary", "expected_sign": "S"})
    rec = set(range(2 * n_fam)) if rec is None else rec
    md = pd.DataFrame({"gene": [f"g{i}" for i in rec], "feature": [f"HOT:{i}" for i in rec], "beta": -1.0})
    return pc, ctl, {"MD-strict": md}


def test_gate_pass_and_fail_branches():
    pc, ctl, L = _gate_inputs()
    assert EP.gate_attempt(pc, ctl, L)["summary"]["pass"]
    pc, ctl, L = _gate_inputs(rec={0, 2, 4, 6, 7})                     # overall 5/8 >= 0.6 but every family 1/2 = 0.5
    s = EP.gate_attempt(pc, ctl, L)["summary"]
    assert s["pass"]
    pc, ctl, L = _gate_inputs(rec={0, 1, 2, 3, 4, 5})                  # 6/8 overall, family F3 0/2 -> fail
    s = EP.gate_attempt(pc, ctl, L)["summary"]
    assert not s["pass"] and s["family_pass"]["F3"] is False
    pc, ctl, L = _gate_inputs(rec={0, 2, 4})                            # 3/8 < 0.6 -> fail
    assert not EP.gate_attempt(pc, ctl, L)["summary"]["pass"]
    L["MD-strict"]["beta"] = 1.0                                        # wrong sign is not recovered
    assert EP.gate_attempt(pc, ctl, L)["summary"]["overall"]["k"] == 0


# ------------------------------------------------------------------------------------------------
# Refusals
# ------------------------------------------------------------------------------------------------


def test_second_r_open_refused_and_reproduction_rules():
    empty = pd.DataFrame(columns=["utc_time", "commit", "purpose", "scope"])
    EP.require_single_r_read(empty)
    log = pd.DataFrame({"utc_time": ["t"], "commit": ["abc"], "purpose": ["replicate.py single R run cfg=XYZ"],
                        "scope": ["replication"]})
    with pytest.raises(EP.EndpointRefusal):
        EP.require_single_r_read(log)
    EP.require_single_r_read(log, {"commit": "abc", "config_sha256": "XYZ"})
    with pytest.raises(EP.EndpointRefusal):
        EP.require_single_r_read(log, {"commit": "other", "config_sha256": "XYZ"})
    with pytest.raises(EP.EndpointRefusal):
        EP.require_single_r_read(log, {"commit": "abc", "config_sha256": "NEW"})
    limits_only = log.assign(scope="limits")
    EP.require_single_r_read(limits_only)                               # set-A access is not an R read


def test_scripts_refuse_before_their_tags(monkeypatch):
    rep = _script("replicate")
    monkeypatch.setattr(EP, "descends_from", lambda tag, repo=None: False)
    with pytest.raises(EP.EndpointRefusal, match="freeze-catalogue"):
        rep.main([])                                                    # pre-tag replication refused
    gate = _script("run_gate")
    with pytest.raises(EP.EndpointRefusal, match="freeze-bounds"):
        gate.run(None, None)                                            # gate before freeze-bounds refused
    pool = _script("build_p1_pool")
    with pytest.raises(EP.EndpointRefusal, match="freeze-bounds"):
        pool.main()


def test_replicate_refuses_second_open_even_after_tag(monkeypatch):
    rep = _script("replicate")
    monkeypatch.setattr(EP, "descends_from", lambda tag, repo=None: True)
    from biomark import guard
    monkeypatch.setattr(guard, "read_log", lambda: pd.DataFrame(
        {"utc_time": ["t"], "commit": ["abc"], "purpose": ["replicate.py single R run cfg=x"], "scope": ["replication"]}))
    called = []
    monkeypatch.setattr(guard, "load_effects", lambda *a, **k: called.append(1))
    with pytest.raises(EP.EndpointRefusal, match="second R open"):
        rep.main([])
    assert not called                                                   # refused before any effect read


def test_catalogue_hash_mismatch_detected(tmp_path):
    f = tmp_path / "a.csv"
    f.write_text("x\n1\n")
    man = EP.catalogue_manifest(["a.csv"], root=tmp_path)
    assert EP.catalogue_mismatches(man, root=tmp_path) == []
    f.write_text("x\n2\n")
    assert EP.catalogue_mismatches(man, root=tmp_path) == ["a.csv"]


def test_pool_power_and_frozen_lists():
    fit, feats = _synthetic_R()
    n = fit.X.shape[0]
    genes = ["g1", "g2", "g3"]
    UG = pd.DataFrame([(g, f) for g in genes for f in feats], columns=["gene", "feature"])
    key = pd.MultiIndex.from_arrays([UG["gene"], UG["feature"]])
    bsh = pd.Series(np.where(UG["feature"] == "HOT:1", -1.0, -0.02), index=key)
    sig = pd.Series(0.4, index=genes)
    expR = pd.DataFrame(True, index=range(n), columns=genes)
    expR.iloc[: n // 2, 2] = False                                        # g3 expected NaN on half of R
    drops = pd.Series("", index=key)
    pool = EP.pool_power(UG, bsh, sig, fit.X, fit.fpos, pd.Series(True, index=feats), fit.cov, fit.patients, expR, drops)
    assert pool["R_evaluable"].all()
    assert (pool.loc[pool["feature"] == "HOT:1", "power"] > 0.9).all()
    assert (pool.loc[pool["feature"] != "HOT:1", "power"] < 0.5).all()
    assert (pool.loc[pool["gene"] == "g3", "n_R"] < pool.loc[pool["gene"] == "g1", "n_R"].max()).all()
    # powered pool + lists: dedup within (g, variant group), unrankable rows excluded, k rule
    arms_tab = pd.DataFrame({"arm": "EN-pred", "gene": UG["gene"], "feature": UG["feature"],
                             "rank": np.arange(1, len(UG) + 1, dtype=float), "rank_stat": 1.0,
                             "variant_group": np.where(UG["feature"].isin(["HOT:1", "HOT:2"]), "VG", None),
                             "beta": -0.1, "maj_sign": -1.0})
    arms_tab.loc[arms_tab["feature"] == "HOT:3", "rank"] = np.nan
    L = EP.build_lists(arms_tab, UG, ["EN-pred"], k=5)["EN-pred"]
    got = [(g, f) for g, f, _ in L["rows"]]
    assert ("g1", "HOT:2") not in got and ("g1", "HOT:1") in got          # group member de-duplicated
    assert all(f != "HOT:3" for _, f in got)                              # unrankable never enters
    assert L["n_list"] == 5 and EP.p1_k(len(UG)) == min(100, len(UG) // 3)


def test_rank_index_slim_equals_full():
    rng = np.random.default_rng(5)
    n = 300
    a = pd.DataFrame({"arm": rng.choice(["S1-obs", "EN-pred"], n), "gene": rng.choice([f"g{i}" for i in range(8)], n),
                      "feature": [f"F:{i}" for i in range(n)], "rank_stat": rng.random(n),
                      "variant_group": rng.choice([None, "LG1", "LG2", "LG3"], n), "beta": rng.normal(size=n)})
    a["rank"] = a.groupby("arm")["rank_stat"].rank(method="first")
    a.loc[a.index[::13], "rank"] = np.nan
    a = A.add_dedup(a)
    full = EP.RankIndex(a)
    qg = ["g1", "g4", "g6"]
    dedup = {arm: s.loc[s["rank_dedup"].notna(), "rank"].to_numpy() for arm, s in a.groupby("arm")}
    slim = EP.RankIndex.slim(dedup, a[a["gene"].isin(qg)])
    for arm in ("S1-obs", "EN-pred"):
        for g, f in zip(a.loc[a["gene"].isin(qg), "gene"], a.loc[a["gene"].isin(qg), "feature"]):
            r1, b1 = full.rank_of(arm, g, f)
            r2, b2 = slim.rank_of(arm, g, f)
            assert (np.isnan(r1) and np.isnan(r2)) or r1 == r2
