"""Post-R D-only completion scripts (DEV-POST-R-COMPLETION): ceiling labels follow the numbers (§10.4), the matched-SD
low-tier draw (§7.3 NC2), the miss decomposition order and NC5 (§10.5, §7.3), NC2/NC3 statistics and pass/fail
labels against the frozen bounds, R power merging (§6.1), and every script refuses before its tag. Synthetic data."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from biomark import endpoints as EP
from biomark import params as PR

ROOT = Path(__file__).resolve().parents[1]


def _script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def ce():
    return _script("run_ceiling")


def test_ceiling_labels_follow_numbers(ce):
    q = float(PR.dp()["clinical_s10"]["ceiling_bh_q"])
    assert ce.label_of(0.01, 0.99, q) == "data-supported"
    assert ce.label_of(0.5, 0.01, q) == "data-contradicted"
    assert ce.label_of(0.5, 0.5, q) == "data-null" and ce.label_of(np.nan, np.nan, q) == "data-null"
    p = np.array([0.001, 0.02, 0.04, 0.5, np.nan])
    qs = ce.bh(p)
    assert np.isnan(qs[-1]) and np.allclose(qs[:4], [0.004, 0.04, 0.0533333, 0.5], atol=1e-6)
    assert ce.family([0.01, 0.2], 2) == 0.02 and ce.family([0.01], 1) == 0.01 and np.isnan(ce.family([np.nan], 2))
    t = pd.DataFrame({"test_id": ["a", "b", "c"], "ceiling_testable": [True, True, False], "q": [0.01, 0.5, np.nan],
                      "q_opposite": [0.9, 0.5, np.nan], "label": ["data-supported", "data-null", "untestable"]})
    ce.check_labels(t, q)
    with pytest.raises(AssertionError):
        ce.check_labels(t.assign(label=["data-null", "data-null", "untestable"]), q)
    rng = np.random.default_rng(0)
    car, non = rng.normal(-1, 1, 20), rng.normal(0, 1, 60)
    assert ce.wilcoxon_one_sided(car, non, -1.0) < 0.01 < ce.wilcoxon_one_sided(car, non, 1.0)


def test_lowtier_matched_sd_draw():
    lt = _script("build_lowtier_genes")
    el = pd.DataFrame({"entrez": [str(i) for i in range(1, 41)], "eligible": [i <= 10 for i in range(1, 41)],
                       "tier": ["high"] * 10 + ["low"] * 30})
    sd = pd.Series(np.r_[np.linspace(0.2, 0.6, 10), np.linspace(0.1, 0.7, 30)], index=el["entrez"])
    a, b = lt.match(el, sd, 5, 7), lt.match(el, sd, 5, 7)
    pd.testing.assert_frame_equal(a, b)                                           # deterministic
    assert len(a) == 5 and a["entrez"].is_unique and not set(a["entrez"]) & set(el["entrez"][:10])
    for _, r in a.iterrows():                                                     # nearest SD among tier-low genes
        assert abs(r["sd"] - r["reference_sd"]) <= 0.6 / 29 + 1e-12


def test_miss_decomposition_order_and_nc5():
    cs = _script("run_cgt_secondary")
    cgt = pd.DataFrame({"test_id": ["ineligible", "nofeat", "ceil", "lost", "below", "hit"],
                        "targets": [["g0"], ["g1"], ["g1"], ["g1"], ["g1"], ["g1"]],
                        "feature_D": ["F:1", None, "F:1", "F:2", "F:1", "F:1"]})
    sc = pd.DataFrame({"arm": "S1-obs", "test_id": cgt["test_id"], "gene": "g1", "rank_of_feature_D": 500.0,
                       "hit": [False] * 5 + [True], "testable": True, "set": "non_proxy", "subset": "all", "k": 200})
    md = cs.decompose(sc, cgt, eligible={"g1"}, gene_class={"g0": "2"}, U={("g1", "F:1"), ("g1", "F:2")},
                      ceil={"ceil": "data-null"}, lost={("S1-obs", "g1", "F:2"): "F:3"})
    got = dict(zip(md["test_id"], md["category"]))
    C = cs.CATS
    assert got == {"ineligible": C[0], "nofeat": C[1], "ceil": C[4], "lost": C[3], "below": C[2]}
    assert "hit" not in got and md.loc[md["test_id"] == "ineligible", "detail"].iloc[0] == "classes 2"
    s5 = pd.DataFrame({"arm": ["S1-obs"] * 3 + ["EN-pred"] * 2, "test_id": ["a", "b", "c", "a", "b"],
                       "rank_of_feature_D": [10.0, 300.0, np.nan, 10.0, 50.0], "beta": [-1.0, -1.0, -1.0, 1.0, -0.5]})
    n5 = cs.nc5(s5, 200)
    assert n5["per_arm"]["S1-obs"] == {"n_ranked": 2, "false_sensitivity": 1, "rate": 0.5, "tests": ["a"]}
    assert n5["per_arm"]["EN-pred"]["rate"] == 0.5 and n5["value"] == 0.5
    flags = cs.testable_flags(pd.DataFrame({"targets": [["g1"], ["g1"], ["g2"]], "feature_D": ["F:1", "F:9", "F:1"]}),
                              {("g1", "F:1"): 7, ("g1", "F:9"): 2}, 5)
    assert flags == [True, False, False]


def test_nc_statistics_and_labels(tmp_path):
    nc = _script("run_nc_real")
    genes = ["1", "2", "3", "4"]
    pd.DataFrame({"arm": "EN-pred", "gene": genes, "r2_oof_beyond_design": [0.0, 0.001, 0.003, np.nan]}).to_csv(
        tmp_path / "lambda.csv", index=False)
    pd.DataFrame({"gene": ["1"], "feature": ["F:1"]}).to_csv(tmp_path / "md_strict.csv", index=False)
    pd.DataFrame({"gene": [], "feature": []}).to_csv(tmp_path / "md_global.csv", index=False)
    n2 = nc.nc2(tmp_path, tmp_path, genes)
    assert n2["NC2a"]["value"] == 0.001 and n2["NC2b"]["value"] == 0.25 and n2["NC2b_md_global"]["value"] == 0.0
    n3 = nc.nc3(pd.DataFrame({"testable": [True, True, False], "frac_dependent": [1.0, 0.0, 1.0]}))
    assert n3 == {"value": 0.5, "n_testable_pairs": 2, "n_dependent": 1, "n_pairs": 3}
    b = nc.bounds()
    sec = {"NC5": {"value": 0.1, "per_arm": {}, "n_tests": 3}}
    null = {"fdr_emp_fixed": {"S1-obs": 0.02, "EN-obs selected": 0.2}}
    lim = {"NC6_lineage_tf_class3": {"frac_class3": 0.0, "n": 18}}
    out = nc.assemble(n2, n3, sec, null, lim, b)
    assert out["NC2a"]["pass"] is True and out["NC2b"]["pass"] is True and out["NC3"]["pass"] is False
    assert out["NC3"]["label"] == "NC3 failed (0.5)" and out["NC5"]["pass"] is True and out["NC6"]["pass"] is False
    assert out["NC1"]["EN-obs selected"]["pass"] is False and "EN-obs selected" in out["NC1"]["EN-obs selected"]["label"]


def test_r_power_merge_keeps_frozen_and_checks_identity():
    rp = _script("run_r_power")
    rt = pd.DataFrame({"gene": ["1", "2"], "feature": ["HOT:1", "PROG:1"], "power": [0.8, np.nan]})
    pw = pd.DataFrame({"gene": ["1", "2"], "feature": ["HOT:1", "PROG:1"], "power": [0.8, 0.3],
                       "R_evaluable": [True, False], "n_R": 1, "n1_R": 1, "n0_R": 1, "sd_R": 1.0, "nu_R": 5.0,
                       "se_R_pred": 0.1, "beta_shrunk": -0.2, "binary": [True, False]})
    m = rp.merge_power(rt, pw).set_index("feature")
    assert m.at["HOT:1", "power"] == 0.8 and m.at["HOT:1", "power_source"].startswith("frozen pool")
    assert m.at["PROG:1", "power"] == 0.3 and m.at["PROG:1", "power_effective"] == 0.0     # not R-evaluable
    with pytest.raises(AssertionError, match="differ from the frozen pool"):
        rp.merge_power(rt.assign(power=[0.7, np.nan]), pw)


@pytest.mark.parametrize("name, tag", [("run_ceiling", "freeze-bounds"), ("run_cgt_secondary", "freeze-bounds"),
                                       ("run_nc_real", "freeze-bounds"), ("run_r_power", "freeze-catalogue")])
def test_scripts_refuse_before_their_tag(name, tag):
    with pytest.raises(EP.EndpointRefusal, match=tag):
        _script(name).run(check=lambda t, repo=None: False)
