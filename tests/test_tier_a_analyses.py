"""XLIB split / discovery eligibility, component readout, VIF, Q4c precision / OR, OncoTree mapping, registry and
timing stems (synthetic data only; scripts checked for refusal before freeze-bounds)."""
import ast
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats as sps

from biomark import analyses as AN
from biomark import arms as A
from biomark import endpoints as EP
from biomark import params as PR
from biomark import stats as ST

from synth_data import make_data

ROOT = Path(__file__).resolve().parents[1]


def _script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_xlib_split_by_library_string():
    lib = pd.Series(["Avana", "KY", "Avana+KY", "Humagne-CD", "", "Brunello+KY", "Avana+Humagne-CD"])
    s = AN.xlib_split(lib).tolist()
    assert s == ["discovery", "replication", "", "replication", "", "replication", ""]


def test_discovery_eligibility_recomputes_r_only():
    rng = np.random.default_rng(0)
    n, G = 200, 3
    obs = rng.normal(size=(n, G))
    pred = np.column_stack([obs[:, 0] + rng.normal(0, 0.3, n), rng.normal(size=n), obs[:, 2] + rng.normal(0, 0.3, n)])
    strata = np.where(rng.random(n) < 0.5, "L1|Avana-only", "L2|Avana-only")
    el = pd.DataFrame({"entrez": ["a", "b", "c"], "class_ok": [True, True, False], "tier_ok": [True, True, True]})
    ok = AN.discovery_eligibility(el, pred, obs, strata, ["a", "b", "c"])
    assert ok.tolist() == [True, False, False]          # b: no r; c: class criterion kept


def test_component_selection_layers_and_carrier_input():
    md = pd.DataFrame({"gene": ["g", "g", "g", "h"], "feature": ["E:1", "M:2", "B:3", "N:4"]})

    def tab(nz):
        return pd.DataFrame({"gene": ["g", "g", "g", "h"], "feature": ["E:1", "M:2", "B:3", "N:4"],
                             "beta_full": nz, "selection_freq": np.abs(nz), "selected": np.array(nz) != 0})
    comp = {"M1": tab([0.2, 0.0, 0.1, 0.0]), "M2n-M1": tab([0.0, -0.3, 0.0, 0.0]), "M3n-M1": tab([0.0, 0.0, 0.2, 0.0])}
    acc = pd.Series(["expression-derived", "proxy", "proxy", "proxy"],
                    index=pd.MultiIndex.from_arrays([md["gene"], md["feature"]]))
    s = AN.component_selection(md, comp, acc).set_index("feature")
    assert s.loc["E:1", "layer"] == "expression-mediated" and s.loc["E:1", "carrier_input"] == ""
    assert s.loc["M:2", "layer"] == "genomics-layer" and "M2 local" in s.loc["M:2", "carrier_input"]
    assert s.loc["B:3", "layer"] == "both" and "M3 pool" in s.loc["B:3", "carrier_input"]
    assert s.loc["N:4", "layer"] == "none"


def test_component_readout_uses_the_component_target_with_the_unchanged_engine():
    d = make_data(seed=11, plant_beta=-1.5)
    U = A.universe(d)
    base, _, _ = A.run_en(d, ("EN-pred",), U=U)
    flat = AN.data_with_pred(d, np.zeros_like(d.pred))           # a component with no signal
    t0, _, _ = A.run_en(flat, ("EN-pred",), U=U)
    assert t0["selected"].sum() == 0 and base["selected"].sum() > 0
    assert d.pred.any()                                            # original Data untouched


def test_design_vif_flags_collinearity():
    rng = np.random.default_rng(1)
    n = 300
    cov = pd.DataFrame({"lineage": rng.choice(["A", "B", "C"], n), "library_group": rng.choice(["x", "y"], n),
                        "doubling_time": rng.normal(size=n), "DT_unknown": 0.0, "DMG_BURDEN": rng.normal(size=n),
                        "MSI": "MSS", "tp53_status": rng.choice(["mut", "WT"], n)})
    v, cond, al = AN.design_vif(ST.CovDesign(cov))
    assert not al and v.max() < PR.dp()["confounders_s5"]["vif_max"] and np.isfinite(cond)
    cov["DMG_BURDEN"] = cov["doubling_time"] * 2 + rng.normal(0, 0.01, n)
    v2, cond2, _ = AN.design_vif(ST.CovDesign(cov))
    assert v2.max() > PR.dp()["confounders_s5"]["vif_max"] and cond2 > cond
    cov["serumfree"] = np.where(np.arange(n) < 10, "unknown", "F")       # aliased culture dummies (as on real D)
    cov["cellformat"] = np.where(np.arange(n) < 10, "x", "adherent")
    _, _, al3 = AN.design_vif(ST.CovDesign(cov))
    assert len(al3) == 1


def test_population_precision_and_carriers_only_or():
    x = np.array([1] * 20 + [0] * 10, float)
    ind = np.array([True] * 10 + [False] * 10 + [True] * 10)
    acc = np.array([1] * 8 + [0] * 2 + [1] * 3 + [0] * 7 + [1] * 10, float)
    r = AN.population_precision(x, acc, ind)
    assert r["carriers_indicated"] == 10 and r["carriers_not_indicated"] == 10
    assert r["precision_indicated"] == 0.8 and r["precision_not_indicated"] == 0.3
    ref = sps.contingency.odds_ratio(np.array([[8, 2], [3, 7]]))
    assert r["or_carriers"] == pytest.approx(ref.statistic)
    assert r["or_fisher_p"] == pytest.approx(sps.fisher_exact([[8, 2], [3, 7]])[1])
    r0 = AN.population_precision(x, acc, np.ones(30, bool))
    assert np.isnan(r0["or_carriers"])                              # no non-indicated carriers


def test_indicated_lines_parent_codes_from_terms():
    model = pd.DataFrame({"OncotreeCode": ["LUAD", "LUSC", "COAD", "SKCM", "XX"],
                          "OncotreeSubtype": ["Lung Adenocarcinoma", "Lung Squamous", "Colon Adenocarcinoma", "Cutaneous Melanoma", ""],
                          "OncotreePrimaryDisease": ["Non-Small Cell Lung Cancer", "Non-Small Cell Lung Cancer",
                                                     "Colorectal Adenocarcinoma", "Melanoma", ""],
                          "OncotreeLineage": ["Lung", "Lung", "Bowel", "Skin", ""]})
    terms = {"NSCLC": {"Non-Small Cell Lung Cancer"}, "MEL": {"Melanoma"}}
    m, um = AN.indicated_lines(["NSCLC", "COAD", "SOLID"], model, terms)
    assert m.tolist() == [True, True, True, False, False]
    assert um == ["SOLID"]


def test_registry_and_timing_writers(tmp_path):
    reg = tmp_path / "registry.csv"
    AN.append_registry("scripts/x.py", "ok", ["a.csv", "b.json"], registry=reg)
    AN.append_registry("scripts/y.py", "dryrun", [], registry=reg)
    r = pd.read_csv(reg)
    assert list(r.columns) == AN.REGISTRY_COLS and r["status"].tolist() == ["ok", "dryrun"]
    AN.write_timing("vif", {"wall_s": 1.5}, timing_dir=tmp_path, workers=1)
    assert (tmp_path / "vif.json").exists()


def test_every_timing_stem_written_by_scripts_is_mapped():
    mapped = {s for v in PR.tiers()["stage_job_map"].values() for s in v}
    found = set()
    for f in (ROOT / "scripts").glob("*.py"):
        for n in ast.walk(ast.parse(f.read_text(encoding="utf-8"))):
            if isinstance(n, ast.Call) and getattr(n.func, "attr", getattr(n.func, "id", "")) == "write_timing":
                if n.args and isinstance(n.args[0], ast.Constant):
                    kw = {k.arg for k in n.keywords}
                    if "timing_dir" not in kw:
                        found.add(n.args[0].value)
    assert found and found <= mapped, found - mapped


@pytest.mark.parametrize("name", ["run_xlib", "run_component_readout", "run_confounder_audit", "run_q4"])
def test_scripts_refuse_before_freeze_bounds(name, monkeypatch):
    mod = _script(name)
    monkeypatch.setattr(EP, "descends_from", lambda tag, repo=None: False)
    with pytest.raises(EP.EndpointRefusal, match="freeze-bounds"):
        mod.main([]) if name != "run_q4" else mod.main()


def test_timing_schema_fields(tmp_path):
    import datetime as dt
    import json
    s = AN.utc_now() - dt.timedelta(seconds=30)
    AN.write_timing("vif", {"note": "x"}, timing_dir=tmp_path, start_utc=s, workers=3)
    r = json.loads((tmp_path / "vif.json").read_text(encoding="utf-8"))
    assert r["job"] == "vif" and r["workers"] == 3 and isinstance(r["wall_s"], float) and r["wall_s"] >= 29
    assert r["start_utc"] < r["utc"] and r["run_id"]
    AN.write_timing("gate_pc", {}, timing_dir=tmp_path, wall_s=1.0, workers=1)
    r2 = json.loads((tmp_path / "gate_pc.json").read_text(encoding="utf-8"))
    assert r2["workers"] == 1 and r2["run_id"] == r["run_id"]          # one run id per process
    with pytest.raises(ValueError):
        AN.write_timing("vif", {}, timing_dir=tmp_path, workers=1)     # numeric wall_s required
    for bad in (None, 0, -1, 1.5, "2", True):                          # workers: no silent default
        with pytest.raises(ValueError):
            AN.write_timing("vif", {}, timing_dir=tmp_path, wall_s=1.0, workers=bad)
    assert AN.write_timing("vif", {"workers": 2}, timing_dir=tmp_path, wall_s=1.0)["workers"] == 2
