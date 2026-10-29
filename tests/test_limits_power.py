"""DESIGN v3.2 §12.2 after R: a catalogue pair tested in R and not replicated must leave class 9 (class 5 when R power
< power_min, class 7 when powered), also when the frozen pool held no power for it (stage-X pairs: r_tests.csv power
NaN). Synthetic one-family tables; no data is read."""
import numpy as np
import pandas as pd
import pytest

from biomark import limits as LM
from biomark import params as PR

PC = PR.pp()["pair_classes"]
PMIN = float(PR.dp()["replication_s6"]["power_min"])


def _inputs(p_one=0.4, replicated=False, power=np.nan):
    fm = pd.DataFrame({"family": ["PROG"], "d_type": ["continuous"], "stage": ["X"], "genes": [()],
                       "locus_group": [None], "variant_group": [None]}, index=["PROG:1"])
    fd = LM.FeatureData(fm, {"g1": "G1"})
    cat = pd.DataFrame({"gene": ["g1"], "feature": ["PROG:1"], "beta": [-0.3], "lists": ["MD-global|MD-strict"],
                        "in_md_global": [True], "in_md_strict": [True], "n_carriers": [np.nan],
                        "collinear_with_design": [False]})
    empty = pd.DataFrame(columns=["gene", "feature", "block", "flagged", "level", "label"])
    rt = pd.DataFrame({"gene": ["g1"], "feature": ["PROG:1"], "replicated": [replicated], "p_one": [p_one],
                       "power": [power]})
    return fd, cat, empty, rt


def _cls(r_power=None, **kw):
    fd, cat, empty, rt = _inputs(**kw)
    out = LM.pair_classes(fd, cat, empty, empty, empty, None, None, {}, rt, None, None, r_power=r_power) \
        if r_power is not None else LM.pair_classes(fd, cat, empty, empty, empty, None, None, {}, rt, None, None)
    return out.set_index(["gene", "feature"]).loc[("g1", "PROG:1")]


def test_not_replicated_without_power_leaves_class_9():
    r = _cls()
    assert int(r["pair_class"]) != PC["valid"]
    assert r["class_status"] == "not replicated, power not computed"


@pytest.mark.parametrize("power, want", [(0.9, "powered_not_replicated"), (0.2, "too_rare")])
def test_r_power_gives_class_7_or_5(power, want):
    rp = pd.DataFrame({"gene": ["g1"], "feature": ["PROG:1"], "power": [power], "power_effective": [power]})
    r = _cls(r_power=rp)
    code = {"powered_not_replicated": 7, "too_rare": PC["too_rare"]}[want]
    assert int(r["pair_class"]) == code and float(r["power"]) == power


def test_replicated_powered_pair_stays_class_9_and_untested_rule():
    rp = pd.DataFrame({"gene": ["g1"], "feature": ["PROG:1"], "power": [0.9], "power_effective": [0.9]})
    assert int(_cls(r_power=rp, replicated=True, p_one=1e-4)["pair_class"]) == PC["valid"]
    r = _cls(r_power=rp, replicated=False, p_one=np.nan)              # R fit without p: rule for untested pairs
    assert int(r["pair_class"]) == PC["valid"] and r["class_status"] == "no R test"
    assert PMIN == 0.5
