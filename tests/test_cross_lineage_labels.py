"""DESIGN §13.2 / §8.3 tests/test_cross_lineage_labels.py: each label of ``audit.cross_lineage`` is reachable
("a 2-lineage, 2-group effect β = −1 with 15 carriers each is cross-lineage") and every case gets exactly one label.
Synthetic data only (tests/synth_data.py; synthetic lineages are not in subtypes_s8.lineage_groups, so each is its
own organ-system group). The toy-pipeline version (gated by BIOMARK_TOY_PIPELINE) is
tests/test_toy_answers.py::test_pipeline_cross_lineage_labels.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import synth_data
from biomark import audit as AU

LABELS = {"cross-lineage", "lineage-restricted", "lineage-unresolvable",
          "lineage-unresolvable (single-lineage carriers)", "undetermined"}
HOT = "HOT:9200001"


def _data(monkeypatch, lineages, n, seed, plant):
    """make_data with the given lineages; ``plant(Y, hot, lin)`` edits the observed matrix."""
    monkeypatch.setattr(synth_data, "LIN", lineages)
    base = synth_data.make_data(n=n, seed=seed)
    Y = base.target("obs").copy()
    hot = base.X[:, base.fmeta.index.get_loc(HOT)]
    lin = base.cov_df["lineage"].to_numpy()
    plant(Y, hot, lin)
    Y[~base.rows_ok] = np.nan
    data = synth_data.make_data(n=n, seed=seed, obs_loader=lambda: Y)
    np.testing.assert_array_equal(data.X, base.X)
    return data, hot, lin


def _label(data, g, f, bsh=None):
    out = AU.cross_lineage(data, pd.DataFrame({"gene": [g], "feature": [f]}), bsh or {}, {})
    assert len(out) == 1 and out["label"].iloc[0] in LABELS            # exactly one label per case
    return out.iloc[0]


def test_two_lineage_two_group_effect_is_cross_lineage(monkeypatch):
    def plant(Y, hot, lin):
        Y[:, 2] += -1.0 * hot
    data, hot, lin = _data(monkeypatch, ("LinA", "LinB"), 300, 31, plant)
    assert min(int(hot[lin == L].sum()) for L in ("LinA", "LinB")) >= 15
    r = _label(data, data.genes[2], HOT)
    assert r["n_powered"] == 2 and r["n_organ_groups"] == 2
    assert r["label"] == "cross-lineage", r.to_dict()


def test_opposite_sign_lineage_effects_are_lineage_restricted(monkeypatch):
    def plant(Y, hot, lin):
        Y[:, 2] += np.where(lin == "LinA", -1.5, np.where(lin == "LinB", 1.0, 0.0)) * hot
    data, _, _ = _data(monkeypatch, ("LinA", "LinB", "LinC", "LinD"), 480, 22, plant)
    r = _label(data, data.genes[2], HOT)
    assert not (r["same_sign_frac"] >= 0.8)                              # NaN when no lineage is powered
    assert r["label"] == "lineage-restricted", r.to_dict()


def test_carriers_in_one_lineage_are_unresolvable_single_lineage(monkeypatch):
    data, hot, lin = _data(monkeypatch, ("LinA", "LinB", "LinC", "LinD"), 480, 23, lambda Y, h, l: None)
    j = data.fmeta.index.get_loc("LOFBI:9300002")
    data.X[:, j] = hot * (lin == "LinA")                                 # carriers only in LinA
    r = _label(data, data.genes[2], "LOFBI:9300002")
    assert r["n_testable"] <= 1 and bool(r["single_lineage_carriers"])
    assert r["label"] == "lineage-unresolvable (single-lineage carriers)", r.to_dict()


def test_powered_null_pair_is_undetermined(monkeypatch):
    data, _, _ = _data(monkeypatch, ("LinA", "LinB", "LinC", "LinD"), 480, 24, lambda Y, h, l: None)
    g = data.genes[2]
    r = _label(data, g, HOT, bsh={(g, HOT): -1.0})                       # prior effect large: lineages powered
    assert r["n_powered"] >= 2 and not (r["q_pooled"] < 0.05)
    assert r["label"] == "undetermined", r.to_dict()
