"""Toy answer key for the MD-global catalogue (DESIGN amendment A7), written by the verification agent.

MD-global = EN-pred selected ∩ S1-obs F1 BH hits (full-family BH), same sign. A7 applies the existing frozen toy error
bounds to it unchanged (configs/toy_answer_key.json:pass; no new numbers):
- FDP_true <= planted_fdp_true_max on the smooth truth toy and on the ridge truth toy;
- P(any non-neutral hit) <= pure_null_bh_lists_p_any_hit_max over the pure-null toys.

Runs only with BIOMARK_TOY_PIPELINE=1, like part B of tests/test_toy_answers.py. The suite's pipeline runs are reused:
this module wraps that module's ``PIPE`` in a memo keyed by (toy seed, scenario, predictor, call arguments), and
calls the pipeline with exactly the suite's arguments. In one pytest session each toy is therefore run once. Run alone,
the module computes the same runs itself.

Until the engine exposes an "MD-global" list, the list is computed here from the S1-obs and EN-pred arms. Once the
engine provides it (``out["md_global"]`` or arm "MD-global"), the two must match exactly.
"""
import os

import numpy as np
import pandas as pd
import pytest

import test_toy_answers as TA
from biomark import paths as P
from biomark import toysim as ts

DP = P.load_config("design_params")
BOUNDS = ts.load_toy_config()["pass"]
NAME = "MD-global"


def _memo(fn):
    if fn is None or getattr(fn, "_md_global_memo", False):
        return fn
    cache = {}

    def wrapped(toy, **kw):
        key = (int(toy.seed), toy.scenario, toy.predictor, tuple(sorted(kw.items())))
        if key not in cache:
            cache[key] = fn(toy, **kw)
        return cache[key]

    wrapped._md_global_memo = True
    wrapped.cache = cache
    return wrapped


TA.PIPE = _memo(TA.PIPE)          # the suite's fixtures look PIPE up at call time, so they share the memo
PIPE = TA.PIPE

needs_pipeline = pytest.mark.skipif(PIPE is None or not TA.RUN_PIPE, reason=TA.needs_pipeline.kwargs["reason"])


def md_global_from_arms(out) -> pd.DataFrame:
    """EN-pred selected ∩ S1-obs selected (F1 BH, not collinear), same sign of β."""
    a = out["arms"]
    en = a[(a["arm"] == "EN-pred") & a["selected"].astype(bool)][["gene", "feature", "beta"]]
    s1 = a[(a["arm"] == "S1-obs") & a["selected"].astype(bool)][["gene", "feature", "beta"]]
    m = en.astype({"gene": str}).merge(s1.astype({"gene": str}), on=["gene", "feature"], suffixes=("_en", ""))
    return m[np.sign(m["beta"]) == np.sign(m["beta_en"])].reset_index(drop=True)


def engine_md_global(out):
    """The engine's MD-global list if it exists, else None."""
    if "md_global" in out:
        return out["md_global"]
    if NAME in set(out["arms"]["arm"]):
        return ts.list_hits(out, NAME)
    return None


def md_global(out) -> pd.DataFrame:
    mine = md_global_from_arms(out)
    eng = engine_md_global(out)
    if eng is not None:
        a = set(zip(mine["gene"].astype(str), mine["feature"].astype(str)))
        b = set(zip(eng["gene"].astype(str), eng["feature"].astype(str)))
        assert a == b, f"engine MD-global differs from EN-pred ∩ S1-obs: only engine {sorted(b - a)[:5]}, " \
                       f"only verifier {sorted(a - b)[:5]}"
    return mine


@pytest.fixture(scope="module")
def toy():
    return ts.make_toy()


@pytest.fixture(scope="module")
def out_truth(toy):
    return PIPE(toy)                                   # same call as test_toy_answers.out_truth


@pytest.fixture(scope="module")
def ridge_toy():
    return ts.make_toy(predictor="ridge")


@pytest.fixture(scope="module")
def out_ridge(ridge_toy):
    return PIPE(ridge_toy)                             # same call as test_toy_answers.out_ridge


@pytest.fixture(scope="module")
def pure_null_outs():
    B = int(os.environ.get("BIOMARK_TOY_PURE_NULL_B", BOUNDS["pure_null_B"]))
    base = P.seeds()[DP["toy_s13"]["null_seed_key"]]
    outs = []
    for b in range(B):
        t = ts.make_toy(seed=base + b, scenario="null")
        outs.append((t, PIPE(t, null_B=0)))            # same call as test_toy_answers.pure_null_outs
    return outs


def test_suite_runs_are_shared():
    """The suite module pytest imported is the one whose PIPE carries the memo (no duplicate pipeline runs)."""
    import sys
    mod = sys.modules.get("test_toy_answers")
    assert mod is TA
    if PIPE is not None:
        assert mod.PIPE is PIPE and PIPE._md_global_memo


def test_md_global_helper_on_synthetic_arms():
    """Intersection, sign rule and collinear exclusion (via S1 'selected') on a hand-made arms table."""
    a = pd.DataFrame({
        "arm": ["EN-pred", "EN-pred", "EN-pred", "S1-obs", "S1-obs", "S1-obs", "EN-pred", "S1-obs"],
        "gene": ["1", "1", "2", "1", "1", "2", "3", "3"],
        "feature": ["A", "B", "C", "A", "B", "C", "D", "D"],
        "beta": [-1.0, 0.5, -0.2, -0.3, -0.4, -0.1, 1.0, 1.0],
        "selected": [True, True, True, True, True, False, False, True]})
    m = md_global_from_arms({"arms": a})
    assert set(zip(m["gene"], m["feature"])) == {("1", "A")}   # B: sign differs; C: not S1-selected; D: not EN-selected
    assert engine_md_global({"arms": a}) is None


@needs_pipeline
def test_pipeline_md_global_fdp_smooth_toy(out_truth, toy):
    f = ts.fdp_true(md_global(out_truth), toy)
    print(f"{NAME} smooth toy:", f)
    assert f["fdp"] <= BOUNDS["planted_fdp_true_max"], f


@needs_pipeline
def test_pipeline_md_global_fdp_ridge_toy(out_ridge, ridge_toy):
    f = ts.fdp_true(md_global(out_ridge), ridge_toy)
    print(f"{NAME} ridge toy:", f)
    assert f["fdp"] <= BOUNDS["planted_fdp_true_max"], f


@needs_pipeline
def test_pipeline_md_global_pure_null(pure_null_outs):
    hits = [ts.fdp_true(md_global(out), t)["n_hits"] for t, out in pure_null_outs]   # neutral pairs excluded
    p_any = float(np.mean([h > 0 for h in hits]))
    print(f"{NAME} pure null: P(any non-neutral hit) = {p_any}, hits per toy = {hits}")
    assert p_any <= BOUNDS["pure_null_bh_lists_p_any_hit_max"], (p_any, hits)
