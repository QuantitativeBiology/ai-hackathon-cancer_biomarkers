"""Real dev loader of the M5 audits returns NaN columns for genes absent from the effect matrix (amplicon neighbours
that were never screened) instead of failing in the parquet reader. Reads two genes for five dev lines."""
import numpy as np

from biomark import audit as AU
from biomark import guard
from biomark import p1 as P1


class _D:
    def __init__(self, lines):
        self.lines = lines


def test_effects_of_returns_nan_for_unscreened_gene():
    cols = sorted(guard.effect_columns() - {"ModelID"})
    assert cols, "effect matrix has no gene columns"
    ids = [str(i) for i in P1.dev_ids()][:5]
    have = cols[0]
    missing = "999999999"
    assert missing not in guard.effect_columns()
    effects_of, _ = AU._real_loaders(_D(ids))
    E = effects_of([have, missing])
    assert list(E.columns) == [have, missing] and list(E.index) == ids
    assert E[missing].isna().all() and np.isfinite(E[have].to_numpy(dtype=float)).any()
