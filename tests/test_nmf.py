"""DESIGN §13.2 "NMF" with §1.7 (D-only fitting: NMF is fit on D and frozen, other lines transformed; "a test rebuilds
them with R rows deleted and asserts identical hashes") and §3 (pinned settings, nmf_seed): the PROG family is
deterministic, reproduces the cached dictionary family, and its D rows do not change when R and set-A rows are
deleted. Expression only (no effects). ~2 NMF fits of 1,022 x 5,000; a few minutes on 2 BLAS threads.
"""
from __future__ import annotations

import numpy as np
import pytest

from biomark import dictionary as Dct
from biomark import p1
from biomark import paths as P

pytestmark = pytest.mark.skipif(not (Dct.CACHE_DIR / "PROG.parquet").exists(), reason="dictionary cache not built")


@pytest.fixture(scope="module")
def fits():
    g = [str(x) for x in p1.entrez_of(["MLH1", "UBB"]).dropna()]
    seed = int(P.seeds()[P.load_config("design_params")["dictionary_s3"]["nmf"]["seed_key"]])
    full = Dct.load_inputs(g, g)
    dro = Dct.load_inputs(g, g, drop_roles=("test", "noncancerous"))
    a = Dct.fam_prog(full, seed)
    b = Dct.fam_prog(dro, seed)
    return full, dro, a, b


def test_d_only_fit_components_identical(fits):
    _, _, a, b = fits
    np.testing.assert_array_equal(a[2].to_numpy(), b[2].to_numpy())          # H (components) fit on D only
    assert a[3]["n_iter"] == b[3]["n_iter"]


def test_d_rows_scores_invariant_to_deleting_r_and_set_a(fits):
    full, dro, a, b = fits
    xa = a[0][full.is_D].to_numpy(dtype=np.float64)
    xb = b[0][dro.is_D].to_numpy(dtype=np.float64)
    assert a[0][full.is_D].index.equals(b[0][dro.is_D].index)
    np.testing.assert_allclose(xa, xb, rtol=0, atol=1e-6 * max(1.0, np.abs(xa).max()))


def test_reproduces_cached_prog_family(fits):
    full, _, a, _ = fits
    cached = Dct.load_family("PROG")
    cols = [c for c in a[0].columns if c in cached.columns]
    assert cols == list(a[0].columns)
    ids = a[0].index[full.is_D]
    np.testing.assert_allclose(a[0].loc[ids, cols].to_numpy(dtype=np.float64),
                               cached.reindex(ids)[cols].to_numpy(dtype=np.float64), rtol=0, atol=1e-5)
