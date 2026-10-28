"""DESIGN §13.2 "covariate drops (ERBB2 × AMP:ERBB2 drops own CN; WRN × SIG:MSI drops MSI; a lineage-dummy d is
flagged)" on the REAL D covariate design, dictionary features and own-CN / EXPGATE caches (no effects or predictions:
the universe pass of ``stats.s1_scan(compute=False)`` needs only the row mask, set to all D lines here).
The synthetic version of the same rules is tests/test_stats.py::test_covariate_drops_rules.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from biomark import dictionary as DCT
from biomark import p1 as P1
from biomark import paths as P
from biomark import stats as ST

pytestmark = pytest.mark.skipif(not (P.CACHE / "covariates.parquet").exists() or not DCT.CACHE_DIR.exists(),
                                reason="real covariate / dictionary caches not built")


def _cov(dev) -> pd.DataFrame:
    c = pd.read_parquet(P.CACHE / "covariates.parquet").reindex(dev)
    out = pd.DataFrame(index=dev)
    for k in ("lineage", "library_group", "cellformat", "serumfree", "medium_family", "doubling_time",
              "DT_unknown", "MSI", "DMG_BURDEN", "PatientID"):
        out[k] = c[k]
    out["tp53_status"] = c["TP53"]
    dropped = set(c.attrs.get("dropped", []))
    out["screen_days"] = 0.0 if "screen_days_dev_humagne" in dropped else c["screen_days_dev_humagne"]
    return out


@pytest.fixture(scope="module")
def universe():
    dev = P1.dev_ids()
    ez = P1.entrez_of(["ERBB2", "WRN", "TP53", "MDM2"]).astype(str).to_dict()
    erbb2, wrn = ez["ERBB2"], ez["WRN"]
    meta = pd.DataFrame(DCT.load_meta()["features"]).set_index("id")
    want = [f"AMP:{erbb2}", "SIG:MSI", f"HOT:{ez['TP53']}"]
    want = [f for f in want if f in meta.index]
    cols = {}
    for f in want:
        cols[f] = DCT.load_family(meta.at[f, "family"], columns=[f]).reindex(dev)[f].to_numpy(dtype=float)
    cov = _cov(dev)
    big = cov["lineage"].value_counts().index[0]
    cols["SIL:LINEAGE_DUMMY_TEST"] = (cov["lineage"] == big).astype(float).to_numpy()   # a lineage-dummy d
    fids = list(cols)
    fm = pd.DataFrame(index=fids)
    fm["family"] = [meta.at[f, "family"] if f in meta.index else "SIL-other" for f in fids]
    fm["stage_family"] = fm["family"]
    fm["d_type"] = [meta.at[f, "d_type"] if f in meta.index else "binary" for f in fids]
    fm["scope"] = "shared"
    fm["genes"] = [tuple(map(str, meta.at[f, "genes"])) if f in meta.index and isinstance(meta.at[f, "genes"], list)
                   else () for f in fids]
    fm["variant_group"] = [meta.at[f, "variant_group"] if f in meta.index else None for f in fids]
    genes = [erbb2, wrn]
    X = np.column_stack([cols[f] for f in fids])
    own = pd.read_parquet(DCT.CACHE_DIR / "own_cn_log2.parquet", columns=genes).reindex(dev).to_numpy(dtype=float)
    eg = pd.read_parquet(DCT.CACHE_DIR / "expgate.parquet", columns=genes).reindex(dev).to_numpy(dtype=float)
    drop = ST.DropInfo.build(fm, tp53_genes=(ez["TP53"],), mdm2_genes=(ez["MDM2"],))
    inp = ST.ScanInputs(X=X, fmeta=fm, Dg={g: np.arange(len(fids)) for g in genes},
                        rows_ok=np.ones((len(dev), len(genes)), dtype=bool), genes=genes, cov=ST.CovDesign(cov),
                        own_cn=own, expgate=eg, patients=cov["PatientID"].astype(str).to_numpy(), drop=drop)
    U = ST.s1_scan(inp, {}, compute=False)
    return U.set_index(["gene", "feature"]), ez


def _drops(U, g, f) -> set:
    return set(filter(None, str(U.loc[(g, f), "dropped_cov"]).split("|")))


def test_erbb2_amp_erbb2_drops_own_cn(universe):
    U, ez = universe
    f = f"AMP:{ez['ERBB2']}"
    assert (ez["ERBB2"], f) in U.index, "AMP:ERBB2 not in U for ERBB2"
    assert "own_CN" in _drops(U, ez["ERBB2"], f)
    assert "own_CN" not in _drops(U, ez["WRN"], f)                 # not WRN's own locus


def test_wrn_sig_msi_drops_msi(universe):
    U, ez = universe
    assert "MSI" in _drops(U, ez["WRN"], "SIG:MSI")
    assert "MSI" in _drops(U, ez["ERBB2"], "SIG:MSI")              # identity rule, gene independent
    assert "own_CN" not in _drops(U, ez["WRN"], "SIG:MSI")


def test_lineage_dummy_flagged_never_dropped(universe):
    U, ez = universe
    for g in (ez["ERBB2"], ez["WRN"]):
        r = U.loc[(g, "SIL:LINEAGE_DUMMY_TEST")]
        assert bool(r["collinear_with_design"]) and r["r2_design"] > 0.49
        assert "lineage" not in _drops(U, g, "SIL:LINEAGE_DUMMY_TEST")
    assert not U.loc[(ez["WRN"], "SIG:MSI"), "collinear_with_design"]


def test_tp53_feature_drops_tp53_and_dropped_cov_within_droppable_set(universe):
    U, ez = universe
    f = f"HOT:{ez['TP53']}"
    if (ez["WRN"], f) not in U.index:
        pytest.skip("HOT:TP53 not in U")
    assert "TP53" in _drops(U, ez["WRN"], f)
    allowed = {"MSI", "DMG_BURDEN", "TP53", "screen_days", "own_CN", "EXPGATE"}
    for s in U["dropped_cov"]:
        assert set(filter(None, str(s).split("|"))) <= allowed
