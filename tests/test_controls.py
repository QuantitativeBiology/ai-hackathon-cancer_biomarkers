"""Positive controls (DESIGN v3 9, review R3-2): schema, spec encodings, id resolution, carrier counts."""
import numpy as np
import pandas as pd
import pytest

from biomark import controls, genesets, p1
from biomark import dictionary as Dct
from biomark import paths as P

SPEC = P.CONFIGS / "positive_controls_spec.csv"
OUT = P.CONFIGS / "positive_controls.csv"
SCHEMA = ["control_id", "family", "dependency_symbol", "dependency_entrez", "feature_primary", "feature_secondary",
          "expected_sign", "scope", "citation", "in_cgt", "n_carriers_D", "n_carriers_R", "access_class",
          "en_max_attainable_freq", "gate_eligible", "note"]
FAMILIES = {"ADD-HOT", "ADD-CN", "ADD-EXP", "PAR", "LOF-SL", "SIL", "SIG", "PROG", "P53", "W", "LINEAGE-CHECK"}


@pytest.fixture(scope="module")
def spec():
    return pd.read_csv(SPEC, dtype=str, keep_default_na=False)


@pytest.fixture(scope="module")
def pc():
    if not OUT.exists():
        pytest.skip("positive_controls.csv not built")
    return pd.read_csv(OUT, dtype={"dependency_entrez": str}, keep_default_na=False)


def test_spec_allele_primaries(spec):
    """R3-2: BRAF and EGFR primaries are ALLELE features; HOT entries of allele genes are secondaries only."""
    # Only BRAF and EGFR: DESIGN 9.2 keeps KRAS/NRAS/PIK3CA on HOT primaries although allele classes exist for
    # them (R3-2's "every HOT entry of an allele gene is secondary" contradicts 9.2; logged in impl notes).
    for dep in ("BRAF", "EGFR"):
        r = spec[(spec["dependency_symbol"] == dep) & (spec["family"] == "ADD-HOT")].iloc[0]
        assert r["feature_primary"].startswith("ALLELE:")
        assert r["feature_secondary"] == f"HOT:{dep}"
    for dep in ("KRAS", "NRAS", "PIK3CA"):          # DESIGN 9.2: HOT primaries kept (PI ruling 3)
        r = spec[(spec["dependency_symbol"] == dep) & (spec["family"] == "ADD-HOT")].iloc[0]
        assert r["feature_primary"] == f"HOT:{dep}"
    s3 = P.load_config("design_params")["dictionary_s3"]
    for _, r in spec.iterrows():
        if r["feature_primary"].startswith("ALLELE:"):
            assert r["feature_primary"].split(":", 1)[1] in s3["allele_classes"]
    assert set(spec["family"]) <= FAMILIES
    assert set(spec["expected_sign"]) <= {"S", "W"}


def test_spec_resolves_or_fails():
    meta = Dct.load_meta() if (P.CONFIGS / "dictionary_v1.json").exists() else {"features": []}
    ids = {f["id"] for f in meta["features"]}
    spec = genesets.control_spec()
    for _, r in spec.iterrows():
        for c in ("feature_primary", "feature_secondary", "fallback_feature"):
            controls.resolve(r[c], ids)            # raises if unresolvable
    with pytest.raises(controls.UnresolvableFeature):
        controls.resolve("SIG:CIN_PROXY", ids)
    with pytest.raises(controls.UnresolvableFeature):
        controls.resolve("FOO:BAR", ids)
    with pytest.raises(controls.UnresolvableFeature):
        controls.resolve("LOSS:NOT_A_GENE_SYMBOL", ids)
    sil = controls.resolve("SIL:MLH1", ids)
    assert sil == [f"{f}:{p1.entrez_of(['MLH1']).iloc[0]}" for f in Dct.SIL_FAMILIES]


def test_schema(pc, spec):
    assert set(SCHEMA) <= set(pc.columns)
    assert pc["control_id"].is_unique
    assert set(pc["scope"]) <= {"gate", "class_check", "untestable"}
    assert set(pc["family"]) <= FAMILIES
    lc = pc[pc["family"] == "LINEAGE-CHECK"]
    assert len(lc) == 18 and (lc["scope"] == "class_check").all()
    assert len(pc) == len(spec) + 18
    ge = pc["gate_eligible"].astype(str) == "True"
    assert (pd.to_numeric(pc.loc[ge, "n_obs_D"]) >= P.load_config("build_params")["controls"]["gate_min_nonnan_D_lines"]).all()
    assert (pc.loc[ge, "scope"] == "gate").all()
    assert set(pc.loc[pc["family"] != "LINEAGE-CHECK", "access_class"]) <= \
        {"direct", "potentially-direct", "expression-derived", "proxy"}


def test_carrier_counts_from_caches(pc):
    roles = pd.read_parquet(P.CACHE / "covariates.parquet", columns=["role"])["role"]
    meta = Dct.load_meta()
    fam_of = {f["id"]: f["family"] for f in meta["features"]}
    for _, r in pc.iterrows():
        ids = [i for i in str(r.get("feature_primary_ids", "")).split(";") if i in fam_of]
        if not ids or r["n_carriers_D"] == "":
            continue
        X = pd.concat([Dct.load_family(fam_of[i])[i] for i in ids], axis=1)
        u = X.max(axis=1, skipna=True).where(X.notna().any(axis=1))
        assert int((u[roles == "dev"] == 1).sum()) == int(float(r["n_carriers_D"])), r["control_id"]
        assert int((u[roles == "test"] == 1).sum()) == int(float(r["n_carriers_R"])), r["control_id"]


def test_mtap_counts_vs_raw_cn(pc):
    s3 = P.load_config("design_params")["dictionary_s3"]
    roles = pd.read_parquet(P.CACHE / "covariates.parquet", columns=["role"])["role"]
    mtap = p1.entrez_of(["MTAP"]).iloc[0]
    cn = p1.raw_matrix("cn", columns=[mtap])[mtap]
    r = pc[pc["control_id"] == "PAR_PRMT5_MTAP"].iloc[0]
    for role, col in (("dev", "n_carriers_D"), ("test", "n_carriers_R")):
        ids = roles.index[roles == role]
        assert int((cn.reindex(ids) < s3["cn_del"]).sum()) == int(float(r[col]))


def test_untestable_means_below_floor(pc):
    g = pc[pc["scope"] == "untestable"]
    for _, r in g.iterrows():
        vals = [r.get(f"n_carriers_D_nonnan_{k}", "") for k in ("primary", "secondary", "fallback")]
        nums = [float(v) for v in vals if v not in ("", None) and not (isinstance(v, float) and np.isnan(v))]
        floor = P.load_config("design_params")["dictionary_s3"]["min_carriers"]
        assert all(n < floor for n in nums) or "pending" in r["note"] or float(r["n_obs_D"]) < 100 or \
            str(r.get("in_D_g_primary")) == "False"
