"""Dictionary D (DESIGN v3.1 3): carrier counts vs raw caches, locus groups, self-null exclusions, cap,
continuous minima, data-status anchor, D-only fitting and no replication-effect reads."""
import ast
import hashlib
import json

import numpy as np
import pandas as pd
import pytest

from biomark import dictionary as Dct
from biomark import p1
from biomark import paths as P

pytestmark = pytest.mark.skipif(not (P.CONFIGS / "dictionary_v1.json").exists(), reason="dictionary not built")


@pytest.fixture(scope="module")
def meta():
    return Dct.load_meta()


@pytest.fixture(scope="module")
def feats(meta):
    return pd.DataFrame(meta["features"]).set_index("id", drop=False)


@pytest.fixture(scope="module")
def roles():
    return pd.read_parquet(P.CACHE / "covariates.parquet", columns=["role"])["role"]


@pytest.fixture(scope="module")
def gf():
    return pd.read_csv(P.CONFIGS / "dictionary_gene_features.csv", dtype={"entrez": str})


def _ez(sym):
    return p1.entrez_of([sym]).iloc[0]


def _raw_cn_count(gene, op, thr, ids):
    cn = p1.raw_matrix("cn", columns=[gene], ids=list(ids), strict=False)[gene]
    return int((op(cn, thr)).sum())


@pytest.mark.parametrize("sym,fam", [("CDKN2A", "DEL"), ("MTAP", "DEL"), ("PTEN", "DEL"), ("ERBB2", "AMP"),
                                     ("MYCN", "AMP"), ("CCNE1", "AMP")])
def test_cn_carrier_counts_match_raw(feats, roles, sym, fam):
    s3 = P.load_config("design_params")["dictionary_s3"]
    fid = f"{fam}:{_ez(sym)}"
    assert fid in feats.index
    op, thr = (np.less, s3["cn_del"]) if fam == "DEL" else (np.greater, s3["cn_amp"])
    for role, key in (("dev", "n_carriers_D"), ("test", "n_carriers_R")):
        ids = roles.index[roles == role]
        assert _raw_cn_count(_ez(sym), op, thr, ids) == feats.at[fid, key]


@pytest.mark.parametrize("sym", ["KRAS", "NRAS", "BRAF", "PIK3CA", "CTNNB1"])
def test_hot_carrier_counts_match_raw(feats, roles, sym):
    g = _ez(sym)
    h = p1.raw_matrix("hotspot", columns=[g])[g]
    for role, key in (("dev", "n_carriers_D"), ("test", "n_carriers_R")):
        ids = roles.index[roles == role]
        assert int((h.reindex(ids) > 0).sum()) == feats.at[f"HOT:{g}", key]


def test_loss_counts_match_raw(feats, roles):
    s3 = P.load_config("design_params")["dictionary_s3"]
    loss = feats[feats["family"] == "LOSS"].sample(20, random_state=0)
    ids = roles.index[roles == "dev"]
    genes = [i.split(":")[1] for i in loss.index]
    cn = p1.raw_matrix("cn", columns=genes, ids=list(ids), strict=False)
    for fid, g in zip(loss.index, genes):
        assert int((cn[g] < s3["cn_loss"]).sum()) == feats.at[fid, "n_carriers_D"]


def test_cached_matrices_match_json_counts(feats, roles):
    isD = (roles == "dev").to_numpy()
    for fam in ["HOT", "DEL", "AMP", "LOFBI", "SIL-METH", "SIL-other", "METH", "FUS", "YABS", "COMP"]:
        X = Dct.load_family(fam)
        got = np.nansum(X.to_numpy()[isD] == 1, axis=0)
        exp = feats.loc[X.columns, "n_carriers_D"].to_numpy()
        np.testing.assert_array_equal(got, exp)
        assert (exp >= P.load_config("design_params")["dictionary_s3"]["min_carriers"]).all()


def test_locus_cdkn2a_cdkn2b_mtap(feats):
    ids = {s: f"DEL:{_ez(s)}" for s in ("CDKN2A", "CDKN2B", "MTAP")}
    groups = {s: feats.at[i, "locus_group"] for s, i in ids.items()}
    assert groups["CDKN2A"] == groups["CDKN2B"] == groups["MTAP"] is not None
    assert len(set(ids.values())) == 3                      # separate features
    assert feats.at[ids["CDKN2A"], "hash"] != feats.at[ids["MTAP"], "hash"]
    assert feats.at[ids["CDKN2A"], "variant_group"] == feats.at[ids["MTAP"], "variant_group"]


def test_locus_rule_jaccard_window():
    lg = pd.read_csv(P.RESULTS / "dictionary" / "locus_groups.csv", dtype={"entrez": str})
    multi = lg[lg["group_size"] > 1]
    assert (multi.groupby("locus_group")["chrom"].nunique() == 1).all()


def test_self_null_exclusions(meta, feats, gf):
    lg = pd.read_csv(P.RESULTS / "dictionary" / "locus_groups.csv", dtype={"entrez": str})
    dl = lg[lg["family"] == "DEL"].set_index("entrez")["locus_group"]
    members = dl.groupby(dl).groups
    fm = feats.rename(columns={"d_type": "type"})
    targets = sorted(set(gf["entrez"]))
    rng = np.random.default_rng(0)
    pick = list(rng.choice(targets, 40, replace=False))
    # every target in the 9p21 DEL group, plus a random sample
    g9 = dl.get(_ez("CDKN2A"))
    pick += [t for t in targets if dl.get(t) == g9]
    for g in pick:
        d_g = set(Dct.features_for_gene(g, fm, gf))
        sn = {f"DEL:{g}", f"LOFBI:{g}", f"METH:{g}"} | {f"{f}:{g}" for f in Dct.SIL_FAMILIES}
        if g in dl.index:
            sn |= {f"DEL:{h}" for h in members[dl[g]]}
        assert not (d_g & sn), (g, d_g & sn)
        # own LOSS and AMP stay
        for f in (f"LOSS:{g}", f"AMP:{g}"):
            if f in feats.index:
                assert f in d_g
    # PRMT5 keeps DEL:MTAP (control) and its own LOSS
    pr = _ez("PRMT5")
    assert f"DEL:{_ez('MTAP')}" in set(Dct.features_for_gene(pr, fm, gf))


def test_feature_cap_respected(meta, feats):
    cap = P.load_config("design_params")["dictionary_s3"]["shared_feature_cap"]
    assert int((feats["scope"] == "shared").sum()) <= cap
    fc = pd.read_csv(P.RESULTS / "dictionary" / "family_counts.csv")
    assert fc["total"].iloc[-1] == int((feats["scope"] == "shared").sum())
    levels = P.load_config("design_params")["dictionary_s3"]["cap_escalation_min_carriers"]
    for fam, m in meta["cap"]["family_min_carriers"].items():
        sh = feats[(feats["family"] == fam) & (feats["scope"] == "shared") & (feats["d_type"] == "binary")]
        assert (sh["n_carriers_D"] >= m).all()
        assert m in [P.load_config("design_params")["dictionary_s3"]["min_carriers"]] + levels


def test_continuous_minima(feats):
    s3 = P.load_config("design_params")["dictionary_s3"]
    c = feats[feats["d_type"] == "continuous"]
    assert (c["n_nonnan_D"] >= s3["continuous_min_nonnan"]).all()
    assert (c["sd_D"] > 0).all()
    e = c[c["family"] == "EXPDOSE"]
    assert (e["n_expressed_D"] >= s3["expdose_min_expressed_lines"]).all()


def test_sil_mechanisms_exclusive(roles):
    X = [Dct.load_family(f) for f in Dct.SIL_FAMILIES if (Dct.CACHE_DIR / f"{f}.parquet").exists()]
    genes = {}
    for x in X:
        for c in x.columns:
            genes.setdefault(c.split(":")[1], []).append(x[c])
    for g, cols in genes.items():
        if len(cols) > 1:
            s = pd.concat(cols, axis=1).fillna(0).sum(axis=1)
            assert s.max() <= 1, g


def test_data_status_anchor(meta):
    sha = hashlib.sha256((P.RESULTS / "data_status.json").read_bytes()).hexdigest()
    assert meta["data_status_sha256"] == sha


def test_ids_and_families(feats):
    assert feats.index.is_unique
    assert feats["id"].str.match(r"^[A-Za-z-]+:[A-Za-z0-9_]+$").all()
    st = P.load_config("design_params")["surrogate_s4"]["stages"]
    allowed = set(st["G"]) | set(st["X"]) | {"METH", "SIL-METH"}
    assert set(feats["family"]) <= allowed


SRC = P.SRC / "biomark"
NO_EFFECTS = [SRC / "dictionary.py", SRC / "molecular.py", SRC / "genesets.py", SRC / "covariates.py",
              SRC / "datastatus.py", P.SCRIPTS / "build_dictionary.py", P.SCRIPTS / "build_gene_lists.py",
              P.SCRIPTS / "build_screen_covariates.py", P.SCRIPTS / "data_status.py"]


@pytest.mark.parametrize("path", NO_EFFECTS, ids=lambda p: p.name)
def test_dictionary_code_reads_no_effects(path):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    names = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name)} | \
            {n.attr for n in ast.walk(tree) if isinstance(n, ast.Attribute)}
    mods = {a.name for n in ast.walk(tree) if isinstance(n, (ast.Import, ast.ImportFrom)) for a in n.names}
    assert "load_effects" not in names and "guard" not in mods and "CRISPR_PARQUET" not in names
    consts = [n.value for n in ast.walk(tree) if isinstance(n, ast.Constant) and isinstance(n.value, str)]
    bad = ("crispr.parquet", "crisprgene", "screengeneeffect", "results/test", "cache/test")
    assert not any(b_ in c.lower() for c in consts for b_ in bad)


def test_effect_readers_use_dev_scope_only():
    for path in [SRC / "eligibility.py", P.SCRIPTS / "build_effect_nan_mask.py"]:
        tree = ast.parse(path.read_text(encoding="utf-8"))
        calls = [n for n in ast.walk(tree) if isinstance(n, ast.Call) and getattr(n.func, "attr", "") == "load_effects"]
        assert calls
        for c in calls:
            kw = {k.arg: k.value for k in c.keywords}
            assert isinstance(kw.get("scope"), ast.Constant) and kw["scope"].value == "dev"


def test_eligibility_runtime_reads_dev_only(monkeypatch):
    from biomark import eligibility, guard
    calls = []
    orig = guard.load_effects

    def rec(ids=None, scope="dev", purpose=None, genes=None):
        calls.append((scope, set(p1.roles_of(ids))))
        return orig(ids, scope=scope, purpose=purpose, genes=genes)

    monkeypatch.setattr(guard, "load_effects", rec)
    tiers = p1.gene_tiers("m5").iloc[:30]
    monkeypatch.setattr(p1, "gene_tiers", lambda model="m5": tiers)
    cov = pd.read_parquet(P.CACHE / "covariates.parquet")
    out = eligibility.compute(cov)
    assert len(out) == 30 and calls
    assert all(s == "dev" and r == {"dev"} for s, r in calls)


def test_d_only_fitting_invariance():
    """SIL, EXPDOSE and PSCORE fitted quantities do not change when R and set-A rows are deleted (DESIGN 1.7)."""
    genes = [_ez(s) for s in ("MLH1", "UBB", "MTAP", "VHL", "MGMT")]
    full = Dct.load_inputs(genes, genes)
    dro = Dct.load_inputs(genes, genes, drop_roles=("test", "noncancerous"))
    a, _ = Dct.sil_matrix(full, genes)
    b, _ = Dct.sil_matrix(dro, genes)
    np.testing.assert_array_equal(a[full.is_D], b[dro.is_D])
    xa = Dct.fam_expdose(full, genes)[0]
    xb = Dct.fam_expdose(dro, genes)[0]
    np.testing.assert_allclose(xa[full.is_D].to_numpy(), xb[dro.is_D].to_numpy(), atol=1e-6)
    pa = Dct.fam_pscore(full)[0]
    pb = Dct.fam_pscore(dro)[0]
    np.testing.assert_allclose(pa[full.is_D].to_numpy(), pb[dro.is_D].to_numpy(), atol=1e-5)


def test_gene_features_reference_existing(feats, gf):
    assert gf["feature_id"].isin(feats.index).all()
    assert set(gf["relation"]) >= {"own_EXPDOSE", "self_null"}
    assert gf["access_class"].isin(["direct", "potentially-direct", "expression-derived", "proxy"]).all()


def test_failed_sanity_spectra_excluded(meta, feats):
    """A spectrum failing its pre-declared sanity check is excluded from D and recorded (PI ruling 5)."""
    if "failed_sanity" not in meta:
        pytest.skip("dictionary built before the failed-sanity rule; rebuild pending")
    chk = json.loads((P.RESULTS / "dictionary" / "checks.json").read_text(encoding="utf-8"))["sig"]["spectra_checks"]
    m = P.load_config("build_params")["dictionary"]["spectra_sanity_feature"]
    failed = {fid for k, fid in m.items() if chk.get(k) is False}
    assert failed == {r["id"] for r in meta["failed_sanity"]}
    assert not (failed & set(feats.index))


def test_fusion_driver_orientation():
    """3' drivers except the 5' list (design FGFR2/3 + build_params extra, DEV-DICT-1 KMT2A)."""
    d = pd.read_csv(P.CONFIGS / "fusion_drivers.csv", dtype=str)
    five = set(P.load_config("design_params")["dictionary_s3"]["fusion"]["kinase_3prime_except_5prime"]) | \
        set(P.load_config("build_params")["dictionary"]["fusion_5prime_drivers_extra"])
    assert set(d.loc[d["orientation"] == "5prime", "symbol"]) == five
    assert set(d["orientation"]) <= {"3prime", "5prime"}


def test_cgt_mut_class():
    from biomark import cgt_features as CF
    m = pd.DataFrame({"ProteinChange": ["p.V600E", "p.E746_A750del", "p.L747FfsTer14", "", "p.D770_N771insG"],
                      "VariantInfo": ["missense_variant", "inframe_deletion", "frameshift_variant",
                                      "splice_acceptor_variant", "inframe_insertion"]})
    assert list(CF.mut_class(m)) == ["missense", "inframe", "frameshift", "splice", "inframe"]
