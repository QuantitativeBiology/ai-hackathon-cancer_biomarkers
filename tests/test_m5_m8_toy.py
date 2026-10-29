"""M5 / M6 / M8 scripts end to end on the toy (DESIGN v3.2 §8, §11, §12; amendment A7).

One module-scoped toy workspace: ``m5.prepare_toy`` runs the unchanged engine (four arms, no null, 1 worker) on
``toysim.make_toy()`` and writes the M4 outputs in the real layout; then every script's ``main(['--toy', DIR])``
runs in the real order. The planted labels must land where the answer key (configs/toy_answer_key.json:pass,
docs/toy_answer_key.md) states a class: co-location, cutting, culture-confounded, cross-lineage, lineage-only. The
bounds are the frozen answer-key bounds, applied to the catalogue (MD-global, A7) exactly as the toy suite applies
them; nothing is relaxed. Synthetic unit checks cover pieces the toy cannot reach (combinatorial labels and matched
background, shortlist ⊆ class 9 ∧ R replicated, chunked cross-lineage = one call, streamed EB = stats EB, rank lookup
= endpoints.RankIndex). Runs in a few minutes with 1 worker. Never touches the project's results/.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from biomark import audit as AU
from biomark import endpoints as EP
from biomark import m5 as M5
from biomark import paths as P
from biomark import stats as ST
from biomark import toysim as ts

ROOT = Path(__file__).resolve().parents[1]
BOUNDS = ts.load_toy_config()["pass"]
ORDER = ["run_subtypes", "run_cross_lineage", "run_m5_audits", "run_combinatorial_lookup", "run_limits",
         "build_candidate_cards"]
WORKER_ARGS = {"run_cross_lineage", "run_m5_audits"}


def _script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _synthetic_external(root: Path, toy) -> dict:
    """Esmaeili- and Flister-format files for three toy paralog pairs (planted PAR-like plants when available):
    one hit in two studies (replicated), one tested without hits (tested-not-replicated), the rest untested."""
    c = M5.cfg()["combinatorial"]
    pl = toy.plants("pan_lineage")
    sym = toy.genes["symbol"]
    pairs = []
    for _, r in pl.iterrows():
        f = r["feature"]
        fam = f.split(":")[0]
        if fam in c["loss_families"] and f.split(":")[1] != r["gene"]:
            pairs.append((sym[r["gene"]], sym[f.split(":")[1]]))
    assert len(pairs) >= 2, "toy has fewer than two paralog-loss plants"
    e = c["esmaeili"]
    base = root / "external" / e["dir"]
    base.mkdir(parents=True, exist_ok=True)
    sep = e["pair_sep"]
    lab = [f"{a}{sep}{b}" for a, b in pairs]
    s1, s2 = e["studies"][0], e["studies"][1]
    screened = pd.DataFrame({s1: pd.Series(lab[:2]), s2: pd.Series(lab[:2])})
    hits = pd.DataFrame({s1: pd.Series(lab[:1]), s2: pd.Series(lab[:1])})
    with pd.ExcelWriter(base / e["summary_file"]) as w:
        screened.to_excel(w, sheet_name=e["screened_sheet"], index=False)
        hits.to_excel(w, sheet_name=e["study_hits_sheet"], index=False)
    big, small = e["line_hit_size_min"] + 1, e["line_hit_effect_max"] - 1
    line = pd.DataFrame({"pair": lab[:2], e["line_size_col"]: [big, 0.0], e["line_effect_col"]: [small, 0.0]})
    with pd.ExcelWriter(base / e["per_line_file"]) as w:
        line.to_excel(w, sheet_name=f"{s1}_L1", index=False)
        line.to_excel(w, sheet_name=f"{s2}_L2", index=False)
    return {"replicated": pairs[0], "tested": pairs[1], "untested": pairs[2:]}


@pytest.fixture(scope="module")
def ws(tmp_path_factory):
    root = tmp_path_factory.mktemp("toy_m5_m8")
    w = M5.prepare_toy(root, workers=1)
    src = M5.source_for(w)
    ext = _synthetic_external(root, src.toy)
    for name in ORDER:
        args = ["--toy", str(root)] + (["--workers", "1"] if name in WORKER_ARGS else [])
        assert _script(name).main(args) == 0, name
    return {"ws": w, "src": src, "toy": src.toy, "root": root, "ext": ext}


def _read(ws, rel, **kw):
    return ws["ws"].read_csv(rel, **kw)


# ------------------------------------------------------------------------------------------------
# Plumbing
# ------------------------------------------------------------------------------------------------


def test_outputs_written_only_under_the_toy_root(ws):
    root = ws["root"]
    for rel in ["results/subtypes/subtypes.csv", "results/subtypes/modules.csv", "results/lineage/cross_lineage.csv",
                "results/locus/resolution.csv", "results/locus/amplicon.csv", "results/subtypes/culture_audit.csv",
                "results/subtypes/hypermut_tp53.csv", "results/combinatorial/lookup.csv",
                "results/combinatorial/validation.json", "results/limits/pair_classes.csv",
                "results/limits/gene_classes.csv", "results/candidates/known_recovery.csv",
                "results/candidates/shortlist.csv", "results/registry.csv"]:
        assert (root / rel).exists(), rel
    stems = {f.stem for f in (root / "results" / "timing").glob("*.json")}
    want = {"subtypes", "modules", "cross_lineage", "locus_resolution", "amplicon_test", "culture_audit",
            "hypermut_tp53", "combinatorial_lookup", "limits", "candidate_cards"}
    assert want <= stems
    mapped = {s for v in P.load_config("tiers")["stage_job_map"].values() for s in v}
    assert want <= mapped
    reg = pd.read_csv(root / "results" / "registry.csv")
    assert set(reg["status"]) == {"toy"}


@pytest.mark.parametrize("name", ORDER)
def test_scripts_refuse_before_freeze_bounds(name, monkeypatch, tmp_path):
    mod = _script(name)
    monkeypatch.setattr(EP, "descends_from", lambda tag, repo=None: False)
    with pytest.raises(EP.EndpointRefusal, match="freeze-bounds"):
        mod.main(["--toy", str(tmp_path)])
    assert not any(tmp_path.iterdir()), "nothing may be written before the tag check"


def test_catalogue_is_md_global_with_md_strict_alongside(ws):
    cat = M5.catalogue(ws["ws"])
    mg = ws["ws"].md("MD-global")
    assert set(zip(cat.loc[cat["in_md_global"], "gene"], cat.loc[cat["in_md_global"], "feature"])) == \
        set(zip(mg["gene"], mg["feature"]))
    sub = _read(ws, "results/subtypes/subtypes.csv")
    assert set(zip(sub["gene"], sub["feature"])) == set(zip(cat["gene"], cat["feature"]))


# ------------------------------------------------------------------------------------------------
# Planted labels (answer-key bounds, unchanged)
# ------------------------------------------------------------------------------------------------


def test_co_location(ws):
    lr = _read(ws, "results/locus/resolution.csv").astype({"feature": str})
    co = ws["toy"].answer_key["colocated"]
    ok = pass_gene = true_gene = tight_gene = 0
    for _, r in co.iterrows():
        sub = lr[(lr["gene"] == r["gene"]) & (lr["locus_group"] == r["locus_group"])]
        pas_gene = ((sub["feature"] == r["passenger"]) & (sub["level"] == "gene")).any()
        tr_gene = ((sub["feature"] == r["true_feature"]) & (sub["level"] == "gene")).any()
        if r["kind"] == "co_located":
            ok += len(sub) > 0 and not pas_gene
            pass_gene += pas_gene
            true_gene += tr_gene
        else:
            tight_gene += int((sub["level"] == "gene").sum())
    assert ok >= BOUNDS["co_located_locus_level_min"]
    assert pass_gene <= BOUNDS["co_located_passenger_gene_level_max"]
    assert true_gene >= BOUNDS["co_located_true_gene_level_min"]
    assert tight_gene <= BOUNDS["co_located_tight_gene_level_max"]


def test_cutting_and_addiction(ws):
    am = _read(ws, "results/locus/amplicon.csv").set_index(["gene", "feature"])["label"]
    toy = ws["toy"]
    cut = toy.plants("cutting_artefact").query("subkind == 'core'")
    assert sum(am.get((g, f)) == "cutting artefact" for g, f in zip(cut["gene"], cut["feature"])) \
        >= BOUNDS["cutting_core_label_min"]
    own = toy.plants("own_locus_cn")
    assert sum(am.get((g, f)) == "amplified-driver addiction" for g, f in zip(own["gene"], own["feature"])) \
        >= BOUNDS["own_locus_addiction_label_min"]


def test_culture_confounded(ws):
    traps = ws["toy"].answer_key["culture"]
    key = set(zip(traps["gene"], traps["feature"]))
    cat = M5.catalogue(ws["ws"])
    for col in ("in_md_global", "in_md_strict"):
        c = cat[cat[col]]
        assert sum(k in key for k in zip(c["gene"], c["feature"])) <= BOUNDS["culture_confounded_in_md_strict_max"]
    pc = _read(ws, "results/limits/pair_classes.csv")
    hit = pc[[k in key for k in zip(pc["gene"], pc["feature"])]]
    confounded = M5.PR.pp()["pair_classes"]["confounded"]
    assert (hit["pair_class"] <= confounded).all(), "a culture trap in the catalogue must not be a valid pair"
    de = _read(ws, "results/confounders/delta.csv")
    rep = de[np.array([k in key for k in zip(de["gene"], de["feature"])]) & (de["block"] == "culture").to_numpy()]
    if len(rep):
        assert rep["flagged"].astype(bool).mean() >= BOUNDS["culture_flagged_if_reported_frac_min"]


def test_cross_lineage_labels(ws):
    cl = _read(ws, "results/lineage/cross_lineage.csv").set_index(["gene", "feature"])["label"]
    toy = ws["toy"]
    pan = toy.plants("pan_lineage")
    lab = [cl.get((g, f)) for g, f in zip(pan["gene"], pan["feature"])]
    lab = [x for x in lab if x is not None]
    assert lab and np.mean([x == "cross-lineage" for x in lab]) >= BOUNDS["cross_lineage_pan_frac_min"]
    rs = toy.plants("lineage_restricted")
    assert sum(cl.get((g, f)) == "cross-lineage" for g, f in zip(rs["gene"], rs["feature"])) \
        <= BOUNDS["cross_lineage_restricted_max"]
    full = _read(ws, "results/lineage/cross_lineage.csv")
    assert "single_lineage_driven" in full and full["source"].notna().all()


def test_lineage_only(ws):
    lo = set(ws["toy"].plants("lineage_only")["gene"])
    cat = M5.catalogue(ws["ws"])
    for col in ("in_md_global", "in_md_strict"):
        assert cat.loc[cat[col], "gene"].isin(lo).sum() <= BOUNDS["lineage_only_in_md_strict_max"]
    gc = _read(ws, "results/limits/gene_classes.csv").set_index("gene")["failure_class"].astype(str)
    assert gc.reindex(sorted(lo)).eq("3").mean() >= BOUNDS["lineage_only_class3_frac_min"]
    ld = ws["toy"].plants("lineage_dummy_feature").iloc[0]
    assert (cat["feature"] == ld["feature"]).sum() <= BOUNDS["lineage_dummy_in_md_strict_max"]


# ------------------------------------------------------------------------------------------------
# Subtypes, modules, audits, limits, cards on the toy
# ------------------------------------------------------------------------------------------------


def test_subtype_members_are_accurate_dependent_carriers(ws):
    sub = _read(ws, "results/subtypes/subtypes.csv")
    lines = _read(ws, "results/subtypes/subtype_lines.csv")
    assert (sub["n_members"] > 0).any()
    data = ws["src"].data()
    dep, acc = M5.dependency_status(data)
    smin = M5.PR.dp()["subtypes_s8"]["line_share_min"]
    assert (lines["share"] >= smin).all()
    for (g, f), s in lines.groupby(["gene", "feature"]):
        assert (acc.loc[s["ModelID"], g] == 1).all()
        if data.fmeta.at[f, "d_type"] == "binary":
            x = pd.Series(data.X[:, data.fmeta.index.get_loc(f)], index=data.lines)
            assert (x.loc[s["ModelID"]] == 1).all()
    pan = ws["toy"].plants("pan_lineage")
    hit = sub[[k in set(zip(pan["gene"], pan["feature"])) for k in zip(sub["gene"], sub["feature"])]]
    assert len(hit) and (hit["n_members"] > 0).mean() >= BOUNDS["pan_lineage_recall_min"]


def test_modules_and_clusters(ws):
    m = _read(ws, "results/subtypes/modules.csv", genes=str)
    assert set(m["list"]) <= {"MD-global", "MD-strict"} and len(m)
    assert ((m["n_effective_clusters"] >= 1) & (m["n_effective_clusters"] <= m["n_genes"])).all()


def test_hypermut_tp53_and_culture_tables(ws):
    h = _read(ws, "results/subtypes/hypermut_tp53.csv")
    fams = set(M5.cfg()["hypermutation"]["persist_families"]) | set(M5.cfg()["hypermutation"]["tp53_families"])
    assert set(h["family"]) <= fams and len(h)
    msi = ws["toy"].plants("msi_like")
    sig = h[[k in set(zip(msi["gene"], msi["feature"])) for k in zip(h["gene"], h["feature"])]]
    if len(sig):
        assert (sig["test"] == "TP53-mutant").all()
    cu = _read(ws, "results/subtypes/culture_audit.csv")
    assert {"V_max", "culture_flag", "rescued", "culture_confounded"} <= set(cu.columns)
    assert (cu["culture_confounded"] <= cu["culture_flag"]).all()


def test_pair_classes_pending_R_and_gene_classes(ws):
    pc = _read(ws, "results/limits/pair_classes.csv")
    codes = M5.PR.pp()["pair_classes"]
    pend = pc[pc["pair_class"].isin([codes["valid"], M5.cfg()["limits"]["pair_class_codes"]["ko_drug_discordance"]])]
    assert len(pend) and (pend["class_status"] == "pending R").all()
    gc = _read(ws, "results/limits/gene_classes.csv")
    toy = ws["toy"]
    low = set(toy.plants("low_tier")["gene"])
    assert gc.set_index("gene")["failure_class"].astype(str).reindex(sorted(low)).eq("2").all()
    assert set(gc["gene"]) >= set(toy.genes.index[toy.genes["variable"] | toy.genes["selective"]])


def test_known_recovery_and_empty_shortlist_before_R(ws):
    kr = _read(ws, "results/candidates/known_recovery.csv")
    toy = ws["toy"]
    assert set(toy.controls["control_id"]) <= set(kr.loc[kr["kind"] == "control", "id"])
    assert (kr["R_status"] == "pending R").all()
    gate = kr[(kr["kind"] == "control") & (kr["scope"] == "gate")]
    assert gate["recovered_MD-global"].mean() >= BOUNDS["gate_overall_recall_min"]
    sl = _read(ws, "results/candidates/shortlist.csv")
    assert len(sl) == 0
    summ = ws["ws"].read_json("results/limits/limits_summary.json")
    assert summ["R_status"] == "pending R"
    assert summ["NC6_lineage_tf_class3"]["n"] == int((toy.controls["scope"] == "class_check").sum())


def test_combinatorial_labels_on_synthetic_tables(ws):
    lk = _read(ws, "results/combinatorial/lookup.csv", partner=str)
    ext = ws["ext"]
    lab = {(a, b): x for a, b, x in zip(lk["symbol"], lk["partner_symbol"], lk["label"])}
    lab.update({(b, a): x for (a, b), x in list(lab.items())})
    labels = M5.cfg()["combinatorial"]["labels"]
    assert lab.get(ext["replicated"]) == labels[0]          # hit in two studies
    assert lab.get(ext["tested"]) == labels[1]              # screened, never a hit
    for k in ext["untested"]:
        if k in lab:
            assert lab[k] == labels[2]
    val = ws["ws"].read_json("results/combinatorial/validation.json")
    assert val["_status"] == "ok" and set(M5.cfg()["combinatorial"]["lists"]) <= set(val)
    raw = (ws["root"] / "results" / "combinatorial" / "lookup.csv").read_text(encoding="utf-8")
    e = M5.cfg()["combinatorial"]["esmaeili"]
    assert e["line_size_col"] not in raw and e["line_effect_col"] not in raw      # no external values written


# ------------------------------------------------------------------------------------------------
# Synthetic unit checks
# ------------------------------------------------------------------------------------------------


def test_combinatorial_rules_and_matched_background():
    from biomark import combinatorial as CB
    c = M5.cfg()["combinatorial"]
    rep, nrep, unt = c["labels"]
    k1, k2, k3 = CB.pair_key("1", "2"), CB.pair_key("3", "4"), CB.pair_key("5", "6")
    idx = {"screens": {"A": {"tested": {k1, k2}, "hit": {k1}}, "B": {"tested": {k1}, "hit": {k1}}},
           "lines": [{"study": "A", "line": "x", "tested": {k2}, "hit": {k2}}]}
    lab = CB.label_of([k1, k2, k3], idx)
    assert list(lab["label"]) == [rep, nrep, unt]
    idx["lines"].append({"study": "A", "line": "y", "tested": {k2}, "hit": {k2}})
    assert CB.label_of([k2], idx)["label"].iloc[0] == rep              # >= 2 cell lines
    rng = np.random.default_rng(0)
    n = 400
    bg = pd.DataFrame({"key": [("b", str(i)) for i in range(n)], "identity": rng.uniform(20, 90, n),
                       "family_size": rng.integers(1, 12, n), "label": np.where(rng.random(n) < 0.2, rep, nrep)})
    cand = pd.DataFrame({"list": "MD-global", "relation": "paralog", "sl_direction": True,
                         "key": [("c", str(i)) for i in range(30)], "identity": rng.uniform(20, 90, 30),
                         "family_size": rng.integers(1, 12, 30), "label": [rep] * 15 + [nrep] * 15})
    v = CB.validation(cand, bg, ["MD-global"])["MD-global"]
    assert v["n_tested"] == 30 and v["n_replicated"] == 15 and v["rate"] == 0.5
    assert 0 < v["matched_expected_rate"] < v["rate"] and v["vs_background_odds_ratio"] > 1


def test_streamed_eb_equals_stats_eb():
    rng = np.random.default_rng(1)
    s = rng.uniform(0.1, 1.0, 3000)
    b = rng.normal(0, 1, 3000) * (rng.random(3000) < 0.2) + rng.normal(0, s)
    ref = ST.eb_scale_mixture(b, s)
    sig, pi = M5.eb_fit_streamed(b, s, chunk=257)
    assert np.allclose(sig, ref["sigma"]) and np.allclose(pi, ref["pi"], atol=1e-10)
    assert np.allclose(M5.eb_posterior(b, s, sig, pi), ref["posterior_mean"], atol=1e-10)


def test_chunked_cross_lineage_equals_single_call(ws):
    data = ws["src"].data()
    cat, _, pairs = M5.analysis_pairs(ws["ws"], ws["src"])
    pairs = pairs.head(24)
    bsh = M5.shrunk_beta_ws(ws["ws"], set(zip(pairs["gene"], pairs["feature"])))
    one = AU.cross_lineage(data, pairs, bsh, M5.sign_map(cat))
    R, P5 = [], []
    for i in range(0, len(pairs), 7):
        r, p = AU.cross_lineage_stats(data, pairs.iloc[i:i + 7], bsh, M5.sign_map(cat))
        R.append(r), P5.append(p)
    many = AU.cross_lineage_labels(pd.concat(R), pd.concat(P5))
    assert list(one["label"]) == list(many["label"])
    assert np.allclose(one["q_pooled"].to_numpy(float), many["q_pooled"].to_numpy(float), equal_nan=True)


def test_rank_lookup_equals_rank_index(ws):
    w = ws["ws"]
    a = w.arms()
    ridx = EP.RankIndex(a)
    genes = sorted(a["gene"].unique())[:15]
    rl = M5.RankLookup(w, ["S1-obs", "EN-pred"], genes)
    sub = a[a["gene"].isin(genes) & a["arm"].isin(["S1-obs", "EN-pred"])].sample(200, random_state=0)
    for arm, g, f in zip(sub["arm"], sub["gene"], sub["feature"]):
        x, y = ridx.rank_of(arm, g, f), rl.rank_of(arm, g, f)
        assert (np.isnan(x[0]) and np.isnan(y[0])) or x[0] == y[0]


def test_shortlist_is_class9_and_replicated(ws):
    """Pseudo r_tests (half of the class-9 pairs replicated): the shortlist is exactly those, ordered by label."""
    from biomark import limits as LM
    w = ws["ws"]
    pc = _read(ws, "results/limits/pair_classes.csv")
    for col in ("in_md_global", "in_md_strict"):
        pc[col] = pc[col].astype(str).str.lower().isin(["true", "1"])
    valid = M5.PR.pp()["pair_classes"]["valid"]
    v = pc[(pc["pair_class"] == valid) & pc["in_md_global"]]
    assert len(v) >= 2
    rt = pd.DataFrame({"gene": v["gene"], "feature": v["feature"], "power": 1.0,
                       "replicated": [i % 2 == 0 for i in range(len(v))], "p_one": np.linspace(1e-4, 1e-2, len(v)),
                       "n_carriers": 5})
    rt = pd.concat([rt, pd.DataFrame({"gene": ["x"], "feature": ["y"], "power": [1.0], "replicated": [True],
                                      "p_one": [1e-9], "n_carriers": [5]})])
    cl = _read(ws, "results/lineage/cross_lineage.csv")
    cat = M5.catalogue(w)
    sl = LM.shortlist(pc, rt, cl, None, cat, {}, ws["src"].symbols(), {}, {}, {}, None, {}, set())
    want = set(zip(v["gene"].iloc[::2], v["feature"].iloc[::2]))
    assert set(zip(sl["gene"], sl["feature"])) == want
    prec = M5.cfg()["cards"]["cross_lineage_precedence"]
    order = [prec.index(x) if x in prec else len(prec) for x in sl["cross_lineage_label"]]
    assert order == sorted(order)
    mod = _script("build_candidate_cards")
    mod.check_shortlist(sl, pc, rt)
    bad = pd.concat([sl, pd.DataFrame({"gene": ["x"], "feature": ["y"]})])
    with pytest.raises(AssertionError):
        mod.check_shortlist(bad, pc, rt)
