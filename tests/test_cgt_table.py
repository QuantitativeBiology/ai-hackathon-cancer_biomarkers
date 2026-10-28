"""Checks on the clinical ground truth (DESIGN v3 §10; R1-CL1, CL3, CL6, CL7, CL8; PI D5).

Run after `python -X utf8 scripts/build_cgt.py`. Reads only the CGT outputs, configs and the frozen clinical
sources (no molecular or dependency data).
"""
from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
CGT = ROOT / "configs" / "clinical_ground_truth.csv"
MAPPING = ROOT / "configs" / "cgt_curated_mapping.csv"
RULES = ROOT / "configs" / "cgt_allele_rules.json"
SUMMARY = ROOT / "results" / "cgt" / "summary.json"
DROPPED = ROOT / "results" / "clinical" / "cgt_dropped.csv"
BUILD = ROOT / "scripts" / "build_cgt.py"
S10 = json.loads((ROOT / "configs" / "design_params.json").read_text(encoding="utf-8"))["clinical_s10"]

REQUIRED = ["test_id", "feature", "feature_D", "target_set", "comparator", "indication_codes", "sign", "partition",
            "ko_class", "tier", "drugs", "source_ids", "source_type", "in_gate", "tumour_agnostic",
            "excluded_indications", "allele_level_required", "proxy", "citation", "n_carriers_D", "n_carriers_R",
            "testable", "untestable_reason"]
NEVER_PROXIED = ("ALLELE:KRAS_G12C", "ALLELE:EGFR_T790M", "ALLELE:EGFR_C797S", "ALLELE:KIT_D816V", "ALLELE:KIT_EX11",
                 "ALLELE:ERBB2_TKD")

pytestmark = pytest.mark.skipif(not CGT.exists(), reason="run scripts/build_cgt.py first")


@pytest.fixture(scope="module")
def cgt() -> pd.DataFrame:
    return pd.read_csv(CGT, dtype=str, keep_default_na=False)


@pytest.fixture(scope="module")
def build():
    spec = importlib.util.spec_from_file_location("build_cgt", BUILD)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_schema(cgt):
    miss = [c for c in REQUIRED if c not in cgt.columns]
    assert not miss, miss
    assert len(cgt) > 0


def test_vocabularies(cgt):
    assert set(cgt["sign"]) <= {"S", "W", "N"}
    assert set(cgt["partition"]) <= set(S10["partitions"])
    assert set(cgt["ko_class"]) <= set(S10["ko_classes"])
    assert set(cgt["tier"]) <= {"T1", "T2", "T3", "refuted"}
    assert set(cgt["proxy"]) <= {"0", "1"}
    for c in ("in_gate", "tumour_agnostic", "allele_level_required", "testable"):
        assert set(cgt[c]) <= {"0", "1"}, c
    comp = cgt["comparator"]
    assert comp.str.match(r"^(within-indication(:.+)?|pan-cancer\+lineage|listed-only)$").all()


def test_sign_partition_consistency(cgt):
    exp = {"sensitivity": "S", "bypass-resistance": "W", "clinically-refuted": "S", "expected-null": "N",
           "acquired-resistance": "N"}
    bad = cgt[cgt["partition"].map(exp) != cgt["sign"]]
    assert bad.empty, bad[["test_id", "partition", "sign"]].to_dict("records")
    acq = cgt[cgt["partition"] == "acquired-resistance"]
    assert (acq["comparator"] == "listed-only").all() and (acq["testable"] == "0").all()
    assert ((cgt["tier"] == "refuted") == (cgt["partition"] == "clinically-refuted")).all()


def test_unique_tests(cgt):
    assert not cgt["test_id"].duplicated().any()
    assert not cgt.duplicated(["feature", "target_set", "comparator"]).any()


def test_every_record_id_resolvable(cgt, build):
    src = build.load_sources()
    regex = pd.read_csv(MAPPING, dtype=str, keep_default_na=False).set_index("test_id")["ctg_biomarker_regex"]
    for _, r in cgt.iterrows():
        ids = [t for t in r["source_ids"].split(";") if t]
        assert ids, f"{r['test_id']} has no source id"
        primary = [t for t in ids if not (t.startswith("PMID:") and r["source_type"] != "pmid_only")]
        assert primary, f"{r['test_id']} has only supplementary PMIDs"
        for t in ids:
            ok, why, _ = build.resolve_id(t, src, regex.get(r["test_id"], ""))
            assert ok, f"{r['test_id']}: {t} unresolved ({why})"


def test_pmid_only_only_in_refuted(cgt):
    po = cgt[cgt["source_type"].str.contains("pmid_only")]
    assert (po["partition"] == "clinically-refuted").all()
    for _, r in po.iterrows():
        assert all(t.startswith("PMID:") for t in r["source_ids"].split(";"))
        for t in r["source_ids"].split(";"):
            assert (ROOT / "external" / "clinical" / "pubmed" / f"{t[5:]}.json").exists()


def test_t3_rows_have_matching_ctg(cgt, build):
    regex = pd.read_csv(MAPPING, dtype=str, keep_default_na=False).set_index("test_id")["ctg_biomarker_regex"]
    for _, r in cgt[cgt["tier"] == "T3"].iterrows():
        ncts = [t[4:] for t in r["source_ids"].split(";") if t.startswith("CTG:")]
        assert ncts, r["test_id"]
        assert any(re.search(regex[r["test_id"]], build.ctg_text(n), flags=re.I) for n in ncts), r["test_id"]


def test_cl1_corrections(cgt):
    imat = cgt[cgt["feature"].eq("ALLELE:KIT_D816V") & cgt["drugs"].str.contains("imatinib")]
    assert (imat["sign"] != "S").all() and (imat["partition"] == "acquired-resistance").all()          # CL1a
    assert cgt[cgt["feature"].eq("ALLELE:KIT_D816V") & (cgt["sign"] == "S")]["drugs"].str.contains("midostaurin").all()
    d17 = cgt[cgt["feature"].str.contains("DEL17P|TP53") & cgt["target_set"].str.contains("BCL2")]
    assert len(d17) == 1 and (d17["sign"] == "N").all() and (d17["ko_class"] == "population-only").all()  # CL1b
    braf_s = cgt[cgt["target_set"].eq("BRAF") & (cgt["sign"] == "S") & (cgt["proxy"] == "0") &
                 ~cgt["feature"].str.contains("FUS")]
    assert len(braf_s) >= 1 and braf_s["feature"].str.contains("V600").all()                               # CL1c
    egfr_amp = cgt[cgt["feature"].isin(["AMP:EGFR"]) | cgt["feature_D"].eq("AMP:1956")]
    assert (egfr_amp["partition"] == "clinically-refuted").all()                                           # CL1d
    assert cgt[cgt["test_id"].str.startswith("EGFR_ACT") & (cgt["proxy"] == "0")]["feature"].eq("ALLELE:EGFR_ACT").all()
    akt = cgt[cgt["target_set"].str.contains("AKT")]
    assert len(akt) >= 1 and akt["target_set"].map(lambda s: sorted(s.split(";")) == ["AKT1", "AKT2", "AKT3"]).all()  # CL1e
    ras = cgt[cgt["target_set"].eq("EGFR") & cgt["feature"].str.contains("RAS")]
    assert len(ras) == 1 and ras["sign"].iloc[0] == "W" and ras["partition"].iloc[0] == "bypass-resistance"   # CL1f
    esr = cgt[cgt["feature"].eq("ALLELE:ESR1_LBD")]
    assert len(esr) == 1 and esr["sign"].iloc[0] == "N"                                                     # CL1g
    assert not (cgt["feature"].str.contains("CDKN2A") & cgt["target_set"].str.contains("PRMT5|MAT2A")).any()  # CL1h


def test_ko_classes(cgt):
    fam = cgt[cgt["ko_class"] == "concordant-family"]
    assert (fam["target_set"].str.count(";") >= 1).all()
    assert fam["scoring_rule"].str.contains("any member").all()
    assert not ((cgt["feature"] == "ALLELE:KRAS_G12C") & (cgt["ko_class"] == "discordant-mechanism")).any()


def test_refuted_set(cgt):
    ref = cgt[cgt["partition"] == "clinically-refuted"]
    assert len(ref) >= 3


def test_proxy_rows(cgt):
    px = cgt[cgt["proxy"] == "1"]
    base = cgt.set_index("test_id")
    for _, r in px.iterrows():
        assert r["feature"].startswith("HOT:")
        assert r["proxy_of"] in base.index
        assert base.loc[r["proxy_of"], "feature"] not in NEVER_PROXIED


def test_testability_rule(cgt):
    tmin, cmin = int(S10["testable_min_carriers"]), int(S10["ceiling_min_carriers"])
    t = cgt[cgt["testable"] == "1"]
    assert (t["feature_D"] != "").all()
    assert (t["n_carriers_D"].astype(int) >= tmin).all()
    assert (t["n_carriers_comparator_D"].astype(int) >= cmin).all()
    u = cgt[cgt["testable"] == "0"]
    assert (u["untestable_reason"] != "").all()


def test_summary_consistent(cgt):
    s = json.loads(SUMMARY.read_text(encoding="utf-8"))
    head = cgt[cgt["proxy"] == "0"]
    assert s["n_tests"] == len(cgt)
    assert s["testable"] == int((head["testable"] == "1").sum())
    assert sum(s["by_partition"].values()) == len(head)
    assert s["checks"]["fusion_known_positives_pass"] is True


def test_alk_carriers_have_rule_partners():
    """PI ruling: the ALK <= 10 cap is a logged deviation; every ALK carrier must use a frozen rule partner."""
    s = json.loads(SUMMARY.read_text(encoding="utf-8"))
    rule = json.loads(RULES.read_text(encoding="utf-8"))["features"]["CLIN:ALK_FUS"]["partners"]
    car = s["checks"]["alk_carriers"]
    assert len(car) >= s["checks"]["alk_carriers_D"] > 0
    assert all(c["partners"] and set(c["partners"]) <= set(rule) for c in car), car
    if s["checks"]["alk_carriers_D"] > int(S10.get("alk_carriers_max", 10)):
        assert any(d["id"] == "DEV-CGT-1" for d in s["deviations"])


def test_kmt2a_rule_is_mechanistic():
    r = json.loads(RULES.read_text(encoding="utf-8"))["features"]["CLIN:KMT2A_R"]
    assert r["partners"] == "*" and r["driver"] == "KMT2A" and r["driver_end"] == "5p"
    s = json.loads(SUMMARY.read_text(encoding="utf-8"))
    assert all("KMT2A" not in c["partners"] for c in s["checks"]["kmt2a_carriers"])


def test_indications_from_config(cgt):
    im = pd.read_csv(ROOT / "configs" / "indication_map.csv", dtype=str, keep_default_na=False)
    assert {"indication_code", "indication_label", "match_field", "match_value", "source", "citation"} <= set(im.columns)
    assert (im["source"] != "").all() and (im["citation"] != "").all()
    codes = {c for s in cgt["indication_codes"] for c in s.split(";") if c}
    codes |= {t[4:] for s in cgt["comparator_exclude"] for t in s.split(";") if t.startswith("IND:")}
    assert codes <= set(im["indication_code"]), codes - set(im["indication_code"])
    assert "IND = {" not in BUILD.read_text(encoding="utf-8")


def test_er_concordance_reported():
    s = json.loads(SUMMARY.read_text(encoding="utf-8"))
    ec = s["checks"]["er_concordance"]
    assert "er_pos_vs_pam50_luminal" in ec and "er_pos_vs_subtype_ERplus" in ec
    assert len(s["recorded_gaps"]) == 3


def test_dropped_rows_have_reasons():
    if not DROPPED.exists():
        pytest.skip("no dropped file")
    d = pd.read_csv(DROPPED, dtype=str, keep_default_na=False)
    assert (d["drop_reason"] != "").all()


def test_build_reads_no_effects():
    src = BUILD.read_text(encoding="utf-8")
    forbidden = [r"load_effects", r"oof_predictions", r"test_predictions", r"component_oof", r"p1_test_metrics",
                 r"CRISPRGeneEffect", r"crispr\.parquet", r"raw_matrix\(\s*['\"]crispr", r"import_essml"]
    hits = [f for f in forbidden if re.search(f, src)]
    assert not hits, hits
