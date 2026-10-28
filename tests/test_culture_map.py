"""DESIGN §13.2 "culture map" (§5.3 culture block: screened CellFormat {adherent, suspension, 3D, mixed/unknown},
SerumFree {T, F, unknown}, base-medium family; the formulation column is picked by content because Model.csv
OnboardedMedia / FormulationID appear swapped; lines without media -> unknown). Maps live in
configs/build_params.json:covariates; code in src/biomark/covariates.py. No effects are read.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from biomark import covariates as CV
from biomark import paths as P

CFG = P.load_config("build_params")["covariates"]


@pytest.mark.parametrize("text, fam", [
    ("RPMI 1640 + 10% FBS", "RPMI"), ("RPMI-1640", "RPMI"), ("RPMI", "RPMI"), ("rpmi 1640 + 2mM glutamine", "RPMI"),
    ("DMEM:F12 + 10% FBS", "DMEM-F12"), ("DMEM/F-12", "DMEM-F12"), ("DMEM-F12", "DMEM-F12"),
    ("DMEM + 10% FBS", "DMEM"), ("DMEM high glucose", "DMEM"), ("IMDM + 20% FBS", "IMDM"),
    ("McCoy's 5A + 10% FBS", "other"), ("L-15", "other"), ("RPMI + DMEM", "RPMI"),
])
def test_medium_family_of(text, fam):
    assert CV.medium_family_of(text) == fam


@pytest.mark.parametrize("empty", [None, np.nan, "", "   "])
def test_missing_medium_is_none(empty):
    assert CV.medium_family_of(empty) is None


def test_combine_conditions():
    assert CV._combine([None, None]) is None
    assert CV._combine(["RPMI", None, "RPMI"]) == "RPMI"
    assert CV._combine(["RPMI", "DMEM"]) == CFG["medium_multi_family"]


def test_formulation_column_picked_by_content():
    df = pd.DataFrame({"OnboardedMedia": ["MF-001-041", "MF-002-001", "MF-003-001", None],
                       "FormulationID": ["RPMI + 10% FBS", "DMEM + 10% FBS", "IMDM + 20% FBS", None]})
    assert CV._detect_formulation_col(df, ["OnboardedMedia", "FormulationID"]) == "FormulationID"
    df2 = df.rename(columns={"OnboardedMedia": "FormulationID", "FormulationID": "OnboardedMedia"})
    assert CV._detect_formulation_col(df2, ["OnboardedMedia", "FormulationID"]) == "OnboardedMedia"
    assert CV._detect_formulation_col(df[["OnboardedMedia"]], ["OnboardedMedia"]) is None


def test_maps_cover_design_levels():
    assert set(CFG["cellformat_map"].values()) | {CFG["cellformat_missing"], CFG["cellformat_multi"]} == \
        {"adherent", "suspension", "3D", "mixed/unknown"}
    assert set(CFG["serumfree_map"].values()) | {CFG["serumfree_missing"]} == {"T", "F", "unknown"}


@pytest.mark.skipif(not (P.CACHE / "covariates.parquet").exists(), reason="covariate cache not built")
def test_built_covariate_levels_and_sources():
    c = pd.read_parquet(P.CACHE / "covariates.parquet")
    fams = {f for f, _ in CFG["medium_family_regex"]} | {CFG["medium_other"], CFG["medium_unknown"]}
    assert set(c["cellformat"]) <= {"adherent", "suspension", "3D", "mixed/unknown"}
    assert set(c["serumfree"]) <= {"T", "F", "unknown"}
    assert set(c["medium_family"]) <= fams
    assert c[["cellformat", "serumfree", "medium_family"]].notna().all().all()
    # lines without media -> unknown, and only those
    assert ((c["medium_source"] == "none") == (c["medium_family"] == CFG["medium_unknown"])).all()
    assert ((c["serumfree_source"] == "none") == (c["serumfree"] == CFG["serumfree_missing"])).all()
    assert set(c["role"]) == {"dev", "test", "noncancerous"}       # D ∪ R ∪ set A (§5.3)
