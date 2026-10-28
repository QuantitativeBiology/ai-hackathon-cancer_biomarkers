"""No raw/restricted data or large binaries are tracked in git (DepMap terms; DATA.md)."""
import os

from biomark import gitutil
from biomark import paths as P

MAX_BYTES = 5 * 1024 * 1024
DATA_SUFFIXES = (".parquet", ".feather", ".npy", ".npz", ".pkl", ".joblib", ".h5", ".hdf5", ".csv.gz",
                 ".tsv.gz", ".gz", ".zip", ".maf", ".xlsx")
RAW_NAMES = {"CRISPRGeneEffect.csv", "Model.csv", "ModelCondition.csv", "OmicsCNGeneWGS.csv",
             "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", "OmicsFusionFiltered.csv",
             "OmicsSomaticMutationsMatrixDamaging.csv", "OmicsSomaticMutationsMatrixHotspot.csv",
             "CRISPRScreenMap.csv", "ScreenSequenceMap.csv"}


def test_no_large_tracked_files():
    big = [f for f in gitutil.ls_files(P.PROJECT)
           if (P.PROJECT / f).exists() and os.path.getsize(P.PROJECT / f) > MAX_BYTES]
    assert big == []


def test_no_data_files_tracked():
    files = gitutil.ls_files(P.PROJECT)
    bad = [f for f in files if f.lower().endswith(DATA_SUFFIXES) or os.path.basename(f) in RAW_NAMES
           or f.startswith(("data/", "cache/", "external/"))]
    assert bad == []
