"""run_null's observed inputs must carry MD-global, so observed_hits counts it (the first real null summary had
MD-global hits_obs = 0 because only md_strict was loaded). Uses a synthetic observed-input dict for the count and a
static check of the script for the loader."""
import ast
from pathlib import Path

import pandas as pd

from biomark import null as NL

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_null.py"


def test_observed_hits_counts_md_global_on_the_null_genes():
    arms = pd.DataFrame({"arm": ["S1-obs", "S1-pred", "EN-obs", "EN-pred"], "gene": ["1"] * 4, "feature": ["f"] * 4,
                         "selected": [True] * 4})
    md = pd.DataFrame({"gene": ["1", "1", "2"], "feature": ["a", "b", "c"]})
    obs = NL.observed_hits({"arms": arms, "md_strict": md, "md_global": md.iloc[:2]}, ["1"])
    assert obs["MD-strict"] == 2 and obs["MD-global"] == 2


def test_run_null_observed_inputs_loads_md_global():
    tree = ast.parse(SCRIPT.read_text(encoding="utf-8"))
    fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "observed_inputs")
    keys = {k.value for n in ast.walk(fn) if isinstance(n, ast.Dict) for k in n.keys if isinstance(k, ast.Constant)}
    assert {"md_strict", "md_global"} <= keys
