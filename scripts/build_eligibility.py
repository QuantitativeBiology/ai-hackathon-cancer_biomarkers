"""Eligibility (DESIGN v3 2): (variable OR selective) AND m5 tier in {high, medium} AND r_within >= 0.10.
Reads D (dev) effects only, via guard.load_effects(scope="dev"). Writes configs/eligible_genes.csv (entrez,
symbol, tier, r_beyond, eligible, reason, ...), results/eligibility.csv (same table) and
results/dictionary/eligibility_summary.json. Usage: python scripts/build_eligibility.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd  # noqa: E402

from biomark import eligibility, genesets  # noqa: E402
from biomark import paths as P  # noqa: E402


def main() -> int:
    cov = pd.read_parquet(P.CACHE / "covariates.parquet")
    el = eligibility.compute(cov)
    ctrl = set(genesets.control_dependency_genes()["entrez"].dropna())
    el["control_gene"] = el.index.isin(ctrl)
    cols = ["symbol", "tier", "r_beyond", "eligible", "reason", "variable", "selective", "n_lines_r", "control_gene"]
    out = el[cols].reset_index()
    out.to_csv(P.CONFIGS / "eligible_genes.csv", index=False, float_format="%.6g")
    out.to_csv(P.RESULTS / "eligibility.csv", index=False, float_format="%.6g")
    e = P.load_config("design_params")["eligibility_s2"]
    summ = {
        "rule": "(variable OR selective) AND m5 tier in {high, medium} AND r_within >= r_within_lin_lib_min",
        "r_within_lin_lib_min": e["r_within_lin_lib_min"], "stratum_min_lines": P.load_config("build_params")["eligibility"]["eligibility_min_stratum_lines"],
        "n_genes_evaluated": int(len(el)), "n_class_ok": int(el["class_ok"].sum()),
        "n_class_and_tier_ok": int((el["class_ok"] & el["tier_ok"]).sum()),
        "n_eligible": int(el["eligible"].sum()),
        "n_eligible_by_tier": {k: int(v) for k, v in el.loc[el["eligible"], "tier"].value_counts().items()},
        "reasons": {k: int(v) for k, v in el["reason"].value_counts().items()},
        "control_genes": int(el["control_gene"].sum()),
        "control_genes_eligible": int((el["control_gene"] & el["eligible"]).sum()),
        "median_r_beyond_class_tier_ok": float(el.loc[el["class_ok"] & el["tier_ok"], "r_beyond"].median()),
    }
    (P.RESULTS / "dictionary").mkdir(parents=True, exist_ok=True)
    (P.RESULTS / "dictionary" / "eligibility_summary.json").write_text(json.dumps(summ, indent=2), encoding="utf-8")
    print(json.dumps(summ, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
