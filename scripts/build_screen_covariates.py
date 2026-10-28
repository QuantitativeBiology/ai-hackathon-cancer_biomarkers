"""Screen covariates of the confounder model (DESIGN v3 5.3) for D, R and set A from the screened condition.
Reads no effects. Writes cache/covariates.parquet, results/dictionary/covariate_summary.csv,
results/covariates/vif.csv and results/dictionary/covariate_checks.json.
Usage: python scripts/build_screen_covariates.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from biomark import covariates  # noqa: E402
from biomark import paths as P  # noqa: E402


def main() -> int:
    cov, info = covariates.build()
    P.CACHE.mkdir(parents=True, exist_ok=True)
    cov.to_parquet(P.CACHE / "covariates.parquet")
    out = P.RESULTS / "dictionary"
    out.mkdir(parents=True, exist_ok=True)
    covariates.summary(cov).to_csv(out / "covariate_summary.csv", index=False)
    (P.RESULTS / "covariates").mkdir(parents=True, exist_ok=True)
    info["vif"].rename("vif").to_csv(P.RESULTS / "covariates" / "vif.csv", index_label="column")
    chk = info["checks"]
    chk["n_lines"] = {k: int(v) for k, v in cov["role"].value_counts().items()}
    chk["columns_categorical"] = covariates.CATEGORICAL
    chk["columns_continuous"] = covariates.CONTINUOUS
    chk["columns_flags"] = covariates.FLAGS
    (out / "covariate_checks.json").write_text(json.dumps(chk, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in chk.items() if not isinstance(v, list)}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
