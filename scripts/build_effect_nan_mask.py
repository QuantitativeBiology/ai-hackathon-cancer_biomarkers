"""D effect NaN mask (DESIGN v3 1.6): cache/effect_nan_mask.parquet (D lines x genes, True = effect is NaN),
the only effect-derived input EN-pred may read. Reads D (dev) effects via guard.load_effects(scope="dev");
records sha256 in results/bulk_manifest.csv. Usage: python scripts/build_effect_nan_mask.py
"""
from __future__ import annotations

import datetime as dt
import hashlib
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd  # noqa: E402

from biomark import guard, p1  # noqa: E402
from biomark import paths as P  # noqa: E402


def main() -> int:
    dev = list(p1.dev_ids())
    y = guard.load_effects(dev, scope="dev")
    mask = y.isna()
    out = P.CACHE / "effect_nan_mask.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    mask.to_parquet(out)
    sha = hashlib.sha256(out.read_bytes()).hexdigest()
    man = P.RESULTS / "bulk_manifest.csv"
    row = pd.DataFrame([{"path": out.relative_to(P.PROJECT).as_posix(), "sha256": sha, "rows": mask.shape[0],
                         "cols": mask.shape[1], "written_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                         "script": "scripts/build_effect_nan_mask.py"}])
    if man.exists():
        old = pd.read_csv(man)
        row = pd.concat([old[old["path"] != row["path"].iloc[0]], row], ignore_index=True)
    row.to_csv(man, index=False)
    print(f"mask {mask.shape}, NaN cells {int(mask.to_numpy().sum())}, sha256 {sha[:12]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
