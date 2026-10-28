"""configs/positive_controls.csv and results/controls/pc_eligibility.csv (DESIGN v3 9) from
configs/positive_controls_spec.csv, configs/dictionary_v1.json and the dictionary caches. Fails on any
unresolvable feature id (R3-2). Usage: python scripts/build_positive_controls.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from biomark import controls  # noqa: E402
from biomark import p1  # noqa: E402,F401
from biomark import paths as P  # noqa: E402

FIRST = ["control_id", "family", "dependency_symbol", "dependency_entrez", "feature_primary", "feature_secondary",
         "fallback_feature", "expected_sign", "scope", "citation", "in_cgt", "n_carriers_D", "n_carriers_R",
         "access_class", "en_max_attainable_freq", "gate_eligible", "note"]


def main() -> int:
    pc = controls.build()
    cols = FIRST + [c for c in pc.columns if c not in FIRST]
    pc = pc[cols]
    pc.to_csv(P.CONFIGS / "positive_controls.csv", index=False)
    out = P.RESULTS / "controls"
    out.mkdir(parents=True, exist_ok=True)
    pc[["control_id", "family", "dependency_symbol", "gate_feature", "scope", "n_obs_D", "n_carriers_D",
        "n_carriers_R", "gate_eligible", "access_class", "en_max_attainable_freq"]].to_csv(out / "pc_eligibility.csv", index=False)
    g = pc[pc["scope"] != "class_check"]
    summ = {"n_rows": int(len(pc)), "n_gate_eligible": int(g["gate_eligible"].sum()),
            "n_untestable": int((g["scope"] == "untestable").sum()),
            "gate_eligible_by_family": {k: int(v) for k, v in g.groupby("family")["gate_eligible"].sum().items()},
            "untestable": list(g.loc[g["scope"] == "untestable", "control_id"])}
    (out / "pc_summary.json").write_text(json.dumps(summ, indent=1), encoding="utf-8")
    print(json.dumps(summ, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
