"""P1 secondaries from the single R run, without re-reading R (DESIGN §6.1, primary_endpoints.P1_primary secondary_k).

    python -B scripts/p1_secondary.py

Reads results/replication/r_tests.csv (the tests of the frozen R requests; `replicated` = same sign and one-sided
p < 0.05) and results/replication/bootstrap_p1.csv. Writes results/replication/p1_secondary.json:
  - secondary_k: replication rate per arm for the frozen k = 50 / 200 / 500 lists, with exact binomial 95% CIs;
  - bootstrap_diagnostic: per list, the observed rate against the bootstrap mean and 2.5-97.5% range, and the
    headline contrasts with the pre-declared percentile CI alongside a basic (reverse-percentile) CI. This is a
    reported diagnostic only; the P1 verdicts stay those of the pre-declared percentile rule.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from biomark import analyses as AN  # noqa: E402
from biomark import endpoints as EP  # noqa: E402
from biomark import paths as P  # noqa: E402

ARMS = ("S1-obs", "S1-pred", "EN-obs", "EN-pred")
CONTRASTS = (("EN-pred", "EN-obs"), ("EN-obs", "S1-obs"))


def exact_ci(k: int, n: int) -> list:
    lo = stats.beta.ppf(0.025, k, n - k + 1) if k > 0 else 0.0
    hi = stats.beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return [float(lo), float(hi)]


def main() -> int:
    EP.require_descends("freeze-catalogue")
    t_start = AN.utc_now()
    rt = pd.read_csv(P.RESULTS / "replication" / "r_tests.csv", dtype={"gene": str})
    rep = rt["replicated"].astype(str).str.lower().eq("true")
    reasons = rt["reasons"].fillna("").str.split("|")
    out = {"source": "results/replication/r_tests.csv (single R run; no R re-read)", "secondary_k": {}}
    for k in ("k50", "k200", "k500"):
        out["secondary_k"][k] = {}
        for arm in ARMS:
            m = reasons.apply(lambda r, tag=f"P1:{k}:{arm}": tag in r)
            n, h = int(m.sum()), int(rep[m].sum())
            if n:
                out["secondary_k"][k][arm] = {"n": n, "replicated": h, "rate": h / n, "ci95": exact_ci(h, n)}
    bt = pd.read_csv(P.RESULTS / "replication" / "bootstrap_p1.csv")
    obs = json.loads((P.RESULTS / "replication" / "p1_summary.json").read_text(encoding="utf-8"))["strata"]
    diag = {}
    for stratum, v in obs.items():
        d = {"lists": {}, "contrasts": {}}
        for arm in ARMS:
            col = f"{stratum}|{arm}"
            if col in bt:
                d["lists"][arm] = {"observed": v["rates"][arm], "boot_mean": float(bt[col].mean()),
                                   "boot_q025": float(bt[col].quantile(0.025)), "boot_q975": float(bt[col].quantile(0.975))}
        for a, b in CONTRASTS:
            ca, cb = f"{stratum}|{a}", f"{stratum}|{b}"
            if ca in bt and cb in bt:
                diff_obs = v["rates"][a] - v["rates"][b]
                dd = bt[ca] - bt[cb]
                q_lo, q_hi = float(dd.quantile(0.025)), float(dd.quantile(0.975))
                d["contrasts"][f"{a} - {b}"] = {"observed": diff_obs, "percentile_ci95": [q_lo, q_hi],
                                               "basic_ci95": [2 * diff_obs - q_hi, 2 * diff_obs - q_lo],
                                               "boot_mean": float(dd.mean())}
        diag[stratum] = d
    out["bootstrap_diagnostic"] = diag
    out["note"] = ("bootstrap rates sit below the observed rates (selection of strong effects into the frozen lists); "
                   "basic CIs are reported as a diagnostic only, verdicts follow the pre-declared percentile rule")
    f = P.RESULTS / "replication" / "p1_secondary.json"
    f.write_text(json.dumps(out, indent=1), encoding="utf-8")
    # timing kept next to the outputs: p1_secondary is not a tiers.json stage (stage_job_map is frozen)
    (P.RESULTS / "replication" / "timing_p1_secondary.json").write_text(json.dumps(
        {"job": "p1_secondary", "start_utc": str(t_start), "utc": str(AN.utc_now()), "workers": 1}), encoding="utf-8")
    AN.append_registry("scripts/p1_secondary.py", "ok", ["results/replication/p1_secondary.json"])
    print(json.dumps(out["secondary_k"], indent=1))
    print(json.dumps({s: d["contrasts"] for s, d in diag.items()}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
