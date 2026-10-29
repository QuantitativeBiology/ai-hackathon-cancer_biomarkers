"""Negative controls on real D data (DESIGN v3.2 §7.3) -> results/controls/nc_real.json (+ nc3_self_null.csv).

    python -B scripts/run_nc_real.py        # after run_arms.py --universe lowtier, run_cgt_secondary.py, run_limits.py

Statistics as defined before any value was seen (docs/toy_nc_bounds.md, real column of its table):
- NC1 FDR_emp_fixed per list (results/null/summary.json; bound NC1_fdr_emp_max).
- NC2a median EN-pred fidelity (lambda r2_oof_beyond_design) over the 100 matched-SD low-tier genes
  (configs/lowtier_genes.txt; the separate low-tier universe of run_arms.py --universe lowtier, results/benchmark_lowtier/);
  NC2b MD-strict rows of that universe / number of genes; MD-global rows / genes reported against the provisional
  0.1 for information (outcome_actions nc2b_lowtier_md_global).
- NC3 mean frac_dependent (audit.self_null: carrier median D effect < dependency threshold) over the testable
  self-null pairs (configs/dictionary_gene_features.csv relation self_null) of the eligible genes; D effects (guard
  scope dev). Per-pair rows in results/controls/nc3_self_null.csv.
- NC5 from results/clinical/cgt_secondary.json (run_cgt_secondary.py).
- NC6 from results/limits/limits_summary.json.
Bounds: configs/frozen_hashes.json bounds_values ('final'); labels: configs/outcome_actions.csv report_wording.
No threshold, bound or list is changed. Registry row. Refuses before freeze-bounds.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from biomark import analyses as AN  # noqa: E402
from biomark import endpoint_io as IO  # noqa: E402
from biomark import endpoints as EP  # noqa: E402
from biomark import paths as P  # noqa: E402

OUT = P.RESULTS / "controls" / "nc_real.json"
NC3_CSV = P.RESULTS / "controls" / "nc3_self_null.csv"


def bounds() -> dict:
    return json.loads((P.CONFIGS / "frozen_hashes.json").read_text(encoding="utf-8"))["bounds_values"]


def wording(check: str, failed: bool) -> str:
    oa = pd.read_csv(P.CONFIGS / "outcome_actions.csv", dtype=str)
    r = oa[oa["check"] == check]
    if r.empty:
        return ""
    strict1 = r["outcome"].str.match(r"^\s*(>|<)(?!=)")
    if len(r) == 1 and bool(strict1.iloc[0]) and not failed:
        return "passed"                                         # only a failure wording is pre-declared
    if len(r) > 1:
        strict = r["outcome"].str.match(r"^\s*(>|<)(?!=)")          # '> bound' / '< bound' rows = failed control
        r = r[strict] if failed else r[~strict]
    return str(r["report_wording"].iloc[0]) if len(r) else ""


def nc2(lowtier_dir: Path, md_dir: Path, genes: list) -> dict:
    lam = pd.read_csv(lowtier_dir / "lambda.csv", dtype={"gene": str})
    fid = lam[lam["arm"] == "EN-pred"].set_index("gene")["r2_oof_beyond_design"].reindex(genes)
    ms = pd.read_csv(md_dir / "md_strict.csv", dtype={"gene": str})
    mg = pd.read_csv(md_dir / "md_global.csv", dtype={"gene": str})
    n = len(genes)
    return {"NC2a": {"value": float(np.nanmedian(fid.to_numpy())) if fid.notna().any() else None,
                     "n_genes_with_fidelity": int(fid.notna().sum())},
            "NC2b": {"value": len(ms) / n, "md_strict_rows": int(len(ms)), "n_genes": n},
            "NC2b_md_global": {"value": len(mg) / n, "md_global_rows": int(len(mg)), "n_genes": n}}


def nc3(sn: pd.DataFrame) -> dict:
    t = sn[sn["testable"].astype(bool)]
    return {"value": float(t["frac_dependent"].mean()) if len(t) else None, "n_testable_pairs": int(len(t)),
            "n_dependent": int(np.nansum(t["frac_dependent"].to_numpy())), "n_pairs": int(len(sn))}


def self_null_table() -> pd.DataFrame:
    from biomark import audit as AU
    from biomark import guard, p1
    el = pd.read_csv(P.CONFIGS / "eligible_genes.csv", dtype={"entrez": str})
    elig = set(el.loc[el["eligible"].astype(bool), "entrez"])
    gf = pd.read_csv(P.CONFIGS / "dictionary_gene_features.csv", dtype={"entrez": str})
    sn = gf[(gf["relation"] == "self_null") & gf["entrez"].isin(elig)][["entrez", "feature_id"]]
    sn = sn.rename(columns={"entrez": "gene", "feature_id": "feature"}).drop_duplicates().reset_index(drop=True)
    dev = list(p1.dev_ids())
    X = IO.line_values(dev, sorted(set(sn["feature"])))
    return AU.self_null(sn, lambda gs: guard.load_effects(dev, scope="dev", genes=list(gs)),
                        lambda f: X[f].to_numpy(dtype=np.float64) if f in X else np.full(len(dev), np.nan))


def assemble(n2: dict, n3: dict, sec: dict, null: dict, lim: dict, b: dict) -> dict:
    import re

    def entry(v, bound, kind, check, name=""):
        fail = None if v is None else bool(v > bound if kind == "max" else v < bound)
        lab = wording(check, bool(fail)) if fail is not None else ""
        lab = re.sub(r"\bv\b", f"{v:.4g}", lab.replace("<list>", name)) if (lab and v is not None) else lab
        return {"value": v, "bound": bound, "bound_kind": kind, "pass": None if fail is None else not fail,
                "label": lab}

    uncal = ("S1-pred", "EN-pred selected", "MD-strict")             # outcome_actions: 'not calibrated (A6/A7)'
    nc1 = {k: entry(v, b["NC1_fdr_emp_max"]["final"], "max", "nc1_fdr_emp_fixed", k)
           for k, v in null["fdr_emp_fixed"].items()}
    for k in nc1:
        if k in uncal:
            oa = pd.read_csv(P.CONFIGS / "outcome_actions.csv", dtype=str)
            nc1[k]["label"] = oa.loc[oa["outcome"].str.startswith("any value"), "report_wording"].iloc[0]
    out = {"NC1": nc1,
           "NC2a": dict(entry(n2["NC2a"]["value"], b["NC2a_lowtier_median_fidelity_max"]["final"], "max",
                              "nc2a_lowtier_fidelity"), **{k: v for k, v in n2["NC2a"].items() if k != "value"}),
           "NC2b": dict(entry(n2["NC2b"]["value"], b["NC2b_lowtier_md_strict_per_gene_max"]["final"], "max",
                              "nc2b_lowtier_md_strict"), **{k: v for k, v in n2["NC2b"].items() if k != "value"}),
           "NC2b_md_global_information": dict(n2["NC2b_md_global"],
                                              provisional=b["NC2b_lowtier_md_strict_per_gene_max"]["provisional"],
                                              label=wording("nc2b_lowtier_md_global", False)),
           "NC3": dict(entry(n3["value"], b["NC3_self_null_dependent_frac_max"]["final"], "max", "nc3_self_null"),
                       **{k: v for k, v in n3.items() if k != "value"}),
           "NC5": dict(entry(sec["NC5"]["value"], b["NC5_false_sensitivity_rate_max"]["final"], "max",
                             "nc5_false_sensitivity"), per_arm=sec["NC5"]["per_arm"], n_tests=sec["NC5"]["n_tests"]),
           "NC6": dict(entry(lim["NC6_lineage_tf_class3"]["frac_class3"], b["NC6_lineage_tf_class3_frac_min"]["final"],
                             "min", "nc6_lineage_tf_class3"), n=lim["NC6_lineage_tf_class3"]["n"]),
           "NC4": {"value": None, "label": wording("nc4_fold_artefact", False), "note": "toy-only control (A7-iii)"}}
    return out


def run(check=None) -> dict:
    EP.require_descends("freeze-bounds", check=check)
    t0 = AN.utc_now()
    genes = [x.strip() for x in (P.CONFIGS / "lowtier_genes.txt").read_text(encoding="utf-8").splitlines() if x.strip()]
    n2 = nc2(IO.benchmark_dir("lowtier"), IO.md_dir("lowtier"), genes)
    sn = self_null_table()
    NC3_CSV.parent.mkdir(parents=True, exist_ok=True)
    sn.to_csv(NC3_CSV, index=False)
    sec = IO.read_json(P.RESULTS / "clinical" / "cgt_secondary.json")
    null = IO.read_json(P.RESULTS / "null" / "summary.json")
    lim = IO.read_json(P.RESULTS / "limits" / "limits_summary.json")
    out = assemble(n2, nc3(sn), sec, null, lim, bounds())
    out["_meta"] = {"bounds": "configs/frozen_hashes.json bounds_values (final)", "labels": "configs/outcome_actions.csv",
                    "deviation": "DEV-POST-R-COMPLETION (run after the R read on D data only)",
                    "lowtier_universe": "results/benchmark_lowtier/, results/md_lowtier/ (run_arms.py --universe lowtier)",
                    "start_utc": t0.isoformat(timespec="seconds"), "commit": __import__("biomark.gitutil", fromlist=["h"]).head(P.PROJECT)}
    IO.write_json(OUT, out)
    AN.append_registry("scripts/run_nc_real.py", "ok", ["results/controls/nc_real.json", "results/controls/nc3_self_null.csv"])
    return out


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    out = run()
    for k in ("NC2a", "NC2b", "NC3", "NC5", "NC6"):
        print(k, out[k]["value"], "bound", out[k]["bound"], "pass", out[k]["pass"], out[k]["label"])
    print("NC2b MD-global (information):", out["NC2b_md_global_information"]["value"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
