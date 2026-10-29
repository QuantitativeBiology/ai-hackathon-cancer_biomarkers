"""Limits: gene classes (DESIGN v3.2 §12.1) and pair classes (§12.2) on the catalogue (amendment A7: MD-global;
MD-strict alongside), extra labels, Q1b (§8.3) and the NC6 check.

    python -B scripts/run_limits.py              # real (after run_confounder_audit, run_m5_audits, run_cross_lineage)
    python -B scripts/run_limits.py --toy DIR    # toy workspace

Inputs (frozen outputs): results/md/*.csv, results/benchmark/lambda.csv, results/confounders/delta.csv,
results/locus/{resolution,amplicon}.csv, results/subtypes/{culture_audit,hypermut_tp53}.csv,
results/lineage/cross_lineage.csv, configs/catalogue/r_requests.csv (pre-tag R power, optional),
results/replication/r_tests.csv (optional: R-dependent classes are 'pending R' until it exists),
results/clinical/cgt_ceiling.csv (optional), configs/clinical_ground_truth.csv, configs/eligible_genes.csv, Part 1
gene_classes (D SD, common-essential flag). No effects are read.
Pair classes: 4 expression gating only; 5 too rare or R power < power_min; 6 confounded (confounding delta,
subtype-confounded, unresolved co-location, cutting artefact, collinear_with_design, culture flag, hypermutation,
TP53); 7 powered but not replicated in R; 8 KO-drug discordance; 9 valid. Gene classes 1 / 2 / 3 / 3b, then 9 if any
pair is 9, else the furthest pair class.
Writes results/limits/{pair_classes.csv, gene_classes.csv, pair_labels.csv, limits_summary.json}; timing stem limits.
Refuses unless HEAD descends from freeze-bounds.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from biomark import analyses as AN  # noqa: E402
from biomark import limits as LM  # noqa: E402
from biomark import m5 as M5  # noqa: E402
from biomark import params as PR  # noqa: E402
from biomark import paths as P  # noqa: E402

SCRIPT = "scripts/run_limits.py"
OUT = "results/limits"


def nc6_bound() -> float:
    import json
    fh = json.loads((P.CONFIGS / "frozen_hashes.json").read_text(encoding="utf-8"))
    v = fh.get("bounds_values", {}).get("NC6_lineage_tf_class3_frac_min", {})
    return float(v.get("final", PR.dp()["negative_controls_s7"]["NC6_lineage_tf_class3_frac_min"]))


def run(ws, src) -> dict:
    t0 = AN.utc_now()
    L = M5.cfg()["limits"]
    cat = M5.catalogue(ws)
    fd = LM.FeatureData(src.fmeta(), src.symbols())
    delta = ws.read_csv("results/confounders/delta.csv")
    locus = ws.read_csv("results/locus/resolution.csv")
    amp = ws.read_csv("results/locus/amplicon.csv")
    culture = ws.read_csv("results/subtypes/culture_audit.csv")
    hyper = ws.read_csv("results/subtypes/hypermut_tp53.csv")
    cl = ws.read_csv("results/lineage/cross_lineage.csv", required=False)
    req = ws.read_csv("configs/catalogue/r_requests.csv", required=False)
    rt = ws.read_csv("results/replication/r_tests.csv", required=False)
    ceil = ws.read_csv(L["ceiling_file"], required=False, test_id=str)
    attr = ws.read_csv(L["attribution_file"], required=False)
    power = LM._key(req, "power")
    cgt = src.cgt()
    gi = src.gene_info()
    ctl = src.controls()
    universe = LM.universe_genes(gi, ctl)
    rpow = ws.read_csv("results/replication/r_power.csv", required=False)
    pc = LM.pair_classes(fd, cat, delta, locus, amp, culture, hyper, power, rt, cgt, ceil, gene_universe=universe,
                         r_power=rpow)
    gc = LM.gene_classes(universe, gi, ws.lam(), pc, attr, ctl)
    s1h = ws.arms(arms=["S1-obs"], columns=["gene", "feature"], selected_only=True)
    en_sel = ws.arms(arms=["EN-pred"], columns=["gene", "feature"], selected_only=True)
    labels = LM.pair_labels(fd, cat, s1h, en_sel, gi, rt)
    outs = [ws.write_csv(pc, f"{OUT}/pair_classes.csv"), ws.write_csv(gc, f"{OUT}/gene_classes.csv"),
            ws.write_csv(labels, f"{OUT}/pair_labels.csv")]
    lc = ctl[ctl["scope"] == "class_check"]["gene"].unique() if len(ctl) else []
    c3 = gc.set_index("gene")["failure_class"].reindex(lc)
    pg = pc[pc["in_md_global"].astype(bool)]
    summ = {"catalogue": "MD-global (A7); MD-strict alongside (failure_class_md_strict, in_md_strict)",
            "R_status": "available" if rt is not None else L["pending"],
            "power_source": "r_tests.csv" if rt is not None else ("r_requests.csv (pre-tag)" if req is not None
                                                                   else "missing: class 5 by power pending"),
            "ceiling": "available" if ceil is not None else "cgt_ceiling.csv missing: class 8 from KO class / sign only",
            "m1_attribution": "available" if attr is not None else "not run (Tier B); class 3 uses r_within and fidelity",
            "pair_classes_md_global": pg["pair_class"].value_counts().sort_index().to_dict() if len(pg) else {},
            "pair_classes_md_strict": pc[pc["in_md_strict"].astype(bool)]["pair_class"].value_counts().sort_index().to_dict(),
            "pair_class_status": pc["class_status"].value_counts().to_dict() if len(pc) else {},
            "gene_classes": gc["failure_class"].value_counts().to_dict(),
            "gene_classes_md_strict": gc["failure_class_md_strict"].value_counts().to_dict(),
            "n_universe_genes": len(universe),
            "NC6_lineage_tf_class3": {"n": int(len(lc)), "frac_class3": float((c3 == "3").mean()) if len(lc) else np.nan,
                                     "bound_min": nc6_bound(),
                                     "pass": bool(len(lc) and (c3 == "3").mean() >= nc6_bound())},
            "Q1b": LM.q1b(cl, cat, cgt, rt),
            "methylation_coverage_D_lines": src.meth_profiled(),
            "extra_labels": {"processing_limited": int(pd.to_numeric(labels["processing_limited"], errors="coerce")
                                                        .fillna(0).astype(bool).sum()),
                             "predictor_blind_rows": int((labels["predictor_blind"] != "").sum())}}
    outs.append(ws.write_json(summ, f"{OUT}/limits_summary.json"))
    ws.timing("limits", {"n_pairs": int(len(pc)), "n_genes": int(len(gc))}, t0, workers=1)
    ws.register(SCRIPT, outs)
    return {"pair_classes": pc, "gene_classes": gc, "labels": labels, "summary": summ}


def main(argv=None) -> int:
    a = M5.parse_args(argv, __doc__)
    ws, src = M5.open_workspace(a)
    r = run(ws, src)
    print(f"limits: pair classes {r['summary']['pair_classes_md_global']}; gene classes {r['summary']['gene_classes']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
