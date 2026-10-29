"""CGT secondary metrics on D (DESIGN v3.2 §10.5, §10.6, §7.3 NC5) after the data ceiling.

    python -B scripts/run_cgt_secondary.py      # after run_q4.py and run_ceiling.py

Writes results/clinical/miss_decomposition.csv (one row per missed (arm, test, set, subset, k) of p2_scores.csv) and
results/clinical/cgt_secondary.json:
- **miss_decomposition** per arm and set|subset|k: every miss gets the first applicable category of
  clinical_s10.miss_decomposition, checked in this order: (1) ineligible gene (class 1/2/3; no target gene eligible;
  the classes of the targets from results/limits/gene_classes.csv); (2) no feature_D or feature not in 𝒟_g (no
  eligible target has (g, feature_D) in U); (3) ceiling data-null (ceiling label data-null or data-contradicted,
  results/clinical/cgt_ceiling.csv; the label is kept in a column); (4) lost_to_group_member (the best target's
  feature_D row is not its variant group's de-duplicated representative); (5) eligible, ranked below k. Order (3)
  before (4)-(5) puts the method-independent reason first (agent ordering; the design lists the categories only).
- **recall over data-supported tests** per arm and set|subset|k (testable tests with ceiling label data-supported).
- **NC5** false-sensitivity rate: per arm, the CGT tests with partition bypass-resistance or sign N whose feature_D is
  ranked in the arm (docs/toy_nc_bounds.md); false sensitivity = hit (rank_of_feature_D <= k = 200) with beta < 0;
  statistic = max over the four arms; bound = configs/frozen_hashes.json bounds_values NC5 (final).
- **recall_refuted vs recall_approved** (§10.6) per arm at k = 200 with exact CIs: refuted = clinically-refuted tests
  (expected KO sign S); approved = T1 sensitivity non-proxy tests; both restricted to testable tests (feature_D in U
  with >= testable_min_carriers carriers for some target); hit = rank_of_feature_D <= k, any sign.
D lists only (no R). Registry row; no timing stem. Refuses before freeze-bounds.
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
from biomark import arms as A  # noqa: E402
from biomark import endpoint_io as IO  # noqa: E402
from biomark import endpoints as EP  # noqa: E402
from biomark import params as PR  # noqa: E402
from biomark import paths as P  # noqa: E402

OUT_CSV = P.RESULTS / "clinical" / "miss_decomposition.csv"
OUT_JSON = P.RESULTS / "clinical" / "cgt_secondary.json"
CATS = list(PR.dp()["clinical_s10"]["miss_decomposition"])


def decompose(scores: pd.DataFrame, cgt: pd.DataFrame, eligible: set, gene_class: dict, U: set, ceil: dict,
              lost: dict) -> pd.DataFrame:
    """Category of every missed (arm, test) row of ``scores`` (p2_scores layout; ``cgt`` has a ``targets`` column)."""
    T = cgt.set_index("test_id")
    rows = []
    for _, s in scores[~scores["hit"].astype(bool)].iterrows():
        t = T.loc[s["test_id"]]
        targets = [str(g) for g in t["targets"]]
        el = [g for g in targets if g in eligible]
        fD = t["feature_D"] if isinstance(t["feature_D"], str) and t["feature_D"] else None
        lab = ceil.get(s["test_id"], "")
        if not el:
            cat = CATS[0]
            detail = "classes " + "/".join(str(gene_class.get(g, "n/a")) for g in targets)
        elif fD is None or not any((g, fD) in U for g in el):
            cat, detail = CATS[1], ("no feature_D" if fD is None else "feature not in U for an eligible target")
        elif lab in ("data-null", "data-contradicted"):
            cat, detail = CATS[4], f"ceiling {lab}"
        elif str(s.get("gene", "")) and lost.get((s["arm"], str(s["gene"]), fD)):
            cat, detail = CATS[3], f"lost to {lost[(s['arm'], str(s['gene']), fD)]}"
        else:
            cat, detail = CATS[2], f"rank {s['rank_of_feature_D']}"
        rows.append({"arm": s["arm"], "set": s["set"], "subset": s["subset"], "k": int(s["k"]), "test_id": s["test_id"],
                     "testable": bool(s["testable"]), "category": cat, "detail": detail, "ceiling_label": lab})
    return pd.DataFrame(rows, columns=["arm", "set", "subset", "k", "test_id", "testable", "category", "detail",
                                       "ceiling_label"])


def nc5(scores: pd.DataFrame, k: int) -> dict:
    """``scores``: p2_score rows of the NC5 tests; false sensitivity = hit at k with beta < 0, over ranked tests."""
    per = {}
    for arm, s in scores.groupby("arm"):
        ranked = s[np.isfinite(s["rank_of_feature_D"].astype(float))]
        fs = ranked[(ranked["rank_of_feature_D"] <= k) & (ranked["beta"] < 0)]
        per[arm] = {"n_ranked": int(len(ranked)), "false_sensitivity": int(len(fs)),
                    "rate": (len(fs) / len(ranked)) if len(ranked) else None, "tests": sorted(fs["test_id"])}
    rates = [v["rate"] for v in per.values() if v["rate"] is not None]
    return {"value": max(rates) if rates else None, "per_arm": per}


def recall(scores: pd.DataFrame) -> dict:
    out = {}
    for arm, s in scores.groupby("arm"):
        n = int(s["test_id"].nunique())
        h = int(s.loc[s["hit"].astype(bool), "test_id"].nunique())
        out[arm] = {"hits": h, "n": n, "recall": h / n if n else None, "ci95": list(EP.clopper_pearson(h, n))}
    return out


def testable_flags(cgt: pd.DataFrame, Uc: dict, nmin: int) -> list:
    ok = []
    for ts, f in zip(cgt["targets"], cgt["feature_D"]):
        good = isinstance(f, str) and bool(f) and any((nc := Uc.get((str(g), f))) is not None and
                                                     (pd.isna(nc) or nc >= nmin) for g in ts)
        ok.append(bool(good))
    return ok


def run(check=None) -> dict:
    EP.require_descends("freeze-bounds", check=check)
    from biomark import genesets as GS
    pc = PR.endpoints()["P2_primary_clinical"]
    k = int(pc["k"])
    cgt = pd.read_csv(P.CONFIGS / "clinical_ground_truth.csv")
    ez = GS.symbol_to_entrez(sorted({x for v in cgt["target_set"] for x in EP._targets(v)})).dropna().astype(str).to_dict()
    cgt["targets"] = [[str(ez.get(x, x)) for x in EP._targets(v)] for v in cgt["target_set"]]
    tgenes = sorted({g for ts in cgt["targets"] for g in ts})
    Udf = IO.read_arms(arms=["S1-obs"], columns=["gene", "feature", "n_carriers"], genes=tgenes)
    Uc = dict(zip(zip(Udf["gene"].astype(str), Udf["feature"]), Udf["n_carriers"]))
    rows = IO.read_arms(columns=["arm", "gene", "feature", "rank", "rank_dedup", "variant_group", "beta",
                                 "lost_to_group_member"], genes=tgenes)
    dedup = {}
    for arm in A.ARMS:
        r = IO.read_arms(arms=[arm], columns=["rank", "rank_dedup"])
        dedup[arm] = r.loc[r["rank_dedup"].notna(), "rank"].to_numpy()
        del r
    ridx = EP.RankIndex.slim(dedup, rows.drop(columns=["lost_to_group_member"]))
    lost = {(a, str(g), f): x for a, g, f, x in zip(rows["arm"], rows["gene"], rows["feature"], rows["lost_to_group_member"])
            if isinstance(x, str) and x}
    el = pd.read_csv(P.CONFIGS / "eligible_genes.csv", dtype={"entrez": str})
    eligible = set(el.loc[el["eligible"].astype(bool), "entrez"])
    gc = pd.read_csv(P.RESULTS / "limits" / "gene_classes.csv", dtype={"gene": str, "failure_class": str})
    gene_class = dict(zip(gc["gene"], gc["failure_class"]))
    ce = pd.read_csv(P.RESULTS / "clinical" / "cgt_ceiling.csv", dtype={"test_id": str})
    ceil = dict(zip(ce["test_id"], ce["label"]))
    scores = pd.read_csv(P.RESULTS / "clinical" / "p2_scores.csv", dtype={"gene": str})
    scores["gene"] = scores["gene"].fillna("").astype(str)
    md = decompose(scores, cgt, eligible, gene_class, set(Uc), ceil, lost)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    md.to_csv(OUT_CSV, index=False)
    summ = {"miss_decomposition": {}, "recall_data_supported": {}, "k": k}
    for (st, sub, kk), g in scores.groupby(["set", "subset", "k"]):
        key = f"{st}|{sub}|k{int(kk)}"
        m = md[(md["set"] == st) & (md["subset"] == sub) & (md["k"] == int(kk))]
        summ["miss_decomposition"][key] = {
            arm: {"n_misses": int((~x["hit"].astype(bool)).sum()),
                  "by_category": m[m["arm"] == arm]["category"].value_counts().reindex(CATS, fill_value=0).to_dict()}
            for arm, x in g.groupby("arm")}
        ds = g[g["testable"].astype(bool) & g["test_id"].map(lambda t: ceil.get(t) == "data-supported")]
        summ["recall_data_supported"][key] = recall(ds) if len(ds) else {}
    for key, v in summ["miss_decomposition"].items():
        for arm, d in v.items():
            assert sum(d["by_category"].values()) == d["n_misses"], (key, arm)
    nmin = int(PR.dp()["clinical_s10"]["testable_min_carriers"])
    cgt["testable_U"] = testable_flags(cgt, Uc, nmin)
    nc5_tests = cgt[(cgt["partition"] == "bypass-resistance") | (cgt["sign"] == "N")]
    sc5 = EP.p2_score(None, nc5_tests.assign(testable=nc5_tests["testable_U"]), k=k, ridx=ridx, arms=A.ARMS)
    bounds = json.loads((P.CONFIGS / "frozen_hashes.json").read_text(encoding="utf-8"))["bounds_values"]
    n5 = nc5(sc5, k)
    n5.update({"bound_max": bounds["NC5_false_sensitivity_rate_max"]["final"], "n_tests": int(len(nc5_tests))})
    n5["pass"] = None if n5["value"] is None else bool(n5["value"] <= n5["bound_max"])
    summ["NC5"] = n5
    ref = cgt[(cgt["partition"] == "clinically-refuted") & cgt["testable_U"]]
    appr = cgt[(cgt["tier"] == "T1") & (cgt["partition"] == "sensitivity") &
               (pd.to_numeric(cgt["proxy"], errors="coerce").fillna(0) == 0) & cgt["testable_U"]]
    sr = EP.p2_score(None, ref.assign(testable=True), k=k, ridx=ridx, arms=A.ARMS)
    sa = EP.p2_score(None, appr.assign(testable=True), k=k, ridx=ridx, arms=A.ARMS)
    summ["refuted_vs_approved"] = {"k": k, "n_refuted": int(len(ref)), "n_approved": int(len(appr)),
                                   "refuted_tests": sorted(ref["test_id"]), "recall_refuted": recall(sr),
                                   "recall_approved": recall(sa)}
    summ["rules"] = {"miss_order": "ineligible > no feature > ceiling data-null/contradicted > lost_to_group_member > "
                                   "ranked below k (agent ordering of the design's categories)",
                     "nc5": "docs/toy_nc_bounds.md NC5 analogue on the real CGT table and D lists",
                     "lists": "D only (no R)"}
    IO.write_json(OUT_JSON, summ)
    AN.append_registry("scripts/run_cgt_secondary.py", "ok",
                       ["results/clinical/miss_decomposition.csv", "results/clinical/cgt_secondary.json"])
    return summ


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    s = run()
    print(json.dumps({"NC5": {k: s["NC5"][k] for k in ("value", "bound_max", "pass")},
                      "refuted_vs_approved": {a: (s["refuted_vs_approved"]["recall_refuted"].get(a, {}).get("recall"),
                                                  s["refuted_vs_approved"]["recall_approved"].get(a, {}).get("recall"))
                                              for a in A.ARMS}}, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
