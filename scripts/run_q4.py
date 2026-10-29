"""P2 / Q4 clinical scoring on D (DESIGN v3.2 §10.5; primary_endpoints.P2_primary_clinical; amendment A4).

    python -B scripts/run_q4.py

Q4a / P2: recall of headline CGT tests per arm at k = 200 (secondary 100 / 500) over U_all, sets non_proxy (headline)
and incl_proxy, subsets all and gate_disjoint (headline when results/gate_iterations.csv has > 1 row); exact McNemar
EN-pred vs S1-obs; descriptive below min_testable_tests; per-tier recall; direction accuracy and sign errors.
Q4b actionability (MD-strict and S1-obs vs within-list gene-label permutations). Q4c population: within-indication
and per non-indicated lineage β, test labels, population precision and the carriers-only OR of accurate-dependent
membership (indicated vs not). Indications map to lines through OncotreeCode and, for parent codes, the OncoTree
term paired with the code in the frozen MOAlmanac release matched against Model.csv Subtype / PrimaryDisease /
Lineage names; unmapped codes are listed. D effects only. Writes results/clinical/{q4_summary.json, p2_scores.csv,
actionability.csv, q4c_population.csv, q4c_labels.csv}; timing stem q4_summary. Refuses before freeze-bounds.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from biomark import analyses as AN  # noqa: E402
from biomark import eligibility as ELG  # noqa: E402
from biomark import endpoint_io as IO  # noqa: E402
from biomark import endpoints as EP  # noqa: E402
from biomark import genesets as GS  # noqa: E402
from biomark import guard  # noqa: E402
from biomark import p1  # noqa: E402
from biomark import params as PR  # noqa: E402
from biomark import paths as P  # noqa: E402


def main() -> int:
    EP.require_descends("freeze-bounds")
    t0 = time.perf_counter()
    pc = PR.endpoints()["P2_primary_clinical"]
    cgt = pd.read_csv(P.CONFIGS / "clinical_ground_truth.csv")
    ez = GS.symbol_to_entrez(sorted({x for v in cgt["target_set"] for x in EP._targets(v)})).dropna().astype(str).to_dict()
    tgenes = sorted({str(ez.get(x, x)) for v in cgt["target_set"] for x in EP._targets(v)})
    # slim reads (the arms table is ~26 M rows): U and drops for the CGT target genes, ranks via RankIndex.slim
    U = IO.read_arms(arms=["S1-obs"], columns=["gene", "feature", "n_carriers", "dropped_cov"], genes=tgenes)
    rows = IO.read_arms(columns=["arm", "gene", "feature", "rank", "rank_dedup", "variant_group", "beta"], genes=tgenes)
    dedup = {}
    for arm in ("S1-obs", "S1-pred", "EN-obs", "EN-pred"):
        r = IO.read_arms(arms=[arm], columns=["rank", "rank_dedup"])
        dedup[arm] = r.loc[r["rank_dedup"].notna(), "rank"].to_numpy()
        del r
    it = P.RESULTS / "gate_iterations.csv"
    gd_headline = it.exists() and len(pd.read_csv(it)) > 1
    ridx = EP.RankIndex.slim(dedup, rows)
    summ = {"headline_subset": "gate_disjoint" if gd_headline else "all", "results": {}}
    scores_all = []
    for set_ in ("non_proxy", "incl_proxy"):
        for subset in ("all", "gate_disjoint"):
            T = EP.p2_test_set(cgt, U, set_, gate_disjoint=subset == "gate_disjoint", entrez_of=ez)
            for k in [pc["k"]] + list(pc["secondary_k"]):
                sc = EP.p2_score(None, T, k=k, ridx=ridx)
                s = EP.p2_summary(sc)
                s["per_tier"] = {tier: {arm: float(g["hit"].mean()) for arm, g in sub.groupby("arm")}
                                 for tier, sub in sc[sc["testable"]].groupby("tier")}
                s["n_per_tier"] = sc[sc["testable"] & (sc["arm"] == "S1-obs")].groupby("tier").size().to_dict()
                s["trend_test"] = ("not run: a tier has < trend_test_min_per_tier tests"
                                   if min(s["n_per_tier"].values(), default=0) < PR.dp()["clinical_s10"]["trend_test_min_per_tier"]
                                   else "eligible (not implemented)")
                headline = set_ == "non_proxy" and k == pc["k"] and subset == summ["headline_subset"]
                s["primary"] = bool(headline)
                summ["results"][f"{set_}|{subset}|k{k}"] = s
                scores_all.append(sc.assign(set=set_, subset=subset, k=k))
    out = P.RESULTS / "clinical"
    out.mkdir(parents=True, exist_ok=True)
    pd.concat(scores_all, ignore_index=True).to_csv(out / "p2_scores.csv", index=False)
    # Q4b actionability
    cg = cgt.assign(target_set=[";".join(ez.get(x, x) for x in EP._targets(v)) for v in cgt["target_set"]])
    s1h = IO.read_arms(arms=["S1-obs"], columns=["gene", "feature", "beta"], selected_only=True)
    act = EP.actionability(IO.catalogue(), s1h, cg)                    # A7: catalogue = MD-global
    act["list"] = act["list"].replace({"MD-strict": "MD-global"})
    act2 = EP.actionability(IO.md_strict(), s1h.iloc[:0], cg)          # MD-strict reported alongside
    act2 = act2[act2["list"] == "MD-strict"]
    pd.concat([act, act2], ignore_index=True).to_csv(out / "actionability.csv", index=False)
    # Q4c population on D
    T = EP.p2_test_set(cgt, U, "non_proxy", entrez_of=ez)
    T = T[T["testable"]]
    dev = list(p1.dev_ids())
    genes = sorted({g for ts in T["targets"] for g in ts if g in set(U["gene"])})
    feats = sorted(set(T["feature_D"]))
    Y = guard.load_effects(dev, scope="dev", genes=genes)
    fitter = IO.fitter_for(dev, genes, feats, Y)
    drops = {(str(g), f): (d if isinstance(d, str) else "") for g, f, d in zip(U["gene"], U["feature"], U["dropped_cov"])}
    model = p1.model_table(["OncotreeCode", "OncotreeSubtype", "OncotreePrimaryDisease", "OncotreeLineage"]).reindex(dev)
    terms = AN.oncotree_terms()
    unmapped = set()

    def ind_of(t):
        m, um = AN.indicated_lines(EP._targets(str(t.get("indication_codes", ""))), model, terms)
        unmapped.update(um)
        return m

    per, lab = EP.population_by_lineage(fitter, T, drops, model["OncotreeLineage"].fillna("").to_numpy(),
                                        model["OncotreeCode"].fillna("").to_numpy(), indicated_of=ind_of)
    rule = ELG.dependent_rule(Y)
    pc_ = p1.fold_centre(p1.oof_predictions("m5", columns=genes, ids=dev))
    acc = ELG.accurate_dependent_D(ELG.dependent(Y, rule), pc_)
    prec = []
    for _, t in T.iterrows():
        agn = str(t.get("tumour_agnostic", "0")) in ("1", "True", "true")
        ind = np.ones(len(dev), bool) if agn else ind_of(t)
        x = fitter.X[:, int(fitter.fpos[t["feature_D"]])]
        for g in t["targets"]:
            if g in acc.columns:
                prec.append({"test_id": t["test_id"], "gene": g, **AN.population_precision(x, acc[g].to_numpy(), ind)})
    lab = lab.merge(pd.DataFrame(prec), on=["test_id", "gene"], how="left") if len(lab) else pd.DataFrame(prec)
    per.to_csv(out / "q4c_population.csv", index=False)
    lab.to_csv(out / "q4c_labels.csv", index=False)
    summ["q4c"] = {"labels": lab["label"].value_counts().to_dict() if "label" in lab else {},
                   "unmapped_indication_codes": sorted(unmapped),
                   "indication_mapping": "OncotreeCode, or OncoTree term (frozen MOAlmanac code/term pairs) equal to the "
                                         "line's Model.csv Subtype / PrimaryDisease / Lineage name; no full OncoTree file"}
    IO.write_json(out / "q4_summary.json", summ)
    AN.write_timing("q4_summary", {"n_tests": int(len(T))}, wall_s=time.perf_counter() - t0, workers=1)
    AN.append_registry("scripts/run_q4.py", "ok", ["results/clinical/q4_summary.json", "results/clinical/p2_scores.csv",
                                                   "results/clinical/actionability.csv",
                                                   "results/clinical/q4c_population.csv",
                                                   "results/clinical/q4c_labels.csv"])
    print(summ["results"][f"non_proxy|{summ['headline_subset']}|k{pc['k']}"]["verdict"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
