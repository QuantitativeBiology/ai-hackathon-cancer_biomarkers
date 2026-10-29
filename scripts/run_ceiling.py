"""Data ceiling per CGT test (DESIGN v3.2 §10.4; F7; design_params clinical_s10) on D only.

    python -B scripts/run_ceiling.py          # -> results/clinical/cgt_ceiling.csv (+ ceiling_summary.json)

For every ceiling-testable test of configs/clinical_ground_truth.csv (frozen ``ceiling_testable``: an S/W test whose
feature is assayable, target in the library, >= testable_min_carriers D carriers and >= ceiling_min_carriers carriers
and non-carriers in its comparator), S1-obs on that test itself:
- lines = the test's comparator among D lines (build_cgt.comparator_mask, the frozen comparator rule) with a value of
  the test's feature (build_cgt.Molecular, the frozen CGT feature rule; the carrier counts are asserted equal to the
  frozen n_carriers_comparator_D / n_noncarriers_comparator_D);
- model = the §5.3 confounder model (endpoints.RFitter: design blocks, droppable covariates with the target's D drops
  from the S1-obs arm, own CN / EXPGATE; CR2 + the df rule), on D effects (guard scope dev);
- p one-sided in the expected direction (S: beta < 0, W: beta > 0), and in the opposite direction; target sets with
  |T| >= 2 use the family rule (min over members x |T|, capped at 1);
- robustness: one-sided Wilcoxon rank-sum of the effects, carriers vs non-carriers in the comparator (same family rule);
- BH over the ceiling-testable tests (F7), separately for the expected and the opposite direction; label
  ``data-supported`` (q < ceiling_bh_q, expected direction), ``data-contradicted`` (q_opposite < ceiling_bh_q),
  ``data-null`` otherwise (a test whose carrier indicator is collinear with the design blocks in its comparator,
R² > collinear_with_design_r2_gt, is flagged and never data-supported, as in S1); every other test is ``untestable`` (frozen untestable_reason, or "sign not S/W").
The labels are asserted to follow from the numbers. Timing stem cgt_ceiling; registry row. Refuses before freeze-bounds.
Comparators of type within-indication:stratified get no separate indication fixed effect: the §5.3 lineage block of
the model is used as built (stated in the summary).
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats as sps  # noqa: E402

from biomark import analyses as AN  # noqa: E402
from biomark import endpoint_io as IO  # noqa: E402
from biomark import endpoints as EP  # noqa: E402
from biomark import params as PR  # noqa: E402
from biomark import paths as P  # noqa: E402

OUT = P.RESULTS / "clinical" / "cgt_ceiling.csv"
SUMMARY = P.RESULTS / "clinical" / "ceiling_summary.json"
LABELS = ("data-supported", "data-contradicted", "data-null", "untestable")


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    n = int(ok.sum())
    if n:
        o = np.argsort(p[ok])
        r = p[ok][o] * n / np.arange(1, n + 1)
        r = np.minimum.accumulate(r[::-1])[::-1]
        out = np.empty(n)
        out[o] = np.minimum(r, 1.0)
        q[ok] = out
    return q


def family(ps: list, n_targets: int) -> float:
    """§10.2 family rule: min over members x |T| (capped at 1); single target: the member's p."""
    ps = [p for p in ps if np.isfinite(p)]
    if not ps:
        return np.nan
    return float(min(1.0, min(ps) * (n_targets if n_targets >= 2 else 1)))


def label_of(q: float, q_opp: float, qmax: float, collinear: bool = False) -> str:
    """S1 rule (stats: collinear_with_design is tested but never a hit): a collinear test is not data-supported."""
    if np.isfinite(q) and q < qmax and not collinear:
        return "data-supported"
    if np.isfinite(q_opp) and q_opp < qmax:
        return "data-contradicted"
    return "data-null"


def check_labels(t: pd.DataFrame, qmax: float) -> None:
    """The labels follow from the numbers (DESIGN §10.4 'recomputed from numbers')."""
    for _, r in t.iterrows():
        want = (label_of(r["q"], r["q_opposite"], qmax, bool(r.get("collinear_with_design", False)))
                if r["ceiling_testable"] else "untestable")
        if r["label"] != want:
            raise AssertionError(f"{r['test_id']}: label {r['label']} != {want} from the numbers")


def wilcoxon_one_sided(y_car, y_non, sign: float) -> float:
    if len(y_car) < 1 or len(y_non) < 1:
        return np.nan
    alt = "less" if sign < 0 else "greater"
    return float(sps.mannwhitneyu(y_car, y_non, alternative=alt).pvalue)


def _build_cgt():
    spec = importlib.util.spec_from_file_location("build_cgt", Path(__file__).resolve().parent / "build_cgt.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def drops_of(genes) -> dict:
    """(g, feature) -> D drops of the S1-obs row (eligible universe, then U_ctrl)."""
    out = {}
    for u in IO.ALL_UNIVERSES:
        d = IO.arms_dir(u)
        if not any(d.glob("part-*.parquet")):
            continue
        a = IO.read_arms(arms=["S1-obs"], columns=["gene", "feature", "dropped_cov"], genes=sorted(genes), dataset_dir=d)
        for g, f, x in zip(a["gene"].astype(str), a["feature"], a["dropped_cov"]):
            out.setdefault((g, f), x if isinstance(x, str) else "")
    return out


def run(check=None) -> dict:
    EP.require_descends("freeze-bounds", check=check)
    from biomark import genesets as GS
    from biomark import guard, p1
    t_start = AN.utc_now()
    t0 = time.perf_counter()
    s10 = PR.dp()["clinical_s10"]
    qmax = float(s10["ceiling_bh_q"])
    cgt = pd.read_csv(P.CONFIGS / "clinical_ground_truth.csv", keep_default_na=False)
    syms = sorted({x for v in cgt["target_set"] for x in EP._targets(v)})
    ez = GS.symbol_to_entrez(syms).dropna().astype(str).to_dict()
    bc = _build_cgt()
    import json
    rules = json.loads(bc.RULES.read_text(encoding="utf-8"))
    sp = p1.splits()
    mt = p1.model_table(["OncotreeLineage", "OncotreePrimaryDisease", "OncotreeCode", "ModelSubtypeFeatures",
                         "CellLineName", "StrippedCellLineName"])
    dev_ids = list(p1.dev_ids())
    pop = sp[sp["role"] == "dev"].copy().join(mt.drop(columns=["OncotreeLineage"]), how="left").reindex(dev_ids)
    mol = bc.Molecular(rules, pop, mt)
    ct = cgt[cgt["ceiling_testable"].astype(str) == "1"]
    tgenes = sorted({ez[x] for v in ct["target_set"] for x in EP._targets(v) if x in ez})
    Y = guard.load_effects(dev_ids, scope="dev", genes=tgenes)
    drops = drops_of(tgenes)
    cov_all = IO.line_covariates(dev_ids)
    oc_all, eg_all = IO.line_gene_covs(dev_ids, tgenes)
    rows = []
    for _, r in cgt.iterrows():
        rec = {"test_id": r["test_id"], "feature": r["feature"], "feature_D": r["feature_D"], "target_set": r["target_set"],
               "sign": r["sign"], "partition": r["partition"], "ko_class": r["ko_class"], "tier": r["tier"],
               "proxy": r["proxy"], "comparator": r["comparator"], "ceiling_testable": str(r["ceiling_testable"]) == "1"}
        if not rec["ceiling_testable"]:
            why = r["untestable_reason"] or ("sign not S/W" if r["sign"] not in ("S", "W") else "not ceiling-testable")
            rows.append(dict(rec, label="untestable", untestable_reason=why))
            continue
        sign = -1.0 if r["sign"] == "S" else 1.0
        x = mol.feature(r["feature"]).reindex(dev_ids)
        cm = bc.comparator_mask(r, pop, mol).reindex(dev_ids).fillna(False).astype(bool)
        keep = cm.to_numpy() & np.isfinite(x.to_numpy(dtype=float))
        n1, n0 = int((x[keep] == 1).sum()), int((x[keep] == 0).sum())
        if (n1, n0) != (int(r["n_carriers_comparator_D"]), int(r["n_noncarriers_comparator_D"])):
            raise AssertionError(f"{r['test_id']}: comparator counts {(n1, n0)} != frozen "
                                 f"{(r['n_carriers_comparator_D'], r['n_noncarriers_comparator_D'])}")
        ids = [i for i, k in zip(dev_ids, keep) if k]
        T = [ez[t] for t in EP._targets(r["target_set"]) if t in ez and ez[t] in Y.columns]
        X = x.reindex(ids).to_numpy(dtype=np.float64)[:, None]
        fit = EP.RFitter(X, pd.Series([0], index=[r["feature"]]), cov_all.loc[ids], oc_all.reindex(index=ids),
                         eg_all.reindex(index=ids), Y.reindex(index=ids), cov_all.loc[ids, "PatientID"].astype(str).to_numpy(),
                         binary=pd.Series([True], index=[r["feature"]]))
        from biomark import stats as ST
        rows_ok = np.ones(len(ids), dtype=bool)
        Qd = ST.orth(fit.cov.Z(rows_ok, fit.cov.design_names))
        r2d = float(ST.r2_on(Qd, X, has_intercept=True)[0])
        per = []
        for g in T:
            f = fit.fit(g, r["feature"], drops.get((g, r["feature_D"]), ""), sign=sign)
            y = Y.loc[ids, g].to_numpy(dtype=float)
            xv = X[:, 0]
            ok = np.isfinite(y)
            pw = wilcoxon_one_sided(y[ok & (xv == 1)], y[ok & (xv == 0)], sign)
            pw_opp = wilcoxon_one_sided(y[ok & (xv == 1)], y[ok & (xv == 0)], -sign)
            p_opp = 1.0 - f["p_one"] if np.isfinite(f["p_one"]) else np.nan
            per.append((g, f, pw, pw_opp, p_opp))
        best = min(per, key=lambda z: (np.inf if not np.isfinite(z[1]["p_one"]) else z[1]["p_one"])) if per else None
        nT = len(EP._targets(r["target_set"]))
        rec.update({"targets": ";".join(T), "gene_best": best[0] if best else "", "n_lines": len(ids),
                    "n_carriers": n1, "n_noncarriers": n0,
                    "beta": best[1]["beta"] if best else np.nan, "se": best[1]["se"] if best else np.nan,
                    "df": best[1]["df"] if best else np.nan, "df_rule": best[1].get("df_rule", "") if best else "",
                    "p_one_sided": family([z[1]["p_one"] for z in per], nT),
                    "p_opposite": family([z[4] for z in per], nT),
                    "wilcoxon_p_one_sided": family([z[2] for z in per], nT),
                    "wilcoxon_p_opposite": family([z[3] for z in per], nT),
                    "drops": drops.get((best[0], r["feature_D"]), "") if best else "",
                    "r2_design": r2d, "collinear_with_design": bool(r2d > PR.stats_params()["collinear_r2"]),
                    "untestable_reason": ""})
        rows.append(rec)
    t = pd.DataFrame(rows)
    m = t["ceiling_testable"]
    t["q"] = np.nan
    t["q_opposite"] = np.nan
    t["wilcoxon_q"] = np.nan
    t.loc[m, "q"] = bh(t.loc[m, "p_one_sided"].to_numpy())
    t.loc[m, "q_opposite"] = bh(t.loc[m, "p_opposite"].to_numpy())
    t.loc[m, "wilcoxon_q"] = bh(t.loc[m, "wilcoxon_p_one_sided"].to_numpy())
    t["collinear_with_design"] = t["collinear_with_design"].fillna(False).astype(bool)
    t.loc[m, "label"] = [label_of(a, b, qmax, c) for a, b, c in
                         zip(t.loc[m, "q"], t.loc[m, "q_opposite"], t.loc[m, "collinear_with_design"])]
    t["wilcoxon_agrees"] = np.where(m, (t["wilcoxon_q"] < qmax) == (t["label"] == "data-supported"), np.nan)
    check_labels(t, qmax)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    t.to_csv(OUT, index=False)
    summ = {"n_tests": int(len(t)), "n_ceiling_testable": int(m.sum()),
            "labels": t["label"].value_counts().to_dict(),
            "labels_ceiling_testable_by_ko_class": t[m].groupby(["ko_class", "label"]).size().unstack(fill_value=0).to_dict("index"),
            "wilcoxon_agrees": int(np.nansum(t.loc[m, "wilcoxon_agrees"].astype(float))),
            "collinear_with_design": sorted(t.loc[m & t["collinear_with_design"], "test_id"]),
            "rules": {"bh_q": qmax, "family": s10["family_rule"], "robustness": s10["ceiling_robustness"],
                      "stratified_comparators": "no separate indication fixed effect; the §5.3 lineage block is used",
                      "lines": "D lines only (guard scope dev)"},
            "wall_s": time.perf_counter() - t0, "commit": __import__("biomark.gitutil", fromlist=["head"]).head(P.PROJECT)}
    IO.write_json(SUMMARY, summ)
    AN.write_timing("cgt_ceiling", {"n_tests": int(m.sum())}, start_utc=t_start, workers=1)
    AN.append_registry("scripts/run_ceiling.py", "ok",
                       ["results/clinical/cgt_ceiling.csv", "results/clinical/ceiling_summary.json"])
    return summ


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    s = run()
    print(f"ceiling: {s['n_ceiling_testable']} testable of {s['n_tests']}; labels {s['labels']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
