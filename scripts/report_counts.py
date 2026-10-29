"""Derived counts and rates quoted in REPORT.md, computed from committed result files only (DESIGN §13.5).

    python -B scripts/report_counts.py                 # -> results/report_counts.json (+ a results/registry.csv row)
    python -B scripts/report_counts.py --no-register   # draft: write the file only

Every value is counted or derived (percentage, ratio, total) from tracked files under results/ or configs/; no effects,
no R data and no cache file are read. Each top-level key names what it counts; ``_sources`` lists the input files
and ``_rules`` how each block is derived. Values are stored unrounded: the report rounds them to its printed
precision and scripts/check_report_numbers.py checks that.
"""
from __future__ import annotations

import json
import statistics
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd  # noqa: E402

from biomark import paths as P  # noqa: E402

OUT = "results/report_counts.json"


def _csv(rel: str, **kw) -> pd.DataFrame:
    return pd.read_csv(P.PROJECT / rel, **kw)


def _json(rel: str):
    return json.loads((P.PROJECT / rel).read_text(encoding="utf-8"))


def _true(s: pd.Series) -> pd.Series:
    return s.astype(str).str.lower().isin(["true", "1"])


def _in_list(col: pd.Series, name: str) -> pd.Series:
    return col.fillna("").astype(str).str.split("|").map(lambda v: name in v)


def confounders() -> dict:
    d = _csv("results/confounders/delta_flags.csv", dtype={"gene": str})
    df, sc = _true(d["delta_flagged"]), _true(d["subtype_confounded"])
    out = {"n_pairs": len(d), "n_pairs_by_universe": d["universe"].value_counts().to_dict(),
           "n_delta_flagged": int(df.sum()), "n_subtype_confounded": int(sc.sum()), "n_flagged_either": int((df | sc).sum()),
           "by_list": {}}
    for name in ("MD-global", "MD-strict", "S1-obs BH", "control"):
        m = _in_list(d["source"], name)
        out["by_list"][name] = {"n_pairs": int(m.sum()), "n_delta_flagged": int((m & df).sum()),
                                "n_subtype_confounded": int((m & sc).sum())}
    bf = d.loc[df, "blocks_flagged"].fillna("").astype(str)
    out["n_delta_flagged_by_blocks"] = bf.value_counts().to_dict()
    return out


def components() -> dict:
    c = _csv("results/components/component_selection.csv", dtype={"gene": str})
    g = c[_true(c["in_md_global"])]
    acc = g["access_class"].value_counts()
    return {"md_global_n_pairs": len(g), "md_global_by_layer": g["layer"].value_counts().to_dict(),
            "md_global_by_access_class": acc.to_dict(),
            "md_global_expression_derived_pct": 100.0 * acc.get("expression-derived", 0) / len(g)}


def locus() -> dict:
    r = _csv("results/locus/resolution.csv", dtype={"gene": str})
    mg = r[_in_list(r["source"], "MD-global")]
    lv = mg["level"].value_counts()
    a = _csv("results/locus/amplicon.csv", dtype={"gene": str})
    ma = a[_in_list(a["source"], "MD-global")]
    return {"md_global_cn_sil_pairs": len(mg), "md_global_by_level": lv.to_dict(),
            "md_global_locus_level_pct": 100.0 * lv.get("locus", 0) / len(mg),
            "md_global_amplicon_pairs": len(ma), "md_global_amplicon_by_label": ma["label"].value_counts().to_dict()}


def shortlist() -> dict:
    s = _csv("results/candidates/shortlist.csv", dtype={"gene": str})
    return {"n_rows": len(s), "by_feature_family": s["feature"].str.split(":").str[0].value_counts().to_dict()}


def set_a() -> dict:
    a = _csv("results/replication/setA_lines.csv", dtype={"gene": str})
    f = a["frac_dependent"].dropna()
    return {"n_genes": len(a), "n_lines": sorted(a["n_setA"].dropna().astype(int).unique().tolist()),
            "median_frac_dependent": float(statistics.median(f)), "n_genes_frac_dependent_ge_0.5": int((f >= 0.5).sum())}


def known_recovery() -> dict:
    k = _csv("results/candidates/known_recovery.csv", dtype={"genes": str})
    out = {"n_rows": len(k), "by_kind": k["kind"].value_counts().to_dict(), "rows": {}}
    for nm, m in (("gate_controls", (k["kind"] == "control") & (k["scope"] == "gate")),
                  ("class_check_controls", (k["kind"] == "control") & (k["scope"] == "class_check")),
                  ("cgt_tests", k["kind"] == "cgt")):
        x = k[m]
        out["rows"][nm] = {"n": len(x), "recovered_MD-global": int(_true(x["recovered_MD-global"]).sum()),
                           "recovered_MD-strict": int(_true(x["recovered_MD-strict"]).sum())}
    return out


def combinatorial() -> dict:
    v = _json("results/combinatorial/validation.json")
    out = {}
    for name in ("MD-global", "MD-strict", "S1-obs"):
        x = v[name]
        out[name] = {"rate_pct": 100.0 * x["rate"], "rate_ci_pct": [100.0 * y for y in x["rate_ci"]],
                     "background_rate_pct": 100.0 * x["background_rate"],
                     "matched_expected_rate_pct": 100.0 * x["matched_expected_rate"],
                     "rate_over_matched_expected": x["rate"] / x["matched_expected_rate"]}
    return out


def md_lists() -> dict:
    out = {}
    for name, rel in (("MD-global", "results/md/md_global.csv"), ("MD-strict", "results/md/md_strict.csv")):
        m = _csv(rel, dtype={"gene": str})
        out[name] = {"n_pairs": len(m), "n_genes": int(m["gene"].nunique())}
    return out


def q4c() -> dict:
    lab = _csv("results/clinical/q4c_labels.csv")
    return {"n_tests": len(lab), "by_label": lab["label"].value_counts().to_dict()}


def null() -> dict:
    s = _json("results/null/summary.json")
    b_raised = int(_json("configs/design_params.json")["null_s7"]["B_raised"])
    return {"wall_draw0_s": s["wall_draw0_s"], "B_raised": b_raised,
            "projected_B_raised_min": s["wall_draw0_s"] * b_raised / 60.0}


def lineage() -> dict:
    s = _json("results/lineage/summary.json")["by_source"]
    out = {}
    for name in ("MD-global", "MD-strict"):
        n, lab = s[name]["n"], s[name]["labels"]
        out[name] = {"n": n, "pct_by_label": {k: 100.0 * v / n for k, v in lab.items()}}
    return out


def eligibility() -> dict:
    """From the frozen configs/eligible_genes.csv."""
    e = _csv("configs/eligible_genes.csv", dtype={"entrez": str})
    el = e[_true(e["eligible"])]
    return {"n_genes_evaluated": len(e), "n_eligible": len(el), "n_eligible_by_tier": el["tier"].value_counts().to_dict()}


def controls_spec() -> dict:
    """From the frozen configs/positive_controls.csv: rows specified for the gate (scope gate or untestable), gate rows
    and their families."""
    c = _csv("configs/positive_controls.csv", dtype=str)
    g = c[_true(c["gate_eligible"])]
    return {"n_specified_for_gate": int(c["scope"].isin(["gate", "untestable"]).sum()), "n_gate_rows": len(g),
            "n_gate_families": int(g["family"].nunique())}


def independent() -> dict:
    """Derived from results/gates/G_indep.json (task 3 beta tolerance)."""
    d = _json("results/gates/G_indep.json")
    b = d["task3_en_obs"]["beta"]
    tol = float(d["tolerances"]["en_beta_rel_tol"])
    return {"mean_beta_max_rel": b["mean_beta_max_rel"], "en_beta_rel_tol": tol,
            "mean_beta_excess_over_tol_pct": 100.0 * (b["mean_beta_max_rel"] / tol - 1.0)}


def bh_family() -> dict:
    s = _json("results/benchmark/summary.json")
    v = _json("results/gates/G_indep.json")
    n = v["task2_bh_and_gate"]["S1-obs"]["F1"]["n_tests"]
    return {"n_rows_arms": s["n_rows_arms"], "n_rows_arms_millions": s["n_rows_arms"] / 1e6,
            "f1_n_tests": n, "f1_n_tests_millions": n / 1e6}


def m5_audits() -> dict:
    s = _json("results/m5_audits_summary.json")
    return {"hypermut_tp53_n_failed": int(sum(s["hypermut_tp53"]["failed"].values()))}


def lineage_tf() -> dict:
    """Lineage-TF class-check control genes (configs/positive_controls.csv scope class_check) in gene_classes.csv."""
    c = _csv("configs/positive_controls.csv", dtype=str)
    g = _csv("results/limits/gene_classes.csv", dtype={"gene": str, "failure_class": str})
    x = g[g["gene"].isin(set(c.loc[c["scope"] == "class_check", "dependency_entrez"]))]
    hi = x[x["r_within"] >= 0.10]
    return {"n_genes": len(x), "n_r_within_ge_0.10": len(hi), "classes_r_within_ge_0.10": hi["failure_class"].value_counts().to_dict(),
            "n_class3": int((x["failure_class"] == "3").sum())}


def lineage_by_class() -> dict:
    """MD-global pairs by cross-lineage label (base label; the '(single-lineage carriers)' variant counts as
    lineage-unresolvable), all pairs and class-9 pairs (results/lineage/cross_lineage.csv x pair_classes.csv)."""
    pc = _csv("results/limits/pair_classes.csv", dtype={"gene": str})
    cl = _csv("results/lineage/cross_lineage.csv", dtype={"gene": str})
    g = pc[_true(pc["in_md_global"])].merge(cl[["gene", "feature", "label"]], on=["gene", "feature"], how="left")
    g["base"] = g["label"].astype(str).str.split(" (", regex=False).str[0]
    allc = g["base"].value_counts()
    c9 = g.loc[pd.to_numeric(g["pair_class"], errors="coerce") == 9, "base"].value_counts()
    return {"all": allc.to_dict(), "class_9": c9.to_dict(),
            "class_9_share_pct": {k: 100.0 * c9.get(k, 0) / v for k, v in allc.items()}}


def ceiling_class8() -> dict:
    pc = _csv("results/limits/pair_classes.csv", dtype={"gene": str})
    gc = _csv("results/limits/gene_classes.csv", dtype={"gene": str, "failure_class": str})
    r = pc[pc["source"] == "cgt_ceiling"]
    g8 = gc[gc["failure_class"] == "8"]
    return {"n_ceiling_class8_rows": len(r), "n_ceiling_class8_genes": int(r["gene"].nunique()),
            "n_gene_class8": len(g8), "gene_class8_symbols": sorted(g8["symbol"].astype(str))}


def miss_testable() -> dict:
    """Misses among testable tests per arm and category (results/clinical/miss_decomposition.csv)."""
    m = _csv("results/clinical/miss_decomposition.csv")
    m = m[_true(m["testable"])]
    out = {}
    for (st, sub, k), g in m.groupby(["set", "subset", "k"]):
        out[f"{st}|{sub}|k{int(k)}"] = {a: x["category"].value_counts().to_dict() | {"n_misses": len(x)}
                                         for a, x in g.groupby("arm")}
    return out


def nc_detail() -> dict:
    lam = _csv("results/benchmark/lambda.csv", dtype={"gene": str})
    s = _csv("results/controls/nc3_self_null.csv", dtype={"gene": str})
    t = s[_true(s["testable"])].assign(fam=lambda d: d["feature"].str.split(":").str[0])
    return {"eligible_en_pred_fidelity_median": float(lam.loc[lam["arm"] == "EN-pred", "r2_oof_beyond_design"].median()),
            "nc3_by_family": {f: {"n_testable": len(x), "dependent_frac": float(x["frac_dependent"].mean())}
                              for f, x in t.groupby("fam")}}


def subtype_examples() -> dict:
    """Accurate-dependent lines outside the subtype (n_accurate_dependent - n_members) for the §6.5 examples."""
    s = _csv("results/subtypes/subtypes.csv", dtype={"gene": str})
    out = {}
    ex = {tuple(x) for x in json.loads((P.CONFIGS / "report_examples.json").read_text(encoding="utf-8"))["examples"]}
    for _, r in s[s["lists"].fillna("").str.contains("MD-global")].iterrows():
        if (str(r["gene"]), r["feature"]) in ex:
            out[f"{r['symbol']}|{r['feature']}"] = {"n_accurate_dependent": int(r["n_accurate_dependent"]),
                                                    "n_members": int(r["n_members"]),
                                                    "n_accurate_non_members": int(r["n_accurate_dependent"] - r["n_members"])}
    return out


PRE_POWER_FIX_COMMIT = "5de3c19"     # last commit before the R-power fix (50d40aa); committed files only


def pre_power_fix() -> dict:
    """Class-9 / shortlist counts as committed before the R-power fix (git show <commit>:<file>)."""
    import subprocess
    def show(rel):
        r = subprocess.run(["git", "-C", str(P.PROJECT), "show", f"{PRE_POWER_FIX_COMMIT}:{rel}"], capture_output=True,
                           text=True, encoding="utf-8", check=True)
        return r.stdout
    lim = json.loads(show("results/limits/limits_summary.json"))
    cards = json.loads(show("results/candidates/cards_summary.json"))
    return {"commit": PRE_POWER_FIX_COMMIT, "pair_classes_md_global": lim["pair_classes_md_global"],
            "n_shortlist": cards["n_shortlist"]}


def controls() -> dict:
    c = _json("results/controls/pc_summary.json")
    fam = c["gate_eligible_by_family"]
    return {"n_gate_eligible": c["n_gate_eligible"], "n_gate_families": len(fam),
            "n_families_with_ge_2_rows": sum(1 for v in fam.values() if v >= 2)}


def p2() -> dict:
    """Recovered (hit) and testable CGT tests per set|subset|k and arm, counted from p2_scores.csv (unique test ids)."""
    s = _csv("results/clinical/p2_scores.csv", dtype={"gene": str})
    s = s[_true(s["testable"])]
    out = {}
    for (st, sub, k), g in s.groupby(["set", "subset", "k"]):
        key = f"{st}|{sub}|k{int(k)}"
        tests = g.drop_duplicates("test_id")
        out[key] = {"n_testable": int(tests["test_id"].nunique()),
                    "testable_by_tier": tests["tier"].value_counts().to_dict(),
                    "hits": {a: int(x.loc[_true(x["hit"]), "test_id"].nunique()) for a, x in g.groupby("arm")}}
    return out


def main(argv=None) -> int:
    register = "--no-register" not in (sys.argv[1:] if argv is None else argv)
    out = {"confounders": confounders(), "components": components(), "locus": locus(), "shortlist": shortlist(),
           "set_a": set_a(), "known_recovery": known_recovery(), "combinatorial": combinatorial(), "p2": p2(), "controls": controls(), "md_lists": md_lists(), "q4c": q4c(), "null": null(), "lineage": lineage(), "eligibility": eligibility(),
           "controls_spec": controls_spec(), "independent": independent(), "bh_family": bh_family(),
           "m5_audits": m5_audits(), "lineage_tf": lineage_tf(),
           "lineage_by_class": lineage_by_class(), "ceiling_class8": ceiling_class8(), "miss_testable": miss_testable(),
           "nc_detail": nc_detail(), "pre_power_fix": pre_power_fix(),
           "subtype_examples": subtype_examples()}
    out["_sources"] = ["results/confounders/delta_flags.csv", "results/components/component_selection.csv",
                       "results/locus/resolution.csv", "results/locus/amplicon.csv", "results/candidates/shortlist.csv",
                       "results/replication/setA_lines.csv", "results/candidates/known_recovery.csv",
                       "results/combinatorial/validation.json", "results/clinical/p2_scores.csv",
                       "results/controls/pc_summary.json", "results/md/md_global.csv", "results/md/md_strict.csv",
                       "results/clinical/q4c_labels.csv", "results/null/summary.json", "configs/design_params.json",
                       "results/lineage/summary.json", "configs/eligible_genes.csv", "configs/positive_controls.csv",
                       "results/gates/G_indep.json", "results/benchmark/summary.json", "results/m5_audits_summary.json",
                       "results/limits/gene_classes.csv", "results/limits/pair_classes.csv",
                       "results/lineage/cross_lineage.csv", "results/clinical/miss_decomposition.csv",
                       "results/benchmark/lambda.csv", "results/controls/nc3_self_null.csv",
                       "results/subtypes/subtypes.csv", "git show 5de3c19:results/limits/limits_summary.json",
                       "git show 5de3c19:results/candidates/cards_summary.json"]
    out["_rules"] = {"by_list": "a row counts for every list named in its '|'-joined source column",
                     "pct": "100 x the fraction in the source file", "md_global": "rows whose source / flag marks MD-global"}
    f = P.PROJECT / OUT
    f.write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    if register:
        from biomark import analyses as AN
        AN.append_registry("scripts/report_counts.py", "ok", [OUT])
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
