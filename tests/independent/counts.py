"""Task 2: BH counts from the engine's p values with an independently written family definition, and the gate
recall from the frozen inputs of results/gates/G_pc_attempt1.json.

Family F1 (design 5.6, test_families_s5.F1_S1_main): every (g, d) in U (g eligible; d in D_g minus self-null) whose
d is not g's own EXPDOSE, with a finite p (tests the engine dropped for < 5 carriers or zero residual variance have no
p and are not tests). F2 = the own-EXPDOSE tests. BH (step-up) at q = s1_s5.bh_q over the whole family.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from . import data as DA


def bh_count(p: np.ndarray, q: float) -> tuple[int, float]:
    p = np.asarray(p, float)
    p = p[np.isfinite(p)]
    m = len(p)
    if m == 0:
        return 0, np.nan
    s = np.sort(p)
    ok = np.where(s <= q * np.arange(1, m + 1) / m)[0]
    if not len(ok):
        return 0, np.nan
    k = int(ok[-1] + 1)
    return k, float(s[k - 1])


def bh_reject(p: np.ndarray, q: float) -> np.ndarray:
    k, cut = bh_count(p, q)
    p = np.asarray(p, float)
    return np.isfinite(p) & (p <= cut) if k else np.zeros(len(p), bool)


def read_arm(arm: str, columns):
    from biomark import endpoint_io as IO
    return IO.read_arms(arms=[arm], columns=columns)


def universe_check(s1: pd.DataFrame) -> dict:
    """Engine S1 rows vs my D_g for every eligible gene."""
    elig = DA.eligible()
    eng = s1.groupby("gene")["feature"].apply(set).to_dict()
    n_mine = n_eng = n_only_mine = n_only_eng = 0
    genes_diff = []
    ex = []
    for g in elig:
        u = set(DA.universe(g))
        e = eng.get(g, set())
        n_mine += len(u); n_eng += len(e)
        a, b = u - e, e - u
        n_only_mine += len(a); n_only_eng += len(b)
        if a or b:
            genes_diff.append(g)
            if len(ex) < 10:
                ex.append({"gene": g, "only_mine": sorted(a)[:5], "only_engine": sorted(b)[:5]})
    return {"n_pairs_mine": n_mine, "n_pairs_engine": n_eng, "only_mine": n_only_mine, "only_engine": n_only_eng,
            "n_genes_differing": len(genes_diff), "genes_differing": genes_diff, "examples": ex,
            "engine_genes_not_eligible": sorted(set(eng) - set(elig))[:10]}


def family_bh(s1: pd.DataFrame, q: float) -> dict:
    """BH over the whole family (collinear_with_design tests included: they are tests of F1); list membership
    then excludes collinear_with_design pairs (5.3: 'never counted as determinants') - fix 1 in README."""
    own = s1["feature"] == ("EXPDOSE:" + s1["gene"].astype(str))
    fin = np.isfinite(s1["p"].to_numpy(float))
    col = s1["collinear_with_design"].fillna(False).astype(bool).to_numpy()
    out = {"n_rows_no_p": int((~fin).sum())}
    idx = {}
    for name, m in (("F1", ~own & fin), ("F2", own & fin)):
        F = s1[m]
        k, cut = bh_count(F["p"].to_numpy(), q)
        rej = (F["p"].to_numpy() <= cut) if k else np.zeros(len(F), bool)
        lst = rej & ~col[m.to_numpy()]
        out[name] = {"n_tests": int(len(F)), "bh_rejections": k, "p_cut": cut,
                     "collinear_among_rejections": int((rej & col[m.to_numpy()]).sum()),
                     "list_hits": int(lst.sum())}
        idx[name] = F.index[lst]
    return out, idx


def run(q: float | None = None) -> dict:
    q = float(DA.design_params()["s1_s5"]["bh_q"]) if q is None else q
    res = {"q": q}
    cols = ["gene", "feature", "p", "beta", "test_family", "selected", "selected_F2", "q", "collinear_with_design"]
    hits_by_arm = {}
    for arm in ("S1-obs", "S1-pred"):
        s1 = read_arm(arm, cols)
        s1["gene"] = s1["gene"].astype(str)
        fam, hidx = family_bh(s1, q)
        hit_idx = hidx["F1"]
        eng_F1 = int((s1["selected"].fillna(False).astype(bool) & (s1["test_family"] == "F1")).sum())
        eng_F2 = int(s1["selected_F2"].fillna(False).astype(bool).sum()) if "selected_F2" in s1 else None
        own = s1["feature"] == ("EXPDOSE:" + s1["gene"])
        fam["engine"] = {"F1_selected": eng_F1, "F2_selected": eng_F2,
                         "F1_rows_engine_label": int((s1["test_family"] == "F1").sum()),
                         "F1_rows_engine_label_with_p": int(((s1["test_family"] == "F1") & s1["p"].notna()).sum()),
                         "label_vs_mine_mismatch": int(((s1["test_family"] == "F1") != ~own).sum())}
        mine = pd.Series(False, index=s1.index)
        mine.loc[hit_idx] = True
        eng = s1["selected"].fillna(False).astype(bool) & (s1["test_family"] == "F1")
        fam["F1_hit_set_symmetric_difference"] = int((mine != eng).sum())
        mine2 = pd.Series(False, index=s1.index)
        mine2.loc[hidx["F2"]] = True
        fam["F2_hit_set_symmetric_difference"] = int((mine2 != s1["selected_F2"].fillna(False).astype(bool)).sum())
        fam["F1_match_exact"] = fam["F1"]["list_hits"] == eng_F1 and fam["F1_hit_set_symmetric_difference"] == 0
        fam["F2_match_exact"] = fam["F2"]["list_hits"] == eng_F2 and fam["F2_hit_set_symmetric_difference"] == 0
        if arm == "S1-obs":
            fam["universe"] = universe_check(s1)
            hits = s1.loc[hit_idx, ["gene", "feature", "beta"]].copy()
            hits_by_arm[arm] = hits
            res["_s1obs_hit_genes"] = sorted(set(hits["gene"]))
        res[arm] = fam
        del s1
    res["gate"], res["md"] = gate(hits_by_arm["S1-obs"], q)
    return res


def _cp(k, n):
    lo = 0.0 if k == 0 else stats.beta.ppf(0.025, k, n - k + 1)
    hi = 1.0 if k == n else stats.beta.ppf(0.975, k + 1, n - k)
    return [float(lo), float(hi)]


def gate(s1_hits: pd.DataFrame, q: float):
    """Gate recall (9.3) on MD-strict from the attempt-1 inputs, plus my recomputation of the MD-strict and
    MD-global lists from md_strict_all.csv (EN-pred selected + D-stage tests) and my S1-obs BH hits."""
    gp = DA.design_params()["gate_s9"]
    elig = pd.read_csv(DA.ROOT / "results" / "controls" / "pc_eligibility.csv")
    pc = pd.read_csv(DA.ROOT / "configs" / "positive_controls.csv", dtype={"dependency_entrez": str})
    den = elig[(elig["gate_eligible"].astype(bool)) & (elig["scope"] == "gate")]["control_id"]
    pc = pc.set_index("control_id").loc[den]
    # MD-strict: engine catalogue file (the G_pc input) and my recomputation from md_strict_all.csv
    ms = pd.read_csv(DA.ROOT / "results" / "md" / "md_strict.csv", dtype={"gene": str})
    ma = pd.read_csv(DA.ROOT / "results" / "md" / "md_strict_all.csv", dtype={"gene": str})
    md_check = {}
    testable = ma["d_stage_testable"].fillna(False).astype(bool) & ma["p_sc"].notna()
    rej = np.zeros(len(ma), bool)
    rej[testable.to_numpy()] = bh_reject(ma.loc[testable, "p_sc"].to_numpy(), float(DA.design_params()["md_strict_s5"]["q_sc_bh"]))
    same = np.sign(ma["beta"].to_numpy()) == np.sign(ma["maj_sign"].to_numpy())
    col = ma["collinear_with_design"].fillna(False).astype(bool).to_numpy()
    my_ms = ma[rej & same & ~col]
    md_check["MD-strict"] = {"n_EN_pred_selected": int(len(ma)), "n_testable": int(testable.sum()),
                             "bh_rejections_q_sc": int(rej.sum()),
                             "p_cut_sc": float(ma.loc[rej, "p_sc"].max()) if rej.any() else None,
                             "mine_same_sign_noncollinear": int(len(my_ms)),
                             "engine": int(len(ms)),
                             "set_equal": set(zip(my_ms.gene, my_ms.feature)) == set(zip(ms.gene, ms.feature))}
    hs = s1_hits.assign(gene=s1_hits["gene"].astype(str))
    hk = {(g, f): b for g, f, b in zip(hs.gene, hs.feature, hs.beta)}
    mg = ma[[((g, f) in hk) and np.sign(hk[(g, f)]) == np.sign(s) for g, f, s in zip(ma.gene, ma.feature, ma.maj_sign)]]
    eng_mg = pd.read_csv(DA.ROOT / "results" / "md" / "md_global.csv", dtype={"gene": str})
    md_check["MD-global"] = {"mine": int(len(mg)), "engine": int(len(eng_mg)),
                             "set_equal": set(zip(mg.gene, mg.feature)) == set(zip(eng_mg.gene, eng_mg.feature))}

    def recovered(listdf, beta_col):
        key = {(g, f): b for g, f, b in zip(listdf["gene"].astype(str), listdf["feature"], listdf[beta_col])}
        rows = []
        for cid, r in pc.iterrows():
            feats = []
            for col in ("feature_primary_ids", "feature_secondary_ids"):
                v = r[col]
                if isinstance(v, str) and v.strip():
                    feats += [x for x in v.split(";") if x]
            want = -1 if r["expected_sign"] == "S" else 1
            g = str(r["dependency_entrez"])
            ok = any((g, f) in key and np.sign(key[(g, f)]) == want for f in feats)
            rows.append({"control_id": cid, "family": r["family"], "recovered": ok})
        return pd.DataFrame(rows)

    def summarise(rec):
        k, n = int(rec["recovered"].sum()), len(rec)
        fam = {}
        for f, s in rec.groupby("family"):
            fam[f] = {"k": int(s["recovered"].sum()), "n": int(len(s))}
        fam_pass = {f: v["k"] / v["n"] >= gp["family_recall_min"] for f, v in fam.items()
                    if v["n"] >= gp["family_min_eligible_rows"]}
        return {"k": k, "n": n, "recall": k / n, "ci95": _cp(k, n), "by_family": fam, "family_pass": fam_pass,
                "pass": bool(k / n >= gp["overall_recall_min"] and all(fam_pass.values()))}

    import json
    att = json.load(open(DA.ROOT / "results" / "gates" / "G_pc_attempt1.json", encoding="utf-8"))
    out = {"MD-strict_engine_file": summarise(recovered(ms, "beta")),
           "MD-strict_my_recompute": summarise(recovered(my_ms, "beta")),
           "S1-obs_my_bh": summarise(recovered(hs, "beta")),
           "MD-global_my_recompute": summarise(recovered(mg, "beta"))}
    a = att["summary"]
    out["attempt1"] = {"k": a["overall"]["k"], "n": a["overall"]["n"],
                       "by_family": {f: {"k": v["k"], "n": v["n"]} for f, v in a["by_family"].items()},
                       "S1-obs": {"k": a["by_list"]["S1-obs"]["k"]}, "MD-global": {"k": a["by_list"]["MD-global"]["k"]}}
    per = {r["control_id"]: r["recovered_MD-strict"] for r in att["rows"]}
    mine = recovered(ms, "beta").set_index("control_id")["recovered"].to_dict()
    out["row_disagreements_vs_attempt1"] = sorted(c for c in per if per[c] != mine.get(c))
    e, m = out["attempt1"], out["MD-strict_engine_file"]
    out["exact_match"] = (m["k"] == e["k"] and m["n"] == e["n"] and not out["row_disagreements_vs_attempt1"]
                          and all(m["by_family"][f] == e["by_family"][f] for f in e["by_family"]))
    out["exact_match_my_recompute"] = out["MD-strict_my_recompute"]["k"] == e["k"]
    return out, md_check
