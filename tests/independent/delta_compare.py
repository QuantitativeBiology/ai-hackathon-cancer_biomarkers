"""Task 4 follow-up: my confounding deltas (out/task4_delta.csv, seeds[2] pairs) vs the engine's
results/confounders/delta.csv and delta_flags.csv (sha256 checked against results/bulk_manifest.csv).

Tolerances (design_params.indep_s13 has no delta-specific key): beta_full, beta_minus, delta beta and SE_minus
relative <= s1_p_tol_abs_log10 (1e-6; the S1 numeric tolerance, since the delta refits are S1 fits), also reported
against en_beta_rel_tol (1e-3); flags per block and per pair exact (indep_s13.counts = exact).

    python -X utf8 -m tests.independent.delta_compare
"""
from __future__ import annotations

import hashlib
import json

import numpy as np
import pandas as pd

from . import data as DA

MAP = {"lineage": "lineage", "library": "library", "culture": "culture", "growth": "growth", "burden": "burden",
       "TP53": "TP53", "own CN": "own_CN"}


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    tol = DA.design_params()["indep_s13"]
    rtol, ltol = float(tol["s1_p_tol_abs_log10"]), float(tol["en_beta_rel_tol"])
    cdir = DA.ROOT / "results" / "confounders"
    man = pd.read_csv(DA.ROOT / "results" / "bulk_manifest.csv", header=None)
    shas = {}
    for name in ("delta.csv", "delta_flags.csv"):
        row = man[man[0].astype(str).str.endswith("confounders/" + name)]
        shas[name] = {"manifest": row[1].iloc[0] if len(row) else None, "actual": sha256(cdir / name)}
        shas[name]["match"] = shas[name]["manifest"] == shas[name]["actual"]
    mine = pd.read_csv(DA.OUT / "task4_delta.csv", dtype={"gene": str})
    eng = pd.read_csv(cdir / "delta.csv", dtype={"gene": str})
    fl = pd.read_csv(cdir / "delta_flags.csv", dtype={"gene": str})
    rows, pair_rows = [], []
    for _, m in mine.iterrows():
        e = eng[(eng["gene"] == m["gene"]) & (eng["feature"] == m["feature"])]
        for blk_e in e["block"]:
            pass
        for blk, key in list(MAP.items()) + [("primary_disease", "primary_disease")]:
            if blk == "primary_disease":
                er = e[e["block"].str.startswith("primary_disease")]
            else:
                er = e[e["block"] == blk]
            rec = {"gene": m["gene"], "feature": m["feature"], "block": blk, "n_engine_rows": int(len(er)),
                   "beta_full_mine": m["beta_full"], "beta_minus_mine": m[f"beta_minus_{key}"],
                   "se_minus_mine": m[f"se_minus_{key}"], "flag_mine": bool(m[f"flag_{key}"])}
            if len(er) == 1:
                r = er.iloc[0]
                rec.update({"beta_full_eng": r["beta_full"], "beta_minus_eng": r["beta_minus"],
                            "se_minus_eng": r["se_minus"], "flag_eng": bool(r["flagged"])})
                dm = rec["beta_full_mine"] - rec["beta_minus_mine"]
                de = r["beta_full"] - r["beta_minus"]
                rec["dbeta_mine"], rec["dbeta_eng"] = dm, de
                rel = lambda a, b: abs(a - b) / max(abs(b), 1e-300)
                rec["rel_beta_full"] = rel(rec["beta_full_mine"], r["beta_full"])
                rec["rel_beta_minus"] = rel(rec["beta_minus_mine"], r["beta_minus"])
                rec["rel_dbeta"] = rel(dm, de)
                rec["rel_se_minus"] = rel(rec["se_minus_mine"], r["se_minus"])
                rec["max_rel"] = max(rec["rel_beta_full"], rec["rel_beta_minus"], rec["rel_dbeta"], rec["rel_se_minus"])
            rows.append(rec)
        f = fl[(fl["gene"] == m["gene"]) & (fl["feature"] == m["feature"])]
        mb = [b for b in MAP if m[f"flag_{MAP[b]}"]]
        pr = {"gene": m["gene"], "feature": m["feature"], "blocks_flagged_mine": ";".join(mb),
              "delta_flagged_mine": bool(mb), "subtype_confounded_mine": bool(m["flag_primary_disease"]),
              "n_engine_rows": int(len(f))}
        if len(f) == 1:
            r = f.iloc[0]
            bf = r["blocks_flagged"]
            eb = sorted(x for x in (bf.replace("|", ";").split(";") if isinstance(bf, str) else []) if x)
            pr.update({"blocks_flagged_eng": ";".join(eb), "delta_flagged_eng": bool(r["delta_flagged"]),
                       "subtype_confounded_eng": bool(r["subtype_confounded"])})
            pr["agree"] = (sorted(mb) == eb and pr["delta_flagged_mine"] == pr["delta_flagged_eng"]
                           and pr["subtype_confounded_mine"] == pr["subtype_confounded_eng"])
        pair_rows.append(pr)
    D = pd.DataFrame(rows)
    P = pd.DataFrame(pair_rows)
    D.to_csv(DA.OUT / "task4_engine_comparison.csv", index=False)
    have = D["n_engine_rows"] == 1
    Dh = D[have]
    worst = Dh.sort_values("max_rel", ascending=False).head(5)
    flag_dis = Dh[Dh["flag_mine"] != Dh["flag_eng"]]
    res = {"tolerances": {"numeric_rel": rtol, "numeric_rel_loose": ltol, "flags": "exact",
                          "note": "indep_s13 has no delta key; S1 tolerance used for the S1 refits"},
           "sha256": shas, "n_pairs": int(len(P)), "n_block_rows": int(len(D)),
           "n_block_rows_found_in_engine": int(have.sum()),
           "n_pairs_found_in_flags": int((P["n_engine_rows"] == 1).sum()),
           "max_rel_beta_full": float(Dh["rel_beta_full"].max()), "max_rel_beta_minus": float(Dh["rel_beta_minus"].max()),
           "max_rel_dbeta": float(Dh["rel_dbeta"].max()), "max_rel_se_minus": float(Dh["rel_se_minus"].max()),
           "n_rows_within_rel": int((Dh["max_rel"] <= rtol).sum()),
           "n_rows_within_rel_loose": int((Dh["max_rel"] <= ltol).sum()),
           "n_block_flag_disagreements": int(len(flag_dis)),
           "n_pair_flag_disagreements": int((~P.get("agree", pd.Series(False, index=P.index)).fillna(False).astype(bool)).sum()),
           "worst_rows": worst.to_dict("records"),
           "flag_disagreements": flag_dis.to_dict("records"),
           "pairs": P.to_dict("records")}
    # diagnostic (not my reading): library block = library_group only, screen days kept
    from .confounder import s1_gene
    E = DA.effects(sorted(set(mine["gene"])))
    el = eng[eng["block"] == "library"]
    diag = []
    for _, m in mine.iterrows():
        x = s1_gene(m["gene"], E[m["gene"]].to_numpy(float), [m["feature"]], drop_extra={m["feature"]: {"library"}})[0]
        e = el[(el["gene"] == m["gene"]) & (el["feature"] == m["feature"])].iloc[0]
        diag.append(max(abs(x["beta"] - e["beta_minus"]) / abs(e["beta_minus"]), abs(x["se"] - e["se_minus"]) / e["se_minus"]))
    nonlib = Dh[Dh["block"] != "library"]
    res["by_block_max_rel"] = Dh.groupby("block")["max_rel"].max().to_dict()
    res["non_library_rows"] = {"n": int(len(nonlib)), "n_within_rel": int((nonlib["max_rel"] <= rtol).sum()),
                               "max_rel": float(nonlib["max_rel"].max())}
    res["diagnostic_library_group_only"] = {
        "max_rel_beta_or_se_minus_vs_engine": float(max(diag)),
        "reading": "engine's library delta removes library_group only (screen days stay, as a droppable covariate); "
                   "mine removes library_group + screen days because DESIGN 5.3 defines the library block as "
                   "'library group (...) + screen days within Humagne'. Not adopted (pre-declared reading R0 #8)."}
    res["sha_note"] = "delta.csv matches results/bulk_manifest.csv; delta_flags.csv has no manifest row"
    res["pass_flags"] = bool(res["n_block_flag_disagreements"] == 0 and res["n_pair_flag_disagreements"] == 0)
    res["pass_numeric"] = bool(res["n_rows_within_rel"] == len(D))
    res["pass"] = bool(shas["delta.csv"]["match"] and have.all()
                       and res["n_rows_within_rel"] == len(D) and res["n_block_flag_disagreements"] == 0
                       and res["n_pair_flag_disagreements"] == 0)
    p = DA.ROOT / "results" / "gates" / "G_indep.json"
    g = json.load(open(p, encoding="utf-8"))
    g["task4_engine_comparison"] = json.loads(json.dumps(res, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    g.setdefault("overall", {})["task4_engine_comparison"] = res["pass"]
    with open(p, "w", encoding="utf-8") as fh:
        json.dump(g, fh, indent=1)
    print(json.dumps({k: v for k, v in res.items() if k not in ("pairs", "worst_rows", "flag_disagreements")}, indent=1))
    print(pd.DataFrame(res["worst_rows"])[["gene", "feature", "block", "beta_minus_mine", "beta_minus_eng", "se_minus_mine",
                                           "se_minus_eng", "max_rel"]].to_string())
    print(P.to_string())


if __name__ == "__main__":
    main()
