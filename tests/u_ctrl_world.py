"""Toy world with a non-eligible positive control (DESIGN v3.2 §2 U_ctrl), in the real on-disk layout under ``root``.

A planted, gate-eligible ADD-HOT control gene of biomark.toysim is declared non-eligible: the eligible arms run never
sees it, and the separate U_ctrl run (that gene + the toy's non-eligible LINEAGE-CHECK genes; BH within U_ctrl)
recovers it. Written: configs/{eligible_genes.csv, u_ctrl_genes.txt, positive_controls.csv},
results/controls/pc_eligibility.csv, cache/benchmark/arms{,_ctrl}/part-00000.parquet,
results/md{,_ctrl}/{md_strict,md_global,md_strict_all,md_global_all}.csv, results/benchmark/lambda.csv,
results/benchmark_ctrl/summary.json and cache/toy/meta.json (biomark.m5 toy workspace). No real data."""
import hashlib
import json
from pathlib import Path

import pandas as pd

from biomark import arms as A
from biomark import pipeline as PL
from biomark import toysim as TS

N_ELIGIBLE_CONTROLS = 6          # eligible-universe controls kept (speed)
MD_FILES = ("md_strict", "md_global", "md_strict_all", "md_global_all")


def _sha(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def build(root: Path) -> dict:
    toy = TS.make_toy()
    el = PL.toy_eligibility(toy)
    c = toy.controls.assign(dependency_entrez=toy.controls["dependency_entrez"].astype(str))
    gate = c[c["gate_eligible"].astype(bool) & c["dependency_entrez"].map(el["eligible"]).astype(bool)]
    moved = gate[gate["family"] == "ADD-HOT"].iloc[0]                    # declared non-eligible -> U_ctrl
    kept = gate[gate["dependency_entrez"] != moved["dependency_entrez"]].iloc[:N_ELIGIBLE_CONTROLS]
    lineage = sorted(g for g in c["dependency_entrez"] if not el.loc[g, "eligible"])
    elig_genes = sorted(set(kept["dependency_entrez"]))
    ctrl_genes = sorted({moved["dependency_entrez"], *lineage})
    cfg, res, cache = root / "configs", root / "results", root / "cache"
    for d in (cfg, res / "controls", res / "benchmark", res / "benchmark_ctrl", cache / "toy"):
        d.mkdir(parents=True, exist_ok=True)
    all_genes = list(map(str, toy.genes.index))
    eligible = set(el.index[el["eligible"]].astype(str)) - {moved["dependency_entrez"]}
    pd.DataFrame({"entrez": all_genes, "eligible": [g in eligible for g in all_genes]}).to_csv(
        cfg / "eligible_genes.csv", index=False)
    (cfg / "u_ctrl_genes.txt").write_text("".join(f"{g}\n" for g in ctrl_genes), encoding="utf-8")
    ctl = pd.concat([pd.DataFrame([moved]), kept])
    pd.DataFrame({"control_id": ctl["control_id"], "family": ctl["family"], "dependency_entrez": ctl["dependency_entrez"],
                  "feature_primary_ids": ctl["feature_primary"], "feature_secondary_ids": ctl["feature_secondary"],
                  "fallback_feature_ids": "", "gate_feature": "primary", "expected_sign": ctl["expected_sign"]}).to_csv(
        cfg / "positive_controls.csv", index=False)
    pd.DataFrame({"control_id": ctl["control_id"], "family": ctl["family"], "gate_eligible": "True",
                  "access_class": "direct", "en_max_attainable_freq": "1.0"}).to_csv(
        res / "controls" / "pc_eligibility.csv", index=False)
    runs = {}
    for u, genes in (("eligible", elig_genes), ("ctrl", ctrl_genes)):
        out = A.run_arms(PL.Data.from_toy(toy, genes=genes), arms=A.ARMS, workers=1)
        sfx = "" if u == "eligible" else "_ctrl"
        ds, md = cache / "benchmark" / f"arms{sfx}", res / f"md{sfx}"
        ds.mkdir(parents=True, exist_ok=True)
        md.mkdir(parents=True, exist_ok=True)
        a = out["arms"].copy()
        a["gene"] = a["gene"].astype(str)
        a.to_parquet(ds / "part-00000.parquet")                            # layout of endpoint_io.read_arms
        for k in MD_FILES:
            out[k].to_csv(md / f"{k}.csv", index=False)
        runs[u] = out
    runs["eligible"]["lam"].to_csv(res / "benchmark" / "lambda.csv", index=False)
    (res / "benchmark_ctrl" / "summary.json").write_text(
        json.dumps({"universe": "ctrl", "u_ctrl_genes_sha256": _sha(cfg / "u_ctrl_genes.txt")}), encoding="utf-8")
    (cache / "toy" / "meta.json").write_text(json.dumps({"seed": int(toy.seed), "predictor": toy.predictor,
                                                         "scenario": toy.scenario}), encoding="utf-8")
    m = runs["ctrl"]["md_strict"]
    assert ((m["gene"].astype(str) == moved["dependency_entrez"]) & (m["feature"] == moved["feature_primary"])
            & (m["beta"] < 0)).any(), "precondition: the moved control is recovered on MD-strict within U_ctrl"
    assert not (runs["eligible"]["arms"]["gene"].astype(str) == moved["dependency_entrez"]).any()
    return {"root": root, "toy": toy, "moved": moved, "kept": kept, "runs": runs, "ctrl_genes": ctrl_genes}
