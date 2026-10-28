"""Build the interpretable feature dictionary D (DESIGN v3.1 section 3).

Inputs (all must exist; run in this order): results/data_status.json (scripts/data_status.py; its sha256 is
stored here = data deadline), cache/covariates.parquet, configs/eligible_genes.csv, configs/gamma_genes.csv,
configs/fusion_drivers.csv, configs/partner_priors.csv (scripts/build_gene_lists.py) and, if present,
configs/cgt_features.csv (ALLELE / CLIN, CGT track).

Outputs: cache/dictionary/<FAMILY>.parquet (population rows x kept features, float32), sil_mechanism.parquet,
expgate.parquet, own_cn_log2.parquet, nmf_components.parquet; configs/dictionary_v1.json;
configs/dictionary_gene_features.csv; results/dictionary/{summary.csv, family_counts.csv, locus_groups.csv,
prog_annotation.csv, checks.json}. No gene-effect value is read.
Usage: python scripts/build_dictionary.py
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import re
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from biomark import dictionary as Dct  # noqa: E402
from biomark import genesets, gitutil, p1  # noqa: E402
from biomark import paths as P  # noqa: E402


def sha_file(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def known_positive_check(fus: pd.DataFrame) -> list[dict]:
    """configs/fusion_known_positives.csv lines (names resolved via Model.csv) vs FUS features (any FUS feature of
    the line's fusion partners)."""
    kp = pd.read_csv(P.PROJECT / P.load_config("build_params")["controls"]["fusion_known_positives"], dtype=str)
    mt = p1.model_table(["StrippedCellLineName"])["StrippedCellLineName"].dropna()
    mt = mt[mt.index.isin(fus.index)]          # population lines only (derived models share names)
    norm = lambda s: re.sub(r"[^A-Z0-9]", "", str(s).upper())  # noqa: E731
    by = {}
    for k, v in sorted(mt.items()):
        by.setdefault(norm(v), k)
    out = []
    for _, r in kp.iterrows():
        mid = by.get(norm(r["cell_line_name"]))
        genes = [g for g in re.split(r"::", r["fusion"])]
        ez = genesets.symbol_to_entrez(genes)
        feats = [f"FUS:{e}" for e in ez.dropna() if f"FUS:{e}" in fus.columns]
        car = [f for f in feats if mid in fus.index and fus.at[mid, f] == 1]
        out.append({"cell_line_name": r["cell_line_name"], "ModelID": mid, "fusion": r["fusion"],
                    "fus_features_available": feats, "carrier_of": car,
                    "pass": None if mid is None else bool(car),
                    "status": "line not in the Part 1 population" if mid is None else
                    ("carrier" if car else ("no FUS driver for this fusion" if not feats else "not a carrier"))})
    return out


def main() -> int:
    t0 = dt.datetime.now(dt.timezone.utc)
    cfg = P.load_config("build_params")["dictionary"]
    ds = P.RESULTS / "data_status.json"
    if not ds.exists():
        raise SystemExit("results/data_status.json missing: run scripts/data_status.py first (data deadline)")
    data_status_sha = sha_file(ds)
    prev = P.CONFIGS / "dictionary_v1.json"
    build_count = 1
    if prev.exists():
        try:
            build_count = int(json.loads(prev.read_text(encoding="utf-8")).get("build_count", 3)) + 1
        except (ValueError, json.JSONDecodeError):
            pass
    el = pd.read_csv(P.CONFIGS / "eligible_genes.csv", dtype={"entrez": str})
    elig = list(el.loc[el["eligible"], "entrez"])
    ctrl = list(genesets.control_dependency_genes()["entrez"].dropna())
    targets = list(dict.fromkeys(elig + ctrl))
    gamma = list(pd.read_csv(P.CONFIGS / "gamma_genes.csv", dtype={"entrez": str})["entrez"])
    drivers = pd.read_csv(P.CONFIGS / "fusion_drivers.csv", dtype={"entrez": str})
    partners = pd.read_csv(P.CONFIGS / "partner_priors.csv", dtype={"entrez": str, "partner_entrez": str})
    seed = int(P.seeds()[P.load_config("design_params")["dictionary_s3"]["nmf"]["seed_key"]])
    B = Dct.build(targets, gamma, drivers, partners, seed)
    inp = B.extras["inp"]
    meta = B.meta.copy()
    meta["variant_group"] = Dct.variant_groups(meta)
    nh = int(cfg["feature_hash_chars"])

    # ---- caches
    cdir = Dct.CACHE_DIR
    cdir.mkdir(parents=True, exist_ok=True)
    cache_sha = {}
    hashes, hashes_D = {}, {}
    for fam, X in B.values.items():
        X = X.astype(np.float32)
        f = cdir / f"{fam}.parquet"
        X.to_parquet(f)
        cache_sha[f.name] = sha_file(f)
        A = X.to_numpy()
        for j, c in enumerate(X.columns):
            hashes[c] = Dct.feature_hash(A[:, j], nh)
            hashes_D[c] = Dct.feature_hash(A[inp.is_D, j], nh)
    for name, df in (("sil_mechanism", B.extras["mech"]), ("expgate", B.extras["gate"]),
                     ("own_cn_log2", Dct.own_cn_log2(inp, targets))):
        f = cdir / f"{name}.parquet"
        df.to_parquet(f)
        cache_sha[f.name] = sha_file(f)
    if "nmf_H" in B.extras:
        f = cdir / "nmf_components.parquet"
        B.extras["nmf_H"].to_parquet(f)
        cache_sha[f.name] = sha_file(f)

    # ---- gene-specific map with access class
    gf = B.gene_features.copy()
    par = partners[partners["relation"].str.startswith("paralog_rank")].groupby("entrez")["partner_entrez"].apply(set)
    fam_of = dict(zip(meta["id"], meta["family"]))
    genes_of = dict(zip(meta["id"], meta["genes"]))
    gf["access_class"] = [Dct.access_class(g, f, fam_of.get(f, ""), genes_of.get(f, []), par.get(g, set()))
                          for g, f in zip(gf["entrez"], gf["feature_id"])]
    gf.to_csv(P.CONFIGS / "dictionary_gene_features.csv", index=False)

    # ---- results
    rdir = P.RESULTS / "dictionary"
    rdir.mkdir(parents=True, exist_ok=True)
    allm = B.meta_all
    summ = []
    for fam in dict.fromkeys(allm["family"]):
        a = allm[allm["family"] == fam]
        k = meta[meta["family"] == fam]
        summ.append({"family": fam, "type": ";".join(sorted(set(a["type"]))),
                     "scope": ";".join(sorted(set(k["scope"]) or set(a["scope"]))),
                     "stage": ";".join(sorted(set(k["stage"].dropna()))),
                     "n_candidates": len(a), "n_kept": len(k),
                     "n_kept_shared": int((k["scope"] == "shared").sum()),
                     "n_kept_gene_specific": int((k["scope"] == "gene_specific").sum()),
                     "min_carriers_D": float(B.fam_min.get(fam, P.load_config("design_params")["dictionary_s3"]["min_carriers"])),
                     "median_carriers_D_kept": float(k["n_carriers_D"].median()) if k["n_carriers_D"].notna().any() else np.nan,
                     "n_lineage_proxy": int(k["lineage_proxy"].sum())})
    summ = pd.DataFrame(summ)
    summ.loc[len(summ)] = {"family": "TOTAL", "n_candidates": len(allm), "n_kept": len(meta),
                           "n_kept_shared": int((meta["scope"] == "shared").sum()),
                           "n_kept_gene_specific": int((meta["scope"] == "gene_specific").sum())}
    summ["note"] = ""
    cov = B.extras.get("meth_info", {}).get("rrbs_profiled", {})
    nD = int(inp.is_D.sum())
    if cov:
        txt = f"RRBS coverage D {cov.get('dev', 0)}/{nD}, R {cov.get('test', 0)}/{int(inp.role_mask('test').sum())}; NA without RRBS (METH_absent)"
        summ.loc[summ["family"].isin(["METH", "SIL-METH"]), "note"] = txt
    fs = B.extras.get("failed_sanity", [])
    if fs:
        summ.loc[summ["family"] == "SIG", "note"] = "failed pre-declared sanity check, excluded: " + ", ".join(fs)
    summ.to_csv(rdir / "summary.csv", index=False)
    # configs/cgt_features.csv: the built ALLELE/CLIN list with carrier counts and D membership
    lst = B.extras["cgt_info"].get("list")
    if lst is not None and len(lst):
        allc = B.meta_all.set_index("id")
        lst = lst.copy()
        for c in ("n_carriers_D", "n_carriers_R", "n_NA_D", "n_NA_R"):
            lst[c] = allc[c].reindex(lst["feature_id"]).to_numpy()
        lst["in_D"] = lst["feature_id"].isin(set(meta["id"]))
        lst["reason"] = np.where(lst["in_D"], "", "fewer than min_carriers D carriers")
        nb = pd.DataFrame(B.extras["cgt_info"].get("not_built", []))
        if len(nb):
            nb["in_D"] = False
            lst = pd.concat([lst, nb], ignore_index=True)
        lst.to_csv(P.CONFIGS / "cgt_features.csv", index=False)
    B.cap_log.to_csv(rdir / "family_counts.csv", index=False)
    loc = B.locus.reset_index()
    loc.insert(2, "symbol", genesets.symbol_of().reindex(loc["entrez"]).to_numpy())
    loc.to_csv(rdir / "locus_groups.csv", index=False)
    if "nmf_H" in B.extras:
        ann = Dct.prog_annotation(B.values["PROG"], B.extras["nmf_H"], inp)
        ann.to_csv(rdir / "prog_annotation.csv", index=False)
        meta.loc[meta["family"] == "PROG", "annotation"] = meta.loc[meta["family"] == "PROG", "id"].map(
            dict(zip(ann["id"], ann["annotation"])))
    fusion_kp = known_positive_check(B.extras["fus_all"])
    kept_fus = set(B.values["FUS"].columns)  # literal-ok: family id
    for k in fusion_kp:
        k["carrier_feature_in_D"] = [f for f in k["carrier_of"] if f in kept_fus]
    checks = {"sig": B.extras["sig_info"], "meth": B.extras.get("meth_info", {}), "sil": B.extras["sil_info"],
              "cgt_features": B.extras["cgt_info"], "nmf": B.extras.get("nmf_info", {}),
              "fusion_known_positives": fusion_kp,
              "fus_carriers_D": {c: int(np.nansum(B.values["FUS"][c].to_numpy()[inp.is_D])) for c in B.values["FUS"].columns}}  # literal-ok: family id
    (rdir / "checks.json").write_text(json.dumps(checks, indent=1, default=str), encoding="utf-8")

    # ---- dictionary_v1.json
    defs = {
        "HOT": "hotspot > 0 (gene-level)", "ALLELE": "allele class (configs/cgt_allele_rules.json)",
        "LOFBI": P.load_config("design_params")["dictionary_s3"]["lof_bi_rule"],
        "AMP": "relative CN > cn_amp", "DEL": "relative CN < cn_del", "LOSS": "relative CN < cn_loss",
        "FUS": "filtered fusion, FFPM >= ffpm_min, driver (configs/fusion_drivers.csv) in its required orientation",  # literal-ok: family id
        "SIL-DEL": "SIL with CN < cn_del", "SIL-MUT": "SIL with damaging dosage >= 1 (not SIL-DEL)",
        "SIL-METH": "SIL with promoter methylation >= meth_high (not SIL-DEL/MUT)", "SIL-other": "SIL, other mechanism",
        "METH": P.load_config("design_params")["dictionary_s3"]["methylation"]["METH_rule"],
        "YABS": "Y paralog absent: female, or male with SIL of the Y gene",
        "CLIN": "clinical composite (configs/cgt_features.csv)",
        "EXPDOSE": "own log2(TPM+1) - D mean among expressed lines (>= expressed_log2tpm), 0 elsewhere",
        "EXPPAR": "log2(TPM+1) of a top-3 paralog", "PROG": "NMF programme score (D-fitted basis)",
        "PSCORE": "mean D-fitted z of member genes", "SIG": "global signature / burden / spectrum",
        "COMP": "P1 composite flag",  # literal-ok: family id
    }
    keep_cols = ["id", "family", "stage", "scope", "type", "genes", "n_carriers_D", "n_carriers_R", "n_carriers_A",
                 "n_nonnan_D", "n_NA_D", "n_NA_R", "cramers_v_lineage", "lineage_proxy", "locus_group",
                 "variant_group", "source"]
    feats = []
    src_by_family: dict = {}
    sym = genesets.symbol_of()
    for _, r in meta.iterrows():
        d = {c: r[c] for c in keep_cols if c in meta.columns}
        d["symbols"] = [sym.get(g) for g in r["genes"]]
        d["d_type"] = d.pop("type")
        d["hash"] = hashes[r["id"]]
        d["hash_D"] = hashes_D[r["id"]]
        for c in ("n_expressed_D", "annotation", "sd_D"):
            if c in meta.columns and pd.notna(r.get(c)):
                d[c] = r[c]
        for k_, v in list(d.items()):
            if isinstance(v, np.generic):
                v = v.item()
                d[k_] = v
            if v is None or (isinstance(v, float) and np.isnan(v)) or v is False or v == []:
                del d[k_]
            elif isinstance(v, float):
                d[k_] = round(v, 6)
        src_by_family.setdefault(r["family"], set()).add(d.pop("source", ""))
        feats.append(d)
    not_assayed = [{"family": r["feature_id"].split(":")[0], "feature_id": r["feature_id"], "reason": r["reason"]}
                   for r in B.extras["cgt_info"].get("not_built", [])]
    if str(B.extras["cgt_info"].get("status", "")).startswith("pending"):
        not_assayed.append({"family": "ALLELE", "reason": B.extras["cgt_info"]["status"]})
        not_assayed.append({"family": "CLIN", "reason": B.extras["cgt_info"]["status"]})
    if inp.meth is None:
        not_assayed.append({"family": "METH/SIL-METH", "reason": "no RRBS file"})
    else:
        not_assayed.append({"family": "METH/SIL-METH (partial)", "reason": "RRBS covers only CCLE lines: "
                            + json.dumps(B.extras.get("meth_info", {}).get("rrbs_profiled", {}))})
    try:
        commit = gitutil.head(P.PROJECT)
    except Exception:  # noqa: BLE001
        commit = ""
    dirty = bool(gitutil.status_porcelain(P.PROJECT, untracked=False))
    out = {
        "version": "dictionary_v1", "built_utc": t0.isoformat(timespec="seconds"), "code_commit": commit,
        "dirty_tree": dirty,
        "data_status_sha256": data_status_sha,
        "design_params_sha256": sha_file(P.CONFIGS / "design_params.json"),
        "build_params_sha256": sha_file(P.CONFIGS / "build_params.json"),
        "input_configs_sha256": {n: sha_file(P.CONFIGS / n) for n in
                                 ["eligible_genes.csv", "gamma_genes.csv", "fusion_drivers.csv", "partner_priors.csv",
                                  "positive_controls_spec.csv"] + (["cgt_features.csv"] if (P.CONFIGS / "cgt_features.csv").exists() else [])},
        "gene_features_file": {"path": "configs/dictionary_gene_features.csv",
                               "sha256": sha_file(P.CONFIGS / "dictionary_gene_features.csv")},
        "cache_dir": "cache/dictionary", "cache_sha256": cache_sha,
        "rows": {"roles": {k: int(v) for k, v in inp.pop["role"].value_counts().items()},
                 "fit_role": cfg["fit_role"]},
        "n_targets": len(targets), "n_eligible": len(elig), "n_gamma": len(gamma),
        "families": summ.replace({np.nan: None}).to_dict("records"),
        "failed_sanity": [{"id": f, "check": "spectra sanity (design_params dictionary_s3.spectra.sanity)",
                           "values": B.extras["sig_info"]["spectra_checks"]} for f in B.extras.get("failed_sanity", [])],
        "rrbs_coverage": B.extras.get("meth_info", {}).get("rrbs_profiled", {}),
        "build_count": build_count,
        "cap": {"shared_cap": P.load_config("design_params")["dictionary_s3"]["shared_feature_cap"],
                "family_min_carriers": {k: int(v) for k, v in B.fam_min.items()},
                "log": B.cap_log.replace({np.nan: None}).to_dict("records")},
        "definitions": defs,
        "feature_record_note": "keys with null/false/empty values are omitted; source is per family (sources_by_family)",
        "sources_by_family": {k: sorted(v - {""}) for k, v in src_by_family.items()},
        "sources": B.extras["sig_info"]["sources"],
        "not_assayed": not_assayed,
        "covariate_only": {"EXPGATE": "cache/dictionary/expgate.parquet (own expressed vs not; never a determinant)",
                           "own_CN": "cache/dictionary/own_cn_log2.parquet (NaN = not profiled -> 0 + CN_absent)"},
        "features": feats,
    }
    (P.CONFIGS / "dictionary_v1.json").write_text(json.dumps(out, separators=(",", ":"), default=str), encoding="utf-8")
    print(summ.to_string(index=False))
    print("fusion known positives:", [(k["cell_line_name"], k["pass"]) for k in fusion_kp])
    print("spectra checks:", B.extras["sig_info"]["spectra_checks"])
    print("cgt:", B.extras["cgt_info"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
