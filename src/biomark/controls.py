"""Positive controls (DESIGN v3 9): resolve configs/positive_controls_spec.csv to dictionary ids, count carriers
from the dictionary caches (molecular data only; D and R rows), gate eligibility, access class and the
maximum attainable EN selection frequency over configs/subsamples.json.

The only effect-derived input is the D NaN mask (cache/effect_nan_mask.parquet).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from . import dictionary as Dct
from . import genesets
from . import paths as P


class UnresolvableFeature(ValueError):
    """A spec feature id that is neither a dictionary id nor a declared family pattern (R3-2)."""


def declared_named() -> dict[str, set]:
    s3 = P.load_config("design_params")["dictionary_s3"]
    bp = P.load_config("build_params")["dictionary"]
    return {"ALLELE": {f"ALLELE:{x}" for x in s3["allele_classes"]},
            "CLIN": set(s3["clin_features"]),
            "SIG": set(s3["sig_features"]) | set(s3["spectra"]) - {"sanity", "input", "min_snvs"},
            "PSCORE": set(s3["programme_scores"]),
            "COMP": set(bp["composite_features"])}  # literal-ok: family id (DESIGN 3)


def resolve(ref: str, dict_ids: set) -> list[str]:
    """Spec reference -> list of dictionary ids (SIL:<x> expands to its mechanism features). Raises
    UnresolvableFeature for unknown families, unresolvable symbols or undeclared named features."""
    ref = ref.strip()
    if not ref:
        return []
    fam, name = genesets.parse_feature_ref(ref)
    if fam in genesets.GENE_FAMILIES:
        ez = genesets.symbol_to_entrez([name]).iloc[0]
        if not isinstance(ez, str):
            raise UnresolvableFeature(f"{ref}: symbol {name!r} has no Entrez id")
        return Dct.expand_ref(f"{fam}:{ez}")
    if ref in dict_ids:
        return [ref]
    dec = declared_named()
    if fam in dec and ref in dec[fam]:
        return [ref]
    raise UnresolvableFeature(f"{ref}: not a dictionary id or declared family pattern")


def _column(fid: str, fam_of: dict, cache: dict) -> pd.Series | None:
    fam = fam_of.get(fid)
    if fam is None:
        return None
    if fam not in cache:
        cache[fam] = Dct.load_family(fam)
    return cache[fam][fid]


def citation_resolved(cit: str) -> bool:
    text = (P.PROJECT / "LITERATURE.md").read_text(encoding="utf-8")
    for p in sorted((P.DOCS / "lit").glob("*.md")):
        text += p.read_text(encoding="utf-8")
    keys = [k.strip() for k in cit.split(";") if k.strip()]
    return bool(keys) and all(k in text for k in keys)


def build() -> pd.DataFrame:
    meta = Dct.load_meta()
    feats = pd.DataFrame(meta["features"])
    dict_ids = set(feats["id"])
    fam_of = dict(zip(feats["id"], feats["family"]))
    dtype_of = dict(zip(feats["id"], feats["d_type"]))
    genes_of = {i: (g if isinstance(g, list) else []) for i, g in zip(feats["id"], feats["genes"])}
    gf = pd.read_csv(P.CONFIGS / "dictionary_gene_features.csv", dtype={"entrez": str})
    fmeta = feats.rename(columns={"d_type": "type"})
    partners = pd.read_csv(P.CONFIGS / "partner_priors.csv", dtype={"entrez": str, "partner_entrez": str})
    par = partners[partners["relation"].str.startswith("paralog_rank")].groupby("entrez")["partner_entrez"].apply(set)
    mask = pd.read_parquet(P.CACHE / "effect_nan_mask.parquet")
    subs = json.loads((P.CONFIGS / "subsamples.json").read_text(encoding="utf-8"))["subsamples"]
    pop = pd.read_parquet(P.CACHE / "covariates.parquet", columns=["role"])
    isD = pop["role"] == "dev"
    isR = pop["role"] == "test"
    s3 = P.load_config("design_params")["dictionary_s3"]
    min_c = int(s3["min_carriers"])
    min_fit = int(P.load_config("design_params")["surrogate_s4"]["min_carriers_in_fit"])
    cc = P.load_config("build_params")["controls"]
    gate_n = int(cc["gate_min_nonnan_D_lines"])
    spec = genesets.control_spec()
    cache: dict = {}
    rows = []
    for _, r in spec.iterrows():
        dep = r["dependency_symbol"]
        ez = genesets.symbol_to_entrez([dep]).iloc[0]
        if not isinstance(ez, str):
            raise UnresolvableFeature(f"dependency {dep!r} has no Entrez id")
        d_g = set(Dct.features_for_gene(ez, fmeta, gf))
        nonnan = ~mask[ez] if ez in mask.columns else pd.Series(False, index=mask.index)
        n_obs = int(nonnan.sum())
        res = {}
        for col in ("feature_primary", "feature_secondary", "fallback_feature"):
            res[col] = resolve(r[col], dict_ids)
        rec = {"family": r["family"], "dependency_symbol": dep, "dependency_entrez": ez,
               "feature_primary": r["feature_primary"], "feature_secondary": r["feature_secondary"],
               "fallback_feature": r["fallback_feature"],
               "feature_primary_ids": ";".join(res["feature_primary"]),
               "feature_secondary_ids": ";".join(res["feature_secondary"]),
               "fallback_feature_ids": ";".join(res["fallback_feature"]),
               "expected_sign": r["expected_sign"], "citation": r["citation"],
               "citation_resolved": citation_resolved(r["citation"]), "in_cgt": "", "n_obs_D": n_obs}
        best = None
        for col in ("feature_primary", "feature_secondary", "fallback_feature"):
            ids = [f for f in res[col] if f in dict_ids]
            in_dg = [f for f in ids if f in d_g]
            cols = [c for c in (_column(f, fam_of, cache) for f in ids) if c is not None]
            binary = all(dtype_of.get(f) == "binary" for f in ids) if ids else True
            if cols and binary:
                X = pd.concat(cols, axis=1)
                union = X.max(axis=1, skipna=True).where(X.notna().any(axis=1))
                nD = int((union[isD] == 1).sum())
                nR = int((union[isR] == 1).sum())
                u_dev = union.reindex(mask.index)
                n_nn = int(((u_dev == 1) & nonnan).sum())
                fit_ok = [(int(((u_dev.loc[s] == 1) & nonnan.loc[s]).sum()) >= min_fit) for s in subs]
                maxfreq = float(np.mean(fit_ok))
                ok = bool(in_dg) and n_nn >= min_c
            elif cols:
                nD = nR = None
                n_nn = None
                maxfreq = 1.0
                ok = bool(in_dg)
            else:
                nD = nR = n_nn = None
                maxfreq = 0.0
                ok = False
            short = col.replace("feature_", "").replace("fallback_feature", "fallback")
            rec[f"n_carriers_D_{short}"] = nD
            rec[f"n_carriers_R_{short}"] = nR
            rec[f"n_carriers_D_nonnan_{short}"] = n_nn
            rec[f"in_D_g_{short}"] = bool(in_dg)
            rec[f"en_max_attainable_freq_{short}"] = maxfreq
            if ok and best is None and res[col]:
                best = (short, in_dg[0])
            if col == "feature_primary":
                f0 = (in_dg or ids or res[col] or [""])[0]
                rec["access_class"] = Dct.access_class(ez, f0, fam_of.get(f0, genesets.parse_feature_ref(f0)[0]),
                                                       genes_of.get(f0, []), par.get(ez, set())) if f0 else ""
        rec["n_carriers_D"] = rec["n_carriers_D_primary"]
        rec["n_carriers_R"] = rec["n_carriers_R_primary"]
        rec["en_max_attainable_freq"] = rec["en_max_attainable_freq_primary"]
        rec["gate_feature"] = "" if best is None else best[0]
        rec["gate_eligible"] = bool(best is not None and n_obs >= gate_n)
        rec["scope"] = "gate" if rec["gate_eligible"] else "untestable"
        notes = [r["note"]] if r["note"] else []
        if not rec["citation_resolved"]:
            notes.append("citation not in LITERATURE.md/docs/lit")
        pend = [f for f in res["feature_primary"] if f not in dict_ids and fam_of.get(f) is None
                and genesets.parse_feature_ref(f)[0] in ("ALLELE", "CLIN")]
        if pend:
            notes.append("primary pending configs/cgt_features.csv")
        rec["note"] = "; ".join(notes)
        rows.append(rec)
    # LINEAGE-CHECK rows (class_check)
    for _, r in genesets.lineage_check_rows().iterrows():
        ez = genesets.symbol_to_entrez([r["dependency"]]).iloc[0]
        nonnan = ~mask[ez] if isinstance(ez, str) and ez in mask.columns else pd.Series(False, index=mask.index)
        rows.append({"family": "LINEAGE-CHECK", "dependency_symbol": r["dependency"], "dependency_entrez": ez,
                     "feature_primary": "", "feature_secondary": "", "fallback_feature": "",
                     "expected_sign": "", "citation": r.get("citation", ""), "citation_resolved": None,
                     "in_cgt": "", "n_obs_D": int(nonnan.sum()), "scope": "class_check", "gate_eligible": False,
                     "access_class": "", "note": f"P1 biomarker_pairs {r['pair_id']}: lineage covariate; expected class 3"})
    out = pd.DataFrame(rows)
    def cid(r):
        f = (r["feature_primary"] or "").split(":")[-1]
        return f"{r['family']}_{r['dependency_symbol']}" + (f"_{f}" if f else "")
    out.insert(0, "control_id", out.apply(cid, axis=1))
    if out["control_id"].duplicated().any():
        raise ValueError(f"duplicate control ids: {list(out.loc[out['control_id'].duplicated(), 'control_id'])}")
    return out
