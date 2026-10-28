"""Line-level molecular annotations from the DepMap 26Q3 extra downloads (DESIGN 1.5, 3): global signatures,
somatic mutations (TMB, spectra), MSI status with the design's source order. Molecular data only: no
gene-effect values are read here, so D, R and set A rows are all allowed.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from . import datastatus, p1
from . import paths as P

CACHE_DIR = P.CACHE / "dictionary"


def _bp() -> dict:
    return P.load_config("build_params")["dictionary"]


def population() -> pd.DataFrame:
    """Rows of the dictionary: Part 1 population (dev, test, noncancerous) with role, cv_fold, lineage,
    PatientID, Sex. Sorted by ModelID."""
    sp = p1.splits()
    sp = sp[sp["role"].isin(_bp()["population_roles"])].sort_index()
    sex = p1.model_table(["Sex"])["Sex"].reindex(sp.index)
    return sp.assign(Sex=sex)


# ---------------------------------------------------------------------------------------------
# Global signatures
# ---------------------------------------------------------------------------------------------


def global_signatures() -> pd.DataFrame | None:
    """OmicsGlobalSignatures, default entries, one row per ModelID (first non-null value per column across
    the model's default rows: RNA rows hold PAM50, DNA rows MSI/CIN/...). None if the file is absent."""
    p = datastatus.find("OmicsGlobalSignatures.csv")
    if p is None:
        return None
    gs = pd.read_csv(p, index_col=0)
    gs = gs[gs["IsDefaultEntryForModel"] == "Yes"]
    num = [c for c in gs.columns if c not in ("SequencingID", "ModelID", "ModelConditionID", "IsDefaultEntryForModel",
                                              "IsDefaultEntryForMC")]
    return gs.groupby("ModelID")[num].first()


def crispr_confounders() -> pd.DataFrame | None:
    p = datastatus.find("CRISPRConfounders.csv")
    if p is None:
        return None
    cf = pd.read_csv(p, index_col=0)
    cf.index = cf.index.astype(str)
    cf.index.name = "ModelID"
    return cf


# ---------------------------------------------------------------------------------------------
# MSI (DESIGN v3 1.5: MSI = MSIScore >= data_s1.msi_threshold; Cell Model Passports = concordance check only)
# ---------------------------------------------------------------------------------------------


def msi_threshold() -> float:
    return float(P.load_config("design_params")["data_s1"]["msi_threshold"])


def cmp_msi() -> pd.Series:
    p = P.PROJECT / _bp()["cmp_models"]
    cmp = pd.read_csv(p, usecols=["BROAD_ID", "msi_status"]).dropna(subset=["BROAD_ID"])
    cmp = cmp.assign(BROAD_ID=cmp["BROAD_ID"].astype(str).str.split(";")).explode("BROAD_ID")
    cmp["BROAD_ID"] = cmp["BROAD_ID"].str.strip()
    s = cmp.dropna(subset=["msi_status"]).drop_duplicates("BROAD_ID").set_index("BROAD_ID")["msi_status"]
    return s.where(s.isin(["MSI", "MSS"]))


def subtype_tokens() -> pd.Series:
    """ModelSubtypeFeatures split into tokens (',' and 'and' separated), per ModelID."""
    m = p1.model_table(["ModelSubtypeFeatures"])["ModelSubtypeFeatures"].dropna()
    return m.str.split(r"\s*(?:,|\band\b)\s*", regex=True).map(lambda t: [x.strip() for x in t if x.strip()])


def subtype_msi() -> pd.Series:
    tok = _bp()["subtype_msi_token"]
    t = subtype_tokens()
    hit = t.map(lambda x: tok in x)
    return hit[hit]


def msi_status(ids) -> pd.DataFrame:
    """[MSI (1/0/NaN), MSIScore, msi_source]; unknown (no MSIScore) = NaN, never MSS."""
    ids = pd.Index(ids, name="ModelID")
    out = pd.DataFrame(index=ids)
    out["MSI"] = np.nan
    out["MSIScore"] = np.nan
    out["msi_source"] = None
    gs = global_signatures()
    if gs is not None:
        sc = gs["MSIScore"].reindex(ids).astype(float)
        ok = sc.notna()
        out.loc[ok, "MSI"] = (sc[ok] >= msi_threshold()).astype(float)
        out.loc[ok, "MSIScore"] = sc[ok]
        out.loc[ok, "msi_source"] = "OmicsGlobalSignatures.MSIScore"
    return out


def msi_crosscheck(ids) -> dict:
    """Concordance of the MSIScore call with Cell Model Passports and ModelSubtypeFeatures (impl notes)."""
    gs = global_signatures()
    res = {}
    if gs is None:
        return res
    thr = msi_threshold()
    sc = gs["MSIScore"].reindex(ids)
    call = (sc >= thr).where(sc.notna())
    c = cmp_msi().reindex(ids)
    j = pd.DataFrame({"call": call, "cmp": c}).dropna()
    res["vs_cmp"] = {f"score_{'MSI' if a else 'MSS'}_cmp_{b}": int(n)
                     for (a, b), n in j.groupby(["call", "cmp"]).size().items()}
    st = subtype_msi()
    st = st[st.index.isin(ids)]
    res["subtype_MSI_lines_called_MSI"] = int(call.reindex(st.index).fillna(False).astype(bool).sum())
    res["subtype_MSI_lines"] = int(len(st))
    res["threshold"] = thr
    return res


# ---------------------------------------------------------------------------------------------
# Somatic mutations: TMB and spectra (DESIGN 3 SIG, SIG-SPEC)
# ---------------------------------------------------------------------------------------------


def mutations(rebuild: bool = False) -> pd.DataFrame | None:
    """Default-entry rows of OmicsSomaticMutations (column subset, cached to cache/dictionary)."""
    p = datastatus.find("OmicsSomaticMutations.csv")
    if p is None:
        return None
    cache = CACHE_DIR / "mutations_default.parquet"
    cols = _bp()["mutation_columns"]
    if cache.exists() and not rebuild and cache.stat().st_mtime >= p.stat().st_mtime:
        have = set(pq.read_schema(cache).names)
        if set(cols) - {"IsDefaultEntryForModel"} <= have:
            return pd.read_parquet(cache)
    import pyarrow.csv as pc

    tb = pc.read_csv(p, convert_options=pc.ConvertOptions(include_columns=cols))
    m = tb.to_pandas()
    m = m[m["IsDefaultEntryForModel"] == "Yes"].drop(columns="IsDefaultEntryForModel").reset_index(drop=True)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    m.to_parquet(cache)
    return m


_COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def _rc(s: pd.Series) -> pd.Series:
    return s.map(lambda x: "".join(_COMP.get(c, "N") for c in reversed(x)))


def mutation_summaries(ids, profiled) -> pd.DataFrame:
    """Per line: n_variants, n_snv, n_indel1, n_indel, C:G SNV counts by substitution class, CC>TT DNVs.
    ``profiled``: lines with mutation data (others NaN)."""
    m = mutations()
    ids = pd.Index(ids, name="ModelID")
    if m is None:
        return pd.DataFrame(index=ids)
    m = m[m["ModelID"].isin(ids)]
    ref, alt = m["Ref"].astype(str), m["Alt"].astype(str)
    vt = m["VariantType"].astype(str)
    snv = (vt == "SNV") & (ref.str.len() == 1) & (alt.str.len() == 1)
    indel = vt.isin(["insertion", "deletion"])
    indel1 = indel & ((ref.str.len() - alt.str.len()).abs() == 1)
    # pyrimidine-collapsed SNV classes
    r1 = ref.where(~ref.isin(["G", "A"]), ref.map(_COMP))
    a1 = alt.where(~ref.isin(["G", "A"]), alt.map(_COMP))
    cls = (r1 + ">" + a1).where(snv)
    dnv = (vt == "substitution") & (ref.str.len() == 2) & (alt.str.len() == 2)
    cc_tt = dnv & (((ref == "CC") & (alt == "TT")) | ((ref == "GG") & (alt == "AA")))
    df = pd.DataFrame({"ModelID": m["ModelID"], "snv": snv, "indel": indel, "indel1": indel1,
                       "c_snv": cls.str.startswith("C>", na=False), "c_to_g": cls == "C>G",
                       "c_to_a": cls == "C>A", "dnv": dnv, "cc_tt": cc_tt, "one": True})
    g = df.groupby("ModelID").sum(numeric_only=True)
    g = g.reindex(ids).fillna(0.0)
    g = g.rename(columns={"one": "n_variants"})
    prof = ids.isin(pd.Index(profiled))
    g.loc[~prof, :] = np.nan
    return g.astype(float)


def spectra(summ: pd.DataFrame, fit_ids) -> tuple[pd.DataFrame, dict]:
    """SIG spectra (reference-free definitions of design_params dictionary_s3.spectra); lines with < min_snvs
    SNVs are NA for every spectrum. The SPEC_POLE SNV percentile is fit on ``fit_ids`` (D)."""
    sp = P.load_config("design_params")["dictionary_s3"]["spectra"]
    cfg = _bp()
    nmin = float(sp["min_snvs"])
    ok = summ["snv"] >= nmin
    out = pd.DataFrame(index=summ.index)
    out["SIG:SPEC_MMR"] = np.log1p(summ["indel1"]).where(ok)
    q = float(np.nanpercentile(summ.loc[summ.index.isin(fit_ids), "snv"], float(cfg["spec_pole_snv_percentile"])))
    ratio = summ["indel"] / summ["snv"].where(summ["snv"] > 0)
    pole = ((summ["snv"] >= q) & (ratio < float(cfg["spec_pole_indel_snv_ratio_max"]))).astype(float)
    out["SIG:SPEC_POLE"] = pole.where(ok)
    out["SIG:SPEC_APOBEC"] = (summ["c_to_g"] / summ["c_snv"].where(summ["c_snv"] > 0)).where(ok)
    tob = summ["c_to_a"] / summ["snv"].where(ok)
    out["SIG:SPEC_TOBACCO"] = tob.where(out["SIG:SPEC_POLE"] != 1)
    out["SIG:SPEC_UV"] = summ["cc_tt"] / (summ["snv"] + summ["dnv"]).where(ok)
    return out, {"spec_pole_snv_threshold_D": q}


def tmb(summ: pd.DataFrame) -> pd.Series:
    return np.log1p(summ["n_variants"]).rename(_bp()["tmb_feature"])


# ---------------------------------------------------------------------------------------------
# RRBS promoter methylation (DESIGN v3.1 amendment D7)
# ---------------------------------------------------------------------------------------------


def rrbs_gene_matrix(rebuild: bool = False) -> tuple[pd.DataFrame | None, dict]:
    """Gene-level promoter methylation: max over the gene's TSS+-1kb loci (locus_id prefix = HGNC symbol) of
    the methylation fraction; columns CCLE names mapped to ModelID via Model.csv CCLEName (unmapped dropped and
    counted). Returns (ModelID x entrez float32 or None, info)."""
    from . import genesets

    p = datastatus.find("RRBS")
    if p is None:
        return None, {"present": False}
    cache = CACHE_DIR / "rrbs_gene.parquet"
    info_p = CACHE_DIR / "rrbs_gene_info.json"
    import json
    if cache.exists() and info_p.exists() and not rebuild and cache.stat().st_mtime >= p.stat().st_mtime:
        return pd.read_parquet(cache), json.loads(info_p.read_text(encoding="utf-8"))
    r = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
    meta_cols = [c for c in r.columns[:3]]  # literal-ok: RRBS file layout: 3 leading metadata columns
    loc = r[meta_cols[0]]
    sym = loc.str.rsplit("_", n=3).str[0]  # literal-ok: locus_id = SYMBOL_chr_start_end
    vals = r.drop(columns=meta_cols).apply(pd.to_numeric, errors="coerce").astype(np.float32)
    vals["__sym"] = sym.to_numpy()
    g = vals.groupby("__sym").max()
    ez = genesets.symbol_to_entrez(list(g.index))
    keep = ez.notna().to_numpy()
    g = g[keep]
    g.index = ez[keep].to_numpy()
    g = g.groupby(level=0).max()
    ccle = p1.model_table(["CCLEName"])["CCLEName"].dropna()
    c2m = pd.Series(ccle.index, index=ccle.to_numpy())
    c2m = c2m[~c2m.index.duplicated(keep=False)]
    cols = list(g.columns)
    mapped = [c for c in cols if c in c2m.index]
    unmapped = [c for c in cols if c not in c2m.index]
    M = g[mapped].T
    M.index = c2m.reindex(mapped).to_numpy()
    M = M[~M.index.duplicated(keep="first")]
    M.index.name = "ModelID"
    M.columns = M.columns.astype(str)
    info = {"present": True, "file": p.name, "n_loci": int(len(r)), "n_columns": len(cols),
            "n_mapped": len(mapped), "n_unmapped": len(unmapped), "unmapped_columns": unmapped,
            "n_symbols": int(len(keep)), "n_symbols_unmapped": int((~keep).sum()), "n_genes": int(M.shape[1])}
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    M.astype(np.float32).to_parquet(cache)
    info_p.write_text(json.dumps(info, indent=1), encoding="utf-8")
    return M.astype(np.float32), info


def likely_lof_af(ids, genes, af_min: float) -> pd.DataFrame | None:
    """Boolean (ModelID x entrez): >= 1 LikelyLoF call with AF >= af_min (default-entry mutation calls)."""
    m = mutations()
    if m is None:
        return None
    m = m[(m["LikelyLoF"] == True) & (m["AF"] >= af_min) & m["ModelID"].isin(ids)]  # noqa: E712
    ez = pd.to_numeric(m["EntrezGeneID"], errors="coerce")
    m = m.assign(ez=ez).dropna(subset=["ez"])
    m["ez"] = m["ez"].astype(np.int64).astype(str)
    m = m[m["ez"].isin(set(genes))]
    out = pd.DataFrame(False, index=pd.Index(ids, name="ModelID"), columns=list(genes))
    for (line, g), _ in m.groupby(["ModelID", "ez"]):
        out.at[line, g] = True
    return out
