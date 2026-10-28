"""Interpretable feature dictionary D (DESIGN v2 section 3).

Rows: Part 1 population (D = dev, R = test, set A = noncancerous), molecular data only. Every fitted quantity
(SIL lineage medians and expressed fractions, EXPDOSE centres, PSCORE z-scalers, NMF basis, SPEC_POLE
threshold, locus groups, carrier minima, cap escalation) is fitted on D rows only; R and set-A rows are
transformed. No gene-effect value is read in this module.

Feature ids are ``FAMILY:<entrez>`` for gene-level features and ``FAMILY:<NAME>`` for named ones.
Shared features apply to every gene; gene-specific ones (own AMP/LOSS, partner LOSS, EXPDOSE, EXPPAR,
YABS) are attached per target gene in ``configs/dictionary_gene_features.csv`` together with the
self-null exclusions (DEL_g, SIL_g, LOFBI_g and DEL_h for h in g's DEL locus group).

All thresholds come from ``configs/design_params.json:dictionary_s3`` and implementation details from
``configs/build_params.json:dictionary``.
"""
from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from . import genesets, molecular, p1
from . import paths as P

CACHE_DIR = P.CACHE / "dictionary"
BINARY, CONTINUOUS = "binary", "continuous"


def dp() -> dict:
    return P.load_config("design_params")


def ds3() -> dict:
    return dp()["dictionary_s3"]


def bpd() -> dict:
    return P.load_config("build_params")["dictionary"]


# ---------------------------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------------------------


@dataclass
class Inputs:
    pop: pd.DataFrame                    # index ModelID: role, cv_fold, OncotreeLineage, PatientID, Sex
    expr: pd.DataFrame                   # log2(TPM+1), float32
    cn: pd.DataFrame                     # relative CN (linear), NaN where not profiled
    dmg: pd.DataFrame                    # damaging dosage 0/1/2, NaN where no mutation data
    hot: pd.DataFrame                    # hotspot counts
    coords: pd.DataFrame                 # entrez -> chrom, start, end
    meth: pd.DataFrame | None = None     # RRBS gene-level methylation fraction (NaN rows = not profiled)
    targets: list = field(default_factory=list)   # gene-specific target genes G* (entrez)
    gamma: list = field(default_factory=list)

    @property
    def ids(self) -> pd.Index:
        return self.pop.index

    @property
    def is_D(self) -> np.ndarray:
        return (self.pop["role"] == bpd()["fit_role"]).to_numpy()

    def role_mask(self, role) -> np.ndarray:
        return (self.pop["role"] == role).to_numpy()

    @property
    def male(self) -> np.ndarray:
        return (self.pop["Sex"] == "Male").to_numpy()

    @property
    def female(self) -> np.ndarray:
        return (self.pop["Sex"] == "Female").to_numpy()

    @property
    def sex_known(self) -> np.ndarray:
        return self.male | self.female

    @property
    def cn_profiled(self) -> np.ndarray:
        return self.cn.notna().any(axis=1).to_numpy()


def load_inputs(targets, gamma, columns_needed=None, drop_roles=()) -> Inputs:
    """Raw P1 caches reindexed to the population. ``drop_roles``: roles removed from every input (used by
    the D-only fitting invariance test)."""
    pop = molecular.population()
    if drop_roles:
        pop = pop[~pop["role"].isin(drop_roles)]
    ids = list(pop.index)

    def rd(name):
        df = p1.raw_matrix(name, ids=ids, strict=False)
        return df.reindex(ids).astype(np.float32)

    expr, cn, dmg, hot = rd("expression"), rd("cn"), rd("damaging"), rd("hotspot")
    meth, _ = molecular.rrbs_gene_matrix()
    if meth is not None:
        meth = meth.reindex(ids).astype(np.float32)
    return Inputs(pop=pop, expr=expr, cn=cn, dmg=dmg, hot=hot, coords=genesets.gene_coords(), meth=meth,
                  targets=[str(t) for t in targets], gamma=[str(g) for g in gamma])


# ---------------------------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------------------------


def _frame(arr, inp: Inputs, cols) -> pd.DataFrame:
    return pd.DataFrame(np.asarray(arr, dtype=np.float32), index=inp.ids, columns=list(cols))


def _bin(cond: np.ndarray, na: np.ndarray) -> np.ndarray:
    out = cond.astype(np.float32)
    out[na] = np.nan
    return out


def carriers(X: pd.DataFrame, mask: np.ndarray) -> np.ndarray:
    return np.nansum(X.to_numpy()[mask] == 1, axis=0).astype(int)


def n_na(X: pd.DataFrame, mask: np.ndarray) -> np.ndarray:
    return np.isnan(X.to_numpy()[mask]).sum(0).astype(int)


def chrom_of(inp: Inputs, genes) -> pd.Series:
    return inp.coords["chrom"].reindex([str(g) for g in genes])


def feature_hash(col: np.ndarray, n: int) -> str:
    return hashlib.sha256(np.ascontiguousarray(col, dtype=np.float32).tobytes()).hexdigest()[:n]


def cramers_v(B: np.ndarray, lineage: np.ndarray) -> np.ndarray:
    """Cramer's V of each binary column (NaN = missing) with a categorical vector (2 x K table)."""
    _, code = np.unique(lineage, return_inverse=True)
    L = np.zeros((len(code), code.max() + 1))
    L[np.arange(len(code)), code] = 1
    M = ~np.isnan(B)
    Bv = np.where(M, B, 0.0)
    Nk = L.T @ M
    O1 = L.T @ Bv
    O0 = Nk - O1
    n = Nk.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        p = O1.sum(0) / n
        E1, E0 = Nk * p, Nk * (1 - p)
        chi = np.where(E1 > 0, (O1 - E1) ** 2 / E1, 0).sum(0) + np.where(E0 > 0, (O0 - E0) ** 2 / E0, 0).sum(0)
        v = np.sqrt(chi / n)
    v[(p <= 0) | (p >= 1)] = np.nan
    return v


# ---------------------------------------------------------------------------------------------
# Families (each returns (values DataFrame, meta DataFrame))
# ---------------------------------------------------------------------------------------------


def _meta(cols, family, typ, genes=None, **extra) -> pd.DataFrame:
    m = pd.DataFrame({"id": list(cols), "family": family, "type": typ})
    m["genes"] = genes if genes is not None else [[] for _ in cols]
    for k, v in extra.items():
        m[k] = v
    return m


def fam_hot(inp: Inputs):
    H = inp.hot
    X = _bin(H.to_numpy() > 0, np.isnan(H.to_numpy()))
    cols = [f"HOT:{g}" for g in H.columns]
    return _frame(X, inp, cols), _meta(cols, "HOT", BINARY, [[g] for g in H.columns],
                                       source="OmicsSomaticMutationsMatrixHotspot (P1 cache hotspot.parquet)")


def fam_cn(inp: Inputs, family: str, genes) -> tuple[pd.DataFrame, pd.DataFrame]:
    s = ds3()
    genes = [g for g in genes if g in inp.cn.columns]
    C = inp.cn[genes].to_numpy()
    na = np.isnan(C)
    if family == "AMP":
        cond = C > float(s["cn_amp"])
    elif family == "DEL":
        cond = C < float(s["cn_del"])
    elif family == "LOSS":
        cond = C < float(s["cn_loss"])
    else:
        raise ValueError(family)
    cols = [f"{family}:{g}" for g in genes]
    return _frame(_bin(cond, na), inp, cols), _meta(cols, family, BINARY, [[g] for g in genes],
                                                   source="OmicsCNGeneWGS relative CN (P1 cache cn.parquet)")


def sil_stats(E_fit: np.ndarray, lineage_fit: np.ndarray, min_lines: int, expressed: float):
    """Per lineage with >= min_lines fit rows: (expressed fraction, median) per gene; plus pan-fit stats."""
    stats = {}
    for lin in np.unique(lineage_fit):
        m = lineage_fit == lin
        if m.sum() >= min_lines:
            stats[lin] = (np.nanmean(E_fit[m] >= expressed, 0), np.nanmedian(E_fit[m], 0))
    pan = (np.nanmean(E_fit >= expressed, 0), np.nanmedian(E_fit, 0))
    return stats, pan


def sil_matrix(inp: Inputs, genes, rows_ok: np.ndarray | None = None, min_lines: int | None = None):
    """SIL calls (DESIGN 3) for ``genes`` on all rows; stats fit on D rows (restricted to ``rows_ok``, e.g.
    males for Y-linked genes). Rows outside ``rows_ok`` are NaN. Returns (values, fitted stats summary)."""
    s = ds3()
    expressed = float(s["expressed_log2tpm"])
    frac_min = float(s["sil_lineage_expressed_frac"])
    delta = float(s["sil_below_lineage_median_by"])
    min_lines = int(dp()["eligibility_s2"]["lineage_min_lines"]) if min_lines is None else min_lines
    E = inp.expr[genes].to_numpy(dtype=np.float32)
    ok = np.ones(len(inp.ids), bool) if rows_ok is None else rows_ok
    fit = inp.is_D & ok
    lin = inp.pop["OncotreeLineage"].fillna("_none").to_numpy()
    stats, pan = sil_stats(E[fit], lin[fit], min_lines, expressed)
    F = np.empty_like(E)
    Md = np.empty_like(E)
    for i, l in enumerate(lin):
        f, m = stats.get(l, pan)
        F[i], Md[i] = f, m
    cond = (E < expressed) & (F >= frac_min) & (E < Md - delta)
    out = _bin(cond, np.isnan(E) | ~ok[:, None])
    return out, {"n_lineages_with_own_stats": len(stats)}


def sil_mechanism(inp: Inputs, sil: np.ndarray, genes) -> np.ndarray:
    """Mechanism code per SIL carrier, in order SIL-DEL (CN < cn_del), SIL-MUT (dosage >= 1), SIL-METH
    (promoter methylation >= meth_high; RRBS lines only, D7), SIL-other. 0 = not SIL / NaN."""
    codes = bpd()["sil_mechanism_codes"]
    C = inp.cn.reindex(columns=genes).to_numpy()
    Dm = inp.dmg.reindex(columns=genes).to_numpy()
    out = np.zeros(sil.shape, dtype=np.int8)
    car = sil == 1
    is_del = car & (C < float(ds3()["cn_del"]))
    is_mut = car & ~is_del & (Dm >= 1)
    is_meth = np.zeros_like(car)
    if inp.meth is not None:
        Mv = inp.meth.reindex(columns=genes).to_numpy()
        is_meth = car & ~is_del & ~is_mut & (Mv >= float(ds3()["methylation"]["meth_high"]))
    out[car] = codes["SIL-other"]
    out[is_meth] = codes["SIL-METH"]
    out[is_mut] = codes["SIL-MUT"]
    out[is_del] = codes["SIL-DEL"]
    return out


def fam_sil(inp: Inputs, genes):
    genes = [g for g in genes if g in inp.expr.columns]
    ch = chrom_of(inp, genes)
    ygenes = [g for g in genes if ch.get(g) == "Y"]
    agenes = [g for g in genes if ch.get(g) != "Y"]
    Xa, info = sil_matrix(inp, agenes)
    parts = [pd.DataFrame(Xa, index=inp.ids, columns=agenes)]
    if ygenes:
        Xy, _ = sil_matrix(inp, ygenes, rows_ok=inp.male, min_lines=int(bpd()["sil_y_min_male_lines_in_lineage"]))
        parts.append(pd.DataFrame(Xy, index=inp.ids, columns=ygenes))
    S = pd.concat(parts, axis=1)[genes]
    mech = sil_mechanism(inp, S.to_numpy(), genes)
    cols = [f"SIL:{g}" for g in genes]
    vals = _frame(S.to_numpy(), inp, cols)
    mech_df = pd.DataFrame(mech, index=inp.ids, columns=cols)
    meta = _meta(cols, "SIL", BINARY, [[g] for g in genes],
                 source="expression (P1 cache expression.parquet) + CN/damaging for mechanism",
                 y_linked=[g in ygenes for g in genes])
    return vals, meta, mech_df, info


SIL_FAMILIES = ["SIL-DEL", "SIL-MUT", "SIL-METH", "SIL-other"]


def sil_split(inp: Inputs, S: pd.DataFrame, meta: pd.DataFrame, mech: pd.DataFrame) -> dict:
    """Mutually exclusive SIL mechanism features (R3-2): SIL-<mech>:<g> = SIL carrier whose mechanism (precedence
    DEL > MUT > METH > other) is <mech>; NaN where SIL is NaN. One variant group per gene (VG_SIL:<g>)."""
    codes = bpd()["sil_mechanism_codes"]
    Sv = S.to_numpy()
    Mv = mech.to_numpy()
    out = {}
    genes = list(meta["genes"])
    for fam in SIL_FAMILIES:
        X = np.where(np.isnan(Sv), np.nan, (Mv == codes[fam]).astype(np.float32))
        cols = [f"{fam}:{g[0]}" for g in genes]
        m = _meta(cols, fam, BINARY, genes, source=meta["source"].iloc[0] if len(meta) else "",
                  y_linked=list(meta["y_linked"]) if len(meta) else [])
        out[fam] = (_frame(X, inp, cols), m)
    return out


def fam_lofbi(inp: Inputs, genes, sil: pd.DataFrame):
    """Biallelic LOF: dosage 2, or >= 1 with own CN < cn_loss, SIL, or a LikelyLoF call with AF >= lof_bi_af_min;
    or X-linked dosage >= 1 in males."""
    s = ds3()
    genes = [g for g in genes if g in inp.dmg.columns]
    Dm = inp.dmg[genes].to_numpy()
    C = inp.cn.reindex(columns=genes).to_numpy()
    Sl = sil.reindex(columns=[f"SIL:{g}" for g in genes]).to_numpy()
    sil_absent = ~np.isin([f"SIL:{g}" for g in genes], list(sil.columns))[None, :]
    Sl = np.where(sil_absent, 0.0, Sl)
    lof = molecular.likely_lof_af(list(inp.ids), genes, float(s["lof_bi_af_min"]))
    Lf = np.zeros(Dm.shape, bool) if lof is None else lof.reindex(index=inp.ids, columns=genes).fillna(False).to_numpy()
    xl = (chrom_of(inp, genes) == "X").to_numpy()[None, :]
    male = inp.male[:, None]
    known_sex = inp.sex_known[:, None]
    het = Dm >= 1
    c1 = Dm == 2
    c2 = het & (C < float(s["cn_loss"]))
    c3 = het & (Sl == 1)
    c4 = het & xl & male
    c5 = het & Lf
    true = c1 | c2 | c3 | c4 | c5
    unknown = het & ~true & (np.isnan(C) | np.isnan(Sl) | (xl & ~known_sex))
    X = true.astype(np.float32)
    X[unknown | np.isnan(Dm)] = np.nan
    cols = [f"LOFBI:{g}" for g in genes]
    return _frame(X, inp, cols), _meta(cols, "LOFBI", BINARY, [[g] for g in genes],
                                       source="damaging dosage + CN + SIL (P1 caches) + OmicsSomaticMutations LikelyLoF/AF")


def y_paralogs(inp: Inputs, targets) -> pd.DataFrame:
    k = int(ds3()["n_paralogs"])
    par = genesets.top_paralogs_long(targets, k)
    ch = inp.coords["chrom"]
    par = par[(par["partner_entrez"].map(ch) == "Y") & (par["entrez"].map(ch) != "Y")]
    return par


def fam_yabs(inp: Inputs, targets):
    par = y_paralogs(inp, targets)
    ys = [y for y in dict.fromkeys(par["partner_entrez"]) if y in inp.expr.columns]
    if not ys:
        return _frame(np.zeros((len(inp.ids), 0)), inp, []), _meta([], "YABS", BINARY), par
    Sy, _ = sil_matrix(inp, ys, rows_ok=inp.male, min_lines=int(bpd()["sil_y_min_male_lines_in_lineage"]))
    X = np.where(inp.female[:, None], 1.0, np.where(inp.male[:, None], Sy, np.nan)).astype(np.float32)
    cols = [f"YABS:{y}" for y in ys]
    return _frame(X, inp, cols), _meta(cols, "YABS", BINARY, [[y] for y in ys],
                                       source="Model.csv Sex + SIL of the Y paralog within males"), par


def fam_expdose(inp: Inputs, targets):
    s = ds3()
    thr = float(s["expressed_log2tpm"])
    genes = [g for g in targets if g in inp.expr.columns]
    E = inp.expr[genes].to_numpy(dtype=np.float64)
    expd = E >= thr
    ED = np.where(expd & inp.is_D[:, None], E, np.nan)
    with np.errstate(invalid="ignore"):
        mu = np.nanmean(ED, 0)
    keep = ~np.isnan(mu)
    X = np.where(expd, E - mu, 0.0)
    X[np.isnan(E)] = np.nan
    gate = expd.astype(np.float32)
    gate[np.isnan(E)] = np.nan
    genes_k = [g for g, k in zip(genes, keep) if k]
    cols = [f"EXPDOSE:{g}" for g in genes_k]
    vals = _frame(X[:, keep], inp, cols)
    gate_df = pd.DataFrame(gate, index=inp.ids, columns=genes)
    return vals, _meta(cols, "EXPDOSE", CONTINUOUS, [[g] for g in genes_k], source="expression (P1 cache)",
                       centre_D=list(mu[keep])), gate_df


def fam_exppar(inp: Inputs, paralogs):
    genes = [g for g in dict.fromkeys(paralogs) if g in inp.expr.columns]
    E = inp.expr[genes]
    sd = E[inp.is_D].std(ddof=0)
    genes = [g for g in genes if sd[g] > float(bpd()["exppar_drop_if_d_sd_le"])]
    cols = [f"EXPPAR:{g}" for g in genes]
    return _frame(inp.expr[genes].to_numpy(), inp, cols), _meta(cols, "EXPPAR", CONTINUOUS, [[g] for g in genes],
                                                                source="expression (P1 cache)")


def _parse_fusion_gene(s: pd.Series) -> pd.Series:
    ens = s.astype(str).str.extract(r"\((ENSG\d{11})")[0]
    return ens.map(genesets.ensembl_to_entrez())


def fam_fus(inp: Inputs, drivers: pd.DataFrame):
    s = ds3()["fusion"]
    cfg = bpd()
    f = p1.raw_matrix("fusion")
    f = f[f["IsDefaultEntryForModel"] == "Yes"]
    f = f[f["FFPM"] >= float(s["ffpm_min"])]
    g5 = _parse_fusion_gene(f[cfg["fusion_5prime_column"]])
    g3 = _parse_fusion_gene(f[cfg["fusion_3prime_column"]])
    f = f.assign(g5=g5.to_numpy(), g3=g3.to_numpy())
    f = f[f["g5"] != f["g3"]]
    profiled = inp.expr.notna().any(axis=1).to_numpy()
    cols, X, genes = [], [], []
    for _, d in drivers.iterrows():
        ez, orient = str(d["entrez"]), d["orientation"]
        if orient == "3prime":
            hit = f["g3"] == ez
        elif orient == "5prime":
            hit = f["g5"] == ez
        else:
            hit = (f["g3"] == ez) | (f["g5"] == ez)
        lines = set(f.loc[hit, "ModelID"])
        x = inp.ids.isin(lines).astype(np.float32)
        x[~profiled] = np.nan
        cols.append(f"FUS:{ez}")
        X.append(x)
        genes.append([ez])
    Xa = np.column_stack(X) if X else np.zeros((len(inp.ids), 0))
    return _frame(Xa, inp, cols), _meta(cols, "FUS", BINARY, genes,  # literal-ok: family id (DESIGN 3)
                                        source="OmicsFusionFiltered (P1 cache fusion.parquet) x configs/fusion_drivers.csv")


def fam_prog(inp: Inputs, seed: int):
    """NMF programmes (pinned settings) fitted on D, every row scored with model.transform."""
    from sklearn.decomposition import NMF

    n = ds3()["nmf"]
    E = inp.expr
    ED = E[inp.is_D]
    var = ED.var(ddof=0).sort_values(ascending=False, kind="stable")
    genes = list(var.index[: int(n["n_top_var_genes"])])
    model = NMF(n_components=int(n["n_components"]), init=n["init"], solver=n["solver"],
                max_iter=int(n["max_iter"]), random_state=int(seed))
    XD = ED[genes].to_numpy(dtype=np.float32)
    model.fit(XD)
    W = model.transform(np.nan_to_num(E[genes].to_numpy(dtype=np.float32)))
    k = W.shape[1]
    cols = [f"PROG:NMF{i + 1:02d}" for i in range(k)]
    H = pd.DataFrame(model.components_.astype(np.float32), index=cols, columns=genes)
    info = {"n_iter": int(model.n_iter_), "reconstruction_err": float(model.reconstruction_err_),
            "n_genes": len(genes), "seed": int(seed)}
    return _frame(W, inp, cols), _meta(cols, "PROG", CONTINUOUS, source="expression NMF (D-fitted)"), H, info


def prog_annotation(W: pd.DataFrame, H: pd.DataFrame, inp: Inputs) -> pd.DataFrame:
    """Hallmark enrichment (hypergeometric, top-N loading genes vs the NMF universe, BH within component) and
    lineage association (eta^2 on D, top lineage by mean score)."""
    from scipy.stats import hypergeom

    cfg = bpd()
    top_n = int(cfg["nmf_annotation_top_genes"])
    sets = read_gmt(genesets.resolve(cfg["hallmark_gmt"]))
    universe = set(H.columns)
    sym = genesets.symbol_to_entrez
    set_ez = {k: set(sym(v).dropna()) & universe for k, v in sets.items()}
    lin = inp.pop["OncotreeLineage"].fillna("_none").to_numpy()[inp.is_D]
    rows = []
    for c in H.index:
        top = set(H.loc[c].sort_values(ascending=False).index[:top_n])
        pv = {k: hypergeom.sf(len(top & s) - 1, len(universe), len(s), top_n) for k, s in set_ez.items() if s}
        ks = list(pv)
        p = np.array([pv[k] for k in ks])
        q = bh(p)
        j = int(np.argmin(p))
        w = W.loc[inp.is_D, c].to_numpy(dtype=np.float64)
        dfw = pd.DataFrame({"w": w, "lin": lin})
        means = dfw.groupby("lin")["w"].mean()
        ss_b = (dfw.groupby("lin")["w"].count() * (means - w.mean()) ** 2).sum()
        ss_t = ((w - w.mean()) ** 2).sum()
        rows.append({"id": c, "top_hallmark": ks[j], "hallmark_overlap": len(top & set_ez[ks[j]]),
                     "hallmark_p": float(p[j]), "hallmark_q": float(q[j]),
                     "top_lineage": str(means.idxmax()), "lineage_eta2": float(ss_b / ss_t) if ss_t > 0 else np.nan,
                     "top_genes": ";".join(genesets.symbol_of().reindex(
                         list(H.loc[c].sort_values(ascending=False).index[:10])).fillna("?"))})  # literal-ok: number of genes shown in annotation text
    out = pd.DataFrame(rows)
    out["annotation"] = np.where(out["hallmark_q"] < float(cfg["nmf_annotation_bh_q"]),
                                 out["top_hallmark"], "lineage:" + out["top_lineage"])
    return out


def bh(p) -> np.ndarray:
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1)
    return out


def read_gmt(path) -> dict:
    out = {}
    with open(path, encoding="utf-8") as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            out[f[0]] = [x for x in f[2:] if x]
    return out


def fam_pscore(inp: Inputs):
    """Programme scores: mean D-fitted z of member genes."""
    cfg = bpd()
    ps = ds3()["programme_scores"]
    gmt = read_gmt(genesets.resolve(cfg["hallmark_gmt"]))
    cols, X, genes, n_mem = [], [], [], []
    for fid, spec in ps.items():
        syms = spec["genes"] if "genes" in spec else gmt[cfg["isg_set"]]
        ez = [g for g in genesets.symbol_to_entrez(syms).dropna() if g in inp.expr.columns]
        E = inp.expr[ez].to_numpy(dtype=np.float64)
        mu = E[inp.is_D].mean(0)
        sd = E[inp.is_D].std(0)
        ok = sd > 0
        Z = (E[:, ok] - mu[ok]) / sd[ok]
        X.append(Z.mean(1))
        cols.append(fid)
        genes.append([g for g, k in zip(ez, ok) if k])
        n_mem.append(int(ok.sum()))
    return _frame(np.column_stack(X), inp, cols), _meta(cols, "PSCORE", CONTINUOUS, genes,
                                                        source="expression; members from design_params / P1 Hallmark GMT",
                                                        n_members=n_mem)


def fam_sig(inp: Inputs) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """SIG (DESIGN v3: MSI, MSIScore, ploidy, WGD, CIN, aneuploidy, LoH, TMB, DMG_BURDEN, spectra)."""
    cfg = bpd()
    ids = inp.ids
    vals, metas, sources = [], [], {}
    msi = molecular.msi_status(ids)
    fb, fs = cfg["msi_feature_ids"]["binary"], cfg["msi_feature_ids"]["score"]
    vals.append(pd.DataFrame({fb: msi["MSI"].astype(np.float32), fs: msi["MSIScore"].astype(np.float32)}, index=ids))
    metas.append(_meta([fb], "SIG", BINARY, source=f"OmicsGlobalSignatures.MSIScore >= {molecular.msi_threshold()}"))
    metas.append(_meta([fs], "SIG", CONTINUOUS, source="OmicsGlobalSignatures.MSIScore"))
    sources[fb] = {"rule": dp()["data_s1"]["msi_rule"], "threshold": molecular.msi_threshold(),
                   "lines_with_score": int(msi["MSIScore"].notna().sum()),
                   "cmp_concordance": molecular.msi_crosscheck(ids)}
    gs = molecular.global_signatures()
    if gs is not None:
        for fid, col in cfg["global_signatures_columns"].items():
            vals.append(pd.DataFrame({fid: gs[col].reindex(ids).astype(np.float32)}, index=ids))
            typ = BINARY if fid in cfg["global_signatures_binary"] else CONTINUOUS
            metas.append(_meta([fid], "SIG", typ, source=f"OmicsGlobalSignatures.{col}"))
            sources[fid] = f"OmicsGlobalSignatures.{col}"
    mb = p1.feature_block("M_burden")["M_burden"].reindex(ids).astype(np.float32)
    vals.append(pd.DataFrame({cfg["dmg_burden_feature"]: mb}, index=ids))
    metas.append(_meta([cfg["dmg_burden_feature"]], "SIG", CONTINUOUS, source="P1 M_burden = log1p damaged-gene count"))
    prof = p1.raw_matrix("damaging", columns=[]).index
    summ = molecular.mutation_summaries(ids, prof)
    checks = {}
    if len(summ.columns):
        vals.append(molecular.tmb(summ).astype(np.float32).to_frame())
        metas.append(_meta([cfg["tmb_feature"]], "SIG", CONTINUOUS, source="OmicsSomaticMutations: " + ds3()["tmb_definition"]))
        sources[cfg["tmb_feature"]] = "OmicsSomaticMutations"
        spec, fitted = molecular.spectra(summ, ids[inp.is_D])
        vals.append(spec.astype(np.float32))
        for c in spec.columns:
            metas.append(_meta([c], "SIG", BINARY if c.endswith("POLE") else CONTINUOUS,  # literal-ok: spectrum id SIG:SPEC_POLE (DESIGN 3)
                               source="OmicsSomaticMutations (reference-free; coding protein-altering calls only)"))
        checks = spectra_sanity(spec, msi, inp)
        checks["fitted"] = fitted
        sources["spectra"] = "OmicsSomaticMutations"
    return pd.concat(vals, axis=1), pd.concat(metas, ignore_index=True), {"sources": sources, "spectra_checks": checks}


def fam_meth(inp: Inputs, genes):
    """METH (D7): gene_value >= meth_high AND gene_value - D lineage median >= meth_above_lineage_median; NA for
    lines without RRBS. Lineage medians over RRBS-profiled D lines (lineages with >= lineage_min_lines such lines,
    else pan-D)."""
    m = ds3()["methylation"]
    if inp.meth is None:
        return _frame(np.zeros((len(inp.ids), 0)), inp, []), _meta([], "METH", BINARY), {}
    genes = [g for g in genes if g in inp.meth.columns]
    V = inp.meth[genes].to_numpy(dtype=np.float64)
    prof = inp.meth.notna().any(axis=1).to_numpy()
    fit = inp.is_D & prof
    lin = inp.pop["OncotreeLineage"].fillna("_none").to_numpy()
    nmin = int(dp()["eligibility_s2"]["lineage_min_lines"])
    pan = np.nanmedian(V[fit], 0)
    med = np.tile(pan, (len(lin), 1))
    n_own = 0
    for l_ in np.unique(lin[fit]):
        mm = fit & (lin == l_)
        if mm.sum() >= nmin:
            med[lin == l_] = np.nanmedian(V[mm], 0)
            n_own += 1
    with np.errstate(invalid="ignore"):
        cond = (V >= float(m["meth_high"])) & (V - med >= float(bpd()["meth_above_lineage_median"]))
    X = _bin(cond, np.isnan(V))
    cols = [f"METH:{g}" for g in genes]
    info = {"rrbs_profiled": {r: int((prof & inp.role_mask(r)).sum()) for r in ("dev", "test", "noncancerous")},
            "n_lineages_with_own_median": n_own, "n_genes_gamma_with_rrbs": len(genes)}
    return _frame(X, inp, cols), _meta(cols, "METH", BINARY, [[g] for g in genes],
                                       source="CCLE RRBS TSS1kb (max over loci)"), info


def read_cgt_features(ids) -> tuple[pd.DataFrame | None, dict]:
    """ALLELE:* and CLIN:* columns from configs/cgt_features.csv (written by the CGT track's
    scripts/build_cgt_features.py). Accepts a long table (feature_id, ModelID, value with 1/0/NA) or a wide
    table (ModelID + feature columns). Lines absent from a long table are 0; explicit NA stays NaN."""
    path = P.PROJECT / bpd()["cgt_features_file"]
    if not path.exists():
        return None, {"status": "pending: configs/cgt_features.csv absent at build time"}
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    ids = pd.Index(ids, name="ModelID")
    if {"feature_id", "ModelID", "value"} <= set(df.columns):
        feats = list(dict.fromkeys(df["feature_id"]))
        X = pd.DataFrame(0.0, index=ids, columns=feats)
        for f, grp in df.groupby("feature_id"):
            v = pd.to_numeric(grp["value"].replace({"NA": np.nan, "": np.nan}), errors="coerce")
            sub = pd.Series(v.to_numpy(), index=grp["ModelID"].to_numpy())
            sub = sub[sub.index.isin(ids)]
            X.loc[sub.index, f] = sub.to_numpy()
    else:
        key = "ModelID" if "ModelID" in df.columns else df.columns[0]
        X = df.set_index(key).replace({"NA": np.nan, "": np.nan}).apply(pd.to_numeric, errors="coerce").reindex(ids)
    X = X[[c for c in X.columns if c.split(":", 1)[0] in ("ALLELE", "CLIN")]]
    X = X.loc[:, X.notna().any()]
    return X.astype(np.float32), {"status": "read", "file": str(path.relative_to(P.PROJECT)), "n_features": X.shape[1]}


def fam_cgt(inp: Inputs):
    """ALLELE and CLIN families built from the CGT track's configs/cgt_allele_rules.json (src/biomark/cgt_features)."""
    from . import cgt_features as CF

    X, lst, info = CF.build(inp)
    out = {}
    for fam in ("ALLELE", "CLIN"):
        keys = list(lst.loc[lst["family"] == fam, "feature_id"])
        genes = [[e for e in genesets.symbol_to_entrez([g for g in s.split(";") if g]).dropna()]
                 for s in lst.loc[lst["family"] == fam, "gene"]]
        out[fam] = (X[keys], _meta(keys, fam, BINARY, genes, source="configs/cgt_allele_rules.json"))
    info["list"] = lst
    return out, info


def spectra_sanity(spec: pd.DataFrame, msi: pd.DataFrame, inp: Inputs) -> dict:
    from scipy.stats import mannwhitneyu, spearmanr

    san = ds3()["spectra"]["sanity"]
    D = inp.is_D
    lin = inp.pop["OncotreeLineage"].to_numpy()
    out = {}
    x = spec["SIG:SPEC_MMR"].to_numpy()[D]
    ref = msi["MSIScore"].to_numpy()[D]
    ok = ~np.isnan(x) & ~np.isnan(ref)
    rho = float(spearmanr(x[ok], ref[ok]).statistic)
    out["mmr_vs_msiscore_spearman"] = rho
    out["mmr_pass"] = rho > float(san["mmr_vs_msi_spearman_min"])
    uv = pd.Series(spec["SIG:SPEC_UV"].to_numpy()[D]).groupby(lin[D]).mean()
    out["uv_top_lineage"] = str(uv.idxmax())
    out["uv_pass"] = out["uv_top_lineage"] == san["uv_top_lineage"]
    tob = spec["SIG:SPEC_TOBACCO"].to_numpy()[D]
    lung = lin[D] == bpd()["tobacco_lineage"]
    a, b = tob[lung & ~np.isnan(tob)], tob[~lung & ~np.isnan(tob)]
    p = float(mannwhitneyu(a, b, alternative="greater").pvalue)
    out["tobacco_lung_vs_rest_p"] = p
    out["tobacco_pass"] = p < float(san["tobacco_lung_vs_rest_wilcoxon_p_max"])
    return out


def fam_comp(inp: Inputs):
    comp = p1.feature_block("composite")
    cols = bpd()["composite_features"]
    X = pd.DataFrame({fid: comp[c].reindex(inp.ids).astype(np.float32) for fid, c in cols.items()})
    X = (X > 0).astype(np.float32).where(X.notna())
    return X, _meta(list(cols), "COMP", BINARY, source="P1 composite block (" + ", ".join(cols.values()) + ")")  # literal-ok: family id (DESIGN 3)


# ---------------------------------------------------------------------------------------------
# Locus groups (DESIGN 3)
# ---------------------------------------------------------------------------------------------


def locus_groups(inp: Inputs, family: str, calls: pd.DataFrame) -> pd.DataFrame:
    """Connected components of genes within window_mb whose binary calls have Jaccard >= min among
    CN-profiled D lines. ``calls``: lines x genes (entrez columns). Returns entrez -> group, size, chrom, span."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    lg = ds3()["locus_groups"]
    bpm = float(bpd()["window_bp_per_mb"])
    win = float(lg["window_mb"]) * bpm
    jmin = float(lg["binary_jaccard_min"])
    rows = inp.is_D & inp.cn_profiled
    co = inp.coords.reindex(calls.columns)
    out = []
    for ch, grp in co.dropna(subset=["chrom"]).groupby("chrom"):
        grp = grp.sort_values("start")
        genes = list(grp.index)
        B = np.nan_to_num(calls.loc[rows, genes].to_numpy(dtype=np.float32)) == 1
        Bf = B.astype(np.float32)
        inter = Bf.T @ Bf
        cnt = Bf.sum(0)
        union = cnt[:, None] + cnt[None, :] - inter
        with np.errstate(invalid="ignore", divide="ignore"):
            J = np.where(union > 0, inter / union, 0)
        st, en = grp["start"].to_numpy(), grp["end"].to_numpy()
        gap = np.maximum(0, np.maximum(st[None, :] - en[:, None], st[:, None] - en[None, :]))
        A = (J >= jmin) & (gap <= win)
        np.fill_diagonal(A, False)
        i, j = np.nonzero(A)
        n = len(genes)
        _, lab = connected_components(coo_matrix((np.ones(len(i)), (i, j)), shape=(n, n)), directed=False)
        for k in np.unique(lab):
            idx = np.where(lab == k)[0]
            mb0, mb1 = st[idx].min() / bpm, en[idx].max() / bpm
            gid = f"{family}_LG:{ch}:{mb0:.2f}-{mb1:.2f}"
            for t in idx:
                out.append((genes[t], gid, len(idx), ch, mb0, mb1))
    res = pd.DataFrame(out, columns=["entrez", "locus_group", "group_size", "chrom", "start_mb", "end_mb"])
    nocoord = [g for g in calls.columns if g not in set(res["entrez"])]
    if nocoord:
        res = pd.concat([res, pd.DataFrame({"entrez": nocoord, "locus_group": [f"{family}_LG:nocoord:{g}" for g in nocoord],
                                            "group_size": 1, "chrom": None, "start_mb": np.nan, "end_mb": np.nan})])
    res.insert(0, "family", family)
    return res.set_index("entrez")


# ---------------------------------------------------------------------------------------------
# Cap (DESIGN 3)
# ---------------------------------------------------------------------------------------------


def apply_cap(meta: pd.DataFrame, min_c: int) -> tuple[pd.Series, pd.DataFrame]:
    """Per-feature minimum carriers after the cap rule on the shared dictionary. Returns (min per family,
    log of steps). Only binary shared families take part in escalation."""
    s = ds3()
    cap = int(s["shared_feature_cap"])
    levels = [int(x) for x in s["cap_escalation_min_carriers"]]
    sh = meta[meta["scope"] == "shared"]
    fam_min = {f: min_c for f in sh["family"].unique()}
    binary_fams = set(sh.loc[sh["type"] == BINARY, "family"])
    log = []

    def counts():
        keep = sh[(((sh["type"] != BINARY) & sh["cont_ok"]) | (sh["n_carriers_D"] >= sh["family"].map(fam_min)))
                  & ~sh["failed_sanity"]]
        return keep["family"].value_counts()

    step = 0
    c = counts()
    log.append({"step": step, "action": "initial", **{f"n_{k}": int(v) for k, v in c.items()}, "total": int(c.sum())})
    while c.sum() > cap:
        cand = [f for f in c.index if f in binary_fams and fam_min[f] < levels[-1]]
        if not cand:
            raise RuntimeError(f"shared dictionary {int(c.sum())} > cap {cap} with every family at the last level")
        big = max(cand, key=lambda f: (c[f], f))
        fam_min[big] = next(x for x in levels if x > fam_min[big])
        step += 1
        c = counts()
        log.append({"step": step, "action": f"{big} min carriers -> {fam_min[big]}",
                    **{f"n_{k}": int(v) for k, v in c.items()}, "total": int(c.sum())})
    return pd.Series(fam_min), pd.DataFrame(log)


# ---------------------------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------------------------


@dataclass
class Built:
    values: dict            # family -> DataFrame (kept features only)
    meta: pd.DataFrame      # one row per kept feature
    meta_all: pd.DataFrame  # every candidate feature (incl. dropped by carriers/cap)
    gene_features: pd.DataFrame
    locus: pd.DataFrame
    cap_log: pd.DataFrame
    fam_min: pd.Series
    extras: dict


def _stage_of(family: str) -> str:
    st = dp()["surrogate_s4"]["stages"]
    if family in st["G"] or family in ("METH", "SIL-METH"):
        return "G"
    if family in st["X"]:
        return "X"
    return ""


def build(targets, gamma, drivers: pd.DataFrame, partners: pd.DataFrame, seed: int,
          drop_roles=(), with_prog: bool = True) -> Built:
    """Build D. ``targets`` = gene-specific target genes (eligible U control genes), ``gamma`` = restricted
    gene set, ``drivers`` = configs/fusion_drivers.csv, ``partners`` = configs/partner_priors.csv."""
    cfg = bpd()
    s = ds3()
    min_c = int(s["min_carriers"])
    inp = load_inputs(targets, gamma, drop_roles=drop_roles)
    T = list(inp.targets)
    G = set(inp.gamma)
    tset = set(T)
    part = partners[partners["entrez"].astype(str).isin(tset)]
    partner_genes = list(dict.fromkeys(part["partner_entrez"].astype(str)))
    paralogs = list(dict.fromkeys(part.loc[part["relation"].str.startswith("paralog_rank"),
                                           "partner_entrez"].astype(str)))

    vals, metas, extras = {}, [], {}

    def add(fam, X, M):
        vals[fam] = X
        metas.append(M)

    add("HOT", *fam_hot(inp))
    add("DEL", *fam_cn(inp, "DEL", list(inp.cn.columns)))
    add("AMP", *fam_cn(inp, "AMP", [g for g in inp.cn.columns if g in G or g in tset]))
    loss_genes = list(dict.fromkeys([g for g in T + partner_genes if g in inp.cn.columns]))
    add("LOSS", *fam_cn(inp, "LOSS", loss_genes))
    Xs, Ms, mech, sil_info = fam_sil(inp, [g for g in inp.expr.columns if g in G])
    for fam, (Xf, Mf) in sil_split(inp, Xs, Ms, mech).items():
        add(fam, Xf, Mf)
    extras["sil_info"] = sil_info
    add("LOFBI", *fam_lofbi(inp, [g for g in inp.dmg.columns if g in G], Xs))
    Xm, Mm_, meth_info = fam_meth(inp, [g for g in inp.gamma])
    add("METH", Xm, Mm_)
    extras["meth_info"] = meth_info
    Xy, My, ypar = fam_yabs(inp, T)
    add("YABS", Xy, My)
    Xe, Me, gate = fam_expdose(inp, T)
    add("EXPDOSE", Xe, Me)
    add("EXPPAR", *fam_exppar(inp, paralogs))
    add("FUS", *fam_fus(inp, drivers))  # literal-ok: family id (DESIGN 3)
    if with_prog:
        Wp, Mp, H, nmf_info = fam_prog(inp, seed)
        add("PROG", Wp, Mp)
        extras.update(nmf_H=H, nmf_info=nmf_info)
    add("PSCORE", *fam_pscore(inp))
    Xg, Mg, sig_info = fam_sig(inp)
    add("SIG", Xg, Mg)
    extras["sig_info"] = sig_info
    add("COMP", *fam_comp(inp))  # literal-ok: family id (DESIGN 3)
    cgt, cgt_info = fam_cgt(inp)
    for fam, (Xc, Mc) in cgt.items():
        add(fam, Xc, Mc)
    extras["cgt_info"] = cgt_info
    extras.update(gate=gate, mech=mech, ypar=ypar, inp=inp, fus_all=vals["FUS"])  # literal-ok: family id (DESIGN 3)

    meta = pd.concat(metas, ignore_index=True)
    shared_f = set(cfg["scope_shared_families"])
    meta["scope"] = np.where(meta["family"].isin(shared_f), "shared", "gene_specific")
    is_amp = (meta["family"] == "AMP").to_numpy()
    meta.loc[is_amp, "scope"] = np.where(meta.loc[is_amp, "genes"].map(lambda g: g[0] in G),
                                         "shared", "gene_specific")
    # carrier and NA counts per role (concatenation order of metas == order of vals)
    roles = {"D": inp.is_D, "R": inp.role_mask("test"), "A": inp.role_mask("noncancerous")}
    cols = {f"n_carriers_{r}": [] for r in roles} | {f"n_NA_{r}": [] for r in roles}
    for fam, X in vals.items():
        typ = meta.loc[meta["family"] == fam, "type"].to_numpy()
        assert list(meta.loc[meta["family"] == fam, "id"]) == list(X.columns), fam
        for r, m in roles.items():
            c = carriers(X, m).astype(float)
            c[typ != BINARY] = np.nan
            cols[f"n_carriers_{r}"].extend(c)
            cols[f"n_NA_{r}"].extend(n_na(X, m))
    for k_, v in cols.items():
        meta[k_] = v
    # continuous features: >= continuous_min_nonnan non-NaN D lines and SD > 0 (EXPDOSE: >= 5 expressed D lines)
    sd, nn = [], []
    for fam, X in vals.items():
        XD = X.to_numpy(dtype=np.float64)[inp.is_D]
        nn.extend((~np.isnan(XD)).sum(0))
        with np.errstate(invalid="ignore"):
            sd.extend(np.nanstd(XD, 0) if XD.shape[1] else [])
    meta["n_nonnan_D"] = nn
    meta["sd_D"] = sd
    n_expr = gate[inp.is_D].sum(0)
    meta["n_expressed_D"] = [float(n_expr.get(g[0], np.nan)) if f == "EXPDOSE" else np.nan
                             for f, g in zip(meta["family"], meta["genes"])]
    cont_ok = (meta["n_nonnan_D"] >= int(s["continuous_min_nonnan"])) & (meta["sd_D"] > 0)
    cont_ok &= (meta["family"] != "EXPDOSE") | (meta["n_expressed_D"] >= int(s["expdose_min_expressed_lines"]))
    meta["cont_ok"] = cont_ok | (meta["type"] == BINARY)
    # pre-declared sanity checks of spectra: a failing spectrum is excluded from D (PI ruling) and recorded
    chk = sig_info.get("spectra_checks", {})
    failed = [fid for k, fid in cfg["spectra_sanity_feature"].items() if chk.get(k) is False]
    meta["failed_sanity"] = meta["id"].isin(failed)
    extras["failed_sanity"] = failed
    # minimum carriers and cap
    fam_min, cap_log = apply_cap(meta, min_c)
    need = meta["family"].map(fam_min).fillna(min_c).astype(float)
    need[meta["scope"] != "shared"] = float(min_c)          # the cap escalates shared families only
    keep = np.where(meta["type"] == BINARY, meta["n_carriers_D"] >= need, meta["cont_ok"]) & ~meta["failed_sanity"].to_numpy()
    meta["min_carriers_applied"] = np.where(meta["type"] == BINARY, need, np.nan)
    meta["kept"] = keep
    meta_all = meta.copy()
    meta = meta[keep].reset_index(drop=True)
    kept_vals = {fam: X[list(meta.loc[meta["family"] == fam, "id"])] for fam, X in vals.items()}
    # Cramer's V with (pooled) lineage on D
    from .covariates import pooled_lineage
    lin = pooled_lineage(inp.pop["OncotreeLineage"], inp.pop["role"]).to_numpy()[inp.is_D]
    v = np.full(len(meta), np.nan)
    for fam, X in kept_vals.items():
        idx = np.where((meta["family"] == fam) & (meta["type"] == BINARY))[0]
        if len(idx):
            v[idx] = cramers_v(X[list(meta.loc[idx, "id"])].to_numpy(dtype=np.float64)[inp.is_D], lin)
    meta["cramers_v_lineage"] = v
    meta["lineage_proxy"] = meta["family"].isin(SIL_FAMILIES + ["METH"]) & (meta["cramers_v_lineage"] > float(s["sil_lineage_proxy_cramers_v"]))
    meta["stage"] = meta["family"].map(_stage_of)
    # locus groups (AMP, DEL, LOSS)
    loc_parts = []
    for fam in ("AMP", "DEL", "LOSS"):
        genes = [g[0] for g in meta.loc[meta["family"] == fam, "genes"]]
        if fam == "DEL":
            Xall = vals["DEL"]
            dcar = carriers(Xall, inp.is_D & inp.cn_profiled)
            extra = [c.split(":", 1)[1] for c, n in zip(Xall.columns, dcar)
                     if n >= 1 and c.split(":", 1)[1] in tset]
            genes = list(dict.fromkeys(genes + extra))
            calls = Xall[[f"DEL:{g}" for g in genes]]
        else:
            calls = kept_vals[fam]
        loc_parts.append(locus_groups(inp, fam, calls.set_axis(genes, axis=1)))
    locus = pd.concat(loc_parts)
    lg_map = {(r.family, e): r.locus_group for e, r in locus.iterrows()}
    meta["locus_group"] = [lg_map.get((f, g[0])) if f in ("AMP", "DEL", "LOSS") else None
                           for f, g in zip(meta["family"], meta["genes"])]
    gf = gene_feature_map(meta, T, part, ypar, locus)
    return Built(values=kept_vals, meta=meta, meta_all=meta_all, gene_features=gf, locus=locus,
                 cap_log=cap_log, fam_min=fam_min, extras=extras)


def expand_ref(fid: str) -> list[str]:
    """SIL:<x> expands to its mechanism features SIL-DEL/SIL-MUT/SIL-METH/SIL-other:<x> (R3-2)."""
    fam, name = genesets.parse_feature_ref(fid)
    if fam == "SIL":
        return [f"{f}:{name}" for f in SIL_FAMILIES]
    return [fid]


def resolve_ref(ref: str) -> str:
    """'LOSS:VPS4B' -> 'LOSS:<entrez>'; named features unchanged."""
    fam, name = genesets.parse_feature_ref(ref)
    if fam in genesets.GENE_FAMILIES:
        ez = genesets.symbol_to_entrez([name]).iloc[0]
        return f"{fam}:{ez}" if isinstance(ez, str) else ref
    return ref


def gene_feature_map(meta: pd.DataFrame, targets, part: pd.DataFrame, ypar: pd.DataFrame,
                     locus: pd.DataFrame) -> pd.DataFrame:
    """Per target gene: gene-specific features and self-null exclusions (DESIGN 3)."""
    have = set(meta["id"])
    shared = set(meta.loc[meta["scope"] == "shared", "id"])
    del_lg = locus.loc[locus["family"] == "DEL", "locus_group"]
    members = {k: list(v) for k, v in del_lg.groupby(del_lg).groups.items()}
    spec = genesets.control_spec()
    c_ez = genesets.symbol_to_entrez(spec["dependency_symbol"]).to_numpy()
    ctrl: dict[str, list] = {}
    for (_, r), dep in zip(spec.iterrows(), c_ez):
        for col in ("feature_primary", "feature_secondary", "fallback_feature"):
            ref = r[col].strip()
            if ref:
                ctrl.setdefault(dep, []).extend(expand_ref(resolve_ref(ref)))
    part = part.assign(entrez=part["entrez"].astype(str), partner_entrez=part["partner_entrez"].astype(str))
    pgrp = dict(tuple(part.groupby("entrez")))
    ygrp = ypar.groupby("entrez")["partner_entrez"].apply(list).to_dict() if len(ypar) else {}
    rows = []
    for g in targets:
        def add(fid, rel, src):
            if fid in have:
                rows.append((g, fid, rel, src))
        add(f"EXPDOSE:{g}", "own_EXPDOSE", "DESIGN 3")
        add(f"AMP:{g}", "own_AMP", "DESIGN 3")
        add(f"LOSS:{g}", "own_LOSS", "DESIGN 3")
        if g in pgrp:
            for _, r in pgrp[g].iterrows():
                p_ = r["partner_entrez"]
                add(f"LOSS:{p_}", f"partner_LOSS:{r['relation']}", "configs/partner_priors.csv")
                if r["relation"].startswith("paralog_rank"):
                    add(f"EXPPAR:{p_}", f"EXPPAR:{r['relation']}", "configs/partner_priors.csv")
        for y in ygrp.get(g, []):
            add(f"YABS:{y}", "YABS", "top-3 paralog on chrY")
        for fid in ctrl.get(g, []):
            add(fid, "control_feature", "configs/positive_controls_spec.csv")
        sn = {f"DEL:{g}", f"LOFBI:{g}", f"METH:{g}"} | {f"{fam}:{g}" for fam in SIL_FAMILIES}
        if g in del_lg.index:
            sn |= {f"DEL:{h}" for h in members[del_lg[g]]}
        for fid in sorted(sn & shared):
            rows.append((g, fid, "self_null", "DESIGN 3 self-null rule"))
    gf = pd.DataFrame(rows, columns=["entrez", "feature_id", "relation", "source"])
    gf = gf.drop_duplicates(["entrez", "feature_id", "relation"]).reset_index(drop=True)
    gf.insert(1, "symbol", genesets.symbol_of().reindex(gf["entrez"]).to_numpy())
    return gf


def features_for_gene(gene: str, meta: pd.DataFrame, gene_features: pd.DataFrame) -> list[str]:
    """D_g = shared + gene-specific - self-null (control features that are self-null stay excluded)."""
    gf = gene_features[gene_features["entrez"].astype(str) == str(gene)]
    sn = set(gf.loc[gf["relation"] == "self_null", "feature_id"])
    shared = [f for f in meta.loc[meta["scope"] == "shared", "id"] if f not in sn]
    sset = set(shared)
    spec = [f for f in dict.fromkeys(gf.loc[gf["relation"] != "self_null", "feature_id"]) if f not in sset and f not in sn]
    return shared + spec


def own_cn_log2(inp: Inputs, genes) -> pd.DataFrame:
    """log2((CN + off) / (1 + off)) clipped (P1 C-block transform without the NaN fill); NaN = not profiled."""
    c = P.load_config("build_params")["covariates"]
    off, clip = float(c["cn_log2_offset"]), float(c["cn_log2_clip"])
    C = inp.cn.reindex(columns=[g for g in genes if g in inp.cn.columns])
    X = np.clip(np.log2((C.to_numpy(dtype=np.float64) + off) / (1 + off)), -clip, clip)
    return pd.DataFrame(X.astype(np.float32), index=C.index, columns=C.columns)


def load_family(family: str, columns=None) -> pd.DataFrame:
    """Cached dictionary family matrix (population rows)."""
    return pd.read_parquet(CACHE_DIR / f"{family}.parquet", columns=columns)


def load_meta() -> dict:
    import json
    with open(P.CONFIGS / "dictionary_v1.json", encoding="utf-8") as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------------------------
# Access class (DESIGN v3 4.4) and variant groups (3)
# ---------------------------------------------------------------------------------------------

_DIRECT_OWN = {"EXPDOSE", "EXPGATE", "AMP", "DEL", "LOSS", "LOFBI"}   # M2: own expression, CN, deep deletion, damaging
_DIRECT_PARALOG = {"EXPPAR", "DEL", "LOFBI"}
_POT_DIRECT = {"HOT", "LOFBI", "AMP", "DEL", "LOSS"}
_EXPR = {"EXPDOSE", "EXPPAR", "PROG", "PSCORE", "YABS", *SIL_FAMILIES}


def access_class(gene: str, feature_id: str, family: str, feature_genes, paralogs_of_gene) -> str:
    """direct / potentially-direct / expression-derived / proxy for (g, d); rules in build_params
    dictionary.access_class_rules (M2's 13 inputs: own expression, CN, deep deletion, damaging mutation; top-3
    paralog expression, deep deletion, damaging mutation)."""
    fg = set(str(x) for x in (feature_genes or []))
    g = str(gene)
    if family in _DIRECT_OWN and fg == {g}:
        return "direct"
    if family in _DIRECT_PARALOG and fg and fg <= set(paralogs_of_gene):
        return "direct"
    if family in _POT_DIRECT:
        return "potentially-direct"
    if family in _EXPR:
        return "expression-derived"
    return "proxy"


def variant_groups(meta: pd.DataFrame) -> list:
    """Variant group per feature: CN locus group; allele group 'VG_ALLELE:<entrez>' for HOT:x, ALLELE:x_* and
    CLIN features built from them (';'-joined when several genes)."""
    out = []
    for f, genes, lg in zip(meta["family"], meta["genes"], meta.get("locus_group", [None] * len(meta))):
        if f in ("AMP", "DEL", "LOSS"):
            out.append(lg)
        elif f in SIL_FAMILIES and genes:
            out.append(f"VG_SIL:{genes[0]}")
        elif f in ("HOT", "ALLELE", "CLIN") and genes:
            out.append(";".join(f"VG_ALLELE:{g}" for g in genes))
        else:
            out.append(None)
    return out
