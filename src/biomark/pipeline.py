"""Pipeline adapter: one ``Data`` container for the toy (``biomark.toysim.Toy``) and the real D data, and
``run_pipeline(data)`` producing the arms, MD-strict, audits and null calibration.

Real data (``Data.from_real``): D lines = Part 1 dev lines; genes = configs/eligible_genes.csv (eligible) or a given
list (U_ctrl); 𝒟 from configs/dictionary_v1.json + cache/dictionary/<FAMILY>.parquet + configs/dictionary_gene_
features.csv; covariates from cache/covariates.parquet; own CN / EXPGATE from cache/dictionary; row masks from
cache/effect_nan_mask.parquet (the only effect-derived input EN-pred reads); predictions = fold-centred m5 OOF;
observed effects only through ``guard.load_effects`` (scope dev), loaded lazily and only for obs arms.
"""
from __future__ import annotations

import json
import time
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from . import arms as A
from . import params as PR
from . import paths as P
from . import stats as ST
from . import en as EN

ARMS = A.ARMS


@dataclass
class Data:
    lines: pd.Index
    genes: list
    X: np.ndarray                    # n x F raw feature values
    fmeta: pd.DataFrame              # F rows: family, stage, d_type, scope, genes (tuple), variant_group, locus_group
    Dg: dict                         # gene -> positions of 𝒟_g (minus self-null)
    forced: dict                     # gene -> set of forced feature ids
    cov_df: pd.DataFrame             # canonical covariates (n rows)
    na_flags: np.ndarray             # n x k (CN_absent, MSI_unknown, METH_absent)
    own_cn: np.ndarray | None        # n x G
    expgate: np.ndarray | None       # n x G
    rows_ok: np.ndarray              # n x G (effect non-NaN)
    pred: np.ndarray                 # n x G fold-centred predictions
    folds: np.ndarray
    sub_masks: list
    patients: np.ndarray
    access_map: dict = field(default_factory=dict)      # (gene, feature) -> access class (gene-specific rows)
    access_fn: object = None
    tp53_genes: tuple = ()
    mdm2_genes: tuple = ()
    obs_loader: object = None        # callable() -> n x G observed effects (guard.load_effects for real data)
    _obs: np.ndarray | None = None
    extras: dict = field(default_factory=dict)

    def __post_init__(self):
        self.cov = ST.CovDesign(self.cov_df)
        self.drop = ST.DropInfo.build(self.fmeta, self.tp53_genes, self.mdm2_genes)
        self.Xstd = EN.standardise(self.X)
        self._gpos = {g: i for i, g in enumerate(self.genes)}

    def target(self, kind: str) -> np.ndarray:
        if kind == "pred":
            return self.pred
        if kind == "obs":
            if self._obs is None:
                self._obs = np.asarray(self.obs_loader(), dtype=np.float64)
            return self._obs
        raise ValueError(kind)

    def access(self, g, f) -> str:
        k = (g, f)
        if k in self.access_map:
            return self.access_map[k]
        return self.access_fn(g, f) if self.access_fn else ""

    def for_workers(self, include_obs: bool) -> "Data":
        """Shallow copy for worker processes: no loaders/closures; observed effects only if ``include_obs``
        (EN-pred / S1-pred workers never receive them)."""
        import copy
        if include_obs:
            self.target("obs")
        d = copy.copy(self)            # goes through __getstate__: loaders/closures dropped
        d._obs = self._obs if include_obs else None
        return d

    def __getstate__(self):
        d = dict(self.__dict__)
        d["obs_loader"] = None
        d["access_fn"] = None
        return d

    # --------------------------------------------------------------------------------------------
    # Toy adapter
    # --------------------------------------------------------------------------------------------
    @staticmethod
    def from_toy(toy, fold_centre: bool = True, genes=None) -> "Data":
        from . import dictionary as DCT
        from . import p1 as P1
        dp = PR.dp()
        lines = toy.lines
        fm0 = toy.feature_meta
        fam_sil = PR.pp().get("toy_sil_family", None)
        fmeta = pd.DataFrame(index=fm0.index)
        fmeta["family"] = fm0["family"].to_numpy()
        fmeta["stage_family"] = fm0["stage_family"].to_numpy()
        fmeta["stage"] = fm0["stage"].to_numpy()
        fmeta["d_type"] = fm0["d_type"].to_numpy()
        fmeta["scope"] = np.where(fm0["scope"] == "shared", "shared", "gene_specific")
        fmeta["genes"] = [(g,) if g else () for g in fm0["gene"]]
        lg = fm0["locus_group"].fillna("").astype(str)
        fmeta["locus_group"] = lg.where(lg != "", None).to_numpy()
        fmeta["variant_group"] = fmeta["locus_group"]
        X = toy.features[fm0.index].to_numpy(dtype=np.float64)
        pos = pd.Series(np.arange(len(fm0)), index=fm0.index)
        all_genes = list(toy.genes.index)
        el = toy_eligibility(toy)
        gl = list(el.index[el["eligible"]]) if genes is None else list(genes)
        Dg = {g: np.sort(pos.reindex(toy.features_of(g)).dropna().to_numpy().astype(np.int64)) for g in gl}
        gf = toy.gene_features
        forced = {g: set(s["feature"]) for g, s in gf.groupby("gene")}
        cv = toy.covariates
        el2 = dp["eligibility_s2"]
        vc = cv["lineage"].value_counts()
        lin = cv["lineage"].where(cv["lineage"].map(vc) >= el2["lineage_min_lines"], el2["rare_lineage_label"])
        dt = cv["ScreenDoublingTime"]
        cov = pd.DataFrame({"lineage": lin, "library_group": cv["library_group"], "cellformat": cv["CellFormat"],
                            "serumfree": cv["SerumFree"], "medium_family": cv["medium_family"],
                            "doubling_time": dt.fillna(dt.median()), "DT_unknown": cv["DT_unknown"].astype(float),
                            "MSI": cv["MSI"], "DMG_BURDEN": cv["DMG_BURDEN"], "tp53_status": cv["TP53_status"],
                            "PatientID": cv["PatientID"], "primary_disease": _pool_disease(cv["OncotreePrimaryDisease"]),
                            "lineage_raw": cv["lineage"]}, index=lines)
        na = np.column_stack([cv["CN_absent"].astype(float), cv["MSI_unknown"].astype(float), np.zeros(len(lines))])
        thr = dp["dictionary_s3"]["expressed_log2tpm"]
        own_cn = toy.cn_log2[gl].to_numpy(dtype=np.float64)
        expg = (toy.expression[gl] >= thr).astype(float).to_numpy()
        rows_ok = toy.effects[gl].notna().to_numpy()
        pr = toy.predictions[gl]
        if fold_centre:
            pr = P1.fold_centre(pr, toy.folds)
        subs = EN.make_subsamples(lines, cov["lineage"], cov["PatientID"], int(PR.seeds()["subsample_seed"]),
                                  PR.en_params()["n_pairs"])
        sub_masks = [lines.isin(s) for s in subs]
        par = {g: set(str(p).split("|")) - {""} for g, p in toy.genes["paralogs"].items()}
        fam_of = dict(zip(fmeta.index, fmeta["stage_family"]))
        gen_of = dict(zip(fmeta.index, fmeta["genes"]))

        def acc(g, f):
            return DCT.access_class(g, f, fam_of.get(f, ""), list(gen_of.get(f, ())), par.get(g, set()))

        eff = toy.effects[gl].to_numpy(dtype=np.float64)
        d = Data(lines=lines, genes=gl, X=X, fmeta=fmeta, Dg=Dg, forced=forced, cov_df=cov, na_flags=na,
                 own_cn=own_cn, expgate=expg, rows_ok=rows_ok, pred=pr.to_numpy(dtype=np.float64),
                 folds=toy.folds.reindex(lines).to_numpy(), sub_masks=sub_masks,
                 patients=cv["PatientID"].astype(str).to_numpy(), access_fn=acc, obs_loader=lambda: eff,
                 extras={"eligibility": el, "toy": True, "all_genes": all_genes})
        return d

    # --------------------------------------------------------------------------------------------
    # Real-data adapter
    # --------------------------------------------------------------------------------------------
    @staticmethod
    def from_real(genes=None, nan_mask_path=None, subsamples_path=None, save_subsamples: bool = True) -> "Data":
        from . import dictionary as DCT
        from . import genesets as GS
        from . import guard
        from . import p1 as P1
        dev = P1.dev_ids()
        meta = DCT.load_meta()
        fl = pd.DataFrame(meta["features"]).set_index("id")
        gfe = pd.read_csv(P.CONFIGS / "dictionary_gene_features.csv", dtype={"entrez": str})
        if genes is None:
            el = pd.read_csv(P.CONFIGS / "eligible_genes.csv", dtype={"entrez": str})
            genes = list(el.loc[el["eligible"].astype(bool), "entrez"])
        genes = [str(g) for g in genes]
        # feature values (D rows)
        vals = []
        for fam in dict.fromkeys(fl["family"]):
            ids = list(fl.index[fl["family"] == fam])
            v = DCT.load_family(fam, columns=ids).reindex(dev)
            vals.append(v[ids])
        Xdf = pd.concat(vals, axis=1)
        fmeta = pd.DataFrame(index=Xdf.columns)
        fmeta["family"] = fl["family"].reindex(fmeta.index).to_numpy()
        fmeta["stage_family"] = fmeta["family"]
        fmeta["stage"] = fl["stage"].reindex(fmeta.index).to_numpy()
        fmeta["d_type"] = fl["d_type"].reindex(fmeta.index).to_numpy()
        fmeta["scope"] = fl["scope"].reindex(fmeta.index).to_numpy()
        fmeta["genes"] = [tuple(map(str, g)) if isinstance(g, (list, tuple)) else () for g in fl["genes"].reindex(fmeta.index)]
        fmeta["locus_group"] = fl["locus_group"].reindex(fmeta.index).to_numpy() if "locus_group" in fl else None
        fmeta["variant_group"] = fl["variant_group"].reindex(fmeta.index).to_numpy() if "variant_group" in fl else None
        pos = pd.Series(np.arange(len(fmeta)), index=fmeta.index)
        mlike = pd.DataFrame({"id": fmeta.index, "scope": fmeta["scope"].to_numpy()})
        Dg, forced, amap = {}, {}, {}
        excl = set(PR.en_params()["forced_exclude"])
        for g, sub in gfe.groupby("entrez"):
            if g not in set(genes):
                continue
            fs = DCT.features_for_gene(g, mlike, sub)
            Dg[g] = np.sort(pos.reindex(fs).dropna().to_numpy().astype(np.int64))
            forced[g] = set(sub.loc[~sub["relation"].isin(excl), "feature_id"])
            if "access_class" in sub:
                amap.update({(g, f): a for f, a in zip(sub["feature_id"], sub["access_class"])})
        shared_pos = np.flatnonzero(fmeta["scope"].to_numpy() == "shared")
        for g in genes:
            Dg.setdefault(g, shared_pos)
        cov = real_covariates(dev)
        na = cov[PR.pp()["covariates"]["na_flags"]].astype(float).to_numpy()
        oc = pd.read_parquet(DCT.CACHE_DIR / "own_cn_log2.parquet").reindex(index=dev, columns=genes)
        eg = pd.read_parquet(DCT.CACHE_DIR / "expgate.parquet").reindex(index=dev, columns=genes)
        mpath = nan_mask_path or (P.CACHE / "effect_nan_mask.parquet")
        if not mpath.exists():
            raise FileNotFoundError(f"{mpath} missing (written by the NaN-mask job, DESIGN §1.6)")
        nm = pd.read_parquet(mpath).reindex(index=dev, columns=genes)
        rows_ok = ~nm.fillna(True).astype(bool).to_numpy()
        pr = P1.fold_centre(P1.oof_predictions("m5", columns=genes, ids=list(dev)))
        sp = subsamples_path or (P.CONFIGS / "subsamples.json")
        if sp.exists():
            subs = EN.read_subsamples(sp)
        else:
            subs = EN.make_subsamples(dev, cov["lineage"], cov["PatientID"], int(PR.seeds()["subsample_seed"]),
                                      PR.en_params()["n_pairs"])
            if save_subsamples:
                sp.write_text(json.dumps(EN.subsample_manifest(subs, PR.seeds()["subsample_seed"])), encoding="utf-8")
        sub_masks = [dev.isin(s) for s in subs]
        part = pd.read_csv(P.CONFIGS / "partner_priors.csv", dtype={"entrez": str, "partner_entrez": str})
        par = part[part["relation"].str.startswith("paralog_rank")].groupby("entrez")["partner_entrez"].apply(set)
        fam_of, gen_of = dict(zip(fmeta.index, fmeta["family"])), dict(zip(fmeta.index, fmeta["genes"]))

        def acc(g, f):
            return DCT.access_class(g, f, fam_of.get(f, ""), list(gen_of.get(f, ())), par.get(g, set()))

        tdv = PR.pp()["covariates"]["tp53_derivation"]
        tp = tuple(GS.symbol_to_entrez(tdv["gene_symbols"]).dropna().astype(str))
        md = tuple(GS.symbol_to_entrez(tdv["amp_symbols"]).dropna().astype(str))
        dev_list = list(dev)

        def obs():
            return guard.load_effects(dev_list, scope="dev", genes=genes).to_numpy(dtype=np.float64)

        return Data(lines=dev, genes=genes, X=Xdf.to_numpy(dtype=np.float64), fmeta=fmeta, Dg=Dg, forced=forced,
                    cov_df=cov, na_flags=na, own_cn=oc.to_numpy(dtype=np.float64), expgate=eg.to_numpy(dtype=np.float64),
                    rows_ok=rows_ok, pred=pr.to_numpy(dtype=np.float64), folds=P1.folds().reindex(dev).to_numpy(),
                    sub_masks=sub_masks, patients=cov["PatientID"].astype(str).to_numpy(), access_map=amap,
                    access_fn=acc, tp53_genes=tp, mdm2_genes=md, obs_loader=obs,
                    extras={"dictionary_built_utc": meta.get("built_utc"), "toy": False})


def _pool_disease(s: pd.Series) -> pd.Series:
    nmin = int(PR.dp()["confounders_s5"]["primary_disease_sensitivity_min_lines"])
    vc = s.value_counts()
    return s.where(s.map(vc) >= nmin, "_other")


def real_covariates(dev) -> pd.DataFrame:
    """Canonical covariates for D lines from cache/covariates.parquet (+ METH_absent when absent there)."""
    from . import p1 as P1
    c = pd.read_parquet(P.CACHE / "covariates.parquet").reindex(dev)
    out = pd.DataFrame(index=dev)
    for k in ("lineage", "library_group", "cellformat", "serumfree", "medium_family", "doubling_time",
              "DT_unknown", "MSI", "DMG_BURDEN", "PatientID", "CN_absent", "MSI_unknown"):
        out[k] = c[k]
    out["tp53_status"] = c[PR.pp()["covariates"]["tp53_derivation"]["covariate"]]
    vif_dropped = set(c.attrs.get("dropped", [])) if hasattr(c, "attrs") else set()
    out["screen_days"] = 0.0 if "screen_days_dev_humagne" in vif_dropped else c["screen_days_dev_humagne"]
    if "METH_absent" in c:
        out["METH_absent"] = c["METH_absent"].astype(float)
    else:
        from . import molecular
        meth, _ = molecular.rrbs_gene_matrix()
        out["METH_absent"] = (~dev.isin(meth.index if meth is not None else [])).astype(float)
    mt = P1.model_table(["OncotreePrimaryDisease"]).reindex(dev)
    out["primary_disease"] = _pool_disease(mt["OncotreePrimaryDisease"].fillna("_other"))
    out["lineage_raw"] = c["lineage_raw"] if "lineage_raw" in c else c["lineage"]
    return out


# ------------------------------------------------------------------------------------------------
# Toy eligibility (DESIGN §2 rule applied to toy columns)
# ------------------------------------------------------------------------------------------------


def toy_eligibility(toy) -> pd.DataFrame:
    from . import eligibility as ELG
    from . import p1 as P1
    dp = PR.dp()
    e = dp["eligibility_s2"]
    cv = toy.covariates
    vc = cv["lineage"].value_counts()
    lin = cv["lineage"].where(cv["lineage"].map(vc) >= e["lineage_min_lines"], e["rare_lineage_label"])
    strata = (lin.astype(str) + "|" + cv["library_group"].astype(str)).to_numpy()
    pc = P1.fold_centre(toy.predictions, toy.folds)
    min_total = int(P.load_config("build_params")["eligibility"]["min_lines_for_r"])
    r, n = ELG.r_within_strata(pc.to_numpy(), toy.effects.to_numpy(), strata, int(e["r_within_stratum_min_lines"]),
                               min_total)
    gi = toy.genes
    out = pd.DataFrame({"gene": gi.index, "r_within": r, "n_lines_r": n}, index=gi.index)
    cls = gi["variable"].astype(bool) | gi["selective"].astype(bool)
    out["tier"] = gi["tier"].to_numpy()
    out["eligible"] = cls.to_numpy() & gi["tier"].isin(e["tiers"]).to_numpy() & (out["r_within"] >= e["r_within_lin_lib_min"])
    out["eligible"] = out["eligible"].fillna(False).astype(bool)
    return out


# ------------------------------------------------------------------------------------------------
# run_pipeline
# ------------------------------------------------------------------------------------------------


def run_pipeline(data, *, fold_centre: bool = True, null_B=None, arms=ARMS, workers: int | None = None,
                 null_workers: int | None = None, audits: bool = True) -> dict:
    """Run the arms, MD-strict, audits and fixed-threshold null on a ``Data`` object or a ``toysim.Toy``.
    ``null_B=None`` -> design B (null_s7.B_base), 0 -> skip. Returns the dict of toysim.REQUIRED_OUTPUTS plus
    engine tables (s1, en, lam, md_strict_all, timings)."""
    from . import audit as AU
    from . import null as NL
    npar = PR.null_params()
    workers = npar["fg_workers"] if workers is None else workers
    null_workers = npar["workers"] if null_workers is None else null_workers
    t_all = time.perf_counter()
    is_toy = not isinstance(data, Data)
    toy = data if is_toy else None
    if is_toy:
        data = Data.from_toy(toy, fold_centre=fold_centre)
    timings = {}
    res = A.run_arms(data, arms=arms, workers=workers, timings=timings)
    out = {"meta": {"design_params_sha256": PR.sha256("design_params.json"),
                    "pipeline_params_sha256": PR.sha256("pipeline_params.json"),
                    "n_genes": len(data.genes), "n_features": int(data.X.shape[1]), "arms": list(arms),
                    "fold_centre": fold_centre}}
    el = data.extras.get("eligibility")
    if el is None:
        el = pd.DataFrame({"gene": data.genes, "eligible": True, "r_within": np.nan})
    out["eligibility"] = el.reset_index(drop=True)
    out["arms"] = res["arms"]
    out["s1"], out["en"], out["lam"] = res["s1"], res["en"], res["lam"]
    md = res.get("md_strict", pd.DataFrame(columns=["gene", "feature", "beta", "q_sc", "dropped_cov"]))
    out["md_strict"] = md
    out["md_strict_all"] = res.get("md_strict_all")
    out["md_global"] = res.get("md_global", pd.DataFrame(columns=["gene", "feature", "beta", "q", "dropped_cov"]))
    out["md_global_all"] = res.get("md_global_all")
    out["en_obs_strict"] = res.get("en_obs_strict")
    out["p_cut"] = dict(res.get("p_cut", {}), **({"MD-strict": res["p_cut_sc"]} if "p_cut_sc" in res else {}))
    s1 = res["s1"]
    out["collinear_with_design"] = (s1.loc[s1["collinear_with_design"].astype(bool), ["gene", "feature", "r2_design"]]
                                    .drop_duplicates().reset_index(drop=True) if len(s1)
                                    else pd.DataFrame(columns=["gene", "feature"]))
    if audits:
        out.update(AU.run_audits(data, res, toy=toy))
    B = PR.null_params()["B_base"] if null_B is None else int(null_B)
    t0 = time.perf_counter()
    if B > 0 and "EN-pred" in arms:
        nul = NL.run_null(data, res, B=B, workers=null_workers, genes=NL.null_genes(data))
        out["null_calibration"] = nul["summary"]
        out["null_draws"] = nul["draws"]
    else:
        out["null_calibration"] = pd.DataFrame(columns=["list", "fdr_emp_fixed", "hits_obs", "mean_hits_perm"])
    timings["null_fixed"] = time.perf_counter() - t0
    timings["total"] = time.perf_counter() - t_all
    out["timings"] = timings
    return out
