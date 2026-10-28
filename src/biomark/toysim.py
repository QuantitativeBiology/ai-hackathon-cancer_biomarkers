"""Toy simulation with planted truths: the verification answer key (DESIGN v3 Â§13.1).

Written by the verification agent, deliberately not by the pipeline author. It reads no real data: every
number comes from ``configs/toy_answer_key.json`` (generator + pass bounds), ``configs/design_params.json``
(dimensions, thresholds) and ``configs/seeds.json`` (``toy_seed``, ``toy_null_seed``).

Public API
----------
make_toy(seed, scenario)      -> Toy  (scenario "truth" = planted determinants; "null" = none)
oracle_s1 / oracle_pair / oracle_r_within / oracle_discordant / oracle_cross_lineage / oracle_amplicon
                              -> simple oracle analyses used to derive and sanity-check the pass bounds
oracle_summary, derive        -> bound derivation (``python -m biomark.toysim derive --n-seeds 20 --write``)
recovered, fdp_true, list_hits -> scoring helpers shared by the oracle and tests/test_toy_answers.py
get_pipeline, validate_outputs, REQUIRED_OUTPUTS -> the adapter interface the pipeline must expose

Pipeline adapter (``biomark.pipeline.run_pipeline``)
----------------------------------------------------
``run_pipeline(toy, *, fold_centre=True, null_B=None, arms=ARMS) -> dict[str, DataFrame | dict]`` with the keys
and columns of ``REQUIRED_OUTPUTS``. ``null_B=None`` means the design value (null_s7.B_base); ``null_B=0`` skips
the fixed-threshold null; ``arms`` restricts the arms (MD-strict needs EN-pred). Genes are Entrez-like strings,
features ``FAMILY:<entrez or name>``, lines ``TOYL####``; see ``Toy`` for every field.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.util
import json
import sys
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats

from . import paths as P

CONFIG = "toy_answer_key.json"
SCENARIOS = ("truth", "null")
ARMS = ("S1-obs", "S1-pred", "EN-obs", "EN-pred")
LISTS = ("S1-obs", "S1-pred", "EN-obs selected", "EN-pred selected", "MD-strict")
BH_LISTS = ("S1-obs", "S1-pred", "MD-strict")
CN_KEYS = ("amp", "del", "loss")
PIPELINE_MODULE = "biomark.pipeline"
PIPELINE_ENTRY = "run_pipeline"
TRUE_ROLES = frozenset({"true", "proxy_ok"})
NEUTRAL_ROLES = frozenset({"neutral"})

REQUIRED_OUTPUTS = {
    "meta": ("design_params_sha256",),
    "eligibility": ("gene", "eligible", "r_within"),
    # one row per (arm, g, d) in U; EN: rank_stat = selection frequency, selected = freq/sign rule, stage G/X;
    # S1: rank_stat = p, selected = BH q < s1_s5.bh_q (F1), stage = "".
    "arms": ("arm", "gene", "feature", "rank", "rank_stat", "selected", "beta", "stage"),
    "md_strict": ("gene", "feature", "beta", "q_sc", "dropped_cov"),
    "collinear_with_design": ("feature",),
    "locus_resolution": ("gene", "feature", "locus_group", "level", "n_discordant_carriers"),
    "confounding_delta": ("gene", "feature", "block", "flagged"),
    "amplicon_test": ("gene", "feature", "label"),
    "cross_lineage": ("gene", "feature", "label"),
    "gene_classes": ("gene", "failure_class"),
    "self_null": ("gene", "feature", "frac_dependent"),
    "null_calibration": ("list", "fdr_emp_fixed", "hits_obs", "mean_hits_perm"),
    "cgt_scores": ("arm", "test_id", "hit", "sign_correct"),
    "gate": ("control_id", "recovered"),
}


# ---------------------------------------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------------------------------------


def load_toy_config() -> dict:
    return P.load_config(CONFIG)


def _sha(path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def design_params_sha256() -> str:
    return _sha(P.CONFIGS / "design_params.json")


def toy_config_sha256() -> str:
    return _sha(P.CONFIGS / CONFIG)


def pass_bounds() -> dict:
    """Pre-declared pass bounds (configs/toy_answer_key.json: pass)."""
    return load_toy_config()["pass"]


# ---------------------------------------------------------------------------------------------------------
# The toy object
# ---------------------------------------------------------------------------------------------------------


@dataclass
class Toy:
    """Toy dataset in the real pipeline's conventions (rows = lines ``TOYL####``, gene columns = Entrez strings).

    genes          gene_info indexed by entrez: symbol, chrom, start_mb, tier (m5-like), variable, selective,
                   paralogs ('|'-joined entrez), locus_AMP / locus_DEL / locus_LOSS (group id or '').
    features       lines x feature ids ``FAMILY:<entrez or name>``; binary 0/1 or continuous; NaN where the
                   source is missing (CN families on CN_absent lines, SIG:MSI on MSI unknown, COMP:TP53MUT on
                   TP53 unknown). Shared ð’Ÿ (scope 'shared') + gene-specific columns (EXPDOSE, LOSS, EXPPAR).
    feature_meta   indexed by feature id: family, gene, d_type (binary/continuous), stage (G/X), stage_family,
                   scope (shared/gene), locus_group, n_carriers, n_na.
    gene_features  rows (gene, feature, relation in own/paralog): gene-specific membership of ð’Ÿ_g.
    covariates     lineage, library_group, CellFormat, SerumFree, medium_family, ScreenDoublingTime (NaN =
                   unknown), DT_unknown, MSI {MSI, MSS, unknown}, DMG_BURDEN, TMB, TP53_status {mut, WT,
                   unknown}, CN_absent, MSI_unknown, PatientID, OncotreePrimaryDisease, role.
    expression     lines x genes log2(TPM+1) (EXPGATE = expression >= dictionary_s3.expressed_log2tpm).
    cn_log2        lines x genes log2 relative CN (NaN on CN_absent lines); own-CN covariate source.
    effects        observed effects (truth + library offset + noise; NaN for one library's missing genes).
    predictions    OOF-like predictions WITH per-fold offsets (not fold-centred: the pipeline must centre).
    folds          cv_fold per line (patient-grouped, lineage-stratified).
    patients       PatientID per line.
    controls       positive-control table (DESIGN Â§9.1 subset).
    cgt            toy clinical ground truth (DESIGN Â§10.2 subset); the true signs are only in answer_key.
    truth          noiseless true effects (never to be read by the pipeline).
    answer_key     dict: plants, pairs, colocated, fold_artefact, self_null, culture, cgt_truth, bounds, ...
    """

    seed: int
    scenario: str
    predictor: str
    genes: pd.DataFrame
    features: pd.DataFrame
    feature_meta: pd.DataFrame
    gene_features: pd.DataFrame
    covariates: pd.DataFrame
    expression: pd.DataFrame
    cn_log2: pd.DataFrame
    effects: pd.DataFrame
    predictions: pd.DataFrame
    folds: pd.Series
    patients: pd.Series
    controls: pd.DataFrame
    cgt: pd.DataFrame
    truth: pd.DataFrame
    answer_key: dict = field(default_factory=dict)
    params: dict = field(default_factory=dict)

    @property
    def lines(self) -> pd.Index:
        return self.features.index

    def plants(self, *kinds) -> pd.DataFrame:
        pl = self.answer_key["plants"]
        return pl if not kinds else pl[pl["kind"].isin(kinds)]

    def features_of(self, gene: str) -> list[str]:
        """ð’Ÿ_g before carrier filters: shared minus self-null, plus gene-specific."""
        fm = self.feature_meta
        sn = set(self_null_features(self, gene))
        shared = [f for f in fm.index[fm["scope"] == "shared"] if f not in sn]
        own = self.gene_features.loc[self.gene_features["gene"] == gene, "feature"].tolist()
        return list(dict.fromkeys(shared + own))


def self_null_features(toy: Toy, gene: str) -> list[str]:
    """DESIGN Â§3 self-null set of g: DEL_g, SIL_g, LOFBI_g and DEL_h for h in g's DEL locus group."""
    fam = toy.params["families"]
    out = [f"{fam['del']}:{gene}", f"{fam['sil']}:{gene}", f"{fam['lof']}:{gene}"]
    grp = toy.genes.at[gene, "locus_DEL"]
    if grp:
        out += [f"{fam['del']}:{h}" for h in toy.genes.index[toy.genes["locus_DEL"] == grp] if h != gene]
    return [f for f in out if f in toy.feature_meta.index]


# ---------------------------------------------------------------------------------------------------------
# Generator
# ---------------------------------------------------------------------------------------------------------


class _Ctx:
    """Mutable bag of generator state."""


def _u(rng, x) -> float:
    """Uniform draw from [lo, hi] when x is a pair, else x itself (no draw)."""
    if isinstance(x, (list, tuple)):
        return float(rng.uniform(x[0], x[1]))
    return float(x)


def _largest_remainder(w, n) -> np.ndarray:
    q = np.asarray(w, float) / np.sum(w) * n
    base = np.floor(q).astype(int)
    base[np.argsort(-(q - base))[: n - base.sum()]] += 1
    return base


def _lineages(dp: dict, L: int):
    """12 toy lineages, round-robin over the design's organ-system groups (so every group is represented)."""
    groups = dp["subtypes_s8"]["lineage_groups"]
    out, depth = [], 0
    while len(out) < L and depth < max(len(v) for v in groups.values()):
        for gname, members in groups.items():
            if depth < len(members) and len(out) < L:
                out.append((members[depth], gname))
        depth += 1
    if len(out) < L:
        raise ValueError(f"design lineage_groups name only {len(out)} lineages < {L}")
    return [o[0] for o in out], {o[0]: o[1] for o in out}


def _choice(rng, spec: dict, n: int) -> np.ndarray:
    return rng.choice(np.asarray(spec["levels"], dtype=object), size=n, p=spec["p"])


def make_toy(seed: int | None = None, scenario: str = "truth", n_lines: int | None = None,
             predictor: str | None = None) -> Toy:
    """Build a toy dataset. ``scenario='truth'`` plants every determinant of the answer key; ``'null'`` uses the
    identical world (same draws: lines, covariates, features, noise, lineage/library/latent effects, fold offsets,
    self-null and lineage-only genes) with every molecular/culture determinant set to zero. Default seed:
    seeds.json toy_seed (truth) or toy_null_seed (null). ``n_lines`` overrides the line count (e.g. 1022 for
    real-D-dimension NC bounds); ``predictor`` in {'ridge' (DESIGN Â§13.1, default), 'smooth'}."""
    if scenario not in SCENARIOS:
        raise ValueError(f"scenario must be one of {SCENARIOS}")
    c = _Ctx()
    c.cfg = load_toy_config()
    c.g, c.fam, c.nf = c.cfg["generator"], c.cfg["families"], c.cfg["named_features"]
    c.dp = P.load_config("design_params")
    t13 = c.dp["toy_s13"]
    sd = P.seeds()
    if seed is None:
        seed = sd[t13["seed_key"]] if scenario == "truth" else sd[t13["null_seed_key"]]
    c.seed, c.scenario = int(seed), scenario
    c.on = 1.0 if scenario == "truth" else 0.0
    c.n = int(n_lines or t13["n_lines"])
    c.G, c.L, c.K = int(t13["n_genes"]), int(t13["n_lineages"]), int(t13["n_folds"])
    c.n_shared = int(t13["n_features"])
    c.predictor = predictor or c.g["predictor"]["default"]
    c.rw, c.rp, c.rn = (np.random.default_rng(s) for s in np.random.SeedSequence(c.seed).spawn(3))  # literal-ok: 3 streams
    _lines(c)
    _genome(c)
    _plan(c)
    _features(c)
    _effects(c)
    _predict(c)
    return _assemble(c)


def _lines(c) -> None:
    g, n, L, K, rw = c.g, c.n, c.L, c.K, c.rw
    c.lin_names, c.lin_group = _lineages(c.dp, L)
    w = np.asarray(g["lineage_weights"], float)
    if len(w) < L:
        w = np.r_[w, np.full(L - len(w), w.mean())]
    c.lin = rw.permutation(np.repeat(np.arange(L), _largest_remainder(w[:L], n)))
    c.ids = pd.Index([f"TOYL{i + 1:04d}" for i in range(n)], name="ModelID")
    c.lin_size = np.bincount(c.lin, minlength=L)
    # patients: a few same-lineage pairs share a patient
    patient = np.arange(n)
    pools = {l: list(rw.permutation(np.flatnonzero(c.lin == l))) for l in range(L)}
    target, made = int(round(g["patient_pair_frac"] * n / 2)), 0
    order = list(rw.permutation(L))
    while made < target and any(len(v) >= 2 for v in pools.values()):
        for l in order:
            if made < target and len(pools[l]) >= 2:
                a, b = pools[l].pop(), pools[l].pop()
                patient[b] = patient[a]
                made += 1
    c.patient = patient
    # folds: patient-grouped, lineage-stratified
    fold = np.zeros(n, int)
    for l in range(L):
        pats = rw.permutation(np.unique(patient[c.lin == l]))
        off = int(rw.integers(K))
        for i, p in enumerate(pats):
            fold[patient == p] = (off + i) % K
    c.fold = fold
    # covariates
    grp = np.array([c.lin_group[c.lin_names[l]] for l in c.lin], dtype=object)
    susp = np.isin(grp, g["suspension_groups"])
    libs = np.asarray(g["library_groups"], dtype=object)
    lib_other = rw.choice(libs, size=n, p=g["library_other_p"])
    lib_nested = np.where(rw.random(n) < g["library_nested_p"], g["library_nested_group"], lib_other)
    c.lib = np.where(np.isin(c.lin, g["library_nested_lineage_idx"]), lib_nested, lib_other)
    c.cellformat = np.where(susp, _choice(rw, g["cellformat_suspension_lineages"], n), _choice(rw, g["cellformat_other"], n))
    pT = np.array([g["serumfree_T_p_by_group"].get(x, g["serumfree_T_p_default"]) for x in grp])
    c.serum_true = rw.random(n) < pT
    sf_unknown = rw.random(n) < g["serumfree_unknown_p"]
    c.serumfree = np.where(sf_unknown, "unknown", np.where(c.serum_true, "T", "F")).astype(object)
    c.medium = np.where(susp, _choice(rw, g["medium_suspension"], n), _choice(rw, g["medium_other"], n))
    dt = g["doubling_time"]
    c.dt = rw.normal(dt["mean_h"], dt["sd_h"], n) + dt["suspension_shift_h"] * (c.cellformat == "suspension")
    c.dt_unknown = rw.random(n) < dt["unknown_p"]
    rate = np.full(n, g["msi_rate_default"])
    for k, v in g["msi_rate_by_lineage_idx"].items():
        rate[c.lin == int(k)] = v
    c.msi_true = rw.random(n) < rate
    c.msi_unknown = rw.random(n) < g["msi_unknown_p"]
    c.msi_label = np.where(c.msi_unknown, "unknown", np.where(c.msi_true, "MSI", "MSS")).astype(object)
    db = g["dmg_burden"]
    c.dmg = rw.normal(db["mean"], db["sd"], n) + db["msi_shift"] * c.msi_true
    c.tmb = np.exp(c.dmg - g["tmb_offset"] + rw.normal(0, g["tmb_log_sd"], n))
    pt = np.array([g["tp53_rate_by_group"].get(x, g["tp53_rate_default"]) for x in grp])
    c.tp53_true = rw.random(n) < pt
    c.tp53_unknown = rw.random(n) < g["tp53_unknown_p"]
    c.tp53_label = np.where(c.tp53_unknown, "unknown", np.where(c.tp53_true, "mut", "WT")).astype(object)
    c.cn_absent = rw.random(n) < g["cn_absent_p"]
    sub = np.where(rw.random(n) < g["subtype_A_p"], "A", "B")
    c.disease = np.array([f"{c.lin_names[l]}:{s}" for l, s in zip(c.lin, sub)], dtype=object)
    c.culture_var = {"suspension": c.cellformat == "suspension", "serum_free": c.serum_true}


def _genome(c) -> None:
    g, G, n, rw = c.g, c.G, c.n, c.rw
    nc = int(g["n_chrom"])
    gpc = int(np.ceil(G / nc))
    idx = np.arange(G)
    c.chrom, c.jpos = idx // gpc + 1, idx % gpc
    c.pos = g["first_gene_mb"] + g["gene_spacing_mb"] * c.jpos
    c.entrez = np.array([str(int(g["entrez_base"]) + i) for i in idx], dtype=object)
    c.symbol = np.array([f"TG{i + 1:03d}" for i in idx], dtype=object)
    cn = g["cn"]
    arm_id = (c.chrom - 1) * 2 + (c.jpos >= gpc // 2)
    A = 2 * nc
    arm_off = (rw.normal(0, cn["arm_sd"], (n, A)) + (rw.random((n, A)) < cn["arm_loss_p"]) * cn["arm_loss_shift"]
               + (rw.random((n, A)) < cn["arm_gain_p"]) * cn["arm_gain_shift"])
    log2 = arm_off[:, arm_id] + rw.normal(0, cn["gene_sd"], (n, G))
    first, step = int(g["locus_centre_first"]), int(g["locus_centre_step"])
    cands = [i for i in idx if c.jpos[i] >= first and (c.jpos[i] - first) % step == 0
             and c.jpos[i] + 1 < gpc and i + 1 < G]
    amp_kinds = [k for k, m in g["amp_loci"] for _ in range(m)]
    del_kinds = [k for k, m in g["del_loci"] for _ in range(m)]
    chosen = c.rp.choice(cands, size=len(amp_kinds) + len(del_kinds), replace=False)
    c.loci = []
    for k, (core, kind) in enumerate(zip(chosen, amp_kinds + del_kinds)):
        is_amp = k < len(amp_kinds)
        prm = g["loci"][kind]
        p, q, nbo = _u(rw, prm["p"]), _u(rw, prm["q"]), float(prm["nb_only"])
        lo, hi = cn["amp_shift"] if is_amp else cn["del_shift"]
        sgn = 1.0 if is_amp else -1.0
        ev = rw.random(n) < p
        log2[:, core] += sgn * ev * rw.uniform(lo, hi, n)
        for nb in (core - 1, core + 1):
            hit = (ev & (rw.random(n) < q)) | (~ev & (rw.random(n) < nbo))
            log2[:, nb] += sgn * hit * rw.uniform(lo, hi, n)
        c.loci.append({"kind": kind, "key": "amp" if is_amp else "del", "core": int(core),
                       "left": int(core - 1), "right": int(core + 1)})
    c.loci_by_kind = {}
    for loc in c.loci:
        c.loci_by_kind.setdefault(loc["kind"], []).append(loc)
    d3 = c.dp["dictionary_s3"]
    c.log2 = log2
    c.call = {"amp": log2 > np.log2(d3["cn_amp"]), "del": log2 < np.log2(d3["cn_del"]),
              "loss": log2 < np.log2(d3["cn_loss"])}
    lg = d3["locus_groups"]
    c.locus = {k: _locus_groups(c.call[k][~c.cn_absent], c.chrom, c.pos, lg["window_mb"], lg["binary_jaccard_min"],
                                c.symbol, c.fam[k]) for k in CN_KEYS}


def _locus_groups(call, chrom, pos, window, jmin, symbol, fam) -> np.ndarray:
    """DESIGN Â§3 locus groups: connected components of genes within ``window`` Mb with Jaccard >= jmin."""
    X = call.astype(np.int64)
    inter = X.T @ X
    cnt = X.sum(0)
    union = cnt[:, None] + cnt[None, :] - inter
    with np.errstate(divide="ignore", invalid="ignore"):
        J = np.where(union > 0, inter / np.maximum(union, 1), 0.0)
    near = (chrom[:, None] == chrom[None, :]) & (np.abs(pos[:, None] - pos[None, :]) <= window)
    G = len(chrom)
    parent = list(range(G))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for i, j in np.argwhere(np.triu(near & (J >= jmin), 1)):
        ri, rj = find(int(i)), find(int(j))
        if ri != rj:
            parent[max(ri, rj)] = min(ri, rj)
    roots = np.array([find(i) for i in range(G)])
    labels = np.full(G, "", dtype=object)
    for r in np.unique(roots):
        mem = np.flatnonzero(roots == r)
        if len(mem) >= 2:
            labels[mem] = f"{fam}:{symbol[mem[0]]}-{symbol[mem[-1]]}"
    return labels


def _plan(c) -> None:
    pl, e, rp = c.g["plants"], c.g["effects"], c.rp
    locus_genes = {x for loc in c.loci for x in (loc["left"], loc["core"], loc["right"])}
    pool = [int(i) for i in rp.permutation(c.G) if int(i) not in locus_genes]

    def take(k):
        out = pool[: int(k)]
        del pool[: int(k)]
        if len(out) < int(k):
            raise ValueError("toy has too few genes for the planted roles")
        return out

    pan = pl["pan_lineage"]
    R = {}
    R["hot_own"] = take(pan["hot_own"])
    R["hot_path_dep"], R["hot_path_x"] = take(pan["hot_pathway"]), take(pan["hot_pathway"])
    R["lof_B"], R["lof_A"] = take(pan["lofbi_paralog"]), take(pan["lofbi_paralog"])
    R["loss_B"], R["loss_A"] = take(pan["loss_paralog"]), take(pan["loss_paralog"])
    R["sil_B"], R["sil_A"] = take(pan["sil_partner"]), take(pan["sil_partner"])
    R["amp_partner_dep"] = take(pan["amp_partner"])
    R["restr_dep"], R["restr_x"] = take(pl["lineage_restricted"]), take(pl["lineage_restricted"])
    R["msi_dep"], R["msi_pass"] = take(pl["msi_like"]), take(e["n_msi_passengers"])
    R["cult_dep"], R["cult_x"] = take(pl["culture_confounded"]), take(pl["culture_confounded"])
    R["coloc_B"], R["tight_B"] = take(pl["co_located"]), take(pl["co_located_tight"])
    R["fold_dep"], R["fold_x"] = take(pl["fold_offset_artefact"]), take(pl["fold_offset_artefact"])
    R["proxy_B"], R["proxy_A"] = take(pl["genomic_with_expression_proxy"]), take(pl["genomic_with_expression_proxy"])
    R["prog_dep"] = take(pl["prog_driven"])
    R["p53_dep"] = take(pl["p53_w"])
    R["blind_dep"], R["blind_x"] = take(pl["predictor_blind"]), take(pl["predictor_blind"])
    R["wrong_dep"], R["wrong_x"] = take(pl["wrong_sign_cgt"]), take(pl["wrong_sign_cgt"])
    R["six_dep"], R["six_x"] = take(pl["six_carrier_strong"]), take(pl["six_carrier_strong"])
    R["ldum_dep"], R["ldum_x"] = take(pl["lineage_dummy_feature"]), take(pl["lineage_dummy_feature"])
    R["lin_only"] = take(pl["lineage_only"])
    R["low_tier"] = take(pl["low_tier"])
    R["nan_genes"] = take(c.g["noise"]["n_library_nan_genes"])
    c.R = R
    c.background = pool
    feat_keys = ("hot_own", "hot_path_x", "lof_A", "loss_A", "sil_A", "restr_x", "msi_pass", "cult_x", "fold_x",
                 "proxy_A", "blind_x", "wrong_x", "six_x", "ldum_x")
    c.planted_feature_genes = {x for k in feat_keys for x in R[k]}
    c.bg_feature_pool = [i for i in range(c.G) if i not in locus_genes and i not in c.planted_feature_genes]
    # paralogs (Ensembl-top-3-like); planted paralog pairs put A first
    lo, hi = c.g["n_paralogs"]
    par = []
    for i in range(c.G):
        k = int(rp.integers(lo, hi + 1))
        others = rp.choice(np.delete(np.arange(c.G), i), size=k, replace=False)
        par.append([int(x) for x in others])
    forced = list(zip(R["lof_B"], R["lof_A"])) + list(zip(R["loss_B"], R["loss_A"])) + \
        list(zip(R["sil_B"], R["sil_A"])) + list(zip(R["proxy_B"][:-1], R["proxy_A"][:-1]))
    for b, a in forced:
        par[b] = ([a] + [x for x in par[b] if x != a])[:hi]
    c.paralogs = par
    # lineages usable for lineage-specific plants (>= lineage_min_lines_for_plants lines)
    c.big_lineages = [l for l in range(c.L) if c.lin_size[l] >= e["lineage_min_lines_for_plants"]]


def _features(c) -> None:
    g, e, fam, nf, rw, n, L = c.g, c.g["effects"], c.fam, c.nf, c.rw, c.n, c.L
    F, M = {}, {}
    stagesG = set(c.dp["surrogate_s4"]["stages"]["G"])
    E_ = c.entrez

    def add(name, values, family, gene=None, d_type="binary", scope="shared", stage_family=None):
        if name in F:
            raise KeyError(f"duplicate feature {name}")
        F[name] = np.asarray(values, float)
        sf = stage_family or family
        M[name] = {"family": family, "gene": "" if gene is None else E_[gene], "d_type": d_type,
                   "stage_family": sf, "stage": "G" if sf in stagesG else "X", "scope": scope}

    def bern(p):
        return (rw.random(n) < p).astype(float)

    def bg_p():
        p = rw.uniform(*g["background_binary_p"])
        mult = np.ones(L)
        if rw.random() < g["background_lineage_enriched_frac"]:
            k = int(rw.integers(1, 3))  # literal-ok: 1-2 enriched lineages
            mult[rw.choice(L, size=k, replace=False)] = g["background_lineage_enrichment"]
        return np.minimum(p * mult[c.lin], g["background_lineage_enriched_cap"])

    def exact(rows, k):
        v = np.zeros(n)
        v[rw.choice(rows, size=int(k), replace=False)] = 1.0
        return v

    R = c.R
    # --- HOT -------------------------------------------------------------------------------------------
    hot = fam["hot"]
    for x in R["hot_own"] + R["hot_path_x"]:
        add(f"{hot}:{E_[x]}", bern(_u(rw, g["plant_binary_p"])), hot, x)
    for x in R["restr_x"]:
        add(f"{hot}:{E_[x]}", bern(_u(rw, e["restricted_feature_p"])), hot, x)
    for x in R["blind_x"]:
        add(f"{hot}:{E_[x]}", bern(_u(rw, e["blind_p"])), hot, x)
    for x in R["wrong_x"]:
        add(f"{hot}:{E_[x]}", bern(e["wrong_sign_p"]), hot, x)
    for x in R["six_x"]:
        lins = rw.choice(c.big_lineages, size=int(e["six_carrier_n"]), replace=False)
        rows = [int(rw.choice(np.flatnonzero(c.lin == l))) for l in lins]
        v = np.zeros(n)
        v[rows] = 1.0
        add(f"{hot}:{E_[x]}", v, hot, x)
    pin, pout = e["culture_trap_p"]
    c.trap_feature = []
    for i, x in enumerate(R["cult_x"]):
        key = e["culture_trap_families"][i]
        var = c.culture_var[e["culture_vars"][i]]
        name = f"{fam[key]}:{E_[x]}"
        add(name, bern(np.where(var, pin, pout)), fam[key], x,
            stage_family=c.cfg["stage_family_of_sil"] if key == "sil" else None)
        c.trap_feature.append(name)
    c.fold_k = []
    nc_ = int(round(e["fold_artefact_carrier_frac"] * n))
    nin = int(round(e["fold_artefact_in_fold_frac"] * nc_))
    for i, x in enumerate(R["fold_x"]):
        k = i % c.K
        v = exact(np.flatnonzero(c.fold == k), nin) + exact(np.flatnonzero(c.fold != k), nc_ - nin)
        add(f"{hot}:{E_[x]}", v, hot, x)
        c.fold_k.append(k)
    n_hot_bg = int(g["n_hot"]) - sum(1 for m in M.values() if m["family"] == hot)
    for x in rw.choice(c.bg_feature_pool, size=n_hot_bg, replace=False):
        add(f"{hot}:{E_[x]}", bern(bg_p()), hot, int(x))
    # --- LOFBI -----------------------------------------------------------------------------------------
    lof = fam["lof"]
    for x in R["lof_A"] + R["proxy_A"]:
        add(f"{lof}:{E_[x]}", bern(_u(rw, g["plant_binary_p"])), lof, x)
    pm, ps = e["msi_passenger_p"]
    for x in R["msi_pass"]:
        add(f"{lof}:{E_[x]}", bern(np.where(c.msi_true, pm, ps)), lof, x)
    n_amp = 3 * len([l for l in c.loci if l["key"] == "amp"])  # literal-ok: 3 genes per locus
    n_del = 3 * len([l for l in c.loci if l["key"] == "del"])  # literal-ok: 3 genes per locus
    n_prog = 2 * int(g["n_prog_base"])
    n_fixed = 5 + 2 + 2  # literal-ok: SIG 5, COMP 2, PSCORE 2 named features
    n_lof_total = c.n_shared - (int(g["n_hot"]) + int(g["n_sil"]) + n_amp + n_del + n_prog + n_fixed)
    n_lof_bg = n_lof_total - sum(1 for m in M.values() if m["family"] == lof)
    if n_lof_bg < 0:
        raise ValueError("n_features too small for the planted features")
    for x in rw.choice(c.bg_feature_pool, size=n_lof_bg, replace=False):
        add(f"{lof}:{E_[x]}", bern(bg_p()), lof, int(x))
    # --- SIL (mechanism SIL-other -> stage X) -----------------------------------------------------------
    sil, silsf = fam["sil"], c.cfg["stage_family_of_sil"]
    for x in R["sil_A"]:
        add(f"{sil}:{E_[x]}", bern(_u(rw, g["plant_binary_p"])), sil, x, stage_family=silsf)
    c.ldum_lineage = []
    for x in R["ldum_x"]:
        l = int(rw.choice(c.big_lineages))
        add(f"{sil}:{E_[x]}", (c.lin == l).astype(float), sil, x, stage_family=silsf)
        c.ldum_lineage.append(l)
    have = {M[k]["gene"] for k in M if M[k]["family"] == sil}
    sil_pool = [i for i in c.bg_feature_pool if E_[i] not in have]
    sil_genes = list(rw.choice(sil_pool, size=int(g["n_sil"]) - sum(1 for m in M.values() if m["family"] == sil),
                               replace=False))
    c.lineage_proxy_features, c.lineage_proxy_lineage = [], []
    for j, x in enumerate(sil_genes):
        x = int(x)
        if j < int(g["n_sil_lineage_proxy"]):
            l = int(rw.choice(c.big_lineages))
            p = np.where(c.lin == l, _u(rw, g["sil_lineage_proxy_p_in"]), g["sil_lineage_proxy_p_out"])
            c.lineage_proxy_features.append(f"{sil}:{E_[x]}")
            c.lineage_proxy_lineage.append(l)
        else:
            p = _u(rw, g["sil_regular_p"])
        add(f"{sil}:{E_[x]}", bern(p), sil, x, stage_family=silsf)
    # --- AMP / DEL (CN-derived; NaN on CN_absent lines) --------------------------------------------------
    for loc in c.loci:
        key = loc["key"]
        for x in (loc["left"], loc["core"], loc["right"]):
            v = c.call[key][:, x].astype(float)
            v[c.cn_absent] = np.nan
            add(f"{fam[key]}:{E_[x]}", v, fam[key], x)
    # --- SIG / COMP / PSCORE / PROG -------------------------------------------------------------------
    sig, comp = fam["sig"], fam["comp"]
    msi = np.where(c.msi_unknown, np.nan, c.msi_true.astype(float))
    add(f"{sig}:{nf['msi']}", msi, sig)
    add(f"{sig}:{nf['dmg']}", c.dmg, sig, d_type="continuous")
    wgd = bern(rw.uniform(*g["background_binary_p"]) + 0.2)  # literal-ok: WGD prevalence ~0.2-0.3
    add(f"{sig}:{nf['wgd']}", wgd, sig)
    add(f"{sig}:{nf['ploidy']}", 2.1 + 1.6 * wgd + rw.normal(0, 0.25, n), sig, d_type="continuous")  # literal-ok: ploidy model
    add(f"{sig}:{nf['cin']}", 0.6 * c.tp53_true + rw.normal(0, 1, n), sig, d_type="continuous")  # literal-ok: CIN model
    add(f"{comp}:{nf['tp53']}", np.where(c.tp53_unknown, np.nan, c.tp53_true.astype(float)), comp)
    add(f"{comp}:{nf['rb']}", bern(0.07), comp)  # literal-ok: RB_ALT prevalence
    ps_ = fam["pscore"]
    lin_load = rw.normal(0, 1, L)
    add(f"{ps_}:{nf['isg']}", lin_load[c.lin] + rw.normal(0, 1, n), ps_, d_type="continuous")
    k_lof = F[f"{lof}:{E_[R['proxy_A'][-1]]}"]
    add(f"{ps_}:{nf['nrf2']}", rw.normal(0, 1, n) + e["nrf2_score_shift"] * k_lof, ps_, d_type="continuous")
    prog, nb = fam["prog"], int(g["n_prog_base"])
    r = float(g["prog_partner_r"])
    c.prog_nolineage = []
    base = []
    for k in range(nb):
        lw = rw.normal(0, 1, L) if rw.random() < g["prog_lineage_frac"] else np.zeros(L)
        if not lw.any():
            c.prog_nolineage.append(k)
        z = lw[c.lin] + rw.normal(0, 1, n)
        base.append((z - z.mean()) / z.std())
    for k in range(nb):
        add(f"{prog}:{k + 1}", base[k], prog, d_type="continuous")
    for k in range(nb):
        z = r * base[k] + np.sqrt(1 - r * r) * rw.normal(0, 1, n)
        add(f"{prog}:{nb + k + 1}", (z - z.mean()) / z.std(), prog, d_type="continuous")
    n_sh = len(F)
    if n_sh != c.n_shared:
        raise AssertionError(f"shared dictionary has {n_sh} features, expected {c.n_shared}")
    # --- expression ------------------------------------------------------------------------------------
    ex = g["expression"]
    mu_e = rw.uniform(*ex["mean"], c.G)
    lin_e = rw.normal(0, ex["lineage_sd"], (c.G, L))
    Ex = mu_e + lin_e[:, c.lin].T + ex["cn_slope"] * c.log2 + rw.normal(0, ex["sd"], (n, c.G))
    proxyA = set(R["proxy_A"])
    for name, m in M.items():
        if m["family"] == lof and m["gene"]:
            x = int(np.flatnonzero(E_ == m["gene"])[0])
            Ex[F[name] == 1, x] += ex["proxy_lof_shift"] if x in proxyA else ex["lof_shift"]
    for name, m in M.items():
        if m["family"] == sil and m["gene"]:
            x = int(np.flatnonzero(E_ == m["gene"])[0])
            on = F[name] == 1
            Ex[on, x] = rw.uniform(0, ex["sil_max"], int(on.sum()))
    c.unexpressed = []
    for loc in c.loci_by_kind.get("cutting_artefact", []):
        c.unexpressed += [loc["left"], loc["right"]]
    for loc in c.loci_by_kind.get("own_locus_cn", []):
        c.unexpressed.append(loc["left"])
    for x in c.unexpressed:
        Ex[:, x] = rw.uniform(0, ex["unexpressed_max"], n)
    c.expr = np.clip(Ex, 0, None)
    # --- gene-specific features ------------------------------------------------------------------------
    thr = c.dp["dictionary_s3"]["expressed_log2tpm"]
    gf = []
    for x in range(c.G):
        v = np.where(c.expr[:, x] >= thr, c.expr[:, x], 0.0)
        add(f"{fam['expdose']}:{E_[x]}", v - v.mean(), fam["expdose"], x, d_type="continuous", scope="gene")
        gf.append((E_[x], f"{fam['expdose']}:{E_[x]}", "own"))
    for x in range(c.G):
        v = c.call["loss"][:, x].astype(float)
        v[c.cn_absent] = np.nan
        add(f"{fam['loss']}:{E_[x]}", v, fam["loss"], x, scope="gene")
        gf.append((E_[x], f"{fam['loss']}:{E_[x]}", "own"))
        amp_name = f"{fam['amp']}:{E_[x]}"
        if amp_name in F:
            gf.append((E_[x], amp_name, "own"))
    par_union = sorted({p for ps in c.paralogs for p in ps})
    for x in par_union:
        add(f"{fam['exppar']}:{E_[x]}", c.expr[:, x], fam["exppar"], x, d_type="continuous", scope="gene")
    for b in range(c.G):
        for a in c.paralogs[b]:
            gf.append((E_[b], f"{fam['loss']}:{E_[a]}", "paralog"))
            gf.append((E_[b], f"{fam['exppar']}:{E_[a]}", "paralog"))
    c.F, c.M = F, M
    c.gene_features = pd.DataFrame(gf, columns=["gene", "feature", "relation"]).drop_duplicates(["gene", "feature"])


def _effects(c) -> None:
    e, rw, n, G, L, on = c.g["effects"], c.rw, c.n, c.G, c.L, c.on
    F, fam, nf, E_, R = c.F, c.fam, c.nf, c.entrez, c.R
    mu = rw.uniform(*e["gene_mean"], G)
    lin_eff = rw.normal(0, e["lineage_sd"], (G, L))
    nfac = int(e["latent_factors"])
    Lat = rw.normal(0, 1, (n, nfac))
    W = rw.normal(0, e["latent_loading_sd"], (nfac, G))
    planted = set()
    for k, v in R.items():
        if k.endswith("_dep") or k.endswith("_B") or k in ("hot_own",):
            planted |= set(v)
    for loc in c.loci:
        if loc["kind"] in ("own_locus_cn", "self_null", "cutting_artefact"):
            planted |= {loc["core"]} if loc["kind"] != "cutting_artefact" else {loc["left"], loc["core"], loc["right"]}
    u = rw.random(G)
    has_lat = np.where([i in planted for i in range(G)], u < e["latent_planted_frac"], u < e["latent_gene_frac"])
    for k in ("fold_dep", "blind_dep", "six_dep"):
        has_lat[R[k]] = True
    has_lat[R["low_tier"]] = True
    has_lat[R["lin_only"]] = False
    W = W * has_lat
    c.has_latent, c.Lat = has_lat, Lat
    D = np.zeros((n, G))      # molecular/culture determinants (truth only)
    Bl = np.zeros((n, G))     # predictor-blind determinants (truth only; absent from predictor inputs)
    c.beta = {}

    def setmu(x, rng_):
        mu[x] = rw.uniform(*rng_) if isinstance(rng_, list) else rng_

    def plant(x, vec, beta, blind=False):
        (Bl if blind else D)[:, x] += beta * vec
        c.beta[(x, len(c.beta))] = beta

    def col(name):
        return np.nan_to_num(F[name], nan=0.0)

    hot, lof, sil = fam["hot"], fam["lof"], fam["sil"]
    c.pl_rows = []

    def row(kind, sub, x, feature, sign, beta, stage=None, **kw):
        c.pl_rows.append({"kind": kind, "subkind": sub, "gene": E_[x], "feature": feature, "sign": sign,
                          "beta": float(beta), **kw})

    # pan-lineage --------------------------------------------------------------------------------------
    for x in R["hot_own"]:
        setmu(x, e["plant_gene_mean"])
        b = -_u(rw, e["pan_beta"])
        plant(x, col(f"{hot}:{E_[x]}"), b)
        row("pan_lineage", "hot_own", x, f"{hot}:{E_[x]}", "S", b)
    for x, y in zip(R["hot_path_dep"], R["hot_path_x"]):
        setmu(x, e["plant_gene_mean"])
        b = -_u(rw, e["pan_beta"])
        plant(x, col(f"{hot}:{E_[y]}"), b)
        row("pan_lineage", "hot_pathway", x, f"{hot}:{E_[y]}", "S", b)
    for x, y in zip(R["lof_B"], R["lof_A"]):
        setmu(x, e["plant_gene_mean"])
        b = -_u(rw, e["pan_beta"])
        plant(x, col(f"{lof}:{E_[y]}"), b)
        row("pan_lineage", "lofbi_paralog", x, f"{lof}:{E_[y]}", "S", b)
    for x, y in zip(R["loss_B"], R["loss_A"]):
        setmu(x, e["plant_gene_mean"])
        b = -_u(rw, e["loss_beta"])
        plant(x, c.call["loss"][:, y].astype(float), b)
        row("pan_lineage", "loss_paralog", x, f"{fam['loss']}:{E_[y]}", "S", b)
    for x, y in zip(R["sil_B"], R["sil_A"]):
        setmu(x, e["plant_gene_mean"])
        b = -_u(rw, e["pan_beta"])
        plant(x, col(f"{sil}:{E_[y]}"), b)
        row("pan_lineage", "sil_partner", x, f"{sil}:{E_[y]}", "S", b)
    for x, loc in zip(R["amp_partner_dep"], c.loci_by_kind.get("pan_amp", [])):
        setmu(x, e["plant_gene_mean"])
        b = -_u(rw, e["pan_beta"])
        plant(x, c.call["amp"][:, loc["core"]].astype(float), b)
        row("pan_lineage", "amp_partner", x, f"{fam['amp']}:{E_[loc['core']]}", "S", b)
    # lineage-restricted ----------------------------------------------------------------------------------
    sets = [[c.lin_names.index(s) for s in st if s in c.lin_names] for st in e["restricted_lineage_sets"]]
    sets = [s for s in sets if s] or [[c.big_lineages[0]]]
    for i, (x, y) in enumerate(zip(R["restr_dep"], R["restr_x"])):
        setmu(x, e["plant_gene_mean"])
        lins = sets[i % len(sets)]
        b = -_u(rw, e["restricted_beta"])
        plant(x, col(f"{hot}:{E_[y]}") * np.isin(c.lin, lins), b)
        row("lineage_restricted", "hot", x, f"{hot}:{E_[y]}", "S", b,
            lineages="|".join(c.lin_names[l] for l in lins))
    # own-locus CN (AMP of own locus) ---------------------------------------------------------------------
    for loc in c.loci_by_kind.get("own_locus_cn", []):
        x = loc["core"]
        setmu(x, e["plant_gene_mean"])
        b = -_u(rw, e["own_locus_beta"])
        plant(x, c.call["amp"][:, x].astype(float), b)
        row("own_locus_cn", "amp_own", x, f"{fam['amp']}:{E_[x]}", "S", b)
    # MSI-like signature -----------------------------------------------------------------------------------
    for x in R["msi_dep"]:
        setmu(x, e["plant_gene_mean"])
        b = -_u(rw, e["msi_beta"])
        plant(x, c.msi_true.astype(float), b)
        row("msi_like", "sig", x, f"{fam['sig']}:{nf['msi']}", "S", b)
    # culture-confounded: dependency on culture; the trap feature is correlated with it -------------------
    for i, x in enumerate(R["cult_dep"]):
        setmu(x, e["plant_gene_mean"])
        b = -_u(rw, e["culture_beta"])
        var = e["culture_vars"][i]
        plant(x, c.culture_var[var].astype(float), b)
        row("culture_confounded", var, x, c.trap_feature[i], "S", b, culture_var=var)
    # co-located (resolvable) and tight (unresolvable) DEL loci: DEL of the core -> dependency of B ---------
    for kind, Bs in (("co_located", R["coloc_B"]), ("co_located_tight", R["tight_B"])):
        for x, loc in zip(Bs, c.loci_by_kind.get(kind, [])):
            setmu(x, e["plant_gene_mean"])
            b = -_u(rw, e["co_located_beta"])
            plant(x, c.call["del"][:, loc["core"]].astype(float), b)
            row(kind, "del", x, f"{fam['del']}:{E_[loc['core']]}", "S", b,
                passenger=f"{fam['del']}:{E_[loc['left']]}")
    # genomic determinant with an expression proxy (stage G must win) --------------------------------------
    for i, (x, y) in enumerate(zip(R["proxy_B"], R["proxy_A"])):
        setmu(x, e["plant_gene_mean"])
        b = -_u(rw, e["proxy_beta"])
        plant(x, col(f"{lof}:{E_[y]}"), b)
        proxy = (f"{fam['pscore']}:{nf['nrf2']}" if i == len(R["proxy_B"]) - 1 else f"{fam['exppar']}:{E_[y]}")
        row("genomic_with_expression_proxy", "lofbi", x, f"{lof}:{E_[y]}", "S", b, expression_proxy=proxy)
    # expression programme (stage X) ------------------------------------------------------------------------
    nb = int(c.g["n_prog_base"])
    progs = c.prog_nolineage or list(range(nb))
    for i, x in enumerate(R["prog_dep"]):
        setmu(x, e["plant_gene_mean"])
        k = progs[i % len(progs)]
        b = -_u(rw, e["prog_beta_per_sd"])
        plant(x, F[f"{fam['prog']}:{k + 1}"], b)
        row("prog_driven", "prog", x, f"{fam['prog']}:{k + 1}", "S", b, partner=f"{fam['prog']}:{nb + k + 1}")
    # TP53 (W): TP53-mutant lines are less dependent ----------------------------------------------------------
    for x in R["p53_dep"]:
        setmu(x, e["w_sign_gene_mean"])
        b = _u(rw, e["p53_beta"])
        plant(x, c.tp53_true.astype(float), b)
        row("p53_w", "comp", x, f"{fam['comp']}:{nf['tp53']}", "W", b)
    # predictor-blind: in the truth/observations, absent from every predictor input -------------------------
    for x, y in zip(R["blind_dep"], R["blind_x"]):
        setmu(x, e["plant_gene_mean"])
        b = -_u(rw, e["blind_beta"])
        plant(x, col(f"{hot}:{E_[y]}"), b, blind=True)
        row("predictor_blind", "hot", x, f"{hot}:{E_[y]}", "S", b)
    # wrong-sign CGT pair: truth is W, the toy CGT table declares S ------------------------------------------
    for x, y in zip(R["wrong_dep"], R["wrong_x"]):
        setmu(x, e["w_sign_gene_mean"])
        b = float(e["wrong_sign_beta"])
        plant(x, col(f"{hot}:{E_[y]}"), b)
        row("wrong_sign_cgt", "hot", x, f"{hot}:{E_[y]}", "W", b)
    # six-carrier strong plant -------------------------------------------------------------------------------
    for x, y in zip(R["six_dep"], R["six_x"]):
        setmu(x, e["plant_gene_mean"])
        b = -float(e["six_carrier_beta"])
        plant(x, col(f"{hot}:{E_[y]}"), b)
        row("six_carrier_strong", "hot", x, f"{hot}:{E_[y]}", "S", b)
    # lineage-dummy feature: SIL feature == lineage indicator, gene dependent in that lineage (both scenarios)
    for x, y, l in zip(R["ldum_dep"], R["ldum_x"], c.ldum_lineage):
        lin_eff[x, l] = -float(e["lineage_dummy_beta"])
        row("lineage_dummy_feature", "sil", x, f"{sil}:{E_[y]}", "S", -float(e["lineage_dummy_beta"]),
            lineages=c.lin_names[l])
    # lineage-only genes (both scenarios): strong lineage effects, no molecular determinant ------------------
    lo, hi = e["lineage_only_n_lineages"]
    prox_l = list(dict.fromkeys(c.lineage_proxy_lineage)) or c.big_lineages
    for i, x in enumerate(R["lin_only"]):
        k = int(rw.integers(lo, hi + 1))
        first = prox_l[i % len(prox_l)]
        rest = [l for l in c.big_lineages if l != first]
        lins = [first] + list(rw.choice(rest, size=k - 1, replace=False))
        for l in lins:
            lin_eff[x, l] = -_u(rw, e["lineage_only_dep"])
        row("lineage_only", "", x, "", "", 0.0, lineages="|".join(c.lin_names[l] for l in lins))
    for x in R["low_tier"]:
        row("low_tier", "", x, "", "", 0.0)
    # cutting artefact (truth only): every gene of the amplicon is depleted in carriers ----------------------
    for loc in c.loci_by_kind.get("cutting_artefact", []):
        for x in (loc["left"], loc["core"], loc["right"]):
            plant(x, c.call["amp"][:, x].astype(float), -float(e["cutting_beta"]))
            row("cutting_artefact", "core" if x == loc["core"] else "neighbour", x, f"{fam['amp']}:{E_[x]}", "S",
                -float(e["cutting_beta"]), unexpressed=x != loc["core"])
    # truth -----------------------------------------------------------------------------------------------
    T = mu + lin_eff[:, c.lin].T + Lat @ W + on * D + on * Bl
    # self-null (both scenarios): carriers of DEL of the gene have no target -> effect ~ 0
    for loc in c.loci_by_kind.get("self_null", []):
        x = loc["core"]
        mu[x] = e["self_null_gene_mean"]
        T[:, x] = mu[x] + lin_eff[x, c.lin] + Lat @ W[:, x]
        car = c.call["del"][:, x]
        T[car, x] = rw.normal(0, e["self_null_carrier_sd"], int(car.sum()))
        row("self_null", "del", x, f"{fam['del']}:{E_[x]}", "W", 0.0)
    c.mu, c.lin_eff, c.T, c.Tblind = mu, lin_eff, T, on * Bl
    # observations ------------------------------------------------------------------------------------------
    nz = c.g["noise"]
    libs = list(c.g["library_groups"])
    lib_i = np.array([libs.index(v) for v in c.lib])
    off = rw.normal(0, nz["library_offset_sd"], (G, len(libs)))
    c.lib_off = off[:, lib_i].T
    sd_obs = rw.uniform(*nz["obs_sd"], G)
    Y = T + c.lib_off + rw.normal(0, 1, (n, G)) * sd_obs
    Y[np.ix_(c.lib == nz["library_nan_group"], np.asarray(R["nan_genes"], int))] = np.nan
    c.Y = Y
    fo = rw.normal(0, nz["fold_offset_sd"], (c.K, G))
    for x, k in zip(R["fold_dep"], c.fold_k):
        fo[k, x] = e["fold_artefact_offset"]
    c.fold_off = fo
    c.eps_pred = rw.normal(0, 1, (n, G))


def _predict(c) -> None:
    """OOF-like predictions + per-fold offsets. 'ridge': 5-fold CV ridge on the noisy observations (DESIGN Â§13.1);
    'smooth': shrunk tanh of the (predictor-visible) truth + noise."""
    pc = c.g["predictor"]
    n, G = c.n, c.G
    if c.predictor == "smooth":
        s = pc["smooth"]
        Tp = c.T - c.Tblind
        a = s["pred_tanh_scale"]
        Pr = (c.mu + s["pred_shrink"] * a * np.tanh((Tp - c.mu) / a) + s["pred_library_share"] * c.lib_off
              + s["pred_sd"] * c.eps_pred)
    elif c.predictor == "ridge":
        from sklearn.linear_model import Ridge, RidgeCV
        blind = {f"{c.fam['hot']}:{c.entrez[y]}" for y in c.R["blind_x"]}
        cols = [k for k in c.F if k not in blind]
        X = np.column_stack([np.nan_to_num(c.F[k], nan=0.0) for k in cols])
        oh = [(c.lin == l).astype(float) for l in range(c.L)]
        oh += [(c.lib == v).astype(float) for v in c.g["library_groups"]]
        X = np.column_stack([X, np.column_stack(oh), c.Lat])
        sdv = X.std(0)
        X = (X[:, sdv > 0] - X[:, sdv > 0].mean(0)) / sdv[sdv > 0]
        Pr = np.zeros((n, G))
        for k in range(c.K):
            tr, te = c.fold != k, c.fold == k
            Yt = c.Y[tr].copy()
            m = np.nanmean(Yt, axis=0)
            Yt = np.where(np.isnan(Yt), m, Yt)
            rc = pc["ridge"]
            if "alphas" in rc:
                mdl = RidgeCV(alphas=np.asarray(rc["alphas"], float), alpha_per_target=True).fit(X[tr], Yt)
            else:
                mdl = Ridge(alpha=float(rc["alpha"])).fit(X[tr], Yt)
            Pr[te] = mdl.predict(X[te])
    else:
        raise ValueError(f"unknown predictor {c.predictor!r}")
    c.Pred = Pr + c.fold_off[c.fold]


def _assemble(c) -> Toy:
    E_, fam, nf, R = c.entrez, c.fam, c.nf, c.R
    ids = c.ids
    genes = pd.Index(E_, name="entrez")
    # gene info ----------------------------------------------------------------------------------------
    tier = np.where(c.has_latent, "high", "medium").astype(object)
    planted_dep = {int(np.flatnonzero(E_ == r["gene"])[0]) for r in c.pl_rows if r["kind"] not in ("low_tier",)}
    for x in planted_dep:
        tier[x] = "high"
    tier[R["low_tier"]] = "low"
    sd_obs = np.nanstd(c.Y, axis=0)
    gi = pd.DataFrame({"symbol": c.symbol, "chrom": c.chrom, "start_mb": c.pos, "tier": tier,
                       "variable": sd_obs >= c.dp["dependency_s2"]["no_contrast_sd_max"], "selective": False,
                       "paralogs": ["|".join(E_[p] for p in ps) for ps in c.paralogs],
                       "locus_AMP": c.locus["amp"], "locus_DEL": c.locus["del"], "locus_LOSS": c.locus["loss"]},
                      index=genes)
    # features + meta -----------------------------------------------------------------------------------
    feats = pd.DataFrame({k: v for k, v in c.F.items()}, index=ids)
    fm = pd.DataFrame.from_dict(c.M, orient="index")
    fm.index.name = "feature"
    key_of = {fam[k]: k for k in CN_KEYS}
    lg = []
    for f, m in fm.iterrows():
        k = key_of.get(m["family"])
        lg.append(c.locus[k][int(np.flatnonzero(E_ == m["gene"])[0])] if k and m["gene"] else "")
    fm["locus_group"] = lg
    fm["n_carriers"] = np.where(fm["d_type"] == "binary", (feats == 1).sum(0).reindex(fm.index), np.nan)
    fm["n_na"] = feats.isna().sum(0).reindex(fm.index).to_numpy()
    cov = pd.DataFrame({
        "lineage": [c.lin_names[l] for l in c.lin], "library_group": c.lib, "CellFormat": c.cellformat,
        "SerumFree": c.serumfree, "medium_family": c.medium,
        "ScreenDoublingTime": np.where(c.dt_unknown, np.nan, c.dt), "DT_unknown": c.dt_unknown,
        "MSI": c.msi_label, "DMG_BURDEN": c.dmg, "TMB": c.tmb, "TP53_status": c.tp53_label,
        "CN_absent": c.cn_absent, "MSI_unknown": c.msi_unknown,
        "PatientID": [f"TOYP{p + 1:04d}" for p in c.patient], "OncotreePrimaryDisease": c.disease,
        "role": "dev"}, index=ids)
    cn = pd.DataFrame(np.where(c.cn_absent[:, None], np.nan, c.log2), index=ids, columns=genes)
    toy = Toy(seed=c.seed, scenario=c.scenario, predictor=c.predictor, genes=gi, features=feats, feature_meta=fm,
              gene_features=c.gene_features.reset_index(drop=True), covariates=cov,
              expression=pd.DataFrame(c.expr, index=ids, columns=genes), cn_log2=cn,
              effects=pd.DataFrame(c.Y, index=ids, columns=genes),
              predictions=pd.DataFrame(c.Pred, index=ids, columns=genes),
              folds=pd.Series(c.fold, index=ids, name="cv_fold"),
              patients=pd.Series(cov["PatientID"].to_numpy(), index=ids, name="PatientID"),
              controls=pd.DataFrame(), cgt=pd.DataFrame(),
              truth=pd.DataFrame(c.T, index=ids, columns=genes),
              params={"families": fam, "named_features": nf, "design_params_sha256": design_params_sha256(),
                      "toy_config_sha256": toy_config_sha256(), "lineage_groups": c.lin_group})
    _answer_key(c, toy)
    return toy


def _answer_key(c, toy: Toy) -> None:
    fam, nf, E_ = c.fam, c.nf, c.entrez
    fm, F = toy.feature_meta, toy.features
    pl = pd.DataFrame(c.pl_rows)
    pl.insert(0, "plant_id", [f"P{i + 1:03d}" for i in range(len(pl))])
    for col in ("lineages", "passenger", "expression_proxy", "partner", "culture_var", "unexpressed"):
        if col not in pl:
            pl[col] = np.nan
    pl["stage"] = [fm.at[f, "stage"] if f in fm.index else "" for f in pl["feature"]]
    pl["active"] = (c.scenario == "truth") | pl["kind"].isin(["lineage_only", "low_tier", "self_null",
                                                              "lineage_dummy_feature"])
    key_of = {fam[k]: k for k in CN_KEYS}

    def mates(f):
        if f not in fm.index or not fm.at[f, "locus_group"]:
            return []
        m = fm[(fm["family"] == fm.at[f, "family"]) & (fm["locus_group"] == fm.at[f, "locus_group"])]
        return [x for x in m.index if x != f]

    proxies, neutral = [], []
    exppar, expdose = fam["exppar"], fam["expdose"]
    for _, r in pl.iterrows():
        f, g = r["feature"], r["gene"]
        px = mates(f) if f else []
        if r["kind"] == "prog_driven":
            px.append(r["partner"])
        nu = []
        if f and f in fm.index:
            x = fm.at[f, "gene"]
            if x:
                nu += [f"{exppar}:{x}"] + ([f"{expdose}:{x}"] if x == g else [])
                for m in mates(f):
                    mx = fm.at[m, "gene"]
                    nu += [f"{exppar}:{mx}"] + ([f"{expdose}:{mx}"] if mx == g else [])
        if r["kind"] == "msi_like":
            nu.append(f"{fam['sig']}:{nf['dmg']}")
        if r["kind"] == "p53_w":
            nu.append(f"{fam['sig']}:{nf['cin']}")
        if r["kind"] == "genomic_with_expression_proxy":
            nu.append(r["expression_proxy"])
        if r["kind"] == "self_null":  # LOSS_g stays in U (CYCLOPS) and contains the DEL_g carriers (effect ~ 0)
            nu.append(f"{fam['loss']}:{g}")
        proxies.append("|".join(dict.fromkeys(px)))
        neutral.append("|".join(dict.fromkeys(nu)))
    pl["proxies_ok"], pl["neutral"] = proxies, neutral
    # pairs table: (gene, feature) -> role ------------------------------------------------------------------
    role_of_kind = {"culture_confounded": "culture_trap", "fold_offset_artefact": "fold_artefact",
                    "self_null": "self_null", "cutting_artefact": "neutral", "lineage_dummy_feature": "lineage_collinear"}
    rows = []
    for _, r in pl.iterrows():
        if not r["feature"]:
            continue
        rows.append((r["gene"], r["feature"], role_of_kind.get(r["kind"], "true"), r["plant_id"]))
        rows += [(r["gene"], p, "proxy_ok", r["plant_id"]) for p in r["proxies_ok"].split("|") if p]
        rows += [(r["gene"], p, "neutral", r["plant_id"]) for p in r["neutral"].split("|") if p]
    for x in c.R["msi_dep"]:
        rows += [(E_[x], f"{fam['lof']}:{E_[y]}", "passenger", "") for y in c.R["msi_pass"]]
    for i, (x, y) in enumerate(zip(c.R["fold_dep"], c.R["fold_x"])):
        rows.append((E_[x], f"{fam['hot']}:{E_[y]}", "fold_artefact", ""))
    prio = {"true": 0, "proxy_ok": 1, "neutral": 2}
    pairs = pd.DataFrame(rows, columns=["gene", "feature", "role", "plant_id"])
    pairs["_p"] = pairs["role"].map(prio).fillna(len(prio))
    pairs = pairs.sort_values("_p", kind="stable").drop_duplicates(["gene", "feature"]).drop(columns="_p")
    # fold artefact table
    fold_tab = pd.DataFrame({"gene": [E_[x] for x in c.R["fold_dep"]],
                             "feature": [f"{fam['hot']}:{E_[y]}" for y in c.R["fold_x"]],
                             "fold": c.fold_k, "offset": [c.fold_off[k, x] for x, k in zip(c.R["fold_dep"], c.fold_k)]})
    # co-located table (realised discordance on CN-profiled lines)
    md = c.dp["dictionary_s3"]["locus_resolution"]["min_discordant_carriers"]
    co = []
    prof = ~c.cn_absent
    for kind in ("co_located", "co_located_tight"):
        for x, loc in zip(c.R["coloc_B" if kind == "co_located" else "tight_B"], c.loci_by_kind.get(kind, [])):
            core = c.call["del"][prof, loc["core"]]
            nbs = {s: c.call["del"][prof, loc[s]] for s in ("left", "right")}
            only = {s: int((core & ~v).sum()) for s, v in nbs.items()}
            pas_only = int((nbs["left"] & ~core).sum())
            jac = float((core & nbs["left"]).sum() / max(1, (core | nbs["left"]).sum()))
            grp = c.locus["del"][loc["core"]]
            same = grp != "" and grp == c.locus["del"][loc["left"]]
            co.append({"kind": kind, "gene": E_[x], "true_feature": f"{fam['del']}:{E_[loc['core']]}",
                       "passenger": f"{fam['del']}:{E_[loc['left']]}", "locus_group": grp, "same_group": same,
                       "jaccard_true_passenger": jac, "n_true_only_vs_left": only["left"],
                       "n_true_only_vs_right": only["right"], "n_passenger_only": pas_only,
                       "expected_true_level": "gene" if min(only.values()) >= md else "locus",
                       "expected_passenger_level": "locus"})
    co = pd.DataFrame(co)
    sn_rows = []
    for loc in c.loci_by_kind.get("self_null", []):
        gnm = E_[loc["core"]]
        sn_rows += [(gnm, f) for f in self_null_features(toy, gnm)]
    sn = pd.DataFrame(sn_rows, columns=["gene", "feature"])
    # controls (DESIGN Â§9.1 subset) -------------------------------------------------------------------------
    fam_of = {"hot_own": "ADD-HOT", "lofbi_paralog": "PAR", "loss_paralog": "PAR"}
    ctl = []
    for _, r in pl.iterrows():
        cf = None
        if r["kind"] == "pan_lineage" and r["subkind"] in fam_of:
            cf = fam_of[r["subkind"]]
        elif r["kind"] == "own_locus_cn":
            cf = "ADD-CN"
        elif r["kind"] == "msi_like":
            cf = "SIG"
        elif r["kind"] == "p53_w":
            cf = "P53"
        elif r["kind"] == "lineage_only":
            cf = "LINEAGE-CHECK"
        if cf:
            ctl.append({"family": cf, "dependency_entrez": r["gene"], "feature_primary": r["feature"],
                        "feature_secondary": "", "expected_sign": r["sign"] or "",
                        "scope": "class_check" if cf == "LINEAGE-CHECK" else "gate", "plant_id": r["plant_id"]})
    rare = fm[(fm["family"] == fam["hot"]) & (fm["n_carriers"] < c.dp["dictionary_s3"]["min_carriers"])
              & (fm["n_carriers"] > 0)]
    if len(rare):
        f = rare.index[0]
        ctl.append({"family": "ADD-HOT", "dependency_entrez": fm.at[f, "gene"], "feature_primary": f,
                    "feature_secondary": "", "expected_sign": "S", "scope": "untestable", "plant_id": ""})
    ctl = pd.DataFrame(ctl)
    ctl.insert(0, "control_id", [f"TOYC{i + 1:02d}" for i in range(len(ctl))])
    nmin = c.dp["dictionary_s3"]["min_carriers"]
    ctl["n_carriers_D"] = [int((F[f] == 1).sum()) if f in F else 0 for f in ctl["feature_primary"]]
    nonnan = toy.effects.notna().sum(0)
    ctl["gate_eligible"] = [(s == "gate") and (nonnan[g] >= 100) and (k >= nmin)
                            for s, g, k in zip(ctl["scope"], ctl["dependency_entrez"], ctl["n_carriers_D"])]
    toy.controls = ctl
    # toy CGT --------------------------------------------------------------------------------------------
    cg = []

    def cgt_row(tid, r, sign, part, tier, truth):
        cg.append({"test_id": tid, "feature": r["feature"], "feature_D": r["feature"], "target_set": r["gene"],
                   "comparator": "pan-cancer + lineage covariate", "indication_codes": "all", "sign": sign,
                   "partition": part, "ko_class": "concordant", "tier": tier, "true_sign": truth})

    hot_own = pl[(pl["kind"] == "pan_lineage") & (pl["subkind"] == "hot_own")]
    for i, (_, r) in enumerate(hot_own.head(4).iterrows()):  # literal-ok: four sensitivity tests
        cgt_row(f"TOYCGT{i + 1:02d}", r, "S", "sensitivity", "T1", "S")
    for kind, sign, part, tier in (("msi_like", "S", "sensitivity", "T2"), ("p53_w", "W", "bypass-resistance", "T2"),
                                   ("wrong_sign_cgt", "S", "sensitivity", "T1")):
        r = pl[pl["kind"] == kind].iloc[0]
        cgt_row(f"TOYCGT{len(cg) + 1:02d}", r, sign, part, tier, r["sign"])
    bgf = fm[(fm["family"] == fam["hot"]) & (fm["n_carriers"] >= nmin)
             & fm["gene"].isin([E_[x] for x in c.background])]
    if len(bgf):
        f = bgf.index[0]
        cgt_row(f"TOYCGT{len(cg) + 1:02d}", {"feature": f, "gene": fm.at[f, "gene"]}, "N", "expected-null", "T2", "N")
    cg = pd.DataFrame(cg)
    toy.cgt = cg.drop(columns="true_sign")
    toy.answer_key = {
        "scenario": c.scenario, "seed": c.seed, "predictor": c.predictor,
        "design_params_sha256": design_params_sha256(), "toy_config_sha256": toy_config_sha256(),
        "plants": pl, "pairs": pairs.reset_index(drop=True), "colocated": co, "fold_artefact": fold_tab,
        "self_null": sn, "cgt_truth": cg[["test_id", "sign", "true_sign"]],
        "culture": pl[pl["kind"] == "culture_confounded"][["gene", "feature", "culture_var"]].reset_index(drop=True),
        "lineage_proxy_features": list(c.lineage_proxy_features),
        "msi_passengers": [f"{fam['lof']}:{E_[y]}" for y in c.R["msi_pass"]],
        "unexpressed_genes": [E_[x] for x in c.unexpressed],
        "bounds": load_toy_config()["pass"],
    }


# ---------------------------------------------------------------------------------------------------------
# Oracle analyses (simple, homoskedastic OLS; used to derive bounds, never a stand-in for the pipeline)
# ---------------------------------------------------------------------------------------------------------


def _dummies(s: pd.Series, prefix: str) -> pd.DataFrame:
    d = pd.get_dummies(s.astype(str), prefix=prefix, dtype=float)
    return d.iloc[:, 1:] if d.shape[1] > 1 else d.iloc[:, :0]


def oracle_design(toy: Toy, design: str = "full"):
    """Return (Z, droppable) with Z the design matrix and droppable {group: [columns]}.
    design: 'lineage_library' (lineage + library + NA flags), 'design' (DESIGN v3 never-dropped blocks: + culture +
    growth; = EN FWL block), 'full' (design + droppable MSI, DMG_BURDEN, TP53; own CN/EXPGATE omitted)."""
    cov = toy.covariates
    dp = P.load_config("design_params")
    lmin = dp["eligibility_s2"]["lineage_min_lines"]
    lin = cov["lineage"].where(cov["lineage"].map(cov["lineage"].value_counts()) >= lmin,
                               dp["eligibility_s2"]["rare_lineage_label"])
    parts = [pd.DataFrame({"intercept": 1.0}, index=cov.index), _dummies(lin, "lin"),
             _dummies(cov["library_group"], "lib"),
             cov[["CN_absent", "MSI_unknown"]].astype(float)]
    droppable = {}
    if design in ("design", "full"):
        dt = cov["ScreenDoublingTime"]
        dtf = dt.fillna(dt.median())
        parts += [_dummies(cov["CellFormat"], "cf"), _dummies(cov["SerumFree"], "sf"),
                  _dummies(cov["medium_family"], "med"),
                  pd.DataFrame({"DT": (dtf - dtf.mean()) / dtf.std(), "DT_unknown": cov["DT_unknown"].astype(float)})]
    if design == "full":
        msi = pd.DataFrame({"MSI_MSI": (cov["MSI"] == "MSI").astype(float)})
        dmg = pd.DataFrame({"DMG_BURDEN": (cov["DMG_BURDEN"] - cov["DMG_BURDEN"].mean()) / cov["DMG_BURDEN"].std()})
        tp = pd.DataFrame({"TP53_mut": (cov["TP53_status"] == "mut").astype(float),
                           "TP53_unknown": (cov["TP53_status"] == "unknown").astype(float)})
        parts += [msi, dmg, tp]
        droppable = {"MSI": list(msi.columns), "DMG_BURDEN": list(dmg.columns), "TP53_status": list(tp.columns)}
    elif design not in ("lineage_library", "design"):
        raise ValueError(design)
    Z = pd.concat(parts, axis=1)
    keep = [col for col in Z.columns if col == "intercept" or Z[col].std() > 0]
    return Z[keep], {k: [x for x in v if x in keep] for k, v in droppable.items()}


def _orth(Z: np.ndarray) -> np.ndarray:
    U, s, _ = np.linalg.svd(Z, full_matrices=False)
    return U[:, s > s[0] * 1e-10]  # literal-ok: numerical rank tolerance


def _r2_on(Q: np.ndarray, X: np.ndarray) -> np.ndarray:
    """RÂ² of each column of X on span(Q) (Q contains the intercept)."""
    Xc = X - X.mean(0)
    tot = (Xc ** 2).sum(0)
    fit = (Q.T @ X)
    ss = (fit ** 2).sum(0) - X.shape[0] * X.mean(0) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(tot > 0, ss / tot, 0.0)


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    if pv.size == 0:
        return q
    o = np.argsort(pv)
    r = pv[o] * pv.size / np.arange(1, pv.size + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty_like(r)
    out[o] = np.minimum(r, 1.0)
    q[ok] = out
    return q


def fold_centred(toy: Toy) -> pd.DataFrame:
    """Independent re-implementation of per-gene within-fold centring (same contract as p1.fold_centre)."""
    return toy.predictions - toy.predictions.groupby(toy.folds).transform("mean")


def _target(toy: Toy, target: str) -> pd.DataFrame:
    if target == "obs":
        return toy.effects
    if target == "pred":
        return fold_centred(toy)
    if target == "pred_raw":
        return toy.predictions
    raise ValueError(target)


def universe_mask(toy: Toy, genes=None) -> pd.DataFrame:
    """features x genes boolean: d in ð’Ÿ_g (shared minus self-null, plus gene-specific); no carrier filter."""
    fm = toy.feature_meta
    genes = list(toy.genes.index if genes is None else genes)
    U = pd.DataFrame(False, index=fm.index, columns=genes)
    U.loc[fm.index[fm["scope"] == "shared"], :] = True
    for g in genes:
        sn = self_null_features(toy, g)
        if sn:
            U.loc[sn, g] = False
    gf = toy.gene_features[toy.gene_features["gene"].isin(genes)]
    for g, sub in gf.groupby("gene"):
        U.loc[sub["feature"].to_numpy(), g] = True
    return U


def oracle_r_within(toy: Toy) -> pd.Series:
    """DESIGN Â§2 r_within: Pearson r between fold-centred prediction and observed effect after centring both within
    lineage x library strata. Lineages < lineage_min_lines are pooled to _rare first; lines in strata with < 5 lines
    are EXCLUDED (pooling them across libraries leaks library offsets into r_within: verified on the toy, where it
    makes ~30% of lineage-only genes 'eligible')."""
    cov = toy.covariates
    el = P.load_config("design_params")["eligibility_s2"]
    lin = cov["lineage"].where(cov["lineage"].map(cov["lineage"].value_counts()) >= el["lineage_min_lines"],
                               el["rare_lineage_label"])
    strat = lin + "|" + cov["library_group"]
    keep = strat.map(strat.value_counts()) >= el["r_within_stratum_min_lines"]
    Pc = fold_centred(toy)[keep]
    Y = toy.effects[keep]
    strat = strat[keep]
    Pc = Pc - Pc.groupby(strat).transform("mean")
    Yc = Y - Y.groupby(strat).transform("mean")
    Pm = Pc.where(Y.notna())
    num = (Pm * Yc).sum()
    den = np.sqrt((Pm ** 2).sum() * (Yc ** 2).sum())
    return (num / den).rename("r_within")


def oracle_eligible(toy: Toy) -> pd.Series:
    dp = P.load_config("design_params")["eligibility_s2"]
    r = oracle_r_within(toy)
    gi = toy.genes
    return ((gi["variable"] | gi["selective"]) & gi["tier"].isin(dp["tiers"]) & (r >= dp["r_within_lin_lib_min"]))


def oracle_s1(toy: Toy, target: str = "obs", design: str = "full", eligible_only: bool = True) -> pd.DataFrame:
    """Vectorised univariate partial-regression test of every (g, d) in U (S1-like oracle).

    Drop rules (DESIGN v3 Â§5.3): only droppable groups (MSI, DMG_BURDEN, TP53) are dropped, when the group is d
    itself or RÂ²(d ~ group) > drop_assoc_r2_gt; design blocks are never dropped; RÂ²(d ~ design blocks) >
    collinear_with_design_r2_gt labels the test ``collinear`` (never a hit). Own CN and EXPGATE are not in the
    oracle model (equivalent to the own-CN drop for own-locus plants). SEs are homoskedastic OLS (oracle only).
    Returns long rows: gene, feature, family, stage, d_type, beta, t, df, p, n_carriers, test_family (F1/F2),
    collinear, dropped, q (BH within test family over non-collinear rows), hit (F1), hit_f2."""
    dp = P.load_config("design_params")
    c5, s1, d3 = dp["confounders_s5"], dp["s1_s5"], dp["dictionary_s3"]
    fam, nf = toy.params["families"], toy.params["named_features"]
    Y = _target(toy, target)
    genes = list(Y.columns)
    if eligible_only:
        el = oracle_eligible(toy)
        genes = [g for g in genes if el[g]]
    Yv = Y[genes].to_numpy(float)
    fm = toy.feature_meta
    X = toy.features[fm.index].to_numpy(float)
    Xi = np.nan_to_num(X, nan=0.0)
    Z, drop_groups = oracle_design(toy, design)
    Zv = Z.to_numpy(float)
    n = Zv.shape[0]
    # per-feature drop sets and collinearity with the design blocks
    Zd, _ = oracle_design(toy, "design")
    r2_design = _r2_on(_orth(Zd.to_numpy(float)), Xi)
    collinear = r2_design > c5["collinear_with_design_r2_gt"]
    ident = {f"{fam['sig']}:{nf['msi']}": "MSI", f"{fam['sig']}:{nf['dmg']}": "DMG_BURDEN",
             f"{fam['comp']}:{nf['tp53']}": "TP53_status"}
    drops = [set() for _ in range(len(fm))]
    for grp, cols in drop_groups.items():
        Qg = _orth(np.column_stack([np.ones(n), Z[cols].to_numpy(float)]))
        r2 = _r2_on(Qg, Xi)
        for j in np.flatnonzero(r2 > c5["drop_assoc_r2_gt"]):
            drops[j].add(grp)
    for j, f in enumerate(fm.index):
        if f in ident and ident[f] in drop_groups:
            drops[j].add(ident[f])
    U = universe_mask(toy, genes).to_numpy()
    binary = (fm["d_type"] == "binary").to_numpy()
    isnanY = np.isnan(Yv)
    pat_keys = {}
    for gi_ in range(len(genes)):
        pat_keys.setdefault(isnanY[:, gi_].tobytes(), []).append(gi_)
    Fn, Gn = X.shape[1], len(genes)
    B = np.full((Fn, Gn), np.nan)
    T = np.full((Fn, Gn), np.nan)
    DF = np.full((Fn, Gn), np.nan)
    NC = np.full((Fn, Gn), np.nan)
    by_drop = {}
    for j, ds in enumerate(drops):
        by_drop.setdefault(frozenset(ds), []).append(j)
    for ds, fidx in by_drop.items():
        keep = [col for col in Z.columns if not any(col in drop_groups[g] for g in ds)]
        Zk = Z[keep].to_numpy(float)
        fidx = np.asarray(fidx)
        for key, gidx in pat_keys.items():
            rows = ~np.frombuffer(key, dtype=bool)
            Q = _orth(Zk[rows])
            Yr = Yv[np.ix_(rows, gidx)]
            Yr = Yr - Q @ (Q.T @ Yr)
            Xr = Xi[np.ix_(rows, fidx)]
            Xr = Xr - Q @ (Q.T @ Xr)
            df = rows.sum() - Q.shape[1] - 1
            xx = (Xr ** 2).sum(0)
            XY = Xr.T @ Yr
            yy = (Yr ** 2).sum(0)
            with np.errstate(divide="ignore", invalid="ignore"):
                b = XY / xx[:, None]
                rss = np.maximum(yy[None, :] - XY ** 2 / xx[:, None], 0)
                t = b / np.sqrt(rss / df / xx[:, None])
            bad = xx < 1e-8  # literal-ok: residual variance floor
            b[bad], t[bad] = np.nan, np.nan
            B[np.ix_(fidx, gidx)] = b
            T[np.ix_(fidx, gidx)] = t
            DF[np.ix_(fidx, gidx)] = df
            nc = np.where(binary[fidx], (X[np.ix_(rows, fidx)] == 1).sum(0), np.sum(~np.isnan(X[np.ix_(rows, fidx)]), 0))
            NC[np.ix_(fidx, gidx)] = nc[:, None]
    cont_min = d3.get("min_continuous_nonnan", 20)  # literal-ok: DESIGN Â§3 continuous >= 20 non-NaN lines
    okc = np.where(binary[:, None], NC >= s1["min_carriers"], NC >= cont_min)
    keep = U & okc & (~np.isnan(T) | collinear[:, None])
    fi, gi_ = np.nonzero(keep)
    feat = fm.index.to_numpy()[fi]
    out = pd.DataFrame({"gene": np.asarray(genes, dtype=object)[gi_], "feature": feat,
                        "family": fm["family"].to_numpy()[fi], "stage": fm["stage"].to_numpy()[fi],
                        "d_type": fm["d_type"].to_numpy()[fi], "beta": B[fi, gi_], "t": T[fi, gi_],
                        "df": DF[fi, gi_], "n_carriers": NC[fi, gi_], "collinear": collinear[fi],
                        "dropped": ["|".join(sorted(drops[j])) for j in fi]})
    out["p"] = 2 * stats.t.sf(np.abs(out["t"].to_numpy()), out["df"].to_numpy())
    own = (out["family"] == fam["expdose"]) & (out["feature"].str.split(":").str[1] == out["gene"])
    out["test_family"] = np.where(own, "F2", "F1")
    out["q"] = np.nan
    for tf in ("F1", "F2"):
        m = (out["test_family"] == tf) & ~out["collinear"]
        out.loc[m, "q"] = bh(out.loc[m, "p"].to_numpy())
    sig = (out["q"] < s1["bh_q"]) & ~out["collinear"]
    out["hit"] = sig & (out["test_family"] == "F1")      # S1 lists = F1 BH (DESIGN ?5.6)
    out["hit_f2"] = sig & (out["test_family"] == "F2")   # own-EXPDOSE family, reported separately
    return out


def oracle_pair(toy: Toy, gene: str, feature: str, target: str = "obs", design: str = "full",
                extra: pd.DataFrame | None = None, rows=None) -> dict:
    """Single OLS test of y_g on d + design (+ extra columns), droppable groups dropped by the same rules."""
    Z, drop = oracle_design(toy, design)
    y = _target(toy, target)[gene]
    x = toy.features[feature].fillna(0.0)
    if design == "full":
        c5 = P.load_config("design_params")["confounders_s5"]
        for grp, cols in drop.items():
            Q = _orth(np.column_stack([np.ones(len(Z)), Z[cols].to_numpy(float)]))
            if _r2_on(Q, x.to_numpy()[:, None])[0] > c5["drop_assoc_r2_gt"]:
                Z = Z.drop(columns=cols)
    if extra is not None:
        Z = pd.concat([Z, extra], axis=1)
    ok = y.notna() if rows is None else (y.notna() & rows)
    Zv, xv, yv = Z[ok].to_numpy(float), x[ok].to_numpy(float), y[ok].to_numpy(float)
    Q = _orth(Zv)
    xr, yr = xv - Q @ (Q.T @ xv), yv - Q @ (Q.T @ yv)
    df = len(yv) - Q.shape[1] - 1
    b = xr @ yr / (xr @ xr)
    se = np.sqrt(((yr - b * xr) ** 2).sum() / df / (xr @ xr))
    return {"beta": float(b), "se": float(se), "t": float(b / se), "df": int(df),
            "p": float(2 * stats.t.sf(abs(b / se), df))}


def oracle_discordant(toy: Toy, gene: str, feature: str, other: str) -> dict:
    """DESIGN Â§8.4 resolution oracle: two-sample test of y_g on d among CN-profiled lines discordant for d and other."""
    y = toy.effects[gene]
    a, b = toy.features[feature], toy.features[other]
    m = a.notna() & b.notna() & y.notna() & (a != b)
    yd, ad = y[m], a[m] == 1
    n1 = int(ad.sum())
    if n1 < 2 or (~ad).sum() < 2:
        return {"n_discordant_carriers": n1, "beta": np.nan, "p": np.nan}
    r = stats.ttest_ind(yd[ad], yd[~ad])
    return {"n_discordant_carriers": n1, "beta": float(yd[ad].mean() - yd[~ad].mean()), "p": float(r.pvalue)}


def oracle_resolution_level(toy: Toy, gene: str, feature: str, sign: str) -> str:
    """'gene' iff, against every locus-group neighbour, >= min discordant carriers and same sign at p < p_max."""
    rs = P.load_config("design_params")["dictionary_s3"]["locus_resolution"]
    fm = toy.feature_meta
    grp = fm.at[feature, "locus_group"]
    if not grp:
        return "gene"
    nbs = [f for f in fm.index[(fm["locus_group"] == grp) & (fm["family"] == fm.at[feature, "family"])] if f != feature]
    for nb in nbs:
        r = oracle_discordant(toy, gene, feature, nb)
        same = (r["beta"] < 0) if sign == "S" else (r["beta"] > 0)
        if not (r["n_discordant_carriers"] >= rs["min_discordant_carriers"] and same and r["p"] < rs["p_max"]):
            return "locus"
    return "gene"


def oracle_cross_lineage(toy: Toy, gene: str, feature: str) -> str:
    """Simplified DESIGN Â§8.3 label: per-lineage OLS (library-adjusted) in testable lineages, fixed-effect pooling,
    powered = |beta_pool| >= z_0.95 * SE_L, cross-lineage = pooled p < 0.05 / n_planted-scale (p < 1e-3), >= 2
    powered lineages in >= 2 organ groups, >= 80% same sign, leave-one-lineage-out sign kept."""
    dp = P.load_config("design_params")
    cl = dp["subtypes_s8"]["cross_lineage"]
    groups = toy.params["lineage_groups"]
    y = toy.effects[gene]
    x = toy.features[feature]
    cov = toy.covariates
    res = []
    for lin, sub in cov.groupby("lineage"):
        m = y.index.isin(sub.index) & y.notna().to_numpy() & x.notna().to_numpy()
        xs, ys = x[m], y[m]
        nc = int((xs == 1).sum())
        if m.sum() < dp["eligibility_s2"]["lineage_min_lines"] or nc < 3 or (len(xs) - nc) < 3:  # literal-ok: DESIGN Â§8.3 >= 3/3
            continue
        Zl = _dummies(cov.loc[m, "library_group"], "lib")
        Q = _orth(np.column_stack([np.ones(int(m.sum())), Zl.to_numpy(float)]))
        xr = xs.to_numpy() - Q @ (Q.T @ xs.to_numpy())
        yr = ys.to_numpy() - Q @ (Q.T @ ys.to_numpy())
        df = len(ys) - Q.shape[1] - 1
        if df < 2 or xr @ xr < 1e-8:  # literal-ok: tolerance
            continue
        b = xr @ yr / (xr @ xr)
        se = np.sqrt(((yr - b * xr) ** 2).sum() / df / (xr @ xr))
        res.append((lin, b, se))
    if len(res) < 2:
        return "lineage-unresolvable"
    lins, b, se = zip(*res)
    b, se = np.array(b), np.array(se)
    w = 1 / se ** 2
    bp = (w * b).sum() / w.sum()
    p = 2 * stats.norm.sf(abs(bp) * np.sqrt(w.sum()))
    z = stats.norm.ppf(0.95)  # literal-ok: one-sided 0.05
    powered = np.abs(bp) >= z * se
    if powered.sum() < cl["min_powered_lineages"]:
        return "lineage-unresolvable"
    same = np.sign(b[powered]) == np.sign(bp)
    ngroups = len({groups.get(l, l) for l, pw in zip(lins, powered) if pw})
    lolo = all(np.sign(((w * b).sum() - w[i] * b[i]) / (w.sum() - w[i])) == np.sign(bp) for i in range(len(b)))
    if p < 1e-3 and ngroups >= cl["min_organ_groups"] and same.mean() >= cl["same_sign_frac_min"] and lolo:  # literal-ok: oracle pooled alpha
        return "cross-lineage"
    if (p < 1e-3) and (np.abs(b / se) > stats.norm.ppf(0.975)).any():  # literal-ok: oracle per-lineage alpha
        return "lineage-restricted"
    return "undetermined"


def oracle_amplicon(toy: Toy, gene: str, feature: str) -> str:
    """DESIGN Â§8.5 rule applied with oracle statistics."""
    s8 = P.load_config("design_params")["subtypes_s8"]["cutting"]
    fm = toy.feature_meta
    grp = fm.at[feature, "locus_group"]
    mem = fm.index[(fm["locus_group"] == grp) & (fm["family"] == fm.at[feature, "family"])] if grp else [feature]
    car = toy.features[feature] == 1
    prof = toy.features[feature].notna()
    diffs = []
    for f in mem:
        h = fm.at[f, "gene"]
        if h == gene:
            continue
        unexp = toy.expression[h] < s8["unexpressed_log2tpm"]
        y = toy.effects[h]
        a, b = y[car & unexp & prof], y[~car & unexp & prof]
        if len(a) >= 3 and len(b) >= 3:  # literal-ok: minimal group size
            diffs.append(a.median() - b.median())
    yg, eg = toy.effects[gene][car], toy.expression[gene][car]
    ok = yg.notna()
    r, p = stats.pearsonr(yg[ok], eg[ok]) if ok.sum() > 3 else (0.0, 1.0)  # literal-ok: minimal n
    neg_sig = (r < 0) and (p < 0.05)  # literal-ok: DESIGN Â§8.5 alpha
    if not diffs:
        return "unresolved"
    d = float(np.median(diffs))
    thr = s8["neighbour_effect_delta_max"]
    dep_thr = P.load_config("design_params")["dependency_s2"]["threshold"]
    if d < thr and not neg_sig:
        return "cutting artefact"
    if d >= thr and toy.effects[gene][car].median() < dep_thr:
        return "amplified-driver addiction"
    return "unresolved"


# ---------------------------------------------------------------------------------------------------------
# Scoring helpers (shared by the oracle derivation and the pipeline answer-key tests)
# ---------------------------------------------------------------------------------------------------------


def _norm_hits(hits: pd.DataFrame) -> pd.DataFrame:
    h = hits.copy()
    h["gene"] = h["gene"].astype(str)
    h["feature"] = h["feature"].astype(str)
    if "beta" not in h:
        h["beta"] = np.nan
    return h


def recovered(hits: pd.DataFrame, toy: Toy, kind: str, require_sign: bool = True, accept_proxies: bool = True,
              stage: str | None = None, subkind: str | None = None) -> pd.Series:
    """Per plant of ``kind``: True if (gene, primary or proxy_ok feature) is among ``hits`` with the planted sign
    (S: beta < 0, W: beta > 0). ``stage`` additionally requires hits['stage'] == stage."""
    h = _norm_hits(hits)
    if stage is not None:
        h = h[h["stage"] == stage]
    pl = toy.plants(kind)
    if subkind is not None:
        pl = pl[pl["subkind"] == subkind]
    out = {}
    for _, r in pl.iterrows():
        feats = [r["feature"]] + ([p for p in str(r["proxies_ok"]).split("|") if p] if accept_proxies else [])
        m = h[(h["gene"] == r["gene"]) & h["feature"].isin(feats)]
        if require_sign and r["sign"] in ("S", "W"):
            m = m[(m["beta"] < 0) if r["sign"] == "S" else (m["beta"] > 0)]
        out[r["plant_id"]] = len(m) > 0
    return pd.Series(out, dtype=bool)


def fdp_true(hits: pd.DataFrame, toy: Toy) -> dict:
    """True false-discovery proportion of a hit list against the answer key: neutral pairs (expression proxies of a
    true determinant, cutting-artefact pairs, DMG_BURDEN for MSI genes, ...) are excluded; every other hit that is not
    a true/proxy_ok pair is false (culture traps, MSI passengers, fold artefacts, lineage proxies, background)."""
    h = _norm_hits(hits)[["gene", "feature"]].drop_duplicates()
    roles = toy.answer_key["pairs"].set_index(["gene", "feature"])["role"]
    r = pd.Series([roles.get((g, f), "none") for g, f in zip(h["gene"], h["feature"])], dtype=object)
    counted = ~r.isin(NEUTRAL_ROLES)
    n = int(counted.sum())
    nf_ = int((counted & ~r.isin(TRUE_ROLES)).sum())
    return {"n_hits": n, "n_false": nf_, "fdp": nf_ / max(1, n), "n_neutral": int((~counted).sum())}


def list_hits(out: dict, name: str) -> pd.DataFrame:
    """Hit list ``name`` (one of LISTS) from pipeline outputs."""
    if name == "MD-strict":
        return out["md_strict"]
    arm = name.replace(" selected", "")
    a = out["arms"]
    return a[(a["arm"] == arm) & a["selected"].astype(bool)]


# ---------------------------------------------------------------------------------------------------------
# Bound derivation
# ---------------------------------------------------------------------------------------------------------


def oracle_md_strict(toy: Toy, s1_pred: pd.DataFrame, s1_obs: pd.DataFrame) -> pd.DataFrame:
    """MD-strict analogue: S1-pred hits (stand-in for EN-pred selection) that pass the observed-effect screen with the
    same sign at BH q < md_strict_s5.q_sc_bh within the selected set."""
    qmax = P.load_config("design_params")["md_strict_s5"]["q_sc_bh"]
    sel = s1_pred[s1_pred["hit"]][["gene", "feature", "beta"]].rename(columns={"beta": "beta_sel"})
    m = sel.merge(s1_obs[["gene", "feature", "beta", "p", "collinear"]], on=["gene", "feature"], how="inner")
    m = m[~m["collinear"]].copy()
    m["q_sc"] = bh(m["p"].to_numpy())
    return m[(np.sign(m["beta"]) == np.sign(m["beta_sel"])) & (m["q_sc"] < qmax)].reset_index(drop=True)


def oracle_summary(truth: Toy, null: Toy | None = None) -> dict:
    """Oracle statistics of one (truth, pure-null) toy pair; the inputs of ``derive``."""
    s = {}
    so = oracle_s1(truth, "obs", "full")
    sp = oracle_s1(truth, "pred", "full")
    ho, hp = so[so["hit"]], sp[sp["hit"]]
    el = oracle_eligible(truth)
    pl = truth.plants()
    det = pl[pl["kind"].isin(["pan_lineage", "lineage_restricted", "own_locus_cn", "msi_like", "co_located",
                              "co_located_tight", "genomic_with_expression_proxy", "prog_driven", "p53_w",
                              "wrong_sign_cgt", "six_carrier_strong"])]
    s["planted_eligible_frac"] = float(el.reindex(det["gene"]).mean())
    s["lineage_only_ineligible_frac"] = float((~el.reindex(truth.plants("lineage_only")["gene"])).mean())
    s["pan_recall_s1_obs"] = float(recovered(ho, truth, "pan_lineage").mean())
    s["pan_recall_s1_pred"] = float(recovered(hp, truth, "pan_lineage").mean())
    s["own_locus_recall_s1_obs"] = float(recovered(ho, truth, "own_locus_cn").mean())
    s["msi_recovered_s1_obs"] = int(recovered(ho, truth, "msi_like", accept_proxies=False).sum())
    msi_genes = set(truth.plants("msi_like")["gene"])
    pas = set(truth.answer_key["msi_passengers"])
    s["msi_passenger_genes_s1_obs"] = len({g for g, f in zip(ho["gene"], ho["feature"]) if g in msi_genes and f in pas})
    s["proxy_stage_G_s1_pred"] = int(recovered(hp, truth, "genomic_with_expression_proxy", stage="G").sum())
    s["prog_recovered_s1_pred"] = int(recovered(hp, truth, "prog_driven").sum())
    s["p53_recovered_s1_obs"] = int(recovered(ho, truth, "p53_w").sum())
    s["blind_recall_s1_obs"] = float(recovered(ho, truth, "predictor_blind").mean())
    s["blind_hits_s1_pred"] = int(recovered(hp, truth, "predictor_blind").sum())
    s["six_carrier_s1_obs"] = int(recovered(ho, truth, "six_carrier_strong").sum())
    ws = truth.plants("wrong_sign_cgt").iloc[0]
    s["wrong_sign_hit_w_s1_obs"] = int(((ho["gene"] == ws["gene"]) & (ho["feature"] == ws["feature"]) & (ho["beta"] > 0)).any())
    traps = truth.answer_key["culture"]
    trap_set = set(zip(traps["gene"], traps["feature"]))
    s["culture_trap_hits_s1_obs"] = sum((g, f) in trap_set for g, f in zip(ho["gene"], ho["feature"]))
    ll = oracle_s1(truth, "pred", "lineage_library")
    s["culture_trap_live_pred_no_culture"] = sum((g, f) in trap_set for g, f in zip(ll.loc[ll["hit"], "gene"], ll.loc[ll["hit"], "feature"]))
    lo = set(truth.plants("lineage_only")["gene"]) | set(truth.plants("lineage_dummy_feature")["gene"])
    s["lineage_only_hits_s1_obs"] = int(ho["gene"].isin(lo).sum())
    ld = truth.plants("lineage_dummy_feature").iloc[0]
    s["lineage_dummy_collinear"] = int(so.loc[so["feature"] == ld["feature"], "collinear"].all()
                                       and (so["feature"] == ld["feature"]).any())
    fa = truth.answer_key["fold_artefact"]
    raw = oracle_s1(truth, "pred_raw", "full")
    fa_set = set(zip(fa["gene"], fa["feature"]))
    s["fold_live_raw_bh"] = sum((g, f) in fa_set for g, f in zip(raw.loc[raw["hit"], "gene"], raw.loc[raw["hit"], "feature"]))
    s["fold_hits_centred"] = sum((g, f) in fa_set for g, f in zip(hp["gene"], hp["feature"]))
    # co-location
    co = truth.answer_key["colocated"]
    lv = {"co_located_true_gene": 0, "co_located_passenger_gene": 0, "tight_gene": 0}
    for _, r in co.iterrows():
        tl = oracle_resolution_level(truth, r["gene"], r["true_feature"], "S")
        pl_ = oracle_resolution_level(truth, r["gene"], r["passenger"], "S")
        if r["kind"] == "co_located":
            lv["co_located_true_gene"] += tl == "gene"
            lv["co_located_passenger_gene"] += pl_ == "gene"
        else:
            lv["tight_gene"] += (tl == "gene") + (pl_ == "gene")
    s.update(lv)
    s["co_located_same_group"] = int(co["same_group"].sum())
    s["co_located_min_jaccard"] = float(co["jaccard_true_passenger"].min())
    # cross-lineage
    pan = truth.plants("pan_lineage")
    s["cross_lineage_pan_frac"] = float(np.mean([oracle_cross_lineage(truth, g, f) == "cross-lineage"
                                                 for g, f in zip(pan["gene"], pan["feature"])]))
    rs = truth.plants("lineage_restricted")
    s["cross_lineage_restricted_n"] = sum(oracle_cross_lineage(truth, g, f) == "cross-lineage"
                                          for g, f in zip(rs["gene"], rs["feature"]))
    # cutting / addiction
    cut = truth.plants("cutting_artefact")
    cut = cut[cut["subkind"] == "core"]
    s["cutting_core_labelled"] = sum(oracle_amplicon(truth, g, f) == "cutting artefact" for g, f in zip(cut["gene"], cut["feature"]))
    own = truth.plants("own_locus_cn")
    s["own_locus_addiction"] = sum(oracle_amplicon(truth, g, f) == "amplified-driver addiction" for g, f in zip(own["gene"], own["feature"]))
    # self-null (NC3)
    sn = truth.answer_key["self_null"]
    thr = P.load_config("design_params")["dependency_s2"]["threshold"]
    dep = [truth.effects.loc[truth.features[f] == 1, g].median() < thr for g, f in zip(sn["gene"], sn["feature"])]
    s["nc3_frac_dependent"] = float(np.mean(dep)) if dep else 0.0
    # gate (S1-obs oracle)
    ctl = truth.controls[truth.controls["gate_eligible"]]
    rec = []
    for _, r in ctl.iterrows():
        m = ho[(ho["gene"] == r["dependency_entrez"]) & (ho["feature"] == r["feature_primary"])]
        m = m[(m["beta"] < 0) if r["expected_sign"] == "S" else (m["beta"] > 0)]
        rec.append(len(m) > 0)
    s["gate_recall_s1_obs"] = float(np.mean(rec)) if rec else 0.0
    # error rates on the planted toy
    s["fdp_s1_obs"] = fdp_true(ho, truth)["fdp"]
    s["fdp_s1_pred"] = fdp_true(hp, truth)["fdp"]
    s["n_hits_s1_obs"] = int(len(ho))
    s["n_hits_s1_pred"] = int(len(hp))
    # oracle MD-strict: S1-pred hits as the "selected" set, screened on observed effects (same sign, BH within)
    md = oracle_md_strict(truth, sp, so)
    s["fdp_md_strict_oracle"] = fdp_true(md, truth)["fdp"]
    s["n_md_strict_oracle"] = int(len(md))
    s["pan_recall_md_strict_oracle"] = float(recovered(md, truth, "pan_lineage").mean())
    s["culture_trap_md_strict_oracle"] = sum((g, f) in trap_set for g, f in zip(md["gene"], md["feature"]))
    if null is not None:
        no = oracle_s1(null, "obs", "full")
        npr = oracle_s1(null, "pred", "full")
        # pure-null hits exclude neutral pairs (the self-null gene's own LOSS/EXPDOSE are real, kept on purpose)
        s["null_hits_s1_obs"] = fdp_true(no[no["hit"]], null)["n_hits"]
        s["null_hits_s1_pred"] = fdp_true(npr[npr["hit"]], null)["n_hits"]
        f1 = no[(no["test_family"] == "F1") & ~no["collinear"]]
        s["null_frac_p_lt_0.05"] = float((f1["p"] < 0.05).mean())  # literal-ok: nominal alpha
        s["null_hits_md_strict_oracle"] = fdp_true(oracle_md_strict(null, npr, no), null)["n_hits"]
    return s


def derive(n_seeds: int = 20, predictors=("ridge", "smooth")) -> dict:  # literal-ok: default seed count
    """Simulate truth toys (toy_seed + i) and pure-null toys (toy_null_seed + i) per predictor; summarise the oracle."""
    sd = P.seeds()
    out = {"_n_seeds": n_seeds}
    for pr in predictors:
        rows = []
        for i in range(n_seeds):
            t = make_toy(sd["toy_seed"] + i, "truth", predictor=pr)
            nl = make_toy(sd["toy_null_seed"] + i, "null", predictor=pr)
            rows.append(oracle_summary(t, nl))
        df = pd.DataFrame(rows)
        summ = {k: {"mean": round(float(df[k].mean()), 4), "p05": round(float(df[k].quantile(0.05)), 4),  # literal-ok: quantiles
                    "p95": round(float(df[k].quantile(0.95)), 4), "min": round(float(df[k].min()), 4),  # literal-ok: quantiles
                    "max": round(float(df[k].max()), 4)} for k in df.columns}  # literal-ok: rounding
        for lst in ("s1_obs", "s1_pred", "md_strict_oracle"):
            summ[f"null_p_any_hit_{lst}"] = round(float((df[f"null_hits_{lst}"] > 0).mean()), 4)  # literal-ok: rounding
        out[pr] = summ
    return out


def _floor_step(x: float, step: float) -> float:
    return round(float(np.floor(x / step + 1e-9) * step), 4)  # literal-ok: float tolerance, rounding


def propose_bounds(deriv: dict) -> dict:
    """Pass bounds from the oracle derivation (rules documented in docs/toy_answer_key.md):
    fraction 'min' bounds = min(DESIGN value, floor_0.05(oracle p05 - headroom)); count 'min' bounds over n plants =
    min(DESIGN value, floor(n * (p05/n - headroom))); 'max' bounds = max(DESIGN value, ceil(oracle p95)); the error-rate
    bounds are the DESIGN values (the oracle only checks feasibility). Headroom = 0.10 (the pipeline's EN and CR2 tests
    are expected to be somewhat less powerful than the homoskedastic oracle)."""
    dp = P.load_config("design_params")
    dpass = dp["toy_s13"]["pass"]
    gate = dp["gate_s9"]
    cfg = load_toy_config()
    pl = cfg["generator"]["plants"]
    H, step = 0.1, 0.05  # literal-ok: headroom and rounding step of the bound rule
    sm = deriv["smooth"]

    def p05(k):
        return sm[k]["p05"]

    def p95(k):
        return sm[k]["p95"]

    def fmin(k, design=None):
        v = max(0.0, _floor_step(p05(k) - H, step))
        return v if design is None else min(design, v)

    def cmin(k, n, design=None):
        v = int(max(0, np.floor(n * (p05(k) / n - H) + 1e-9)))  # literal-ok: float tolerance
        return v if design is None else min(design, v)

    def cmax(k, design=None):
        v = int(np.ceil(p95(k) - 1e-9))  # literal-ok: float tolerance
        return v if design is None else max(design, v)

    b = {
        "_rule": "see propose_bounds docstring; derived from derivation.smooth over derivation._n_seeds seeds",
        "_lists_hard_on_smooth_toy": list(LISTS),
        "_lists_hard_on_ridge_toy": ["S1-obs", "EN-obs selected", "MD-strict"],
        "_lists_reported_on_ridge_toy": ["S1-pred", "EN-pred selected"],
        "pan_lineage_recall_min": fmin("pan_recall_md_strict_oracle", dpass["pan_lineage_recall_min"]),
        "pan_lineage_recall_s1_obs_min": fmin("pan_recall_s1_obs"),
        "own_locus_cn_recall_min": fmin("own_locus_recall_s1_obs", dpass["own_locus_cn_recall_min"]),
        "msi_like_signature_recovered_min": cmin("msi_recovered_s1_obs", pl["msi_like"],
                                                 dpass["msi_like_signature_recovered_min"]),
        "msi_like_passenger_lof_selected_max": cmax("msi_passenger_genes_s1_obs",
                                                    dpass["msi_like_passenger_lof_selected_max"]),
        "co_located_locus_level_min": dpass["co_located_locus_level_min"],
        "co_located_true_gene_level_min": cmin("co_located_true_gene", pl["co_located"]),
        "co_located_passenger_gene_level_max": cmax("co_located_passenger_gene",
                                                    pl["co_located"] - dpass["co_located_locus_level_min"]),
        "co_located_tight_gene_level_max": cmax("tight_gene"),
        "culture_confounded_in_md_strict_max": cmax("culture_trap_md_strict_oracle",
                                                    dpass["culture_confounded_in_md_strict_max"]),
        "culture_flagged_if_reported_frac_min": 1.0 - 0.2,  # literal-ok: >= 80% of reported trap rows flagged
        "lineage_only_in_md_strict_max": cmax("lineage_only_hits_s1_obs", dpass["lineage_only_in_md_strict_max"]),
        "lineage_only_class3_frac_min": fmin("lineage_only_ineligible_frac", dpass["lineage_only_class3_frac_min"]),
        "lineage_dummy_in_md_strict_max": 0,
        "lineage_dummy_flagged_collinear": True,
        "fold_offset_selected_after_centring_max": cmax("fold_hits_centred",
                                                        dpass["fold_offset_selected_after_centring_max"]),
        "fold_offset_selected_without_centring_min": cmin("fold_live_raw_bh", pl["fold_offset_artefact"],
                                                          dpass["fold_offset_selected_without_centring_min"]),
        "genomic_stage_G_recovered_min": cmin("proxy_stage_G_s1_pred", pl["genomic_with_expression_proxy"],
                                              dpass["genomic_stage_G_recovered_min"]),
        "prog_driven_en_pred_min": cmin("prog_recovered_s1_pred", pl["prog_driven"]),
        "p53_md_strict_min": cmin("p53_recovered_s1_obs", pl["p53_w"]),
        "six_carrier_selected_subsample_frac_min": dpass["six_carrier_selected_subsample_frac_min"],
        "predictor_blind_en_pred_selected_max": cmax("blind_hits_s1_pred"),
        "predictor_blind_s1_obs_recall_min": fmin("blind_recall_s1_obs"),
        "wrong_sign_counts_as_sign_error": dpass["wrong_sign_counts_as_sign_error"],
        "cross_lineage_pan_frac_min": fmin("cross_lineage_pan_frac"),
        "cross_lineage_restricted_max": cmax("cross_lineage_restricted_n", 1),
        "cutting_core_label_min": cmin("cutting_core_labelled", pl["cutting_artefact_loci"]),
        "own_locus_addiction_label_min": cmin("own_locus_addiction", pl["own_locus_cn"]),
        "self_null_in_universe_max": 0,
        "nc3_frac_dependent_max": dp["negative_controls_s7"]["NC3_self_null_dependent_frac_max"],
        "gate_overall_recall_min": gate["overall_recall_min"],
        "gate_family_recall_min": gate["family_recall_min"],
        "planted_eligible_frac_min": fmin("planted_eligible_frac"),
        "planted_fdp_true_max": dpass["planted_fdp_true_max"],
        "fdr_emp_fixed_minus_fdp_true_min": dpass["fdr_emp_fixed_minus_fdp_true_min"],
        "pure_null_bh_lists": list(BH_LISTS),
        "pure_null_bh_lists_p_any_hit_max": dpass["pure_null_bh_lists_p_any_hit_max"],
        "pure_null_B": dpass["pure_null_B"],
    }
    b["_feasibility"] = {
        "smooth_fdp_p95": {k: sm[k]["p95"] for k in ("fdp_s1_obs", "fdp_s1_pred", "fdp_md_strict_oracle")},
        "smooth_null_p_any_hit": {k: sm[f"null_p_any_hit_{k}"] for k in ("s1_obs", "s1_pred", "md_strict_oracle")},
        "ridge_fdp_p05": {k: deriv["ridge"][k]["p05"] for k in ("fdp_s1_obs", "fdp_s1_pred", "fdp_md_strict_oracle")},
        "ridge_null_p_any_hit": {k: deriv["ridge"][f"null_p_any_hit_{k}"] for k in ("s1_obs", "s1_pred", "md_strict_oracle")},
    }
    return b


# ---------------------------------------------------------------------------------------------------------
# Pipeline adapter
# ---------------------------------------------------------------------------------------------------------


def get_pipeline():
    """``biomark.pipeline.run_pipeline`` if the module exists, else None (the answer-key tests then skip)."""
    if importlib.util.find_spec(PIPELINE_MODULE) is None:
        return None
    mod = importlib.import_module(PIPELINE_MODULE)
    return getattr(mod, PIPELINE_ENTRY, None)


def validate_outputs(out: dict) -> list[str]:
    """Problems with a pipeline output dict versus REQUIRED_OUTPUTS (empty list = conforming)."""
    probs = []
    for key, cols in REQUIRED_OUTPUTS.items():
        if key not in out:
            probs.append(f"missing output {key!r}")
            continue
        obj = out[key]
        have = set(obj.keys()) if isinstance(obj, dict) else set(getattr(obj, "columns", []))
        miss = [c for c in cols if c not in have]
        if miss:
            probs.append(f"output {key!r} lacks {miss}")
    return probs


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m biomark.toysim")
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("derive", help="oracle statistics over seeds (bound derivation)")
    d.add_argument("--n-seeds", type=int, default=20)  # literal-ok: CLI default
    d.add_argument("--write", action="store_true", help="store the statistics in configs/toy_answer_key.json:derivation")
    a = ap.parse_args(argv)
    if a.cmd == "derive":
        summ = derive(a.n_seeds)
        print(json.dumps(summ, indent=1))
        if a.write:
            path = P.CONFIGS / CONFIG
            cfg = json.loads(path.read_text(encoding="utf-8"))
            cfg["derivation"] = summ
            cfg["pass"] = propose_bounds(summ)
            path.write_text(json.dumps(cfg, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
