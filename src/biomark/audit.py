"""Confounder audits and downstream labels on the arms' outputs (first implementation, used by the toy adapter
and later by the M5/M6 scripts): confounding delta (§5.5), co-location resolution (§8.4), CN-cutting artefact (§8.5),
within vs across lineage (§8.3), self-null NC3 (§7.3), gene classes (§12.1, pre-R classes only), CGT recovery
(§10.5 hit rule) and the positive-control gate (§9.3). Every (g, d) model is the §5.3 confounder model with the same
drop rules as the S1 scan (stats.pair_test, CR2 + ν).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as sps

from . import arms as A
from . import params as PR
from . import stats as ST


# ------------------------------------------------------------------------------------------------
# One (g, d) confounder model with options
# ------------------------------------------------------------------------------------------------


def pair_fit(data, g, f, extra_rows=None, drop_blocks=(), add_blocks=(), add_cols=None, target="obs") -> dict:
    sp = PR.stats_params()
    gj = data._gpos[g]
    fj = data.fmeta.index.get_loc(f)
    y = data.target(target)[:, gj]
    x = data.X[:, fj]
    rows = data.rows_ok[:, gj] & np.isfinite(x)
    if extra_rows is not None:
        rows = rows & extra_rows
    if rows.sum() < 2:
        return {"beta": np.nan, "se": np.nan, "p": np.nan, "n": int(rows.sum()), "df": np.nan, "nu": np.nan,
                "n_carriers": 0, "sigma_resid": np.nan}
    xr = x[rows]
    S = set(data.drop.ident[fj])
    for nm in data.cov.droppable_names:
        Bk = data.cov.block(nm, rows)
        if Bk.shape[1] and ST.r2_on(ST.orth(Bk), xr)[0] > sp["drop_r2"]:
            S.add(nm)
    expand = set()
    for b in drop_blocks:
        expand |= set(PR.pp()["covariates"]["delta_blocks"].get(b, [b]))
    names = [b for b in data.cov.design_names if b not in expand]
    names += [nm for nm in data.cov.droppable_names if nm not in S and nm not in expand]
    names += list(add_blocks)
    Z = data.cov.Z(rows, names)
    inp = A.scan_inputs(data)
    gc = ST._gene_covs(inp, gj, rows)
    if gc["cols"]:
        codes, _ = ST._gene_cov_codes(gc, xr[:, None], np.array([fj]), g, data.drop, sp)
        keep = ST._code_keep(int(codes[0]), gc)
        cols = [c for c, k, kd in zip(gc["cols"], keep, gc["kind"]) if k and kd not in expand]
        if cols:
            Z = np.hstack([Z, np.column_stack(cols)])
    if add_cols is not None:
        Z = np.hstack([Z, np.asarray(add_cols, float)[rows]])
    r = ST.pair_test(y[rows], xr, Z, data.patients[rows])
    r["rows"] = rows
    return r


def _pairs_of(res: dict, data, controls: pd.DataFrame | None, top_s1: int | None = None) -> pd.DataFrame:
    a = res["arms"]
    parts = []
    md = res.get("md_strict")
    if md is not None and len(md):
        parts.append(md[["gene", "feature"]].assign(source="MD-strict"))
    s1 = a[(a["arm"] == "S1-obs") & a["selected"]]
    if top_s1 is not None:
        s1 = s1.sort_values("rank").head(top_s1)
    parts.append(s1[["gene", "feature"]].assign(source="S1-obs"))
    if controls is not None and len(controls):
        U = res["universe"]
        uk = set(zip(U["gene"], U["feature"]))
        c = [(g, f) for g, f in zip(controls["dependency_entrez"], controls["feature_primary"]) if (g, f) in uk]
        if c:
            parts.append(pd.DataFrame(c, columns=["gene", "feature"]).assign(source="control"))
    p = pd.concat(parts, ignore_index=True)
    return p.groupby(["gene", "feature"], as_index=False)["source"].agg(lambda s: "|".join(sorted(set(s))))


# ------------------------------------------------------------------------------------------------
# §5.5 confounding delta
# ------------------------------------------------------------------------------------------------


def confounding_delta(data, pairs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    blocks = list(PR.pp()["covariates"]["delta_blocks"])
    for g, f in zip(pairs["gene"], pairs["feature"]):
        full = pair_fit(data, g, f)
        for b in blocks:
            m = pair_fit(data, g, f, drop_blocks=(b,))
            dlt = abs(full["beta"] - m["beta"])
            flag = bool(dlt > 0.5 * abs(m["beta"]) and dlt > m["se"])
            rows.append({"gene": g, "feature": f, "block": b, "beta_full": full["beta"], "beta_minus": m["beta"],
                         "se_minus": m["se"], "flagged": flag})
        s = pair_fit(data, g, f, add_blocks=("primary_disease",))
        dlt = abs(full["beta"] - s["beta"])
        rows.append({"gene": g, "feature": f, "block": "primary_disease (sensitivity: subtype-confounded)",
                     "beta_full": full["beta"], "beta_minus": s["beta"], "se_minus": s["se"],
                     "flagged": bool(dlt > 0.5 * abs(s["beta"]) and dlt > s["se"])})
    return pd.DataFrame(rows, columns=["gene", "feature", "block", "beta_full", "beta_minus", "se_minus", "flagged"])


# ------------------------------------------------------------------------------------------------
# §8.4 co-location
# ------------------------------------------------------------------------------------------------


def _group_members(data, f) -> list:
    fm = data.fmeta
    grp = fm.at[f, "locus_group"]
    if grp is None or (isinstance(grp, float) and np.isnan(grp)) or grp == "":
        return []
    m = fm.index[(fm["locus_group"] == grp) & (fm["family"] == fm.at[f, "family"])]
    return [x for x in m if x != f]


def locus_resolution(data, pairs: pd.DataFrame, sign_of: dict) -> pd.DataFrame:
    """Discordance test per neighbour: rows without concordant carriers (d = nb = 1); y ~ model + nb + d; gene-level
    iff for every neighbour >= min discordant carriers, same sign and p < p_max."""
    lr = PR.dp()["dictionary_s3"]["locus_resolution"]
    cnf = set(PR.pp()["covariates"]["cn_derived_families"]) | {"SIL", "SIL-DEL", "SIL-MUT", "SIL-METH", "SIL-other"}
    out = []
    for g, f in zip(pairs["gene"], pairs["feature"]):
        fam = data.fmeta.at[f, "family"]
        nbs = _group_members(data, f)
        if fam not in cnf or not nbs:
            continue
        x = data.X[:, data.fmeta.index.get_loc(f)]
        level, nmin, det = "gene", np.inf, []
        for nb in nbs:
            v = data.X[:, data.fmeta.index.get_loc(nb)]
            conc = (x == 1) & (v == 1)
            ok = np.isfinite(v) & ~conc
            r = pair_fit(data, g, f, extra_rows=ok, add_cols=np.nan_to_num(v)[:, None])
            ndis = int(((x == 1) & (v == 0) & data.rows_ok[:, data._gpos[g]]).sum())
            nmin = min(nmin, ndis)
            s = sign_of.get((g, f), np.sign(r["beta"]))
            good = ndis >= lr["min_discordant_carriers"] and np.sign(r["beta"]) == s and r["p"] < lr["p_max"]
            det.append(f"{nb}:{ndis}:{r['p']:.3g}")
            if not good:
                level = "locus"
        out.append({"gene": g, "feature": f, "locus_group": data.fmeta.at[f, "locus_group"], "level": level,
                    "n_discordant_carriers": int(nmin), "neighbours": ";".join(det)})
    return pd.DataFrame(out, columns=["gene", "feature", "locus_group", "level", "n_discordant_carriers", "neighbours"])


# ------------------------------------------------------------------------------------------------
# §8.5 CN-cutting artefact
# ------------------------------------------------------------------------------------------------


def amplicon_test(data, pairs: pd.DataFrame, effects_of, expression_of) -> pd.DataFrame:
    """AMP determinants whose group contains g. ``effects_of(genes)`` / ``expression_of(genes)`` -> D x genes."""
    cut = PR.dp()["subtypes_s8"]["cutting"]
    dep_thr = PR.dp()["dependency_s2"]["threshold"]
    amin = int(PR.pp()["audit"]["min_group_size_amplicon"])
    nmin = int(PR.pp()["audit"]["min_n_corr"])
    out = []
    for g, f in zip(pairs["gene"], pairs["feature"]):
        if data.fmeta.at[f, "family"] != "AMP":
            continue
        mem = [f] + _group_members(data, f)
        mg = [data.fmeta.at[m, "genes"][0] for m in mem if data.fmeta.at[m, "genes"]]
        if g not in mg:
            continue
        x = data.X[:, data.fmeta.index.get_loc(f)]
        car, prof = x == 1, np.isfinite(x)
        others = [h for h in mg if h != g]
        diffs = []
        if others:
            E, X_ = effects_of(others), expression_of(others)
            for h in others:
                une = X_[h].to_numpy() < cut["unexpressed_log2tpm"]
                yh = E[h].to_numpy()
                a_, b_ = yh[car & une & prof], yh[~car & une & prof]
                a_, b_ = a_[np.isfinite(a_)], b_[np.isfinite(b_)]
                if a_.size >= amin and b_.size >= amin:
                    diffs.append(float(np.median(a_) - np.median(b_)))
        yg = data.target("obs")[:, data._gpos[g]]
        eg = expression_of([g])[g].to_numpy()
        m = car & np.isfinite(yg) & np.isfinite(eg)
        neg = False
        if m.sum() >= nmin and np.ptp(eg[m]) > 0:
            r, p = sps.pearsonr(yg[m], eg[m])
            neg = bool(r < 0 and p / 2 < cut["alpha"])
        dmed = float(np.median(diffs)) if diffs else np.nan
        cmed = float(np.nanmedian(yg[car])) if car.any() else np.nan
        if not diffs:
            lab = "unresolved"
        elif dmed < cut["neighbour_effect_delta_max"] and not neg:
            lab = "cutting artefact"
        elif dmed >= cut["neighbour_effect_delta_max"] and cmed < dep_thr:
            lab = "amplified-driver addiction"
        else:
            lab = "unresolved"
        out.append({"gene": g, "feature": f, "label": lab, "neighbour_delta": dmed, "n_unexpressed_neighbours":
                    len(diffs), "carrier_median_effect": cmed, "expr_corr_neg_sig": neg})
    return pd.DataFrame(out, columns=["gene", "feature", "label", "neighbour_delta", "n_unexpressed_neighbours",
                                      "carrier_median_effect", "expr_corr_neg_sig"])


# ------------------------------------------------------------------------------------------------
# §8.3 within vs across lineage
# ------------------------------------------------------------------------------------------------


def cross_lineage(data, pairs: pd.DataFrame, beta_shrunk: dict, sign_of: dict) -> pd.DataFrame:
    s8 = PR.dp()["subtypes_s8"]["cross_lineage"]
    groups = {L: k for k, v in PR.dp()["subtypes_s8"]["lineage_groups"].items() for L in v}
    lmin = int(PR.dp()["eligibility_s2"]["lineage_min_lines"])
    rmin = int(PR.dp()["replication_s6"]["R_min_carriers"])
    q_thr = float(s8["meta_bh_q"])
    rare = PR.dp()["eligibility_s2"]["rare_lineage_label"]
    lin_raw = data.cov_df.get("lineage_raw", data.cov_df["lineage"]).astype(str).to_numpy()
    recs, per = [], []
    for g, f in zip(pairs["gene"], pairs["feature"]):
        gj, fj = data._gpos[g], data.fmeta.index.get_loc(f)
        x = data.X[:, fj]
        base = data.rows_ok[:, gj] & np.isfinite(x)
        binary = data.fmeta.at[f, "d_type"] == "binary"
        testable, powered, lres = [], [], {}
        lin_car = set()
        full = pair_fit(data, g, f)
        bsh = beta_shrunk.get((g, f), full["beta"])
        for L in sorted(set(lin_raw[base])):
            if L == rare:
                continue
            m = base & (lin_raw == L)
            nL = int(m.sum())
            n1 = int((x[m] == 1).sum()) if binary else nL
            n0 = nL - n1 if binary else nL
            if binary and n1 > 0:
                lin_car.add(L)
            if nL < lmin or (binary and (n1 < rmin or n0 < rmin)) or (not binary and np.nanstd(x[m]) <= 0):
                continue
            r = pair_fit(data, g, f, extra_rows=m)
            if not np.isfinite(r.get("beta", np.nan)):
                continue
            testable.append(L)
            sep = ST.se_r_pred(full["sigma_resid"], binary, n1, n0, r.get("sd_dtilde"), nL)
            pw = float(ST.power_t(bsh, sep, r["nu"], s8["powered_alpha_one_sided"]))
            lres[L] = (r["beta"], r["p"], pw)
            if pw >= s8["powered_lineage_power_min"]:
                powered.append(L)
            per.append({"gene": g, "feature": f, "lineage": L, "beta": r["beta"], "p": r["p"], "power": pw})
        rec = {"gene": g, "feature": f, "n_testable": len(testable), "n_powered": len(powered),
               "lineages_powered": "|".join(powered)}
        if len(testable) >= 1:
            mt = base & np.isin(lin_raw, testable)
            pool = pair_fit(data, g, f, extra_rows=mt)
            rec.update(beta_pooled=pool["beta"], p_pooled=pool["p"])
            lolo = True
            for L in testable:
                if len(testable) < 2:
                    break
                r = pair_fit(data, g, f, extra_rows=mt & (lin_raw != L))
                lolo &= np.sign(r["beta"]) == np.sign(pool["beta"])
            rec["lolo_sign_kept"] = bool(lolo)
            sgn = np.sign(pool["beta"])
            same = [np.sign(lres[L][0]) == sgn for L in powered]
            rec["same_sign_frac"] = float(np.mean(same)) if same else np.nan
            rec["n_organ_groups"] = len({groups.get(L, L) for L in powered})
            rec["p_interaction"] = _interaction_p(data, g, f, mt, lin_raw) if len(testable) >= 2 else np.nan
        rec["single_lineage_carriers"] = len(lin_car) == 1
        recs.append(rec)
    R = pd.DataFrame(recs)
    P5 = pd.DataFrame(per)
    if R.empty:
        return pd.DataFrame(columns=["gene", "feature", "label"])
    R["q_pooled"] = ST.bh(R.get("p_pooled", pd.Series(np.nan, index=R.index)).to_numpy())
    if len(P5):
        P5["q"] = ST.bh(P5["p"].to_numpy())
        minq = P5.groupby(["gene", "feature"])["q"].min()
        R["min_q_lineage"] = minq.reindex(pd.MultiIndex.from_arrays([R["gene"], R["feature"]])).to_numpy()
    else:
        R["min_q_lineage"] = np.nan
    labels = []
    for _, r in R.iterrows():
        if (r.get("q_pooled", np.nan) < q_thr and r["n_powered"] >= s8["min_powered_lineages"]
                and r.get("n_organ_groups", 0) >= s8["min_organ_groups"]
                and r.get("same_sign_frac", 0) >= s8["same_sign_frac_min"] and bool(r.get("lolo_sign_kept", False))):
            labels.append("cross-lineage")
        elif (r["n_testable"] >= 2 and r.get("min_q_lineage", np.nan) < q_thr
              and r.get("p_interaction", np.nan) < PR.dp()["dictionary_s3"]["locus_resolution"]["p_max"]):
            labels.append("lineage-restricted")
        elif r["n_powered"] < s8["min_powered_lineages"]:
            labels.append("lineage-unresolvable" + (" (single-lineage carriers)" if r["single_lineage_carriers"] else ""))
        else:
            labels.append("undetermined")
    R["label"] = labels
    return R


def _interaction_p(data, g, f, rows, lin_raw) -> float:
    """Homoskedastic OLS F-test of d x lineage interactions (testable lineages) in the confounder model."""
    base = pair_fit(data, g, f, extra_rows=rows)
    r = base["rows"]
    gj, fj = data._gpos[g], data.fmeta.index.get_loc(f)
    y = data.target("obs")[r, gj]
    x = data.X[r, fj]
    Z = data.cov.Z(r, data.cov.design_names)
    L = lin_raw[r]
    levs = sorted(set(L))[1:]
    I = np.column_stack([x * (L == l) for l in levs]) if levs else np.zeros((len(x), 0))
    Q0 = ST.orth(np.column_stack([Z, x]))
    Q1 = ST.orth(np.column_stack([Z, x, I]))
    e0 = ST.resid(Q0, y[:, None])
    e1 = ST.resid(Q1, y[:, None])
    k = Q1.shape[1] - Q0.shape[1]
    dfd = len(y) - Q1.shape[1]
    if k <= 0 or dfd <= 0:
        return np.nan
    F = (((e0 ** 2).sum() - (e1 ** 2).sum()) / k) / ((e1 ** 2).sum() / dfd)
    return float(sps.f.sf(F, k, dfd))


# ------------------------------------------------------------------------------------------------
# §7.3 NC3 self-null, §12.1 gene classes, §10.5 CGT hits, §9.3 gate
# ------------------------------------------------------------------------------------------------


def self_null(sn_pairs: pd.DataFrame, effects_of, X_of) -> pd.DataFrame:
    """sn_pairs: (gene, feature). frac_dependent = 1 if the carrier median effect < NC3 threshold (so its mean over
    testable pairs is the NC3 statistic); frac_carriers_dependent = share of carriers below the threshold."""
    thr = PR.dp()["dependency_s2"]["threshold"]
    nmin = PR.stats_params()["min_carriers"]
    out = []
    if sn_pairs.empty:
        return pd.DataFrame(columns=["gene", "feature", "frac_dependent"])
    E = effects_of(sorted(set(sn_pairs["gene"])))
    for g, f in zip(sn_pairs["gene"], sn_pairs["feature"]):
        x = X_of(f)
        y = E[g].to_numpy()
        c = (x == 1) & np.isfinite(y)
        n = int(c.sum())
        med = float(np.median(y[c])) if n else np.nan
        out.append({"gene": g, "feature": f, "n_carriers": n, "testable": n >= nmin, "carrier_median": med,
                    "frac_dependent": float(med < thr) if n else np.nan,
                    "frac_carriers_dependent": float((y[c] < thr).mean()) if n else np.nan})
    return pd.DataFrame(out)


def gene_classes(all_genes, gene_info: pd.DataFrame, eligibility: pd.DataFrame, lam: pd.DataFrame,
                 md: pd.DataFrame, pair_flags: pd.DataFrame, sd_obs: pd.Series) -> pd.DataFrame:
    dp = PR.dp()
    sd_max = dp["dependency_s2"]["no_contrast_sd_max"]
    rmin = dp["eligibility_s2"]["r_within_lin_lib_min"]
    fmin = PR.en_params()["fidelity_class3"]
    tiers = set(dp["eligibility_s2"]["tiers"])
    fid = lam[lam["arm"] == "EN-pred"].set_index("gene")["r2_oof_beyond_design"] if len(lam) else pd.Series(dtype=float)
    el = eligibility.set_index("gene") if "gene" in eligibility else eligibility
    rows = []
    for g in all_genes:
        rw = el["r_within"].get(g, np.nan) if "r_within" in el else np.nan
        if sd_obs.get(g, np.nan) < sd_max:
            c = "1"
        elif gene_info["tier"].get(g) not in tiers:
            c = "2"
        elif not (rw >= rmin) or (g in fid.index and fid[g] < fmin):
            c = "3"
        else:
            pf = pair_flags[pair_flags["gene"] == g] if len(pair_flags) else pair_flags
            if len(pf) == 0:
                c = "3b"
            else:
                pc = pf["pair_class"].astype(int)
                c = "9" if (pc == 9).any() else str(int(pc.max()))
        rows.append({"gene": g, "failure_class": c, "r_within": rw, "fidelity": fid.get(g, np.nan)})
    return pd.DataFrame(rows)


def pair_classes(data, md: pd.DataFrame, delta: pd.DataFrame, locus: pd.DataFrame, amp: pd.DataFrame) -> pd.DataFrame:
    """Pre-R pair classes of MD-strict pairs: 4 expression gating only, 5 too rare, 6 confounded, else 9."""
    if md is None or md.empty:
        return pd.DataFrame(columns=["gene", "feature", "pair_class", "reasons"])
    nmin = PR.dp()["dictionary_s3"]["min_carriers"]
    exempt = set(PR.dp()["clinical_s10"]["expression_biomarker_exempt"])
    dflag = delta.groupby(["gene", "feature"])["flagged"].any() if len(delta) else pd.Series(dtype=bool)
    loc = locus.set_index(["gene", "feature"])["level"] if len(locus) else pd.Series(dtype=object)
    amp_l = amp.set_index(["gene", "feature"])["label"] if len(amp) else pd.Series(dtype=object)
    gsym = data.extras.get("symbols", {})
    out = []
    by_gene = md.groupby("gene")["feature"].apply(set).to_dict()
    for _, r in md.iterrows():
        g, f = r["gene"], r["feature"]
        fam = data.fmeta.at[f, "family"]
        reasons = []
        cls = 9
        if fam == "EXPDOSE" and g in set(data.fmeta.at[f, "genes"]) and gsym.get(g) not in exempt:
            backing = [x for x in by_gene.get(g, set()) if data.fmeta.at[x, "family"] in
                       ("AMP", "SIL-DEL", "SIL-MUT", "SIL-METH", "SIL-other", "SIL")]
            if not backing:
                cls, reasons = 4, ["expression gating only"]
        if cls == 9 and data.fmeta.at[f, "d_type"] == "binary" and r.get("n_carriers", nmin) < nmin:
            cls, reasons = 5, ["too rare"]
        if cls == 9:
            if bool(r.get("collinear_with_design", False)):
                reasons.append("collinear_with_design")
            if bool(dflag.get((g, f), False)):
                reasons.append("confounding delta")
            if loc.get((g, f)) == "locus":
                reasons.append("unresolved co-location")
            if amp_l.get((g, f)) == "cutting artefact":
                reasons.append("cutting artefact")
            if reasons:
                cls = 6
        out.append({"gene": g, "feature": f, "pair_class": cls, "reasons": ";".join(reasons)})
    return pd.DataFrame(out)


def cgt_scores(arms_tab: pd.DataFrame, cgt: pd.DataFrame, k: int) -> pd.DataFrame:
    rows = []
    for arm in arms_tab["arm"].unique():
        s = arms_tab[arms_tab["arm"] == arm]
        for _, t in cgt.iterrows():
            T = [x for x in str(t["target_set"]).split("|") if x]
            f = t["feature_D"]
            best, beta = np.nan, np.nan
            for g in T:
                rk = A.rank_of_feature(s, arm, g, f) if isinstance(f, str) and f else np.nan
                if np.isfinite(rk) and (not np.isfinite(best) or rk < best):
                    best = rk
                    beta = float(s[(s["gene"] == g) & (s["feature"] == f)]["beta"].iloc[0])
            hit = bool(np.isfinite(best) and best <= k)
            sc = np.nan
            if hit and t["sign"] in ("S", "W"):
                sc = bool(beta < 0) if t["sign"] == "S" else bool(beta > 0)
            rows.append({"arm": arm, "test_id": t["test_id"], "rank_of_feature_D": best, "hit": hit,
                         "sign_correct": sc, "sign": t["sign"], "partition": t.get("partition", "")})
    return pd.DataFrame(rows)


def gate(controls: pd.DataFrame, lists: dict) -> pd.DataFrame:
    """recovered (decided on MD-strict) and recovered_<list> per control: (g, primary or secondary feature) in the
    list with the expected sign (S: β < 0, W: β > 0)."""
    out = []
    for _, c in controls.iterrows():
        feats = [x for x in (c.get("feature_primary", ""), c.get("feature_secondary", "")) if isinstance(x, str) and x]
        rec = {"control_id": c["control_id"], "family": c.get("family", ""), "gate_eligible": bool(c.get("gate_eligible", False)),
               "scope": c.get("scope", "")}
        for name, L in lists.items():
            m = L[(L["gene"] == c["dependency_entrez"]) & L["feature"].isin(feats)] if len(L) else L
            if len(m) and c.get("expected_sign") in ("S", "W"):
                m = m[(m["beta"] < 0) if c["expected_sign"] == "S" else (m["beta"] > 0)]
            rec[f"recovered_{name}"] = bool(len(m) > 0)
        rec["recovered"] = rec.get("recovered_MD-strict", False)
        out.append(rec)
    return pd.DataFrame(out)


# ------------------------------------------------------------------------------------------------
# Orchestration
# ------------------------------------------------------------------------------------------------


def run_audits(data, res: dict, toy=None) -> dict:
    a = res["arms"]
    md = res.get("md_strict", pd.DataFrame(columns=["gene", "feature", "beta"]))
    controls = toy.controls if toy is not None else _real_controls()
    cgt = toy.cgt if toy is not None else None
    if toy is not None:
        def effects_of(gs):
            return toy.effects[list(gs)]

        def expression_of(gs):
            return toy.expression[list(gs)]
    else:
        effects_of, expression_of = _real_loaders(data)
    has_obs = (a["arm"] == "S1-obs").any()
    pairs = _pairs_of(res, data, controls) if has_obs else pd.DataFrame(columns=["gene", "feature"])
    sign_of = {(g, f): s for g, f, s in zip(md.get("gene", []), md.get("feature", []), np.sign(md.get("beta", [])))}
    out = {}
    out["confounding_delta"] = confounding_delta(data, pairs) if has_obs else pd.DataFrame(
        columns=["gene", "feature", "block", "flagged"])
    out["locus_resolution"] = locus_resolution(data, pairs, sign_of) if has_obs else pd.DataFrame(
        columns=["gene", "feature", "locus_group", "level", "n_discordant_carriers"])
    out["amplicon_test"] = amplicon_test(data, pairs, effects_of, expression_of) if has_obs else pd.DataFrame(
        columns=["gene", "feature", "label"])
    # shrunk β for the power rule: one EB fit over all S1-obs U rows
    bsh = {}
    if has_obs:
        s = a[a["arm"] == "S1-obs"]
        ok = s["beta"].notna() & s["se"].notna()
        sd = s.loc[ok, "sd_dtilde"].to_numpy()
        b = ST.shrunk_beta(s.loc[ok, "beta"].to_numpy(), s.loc[ok, "se"].to_numpy(), sd)
        bsh = dict(zip(zip(s.loc[ok, "gene"], s.loc[ok, "feature"]), b))
        cl_pairs = _pairs_of(res, data, controls, top_s1=int(PR.pp()["audit"]["cross_lineage_top_s1"]))
        out["cross_lineage"] = cross_lineage(data, cl_pairs, bsh, sign_of)
    else:
        out["cross_lineage"] = pd.DataFrame(columns=["gene", "feature", "label"])
    # self-null
    if toy is not None:
        from . import toysim as TS
        gl = data.extras.get("all_genes", data.genes)
        sn = pd.DataFrame([(g, f) for g in gl for f in TS.self_null_features(toy, g)], columns=["gene", "feature"])
        out["self_null"] = self_null(sn, effects_of, lambda f: toy.features[f].to_numpy())
    else:
        out["self_null"] = pd.DataFrame(columns=["gene", "feature", "frac_dependent"])
    pc = pair_classes(data, md, out["confounding_delta"], out["locus_resolution"], out["amplicon_test"])
    out["pair_classes"] = pc
    if toy is not None:
        sd_obs = toy.effects.std()
        out["gene_classes"] = gene_classes(data.extras.get("all_genes", data.genes), toy.genes,
                                           data.extras["eligibility"], res["lam"], md, pc, sd_obs)
        k = int(PR.endpoints()["P2_primary_clinical"]["k"])
        out["cgt_scores"] = cgt_scores(a, cgt, k) if cgt is not None and len(cgt) else pd.DataFrame(
            columns=["arm", "test_id", "hit", "sign_correct"])
    else:
        out["gene_classes"] = pd.DataFrame(columns=["gene", "failure_class"])
        out["cgt_scores"] = pd.DataFrame(columns=["arm", "test_id", "hit", "sign_correct"])
    lists = {"MD-strict": md}
    for arm in A.ARMS:
        lists[arm] = a[(a["arm"] == arm) & a["selected"]]
    out["gate"] = gate(controls, lists) if controls is not None and len(controls) else pd.DataFrame(
        columns=["control_id", "recovered"])
    return out


def _real_controls():
    from . import paths as P
    f = P.CONFIGS / "positive_controls.csv"
    if not f.exists():
        return None
    c = pd.read_csv(f, dtype=str)
    if "gate_eligible" in c:
        c["gate_eligible"] = c["gate_eligible"].astype(str).str.lower().isin(["true", "1"])
    return c


def _real_loaders(data):
    from . import guard
    from . import p1 as P1
    ids = list(data.lines)

    def effects_of(gs):
        return guard.load_effects(ids, scope="dev", genes=list(gs))

    def expression_of(gs):
        return P1.raw_matrix("expression", columns=list(gs), ids=ids, strict=False).reindex(index=ids, columns=list(gs))
    return effects_of, expression_of
