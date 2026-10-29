"""§12 limits and candidate cards (DESIGN v3.2 §12; amendment A7: classes on the MD-global catalogue, MD-strict
alongside).

Pair classes (§12.2) start from ``audit.pair_classes`` (4 expression gating, 5 too rare, 6 confounded by delta /
co-location / cutting / collinear_with_design, 9 valid) and add the remaining inputs: class 6 culture flag (§8.6),
hypermutation and TP53 (§8.7) and subtype-confounded (§5.5 primary-disease block); class 5 R power < power_min
(pre-tag power from configs/catalogue/r_requests.csv, or r_tests.csv); class 7 powered but not replicated in R;
class 8 KO-drug discordance (catalogue pair matching a CGT test of a discordant KO class, a discordant ceiling label,
or the opposite sign; plus gene-level rows for concordant CGT tests whose ceiling is data-null / data-contradicted).
Classes that R can still change are marked 'pending R' until results/replication/r_tests.csv exists.
Gene classes (§12.1) reuse ``audit.gene_classes``. Candidate cards (§12.3) read frozen outputs only.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import audit as AU
from . import params as PR


def lcfg() -> dict:
    from . import m5
    return m5.cfg()["limits"]


class FeatureData:
    """The two attributes of ``pipeline.Data`` that ``audit.pair_classes`` reads (fmeta, extras['symbols'])."""

    def __init__(self, fmeta: pd.DataFrame, symbols):
        self.fmeta = fmeta
        self.extras = {"symbols": dict(symbols)}


def fmeta_from_dictionary() -> pd.DataFrame:
    """Feature metadata of configs/dictionary_v1.json in the ``pipeline.Data.fmeta`` conventions."""
    from . import endpoint_io as IO
    fl = IO.feature_meta()
    fm = pd.DataFrame(index=fl.index)
    fm["family"] = fl["family"]
    fm["d_type"] = fl["d_type"]
    fm["stage"] = fl.get("stage")
    fm["genes"] = [tuple(map(str, g)) if isinstance(g, (list, tuple)) else () for g in fl["genes"]]
    fm["locus_group"] = fl.get("locus_group")
    return fm


# ------------------------------------------------------------------------------------------------
# Inputs keyed by (gene, feature)
# ------------------------------------------------------------------------------------------------


def _key(df, col):
    if df is None or not len(df) or col not in df:
        return {}
    return dict(zip(zip(df["gene"].astype(str), df["feature"]), df[col]))


def r_status(rt: pd.DataFrame | None) -> dict:
    """{(g, d): (replicated, p_one, power)} from r_tests.csv (empty before the single R run)."""
    if rt is None or not len(rt):
        return {}
    rep = rt["replicated"].astype(str).str.lower().isin(["true", "1"])
    pw = rt["power"] if "power" in rt else pd.Series(np.nan, index=rt.index)
    p1 = rt["p_one"] if "p_one" in rt else pd.Series(np.nan, index=rt.index)
    return {(str(g), f): (bool(r), float(p), float(w)) for g, f, r, p, w in
            zip(rt["gene"], rt["feature"], rep, p1, pw)}


def ceiling_labels(ceil: pd.DataFrame | None) -> dict:
    if ceil is None or not len(ceil):
        return {}
    col = next((c for c in lcfg()["ceiling_label_cols"] if c in ceil), None)
    return dict(zip(ceil["test_id"], ceil[col])) if col else {}


def cgt_pairs(cgt: pd.DataFrame) -> dict:
    """{(target gene, feature_D): [test rows]}."""
    out = {}
    for _, t in cgt.iterrows():
        f = t.get("feature_D")
        if not isinstance(f, str) or not f:
            continue
        for g in t["targets"]:
            out.setdefault((str(g), f), []).append(t)
    return out


def discordance(g, f, beta, cg: dict, ceil: dict) -> str:
    """KO-drug discordance reasons for a catalogue pair (empty string if none)."""
    L = lcfg()
    why = []
    for t in cg.get((g, f), []):
        lab = ceil.get(t["test_id"], "")
        if t.get("ko_class") in L["ko_discordant_classes"]:
            why.append(f"{t['test_id']}: KO class {t['ko_class']}")
        elif lab in L["ceiling_discordant_labels"]:
            why.append(f"{t['test_id']}: ceiling {lab}")
        elif t.get("sign") in ("S", "W") and np.isfinite(beta) and np.sign(beta) != (-1 if t["sign"] == "S" else 1):
            why.append(f"{t['test_id']}: sign opposite to the CGT test")
    return "; ".join(why)


# ------------------------------------------------------------------------------------------------
# §12.2 pair classes
# ------------------------------------------------------------------------------------------------


def pair_classes(fd: FeatureData, cat: pd.DataFrame, delta: pd.DataFrame, locus: pd.DataFrame, amp: pd.DataFrame,
                 culture: pd.DataFrame | None, hyper: pd.DataFrame | None, power: dict, rt: pd.DataFrame | None,
                 cgt: pd.DataFrame | None, ceil: pd.DataFrame | None, gene_universe=None) -> pd.DataFrame:
    PC = dict(PR.pp()["pair_classes"], **lcfg()["pair_class_codes"])
    pmin = float(PR.dp()["replication_s6"]["power_min"])
    pending = lcfg()["pending"]
    base = AU.pair_classes(fd, cat, delta, locus, amp)
    base = base.merge(cat[["gene", "feature", "beta", "lists", "in_md_global", "in_md_strict"]], on=["gene", "feature"],
                      how="left")
    cul = _key(culture, "culture_confounded")
    culv = _key(culture, "V_max_variable")
    hy = {}
    if hyper is not None and len(hyper):
        for g, f, t, fl in zip(hyper["gene"], hyper["feature"], hyper["test"], hyper["failed"]):
            if bool(fl):
                hy.setdefault((g, f), []).append(t)
    sub = set()
    if delta is not None and len(delta):
        d = delta[delta["block"].astype(str).str.startswith("primary_disease") & delta["flagged"].astype(bool)]
        sub = set(zip(d["gene"].astype(str), d["feature"]))
    R = r_status(rt)
    have_R = rt is not None
    cg = cgt_pairs(cgt) if cgt is not None else {}
    cl = ceiling_labels(ceil)
    rows = []
    for _, r in base.iterrows():
        g, f = r["gene"], r["feature"]
        k = (g, f)
        cls = int(r["pair_class"])
        reasons = [x for x in str(r["reasons"]).split(";") if x and x != "nan"]
        extra = []
        if bool(cul.get(k, False)):
            extra.append(f"culture flag ({culv.get(k, '')})")
        extra += [f"{t} subset sign/size not kept" for t in hy.get(k, [])]
        if k in sub:
            extra.append("subtype-confounded")
        if cls in (PC["valid"], PC["confounded"]) and extra:
            cls = PC["confounded"]
            reasons += extra
        pw = R[k][2] if k in R and np.isfinite(R[k][2]) else power.get(k, np.nan)
        if cls in (PC["valid"], PC["confounded"]) and np.isfinite(pw) and pw < pmin:
            cls = PC["too_rare"]
            reasons.append("R power < power_min")
        disc = discordance(g, f, r["beta"], cg, cl)
        if cls == PC["valid"] and disc:
            cls = PC["ko_drug_discordance"]
            reasons.append(f"KO-drug discordance ({disc})")
        pre = cls
        status = "final"
        rep, p_one = (R[k][0], R[k][1]) if k in R else (np.nan, np.nan)
        if cls in (PC["valid"], PC["ko_drug_discordance"]):
            if not have_R:
                status = pending
            elif k in R and np.isfinite(pw) and pw >= pmin and not R[k][0]:
                cls = PC["powered_not_replicated"]
                reasons.append("powered, not replicated in R")
            elif k not in R:
                status = "no R test"
        elif cls in (PC["confounded"],) and not np.isfinite(pw):
            status = pending if not have_R else status
        rows.append({"gene": g, "feature": f, "family": fd.fmeta.at[f, "family"], "lists": r["lists"],
                     "in_md_global": bool(r["in_md_global"]), "in_md_strict": bool(r["in_md_strict"]),
                     "source": "catalogue", "pair_class_preR": pre, "pair_class": cls, "class_status": status,
                     "power": pw, "R_replicated": rep, "R_p_one": p_one, "reasons": ";".join(reasons)})
    out = pd.DataFrame(rows)
    # gene-level class-8 rows from concordant CGT tests whose ceiling is discordant (not in the catalogue)
    have = set(zip(out["gene"], out["feature"])) if len(out) else set()
    extra_rows = []
    U = None if gene_universe is None else set(map(str, gene_universe))
    for (g, f), tests in cg.items():
        for t in tests:
            lab = cl.get(t["test_id"], "")
            if str(t.get("ko_class", "")).startswith("concordant") and lab in lcfg()["ceiling_discordant_labels"] \
                    and (g, f) not in have and (U is None or g in U):
                extra_rows.append({"gene": g, "feature": f, "family": fd.fmeta.at[f, "family"] if f in fd.fmeta.index
                                   else "", "lists": "", "in_md_global": False, "in_md_strict": False,
                                   "source": "cgt_ceiling", "pair_class_preR": PC["ko_drug_discordance"],
                                   "pair_class": PC["ko_drug_discordance"], "class_status": "final", "power": np.nan,
                                   "R_replicated": np.nan, "R_p_one": np.nan,
                                   "reasons": f"KO-drug discordance ({t['test_id']}: ceiling {lab})"})
                have.add((g, f))
    if extra_rows:
        out = pd.concat([out, pd.DataFrame(extra_rows)], ignore_index=True)
    return out


# ------------------------------------------------------------------------------------------------
# §12.1 gene classes
# ------------------------------------------------------------------------------------------------


def gene_classes(universe, gene_info: pd.DataFrame, lam: pd.DataFrame, pc: pd.DataFrame, attribution=None,
                 controls: pd.DataFrame | None = None) -> pd.DataFrame:
    """audit.gene_classes on the MD-global pairs (and CGT-ceiling class-8 rows); the MD-strict aggregation is
    reported alongside. M1 named share < attribution_s4 threshold -> class 3 when the Tier B attribution exists."""
    pending = lcfg()["pending"]
    gi = gene_info.copy()
    elig = pd.DataFrame({"gene": gi.index.astype(str), "r_within": gi["r_within"].to_numpy()})
    sd = gi["sd"]
    universe = [str(g) for g in universe]
    out = {}
    for name, m in (("", pc["in_md_global"] | (pc["source"] == "cgt_ceiling")),
                    ("_md_strict", pc["in_md_strict"] | (pc["source"] == "cgt_ceiling"))):
        sub = pc[m.astype(bool)][["gene", "feature", "pair_class"]]
        out[name] = AU.gene_classes(universe, gi, elig, lam, sub, sub, sd).set_index("gene")
    gc = out[""].copy()
    gc["failure_class_md_strict"] = out["_md_strict"]["failure_class"].reindex(gc.index)
    share_min = float(PR.dp()["attribution_s4"]["m1_named_share_min"])
    if attribution is not None and len(attribution):
        a = attribution.set_index(attribution["gene"].astype(str))[lcfg()["attribution_col"]]
        gc["m1_named_share"] = a.reindex(gc.index)
        low = gc["m1_named_share"] < share_min
        for col in ("failure_class", "failure_class_md_strict"):
            gc.loc[low & ~gc[col].isin(["1", "2"]), col] = "3"
        gc["m1_attribution"] = "run"
    else:
        gc["m1_named_share"] = np.nan
        gc["m1_attribution"] = "not run (Tier B)"
    pg = pc[pc["in_md_global"].astype(bool) | (pc["source"] == "cgt_ceiling")]
    summ = pg.groupby("gene").apply(lambda s: ";".join(f"{c}:{n}" for c, n in s["pair_class"].value_counts()
                                                       .sort_index().items()), include_groups=False) \
        if len(pg) else pd.Series(dtype=object)
    pend = pg[pg["class_status"] == pending].groupby("gene").size() if len(pg) else pd.Series(dtype=int)
    gc["pair_classes"] = summ.reindex(gc.index).fillna("")
    gc["n_pairs"] = pg.groupby("gene").size().reindex(gc.index).fillna(0).astype(int) if len(pg) else 0
    gc["class_status"] = np.where(pend.reindex(gc.index).fillna(0).to_numpy() > 0, pending, "final")
    gc["symbol"] = gi["symbol"].reindex(gc.index)
    gc["tier"] = gi["tier"].reindex(gc.index)
    gc["eligible"] = gi["eligible"].reindex(gc.index)
    gc["sd_D"] = sd.reindex(gc.index)
    gc["pan_lethal"] = gi["common_essential"].reindex(gc.index).fillna(False).astype(bool)
    if controls is not None and len(controls):
        cf = controls.groupby("gene")["family"].agg(lambda v: "|".join(sorted(set(v))))
        cs = controls.groupby("gene")["scope"].agg(lambda v: "|".join(sorted(set(v))))
        gc["control_families"] = cf.reindex(gc.index).fillna("")
        gc["control_scope"] = cs.reindex(gc.index).fillna("")
    return gc.reset_index()


def universe_genes(gene_info: pd.DataFrame, controls: pd.DataFrame) -> list:
    """§12.1 universe: variable ∪ selective genes ∪ control dependency genes."""
    g = set(gene_info.index[gene_info["variable"].astype(bool) | gene_info["selective"].astype(bool)].astype(str))
    g |= set(controls["gene"].astype(str)) if controls is not None and len(controls) else set()
    return sorted(g)


# ------------------------------------------------------------------------------------------------
# §12.2 extra labels, Q1b
# ------------------------------------------------------------------------------------------------


def pair_labels(fd: FeatureData, cat: pd.DataFrame, s1_hits: pd.DataFrame, en_sel: pd.DataFrame,
                gene_info: pd.DataFrame, rt: pd.DataFrame | None) -> pd.DataFrame:
    """processing-limited (own-gene CN family: own-arm CN / CYCLOPS), pan-lethal (DepMap common essential),
    predictor-blind (S1-obs BH, not EN-pred selected, replicated in R), predictor-only not replicated (EN-pred
    selected, not S1-obs BH, tested and not replicated in R)."""
    pending = lcfg()["pending"]
    plf = set(lcfg()["processing_limited_families"])
    ce = gene_info["common_essential"].astype(bool).to_dict()
    R = r_status(rt)
    have_R = rt is not None
    en_keys = set(zip(en_sel["gene"].astype(str), en_sel["feature"]))
    s1_keys = set(zip(s1_hits["gene"].astype(str), s1_hits["feature"]))
    rows = []
    for g, f, ls in zip(cat["gene"], cat["feature"], cat["lists"]):
        fam = fd.fmeta.at[f, "family"]
        rows.append({"gene": g, "feature": f, "set": ls, "processing_limited": bool(fam in plf and g in
                                                                                   fd.fmeta.at[f, "genes"]),
                     "pan_lethal": bool(ce.get(g, False)), "predictor_blind": "", "predictor_only_not_replicated": ""})
    for g, f in sorted(s1_keys - en_keys):
        lab = pending if not have_R else (bool(R[(g, f)][0]) if (g, f) in R else "no R test")
        rows.append({"gene": g, "feature": f, "set": "S1-obs BH, not EN-pred selected", "processing_limited": "",
                     "pan_lethal": bool(ce.get(g, False)), "predictor_blind": lab, "predictor_only_not_replicated": ""})
    if have_R:
        for g, f in sorted((en_keys - s1_keys) & set(R)):
            rows.append({"gene": g, "feature": f, "set": "EN-pred selected, not S1-obs BH", "processing_limited": "",
                         "pan_lethal": bool(ce.get(g, False)), "predictor_blind": "",
                         "predictor_only_not_replicated": not R[(g, f)][0]})
    return pd.DataFrame(rows)


def q1b(cl: pd.DataFrame | None, cat: pd.DataFrame, cgt: pd.DataFrame | None, rt: pd.DataFrame | None) -> dict:
    """§8.3 Q1b on MD-global determinants: Fisher OR (CI) of CGT membership and of R replication, cross-lineage vs
    lineage-restricted. A null result is a finding."""
    from . import combinatorial as CB
    pending = lcfg()["pending"]
    level = float(PR.endpoints()["P1_primary"]["ci"]["level"])
    if cl is None or not len(cl):
        return {"status": "cross_lineage.csv missing"}
    mg = set(zip(cat.loc[cat["in_md_global"], "gene"], cat.loc[cat["in_md_global"], "feature"]))
    lab = {(g, f): x for g, f, x in zip(cl["gene"].astype(str), cl["feature"], cl["label"]) if (g, f) in mg}
    groups = {"cross-lineage": [k for k, v in lab.items() if v == "cross-lineage"],
              "lineage-restricted": [k for k, v in lab.items() if v == "lineage-restricted"]}
    cg = set(cgt_pairs(cgt)) if cgt is not None else set()
    out = {"n_cross_lineage": len(groups["cross-lineage"]), "n_lineage_restricted": len(groups["lineage-restricted"])}

    def tab(fn):
        a = sum(fn(k) for k in groups["cross-lineage"])
        c = sum(fn(k) for k in groups["lineage-restricted"])
        return CB._fisher(a, len(groups["cross-lineage"]) - a, c, len(groups["lineage-restricted"]) - c, level, "") | \
            {"k_cross_lineage": int(a), "k_lineage_restricted": int(c)}

    out["cgt_membership"] = tab(lambda k: k in cg)
    if rt is None:
        out["r_replication"] = {"status": pending}
    else:
        R = r_status(rt)
        out["r_replication"] = tab(lambda k: bool(R.get(k, (False,))[0]))
    return out


# ------------------------------------------------------------------------------------------------
# §12.3 candidate cards
# ------------------------------------------------------------------------------------------------


def known_recovery(ridx, arms_list, controls: pd.DataFrame, cgt_testable: pd.DataFrame, cat: pd.DataFrame,
                   rt: pd.DataFrame | None, ceil: dict, access: dict, symbols) -> pd.DataFrame:
    """One row per positive control and per testable CGT test: best rank_of_feature_D and its sign per arm,
    recovery in MD-global / MD-strict (expected sign), R result, ceiling label, access class."""
    pending = lcfg()["pending"]
    R = r_status(rt)
    lists_ = {n: cat[cat[f"in_{n.lower().replace('-', '_')}"]] for n in ("MD-global", "MD-strict")}
    rows = []

    def one(kind, ident, family, pairs, exp_sign, extra):
        rec = {"kind": kind, "id": ident, "family_or_tier": family, "expected_sign": exp_sign,
               "genes": "|".join(dict.fromkeys(g for g, _ in pairs)),
               "symbols": "|".join(dict.fromkeys(str(symbols.get(g, "")) for g, _ in pairs)),
               "features": "|".join(dict.fromkeys(f for _, f in pairs if isinstance(f, str))),
               "access_class": "|".join(sorted({str(access.get(k, "")) for k in pairs if access.get(k)}))}
        for arm in arms_list:
            best, b = np.nan, np.nan
            for g, f in pairs:
                rk, bt = ridx.rank_of(arm, g, f)
                if np.isfinite(rk) and (not np.isfinite(best) or rk < best):
                    best, b = rk, bt
            rec[f"rank_{arm}"] = best
            rec[f"sign_{arm}"] = np.sign(b) if np.isfinite(b) else np.nan
        for n, L in lists_.items():
            m = L[[k in set(pairs) for k in zip(L["gene"], L["feature"])]]
            if exp_sign in ("S", "W"):
                m = m[(m["beta"] < 0) if exp_sign == "S" else (m["beta"] > 0)]
            rec[f"recovered_{n}"] = bool(len(m))
        tested = [R[k] for k in pairs if k in R]
        rec["R_status"] = pending if rt is None else ("tested" if tested else "no R test")
        rec["R_replicated"] = any(t[0] for t in tested) if tested else np.nan
        rec["R_p_one_min"] = min(t[1] for t in tested) if tested else np.nan
        rec.update(extra)
        rows.append(rec)

    for cid, s in controls.groupby("control_id", sort=False):
        one("control", cid, s["family"].iloc[0], list(zip(s["gene"], s["feature"])), s["expected_sign"].iloc[0],
            {"scope": s["scope"].iloc[0], "ceiling_label": ""})
    for _, t in cgt_testable.iterrows():
        pairs = [(str(g), t["feature_D"]) for g in t["targets"]]
        one("cgt", t["test_id"], t.get("tier", ""), pairs, t.get("sign", ""),
            {"scope": t.get("partition", ""), "ceiling_label": ceil.get(t["test_id"], ""),
             "proxy": t.get("proxy", "")})
    return pd.DataFrame(rows)


def shortlist(pc: pd.DataFrame, rt: pd.DataFrame | None, cl: pd.DataFrame | None, per_lineage: pd.DataFrame | None,
              cat: pd.DataFrame, access: dict, symbols, delta_flags: dict, locus: dict, comb: dict, cgt: pd.DataFrame | None,
              setA: dict, curated: set) -> pd.DataFrame:
    """Class-9 MD-global pairs replicated in R, ranked by cross-lineage label precedence, R one-sided p, number of
    powered lineages (limits_s12.cards.shortlist_rank)."""
    from . import m5
    prec = m5.cfg()["cards"]["cross_lineage_precedence"]
    PC = PR.pp()["pair_classes"]
    R = r_status(rt)
    cols = ["gene", "symbol", "feature", "access_class", "cross_lineage_label", "n_powered_lineages",
            "lineages_beta", "n_carriers_D", "n_carriers_R", "R_p_one", "delta_flags", "locus_resolution",
            "combinatorial_label", "cgt_tests", "setA_dependent_frac", "not_in_curated_lists"]
    if rt is None or pc is None or not len(pc):
        return pd.DataFrame(columns=cols)
    ok = pc[(pc["pair_class"] == PC["valid"]) & pc["in_md_global"].astype(bool)]
    ok = ok[[R.get((g, f), (False,))[0] for g, f in zip(ok["gene"], ok["feature"])]]
    lab = {(str(g), f): (x, n) for g, f, x, n in zip(cl["gene"], cl["feature"], cl["label"], cl["n_powered"])} \
        if cl is not None and len(cl) else {}
    plin = {}
    if per_lineage is not None and len(per_lineage):
        for (g, f), s in per_lineage.groupby(["gene", "feature"]):
            plin[(str(g), f)] = s
    cg = cgt_pairs(cgt) if cgt is not None else {}
    nc = dict(zip(zip(cat["gene"], cat["feature"]), cat["n_carriers"]))
    rtc = _key(rt, "n_carriers")
    rows = []
    for g, f in zip(ok["gene"], ok["feature"]):
        k = (g, f)
        L, npw = lab.get(k, ("undetermined", 0))
        s = plin.get(k)
        pw = s[s["power"] >= float(PR.dp()["subtypes_s8"]["cross_lineage"]["powered_lineage_power_min"])] \
            if s is not None else None
        lb = ";".join(f"{a}:{b:.3g}" for a, b in zip(pw["lineage"], pw["beta"])) if pw is not None else ""
        genes_d = set(m5.feature_genes_of_id(f))
        rows.append({"gene": g, "symbol": symbols.get(g, ""), "feature": f, "access_class": access.get(k, ""),
                     "cross_lineage_label": L.split(" (")[0] if isinstance(L, str) else L,
                     "n_powered_lineages": int(npw) if np.isfinite(npw) else 0, "lineages_beta": lb,
                     "n_carriers_D": nc.get(k, np.nan), "n_carriers_R": rtc.get(k, np.nan), "R_p_one": R[k][1],
                     "delta_flags": delta_flags.get(k, ""), "locus_resolution": locus.get(k, ""),
                     "combinatorial_label": comb.get(k, ""), "cgt_tests": "|".join(t["test_id"] for t in cg.get(k, [])),
                     "setA_dependent_frac": setA.get(g, np.nan),
                     "not_in_curated_lists": not (k in curated or any((g, h) in curated for h in genes_d))})
    out = pd.DataFrame(rows, columns=cols)
    if len(out):
        out["_p"] = out["cross_lineage_label"].map(lambda x: prec.index(x) if x in prec else len(prec))
        out = out.sort_values(["_p", "R_p_one", "n_powered_lineages"], ascending=[True, True, False],
                              kind="mergesort").drop(columns="_p").reset_index(drop=True)
        out.insert(0, "rank", np.arange(1, len(out) + 1))
    return out
