"""§11 combinatorial screens: catalogue "loss of A -> dependency on partner B" candidates looked up in the downloaded
combinatorial-knockout tables (DESIGN v3.2 §11; amendment A7: MD-global, MD-strict alongside, vs S1-obs).

Sources (paths and columns in configs/m5_m8_params.json:combinatorial):
- Esmaeili Anvar 2024 harmonised re-analysis of five paralog screens (Dede 2020, Parrish 2021, Thompson 2021,
  CHyMErA / Gonatopoulos-Pournatzis 2020, Ito 2021): per-study screened pairs and the authors' per-study hit lists
  (Data_for_Figure_01), per-cell-line effects (Supp_table1, one sheet per study_line).
- Flister et al. 2024 bioRxiv supplement (analysis only, never redistributed): per-pair calls of the preprint's own
  screens.
Labels come from file rows only: replicated (hit in >= 2 screens or >= 2 cell lines), tested-not-replicated,
untested. Outputs carry identifiers and statistics computed here, never a row or value of the external tables.
Validation rate of the paralog candidates of each list vs a screened background matched on sequence-identity decile
and family size (Fisher exact, CI); MD-global vs S1-obs.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

from . import params as PR


def ccfg() -> dict:
    from . import m5
    return m5.cfg()["combinatorial"]


def pair_key(a, b) -> tuple:
    a, b = str(a), str(b)
    return (a, b) if a <= b else (b, a)


# ------------------------------------------------------------------------------------------------
# Readers (external tables -> sets of entrez pair keys; nothing else is kept)
# ------------------------------------------------------------------------------------------------


def _keys(labels, sep: str, sym2ez: dict) -> pd.Series:
    """Pair label 'A<sep>B' -> entrez pair key (None if either symbol is unmapped)."""
    out = []
    for lab in labels:
        parts = str(lab).split(sep)
        if len(parts) != 2:
            out.append(None)
            continue
        a, b = (sym2ez.get(p.strip()) for p in parts)
        out.append(pair_key(a, b) if a and b and a != b else None)
    return pd.Series(out, dtype=object)


def symbols_in(labels, sep: str) -> set:
    return {p.strip() for lab in labels for p in str(lab).split(sep) if p.strip()}


def read_esmaeili(external: Path, sym2ez_fn) -> dict | None:
    """{'screens': {study: {'tested': set, 'hit': set}}, 'lines': [{'study', 'line', 'tested', 'hit'}]} or None."""
    c = ccfg()["esmaeili"]
    base = Path(external) / c["dir"]
    fs, fl = base / c["summary_file"], base / c["per_line_file"]
    if not (fs.exists() and fl.exists()):
        return None
    scr = pd.read_excel(fs, sheet_name=c["screened_sheet"])
    hit = pd.read_excel(fs, sheet_name=c["study_hits_sheet"])
    xl = pd.ExcelFile(fl)
    per = {s: xl.parse(s) for s in xl.sheet_names if s.split("_", 1)[0] in c["studies"]}
    sep = c["pair_sep"]
    syms = set()
    for d in [scr, hit]:
        for col in d.columns:
            syms |= symbols_in(d[col].dropna(), sep)
    for d in per.values():
        syms |= symbols_in(d.iloc[:, 0].dropna(), sep)
    sym2ez = sym2ez_fn(sorted(syms))
    screens = {}
    for s in c["studies"]:
        if s not in scr.columns:
            continue
        t = {k for k in _keys(scr[s].dropna(), sep, sym2ez) if k}
        h = {k for k in _keys(hit[s].dropna(), sep, sym2ez) if k} if s in hit.columns else set()
        screens[c["study_labels"].get(s, s)] = {"tested": t, "hit": h, "raw": s}
    lines = []
    for sheet, d in per.items():
        study = sheet.split("_", 1)[0]
        lab = c["study_labels"].get(study, study)
        k = _keys(d.iloc[:, 0], sep, sym2ez)
        ok = k.notna().to_numpy()
        hitm = (pd.to_numeric(d[c["line_size_col"]], errors="coerce") >= float(c["line_hit_size_min"])) & \
               (pd.to_numeric(d[c["line_effect_col"]], errors="coerce") <= float(c["line_hit_effect_max"]))
        hk = set(k[ok & hitm.to_numpy()])
        if c["line_hits_only_in_study_hits"]:
            hk &= screens.get(lab, {}).get("hit", set())
        lines.append({"study": lab, "line": sheet.split("_", 1)[1] if "_" in sheet else sheet,
                      "tested": set(k[ok]), "hit": hk})
    return {"screens": screens, "lines": lines}


def read_flister(external: Path, sym2ez_fn) -> dict | None:
    c = ccfg()["flister"]
    f = Path(external) / c["dir"] / c["file"]
    if not f.exists():
        return None
    d = pd.read_excel(f, sheet_name=c["sheet"])
    sep = c["pair_sep"]
    sym2ez = sym2ez_fn(sorted(symbols_in(d[c["label_col"]].dropna(), sep)))
    k = _keys(d[c["label_col"]], sep, sym2ez)
    ok = k.notna().to_numpy()
    tested = set(k[ok])
    lines, anyhit = [], set()
    for ln in c["own_lines"]:
        if ln not in d.columns:
            continue
        h = set(k[ok & (d[ln].astype(str).str.strip() == c["hit_value"]).to_numpy()])
        anyhit |= h
        lines.append({"study": c["study"], "line": ln, "tested": tested, "hit": h})
    return {"screens": {c["study"]: {"tested": tested, "hit": anyhit, "raw": c["study"]}}, "lines": lines}


def merge_sources(*srcs) -> dict:
    out = {"screens": {}, "lines": []}
    for s in srcs:
        if s:
            out["screens"].update(s["screens"])
            out["lines"] += s["lines"]
    return out


def label_of(keys, idx: dict) -> pd.DataFrame:
    """Per pair key: screens / lines tested and hit, studies hit, label."""
    c = ccfg()
    rows = []
    for k in keys:
        st = [n for n, s in idx["screens"].items() if k in s["tested"]]
        sh = [n for n, s in idx["screens"].items() if k in s["hit"]]
        lt = sum(k in ln["tested"] for ln in idx["lines"])
        lh = sum(k in ln["hit"] for ln in idx["lines"])
        if not st and not lt:
            lab = c["labels"][2]
        elif len(sh) >= int(c["replicated_min_screens"]) or lh >= int(c["replicated_min_lines"]):
            lab = c["labels"][0]
        else:
            lab = c["labels"][1]
        rows.append({"n_screens_tested": len(st), "n_screens_hit": len(sh), "n_lines_tested": int(lt),
                     "n_lines_hit": int(lh), "screens_hit": "|".join(sh), "label": lab})
    return pd.DataFrame(rows, columns=["n_screens_tested", "n_screens_hit", "n_lines_tested", "n_lines_hit",
                                       "screens_hit", "label"])


# ------------------------------------------------------------------------------------------------
# Candidates, identity, background
# ------------------------------------------------------------------------------------------------


def identity_table(paralogs: pd.DataFrame) -> tuple:
    """({pair key: % identity}, {gene: family size}) from a paralog_top-like table (par<r>, par<r>_score)."""
    k = int(PR.dp()["dictionary_s3"]["n_paralogs"])
    ident = {}
    for g, row in paralogs.iterrows():
        for r in range(1, k + 1):
            p = row.get(f"par{r}")
            if p is None or (isinstance(p, float) and np.isnan(p)):
                continue
            key = pair_key(g, p)
            s = pd.to_numeric(row.get(f"par{r}_score", np.nan), errors="coerce")
            prev = ident.get(key, np.nan)
            ident[key] = float(max(v for v in (prev, s) if np.isfinite(v))) if np.isfinite(prev) or np.isfinite(s) \
                else np.nan
    fam = paralogs["family_size"].to_dict() if "family_size" in paralogs else {}
    return ident, {str(g): v for g, v in fam.items()}


def candidates(named_lists: dict, fmeta: pd.DataFrame, ident: dict, spec_partners: set) -> pd.DataFrame:
    """Rows (list, gene, feature, partner) for list pairs whose feature encodes loss of another gene A (families of
    combinatorial.loss_families); relation paralog (A among the top paralogs of g or vice versa), spec_partner,
    or other; sl_direction = the list's β sign is the loss -> dependency direction."""
    lf = ccfg()["loss_families"]
    rows = []
    for name, L in named_lists.items():
        for g, f, b in zip(L["gene"].astype(str), L["feature"], L["beta"]):
            if f not in fmeta.index:
                continue
            fam = fmeta.at[f, "family"]
            if fam not in lf:
                continue
            for a in [str(x) for x in fmeta.at[f, "genes"] if str(x) != g]:
                key = pair_key(g, a)
                rel = "paralog" if key in ident else ("spec_partner" if (g, a) in spec_partners else "other")
                want = -1.0 if lf[fam] == "neg" else 1.0
                rows.append({"list": name, "gene": g, "feature": f, "family": fam, "partner": a, "relation": rel,
                             "beta": float(b), "sl_direction": bool(np.sign(b) == want), "key": key})
    return pd.DataFrame(rows, columns=["list", "gene", "feature", "family", "partner", "relation", "beta",
                                       "sl_direction", "key"])


def strata(identity: np.ndarray, family: np.ndarray, edges_id: np.ndarray) -> np.ndarray:
    fe = np.asarray(ccfg()["family_size_edges"], dtype=float)
    i = np.digitize(identity, edges_id)
    j = np.digitize(np.nan_to_num(family, nan=-1), fe)
    s = np.char.add(np.char.add(i.astype(str), "|"), j.astype(str)).astype(object)
    s[~np.isfinite(identity)] = None
    return s


def background(idx: dict, ident: dict, famsize: dict) -> pd.DataFrame:
    """Every screened pair with a known sequence identity: key, identity, family size (max of the two genes), label."""
    keys = set().union(*[s["tested"] for s in idx["screens"].values()]) if idx["screens"] else set()
    keys |= set().union(*[ln["tested"] for ln in idx["lines"]]) if idx["lines"] else set()
    keys = sorted(k for k in keys if k in ident and np.isfinite(ident[k]))
    lab = label_of(keys, idx)
    fs = [np.nanmax([famsize.get(a, np.nan), famsize.get(b, np.nan)]) for a, b in keys]
    return pd.concat([pd.DataFrame({"key": keys, "identity": [ident[k] for k in keys], "family_size": fs}), lab],
                     axis=1)


def validation(cands: pd.DataFrame, bg: pd.DataFrame, lists_: list) -> dict:
    """Per list: validation rate of tested SL-direction paralog candidates (Clopper-Pearson CI), the matched
    background rate (strata = identity decile x family-size bin; background minus the list's own keys, strata
    occupied by the list) and Fisher exact OR / CI / p vs that background; plus MD-global vs S1-obs."""
    from . import endpoints as EP
    c = ccfg()
    rep, nrep = c["labels"][0], c["labels"][1]
    level = float(PR.endpoints()["P1_primary"]["ci"]["level"])
    nb = int(c["identity_bins"])
    edges = np.nanquantile(bg["identity"], np.arange(1, nb) / nb) if len(bg) else np.array([])
    bg = bg.assign(stratum=strata(bg["identity"].to_numpy(float), bg["family_size"].to_numpy(float), edges)) \
        if len(bg) else bg.assign(stratum=None)
    out, counts = {}, {}
    for name in lists_:
        cl = cands[(cands["list"] == name) & (cands["relation"] == "paralog") & cands["sl_direction"]]
        cl = cl.drop_duplicates("key")
        t = cl[cl["label"].isin([rep, nrep])]
        k, n = int((t["label"] == rep).sum()), int(len(t))
        lo, hi = EP.clopper_pearson(k, n, level)
        rec = {"n_candidates": int(len(cl)), "n_tested": n, "n_replicated": k, "rate": k / n if n else np.nan,
               "rate_ci": [lo, hi]}
        counts[name] = (k, n)
        if n and len(bg):
            st = strata(t["identity"].to_numpy(float), t["family_size"].to_numpy(float), edges)
            b = bg[~bg["key"].isin(set(cl["key"])) & bg["label"].isin([rep, nrep])]
            occ = set(s for s in st if s is not None)
            bm = b[b["stratum"].isin(occ)]
            br = bm.groupby("stratum")["label"].apply(lambda v: float((v == rep).mean()))
            w = pd.Series([s for s in st if s is not None]).value_counts()
            w = w[w.index.isin(br.index)]
            exp = float((w * br.reindex(w.index)).sum() / w.sum()) if w.sum() else np.nan
            kb, nb_ = int((bm["label"] == rep).sum()), int(len(bm))
            rec.update(n_background=nb_, n_background_replicated=kb, background_rate=kb / nb_ if nb_ else np.nan,
                       matched_expected_rate=exp, n_candidates_in_strata=int(w.sum()))
            rec.update(_fisher(k, n - k, kb, nb_ - kb, level, "vs_background"))
        out[name] = rec
    if "MD-global" in counts and "S1-obs" in counts:
        (k1, n1), (k2, n2) = counts["MD-global"], counts["S1-obs"]
        out["MD-global_vs_S1-obs"] = _fisher(k1, n1 - k1, k2, n2 - k2, level, "")
    if "MD-strict" in counts and "S1-obs" in counts:
        (k1, n1), (k2, n2) = counts["MD-strict"], counts["S1-obs"]
        out["MD-strict_vs_S1-obs"] = _fisher(k1, n1 - k1, k2, n2 - k2, level, "")
    return out


def _fisher(a, b, c, d, level, prefix) -> dict:
    p = f"{prefix}_" if prefix else ""
    t = np.array([[a, b], [c, d]])
    if t.sum(1).min() == 0:
        return {f"{p}odds_ratio": np.nan, f"{p}or_ci": [np.nan, np.nan], f"{p}fisher_p": np.nan}
    r = sps.contingency.odds_ratio(t)
    ci = r.confidence_interval(confidence_level=level)
    return {f"{p}odds_ratio": float(r.statistic), f"{p}or_ci": [float(ci.low), float(ci.high)],
            f"{p}fisher_p": float(sps.fisher_exact(t)[1])}


def lookup(named_lists: dict, fmeta: pd.DataFrame, paralogs: pd.DataFrame, spec_partners: set, idx: dict) -> tuple:
    """(candidate table with labels, validation dict)."""
    ident, famsize = identity_table(paralogs)
    cand = candidates(named_lists, fmeta, ident, spec_partners)
    cand = pd.concat([cand.reset_index(drop=True), label_of(list(cand["key"]), idx)], axis=1)
    cand["identity"] = [ident.get(k, np.nan) for k in cand["key"]]
    cand["family_size"] = [np.nanmax([famsize.get(a, np.nan), famsize.get(b, np.nan)]) for a, b in cand["key"]]
    bg = background(idx, ident, famsize)
    val = validation(cand, bg, list(named_lists))
    val["_n_background_pairs"] = int(len(bg))
    val["_sources"] = {n: {"n_pairs_tested": len(s["tested"]), "n_pairs_hit": len(s["hit"])}
                       for n, s in idx["screens"].items()}
    val["_n_cell_lines"] = len(idx["lines"])
    return cand.drop(columns="key"), val
