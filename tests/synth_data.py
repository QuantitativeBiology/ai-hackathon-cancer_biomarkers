"""Small synthetic ``pipeline.Data`` objects for the engine unit tests (no real data)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from biomark import en as EN
from biomark import params as PR
from biomark import pipeline as PL

LIN = ("LinA", "LinB", "LinC", "LinD")
LIB = ("Avana-only", "KY-only")


def make_data(n=240, G=4, seed=0, plant_beta=-1.0, obs_loader=None, n_bg=12) -> PL.Data:
    rng = np.random.default_rng(seed)
    lines = pd.Index([f"SYN{i:04d}" for i in range(n)], name="ModelID")
    lin = rng.choice(LIN, size=n)
    lib = rng.choice(LIB, size=n, p=[0.7, 0.3])
    pat = np.array([f"P{i:04d}" for i in range(n)], dtype=object)
    for i in range(0, 20, 2):           # 10 two-line patients
        pat[i + 1] = pat[i]
    pat[30] = pat[31] = pat[32] = "P_TRIPLE"
    cov = pd.DataFrame({
        "lineage": lin, "library_group": lib,
        "cellformat": rng.choice(["adherent", "suspension"], size=n, p=[0.8, 0.2]),
        "serumfree": rng.choice(["F", "T"], size=n, p=[0.9, 0.1]),
        "medium_family": rng.choice(["RPMI", "DMEM"], size=n),
        "doubling_time": rng.normal(40, 8, n), "DT_unknown": (rng.random(n) < 0.1).astype(float),
        "MSI": np.where(rng.random(n) < 0.1, "MSI", "MSS"), "DMG_BURDEN": rng.normal(5, 0.5, n),
        "TP53": np.where(rng.random(n) < 0.5, "mut", "WT"), "PatientID": pat,
        "primary_disease": rng.choice(["d1", "d2"], size=n), "lineage_raw": lin}, index=lines)
    genes = [f"{9100001 + j}" for j in range(G)]
    cn_absent = rng.random(n) < 0.05
    feats, meta = {}, {}

    def add(fid, v, fam, stage, dtype="binary", scope="shared", gene=(), vg=None):
        feats[fid] = np.asarray(v, float)
        meta[fid] = {"family": fam, "stage_family": fam, "stage": stage, "d_type": dtype, "scope": scope,
                     "genes": tuple(gene), "locus_group": vg, "variant_group": vg}

    hot = (rng.random(n) < 0.12).astype(float)
    add("HOT:9200001", hot, "HOT", "G", gene=("9200001",))
    for k in range(n_bg):
        add(f"LOFBI:{9300001 + k}", (rng.random(n) < rng.uniform(0.05, 0.2)).astype(float), "LOFBI", "G",
            gene=(f"{9300001 + k}",))
    add("SIG:MSI", (cov["MSI"] == "MSI").astype(float), "SIG", "G")
    add("SIG:DMG_BURDEN", cov["DMG_BURDEN"].to_numpy(), "SIG", "G", dtype="continuous")
    add("SIL:LINDUMMY", (lin == "LinB").astype(float), "SIL-other", "X")
    amp = (rng.random(n) < 0.1).astype(float)
    amp[cn_absent] = np.nan
    add(f"AMP:{genes[1]}", amp, "AMP", "G", scope="gene_specific", gene=(genes[1],), vg="AMP:LG1")
    add("PROG:1", rng.normal(size=n), "PROG", "X", dtype="continuous")
    for g in genes:
        add(f"EXPDOSE:{g}", rng.normal(size=n), "EXPDOSE", "X", dtype="continuous", scope="gene_specific", gene=(g,))
    fmeta = pd.DataFrame.from_dict(meta, orient="index")
    X = np.column_stack([feats[f] for f in fmeta.index])
    pos = pd.Series(np.arange(len(fmeta)), index=fmeta.index)
    shared = [f for f in fmeta.index if fmeta.at[f, "scope"] == "shared"]
    Dg, forced = {}, {}
    for g in genes:
        own = [f"EXPDOSE:{g}"] + ([f"AMP:{g}"] if f"AMP:{g}" in pos else [])
        Dg[g] = np.sort(pos[shared + own].to_numpy())
        forced[g] = set(own)
    lin_eff = {L: rng.normal(0, 0.2) for L in LIN}
    Y = np.column_stack([np.array([lin_eff[L] for L in lin]) + rng.normal(0, 0.5, n) for _ in genes])
    Y[:, 0] += plant_beta * hot
    Y[:, 1] += -0.8 * np.nan_to_num(amp)
    pred = 0.6 * Y + rng.normal(0, 0.2, Y.shape)
    rows_ok = np.ones_like(Y, dtype=bool)
    rows_ok[lib == "KY-only", G - 1] = False
    Y[~rows_ok] = np.nan
    own_cn = rng.normal(0, 0.3, (n, G))
    own_cn[:, 1] += 1.5 * np.nan_to_num(amp)
    own_cn[cn_absent] = np.nan
    expg = (rng.random((n, G)) < 0.9).astype(float)
    folds = rng.integers(0, 5, n)
    subs = EN.make_subsamples(lines, cov["lineage"], cov["PatientID"], 123, PR.en_params()["n_pairs"])
    sub_masks = [lines.isin(s) for s in subs]
    na = np.column_stack([cn_absent.astype(float), np.zeros(n), np.zeros(n)])
    Yc = Y.copy()
    return PL.Data(lines=lines, genes=genes, X=X, fmeta=fmeta, Dg=Dg, forced=forced, cov_df=cov, na_flags=na,
                   own_cn=own_cn, expgate=expg, rows_ok=rows_ok, pred=np.where(rows_ok, pred, 0.0), folds=folds,
                   sub_masks=sub_masks, patients=pat.astype(str), access_fn=lambda g, f: "proxy",
                   obs_loader=obs_loader or (lambda: Yc), extras={"toy": True})
