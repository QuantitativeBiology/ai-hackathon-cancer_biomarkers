"""Streamed ("blocked") orchestration of the four arms for the full real-D universe (orchestration only: results equal
``arms.run_arms``; tests/test_blocked.py).

Why: one in-memory arms table for 2,193 genes is ~26 M rows (tens of GB). Here the parent never holds it:
1. Genes are cut into chunks at NaN-pattern boundaries (≤ chunk genes). One S1 pool scans all chunks; every chunk's
   raw scan is written at once to a stage part (``arms.finish_s1_local`` columns) and only slim categorical columns
   (arm, gene, feature, p, beta, test family, collinearity) stay in memory. One pool per EN arm does the same
   (``arms.finish_en_local``), with workers evicting FWL caches to one row pattern (cache only).
2. Global steps on the slim columns, exactly as ``finish_s1`` / ``finish_en`` / ``add_dedup`` / ``run_arms``: BH per
   (arm, test family), S1 ranks (p, -|beta|, feature; stable), EN ranks (``en.rank_en``), de-duplicated ranks,
   observed p-cutoffs; MD-strict / MD-global / EN-obs-strict with ``arms.strict`` / ``arms.md_global`` on the
   EN-selected candidates and their S1-obs rows (read back from the stage parts).
3. Each final part = ``arms.annotate`` of the part's S1 and EN rows with the global columns attached, written to
   ``<out_dir>/part-NNNNN.parquet`` with categorical string columns (read with ``endpoint_io.read_arms``).
"""
from __future__ import annotations

import hashlib
import shutil
import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import arms as A
from . import en as EN
from . import params as PR
from . import paths as P
from . import stats as ST


def pattern_chunks(data, genes, chunk: int) -> list:
    """Genes sorted by NaN pattern (as the parallel arms paths sort them), cut at pattern changes and every chunk."""
    gi = {g: i for i, g in enumerate(data.genes)}
    order = sorted(genes, key=lambda g: data.rows_ok[:, gi[g]].tobytes())
    out, cur, key = [], [], None
    for g in order:
        k = data.rows_ok[:, gi[g]].tobytes()
        if cur and (k != key or len(cur) >= chunk):
            out.append(cur)
            cur = []
        cur.append(g)
        key = k
    if cur:
        out.append(cur)
    return out


def to_categorical(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for c in df.columns:
        if df[c].dtype == object:
            vals = df[c]
            if vals.map(lambda v: v is None or isinstance(v, str) or (isinstance(v, float) and np.isnan(v))).all():
                df[c] = vals.astype("category")
    return df


def _write(df: pd.DataFrame, path: Path) -> None:
    to_categorical(df).to_parquet(path, index=False)


def _slim(df: pd.DataFrame, cols, part: int) -> pd.DataFrame:
    s = df[cols].copy()
    for c in ("arm", "gene", "feature", "test_family"):
        if c in s:
            s[c] = s[c].astype("category")
    s["part"] = np.int32(part)
    return s


def run_streamed(data, arms=A.ARMS, out_dir: Path | None = None, stage_dir: Path | None = None, workers_s1: int = 1,
                 workers_en: int = 1, chunk: int | None = None, genes=None, log=print) -> dict:
    rt = PR.runtime()["blocked"]
    sp = PR.stats_params()
    arms = tuple(arms)
    s1_arms = tuple(a for a in A.S1_ARMS if a in arms)
    if "EN-pred" in arms and "S1-obs" not in s1_arms:
        s1_arms = s1_arms + ("S1-obs",)
    en_arms = tuple(a for a in A.EN_ARMS if a in arms)
    s1_out_arms = [a for a in arms if a in A.S1_ARMS]
    out_dir = Path(out_dir or (P.PROJECT / rt["dataset_dir"]))
    stage_dir = Path(stage_dir or (P.PROJECT / rt["stage_dir"]))
    for d in (out_dir, stage_dir):
        shutil.rmtree(d, ignore_errors=True)
        d.mkdir(parents=True, exist_ok=True)
    chunks = pattern_chunks(data, list(data.genes if genes is None else genes), int(chunk or rt["chunk_genes"]))
    tm = {"n_chunks": len(chunks)}
    # ---------------- S1 ----------------------------------------------------------------------------------
    slim_s1, U_parts, n_s1 = [], [], []

    def on_s1(i, df):
        df = A.finish_s1_local(data, df)
        _write(df, stage_dir / f"s1-{i:05d}.parquet")
        n_s1.append(len(df))
        if len(df):
            slim_s1.append(_slim(df, ["arm", "gene", "feature", "p", "beta", "test_family", "collinear_with_design"], i))
            U_parts.append(df.loc[df["arm"] == s1_arms[0], ["gene", "feature"]])

    t0 = time.time()
    if s1_arms:
        A.run_s1(data, s1_arms, workers=workers_s1, chunks=chunks, on_chunk=on_s1, finish=False)
    tm["s1_window"] = (t0, time.time())
    log(f"S1: {sum(n_s1)} rows in {len(chunks)} chunks, {tm['s1_window'][1] - t0:.0f} s")
    U = pd.concat(U_parts, ignore_index=True) if U_parts else A.universe(data, genes)
    gf = A.gene_feats(data, U)
    del U_parts
    # ---------------- EN ----------------------------------------------------------------------------------
    slim_en, cand, lam_parts, n_en = [], [], [], {}

    def on_en(i, arm, tab, lt, tim):
        if len(tab):
            tab = A.finish_en_local(tab)
        _write(tab if len(tab) else pd.DataFrame({"arm": pd.Series(dtype=object)}), stage_dir / f"en-{arm}-{i:05d}.parquet")
        n_en[(arm, i)] = len(tab)
        if len(tab):
            slim_en.append(_slim(tab, ["arm", "gene", "feature", "selection_freq", "beta"], i))
            cand.append(tab[tab["selected"].astype(bool)])
        if len(lt):
            lam_parts.append(lt)

    walls = A.run_en(data, en_arms, U=U, gf=gf, workers=workers_en, chunks=chunks, on_chunk=on_en, finish=False,
                     max_patterns=int(rt["en_max_patterns"])) if en_arms else []
    tm["en_walls"] = walls
    log(f"EN: {sum(n_en.values())} rows")
    # ---------------- global steps on slim columns -------------------------------------------------------------
    t_g = time.time()
    S = pd.concat(slim_s1, ignore_index=True) if slim_s1 else pd.DataFrame(columns=["arm", "part"])
    E = pd.concat(slim_en, ignore_index=True) if slim_en else pd.DataFrame(columns=["arm", "part"])
    del slim_s1, slim_en
    for X in (S, E):
        for c in ("arm", "gene", "feature", "test_family"):
            if c in X and X[c].dtype != "category":
                X[c] = X[c].astype("category")
    gs = {}
    if len(S):
        S_arm, S_tf = S["arm"].astype(str).to_numpy(), S["test_family"].astype(str).to_numpy()
        p = S["p"].to_numpy(dtype=np.float64)
        q = np.full(len(S), np.nan)
        for key in pd.unique(pd.Series(S_arm + "|" + S_tf)):
            a_, t_ = key.split("|")
            m = (S_arm == a_) & (S_tf == t_)
            q[m] = ST.bh(p[m])
        coll = S["collinear_with_design"].astype(bool).to_numpy()
        ok = (q < sp["bh_q"]) & ~coll
        gs["q"], gs["selected"], gs["selected_F2"] = q, ok & (S_tf == "F1"), ok & (S_tf == "F2")
        fstr = S["feature"].astype(str).to_numpy().astype(object)
        rank = np.full(len(S), np.nan)
        beta = S["beta"].to_numpy(dtype=np.float64)
        for a_ in pd.unique(S_arm):
            idx = np.flatnonzero((S_arm == a_) & ~np.isnan(p))
            o = idx[A.s1_order(p[idx], beta[idx], fstr[idx])]
            rank[o] = np.arange(1, len(o) + 1)
        gs["rank"] = rank
        cut = {}
        for key in pd.unique(pd.Series(S_arm + "|" + S_tf)):
            a_, t_ = key.split("|")
            cut[key] = ST.bh_cutoff(p[(S_arm == a_) & (S_tf == t_)], sp["bh_q"])
    else:
        cut = {}
    ge = {}
    if len(E):
        E_arm = E["arm"].astype(str).to_numpy()
        erank = np.full(len(E), np.nan)
        for a_ in pd.unique(E_arm):
            idx = np.flatnonzero(E_arm == a_)
            ids = (E["gene"].astype(str).to_numpy()[idx].astype(object) + "|" + E["feature"].astype(str).to_numpy()[idx].astype(object))
            erank[idx] = EN.rank_en(E["selection_freq"].to_numpy()[idx], E["beta"].to_numpy()[idx], ids)
        ge["rank"] = erank
    # de-duplicated ranks per output arm
    vg_of = data.fmeta["variant_group"]

    def dedup(gene, feat, rank):
        rd = np.full(len(rank), np.nan)
        lost = np.full(len(rank), "", dtype=object)
        idx = np.flatnonzero(~np.isnan(rank))
        idx = idx[np.argsort(rank[idx], kind="stable")]
        vg = vg_of.reindex(feat[idx]).to_numpy()
        has = np.array([isinstance(v, str) and v != "" for v in vg], dtype=bool)
        key = np.where(has, gene[idx].astype(object) + "|" + vg.astype(str).astype(object),
                       gene[idx].astype(object) + "|" + feat[idx].astype(object))
        first = ~pd.Index(key).duplicated()
        rd[idx[first]] = np.arange(1, int(first.sum()) + 1)
        best = pd.Series(feat[idx][first], index=key[first])
        lost[idx[~first]] = best.reindex(key[~first]).to_numpy()
        return rd, lost

    for X, g_, rk in ((S, gs, gs.get("rank")), (E, ge, ge.get("rank"))):
        if not len(X):
            continue
        Xa = X["arm"].astype(str).to_numpy()
        rd = np.full(len(X), np.nan)
        lost = np.full(len(X), "", dtype=object)
        genes_ = X["gene"].astype(str).to_numpy().astype(object)
        feats_ = X["feature"].astype(str).to_numpy().astype(object)
        for a_ in pd.unique(Xa):
            if X is S and a_ not in s1_out_arms:
                continue
            m = np.flatnonzero(Xa == a_)
            r1, l1 = dedup(genes_[m], feats_[m], rk[m])
            rd[m], lost[m] = r1, l1
        g_["rank_dedup"], g_["lost_to_group_member"] = rd, lost
    # MD-strict / MD-global / EN-obs-strict on the candidates and their S1-obs rows
    out = {"p_cut": cut, "timings": tm}
    en_c = pd.concat(cand, ignore_index=True) if cand else pd.DataFrame(columns=["arm", "gene", "feature", "selected"])
    out["lam"] = pd.concat(lam_parts, ignore_index=True) if lam_parts else pd.DataFrame()
    if len(en_c) and len(S):
        part_of = S[["gene", "part"]].drop_duplicates("gene").set_index("gene")["part"].astype(int)
        need = sorted(set(part_of.reindex(en_c["gene"].astype(str).unique()).dropna().astype(int)))
        sel_keys = set(zip(en_c["gene"].astype(str), en_c["feature"].astype(str)))
        rows = []
        for i in need:
            st = pd.read_parquet(stage_dir / f"s1-{i:05d}.parquet")
            off = _offset(n_s1, i)
            for c in ("q", "selected", "selected_F2", "rank"):
                st[c] = gs[c][off:off + len(st)]
            st = st[st["arm"].astype(str) == "S1-obs"]
            k = list(zip(st["gene"].astype(str), st["feature"].astype(str)))
            rows.append(st[[x in sel_keys for x in k]])
        s1c = _decat(pd.concat(rows, ignore_index=True)) if rows else pd.DataFrame()
        en_c = _decat(en_c)
        if "EN-pred" in arms:
            md = A.strict(en_c, s1c, "EN-pred")
            out["md_strict_all"] = md
            out["md_strict"] = md[md["pass"]].reset_index(drop=True)
            out["p_cut_sc"] = A.p_cut_sc(md)
            mg = A.md_global(en_c, s1c, "EN-pred")
            out["md_global_all"] = mg
            out["md_global"] = mg[mg["pass"]].reset_index(drop=True)
            mds = set(zip(out["md_strict"]["gene"], out["md_strict"]["feature"]))
            out["md_global"]["in_md_strict"] = [(g, f) in mds for g, f in zip(out["md_global"]["gene"],
                                                                               out["md_global"]["feature"])]
        if "EN-obs" in arms:
            mo = A.strict(en_c, s1c, "EN-obs")
            out["en_obs_strict_all"] = mo
            out["en_obs_strict"] = mo[mo["pass"]].reset_index(drop=True)
    tm["global_s"] = time.time() - t_g
    # ---------------- final parts -------------------------------------------------------------------------------
    t_f = time.time()
    sha, n_rows, sel = {}, 0, {}
    for i in range(len(chunks)):
        s1p = pd.DataFrame()
        if n_s1 and n_s1[i]:
            s1p = pd.read_parquet(stage_dir / f"s1-{i:05d}.parquet")
            off = _offset(n_s1, i)
            for c in ("q", "selected", "selected_F2", "rank", "rank_dedup", "lost_to_group_member"):
                s1p[c] = gs[c][off:off + len(s1p)]
            s1p = s1p[s1p["arm"].astype(str).isin(s1_out_arms)]
        enp = []
        for arm in en_arms:
            if n_en.get((arm, i)):
                t = pd.read_parquet(stage_dir / f"en-{arm}-{i:05d}.parquet")
                off = _offset_en(n_en, en_arms, chunks, arm, i)
                for c in ("rank", "rank_dedup", "lost_to_group_member"):
                    t[c] = ge[c][off:off + len(t)]
                enp.append(t)
        enp = pd.concat(enp, ignore_index=True) if enp else pd.DataFrame()
        s1p, enp = _decat(s1p), _decat(enp)
        if not len(s1p) and not len(enp):
            continue
        a = A.annotate(data, s1p, enp)
        a["rank_dedup"] = a["rank_dedup"].astype(float)
        a["lost_to_group_member"] = a["lost_to_group_member"].fillna("").astype(object)
        f = out_dir / f"part-{i:05d}.parquet"
        _write(a, f)
        sha[f.name] = hashlib.sha256(f.read_bytes()).hexdigest()
        n_rows += len(a)
        for arm, v in a.groupby("arm", observed=True)["selected"].sum().items():
            sel[arm] = sel.get(arm, 0) + int(v)
    shutil.rmtree(stage_dir, ignore_errors=True)
    tm["final_s"] = time.time() - t_f
    tm["global_window"] = (t_g, time.time())
    out.update(parts=sorted(sha), sha256=sha, n_rows=n_rows, selected=sel, out_dir=str(out_dir))
    write_manifest(out_dir, data, arms, s1_out_arms + list(en_arms), out, tm)
    log(f"final: {n_rows} rows in {len(sha)} parts; global {tm['global_s']:.0f} s, write {tm['final_s']:.0f} s")
    return out


# ------------------------------------------------------------------------------------------------
# Completed-dataset manifest (reuse of a finished streamed dataset, e.g. run_xlib.py --reuse-arms)
# ------------------------------------------------------------------------------------------------

MANIFEST = "_manifest.json"


def _ids_sha(ids) -> str:
    return hashlib.sha256("\n".join(sorted(map(str, ids))).encode()).hexdigest()


def write_manifest(out_dir: Path, data, arms, out_arms, out: dict, tm: dict, extra: dict | None = None) -> Path:
    """``<out_dir>/_manifest.json``, written last: the dataset is complete iff it exists and matches the parts."""
    import datetime as dt
    import json
    rec = {"complete": True, "utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
           "arms": list(arms), "out_arms": list(out_arms), "parts": dict(sorted(out["sha256"].items())),
           "n_rows": int(out["n_rows"]), "n_chunks": int(tm.get("n_chunks", 0)), "selected": out.get("selected", {}),
           "genes_sha256": _ids_sha(data.genes), "lines_sha256": _ids_sha(data.lines), "n_genes": len(data.genes),
           "n_lines": len(data.lines),
           "timings": {k: tm[k] for k in ("s1_window", "global_window", "en_walls") if k in tm}, **(extra or {})}
    f = Path(out_dir) / MANIFEST
    f.write_text(json.dumps(rec, indent=1, default=float), encoding="utf-8")
    return f


def verify_dataset(out_dir: Path, arms, genes, lines=None, stage_dir: Path | None = None,
                   expect_parts: int | None = None, expect_rows: int | None = None) -> dict:
    """Check that ``out_dir`` holds a complete streamed dataset for ``arms`` / ``genes`` (/ ``lines``). With a
    manifest: same part set, every part's sha256, arms, gene and line sets. Without one (datasets written before the
    manifest existed): the stage dir must be gone (removed only after the last part), and the part count and total
    row count must equal ``expect_parts`` / ``expect_rows`` (from the run's log), the arm values must be the output
    arms and the genes on disk a subset of ``genes`` (genes without universe rows have none). Returns the manifest (legacy: a reconstructed one); raises ValueError."""
    import json

    import pyarrow.parquet as pq
    out_dir = Path(out_dir)
    parts = sorted(p.name for p in out_dir.glob("part-*.parquet"))
    if not parts:
        raise ValueError(f"no parts under {out_dir}")
    if stage_dir is not None and Path(stage_dir).exists():
        raise ValueError(f"stage dir {stage_dir} still exists: the streamed run did not finish")
    f = out_dir / MANIFEST
    if f.exists():
        m = json.loads(f.read_text(encoding="utf-8"))
        if not m.get("complete") or sorted(m["parts"]) != parts:
            raise ValueError("parts on disk differ from the manifest")
        bad = [p for p in parts if hashlib.sha256((out_dir / p).read_bytes()).hexdigest() != m["parts"][p]]
        if bad:
            raise ValueError(f"sha256 mismatch: {bad}")
        if list(m["arms"]) != list(arms):
            raise ValueError(f"dataset arms {m['arms']} != {list(arms)}")
        if m["genes_sha256"] != _ids_sha(genes) or (lines is not None and m["lines_sha256"] != _ids_sha(lines)):
            raise ValueError("dataset genes / lines differ from this run")
        return m
    if expect_parts is None or expect_rows is None:
        raise ValueError(f"{out_dir} has no {MANIFEST}: pass the part and row counts of the run's log")
    n = sum(pq.ParquetFile(out_dir / p).metadata.num_rows for p in parts)
    if len(parts) != int(expect_parts) or n != int(expect_rows):
        raise ValueError(f"{len(parts)} parts / {n} rows on disk, expected {expect_parts} / {expect_rows}")
    seen_arms, seen_genes = set(), set()
    for p in parts:
        t = pq.read_table(out_dir / p, columns=["arm", "gene"])
        seen_arms |= set(map(str, t.column("arm").unique().to_pylist()))
        seen_genes |= set(map(str, t.column("gene").unique().to_pylist()))
    want = [a for a in A.S1_ARMS if a in arms] + [a for a in A.EN_ARMS if a in arms]
    if seen_arms != set(want):
        raise ValueError(f"arms on disk {sorted(seen_arms)} != {want}")
    extra_g = seen_genes - set(map(str, genes))
    if extra_g:
        raise ValueError(f"{len(extra_g)} genes on disk are not genes of this run")
    sha = {p: hashlib.sha256((out_dir / p).read_bytes()).hexdigest() for p in parts}
    return {"complete": True, "arms": list(arms), "parts": sha, "n_rows": int(n), "genes_sha256": _ids_sha(genes),
            "lines_sha256": None if lines is None else _ids_sha(lines), "n_genes_with_rows": len(seen_genes),
            "n_genes": len(set(map(str, genes))), "verified_from": "legacy part / row counts (lines not verifiable)"}


def _offset(counts: list, i: int) -> int:
    return int(sum(counts[:i]))


def _offset_en(n_en: dict, en_arms, chunks, arm: str, i: int) -> int:
    off = 0
    for a in en_arms:
        for j in range(len(chunks)):
            if a == arm and j == i:
                return off
            off += n_en.get((a, j), 0)
    return off


def _decat(df: pd.DataFrame) -> pd.DataFrame:
    if not len(df):
        return df
    df = df.copy()
    for c in df.columns:
        if isinstance(df[c].dtype, pd.CategoricalDtype):
            df[c] = df[c].astype(object)
    return df


def read_dataset(out_dir, columns=None) -> pd.DataFrame:
    parts = sorted(Path(out_dir).glob("part-*.parquet"))
    return pd.concat([_decat(pd.read_parquet(p, columns=columns)) for p in parts], ignore_index=True) if parts \
        else pd.DataFrame()
