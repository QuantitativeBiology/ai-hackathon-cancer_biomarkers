"""Confounder audit (DESIGN v3.2 §5.5 confounding delta; §5.3 VIF / condition number).

    python -B scripts/run_confounder_audit.py [--workers N]

Pairs = MD-global ∪ MD-strict ∪ S1-obs BH hits ∪ gate-scope positive controls (in U), plus the gate-scope controls
whose gene is in U_ctrl (DESIGN §2; configs/u_ctrl_genes.txt), taken where they are rows of U_ctrl (the ctrl dataset
cache/benchmark/arms_ctrl/ of run_arms.py --universe ctrl) and marked universe = ctrl in delta.csv / delta_flags.csv
(eligible pairs: universe = eligible, unchanged). For each pair: drop-one-block refits of the
§5.3 model (lineage, library, culture, growth, burden, TP53, own CN) and the OncotreePrimaryDisease sensitivity block
(-> subtype-confounded); flag if |β_full − β_minus| > 0.5·|β_minus| and > SE_minus (audit.confounding_delta, unchanged).
Writes results/confounders/delta.csv (+ delta_flags.csv: one row per pair), results/covariates/vif_check.json (max VIF
and condition number of the D design as the S1 scan codes it, vs confounders_s5.vif_max); results/covariates/vif.csv
is written only if absent (the covariate builder owns it). Timing stems confounder_delta and vif.
Refuses unless HEAD descends from freeze-bounds.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd  # noqa: E402

from biomark import analyses as AN  # noqa: E402
from biomark import audit as AU  # noqa: E402
from biomark import endpoint_io as IO  # noqa: E402
from biomark import endpoints as EP  # noqa: E402
from biomark import params as PR  # noqa: E402
from biomark import paths as P  # noqa: E402
from biomark import pipeline as PL  # noqa: E402

_W = {}


def _init(path):
    from biomark import workerstore as WS
    _W["data"] = WS.load(path)["data"]


def _chunk(pairs):
    return AU.confounding_delta(_W["data"], pairs)


def _in_u(p: pd.DataFrame, universe: str) -> pd.DataFrame:
    u = IO.read_arms(arms=["S1-obs"], columns=["gene", "feature"], genes=sorted(set(p["gene"])),
                     dataset_dir=IO.arms_dir(universe))
    U = set(zip(u["gene"].astype(str), u["feature"]))
    p = p[[k in U for k in zip(p["gene"], p["feature"])]]
    return p.groupby(["gene", "feature"], as_index=False)["source"].agg(lambda s: "|".join(sorted(set(s))))


def _ctrl_ready() -> bool:
    """U_ctrl outputs exist and were built from the current configs/u_ctrl_genes.txt (as run_gate.py checks)."""
    import hashlib
    f = IO.benchmark_dir("ctrl") / "summary.json"
    if not (any(IO.arms_dir("ctrl").glob("part-*.parquet")) and f.exists()):
        return False
    got = IO.read_json(f).get("u_ctrl_genes_sha256")
    if got != hashlib.sha256((P.CONFIGS / "u_ctrl_genes.txt").read_bytes()).hexdigest():
        raise SystemExit("results/benchmark_ctrl/summary.json was built from another u_ctrl_genes.txt")
    return True


def audit_pairs(include_ctrl: bool = True) -> pd.DataFrame:
    """Audit pairs with a ``universe`` column: the eligible pairs (catalogue, MD-strict, S1-obs BH, gate controls in
    U), then the gate controls of U_ctrl genes in U_ctrl (``include_ctrl``; refuses if U_ctrl genes exist but the
    ctrl run is missing)."""
    md, mg = IO.md_strict(), IO.catalogue()                            # A7: MD-global (catalogue) and MD-strict
    s1h = IO.read_arms(arms=["S1-obs"], columns=["gene", "feature"], selected_only=True,
                       dataset_dir=IO.arms_dir("eligible"))
    ctl = pd.read_csv(P.CONFIGS / "positive_controls.csv", dtype=str)
    ctl = ctl[ctl["scope"] == "gate"]
    rows = [(g, f, "MD-global") for g, f in zip(mg["gene"], mg["feature"])]
    rows += [(g, f, "MD-strict") for g, f in zip(md["gene"], md["feature"])]
    rows += [(g, f, "S1-obs BH") for g, f in zip(s1h["gene"], s1h["feature"])]
    crow = []
    for _, c in ctl.iterrows():
        crow += [(str(c["dependency_entrez"]), f, "control") for f in EP.control_features(c)]
    cols = ["gene", "feature", "source"]
    p = _in_u(pd.DataFrame(rows + crow, columns=cols), "eligible").assign(universe="eligible")
    cg = IO.universe_genes()["ctrl"] if include_ctrl else set()
    pc = pd.DataFrame([r for r in crow if r[0] in cg], columns=cols)
    if not len(pc):
        return p
    if not _ctrl_ready():
        raise SystemExit("gate controls in U_ctrl but no U_ctrl run: run scripts/run_arms.py --universe ctrl first")
    pc = _in_u(pc, "ctrl").assign(universe="ctrl")
    return pd.concat([p, pc], ignore_index=True) if len(pc) else p


def delta_tables(data, pairs: pd.DataFrame, w: int = 1) -> tuple:
    """(delta, delta_flags) for ``pairs`` (audit.confounding_delta per pair chunk; pair columns incl. universe merged)."""
    chunks = [pairs.iloc[i:i + PR.null_params()["chunk_genes"]] for i in range(0, len(pairs), PR.null_params()["chunk_genes"])]
    if w <= 1:
        parts = [AU.confounding_delta(data, c) for c in chunks]
    else:
        from concurrent.futures import ProcessPoolExecutor

        from biomark import workerstore as WS
        with WS.staged({"data": data.for_workers(include_obs=True)}, "delta") as path, \
                ProcessPoolExecutor(max_workers=w, initializer=_init, initargs=(path,)) as ex:
            parts = list(ex.map(_chunk, chunks))
    delta = pd.concat(parts, ignore_index=True).merge(pairs, on=["gene", "feature"], how="left")
    sens = delta["block"].str.startswith("primary_disease")
    flags = delta[~sens].groupby(["gene", "feature"])["flagged"].any().rename("delta_flagged").to_frame()
    flags["subtype_confounded"] = delta[sens].set_index(["gene", "feature"])["flagged"].reindex(flags.index).fillna(False)
    flags["blocks_flagged"] = delta[~sens & delta["flagged"]].groupby(["gene", "feature"])["block"].agg("|".join).reindex(
        flags.index).fillna("")
    return delta, flags.reset_index().merge(pairs, on=["gene", "feature"])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=None)
    a = ap.parse_args(argv)
    EP.require_descends("freeze-bounds")
    w = a.workers or PR.null_params()["fg_workers"]
    t0 = time.perf_counter()
    pairs = audit_pairs()
    data = PL.Data.from_real(genes=sorted(pairs["gene"].unique()))
    data.target("obs")
    delta, flags = delta_tables(data, pairs, w)
    out = P.RESULTS / "confounders"
    out.mkdir(parents=True, exist_ok=True)
    delta.to_csv(out / "delta.csv", index=False)
    flags.to_csv(out / "delta_flags.csv", index=False)
    AN.write_timing("confounder_delta", {"n_pairs": int(len(pairs)),
                                         "n_pairs_u_ctrl": int((pairs["universe"] == "ctrl").sum())},
                    wall_s=time.perf_counter() - t0, workers=w)
    # VIF / condition number of the D design
    t1 = time.perf_counter()
    v, cond, aliases = AN.design_vif(data.cov)
    vmax = float(PR.dp()["confounders_s5"]["vif_max"])
    cdir = P.RESULTS / "covariates"
    cdir.mkdir(parents=True, exist_ok=True)
    if not (cdir / "vif.csv").exists():
        v.rename("vif").rename_axis("column").reset_index().to_csv(cdir / "vif.csv", index=False)
    IO.write_json(cdir / "vif_check.json", {"vif_max_observed": float(v.max()), "vif_max_column": str(v.idxmax()),
                                            "vif_max_allowed": vmax, "pass": bool(v.max() <= vmax),
                                            "condition_number": cond, "n_columns": int(len(v)),
                                            "aliased_columns": aliases,
                                            "top": {k: float(x) for k, x in v.head(10).items()}})
    AN.write_timing("vif", {}, wall_s=time.perf_counter() - t1, workers=1)
    AN.append_registry("scripts/run_confounder_audit.py", "ok",
                       ["results/confounders/delta.csv", "results/confounders/delta_flags.csv",
                        "results/covariates/vif_check.json"])
    print(f"pairs {len(pairs)}; flagged {int(flags['delta_flagged'].sum())}; VIF max {v.max():.2f} (<= {vmax}: {v.max() <= vmax})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
