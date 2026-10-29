"""scripts/run_arms.py --universe ctrl writes only ctrl-specific paths (cache/benchmark/arms_ctrl/, results/benchmark_ctrl/,
results/md_ctrl/; results/bulk_manifest.csv rows keyed by ctrl paths), so it can never overwrite an eligible output; the
U_ctrl gene list is disjoint from the eligible set. Toy data in a temporary project root; no real data is read."""
import importlib.util
import time
from pathlib import Path

import pandas as pd
import pytest

from biomark import arms as A
from biomark import paths as P
from biomark import pipeline as PL
from biomark import toysim as TS

ROOT = Path(__file__).resolve().parents[1]


def _script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def toy_out():
    toy = TS.make_toy()
    el = PL.toy_eligibility(toy)
    genes = sorted(el.index[el["eligible"]].astype(str))[:3]
    data = PL.Data.from_toy(toy, genes=genes)
    return data, A.run_arms(data, arms=A.ARMS, workers=1)


def _files(root: Path) -> set:
    return {p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file()}


def _write(ra, data, out, universe):
    paths = ra.output_paths(universe)
    ds = paths["dataset"]
    ds.mkdir(parents=True, exist_ok=True)
    out["arms"].to_parquet(ds / "part-00000.parquet")
    now = time.time()
    res = dict(out, parts=["part-00000.parquet"], n_rows=len(out["arms"]), selected={},
               timings={"s1_window": (now - 2, now - 1), "global_window": (now - 1, now),
                        "en_walls": [{"arm": a, "workers": 1, "start_epoch": now - 2, "end_epoch": now - 1}
                                     for a in A.EN_ARMS]})
    return ra.write_outputs(res, data, universe, {"s1": {}, "en": {}}, 1)


def test_ctrl_outputs_are_ctrl_specific(toy_out, tmp_path, monkeypatch):
    data, out = toy_out
    ra = _script("run_arms")
    monkeypatch.setattr(P, "PROJECT", tmp_path)
    monkeypatch.setattr(P, "RESULTS", tmp_path / "results")
    monkeypatch.setattr(ra.gitutil, "head", lambda repo=None: "test")
    pe, pc = ra.output_paths("eligible"), ra.output_paths("ctrl")
    for k in ("dataset", "stage", "benchmark", "md"):
        assert pe[k] != pc[k] and "ctrl" in pc[k].name, k
    assert pc["timing"].is_relative_to(pc["benchmark"])
    _write(ra, data, out, "eligible")
    elig = _files(tmp_path)
    man_e = pd.read_csv(tmp_path / "results" / "bulk_manifest.csv")
    _write(ra, data, out, "ctrl")
    both = _files(tmp_path)
    ctrl = both - elig
    allowed = ("cache/benchmark/arms_ctrl/", "results/benchmark_ctrl/", "results/md_ctrl/")
    assert ctrl and all(f.startswith(allowed) for f in ctrl), sorted(ctrl)
    assert elig & {f for f in both if f.startswith(allowed)} == set()
    man = pd.read_csv(tmp_path / "results" / "bulk_manifest.csv")
    assert set(man_e["path"]) <= set(man["path"])                    # eligible rows kept
    new = set(man["path"]) - set(man_e["path"])
    assert new and all(p.startswith("cache/benchmark/arms_ctrl/") for p in new)
    assert {"s1_obs.json", "en_pred.json", "md_strict.json"} <= {Path(f).name for f in ctrl
                                                                  if f.startswith("results/benchmark_ctrl/timing/")}


def test_u_ctrl_list_disjoint_from_eligible_and_built_from_configs():
    bu = _script("build_u_ctrl_genes")
    t = bu.u_ctrl_table()
    el = pd.read_csv(P.CONFIGS / "eligible_genes.csv", dtype={"entrez": str})
    assert not set(t["entrez"]) & set(el.loc[el["eligible"].astype(bool), "entrez"])
    ctl = pd.read_csv(P.CONFIGS / "positive_controls.csv", dtype=str)
    nonelig = set(ctl["dependency_entrez"]) - set(el.loc[el["eligible"].astype(bool), "entrez"])
    assert nonelig <= set(t["entrez"])
    f = P.CONFIGS / "u_ctrl_genes.txt"
    if f.exists():
        assert bu.read_u_ctrl(f) == t["entrez"].tolist(), "configs/u_ctrl_genes.txt is stale"
