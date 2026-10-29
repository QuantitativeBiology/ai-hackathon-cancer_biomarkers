"""Positive-control gate attempt (DESIGN v3.2 §9.3).

    python -B scripts/run_gate.py                                   # attempt 1
    python -B scripts/run_gate.py --reason "<bug fix>" --toy-test tests/test_x.py::test_y   # attempt n > 1

Denominator: gate-eligible rows of results/controls/pc_eligibility.csv (frozen at freeze-bounds). Recovered = (g,
primary or secondary feature) in the list with the expected sign; decided on MD-strict, MD-global and every arm
reported (selected rows). Each control is read from the universe of its gene (DESIGN §2): eligible genes from the
eligible lists (results/md/, cache/benchmark/arms/), non-eligible control genes from U_ctrl (configs/u_ctrl_genes.txt;
results/md_ctrl/, cache/benchmark/arms_ctrl/, written by run_arms.py --universe ctrl with BH within U_ctrl). The
universe is recorded per row and the input hashes of both datasets are logged. Refuses if a gate-eligible control gene
is in neither universe, or if U_ctrl is needed and its outputs are missing or were built from another
u_ctrl_genes.txt. Writes results/gates/G_pc_attempt{n}.json (never overwritten) and appends a row to
results/gate_iterations.csv. Refuses unless HEAD descends from tag freeze-bounds.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd  # noqa: E402

from biomark import arms as A  # noqa: E402
from biomark import endpoint_io as IO  # noqa: E402
from biomark import endpoints as EP  # noqa: E402
from biomark import gitutil  # noqa: E402
from biomark import paths as P  # noqa: E402

TAG = "freeze-bounds"
GATES = P.RESULTS / "gates"
ITER = P.RESULTS / "gate_iterations.csv"


def sha(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def next_attempt() -> int:
    GATES.mkdir(parents=True, exist_ok=True)
    have = [int(f.stem.replace("G_pc_attempt", "")) for f in GATES.glob("G_pc_attempt*.json")]
    return max(have, default=0) + 1


def arms_dataset_sha(d: Path | None = None) -> str:
    """sha256 over the sorted (name, sha256) of the blocked arms parts in ``d`` (default: the eligible dataset; legacy
    single file cache/benchmark/arms.parquet if there are no parts)."""
    d = Path(d or IO.arms_dir("eligible"))
    parts = sorted(d.glob("part-*.parquet"))
    if not parts:
        return sha(d.parent / f"{d.name}.parquet")
    return hashlib.sha256("\n".join(f"{p.name} {sha(p)}" for p in parts).encode()).hexdigest()


def gate_lists(genes_by_universe: dict) -> dict:
    """{list name: DataFrame(gene, feature, beta, universe)} for MD-strict, MD-global and every arm, each control
    gene's rows read only from its own universe ({'eligible': [...], 'ctrl': [...]})."""
    names = ["MD-strict", "MD-global", *A.ARMS]
    out = {k: [] for k in names}
    for u, genes in genes_by_universe.items():
        genes = sorted(set(map(str, genes)))
        if not genes:
            continue
        md = {"MD-strict": IO.md_strict(u), "MD-global": IO.md_global(u)}
        for k, L in md.items():
            L = L[L["gene"].astype(str).isin(genes)]
            out[k].append(L[["gene", "feature", "beta"]].assign(universe=u))
        for arm in A.ARMS:
            L = IO.read_arms(arms=[arm], columns=["gene", "feature", "beta"], genes=genes, selected_only=True,
                             dataset_dir=IO.arms_dir(u))
            out[arm].append(L[["gene", "feature", "beta"]].assign(universe=u))
    cols = ["gene", "feature", "beta", "universe"]
    return {k: (pd.concat(v, ignore_index=True) if v else pd.DataFrame(columns=cols)) for k, v in out.items()}


def ctrl_inputs_ok() -> None:
    """U_ctrl outputs exist and were built from the current configs/u_ctrl_genes.txt."""
    need = [IO.md_dir("ctrl") / "md_strict.csv", IO.md_dir("ctrl") / "md_global.csv",
            IO.benchmark_dir("ctrl") / "summary.json"]
    miss = [str(p) for p in need if not p.exists()]
    if not any(IO.arms_dir("ctrl").glob("part-*.parquet")):
        miss.append(str(IO.arms_dir("ctrl")))
    if miss:
        raise EP.EndpointRefusal(f"U_ctrl outputs missing (run scripts/run_arms.py --universe ctrl): {miss}")
    s = IO.read_json(IO.benchmark_dir("ctrl") / "summary.json")
    if s.get("u_ctrl_genes_sha256") != sha(P.CONFIGS / "u_ctrl_genes.txt"):
        raise EP.EndpointRefusal("results/benchmark_ctrl/summary.json was built from another u_ctrl_genes.txt")


def inputs_sha256(with_ctrl: bool) -> dict:
    h = {"pc_eligibility.csv": sha(P.RESULTS / "controls" / "pc_eligibility.csv"),
         "positive_controls.csv": sha(P.CONFIGS / "positive_controls.csv"),
         "eligible_genes.csv": sha(P.CONFIGS / "eligible_genes.csv"),
         "md_strict.csv": sha(IO.md_dir("eligible") / "md_strict.csv"),
         "md_global.csv": sha(IO.md_dir("eligible") / "md_global.csv"),
         "arms_dataset": arms_dataset_sha(IO.arms_dir("eligible"))}
    if with_ctrl:
        h.update({"u_ctrl_genes.txt": sha(P.CONFIGS / "u_ctrl_genes.txt"),
                  "md_ctrl/md_strict.csv": sha(IO.md_dir("ctrl") / "md_strict.csv"),
                  "md_ctrl/md_global.csv": sha(IO.md_dir("ctrl") / "md_global.csv"),
                  "arms_dataset_ctrl": arms_dataset_sha(IO.arms_dir("ctrl"))})
    return h


def run(reason: str | None, toy_test: str | None, check=None) -> dict:
    from biomark import analyses as AN
    EP.require_descends(TAG, check=check)
    t_start = AN.utc_now()
    n = next_attempt()
    if n > 1 and not (reason and toy_test):
        raise EP.EndpointRefusal("attempt > 1 needs --reason and --toy-test (a failing-then-passing toy test)")
    pc = pd.read_csv(P.RESULTS / "controls" / "pc_eligibility.csv", dtype=str)
    ctl = pd.read_csv(P.CONFIGS / "positive_controls.csv", dtype=str)
    e = pc[pc["gate_eligible"].astype(str).str.lower().isin(["true", "1"])]
    gene_of = ctl.set_index("control_id")["dependency_entrez"].astype(str)
    genes = sorted(set(gene_of.reindex(e["control_id"])))
    uof = IO.universe_of(genes)
    lost = sorted(g for g, u in uof.items() if u == "none")
    if lost:
        raise EP.EndpointRefusal(f"gate-eligible control genes in neither the eligible set nor U_ctrl: {lost}")
    by_u = {u: [g for g in genes if uof[g] == u] for u in IO.UNIVERSES}
    if by_u["ctrl"]:
        ctrl_inputs_ok()
    lists = gate_lists(by_u)             # decided on MD-strict (A7); MD-global and every arm reported
    res = EP.gate_attempt(pc, ctl, lists, universe_of=uof)
    commit = gitutil.head(P.PROJECT)
    out = {"attempt": n, "utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"), "commit": commit,
           "reason": reason or "attempt 1", "summary": res["summary"],
           "rows": res["rows"].to_dict("records"),
           "universe_genes": by_u,
           "inputs_sha256": inputs_sha256(bool(by_u["ctrl"]))}
    f = GATES / f"G_pc_attempt{n}.json"
    if f.exists():
        raise EP.EndpointRefusal(f"{f.name} exists and is never overwritten")
    IO.write_json(f, out)
    prev = None
    if n > 1:
        prev = IO.read_json(GATES / f"G_pc_attempt{n - 1}.json")["summary"]
    row = {"attempt": n, "utc": out["utc"], "commit": commit, "reason": out["reason"], "toy_test": toy_test or "",
           "gate_before": "" if prev is None else f"{prev['pass']} ({prev['overall']['k']}/{prev['overall']['n']})",
           "gate_after": f"{res['summary']['pass']} ({res['summary']['overall']['k']}/{res['summary']['overall']['n']})",
           "attempt_file_sha256": sha(f)}
    it = pd.read_csv(ITER) if ITER.exists() else pd.DataFrame(columns=list(row))
    it = pd.concat([it, pd.DataFrame([row])], ignore_index=True)
    it.to_csv(ITER, index=False)
    AN.write_timing("gate_pc", {"attempt": n}, start_utc=t_start, workers=1)
    AN.append_registry("scripts/run_gate.py", "ok", [f"results/gates/{f.name}", "results/gate_iterations.csv"])
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reason", default=None)
    ap.add_argument("--toy-test", default=None)
    a = ap.parse_args(argv)
    out = run(a.reason, a.toy_test)
    s = out["summary"]
    print(f"gate attempt {out['attempt']}: pass={s['pass']} overall {s['overall']['k']}/{s['overall']['n']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
