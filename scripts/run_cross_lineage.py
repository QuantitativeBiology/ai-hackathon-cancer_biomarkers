"""Within vs across lineage (DESIGN v3.2 §8.3) on MD-global ∪ MD-strict ∪ positive controls ∪ the S1-obs BH top 200
(amendment A7: MD-global is the catalogue, MD-strict reported alongside).

    python -B scripts/run_cross_lineage.py [--workers N]          # real D
    python -B scripts/run_cross_lineage.py --toy DIR --workers 1  # toy workspace

Per pair (audit.cross_lineage_stats, unchanged rules): testable lineages (>= lineage_min_lines D lines, >= 3
carriers and >= 3 non-carriers), per-lineage power under the §6.1 rule (EB-shrunk β from one fit over every S1-obs
U row), pooled confounder-model test on testable lineages, LOLO sign, organ-system groups, d x lineage interaction.
audit.cross_lineage_labels then applies BH over F4 (pooled) and F5 (per lineage) across ALL pairs of this run and
the labels in precedence order: cross-lineage > lineage-restricted > lineage-unresolvable > undetermined.
Known limitation (A7-iii): a one-lineage effect can be labelled cross-lineage when the null lineages share its sign
by chance; ``single_lineage_driven`` flags cross-lineage labels with at most one powered lineage that is nominally
significant (p < cross_lineage.single_lineage_p_max) with the pooled sign.
Writes results/lineage/{cross_lineage.csv, cross_lineage_per_lineage.csv, summary.json}; timing stem cross_lineage.
Refuses unless HEAD descends from freeze-bounds.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from biomark import analyses as AN  # noqa: E402
from biomark import m5 as M5  # noqa: E402
from biomark import params as PR  # noqa: E402

SCRIPT = "scripts/run_cross_lineage.py"
OUT = "results/lineage"
LIMITATION = ("A7-iii: the cross-lineage rule can mislabel a one-lineage effect when the null lineages share its "
              "sign by chance; single_lineage_driven marks cross-lineage labels with <= 1 nominally significant "
              "powered lineage of the pooled sign")


def run(ws, src, workers: int) -> dict:
    t0 = AN.utc_now()
    cat, _, pairs = M5.analysis_pairs(ws, src, top_s1=int(PR.pp()["audit"]["cross_lineage_top_s1"]))
    bsh = M5.shrunk_beta_ws(ws, keep=set(zip(pairs["gene"], pairs["feature"])))
    data = src.data(pairs["gene"].unique())
    lab, per = M5.cross_lineage(data, pairs, bsh, M5.sign_map(cat), workers=workers)
    outs = [ws.write_csv(lab, f"{OUT}/cross_lineage.csv"), ws.write_csv(per, f"{OUT}/cross_lineage_per_lineage.csv")]
    base = lab["label"].astype(str).str.split(" \\(").str[0]
    by = {}
    for src_name in ["MD-global", "MD-strict", "control", "S1-obs BH"]:
        m = lab["source"].astype(str).str.split("|").map(lambda v, s=src_name: s in v)
        by[src_name] = {"n": int(m.sum()), "labels": base[m].value_counts().to_dict(),
                        "single_lineage_driven": int((lab["single_lineage_driven"] & m).sum())}
    summ = {"n_pairs": int(len(lab)), "labels": base.value_counts().to_dict(), "by_source": by,
            "bh_families": "F4 pooled and F5 per-lineage over every pair of this run", "limitation": LIMITATION,
            "workers": workers}
    outs.append(ws.write_json(summ, f"{OUT}/summary.json"))
    ws.timing("cross_lineage", {"n_pairs": int(len(lab))}, t0, workers=workers)
    ws.register(SCRIPT, outs)
    return {"labels": lab, "per_lineage": per, "summary": summ}


def main(argv=None) -> int:
    a = M5.parse_args(argv, __doc__, workers=True)
    ws, src = M5.open_workspace(a)
    r = run(ws, src, a.workers)
    print(f"cross-lineage: {r['summary']['n_pairs']} pairs; {r['summary']['labels']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
