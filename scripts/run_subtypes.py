"""Per-gene subtypes (DESIGN v3.2 §8.1) and biomarker modules (§8.2) on the catalogue (amendment A7: MD-global;
MD-strict reported alongside).

    python -B scripts/run_subtypes.py              # real D (after run_arms.py)
    python -B scripts/run_subtypes.py --toy DIR    # toy workspace (engine unchanged; never the project results/)

Inputs: results/md/md_global.csv, md_strict.csv; cache/benchmark/arms/ via endpoint_io.read_arms (EN-pred selected rows: selected, full-D
coefficient beta_full); D effects (guard.load_effects scope dev), fold-centred m5 OOF predictions, 𝒟 and covariates
(pipeline.Data.from_real).
S(g, d): accurate-dependent lines of g whose line share s_ld = β_gd(x̃_ld − x̄_d) / Σ_{selected d'} β_gd'(x̃_ld' − x̄_d')
(full-D EN-pred fit, both stages, FWL-residualised 𝒟; lines with a negative sum only) is >= line_share_min; binary
d also needs carrier status. M(d): genes with (g, d, same sign) in the list; member lines = carriers
accurate-dependent on >= module_line_frac_min of them; size as genes and effective gene clusters (average linkage on
1 - |r| of D effects, cut at gene_cluster_abs_r).
Writes results/subtypes/{subtypes.csv, subtype_lines.csv, subtype_genes.csv, modules.csv, module_lines.csv,
summary.json}; timing stems subtypes, modules. Refuses unless HEAD descends from freeze-bounds.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from biomark import analyses as AN  # noqa: E402
from biomark import m5 as M5  # noqa: E402

SCRIPT = "scripts/run_subtypes.py"
OUT = "results/subtypes"


def run(ws, src) -> dict:
    t0 = AN.utc_now()
    cat = M5.catalogue(ws)
    genes = sorted(cat["gene"].unique())
    data = src.data(genes)
    en = ws.arms(arms=["EN-pred"], columns=["gene", "feature", "selected", "beta_full", "stage"], genes=genes,
                 selected_only=True)
    dep, acc = M5.dependency_status(data)
    sym = src.symbols()
    sub = M5.subtypes(data, en, cat, dep, acc, sym)
    outs = [ws.write_csv(sub["summary"], f"{OUT}/subtypes.csv"), ws.write_csv(sub["lines"], f"{OUT}/subtype_lines.csv"),
            ws.write_csv(sub["genes"], f"{OUT}/subtype_genes.csv")]
    s = sub["summary"]
    ws.timing("subtypes", {"n_pairs": int(len(s)), "n_genes": len(genes)}, t0, workers=1)
    t1 = AN.utc_now()
    mods = M5.modules(data, cat, acc, sym)
    outs += [ws.write_csv(mods["modules"], f"{OUT}/modules.csv"), ws.write_csv(mods["lines"], f"{OUT}/module_lines.csv")]
    m = mods["modules"]
    summ = {"catalogue": "MD-global (A7); MD-strict alongside", "n_pairs": int(len(s)),
            "n_pairs_with_members": int((s["n_members"] > 0).sum()) if len(s) else 0,
            "by_list": {n: {"n_pairs": int(s["lists"].str.contains(n).sum()),
                            "n_with_members": int((s["lists"].str.contains(n) & (s["n_members"] > 0)).sum()),
                            "n_modules": int((m["list"] == n).sum()) if len(m) else 0,
                            "n_multi_gene_modules": int(((m["list"] == n) & (m["n_genes"] > 1)).sum()) if len(m) else 0}
                        for n in M5.lists()},
            "beta_full_zero": int(s["note"].astype(str).str.startswith("full-D").sum()) if len(s) else 0,
            "genes": {"n": int(len(sub["genes"])),
                      "unexplained_lines": int(sub["genes"]["n_unexplained"].sum()) if len(sub["genes"]) else 0,
                      "accurate_unassigned_lines": int(sub["genes"]["n_accurate_unassigned"].sum()) if len(sub["genes"]) else 0}}
    outs.append(ws.write_json(summ, f"{OUT}/summary.json"))
    ws.timing("modules", {"n_modules": int(len(m))}, t1, workers=1)
    ws.register(SCRIPT, outs)
    return {"subtypes": sub, "modules": mods, "summary": summ}


def main(argv=None) -> int:
    a = M5.parse_args(argv, __doc__)
    ws, src = M5.open_workspace(a)
    r = run(ws, src)
    print(f"subtypes: {r['summary']['n_pairs']} catalogue pairs, {r['summary']['n_pairs_with_members']} with members; "
          f"modules: {len(r['modules']['modules'])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
