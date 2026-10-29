"""M5 audits (DESIGN v3.2 §8.4 co-location, §8.5 CN-cutting artefact, §8.6 culture, §8.7 hypermutation / TP53) on
the catalogue (amendment A7: MD-global, MD-strict alongside).

    python -B scripts/run_m5_audits.py [--workers N]          # real D (after run_subtypes.py)
    python -B scripts/run_m5_audits.py --toy DIR --workers 1  # toy workspace

- Co-location (audit.locus_resolution, unchanged): every CN/SIL determinant in a locus group among MD-global ∪
  MD-strict ∪ S1-obs BH ∪ positive controls, refitted per neighbour on lines without concordant carriers; >= 5
  discordant carriers, same sign, p < 0.05 for every neighbour => gene level, else locus level.
- Cutting (audit.amplicon_test, unchanged): AMP determinants whose group contains g -> cutting artefact /
  amplified-driver addiction / unresolved.
- Culture: Cramér's V of each subtype (results/subtypes/subtypes.csv, subtype_lines.csv) with CellFormat, 3D,
  serum-free, medium family and library; V >= culture_flag_cramers_v flags it unless the confounder model within the
  modal member stratum keeps the sign at p < 0.05; module genes vs the culture gene sets (hypergeometric, BH).
- Hypermutation / TP53: HOT/ALLELE/LOFBI determinants refitted on non-hypermutated lines (MSI ∪ top TMB), SIG
  determinants within TP53-mutant lines.
Writes results/locus/{resolution.csv, amplicon.csv}, results/subtypes/{culture_audit.csv, culture_overlap.csv,
hypermut_tp53.csv}, results/m5_audits_summary.json; timing stems locus_resolution, amplicon_test, culture_audit,
hypermut_tp53. Refuses unless HEAD descends from freeze-bounds.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd  # noqa: E402

from biomark import analyses as AN  # noqa: E402
from biomark import audit as AU  # noqa: E402
from biomark import m5 as M5  # noqa: E402
from biomark import params as PR  # noqa: E402

SCRIPT = "scripts/run_m5_audits.py"


def run(ws, src, workers: int) -> dict:
    t0 = AN.utc_now()
    cat, _, pairs = M5.analysis_pairs(ws, src)
    data = src.data(pairs["gene"].unique())
    sign_of = M5.sign_map(cat)
    effects_of, expression_of = src.loaders(data)
    outs = []
    # §8.4 co-location
    fm = data.fmeta
    cnf = set(PR.pp()["covariates"]["cn_derived_families"]) | {f for f in fm["family"].unique()
                                                                   if str(f).startswith("SIL")}
    grp = fm["locus_group"].reindex(pairs["feature"]).astype(object)
    cn = pairs[fm["family"].reindex(pairs["feature"]).isin(cnf).to_numpy() &
               grp.notna().to_numpy() & (grp.astype(str) != "").to_numpy()]
    parts = M5.map_pairs("locus", data, cn[["gene", "feature"]], workers, sign_of=sign_of) if len(cn) else []
    locus = pd.concat(parts, ignore_index=True) if parts else AU.locus_resolution(data, cn, sign_of)
    locus = locus.merge(pairs, on=["gene", "feature"], how="left")
    outs.append(ws.write_csv(locus, "results/locus/resolution.csv"))
    ws.timing("locus_resolution", {"n_pairs": int(len(locus))}, t0, workers=workers)
    # §8.5 cutting
    t1 = AN.utc_now()
    amp = AU.amplicon_test(data, pairs, effects_of, expression_of).merge(pairs, on=["gene", "feature"], how="left")
    outs.append(ws.write_csv(amp, "results/locus/amplicon.csv"))
    ws.timing("amplicon_test", {"n_pairs": int(len(amp))}, t1, workers=1)
    # §8.6 culture
    t2 = AN.utc_now()
    sub = ws.read_csv("results/subtypes/subtypes.csv")
    lines = ws.read_csv("results/subtypes/subtype_lines.csv")
    cul = M5.culture_audit(data, sub, lines, sign_of)
    outs.append(ws.write_csv(cul, "results/subtypes/culture_audit.csv"))
    mods = ws.read_csv("results/subtypes/modules.csv", genes=str)
    gi = src.gene_info()
    universe = sorted(gi.index[gi["eligible"].astype(bool)].astype(str))          # the arms' gene universe
    sym = src.symbols()
    set_syms = sorted({x for v in PR.dp()["subtypes_s8"]["culture_genesets"].values() for x in v})
    ov = M5.culture_overlap(mods, universe, src.symbol_to_entrez(set_syms)) if len(mods) else pd.DataFrame()
    outs.append(ws.write_csv(ov, "results/subtypes/culture_overlap.csv"))
    ws.timing("culture_audit", {"n_subtypes": int(len(cul)), "n_overlap_tests": int(len(ov))}, t2, workers=1)
    # §8.7 hypermutation / TP53
    t3 = AN.utc_now()
    hyp = M5.hypermut_tp53(data, cat, src.tmb(data))
    hyp = hyp.merge(cat[["gene", "feature", "lists"]], on=["gene", "feature"], how="left")
    outs.append(ws.write_csv(hyp, "results/subtypes/hypermut_tp53.csv"))
    ws.timing("hypermut_tp53", {"n_tests": int(len(hyp))}, t3, workers=1)
    summ = {"locus_resolution": {"n": int(len(locus)), "levels": locus["level"].value_counts().to_dict() if len(locus) else {}},
            "amplicon": {"n": int(len(amp)), "labels": amp["label"].value_counts().to_dict() if len(amp) else {}},
            "culture": {"n_subtypes": int(len(cul)), "flagged": int(cul["culture_flag"].sum()) if len(cul) else 0,
                        "rescued": int(cul["rescued"].sum()) if len(cul) else 0,
                        "confounded": int(cul["culture_confounded"].sum()) if len(cul) else 0,
                        "module_overlap_flags": int(ov["flag"].sum()) if len(ov) else 0},
            "hypermut_tp53": {"n": int(len(hyp)), "failed": hyp.groupby("test")["failed"].sum().astype(int).to_dict()
                              if len(hyp) else {}},
            "symbols_known": int(len(sym))}
    outs.append(ws.write_json(summ, "results/m5_audits_summary.json"))
    ws.register(SCRIPT, outs)
    return {"locus": locus, "amplicon": amp, "culture": cul, "culture_overlap": ov, "hypermut": hyp, "summary": summ}


def main(argv=None) -> int:
    a = M5.parse_args(argv, __doc__, workers=True)
    ws, src = M5.open_workspace(a)
    r = run(ws, src, a.workers)
    print(r["summary"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
