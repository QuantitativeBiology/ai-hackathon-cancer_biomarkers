"""Combinatorial-screen lookup (DESIGN v3.2 §11; Tier A, never cut) for the catalogue (amendment A7: MD-global;
MD-strict alongside) and S1-obs candidates "loss of A -> dependency on partner B".

    python -B scripts/run_combinatorial_lookup.py              # real
    python -B scripts/run_combinatorial_lookup.py --toy DIR    # toy workspace (external tables under DIR/external/)

Candidates: list pairs (g, d) whose feature d encodes loss of another gene A (families in
m5_m8_params.combinatorial.loss_families; EXPPAR in the low-expression direction); A is a paralog of g (Part 1 top
Ensembl Compara paralogs, either direction), a positive-control spec partner, or other. The pair {A, g} is looked
up in the Esmaeili Anvar 2024 harmonised table (Dede, Parrish, Thompson, CHyMErA, Ito: per-study screened pairs and
the authors' hit lists, per-line effects) and the Flister 2024 bioRxiv supplement (analysis only; its own screens
only). Labels from file rows only: replicated (hit in >= 2 screens or >= 2 cell lines), tested-not-replicated,
untested. Validation rate of SL-direction paralog candidates per list vs a screened background matched on
sequence-identity decile and family size (Fisher exact, CI); MD-global vs S1-obs.
Writes results/combinatorial/{lookup.csv, validation.json}: identifiers and statistics computed here only, never a
row or value of the external tables. Timing stem combinatorial_lookup. Refuses unless HEAD descends from
freeze-bounds.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from biomark import analyses as AN  # noqa: E402
from biomark import combinatorial as CB  # noqa: E402
from biomark import m5 as M5  # noqa: E402

SCRIPT = "scripts/run_combinatorial_lookup.py"
OUT = "results/combinatorial"
COLS = ["list", "gene", "symbol", "feature", "family", "partner", "partner_symbol", "relation", "beta", "sl_direction",
        "identity", "family_size", "n_screens_tested", "n_screens_hit", "n_lines_tested", "n_lines_hit", "screens_hit",
        "label"]


def run(ws, src) -> dict:
    t0 = AN.utc_now()
    cat = M5.catalogue(ws)
    s1h = ws.arms(arms=["S1-obs"], columns=["gene", "feature", "beta"], selected_only=True)
    named = {"MD-global": cat[cat["in_md_global"]], "MD-strict": cat[cat["in_md_strict"]], "S1-obs": s1h}
    named = {k: v for k, v in named.items() if k in M5.cfg()["combinatorial"]["lists"]}
    idx = CB.merge_sources(CB.read_esmaeili(ws.external, src.symbol_to_entrez),
                           CB.read_flister(ws.external, src.symbol_to_entrez))
    cand, val = CB.lookup(named, src.fmeta(), src.paralogs(), src.spec_partners(), idx)
    sym = src.symbols()
    cand["symbol"] = cand["gene"].map(sym)
    cand["partner_symbol"] = cand["partner"].map(sym)
    cand = cand[COLS]
    val["_status"] = "ok" if idx["screens"] else "no combinatorial source file found: every candidate untested"
    val["_labels_rule"] = "replicated = hit in >= replicated_min_screens screens or >= replicated_min_lines cell lines"
    val["_licence_note"] = "Flister 2024 bioRxiv file used for analysis only; no external row or value is written"
    outs = [ws.write_csv(cand, f"{OUT}/lookup.csv"), ws.write_json(val, f"{OUT}/validation.json")]
    ws.timing("combinatorial_lookup", {"n_candidates": int(len(cand))}, t0, workers=1)
    ws.register(SCRIPT, outs)
    return {"lookup": cand, "validation": val}


def main(argv=None) -> int:
    a = M5.parse_args(argv, __doc__)
    ws, src = M5.open_workspace(a)
    r = run(ws, src)
    c = r["lookup"]
    print(f"combinatorial: {len(c)} candidate rows; labels {c['label'].value_counts().to_dict()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
