"""Candidate cards (DESIGN v3.2 §12.3) from frozen outputs only (no effects, no predictions, no refits).

    python -B scripts/build_candidate_cards.py              # real (after run_limits and, for the shortlist, replicate)
    python -B scripts/build_candidate_cards.py --toy DIR    # toy workspace

- results/candidates/known_recovery.csv: one row per positive control and per testable CGT test (headline-eligible
  set incl. proxy copies): best rank_of_feature_D and sign per arm, recovery in MD-global / MD-strict with the
  expected sign, R result (pending until results/replication/r_tests.csv exists), ceiling label, access class.
- results/candidates/shortlist.csv: class-9 MD-global pairs replicated in R, ranked by cross-lineage label
  precedence, R one-sided p, number of powered lineages; columns g, d, access class, powered lineages with
  per-lineage β, carriers D / R, delta flags, locus resolution, combinatorial label, CGT tests, set-A dependent
  fraction (normal-line flag), not_in_curated_lists. Empty (status 'pending R') before the single R run. Asserted:
  shortlist ⊆ (class 9 ∧ R replicated).
- results/candidates/cards_summary.json (top report rows and the three case-study pairs).
Timing stem candidate_cards. Refuses unless HEAD descends from freeze-bounds.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from biomark import analyses as AN  # noqa: E402
from biomark import arms as A  # noqa: E402
from biomark import endpoints as EP  # noqa: E402
from biomark import limits as LM  # noqa: E402
from biomark import m5 as M5  # noqa: E402
from biomark import params as PR  # noqa: E402

SCRIPT = "scripts/build_candidate_cards.py"
OUT = "results/candidates"


def check_shortlist(sl, pc, rt) -> None:
    """DESIGN §12.3: shortlist ⊆ (class 9 ∧ R replicated)."""
    if not len(sl):
        return
    valid = PR.pp()["pair_classes"]["valid"]
    ok9 = set(zip(pc.loc[pc["pair_class"] == valid, "gene"], pc.loc[pc["pair_class"] == valid, "feature"]))
    R = LM.r_status(rt)
    bad = [k for k in zip(sl["gene"], sl["feature"]) if k not in ok9 or not R.get(k, (False,))[0]]
    if bad:
        raise AssertionError(f"shortlist rows outside class 9 ∧ R replicated: {bad[:5]}")


def run(ws, src) -> dict:
    t0 = AN.utc_now()
    L = M5.cfg()["limits"]
    cards = PR.dp()["limits_s12"]["cards"]
    cat = M5.catalogue(ws)
    ctl = src.controls()
    cgt = src.cgt()
    sym = src.symbols()
    rt = ws.read_csv("results/replication/r_tests.csv", required=False)
    ceil = LM.ceiling_labels(ws.read_csv(L["ceiling_file"], required=False, test_id=str))
    # testable CGT tests (headline-eligible incl. proxy copies); target genes already entrez in src.cgt()
    tg_genes = sorted({g for ts in cgt["targets"] for g in ts})
    Ut = ws.arms(arms=["S1-obs"], columns=["gene", "feature", "n_carriers"], genes=tg_genes)
    cg = cgt.assign(target_set=[";".join(ts) for ts in cgt["targets"]])
    T = EP.p2_test_set(cg, Ut, "incl_proxy")
    T = T[T["testable"]]
    genes = set(ctl["gene"]) | {g for ts in T["targets"] for g in ts}
    arms_list = [a for a in A.ARMS]
    ridx = M5.RankLookup(ws, arms_list, sorted(genes | set(cat["gene"])))
    kr = LM.known_recovery(ridx, arms_list, ctl, T, cat, rt, ceil, ridx.access, sym)
    pc = ws.read_csv("results/limits/pair_classes.csv")
    for c in ("in_md_global", "in_md_strict"):
        pc[c] = pc[c].astype(str).str.lower().isin(["true", "1"])
    cl = ws.read_csv("results/lineage/cross_lineage.csv", required=False)
    per = ws.read_csv("results/lineage/cross_lineage_per_lineage.csv", required=False)
    delta = ws.read_csv("results/confounders/delta.csv", required=False)
    dflags = {}
    if delta is not None and len(delta):
        d = delta[delta["flagged"].astype(str).str.lower().isin(["true", "1"])]
        dflags = d.groupby(["gene", "feature"])["block"].agg("|".join).to_dict()
    locus = LM._key(ws.read_csv("results/locus/resolution.csv", required=False), "level")
    comb = ws.read_csv("results/combinatorial/lookup.csv", required=False)
    combl = {}
    if comb is not None and len(comb):
        c = comb[comb["list"] == M5.PRIMARY]
        combl = c.groupby(["gene", "feature"])["label"].agg(lambda v: "|".join(sorted(set(v)))).to_dict()
    sa = ws.read_csv("results/replication/setA_lines.csv", required=False)
    setA = dict(zip(sa["gene"], sa["frac_dependent"])) if sa is not None and "frac_dependent" in sa else {}
    sl = LM.shortlist(pc, rt, cl, per, cat, ridx.access, sym, dflags, locus, combl, cgt, setA, src.curated_pairs())
    check_shortlist(sl, pc, rt)
    outs = [ws.write_csv(kr, f"{OUT}/known_recovery.csv"), ws.write_csv(sl, f"{OUT}/shortlist.csv")]
    top = int(cards["report_top_rows"])
    ncs = int(cards["case_studies"])
    summ = {"status": "ok" if rt is not None else L["pending"],
            "n_known_recovery": int(len(kr)), "n_controls": int((kr["kind"] == "control").sum()),
            "n_cgt_tests": int((kr["kind"] == "cgt").sum()),
            "recovered_md_global": int(kr["recovered_MD-global"].sum()), "recovered_md_strict": int(kr["recovered_MD-strict"].sum()),
            "n_shortlist": int(len(sl)), "report_top_rows": top,
            "case_studies": [f"{g}|{f}" for g, f in zip(sl["gene"].head(ncs), sl["feature"].head(ncs))],
            "rule": "shortlist = class-9 MD-global pairs replicated in R (asserted); 'novel' is never applied to a gene"}
    outs.append(ws.write_json(summ, f"{OUT}/cards_summary.json"))
    ws.timing("candidate_cards", {"n_shortlist": int(len(sl)), "n_known": int(len(kr))}, t0, workers=1)
    ws.register(SCRIPT, outs)
    return {"known_recovery": kr, "shortlist": sl, "summary": summ}


def main(argv=None) -> int:
    a = M5.parse_args(argv, __doc__)
    ws, src = M5.open_workspace(a)
    r = run(ws, src)
    print(f"cards: known_recovery {len(r['known_recovery'])} rows; shortlist {len(r['shortlist'])} ({r['summary']['status']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
