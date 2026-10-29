"""NC2 gene set (DESIGN v3.2 §7.3): 100 low-tier genes matched in D effect SD to the eligible genes -> configs/lowtier_genes.txt.

    python -B scripts/build_lowtier_genes.py [--check]

Frozen inputs only (no effects are read): configs/eligible_genes.csv (eligible flag, m5 tier) and the D effect SD of
Part 1's gene table (biomark.p1.gene_classes()['sd'], 1,022 dev lines; read-only). The design fixes the count
(negative_controls_s7.NC2_lowtier_n_genes), the seed (seeds.json lowtier_seed) and "matched-SD"; the matching is
operationalised here (agent rule, stated in the output header line of the summary):
  1. reference = NC2_lowtier_n_genes eligible genes drawn without replacement with default_rng(lowtier_seed)
     (eligible genes in file order, as build_null_genes draws the null subset);
  2. candidates = genes with tier 'low' and a finite SD;
  3. for each reference gene in draw order, take the unused candidate with the smallest |SD - SD_ref| (ties: lowest
     Entrez id).
The list is disjoint from the eligible set by construction (eligibility requires tier high/medium); asserted.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from biomark import params as PR  # noqa: E402
from biomark import paths as P  # noqa: E402

OUT = "lowtier_genes.txt"
SUMMARY = P.RESULTS / "controls" / "lowtier_genes_summary.json"


def match(el: pd.DataFrame, sd: pd.Series, n: int, seed: int) -> pd.DataFrame:
    """Rows (entrez, sd, reference_entrez, reference_sd) of the matched low-tier genes (see module doc)."""
    el = el.assign(entrez=el["entrez"].astype(str))
    sd = sd.copy()
    sd.index = sd.index.astype(str)
    elig = [g for g in el.loc[el["eligible"].astype(bool), "entrez"] if np.isfinite(sd.get(g, np.nan))]
    ref = np.random.default_rng(seed).choice(np.asarray(elig, dtype=object), size=min(n, len(elig)), replace=False)
    cand = el.loc[(el["tier"] == "low"), "entrez"]
    cand = pd.Series({g: float(sd[g]) for g in cand if np.isfinite(sd.get(g, np.nan))})
    cand = cand.iloc[np.argsort(cand.index.astype(int), kind="stable")]
    used, rows = set(), []
    for g in ref:
        c = cand[~cand.index.isin(used)]
        if c.empty:
            break
        d = (c - float(sd[g])).abs()
        best = d.index[int(np.argmin(d.to_numpy()))]          # first minimum = lowest Entrez id among ties
        used.add(best)
        rows.append({"entrez": best, "sd": float(c[best]), "reference_entrez": str(g), "reference_sd": float(sd[g])})
    out = pd.DataFrame(rows)
    assert not set(out["entrez"]) & set(el.loc[el["eligible"].astype(bool), "entrez"]), "low-tier set overlaps eligible"
    return out


def main(argv=None) -> int:
    from biomark import p1
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args(argv)
    el = pd.read_csv(P.CONFIGS / "eligible_genes.csv", dtype={"entrez": str})
    n = int(PR.dp()["negative_controls_s7"]["NC2_lowtier_n_genes"])
    seed = int(PR.seeds()["lowtier_seed"])
    m = match(el, p1.gene_classes()["sd"], n, seed)
    genes = sorted(m["entrez"], key=int)
    f = P.CONFIGS / OUT
    text = "".join(f"{g}\n" for g in genes)
    if a.check:
        ok = f.exists() and f.read_text(encoding="utf-8") == text
        print(f"{f.name}: {'matches' if ok else 'DIFFERS'}")
        return 0 if ok else 1
    f.write_text(text, encoding="utf-8", newline="\n")
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY.write_text(json.dumps({
        "n": len(genes), "seed_key": "lowtier_seed", "seed": seed,
        "rule": "reference = n eligible genes drawn with default_rng(lowtier_seed); each matched to the unused tier-low "
                "gene with the nearest D effect SD (Part 1 gene table), ties by lowest Entrez id (agent rule)",
        "sd_median_lowtier": float(m["sd"].median()), "sd_median_reference": float(m["reference_sd"].median()),
        "max_abs_sd_gap": float((m["sd"] - m["reference_sd"]).abs().max()),
        "pairs": m.to_dict("records")}, indent=1), encoding="utf-8")
    from biomark import analyses as AN
    AN.append_registry("scripts/build_lowtier_genes.py", "ok", [f"configs/{OUT}", "results/controls/lowtier_genes_summary.json"])
    print(f"{len(genes)} low-tier genes; SD median {m['sd'].median():.3f} vs reference {m['reference_sd'].median():.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
