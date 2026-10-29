"""Task 4: confounding delta (DESIGN 5.5) for 10 seeded MD-strict pairs.

Full model = the 5.3 S1 model on observed D effects (my confounder.py, same per-test drops). Drop-one-block refits
for the blocks of confounders_s5.delta_blocks: lineage, library (library group + screen days, which 5.3 places in
the library block), culture, growth, burden (MSI + DMG_BURDEN), TP53, own CN (own CN + its CN_absent indicator).
Flag (confounders_s5.delta_flag): |b_full - b_minus| > 0.5*|b_minus| and |b_full - b_minus| > SE_minus.
Sensitivity: OncotreePrimaryDisease block (levels with >= 5 D lines) added to the full model; flagged
'subtype-confounded' with the same rule (my reading; 5.5 does not restate the rule for this block).
A block already removed by the per-test drops gives b_minus = b_full (no flag).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import data as DA
from .confounder import s1_gene

BLOCKS = {"lineage": {"lineage"}, "library": {"library", "screen_days"}, "culture": {"culture"},
          "growth": {"growth"}, "burden": {"MSI", "DMG_BURDEN"}, "TP53": {"TP53"}, "own CN": {"own_CN"}}


def draw_pairs(n=10, seed_key="seeds", seed_index=2):
    ms = pd.read_csv(DA.ROOT / "results" / "md" / "md_strict.csv", dtype={"gene": str})
    seed = DA.seeds()[seed_key][seed_index]
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(len(ms), n, replace=False))
    return ms.iloc[idx][["gene", "feature", "beta", "se", "p_sc"]].reset_index(drop=True), \
        {"seed_key": f"seeds.json:{seed_key}[{seed_index}]", "seed": int(seed)}


def run(pairs: pd.DataFrame) -> pd.DataFrame:
    E = DA.effects(sorted(set(pairs["gene"])))
    out = []
    for _, r in pairs.iterrows():
        g, f = r["gene"], r["feature"]
        y = E[g].to_numpy(float)
        full = s1_gene(g, y, [f])[0]
        rec = {"gene": g, "feature": f, "beta_full": full["beta"], "se_full": full["se"], "p_full": full["p"],
               "dropped_cov": full["dropped_cov"], "beta_md_strict_engine": r["beta"], "se_md_strict_engine": r["se"],
               "p_md_strict_engine": r["p_sc"]}
        flags = []
        for blk, names in list(BLOCKS.items()) + [("primary_disease", {"primary_disease"})]:
            mm = s1_gene(g, y, [f], drop_extra={f: names})[0]
            bm, sm = mm["beta"], mm["se"]
            d = abs(full["beta"] - bm)
            flag = bool(d > 0.5 * abs(bm) and d > sm)
            key = blk.replace(" ", "_")
            rec[f"beta_minus_{key}"] = bm
            rec[f"se_minus_{key}"] = sm
            rec[f"flag_{key}"] = flag
            if flag:
                flags.append(blk)
        rec["flagged_blocks"] = ";".join(flags)
        out.append(rec)
    return pd.DataFrame(out)
