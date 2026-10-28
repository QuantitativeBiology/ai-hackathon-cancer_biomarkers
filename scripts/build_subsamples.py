"""Stability subsamples (DESIGN v3 4.1.5): n_complementary_pairs pairs of complementary halves of D,
patient-grouped and lineage-stratified, drawn once from seeds.json:subsample_seed into configs/subsamples.json.
Usage: python scripts/build_subsamples.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402

from biomark import covariates, p1  # noqa: E402
from biomark import paths as P  # noqa: E402


def draw(seed: int, n_pairs: int) -> list[list[str]]:
    sp = p1.splits()
    sp = sp[sp["role"] == "dev"].sort_index()
    lin = covariates.pooled_lineage(sp["OncotreeLineage"], sp["role"])
    pat = sp["PatientID"].fillna(sp.index.to_series())
    plin = lin.groupby(pat).agg(lambda s: s.value_counts().index[0])     # patient -> modal pooled lineage
    psize = pat.value_counts()
    rng = np.random.default_rng(seed)
    subs = []
    for _ in range(n_pairs):
        half = {}
        for _, pts in plin.groupby(plin):
            order = rng.permutation(np.asarray(sorted(pts.index)))
            n = [0, 0]
            for q in order:
                k = int(np.argmin(n)) if n[0] != n[1] else int(rng.integers(len(n)))
                half[q] = k
                n[k] += int(psize[q])
        a = sorted(sp.index[pat.map(half) == 0])
        b = sorted(sp.index[pat.map(half) == 1])
        subs += [a, b]
    return subs


def main() -> int:
    s4 = P.load_config("design_params")["surrogate_s4"]
    seed = int(P.seeds()["subsample_seed"])
    n_pairs = int(s4["n_complementary_pairs"])
    subs = draw(seed, n_pairs)
    out = {"seed_key": "subsample_seed", "seed": seed, "n_pairs": n_pairs,
           "design": s4["subsample_design"], "pair_of_subsample": [i // 2 for i in range(len(subs))],
           "subsamples": subs}
    (P.CONFIGS / "subsamples.json").write_text(json.dumps(out, separators=(",", ":")), encoding="utf-8")
    print(f"{len(subs)} subsamples; sizes {min(map(len, subs))}-{max(map(len, subs))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
