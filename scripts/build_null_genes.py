"""Pre-draw the 300-gene null subset for real data (DESIGN v3.2 §7.1) -> configs/null_genes.txt.

Exactly the draw of ``biomark.null.null_genes`` on ``Data.from_real().genes`` without building ``Data`` (no
dependency data, predictions or effects are read):
    genes = [str(entrez) for eligible rows of configs/eligible_genes.csv, in file order]   (pipeline.Data.from_real)
    rng = np.random.default_rng(null_params()['gene_subset_seed'])
    sorted(rng.choice(np.asarray(genes, dtype=object), size=null_params()['gene_subset_n'], replace=False).tolist())
One Entrez id per line (the id format of Data.from_real().genes).

    python -X utf8 scripts/build_null_genes.py [--check]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from biomark import params as PR  # noqa: E402
from biomark import paths as P  # noqa: E402


def eligible_genes() -> list:
    """Same expression as pipeline.Data.from_real(genes=None)."""
    el = pd.read_csv(P.CONFIGS / "eligible_genes.csv", dtype={"entrez": str})
    genes = list(el.loc[el["eligible"].astype(bool), "entrez"])
    return [str(g) for g in genes]


def draw(genes: list) -> list:
    """Same expression as null.null_genes for non-toy data without configs/null_genes.txt."""
    np_ = PR.null_params()
    n = np_["gene_subset_n"]
    if len(genes) <= n:
        return list(genes)
    rng = np.random.default_rng(np_["gene_subset_seed"])
    return sorted(rng.choice(np.asarray(genes, dtype=object), size=n, replace=False).tolist())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="verify the existing file instead of writing")
    a = ap.parse_args(argv)
    genes = eligible_genes()
    sub = draw(genes)
    assert len(sub) == len(set(sub)) == min(len(genes), PR.null_params()["gene_subset_n"])
    assert set(sub) <= set(genes), "null subset must be a subset of the eligible genes (DESIGN §2)"
    out = P.CONFIGS / "null_genes.txt"
    text = "\n".join(sub) + "\n"
    if a.check:
        ok = out.exists() and out.read_text(encoding="utf-8") == text
        print(f"{out}: {'matches' if ok else 'DIFFERS'} the draw ({len(sub)} of {len(genes)} eligible genes)")
        return 0 if ok else 1
    out.write_bytes(text.encode("utf-8"))
    print(f"wrote {out}: {len(sub)} of {len(genes)} eligible genes, seed {PR.null_params()['gene_subset_seed']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
