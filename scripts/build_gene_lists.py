"""Gene lists for the dictionary (DESIGN v3 3), each with a source column:
configs/gamma_genes.csv (restricted gene set), configs/fusion_drivers.csv, configs/partner_priors.csv
(top-3 Ensembl paralogs from Part 1 priors + positive-control spec partners, per target gene).
Targets = eligible genes (configs/eligible_genes.csv) U control dependency genes. Usage:
python scripts/build_gene_lists.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd  # noqa: E402

from biomark import genesets  # noqa: E402
from biomark import p1  # noqa: E402,F401
from biomark import paths as P  # noqa: E402


def targets() -> list[str]:
    el = pd.read_csv(P.CONFIGS / "eligible_genes.csv", dtype={"entrez": str})
    elig = list(el.loc[el["eligible"], "entrez"])
    ctrl = list(genesets.control_dependency_genes()["entrez"].dropna())
    return list(dict.fromkeys(elig + ctrl))


def main() -> int:
    el = pd.read_csv(P.CONFIGS / "eligible_genes.csv", dtype={"entrez": str})
    elig = list(el.loc[el["eligible"], "entrez"])
    gam = genesets.build_gamma(elig)
    gam.to_csv(P.CONFIGS / "gamma_genes.csv", index=False)
    drv = genesets.build_fusion_drivers()
    drv.to_csv(P.CONFIGS / "fusion_drivers.csv", index=False)
    pp = genesets.build_partner_priors(targets())
    pp.to_csv(P.CONFIGS / "partner_priors.csv", index=False)
    print(f"gamma {len(gam)}; fusion drivers {len(drv)} ({drv['orientation'].value_counts().to_dict()}); "
          f"partner rows {len(pp)} ({pp['relation'].str.split(':').str[0].value_counts().to_dict()})")
    unresolved = genesets.control_partner_genes()
    miss = unresolved[genesets.symbol_to_entrez(unresolved["partner_symbol"]).isna().to_numpy()]
    if len(miss):
        print("unresolved control partner symbols:", miss.to_dict("records"))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
