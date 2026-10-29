"""U_ctrl gene list (DESIGN v3.2 §2): non-eligible positive-control dependency genes and CGT target genes.

    python -B scripts/build_u_ctrl_genes.py          # -> configs/u_ctrl_genes.txt (one Entrez id per line)

U_ctrl = ({dependency_entrez of configs/positive_controls.csv} ∪ {CGT target genes of configs/clinical_ground_truth.csv
target_set}) minus the eligible genes of configs/eligible_genes.csv. Frozen configs only: no effects, predictions or
results are read. CGT target symbols are mapped to Entrez through the symbol column of eligible_genes.csv (the library
gene table); a symbol that is missing or ambiguous there stops the build. The list is asserted disjoint from the
eligible set and consistent with eligible_genes.csv control_gene. U_ctrl genes run in a separate universe
(scripts/run_arms.py --universe ctrl) and never enter a top-k list, the null subset or P1/P2.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd  # noqa: E402

from biomark import endpoints as EP  # noqa: E402
from biomark import paths as P  # noqa: E402

OUT = "u_ctrl_genes.txt"


def u_ctrl_table(configs: Path | None = None) -> pd.DataFrame:
    """Rows (entrez, symbol, source) of U_ctrl, sorted by numeric Entrez; source = 'control', 'cgt' or 'control|cgt'."""
    c = Path(configs or P.CONFIGS)
    el = pd.read_csv(c / "eligible_genes.csv", dtype={"entrez": str})
    eligible = set(el.loc[el["eligible"].astype(bool), "entrez"])
    pc = pd.read_csv(c / "positive_controls.csv", dtype=str)
    ctrl = set(pc["dependency_entrez"].dropna().astype(str))
    if "control_gene" in el:
        flagged = set(el.loc[el["control_gene"].astype(bool), "entrez"])
        if flagged != ctrl:
            raise SystemExit(f"eligible_genes.csv control_gene != positive_controls.csv genes: {sorted(flagged ^ ctrl)}")
    cgt = pd.read_csv(c / "clinical_ground_truth.csv", dtype=str)
    syms = sorted({x for v in cgt["target_set"].dropna() for x in EP._targets(v)})
    by_sym = el.groupby("symbol")["entrez"].apply(list)
    bad = [s for s in syms if s not in by_sym.index or len(by_sym[s]) != 1]
    if bad:
        raise SystemExit(f"CGT target symbols missing or ambiguous in eligible_genes.csv: {bad}")
    tgt = {by_sym[s][0] for s in syms}
    unknown = (ctrl | tgt) - set(el["entrez"])
    if unknown:
        raise SystemExit(f"genes outside the library gene table: {sorted(unknown)}")
    genes = sorted((ctrl | tgt) - eligible, key=int)
    sym = el.set_index("entrez")["symbol"]
    src = ["|".join(s for s, m in (("control", g in ctrl), ("cgt", g in tgt)) if m) for g in genes]
    out = pd.DataFrame({"entrez": genes, "symbol": sym.reindex(genes).to_numpy(), "source": src})
    assert not set(out["entrez"]) & eligible, "U_ctrl must be disjoint from the eligible set"
    return out


def read_u_ctrl(path: Path | None = None) -> list[str]:
    f = Path(path or P.CONFIGS / OUT)
    return [x.strip() for x in f.read_text(encoding="utf-8").splitlines() if x.strip()]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="only compare with the existing file")
    a = ap.parse_args(argv)
    t = u_ctrl_table()
    f = P.CONFIGS / OUT
    if a.check:
        ok = f.exists() and read_u_ctrl(f) == t["entrez"].tolist()
        print("u_ctrl_genes.txt", "matches" if ok else "DIFFERS")
        return 0 if ok else 1
    f.write_text("".join(f"{g}\n" for g in t["entrez"]), encoding="utf-8", newline="\n")
    print(t.to_string(index=False))
    print(f"{len(t)} U_ctrl genes ({(t['source'].str.contains('control')).sum()} control, "
          f"{(t['source'].str.contains('cgt')).sum()} CGT target) -> {f.relative_to(P.PROJECT).as_posix()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
