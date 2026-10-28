"""Component readout (DESIGN v3.2 §4.2, Tier A): the EN (unchanged engine) on fold-centred ŷ_M1, ŷ_M2n − ŷ_M1 and
ŷ_M3n − ŷ_M1 for the MD-strict genes, on the same U rows, subsamples and configs as EN-pred.

    python -B scripts/run_component_readout.py [--workers N]

Writes results/components/component_selection.csv (per MD-strict pair: full-D β and selection frequency per component,
layer = expression-mediated / genomics-layer / both / none, carrier_input for proxy-class features),
cache/benchmark/component_en.parquet and results/timing/component_readout.json. Reads predictions only (no effects).
Refuses unless HEAD descends from freeze-bounds.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd  # noqa: E402

from biomark import analyses as AN  # noqa: E402
from biomark import arms as A  # noqa: E402
from biomark import endpoint_io as IO  # noqa: E402
from biomark import endpoints as EP  # noqa: E402
from biomark import p1  # noqa: E402
from biomark import params as PR  # noqa: E402
from biomark import paths as P  # noqa: E402
from biomark import pipeline as PL  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=None)
    a = ap.parse_args(argv)
    EP.require_descends("freeze-bounds")
    t0 = time.perf_counter()
    w = a.workers or PR.null_params()["fg_workers"]
    md = IO.md_strict()
    genes = sorted(md["gene"].astype(str).unique())
    data = PL.Data.from_real(genes=genes)
    arms_tab = IO.arms_table()
    U = arms_tab[(arms_tab["arm"] == "EN-pred") & arms_tab["gene"].isin(genes)][["gene", "feature"]]
    dev = list(data.lines)
    oof = {k: p1.oof_predictions(k, columns=genes, ids=dev) for k in ("m1", "m2n", "m3n")}
    targets = AN.component_targets(oof, p1.folds().reindex(dev))
    tabs = []
    comp = {}
    for c, T in targets.items():
        dc = AN.data_with_pred(data, T.reindex(index=dev, columns=genes).to_numpy())
        t, _, _ = A.run_en(dc, ("EN-pred",), genes=genes, U=U, workers=w)
        t = t.assign(arm=c)
        comp[c] = t
        tabs.append(t)
    (P.CACHE / "benchmark").mkdir(parents=True, exist_ok=True)
    pd.concat(tabs, ignore_index=True).to_parquet(P.CACHE / "benchmark" / "component_en.parquet")
    acc = arms_tab[arms_tab["arm"] == "EN-pred"].set_index(["gene", "feature"])["access_class"]
    sel = AN.component_selection(md, comp, acc)
    out = P.RESULTS / "components"
    out.mkdir(parents=True, exist_ok=True)
    sel.to_csv(out / "component_selection.csv", index=False)
    AN.write_timing("component_readout", {"wall_s": time.perf_counter() - t0, "workers": w, "n_genes": len(genes),
                                          "n_pairs": int(len(sel))})
    AN.append_registry("scripts/run_component_readout.py", "ok", ["results/components/component_selection.csv"])
    print(sel["layer"].value_counts().to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
