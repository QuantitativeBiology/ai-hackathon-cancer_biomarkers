"""Write results/data_status.json (DESIGN 1.5): which pending downloads are present in external/depmap_extra,
their sha256 (vs configs/depmap_extra_sha256.txt), which fallbacks are active and which source each
dependent feature uses. Usage: python scripts/data_status.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from biomark import datastatus, gitutil  # noqa: E402
from biomark import p1  # noqa: E402
from biomark import paths as P  # noqa: E402


def main() -> int:
    st = datastatus.status(compute_sha=True)
    pres = {r["file"]: r["present"] for r in st["files"]}
    msi_order = P.load_config("design_params")["dictionary_s3"]["msi_source_order"]
    st["sources_selected"] = {
        "MSI": msi_order[0] if pres.get("OmicsGlobalSignatures.csv") else msi_order[1] + " -> " + msi_order[2],
        "ploidy_WGD_CIN_LoH": "OmicsGlobalSignatures" if pres.get("OmicsGlobalSignatures.csv") else "CIN proxy C_instab; others not built",
        "TMB_spectra": "OmicsSomaticMutations" if pres.get("OmicsSomaticMutations.csv") else "not assayed",
        "doubling_time": "CRISPRConfounders.ScreenDoublingTime" if pres.get("CRISPRConfounders.csv") else "none",
        "methylation": "RRBS (METH, SIL-METH; amendment D7)" if pres.get("RRBS") else "not assayed (SIL-other proxy)",
    }
    s3 = P.load_config("design_params")["dictionary_s3"]
    st["thresholds"] = {"msi_threshold": P.load_config("design_params")["data_s1"]["msi_threshold"],
                        "meth_high": s3["methylation"]["meth_high"],
                        "meth_above_lineage_median": P.load_config("build_params")["dictionary"]["meth_above_lineage_median"],
                        "lof_bi_af_min": s3["lof_bi_af_min"], "spectra_min_snvs": s3["spectra"]["min_snvs"]}
    if pres.get("RRBS"):
        from biomark import molecular
        M, info = molecular.rrbs_gene_matrix()
        roles = p1.splits()["role"]
        info["profiled_by_role"] = {r: int(M.index.isin(roles.index[roles == r]).sum())
                                    for r in ("dev", "test", "noncancerous")}
        info["n_by_role"] = {r: int((roles == r).sum()) for r in ("dev", "test", "noncancerous")}
        st["rrbs"] = info
    try:
        st["commit"] = gitutil.head(P.PROJECT)
    except Exception:  # noqa: BLE001
        st["commit"] = ""
    out = P.RESULTS / "data_status.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(st, indent=2), encoding="utf-8")
    for r in st["files"]:
        print(f"{r['file']:40s} present={r['present']!s:5s} sha_match={r['sha256_matches']}")
    print("fallbacks active:", [f["file"] for f in st["fallbacks_active"]])
    return 0 if st["all_present_match_recorded"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
