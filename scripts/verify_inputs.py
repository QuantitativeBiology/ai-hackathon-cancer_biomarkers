"""Verify Part 1 inputs before any Part 2 analysis.

1. sha256 of every file listed in P1/configs/depmap_local_sha256.txt (files in P1/data/) must match.
2. Part 1 git HEAD must equal configs/part1_commit.txt and its tree must be clean.
3. Records (no reference to compare against) the sha256 of the Part 1 caches Part 2 reads, so later
   runs can detect drift: raw parquet caches, feature blocks, the headline/component OOF predictions
   and the m5_avg test predictions.

Writes results/input_verification.json; exit code 1 on any mismatch. Reads only; never writes into P1.
Usage: python scripts/verify_inputs.py [--no-caches]
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from biomark import gitutil  # noqa: E402
from biomark import p1  # noqa: E402,F401  (installs the Part 1 write blocker)
from biomark import paths as P  # noqa: E402


def sha256(path: Path, chunk: int = 8 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def parse_sha_file(path: Path) -> dict[str, str]:
    out = {}
    for ln in path.read_text(encoding="utf-8").splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        digest, name = ln.split(None, 1)
        out[name.lstrip("*").strip()] = digest.lower()
    return out


def cache_files() -> list[Path]:
    fs = [P.P1_CRISPR_PARQUET, P.P1_GENE_MAP, *P.P1_RAW_CACHES.values()]
    fs += sorted(P.P1_FEATURES.glob("*.parquet"))
    fs += [P.P1_RUNS / r / "pred.parquet" for r in P.P1_RUN_IDS.values()]
    fs += [P.P1_TEST_PREDS / "m5_avg.parquet"]
    fs += [P.P1_SPLITS, P.P1_GENE_CLASSES, *P.P1_PREDICTABLE_GENES.values(), P.P1_LIBRARY_GROUPS]
    return fs


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-caches", action="store_true", help="skip recording cache checksums")
    a = ap.parse_args(argv)
    t0 = time.time()
    ok = True
    ref = parse_sha_file(P.P1_SHA256)
    data = {}
    for name, digest in ref.items():
        f = P.P1_DATA / name
        if not f.exists():
            data[name] = {"expected": digest, "observed": None, "match": False, "bytes": None}
            ok = False
            continue
        got = sha256(f)
        data[name] = {"expected": digest, "observed": got, "match": got == digest, "bytes": f.stat().st_size}
        ok &= got == digest
        print(f"{'OK ' if got == digest else 'BAD'} {name}", flush=True)

    head = gitutil.head(P.P1)
    pinned = P.part1_commit()
    dirty = gitutil.status_porcelain(P.P1)
    ok &= head == pinned and not dirty

    caches = {}
    if not a.no_caches:
        for f in cache_files():
            caches[f.relative_to(P.P1).as_posix()] = {"sha256": sha256(f), "bytes": f.stat().st_size}
            print(f"rec {f.relative_to(P.P1).as_posix()}", flush=True)

    out = {
        "utc_time": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "project_b_commit": gitutil.head(P.PROJECT),
        "part1_root": str(P.P1),
        "part1_head": head,
        "part1_pinned_commit": pinned,
        "part1_head_matches_pin": head == pinned,
        "part1_tree_clean": not dirty,
        "part1_dirty_entries": dirty[:50],
        "sha256_reference": P.P1_SHA256.relative_to(P.P1).as_posix(),
        "data_files": data,
        "n_data_files": len(data),
        "n_data_match": sum(v["match"] for v in data.values()),
        "p1_cache_sha256_recorded": caches,
        "all_ok": bool(ok),
        "seconds": round(time.time() - t0, 1),
    }
    P.RESULTS.mkdir(parents=True, exist_ok=True)
    P.INPUT_VERIFICATION.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    print(f"all_ok={ok}  data {out['n_data_match']}/{out['n_data_files']}  P1 HEAD {head[:10]} "
          f"(pin {pinned[:10]})  clean={not dirty}  -> {P.INPUT_VERIFICATION}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
