"""Frozen-file and freeze-tag check (DESIGN §13.4).

    python -X utf8 scripts/check_frozen.py [--root DIR] [--json out.json]

1. ``configs/frozen_hashes.json``: every file of section ``design`` (tag ``freeze-design``) and, when present,
   ``bounds`` (tag ``freeze-bounds``) must have the recorded sha256 (LF-normalised) now and at its tag. DESIGN.md is
   the one file allowed to differ from ``design``: then ``design_amendments.DESIGN.md`` must match the current file
   and the text above ``## Amendments after freeze`` must equal DESIGN.md at ``freeze-design`` (trailing
   whitespace ignored).
2. If tag ``freeze-bounds`` exists it must be an ancestor of the first commit (in HEAD's history) touching
   ``results/benchmark/``, ``results/md/``, ``results/gates/``, ``results/null/`` or the ctrl-universe
   ``results/benchmark_ctrl/`` / ``results/md_ctrl/``; if it does not exist, no such
   commit may exist yet.
3. Coverage (reported, ``--strict-coverage`` makes it fail): files DESIGN §13.4 names for each freeze tag that are not
   in the hash record, with whether they existed at the tag and changed since.

A mismatch prints any RESEARCH_TREE.md decision line (``- D<n> ...``) naming the file; a decision does not turn a
mismatch into a pass (the record must be updated with the decision). Exit 1 on any FAIL.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HASHES = Path("configs") / "frozen_hashes.json"
SECTION_TAG = {"design": "freeze-design", "bounds": "freeze-bounds"}
AMEND_HEADING = b"## Amendments after freeze"
PROTECTED = ("results/benchmark/", "results/md/", "results/gates/", "results/null/", "results/benchmark_ctrl/",
             "results/md_ctrl/", "results/benchmark_lowtier/", "results/md_lowtier/")
# DESIGN §13.4 lists (file paths as found in this repository; "culture and lineage lists" live in build_params.json,
# toy plants in toy_answer_key.json).
EXPECTED = {
    "design": ["DESIGN.md", "configs/design_params.json", "configs/primary_endpoints.json", "configs/tiers.json",
               "configs/seeds.json", "configs/build_params.json", "configs/positive_controls_spec.csv",
               "configs/cgt_allele_rules.json", "configs/fusion_known_positives.csv", "configs/cgt_corrections.csv",
               "configs/toy_answer_key.json"],
    "bounds": ["results/toy/nc_bounds.json", "results/controls/pc_eligibility.csv", "configs/positive_controls.csv",
               "configs/null_genes.txt", "configs/subsamples.json", "configs/dictionary_v1.json",
               "results/data_status.json"],
}


def sha_lf(data: bytes) -> str:
    return hashlib.sha256(data.replace(b"\r\n", b"\n")).hexdigest()


def git(root: Path, *args) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "--no-optional-locks", *args], cwd=str(root), capture_output=True)


def tag_exists(root: Path, tag: str) -> bool:
    return git(root, "rev-parse", "-q", "--verify", f"refs/tags/{tag}").returncode == 0


def at_tag(root: Path, tag: str, rel: str) -> bytes | None:
    r = git(root, "show", f"{tag}:{rel}")
    return r.stdout if r.returncode == 0 else None


def decisions_naming(root: Path, rel: str) -> list:
    p = root / "RESEARCH_TREE.md"
    if not p.exists():
        return []
    name = Path(rel).name
    return [ln.strip()[:160] for ln in p.read_text(encoding="utf-8").splitlines()
            if re.match(r"\s*-\s*D\d+\b", ln) and (rel in ln or name in ln)]


def check_section(root: Path, name: str, sec: dict, amend: dict | None, fails: list, notes: list) -> None:
    tag = SECTION_TAG.get(name)
    has_tag = bool(tag) and tag_exists(root, tag)
    if tag and not has_tag:
        fails.append(f"[{name}] tag {tag} missing although the section exists")
    for rel, h in sec.items():
        if not isinstance(h, str) or not re.fullmatch(r"[0-9a-f]{64}", h):
            continue                                          # comments / metadata keys
        p = root / rel
        if not p.exists():
            fails.append(f"[{name}] {rel}: file missing")
            continue
        cur = sha_lf(p.read_bytes())
        if has_tag:
            old = at_tag(root, tag, rel)
            if old is None:
                fails.append(f"[{name}] {rel}: not present at tag {tag}")
            elif sha_lf(old) != h:
                fails.append(f"[{name}] {rel}: recorded hash != file at tag {tag}")
        if cur == h:
            notes.append(f"[{name}] {rel}: OK")
            continue
        if name == "design" and rel == "DESIGN.md" and amend and amend.get("DESIGN.md"):
            if cur != amend["DESIGN.md"]:
                fails.append(f"[design_amendments] DESIGN.md: current sha256 != design_amendments record "
                             f"(amendments {amend.get('amendments')}); decisions: {decisions_naming(root, rel)}")
                continue
            data = p.read_bytes().replace(b"\r\n", b"\n")
            i = data.find(AMEND_HEADING)
            old = at_tag(root, tag, rel) if has_tag else None
            if i < 0:
                fails.append("[design_amendments] DESIGN.md: no '## Amendments after freeze' heading")
            elif old is not None and data[:i].rstrip() != old.replace(b"\r\n", b"\n").rstrip():
                fails.append("[design_amendments] DESIGN.md: text above the amendments heading differs from freeze-design")
            else:
                notes.append(f"[design_amendments] DESIGN.md: OK (amended {amend.get('amendments')}, "
                             f"text above heading unchanged from {tag})")
            continue
        fails.append(f"[{name}] {rel}: current sha256 {cur[:12]} != recorded {h[:12]}; "
                     f"decisions naming it: {decisions_naming(root, rel) or 'none'}")


def coverage(root: Path, rec: dict) -> list:
    out = []
    for name, files in EXPECTED.items():
        tag = SECTION_TAG[name]
        has_tag = tag_exists(root, tag)
        recorded = set((rec.get(name) or {}).keys())
        for rel in files:
            if rel in recorded:
                continue
            p = root / rel
            line = f"[{name}] {rel} named in DESIGN §13.4 but not in frozen_hashes.{name}"
            if not p.exists():
                line += "; file absent"
            elif has_tag:
                old = at_tag(root, tag, rel)
                line += ("; absent at tag " + tag) if old is None else \
                    ("; unchanged since " + tag if sha_lf(old) == sha_lf(p.read_bytes()) else "; CHANGED since " + tag)
            else:
                line += f"; tag {tag} not yet set"
            out.append(line)
    return out


def check_tag_order(root: Path, fails: list, notes: list) -> None:
    r = git(root, "log", "--reverse", "--format=%H", "HEAD", "--", *PROTECTED)
    commits = r.stdout.decode().split() if r.returncode == 0 else []
    first = commits[0] if commits else None
    if not tag_exists(root, "freeze-bounds"):
        if first:
            fails.append(f"commit {first[:10]} touches {PROTECTED} but tag freeze-bounds does not exist")
        else:
            notes.append("freeze-bounds: tag not set; no commit touches the protected results folders (OK)")
        untracked = [d for d in PROTECTED if (root / d).exists()]
        if untracked:
            notes.append(f"WARN: protected result folders exist in the working tree before freeze-bounds: {untracked}")
        return
    if not first:
        notes.append("freeze-bounds: tag set; no commit touches the protected results folders yet (OK)")
        return
    anc = git(root, "merge-base", "--is-ancestor", "freeze-bounds", first).returncode == 0
    if anc:
        notes.append(f"freeze-bounds is an ancestor of the first protected-results commit {first[:10]} (OK)")
    else:
        fails.append(f"freeze-bounds is NOT an ancestor of {first[:10]}, the first commit touching {PROTECTED}")


def evaluate(root: Path = ROOT, strict_coverage: bool = False) -> dict:
    root = Path(root)
    fails, notes = [], []
    rec = json.loads((root / HASHES).read_text(encoding="utf-8"))
    if "design" not in rec:
        fails.append("frozen_hashes.json has no 'design' section")
    for name in ("design", "bounds"):
        if name in rec:
            check_section(root, name, rec[name], rec.get("design_amendments") if name == "design" else None,
                          fails, notes)
        elif name == "bounds" and tag_exists(root, "freeze-bounds"):
            fails.append("tag freeze-bounds exists but frozen_hashes.json has no 'bounds' section")
    check_tag_order(root, fails, notes)
    cov = coverage(root, rec)
    if strict_coverage:
        fails += [f"coverage: {c}" for c in cov]
    return {"fails": fails, "notes": notes, "coverage": cov}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="DESIGN §13.4 frozen-file check")
    ap.add_argument("--root", type=Path, default=ROOT)
    ap.add_argument("--strict-coverage", action="store_true")
    ap.add_argument("--json", type=Path, default=None)
    a = ap.parse_args(argv)
    res = evaluate(a.root, a.strict_coverage)
    for n in res["notes"]:
        print(n)
    for c in res["coverage"]:
        print(f"COVERAGE: {c}")
    for f in res["fails"]:
        print(f"FAIL: {f}")
    print("RESULT: " + ("FAIL" if res["fails"] else "PASS"))
    if a.json:
        a.json.write_text(json.dumps(res, indent=1), encoding="utf-8")
    return 1 if res["fails"] else 0


if __name__ == "__main__":
    sys.exit(main())
