#!/usr/bin/env python3
"""Targeted fetch of the ClinicalTrials.gov (API v2) and PubMed (E-utilities esummary) records
named in configs/cgt_curated_mapping.csv (DESIGN §10.1: Tier A, <= 30 CTG records).

Writes
  external/clinical/ctgov/<NCT>.json      raw API v2 study record
  external/clinical/ctgov/manifest.json   {NCT: {url, retrieved_utc, sha256, size, api_version}}
  external/clinical/pubmed/<PMID>.json    raw esummary JSON (db=pubmed, retmode=json)
  external/clinical/pubmed/manifest.json  {PMID: {url, retrieved_utc, sha256, size}}

Network policy: no user identifiers are sent (no email/tool parameters, generic User-Agent).
Idempotent: records already cached are kept unless --refresh.

Usage:
  python -X utf8 scripts/fetch_cgt_records.py                  # ids from the curated mapping
  python -X utf8 scripts/fetch_cgt_records.py --nct NCT05245500 --pmid 19204207
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import re
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
CTG_DIR = ROOT / "external" / "clinical" / "ctgov"
PM_DIR = ROOT / "external" / "clinical" / "pubmed"
MAPPING = ROOT / "configs" / "cgt_curated_mapping.csv"
UA = {"User-Agent": "Mozilla/5.0 (research data fetch; python-requests)"}
CTG_URL = "https://clinicaltrials.gov/api/v2/studies/{}"
CTG_VERSION_URL = "https://clinicaltrials.gov/api/v2/version"
PM_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={}&retmode=json"
MAX_CTG = 30


def utcnow() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def ids_from_mapping() -> tuple[list[str], list[str]]:
    nct, pmid = set(), set()
    with open(MAPPING, encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            for tok in str(row.get("source_record_ids", "")).split(";"):
                tok = tok.strip()
                if tok.startswith("CTG:"):
                    nct.add(tok[4:])
                elif tok.startswith("PMID:") and row.get("source_type") == "pmid_only":
                    pmid.add(tok[5:])
    return sorted(nct), sorted(pmid)


def _get(url: str) -> requests.Response:
    last = None
    for attempt in range(3):
        try:
            r = requests.get(url, headers=UA, timeout=60)
            if r.status_code == 200:
                return r
            last = RuntimeError(f"HTTP {r.status_code} for {url}")
        except requests.RequestException as exc:  # pragma: no cover - network
            last = exc
        time.sleep(1.5 * (attempt + 1))
    raise last  # type: ignore[misc]


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def _save(path: Path, obj: dict) -> None:
    path.write_text(json.dumps(obj, indent=2, sort_keys=True), encoding="utf-8")


def fetch_ctg(ncts: list[str], refresh: bool) -> None:
    if len(ncts) > MAX_CTG:
        raise SystemExit(f"{len(ncts)} NCT ids > Tier-A cap {MAX_CTG}")
    CTG_DIR.mkdir(parents=True, exist_ok=True)
    man_path = CTG_DIR / "manifest.json"
    man = _load(man_path)
    api_version = None
    for nct in ncts:
        if not re.fullmatch(r"NCT\d{8}", nct):
            raise SystemExit(f"bad NCT id {nct!r}")
        out = CTG_DIR / f"{nct}.json"
        if out.exists() and nct in man and not refresh:
            continue
        if api_version is None:
            api_version = _get(CTG_VERSION_URL).json()
        url = CTG_URL.format(nct)
        r = _get(url)
        r.json()  # must parse
        out.write_bytes(r.content)
        man[nct] = {"url": url, "retrieved_utc": utcnow(), "sha256": hashlib.sha256(r.content).hexdigest(),
                    "size": len(r.content), "api_version": api_version}
        print("CTG", nct, len(r.content))
        time.sleep(0.5)
    _save(man_path, man)


def fetch_pubmed(pmids: list[str], refresh: bool) -> None:
    PM_DIR.mkdir(parents=True, exist_ok=True)
    man_path = PM_DIR / "manifest.json"
    man = _load(man_path)
    for pmid in pmids:
        if not pmid.isdigit():
            raise SystemExit(f"bad PMID {pmid!r}")
        out = PM_DIR / f"{pmid}.json"
        if out.exists() and pmid in man and not refresh:
            continue
        url = PM_URL.format(pmid)
        r = _get(url)
        js = r.json()
        rec = js.get("result", {}).get(pmid, {})
        if not rec or "error" in rec:
            raise SystemExit(f"PMID {pmid} not resolvable via esummary: {rec}")
        out.write_bytes(r.content)
        man[pmid] = {"url": url, "retrieved_utc": utcnow(), "sha256": hashlib.sha256(r.content).hexdigest(),
                     "size": len(r.content), "title": rec.get("title", ""), "source": rec.get("source", ""),
                     "pubdate": rec.get("pubdate", "")}
        print("PMID", pmid, rec.get("source"), rec.get("pubdate"), rec.get("title", "")[:90])
        time.sleep(0.4)
    _save(man_path, man)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nct", nargs="*", default=None)
    ap.add_argument("--pmid", nargs="*", default=None)
    ap.add_argument("--refresh", action="store_true")
    a = ap.parse_args()
    if a.nct is None and a.pmid is None:
        nct, pmid = ids_from_mapping()
    else:
        nct, pmid = a.nct or [], a.pmid or []
    fetch_ctg(nct, a.refresh)
    fetch_pubmed(pmid, a.refresh)


if __name__ == "__main__":
    main()
