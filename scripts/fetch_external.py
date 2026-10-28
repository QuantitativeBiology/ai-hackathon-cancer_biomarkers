#!/usr/bin/env python3
"""Fetch external clinical-evidence and combinatorial-screen resources.

Downloads into external/<group>/<source>/ (git-ignored) and records, for every file,
url, retrieval UTC timestamp, size, sha256, declared version/date and licence in
configs/external_manifest.json, plus a sha256sum-style configs/external_sha256.txt
(paths relative to the project root).

Idempotent: a file already on disk whose sha256 matches the manifest is skipped.
Network policy: 60 s timeout, 2 retries; sources that cannot be fetched are recorded
with status "blocked" (manifest key "blocked") so they can be requested from the user.

Usage:
  python scripts/fetch_external.py                 # all sources
  python scripts/fetch_external.py --only civic    # one or more (repeatable / comma list)
  python scripts/fetch_external.py --list
No identifiers of the user are sent anywhere (generic User-Agent only).
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import io
import json
import re
import sys
import time
import zipfile
from html.parser import HTMLParser
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
EXT = ROOT / "external"
MANIFEST = ROOT / "configs" / "external_manifest.json"
SHA_TXT = ROOT / "configs" / "external_sha256.txt"
UA = {"User-Agent": "Mozilla/5.0 (research data fetch; python-requests)"}
TIMEOUT, RETRIES = 60, 2
CAP_BYTES = 200 * 1024 * 1024


def utcnow() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# ----------------------------------------------------------------------------- manifest
def load_manifest() -> dict:
    if MANIFEST.exists():
        return json.loads(MANIFEST.read_text(encoding="utf-8"))
    return {"files": {}, "blocked": {}}


def save_manifest(m: dict) -> None:
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    m["generated_utc"] = utcnow()
    MANIFEST.write_text(json.dumps(m, indent=2, sort_keys=True), encoding="utf-8")
    lines = [f"{e['sha256']}  {p}" for p, e in sorted(m["files"].items())]
    SHA_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class Ctx:
    def __init__(self) -> None:
        self.m = load_manifest()
        self.s = requests.Session()
        self.s.headers.update(UA)

    # ---- http
    def request(self, method: str, url: str, **kw):
        last = None
        for attempt in range(RETRIES + 1):
            try:
                r = self.s.request(method, url, timeout=TIMEOUT, allow_redirects=True, **kw)
                if r.status_code >= 400:
                    raise requests.HTTPError(f"HTTP {r.status_code}", response=r)
                return r
            except Exception as e:  # noqa: BLE001
                last = e
                code = getattr(getattr(e, "response", None), "status_code", None)
                if code in (401, 403, 404, 410):
                    break
                time.sleep(2 * (attempt + 1))
        raise RuntimeError(f"{url}: {last}")

    def get(self, url, **kw):
        return self.request("GET", url, **kw)

    # ---- recording
    def rel(self, p: Path) -> str:
        return p.relative_to(ROOT).as_posix()

    def have(self, p: Path) -> bool:
        e = self.m["files"].get(self.rel(p))
        return bool(e and p.exists() and sha256_file(p) == e["sha256"])

    def record(self, p: Path, url: str, source: str, version: str, licence: str,
               retrieved: str | None = None, note: str = "") -> None:
        self.m["files"][self.rel(p)] = {
            "source": source, "url": url, "retrieved_utc": retrieved or utcnow(),
            "size": p.stat().st_size, "sha256": sha256_file(p),
            "version": version, "licence": licence, "status": "downloaded", "note": note,
        }
        self.m["blocked"].pop(f"{source}|{url}", None)
        print(f"  [ok] {self.rel(p)} ({p.stat().st_size:,} B)")

    def block(self, source: str, url: str, reason: str, licence: str = "") -> None:
        self.m["blocked"][f"{source}|{url}"] = {
            "source": source, "url": url, "status": "blocked", "reason": reason,
            "checked_utc": utcnow(), "licence": licence}
        print(f"  [BLOCKED] {source} {url}: {reason}")

    def download(self, url: str, dest: Path, source: str, version: str, licence: str,
                 note: str = "", headers_version=None) -> Path | None:
        """Download url to dest unless present with matching checksum."""
        if self.have(dest):
            print(f"  [skip] {self.rel(dest)}")
            return dest
        dest.parent.mkdir(parents=True, exist_ok=True)
        try:
            r = self.get(url)
        except Exception as e:  # noqa: BLE001
            self.block(source, url, str(e)[:300], licence)
            return None
        dest.write_bytes(r.content)
        v = version
        if headers_version:
            v = headers_version(r) or version
        self.record(dest, url, source, v, licence, note=note)
        return dest


# ----------------------------------------------------------------------------- clinical
CIVIC_FILES = ["ClinicalEvidenceSummaries", "AssertionSummaries",
               "MolecularProfileSummaries", "VariantSummaries"]


def fetch_civic(c: Ctx) -> None:
    out = EXT / "clinical" / "civic"
    today = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d")
    for name in CIVIC_FILES:
        url = f"https://civicdb.org/downloads/nightly/nightly-{name}.tsv"
        dest = out / f"nightly-{name}.tsv"
        if c.have(dest):
            print(f"  [skip] {c.rel(dest)}")
            continue

        def ver(r, today=today):
            lm = r.headers.get("Last-Modified", "")
            return f"nightly retrieved {today}" + (f" (Last-Modified: {lm})" if lm else "")
        c.download(url, dest, "civic", f"nightly retrieved {today}", "CC0 1.0",
                   headers_version=ver)


def fetch_moalmanac(c: Ctx) -> None:
    out = EXT / "clinical" / "moalmanac"
    try:
        r = c.get("https://zenodo.org/api/records/17089309/versions/latest")
        rec = r.json()
        files = rec["files"]
        version = rec["metadata"].get("version", "?")
        doi = rec.get("doi", "")
        pub = rec["metadata"].get("publication_date", "")
        lic = ("Repo LICENSE and Zenodo record state GPL-2.0 (id %s); the database being CC BY 4.0 "
               "is NOT confirmed in the release files (only in earlier literature notes) - verify before reuse"
               % rec["metadata"].get("license", {}).get("id", ""))
        vtxt = f"{version} (published {pub}); Zenodo DOI {doi}; concept DOI 10.5281/zenodo.17089309"
    except Exception as e:  # noqa: BLE001
        c.block("moalmanac", "https://zenodo.org/api/records/17089309", str(e)[:300], "CC BY 4.0")
        return
    for f in files:
        dest = out / Path(f["key"]).name
        url = f["links"]["self"]
        if c.have(dest):
            print(f"  [skip] {c.rel(dest)}")
        else:
            got = c.download(url, dest, "moalmanac", vtxt, lic)
            if got is None:
                continue
        if dest.suffix == ".zip" and dest.exists():
            with zipfile.ZipFile(dest) as z:
                for info in z.infolist():
                    n = info.filename
                    # repo-root files only (database json, README, LICENSE, LOG); the
                    # per-record dereferenced/referenced trees stay inside the zip
                    if info.is_dir() or n.count("/") != 1 or not re.search(
                            r"\.(json|md)$|LICENSE$", n, re.I):
                        continue
                    if info.file_size > 50 * 1024 * 1024:
                        continue
                    xp = out / "extracted" / n
                    if not xp.exists():
                        xp.parent.mkdir(parents=True, exist_ok=True)
                        xp.write_bytes(z.read(info))
                    if not c.have(xp):
                        c.record(xp, f"{url}#{n}", "moalmanac", vtxt, lic,
                                 note="member extracted from release zip")


class TableParser(HTMLParser):
    """Minimal tolerant HTML table extractor (cells -> text and hrefs)."""

    BLOCK = {"br", "p", "li", "ul", "div"}

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.tables: list[dict] = []
        self._t = None
        self._row = None
        self._cell = None
        self._in_head = False

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag == "table":
            self._t = {"head": [], "rows": [], "row_cells_raw": 0}
            self.tables.append(self._t)
        elif self._t is None:
            return
        elif tag == "thead":
            self._in_head = True
        elif tag == "tr":
            self._row = []
        elif tag in ("td", "th") and self._row is not None:
            self._cell = {"text": [], "links": []}
        elif self._cell is not None:
            if tag in self.BLOCK:
                self._cell["text"].append("\n")
            if tag == "a" and a.get("href"):
                self._cell["links"].append(a["href"])

    def handle_endtag(self, tag):
        if self._t is None:
            return
        if tag in ("td", "th") and self._cell is not None and self._row is not None:
            txt = "".join(self._cell["text"])
            parts = [re.sub(r"\s+", " ", p).strip() for p in txt.split("\n")]
            self._row.append({"text": " | ".join(p for p in parts if p),
                              "links": self._cell["links"]})
            self._cell = None
        elif tag == "tr" and self._row is not None:
            if self._row:
                self._t["head" if self._in_head else "rows"].append(self._row)
            self._row = None
        elif tag == "thead":
            self._in_head = False
        elif tag == "table":
            self._t = None
        elif self._cell is not None and tag in self.BLOCK:
            self._cell["text"].append("\n")

    def handle_data(self, data):
        if self._cell is not None:
            self._cell["text"].append(data)


def decode_html(b: bytes) -> str:
    try:
        return b.decode("utf-8")
    except UnicodeDecodeError:
        return b.decode("cp1252", errors="replace")


DATE_RE = re.compile(r"\(\s*(\d{1,2}/\s*\d{1,2}/\s*\d{4})\s*\)")


def parse_fda_cdx(html: str) -> tuple[list[str], list[dict], dict]:
    """Parse the main FDA CDx table. Returns (columns, rows, report).

    Robust to row misalignment: rows whose cell count differs from the header are not
    silently dropped; they are kept (padded/truncated) and flagged. Row count is
    cross-checked against the number of <tr> elements in the first table.
    """
    p = TableParser()
    p.feed(html)
    main = p.tables[0]
    header = [re.sub(r"\s+", " ", h["text"].replace(" | ", " ")).strip() for h in main["head"][0]]
    ncol = len(header)
    rows, report = [], {"n_header_cols": ncol, "n_rows_parsed": 0, "suspicious": []}
    raw_tr = len(re.findall(r"<tr\b", html.split("</table>")[0]))  # head row + body rows
    report["n_tr_in_first_table"] = raw_tr
    seen = {}
    for i, r in enumerate(main["rows"]):
        flags = []
        cells = [x["text"] for x in r]
        links = [x["links"] for x in r]
        if len(cells) != ncol:
            flags.append(f"cell_count={len(cells)}!={ncol}")
            cells = (cells + [""] * ncol)[:ncol]
            links = (links + [[]] * ncol)[:ncol]
        rec = dict(zip(header, cells))
        last = cells[-1]
        m = DATE_RE.findall(last)
        if not m:
            flags.append("no_approval_date_in_last_col")
        if not all(cells[:4]):
            flags.append("empty_core_field")
        if len(cells[0]) < 4 or len(cells[0]) > 250:
            flags.append("odd_device_name_length")
        if cells[3] and re.search(r"\(\d{1,2}/\d{1,2}/\d{4}\)|NDA|BLA", cells[3]):
            flags.append("biomarker_field_looks_like_drug_or_date")
        if cells[2] and not re.search(r"NDA|BLA|ANDA|N/A|\d{5,}", cells[2] + " ".join(links[2])):
            flags.append("drug_field_without_NDA_BLA_number")
        pma = re.findall(r"\b([PKH]\d{5,7}(?:/S\d+)?|K\d{6}|DEN\d{6})\b", last)
        rec["pma_ids"] = ";".join(pma)
        rec["approval_dates"] = ";".join(m)
        rec["drug_application_numbers"] = ";".join(
            re.findall(r"\b((?:NDA|BLA)\s?\d{5,6})\b", cells[2]))
        key = tuple(cells)
        if key in seen:
            flags.append(f"exact_duplicate_of_row_{seen[key]}")
        seen[key] = i
        rec["row_index"] = i
        rec["suspicious"] = ";".join(flags)
        if flags:
            report["suspicious"].append((i, flags, cells[0][:60]))
        rows.append(rec)
    report["n_rows_parsed"] = len(rows)
    report["row_count_matches_tr"] = (len(rows) == raw_tr - 1)
    return header, rows, report


def parse_fda_group_table(html: str) -> tuple[list[str], list[dict]]:
    p = TableParser()
    p.feed(html)
    if len(p.tables) < 2:
        return [], []
    t = p.tables[1]
    header = [re.sub(r"\s+", " ", h["text"].replace(" | ", " ")).strip() for h in t["head"][0]]
    rows = []
    for r in t["rows"]:
        cells = ([x["text"] for x in r] + [""] * len(header))[: len(header)]
        rows.append(dict(zip(header, cells)))
    return header, rows


def fetch_fda_cdx(c: Ctx) -> None:
    out = EXT / "clinical" / "fda_cdx"
    url = ("https://www.fda.gov/medical-devices/in-vitro-diagnostics/"
           "list-cleared-or-approved-companion-diagnostic-devices-in-vitro-and-imaging-tools")
    dest = out / "fda_cdx.html"
    if not c.have(dest):
        out.mkdir(parents=True, exist_ok=True)
        try:
            r = c.get(url)
        except Exception as e:  # noqa: BLE001
            c.block("fda_cdx", url, str(e)[:300], "US government public domain")
            return
        dest.write_bytes(r.content)
        html = decode_html(r.content)
        m = re.search(r"Content current as of:?</h2>\s*<p[^>]*>\s*<time datetime=\"([^\"]+)\"", html)
        cur = m.group(1)[:10] if m else "unknown"
        c.record(dest, url, "fda_cdx", f"content current as of {cur}", "US government public domain")
    else:
        print(f"  [skip] {c.rel(dest)}")
    html = decode_html(dest.read_bytes())
    entry = c.m["files"][c.rel(dest)]
    header, rows, rep = parse_fda_cdx(html)
    pcsv = out / "fda_cdx_parsed.csv"
    cols = header + ["pma_ids", "approval_dates", "drug_application_numbers", "row_index", "suspicious"]
    with open(pcsv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    h2, rows2 = parse_fda_group_table(html)
    p2 = out / "fda_cdx_group_indications_parsed.csv"
    if rows2:
        with open(p2, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=h2)
            w.writeheader()
            w.writerows(rows2)
    print(f"  parsed {rep['n_rows_parsed']} rows x {rep['n_header_cols']} cols; "
          f"<tr> in table = {rep['n_tr_in_first_table']} (match: {rep['row_count_matches_tr']}); "
          f"suspicious rows: {len(rep['suspicious'])}")
    for i, fl, name in rep["suspicious"][:25]:
        print(f"    row {i}: {fl} | {name}")
    note = (f"derived from {entry['size']} B raw html sha256 {entry['sha256'][:12]}; "
            f"rows={rep['n_rows_parsed']}; suspicious={len(rep['suspicious'])}")
    for pth in (pcsv, p2):
        if pth.exists():
            c.record(pth, entry["url"], "fda_cdx", entry["version"], entry["licence"],
                     retrieved=entry["retrieved_utc"], note=note)
    c.m.setdefault("parse_reports", {})["fda_cdx"] = {
        k: v for k, v in rep.items() if k != "suspicious"} | {"n_suspicious": len(rep["suspicious"])}


def fetch_oncokb(c: Ctx) -> None:
    out = EXT / "clinical" / "oncokb"
    base = "https://www.oncokb.org/api/v1"
    lic = ("OncoKB terms: free academic non-commercial use; academic/commercial terms of use "
           "UNCONFIRMED here (could not read terms page); do not redistribute or commit")
    ver = "unknown"
    try:
        info = c.get(base + "/info").json()
        ver = f"dataVersion={info.get('dataVersion', {})}; oncotree={info.get('oncoTreeVersion', '')}"
    except Exception as e:  # noqa: BLE001
        ver = f"info endpoint failed: {str(e)[:100]}"
    for name in ("cancerGeneList.txt", "allCuratedGenes.txt"):
        c.download(f"{base}/utils/{name}", out / name, "oncokb", ver, lic,
                   note="public endpoint, no token")


# ----------------------------------------------------------------------------- screens
def fetch_parrish(c: Ctx) -> None:
    out = EXT / "combinatorial" / "parrish2021"
    base = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178179/suppl/"
    lic = "Article CC BY-NC-ND (GEO data: public, NCBI GEO terms); check before redistributing derived tables"
    try:
        listing = c.get(base).text
    except Exception as e:  # noqa: BLE001
        c.block("parrish2021", base, str(e)[:300], lic)
        return
    names = sorted(set(re.findall(r'href="(GSE178179_[^"?/]+)"', listing)))
    raw = re.compile(r"\.(fastq|fq|sra|bam|cram)(\.gz)?$|RAW\.tar$", re.I)
    total = 0
    ver = "GEO GSE178179, series last updated 2021-10-27 (public 2021-08-31)"
    for n in names:
        if raw.search(n):
            print(f"  [skip raw] {n}")
            continue
        url = base + n
        try:
            size = int(c.request("HEAD", url).headers.get("Content-Length", 0))
        except Exception:  # noqa: BLE001
            size = 0
        if total + size > CAP_BYTES:
            print(f"  [skip cap] {n} ({size:,} B)")
            continue
        total += size
        c.download(url, out / n, "parrish2021", ver, lic)
    # series + sample metadata (small, useful for column/sample interpretation)
    soft = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178179/miniml/"
            "GSE178179_family.xml.tgz")
    c.download(soft, out / "GSE178179_family.xml.tgz", "parrish2021", ver, lic,
               note="sample metadata (MINiML)")
    # Processed GI scores are in the article supplement (PMC8534300 / Cell Rep). The Europe PMC
    # bundle holds figures only, and PMC serves a bot-challenge page for the bin/ files.
    pbase = "https://pmc.ncbi.nlm.nih.gov/articles/instance/8534300/bin/NIHMS1737332-supplement-{}.xlsx"
    ok = 0
    for i in range(2, 9):
        dest = out / f"NIHMS1737332-supplement-{i}.xlsx"
        if c.have(dest):
            ok += 1
            continue
        try:
            r = c.get(pbase.format(i))
        except Exception:  # noqa: BLE001
            continue
        if r.content[:2] != b"PK":      # HTML challenge page, not an xlsx
            continue
        dest.write_bytes(r.content)
        c.record(dest, pbase.format(i), "parrish2021", "Cell Rep 2021 supplement (PMC8534300)",
                 "Article CC BY-NC-ND")
        ok += 1
    if ok == 0:
        c.block("parrish2021", "https://pmc.ncbi.nlm.nih.gov/articles/PMC8534300/#supplementary-materials",
                "supplementary xlsx (processed GI scores, NIHMS1737332-supplement-2..8) return a bot "
                "challenge page; request manual download or Cell Rep mmc files", lic)
    # authors' analysis code (defines GI score computation)
    zurl = ("https://zenodo.org/api/records/5081113/files/bergerbio/"
            "pgRNA_sequencing_analysis-v1.0.0.zip/content")
    c.download(zurl, out / "pgRNA_sequencing_analysis-v1.0.0.zip", "parrish2021",
               "Zenodo DOI 10.5281/zenodo.5081113 v1.0.0", "MIT (code)")


def fetch_pmc_supp(c: Ctx, pmcid: str, source: str, ver: str, lic: str, subdir: Path,
                   ext_ok=r"\.(xlsx?|csv|tsv|txt|docx|pdf)$") -> None:
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/supplementaryFiles"
    # skip if every wanted file is already present
    existing = [p for p in subdir.glob("*") if p.is_file()] if subdir.exists() else []
    if existing and all(c.have(p) for p in existing):
        for p in existing:
            print(f"  [skip] {c.rel(p)}")
        return
    try:
        r = c.get(url)
        z = zipfile.ZipFile(io.BytesIO(r.content))
    except Exception as e:  # noqa: BLE001
        c.block(source, url, str(e)[:300], lic)
        return
    t = utcnow()
    subdir.mkdir(parents=True, exist_ok=True)
    for info in z.infolist():
        n = info.filename
        if not re.search(ext_ok, n, re.I) or re.search(r"_Fig\d|Equ\d|_Tab\d", n):
            continue
        dest = subdir / n
        dest.write_bytes(z.read(info))
        c.record(dest, f"{url}#{n}", source, ver + f"; Europe PMC supplementaryFiles for {pmcid}",
                 lic, retrieved=t)


def fetch_dede(c: Ctx) -> None:
    fetch_pmc_supp(c, "PMC7558751", "dede2020", "Dede 2020 Genome Biol 21:262",
                   "Article CC BY 4.0; data CC0 waiver", EXT / "combinatorial" / "dede2020")


def fetch_thompson(c: Ctx) -> None:
    fetch_pmc_supp(c, "PMC7910459", "thompson2021", "Thompson 2021 Nat Commun 12:1302",
                   "CC BY 4.0", EXT / "combinatorial" / "thompson2021",
                   ext_ok=r"\.(xlsx?|csv|tsv|txt|pdf)$")


def fetch_esmaeili(c: Ctx) -> None:
    out = EXT / "combinatorial" / "esmaeili2024"
    api = "https://api.figshare.com/v2/articles/24243832"
    try:
        d = c.get(api).json()
    except Exception as e:  # noqa: BLE001
        c.block("esmaeili2024", api, str(e)[:300], "CC BY 4.0")
        return
    ver = f"figshare DOI {d.get('doi')} version {d.get('version')}, modified {d.get('modified_date')}"
    lic = (d.get("license") or {}).get("name", "CC BY 4.0")
    total = 0
    for f in sorted(d["files"], key=lambda x: x["name"]):
        n = f["name"]
        if n.endswith(".ipynb"):
            continue                       # code, not data
        if not re.search(r"\.(xlsx|txt|csv|tsv|zip)$", n, re.I):
            continue
        if total + f["size"] > CAP_BYTES:
            print(f"  [skip cap] {n}")
            continue
        total += f["size"]
        dest = out / n
        if c.have(dest):
            print(f"  [skip] {c.rel(dest)}")
            continue
        got = c.download(f["download_url"], dest, "esmaeili2024", ver, lic,
                         note=f"figshare file id {f['id']}, md5 {f.get('computed_md5')}")
        if got and f.get("computed_md5"):
            import hashlib as _h
            if _h.md5(got.read_bytes()).hexdigest() != f["computed_md5"]:
                print(f"  WARNING md5 mismatch for {n}")


def fetch_flister(c: Ctx) -> None:
    """bioRxiv 10.1101/2024.07.16.603642 (v3 current). The brief called it 'Klingbeil';
    the bioRxiv record lists Flister et al. (Pfizer)."""
    out = EXT / "combinatorial" / "klingbeil2024_biorxiv"
    page = "https://www.biorxiv.org/content/10.1101/2024.07.16.603642v3.supplementary-material"
    lic = "bioRxiv license cc_no (all rights reserved / no reuse licence); personal analysis only, do not redistribute"
    try:
        html = c.get(page).text
    except Exception as e:  # noqa: BLE001
        c.block("flister2024", page, str(e)[:300], lic)
        return
    links = sorted(set(re.findall(
        r'href="(https://www\.biorxiv\.org/content/biorxiv/early/[^"]+/DC\d+/embed/media-\d+\.[a-z]+\?download=true)"',
        html)), key=lambda u: int(re.search(r"media-(\d+)", u).group(1)))
    if not links:
        c.block("flister2024", page, "no supplementary links found on page", lic)
        return
    ver = "bioRxiv 10.1101/2024.07.16.603642 v3 (2025-01-16); title 'Genetic Mapping of All Human Paralog Pairs...'"
    for u in links:
        n = re.search(r"(media-\d+\.[a-z]+)", u).group(1)
        c.download(u, out / n, "flister2024", ver, lic,
                   note="media-N = Supplementary Table S(N); verify against supplement page")


SOURCES = {
    "civic": fetch_civic, "moalmanac": fetch_moalmanac, "fda_cdx": fetch_fda_cdx,
    "oncokb": fetch_oncokb, "parrish2021": fetch_parrish, "dede2020": fetch_dede,
    "esmaeili2024": fetch_esmaeili, "klingbeil2024": fetch_flister,
    "thompson2021": fetch_thompson,
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", action="append", help="source name(s); repeat or comma-separate")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    if a.list:
        print("\n".join(SOURCES))
        return 0
    names = [x for o in (a.only or []) for x in o.split(",")] or list(SOURCES)
    bad = [n for n in names if n not in SOURCES]
    if bad:
        print("unknown source(s):", bad, "choose from", list(SOURCES))
        return 2
    c = Ctx()
    for n in names:
        print(f"== {n}")
        try:
            SOURCES[n](c)
        except Exception as e:  # noqa: BLE001
            c.block(n, "(unhandled)", repr(e)[:300])
        save_manifest(c.m)
    print(f"manifest: {MANIFEST} ({len(c.m['files'])} files, {len(c.m['blocked'])} blocked)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
