"""Check every number in REPORT.md and docs/effective_rules.md against its cited source (DESIGN v3.2 §13.5).

    python -B scripts/check_report_numbers.py                       # REPORT.md and docs/effective_rules.md (those present)
    python -B scripts/check_report_numbers.py REPORT.md --strict    # final gate
    python -B scripts/check_report_numbers.py --json out.json       # also write every checked number

Ported from Part 1 (scripts/check_report_numbers.py there; read-only). Two citation forms:

* **Inline citation** (the writer's form): a repo-relative source in parentheses, e.g. ``... failed 21/39
  (results/gates/G_pc_attempt2.json)`` or ``(results/a.json; configs/design_params.json)``. Sources are
  ``results/``, ``configs/`` or ``docs/`` files ending in ``.json`` or ``.csv``, or the text sources DESIGN.md,
  RESEARCH_TREE.md and results/deviations.csv (a number cited to these must occur there verbatim, as a whole token;
  a qualifier may follow the path, e.g. ``(DESIGN.md §6.1; RESEARCH_TREE.md D15)``; they need no registry row).
  A citation covers every number of its
  **segment**: the sentence it stands in (prose, list items), or the table row; a table row without a citation
  inherits the citations of the table's header row and of the line just above the table (its caption). A covered
  number passes if it appears in one of the segment's cited files at the printed precision:
  |printed − value| <= 0.5 x 10^-decimals (+1e-12), searched over every numeric value of the file (JSON leaves and
  numbers inside JSON strings such as "21/39"; every CSV cell). A printed percentage also matches value x 100 or
  value. Signs are checked (U+2212 is a minus); a hyphen between two numbers is a range, not a sign.
* **Pinned tag** (Part 1 form, exact path): ``<!--n:SOURCE#PATH[|pct|abs|neg]-->0.308`` with optional aliases
  ``<!--n-src:alias=results/x.json-->``; PATH is '/'-separated (dict keys or list indices).

Not numbers: code spans / blocks, HTML comments (incl. ``<!--c:claim_id-->``), link targets, the citation itself and
digits inside identifiers (S1-obs, Q4a, class 3b, §5.5, commit hashes). Structural numbers (section numbers, years,
dates, list / heading numbering, sample sizes and thresholds defined in DESIGN) are exempted by
configs/report_number_exemptions.json: entries {id, regex, reason}; a number inside a match of ``regex`` on its line
is exempt.

Sources (DESIGN §13.5): a cited file must exist, be tracked by git and, under results/, its latest
results/registry.csv row listing it as an output (or a parent directory) must have status ok (``exploratory`` and
``dryrun`` are uncitable). Run ids cited in the text (``YYYYMMDDTHHMMSS_<script>``) must be registry rows with status
ok.

Without ``--strict`` failures are: a cited number not found in its sources, a bad pinned tag, a missing cited file or
run id; source tracking / registry problems are warnings. ``--strict`` also fails on every uncited (non-exempt) number
and on the source warnings. Exit status 1 on any failure.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DOCS = ("REPORT.md", "docs/effective_rules.md")
EXEMPTIONS = "configs/report_number_exemptions.json"
UNCITABLE = ("exploratory", "dryrun")

NUM = r"(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?(?:[eE][+\-]?\d+)?%?"
NUM_RE = re.compile(r"(?<![\w.§#\\])" + NUM + r"(?![\w])")        # sign handled from the preceding characters
NUM_ANY = re.compile(r"[+\-−]?" + NUM)
TAG = re.compile(r"<!--\s*n:(?P<key>[^>]*?)\s*-->[ \t]*(?P<num>[+\-−]?" + NUM + r")?")
SRC = re.compile(r"<!--\s*n-src:(?P<alias>[\w.\-]+)\s*=\s*(?P<path>[^>]*?)\s*-->")
COMMENT = re.compile(r"<!--.*?-->", re.S)
TEXT_SOURCES = ("DESIGN.md", "RESEARCH_TREE.md", "results/deviations.csv")   # number must appear verbatim
CITE_PATH = (r"(?:" + "|".join(re.escape(s) for s in TEXT_SOURCES) +
             r"|(?:results|configs|docs)/[\w./\-]+?\.(?:json|csv))(?![\w./\-])")
CITE_ITEM = CITE_PATH + r"[^;,()]*"                  # a path, then an optional qualifier such as "§6.1" or "D15"
CITE = re.compile(r"\((?P<body>\s*" + CITE_ITEM + r"(?:[;,]\s*" + CITE_ITEM + r")*)\)")
CITE_ONE = re.compile(CITE_PATH)
LINK_TARGET = re.compile(r"\]\([^)]*\)")
RUN_ID = re.compile(r"\b\d{8}T\d{6}_[A-Za-z0-9_]+\b")
SENT_SPLIT = re.compile(r"(?<=[.!?])\s+(?=[A-Z\"'*_\[])")
ABBREV = re.compile(r"(?:\b(?:e\.g|i\.e|vs|cf|al|approx|Fig|Figs|Tab|resp|incl|excl|no|ca)\.|\b[A-Z]\.)$")
LIST_ITEM = re.compile(r"^\s*(?:[-*+]|\d+[.)])\s+")


# ------------------------------------------------------------------------------------------------
# numbers and sources
# ------------------------------------------------------------------------------------------------


def parse_number(s: str):
    """(value, decimals, is_percent) of a printed number (sign, thousands commas, exponent, trailing %)."""
    pct = s.endswith("%")
    t = s.rstrip("%").replace("−", "-").replace(",", "")
    mant, ex = re.split(r"[eE]", t) if re.search(r"[eE]", t) else (t, "0")
    dec = len(mant.split(".")[1]) if "." in mant else 0
    return float(t), max(dec - int(ex), 0), pct


def tolerance(dec: int) -> float:
    return 0.5 * 10 ** (-dec) * (1 + 1e-9) + 1e-12


def resolve(obj, path: str):
    cur = obj
    for part in [p for p in path.split("/") if p != ""]:
        if isinstance(cur, list):
            cur = cur[int(part)]
        elif isinstance(cur, dict):
            if part not in cur:
                raise KeyError(part)
            cur = cur[part]
        else:
            raise KeyError(part)
    return cur


def _numbers_in_str(s: str):
    for m in NUM_ANY.finditer(s):
        try:
            v = parse_number(m.group(0))
        except ValueError:
            continue
        yield v[0] / 100 if v[2] else v[0]
        if v[2]:
            yield v[0]


def _leaves(obj):
    if isinstance(obj, dict):
        for v in obj.values():
            yield from _leaves(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from _leaves(v)
    elif isinstance(obj, bool) or obj is None:
        return
    elif isinstance(obj, (int, float)):
        if math.isfinite(obj):
            yield float(obj)
    elif isinstance(obj, str):
        yield from _numbers_in_str(obj)


class Sources:
    """Cited files (parsed once, numeric values kept as one sorted array per file), git tracking and registry."""

    def __init__(self, root: Path, registry: Path | None = None, tracked: set | None = None):
        self.root = Path(root)
        self.registry = registry or self.root / "results" / "registry.csv"
        self._vals, self._json, self._reg, self._tracked = {}, {}, None, tracked

    def json(self, rel: str):
        if rel not in self._json:
            self._json[rel] = json.loads((self.root / rel).read_text(encoding="utf-8"))
        return self._json[rel]

    def values(self, rel: str) -> np.ndarray:
        if rel not in self._vals:
            p = self.root / rel
            if rel.endswith(".json"):
                vals = list(_leaves(self.json(rel)))
            else:
                vals = []
                with open(p, newline="", encoding="utf-8") as fh:
                    for row in csv.reader(fh):
                        for cell in row:
                            vals.extend(_numbers_in_str(cell))
            self._vals[rel] = np.unique(np.asarray(vals, dtype=np.float64))
        return self._vals[rel]

    def contains_text(self, rel: str, printed: str) -> bool:
        """Text sources (TEXT_SOURCES): the printed number occurs verbatim, as a whole number token (minus written
        as '-' or U+2212; thousands commas optional)."""
        if rel not in self._json:
            raw = (self.root / rel).read_text(encoding="utf-8")
            if rel.endswith(".csv"):                      # cells joined by ' | ' so a delimiter is not a digit group
                raw = "\n".join(" | ".join(row) for row in csv.reader(raw.splitlines()))
            self._json[rel] = raw
        txt = self._json[rel]
        body = printed.lstrip("+-−")
        sign = printed[: len(printed) - len(body)]
        signs = {"-", "−"} if sign in ("-", "−") else {sign}
        bodies = {body, body.replace(",", "")}
        mant = body.rstrip("%").split(".")[0].replace(",", "")
        if "," not in body and mant.isdigit() and len(mant) > 3:
            bodies.add(f"{int(mant):,}" + body[len(mant):])
        alts = "|".join(re.escape(s + b) for s in signs for b in sorted(bodies))
        lead = r"(?<![\d.,])" if not sign else r"(?<![\w.,])"
        return re.search(lead + "(?:" + alts + r")(?![\d]|[.,]\d)", txt) is not None

    def contains(self, rel: str, value: float, tol: float) -> bool:
        a = self.values(rel)
        i = np.searchsorted(a, value - tol, side="left")
        return bool(i < len(a) and a[i] <= value + tol)

    def reg(self) -> list:
        if self._reg is None:
            self._reg = []
            if self.registry.exists():
                with open(self.registry, newline="", encoding="utf-8") as fh:
                    self._reg = list(csv.DictReader(fh))
        return self._reg

    def tracked(self) -> set | None:
        if self._tracked is None:
            try:
                out = subprocess.run(["git", "-C", str(self.root), "ls-files", "-z"], capture_output=True, check=True)
                self._tracked = {x for x in out.stdout.decode("utf-8", "replace").split("\0") if x}
            except (OSError, subprocess.CalledProcessError):
                self._tracked = set()
        return self._tracked

    def source_problems(self, rel: str) -> tuple:
        """(fatal, warning) messages for a cited file: missing is fatal; untracked / registry problems are warnings
        (failures under --strict)."""
        if not (self.root / rel).is_file():
            return f"cited file {rel} does not exist", None
        if rel not in self.tracked():
            return None, f"cited file {rel} is not tracked by git"
        if not rel.startswith("results/") or rel in ("results/registry.csv", *TEXT_SOURCES):
            return None, None
        rows = [r for r in self.reg() if any(rel == o or rel.startswith(o.rstrip("/") + "/")
                                             for o in str(r.get("outputs") or "").split(";") if o)]
        if not rows:
            return None, f"cited file {rel} has no results/registry.csv row"
        last = max(rows, key=lambda r: str(r.get("utc") or ""))
        st = str(last.get("status") or "")
        if st != "ok":
            why = "uncitable" if st in UNCITABLE else "not ok"
            return None, f"cited file {rel}: latest registry row {last.get('run_id')} has status {st!r} ({why})"
        return None, None


def load_exemptions(root: Path, path: str | None = None) -> list:
    f = Path(path) if path else Path(root) / EXEMPTIONS
    if not f.exists():
        f = ROOT / EXEMPTIONS
    if not f.exists():
        return []
    ex = json.loads(f.read_text(encoding="utf-8"))["exemptions"]
    for e in ex:
        if not str(e.get("reason", "")).strip():
            raise ValueError(f"exemption {e.get('id')} has no reason")
        e["_re"] = re.compile(e["regex"], re.M)
    return ex


# ------------------------------------------------------------------------------------------------
# text layout: blanking, segments
# ------------------------------------------------------------------------------------------------


def _blank(text: str, pattern, flags=0) -> str:
    rx = pattern if hasattr(pattern, "sub") else re.compile(pattern, flags)
    return rx.sub(lambda m: re.sub(r"[^\n]", " ", m.group(0)), text)


def clean_text(text: str) -> str:
    """Code blocks, inline code, HTML comments and link targets blanked (positions and line numbers kept)."""
    t = _blank(text, r"```.*?```", re.S)
    t = _blank(t, r"`[^`\n]*`")
    t = _blank(t, COMMENT)
    return _blank(t, LINK_TARGET)


def segments(clean: str) -> list:
    """[(start, end, inherited_citation_spans)] covering the text: one per table row, heading, list item or sentence
    of a paragraph. Table rows inherit the citation-bearing header row and the caption line above the table."""
    lines = clean.split("\n")
    offs, pos = [], 0
    for ln in lines:
        offs.append(pos)
        pos += len(ln) + 1
    segs = []
    i = 0
    while i < len(lines):
        ln = lines[i]
        s = ln.strip()
        if not s:
            i += 1
            continue
        if s.startswith("|"):
            j = i
            while j < len(lines) and lines[j].strip().startswith("|"):
                j += 1
            k = i - 1
            while k >= 0 and not lines[k].strip():
                k -= 1
            inherit = [(offs[i], offs[i] + len(lines[i]))]
            if k >= 0 and not lines[k].strip().startswith("|"):
                inherit.append((offs[k], offs[k] + len(lines[k])))
            for r in range(i, j):
                segs.append((offs[r], offs[r] + len(lines[r]), inherit))
            i = j
            continue
        if s.startswith("#"):
            segs.append((offs[i], offs[i] + len(ln), []))
            i += 1
            continue
        j = i + 1                                    # a paragraph or list item: until blank / table / heading / item
        while j < len(lines):
            t = lines[j].strip()
            if not t or t.startswith(("|", "#")) or LIST_ITEM.match(lines[j]):
                break
            j += 1
        a, b = offs[i], offs[j - 1] + len(lines[j - 1])
        block = clean[a:b]
        cuts = [0] + [m.end() for m in SENT_SPLIT.finditer(block) if not ABBREV.search(block, 0, m.start())]
        cuts += [len(block)]
        for x, y in zip(cuts[:-1], cuts[1:]):
            if block[x:y].strip():
                segs.append((a + x, a + y, []))
        i = j
    return segs


# ------------------------------------------------------------------------------------------------
# checks
# ------------------------------------------------------------------------------------------------


def _signed(clean: str, m) -> str:
    """The printed number with its sign: '-' / '−' / '+' directly before it count as a sign unless they join two
    numbers (a range such as 0.37-0.70) or follow a word character."""
    s = m.start()
    if s > 0 and clean[s - 1] in "+-−":
        prev = clean[s - 2] if s > 1 else " "
        if not (prev.isalnum() or prev in ")%]"):
            return clean[s - 1] + m.group(0)
    return m.group(0)


def check_tags(text: str, src: Sources, line_of) -> list:
    """Part 1 pinned tags <!--n:SOURCE#PATH|flags-->number (aliases <!--n-src:a=path-->)."""
    aliases = {m.group("alias"): m.group("path") for m in SRC.finditer(text)}
    rows = []
    for m in TAG.finditer(text):
        key, num = m.group("key").strip(), m.group("num")
        row = {"line": line_of(m.start()), "kind": "tag", "printed": num, "sources": key, "ok": False, "msg": ""}
        rows.append(row)
        if num is None:
            row["msg"] = "tag not followed by a number"
            continue
        spec, *flags = [x.strip() for x in key.split("|")]
        if "#" not in spec:
            row["msg"] = "key must be SOURCE#PATH"
            continue
        source, path = (x.strip() for x in spec.split("#", 1))
        rel = aliases.get(source, source)
        try:
            v = resolve(src.json(rel), path)
        except (KeyError, IndexError, ValueError, FileNotFoundError, json.JSONDecodeError) as e:
            row["msg"] = f"unresolved: {type(e).__name__}: {e}"
            continue
        if v is None or isinstance(v, (dict, list, str, bool)):
            row["msg"] = f"source value is not a number: {v!r}"[:200]
            continue
        v = float(v)
        bad = [f for f in flags if f not in ("pct", "abs", "neg")]
        if bad:
            row["msg"] = f"unknown flag {bad[0]}"
            continue
        v = v * 100 if "pct" in flags else v
        v = abs(v) if "abs" in flags else v
        v = -v if "neg" in flags else v
        pv, dec, pct = parse_number(num)
        if pct and "pct" not in flags:
            row["msg"] = "printed % but no |pct flag"
            continue
        row["value"] = v
        if not math.isfinite(v):
            row["msg"] = "source value not finite"
        elif abs(pv - v) <= tolerance(dec):
            row["ok"] = True
        else:
            row["msg"] = f"mismatch: source {v:.{dec + 3}f}"
    return rows


def check(report: Path, root: Path = ROOT, strict: bool = False, registry: Path | None = None,
          exemptions: list | None = None, tracked: set | None = None) -> list:
    """One row per checked number / citation / run id: line, kind (cited, uncited, exempt, tag, source, run id),
    printed, sources, ok, msg, warn (source problem that only --strict fails)."""
    text = Path(report).read_text(encoding="utf-8")
    src = Sources(root, registry, tracked)
    ex = load_exemptions(root) if exemptions is None else exemptions

    def line_of(pos):
        return text.count("\n", 0, pos) + 1

    rows = check_tags(text, src, line_of)
    clean = clean_text(text)
    cites = [(m.start(), m.end(), CITE_ONE.findall(m.group("body"))) for m in CITE.finditer(clean)]
    seen = set()
    for a, _, files in cites:
        for rel in files:
            if rel in seen:
                continue
            seen.add(rel)
            fatal, warn = src.source_problems(rel)
            if fatal or warn:
                rows.append({"line": line_of(a), "kind": "source", "printed": "", "sources": rel, "ok": not fatal,
                             "warn": bool(warn) and not fatal, "msg": fatal or warn})
    scan = clean
    for a, b, _ in cites:                                    # the citation text itself holds no report numbers
        scan = scan[:a] + re.sub(r"[^\n]", " ", scan[a:b]) + scan[b:]
    lines = text.split("\n")
    line_start = np.cumsum([0] + [len(x) + 1 for x in lines[:-1]])
    for s0, s1, inherit in segments(clean):
        files = [f for a, b, fs in cites if a >= s0 and b <= s1 for f in fs]
        if not files:
            files = [f for a, b, fs in cites for (x, y) in inherit if a >= x and b <= y for f in fs]
        files = list(dict.fromkeys(files))
        for m in NUM_RE.finditer(scan, s0, s1):
            printed = _signed(scan, m)
            ln = line_of(m.start())
            col0 = m.start() - int(line_start[ln - 1])
            col1 = col0 + len(m.group(0))
            exempt = next((e for e in ex for mm in e["_re"].finditer(lines[ln - 1])
                           if mm.start() <= col0 and col1 <= mm.end()), None)
            row = {"line": ln, "printed": printed, "sources": ";".join(files), "ok": True, "msg": ""}
            if exempt is not None:
                row.update(kind="exempt", msg=f"exempt: {exempt['id']}")
            elif not files:
                row.update(kind="uncited", ok=not strict, msg="uncited number")
            else:
                pv, dec, pct = parse_number(printed)
                tol = tolerance(dec)
                cands = [(pv, tol)] + ([(pv / 100, tol / 100)] if pct else [])
                usable = [f for f in files if (root / f).is_file()]
                hit = next((f for f in usable if f in TEXT_SOURCES and src.contains_text(f, printed)), None)
                hit = hit or next((f for f in usable if f not in TEXT_SOURCES for v, t in cands
                                   if src.contains(f, v, t)), None)
                row.update(kind="cited", ok=hit is not None, matched_in=hit or "",
                           msg="" if hit else ("not found in " + ";".join(files) if usable else "cited file missing"))
            rows.append(row)
    reg = {r.get("run_id"): r for r in src.reg()}
    for m in RUN_ID.finditer(clean):
        r = reg.get(m.group(0))
        msg = ("run id not in results/registry.csv" if r is None else
               "" if str(r.get("status") or "") == "ok" else f"cited run has status {r.get('status')!r}")
        if msg:
            rows.append({"line": line_of(m.start()), "kind": "run id", "printed": m.group(0), "sources": "",
                         "ok": False, "msg": msg})
    if strict:
        for r in rows:
            if r.get("warn"):
                r["ok"] = False
    return rows


def summary(rows: list) -> dict:
    out = {}
    for r in rows:
        k = r["kind"]
        d = out.setdefault(k, {"n": 0, "ok": 0})
        d["n"] += 1
        d["ok"] += int(bool(r["ok"]))
    out["failures"] = sum(1 for r in rows if not r["ok"])
    out["warnings"] = sum(1 for r in rows if r.get("warn") and r["ok"])
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("docs", nargs="*", help="markdown files (default: REPORT.md and docs/effective_rules.md if present)")
    ap.add_argument("--strict", action="store_true")
    ap.add_argument("--root", default=str(ROOT))
    ap.add_argument("--exemptions", default=None, help=f"exemption config (default {EXEMPTIONS})")
    ap.add_argument("--json", default=None, help="write every checked row to this JSON file")
    ap.add_argument("--quiet", action="store_true", help="print the summary only")
    a = ap.parse_args(argv)
    root = Path(a.root)
    docs = [Path(d) for d in a.docs] or [root / d for d in DEFAULT_DOCS if (root / d).exists()]
    if not docs:
        print("no report found (REPORT.md, docs/effective_rules.md)")
        return 1
    ex = load_exemptions(root, a.exemptions)
    status, dump = 0, {}
    for d in docs:
        if not d.exists():
            print(f"{d} not found")
            status = 1
            continue
        rows = check(d, root, strict=a.strict, exemptions=ex)
        s = summary(rows)
        dump[str(d)] = {"summary": s, "rows": rows}
        if not a.quiet:
            for r in rows:
                if not r["ok"] or r.get("warn"):
                    tag = "FAIL" if not r["ok"] else "WARN"
                    print(f"{tag} {d.name}:{r['line']} [{r['kind']}] {r['printed']} {r['sources']}: {r['msg']}")
        counts = ", ".join(f"{k} {v['ok']}/{v['n']}" for k, v in s.items() if isinstance(v, dict))
        print(f"{d.name}: {counts}; {s['failures']} failures, {s['warnings']} warnings"
              f"{' (strict)' if a.strict else ''}")
        status = max(status, int(s["failures"] > 0))
    if a.json:
        Path(a.json).write_text(json.dumps(dump, indent=1, default=str), encoding="utf-8")
    return status


if __name__ == "__main__":
    sys.exit(main())
