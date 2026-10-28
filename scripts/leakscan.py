"""Leakscan (DESIGN §13.4): AST scan of ``src/`` and ``scripts/`` for hard-coded answers.

    python -X utf8 scripts/leakscan.py [--json out.json] [--allow configs/leakscan_allow.json]

Hits:
  (i)  ``gene``: a string constant exactly equal to an HGNC symbol of >= 3 characters from Part 1 ``gene_map``
       (``biomark.p1.gene_map()``, ../project/cache/gene_map.parquet, read-only), anywhere in src/ or scripts/;
       ``ach``: a string constant containing ``ACH-\\d{6}``, anywhere in src/ or scripts/;
  (ii) ``num``: a numeric constant (``-x`` counts as one constant) outside {0, 1, -1, 2, 0.5, 100} in
       ``src/biomark/``, except the params loader ``src/biomark/params.py`` and lines carrying
       ``# literal-ok: <reason>`` (non-empty reason).
The ``# literal-ok`` marker exempts (ii) only. A (i) ``gene`` hit is exempt only through the explicit allowlist
``configs/leakscan_allow.json`` ({"allow": [{file, literal, reason}]}): an entry covers that exact literal in that
exact file (repo-relative posix path) and nowhere else; ``ach`` hits have no exemption. Booleans are not numbers.
Exit 1 on any hit or on an invalid allowlist entry (missing field, empty reason, ``ACH-`` literal, file absent).
Every exemption is listed (``EXEMPT`` lines): each ``literal-ok`` comment with what it suppressed, its reason and
whether it is well formed, and each allowlisted occurrence with its reason; unused allowlist entries are reported.
"""
from __future__ import annotations

import argparse
import ast
import json
import io
import re
import sys
import tokenize
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ALLOWED_NUMBERS = frozenset({0, 1, -1, 2, 0.5, 100})
MIN_SYMBOL_LEN = 3
ACH_RE = re.compile(r"ACH-\d{6}")
LITERAL_OK_RE = re.compile(r"#\s*literal-ok:\s*\S")
LITERAL_OK_ANY_RE = re.compile(r"literal-ok", re.IGNORECASE)
PARAMS_LOADER = Path("src") / "biomark" / "params.py"
NUMERIC_SCOPE = Path("src") / "biomark"
SCAN_DIRS = ("src", "scripts")
ALLOW_FILE = Path("configs") / "leakscan_allow.json"


def load_allow(path: Path) -> tuple[list, list]:
    """(entries, errors) from an allowlist file; a missing file is an empty allowlist."""
    path = Path(path)
    if not path.exists():
        return [], []
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except Exception as e:  # noqa: BLE001
        return [], [f"{path}: unreadable JSON ({e})"]
    items = raw.get("allow", []) if isinstance(raw, dict) else raw
    return validate_allow(items)


def validate_allow(items) -> tuple[list, list]:
    ok, errs = [], []
    for i, e in enumerate(items or []):
        if not isinstance(e, dict) or not all(isinstance(e.get(k), str) for k in ("file", "literal", "reason")):
            errs.append(f"allowlist entry {i}: needs string fields file, literal, reason: {e!r}")
            continue
        if not e["reason"].strip():
            errs.append(f"allowlist entry {i} ({e['file']}, {e['literal']!r}): empty reason")
            continue
        if ACH_RE.search(e["literal"]):
            errs.append(f"allowlist entry {i} ({e['file']}, {e['literal']!r}): ACH- ids cannot be allowlisted")
            continue
        ok.append(e)
    return ok, errs


def gene_symbols() -> set:
    """HGNC symbols of >= 3 characters from Part 1 gene_map (read-only)."""
    sys.path.insert(0, str(ROOT / "src"))
    from biomark import p1
    sym = p1.gene_map()["symbol"].dropna().astype(str)
    return {s for s in sym if len(s) >= MIN_SYMBOL_LEN}


def _is_num(v) -> bool:
    return isinstance(v, (int, float, complex)) and not isinstance(v, bool)


class _Visitor(ast.NodeVisitor):
    def __init__(self, symbols: set, numeric: bool):
        self.symbols, self.numeric = symbols, numeric
        self.hits = []
        self._neg = set()

    def visit_UnaryOp(self, node):
        if isinstance(node.op, (ast.USub, ast.UAdd)) and isinstance(node.operand, ast.Constant) \
                and _is_num(node.operand.value):
            self._neg.add(id(node.operand))
            v = -node.operand.value if isinstance(node.op, ast.USub) else node.operand.value
            self._num(node, v)
            return
        self.generic_visit(node)

    def visit_Constant(self, node):
        v = node.value
        if isinstance(v, str):
            if v in self.symbols:
                self.hits.append(("gene", node.lineno, node.col_offset, v))
            elif ACH_RE.search(v):
                self.hits.append(("ach", node.lineno, node.col_offset, ACH_RE.search(v).group(0)))
        elif _is_num(v) and id(node) not in self._neg:
            self._num(node, v)

    def _num(self, node, v):
        if not self.numeric:
            return
        if isinstance(v, complex) or v not in ALLOWED_NUMBERS:
            self.hits.append(("num", node.lineno, node.col_offset, repr(v)))


def literal_ok_markers(src: str, rel: Path) -> list:
    """Every comment mentioning ``literal-ok`` (well-formed or not), so exemptions stay visible."""
    out = []
    try:
        for tok in tokenize.generate_tokens(io.StringIO(src).readline):
            if tok.type == tokenize.COMMENT and LITERAL_OK_ANY_RE.search(tok.string):
                m = re.search(r"literal-ok:\s*(.*)$", tok.string)
                out.append({"file": rel.as_posix(), "line": tok.start[0], "comment": tok.string.strip(),
                            "reason": m.group(1).strip() if m else "",
                            "well_formed": bool(LITERAL_OK_RE.search(tok.string)), "suppressed": []})
    except (tokenize.TokenError, IndentationError, SyntaxError):
        pass
    return out


def scan_file(path: Path, rel: Path, symbols: set, exemptions: list | None = None) -> list:
    src = path.read_text(encoding="utf-8")
    lines = src.splitlines()
    numeric = rel.parts[: len(NUMERIC_SCOPE.parts)] == NUMERIC_SCOPE.parts and rel != PARAMS_LOADER
    try:
        tree = ast.parse(src, filename=str(rel))
    except SyntaxError as e:
        return [{"file": rel.as_posix(), "line": e.lineno or 0, "col": 0, "kind": "syntax", "value": str(e),
                 "source": ""}]
    v = _Visitor(symbols, numeric)
    v.visit(tree)
    marks = {m["line"]: m for m in literal_ok_markers(src, rel)}
    out = []
    for kind, ln, col, val in v.hits:
        text = lines[ln - 1] if 0 < ln <= len(lines) else ""
        if kind == "num" and ln in marks and marks[ln]["well_formed"]:
            marks[ln]["suppressed"].append(val)
            continue
        out.append({"file": rel.as_posix(), "line": ln, "col": col, "kind": kind, "value": val,
                    "source": text.strip(), "literal_ok_line": ln in marks})
    if exemptions is not None:
        for m in marks.values():
            m["source"] = (lines[m["line"] - 1].strip() if 0 < m["line"] <= len(lines) else "")
            exemptions.append(m)
    return out


def scan(root: Path = ROOT, symbols: set | None = None, dirs=SCAN_DIRS, exemptions: list | None = None,
         allow: list | None = None, allowed: list | None = None, errors: list | None = None) -> list:
    """Hits. If ``exemptions`` (a list) is given, every literal-ok marker is appended to it. ``allow`` = allowlist
    entries (default: ``<root>/configs/leakscan_allow.json``); allowlisted ``gene`` hits go to ``allowed`` (with the
    entry's reason) instead of the hits, followed by unused entries (``unused: True``); invalid entries -> ``errors``."""
    root = Path(root)
    symbols = gene_symbols() if symbols is None else set(symbols)
    if allow is None:
        allow, errs = load_allow(root / ALLOW_FILE)
    else:
        allow, errs = validate_allow(allow)
    errs += [f"allowlist entry ({e['file']}, {e['literal']!r}): file does not exist" for e in allow
             if not (root / e["file"]).exists()]
    if errors is not None:
        errors += errs
    hits = []
    for d in dirs:
        base = root / d
        if not base.exists():
            continue
        for p in sorted(base.rglob("*.py")):
            if "__pycache__" in p.parts:
                continue
            hits += scan_file(p, p.relative_to(root), symbols, exemptions)
    index = {(e["file"], e["literal"]): e for e in allow}
    used, kept = set(), []
    for h in hits:
        e = index.get((h["file"], h["value"])) if h["kind"] == "gene" else None
        if e is None:
            kept.append(h)
            continue
        used.add((e["file"], e["literal"]))
        if allowed is not None:
            allowed.append({**h, "reason": e["reason"]})
    if allowed is not None:
        allowed += [{"file": f, "literal": lit, "reason": e["reason"], "unused": True}
                    for (f, lit), e in index.items() if (f, lit) not in used]
    return kept


def render(hits: list, exemptions: list | None = None, allowed: list | None = None,
           errors: list | None = None) -> str:
    L = []
    for a in allowed or []:
        if a.get("unused"):
            L.append(f"WARN: allowlist entry ({a['file']}, {a['literal']!r}) matched nothing (stale?)")
        else:
            L.append(f"EXEMPT {a['file']}:{a['line']}:{a['col']}: [gene] {a['value']} allowlisted  "
                     f"reason: {a['reason']}    | {a['source']}")
    if allowed is not None:
        L.append(f"allowlisted string hits: {sum(1 for a in allowed if not a.get('unused'))}")
    for m in exemptions or []:
        what = (f"suppressed {m['suppressed']}" if m["suppressed"] else "suppresses nothing")
        form = "" if m["well_formed"] else "  MALFORMED (no reason; exempts nothing)"
        L.append(f"EXEMPT {m['file']}:{m['line']}: {what}{form}  reason: {m['reason'] or '-'}    | {m['source']}")
    if exemptions is not None:
        n_sup = sum(len(m["suppressed"]) for m in exemptions)
        L.append(f"literal-ok markers: {len(exemptions)} ({n_sup} numeric constant(s) suppressed)")
    L += [f"{h['file']}:{h['line']}:{h['col']}: [{h['kind']}] {h['value']}"
         f"{'  (literal-ok on line: not an exemption for (i))' if h.get('literal_ok_line') else ''}    | {h['source']}"
         for h in hits]
    by = {}
    for h in hits:
        by.setdefault((h["file"], h["kind"]), 0)
        by[(h["file"], h["kind"])] += 1
    L.append(f"leakscan: {len(hits)} hit(s)")
    for (f, k), n in sorted(by.items()):
        L.append(f"  {f} [{k}]: {n}")
    for e in errors or []:
        L.append(f"FAIL: {e}")
    L.append("RESULT: " + ("FAIL" if hits or errors else "PASS"))
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="DESIGN §13.4 leakscan")
    ap.add_argument("--root", type=Path, default=ROOT)
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--allow", type=Path, default=None, help="allowlist (default <root>/configs/leakscan_allow.json)")
    a = ap.parse_args(argv)
    ex, allowed, errors = [], [], []
    allow = None
    if a.allow is not None:
        allow, errs = load_allow(a.allow)
        errors += errs
    hits = scan(a.root, exemptions=ex, allow=allow, allowed=allowed, errors=errors)
    print(render(hits, ex, allowed, errors))
    if a.json:
        a.json.write_text(json.dumps({"hits": hits, "literal_ok": ex, "allowlisted": allowed, "errors": errors},
                                     indent=1), encoding="utf-8")
    return 1 if hits or errors else 0


if __name__ == "__main__":
    sys.exit(main())
