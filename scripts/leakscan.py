"""Leakscan (DESIGN §13.4): AST scan of ``src/`` and ``scripts/`` for hard-coded answers.

    python -X utf8 scripts/leakscan.py [--json out.json]

Hits:
  (i)  ``gene``: a string constant exactly equal to an HGNC symbol of >= 3 characters from Part 1 ``gene_map``
       (``biomark.p1.gene_map()``, ../project/cache/gene_map.parquet, read-only), anywhere in src/ or scripts/;
       ``ach``: a string constant containing ``ACH-\\d{6}``, anywhere in src/ or scripts/;
  (ii) ``num``: a numeric constant (``-x`` counts as one constant) outside {0, 1, -1, 2, 0.5, 100} in
       ``src/biomark/``, except the params loader ``src/biomark/params.py`` and lines carrying
       ``# literal-ok: <reason>`` (non-empty reason).
The ``# literal-ok`` marker exempts (ii) only; (i) has no exemption. Booleans are not numbers. Exit 1 on any hit.
Every ``literal-ok`` comment in the scanned tree is listed (``EXEMPT`` lines: what it suppressed, its reason, and
whether it is well formed), so an exemption is never a silent bypass.
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


def scan(root: Path = ROOT, symbols: set | None = None, dirs=SCAN_DIRS, exemptions: list | None = None) -> list:
    """Hits; if ``exemptions`` (a list) is given, every literal-ok marker is appended to it."""
    root = Path(root)
    symbols = gene_symbols() if symbols is None else set(symbols)
    hits = []
    for d in dirs:
        base = root / d
        if not base.exists():
            continue
        for p in sorted(base.rglob("*.py")):
            if "__pycache__" in p.parts:
                continue
            hits += scan_file(p, p.relative_to(root), symbols, exemptions)
    return hits


def render(hits: list, exemptions: list | None = None) -> str:
    L = []
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
    L.append("RESULT: " + ("FAIL" if hits else "PASS"))
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="DESIGN §13.4 leakscan")
    ap.add_argument("--root", type=Path, default=ROOT)
    ap.add_argument("--json", type=Path, default=None)
    a = ap.parse_args(argv)
    ex = []
    hits = scan(a.root, exemptions=ex)
    print(render(hits, ex))
    if a.json:
        a.json.write_text(json.dumps({"hits": hits, "literal_ok": ex}, indent=1), encoding="utf-8")
    return 1 if hits else 0


if __name__ == "__main__":
    sys.exit(main())
