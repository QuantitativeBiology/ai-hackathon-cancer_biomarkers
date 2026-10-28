"""Test setup: put src/ on sys.path and install the Part 1 write blocker for the whole session."""
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from biomark import p1 as _p1  # noqa: E402

_p1.install_p1_write_guard()
