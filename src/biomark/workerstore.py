"""File-backed worker payloads for process pools (infrastructure only; no effect on results).

On Windows the spawn start method pickles every initializer argument through a pipe, which fails for the full real-D
matrices (OSError 22 / "pickle data was truncated"). ``stage(payload)`` writes every large ndarray found in the
payload (recursively, including the fields of ``pipeline.Data``) to ``cache/workers/<run_id>/<pool>/a<k>.npy`` and the
rest to ``payload.pkl``; only the directory path crosses the pipe. ``load(path)`` rebuilds the payload in the worker
with ``np.load(mmap_mode='r')`` (read-only, shared through the OS page cache). ``staged(payload)`` is a context
manager that deletes the directory when the pool is done. Blindness: a Data staged with ``include_obs=False`` has no
``_obs`` array, so no observed-effect file exists for pred-arm pools (``files(path)`` lists what was written).
"""
from __future__ import annotations

import contextlib
import os
import pickle
import shutil
import uuid
from pathlib import Path

import numpy as np

from . import paths as P

MIN_BYTES = 1 << 20  # literal-ok: 1 MiB threshold for arrays written as .npy files


class _ArrRef:
    __slots__ = ("name",)

    def __init__(self, name: str):
        self.name = name


class _DataRef:
    __slots__ = ("state",)

    def __init__(self, state: dict):
        self.state = state


def root() -> Path:
    return P.CACHE / "workers"


def _strip(obj, out_dir: Path, written: dict, keyhint: str, min_bytes: int = MIN_BYTES):
    from .pipeline import Data
    if isinstance(obj, np.ndarray) and obj.dtype != object and obj.nbytes >= min_bytes:
        name = f"a{len(written):03d}_{keyhint}.npy"[:80]  # literal-ok: file-name length cap
        np.save(out_dir / name, np.ascontiguousarray(obj), allow_pickle=False)
        written[name] = obj.nbytes
        return _ArrRef(name)
    if isinstance(obj, Data):
        st = {}
        for k, v in obj.__dict__.items():
            if k in ("obs_loader", "access_fn"):          # closures never travel (see Data.__getstate__)
                st[k] = None
                continue
            st[k] = _strip(v, out_dir, written, k.strip("_"), min_bytes)
        return _DataRef(st)
    if isinstance(obj, dict):
        return {k: _strip(v, out_dir, written, f"{keyhint}-{k}" if isinstance(k, str) else keyhint, min_bytes)
                for k, v in obj.items()}
    if isinstance(obj, list):
        return [_strip(v, out_dir, written, keyhint, min_bytes) for v in obj]
    if isinstance(obj, tuple):
        return tuple(_strip(v, out_dir, written, keyhint, min_bytes) for v in obj)
    return obj


def _restore(obj, d: Path):
    from .pipeline import Data
    if isinstance(obj, _ArrRef):
        return np.load(d / obj.name, mmap_mode="r", allow_pickle=False)
    if isinstance(obj, _DataRef):
        x = Data.__new__(Data)
        x.__dict__.update({k: _restore(v, d) for k, v in obj.state.items()})
        return x
    if isinstance(obj, dict):
        return {k: _restore(v, d) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_restore(v, d) for v in obj]
    if isinstance(obj, tuple):
        return tuple(_restore(v, d) for v in obj)
    return obj


def stage(payload: dict, name: str = "pool", run_dir: Path | None = None, min_bytes: int = MIN_BYTES) -> str:
    """Write ``payload`` under a fresh directory; returns its path (the only thing passed to workers)."""
    from . import analyses as AN
    base = run_dir or (root() / AN.run_id())
    d = base / f"{name}-{uuid.uuid4().hex[:8]}"      # literal-ok: short unique suffix
    d.mkdir(parents=True, exist_ok=False)
    written: dict = {}
    stripped = _strip(payload, d, written, "p", min_bytes)
    with open(d / "payload.pkl", "wb") as fh:
        pickle.dump(stripped, fh, protocol=pickle.HIGHEST_PROTOCOL)
    return str(d)


def load(path: str) -> dict:
    d = Path(path)
    with open(d / "payload.pkl", "rb") as fh:
        return _restore(pickle.load(fh), d)


def files(path: str) -> list:
    return sorted(os.listdir(path))


def cleanup(path: str) -> None:
    p = Path(path)
    shutil.rmtree(p, ignore_errors=True)
    parent = p.parent
    if parent.exists() and parent != root() and not any(parent.iterdir()):
        shutil.rmtree(parent, ignore_errors=True)


@contextlib.contextmanager
def staged(payload: dict, name: str = "pool"):
    path = stage(payload, name)
    try:
        yield path
    finally:
        cleanup(path)
