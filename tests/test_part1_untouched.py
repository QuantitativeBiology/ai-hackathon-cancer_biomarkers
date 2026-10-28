"""Part 1 (../project) is read-only: its git tree is clean, it is at the pinned commit, and no file or
directory outside .git is newer than this project's first commit."""
import os

import pytest

from biomark import gitutil
from biomark import paths as P

pytestmark = pytest.mark.skipif(not P.P1.exists(), reason="Part 1 tree not present")


def test_part1_git_clean():
    assert gitutil.status_porcelain(P.P1) == []


def test_part1_at_pinned_commit():
    assert gitutil.head(P.P1) == P.part1_commit()


def test_no_part1_file_newer_than_project_start():
    t0 = gitutil.first_commit_time(P.PROJECT)
    newer, n = [], 0
    for root, dirs, files in os.walk(P.P1):
        if os.path.samefile(root, P.P1):
            dirs[:] = [d for d in dirs if d != ".git"]
        if os.stat(root).st_mtime > t0:
            newer.append(root + os.sep)
        for f in files:
            n += 1
            p = os.path.join(root, f)
            if os.stat(p).st_mtime > t0:
                newer.append(p)
    assert n > 100  # really scanned the tree (src, scripts, configs, results, docs, cache, data)
    assert newer == [], f"{len(newer)} Part 1 paths modified after project start, e.g. {newer[:5]}"


def test_no_part1_pycache_written_by_import():
    from biomark import p1
    import sys

    p1.import_essml()
    assert sys.dont_write_bytecode
    t0 = gitutil.first_commit_time(P.PROJECT)
    for root, _dirs, files in os.walk(P.P1_SRC):
        for f in files:
            if f.endswith(".pyc"):
                assert os.stat(os.path.join(root, f)).st_mtime <= t0, f


def test_write_blocker_refuses_part1_writes(tmp_path):
    import numpy as np
    import pandas as pd

    from biomark import p1

    target = P.P1 / "cache" / "__biomark_write_probe__"
    with pytest.raises(PermissionError):
        open(str(target) + ".txt", "w")
    with pytest.raises(PermissionError):
        np.save(str(target) + ".npy", np.zeros(1))
    with pytest.raises(PermissionError):
        pd.DataFrame({"a": [1]}).to_parquet(str(target) + ".parquet")
    with pytest.raises(PermissionError):
        os.mkdir(str(target))
    assert not any(P.P1.joinpath("cache").glob("__biomark_write_probe__*"))
    # writes elsewhere are unaffected
    pd.DataFrame({"a": [1]}).to_parquet(tmp_path / "ok.parquet")
    assert p1.Part1WriteError.__mro__[1] is PermissionError


def test_essml_write_hazards_redirected():
    from biomark import p1

    p1.import_essml()
    import essml.guard as g
    import essml.modelstore as ms
    import essml.models.m2_local as m2
    import essml.models.m5_stack as m5

    for mod, attr in [(m2, "INNER_CACHE"), (m5, "INNER_CACHE"), (m5, "OUTER_CACHE"), (ms, "MODEL_DIR"),
                      (g, "ACCESS_LOG")]:
        assert P.is_under(getattr(mod, attr), P.CACHE), (mod.__name__, attr)
    assert P.is_under(ms.save_fold_model.__defaults__[0], P.CACHE)
    assert P.is_under(ms.load_fold_models.__defaults__[0], P.P1)  # reading P1 fold models still works
    assert all(not (hasattr(v, "__fspath__") and P.is_under(v, P.P1)) for v in g._log.__defaults__)
    import essml.io as eio
    with pytest.raises(PermissionError):
        eio.build_cache()
