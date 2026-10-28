"""Read-only git helpers. ``--no-optional-locks`` stops ``git status`` refreshing the index, so querying
Part 1's repository never writes into ``P1/.git``."""
from __future__ import annotations

import subprocess
from pathlib import Path


def git(repo, *args, check: bool = False) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "--no-optional-locks", *args], cwd=str(Path(repo)), capture_output=True,
                          text=True, check=check)


def head(repo) -> str:
    return git(repo, "rev-parse", "HEAD", check=True).stdout.strip()


def status_porcelain(repo, untracked: bool = True) -> list[str]:
    args = ["status", "--porcelain"] + ([] if untracked else ["--untracked-files=no"])
    r = git(repo, *args, check=True)
    return [ln for ln in r.stdout.splitlines() if ln.strip()]


def tag_exists(repo, tag: str) -> bool:
    return git(repo, "rev-parse", "-q", "--verify", f"refs/tags/{tag}").returncode == 0


def first_commit_time(repo) -> int:
    """Unix committer time of the root commit (``git log --reverse``)."""
    out = git(repo, "log", "--reverse", "--format=%ct", check=True).stdout.split()
    return int(out[0])


def ls_files(repo) -> list[str]:
    return [ln for ln in git(repo, "ls-files", check=True).stdout.splitlines() if ln]
