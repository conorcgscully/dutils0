from __future__ import annotations

import sys
from pathlib import Path


def pytest_configure() -> None:
    """
    Ensure the vendored `rdsl` package is importable for dutils0 tests.

    `rdsl` is vendored in this repo under `rdsl/src/rdsl` but is not a declared
    dependency of dutils0, so we add `rdsl/src` to sys.path for test runs.
    """
    repo_root = Path(__file__).resolve().parents[2]
    rdsl_src = repo_root / "rdsl" / "src"
    if rdsl_src.exists():
        sys.path.insert(0, str(rdsl_src))

