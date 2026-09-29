"""Shared fixtures of the precomp test suite.

The package is imported from the source tree (python/), so the tests run
without installing it. `fake_solver` is an executable wrapper around
`fake_sparlab_form.py`, a test double that writes contract-conforming result
directories; the real sparlab_form is only used by test_integration_sparlab.py,
which is skipped when build/bin/sparlab_form is absent.
"""

import stat
import sys
from pathlib import Path

import numpy as np
import pytest

PYTHON_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = PYTHON_DIR.parent
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from precomp.geometry import Grid, TruncatedCone  # noqa: E402


@pytest.fixture(scope="session")
def fake_solver(tmp_path_factory) -> Path:
    """Path of an executable that behaves like sparlab_form (see fake_sparlab_form.py)."""
    d = tmp_path_factory.mktemp("fake_bin")
    exe = d / "sparlab_form"
    script = Path(__file__).with_name("fake_sparlab_form.py")
    exe.write_text(f'#!/bin/sh\nexec "{sys.executable}" "{script}" "$@"\n', encoding="utf-8")
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return exe


@pytest.fixture
def counter(tmp_path, monkeypatch) -> Path:
    """A file the fake solver appends one line to per run."""
    path = tmp_path / "runs.log"
    monkeypatch.setenv("FAKE_SPARLAB_COUNTER", str(path))
    monkeypatch.delenv("FAKE_SPARLAB_FAIL", raising=False)
    return path


def run_count(path: Path) -> int:
    return len(path.read_text().splitlines()) if path.exists() else 0


@pytest.fixture
def small_setup(fake_solver):
    """A 0.12 m blank with 2.5 mm elements, run by the fake solver."""
    from precomp.fea import FormingSetup
    from precomp.materials import get_material

    return FormingSetup(get_material("AA5754-O"), blank_size=0.12, clamp_margin=0.01,
                        element_size=2.5e-3, layers=1, step_down=1e-3,
                        toolpath_spacing=2e-3, executable=str(fake_solver))


@pytest.fixture
def small_cone():
    return TruncatedCone(0.03, 45.0, 0.012, 0.004, 0.004)


@pytest.fixture
def node_grid() -> Grid:
    """The node grid of `small_setup`'s mesh (spacing = element size)."""
    return Grid.centered(0.12, 2.5e-3)


@pytest.fixture
def rng():
    return np.random.default_rng(20240917)
