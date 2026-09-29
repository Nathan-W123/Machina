"""Integration with the real sparlab_form, skipped when it has not been built.

The C++ forming app is developed separately; this test runs the smallest
meaningful deck through it and checks the result against the contract the
loader enforces, not against physics (that is the C++ suite's job).
"""

import os

import numpy as np
import pytest

from conftest import REPO_ROOT
from precomp.fea import EXECUTABLE_ENV, FormingSetup, simulate
from precomp.geometry import Grid, TruncatedCone
from precomp.materials import get_material

EXE = os.environ.get(EXECUTABLE_ENV) or str(REPO_ROOT / "build" / "bin" / "sparlab_form")

pytestmark = pytest.mark.skipif(not (os.path.isfile(EXE) and os.access(EXE, os.X_OK)),
                                reason=f"{EXE} is not built")


def test_small_cone_runs_through_sparlab_form(tmp_path):
    setup = FormingSetup(get_material("AA5754-O"), blank_size=0.08, clamp_margin=0.01,
                         element_size=4e-3, layers=1, step_down=2e-3, toolpath_spacing=4e-3,
                         executable=EXE, timeout=1800.0)
    target = TruncatedCone(0.02, 40.0, 0.004, 0.003, 0.003).heightmap(Grid.centered(0.08, 1e-3))
    res = simulate(setup, target, tmp_path / "w")
    assert res.step_names[0] == "form"
    formed = res.formed_surface(grid=target.grid)
    assert formed.mask.mean() > 0.9
    assert 0.0 < formed.depth < 2 * target.depth
    assert np.all(np.isfinite(res.forming_forces()["f"]))
