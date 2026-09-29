"""Robot compliance: exact inverses of the linear deflection model."""

import numpy as np
import pandas as pd
import pytest

from precomp import PrecompError
from precomp.robot import (CartesianCompliance, compensate_toolpath, forces_on_path)
from precomp.toolpath import Toolpath


def _path(n=50):
    s = np.linspace(0, 2 * np.pi, n)
    return np.column_stack([0.03 * np.cos(s), 0.03 * np.sin(s), -0.01 * s / (2 * np.pi)])


def test_constant_force_precompensation_is_exact():
    C = CartesianCompliance(np.array([[2e-6, 1e-7, 0.0], [1e-7, 3e-6, 0.0], [0.0, 0.0, 1e-6]]))
    p = _path()
    f = np.column_stack([np.zeros(50), 50 * np.ones(50), 800 * np.ones(50)])
    cmd, info = C.precompensate_path(p, f)
    assert np.abs(C.deflected_path(cmd, f) - p).max() < 1e-18
    assert info["iterations"] == 0
    assert np.allclose(C.deflection(f)[0], C.C @ f[0])


def test_fixed_point_with_a_force_that_depends_on_the_command():
    """The indentation, hence the force, grows with how deep the robot is
    commanded: f = f0 + K p_cmd. The exact command solves p + C f(p) = target."""
    C = CartesianCompliance.diagonal(1e-6, 1e-6, 2e-6)
    K = np.diag([0.0, 0.0, -1e5])                 # N/m: ||C K|| = 0.2 < 1
    f0 = np.array([0.0, 0.0, 500.0])

    def force(cmd):
        return f0 + cmd @ K.T

    p = _path()
    cmd, info = C.precompensate_path(p, force_model=force, tolerance=1e-15)
    assert np.abs(cmd + C.deflection(force(cmd)) - p).max() < 1e-15
    assert info["iterations"] >= 3                           # contraction 0.2 per step
    first_order = p - C.deflection(force(p))
    assert np.abs(first_order + C.deflection(force(first_order)) - p).max() > 1e-8
    unstable = CartesianCompliance.diagonal(0.0, 0.0, 2e-5)   # ||C K|| = 2: diverges
    with pytest.raises(PrecompError, match="converge"):
        unstable.precompensate_path(p, force_model=force, iterations=20)


def test_compliance_from_joint_stiffness():
    J = np.array([[0.5, 0.2, 0.0], [0.0, 0.4, 0.1], [0.1, 0.0, 0.3]])
    K = np.array([1e5, 2e5, 5e4])
    C = CartesianCompliance.from_joint_stiffness(J, K)
    assert np.allclose(C.C, J @ np.diag(1 / K) @ J.T)
    assert np.allclose(C.C, C.C.T)
    with pytest.raises(ValueError):
        CartesianCompliance(np.eye(2))


def test_forces_on_path_interpolate_the_history_with_the_stated_sign():
    pts = _path()
    path = Toolpath(pts, np.ones(len(pts), int), 0.005)
    hist = pd.DataFrame({"step": 0, "increment": range(3), "t": [0.0, 0.5, 1.0], "tool": "tool",
                         "fx": [0.0, 10.0, 20.0], "fy": 0.0, "fz": [-100.0, -300.0, -500.0]})
    f = forces_on_path(path, hist, tool="tool")
    assert np.allclose(f[:, 2], 100.0 + 400.0 * path.t)       # on the tool: minus on-sheet
    g = forces_on_path(path, hist, sign="on_tool")
    assert np.allclose(g, -f)
    comp = compensate_toolpath(path, CartesianCompliance.diagonal(1e-6, 1e-6, 1e-6), f)
    assert comp.metadata["robot_compensation"]["max_correction_m"] == pytest.approx(
        1e-6 * np.linalg.norm([20.0, 0.0, 500.0]))
