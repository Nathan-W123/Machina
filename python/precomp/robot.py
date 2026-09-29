"""Robot compliance: how far the tool is pushed off its path, and the pre-correction.

An industrial robot holding a forming tool is compliant: the forming force
deflects the tool-centre point by ``delta = C f`` with C the 3 x 3 Cartesian
compliance at the working pose [m/N] and f the force the sheet exerts on the
tool [N] (the reaction of the forming force). The commanded path is then
pre-deflected against the load:

* `deflected_path(path, forces)` = path + C f(path) - where the tool
  actually goes;
* `precompensate_path(path, forces)` = path - C f, exact when f is the
  force at the target positions (which is where the corrected tool is); with
  a `force_model` giving the force as a function of the *commanded* path, a
  fixed-point iteration solves ``p_cmd + C f(p_cmd) = p_target``,
  converging when ||C df/dp|| < 1.

Forces from a simulation: `forces_on_path` interpolates the tool force
history of a sparlab_form result at the path's pseudo-time. It takes
tool_forces.csv's fx, fy, fz to be the force the tool exerts on the sheet
(fz < 0 while pushing down), so the force on the robot is their negative;
`sign` states the convention explicitly if the solver's differs.

A single constant C is a first-order model: the compliance of a serial robot
changes with its pose. `CartesianCompliance.from_joint_stiffness` computes C
for one pose; a pose-dependent map needs one C per path point (pass an
(n, 3, 3) array to `deflection`).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from ._util import PrecompError
from .toolpath import Toolpath


@dataclass
class CartesianCompliance:
    """Translational compliance C [m/N] of the tool-centre point (3 x 3, or
    (n, 3, 3) for one matrix per path point)."""

    C: np.ndarray

    def __post_init__(self) -> None:
        C = np.asarray(self.C, dtype=float)
        if C.shape[-2:] != (3, 3) or C.ndim not in (2, 3):
            raise ValueError(f"compliance must be 3 x 3 or (n, 3, 3), got {C.shape}")
        if not np.all(np.isfinite(C)):
            raise ValueError("compliance must be finite")
        self.C = C

    @classmethod
    def diagonal(cls, cx: float, cy: float, cz: float) -> "CartesianCompliance":
        """Uncoupled compliance diag(cx, cy, cz) [m/N]."""
        return cls(np.diag([cx, cy, cz]))

    @classmethod
    def from_joint_stiffness(cls, jacobian: np.ndarray,
                             joint_stiffness: np.ndarray) -> "CartesianCompliance":
        """C = J_v K^-1 J_v^T from the translational rows J_v (3 x n) [m/rad]
        of the robot Jacobian at the pose and the joint stiffnesses K (n,)
        [N m/rad] - the classical joint-compliance model (links rigid)."""
        J = np.asarray(jacobian, dtype=float)
        if J.shape[0] == 6:
            J = J[:3]
        K = np.asarray(joint_stiffness, dtype=float)
        if J.ndim != 2 or J.shape[0] != 3 or K.shape != (J.shape[1],):
            raise ValueError("jacobian must be 3 x n (or 6 x n) and joint_stiffness (n,)")
        if np.any(K <= 0):
            raise ValueError("joint stiffnesses must be > 0")
        return cls(J @ np.diag(1.0 / K) @ J.T)

    def deflection(self, forces: np.ndarray) -> np.ndarray:
        """(n, 3) tool-centre deflections C f [m] under forces (n, 3) [N]
        acting on the tool."""
        f = np.asarray(forces, dtype=float)
        if self.C.ndim == 2:
            return f @ self.C.T
        if self.C.shape[0] != len(f):
            raise ValueError("one compliance matrix per force is required")
        return np.einsum("nij,nj->ni", self.C, f)

    def deflected_path(self, path: np.ndarray, forces: np.ndarray) -> np.ndarray:
        """Where the tool goes when commanded along `path` (n, 3) [m] under
        `forces` (n, 3) [N] on the tool: path + C f."""
        return np.asarray(path, float) + self.deflection(forces)

    def precompensate_path(self, path: np.ndarray, forces: Optional[np.ndarray] = None, *,
                           force_model: Optional[Callable[[np.ndarray], np.ndarray]] = None,
                           iterations: int = 50, tolerance: float = 1e-12
                           ) -> Tuple[np.ndarray, Dict[str, object]]:
        """The commanded path that makes the tool follow `path` (n, 3) [m].

        With `forces` (n, 3) [N on the tool] known at the target positions:
        path - C f, exact. This covers forces that depend on where the tool
        actually is: when the correction works, the tool is at the target, so
        the force there is the one to use.

        With `force_model(commanded_path) -> forces` - the forces the process
        produces when the robot is *commanded* along a path (a re-simulation
        along it, or a model of how the indentation grows with the command) -
        the fixed point ``p_cmd = path - C f(p_cmd)`` is iterated from the
        first-order guess until the miss |p_cmd + C f(p_cmd) - path| is below
        `tolerance` [m]; it converges when ||C df/dp|| < 1. Returns
        (commanded path, info) with the iterations and the miss history [m];
        raises PrecompError when it has not converged.
        """
        target = np.asarray(path, dtype=float)
        if force_model is None:
            if forces is None:
                raise ValueError("give forces or a force_model")
            return target - self.deflection(forces), {"iterations": 0, "miss_m": [0.0]}
        cmd = target.copy()
        history: List[float] = []
        for it in range(1, iterations + 1):
            cmd = target - self.deflection(force_model(cmd))
            miss = float(np.abs(cmd + self.deflection(force_model(cmd)) - target).max())
            history.append(miss)
            if not np.isfinite(miss):
                break
            if miss < tolerance:
                return cmd, {"iterations": it, "miss_m": history}
        raise PrecompError(f"precompensate_path did not converge: miss {history[-1]:.3g} m "
                           f"after {len(history)} iterations (is ||C df/dp|| < 1?)")


def forces_on_path(path: Toolpath, tool_forces: pd.DataFrame, *, tool: Optional[str] = None,
                   step: Optional[int] = None, sign: str = "on_sheet") -> np.ndarray:
    """Forces on the tool (n, 3) [N] at the points of `path`, from a force history.

    `tool_forces` is a FormingResult's tool force table (columns t, tool,
    step, fx, fy, fz); rows are filtered by `tool` and `step` (default: the
    first step, the forming one) and interpolated linearly in the pseudo-time
    t at the path's `t`. `sign`: "on_sheet" when fx, fy, fz are the force the
    tool exerts on the sheet (the result is their negative, the load on the
    robot), "on_tool" when they already are the force on the tool.
    """
    frame = tool_forces
    if tool is not None:
        frame = frame[frame["tool"] == tool]
    if step is None:
        step = int(frame["step"].min())
    frame = frame[frame["step"] == step].sort_values("t")
    if len(frame) < 2:
        raise PrecompError("the force history has fewer than two rows for this tool and step")
    t = frame["t"].to_numpy(float)
    if np.any(np.diff(t) < 0):
        raise PrecompError("the force history's t is not monotone")
    f = frame[["fx", "fy", "fz"]].to_numpy(float)
    if sign == "on_sheet":
        f = -f
    elif sign != "on_tool":
        raise ValueError("sign must be 'on_sheet' or 'on_tool'")
    tp = path.t
    return np.column_stack([np.interp(tp, t, f[:, i]) for i in range(3)])


def compensate_toolpath(path: Toolpath, compliance: CartesianCompliance,
                        forces: np.ndarray) -> Toolpath:
    """A copy of `path` pre-deflected against `forces` (n, 3) [N on the
    tool] to first order (path - C f); metadata records the largest
    correction [m]."""
    cmd, _ = compliance.precompensate_path(path.points, forces)
    meta = dict(path.metadata)
    meta["robot_compensation"] = {
        "max_correction_m": float(np.linalg.norm(cmd - path.points, axis=1).max())}
    return Toolpath(cmd, path.level.copy(), path.tool_radius, meta)
