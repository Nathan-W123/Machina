"""Robot compliance at the forming tool: published joint stiffnesses -> C [m/N].

ABB IRB 7600 (500 kg payload) joint stiffnesses identified by Bharti et al.,
Sci. Rep. 14 (2024) 20183, Table 3 (VJM, rigid links):
https://www.nature.com/articles/s41598-024-70746-3/tables/3
  K = 3.73e6, 4.82e6, 4.04e6, 6.38e5, 7.59e5, 4.72e5 N m/rad.

Kinematics: standard DH of an IRB 7600-500/2.55-class arm (axis-2 height
0.780 m, offset 0.410, lower arm 1.075, axis-3 offset 0.165, upper arm 1.056,
wrist to flange 0.250 m; approximate public dimensions - INFERENCE, not a
calibrated model) plus a forming tool of length L_tool (force/torque sensor,
holder, stylus). The translational compliance at the tool centre point is
C = J_v K^-1 J_v^T (precomp.robot.CartesianCompliance.from_joint_stiffness),
J_v by central differences of the forward kinematics.

Two installations:
  * "vertical_sheet": the Machina-style cell - the sheet stands in a frame
    between two rail-mounted robots, each robot reaches forward horizontally
    (tool axis horizontal); the sheet normal is the robot's +x.
  * "horizontal_sheet": a sheet clamped flat below the robot, tool pointing
    down (the Bharti et al. / most laboratory set-ups).

Writes robot_compliance.json (C in sheet axes: z = sheet normal, pointing
out of the tool side, as in the SparLab deck) for a grid of poses.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "python"))
from precomp.robot import CartesianCompliance  # noqa: E402

K_IRB7600 = np.array([3.73e6, 4.82e6, 4.04e6, 6.38e5, 7.59e5, 4.72e5])   # N m/rad
# a, alpha, d, theta offset (standard DH)
DH = [(0.410, -np.pi / 2, 0.780, 0.0),
      (1.075, 0.0, 0.0, -np.pi / 2),
      (0.165, -np.pi / 2, 0.0, 0.0),
      (0.0, np.pi / 2, 1.056, 0.0),
      (0.0, -np.pi / 2, 0.0, 0.0),
      (0.0, 0.0, 0.250, 0.0)]


def fk(q, tool=0.30):
    T = np.eye(4)
    for i, (a, al, d, th0) in enumerate(DH):
        dd = d + (tool if i == 5 else 0.0)
        th = q[i] + th0
        ct, st, ca, sa = np.cos(th), np.sin(th), np.cos(al), np.sin(al)
        A = np.array([[ct, -st * ca, st * sa, a * ct],
                      [st, ct * ca, -ct * sa, a * st],
                      [0, sa, ca, dd],
                      [0, 0, 0, 1]])
        T = T @ A
    return T


def jac_v(q, tool):
    J = np.zeros((3, 6))
    h = 1e-7
    for i in range(6):
        dq = np.zeros(6)
        dq[i] = h
        J[:, i] = (fk(q + dq, tool)[:3, 3] - fk(q - dq, tool)[:3, 3]) / (2 * h)
    return J


def ik(p, axis, tool, q0):
    axis = np.asarray(axis, float) / np.linalg.norm(axis)

    def res(q):
        T = fk(q, tool)
        return np.concatenate([(T[:3, 3] - p) * 10.0, T[:3, 2] - axis])
    r = least_squares(res, q0, xtol=1e-14, ftol=1e-14, gtol=1e-14)
    return r.x, float(np.abs(res(r.x)).max())


def to_sheet(C, install):
    """C (robot axes) -> sheet axes (x, y in the sheet plane, z its normal
    pointing towards the tool, i.e. against the direction the tool pushes)."""
    if install == "vertical_sheet":
        # tool pushes along robot +x; sheet z = -x_r, sheet x = y_r, sheet y = z_r
        R = np.array([[0, 1, 0], [0, 0, 1], [-1, 0, 0]], float)
    else:
        # tool pushes along robot -z; sheet z = z_r
        R = np.eye(3)
    return R @ C @ R.T


def main():
    out = {"source": "Bharti et al., Sci. Rep. 14 (2024), Table 3 (ABB IRB 7600, VJM)",
           "url": "https://www.nature.com/articles/s41598-024-70746-3",
           "joint_stiffness_Nm_per_rad": K_IRB7600.tolist(), "dh": DH, "poses": []}
    q0 = np.array([0.0, 0.2, -0.2, 0.0, 0.0, 0.0])
    rows = []
    for install in ("vertical_sheet", "horizontal_sheet"):
        for tool in (0.20, 0.30, 0.45):
            for reach in (1.4, 1.8, 2.2):
                for h in ((0.8, 1.3, 1.8) if install == "vertical_sheet" else (0.3, 0.6)):
                    if install == "vertical_sheet":
                        p, axis = np.array([reach, 0.0, h]), [1.0, 0.0, 0.0]
                        qs = np.array([0.0, 0.3, 0.0, 0.0, -0.3, 0.0])
                    else:
                        p, axis = np.array([reach, 0.0, h]), [0.0, 0.0, -1.0]
                        qs = np.array([0.0, 0.5, 0.3, 0.0, 0.8, 0.0])
                    best = None
                    for trial in range(8):
                        start = qs + (0 if trial == 0 else np.random.default_rng(trial).normal(0, 0.4, 6))
                        q, err = ik(p, axis, tool, start)
                        if best is None or err < best[1]:
                            best = (q, err)
                    q, err = best
                    if err > 1e-6:
                        continue
                    J = jac_v(q, tool)
                    C = CartesianCompliance.from_joint_stiffness(J, K_IRB7600).C
                    Cs = to_sheet(C, install)
                    ev = np.linalg.eigvalsh(Cs)
                    rec = {"install": install, "tool_m": tool, "reach_m": reach, "height_m": h,
                           "q_rad": q.tolist(), "ik_err": err, "C_sheet_m_per_N": Cs.tolist(),
                           "cz_um_per_N": Cs[2, 2] * 1e6, "cx_um_per_N": Cs[0, 0] * 1e6,
                           "cy_um_per_N": Cs[1, 1] * 1e6, "cmax_um_per_N": ev.max() * 1e6,
                           "kz_N_per_um": 1e-6 / Cs[2, 2]}
                    out["poses"].append(rec)
                    rows.append(rec)
                    print(f"{install:16s} tool {tool:.2f} reach {reach:.1f} h {h:.1f}: "
                          f"cz {Cs[2,2]*1e6:.3f} cx {Cs[0,0]*1e6:.3f} cy {Cs[1,1]*1e6:.3f} "
                          f"cmax {ev.max()*1e6:.3f} um/N  (kz {1e-6/Cs[2,2]:.2f} N/um)")
    (HERE / "robot_compliance.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
