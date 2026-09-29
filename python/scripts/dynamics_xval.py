"""Cross-validation of SparLab's transient and harmonic analyses, used by
cross_validate.py for a run with a `transient` or `frequency_response` block.

The same discrete problem is solved again, independently:

  * by scikit-fem: its own stiffness and mass matrices - consistent, or the
    lumped mass SparLab uses (row sums for the linear elements, HRZ diagonal
    scaling for the Tet10) - and an HHT-alpha integration written in the
    acceleration (predictor-corrector) form, where SparLab solves in the
    displacement form, with the Rayleigh damping and the amplitude that
    summary.json records and the prescribed DOFs moving with the Newmark
    kinematics of their displacement history, as SparLab's do; a non-linear
    run with the independent J2 / Saint Venant-Kirchhoff system of
    SkfemProblem.j2_system, converged by Newton's method at every step with
    its plastic history committed there. A harmonic response is a direct
    complex solve of [K (1 + i eta) - omega^2 M + i omega (a M + b K)] U = f at
    every frequency.
  * by CalculiX, from the run's `calculix_<lc>_dynamic.inp` (*DYNAMIC, DIRECT,
    ALPHA; see CalculixWriter.hpp), for the solid elements. CalculiX's
    expanded plane elements are not used: in *DYNAMIC their step response
    contradicts CalculiX's own *FREQUENCY result for the same mesh (a CPS4
    cantilever whose first period is 12.1 ms peaks after 4.1 ms instead of
    6.1 ms; CPE4 alike - measured). Its *STEADY STATE DYNAMICS is a modal
    superposition, not a direct solve, so the harmonic response is compared
    with scikit-fem only.

Compared are the monitors at every step (or frequency), the displacement field
at every snapshot and the final displacement, velocity and acceleration (the
complex fields at the snapshot frequencies), each as the largest difference
over the largest value of the reference.
"""

from __future__ import annotations

import json
import os
import tempfile
from typing import Dict, List, Optional

import numpy as np

from sparlab_viz.loaders import ResultError, read_legacy_vtk


def _safe(name: str) -> str:
    return "".join(c if (c.isalnum() or c in "-_") else "_" for c in name) or "unnamed"


# ---------------------------------------------------------------------------
# The amplitude of a transient, as SparLab evaluates it (Dynamics.cpp)
# ---------------------------------------------------------------------------
class Amplitude:
    """A(t) from the `amplitude` block of summary.json: a step, a table
    (piecewise linear, constant beyond its ends, the slope to the right of a
    point) or a harmonic scale sin(2 pi f t + phase)."""

    def __init__(self, block: Dict):
        self.kind = block["type"]
        self.scale = float(block.get("scale", 1.0))
        self.times = np.asarray(block.get("times_s", []), dtype=float)
        self.values = np.asarray(block.get("values", []), dtype=float)
        self.frequency = float(block.get("frequency_Hz", 0.0))
        self.phase = float(block.get("phase_rad", 0.0))

    def value(self, t: float) -> float:
        if self.kind == "step":
            return self.scale
        if self.kind == "harmonic":
            return self.scale * np.sin(2.0 * np.pi * self.frequency * t + self.phase)
        if t <= self.times[0]:
            return self.scale * self.values[0]
        if t >= self.times[-1]:
            return self.scale * self.values[-1]
        i = int(np.searchsorted(self.times, t, side="right"))  # std::upper_bound
        s = (t - self.times[i - 1]) / (self.times[i] - self.times[i - 1])
        return self.scale * (self.values[i - 1] + s * (self.values[i] - self.values[i - 1]))

    def rate(self, t: float) -> float:
        if self.kind == "step":
            return 0.0
        if self.kind == "harmonic":
            w = 2.0 * np.pi * self.frequency
            return self.scale * w * np.cos(w * t + self.phase)
        if t < self.times[0] or t >= self.times[-1]:
            return 0.0
        i = int(np.searchsorted(self.times, t, side="right"))
        return self.scale * (self.values[i] - self.values[i - 1]) / (self.times[i] - self.times[i - 1])

    def second_rate(self, t: float) -> float:
        if self.kind != "harmonic":
            return 0.0
        w = 2.0 * np.pi * self.frequency
        return -self.scale * w * w * np.sin(w * t + self.phase)


# ---------------------------------------------------------------------------
# Mass matrix in scikit-fem
# ---------------------------------------------------------------------------
def assemble_mass(problem, materials: List[Dict], element_materials: Optional[np.ndarray],
                  lumped: bool):
    """The mass matrix in scikit-fem's numbering: int rho u . v, with a rule
    exact for the integrand (the 3 x 3 (x 3) Gauss rule SparLab uses for Q4 /
    Hex8, degree 2 for the linear simplices, degree 4 for the straight Tet10);
    lumped, the row sums of the linear elements or the HRZ diagonal scaling of
    the Tet10 (each element's diagonal scaled to its mass), as SparLab lumps."""
    import scipy.sparse
    from skfem import Basis, BilinearForm, asm
    from skfem.helpers import dot

    element = problem.mesh.element_type
    intorder = {"Quad4": 4, "Hex8": 4, "Tri3": 2, "Tet4": 2, "Tet10": 4}[element]
    ne = problem.mesh.num_elements
    index = (np.zeros(ne, dtype=int) if element_materials is None
             else np.asarray(element_materials, dtype=int))
    rho_e = np.asarray([float(m["density_kg_per_m3"]) for m in materials])[index]
    basis = Basis(problem.skfem_mesh, problem.vector_element, intorder=intorder)
    if not np.array_equal(basis.nodal_dofs, problem.basis.nodal_dofs):
        raise ResultError("the mass basis numbers its DOFs differently")
    rho = np.repeat(rho_e[:, None], basis.X.shape[-1], axis=1)

    @BilinearForm
    def mass(u, v, w):
        return w["rho"] * dot(u, v)

    consistent = asm(mass, basis, rho=rho) * problem.thickness
    if not lumped:
        return consistent.tocsr()
    if element != "Tet10":
        # Every component couples only to itself: the row sum over the whole
        # row is the scalar element's.
        return scipy.sparse.diags(np.asarray(consistent.sum(axis=1)).ravel()).tocsr()
    scalar = Basis(problem.skfem_mesh, problem.scalar_element, intorder=intorder)
    phi = np.stack([scalar.basis[i][0].value for i in range(scalar.Nbfun)])  # (nbf, ne, nq)
    dx = scalar.dx * rho_e[:, None]
    diagonal = np.einsum("ieq,ieq,eq->ie", phi, phi, dx)
    element_mass = dx.sum(axis=1)
    hrz = diagonal * (element_mass / diagonal.sum(axis=0))[None, :]
    per_scalar = np.zeros(scalar.N)
    np.add.at(per_scalar, scalar.element_dofs, hrz)
    values = np.zeros(problem.basis.N)
    for node in range(problem.mesh.num_nodes):
        for c in range(problem.dim):
            values[problem.dof_of_node[node, c]] = per_scalar[problem.scalar_dof_of_node[node]]
    return scipy.sparse.diags(values).tocsr()


# ---------------------------------------------------------------------------
# HHT-alpha in the acceleration form
# ---------------------------------------------------------------------------
class History:
    """u, v, a (steps + 1, N) and the reactions at the prescribed DOFs."""

    def __init__(self, n: int, steps: int):
        self.u = np.zeros((steps + 1, n))
        self.v = np.zeros((steps + 1, n))
        self.a = np.zeros((steps + 1, n))
        self.r = np.zeros((steps + 1, n))


def _parameters(alpha: float):
    beta = 0.25 * (1.0 - alpha) ** 2
    gamma = 0.5 - alpha
    return beta, gamma, 1.0 + alpha


def linear_transient(problem, m, c, amp: Amplitude, dt: float, steps: int, alpha: float,
                     static_start: bool) -> History:
    """M a + C v + K u = A(t) f, u_p = A(t) g, by HHT-alpha in the acceleration
    form: the Newmark predictors from the last step, the weighted equation of
    motion solved for the new acceleration of the free DOFs (the prescribed
    DOFs' acceleration follows from their displacement), then the correctors."""
    import scipy.sparse.linalg as sla

    k = problem.k.tocsr()
    f = problem.f
    g = problem.x
    free, fixed = problem.free, problem.fixed
    beta, gamma, w = _parameters(alpha)
    n = problem.basis.N
    h = History(n, steps)
    u = np.zeros(n)
    v = np.zeros(n)
    a = np.zeros(n)
    a0 = amp.value(0.0)
    u[fixed] = a0 * g[fixed]
    v[fixed] = amp.rate(0.0) * g[fixed]
    a[fixed] = amp.second_rate(0.0) * g[fixed]
    kff = k[free][:, free].tocsc()
    if static_start:
        u[free] = sla.spsolve(kff, a0 * f[free] - k[free][:, fixed] @ u[fixed])
    mff = m[free][:, free].tocsc()
    a[free] = sla.spsolve(mff, (a0 * f - m @ a - c @ v - k @ u)[free])
    s = (m + (w * gamma * dt) * c + (w * beta * dt * dt) * k).tocsr()
    solve = sla.factorized(s[free][:, free].tocsc())
    s_fp = s[free][:, fixed]

    def reactions(uu, vv, aa, amplitude):
        out = np.zeros(n)
        out[fixed] = (m @ aa + c @ vv + k @ uu - amplitude * f)[fixed]
        return out

    h.u[0], h.v[0], h.a[0], h.r[0] = u, v, a, reactions(u, v, a, a0)
    for step in range(1, steps + 1):
        t = step * dt
        a1_amp = amp.value(t)
        up = u + dt * v + dt * dt * (0.5 - beta) * a
        vp = v + dt * (1.0 - gamma) * a
        a1 = np.zeros(n)
        a1[fixed] = (a1_amp * g[fixed] - up[fixed]) / (beta * dt * dt)
        rhs = ((w * a1_amp - alpha * a0) * f - w * (c @ vp + k @ up) + alpha * (c @ v + k @ u))
        a1[free] = solve(rhs[free] - s_fp @ a1[fixed])
        u = up + beta * dt * dt * a1
        v = vp + gamma * dt * a1
        a = a1
        a0 = a1_amp
        h.u[step], h.v[step], h.a[step], h.r[step] = u, v, a, reactions(u, v, a, a1_amp)
    return h


def nonlinear_transient(problem, system, m, c, amp: Amplitude, dt: float, steps: int,
                        alpha: float, tolerance: float = 1e-11,
                        max_iterations: int = 40) -> History:
    """The non-linear transient from rest: Newton's method at every step on
    M a1 + (1 + alpha)(C v1 + f_int(u1) - A1 f) - alpha (C v0 + f_int(u0) - A0 f)
    = 0 with a1, v1 the Newmark functions of u1 and the tangent
    M / (beta dt^2) + (1 + alpha)(gamma / (beta dt) C + K_T); the J2 state is
    committed on convergence. Converged to `tolerance` of the largest force
    involved, or at the round-off floor 64 eps || |K_eff| |u| || once Newton
    stops reducing the residual."""
    import scipy.sparse.linalg as sla

    f = problem.f
    g = problem.x
    free, fixed = problem.free, problem.fixed
    beta, gamma, w = _parameters(alpha)
    n = problem.basis.N
    eps = np.finfo(float).eps
    h = History(n, steps)
    state = system.state
    zero = np.zeros(system.dtemp.shape)
    u = np.zeros(n)
    v = np.zeros(n)
    a = np.zeros(n)
    a0 = amp.value(0.0)
    u[fixed] = a0 * g[fixed]
    v[fixed] = amp.rate(0.0) * g[fixed]
    a[fixed] = amp.second_rate(0.0) * g[fixed]
    f_int, _, _ = system.evaluate(u, state, zero, False)
    mff = m[free][:, free].tocsc()
    a[free] = sla.spsolve(mff, (a0 * f - m @ a - c @ v - f_int)[free])

    def reactions(internal, vv, aa, amplitude):
        out = np.zeros(n)
        out[fixed] = (m @ aa + c @ vv + internal - amplitude * f)[fixed]
        return out

    h.u[0], h.v[0], h.a[0], h.r[0] = u, v, a, reactions(f_int, v, a, a0)
    for step in range(1, steps + 1):
        t = step * dt
        a1_amp = amp.value(t)
        up = u + dt * v + dt * dt * (0.5 - beta) * a
        vp = v + dt * (1.0 - gamma) * a
        previous = alpha * (c @ v + f_int - a0 * f)
        u1 = up.copy()                      # constant-acceleration predictor
        u1[fixed] = a1_amp * g[fixed]
        last = np.inf
        for _ in range(max_iterations):
            f1, trial, kt = system.evaluate(u1, state, zero, True)
            a1 = (u1 - up) / (beta * dt * dt)
            v1 = vp + gamma * dt * a1
            inertia = m @ a1
            damping = c @ v1
            r = inertia + w * (damping + f1 - a1_amp * f) - previous
            rf = r[free]
            jac = (m / (beta * dt * dt) + (w * gamma / (beta * dt)) * c + w * kt).tocsr()
            scale = max(np.linalg.norm(inertia[free]), np.linalg.norm(a1_amp * f),
                        np.linalg.norm(damping[free]), np.linalg.norm(r[fixed]), 1e-300)
            floor = 64.0 * eps * np.linalg.norm((abs(jac) @ np.abs(u1))[free])
            norm = np.linalg.norm(rf)
            if norm <= tolerance * scale or (norm <= floor and norm > 0.5 * last):
                break
            last = norm
            du = sla.spsolve(jac[free][:, free].tocsc(), -rf)
            u1[free] += du
        else:
            raise ResultError(f"scikit-fem's transient Newton did not converge at t = {t}")
        state = trial
        u, v, a, f_int, a0 = u1, v1, a1, f1, a1_amp
        h.u[step], h.v[step], h.a[step], h.r[step] = u, v, a, reactions(f_int, v, a, a1_amp)
    return h


def harmonic_response(problem, m, frequencies: List[float], eta: float, a: float,
                      b: float) -> np.ndarray:
    """U (frequencies, N): [K (1 + i eta) - omega^2 M + i omega (a M + b K)] U
    = f on the free DOFs, the prescribed DOFs at g, by a complex sparse LU."""
    import scipy.sparse.linalg as sla

    k = problem.k.tocsr()
    free, fixed = problem.free, problem.fixed
    out = np.zeros((len(frequencies), problem.basis.N), dtype=complex)
    for j, freq in enumerate(frequencies):
        omega = 2.0 * np.pi * freq
        kc = 1.0 + 1j * (eta + omega * b)
        mc = -omega * omega + 1j * omega * a
        dyn = (kc * k + mc * m).tocsr()
        u = np.zeros(problem.basis.N, dtype=complex)
        u[fixed] = problem.x[fixed]
        rhs = problem.f[free] - dyn[free][:, fixed] @ u[fixed]
        u[free] = sla.spsolve(dyn[free][:, free].tocsc(), rhs)
        out[j] = u
    return out


def monitor_values(problem, monitor: Dict, u, v, a, r) -> np.ndarray:
    """A monitor over the rows of u, v, a, r (steps x N or frequencies x N):
    the mean of a kinematic component over its nodes, the sum of a reaction."""
    dofs = problem.dof_of_node[np.asarray(monitor["nodes"], dtype=int), int(monitor["component"])]
    source = {"displacement": u, "velocity": v, "acceleration": a, "reaction": r}[
        monitor["quantity"]]
    values = source[:, dofs].sum(axis=1)
    return values if monitor["quantity"] == "reaction" else values / len(dofs)


def _relative(reference: np.ndarray, ours: np.ndarray) -> float:
    return float(np.abs(reference - ours).max() / max(np.abs(reference).max(), 1e-300))


# ---------------------------------------------------------------------------
# CalculiX
# ---------------------------------------------------------------------------
def parse_frd_series(path: str, num_nodes: int) -> Dict[int, np.ndarray]:
    """Every displacement block of a CalculiX .frd file by its increment."""
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()
    out: Dict[int, np.ndarray] = {}
    increment = None
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.lstrip().startswith("1PSTEP"):
            increment = int(line.split()[2])
        if line.startswith(" -4") and "DISP" in line:
            values = np.full((num_nodes, 3), np.nan)
            i += 1
            while i < len(lines) and lines[i].startswith(" -5"):
                i += 1
            while i < len(lines) and lines[i].startswith(" -1"):
                row = lines[i]
                values[int(row[3:13]) - 1] = [float(row[13 + 12 * k: 25 + 12 * k])
                                              for k in range(3)]
                i += 1
            if increment is None:
                raise ResultError(f"{path}: a displacement block without its increment")
            out[increment] = values
            continue
        i += 1
    if not out:
        raise ResultError(f"{path} holds no displacement block")
    return out


# ---------------------------------------------------------------------------
# The comparisons
# ---------------------------------------------------------------------------
def _snapshots(case, name: str, dt: float):
    """SparLab's snapshot fields: {step: displacement (nodes, 3)}."""
    series_path = case.path(f"transient_{_safe(name)}.vtk.series")
    if not os.path.isfile(series_path):
        raise ResultError(f"{series_path} is missing; rerun sparlab_solve with VTK output")
    with open(series_path, "r", encoding="utf-8") as handle:
        series = json.load(handle)
    out = {}
    for entry in series["files"]:
        grid = read_legacy_vtk(case.path(entry["name"]))
        out[int(round(float(entry["time"]) / dt))] = grid.point_data["displacement"]
    return out


def transient_comparisons(case, summary: Dict, name: str, tr_case: Dict, problem,
                          ccx_type: str, skfem_element: str, tolerances: Dict[str, float],
                          skip_calculix: bool, run_calculix) -> Dict:
    """scikit-fem's independent integration of the run, and CalculiX's
    *DYNAMIC deck of it, against SparLab's monitors, snapshots and final
    state."""
    entry: Dict = {"codes": {}, "notes": []}
    mesh = case.mesh
    block = summary["transient"]
    dt = float(block["time_step_s"])
    steps = int(tr_case["requested_steps"])
    if not tr_case.get("completed", True) or int(tr_case["steps"]) != steps:
        raise ResultError(f"SparLab's transient of '{name}' did not complete; the comparison "
                          "needs the whole run")
    alpha = float(block["alpha"])
    amp = Amplitude(block["amplitude"])
    lumped = block["mass"] == "lumped"
    materials = summary.get("materials") or [summary["material"]]
    m = assemble_mass(problem, materials, mesh.element_materials, lumped)
    a_ray = float(block["rayleigh_damping"]["mass_1_per_s"])
    b_ray = float(block["rayleigh_damping"]["stiffness_s"])
    c = (a_ray * m + b_ray * problem.k).tocsr()
    nonlinear = block.get("nonlinear")
    if nonlinear is None:
        h = linear_transient(problem, m, c, amp, dt, steps, alpha, block["start"] == "static")
        what = ("linear: an HHT-alpha integration in the acceleration form of scikit-fem's "
                "own K and M")
    else:
        plastic = tr_case.get("plasticity")
        mean = bool(plastic and plastic.get("mean_dilatation_applied"))
        kinematics = nonlinear["kinematics"]
        if kinematics == "finite" and nonlinear["material_model"] != "saint_venant_kirchhoff":
            entry["notes"].append("the scikit-fem transient has the Saint Venant-Kirchhoff "
                                  "law only; the unit tests cover the neo-Hookean one")
            return entry
        system = problem.j2_system(materials, mesh.element_materials, mean,
                                   kinematics="finite" if kinematics == "finite"
                                   else "small_strain")
        h = nonlinear_transient(problem, system, m, c, amp, dt, steps, alpha)
        what = ("non-linear: Newton's method at every step on the HHT-alpha residual of the "
                "independent " + ("J2" if plastic else "elastic") + " system of "
                "SkfemProblem.j2_system (" + kinematics + (", mean dilatation" if mean else "")
                + "), its history committed on convergence")

    # The monitors at every step.
    history = case.table(f"transient_{_safe(name)}.csv")
    monitor_diff = 0.0
    for monitor in tr_case.get("monitors", []):
        ref = monitor_values(problem, monitor, h.u, h.v, h.a, h.r)
        ours = history[f"{monitor['name']}[{monitor['unit']}]"].to_numpy()
        monitor_diff = max(monitor_diff, _relative(ref, ours))
    # The snapshots and the final state.
    snapshots = _snapshots(case, name, dt)
    ref_fields = np.stack([problem.nodal(h.u[s]) for s in sorted(snapshots)])
    our_fields = np.stack([snapshots[s][:, : mesh.dim] for s in sorted(snapshots)])
    snapshot_diff = _relative(ref_fields, our_fields)
    state = case.table(f"transient_state_{_safe(name)}.csv")
    axes = "xyz"[: mesh.dim]
    final = {}
    for key, field, unit in (("u", h.u, "m"), ("v", h.v, "m/s"), ("a", h.a, "m/s^2")):
        ours = state[[f"{key}{x}[{unit}]" for x in axes]].to_numpy()
        final[key] = _relative(problem.nodal(field[-1]), ours)
    # Judged on the full-precision data (the monitors and the final state,
    # 17 significant digits); the snapshots carry the VTK files' nine, a
    # rounding of up to 5e-9 of a value, and must agree to it.
    judged = max(monitor_diff, final["u"], final["v"], final["a"])
    stats = {
        "max_rel_diff": judged,
        "monitors_max_rel_diff": monitor_diff,
        "snapshots_max_rel_diff": snapshot_diff,
        "snapshots_rounding_floor_rel": 5.0e-9,
        "final_displacement_rel_diff": final["u"],
        "final_velocity_rel_diff": final["v"],
        "final_acceleration_rel_diff": final["a"],
        "snapshots": len(snapshots),
        "steps": steps,
        "ref_max_abs_m": float(np.abs(ref_fields).max()),
        "sparlab_max_abs_m": float(np.abs(our_fields).max()),
        "tolerance": tolerances["skfem_transient"],
        "element": skfem_element + (" lumped" if lumped else ""),
        "comparison": what + ("; judged on the monitors at every step and the final u, v, a "
                              "(full precision), the snapshot fields read from the VTK series "
                              "(nine significant digits) checked beside them"),
    }
    stats["passed"] = (judged <= tolerances["skfem_transient"]
                       and snapshot_diff <= max(tolerances["skfem_transient"], 1.0e-8))
    entry["codes"]["scikit-fem transient"] = stats

    if skip_calculix:
        return entry
    deck = case.path(f"calculix_{_safe(name)}_dynamic.inp")
    if ccx_type not in ("C3D8", "C3D4", "C3D10"):
        entry["notes"].append(f"no CalculiX comparison: CalculiX's {ccx_type} in *DYNAMIC "
                              "contradicts its own *FREQUENCY result (see dynamics_xval)")
        return entry
    if not os.path.isfile(deck):
        entry["notes"].append("no CalculiX *DYNAMIC deck: " + (
            "lumped mass" if lumped else "the run cannot be exported (sparlab_solve's warning "
            "says why)"))
        return entry
    with tempfile.TemporaryDirectory(prefix="sparlab_ccx_dynamic_") as work:
        frd = run_calculix(deck, work)
        series = parse_frd_series(frd, mesh.num_nodes)
    common = [s for s in sorted(snapshots) if s in series]
    if len(common) < len([s for s in snapshots if s > 0]):
        raise ResultError(f"CalculiX stored {sorted(series)} but SparLab's snapshots are at "
                          f"{sorted(snapshots)}")
    ref = np.stack([series[s][:, : mesh.dim] for s in common])
    ours = np.stack([snapshots[s][:, : mesh.dim] for s in common])
    if np.isnan(ref).any():
        raise ResultError(f"CalculiX returned no displacement for some nodes of {name}")
    stats = {
        "max_rel_diff": _relative(ref, ours),
        "max_abs_diff_m": float(np.abs(ref - ours).max()),
        "snapshots": len(common),
        "ref_max_abs_m": float(np.abs(ref).max()),
        "sparlab_max_abs_m": float(np.abs(ours).max()),
        "tolerance": tolerances["calculix_transient"],
        "element": f"{ccx_type} *DYNAMIC" + (" NLGEOM" if nonlinear is not None and
                                             nonlinear["kinematics"] == "finite" else "")
                   + (" *PLASTIC" if nonlinear is not None and tr_case.get("plasticity")
                      else ""),
        "frd_rounding_floor_rel": 5.0e-6,
        "comparison": ("same discrete problem: CalculiX's *DYNAMIC, DIRECT, ALPHA="
                       f"{alpha:g} with the amplitude tabulated at every step, the snapshot "
                       "fields from its .frd (six significant digits)"),
    }
    plastic = tr_case.get("plasticity")
    if plastic and plastic.get("mean_dilatation_applied"):
        stats["passed"] = None
        stats["comparison"] += (f"; SparLab averages the dilatation (B-bar), CalculiX's "
                                f"{ccx_type} does not: informational")
    elif nonlinear is not None and nonlinear["kinematics"] == "finite" and plastic:
        stats["passed"] = None
        stats["comparison"] += ("; different plasticity models at finite strain: "
                                "informational")
    else:
        stats["passed"] = stats["max_rel_diff"] <= tolerances["calculix_transient"]
    entry["codes"]["calculix *DYNAMIC"] = stats
    return entry


def frequency_comparisons(case, summary: Dict, name: str, fr_case: Dict, problem,
                          skfem_element: str, tolerances: Dict[str, float]) -> Dict:
    """scikit-fem's direct complex solve at every frequency against SparLab's
    complex monitors and snapshot fields."""
    entry: Dict = {"codes": {}, "notes": []}
    mesh = case.mesh
    block = summary["frequency_response"]
    frequencies = [float(f) for f in block["frequencies_Hz"]]
    lumped = block["mass"] == "lumped"
    materials = summary.get("materials") or [summary["material"]]
    m = assemble_mass(problem, materials, mesh.element_materials, lumped)
    damping = block["damping"]
    eta = float(damping["structural_loss_factor"])
    a_ray = float(damping["mass_1_per_s"])
    b_ray = float(damping["stiffness_s"])
    u = harmonic_response(problem, m, frequencies, eta, a_ray, b_ray)
    omega = 2.0 * np.pi * np.asarray(frequencies)[:, None]
    k = problem.k.tocsr()
    kc = 1.0 + 1j * (eta + omega * b_ray)
    mc = -omega * omega + 1j * omega * a_ray
    reactions = np.zeros_like(u)
    fixed = problem.fixed
    for j in range(len(frequencies)):
        full = kc[j, 0] * (k @ u[j]) + mc[j, 0] * (m @ u[j]) - problem.f
        reactions[j, fixed] = full[fixed]
    table = case.table(f"frequency_response_{_safe(name)}.csv")
    monitor_diff = 0.0
    for monitor in fr_case.get("monitors", []):
        ref = monitor_values(problem, monitor, u, 1j * omega * u, -omega * omega * u, reactions)
        label = monitor["name"]
        unit = monitor["unit"]
        ours = (table[f"{label}_re[{unit}]"].to_numpy()
                + 1j * table[f"{label}_im[{unit}]"].to_numpy())
        monitor_diff = max(monitor_diff, _relative(ref, ours))
    field_diff = 0.0
    snapshots = fr_case.get("snapshot_frequencies_Hz", [])
    width = max(4, len(str(max(len(snapshots) - 1, 0))))
    for i, freq in enumerate(snapshots):
        path = case.path(f"frequency_response_field_{_safe(name)}_{i:0{width}d}.csv")
        field = case.table(os.path.basename(path))
        axes = "xyz"[: mesh.dim]
        ours = np.stack([field[f"u{x}_re[m]"].to_numpy() + 1j * field[f"u{x}_im[m]"].to_numpy()
                         for x in axes], axis=1)
        j = int(np.argmin(np.abs(np.asarray(frequencies) - float(freq))))
        ref = np.stack([u[j][problem.dof_of_node[:, c]] for c in range(mesh.dim)], axis=1)
        field_diff = max(field_diff, _relative(ref, ours))
    judged = max(monitor_diff, field_diff)
    stats = {
        "max_rel_diff": judged,
        "monitors_max_rel_diff": monitor_diff,
        "fields_max_rel_diff": field_diff,
        "frequencies": len(frequencies),
        "snapshots": len(snapshots),
        "ref_max_abs_m": float(np.abs(u).max()),
        "tolerance": tolerances["skfem_harmonic"],
        "element": skfem_element + (" lumped" if lumped else ""),
        "comparison": ("a direct complex solve per frequency of scikit-fem's own K and M, "
                       "structural and Rayleigh damping; the monitors at every frequency, the "
                       "complex fields at the snapshot frequencies"),
    }
    stats["passed"] = judged <= tolerances["skfem_harmonic"]
    entry["codes"]["scikit-fem harmonic"] = stats
    return entry
