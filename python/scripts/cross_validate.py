#!/usr/bin/env python3
"""Cross-validate SparLab's static solution against two independent FE codes.

usage:
    python3 python/scripts/cross_validate.py \
        --case results/cantilever_analysis --case results/block_3d_analysis \
        --output results/cross_validation

For every result directory given (a `sparlab_solve` run made with
`--export-calculix`), the same discrete problem - identical mesh, element type,
material, consistent nodal loads and prescribed DOFs, all read from the run's
own `mesh.json` and CalculiX decks - is solved by

  * CalculiX (`ccx`), from the exported `calculix_<load case>.inp`, and
  * scikit-fem, assembled in this script from `mesh.json` with the linear
    isotropic elasticity form on ElementQuad1 / ElementHex1 / ElementTriP1 /
    ElementTetP1, and on ElementTetP2 over an isoparametric MeshTet2 for the
    quadratic tetrahedra (curved cells included),

and the nodal displacements are compared with SparLab's `displacement_<load
case>.csv` node by node. Reported per code and load case:

  * `max_abs_diff_m`, `max_rel_diff` = max |u_ref - u| / max |u|, and the RMS
    of the same ratio,
  * `ref_max_abs_m` and `sparlab_max_abs_m`,
  * pass/fail against the tolerance stated in the output.

Two things are deliberately measured rather than assumed. scikit-fem uses the
same element formulation, so it agrees to the level the two linear solvers
agree (1e-10 and better); the tolerance is 1e-7. CalculiX writes its nodal
results to the .frd file with six significant digits, so even for C3D8 - the
same element as SparLab's Hex8 - the comparison is bounded below by that
rounding: half a unit in the sixth digit of the largest displacement, 5e-6
relative. The same holds for C3D4, the linear tetrahedron, and for C3D10, the
quadratic one (four-point stiffness integration, as in SparLab). The CalculiX
tolerance is therefore 1e-5 for every element family, and the summary records
the rounding floor next to the measured difference. CalculiX's plane elements
(CPS4, CPS3) are *expanded* into solids through the thickness with the
plane-stress condition imposed on the expanded element, which is a different
discretisation of the plane problem. In plane stress it coincides with
SparLab's element only for a Poisson ratio of zero: measured on the Gmsh lug
bracket, CPS3 agrees to the .frd rounding floor at nu = 0 (2.8e-6) and differs
by about 1e-3 at nu = 0.33. A plane-stress comparison with nu != 0 is
therefore reported as an informational comparison between two idealisations,
and scikit-fem - the same plane element - is the verification for those
cases. The plane-strain expansion (CPE4, CPE3) holds the out-of-plane
displacement at zero and is exact: the two-material plate agrees to 2.1e-6
at nu = 0.3 / 0.33, so plane-strain comparisons are judged.

Loads beyond point loads and tractions are exported in CalculiX's own form -
pressure faces (P), self-weight (GRAV), body forces (BX/BY/BZ), rotation
(CENTRIF), and a temperature field through *EXPANSION and *TEMPERATURE - so the
CalculiX comparison of such a case also tests SparLab's integration of those
loads, not only its stiffness and solve. A case whose temperature is conducted
has a second deck, a steady *HEAT TRANSFER job that CalculiX solves itself;
its nodal temperatures (NT) are compared with SparLab's, relative to the
temperature range of the field.

scikit-fem makes two comparisons. It solves with SparLab's assembled load
vector (mesh.json), a check of the stiffness - per element material, for a
model with several - and of the solve; and, for a deck with such loads, it
integrates them itself from the same deck CalculiX reads ("scikit-fem loads"):
body forces and pressures with order-6 rules (exact on straight cells, as
SparLab's are), the thermal load with SparLab's rule (the element's stiffness
rule: four points for the Tet10, exact for the linear elements). Two
formulation choices of CalculiX differ from SparLab's and are reproduced in
scikit-fem so that CalculiX is judged against its own problem, with its
difference to SparLab recorded beside it: CalculiX takes the element-average
temperature for the thermal strain of a first-order hexahedron (C3D8, and the
hexahedra of an expanded CPS4 / CPE4), and integrates the centrifugal load of
a C3D10 with its four-point rule, which is not exact for that cubic integrand
(SparLab's consistent-mass integration is). Both were measured, not assumed:
reproducing them brings CalculiX to within 1e-6 to 4e-6 of scikit-fem. A third
shows on curved geometry only: CalculiX integrates a pressure on the six-node
face of a C3D10 with a three-point rule, exact on a flat face but not for the
degree-4 integrand of a curved one, whose load SparLab integrates exactly.

Nothing here recomputes SparLab's numbers: they are read from the run.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from typing import Dict, List, Optional

import _bootstrap  # noqa: F401

import numpy as np

from sparlab_viz.loaders import Mesh, ResultError, load_case


# ---------------------------------------------------------------------------
# CalculiX
# ---------------------------------------------------------------------------
def parse_frd_displacements(path: str, num_nodes: int) -> np.ndarray:
    """Nodal displacements (num_nodes, 3) from a CalculiX .frd file (last step)."""
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()
    blocks: List[np.ndarray] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith(" -4") and "DISP" in line:
            values = np.full((num_nodes, 3), np.nan)
            i += 1
            # component definition lines start with " -5"
            while i < len(lines) and lines[i].startswith(" -5"):
                i += 1
            while i < len(lines) and lines[i].startswith(" -1"):
                row = lines[i]
                node = int(row[3:13])
                comps = [float(row[13 + 12 * k: 25 + 12 * k]) for k in range(3)]
                values[node - 1] = comps
                i += 1
            blocks.append(values)
            continue
        i += 1
    if not blocks:
        raise ResultError(f"{path} holds no displacement block")
    return blocks[-1]


def parse_frd_temperatures(path: str, num_nodes: int) -> np.ndarray:
    """Nodal temperatures (num_nodes,) from the NDTEMP block of a CalculiX .frd
    file (last step)."""
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()
    blocks: List[np.ndarray] = []
    i = 0
    while i < len(lines):
        if lines[i].startswith(" -4") and "NDTEMP" in lines[i]:
            values = np.full(num_nodes, np.nan)
            i += 1
            while i < len(lines) and lines[i].startswith(" -5"):
                i += 1
            while i < len(lines) and lines[i].startswith(" -1"):
                row = lines[i]
                values[int(row[3:13]) - 1] = float(row[13:25])
                i += 1
            blocks.append(values)
            continue
        i += 1
    if not blocks:
        raise ResultError(f"{path} holds no temperature block")
    return blocks[-1]


def frd_rounding_floor(values: np.ndarray) -> float:
    """Half a unit in the sixth significant digit of the largest |value|: the
    resolution of a .frd field."""
    top = float(np.abs(values).max())
    if top == 0.0:
        return 0.0
    return 0.5 * 10.0 ** (np.floor(np.log10(top)) - 5.0)


def run_calculix(deck: str, workdir: str) -> str:
    """Run ccx on `deck` inside `workdir` and return the .frd path."""
    job = os.path.splitext(os.path.basename(deck))[0]
    shutil.copy(deck, os.path.join(workdir, job + ".inp"))
    env = dict(os.environ, OMP_NUM_THREADS="1")
    proc = subprocess.run(["ccx", "-i", job], cwd=workdir, capture_output=True,
                          text=True, env=env, check=False)
    frd = os.path.join(workdir, job + ".frd")
    if proc.returncode != 0 or not os.path.isfile(frd):
        raise ResultError(f"ccx failed on {deck}:\n{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}")
    return frd


def run_calculix_nlgeom(deck: str, workdir: str) -> str:
    """Run an NLGEOM deck and return the .frd path - only if CalculiX completed
    its step. A run that cuts back too often still writes its "best solution"
    to the .frd file; the .sta file says whether the last increment converged
    and reached the step's end, and nothing else is accepted."""
    frd = run_calculix(deck, workdir)
    sta = os.path.join(workdir, os.path.splitext(os.path.basename(deck))[0] + ".sta")
    with open(sta, "r", encoding="utf-8", errors="replace") as handle:
        rows = [line.split() for line in handle if line.strip()]
    last = rows[-1] if rows else []
    completed = (len(last) >= 7 and last[0].isdigit() and "U" not in last[2]
                 and abs(float(last[5]) - 1.0) <= 1e-9)
    if not completed:
        raise ResultError(f"CalculiX did not complete the NLGEOM step of {deck} (last .sta "
                          f"row: {' '.join(last)})")
    return frd


def parse_dat_buckling_factors(path: str) -> np.ndarray:
    """Buckling factors from the "B U C K L I N G   F A C T O R" block of a
    CalculiX .dat file."""
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()
    factors: List[float] = []
    for i, line in enumerate(lines):
        if "B U C K L I N G" in line:
            for row in lines[i + 1:]:
                parts = row.split()
                if len(parts) == 2 and parts[0].isdigit():
                    factors.append(float(parts[1]))
                elif factors and not row.strip():
                    break
            break
    if not factors:
        raise ResultError(f"{path} holds no buckling factors")
    return np.asarray(factors)


def run_calculix_buckling(static_deck: str, workdir: str, num_modes: int) -> np.ndarray:
    """The static deck turned into a *BUCKLE step (same mesh, supports and
    nodal loads, which define the load pattern) and solved by ccx."""
    with open(static_deck, "r", encoding="utf-8") as handle:
        text = handle.read()
    if "*STATIC\n" not in text:
        raise ResultError(f"{static_deck} has no *STATIC step to turn into *BUCKLE")
    text = text.replace("*STATIC\n", f"*BUCKLE\n{num_modes}, 1.E-12, "
                                      f"{max(4 * num_modes, 40)}, 10000\n", 1)
    deck = os.path.join(workdir, "buckle.inp")
    with open(deck, "w", encoding="utf-8") as handle:
        handle.write(text)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    proc = subprocess.run(["ccx", "-i", "buckle"], cwd=workdir, capture_output=True,
                          text=True, env=env, check=False)
    dat = os.path.join(workdir, "buckle.dat")
    if proc.returncode != 0 or not os.path.isfile(dat):
        raise ResultError(f"ccx *BUCKLE failed on {static_deck}:\n{proc.stdout[-2000:]}")
    return parse_dat_buckling_factors(dat)


def calculix_version() -> str:
    """The ccx version number ("2.21"), or its raw banner if it cannot be parsed."""
    proc = subprocess.run(["ccx", "-v"], capture_output=True, text=True, check=False)
    text = (proc.stdout + proc.stderr).strip()
    match = re.search(r"[Vv]ersion\s+([0-9][0-9.]*)", text)
    if match:
        return match.group(1)
    lines = text.splitlines()
    return lines[-1].strip() if lines else "unknown"


# ---------------------------------------------------------------------------
# scikit-fem
# ---------------------------------------------------------------------------
class SkfemProblem:
    """The discrete problem of a SparLab run rebuilt in scikit-fem: basis,
    stiffness, nodal loads and prescribed DOFs in scikit-fem's numbering."""

    def __init__(self, mesh: Mesh, youngs, poisson, thickness: float,
                 stress_state: str, forces: np.ndarray, prescribed: List[Dict]):
        """`youngs` and `poisson` are numbers, or per-material sequences indexed
        by `mesh.element_materials` for a model with several materials."""
        from skfem import (Basis, ElementHex1, ElementQuad1, ElementTetP1, ElementTetP2,
                           ElementTriP1, ElementVector, MeshHex, MeshQuad, MeshTet, MeshTet2,
                           MeshTri)
        from skfem import asm
        from skfem.models.elasticity import lame_parameters, linear_elasticity

        dim = mesh.dim
        self.mesh = mesh
        self.dim = dim
        self.thickness = thickness if dim == 2 else 1.0
        # scikit-fem wants (dim, nodes) and (corners, elements) arrays, C-ordered.
        p = np.ascontiguousarray(mesh.nodes.T.astype(float))
        # SparLab node -> scikit-fem DOF per component; the identity map of the
        # nodal DOFs except for quadratic meshes, which scikit-fem renumbers.
        dof_of_node = None
        if mesh.element_type == "Tet10":
            # Quadratic tetrahedra, curved when the file's edge nodes are: an
            # isoparametric MeshTet2 built from all ten nodes of every cell.
            # With sort_t=False scikit-fem keeps the connectivity rows, and its
            # edge order - (0,1), (1,2), (2,0), (0,3), (1,3), (2,3) - is
            # SparLab's (VTK's), which the element_dofs check below confirms.
            # intorder=2 is the 4-point rule SparLab (and CalculiX's C3D10)
            # integrate the stiffness with, so curved cells are the same
            # discrete problem too.
            m = MeshTet2(p, np.ascontiguousarray(mesh.elements.T), sort_t=False)
            element = ElementVector(ElementTetP2())
            basis = Basis(m, element, intorder=2)
            scalar = Basis(m, ElementTetP2(), intorder=2)
            located = m.doflocs[:, scalar.element_dofs].transpose(2, 1, 0)
            expected = mesh.nodes[mesh.elements]
            if not np.allclose(located, expected, rtol=0.0,
                               atol=1e-12 * float(np.abs(mesh.nodes).max())):
                raise ResultError("scikit-fem's Tet10 node order does not match SparLab's")
            corners = np.unique(mesh.elements[:, :4])
            dof_of_node = np.full((mesh.num_nodes, dim), -1, dtype=int)
            dof_of_node[corners] = basis.nodal_dofs.T
            scalar_dof = np.full(mesh.num_nodes, -1, dtype=int)
            scalar_dof[corners] = scalar.nodal_dofs[0]
            for j in range(6):
                dof_of_node[mesh.elements[:, 4 + j]] = basis.edge_dofs[:, m.t2e[j]].T
                scalar_dof[mesh.elements[:, 4 + j]] = scalar.edge_dofs[0, m.t2e[j]]
            if (dof_of_node < 0).any() or (scalar_dof < 0).any():
                raise ResultError("some Tet10 nodes have no scikit-fem DOF")
            # scikit-fem numbers the corner vertices in ascending node order.
            vertex_of_node = np.full(mesh.num_nodes, -1, dtype=int)
            vertex_of_node[corners] = np.arange(corners.size)
            scalar_element = ElementTetP2()
        elif mesh.element_type == "Tri3":
            # Linear simplices: any vertex order describes the same element.
            m = MeshTri(p, np.ascontiguousarray(mesh.elements.T))
            element = ElementVector(ElementTriP1())
            scalar_element = ElementTriP1()
        elif mesh.element_type == "Tet4":
            m = MeshTet(p, np.ascontiguousarray(mesh.elements.T))
            element = ElementVector(ElementTetP1())
            scalar_element = ElementTetP1()
        elif dim == 2:
            # scikit-fem's quad corners are ordered (0,0), (0,1), (1,1), (1,0) on
            # its reference square; SparLab stores them counter-clockwise, so
            # the second and fourth corners swap.
            t = np.ascontiguousarray(mesh.elements[:, [0, 3, 2, 1]].T)
            m = MeshQuad(p, t)
            element = ElementVector(ElementQuad1())
            scalar_element = ElementQuad1()
        else:
            # scikit-fem hex corners: (0,0,0), (0,1,0), (1,0,0), (0,0,1),
            # (1,1,0), (0,1,1), (1,0,1), (1,1,1) in reference coordinates;
            # SparLab uses the VTK order (bottom face counter-clockwise, then
            # the top face).
            t = np.ascontiguousarray(mesh.elements[:, [0, 3, 1, 4, 2, 7, 5, 6]].T)
            m = MeshHex(p, t)
            element = ElementVector(ElementHex1())
            scalar_element = ElementHex1()
        if dof_of_node is None:
            basis = Basis(m, element)
            dof_of_node = basis.nodal_dofs.T
            # Linear meshes keep SparLab's node numbering for their vertices.
            scalar_dof = np.arange(mesh.num_nodes)
            vertex_of_node = np.arange(mesh.num_nodes)
        self.basis = basis
        self.dof_of_node = dof_of_node
        self.skfem_mesh = m
        self.vector_element = element
        self.scalar_element = scalar_element
        self.scalar_dof_of_node = scalar_dof
        self.vertex_of_node = vertex_of_node
        self.stress_state = stress_state

        if stress_state not in ("plane_stress", "plane_strain", "three_dimensional"):
            raise ResultError(f"unsupported stress state {stress_state}")

        def lame(e: float, nu: float):
            lam_, mu_ = lame_parameters(e, nu)
            if stress_state == "plane_stress":
                lam_ = 2.0 * lam_ * mu_ / (lam_ + 2.0 * mu_)
            return lam_, mu_

        if np.ndim(youngs) == 0:
            lam, mu = lame(float(youngs), float(poisson))
            self.k = asm(linear_elasticity(lam, mu), basis) * self.thickness
            self.lam_e = np.full(mesh.num_elements, lam)
            self.mu_e = np.full(mesh.num_elements, mu)
        else:
            # Several materials: the Lame parameters as element-wise fields.
            from skfem import BilinearForm
            from skfem.helpers import ddot, sym_grad, trace, eye

            if mesh.element_materials is None:
                raise ResultError("several materials but mesh.json has no element_materials")
            pairs = np.asarray([lame(float(e), float(n)) for e, n in zip(youngs, poisson)])
            per_element = pairs[mesh.element_materials]            # (elements, 2)
            self.lam_e = per_element[:, 0].copy()
            self.mu_e = per_element[:, 1].copy()
            nqp = basis.X.shape[-1]
            lam = np.repeat(per_element[:, :1], nqp, axis=1)       # (elements, qp)
            mu = np.repeat(per_element[:, 1:], nqp, axis=1)
            dim_ = dim

            @BilinearForm
            def elasticity(u, v, w):
                eps_u = sym_grad(u)
                sigma = w["lam"] * eye(trace(eps_u), dim_) + 2.0 * w["mu"] * eps_u
                return ddot(sigma, sym_grad(v))

            self.k = asm(elasticity, basis, lam=lam, mu=mu) * self.thickness
        self.lam, self.mu = lam, mu

        # Nodal loads and prescribed DOFs mapped through skfem's DOF numbering.
        self.f = np.zeros(basis.N)
        for row in forces:
            node = int(row[0])
            for c in range(dim):
                self.f[dof_of_node[node, c]] += row[1 + c]
        self.x = np.zeros(basis.N)
        fixed = []
        for entry in prescribed:
            c = "xyz".index(entry["component"])
            dof = dof_of_node[int(entry["node"]), c]
            fixed.append(dof)
            self.x[dof] = float(entry["value_m"])
        self.fixed = np.unique(np.asarray(fixed, dtype=int))
        self.free = np.setdiff1d(np.arange(basis.N), self.fixed)
        self._factor = None
        self._k_fp = None

    def static(self, f: Optional[np.ndarray] = None) -> np.ndarray:
        """The static solution in scikit-fem's numbering, for the load vector
        of mesh.json or for `f`: K_ff u_f = f_f - K_fp x_p, with K_ff factorised
        once (SuperLU) and reused for every load vector."""
        import scipy.sparse.linalg as sla

        if self._factor is None:
            rows = self.k[self.free]
            self._factor = sla.factorized(rows[:, self.free].tocsc())
            self._k_fp = rows[:, self.fixed]
        load = self.f if f is None else f
        u = self.x.copy()
        u[self.free] = self._factor(load[self.free] - self._k_fp @ self.x[self.fixed])
        return u

    def nodal(self, u: np.ndarray) -> np.ndarray:
        out = np.zeros((self.mesh.num_nodes, self.dim))
        for c in range(self.dim):
            out[:, c] = u[self.dof_of_node[:, c]]
        return out

    def nonlinear(self, law: str, load_factors: List[float], f: Optional[np.ndarray] = None,
                  tolerance: float = 1e-11, max_iterations: int = 60) -> np.ndarray:
        """The geometrically non-linear static solution under dead loads, by
        an independent total Lagrangian implementation: the internal force
        int P : grad(v) with the first Piola-Kirchhoff stress P = F S and its
        consistent tangent, written here with scikit-fem's tensor helpers, on
        SparLab's quadrature (2 x 2 (x 2) Gauss points for Q4 / Hex8, the
        four-point rule for Tet10; exact for the constant integrand of the
        linear simplices). `law` is "saint_venant_kirchhoff"
        (S = lambda tr E I + 2 mu E, with the plane-stress lambda* in plane
        stress) or "neo_hookean" (S = mu (I - C^-1) + lambda ln J C^-1, plane
        strain or 3-D). Newton's method converges each load factor of
        `load_factors` in turn - the external load and the prescribed
        displacements scale with it - to a residual of `tolerance` relative to
        the load, with a line search on the energy (the slack criterion that
        SparLab uses). Returns the solution at the last factor in
        scikit-fem's numbering."""
        import scipy.sparse.linalg as sla
        from skfem import Basis, BilinearForm, LinearForm, asm
        from skfem.helpers import ddot, det, eye, grad, inv, mul, trace, transpose

        dim = self.dim
        element = self.mesh.element_type
        intorder = {"Quad4": 3, "Hex8": 3, "Tri3": 1, "Tet4": 1, "Tet10": 2}[element]
        basis = Basis(self.skfem_mesh, self.vector_element, intorder=intorder)
        if not np.array_equal(basis.nodal_dofs, self.basis.nodal_dofs):
            raise ResultError("the non-linear basis numbers its DOFs differently")
        nqp = basis.X.shape[-1]
        lam = np.repeat(self.lam_e[:, None], nqp, axis=1)
        mu = np.repeat(self.mu_e[:, None], nqp, axis=1)
        if law == "neo_hookean" and self.stress_state == "plane_stress":
            raise ResultError("the neo-Hookean law has no plane-stress form")

        def kinematics(w):
            h = grad(w["u"])
            ident = eye(np.ones_like(h[0, 0]), dim)
            return ident, ident + h

        def stress(ident, f, w):
            c = mul(transpose(f), f)
            if law == "saint_venant_kirchhoff":
                e = 0.5 * (c - ident)
                return w["lam"] * eye(trace(e), dim) + 2.0 * w["mu"] * e, c
            ci = inv(c)
            lnj = np.log(det(f))
            return w["mu"] * (ident - ci) + w["lam"] * lnj * ci, c

        @LinearForm
        def internal(v, w):
            ident, f = kinematics(w)
            s, _ = stress(ident, f, w)
            return ddot(mul(f, s), grad(v))

        @BilinearForm
        def tangent(du, v, w):
            ident, f = kinematics(w)
            s, c = stress(ident, f, w)
            df = grad(du)
            dc = mul(transpose(df), f) + mul(transpose(f), df)
            if law == "saint_venant_kirchhoff":
                de = 0.5 * dc
                ds = w["lam"] * eye(trace(de), dim) + 2.0 * w["mu"] * de
            else:
                ci = inv(c)
                lnj = np.log(det(f))
                ds = ((w["mu"] - w["lam"] * lnj) * mul(mul(ci, dc), ci)
                      + 0.5 * w["lam"] * trace(mul(ci, dc)) * ci)
            return ddot(mul(df, s) + mul(f, ds), grad(v))

        load = self.f if f is None else f
        free, fixed = self.free, self.fixed
        t = self.thickness

        def residual(u, factor):
            return asm(internal, basis, u=basis.interpolate(u), lam=lam, mu=mu) * t - factor * load

        def stiffness(u):
            with np.errstate(invalid="ignore", divide="ignore"):
                return asm(tangent, basis, u=basis.interpolate(u), lam=lam, mu=mu) * t

        def residual_or_nan(u, factor):
            with np.errstate(invalid="ignore", divide="ignore"):
                return residual(u, factor)

        u = np.zeros(basis.N)
        scale = max(np.linalg.norm(load), 1e-300)
        eps = np.finfo(float).eps
        previous = 0.0
        for factor in load_factors:
            # Tangent predictor from the converged state: K du = dl f - R with
            # the prescribed increment on its DOFs.
            k = stiffness(u)
            step = np.zeros(basis.N)
            step[fixed] = (factor - previous) * self.x[fixed]
            rhs = ((factor - previous) * load - residual(u, previous) - k @ step)[free]
            step[free] = sla.spsolve(k[free][:, free].tocsc(), rhs)
            u = u + step
            u[fixed] = factor * self.x[fixed]
            previous = factor
            if not np.linalg.norm(load) > 0.0:
                scale = max(np.linalg.norm(residual(u, factor)[fixed]), 1e-300)
            floor = 0.0
            for _ in range(max_iterations):
                r = residual_or_nan(u, factor)
                if not np.all(np.isfinite(r)):
                    raise ResultError(f"scikit-fem's Newton reached an inverted element at load "
                                      f"factor {factor}")
                rf = r[free]
                # Converged at the tolerance, or at the round-off floor the
                # rounding of u itself leaves: 64 eps || |K| |u| || on the free
                # DOFs, where Newton's residual stagnates on a slender structure.
                if np.linalg.norm(rf) <= max(tolerance * max(scale, np.linalg.norm(r[fixed])),
                                             floor):
                    break
                k = stiffness(u)
                floor = 64.0 * eps * np.linalg.norm((abs(k) @ np.abs(u))[free])
                du = np.zeros(basis.N)
                du[free] = sla.spsolve(k[free][:, free].tocsc(), -rf)
                # A step into an inverted element is halved until it is valid;
                # then the line search on the energy: the full step unless it
                # leaves |g(1)| > 0.8 |g(0)|, then regula falsi on
                # g(a) = du . R(u + a du).
                alpha = 1.0
                trial = residual_or_nan(u + du, factor)
                while not np.all(np.isfinite(trial)) and alpha > 1e-4:
                    alpha *= 0.5
                    trial = residual_or_nan(u + alpha * du, factor)
                du = alpha * du
                alpha = 1.0
                g0 = du[free] @ rf
                if g0 < 0.0:
                    g1 = du[free] @ trial[free]
                    if abs(g1) > 0.8 * abs(g0) and g1 > 0.0:
                        lo, glo, hi, ghi = 0.0, g0, 1.0, g1
                        for _ in range(5):
                            a = min(max(lo - glo * (hi - lo) / (ghi - glo), 0.1), 1.0)
                            ga = du[free] @ residual(u + a * du, factor)[free]
                            alpha = a
                            if abs(ga) <= 0.8 * abs(g0):
                                break
                            if ga < 0.0:
                                lo, glo = a, ga
                            else:
                                hi, ghi = a, ga
                u = u + alpha * du
            else:
                raise ResultError(f"scikit-fem's Newton did not converge at load factor {factor}")
        return u

    def buckling(self, u: np.ndarray, num_modes: int) -> np.ndarray:
        """Lowest positive load factors of (K + lambda K_G(u)) phi = 0 on the
        free DOFs, K_G assembled here from the stress of `u` at the same
        quadrature points, by a dense generalised eigensolve."""
        import scipy.linalg
        from skfem import BilinearForm, asm
        from skfem.helpers import eye, grad, sym_grad, trace

        dim = self.dim
        eps = sym_grad(self.basis.interpolate(u))
        sigma = self.lam * eye(trace(eps), dim) + 2.0 * self.mu * eps

        @BilinearForm
        def geometric(du, dv, w):
            # sum_ijk sigma_ij d(du_k)/dx_i d(dv_k)/dx_j
            return np.einsum("ij...,ki...,kj...->...", w["s"], grad(du), grad(dv))

        kg = asm(geometric, self.basis, s=sigma) * self.thickness
        free = self.free
        kff = self.k[free][:, free].toarray()
        gff = kg[free][:, free].toarray()
        # 1/lambda from (-K_G) phi = (1/lambda) K phi, K positive definite.
        inverse = scipy.linalg.eigh(-gff, kff, eigvals_only=True)
        positive = inverse[inverse > 0.0]
        return np.sort(1.0 / positive)[:num_modes]


# ---------------------------------------------------------------------------
# The loads of a CalculiX deck, integrated by scikit-fem
# ---------------------------------------------------------------------------
# Corner nodes (0-based, local) of CalculiX's faces 1, 2, ... (manual, 6.2).
CALCULIX_FACES = {
    "Quad4": [(0, 1), (1, 2), (2, 3), (3, 0)],
    "Tri3": [(0, 1), (1, 2), (2, 0)],
    "Hex8": [(0, 1, 2, 3), (4, 7, 6, 5), (0, 4, 5, 1), (1, 5, 6, 2), (2, 6, 7, 3),
             (3, 7, 4, 0)],
    "Tet4": [(0, 1, 2), (0, 3, 1), (1, 3, 2), (2, 3, 0)],
    "Tet10": [(0, 1, 2), (0, 3, 1), (1, 3, 2), (2, 3, 0)],
}


class DeckLoads:
    """The loads of an exported CalculiX static deck, as CalculiX reads them:
    *CLOAD nodal forces, *DLOAD pressure faces, GRAV, BX/BY/BZ on element
    sets and CENTRIF, and the nodal *TEMPERATURE field."""

    def __init__(self, path: str, num_nodes: int, num_elements: int):
        self.cload: List[tuple] = []       # (node, component, value), 0-based
        self.pressure: List[tuple] = []    # (element, face 1-based, pressure)
        self.gravity: Optional[np.ndarray] = None
        self.body: List[tuple] = []        # (element set, component, value)
        self.centrifugal: Optional[tuple] = None  # (omega^2, point, axis)
        self.temperature: Optional[np.ndarray] = None
        self.elsets: Dict[str, np.ndarray] = {"EALL": np.arange(num_elements)}
        keyword, name, ids = "", "", []
        with open(path, "r", encoding="utf-8") as handle:
            lines = [line.strip() for line in handle]
        for line in lines + ["*END"]:
            if not line:
                continue
            if line.startswith("*"):
                if keyword == "*ELSET":
                    self.elsets[name] = np.asarray(ids, dtype=int) - 1
                keyword = line.split(",")[0].upper()
                name, ids = "", []
                if keyword == "*ELSET":
                    name = line.split("ELSET=")[1].split(",")[0].strip()
                if keyword == "*TEMPERATURE":
                    self.temperature = np.full(num_nodes, np.nan)
                continue
            parts = [p.strip() for p in line.split(",") if p.strip()]
            if keyword == "*ELSET":
                ids.extend(int(p) for p in parts)
            elif keyword == "*CLOAD":
                self.cload.append((int(parts[0]) - 1, int(parts[1]) - 1, float(parts[2])))
            elif keyword == "*TEMPERATURE":
                self.temperature[int(parts[0]) - 1] = float(parts[1])
            elif keyword == "*DLOAD":
                kind = parts[1].upper()
                values = [float(v) for v in parts[2:]]
                if kind.startswith("P") and kind[1:].isdigit():
                    self.pressure.append((int(parts[0]) - 1, int(kind[1:]), values[0]))
                elif kind == "GRAV":
                    self.gravity = values[0] * np.asarray(values[1:4])
                elif kind in ("BX", "BY", "BZ"):
                    self.body.append((parts[0], "XYZ".index(kind[1]), values[0]))
                elif kind == "CENTRIF":
                    self.centrifugal = (values[0], np.asarray(values[1:4]),
                                        np.asarray(values[4:7]))
                else:
                    raise ResultError(f"{path}: unsupported *DLOAD type {kind}")
        if self.temperature is not None and np.isnan(self.temperature).any():
            raise ResultError(f"{path}: *TEMPERATURE leaves some nodes without a value")

    @property
    def native(self) -> bool:
        """True when the deck carries any load in CalculiX's own form."""
        return bool(self.pressure or self.body or self.gravity is not None
                    or self.centrifugal is not None or self.temperature is not None)


def sparlab_thermal_intorder(element_type: str) -> int:
    """The rule SparLab integrates the thermal load with: the element's
    stiffness rule, which is the four-point rule (order 2) for the Tet10 and
    exact for the thermal integrand of every linear element."""
    return 2 if element_type == "Tet10" else 6


def native_load_vector(problem: "SkfemProblem", loads: DeckLoads, materials: List[Dict],
                       element_materials: Optional[np.ndarray], thermal: str,
                       body_intorder: int = 6, pressure_intorder: int = 6,
                       thermal_intorder: Optional[int] = None) -> np.ndarray:
    """The load vector of `loads` integrated by scikit-fem, in its numbering.

    Point loads and tractions are the deck's nodal forces; everything else is
    integrated here: the body force density
    rho g + rho omega^2 (I - e e^T)(x - c) + b over the volume with a rule of
    order `body_intorder` (6 is exact on straight cells, as SparLab's
    consistent-mass integration is), the pressure as -p n over the faces with
    a rule of order `pressure_intorder` (6 is exact for the degree-4 integrand
    of a curved six-node face, as SparLab's rule is), and the thermal load int eps(v) : sigma_0 with
    sigma_0 = k_th alpha dT I (k_th = 3 lambda + 2 mu in 3-D and plane strain,
    2 lambda* + 2 mu in plane stress) with a rule of order `thermal_intorder`
    (by default SparLab's, `sparlab_thermal_intorder`). `thermal` is
    "interpolated" - the temperature interpolated by the shape functions, the
    consistent load - or "element_average", each cell at the mean of its nodal
    temperatures, which is how CalculiX treats first-order hexahedra.
    """
    from skfem import Basis, FacetBasis, LinearForm, asm
    from skfem.helpers import div, dot

    dim, mesh, t = problem.dim, problem.mesh, problem.thickness
    m = problem.skfem_mesh
    if thermal_intorder is None:
        thermal_intorder = sparlab_thermal_intorder(mesh.element_type)
    basis = Basis(m, problem.vector_element, intorder=body_intorder)
    nqp = basis.X.shape[-1]
    ne = mesh.num_elements
    f = np.zeros(basis.N)
    for node, component, value in loads.cload:
        f[problem.dof_of_node[node, component]] += value

    index = element_materials if element_materials is not None else np.zeros(ne, dtype=int)
    rho = np.asarray([float(mat.get("density_kg_per_m3", 0.0)) for mat in materials])[index]
    b = np.zeros((dim, ne, nqp))
    if loads.gravity is not None:
        b += rho[None, :, None] * loads.gravity[:dim, None, None]
    if loads.centrifugal is not None:
        omega2, point, axis = loads.centrifugal
        axis = axis / np.linalg.norm(axis)
        x = np.zeros((3, ne, nqp))
        x[:dim] = basis.global_coordinates().value
        r = x - point[:, None, None]
        perp = r - axis[:, None, None] * np.einsum("i,ijk->jk", axis, r)[None]
        b += omega2 * rho[None, :, None] * perp[:dim]
    for elset, component, value in loads.body:
        b[component, loads.elsets[elset], :] += value
    if np.any(b != 0.0):
        @LinearForm
        def body(v, w):
            return dot(w["b"], v)

        f += asm(body, basis, b=b) * t

    if loads.pressure:
        faces = CALCULIX_FACES[mesh.element_type]
        facet_of = {tuple(sorted(col)): i for i, col in enumerate(m.facets.T)}
        chosen, values = [], []
        for element, face, pressure in loads.pressure:
            corners = mesh.elements[element][list(faces[face - 1])]
            key = tuple(sorted(problem.vertex_of_node[corners]))
            if key not in facet_of:
                raise ResultError(f"pressure face {face} of element {element + 1} is not a "
                                  "facet of the scikit-fem mesh")
            chosen.append(facet_of[key])
            values.append(pressure)
        chosen = np.asarray(chosen)
        fb = FacetBasis(m, problem.vector_element, facets=chosen, intorder=pressure_intorder)
        pressure = np.repeat(np.asarray(values)[:, None], fb.X.shape[-1], axis=1)

        @LinearForm
        def load(v, w):
            return -w["p"] * dot(w.n, v)

        f += asm(load, fb, p=pressure) * t

    if loads.temperature is not None:
        alpha = np.asarray([float(mat.get("thermal_expansion_per_K", 0.0))
                            for mat in materials])[index]
        references = {float(mat.get("reference_temperature_K", 0.0)) for mat in materials
                      if float(mat.get("thermal_expansion_per_K", 0.0)) != 0.0}
        if len(references) > 1:
            raise ResultError("materials with different reference temperatures")
        t_ref = references.pop() if references else 0.0
        tbasis = Basis(m, problem.vector_element, intorder=thermal_intorder)
        tqp = tbasis.X.shape[-1]
        if thermal == "interpolated":
            scalar = Basis(m, problem.scalar_element, intorder=thermal_intorder)
            nodal = np.zeros(scalar.N)
            nodal[problem.scalar_dof_of_node] = loads.temperature
            dt = scalar.interpolate(nodal).value - t_ref
        elif thermal == "element_average":
            mean = loads.temperature[mesh.elements].mean(axis=1) - t_ref
            dt = np.repeat(mean[:, None], tqp, axis=1)
        else:
            raise ValueError(thermal)
        lam, mu = problem.lam_e, problem.mu_e
        k_th = (2.0 * lam + 2.0 * mu) if problem.stress_state == "plane_stress" \
            else (3.0 * lam + 2.0 * mu)
        coefficient = np.repeat((k_th * alpha)[:, None], tqp, axis=1)

        @LinearForm
        def thermal_load(v, w):
            return w["c"] * w["dt"] * div(v)

        f += asm(thermal_load, tbasis, c=coefficient, dt=dt) * t
    return f


def calculix_formulation(ccx_type: str, loads: Optional[DeckLoads]) -> Optional[Dict]:
    """Where CalculiX's treatment of a deck's loads departs from SparLab's, the
    options that make scikit-fem reproduce CalculiX's problem instead; None
    when the two formulations coincide. Measured (docs/verification.md):
      * CalculiX evaluates the thermal strain of a first-order hexahedron -
        C3D8, and the hexahedra it expands CPS4 / CPE4 into - at the
        element-average temperature; SparLab integrates the interpolated
        temperature, the consistent load. They agree for a uniform field.
      * CalculiX integrates the centrifugal load of a C3D10 with its four-point
        rule, which is not exact for the cubic integrand N_a rho omega^2 r;
        SparLab's consistent-mass integration is exact.
      * CalculiX integrates a pressure on the six-node face of a C3D10 with a
        three-point rule, exact on a flat face but not for the degree-4
        integrand of a curved one; SparLab's rule is exact for it.
    """
    if loads is None:
        return None
    reasons = []
    options = {"thermal": "interpolated", "body_intorder": 6, "pressure_intorder": 6}
    if (loads.temperature is not None and np.ptp(loads.temperature) > 0.0
            and ccx_type in ("C3D8", "CPS4", "CPE4")):
        options["thermal"] = "element_average"
        reasons.append("element-average temperature for the thermal strain of "
                       "first-order hexahedra")
    if loads.centrifugal is not None and ccx_type == "C3D10":
        options["body_intorder"] = 2
        reasons.append("four-point rule for the centrifugal load of the C3D10")
    if loads.pressure and ccx_type == "C3D10":
        options["pressure_intorder"] = 2
        reasons.append("three-point rule for a pressure on the six-node faces of the C3D10")
    if not reasons:
        return None
    options["reason"] = "; ".join(reasons)
    return options


def solve_with_skfem(mesh: Mesh, youngs: float, poisson: float, thickness: float,
                     stress_state: str, forces: np.ndarray, prescribed: List[Dict],
                     ) -> np.ndarray:
    """Solve the same problem with scikit-fem; returns (num_nodes, dim)."""
    problem = SkfemProblem(mesh, youngs, poisson, thickness, stress_state, forces,
                           prescribed)
    return problem.nodal(problem.static())


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------
def compare(reference: np.ndarray, ours: np.ndarray) -> Dict[str, float]:
    diff = reference - ours
    scale = float(np.abs(ours).max())
    per_node = np.linalg.norm(diff, axis=1) / max(scale, 1e-300)
    return {
        "max_abs_diff_m": float(np.abs(diff).max()),
        "max_rel_diff": float(np.abs(diff).max() / max(scale, 1e-300)),
        "rms_rel_diff": float(np.sqrt(np.mean(per_node ** 2))),
        "ref_max_abs_m": float(np.abs(reference).max()),
        "sparlab_max_abs_m": scale,
    }


def sparlab_displacement(case, load_case: str) -> np.ndarray:
    table = case.displacement(load_case)
    cols = ["ux[m]", "uy[m]", "uz[m]"][: case.dim]
    return table[cols].to_numpy()


def case_nonlinear_displacement(case, load_case: str) -> np.ndarray:
    """SparLab's final non-linear displacement (nonlinear_displacement_<lc>.csv)."""
    table = case.table(f"nonlinear_displacement_{_safe(load_case)}.csv")
    cols = ["ux[m]", "uy[m]", "uz[m]"][: case.dim]
    return table[cols].to_numpy()


def cross_validate_case(case_dir: str, tolerances: Dict[str, float],
                        skip_calculix: bool) -> Dict:
    case = load_case(case_dir)
    mesh = case.mesh
    summary = case.summary
    material = summary["material"]
    thickness = float(summary["mesh"].get("thickness_m", 1.0))
    stress_state = material["stress_state"]
    element_type = summary["mesh"]["element_type"]
    report = {
        "case": case.name,
        "directory": case_dir,
        "dim": mesh.dim,
        "element_type": element_type,
        "num_nodes": mesh.num_nodes,
        "num_elements": mesh.num_elements,
        "stress_state": stress_state,
        "load_cases": [],
    }
    plane_strain = stress_state == "plane_strain"
    ccx_type = {"Quad4": "CPE4" if plane_strain else "CPS4",
                "Tri3": "CPE3" if plane_strain else "CPS3",
                "Hex8": "C3D8", "Tet4": "C3D4", "Tet10": "C3D10"}[element_type]

    buckling = {entry["load_case"]: entry
                for entry in (summary.get("buckling") or {}).get("load_cases", [])}
    nonlinear_cases = {entry["load_case"]: entry
                       for entry in (summary.get("nonlinear") or {}).get("load_cases", [])}
    for name in mesh.load_case_names:
        ours = sparlab_displacement(case, name)
        entry = {"load_case": name, "codes": {}}

        # --- scikit-fem ---
        forces = mesh.nodal_forces(name)
        materials = summary.get("materials")
        if materials:
            youngs = [float(m["youngs_modulus_Pa"]) for m in materials]
            poisson = [float(m["poisson_ratio"]) for m in materials]
        else:
            youngs = float(material["youngs_modulus_Pa"])
            poisson = float(material["poisson_ratio"])
        problem = SkfemProblem(mesh, youngs, poisson, thickness, stress_state,
                               forces, mesh.prescribed)
        u_sk = problem.static()
        sk = problem.nodal(u_sk)
        stats = compare(sk, ours)
        stats["tolerance"] = tolerances["skfem"]
        stats["passed"] = stats["max_rel_diff"] <= tolerances["skfem"]
        skfem_element = {"Quad4": "ElementQuad1", "Hex8": "ElementHex1",
                         "Tri3": "ElementTriP1", "Tet4": "ElementTetP1",
                         "Tet10": "ElementTetP2"}[element_type]
        stats["element"] = skfem_element
        stats["comparison"] = "SparLab's assembled load vector (mesh.json)"
        entry["codes"]["scikit-fem"] = stats

        # --- scikit-fem integrating the deck's loads itself ---
        deck = case.path(f"calculix_{_safe(name)}.inp")
        loads = (DeckLoads(deck, mesh.num_nodes, mesh.num_elements)
                 if os.path.isfile(deck) else None)
        material_list = summary.get("materials") or [material]
        if loads is not None and loads.native:
            f_native = native_load_vector(problem, loads, material_list,
                                          mesh.element_materials, "interpolated")
            stats = compare(problem.nodal(problem.static(f_native)), ours)
            stats["tolerance"] = tolerances["skfem_loads"]
            stats["passed"] = stats["max_rel_diff"] <= tolerances["skfem_loads"]
            stats["element"] = skfem_element
            stats["loads"] = native_loads(deck)
            stats["comparison"] = (
                "loads integrated by scikit-fem from the CalculiX deck's own "
                "definitions: body forces and pressures with order-6 rules, the "
                "thermal load (interpolated temperature) with SparLab's stiffness "
                "rule, point loads and tractions as nodal forces")
            entry["codes"]["scikit-fem loads"] = stats

        # --- CalculiX ---
        if not skip_calculix:
            if not os.path.isfile(deck):
                raise ResultError(f"{deck} is missing; rerun sparlab_solve with "
                                  "--export-calculix")
            with tempfile.TemporaryDirectory(prefix="sparlab_ccx_") as work:
                frd = run_calculix(deck, work)
                ref = parse_frd_displacements(frd, mesh.num_nodes)[:, : mesh.dim]
            if np.isnan(ref).any():
                raise ResultError(f"CalculiX returned no displacement for some nodes of {name}")
            key = ("calculix_solid" if ccx_type in ("C3D8", "C3D4", "C3D10")
                   else "calculix_plane")
            stats = compare(ref, ours)
            stats["tolerance"] = tolerances[key]
            stats["element"] = ccx_type
            stats["loads"] = native_loads(deck)
            poisson = max(float(m["poisson_ratio"])
                          for m in (summary.get("materials") or [material]))
            formulation = calculix_formulation(ccx_type, loads)
            if formulation is not None:
                # CalculiX's problem differs from SparLab's in the way stated:
                # CalculiX is judged against scikit-fem solving CalculiX's
                # problem, and its difference to SparLab is recorded.
                f_ccx = native_load_vector(problem, loads, material_list,
                                           mesh.element_materials, formulation["thermal"],
                                           body_intorder=formulation["body_intorder"],
                                           pressure_intorder=formulation["pressure_intorder"])
                vs = compare(ref, problem.nodal(problem.static(f_ccx)))
                stats["calculix_formulation"] = formulation["reason"]
                stats["vs_skfem_calculix_formulation"] = vs
                stats["max_rel_diff_judged"] = vs["max_rel_diff"]
            plane_stress_expansion = ccx_type in ("CPS4", "CPS3") and poisson != 0.0
            if formulation is not None and not plane_stress_expansion:
                stats["passed"] = vs["max_rel_diff"] <= tolerances[key]
                stats["comparison"] = (
                    f"CalculiX's formulation differs ({formulation['reason']}): judged "
                    "against scikit-fem solving CalculiX's problem; max_rel_diff is "
                    "SparLab's difference to CalculiX, which the scikit-fem loads check "
                    "accounts for")
            elif plane_stress_expansion:
                # Not the same discrete problem (see the module docstring):
                # recorded, not judged.
                stats["passed"] = None
                stats["comparison"] = (
                    f"different idealisation: CalculiX expands {ccx_type} into a 3-D "
                    f"layer, which matches plane stress only for nu = 0 (here nu = "
                    f"{poisson:g}); informational")
            else:
                stats["passed"] = stats["max_rel_diff"] <= tolerances[key]
                stats["comparison"] = "same discrete problem"
            # Six significant digits in the .frd file: half a unit in the last
            # digit of the largest value, relative to that value.
            stats["frd_rounding_floor_rel"] = 5.0e-6
            stats["within_frd_rounding"] = stats["max_rel_diff"] <= 5.0e-6
            entry["codes"]["calculix"] = stats

        # --- steady conduction, when the case conducted its temperature ---
        conduction_deck = case.path(f"calculix_{_safe(name)}_conduction.inp")
        if not skip_calculix and os.path.isfile(conduction_deck):
            with tempfile.TemporaryDirectory(prefix="sparlab_ccx_heat_") as work:
                frd = run_calculix(conduction_deck, work)
                ref_t = parse_frd_temperatures(frd, mesh.num_nodes)
            if np.isnan(ref_t).any():
                raise ResultError(f"CalculiX returned no temperature for some nodes of {name}")
            ours_t = case.temperature(name)["temperature[K]"].to_numpy()
            span = float(ours_t.max() - ours_t.min())
            diff = np.abs(ref_t - ours_t)
            floor = frd_rounding_floor(ref_t)
            with open(conduction_deck, "r", encoding="utf-8") as handle:
                heat_type = re.search(r"TYPE=(\w+)", handle.read()).group(1)
            stats = {
                "element": heat_type,
                "quantity": "nodal temperature",
                "max_abs_diff_K": float(diff.max()),
                "max_rel_diff": float(diff.max() / max(span, 1e-300)),
                "normalised_by": "the temperature range of SparLab's field",
                "temperature_range_K": span,
                "tolerance": tolerances["calculix_conduction"],
                "frd_rounding_floor_K": floor,
                "frd_rounding_floor_rel": floor / max(span, 1e-300),
                "comparison": "same discrete conduction problem, solved by CalculiX "
                              "(*HEAT TRANSFER, STEADY STATE)",
            }
            stats["passed"] = stats["max_rel_diff"] <= tolerances["calculix_conduction"]
            entry["codes"]["calculix conduction"] = stats

        # --- the geometrically non-linear run, when the case has one ---
        nl_case = nonlinear_cases.get(name)
        if nl_case is not None:
            # Every comparison below is made along SparLab's own load factors, so a
            # run that stopped short would be compared - and could agree - at a
            # partial load. The decks are built to reach lambda = 1; one that does
            # not is a regression, not a result to compare.
            if not nl_case.get("completed"):
                raise ResultError(f"SparLab's non-linear run of load case '{name}' did not "
                                  f"reach lambda = 1 ({nl_case.get('termination', '')}); the "
                                  "cross-validation compares the full load")
            nl_ours = case_nonlinear_displacement(case, name)
            law = summary["nonlinear"]["material_model"]
            follower = bool(summary["nonlinear"]["follower_pressure"])
            path_table = case.table(f"nonlinear_{_safe(name)}.csv")
            factors = path_table["load_factor[-]"].to_numpy().tolist()
            dead_only = loads is None or (loads.centrifugal is None and loads.temperature is None
                                          and not (loads.pressure and follower))
            if dead_only:
                u_nl = problem.nonlinear(law, factors)
                stats = compare(problem.nodal(u_nl), nl_ours)
                stats["tolerance"] = tolerances["skfem_nonlinear"]
                stats["passed"] = stats["max_rel_diff"] <= tolerances["skfem_nonlinear"]
                stats["element"] = entry["codes"]["scikit-fem"]["element"]
                stats["law"] = law
                stats["load_factors"] = len(factors)
                stats["comparison"] = (
                    "an independent total Lagrangian implementation (P = F S and its "
                    "consistent tangent in scikit-fem's tensor helpers, SparLab's quadrature), "
                    "Newton through SparLab's load factors, dead loads of mesh.json")
                entry["codes"]["scikit-fem non-linear"] = stats
            else:
                entry["skipped_nonlinear_skfem"] = (
                    "follower pressure, rotation and temperature are outside the dead-load "
                    "scikit-fem comparison; CalculiX NLGEOM and the exact solutions of "
                    "sparlab_verify cover them")
            nl_deck = case.path(f"calculix_{_safe(name)}_nlgeom.inp")
            if not skip_calculix and law != "saint_venant_kirchhoff":
                entry["skipped_nonlinear_calculix"] = (
                    "CalculiX's NEO HOOKE is a different strain energy from SparLab's "
                    "neo-Hookean law; the exact tube solution of sparlab_verify verifies it")
            elif not skip_calculix:
                if not os.path.isfile(nl_deck):
                    raise ResultError(f"{nl_deck} is missing; rerun sparlab_solve with "
                                      "--export-calculix")
                with tempfile.TemporaryDirectory(prefix="sparlab_ccx_nlgeom_") as work:
                    frd = run_calculix_nlgeom(nl_deck, work)
                    ref = parse_frd_displacements(frd, mesh.num_nodes)[:, : mesh.dim]
                if np.isnan(ref).any():
                    raise ResultError(f"CalculiX returned no displacement for some nodes of "
                                      f"{name} (NLGEOM)")
                stats = compare(ref, nl_ours)
                stats["tolerance"] = tolerances["calculix_nlgeom"]
                stats["element"] = f"{ccx_type} NLGEOM"
                stats["loads"] = native_loads(nl_deck)
                linear = entry["codes"].get("calculix", {})
                # Judged only where the linear decks already coincide - the same
                # discrete problem, CalculiX's formulation choices included - and
                # nothing but CalculiX's finite-strain thermal model differs.
                same_problem = (linear.get("passed") is True
                                and "vs_skfem_calculix_formulation" not in linear)
                if loads is not None and loads.temperature is not None:
                    stats["passed"] = None
                    stats["comparison"] = (
                        "different thermal model at finite strain: SparLab splits the thermal "
                        "stretch off multiplicatively (S = D (E - E_theta) / theta), CalculiX "
                        "does not (measured: 4.9 % lower thermal stress than the additive "
                        "split on a restrained cube at alpha dT = 0.05, 0.12 % from SparLab's); "
                        "informational")
                elif not same_problem:
                    stats["passed"] = None
                    stats["comparison"] = (
                        "the linear decks already differ in formulation ("
                        + str(linear.get("calculix_formulation", linear.get("comparison", "")))
                        + "); informational")
                else:
                    stats["passed"] = stats["max_rel_diff"] <= tolerances["calculix_nlgeom"]
                    stats["comparison"] = (
                        "same discrete problem, CalculiX's *STEP, NLGEOM (Saint "
                        "Venant-Kirchhoff *ELASTIC, follower *DLOAD pressure), step completed")
                stats["frd_rounding_floor_rel"] = 5.0e-6
                entry["codes"]["calculix NLGEOM"] = stats

        # --- linear buckling, when the run computed it ---
        ours_lf = np.asarray(buckling.get(name, {}).get("load_factors") or [], dtype=float)
        if ours_lf.size:
            count = int(ours_lf.size)
            if problem.free.size <= tolerances["skfem_buckling_max_dofs"]:
                ref = problem.buckling(u_sk, count)
                stats = compare_factors(
                    ref, ours_lf, tolerances["skfem_buckling"],
                    "K_G assembled by scikit-fem from its own static solution at the "
                    "same quadrature points; dense generalised eigensolve")
                stats["element"] = entry["codes"]["scikit-fem"]["element"] + " K_G"
                entry["codes"]["scikit-fem buckling"] = stats
            else:
                entry["skipped"] = (f"scikit-fem buckling skipped: {problem.free.size} free "
                                    "DOFs exceed the dense-eigensolve limit")
            if not skip_calculix and mesh.dim == 3:
                deck = case.path(f"calculix_{_safe(name)}.inp")
                with tempfile.TemporaryDirectory(prefix="sparlab_ccx_buckle_") as work:
                    ref = run_calculix_buckling(deck, work, count)
                stats = compare_factors(
                    ref, ours_lf, tolerances["calculix_buckling"],
                    f"CalculiX *BUCKLE on the exported {ccx_type} deck (same mesh, "
                    "supports and nodal loads); its own stress-stiffness evaluation, so "
                    "an independent implementation rather than the same discrete problem")
                stats["element"] = f"{ccx_type} *BUCKLE"
                entry["codes"]["calculix buckling"] = stats
        report["load_cases"].append(entry)
    return report


def native_loads(deck: str) -> List[str]:
    """The distributed loads a CalculiX deck applies in CalculiX's own form."""
    found = []
    with open(deck, "r", encoding="utf-8") as handle:
        text = handle.read()
    for label, pattern in (("pressure (P)", r"^\d+, P\d, "), ("self-weight (GRAV)", r", GRAV, "),
                           ("body force (BX/BY/BZ)", r", B[XYZ], "),
                           ("rotation (CENTRIF)", r", CENTRIF, "),
                           ("temperature (*TEMPERATURE)", r"^\*TEMPERATURE$")):
        if re.search(pattern, text, flags=re.MULTILINE):
            found.append(label)
    if re.search(r"^\*CLOAD$", text, flags=re.MULTILINE):
        found.append("assembled nodal forces (*CLOAD)")
    return found


def compare_factors(reference: np.ndarray, ours: np.ndarray, tolerance: float,
                    comparison: str) -> Dict:
    """Mode-by-mode relative difference of two sets of load factors."""
    count = min(reference.size, ours.size)
    if count == 0:
        raise ResultError("no load factors to compare")
    rel = np.abs(reference[:count] - ours[:count]) / np.abs(ours[:count])
    return {
        "element": "buckling",
        "reference_load_factors": reference[:count].tolist(),
        "sparlab_load_factors": ours[:count].tolist(),
        "max_rel_diff": float(rel.max()),
        "per_mode_rel_diff": rel.tolist(),
        "tolerance": tolerance,
        "passed": bool(rel.max() <= tolerance),
        "comparison": comparison,
    }


def _safe(name: str) -> str:
    return "".join(c if (c.isalnum() or c in "-_") else "_" for c in name) or "unnamed"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--case", action="append", required=True,
                        help="result directory of a sparlab_solve run (repeatable)")
    parser.add_argument("--output", default="results/cross_validation")
    parser.add_argument("--skip-calculix", action="store_true")
    parser.add_argument("--tol-skfem", type=float, default=1e-7)
    parser.add_argument("--tol-skfem-loads", type=float, default=1e-7,
                        help="loads integrated by scikit-fem vs SparLab")
    parser.add_argument("--tol-calculix-solid", type=float, default=1e-5)
    parser.add_argument("--tol-calculix-plane", type=float, default=1e-5)
    parser.add_argument("--tol-skfem-buckling", type=float, default=1e-7)
    parser.add_argument("--tol-calculix-buckling", type=float, default=1e-4)
    parser.add_argument("--tol-skfem-nonlinear", type=float, default=1e-7,
                        help="scikit-fem's total Lagrangian solve vs SparLab's non-linear one")
    parser.add_argument("--tol-calculix-nlgeom", type=float, default=1e-5,
                        help="CalculiX *STEP, NLGEOM vs SparLab's non-linear solution")
    parser.add_argument("--tol-calculix-conduction", type=float, default=1e-4,
                        help="max nodal temperature difference over the temperature range "
                             "(the .frd rounding of a temperature near 300 K is 5e-4 K)")
    parser.add_argument("--skfem-buckling-max-dofs", type=int, default=6000,
                        help="largest free-DOF count for the dense buckling eigensolve")
    args = parser.parse_args(argv)

    if not args.skip_calculix and shutil.which("ccx") is None:
        print("ccx (CalculiX) is not on PATH; pass --skip-calculix to compare with "
              "scikit-fem only", file=sys.stderr)
        return 2

    tolerances = {"skfem": args.tol_skfem, "skfem_loads": args.tol_skfem_loads,
                  "calculix_solid": args.tol_calculix_solid,
                  "calculix_plane": args.tol_calculix_plane,
                  "skfem_buckling": args.tol_skfem_buckling,
                  "calculix_buckling": args.tol_calculix_buckling,
                  "calculix_conduction": args.tol_calculix_conduction,
                  "skfem_nonlinear": args.tol_skfem_nonlinear,
                  "calculix_nlgeom": args.tol_calculix_nlgeom,
                  "skfem_buckling_max_dofs": args.skfem_buckling_max_dofs}
    import skfem
    summary = {
        "about": ("Node-by-node comparison of SparLab nodal displacements with "
                  "independent solutions of the same discrete problem"),
        "codes": {
            "scikit-fem": {"version": skfem.__version__,
                           "note": "same element formulations (bilinear/trilinear "
                                   "isoparametric with full integration, linear "
                                   "simplices, isoparametric quadratic tetrahedra with "
                                   "the 4-point rule); differences are linear-solver "
                                   "round-off"},
            "calculix": {"version": None if args.skip_calculix else calculix_version(),
                         "note": "C3D8, C3D4 and C3D10 are the same elements as "
                                 "SparLab's Hex8, Tet4 and Tet10; CPS4/CPS3 (CPE4/CPE3) are plane elements "
                                 "CalculiX expands through the thickness, a different "
                                 "discretisation of the plane problem. Nodal results are "
                                 "read from the .frd file, which carries six significant "
                                 "digits, so differences below 5e-6 relative are its "
                                 "rounding, not a disagreement"},
        },
        "tolerances": tolerances,
        "cases": [],
    }
    all_passed = True
    for case_dir in args.case:
        try:
            report = cross_validate_case(case_dir, tolerances, args.skip_calculix)
        except (ResultError, FileNotFoundError, KeyError) as error:
            print(f"cross-validation failed for {case_dir}: {error}", file=sys.stderr)
            return 1
        summary["cases"].append(report)
        for lc in report["load_cases"]:
            for code, stats in lc["codes"].items():
                if stats["passed"] is None:
                    verdict = "INFO (different idealisation)"
                else:
                    all_passed = all_passed and stats["passed"]
                    verdict = "PASS" if stats["passed"] else "FAIL"
                judged = stats.get("max_rel_diff_judged", stats["max_rel_diff"])
                extra = ("" if "max_rel_diff_judged" not in stats else
                         f" [vs CalculiX's own formulation; SparLab "
                         f"{stats['max_rel_diff']:.1e}]")
                print(f"  {report['case']:<28} {lc['load_case']:<14} {code:<18} "
                      f"{stats['element']:<12} max rel diff {judged:.3e} "
                      f"(tol {stats['tolerance']:.0e}) {verdict}{extra}")
    summary["all_passed"] = all_passed
    os.makedirs(args.output, exist_ok=True)
    path = os.path.join(args.output, "summary.json")
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    print(f"  wrote {path}")
    return 0 if all_passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
