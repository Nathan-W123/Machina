"""Figures for the verification studies, the design study and the benchmark.

These read the CSV / JSON files written by `sparlab_verify`, `sparlab_bench` and
the study script. Like `plots`, nothing here recomputes physics; the only derived
values are ratios and least-squares slopes, and each is labelled as such.
"""

from __future__ import annotations

import glob
import os
from typing import Dict, List, Optional, Tuple

import matplotlib.ticker
import numpy as np
import pandas as pd

from . import style as st
from .loaders import load_csv, load_json


# ---------------------------------------------------------------------------
# Verification: sensitivity finite differences
# ---------------------------------------------------------------------------
def plot_sensitivity_check(directory: str, path: str) -> str:
    """Analytical vs finite-difference gradients, and error vs step size."""
    steps = load_csv(os.path.join(directory, "sensitivity_steps.csv"))
    elements = load_csv(os.path.join(directory, "sensitivity_elements.csv"))
    summary = load_json(os.path.join(directory, "summary.json")).get("sensitivity", {})
    tolerance = summary.get("tolerance", 1e-5)
    best_step = summary.get("best_step")

    fig, axes = st.figure(7.4, 6.4, nrows=2, ncols=1)

    # --- panel 1: error vs step size (3 series: within the all-pairs cap) ---
    ax = axes[0]
    st.require_scatter_series(3, "error measures")
    included = steps["num_excluded"].iloc[0]
    for slot, (column, label) in enumerate(
        [
            ("max_relative_error[-]", "max relative error over elements"),
            ("rms_relative_error[-]", "RMS relative error over elements"),
            ("directional_relative_error[-]", "relative error of the directional derivative"),
        ]
    ):
        ax.plot(steps["step[-]"], steps[column], "o-", color=st.series_color(slot),
                label=label)
    ax.axhline(tolerance, color=st.INK_MUTED, linewidth=1.0, linestyle="--")
    ax.text(
        steps["step[-]"].min(), tolerance, f" pass tolerance {tolerance:g}",
        fontsize=7.8, color=st.INK_SECONDARY, va="bottom",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.invert_xaxis()
    ax.set_xlabel("central-difference step h on the design variable [-]")
    ax.set_ylabel("relative error [-]")
    st.title(
        ax, "Topology sensitivity: analytical gradient vs central differences",
        f"{int(steps['num_tested'].iloc[0])} elements tested, {int(included)} "
        "excluded at active bounds; the error floor is where round-off in the "
        "difference of two compliances overtakes the truncation error",
    )
    st.legend(ax, loc="upper left")

    # --- panel 2: element-by-element comparison at the best step ------------
    ax = axes[1]
    subset = elements
    if best_step is not None:
        subset = elements[np.isclose(elements["step[-]"], best_step)]
    tested = subset[subset["excluded[-]"] < 0.5]
    excluded = subset[subset["excluded[-]"] > 0.5]

    ax.plot(tested["analytical[J]"], tested["finite_difference[J]"], "o",
            color=st.series_color(0), markersize=4.5,
            label=f"tested elements ({len(tested)})")
    if not excluded.empty:
        ax.plot(excluded["analytical[J]"],
                np.zeros(len(excluded)), "x", color=st.series_color(1),
                markersize=5, label=f"excluded at bounds ({len(excluded)})")
    limits = [
        float(min(tested["analytical[J]"].min(), tested["finite_difference[J]"].min())),
        float(max(tested["analytical[J]"].max(), tested["finite_difference[J]"].max())),
    ]
    span = limits[1] - limits[0]
    limits = [limits[0] - 0.05 * span, limits[1] + 0.05 * span]
    ax.plot(limits, limits, "-", color=st.INK_MUTED, linewidth=1.0,
            label="exact agreement")
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_aspect("equal", adjustable="box")
    st.limit_ticks(ax, x=4, y=5)
    ax.set_xlabel(r"analytical $\partial c/\partial x_e$ [J]")
    ax.set_ylabel(r"central difference $\partial c/\partial x_e$ [J]")
    worst = float(tested["relative_error[-]"].max()) if not tested.empty else float("nan")
    st.title(
        ax, f"element-by-element agreement at h = {best_step:g}"
        if best_step else "element-by-element agreement",
        f"worst relative error {worst:.3e}",
    )
    st.legend(ax, loc="lower right")

    st.annotate_note(
        fig,
        "Excluded elements are passive variables (lower bound equals upper "
        "bound) or free variables within h of 0 or 1, where a central "
        "difference would leave the feasible box. They are plotted at zero on "
        "the vertical axis only to show how many there are.",
    )
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Verification: mesh convergence and beam theory
# ---------------------------------------------------------------------------
def plot_mesh_convergence(directory: str, path: str) -> str:
    """Tip deflection and discretisation error vs element size."""
    table = load_csv(os.path.join(directory, "mesh_convergence.csv"))
    summary = load_json(os.path.join(directory, "summary.json"))
    block = summary.get("mesh_convergence", {})
    orders = block.get("observed_convergence_order_tip_deflection", {})

    poissons = sorted(table["poisson[-]"].unique())
    fig, axes = st.figure(7.4, 6.6, nrows=2, ncols=1)

    # --- panel 1: tip deflection vs element size, with both references -----
    ax = axes[0]
    for slot, nu in enumerate(poissons):
        sub = table[table["poisson[-]"] == nu].sort_values("h[m]")
        ax.plot(sub["h[m]"], np.abs(sub["tip_mean[m]"]), "o-",
                color=st.series_color(slot), label=f"FEM, nu = {nu:g}")
        ax.axhline(abs(float(sub["timoshenko[m]"].iloc[0])),
                   color=st.series_color(slot), linewidth=0.9, linestyle=":")
    ax.axhline(abs(float(table["euler_bernoulli[m]"].iloc[0])),
               color=st.INK_MUTED, linewidth=1.0, linestyle="--",
               label="Euler-Bernoulli (bending only)")
    ax.text(
        float(table["h[m]"].min()), abs(float(table["timoshenko[m]"].iloc[0])),
        " Timoshenko (bending + shear), dotted per nu", fontsize=7.6,
        color=st.INK_SECONDARY, va="bottom",
    )
    ax.set_xscale("log")
    ax.invert_xaxis()
    ax.set_xlabel("element size h [m]")
    ax.set_ylabel("|tip deflection| [m]")
    st.title(
        ax, "Cantilever tip deflection vs mesh size",
        "L = 1 m, h = 0.1 m (L/h = 10), t = 10 mm, E = 70 GPa, 1 kN tip "
        "resultant; deflection is the mean u_y over the tip edge",
    )
    st.legend(ax, loc="lower left")

    # --- panel 2: self-convergence error, i.e. verification ----------------
    ax = axes[1]
    for slot, nu in enumerate(poissons):
        sub = table[table["poisson[-]"] == nu].sort_values("h[m]", ascending=False)
        finest = float(sub["tip_mean[m]"].iloc[-1])
        error = np.abs(sub["tip_mean[m]"] - finest) / abs(finest)
        mask = error > 0
        ax.plot(sub["h[m]"][mask], error[mask], "o-", color=st.series_color(slot),
                label=f"nu = {nu:g} (observed order "
                      f"{orders.get(f'{nu:g}', float('nan')):.2f})")
    # Second-order reference slope.
    h_ref = np.array([float(table["h[m]"].min()) * 2.0, float(table["h[m]"].max())])
    anchor = 1e-3
    ax.plot(h_ref, anchor * (h_ref / h_ref[0]) ** 2, "-", color=st.INK_MUTED,
            linewidth=1.0, label=r"reference slope $h^2$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.invert_xaxis()
    ax.set_xlabel("element size h [m]")
    ax.set_ylabel("relative error vs the finest mesh [-]")
    st.title(
        ax, "Self-convergence of the discretisation error",
        "measured against the finest mesh, which isolates discretisation error "
        "from the beam-theory modelling gap",
    )
    st.legend(ax, loc="lower left")

    st.annotate_note(
        fig,
        "Comparing against the finest mesh is verification (does the "
        "discretisation converge, and at what order); comparing against beam "
        "theory in the upper panel is validation (are the modelling assumptions "
        "appropriate), where a finite gap is expected because 1-D beam theory "
        "omits the Poisson coupling that plane stress retains.",
    )
    return st.save_figure(fig, path)


def plot_modal_convergence(directory: str, path: str) -> str:
    """Computed cantilever frequencies against Euler-Bernoulli theory."""
    table = load_csv(os.path.join(directory, "modal_convergence.csv"))
    fig, axes = st.stacked_panels(2, width=7.2, panel_height=2.5)

    ax = axes[0]
    series = [
        ("f1_fem[Hz]", "f1_theory[Hz]", "o", "bending mode 1"),
        ("f2_fem[Hz]", "f2_theory[Hz]", "s", "bending mode 2"),
        ("f_axial_fem[Hz]", "f_axial_theory[Hz]", "^", "axial mode 1"),
    ]
    for slot, (column, reference, marker, label) in enumerate(series):
        if column not in table:
            continue
        ax.plot(table["num_dofs"], table[column], marker + "-",
                color=st.series_color(slot), label=f"{label}, FEM")
        ax.axhline(float(table[reference].iloc[0]), color=st.series_color(slot),
                   linewidth=0.9, linestyle="--", label=f"{label}, theory")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("frequency [Hz]")
    st.title(
        ax, "Cantilever natural frequencies vs mesh size",
        "L = 1 m, h = 0.1 m, t = 10 mm, E = 70 GPa, rho = 2700 kg/m^3, nu = 0; "
        "consistent mass matrix. Bending references are Euler-Bernoulli; the "
        "axial reference is fixed-free rod theory",
    )
    st.legend(ax, loc="center right", ncol=2)

    ax = axes[1]
    errors = [
        ("f1_rel_error[-]", "o", "bending mode 1 vs Euler-Bernoulli"),
        ("f2_rel_error[-]", "s", "bending mode 2 vs Euler-Bernoulli"),
        ("f_axial_rel_error[-]", "^", "axial mode 1 vs rod theory"),
    ]
    for slot, (column, marker, label) in enumerate(errors):
        if column not in table:
            continue
        ax.plot(table["num_dofs"], table[column], marker + "-",
                color=st.series_color(slot), label=label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("number of degrees of freedom")
    ax.set_ylabel("relative difference [-]")
    st.title(
        ax, "difference from the analytical references",
        "the axial mode agrees to round-off because plane stress at nu = 0 IS a "
        "uniaxial bar; the bending modes sit a little below Euler-Bernoulli "
        "because the 2-D model includes shear deformation and rotary inertia "
        "that beam theory omits, and the gap grows with mode number",
    )
    st.legend(ax, loc="center right")
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Runtime scaling
# ---------------------------------------------------------------------------
def plot_runtime_scaling(directory: str, path: str, stem: str = "runtime_scaling",
                         dim: int = 2) -> str:
    """Wall-clock cost of each phase against problem size.

    `stem` selects the benchmark files (`runtime_scaling` for the Q4 plate,
    `runtime_scaling_3d` for the Hex8 block) and `dim` the wording.
    """
    table = load_csv(os.path.join(directory, f"{stem}.csv"))
    summary = load_json(os.path.join(directory, f"{stem}.json"))
    exponents = summary.get("scaling_exponent_vs_dofs", {})
    element = "Hex8" if dim == 3 else "Q4"

    fig, axes = st.stacked_panels(2, width=7.2, panel_height=2.6)

    ax = axes[0]
    phases = [
        ("assemble[s]", "assemble K", "assemble"),
        ("factorize[s]", "sparse Cholesky factorisation", "factorize"),
        ("solve[s]", "one back-substitution", "solve"),
        ("objective_gradient[s]", "objective + gradient (one optimiser iteration)",
         "objective_gradient"),
    ]
    for slot, (column, label, key) in enumerate(phases):
        exponent = exponents.get(key)
        suffix = f"  (slope {exponent:.2f})" if exponent is not None else ""
        ax.plot(table["num_dofs"], table[column], "o-", color=st.series_color(slot),
                label=label + suffix)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("wall-clock time [s]")
    st.title(
        ax, f"Runtime scaling of the solver phases ({element})",
        f"{summary.get('compiler', 'unknown compiler')}, "
        f"{summary.get('build_type', 'unknown build')}, "
        f"{summary.get('linear_solver', 'direct solver')}; "
        f"{summary.get('estimator', 'minimum over repeats')}",
    )
    st.legend(ax, loc="upper left")

    ax = axes[1]
    ax.plot(table["num_dofs"], table["stiffness_nonzeros"] / table["num_dofs"],
            "o-", color=st.series_color(2), label="stored entries per free DOF")
    ax.set_xscale("log")
    ax.set_xlabel("number of degrees of freedom")
    ax.set_ylabel("nonzeros / DOF [-]")
    st.title(
        ax, "sparsity of the reduced stiffness matrix",
        f"flat with size, as expected for a fixed-stencil structured {element} mesh",
    )
    st.legend(ax, loc="lower right")

    if dim == 3:
        asymptote = ("a 3-D sparse Cholesky with a fill-reducing ordering is close "
                     "to O(n^2) asymptotically (nested dissection), which is what "
                     "limits the direct solver on solid problems; the multigrid "
                     "solver compared in solver_scaling.png removes that limit")
    else:
        asymptote = ("a 2-D sparse Cholesky with a fill-reducing ordering is close "
                     "to O(n^1.5) asymptotically")
    st.annotate_note(
        fig,
        "Slopes are least-squares fits of log(time) against log(DOFs) over the "
        f"largest three sizes. Assembly is O(n); {asymptote}, and at these sizes "
        "the measured slope also carries cache effects.",
    )
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Linear solvers: sparse Cholesky vs multigrid CG vs Jacobi CG
# ---------------------------------------------------------------------------
#: Colour slot per linear solver. The solver is the entity the colour names, so
#: it keeps its slot in every panel of every solver figure.
SOLVER_SLOTS = {"ldlt": 0, "amg": 1, "jacobi": 2}
SOLVER_NAMES = {
    "ldlt": "sparse Cholesky (LDL^T)",
    "amg": "multigrid-preconditioned CG",
    "jacobi": "Jacobi-preconditioned CG",
}
ELEMENT_MARKERS = {"Quad4": "o", "Tri3": "v", "Hex8": "s", "Tet4": "^", "Tet10": "D"}
ELEMENT_NAMES = {"Quad4": "Q4", "Tri3": "Tri3", "Hex8": "Hex8", "Tet4": "Tet4",
                 "Tet10": "Tet10"}
ELEMENT_DOMAINS = {"Quad4": "Q4 plate", "Hex8": "Hex8 block", "Tet4": "Tet4 block"}

#: sparlab_bench file stem per (element, solver), as scripts/run_scaling.sh
#: writes them.
SCALING_STEMS = {
    ("Quad4", "ldlt"): "runtime_scaling",
    ("Quad4", "amg"): "runtime_scaling_amg_cg",
    ("Hex8", "ldlt"): "runtime_scaling_3d",
    ("Hex8", "amg"): "runtime_scaling_3d_amg_cg",
    ("Hex8", "jacobi"): "runtime_scaling_3d_conjugate_gradient",
    ("Tet4", "ldlt"): "runtime_scaling_3d_tet",
    ("Tet4", "amg"): "runtime_scaling_3d_tet_amg_cg",
}


def load_solver_scaling(directory: str) -> Dict[Tuple[str, str], Tuple[pd.DataFrame, Dict]]:
    """Every solver-scaling table in `directory`, keyed by (element, solver)."""
    tables: Dict[Tuple[str, str], Tuple[pd.DataFrame, Dict]] = {}
    for key, stem in SCALING_STEMS.items():
        csv_path = os.path.join(directory, f"{stem}.csv")
        if not os.path.isfile(csv_path):
            continue
        table = load_csv(csv_path).sort_values("num_dofs").reset_index(drop=True)
        json_path = os.path.join(directory, f"{stem}.json")
        meta = load_json(json_path) if os.path.isfile(json_path) else {}
        tables[key] = (table, meta)
    return tables


def loglog_slope(x, y, tail: int = 3) -> float:
    """Least-squares slope of log(y) on log(x) over the last `tail` points."""
    x = np.asarray(x, dtype=float)[-tail:]
    y = np.asarray(y, dtype=float)[-tail:]
    keep = (x > 0) & (y > 0)
    if keep.sum() < 2:
        return float("nan")
    return float(np.polyfit(np.log(x[keep]), np.log(y[keep]), 1)[0])


def _log_log_axes(ax) -> None:
    """Log scales on both axes, labelling only the decades on an axis that
    spans one or more: minor-tick labels crowd a small panel."""
    from matplotlib.ticker import NullFormatter

    ax.set_xscale("log")
    ax.set_yscale("log")
    for axis, (low, high) in ((ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())):
        if low > 0 and high / low >= 10.0:
            axis.set_minor_formatter(NullFormatter())


def _one_two_five(axis) -> None:
    """Label a log axis at 1, 2 and 5 times each decade, as plain numbers: for
    counts (iterations, entries per DOF) that span one or two decades."""
    from matplotlib.ticker import LogLocator, NullFormatter, StrMethodFormatter

    axis.set_major_locator(LogLocator(base=10.0, subs=(1.0, 2.0, 5.0)))
    axis.set_major_formatter(StrMethodFormatter("{x:g}"))
    axis.set_minor_formatter(NullFormatter())


def _grid_with_legend_row(width: float, height: float, nrows: int, ncols: int):
    """A `nrows` x `ncols` panel grid plus a thin full-width row beneath it for a
    legend shared by every panel (colour means the same thing in all of them).

    A figure-level legend placed "outside" would compete with the footnote for
    the bottom margin; a row of its own is laid out like any other panel."""
    import matplotlib.pyplot as plt

    st.apply_style()
    fig = plt.figure(figsize=(width, height), layout="constrained")
    spec = fig.add_gridspec(nrows + 1, ncols, height_ratios=[1.0] * nrows + [0.16])
    grid = np.array([[fig.add_subplot(spec[r, c]) for c in range(ncols)]
                     for r in range(nrows)])
    legend_ax = fig.add_subplot(spec[nrows, :])
    legend_ax.axis("off")
    return fig, grid, legend_ax


def _solve_time(table: pd.DataFrame) -> np.ndarray:
    """Factorisation (or preconditioner setup) plus one solve, from scratch."""
    return (table["factorize[s]"] + table["solve[s]"]).to_numpy()


def _ratio_label(ratio: float) -> str:
    return f"{ratio:.1f}x" if ratio < 10.0 else f"{ratio:.0f}x"


def plot_solver_scaling(directory: str, path: str) -> str:
    """Time to solve from scratch and solver storage, per solver and element."""
    tables = load_solver_scaling(directory)
    if not tables:
        raise FileNotFoundError(
            f"no runtime_scaling*.csv in {directory}; run scripts/run_scaling.sh")
    elements = [e for e in ("Quad4", "Hex8", "Tet4") if any(k[0] == e for k in tables)]
    meta = next(iter(tables.values()))[1]

    fig, grid, legend_ax = _grid_with_legend_row(7.8, 6.4, 2, len(elements))
    for column, element in enumerate(elements):
        ax_time, ax_store = grid[0, column], grid[1, column]
        marker = ELEMENT_MARKERS[element]
        reference = None
        for solver in ("ldlt", "amg", "jacobi"):
            entry = tables.get((element, solver))
            if entry is None:
                continue
            table = entry[0]
            dofs = table["num_dofs"].to_numpy(dtype=float)
            total = _solve_time(table)
            color = st.series_color(SOLVER_SLOTS[solver])
            ax_time.plot(dofs, total, marker + "-", color=color, markersize=4.2,
                         linewidth=1.6, label=SOLVER_NAMES[solver])
            # The fitted slope as a direct label at the line end.
            ax_time.annotate(f"{loglog_slope(dofs, total):.2f}", (dofs[-1], total[-1]),
                             textcoords="offset points", xytext=(4, 0), fontsize=7.2,
                             color=st.INK_SECONDARY, va="center", ha="left")
            if solver != "jacobi":
                ax_store.plot(dofs, table["solver_nonzeros"].to_numpy(dtype=float) / dofs,
                              marker + "-", color=color, markersize=4.2, linewidth=1.6,
                              label=SOLVER_NAMES[solver])
            if reference is None or len(table) > len(reference):
                reference = table

        # Speed-up at the largest size both the direct and the multigrid solver ran.
        ldlt, amg = tables.get((element, "ldlt")), tables.get((element, "amg"))
        if ldlt is not None and amg is not None:
            both = pd.merge(ldlt[0], amg[0], on="num_dofs", suffixes=("_d", "_m"))
            if not both.empty:
                last = both.iloc[-1]
                t_direct = float(last["factorize[s]_d"] + last["solve[s]_d"])
                t_multigrid = float(last["factorize[s]_m"] + last["solve[s]_m"])
                ax_time.plot([last["num_dofs"]] * 2, [t_multigrid, t_direct], "-",
                             color=st.INK_MUTED, linewidth=0.9)
                ax_time.annotate(_ratio_label(t_direct / t_multigrid),
                                 (last["num_dofs"], np.sqrt(t_direct * t_multigrid)),
                                 textcoords="offset points", xytext=(-4, 0), fontsize=7.6,
                                 color=st.INK_PRIMARY, va="center", ha="right")
        if reference is not None:
            dofs = reference["num_dofs"].to_numpy(dtype=float)
            ax_store.plot(dofs, reference["stiffness_nonzeros"].to_numpy(dtype=float) / dofs,
                          "--", color=st.INK_MUTED, linewidth=1.1,
                          label="the stiffness matrix itself")

        for ax in (ax_time, ax_store):
            _log_log_axes(ax)
        _one_two_five(ax_store.yaxis)
        ax_time.set_title(ELEMENT_DOMAINS[element], loc="left", fontsize=10)
        ax_store.set_xlabel("degrees of freedom")
        if column == 0:
            ax_time.set_ylabel("factorise or set up, then solve [s]")
            ax_store.set_ylabel("stored entries per DOF [-]")

    # One legend per row for the whole figure: the colour is the solver in every
    # panel, and the marker (the element) is named by the panel title.
    from matplotlib.lines import Line2D

    present = [s for s in ("ldlt", "amg", "jacobi") if any(k[1] == s for k in tables)]
    handles = [Line2D([], [], color=st.series_color(SOLVER_SLOTS[s]), linewidth=1.8,
                      label=SOLVER_NAMES[s]) for s in present]
    handles.append(Line2D([], [], color=st.INK_MUTED, linewidth=1.1, linestyle="--",
                          label="entries of K itself (lower row)"))
    legend_ax.legend(handles=handles, loc="center", ncol=2, fontsize=7.8, frameon=False)

    threads = meta.get("openmp_threads")
    st.figure_title(
        fig, "Linear solvers: cost of one solve from scratch, and storage",
        f"{meta.get('compiler', 'unknown compiler')}, {meta.get('build_type', '')}"
        + (f", {threads} OpenMP threads" if threads else "")
        + "; CG to a relative residual of 1e-10; minimum over repeats. Numbers at "
        "the line ends are least-squares slopes of log(time) on log(DOFs) over the "
        "largest three sizes. The vertical bar is labelled with the Cholesky time "
        "over the multigrid time at the largest size both solvers ran.",
    )
    st.annotate_note(
        fig,
        "Upper row: sparse Cholesky is AMD-ordered factorisation plus one "
        "back-substitution; the CG solvers are preconditioner setup plus the "
        "iterations for one load case from a zero initial guess. Lower row: "
        "entries the solver stores beyond K itself - the factor L and D for "
        "Cholesky; the coarse operators, prolongators, restrictions and dense "
        "coarsest factor for multigrid. Jacobi stores only a diagonal. The "
        "benchmark problem is a cantilever with a tip load (sparlab_bench).",
    )
    return st.save_figure(fig, path)


def plot_solver_iterations(directory: str, path: str) -> str:
    """CG iteration counts against problem size for both preconditioners."""
    tables = load_solver_scaling(directory)
    series = [(k, v) for k, v in tables.items() if k[1] in ("amg", "jacobi")]
    if not series:
        raise FileNotFoundError(
            f"no iterative runtime_scaling*.csv in {directory}; run scripts/run_scaling.sh")
    order = {("Quad4", "amg"): 0, ("Hex8", "amg"): 1, ("Tet4", "amg"): 2,
             ("Hex8", "jacobi"): 3}
    series.sort(key=lambda item: order.get(item[0], 9))

    fig, ax = st.figure(7.2, 3.9)
    single_level = False
    for (element, solver), (table, _meta) in series:
        color = st.series_color(SOLVER_SLOTS[solver])
        marker = ELEMENT_MARKERS[element]
        dofs = table["num_dofs"].to_numpy(dtype=float)
        iterations = table["iterations"].to_numpy(dtype=float)
        # A multigrid "hierarchy" of one level is the direct coarse solve: one
        # iteration by construction, drawn hollow and left out of the fit.
        multi = (table["levels"].to_numpy() >= 2) if solver == "amg" else np.ones(len(table), bool)
        slope = loglog_slope(dofs[multi], iterations[multi], tail=len(dofs))
        label = f"{SOLVER_NAMES[solver]}, {ELEMENT_NAMES[element]}"
        if np.isfinite(slope):
            label += f" (slope {slope:.2f})"
        ax.plot(dofs[multi], iterations[multi], marker + "-", color=color, markersize=5,
                linewidth=1.6, label=label)
        if not multi.all():
            single_level = True
            ax.plot(dofs[~multi], iterations[~multi], marker, color=color, markersize=5,
                    markerfacecolor="none", markeredgewidth=1.2)
    _log_log_axes(ax)
    _one_two_five(ax.yaxis)
    ax.set_xlabel("degrees of freedom")
    ax.set_ylabel("CG iterations to 1e-10 [-]")
    st.title(
        ax, "Conjugate-gradient iterations under refinement",
        "slopes are least-squares fits of log(iterations) on log(DOFs) over every "
        "multi-level size; a scalable preconditioner keeps the count flat",
    )
    st.legend(ax, loc="lower right", fontsize=7.8)
    note = ("Colour is the preconditioner, marker the element. Jacobi CG needs "
            "O(sqrt(condition number)) iterations, which grows as 1/h: n^(1/3) in 3-D "
            "and n^(1/2) in 2-D.")
    if single_level:
        note += (" Hollow markers are sizes below the coarse-grid size, where the "
                 "multigrid solver is a direct solve and converges in one iteration.")
    st.annotate_note(fig, note)
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Verification: the 3-D (Hex8) studies
# ---------------------------------------------------------------------------
def plot_sensitivity_check_3d(directory: str, path: str) -> str:
    """Error of the Hex8 topology gradient vs the central-difference step."""
    steps = load_csv(os.path.join(directory, "sensitivity_steps_3d.csv"))
    summary = load_json(os.path.join(directory, "summary.json")).get("sensitivity_3d", {})
    tolerance = summary.get("tolerance", 1e-5)

    fig, ax = st.figure(7.2, 3.6)
    st.require_scatter_series(3, "error measures")
    for slot, (column, label) in enumerate(
        [
            ("max_relative_error[-]", "max relative error over elements"),
            ("rms_relative_error[-]", "RMS relative error over elements"),
            ("directional_relative_error[-]", "relative error of the directional derivative"),
        ]
    ):
        ax.plot(steps["step[-]"], steps[column], "o-", color=st.series_color(slot),
                label=label)
    ax.axhline(tolerance, color=st.INK_MUTED, linewidth=1.0, linestyle="--")
    ax.text(steps["step[-]"].min(), tolerance, f" pass tolerance {tolerance:g}",
            fontsize=7.8, color=st.INK_SECONDARY, va="bottom")
    _log_log_axes(ax)
    ax.invert_xaxis()
    ax.set_xlabel("central-difference step h on the design variable [-]")
    ax.set_ylabel("relative error [-]")
    st.title(
        ax, "Hex8 topology sensitivity: analytical gradient vs central differences",
        f"{summary.get('mesh', 'Hex8 mesh')}, {summary.get('filter', 'density filter')}; "
        f"{summary.get('tested_elements', '?')} elements tested, "
        f"{summary.get('excluded_elements', 0)} excluded at active bounds; best step "
        f"{summary.get('best_step', float('nan')):g} gives "
        f"{summary.get('best_max_relative_error', float('nan')):.2e}",
    )
    st.legend(ax, loc="upper left")
    st.annotate_note(
        fig,
        "The same adjoint gradient as in 2-D, evaluated with the trilinear "
        "element; the error floor is where round-off in the difference of two "
        "compliances overtakes the truncation error of the central difference.",
    )
    return st.save_figure(fig, path)


def plot_mesh_convergence_3d(directory: str, path: str) -> str:
    """Hex8 cantilever tip deflection and discretisation error vs element size."""
    # Coarse to fine, so the last row is the finest mesh the self-convergence
    # error is measured against.
    table = load_csv(os.path.join(directory, "mesh_convergence_3d.csv")).sort_values(
        "h[m]", ascending=False)
    summary = load_json(os.path.join(directory, "summary.json"))
    block = summary.get("mesh_convergence_3d", {})
    order = block.get("observed_convergence_order_tip_deflection")

    fig, axes = st.figure(7.4, 6.4, nrows=2, ncols=1)

    ax = axes[0]
    ax.plot(table["h[m]"], np.abs(table["tip_mean[m]"]), "o-", color=st.series_color(0),
            label="FEM, Hex8, full 2x2x2 integration")
    ax.axhline(abs(float(table["timoshenko[m]"].iloc[0])), color=st.series_color(0),
               linewidth=0.9, linestyle=":", label="Timoshenko (bending + shear)")
    ax.axhline(abs(float(table["euler_bernoulli[m]"].iloc[0])), color=st.INK_MUTED,
               linewidth=1.0, linestyle="--", label="Euler-Bernoulli (bending only)")
    ax.set_xscale("log")
    ax.invert_xaxis()
    ax.set_xlabel("element size h [m]")
    ax.set_ylabel("|tip deflection| [m]")
    st.title(
        ax, "Solid cantilever tip deflection vs mesh size",
        f"{block.get('geometry', '')}; deflection is the mean u_y over the tip face",
    )
    st.legend(ax, loc="lower left")

    ax = axes[1]
    finest = float(table["tip_mean[m]"].iloc[-1])
    error = np.abs(table["tip_mean[m]"] - finest) / abs(finest)
    mask = error > 0
    label = "self-convergence error"
    if order is not None:
        label += f" (observed order {order:.2f})"
    ax.plot(table["h[m]"][mask], error[mask], "o-", color=st.series_color(0), label=label)
    h_ref = np.array([float(table["h[m]"].min()) * 2.0, float(table["h[m]"].max())])
    ax.plot(h_ref, 1e-3 * (h_ref / h_ref[0]) ** 2, "-", color=st.INK_MUTED,
            linewidth=1.0, label=r"reference slope $h^2$")
    _log_log_axes(ax)
    ax.invert_xaxis()
    ax.set_xlabel("element size h [m]")
    ax.set_ylabel("relative error vs the finest mesh [-]")
    st.title(
        ax, "Self-convergence of the discretisation error",
        f"finest mesh within {block.get('finest_mesh_relative_error_vs_timoshenko', float('nan')):.2e} "
        "of Timoshenko",
    )
    st.legend(ax, loc="lower left")
    st.annotate_note(
        fig,
        block.get("note", "") + " Comparing against the finest mesh is verification; "
        "comparing against beam theory is validation, where a finite gap is expected.",
    )
    return st.save_figure(fig, path)


def plot_modal_convergence_3d(directory: str, path: str) -> str:
    """Hex8 cantilever frequencies for the weak and strong bending axes."""
    table = load_csv(os.path.join(directory, "modal_convergence_3d.csv")).sort_values("num_dofs")
    summary = load_json(os.path.join(directory, "summary.json"))
    block = summary.get("modal_validation_3d", {})
    fig, axes = st.stacked_panels(2, width=7.2, panel_height=2.5)

    ax = axes[0]
    series = [
        ("f1_weak_fem[Hz]", "f1_weak_theory[Hz]", "o", "first bending mode, weak axis"),
        ("f1_strong_fem[Hz]", "f1_strong_theory[Hz]", "s", "first bending mode, strong axis"),
    ]
    for slot, (column, reference, marker, label) in enumerate(series):
        ax.plot(table["num_dofs"], table[column], marker + "-", color=st.series_color(slot),
                label=f"{label}, FEM")
        ax.axhline(float(table[reference].iloc[0]), color=st.series_color(slot),
                   linewidth=0.9, linestyle="--", label=f"{label}, Euler-Bernoulli")
    ax.set_xscale("log")
    ax.set_ylabel("frequency [Hz]")
    st.title(
        ax, "Solid cantilever natural frequencies vs mesh size",
        f"{block.get('geometry', '')}; consistent mass matrix; mass conserved to "
        f"{block.get('max_mass_relative_error', float('nan')):.1e} relative",
    )
    st.legend(ax, loc="center right", ncol=2)

    ax = axes[1]
    for slot, (column, marker, label) in enumerate([
        ("f1_weak_rel_error[-]", "o", "weak axis vs Euler-Bernoulli"),
        ("f1_strong_rel_error[-]", "s", "strong axis vs Euler-Bernoulli"),
    ]):
        ax.plot(table["num_dofs"], table[column], marker + "-", color=st.series_color(slot),
                label=label)
    _log_log_axes(ax)
    ax.set_xlabel("number of degrees of freedom")
    ax.set_ylabel("relative difference [-]")
    st.title(ax, "difference from beam theory", block.get("note", ""))
    st.legend(ax, loc="upper right")
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Verification: linear simplices, multigrid, sensitivities through projection
# ---------------------------------------------------------------------------
def plot_mesh_convergence_simplex(directory: str, path: str) -> str:
    """Tri3 / Tet4 cantilever error against the Q4 / Hex8 results on the same beams."""
    simplex = load_csv(os.path.join(directory, "mesh_convergence_simplex.csv"))
    summary = load_json(os.path.join(directory, "summary.json"))
    block = summary.get("mesh_convergence_simplex", {})

    panels = []
    q4_path = os.path.join(directory, "mesh_convergence.csv")
    q4 = load_csv(q4_path) if os.path.isfile(q4_path) else None
    if q4 is not None:
        q4 = q4[np.isclose(q4["poisson[-]"], 0.3)]
    q4_order = (summary.get("mesh_convergence", {})
                .get("observed_convergence_order_tip_deflection", {}).get("0.3"))
    panels.append(("Tri3", q4, "Q4", q4_order,
                   "plane-stress cantilever L = 1 m, h = 0.1 m, t = 10 mm, "
                   "E = 70 GPa, nu = 0.3, 1 kN tip load"))
    hex_path = os.path.join(directory, "mesh_convergence_3d.csv")
    hexa = load_csv(hex_path) if os.path.isfile(hex_path) else None
    hex_order = summary.get("mesh_convergence_3d", {}).get(
        "observed_convergence_order_tip_deflection")
    panels.append(("Tet4", hexa, "Hex8", hex_order,
                   "solid cantilever 1 x 0.1 x 0.05 m, E = 70 GPa, nu = 0, "
                   "1 kN tip load"))

    fig, axes = st.figure(7.4, 6.4, nrows=2, ncols=1)
    for ax, (element, reference, ref_name, ref_order, geometry) in zip(axes, panels):
        sub = simplex[simplex["element"] == element].sort_values("num_dofs")
        if sub.empty:
            ax.axis("off")
            continue
        order = block.get(element, {}).get("observed_convergence_order_tip_deflection")
        label = f"{element} (linear simplex)"
        if order is not None:
            label += f", observed order {order:.2f}"
        ax.plot(sub["num_dofs"], sub["rel_error_timoshenko[-]"], "o-",
                color=st.series_color(0), label=label)
        if reference is not None and not reference.empty:
            reference = reference.sort_values("num_dofs")
            label = f"{ref_name} on the same beam"
            if ref_order is not None:
                label += f", observed order {ref_order:.2f}"
            ax.plot(reference["num_dofs"], reference["rel_error_timoshenko[-]"], "s-",
                    color=st.series_color(1), label=label)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("number of degrees of freedom")
        ax.set_ylabel("relative error vs Timoshenko [-]")
        finest = block.get(element, {}).get("finest_mesh_relative_error_vs_timoshenko")
        st.title(
            ax, f"{element}: tip deflection error vs mesh size",
            geometry + (f"; finest {element} mesh within {100 * finest:.2f} % of "
                        "Timoshenko" if finest is not None else ""),
        )
        st.legend(ax, loc="upper right")
    st.annotate_note(
        fig,
        "The simplex meshes split each structured cell (2 triangles per quadrilateral, "
        "6 tetrahedra per hexahedron), so both curves in a panel describe the same "
        "beam. Observed orders are in element size h, from self-convergence over "
        "the three finest meshes. Constant-strain simplices lock in bending and "
        "approach the beam value from below, more slowly per DOF than the bilinear "
        "and trilinear elements. The error is measured against beam theory, so it "
        "also contains the modelling gap between a beam and plane or solid "
        "elasticity, which is where the Q4 curve levels off.",
    )
    return st.save_figure(fig, path)


def plot_multigrid_study(directory: str, path: str) -> str:
    """Multigrid CG against Cholesky and Jacobi CG as the mesh is refined."""
    table = load_csv(os.path.join(directory, "multigrid_scaling.csv"))
    block = load_json(os.path.join(directory, "summary.json")).get("multigrid", {})
    elements = [e for e in ("Hex8", "Tet4") if (table["element"] == e).any()]
    if not elements:
        raise ValueError("multigrid_scaling.csv has no Hex8 or Tet4 rows")

    fig, grid, legend_ax = _grid_with_legend_row(7.6, 6.1, 2, len(elements))
    single_level = False
    for column, element in enumerate(elements):
        sub = table[table["element"] == element].sort_values("num_dofs")
        dofs = sub["num_dofs"].to_numpy(dtype=float)
        multi = sub["amg_levels"].to_numpy() >= 2
        marker = ELEMENT_MARKERS[element]
        amg, jacobi = (st.series_color(SOLVER_SLOTS[s]) for s in ("amg", "jacobi"))

        ax = grid[0, column]
        ax.plot(dofs, sub["jacobi_iterations"], marker + "-", color=jacobi, markersize=4.5,
                label=SOLVER_NAMES["jacobi"])
        ax.plot(dofs[multi], sub["amg_iterations"].to_numpy()[multi], marker + "-",
                color=amg, markersize=4.5, label=SOLVER_NAMES["amg"])
        if not multi.all():
            single_level = True
            ax.plot(dofs[~multi], sub["amg_iterations"].to_numpy()[~multi], marker,
                    color=amg, markersize=4.5, markerfacecolor="none", markeredgewidth=1.2)
        for x, y, levels in zip(dofs[multi], sub["amg_iterations"].to_numpy()[multi],
                                sub["amg_levels"].to_numpy()[multi]):
            ax.annotate(f"{int(levels)} levels", (x, y), textcoords="offset points",
                        xytext=(0, -11), ha="center", fontsize=6.8, color=st.INK_SECONDARY)
        _log_log_axes(ax)
        _one_two_five(ax.yaxis)
        ax.set_title(f"{element}: CG iterations to 1e-10", loc="left", fontsize=10)
        if column == 0:
            ax.set_ylabel("iterations [-]")

        ax = grid[1, column]
        for solver, key in (("ldlt", "ldlt_seconds[s]"), ("amg", "amg_seconds[s]"),
                            ("jacobi", "jacobi_seconds[s]")):
            ax.plot(dofs, sub[key], marker + "-", markersize=4.5,
                    color=st.series_color(SOLVER_SLOTS[solver]), label=SOLVER_NAMES[solver])
        _log_log_axes(ax)
        ax.set_xlabel("degrees of freedom")
        ax.set_title(f"{element}: wall-clock time of one solve", loc="left", fontsize=10)
        if column == 0:
            ax.set_ylabel("seconds [s]")
    from matplotlib.lines import Line2D

    legend_ax.legend(
        handles=[Line2D([], [], color=st.series_color(SOLVER_SLOTS[s]), linewidth=1.8,
                        label=SOLVER_NAMES[s]) for s in ("ldlt", "amg", "jacobi")],
        loc="center", ncol=3, fontsize=7.8, frameon=False)

    subtitle = ("cantilever block 2 x 1 x 0.5 m (nx x nx/2 x nx/4 cells), E = 70 GPa, "
                "nu = 0.3, tip load")
    difference = block.get("max_relative_difference_vs_ldlt")
    if difference is None:
        difference = float(table["relative_difference_vs_ldlt[-]"].max())
    subtitle += ("; largest relative displacement difference from the Cholesky "
                 f"solution {difference:.1e}")
    growth = block.get("max_iteration_growth_coarsest_to_finest")
    if growth is not None:
        subtitle += ("; iteration growth from the coarsest to the finest multi-level "
                     f"mesh {growth:.2f}x")
    st.figure_title(
        fig, "Multigrid CG: agreement with Cholesky and iterations under refinement",
        subtitle,
    )
    note = ("Times include the preconditioner setup or the factorisation. "
            "The number under each multigrid point is the depth of its hierarchy.")
    if single_level:
        note += (" Hollow markers are meshes below the coarse-grid size, where the "
                 "multigrid solver is a direct solve and converges in one iteration.")
    st.annotate_note(fig, note)
    return st.save_figure(fig, path)


def plot_sensitivity_projection(directory: str, path: str) -> str:
    """Gradient error through the Heaviside projection vs the finite-difference step."""
    table = load_csv(os.path.join(directory, "sensitivity_projection.csv"))
    summary = load_json(os.path.join(directory, "summary.json"))
    tolerance = next((o["tolerance"] for o in summary.get("outcomes", [])
                      if "Heaviside" in o.get("study", "")), None)
    elements = [e for e in ("Quad4", "Tet4") if (table["element"] == e).any()]
    betas = sorted(table["beta[-]"].unique())
    # beta is ordinal, so a sequential ramp; its light end stops short of yellow
    # so the thin lines keep their contrast on the page.
    colors = st.sequence_colors(len(betas), lo=0.08, hi=0.68)

    fig, axes, legend_ax = _grid_with_legend_row(7.4, 6.9, len(elements), 1)
    for row, element in enumerate(elements):
        ax = axes[row, 0]
        for color, beta in zip(colors, betas):
            sub = table[(table["element"] == element)
                        & np.isclose(table["beta[-]"], beta)].sort_values("step[-]")
            ax.plot(sub["step[-]"], sub["max_scaled_error[-]"], "o-", color=color)
            ax.plot(sub["step[-]"], sub["directional_relative_error[-]"], "s--",
                    color=color, markersize=4, linewidth=1.3)
        if tolerance:
            ax.axhline(tolerance, color=st.INK_MUTED, linewidth=1.0, linestyle=":")
            ax.text(float(table["step[-]"].min()), tolerance, f"pass tolerance {tolerance:g}",
                    fontsize=7.6, color=st.INK_SECONDARY, va="bottom", ha="right")
        _log_log_axes(ax)
        ax.invert_xaxis()
        ax.set_xlabel("central-difference step h on the design variable [-]")
        ax.set_ylabel("relative error [-]")
        tested = int(table[table["element"] == element]["num_tested"].iloc[0])
        st.title(
            ax, f"{ELEMENT_NAMES[element]}: compliance gradient through filter and projection",
            f"{tested} design variables tested at eta = 0.5 on a smoothly varying "
            "design; the chain rule runs through the Heaviside projection and the "
            "density filter",
        )
    from matplotlib.lines import Line2D

    handles = [Line2D([], [], color=color, linewidth=1.9, label=f"beta = {beta:g}")
               for color, beta in zip(colors, betas)]
    handles += [
        Line2D([], [], color=st.INK_SECONDARY, marker="o", linewidth=1.9,
               label="largest entry error"),
        Line2D([], [], color=st.INK_SECONDARY, marker="s", markersize=4, linewidth=1.3,
               linestyle="--", label="directional derivative"),
    ]
    legend_ax.legend(handles=handles, loc="center", ncol=len(handles), fontsize=7.6,
                     frameon=False, handlelength=2.2, columnspacing=1.2)
    st.annotate_note(
        fig,
        "Each entry is judged against max(|analytical|, |finite difference|, "
        "1e-3 ||gradient||_inf): at large beta the projection's derivative "
        "vanishes away from the threshold, and there the comparison is with the "
        "round-off floor of a central difference rather than with a relative "
        "error of a number near zero. The study passes on the best step of each "
        "beta, taking the worse of the two measures.",
    )
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Cross-validation against independent codes
# ---------------------------------------------------------------------------
def plot_cross_validation(directory: str, path: str) -> str:
    """Largest relative nodal-displacement difference per code and load case:
    the linear solutions and, where a deck has one, the large-deflection or
    the elastoplastic state."""
    summary = load_json(os.path.join(directory, "summary.json"))
    codes = summary.get("codes", {})
    tolerances = summary.get("tolerances", {})
    rows = []
    for case in summary.get("cases", []):
        for load_case in case.get("load_cases", []):
            rows.append((case["case"], case.get("element_type", ""),
                         load_case["load_case"], load_case["codes"]))
    if not rows:
        raise ValueError("cross-validation summary lists no comparisons")

    # Colour follows the code - the slot does not depend on which codes are
    # present, so scikit-fem keeps its colour on a runner without CalculiX -
    # and the marker the comparison: a circle for the linear solution, a
    # diamond, set just below, for the large-deflection state, a triangle,
    # set just above, for the elastoplastic one.
    code_names = ["scikit-fem", "calculix"]
    st.require_scatter_series(len(code_names), "codes")
    nonlinear_of = {"scikit-fem": "scikit-fem non-linear", "calculix": "calculix NLGEOM"}
    plastic_of = {"scikit-fem": "scikit-fem J2", "calculix": "calculix *PLASTIC"}
    present = [code for code in code_names
               if any(code in r[3] or nonlinear_of[code] in r[3] or plastic_of[code] in r[3]
                      for r in rows)]
    # The legend gets a row of its own beneath the panel: beside it, it took a
    # third of the width from a log axis spanning ten decades. The row has a
    # fixed height: a share of a panel this tall would leave it mostly empty.
    import matplotlib.pyplot as plt

    st.apply_style()
    panel_height, legend_height = 0.36 * len(rows), 1.2
    fig = plt.figure(figsize=(7.4, panel_height + legend_height + 2.4), layout="constrained")
    spec = fig.add_gridspec(2, 1, height_ratios=[panel_height, legend_height])
    ax = fig.add_subplot(spec[0, 0])
    legend_ax = fig.add_subplot(spec[1, 0])
    legend_ax.axis("off")
    y = np.arange(len(rows))[::-1]
    floors = set()
    informational = 0
    # The round-off scale of each linear system, kappa_1 eps: two
    # backward-stable solutions of it can differ by up to about that much.
    round_off = [(results["scikit-fem"]["round_off_scale"], position)
                 for position, (_c, _e, _l, results) in zip(y, rows)
                 if "round_off_scale" in results.get("scikit-fem", {})]
    if round_off:
        ax.plot([r for r, _p in round_off], [p for _r, p in round_off], "|",
                color=st.INK_MUTED, markersize=11, markeredgewidth=1.4,
                label="round-off scale of the linear system, kappa_1(K) eps")
    for slot, code in enumerate(code_names):
        if code not in present:
            continue
        colour = st.series_color(slot)
        version = codes.get(code, {}).get("version", "")
        for key, marker, offset, kind in ((code, "o", 0.0, "linear"),
                                          (nonlinear_of[code], "D", -0.22, "large deflection"),
                                          (plastic_of[code], "^", 0.22, "elastoplastic")):
            judged_x, judged_y, info_x, info_y, own_x, own_y = [], [], [], [], [], []
            for position, (_c, _e, _l, results) in zip(y, rows):
                entry = results.get(key)
                if entry is None:
                    continue
                floor = entry.get("frd_rounding_floor_rel")
                if floor:
                    floors.add(float(floor))
                # passed is None for a comparison between two idealisations
                # (a plane element CalculiX expands through the thickness at
                # nu != 0, CalculiX's own thermal model at finite strain):
                # recorded, not judged.
                if entry.get("passed") is None:
                    info_x.append(entry["max_rel_diff"])
                    info_y.append(position + offset)
                    continue
                # Where CalculiX's formulation of a load differs from SparLab's,
                # the judged number is CalculiX against scikit-fem solving
                # CalculiX's problem; SparLab's own difference is drawn apart.
                judged_x.append(entry.get("max_rel_diff_judged", entry["max_rel_diff"]))
                judged_y.append(position + offset)
                if "max_rel_diff_judged" in entry:
                    own_x.append(entry["max_rel_diff"])
                    own_y.append(position + offset)
            if not (judged_x or info_x):
                continue
            label = f"{code} {version}".strip() + f", {kind}"
            if key == "scikit-fem non-linear":
                label = "scikit-fem, large deflection (its own total Lagrangian solver)"
            elif key == "calculix NLGEOM":
                label = "calculix, large deflection (*STEP, NLGEOM)"
            elif key == "scikit-fem J2":
                label = "scikit-fem, elastoplastic (its own J2 solver)"
            elif key == "calculix *PLASTIC":
                label = "calculix, elastoplastic (*PLASTIC)"
            size = 7 if marker == "o" else 6
            if judged_x:
                ax.plot(judged_x, judged_y, marker, color=colour, markersize=size, label=label)
            if info_x:
                informational += len(info_x)
                ax.plot(info_x, info_y, marker, color=colour, markersize=size,
                        markerfacecolor="none", markeredgewidth=1.5,
                        label=f"{code}, {kind}: different idealisation, not judged")
            if own_x:
                ax.plot(own_x, own_y, "s", color=colour, markersize=5.5,
                        markerfacecolor="none", markeredgewidth=1.2,
                        label=f"{code}, {kind}: SparLab vs CalculiX's own formulation, "
                              "not judged")
    tol_lines = [
        ("skfem", "scikit-fem tolerance", 0, "--", "scikit-fem"),
        ("calculix_solid", "CalculiX tolerance", 1, "--", "calculix"),
    ]
    for key, label, slot, style, code in tol_lines:
        value = tolerances.get(key)
        if value and code in present:
            ax.axvline(value, color=st.series_color(slot), linewidth=1.0, linestyle=style,
                       label=f"{label} {value:g}")
    for floor in sorted(floors):
        ax.axvline(floor, color=st.INK_MUTED, linewidth=1.0, linestyle=":",
                   label=f".frd six-digit rounding floor {floor:g}")
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{c} ({e}), '{l}'" for c, e, l, _r in rows], fontsize=8.0)
    ax.set_xlabel("max |u_SparLab - u_reference| / max |u_reference| over all nodes [-]")
    ax.set_ylim(-0.7, len(rows) - 0.3)
    judged = [entry["passed"] for *_r, results in rows for entry in results.values()
              if entry.get("passed") is not None]
    verdict = ("every judged comparison within its tolerance" if all(judged)
               else "a comparison FAILED")
    if informational:
        verdict += f"; {informational} informational"
    st.figure_title(
        fig, "Cross-validation: nodal displacements vs independent codes",
        f"{len(rows)} load cases, {len(present)} "
        f"code{'s' if len(present) != 1 else ''}, linear, large-deflection and "
        f"elastoplastic; {verdict}",
    )
    handles, labels = ax.get_legend_handles_labels()
    legend_ax.legend(handles, labels, loc="center", ncol=2, fontsize=7.4, frameon=False,
                     columnspacing=1.6)
    st.annotate_note(
        fig,
        "Same mesh, material, supports and nodal loads in every code. scikit-fem "
        "uses the same elements (bilinear, trilinear, linear and quadratic "
        "simplices), so its differences are solver round-off: each linear one lies "
        "below the round-off scale of its system (the bar, the condition number of "
        "the stiffness matrix times eps), which the slender beams' conditioning "
        "raises. CalculiX C3D8, "
        "C3D4 and C3D10 are the same elements as SparLab's Hex8, Tet4 and Tet10; its "
        "plane elements are a layer of solid elements, which matches plane stress "
        "only for nu = 0, so plane-stress comparisons at nu != 0 are hollow and not "
        "judged. Where CalculiX's own formulation of a load differs (the "
        "element-average temperature of a C3D8, the C3D10's rules for a centrifugal "
        "load and a curved face), the filled point is CalculiX against scikit-fem "
        "solving CalculiX's problem - the judged number - and the square SparLab's "
        "difference to CalculiX. The elastoplastic states are compared at the end of "
        "their load paths: scikit-fem with its own J2 solver (small strain and finite "
        "kinematics), CalculiX "
        "with *PLASTIC (isotropic hardening only - its kinematic hardening does not "
        "reproduce Prager's rule - and under NLGEOM a different finite-strain model, "
        "hollow). CalculiX results are read from the .frd file, which carries six "
        "significant digits, so differences below the dotted floor are its output "
        "rounding.",
    )
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Heaviside projection: grey level and the SIMP-vs-structure gap
# ---------------------------------------------------------------------------
#: (label, run without projection, run with projection). Either may be None.
PROJECTION_RUNS = [
    ("MBB beam (2-D, Q4)", "mbb_beam", "mbb_beam_projected"),
    ("bracket (3-D, Hex8)", "bracket_3d", "bracket_3d_projected"),
    ("lug bracket (2-D, Tri3, Gmsh)", None, "lug_bracket_2d"),
    ("engine mount (3-D, Tet4, Gmsh)", None, "engine_mount_3d"),
    ("large bracket (3-D, Hex8, 356k DOFs)", None, "bracket_3d_large"),
]


def collect_projection(results_dir: str) -> pd.DataFrame:
    """One row per run in PROJECTION_RUNS that exists under `results_dir`."""
    rows = []
    for label, plain, projected in PROJECTION_RUNS:
        for variant, name in (("without projection", plain), ("with projection", projected)):
            if name is None:
                continue
            path = os.path.join(results_dir, name, "summary.json")
            if not os.path.isfile(path):
                continue
            doc = load_json(path)
            result = doc.get("optimization_result", {})
            solid = doc.get("interpreted_solid_analysis") or {}
            projection = result.get("projection") or {}
            rows.append({
                "problem": label,
                "run": name,
                "variant": variant,
                "final_beta": projection.get("final_beta"),
                "iterations": result.get("iterations"),
                "stop_reason": result.get("stop_reason"),
                "grey_level": result.get("grey_level"),
                "filtered_grey_level": projection.get("filtered_grey_level"),
                "compliance_J": result.get("compliance_J"),
                "interpreted_compliance_J": solid.get("weighted_compliance_J"),
                "compliance_vs_simp_ratio": solid.get("compliance_vs_simp_ratio"),
                "seconds": result.get("total_seconds"),
            })
    if not rows:
        raise FileNotFoundError(
            f"none of the projection runs exist under {results_dir}; run "
            "scripts/run_all_benchmarks.sh")
    return pd.DataFrame(rows)


def plot_projection_summary(results_dir: str, path: str) -> str:
    """Grey level and thresholded-vs-SIMP compliance, with and without projection."""
    table = collect_projection(results_dir)
    problems = [label for label, *_ in PROJECTION_RUNS if (table["problem"] == label).any()]
    fig, axes = st.figure(7.6, 0.42 * len(problems) + 3.0, nrows=1, ncols=2, sharey=True)
    y = {label: position for position, label in enumerate(problems[::-1])}
    styles = {"without projection": (0, "o"), "with projection": (1, "D")}
    for ax, column, xlabel in (
        (axes[0], "grey_level", "grey level [-]"),
        (axes[1], "compliance_vs_simp_ratio", "thresholded / SIMP compliance [-]"),
    ):
        for label in problems:
            sub = table[table["problem"] == label]
            if len(sub) == 2 and sub[column].notna().all():
                ax.plot(sub[column], [y[label]] * 2, "-", color=st.GRID, linewidth=3.0,
                        zorder=1)
        for variant, (slot, marker) in styles.items():
            sub = table[table["variant"] == variant]
            ax.plot(sub[column], [y[p] for p in sub["problem"]], marker,
                    color=st.series_color(slot), markersize=7, zorder=3,
                    label=variant)
            for _, row in sub.iterrows():
                if row[column] is None or row[column] != row[column]:
                    continue
                ax.annotate(f"{row[column]:.3f}", (row[column], y[row["problem"]]),
                            textcoords="offset points", xytext=(0, 7), ha="center",
                            fontsize=7.0, color=st.INK_SECONDARY)
        ax.set_xlabel(xlabel, fontsize=8.6)
        ax.set_ylim(-0.7, len(problems) - 0.2)
    axes[1].axvline(1.0, color=st.INK_MUTED, linewidth=0.9, linestyle="--")
    ratios = table["compliance_vs_simp_ratio"].dropna()
    if not ratios.empty:
        axes[1].set_xlim(min(0.9, float(ratios.min()) - 0.05),
                         max(1.05, float(ratios.max()) + 0.03))
    greys = table["grey_level"].dropna()
    axes[0].set_xlim(0.0, max(0.05, float(greys.max()) * 1.15) if not greys.empty else 1.0)
    axes[0].set_yticks(list(y.values()))
    axes[0].set_yticklabels(list(y.keys()), fontsize=8.2)
    st.legend(axes[0], loc="lower right", fontsize=7.8)
    st.figure_title(
        fig, "Heaviside projection: how binary the design is, and what that buys",
        "grey level = 4 mean(rho (1 - rho)) over all elements; the ratio "
        "re-solves the design thresholded at 0.5 as solid material, so 1.0 means "
        "the optimiser's objective is the structure's real compliance",
    )
    st.annotate_note(
        fig,
        "The grey rails join the two runs of one problem; the runs differ only in "
        "the projection block and the settings listed in each deck's header "
        "comment. Problems with one point were run with projection only. Values "
        "are read from each run's summary.json.",
    )
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Design study
# ---------------------------------------------------------------------------
def collect_study(study_dir: str) -> pd.DataFrame:
    """Flatten every `summary.json` under a study directory into one table."""
    rows: List[Dict] = []
    pattern = os.path.join(study_dir, "*", "*", "summary.json")
    for path in sorted(glob.glob(pattern)):
        parts = path.split(os.sep)
        arm, tag = parts[-3], parts[-2]
        doc = load_json(path)
        setup = doc.get("optimization_setup", {})
        result = doc.get("optimization_result", {})
        mesh = doc.get("mesh", {})
        material = doc.get("material", {})
        comparison = doc.get("mass_stiffness_comparison", {})
        interp = doc.get("solid_interpretation", {})
        row = {
            "arm": arm,
            "tag": tag,
            "volume_fraction_target": setup.get("volume_fraction_target"),
            "volume_fraction": result.get("volume_fraction"),
            "compliance_J": result.get("compliance_J"),
            "mass_kg": result.get("mass_kg"),
            "grey_level": result.get("grey_level"),
            "iterations": result.get("iterations"),
            "converged": result.get("converged"),
            "stop_reason": result.get("stop_reason"),
            "seconds": result.get("total_seconds"),
            "num_elements": mesh.get("num_elements"),
            "num_dofs": mesh.get("num_dofs"),
            "penalty": setup.get("simp_penalty"),
            "filter_radius_m": setup.get("filter_radius_m"),
            "filter_support": setup.get("filter_average_support"),
            "youngs_modulus_Pa": material.get("youngs_modulus_Pa"),
            "full_solid_compliance_J": comparison.get("full_solid_compliance_J"),
            "equal_mass_plate_compliance_J":
                comparison.get("equal_mass_uniform_plate_compliance_J"),
            "stiffness_gain_over_equal_mass_plate":
                comparison.get("stiffness_gain_over_equal_mass_plate"),
            "retained_elements": interp.get("elements_retained"),
            "connected_groups": interp.get("connected_groups_above_threshold"),
            "islands_discarded_m3": interp.get("volume_discarded_as_islands_m3"),
        }
        # What the thresholded structure actually delivers, which is not the
        # objective the optimiser minimised.
        solid = doc.get("interpreted_solid_analysis") or {}
        row["interpreted_compliance_J"] = solid.get("weighted_compliance_J")
        row["compliance_vs_simp_ratio"] = solid.get("compliance_vs_simp_ratio")
        row["interpreted_mass_kg"] = solid.get("mass_kg")
        row["interpreted_max_von_mises_Pa"] = solid.get("max_von_mises_Pa")
        row["interpreted_analysis_failed"] = bool(solid.get("analysis_failed", False))

        modal = doc.get("modal_optimised_topology")
        if modal and modal.get("frequencies_hz"):
            row["f1_topology_Hz"] = modal["frequencies_hz"][0]
            row["topology_mass_kg"] = modal.get("total_mass_kg")
        solid_modal = doc.get("modal_initial_solid")
        if solid_modal and solid_modal.get("frequencies_hz"):
            row["f1_solid_Hz"] = solid_modal["frequencies_hz"][0]
        # Per-load-case columns. The objective is the weight-normalised mean of
        # these, so both the weight and the case compliance are carried through:
        # compliance_J == sum_l weight_norm_l * compliance_<case_l>_J exactly.
        names = setup.get("load_case_names") or []
        for load_case, value in zip(names, result.get("load_case_compliance_J") or []):
            row[f"compliance_{load_case}_J"] = value
        for load_case, value in zip(names, setup.get("load_case_weights") or []):
            row[f"weight_{load_case}"] = value
        rows.append(row)

    if not rows:
        raise FileNotFoundError(
            f"no summary.json found under {study_dir}/*/*/; run "
            "scripts/run_aerospace_study.sh first"
        )
    return pd.DataFrame(rows)


def _numeric_tag(tag: str) -> float:
    """Extract the varied quantity from a study tag like `vf_0.35` or `p_3.0`.

    The *last* numeric token wins, not the first: a load-weighting tag is
    `w_<down>_<reversal>_<lateral>`, and the quantity the arm varies is the
    last one. A mesh tag's `AxB` token resolves to the element count, which is
    the right ordering key for `n_120x80`.

    Prefer a recorded column where one exists - `order_by` below - and keep this
    for the arms whose varied quantity is not a column of its own.
    """
    value = float("nan")
    for token in tag.replace("-", "_").split("_"):
        try:
            value = float(token)
            continue
        except ValueError:
            pass
        if "x" in token:
            parts = token.split("x")
            try:
                if len(parts) == 2:
                    value = float(parts[0]) * float(parts[1])
            except ValueError:
                pass
    return value


def plot_mass_stiffness_pareto(table: pd.DataFrame, path: str) -> str:
    """Compliance and first frequency against mass for the volume sweep."""
    arm = table[table["arm"] == "volume_fraction"].copy()
    if arm.empty:
        raise ValueError("the study table has no volume_fraction arm")
    arm["vf"] = arm["tag"].map(_numeric_tag)
    arm = arm.sort_values("mass_kg")

    fig, axes = st.stacked_panels(2, width=7.2, panel_height=2.7)

    ax = axes[0]
    st.require_scatter_series(3, "mass-stiffness curves")
    ax.plot(arm["mass_kg"], arm["compliance_J"], "o-", color=st.series_color(0),
            label="SIMP density field (the objective minimised)")
    if "interpreted_compliance_J" in arm and arm["interpreted_compliance_J"].notna().any():
        # The structure a threshold actually produces. It is stiffer than the
        # objective here, because thresholding promotes the filter's grey
        # boundary band to solid material.
        ax.plot(arm["mass_kg"], arm["interpreted_compliance_J"], "^-",
                color=st.series_color(2),
                label="the same design thresholded at 0.5 and re-solved")
    ax.plot(arm["mass_kg"], arm["equal_mass_plate_compliance_J"], "s--",
            color=st.series_color(1),
            label="uniform plate of the same mass")
    for _, row in arm.iterrows():
        ax.annotate(
            f"{row['vf']:.2f}", (row["mass_kg"], row["compliance_J"]),
            textcoords="offset points", xytext=(0, 7), ha="center", fontsize=7.4,
            color=st.INK_SECONDARY,
        )
    ax.set_yscale("log")
    ax.set_ylabel("weighted compliance [J]")
    st.title(
        ax, "Aerospace bracket: mass-stiffness trade",
        "points are labelled with the volume-fraction target; lower is stiffer",
    )
    st.legend(ax, loc="upper right")

    ax = axes[1]
    if "f1_topology_Hz" in arm and arm["f1_topology_Hz"].notna().any():
        ax.plot(arm["mass_kg"], arm["f1_topology_Hz"], "o-",
                color=st.series_color(2),
                label="first natural frequency of the interpreted topology")
        if arm["f1_solid_Hz"].notna().any():
            ax.axhline(float(arm["f1_solid_Hz"].dropna().iloc[0]),
                       color=st.INK_MUTED, linewidth=1.0, linestyle="--",
                       label="full solid domain (also the equal-mass uniform plate)")
        ax.set_ylabel("f1 [Hz]")
    ax.set_xlabel("mass [kg]")
    st.title(
        ax, "vibration metric across the same sweep",
        "in this 2-D idealisation a uniformly thinned plate keeps the full-solid "
        "frequencies exactly, so the dashed line is the equal-mass baseline at "
        "every mass",
    )
    st.legend(ax, loc="lower right")
    st.annotate_note(
        fig,
        "The equal-mass uniform plate is C_solid / volume_fraction, which follows "
        "from compliance scaling as 1/thickness. Both curves come from solver "
        "runs recorded in the study summaries; no fit or extrapolation is used.",
    )
    return st.save_figure(fig, path)


def plot_numerical_settings(table: pd.DataFrame, path: str) -> str:
    """How the numerical settings change compliance and the design."""
    arms = [
        ("penalty", "SIMP penalty p [-]", "penalty"),
        ("filter_radius", "filter radius [cells]", "tag"),
        ("mesh_fixed_r", "elements", "num_elements"),
        ("youngs_modulus", "Young's modulus [Pa]", "youngs_modulus_Pa"),
    ]
    present = [a for a in arms if not table[table["arm"] == a[0]].empty]
    fig, axes = st.figure(7.4, 2.3 * len(present), nrows=len(present), ncols=1)
    axes = np.atleast_1d(axes).ravel()

    for ax, (arm, xlabel, key) in zip(axes, present):
        sub = table[table["arm"] == arm].copy()
        sub["x"] = sub[key] if key != "tag" else sub["tag"].map(_numeric_tag)
        sub = sub.sort_values("x")
        ax.plot(sub["x"], sub["compliance_J"], "o-", color=st.series_color(0),
                label="weighted compliance [J]")
        ax.set_ylabel("compliance [J]")
        ax.set_xlabel(xlabel)
        if arm in ("mesh_fixed_r", "youngs_modulus"):
            ax.set_xscale("log")
        twin_label = "grey level [-]"
        # No second y-axis: the grey level is shown as a direct label per point
        # so the two quantities never share a scale.
        for _, row in sub.iterrows():
            ax.annotate(
                f"grey {row['grey_level']:.2f}", (row["x"], row["compliance_J"]),
                textcoords="offset points", xytext=(0, 8), ha="center",
                fontsize=7.0, color=st.INK_SECONDARY,
            )
        st.title(ax, f"effect of {xlabel.split('[')[0].strip()}",
                 f"point labels give the {twin_label} of the resulting design")
        st.legend(ax, loc="best")

    st.annotate_note(
        fig,
        "Each arm changes one setting against the shared baseline deck. "
        "Compliance across the Young's-modulus arm scales as 1/E exactly, which "
        "is the expected linear-elastic behaviour and a useful consistency check "
        "rather than a design result.",
    )
    return st.save_figure(fig, path)


def plot_mesh_dependence(table: pd.DataFrame, path: str) -> str:
    """Mesh refinement with the filter radius fixed in metres vs in cells."""
    fixed_r = table[table["arm"] == "mesh_fixed_r"].copy()
    fixed_cells = table[table["arm"] == "mesh_fixed_cells"].copy()
    if fixed_r.empty or fixed_cells.empty:
        raise ValueError("the study table is missing a mesh-refinement arm")

    fig, axes = st.stacked_panels(2, width=7.2, panel_height=2.6)

    ax = axes[0]
    for slot, (sub, label) in enumerate(
        [
            (fixed_r, "filter radius fixed at 3.75 mm (length scale held)"),
            (fixed_cells, "filter radius fixed at 1.5 cells (length scale shrinks)"),
        ]
    ):
        sub = sub.sort_values("num_elements")
        ax.plot(sub["num_elements"], sub["compliance_J"], "o-",
                color=st.series_color(slot), label=label)
    ax.set_xscale("log")
    ax.set_ylabel("compliance [J]")
    st.title(
        ax, "Mesh dependence of the optimised design",
        "holding the filter radius in physical units makes the optimum mesh "
        "convergent; holding it in cells does not",
    )
    st.legend(ax, loc="upper right")

    ax = axes[1]
    for slot, (sub, label) in enumerate(
        [(fixed_r, "radius fixed in metres"), (fixed_cells, "radius fixed in cells")]
    ):
        sub = sub.sort_values("num_elements")
        ax.plot(sub["num_elements"], sub["grey_level"], "o-",
                color=st.series_color(slot), label=label)
    ax.set_xscale("log")
    ax.set_xlabel("number of elements")
    ax.set_ylabel("grey level [-]")
    st.title(
        ax, "how binary the design is under refinement",
        "a radius fixed in cells produces ever-finer members, which is the "
        "classical mesh-dependence pathology filtering is meant to remove",
    )
    st.legend(ax, loc="best")
    return st.save_figure(fig, path)


def plot_load_weighting(table: pd.DataFrame, path: str) -> str:
    """Effect of the lateral load-case weight, objective and per-case.

    The sweep holds the down and reversal weights at (1, 0.5) and varies the
    lateral weight, so the x axis is the recorded weight of the lateral case
    rather than anything parsed out of a directory name. The one point that
    drops the reversal case entirely is a different experiment and is excluded
    from the curves; it is reported in the table and the topology montage.
    """
    arm = table[table["arm"] == "load_weighting"].copy()
    if arm.empty:
        raise ValueError("the study table has no load_weighting arm")

    required = ["weight_lateral", "weight_up_reversal", "compliance_down_limit_J",
                "compliance_up_reversal_J", "compliance_lateral_J"]
    missing = [c for c in required if c not in arm.columns]
    if missing:
        raise KeyError(
            "the load_weighting summaries carry no per-load-case data "
            f"(missing {missing}); re-run the study with a build that records "
            "optimization_setup.load_case_names"
        )

    sweep = arm[np.isclose(arm["weight_up_reversal"], 0.5)].sort_values("weight_lateral")
    if sweep.empty:
        raise ValueError("no load_weighting point holds the reversal weight at 0.5")

    fig, axes = st.stacked_panels(2, width=7.2, panel_height=2.6)

    ax = axes[0]
    ax.plot(sweep["weight_lateral"], sweep["compliance_J"], "o-",
            color=st.series_color(0),
            label="objective: weight-normalised mean over the three cases")
    last = len(sweep) - 1
    for position, (_, row) in enumerate(sweep.iterrows()):
        # End labels are aligned inwards so they cannot fall off the axes.
        align = "left" if position == 0 else "right" if position == last else "center"
        ax.annotate(
            f"grey {row['grey_level']:.3f}",
            (row["weight_lateral"], row["compliance_J"]),
            textcoords="offset points", xytext=(0, -13), ha=align, fontsize=7.0,
            color=st.INK_SECONDARY,
        )
    ax.set_ylabel("objective [J]")
    st.title(
        ax, "Aerospace bracket: effect of the lateral load-case weight",
        "weights are (down 1.0, reversal 0.5, lateral w). The objective is a "
        "mean over a weight set that changes along the axis, so it is not a "
        "like-for-like stiffness measure",
    )
    st.legend(ax, loc="upper right")

    # Each case relative to its value on the w = 0 design. The three absolute
    # levels differ by a factor of 19, so on a shared axis - log or not - the
    # changes that matter are invisible; the ratio is what the arm is about,
    # and the absolute values are in the table beside this figure.
    ax = axes[1]
    st.require_scatter_series(3, "load cases")
    reference = sweep[np.isclose(sweep["weight_lateral"], 0.0)]
    # The reversal case is exactly -0.5 times the down case, so their ratio
    # curves coincide; a dashed overlay keeps both readable instead of one
    # hiding the other.
    for slot, (column, label, style) in enumerate(
        [
            ("compliance_down_limit_J", "down limit case (weight 1.0)", "o-"),
            ("compliance_up_reversal_J",
             "up reversal case (weight 0.5) - coincides exactly", "o--"),
            ("compliance_lateral_J", "lateral case (weight w)", "o-"),
        ]
    ):
        base = float(reference[column].iloc[0]) if not reference.empty else 1.0
        ax.plot(sweep["weight_lateral"], sweep[column] / base, style,
                color=st.series_color(slot), label=label,
                markersize=7.5 if slot == 0 else 4.5,
                markerfacecolor="none" if slot == 1 else None,
                markeredgewidth=1.4 if slot == 1 else 0.0)
    ax.axhline(1.0, color=st.INK_MUTED, linewidth=0.9, linestyle="--")
    ax.set_xlabel("lateral load-case weight w [-]")
    ax.set_ylabel("compliance / its value at w = 0 [-]")
    st.title(
        ax, "what the weighting actually buys",
        "each series is one load case solved on the design that came out of "
        "that run. The lateral case gains 39 %; the two vertical cases do not "
        "pay for it, which is the local-optimum result discussed in the text",
    )
    st.legend(ax, loc="lower left")

    other = arm[~np.isclose(arm["weight_up_reversal"], 0.5)]
    note = (
        "Weights are normalised to sum to one inside the objective, so the "
        "objective is a mean and is not comparable in magnitude across "
        "different weight sets. The lower panel is each case relative to its "
        "own value on the w = 0 design; absolute compliances are tabulated in "
        "docs/aerospace_study.md."
    )
    if not other.empty:
        tags = ", ".join(
            f"{row['tag']} (objective {row['compliance_J']:.4g} J)"
            for _, row in other.iterrows()
        )
        note += (
            " Excluded from the curves because it changes a second weight: "
            f"{tags}."
        )
    st.annotate_note(fig, note)
    return st.save_figure(fig, path)


def plot_study_topologies(study_dir: str, table: pd.DataFrame, arm: str, path: str,
                          title_text: str, subtitle: str,
                          order_by: Optional[str] = None) -> str:
    """Small multiples of the optimised density for one study arm.

    `order_by` names the recorded column the panels are ordered by - the
    quantity the arm actually varies. Reading the order out of the directory
    names is the fallback.
    """
    from . import fields as fld
    from .loaders import load_mesh

    sub = table[table["arm"] == arm].copy()
    if sub.empty:
        raise ValueError(f"the study table has no {arm} arm")
    if order_by is not None and order_by in sub.columns and sub[order_by].notna().all():
        sub["x"] = sub[order_by]
    else:
        sub["x"] = sub["tag"].map(_numeric_tag)
    sub = sub.sort_values(["x", "tag"])

    count = len(sub)
    columns = min(4, count)
    rows = (count + columns - 1) // columns
    fig, axes = st.figure(7.6, 1.6 * rows + 0.4, nrows=rows, ncols=columns)
    axes = np.atleast_1d(axes).ravel()

    for slot, (_, row) in enumerate(sub.iterrows()):
        directory = os.path.join(study_dir, arm, row["tag"])
        mesh = load_mesh(directory)
        density = load_csv(os.path.join(directory, "density_final.csv"))
        values = density["physical_density[-]"].to_numpy()
        ax = axes[slot]
        shape = fld.structured_grid_shape(mesh)
        extent = mesh.extent
        if shape is not None:
            ax.imshow(fld.density_image(values, shape), cmap=st.DENSITY_CMAP_SURFACE,
                      vmin=0.0, vmax=1.0, origin="lower",
                      extent=(extent[0], extent[1], extent[2], extent[3]),
                      interpolation="nearest", aspect="equal")
        else:
            fld.element_collection(ax, mesh, values, cmap=st.DENSITY_CMAP_SURFACE,
                                   vmin=0.0, vmax=1.0)
            fld.set_domain_limits(ax, mesh)
        fld.bare_axes(ax)
        ax.set_title(
            f"{row['tag']}\nc = {st.format_si(row['compliance_J'])} J",
            loc="left", fontsize=8.4,
        )
    for slot in range(count, len(axes)):
        axes[slot].axis("off")

    fig.suptitle(title_text, x=0.01, ha="left", fontsize=11, fontweight="bold",
                 color=st.INK_PRIMARY)
    st.annotate_note(fig, subtitle)
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Quadratic tetrahedra and linear buckling
# ---------------------------------------------------------------------------
#: Fixed categorical slot per element for the Tet10 and buckling figures, so
#: an element keeps its colour from one figure to the next.
ELEMENT_SLOTS = {"Hex8": 0, "Tet4": 1, "Tet10": 2, "Quad4": 3}


def plot_tet10_convergence(directory: str, path: str) -> str:
    """Cantilever tip error of Hex8, Tet4 and Tet10 on the same grids."""
    table = load_csv(os.path.join(directory, "mesh_convergence_tet10.csv"))
    summary = load_json(os.path.join(directory, "summary.json")).get(
        "mesh_convergence_tet10", {})
    fig, ax = st.figure(7.4, 4.4)
    for element in ("Hex8", "Tet4", "Tet10"):
        sub = table[table["element"] == element].sort_values("num_dofs")
        if sub.empty:
            continue
        ax.plot(sub["num_dofs"], sub["rel_error_timoshenko[-]"],
                ELEMENT_MARKERS[element] + "-", color=st.series_color(ELEMENT_SLOTS[element]),
                label=element)
        last = sub.iloc[-1]
        ax.annotate(f"{100.0 * float(last['rel_error_timoshenko[-]']):.2g} %",
                    (float(last["num_dofs"]), float(last["rel_error_timoshenko[-]"])),
                    xytext=(6, 0), textcoords="offset points", va="center", fontsize=8,
                    color=st.INK_SECONDARY)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("number of degrees of freedom")
    ax.set_ylabel("relative error vs Timoshenko [-]")
    st.title(ax, "Tip deflection error: linear against quadratic elements",
             "solid cantilever 1 x 0.1 x 0.05 m, nu = 0, clamped root, tip traction with "
             "a 1 kN resultant; grids 10 x 2 x 1, 20 x 4 x 2, 40 x 8 x 4, each cell one "
             "Hex8 or six Kuhn tetrahedra")
    st.legend(ax, loc="lower left")
    st.annotate_note(
        fig,
        "The quadratic tetrahedron contains the linear-stress field of pure bending, so "
        "it does not lock: on the coarsest grid it is already within "
        f"{100.0 * float(summary.get('Tet10', {}).get('records', [{}])[0].get('relative_error_vs_timoshenko', float('nan'))):.2g} % "
        "of beam theory. The error is measured against Timoshenko beam theory, so it "
        "also holds the modelling gap between a beam and solid elasticity.",
    )
    return st.save_figure(fig, path)


def plot_buckling_verification(directory: str, path: str) -> str:
    """Critical load of a clamped column against Euler-Engesser."""
    table = load_csv(os.path.join(directory, "buckling_euler.csv"))
    fig, ax = st.figure(7.4, 4.4)
    for element in ("Quad4", "Hex8", "Tet4", "Tet10"):
        sub = table[table["element"] == element].sort_values("num_dofs")
        if sub.empty:
            continue
        label = "Q4 (plane stress)" if element == "Quad4" else element
        ax.plot(sub["num_dofs"], sub["rel_error_engesser[-]"],
                ELEMENT_MARKERS[element] + "-", color=st.series_color(ELEMENT_SLOTS[element]),
                label=label)
        last = sub.iloc[-1]
        ax.annotate(f"{100.0 * float(last['rel_error_engesser[-]']):.2g} %",
                    (float(last["num_dofs"]), float(last["rel_error_engesser[-]"])),
                    xytext=(6, 0), textcoords="offset points", va="center", fontsize=8,
                    color=st.INK_SECONDARY)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("number of degrees of freedom")
    ax.set_ylabel("(lambda_1 - P_Engesser) / P_Engesser [-]")
    # The Tet10 meshes halve h exactly, so Richardson extrapolation of the
    # three errors separates the discretisation error from the modelling gap.
    gap_text = ""
    tet10 = table[table["element"] == "Tet10"].sort_values("num_dofs")
    if len(tet10) >= 3:
        e1, e2, e3 = tet10["rel_error_engesser[-]"].to_numpy()[-3:]
        if (e1 - e2) * (e2 - e3) > 0 and e2 != e3:
            ratio = (e1 - e2) / (e2 - e3)
            if ratio > 1.0:
                order = np.log2(ratio)
                limit = e3 - (e2 - e3) / (ratio - 1.0)
                gap_text = (f" Richardson extrapolation of the Tet10 errors (order "
                            f"{order:.2f} in h) puts their limit at {100.0 * limit:.2f} %: the "
                            "gap between the solid model and the beam formula, not a "
                            "discretisation error.")
    st.title(ax, "Linear buckling of a clamped column: first load factor vs mesh",
             "1 m steel column, 50 mm square section (3-D) or 50 x 20 mm (plane stress), "
             "unit axial tip traction, so lambda_1 is the critical load in newtons; "
             "reference Euler P = pi^2 E I / (4 L^2) with Engesser's shear correction")
    st.legend(ax, loc="lower left")
    st.annotate_note(
        fig,
        "Every point lies above the reference: a displacement-based element is too "
        "stiff, and its buckling load converges from above. The linear Hex8 and Tet4 "
        "lock in bending and converge slowly; the Q4 and the quadratic Tet10 reach "
        "the 1 % tolerance of sparlab_verify's buckling-euler study." + gap_text,
    )
    return st.save_figure(fig, path)


def plot_tet10_part_study(directory: str, path: str) -> str:
    """Compliance of the engine mount against DOFs, per element and geometry."""
    table = load_csv(os.path.join(directory, "tet10_part_study.csv"))
    summary = load_json(os.path.join(directory, "summary.json"))
    variants = [("Tet4", "^", 1), ("Tet10 (straight)", "d", 2), ("Tet10 (curved)", "D", 3)]
    fig, axes = st.figure(7.4, 6.6, nrows=2, ncols=1, sharex=True)
    for ax, case, label in zip(axes, ("vertical", "lateral"),
                               ("vertical pin load, 12 kN", "lateral pin load, 5 kN")):
        key = f"compliance_{case}_J"
        for variant, marker, slot in variants:
            sub = table[table["variant"] == variant].sort_values("num_dofs")
            if sub.empty:
                continue
            ax.plot(sub["num_dofs"], sub[key], marker + "-", color=st.series_color(slot - 1),
                    label=variant)
            for _, row in sub.iterrows():
                ax.annotate(f"{row['size_mm']:g}", (row["num_dofs"], row[key]),
                            xytext=(0, 6 if variant == "Tet10 (curved)" else -11),
                            textcoords="offset points", ha="center", fontsize=7,
                            color=st.INK_MUTED)
        ref = summary["references"][case]
        ax.axhline(ref["finest_tet10_curved_J"], color=st.INK_MUTED, linewidth=0.9,
                   linestyle="--")
        ax.set_xscale("log")
        ax.set_ylabel("compliance f.u [J]")
        st.title(ax, f"Engine mount, {label}",
                 f"dashed: the finest curved Tet10 run ({ref['finest_size_mm']:g} mm), itself "
                 "a lower bound; numbers at the points are the element size in mm")
        st.legend(ax, loc="lower right")
    axes[-1].set_xlabel("number of degrees of freedom")
    st.annotate_note(
        fig,
        "Each mesh comes from the same CAD model through Gmsh with the curvature "
        "refinement scaled with the element size. A displacement-based solution is "
        "too stiff, so its compliance lies below the exact value and rises towards it. "
        "The straight-sided Tet10 meshes are the linear meshes with an edge node at "
        "every edge midpoint: the same faceted holes as the Tet4 mesh, a slightly "
        "different part from the curved one.",
    )
    return st.save_figure(fig, path)


def plot_buckling_cross_validation(directory: str, path: str) -> str:
    """Buckling load factors against scikit-fem and CalculiX *BUCKLE."""
    summary = load_json(os.path.join(directory, "summary.json"))
    codes = [("scikit-fem buckling", "scikit-fem (same K_G, dense eigensolve)"),
             ("calculix buckling", "CalculiX *BUCKLE")]
    rows = []
    for case in summary.get("cases", []):
        for load_case in case.get("load_cases", []):
            if any(code in load_case["codes"] for code, _ in codes):
                rows.append((case["case"], case.get("element_type", ""),
                             load_case["load_case"], load_case["codes"]))
    if not rows:
        raise ValueError("cross-validation summary has no buckling comparisons")
    fig, ax = st.figure(7.4, 0.5 * len(rows) + 3.3)
    y = np.arange(len(rows))[::-1]
    for slot, (code, label) in enumerate(codes):
        xs, ys = [], []
        for position, (*_rest, results) in zip(y, rows):
            entry = results.get(code)
            if entry is not None:
                xs.append(entry["max_rel_diff"])
                ys.append(position)
        if xs:
            ax.plot(xs, ys, "o", color=st.series_color(slot), markersize=7, label=label)
    tolerances = summary.get("tolerances", {})
    for slot, key, label in ((0, "skfem_buckling", "scikit-fem tolerance"),
                             (1, "calculix_buckling", "CalculiX tolerance")):
        value = tolerances.get(key)
        if value:
            ax.axvline(value, color=st.series_color(slot), linewidth=1.0, linestyle="--",
                       label=f"{label} {value:g}")
    ax.set_xscale("log")
    ax.set_xlim(1.0e-12, 1.0e-2)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{c} ({e}), '{l}'" for c, e, l, _r in rows], fontsize=8.0)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_xlabel("max over the reported modes of |lambda_ref - lambda| / lambda [-]")
    st.title(ax, "Cross-validation: linear buckling load factors",
             "scikit-fem assembles the geometric stiffness from its own static solution at "
             "the same quadrature points (the same discrete problem); CalculiX runs "
             "*BUCKLE on the exported deck with its own stress-stiffness evaluation")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2, fontsize=8,
              frameon=False)
    st.annotate_note(
        fig,
        "Where an element's stress varies inside it (C3D8, C3D10), CalculiX's factors "
        "lie a few 1e-5 above SparLab's and scikit-fem's, which agree to about 1e-9; "
        "for C3D4, whose stress is constant in an element, the gap is several times "
        "smaller. Which detail of CalculiX's stress stiffness accounts for it has not "
        "been identified; the tolerance records the measured size, it does not "
        "explain it.",
    )
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Verification: geometrically non-linear statics
# ---------------------------------------------------------------------------

def plot_arch_snap_through(directory: str, path: str) -> str:
    """The shallow arch's force-deflection path, its unstable stretch, the
    limit load and where load control stopped."""
    table = load_csv(os.path.join(directory, "arch_snap_through.csv"))
    block = load_json(os.path.join(directory, "summary.json")).get("arch_snap_through")
    if not block:
        raise ValueError("summary.json has no arch_snap_through block")
    mm = 1.0e3
    deflection = table["crown_deflection[m]"].to_numpy() * mm
    force = table["force_arc_length[N]"].to_numpy()
    control = table["force_displacement_control[N]"].to_numpy()
    unstable = table["negative_pivots_arc_length"].to_numpy() > 0
    fig, ax = st.figure(7.4, 4.6)
    # The unstable stretch: states whose tangent has a negative pivot.
    if unstable.any():
        first = int(np.argmax(unstable))
        last = len(unstable) - 1 - int(np.argmax(unstable[::-1]))
        lo = 0.5 * (deflection[first - 1] + deflection[first]) if first > 0 else deflection[0]
        hi = (0.5 * (deflection[last] + deflection[last + 1]) if last + 1 < len(deflection)
              else deflection[-1])
        ax.axvspan(lo, hi, color=st.GRID, alpha=0.6, linewidth=0.0,
                   label="unstable: the tangent has a negative pivot")
    ax.plot(deflection, force, "-", color=st.series_color(0), linewidth=1.6)
    ax.plot(deflection[~unstable], force[~unstable], "o", color=st.series_color(0),
            markersize=5.5, label="arc length, stable")
    ax.plot(deflection[unstable], force[unstable], "o", markerfacecolor=st.SURFACE,
            markeredgecolor=st.series_color(0), markeredgewidth=1.3, markersize=5.5,
            label="arc length, unstable")
    ax.plot(deflection, control, "x", color=st.series_color(1), markersize=6.0,
            markeredgewidth=1.3, label="displacement control")
    limit = float(block["limit_force_N"])
    limit_x = float(block["limit_deflection_m"]) * mm
    ax.annotate(f"limit {limit:.2f} N", (limit_x, limit), xytext=(-6, 9),
                textcoords="offset points", ha="right", fontsize=8,
                color=st.INK_SECONDARY)
    # Load control's stop, and the jump it refuses: the far branch at the same
    # force, on the rising stretch after the lowest sample.
    stop = float(block["load_control_stop_force_N"])
    stop_x = float(block["load_control_stop_deflection_m"]) * mm
    valley = int(np.argmin(np.where(deflection > limit_x, force, np.inf)))
    rising_x, rising_f = deflection[valley:], force[valley:]
    ax.plot([stop_x], [stop], "s", color=st.series_color(2), markersize=6.5,
            label="load control, last converged")
    if len(rising_f) > 1 and rising_f.max() > stop:
        far = float(np.interp(stop, rising_f, rising_x))
        ax.annotate("", xy=(far, stop), xytext=(stop_x, stop),
                    arrowprops=dict(arrowstyle="->", color=st.INK_MUTED, linewidth=1.0,
                                    linestyle="--"))
        ax.annotate("the snap that load control refuses", ((stop_x + far) / 2.0, stop),
                    xytext=(0, -13), textcoords="offset points", ha="center", fontsize=8,
                    color=st.INK_SECONDARY)
    ax.set_xlabel("crown deflection [mm]")
    ax.set_ylabel("crown force on the half model [N]")
    top = 1.35 * limit
    ax.set_ylim(0.0, top)
    # The frame ends a little past where the stiffening branch leaves it.
    exit_x = float(np.interp(top, rising_f, rising_x)) if rising_f.max() > top else deflection[-1]
    ax.set_xlim(0.0, 1.08 * exit_x)
    st.title(ax, "Snap-through of a shallow arch: one path, three solution methods",
             "clamped circular arch, half span 1 m, rise 0.1 m, 20 x 20 mm section, "
             "E = 70 GPa, plane stress; half model of 60 x 4 Q4 cells, Saint "
             "Venant-Kirchhoff; the path continues to 20 kN at 224 mm")
    st.legend(ax, loc="lower right")
    st.annotate_note(
        fig,
        "The arc-length method follows the path through both limit points in "
        f"{int(block['arc_length_steps'])} steps. Displacement control, driving the crown "
        "to the same deflections, reproduces the force at every sample to "
        f"{float(block['largest_force_difference_over_limit']):.1e} of the limit load, and "
        "the force-controlled tangent has one more negative pivot exactly where the "
        f"force falls ({int(block['inertia_points_checked']) - int(block['inertia_mismatches'])} "
        f"of {int(block['inertia_points_checked'])} points). Load control stops at "
        f"{stop:.3f} N and rejects {float(block['load_control_bound_force_N']):.3f} N, "
        f"bracketing the limit found by a 60-station displacement sweep ({limit:.4f} N).",
    )
    return st.save_figure(fig, path)


def plot_elastica(directory: str, path: str) -> str:
    """The Tet10 cantilever against Euler's elastica: the tip state, and the
    error set against the continuum-elastica gap it converges to."""
    from matplotlib.lines import Line2D

    table = load_csv(os.path.join(directory, "elastica.csv"))
    orders = load_csv(os.path.join(directory, "elastica_orders.csv")).set_index("k[-]")
    finest = int(table["nx"].max())
    fig, (left, right) = st.figure(9.6, 4.4, ncols=2)
    ref = table[table["nx"] == finest].sort_values("k[-]")
    k = np.concatenate([[0.0], ref["k[-]"].to_numpy()])
    for slot, (column, label) in enumerate((("deflection", "tip deflection v / L"),
                                            ("shortening", "tip shortening / L"))):
        exact = np.concatenate([[0.0], ref[f"{column}_elastica[-]"].to_numpy()])
        left.plot(k, exact, "-", color=st.series_color(slot), linewidth=1.6)
        left.plot(ref["k[-]"], ref[f"{column}[-]"], "o", color=st.series_color(slot),
                  markersize=5.0)
        # Direct label above the curve, in text ink.
        left.annotate(label, (8.0, float(np.interp(8.0, k, exact))), xytext=(0, 8),
                      textcoords="offset points", ha="center", fontsize=8.5,
                      color=st.INK_SECONDARY)
    linear = np.linspace(0.0, 3.0, 20)
    left.plot(linear, linear / 3.0, "--", color=st.INK_MUTED, linewidth=1.0)
    handles = [Line2D([], [], color=st.INK_SECONDARY, linewidth=1.6, label="elastica"),
               Line2D([], [], color=st.INK_SECONDARY, marker="o", linestyle="none",
                      label=f"Tet10, {finest} cells"),
               Line2D([], [], color=st.INK_MUTED, linestyle="--", linewidth=1.0,
                      label="linear theory, kL/3")]
    left.legend(handles=handles, loc="lower right")
    left.set_xlabel("k = P L^2 / EI [-]")
    left.set_ylabel("displacement / L [-]")
    left.set_ylim(0.0, 1.0)
    left.set_xlim(0.0, 10.4)
    st.title(left, "Tip of the cantilever", wrap=46)

    scale = 1.0e4
    for slot, nx in enumerate(sorted(table["nx"].unique())[-2:]):
        sub = table[table["nx"] == nx].sort_values("k[-]")
        right.plot(sub["k[-]"], scale * sub["deflection_error[-]"], "o-",
                   color=st.series_color(slot), markersize=4.5, label=f"error, {nx} cells")
    gap = orders["deflection_gap[-]"]
    right.plot(gap.index, scale * gap.to_numpy(), "--", color=st.INK_SECONDARY, linewidth=1.2,
               label="their Richardson limit: continuum minus elastica")
    right.axhline(0.0, color=st.INK_MUTED, linewidth=0.8)
    right.set_ylim(-0.05, 1.4 * scale * float(gap.max()))
    right.set_xlabel("k = P L^2 / EI [-]")
    right.set_ylabel("tip deflection error [1e-4 L]")
    st.title(right, "Error, and the gap it converges to", wrap=46)
    st.legend(right, loc="upper left")
    st.figure_title(
        fig, "Euler's elastica: a Tet10 cantilever at large rotation",
        "L = 1 m, h = 0.01 m square, E = 210 GPa, nu = 0, clamped, dead tip force; Saint "
        "Venant-Kirchhoff, load control through k = 1 ... 10; 25, 50, 100 x 2 x 2 cells")
    low = float(orders[["deflection_order[-]", "shortening_order[-]"]].min().min())
    high = float(orders[["deflection_order[-]", "shortening_order[-]"]].max().max())
    st.annotate_note(
        fig,
        "At k = 10 the tip has moved 0.81 L and turned 1.43 rad, where the linear analysis "
        "puts it at 3.3 L. The errors of the 25, 50 and 100-cell meshes converge at "
        f"order {low:.2f} to {high:.2f} (Richardson, three meshes; the 25-cell errors are "
        "negative, below the right panel) to a limit that is not zero: the continuum is "
        "not the beam - shear adds 0.6 (h/L)^2 of the deflection and the Saint "
        "Venant-Kirchhoff bending moment softens with the curvature. The finest mesh sits "
        "just below that gap at every k, so what remains on it is the model's gap, not "
        "discretisation.",
    )
    return st.save_figure(fig, path)


def plot_finite_strain_tube(directory: str, path: str) -> str:
    """RMS displacement error of the tube at finite strain against the exact
    solution, per load and element, under refinement."""
    table = load_csv(os.path.join(directory, "hyperelastic_cylinder.csv"))
    summary = load_json(os.path.join(directory, "summary.json")).get("hyperelastic_cylinder", {})
    cases = [("inflation", "Inflation, neo-Hookean", "follower pressure 1.5 MPa, E = 10 MPa"),
             ("spin", "Spin, neo-Hookean", "200 rad/s, deformed-position load"),
             ("heating", "Heating, Saint Venant-Kirchhoff", "bore +100 K, alpha = 5e-4 /K")]
    fig, axes = st.figure(10.4, 4.2, ncols=3, sharey=True)
    for ax, (case, label, detail) in zip(axes, cases):
        sub_case = table[table["case"] == case]
        # The Q4 goes last and dashed: on the meshes both run, the one-cell-deep
        # Hex8 section is exactly plane strain and its errors equal the Q4's.
        for element in ("Tri3", "Hex8", "Tet10", "Quad4"):
            sub = sub_case[sub_case["element"] == element].sort_values("h[m]")
            if sub.empty:
                continue
            slot = ELEMENT_SLOTS.get(element, 4) if element != "Tri3" else 4
            ax.plot(sub["h[m]"], sub["u_rms_error[-]"],
                    ELEMENT_MARKERS[element] + ("--" if element == "Quad4" else "-"),
                    color=st.series_color(slot), markersize=4.5,
                    label=ELEMENT_NAMES[element])
        ax.set_xscale("log")
        ax.set_yscale("log")
        ticks = [0.002, 0.005, 0.01, 0.02, 0.05]
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{t:g}" for t in ticks])
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlabel("radial cell size h [m]")
        stretch = summary.get(case, {}).get("bore_hoop_stretch_exact")
        subtitle = detail + (f"; bore hoop stretch {stretch:.4f}" if stretch else "")
        st.title(ax, label, subtitle, wrap=34)
    axes[0].set_ylabel("RMS displacement error / RMS exact [-]")
    # Slope guides for orders 2 and 3, anchored at the coarsest Q4 and Tet10 points.
    for order, element in ((2.0, "Quad4"), (3.0, "Tet10")):
        sub = table[(table["case"] == "inflation") & (table["element"] == element)]
        if sub.empty:
            continue
        sub = sub.sort_values("h[m]")
        h = sub["h[m]"].to_numpy()
        e0 = float(sub["u_rms_error[-]"].to_numpy()[-1]) * 0.35
        guide = e0 * (h / h[-1]) ** order
        axes[0].plot(h, guide, ":", color=st.INK_MUTED, linewidth=1.0)
        axes[0].annotate(f"order {order:g}", (h[0], guide[0]), xytext=(2, -8),
                         textcoords="offset points", fontsize=7.5, color=st.INK_MUTED)
    st.legend(axes[-1], loc="lower right")
    st.annotate_note(
        fig,
        "Quarter section of a tube a = 0.1 m, b = 0.2 m in plane strain (the Hex8 and "
        "Tet10 sections one cell deep with u_z = 0), in ten load-control steps. The "
        "reference solves the radial equilibrium of the finite deformation exactly (a "
        "two-point boundary-value problem by shooting), so the error must vanish at the "
        "element's order: 2 for the linear elements and 3 for the Tet10, whose edge nodes "
        "lie on the curved surfaces. On the meshes both run, the Q4 errors (dashed) equal "
        "the Hex8 errors: the one-cell-deep Hex8 section held at u_z = 0 is exactly plane "
        "strain.",
    )
    return st.save_figure(fig, path)


#: The variants of the thick-cylinder collapse study, in legend order: the
#: element, whether it averages its dilatation, and whether that is its
#: default. Colour follows the element; the default variant is solid and
#: filled, the other dashed and open.
CYLINDER_VARIANTS = [("Quad4", "yes", True), ("Hex8", "yes", True), ("Tet10", "no", True),
                     ("Quad4", "no", False), ("Tet10", "yes", False), ("Tri3", "no", False)]
CYLINDER_SLOTS = {"Hex8": 0, "Tet10": 2, "Quad4": 3, "Tri3": 4}


def _variant_label(element: str, md: str) -> str:
    return ELEMENT_NAMES[element] + (" mean dilatation" if md == "yes" else " standard")


def plot_plastic_cylinder(directory: str, path: str) -> str:
    """The thick tube to plastic collapse: the pressure paths of the finest
    meshes, the collapse-pressure error under refinement, and the stress
    through the wall on the plateau against the exact fully plastic field."""
    table = load_csv(os.path.join(directory, "plastic_cylinder.csv"))
    paths = load_csv(os.path.join(directory, "plastic_cylinder_paths.csv"))
    profile = load_csv(os.path.join(directory, "plastic_cylinder_stress.csv"))
    block = load_json(os.path.join(directory, "summary.json")).get("plastic_cylinder", {})
    p_limit = float(block.get("collapse_pressure_exact_Pa", np.nan))
    first_yield = float(block.get("first_yield_pressure_exact_Pa", np.nan)) / p_limit
    fig, axes = st.figure(12.0, 4.6, ncols=3)
    mm = 1.0e3

    ax = axes[0]
    ax.axhline(1.0, color=st.INK_MUTED, linewidth=1.0, linestyle=":",
               label="exact collapse pressure p_L")
    for element, md, default in CYLINDER_VARIANTS:
        sub = paths[(paths["element"] == element) & (paths["mean_dilatation"] == md)]
        if sub.empty:
            continue
        finest = sub["n_r"].max()
        sub = sub[sub["n_r"] == finest]
        color = st.series_color(CYLINDER_SLOTS[element])
        ax.plot(sub["bore_displacement[m]"] * mm, sub["pressure_ratio[-]"],
                "-" if default else "--", color=color, linewidth=1.4,
                label=f"{_variant_label(element, md)}, n_r = {int(finest)}")
    ax.set_xlabel("bore displacement [mm]")
    ax.set_ylabel("bore pressure / p_L [-]")
    ax.set_ylim(0.985, 1.016)
    ax.set_xlim(left=0.0)
    st.title(ax, "The collapse plateau",
             f"finest mesh of each variant; first yield at {first_yield:.3f} p_L, below the "
             "frame", wrap=44)
    st.legend(ax, loc="lower right", fontsize=7)

    ax = axes[1]
    for element, md, default in CYLINDER_VARIANTS:
        sub = table[(table["element"] == element) & (table["mean_dilatation"] == md)]
        if sub.empty:
            continue
        sub = sub.sort_values("h[m]")
        color = st.series_color(CYLINDER_SLOTS[element])
        error = sub["collapse_error[-]"].to_numpy()
        below = bool((error < 0).all())
        marker = ELEMENT_MARKERS[element]
        ax.plot(sub["h[m]"], np.abs(error), ("-" if default else "--") + marker,
                color=color, markersize=5,
                markerfacecolor=color if default else st.SURFACE,
                label=_variant_label(element, md) + (" (from below)" if below else ""))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ticks = [0.003, 0.006, 0.0125, 0.025, 0.05]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("radial cell size h [m]")
    ax.set_ylabel("|collapse pressure / p_L - 1| [-]")
    for order, element, md in ((2.0, "Quad4", "yes"), (3.0, "Tet10", "no")):
        sub = table[(table["element"] == element) & (table["mean_dilatation"] == md)]
        if sub.empty:
            continue
        sub = sub.sort_values("h[m]")
        h = sub["h[m]"].to_numpy()
        e0 = abs(float(sub["collapse_error[-]"].to_numpy()[-1])) * 0.3
        guide = e0 * (h / h[-1]) ** order
        ax.plot(h, guide, ":", color=st.INK_MUTED, linewidth=1.0)
        # Labelled below the middle of the guide, clear of the frame.
        middle = len(h) // 2
        ax.annotate(f"order {order:g}", (h[middle], guide[middle]), xytext=(2, -11),
                    textcoords="offset points", fontsize=7.5, color=st.INK_MUTED)
    st.title(ax, "Collapse pressure under refinement",
             "largest load factor of the path against (2/sqrt 3) sigma_y ln(b/a)", wrap=44)
    st.legend(ax, loc="upper left", fontsize=6.8)

    ax = axes[2]
    sy = 250.0e6
    components = (("sigma_r", "radial"), ("sigma_theta", "hoop"), ("sigma_z", "axial"))
    for fill, md in ((True, "yes"), (False, "no")):
        sub = profile[(profile["element"] == "Quad4") & (profile["mean_dilatation"] == md)]
        if sub.empty:
            continue
        # One point per radial cell: the cells of a ring share their radius.
        rings = sub.groupby("r[m]").mean(numeric_only=True).reset_index()
        for slot, (key, _) in enumerate(components):
            color = st.series_color(slot)
            ax.plot(rings["r[m]"] * mm, rings[f"{key}[Pa]"] / sy, "o", color=color,
                    markersize=3.5 if fill else 4.5,
                    markerfacecolor=color if fill else st.SURFACE, markeredgewidth=1.0)
    exact = profile[(profile["element"] == "Quad4") & (profile["mean_dilatation"] == "yes")]
    exact = exact.groupby("r[m]").mean(numeric_only=True).reset_index().sort_values("r[m]")
    for slot, (key, _) in enumerate(components):
        ax.plot(exact["r[m]"] * mm, exact[f"{key}_exact[Pa]"] / sy, "-",
                color=st.series_color(slot), linewidth=1.0)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=st.series_color(slot), linewidth=1.4, label=name)
               for slot, (_, name) in enumerate(components)]
    handles += [Line2D([], [], color=st.INK_SECONDARY, marker="o", linestyle="none",
                       markersize=4, label="Q4 mean dilatation"),
                Line2D([], [], color=st.INK_SECONDARY, marker="o", linestyle="none",
                       markersize=4.5, markerfacecolor=st.SURFACE,
                       markeredgecolor=st.INK_SECONDARY, markeredgewidth=1.0,
                       label="Q4 standard"),
                Line2D([], [], color=st.INK_SECONDARY, linewidth=1.0, label="exact")]
    ax.legend(handles=handles, loc="upper left", fontsize=7, ncol=2, frameon=False)
    ax.set_xlabel("radius [mm]")
    ax.set_ylabel("stress / sigma_y [-]")
    ax.set_ylim(-0.95, 1.55)
    st.title(ax, "Stress through the wall on the plateau",
             "ring averages of the finest meshes against the exact fully plastic field",
             wrap=44)

    q4 = table[(table["element"] == "Quad4") & (table["mean_dilatation"] == "yes")]
    q4 = q4.sort_values("h[m]")
    st.annotate_note(
        fig,
        "Quarter section of a tube a = 0.1 m, b = 0.2 m in plane strain, perfectly plastic "
        "(E = 200 GPa, nu = 0.3, sigma_y = 250 MPa), small strain, loaded by its bore "
        "pressure along the arc-length path. The Hex8 section is one cell deep with "
        "u_z = 0, exactly plane strain: its results equal the Q4's. With mean dilatation "
        f"the Q4 converges to p_L at second order ({float(q4['collapse_error[-]'].iloc[0]):.1e} "
        f"at n_r = {int(q4['n_r'].iloc[0])}) and its plateau is flat; fully integrated it "
        "locks: the collapse load comes out high and the plateau keeps rising. The Tet10 "
        "converges at third order without mean dilatation, and with it at second order "
        "from below, its constant element pressure oscillating. The checkerboard Tri3 mesh "
        "does not lock here.",
    )
    return st.save_figure(fig, path)


def _bending_moment(k, youngs, sy, t, h):
    """Exact moment of an elastic-perfectly plastic rectangle at curvature k."""
    ky = 2.0 * sy / (youngs * h)
    inertia = t * h ** 3 / 12.0
    mp = sy * t * h * h / 4.0
    k = np.asarray(k, dtype=float)
    ratio = ky / np.maximum(np.abs(k), ky)
    return np.where(np.abs(k) <= ky, youngs * inertia * k,
                    np.sign(k) * mp * (1.0 - ratio ** 2 / 3.0))


def plot_plastic_bending(directory: str, path: str) -> str:
    """Pure bending of an elastic-perfectly plastic beam: moment-curvature,
    the moment error under refinement, and the residual stress after
    unloading to zero moment."""
    table = load_csv(os.path.join(directory, "plastic_bending.csv"))
    residual = load_csv(os.path.join(directory, "plastic_bending_residual.csv"))
    block = load_json(os.path.join(directory, "summary.json")).get("plastic_bending", {})
    youngs, sy, t, h = 200.0e9, 250.0e6, 0.01, 0.05
    ky = 2.0 * sy / (youngs * h)
    mp = sy * t * h * h / 4.0
    inertia = t * h ** 3 / 12.0
    k_res = float(block.get("unloaded_curvature_ratio", np.nan))
    fig, axes = st.figure(12.0, 4.4, ncols=3)

    ax = axes[0]
    k = np.linspace(0.0, 3.2, 400)
    ax.plot(k, _bending_moment(k * ky, youngs, sy, t, h) / mp, "-", color=st.INK_MUTED,
            linewidth=1.2, label="exact, loading")
    m1 = float(_bending_moment(3.0 * ky, youngs, sy, t, h))
    ax.plot([3.0, k_res], [m1 / mp, 0.0], "--", color=st.INK_MUTED, linewidth=1.2,
            label="exact, elastic unloading")
    finest = table["ny"].max()
    sub = table[table["ny"] == finest]
    ax.plot(sub["curvature_ratio[-]"], sub["moment[N m]"] / mp, "o", color=st.series_color(0),
            markersize=6, label=f"Q4, {int(finest)} cells deep")
    residual_moment = float(block["meshes"][-1]["residual_moment_over_mp"]) if block.get(
        "meshes") else np.nan
    ax.plot([k_res], [residual_moment], "s", color=st.series_color(1), markersize=6,
            label="Q4, unloaded to k_res")
    ax.axhline(1.0, color=st.GRID, linewidth=1.0)
    ax.annotate("M_p", (0.02, 1.0), xycoords=("axes fraction", "data"), xytext=(0, 3),
                textcoords="offset points", fontsize=7.5, color=st.INK_SECONDARY)
    ax.set_xlabel("curvature / first-yield curvature k_y [-]")
    ax.set_ylabel("moment / plastic moment M_p [-]")
    ax.set_ylim(-0.05, 1.1)
    st.title(ax, "Moment-curvature", "loading to 3 k_y and back to zero moment", wrap=44)
    st.legend(ax, loc="lower right", fontsize=7.5)

    ax = axes[1]
    for slot, ratio in enumerate(sorted(table["curvature_ratio[-]"].unique())):
        sub = table[table["curvature_ratio[-]"] == ratio].sort_values("h[m]")
        ax.plot(sub["h[m]"], sub["moment_error[-]"], "-o", color=st.series_color(slot),
                markersize=4.5, label=f"k = {ratio:g} k_y")
    ax.set_xscale("log")
    ax.set_yscale("log")
    hs = sorted(table["h[m]"].unique())
    ax.set_xticks(hs)
    ax.set_xticklabels([f"{x * 1e3:g}" for x in hs])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("cell size [mm]")
    ax.set_ylabel("|M - M_exact| / M_p [-]")
    # The order-2 guide runs beneath every curve, labelled below its middle;
    # the legend takes the empty upper left.
    coarse = table[table["h[m]"] == max(hs)]
    if not coarse.empty:
        hh = np.asarray(hs)
        e0 = float(coarse["moment_error[-]"].min()) * 0.4
        guide = e0 * (hh / hh[-1]) ** 2
        ax.plot(hh, guide, ":", color=st.INK_MUTED, linewidth=1.0)
        middle = len(hh) // 2
        ax.annotate("order 2", (hh[middle], guide[middle]), xytext=(4, -14),
                    textcoords="offset points", fontsize=7.5, color=st.INK_MUTED)
    st.title(ax, "Moment error under refinement", "square Q4 cells, plane stress", wrap=44)
    st.legend(ax, loc="upper left", fontsize=7.5)

    ax = axes[2]
    for slot, ny in enumerate(sorted(residual["ny"].unique())[-2:]):
        sub = residual[residual["ny"] == ny].sort_values("y[m]")
        ax.plot(sub["residual_stress[Pa]"] / sy, (sub["y[m]"] - h / 2) * 1e3, "o",
                color=st.series_color(slot), markersize=4 if slot else 3,
                label=f"Q4, {int(ny)} cells deep")
    y = np.linspace(-h / 2, h / 2, 401)
    loaded = np.clip(-youngs * 3.0 * ky * y, -sy, sy)
    exact = loaded - youngs * (k_res * ky - 3.0 * ky) * y
    ax.plot(exact / sy, y * 1e3, "-", color=st.INK_MUTED, linewidth=1.2, label="exact")
    ax.axvline(0.0, color=st.GRID, linewidth=1.0)
    ax.set_xlabel("residual stress / sigma_y [-]")
    ax.set_ylabel("distance from the neutral axis [mm]")
    st.title(ax, "Residual stress after unloading",
             f"loaded to 3 k_y, unloaded to k_res = {k_res:.3f} k_y", wrap=48)
    st.legend(ax, loc="lower right", fontsize=7.5)

    st.annotate_note(
        fig,
        "A beam 0.2 m x 0.05 m, 0.01 m thick, in plane stress (E = 200 GPa, sigma_y = 250 "
        "MPa, no hardening), small strain, bent by end displacements u_x = -k (x - L/2)(y - "
        "h/2). The section is in uniaxial stress, so M = E I k up to k_y = 2 sigma_y / (E h) "
        "and M_p (1 - (k_y / k)^2 / 3) beyond; unloaded elastically from 3 k_y to the "
        "curvature at which that moment vanishes, the beam keeps the loaded stress profile "
        "less the elastic unloading. The moment is read from the end reactions; the "
        "residual stress from the elements left of mid-span.",
    )
    return st.save_figure(fig, path)


def plot_plastic_cycle(directory: str, path: str) -> str:
    """A uniaxial strain cycle with combined hardening: the stress-strain loop
    of a distorted Hex8 bar against the exact uniaxial response."""
    table = load_csv(os.path.join(directory, "plastic_cycle.csv"))
    block = load_json(os.path.join(directory, "summary.json")).get("plastic_cycle", {})
    fig, ax = st.figure(7.2, 4.6)
    strain = np.concatenate([[0.0], table["strain[-]"].to_numpy()])
    exact = np.concatenate([[0.0], table["stress_exact[Pa]"].to_numpy()])
    ax.plot(strain * 100, exact / 1e6, "-", color=st.INK_MUTED, linewidth=1.2,
            label="exact uniaxial response")
    ax.plot(table["strain[-]"] * 100, table["stress[Pa]"] / 1e6, "o",
            color=st.series_color(0), markersize=4.5,
            label="Hex8 bar, end force / area")
    ax.axhline(0.0, color=st.GRID, linewidth=1.0)
    ax.axvline(0.0, color=st.GRID, linewidth=1.0)
    ax.set_xlabel("axial strain [%]")
    ax.set_ylabel("axial stress [MPa]")
    worst = float(block.get("max_stress_error_over_yield", np.nan))
    st.title(ax, "A strain cycle with isotropic, Voce and kinematic hardening",
             f"eps = 0 -> 1 % -> -1 % -> 1 %, 20 steps per leg; largest error "
             f"{worst:.1e} sigma_y")
    st.legend(ax, loc="lower right")
    st.annotate_note(
        fig,
        "A bar 1 m x 0.1 m x 0.1 m on a distorted 4 x 2 x 2 Hex8 mesh (mean dilatation), "
        "E = 200 GPa, nu = 0.3, sigma_y = 250 MPa, linear isotropic hardening 1 GPa plus a "
        "Voce saturation of 100 MPa at rate 30, Prager kinematic hardening 4 GPa. The "
        "reversal yields early (the Bauschinger effect of the back stress) and the loop "
        "grows as the isotropic hardening accumulates. A homogeneous state is exact on any "
        "mesh, so the error is round-off.",
    )
    return st.save_figure(fig, path)


def _order_two_guide(ax, xs, ys, anchor: float = 0.4, label_offset=(4, -14)) -> None:
    """A dotted order-2 slope beneath the curves, labelled below its middle."""
    xs = np.asarray(sorted(xs), dtype=float)
    guide = anchor * float(np.min(ys)) * (xs / xs[-1]) ** -2 if xs[0] > 1.0 else None
    if guide is None:
        return
    ax.plot(xs, guide, ":", color=st.INK_MUTED, linewidth=1.0)
    middle = len(xs) // 2
    ax.annotate("order 2", (xs[middle], guide[middle]), xytext=label_offset,
                textcoords="offset points", fontsize=7.5, color=st.INK_MUTED)


def plot_rod_dynamics(directory: str, path: str) -> str:
    """The fixed-free rod: its harmonic response against the exact damped
    continuum solution, the convergence of the harmonic and the transient
    response, and the end displacement under a ramped force."""
    sweep = load_csv(os.path.join(directory, "rod_harmonic_sweep.csv"))
    harmonic = load_csv(os.path.join(directory, "rod_harmonic.csv"))
    history = load_csv(os.path.join(directory, "rod_transient_history.csv"))
    transient = load_csv(os.path.join(directory, "rod_transient.csv"))
    # The rod of the study: steel, 1 m long, f1 = c / (4 L).
    f1 = np.sqrt(200.0e9 / 7850.0) / 4.0
    fig, axes = st.figure(11.0, 8.2, nrows=2, ncols=2)

    ax = axes[0, 0]
    ax.plot(sweep["f_over_f1"], sweep["exact_abs[m]"] * 1e6, "-", color=st.INK_MUTED,
            linewidth=1.6, label="exact continuum")
    ax.plot(sweep["f_over_f1"], sweep["q4_n40_abs[m]"] * 1e6, "-", color=st.series_color(0),
            linewidth=1.0, label="40 Q4 elements")
    ax.plot(sweep["f_over_f1"], sweep["q4_n10_abs[m]"] * 1e6, "--", color=st.series_color(1),
            linewidth=1.0, label="10 Q4 elements")
    ax.set_yscale("log")
    ax.set_xlabel("frequency / f1 [-]")
    ax.set_ylabel("end displacement amplitude [um]")
    st.title(ax, "Harmonic response of a fixed-free rod",
             f"end force 1 kN, structural damping eta = 0.02; f1 = c / (4 L) = {f1:.1f} Hz",
             wrap=48)
    st.legend(ax, loc="upper right", fontsize=7.5)

    ax = axes[0, 1]
    picks = [("end force undamped", 0.5), ("end force undamped", 2.5),
             ("end force eta 0.02 at resonance", 1.0), ("base motion undamped", 1.5)]
    q4 = harmonic[(harmonic["element"] == "Q4") & (harmonic["mass"] == "consistent")]
    lowest = []
    for slot, (case, ratio) in enumerate(picks):
        sub = q4[(q4["case"] == case) & (np.isclose(q4["f_over_f1"], ratio))].sort_values("n")
        if sub.empty:
            continue
        ax.plot(sub["n"], sub["continuum_error[-]"], "-o", color=st.series_color(slot),
                markersize=4.5, label=f"{case}, f = {ratio:g} f1")
        lowest.append(float(sub["continuum_error[-]"].min()))
    ns = sorted(q4["n"].unique())
    if lowest:
        _order_two_guide(ax, ns, [min(lowest)])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(ns)
    ax.set_xticklabels([f"{int(n)}" for n in ns])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("elements along the rod")
    ax.set_ylabel("|U(L) - U_exact(L)| / |U_exact(L)| [-]")
    worst = float(harmonic["discrete_difference[-]"].max())
    st.title(ax, "Convergence to the continuum",
             f"Q4, consistent mass; against the exact discrete solution {worst:.1e}", wrap=48)
    st.legend(ax, loc="lower left", fontsize=7.0)

    ax = axes[1, 0]
    t = history["t[s]"].to_numpy() * 1e3
    ax.plot(t, history["exact[m]"] * 1e6, "-", color=st.INK_MUTED, linewidth=1.6,
            label="exact (modal series)")
    ax.plot(t, history["q4_n160[m]"] * 1e6, "-", color=st.series_color(0), linewidth=1.0,
            label="160 Q4, dt = h / c")
    coarse = history.dropna(subset=["q4_n20[m]"])
    ax.plot(coarse["t[s]"] * 1e3, coarse["q4_n20[m]"] * 1e6, "o", color=st.series_color(1),
            markersize=2.5, label="20 Q4, dt = h / c")
    ax.set_xlabel("time [ms]")
    ax.set_ylabel("end displacement [um]")
    st.title(ax, "End displacement under a ramped force",
             "1 kN ramped as sin^2 over 0.6 T1, trapezoidal rule", wrap=48)
    st.legend(ax, loc="lower right", fontsize=7.5)

    ax = axes[1, 1]
    lowest = []
    for slot, mass in enumerate(["consistent", "lumped"]):
        sub = transient[(transient["element"] == "Q4") & (transient["mass"] == mass)]
        sub = sub.sort_values("n")
        ax.plot(sub["n"], sub["max_error[-]"], "-o", color=st.series_color(slot),
                markersize=4.5, label=f"Q4, {mass} mass")
        lowest.append(float(sub["max_error[-]"].min()))
    ns = sorted(transient["n"].unique())
    _order_two_guide(ax, ns, [min(lowest)], anchor=0.2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(ns)
    ax.set_xticklabels([f"{int(n)}" for n in ns])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("elements along the rod (the step halved with the element)")
    ax.set_ylabel("max |u(L) - u_exact(L)| / max |u_exact(L)| [-]")
    st.title(ax, "Transient convergence", "Courant number 1, 2.5 T1", wrap=48)
    st.legend(ax, loc="lower left", fontsize=7.5)

    st.annotate_note(
        fig,
        "A steel rod 1 m long, 0.05 m square, nu = 0, its lateral displacements held: an "
        "exactly one-dimensional model (Hex8 gives the same numbers). Harmonic: every node "
        "equals the exact solution of the discrete equations (their dispersion relation) to "
        "round-off, and the end amplitude converges at second order to u = F sin(kx) / (E* A "
        "k cos(kL)) with the complex modulus of the damping. Transient: the end displacement "
        "converges at second order, the element length and the time step halved together, to "
        "the modal series of the continuum (4000 modes, the static part summed exactly).",
    )
    return st.save_figure(fig, path)


def plot_nonlinear_oscillator(directory: str, path: str) -> str:
    """One element in uniaxial strain under a sudden load: the motion of the
    finite-strain elastic and of the elastoplastic oscillator against the
    exact one, and the convergence of the motion and of the plastic
    dissipation with the step."""
    history = load_csv(os.path.join(directory, "nonlinear_oscillator_history.csv"))
    table = load_csv(os.path.join(directory, "nonlinear_oscillator.csv"))
    fig, axes = st.figure(11.0, 8.2, nrows=2, ncols=2)
    laws = ["Saint Venant-Kirchhoff finite strain", "J2 linear hardening small strain"]
    titles = ["Saint Venant-Kirchhoff, finite strain", "J2 with linear hardening"]
    subtitles = ["a sudden pull, 9 % peak strain; exact: Runge-Kutta at 1/64 step",
                 "a sudden load of 0.8 N_y yields on the first swing; exact: closed form"]
    for panel, (law, name, sub) in enumerate(zip(laws, titles, subtitles)):
        ax = axes[0, panel]
        rows = history[history["law"] == law]
        t = rows["t[s]"].to_numpy() * 1e6
        # The exact motion as a wide halo: the fine run lies on it.
        ax.plot(t, rows["exact[m]"] * 1e3, "-", color=st.INK_MUTED, linewidth=4.0, alpha=0.45,
                label="exact")
        ax.plot(t, rows["fine[m]"] * 1e3, "-", color=st.series_color(0), linewidth=1.0,
                label="320 steps per period")
        coarse = rows.dropna(subset=["coarse[m]"])
        ax.plot(coarse["t[s]"] * 1e6, coarse["coarse[m]"] * 1e3, "o", color=st.series_color(1),
                markersize=3.5, label="20 steps per period")
        ax.set_xlabel("time [us]")
        ax.set_ylabel("end displacement [mm]")
        st.title(ax, name, sub, wrap=48)
        st.legend(ax, loc="lower right", fontsize=7.5)

    ax = axes[1, 0]
    q4 = table[(table["element"] == "Q4") & (table["mass"] == "lumped")]
    slot = 0
    lowest = []
    for law, name in zip(laws, ["SVK", "J2"]):
        for alpha in sorted(q4["alpha"].unique(), reverse=True):
            sub = q4[(q4["law"] == law) & (q4["alpha"] == alpha)].sort_values("steps_per_period")
            ax.plot(sub["steps_per_period"], sub["error[-]"], "-o", color=st.series_color(slot),
                    markersize=4.5, label=f"{name}, alpha = {alpha:g}")
            lowest.append(float(sub["error[-]"].min()))
            slot += 1
    steps = sorted(q4["steps_per_period"].unique())
    _order_two_guide(ax, steps, [min(lowest)])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(steps)
    ax.set_xticklabels([f"{int(s)}" for s in steps])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("steps per elastic period")
    ax.set_ylabel("max |u - u_exact| / max |u| [-]")
    worst = float(table["discrete_difference[-]"].max())
    st.title(ax, "Convergence of the motion",
             f"against the scalar HHT-alpha recursion of the same equation: {worst:.1e}",
             wrap=48)
    st.legend(ax, loc="lower left", fontsize=7.5)

    ax = axes[1, 1]
    plastic = q4[q4["law"] == laws[1]]
    for slot, alpha in enumerate(sorted(plastic["alpha"].unique(), reverse=True)):
        sub = plastic[plastic["alpha"] == alpha].sort_values("steps_per_period")
        gap = np.abs(sub["final_balance[J]"] - sub["exact_dissipation[J]"]) / \
            sub["exact_dissipation[J]"]
        ax.plot(sub["steps_per_period"], gap, "-o", color=st.series_color(slot),
                markersize=4.5, label=f"alpha = {alpha:g}")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(steps)
    ax.set_xticklabels([f"{int(s)}" for s in steps])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("steps per elastic period")
    ax.set_ylabel("|energy balance - D_p| / D_p [-]")
    exact = float(plastic["exact_dissipation[J]"].iloc[0]) if not plastic.empty else np.nan
    st.title(ax, "The energy balance holds the plastic dissipation",
             f"E_0 + W - T - U against D_p = sigma_y alpha_p V = {exact:.2f} J; the gap closes "
             "irregularly as the yield point moves within a step", wrap=48)
    st.legend(ax, loc="lower left", fontsize=7.5)

    st.annotate_note(
        fig,
        "One element 0.1 m, steel (E = 200 GPa, nu = 0.3), its lateral displacements held: "
        "uniaxial strain, a single degree of freedom (Q4 and Hex8, lumped and consistent mass "
        "give the same relative errors). The trapezoidal rule and HHT-alpha converge at second "
        "order in the step for both laws; the elastoplastic motion's order wanders about 2 as "
        "the yield point falls at a different place within a step. The final energy balance "
        "of the plastic run tends to the plastic dissipation, the plastic work less the stored "
        "hardening energy.",
    )
    return st.save_figure(fig, path)
