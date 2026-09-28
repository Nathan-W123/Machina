#!/usr/bin/env python3
"""Turn the run summaries into the small Markdown/CSV tables the docs quote.

usage:
    python3 python/scripts/write_result_tables.py --results results \
                                                  --output docs/results

Every number written here is read straight out of a `summary.json`, so the
documentation and the solver output cannot drift apart. Re-run this after any
benchmark run and the tables update.
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Optional

import _bootstrap  # noqa: F401

import pandas as pd

from sparlab_viz.loaders import ResultError, load_json


BENCHMARK_CASES = ["cantilever_beam", "mbb_beam", "mbb_beam_projected", "aerospace_bracket",
                   "wing_rib", "l_bracket_stress", "bracket_3d", "bracket_3d_projected",
                   "lug_bracket_2d", "engine_mount_3d", "bracket_3d_large", "column_buckling",
                   "mbb_beam_robust", "mbb_beam_overhang", "bracket_3d_overhang"]
#: The buckling-constrained column and its comparison runs (the same deck with
#: one feature switched off; scripts/run_all_benchmarks.sh).
BUCKLING_CASES = [("column_buckling_unconstrained", "compliance only"),
                  ("column_buckling_nonrobust", "lambda >= 6, plain projection"),
                  ("column_buckling", "lambda >= 6, robust projection")]
ROBUST_CASES = [("mbb_beam_robust_off", "plain projection, erosion check"),
                ("mbb_beam_robust", "robust formulation"),
                ("column_buckling", "robust formulation, buckling constraint")]
OVERHANG_CASES = [("mbb_beam_overhang_off", "no filter"),
                  ("mbb_beam_overhang", "overhang filter"),
                  ("mbb_beam_overhang_down", "overhang filter"),
                  ("bracket_3d_overhang_off", "no filter"),
                  ("bracket_3d_overhang", "overhang filter")]
ANALYSIS_CASES = ["cantilever_analysis", "block_3d_analysis"]
#: Runs on meshes read from files (Gmsh / Abaqus-CalculiX input).
REAL_GEOMETRY_CASES = ["lug_bracket_2d", "engine_mount_3d"]
ELEMENT_LABELS = {"Quad4": "Q4", "Tri3": "Tri3", "Hex8": "Hex8", "Tet4": "Tet4",
                  "Tet10": "Tet10"}
#: The unconstrained run of the stress-constrained deck (sparlab_topopt
#: --no-stress) that the stress table sets beside it.
STRESS_REFERENCE = {"l_bracket_stress": "l_bracket_unconstrained"}


def _fmt(value, digits: int = 4) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, (int,)):
        return str(value)
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if number != number:  # NaN
        return "n/a"
    if number == 0:
        return "0"
    if 1e-3 <= abs(number) < 1e5:
        return f"{number:.{digits}g}"
    return f"{number:.{digits}g}"


def _markdown_table(headers: List[str], rows: List[List[str]]) -> str:
    widths = [len(h) for h in headers]
    for row in rows:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(cell))
    lines = [
        "| " + " | ".join(h.ljust(widths[i]) for i, h in enumerate(headers)) + " |",
        "|" + "|".join("-" * (w + 2) for w in widths) + "|",
    ]
    for row in rows:
        lines.append(
            "| " + " | ".join(cell.ljust(widths[i]) for i, cell in enumerate(row)) + " |"
        )
    return "\n".join(lines)


def benchmark_table(results_dir: str) -> Optional[str]:
    rows = []
    frames = []
    for case in BENCHMARK_CASES:
        path = os.path.join(results_dir, case, "summary.json")
        if not os.path.isfile(path):
            continue
        doc = load_json(path)
        mesh = doc.get("mesh", {})
        setup = doc.get("optimization_setup", {})
        result = doc.get("optimization_result", {})
        comparison = doc.get("mass_stiffness_comparison", {})
        interp = doc.get("solid_interpretation", {})
        record = {
            "case": case,
            "dim": mesh.get("dim", 2),
            "method": setup.get("method", "oc"),
            "elements": mesh.get("num_elements"),
            "dofs": mesh.get("num_dofs"),
            "volume_fraction_target": setup.get("volume_fraction_target"),
            "volume_fraction_achieved": result.get("volume_fraction"),
            "volume_violation": result.get("volume_constraint_relative_violation"),
            "iterations": result.get("iterations"),
            "stop_reason": result.get("stop_reason"),
            "compliance_J": result.get("compliance_J"),
            "initial_compliance_J": result.get("initial_compliance_J"),
            "full_solid_compliance_J": comparison.get("full_solid_compliance_J"),
            "equal_mass_plate_compliance_J":
                comparison.get("equal_mass_uniform_plate_compliance_J"),
            "stiffness_gain_over_equal_mass_plate":
                comparison.get("stiffness_gain_over_equal_mass_plate"),
            "grey_level": result.get("grey_level"),
            "retained_elements": interp.get("elements_retained"),
            "connected_groups": interp.get("connected_groups_above_threshold"),
            "islands_discarded_m3": interp.get("volume_discarded_as_islands_m3"),
            "runtime_s": result.get("total_seconds"),
            "s_per_iteration": result.get("seconds_per_iteration"),
        }
        solid = doc.get("modal_initial_solid")
        topology = doc.get("modal_optimised_topology")
        if solid and solid.get("frequencies_hz"):
            record["f1_solid_Hz"] = solid["frequencies_hz"][0]
            record["mass_solid_kg"] = solid.get("total_mass_kg")
        if topology and topology.get("frequencies_hz"):
            record["f1_topology_Hz"] = topology["frequencies_hz"][0]
            record["mass_topology_kg"] = topology.get("total_mass_kg")
        element = mesh.get("element_type") or ("Hex8" if record["dim"] == 3 else "Quad4")
        record["element_type"] = element
        solver = (result.get("linear_solver") or {}).get("solver")
        record["linear_solver"] = solver
        projection = result.get("projection") or {}
        record["final_beta"] = projection.get("final_beta")
        frames.append(record)
        rows.append([
            case,
            f"{record['elements']} {ELEMENT_LABELS.get(element, element)}",
            str(record["method"]),
            _fmt(record["volume_fraction_target"]),
            _fmt(record["volume_fraction_achieved"], 6),
            _fmt(record["volume_violation"], 2),
            _fmt(record["iterations"]),
            str(record["stop_reason"]),
            _fmt(record["compliance_J"], 6),
            _fmt(record["equal_mass_plate_compliance_J"], 6),
            _fmt(record["stiffness_gain_over_equal_mass_plate"], 4),
            _fmt(record["grey_level"], 3),
            _fmt(record["runtime_s"], 4),
        ])

    if not rows:
        return None
    pd.DataFrame(frames).to_csv(
        os.path.join(_OUTPUT, "benchmarks.csv"), index=False
    )
    headers = [
        "case", "elements", "method", "vf target", "vf achieved", "vf violation",
        "iterations", "stop reason", "compliance [J]",
        "equal-mass plate [J]", "stiffness gain", "grey", "runtime [s]",
    ]
    return _markdown_table(headers, rows)


def stress_table(results_dir: str) -> Optional[str]:
    """Stress-constrained runs beside the unconstrained run of the same deck."""
    rows = []
    for case in BENCHMARK_CASES:
        path = os.path.join(results_dir, case, "summary.json")
        if not os.path.isfile(path):
            continue
        doc = load_json(path)
        result = doc.get("optimization_result", {})
        block = result.get("stress")
        if not block:
            continue
        interpreted = doc.get("interpreted_solid_analysis", {})
        limit = interpreted.get("stress_limit_Pa")
        first = (block.get("load_cases") or [{}])[0]
        rows.append([
            case, "on", _fmt(limit, 4),
            _fmt(result.get("compliance_J"), 6),
            _fmt(block.get("max_relaxed_stress_ratio"), 5),
            _fmt(first.get("p_norm_ratio"), 5),
            _fmt(first.get("scale"), 5),
            "yes" if result.get("feasible") else "no",
            _fmt(interpreted.get("max_von_mises_Pa"), 5),
            _fmt(interpreted.get("max_von_mises_over_limit"), 4),
            _fmt(result.get("iterations")),
        ])
        reference = STRESS_REFERENCE.get(case)
        ref_path = os.path.join(results_dir, reference or "", "summary.json")
        if reference and os.path.isfile(ref_path):
            ref = load_json(ref_path)
            ref_result = ref.get("optimization_result", {})
            ref_max = ref.get("interpreted_solid_analysis", {}).get("max_von_mises_Pa")
            rows.append([
                f"{case} (constraint off)", "off", _fmt(limit, 4),
                _fmt(ref_result.get("compliance_J"), 6),
                "n/a", "n/a", "n/a", "n/a",
                _fmt(ref_max, 5),
                _fmt(ref_max / limit if (ref_max and limit) else None, 4),
                _fmt(ref_result.get("iterations")),
            ])
    if not rows:
        return None
    return _markdown_table(
        ["case", "constraint", "limit [Pa]", "compliance [J]", "max relaxed ratio",
         "p-norm ratio", "p-norm scale", "feasible", "re-solve max vM [Pa]",
         "re-solve / limit", "iterations"],
        rows,
    )


def cross_validation_table(results_dir: str) -> Optional[str]:
    path = os.path.join(results_dir, "cross_validation", "summary.json")
    if not os.path.isfile(path):
        return None
    doc = load_json(path)
    versions = {name: entry.get("version", "") for name, entry in doc.get("codes", {}).items()}
    rows = []
    frames = []
    for case in doc.get("cases", []):
        for load_case in case.get("load_cases", []):
            for code, entry in load_case.get("codes", {}).items():
                # passed is None for a comparison between two idealisations
                # (a CalculiX plane element expanded through the thickness at
                # nu != 0): recorded, not judged.
                verdict = ("INFO" if entry.get("passed") is None
                           else "PASS" if entry["passed"] else "FAIL")
                judged = entry.get("max_rel_diff_judged", entry["max_rel_diff"])
                rows.append([
                    case["case"], case.get("element_type", ""), load_case["load_case"],
                    f"{code} {versions.get(code, '')}".strip(), entry.get("element", ""),
                    _fmt(judged, 3), _fmt(entry.get("rms_rel_diff"), 3),
                    _fmt(entry["tolerance"], 2), verdict,
                ])
                flat = {k: v for k, v in entry.items() if not isinstance(v, (dict, list))
                        or k in ("reference_load_factors", "sparlab_load_factors",
                                 "per_mode_rel_diff")}
                if "max_rel_diff_judged" in entry:
                    flat["sparlab_vs_calculix"] = entry["max_rel_diff"]
                if isinstance(entry.get("loads"), list):
                    flat["loads"] = "; ".join(entry["loads"])
                frames.append({"case": case["case"], "load_case": load_case["load_case"],
                               "code": code, **flat})
    if not rows:
        return None
    pd.DataFrame(frames).to_csv(os.path.join(_OUTPUT, "cross_validation.csv"), index=False)
    return _markdown_table(
        ["case", "SparLab element", "load case", "code", "reference element",
         "max rel diff", "RMS rel diff", "tolerance", "result"],
        rows,
    )


#: Arm name -> (what it varies, the summary field that records it).
STUDY_ARMS = {
    "volume_fraction": ("volume-fraction target", "volume_fraction_target"),
    "load_weighting": ("lateral load-case weight", "weight_lateral"),
    "mesh_fixed_r": ("elements, filter radius fixed in metres", "num_elements"),
    "mesh_fixed_cells": ("elements, filter radius fixed in cells", "num_elements"),
    "penalty": ("SIMP penalty p", "penalty"),
    "filter_radius": ("filter radius [m]", "filter_radius_m"),
    "youngs_modulus": ("Young's modulus [Pa]", "youngs_modulus_Pa"),
}


def study_table(study_dir: str) -> Optional[str]:
    """One row per design-study arm, straight from the study summaries."""
    from sparlab_viz import studies

    try:
        table = studies.collect_study(study_dir)
    except (ResultError, FileNotFoundError):
        return None

    rows = []
    for arm, (what, column) in STUDY_ARMS.items():
        sub = table[table["arm"] == arm]
        if sub.empty:
            continue
        span = "n/a"
        if column in sub.columns and sub[column].notna().any():
            span = f"{_fmt(sub[column].min())} - {_fmt(sub[column].max())}"
        rows.append([
            arm,
            what,
            str(len(sub)),
            span,
            f"{int(sub['converged'].sum())}/{len(sub)}",
            f"{_fmt(sub['compliance_J'].min(), 6)} - {_fmt(sub['compliance_J'].max(), 6)}",
            f"{_fmt(sub['grey_level'].min(), 3)} - {_fmt(sub['grey_level'].max(), 3)}",
            str(int(sub["connected_groups"].max())),
            _fmt(sub["seconds"].sum(), 4),
        ])
    if not rows:
        return None
    headers = ["arm", "varies", "points", "range", "converged", "compliance [J]",
               "grey", "worst group count", "total runtime [s]"]
    return _markdown_table(headers, rows)


def verification_table(results_dir: str) -> Optional[str]:
    path = os.path.join(results_dir, "verification", "summary.json")
    if not os.path.isfile(path):
        return None
    doc = load_json(path)
    rows = []
    frames = []
    for outcome in doc.get("outcomes", []):
        rows.append([
            outcome["study"],
            outcome["kind"],
            outcome["metric"],
            _fmt(outcome["value"], 4),
            _fmt(outcome["tolerance"], 3),
            "PASS" if outcome["passed"] else "FAIL",
        ])
        frames.append(outcome)
    if not rows:
        return None
    pd.DataFrame(frames).to_csv(
        os.path.join(_OUTPUT, "verification.csv"), index=False
    )
    return _markdown_table(
        ["study", "kind", "metric", "value", "tolerance", "result"], rows
    )


def modal_table(results_dir: str) -> Optional[str]:
    rows = []
    for case in BENCHMARK_CASES + ANALYSIS_CASES:
        path = os.path.join(results_dir, case, "summary.json")
        if not os.path.isfile(path):
            continue
        doc = load_json(path)
        for key, label in (
            ("modal", "static model"),
            ("modal_initial_solid", "full solid domain"),
            ("modal_mass_matched_baseline", "equal-mass uniform plate"),
            ("modal_optimised_topology", "optimised topology"),
        ):
            block = doc.get(key)
            if not block or not block.get("frequencies_hz"):
                continue
            freqs = block["frequencies_hz"][:4]
            rows.append([
                case, label,
                _fmt(block.get("total_mass_kg"), 5),
                *[_fmt(f, 5) for f in freqs],
                *["n/a"] * (4 - len(freqs)),
                _fmt(max(block.get("eigenpair_residuals", [0.0])), 2),
                "yes" if block.get("converged") else "no",
            ])
    if not rows:
        return None
    return _markdown_table(
        ["case", "structure", "mass [kg]", "f1 [Hz]", "f2 [Hz]", "f3 [Hz]",
         "f4 [Hz]", "max residual", "converged"],
        rows,
    )


def scaling_table(results_dir: str, stem: str = "runtime_scaling") -> Optional[str]:
    path = os.path.join(results_dir, "benchmark", f"{stem}.json")
    if not os.path.isfile(path):
        return None
    doc = load_json(path)
    rows = []
    for record in doc.get("records", []):
        mesh = f"{record['nx']} x {record['ny']}"
        if record.get("nz"):
            mesh += f" x {record['nz']}"
        rows.append([
            mesh,
            _fmt(record["num_elements"]),
            _fmt(record["num_dofs"]),
            _fmt(record["stiffness_nonzeros"]),
            _fmt(record["assemble_s"], 3),
            _fmt(record["factorize_s"], 3),
            _fmt(record["solve_s"], 3),
            _fmt(record["objective_gradient_s"], 3),
        ])
    if not rows:
        return None
    exponents = doc.get("scaling_exponent_vs_dofs", {})
    table = _markdown_table(
        ["mesh", "elements", "DOFs", "nonzeros", "assemble [s]",
         "factorise [s]", "solve [s]", "obj+grad [s]"],
        rows,
    )
    slope_line = ", ".join(
        f"{k} {v:.2f}" for k, v in exponents.items()
    )
    return table + f"\n\nFitted slopes of log(time) vs log(DOFs): {slope_line}.\n"


def _mesh_label(record) -> str:
    text = f"{int(record['nx'])} x {int(record['ny'])}"
    if record.get("nz"):
        text += f" x {int(record['nz'])}"
    return text


def solver_comparison_table(results_dir: str, element: str) -> Optional[str]:
    """Cholesky, multigrid CG and Jacobi CG side by side at every benchmark size
    of one element type, merged on the number of DOFs."""
    from sparlab_viz import studies

    tables = studies.load_solver_scaling(os.path.join(results_dir, "benchmark"))
    series = {solver: tables[(element, solver)][0]
              for solver in ("ldlt", "amg", "jacobi") if (element, solver) in tables}
    if not series:
        return None

    long_rows = []
    for (elem, solver), (table, _meta) in tables.items():
        for _, record in table.iterrows():
            long_rows.append({"element": elem, "solver": solver, **record.to_dict()})
    pd.DataFrame(long_rows).to_csv(os.path.join(_OUTPUT, "solver_scaling.csv"), index=False)

    def lookup(solver: str, dofs: int):
        table = series.get(solver)
        if table is None:
            return None
        match = table[table["num_dofs"].astype(int) == dofs]
        return None if match.empty else match.iloc[0]

    all_dofs = sorted({int(n) for table in series.values() for n in table["num_dofs"]})
    rows = []
    for dofs in all_dofs:
        direct, multigrid, jacobi = (lookup(s, dofs) for s in ("ldlt", "amg", "jacobi"))
        any_record = next(r for r in (direct, multigrid, jacobi) if r is not None)
        cells = [_fmt(dofs), _mesh_label(any_record)]
        if direct is not None:
            cells += [_fmt(direct["factorize[s]"] + direct["solve[s]"], 3),
                      _fmt(direct["solver_nonzeros"] / dofs, 3)]
        else:
            cells += ["-", "-"]
        if multigrid is not None:
            cells += [_fmt(multigrid["factorize[s]"], 3), _fmt(multigrid["solve[s]"], 3),
                      f"{int(multigrid['iterations'])} ({int(multigrid['levels'])})",
                      _fmt(multigrid["operator_complexity"], 3),
                      _fmt(multigrid["solver_nonzeros"] / dofs, 3)]
        else:
            cells += ["-"] * 5
        if jacobi is not None:
            cells += [_fmt(jacobi["factorize[s]"] + jacobi["solve[s]"], 3),
                      str(int(jacobi["iterations"]))]
        else:
            cells += ["-", "-"]
        if direct is not None and multigrid is not None:
            cells.append(_fmt((direct["factorize[s]"] + direct["solve[s]"])
                              / (multigrid["factorize[s]"] + multigrid["solve[s]"]), 3))
        else:
            cells.append("-")
        rows.append(cells)
    return _markdown_table(
        ["DOFs", "mesh", "Cholesky [s]", "Cholesky entries/DOF", "MG setup [s]",
         "MG solve [s]", "MG iterations (levels)", "MG operator complexity",
         "MG entries/DOF", "Jacobi CG [s]", "Jacobi iterations",
         "Cholesky / MG time"],
        rows,
    )


def multigrid_table(results_dir: str) -> Optional[str]:
    path = os.path.join(results_dir, "verification", "multigrid_scaling.csv")
    if not os.path.isfile(path):
        return None
    table = pd.read_csv(path)
    rows = []
    for _, r in table.iterrows():
        rows.append([
            str(r["element"]), str(int(r["nx"])), _fmt(int(r["num_dofs"])),
            str(int(r["amg_levels"])), _fmt(r["operator_complexity"], 3),
            str(int(r["amg_iterations"])), _fmt(r["amg_seconds[s]"], 3),
            str(int(r["jacobi_iterations"])), _fmt(r["jacobi_seconds[s]"], 3),
            _fmt(r["ldlt_seconds[s]"], 3), _fmt(r["relative_difference_vs_ldlt[-]"], 2),
        ])
    return _markdown_table(
        ["element", "nx", "DOFs", "MG levels", "operator complexity", "MG iterations",
         "MG [s]", "Jacobi iterations", "Jacobi [s]", "Cholesky [s]",
         "difference vs Cholesky"],
        rows,
    )


def simplex_convergence_table(results_dir: str) -> Optional[str]:
    path = os.path.join(results_dir, "verification", "mesh_convergence_simplex.csv")
    if not os.path.isfile(path):
        return None
    table = pd.read_csv(path)
    rows = []
    for _, r in table.iterrows():
        mesh = f"{int(r['nx'])} x {int(r['ny'])}"
        if int(r["nz"]) > 0:
            mesh += f" x {int(r['nz'])}"
        rows.append([
            str(r["element"]), mesh, _fmt(int(r["num_elements"])), _fmt(int(r["num_dofs"])),
            _fmt(r["h[m]"], 4), _fmt(r["tip_mean[m]"], 6), _fmt(r["timoshenko[m]"], 6),
            _fmt(r["rel_error_timoshenko[-]"], 3), str(r["solver"]),
        ])
    return _markdown_table(
        ["element", "cells", "elements", "DOFs", "h [m]", "tip deflection [m]",
         "Timoshenko [m]", "relative error", "solver"],
        rows,
    )


def projection_table(results_dir: str) -> Optional[str]:
    from sparlab_viz import studies

    try:
        table = studies.collect_projection(results_dir)
    except (ResultError, FileNotFoundError):
        return None
    table.to_csv(os.path.join(_OUTPUT, "projection.csv"), index=False)
    rows = []
    for _, r in table.iterrows():
        beta = r["final_beta"]
        rows.append([
            str(r["problem"]), str(r["run"]),
            "off" if beta is None or beta != beta else f"beta {beta:g}",
            _fmt(r["iterations"]), str(r["stop_reason"]),
            _fmt(r["grey_level"], 3), _fmt(r["filtered_grey_level"], 3),
            _fmt(r["compliance_J"], 5), _fmt(r["interpreted_compliance_J"], 5),
            _fmt(r["compliance_vs_simp_ratio"], 4), _fmt(r["seconds"], 4),
        ])
    return _markdown_table(
        ["problem", "run", "projection", "iterations", "stop reason", "grey",
         "grey before projection", "SIMP compliance [J]", "thresholded compliance [J]",
         "thresholded / SIMP", "runtime [s]"],
        rows,
    )


def real_geometry_table(results_dir: str) -> Optional[str]:
    """The runs whose mesh was read from a file: what was read and how it solved."""
    rows = []
    for case in REAL_GEOMETRY_CASES:
        path = os.path.join(results_dir, case, "summary.json")
        if not os.path.isfile(path):
            continue
        doc = load_json(path)
        mesh = doc.get("mesh", {})
        source = mesh.get("file", {})
        quality = mesh.get("quality", {})
        result = doc.get("optimization_result", {})
        solver = result.get("linear_solver") or {}
        solid = doc.get("modal_initial_solid") or {}
        topology = doc.get("modal_optimised_topology") or {}
        f_solid = (solid.get("frequencies_hz") or [None])[0]
        f_topology = (topology.get("frequencies_hz") or [None])[0]
        rows.append([
            case,
            f"{os.path.basename(source.get('path', '?'))} ({source.get('format', '?')} "
            f"{source.get('version', '')})".replace(" )", ")"),
            ELEMENT_LABELS.get(mesh.get("element_type"), str(mesh.get("element_type"))),
            _fmt(mesh.get("num_nodes")), _fmt(mesh.get("num_elements")),
            _fmt(mesh.get("num_dofs")),
            f"{_fmt(quality.get('min'), 3)} / {_fmt(quality.get('mean'), 3)}",
            str(solver.get("solver", "n/a")),
            ("-" if "LDLT" in str(solver.get("solver", ""))
             else _fmt(solver.get("iterative_iterations_per_solve"), 3)),
            _fmt(result.get("iterations")),
            _fmt(result.get("total_seconds"), 4),
            f"{_fmt(f_solid, 5)} / {_fmt(f_topology, 5)}",
        ])
    if not rows:
        return None
    return _markdown_table(
        ["case", "mesh file", "element", "nodes", "elements", "DOFs",
         "quality min / mean", "linear solver", "CG iterations per solve",
         "optimiser iterations", "runtime [s]", "f1 solid / topology [Hz]"],
        rows,
    )


def _first_factor(block) -> Optional[float]:
    if not block:
        return None
    for entry in block.get("load_cases", []):
        factors = entry.get("load_factors") or []
        if factors:
            return float(factors[0])
    return None


def buckling_table(results_dir: str) -> Optional[str]:
    """The buckling-constrained column beside its comparison runs."""
    rows, frames = [], []
    for case, label in BUCKLING_CASES:
        path = os.path.join(results_dir, case, "summary.json")
        if not os.path.isfile(path):
            continue
        doc = load_json(path)
        result = doc.get("optimization_result", {})
        check = doc.get("buckling_check", {})
        simp = result.get("buckling") or {}
        simp_cases = simp.get("load_cases") or [{}]
        fraction = (simp_cases[0].get("solid_energy_fraction") or [None])[0]
        record = {
            "case": case, "run": label,
            "compliance_J": result.get("compliance_J"),
            "simp_lambda_1": simp.get("min_load_factor"),
            "simp_mode_1_solid_energy_fraction": fraction,
            "part_lambda_1": _first_factor(check.get("interpreted_structure")),
            "full_solid_lambda_1": _first_factor(check.get("full_solid")),
            "iterations": result.get("iterations"),
            "stop_reason": result.get("stop_reason"),
            "linear_solves": result.get("linear_solves"),
            "runtime_s": result.get("total_seconds"),
        }
        frames.append(record)
        rows.append([case, label, _fmt(record["compliance_J"], 6),
                     _fmt(record["simp_lambda_1"], 5), _fmt(fraction, 3),
                     _fmt(record["part_lambda_1"], 5), _fmt(record["full_solid_lambda_1"], 5),
                     f"{record['iterations']}, {record['stop_reason']}",
                     _fmt(record["linear_solves"])])
    if not rows:
        return None
    pd.DataFrame(frames).to_csv(os.path.join(_OUTPUT, "buckling.csv"), index=False)
    return _markdown_table(
        ["case", "run", "compliance [J]", "SIMP lambda_1", "mode 1 energy in solid",
         "exported part lambda_1", "full solid lambda_1", "iterations, stop",
         "linear solves"], rows)


def robust_table(results_dir: str) -> Optional[str]:
    """Eroded / blueprint / dilated designs and the measured length scales."""
    rows, frames = [], []
    for case, label in ROBUST_CASES:
        path = os.path.join(results_dir, case, "summary.json")
        if not os.path.isfile(path):
            continue
        doc = load_json(path)
        result = doc.get("optimization_result", {})
        block = result.get("robust") or result.get("erosion_check")
        if not block:
            continue
        designs = {d["design"].split(" ")[0]: d for d in block["designs"]}
        cell = float(doc.get("mesh", {}).get("mean_element_size_m") or 1.0)
        scale = doc.get("manufacturing_checks", {}).get("length_scale") or {}
        record = {"case": case, "run": label}
        for key in ("eroded", "intermediate", "dilated"):
            record[f"{key}_compliance_J"] = designs[key]["compliance_J"]
            record[f"{key}_volume_fraction"] = designs[key]["volume_fraction"]
        record["solid_min_size_cells"] = (scale.get("solid_min_size_m", float("nan")) / cell
                                          if scale else None)
        record["void_min_size_cells"] = (scale.get("void_min_size_m", float("nan")) / cell
                                         if scale else None)
        frames.append(record)
        blueprint = record["intermediate_compliance_J"]
        rows.append([
            case, label,
            f"{_fmt(record['eroded_compliance_J'], 6)} "
            f"({100.0 * (record['eroded_compliance_J'] / blueprint - 1.0):+.1f} %)",
            _fmt(blueprint, 6),
            f"{_fmt(record['dilated_compliance_J'], 6)} "
            f"({100.0 * (record['dilated_compliance_J'] / blueprint - 1.0):+.1f} %)",
            f"{_fmt(record['eroded_volume_fraction'], 4)} / "
            f"{_fmt(record['intermediate_volume_fraction'], 4)} / "
            f"{_fmt(record['dilated_volume_fraction'], 4)}",
            _fmt(record["solid_min_size_cells"], 3), _fmt(record["void_min_size_cells"], 3),
        ])
    if not rows:
        return None
    pd.DataFrame(frames).to_csv(os.path.join(_OUTPUT, "robust.csv"), index=False)
    return _markdown_table(
        ["case", "run", "eroded compliance [J]", "blueprint [J]", "dilated [J]",
         "volume fractions e / b / d", "smallest member [cells]", "narrowest gap [cells]"],
        rows)


def overhang_table(results_dir: str) -> Optional[str]:
    """Overhang-filtered designs beside the unfiltered runs of the same decks."""
    rows, frames = [], []
    for case, label in OVERHANG_CASES:
        path = os.path.join(results_dir, case, "summary.json")
        if not os.path.isfile(path):
            continue
        doc = load_json(path)
        result = doc.get("optimization_result", {})
        block = doc.get("manufacturing_checks", {}).get("overhang")
        if not block:
            continue
        interp = doc.get("solid_interpretation", {})
        record = {
            "case": case, "run": label,
            "build_direction": block.get("build_direction"),
            "compliance_J": result.get("compliance_J"),
            "unsupported_elements": block.get("unsupported_elements"),
            "solid_elements": block.get("solid_elements"),
            "unsupported_fraction": block.get("unsupported_fraction"),
            "connected_groups": interp.get("connected_groups_above_threshold"),
            "iterations": result.get("iterations"),
            "stop_reason": result.get("stop_reason"),
        }
        frames.append(record)
        rows.append([case, label, str(record["build_direction"]),
                     _fmt(record["compliance_J"], 6),
                     f"{record['unsupported_elements']} of {record['solid_elements']}",
                     f"{100.0 * float(record['unsupported_fraction']):.2f} %",
                     _fmt(record["connected_groups"]),
                     f"{record['iterations']}, {record['stop_reason']}"])
    if not rows:
        return None
    pd.DataFrame(frames).to_csv(os.path.join(_OUTPUT, "overhang.csv"), index=False)
    return _markdown_table(
        ["case", "run", "build", "compliance [J]", "unsupported / solid elements",
         "unsupported volume", "groups at 0.5", "iterations, stop"], rows)


def tet10_study_table(results_dir: str) -> Optional[str]:
    path = os.path.join(results_dir, "tet10_part_study", "tet10_part_study.csv")
    if not os.path.isfile(path):
        return None
    table = pd.read_csv(path)
    table.to_csv(os.path.join(_OUTPUT, "tet10_part_study.csv"), index=False)
    rows = []
    for _, r in table.sort_values(["variant", "size_mm"], ascending=[True, False]).iterrows():
        rows.append([r["variant"], _fmt(r["size_mm"]), _fmt(int(r["num_elements"])),
                     _fmt(int(r["num_dofs"])), _fmt(r["compliance_vertical_J"], 6),
                     f"{100.0 * r['rel_diff_vertical_vs_finest_tet10']:+.2f} %",
                     _fmt(r["compliance_lateral_J"], 6),
                     f"{100.0 * r['rel_diff_lateral_vs_finest_tet10']:+.2f} %",
                     f"{r['volume_error']:+.1e}", _fmt(r["wall_seconds"], 3)])
    return _markdown_table(
        ["elements", "size [mm]", "cells", "DOFs", "vertical C [J]", "vs finest Tet10",
         "lateral C [J]", "vs finest Tet10", "volume error", "solve [s]"], rows)


def tet10_verification_table(results_dir: str) -> Optional[str]:
    rows = []
    for name, stem, column in (
            ("cantilever tip vs Timoshenko", "mesh_convergence_tet10",
             "rel_error_timoshenko[-]"),
            ("column buckling vs Euler-Engesser", "buckling_euler", "rel_error_engesser[-]")):
        path = os.path.join(results_dir, "verification", f"{stem}.csv")
        if not os.path.isfile(path):
            continue
        table = pd.read_csv(path)
        for _, r in table.iterrows():
            grid = " x ".join(str(int(r[k])) for k in ("nx", "ny", "nz")
                              if k in r and int(r[k]) > 0)
            rows.append([name, str(r["element"]), grid, _fmt(int(r["num_dofs"])),
                         f"{100.0 * float(r[column]):.4g} %"])
    if not rows:
        return None
    return _markdown_table(["study", "element", "grid", "DOFs", "error"], rows)


def load_verification_table(results_dir: str) -> Optional[str]:
    """The pressure, volume and thermal load studies: one row per element of
    each quarter-section study (its finest mesh and the order measured between
    its two finest meshes), then the bimetal strip and the hanging bar."""
    base = os.path.join(results_dir, "verification")
    rows = []
    frames = []
    labels = {"lame_cylinder": "Lame cylinder, internal pressure",
              "rotating_disk": "rotating disk / cylinder",
              "thermal_cylinder": "conduction + thermal stress"}
    for stem, study in labels.items():
        path = os.path.join(base, f"{stem}.csv")
        if not os.path.isfile(path):
            continue
        table = pd.read_csv(path)
        for element, sub in table.groupby("element", sort=False):
            last = sub.iloc[-1]
            has_t = "T_rms_error[-]" in sub.columns
            row = {
                "study": study, "element": str(element),
                "mesh": f"{int(last['n_r'])} x {int(last['n_theta'])}",
                "DOFs": int(last["num_dofs"]),
                "u RMS error": float(last["u_rms_error[-]"]),
                "u order": float(last["u_rms_order[-]"]),
                "stress error": float(last["stress_error[-]"]),
                "T RMS error": float(last["T_rms_error[-]"]) if has_t else None,
                "T order": float(last["T_rms_order[-]"]) if has_t else None,
            }
            frames.append(row)
            rows.append([study, ELEMENT_LABELS.get(row["element"], row["element"]), row["mesh"],
                         _fmt(row["DOFs"]), _fmt(row["u RMS error"], 3),
                         _fmt(row["u order"], 3), _fmt(row["stress error"], 3),
                         _fmt(row["T RMS error"], 3) if has_t else "-",
                         _fmt(row["T order"], 3) if has_t else "-"])
    path = os.path.join(base, "bimetal_strip.csv")
    if os.path.isfile(path):
        table = pd.read_csv(path)
        last = table.iloc[-1]
        frames.append({"study": "bimetal strip curvature", "element": "Quad4",
                       "mesh": f"{int(last['nx'])} x {int(last['ny'])}",
                       "DOFs": int(last["num_dofs"]), "u RMS error": float(last["error[-]"]),
                       "u order": float(last["order[-]"])})
        rows.append(["bimetal strip (curvature error)", "Q4",
                     f"{int(last['nx'])} x {int(last['ny'])}", _fmt(int(last["num_dofs"])),
                     _fmt(float(last["error[-]"]), 3), _fmt(float(last["order[-]"]), 3),
                     "-", "-", "-"])
    path = os.path.join(base, "self_weight.csv")
    if os.path.isfile(path):
        table = pd.read_csv(path)
        for element, sub in table.groupby("element", sort=False):
            last = sub.iloc[-1]
            grid = " x ".join(str(int(last[k])) for k in ("nx", "ny", "nz") if int(last[k]) > 0)
            order = last["u_rms_order[-]"]
            frames.append({"study": "hanging bar, self-weight", "element": str(element),
                           "mesh": grid, "DOFs": int(last["num_dofs"]),
                           "u RMS error": float(last["u_rms_error[-]"]),
                           "u order": None if pd.isna(order) else float(order)})
            rows.append(["hanging bar, self-weight", ELEMENT_LABELS.get(element, element), grid,
                         _fmt(int(last["num_dofs"])), _fmt(float(last["u_rms_error[-]"]), 3),
                         "exact" if element == "Tet10" else _fmt(float(order), 3),
                         "-", "-", "-"])
    if not rows:
        return None
    pd.DataFrame(frames).to_csv(os.path.join(_OUTPUT, "loads.csv"), index=False)
    return _markdown_table(["study", "element", "finest mesh", "DOFs", "u RMS error", "order",
                            "stress error", "T RMS error", "T order"], rows)


_OUTPUT = "docs/results"


def main(argv=None) -> int:
    global _OUTPUT
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", default="results")
    parser.add_argument("--output", default="docs/results")
    args = parser.parse_args(argv)
    _OUTPUT = args.output
    os.makedirs(args.output, exist_ok=True)

    sections = [
        ("Verification and validation", verification_table(args.results),
         "Generated by `python/scripts/write_result_tables.py` from "
         "`results/verification/summary.json`. Verification compares against "
         "exact answers for the discrete problem; validation compares against an "
         "independent theory, where a finite gap is expected."),
        ("Linear simplices: cantilever convergence", simplex_convergence_table(args.results),
         "From `results/verification/mesh_convergence_simplex.csv`. The Tri3 beam is "
         "the plane-stress cantilever of the Q4 study at nu = 0.3, the Tet4 beam the "
         "solid cantilever of the Hex8 study (nu = 0), each meshed by splitting the "
         "structured cells. The error is against Timoshenko beam theory."),
        ("Multigrid CG against Cholesky and Jacobi CG", multigrid_table(args.results),
         "From `results/verification/multigrid_scaling.csv`: one solve of a "
         "cantilever block to a relative residual of 1e-10, times including the "
         "preconditioner setup or factorisation. A one-level hierarchy is the "
         "direct coarse solve (below the coarse-grid size) and converges in one "
         "iteration."),
        ("Quadratic tetrahedra and linear buckling against theory",
         tet10_verification_table(args.results),
         "From `results/verification/mesh_convergence_tet10.csv` and "
         "`buckling_euler.csv`: the error of the tip deflection of a solid "
         "cantilever against Timoshenko beam theory, and of the first buckling load "
         "factor of a clamped column against Euler with Engesser's shear correction, "
         "on the same grids for every element."),
        ("Pressure, volume and thermal loads against exact solutions",
         load_verification_table(args.results),
         "From the `lame_cylinder`, `rotating_disk`, `thermal_cylinder`, "
         "`bimetal_strip` and `self_weight` CSVs of `results/verification`. Errors "
         "are RMS nodal errors relative to the RMS exact field (displacement; "
         "temperature change for T), on the finest mesh, and the order is measured "
         "between the two finest meshes; the element's order is 2 (3 for the Tet10). "
         "The stress error is the largest error of sigma_r, sigma_theta, sigma_z at "
         "element centroids over the largest exact stress. The Hex8 and Tet10 "
         "sections are one cell deep with u_z = 0 (plane strain); the rotating Q4 and "
         "Tri3 rows are a plane-stress disk, the Hex8 and Tet10 rows a plane-strain "
         "cylinder. The bimetal row is the error of the curvature against "
         "Timoshenko's formula; the Tet10 hanging bar is exact."),
        ("Cross-validation against independent codes", cross_validation_table(args.results),
         "Generated from `results/cross_validation/summary.json` by "
         "`python/scripts/cross_validate.py`: node-by-node comparison of the "
         "nodal displacements with scikit-fem (same element formulations) and "
         "CalculiX (C3D8, C3D4 and C3D10 are the same elements; CPS4 and CPS3 are "
         "plane elements CalculiX expands into a layer of solids, which matches "
         "plane stress only at nu = 0, so those rows at nu != 0 are INFO: recorded, "
         "not judged; the plane-strain expansion CPE4 is exact). CalculiX results "
         "are read from its .frd output, which carries six significant digits, so "
         "differences below 5e-6 relative are its rounding. The `buckling` rows "
         "compare load factors mode by mode: scikit-fem assembles the geometric "
         "stiffness of the same discrete problem, CalculiX runs *BUCKLE with its own "
         "stress-stiffness evaluation. For the decks of the pressure, volume and "
         "thermal loads, `scikit-fem loads` integrates those loads itself from the "
         "exported deck, CalculiX integrates them from its own load cards, and "
         "`calculix conduction` compares the temperatures of CalculiX's "
         "*HEAT TRANSFER solution (relative to the temperature range). Where "
         "CalculiX's own formulation differs from SparLab's - the element-average "
         "temperature of a C3D8, the four-point centrifugal and three-point face "
         "pressure rules of a C3D10 - the max rel diff shown is CalculiX against "
         "scikit-fem solving CalculiX's problem, the judged number; SparLab's "
         "difference to CalculiX is in cross_validation.csv "
         "(`sparlab_vs_calculix`)."),
        ("Benchmark results", benchmark_table(args.results),
         "Generated from each `results/<case>/summary.json`. `stiffness gain` is "
         "the compliance of an equal-mass uniform plate divided by the optimised "
         "compliance, so values above 1 mean the optimised design is stiffer at "
         "the same mass; it is a plane-stress quantity and is n/a for the solid "
         "cases. The L-bracket's plate would fill its passive void quadrant, so "
         "its gain is not a fair comparison (`docs/benchmarks.md`, section 5)."),
        ("Linear buckling: the constrained column", buckling_table(args.results),
         "`SIMP lambda_1` is the lowest load factor of the SIMP model the "
         "constraint acts on (the eroded design in a robust run), `mode 1 energy "
         "in solid` the share of that mode's strain energy in elements at density "
         ">= 0.5, and `exported part lambda_1` the lowest load factor of the "
         "thresholded part re-analysed as solid material - the number that says "
         "whether the part meets the requirement of 6."),
        ("Robust formulation and erosion check", robust_table(args.results),
         "Eroded (eta + 0.1 for the MBB beam, + 0.05 for the column), blueprint "
         "(eta = 0.5) and dilated (eta - delta) designs of the final filtered "
         "field. The robust runs minimised the eroded compliance; the plain run "
         "was evaluated at the same thresholds once at the end "
         "(`projection.erosion_check`). Member and gap sizes are measured on the "
         "thresholded blueprint by morphological opening and closing, to half a "
         "cell."),
        ("Overhang filter", overhang_table(args.results),
         "A solid element (density >= 0.5) off the build plate is unsupported when "
         "nothing solid lies directly or diagonally below it (3 elements in 2-D, a "
         "cross of 5 in 3-D): a 45-degree overhang limit on these square cells. "
         "The unfiltered rows are the same decks run with `--no-overhang-filter`."),
        ("Tet4 against Tet10 on the engine mount", tet10_study_table(args.results),
         "From `results/tet10_part_study` (`python/scripts/tet10_part_study.py`): "
         "the engine mount meshed by Gmsh at each size, solved on linear "
         "tetrahedra, on the same mesh elevated to straight-sided Tet10, and on "
         "curved Tet10 cells. Loads are uniform tractions on the pin hole with 12 "
         "kN (vertical) and 5 kN (lateral) resultants on the exact cylinder. The "
         "reference is the finest curved Tet10 run, itself a lower bound on the "
         "exact compliance."),
        ("Heaviside projection", projection_table(args.results),
         "Every projected run and, where one exists, the unprojected run of the "
         "same problem. `grey` is 4 mean(rho (1 - rho)) of the physical density; "
         "`thresholded / SIMP` re-solves the design cut at 0.5 as solid material "
         "and divides by the compliance the optimiser minimised, so 1.0 means the "
         "objective is what the structure delivers."),
        ("Real geometry (meshes read from files)", real_geometry_table(args.results),
         "Meshes generated by `python/scripts/make_meshes.py` with Gmsh and read by "
         "SparLab's Gmsh and Abaqus/CalculiX readers. Quality is 4 sqrt(3) A / "
         "sum of squared edge lengths for triangles and 6 sqrt(2) V / (rms edge "
         "length)^3 for tetrahedra: 1 for an equilateral cell, 0 for a degenerate "
         "one."),
        ("Stress-constrained optimisation", stress_table(args.results),
         "The constraint bounds the relaxed (rho^q) von Mises stress aggregated "
         "with a scaled p-norm; `max relaxed ratio` is the true relaxed maximum "
         "over the limit at the returned design, and `re-solve max vM` is the "
         "peak von Mises of the thresholded structure analysed as solid material, "
         "which is the number that says whether the structure meets the limit. "
         "The `constraint off` row is the same deck run with the constraint "
         "disabled."),
        ("Modal results", modal_table(args.results),
         "`max residual` is the largest relative eigenpair residual "
         "||K phi - lambda M phi|| / ||lambda M phi|| over the reported modes."),
        ("Runtime scaling (2-D, Q4)", scaling_table(args.results),
         "Minimum over repeats at each size, which is the standard estimator for "
         "wall-clock benchmarks."),
        ("Runtime scaling (3-D, Hex8)", scaling_table(args.results, "runtime_scaling_3d"),
         "Same protocol on a Hex8 block (nx x nx/2 x nx/4 cells). The direct "
         "solver's fill-in grows much faster in 3-D, which is what limits it on "
         "solid problems; the solver comparisons below show the multigrid solver "
         "on the same blocks."),
        ("Linear solvers, Q4 plate", solver_comparison_table(args.results, "Quad4"),
         "One solve from scratch at each size: Cholesky is factorisation plus one "
         "back-substitution, multigrid (MG) and Jacobi are CG to a relative "
         "residual of 1e-10 including the preconditioner setup. `entries/DOF` "
         "counts what the solver stores beyond K itself (factor L and D; coarse "
         "operators, prolongators and the dense coarsest factor). The long-format "
         "table is `solver_scaling.csv` beside this file."),
        ("Linear solvers, Hex8 block", solver_comparison_table(args.results, "Hex8"),
         "Same protocol on the Hex8 block of the 3-D scaling benchmark."),
        ("Linear solvers, Tet4 block", solver_comparison_table(args.results, "Tet4"),
         "Same protocol with each hexahedral cell split into six tetrahedra."),
        ("Design study", study_table(os.path.join(args.results, "study")),
         "One row per arm of the aerospace parametric study; the per-point table "
         "is `aerospace_study.csv` beside this file and the interpretation is in "
         "`docs/aerospace_study.md`. A worst group count above 1 means that arm "
         "contains a point whose thresholded density field is not a single "
         "connected structure."),
    ]

    lines = [
        "# SparLab result tables",
        "",
        "**Generated file - do not edit by hand.** Refresh with `make results` "
        "after a benchmark run.",
        "",
    ]
    missing = []
    for title, table, note in sections:
        lines.append(f"## {title}")
        lines.append("")
        if table is None:
            missing.append(title)
            lines.append("_Not available: the corresponding run has not been made._")
        else:
            lines.append(note)
            lines.append("")
            lines.append(table)
        lines.append("")

    path = os.path.join(args.output, "README.md")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    print(f"  wrote {path}")
    for name in missing:
        print(f"  missing section: {name}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
