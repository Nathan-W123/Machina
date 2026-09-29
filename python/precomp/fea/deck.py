"""The sparlab_form input deck.

`build_deck` writes, into one directory,

* `deck.json` - a SparLab JSON deck (mesh, material, model,
  boundary_conditions and nonlinear blocks as `sparlab_solve` reads them,
  `docs/configuration.md`) plus a `forming` object with the tool and the
  ordered steps;
* `toolpath.csv` - the tool-centre trajectory, columns `t,x,y,z`, t strictly
  increasing in [0, 1];
* `commanded.npz`, `precomp_deck.json` - the commanded surface and the
  setup, for provenance. The solver does not read them and they are not part
  of the content hash.

The `forming` keys follow the contract in `docs/precomp.md` ("The forming
deck"). They are written in one place, `forming_block`, so a change of the
C++ contract is a change here only.

Geometry: the blank is a structured mesh of the square [-L/2, L/2]^2 by
[-t, 0], L the blank size rounded to whole elements; the clamp fixes every node (through the thickness) whose x or y lies
within `clamp_margin` of the blank edge; the tool contacts the top faces
(z = 0) inside the clamp; the 3-2-1 release support fixes x, y, z of the top
node nearest (-a, -a, 0), y and z of the one nearest (a, -a, 0) and z of the
one nearest (-a, a, 0), with a = L/2 - clamp_margin / 2 (three flange points
in the former clamp): six constraints, no rigid-body mode left and no
constraint on the springback.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from .. import __version__
from .._util import (PathLike, PrecompError, canonical_json, read_json, sha256_bytes,
                     to_jsonable, write_json)
from ..geometry.heightmap import HeightMap
from ..toolpath import AIR, Toolpath, contour_toolpath, spiral_toolpath
from .setup import FormingSetup

DECK_FILE = "deck.json"
TOOLPATH_FILE = "toolpath.csv"
COMMANDED_FILE = "commanded.npz"
PROVENANCE_FILE = "precomp_deck.json"
TOOL_NAME = "tool"


def make_toolpath(setup: FormingSetup, commanded: HeightMap) -> Toolpath:
    """The tool path the setup prescribes for `commanded` (style, step-down,
    spacing, direction and tool radius from the setup)."""
    kwargs = dict(direction=setup.toolpath_direction)
    if setup.toolpath_style == "spiral":
        return spiral_toolpath(commanded, setup.tool_radius, setup.step_down,
                               setup.toolpath_spacing, **kwargs)
    return contour_toolpath(commanded, setup.tool_radius, setup.step_down,
                            setup.toolpath_spacing, **kwargs)


def _clamp_region(setup: FormingSetup) -> Dict[str, Any]:
    e = setup.free_half_width
    return {"name": "clamp", "any_of": [
        {"box": {"xmax": -e}}, {"box": {"xmin": e}},
        {"box": {"ymax": -e}}, {"box": {"ymin": e}}]}


def clamp_condition(setup: FormingSetup) -> Dict[str, Any]:
    """The clamped frame: x, y, z fixed at zero."""
    return {"name": "clamp", "fix": ["x", "y", "z"], "value": [0.0, 0.0, 0.0],
            "region": _clamp_region(setup)}


def support_321(setup: FormingSetup) -> List[Dict[str, Any]]:
    """The 3-2-1 release support (see the module docstring)."""
    a = 0.5 * setup.meshed_blank_size - 0.5 * setup.clamp_margin
    return [
        {"name": "support_xyz", "fix": ["x", "y", "z"], "value": [0.0, 0.0, 0.0],
         "region": {"nearest_node": [-a, -a, 0.0]}},
        {"name": "support_yz", "fix": ["y", "z"], "value": [0.0, 0.0, 0.0],
         "region": {"nearest_node": [a, -a, 0.0]}},
        {"name": "support_z", "fix": ["z"], "value": [0.0, 0.0, 0.0],
         "region": {"nearest_node": [-a, a, 0.0]}},
    ]


def forming_block(setup: FormingSetup, toolpath_file: str = TOOLPATH_FILE) -> Dict[str, Any]:
    """The deck's `forming` object: one spherical tool and the ordered steps.

    Steps: "form" (the tool follows the trajectory, clamp on), "unload" (no
    tool, clamp on) and, with release "321", "release" (no tool, 3-2-1
    support instead of the clamp).
    """
    e = setup.free_half_width
    surface = {"name": "tool_side", "box": {"zmin": 0.0, "xmin": -e, "xmax": e,
                                            "ymin": -e, "ymax": e}}
    tool = {"name": TOOL_NAME, "shape": "sphere", "radius": setup.tool_radius,
            "surface": surface, "friction": setup.friction,
            "contact": dict(setup.contact), "trajectory": {"file": toolpath_file}}
    clamp = clamp_condition(setup)
    steps: List[Dict[str, Any]] = [
        {"name": "form", "type": "form", "tools": [TOOL_NAME], "boundary_conditions": [clamp],
         "increments": {"max_tool_travel": setup.max_tool_travel}},
        {"name": "unload", "type": "release", "tools": [], "boundary_conditions": [clamp]},
    ]
    if setup.release == "321":
        steps.append({"name": "release", "type": "release", "tools": [],
                      "boundary_conditions": support_321(setup)})
    return {"tools": [tool], "steps": steps}


def deck_document(setup: FormingSetup) -> Dict[str, Any]:
    """The full deck.json content for `setup` (independent of the commanded
    shape, which enters through the trajectory file only)."""
    n = setup.elements_per_side
    L = setup.meshed_blank_size
    mesh_type = "structured_hex" if setup.element == "hex8" else "structured_tet"
    nonlinear = {"enabled": True, "kinematics": setup.kinematics}
    nonlinear.update(setup.solver)
    return {
        "name": setup.name,
        "description": (f"Single-point incremental forming of a {setup.material.name} blank, "
                        f"written by precomp {__version__}"),
        "units": "SI",
        "mesh": {"type": mesh_type, "nx": n, "ny": n, "nz": setup.layers,
                 "lx": L, "ly": L, "lz": setup.thickness,
                 "x0": -0.5 * L, "y0": -0.5 * L, "z0": -setup.thickness},
        "material": setup.material.to_sparlab(),
        "model": {"stress_state": "three_dimensional"},
        "boundary_conditions": [clamp_condition(setup)],
        "nonlinear": nonlinear,
        "forming": forming_block(setup),
        "output": {"csv": True},
    }


def check_toolpath(setup: FormingSetup, path: Toolpath) -> None:
    """Raise PrecompError if the tool would reach the clamped frame: every
    in-contact tool centre must stay a tool radius inside the free window."""
    contact = path.level != AIR
    xy = np.abs(path.points[contact, :2])
    reach = float(xy.max()) + setup.tool_radius if len(xy) else 0.0
    if reach > setup.free_half_width:
        raise PrecompError(
            f"the tool path reaches {reach:.4g} m from the centre (tool radius included), "
            f"beyond the unclamped half-width {setup.free_half_width:.4g} m; enlarge the "
            "blank or reduce the clamp margin")


def build_deck(setup: FormingSetup, commanded: HeightMap, out_dir: PathLike,
               toolpath: Optional[Toolpath] = None) -> Path:
    """Write the deck for forming `commanded` with `setup` into `out_dir`.

    The tool path is generated from `commanded` with the setup's path
    settings unless given. Returns the directory. Raises PrecompError when
    the path would reach the clamp.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = toolpath if toolpath is not None else make_toolpath(setup, commanded)
    if abs(path.tool_radius - setup.tool_radius) > 1e-12:
        raise PrecompError("the tool path was made for a different tool radius")
    check_toolpath(setup, path)
    write_json(out / DECK_FILE, deck_document(setup))
    path.to_sparlab_csv(out / TOOLPATH_FILE)
    commanded.save(out / COMMANDED_FILE)
    write_json(out / PROVENANCE_FILE, {
        "precomp_version": __version__,
        "setup": setup.to_dict(),
        "commanded_grid": commanded.grid.to_dict(),
        "commanded_metadata": to_jsonable(commanded.metadata),
        "toolpath": path.summary(),
        "toolpath_metadata": to_jsonable(path.metadata),
    })
    return out


def deck_hash(deck_dir: PathLike, solver_version: str) -> str:
    """Content hash of a deck: SHA-256 of the canonical deck JSON, the
    trajectory CSV and the solver's `--version` text. Two decks with the same
    hash are the same run."""
    d = Path(deck_dir)
    try:
        deck = canonical_json(read_json(d / DECK_FILE)).encode("utf-8")
        trajectory = (d / TOOLPATH_FILE).read_bytes()
    except FileNotFoundError as exc:
        raise PrecompError(f"{d} is not a complete deck directory: {exc}") from exc
    return sha256_bytes(deck, trajectory, solver_version.strip().encode("utf-8"))
