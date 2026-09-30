"""Rim support in the forming deck: a backing plate, or a DSIF support tool.

A single tool pushing from above cannot keep the sheet between the clamp and
the first contour from bending down: the released part comes out too deep
near its rim (`benchmarks/springback`). `FormingSetup.support` adds support
from below, with the existing rigid tools of `sparlab_form`
(docs/forming.md):

* ``"backing_plate"`` - a rigid plate under the sheet whose opening is the
  part's outline grown by `clearance` (partial support die). It is a
  `plane` tool (normal +z, through the sheet's bottom z = -t), held still,
  whose slave surface is the bottom faces lying wholly outside the opening:
  their nodes are listed by id (`node_ids`, the structured mesh's documented
  numbering, include/sparlab/mesh/StructuredMesh.hpp), since no geometric
  region primitive describes an arbitrary outline. The sheet outside the
  opening cannot go below the plate; nothing holds it down onto the plate.
  On the mesh the opening is the union of the elements it touches, up to
  one element wider than asked: `plate_nodes` reports the clearance the
  mesh realises. The plate is active in "form" and removed with the tool
  in "unload" (clamp still on), before the 3-2-1 release.
* ``"dsif"`` - double-sided incremental forming: a second ball (the
  support) under the sheet, on the bottom faces of the free window,
  following `precomp.toolpath.dsif_support_points` - opposite the forming
  tool along the surface normal, never above the sheet not yet formed -
  written with the forming tool's pseudo-time, so both tools move together
  in one "form" step; "unload" removes both. With `rim_pass`, after that
  unload the support alone runs `precomp.toolpath.rim_pass_path` (step
  "rim_pass", pseudo-time 2 to 3, clamp on) and is removed in "rim_unload",
  before the release.

What a command above the sheet plane means (`command_upper_bound`): the
forming tool presses from above and its path is made for the part of the
command below the plane (`precomp.fea.deck.make_toolpath`); a backing
plate at the sheet's bottom holds the sheet up to the plane outside the
opening but pushes nothing higher; the DSIF support pushes up only
opposite the forming tool, whose path never rises above the plane. Only
the rim pass pushes the sheet up by itself: in the band it sweeps (signed
distance to the target's outline in [-rim_inside, rim_outside], the flange
strip included, which displacement adjustment then adjusts instead of
holding - `compensation_masks`) the command may rise above the plane, by at
most `rim_max_raise` and no further than the support ball can push the
commanded underside from below (`command_upper_bound`: a raised band
narrower than the ball is out of its reach). Every other command is kept
at z <= 0, and `build_deck` refuses one above the bound.

All lengths are metres.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.spatial import cKDTree

from .._util import PrecompError
from ..geometry.heightmap import HeightMap
from ..toolpath import (AIR, Toolpath, dsif_support_points, reach_from_below, rim_pass_path,
                        signed_outline_distance)
from .setup import CONTACT_KEYS, FormingSetup

#: Names of the support tools in the deck (fields of tool_forces.csv).
PLATE_TOOL = "plate"
SUPPORT_TOOL = "support"
#: The support tool's trajectory file, beside toolpath.csv.
SUPPORT_PATH_FILE = "support_path.csv"
#: The rim pass runs over this pseudo-time window, after "unload" (which
#: ends at most at 2: the forming tool's path spans [0, 1]).
RIM_PASS_WINDOW = (2.0, 3.0)
#: Commands within this of the bound are accepted [m].
BOUND_TOLERANCE = 1e-9


@dataclass
class SupportPlan:
    """What a support strategy adds to the deck (`plan_support`).

    tools : tool objects after the forming tool; form_tools : their names
    active in "form"; rim_pass_end : the pseudo-time the "rim_pass" step
    ends at (None: no rim pass); files : {name: text} written beside the
    deck (trajectory files - part of the content hash); info : provenance.
    """

    tools: List[Dict[str, Any]] = field(default_factory=list)
    form_tools: List[str] = field(default_factory=list)
    rim_pass_end: Optional[float] = None
    files: Dict[str, str] = field(default_factory=dict)
    info: Dict[str, Any] = field(default_factory=dict)


def _contact(setup: FormingSetup, settings: Dict[str, Any]) -> Dict[str, Any]:
    out = {"friction": float(settings["friction"])}
    out.update({k: setup.contact[k] for k in CONTACT_KEYS if k in setup.contact})
    out["penalty"] = float(settings["penalty"])
    return out


def bottom_node_grid(setup: FormingSetup) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(ids, x, y), each (n + 1, n + 1) [row j, column i], of the bottom
    nodes (z = -t) of the structured blank mesh: node (i, j, 0) has id
    j (n + 1) + i and sits at (-L/2 + i dx, -L/2 + j dx), n elements per
    side, dx = L / n (include/sparlab/mesh/StructuredMesh.hpp; the Tet4 mesh
    splits the same cells and keeps the nodes)."""
    n = setup.elements_per_side
    L = setup.meshed_blank_size
    c = -0.5 * L + L * np.arange(n + 1) / n
    X, Y = np.meshgrid(c, c, indexing="xy")
    ids = np.arange((n + 1) * (n + 1)).reshape(n + 1, n + 1)
    return ids, X, Y


def part_distance(reference: HeightMap, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Distance [m] from points (x, y) to the nearest part node of
    `reference` (nodes deeper than 1 um): 0 on the part, to within the grid
    spacing of the distance to its outline outside it."""
    part = reference.mask & (reference.z < -1e-6)
    if not part.any():
        raise PrecompError("the reference surface has no part: a backing plate needs an "
                           "outline to follow")
    X, Y = reference.grid.mesh()
    tree = cKDTree(np.column_stack([X[part], Y[part]]))
    d, _ = tree.query(np.column_stack([np.ravel(x), np.ravel(y)]))
    return d.reshape(np.shape(x))


def plate_nodes(setup: FormingSetup, reference: HeightMap,
                clearance: Optional[float] = None) -> Tuple[List[int], Dict[str, Any]]:
    """The bottom nodes the backing plate carries, and a report.

    A bottom face (cell) is carried when its four corners lie farther than
    `clearance` (default the setup's) from the part of `reference`; the
    nodes are the corners of the carried cells (so every listed node is a
    slave node of a carried face and sparlab_form's face selection - every
    node of a face in the region - picks exactly those faces). The report
    gives the clearance asked, the one the mesh realises (the smallest
    distance of a carried node to the part), and the counts.
    """
    s = setup.resolved_support()
    c = float(s["clearance"] if clearance is None else clearance)
    ids, X, Y = bottom_node_grid(setup)
    d = part_distance(reference, X, Y)
    out = d > c
    cell = out[:-1, :-1] & out[1:, :-1] & out[:-1, 1:] & out[1:, 1:]
    carried = np.zeros_like(out)
    carried[:-1, :-1] |= cell
    carried[1:, :-1] |= cell
    carried[:-1, 1:] |= cell
    carried[1:, 1:] |= cell
    if not cell.any():
        raise PrecompError("the backing plate carries no face: the opening covers the blank")
    info = {"clearance_m": c, "realised_clearance_m": float(d[carried].min()),
            "carried_faces": int(cell.sum()), "carried_nodes": int(carried.sum()),
            "bottom_faces": int(cell.size)}
    return sorted(int(i) for i in ids[carried]), info


def plate_tool(setup: FormingSetup, reference: HeightMap) -> Tuple[Dict[str, Any],
                                                                   Dict[str, Any]]:
    """The backing plate's tool object (a still plane through the sheet's
    bottom, normal +z) and its report (`plate_nodes`)."""
    s = setup.resolved_support()
    nodes, info = plate_nodes(setup, reference)
    t = float(setup.thickness)
    tool: Dict[str, Any] = {"name": PLATE_TOOL, "shape": "plane", "normal": [0.0, 0.0, 1.0],
                            "surface": {"name": "plate_side", "node_ids": nodes}}
    tool.update(_contact(setup, s))
    tool["trajectory"] = {"times": [0.0, 1.0], "points": [[0.0, 0.0, -t], [0.0, 0.0, -t]]}
    return tool, info


def support_tool(setup: FormingSetup) -> Dict[str, Any]:
    """The DSIF support tool object: a ball on the bottom faces of the free
    window, driven by `SUPPORT_PATH_FILE`."""
    s = setup.resolved_support()
    e = setup.free_half_width
    t = float(setup.thickness)
    tool: Dict[str, Any] = {"name": SUPPORT_TOOL, "shape": "sphere", "radius": float(s["radius"]),
                            "surface": {"name": "support_side",
                                        "box": {"zmax": -t, "xmin": -e, "xmax": e,
                                                "ymin": -e, "ymax": e}}}
    tool.update(_contact(setup, s))
    tool["trajectory"] = {"file": SUPPORT_PATH_FILE}
    return tool


def _check_reach(setup: FormingSetup, centres: np.ndarray, normals: np.ndarray,
                 radius: float, what: str) -> None:
    """PrecompError unless every contact point (centre + radius * normal)
    lies inside the free window."""
    if not len(centres):
        return
    tip = centres[:, :2] + radius * normals[:, :2]
    reach = float(np.abs(tip).max())
    if reach > setup.free_half_width + 1e-12:
        raise PrecompError(f"the {what} reaches {reach:.4g} m from the centre, beyond the "
                           f"unclamped half-width {setup.free_half_width:.4g} m; enlarge the "
                           "blank or reduce the clamp margin")


def _csv(t: np.ndarray, pts: np.ndarray) -> str:
    lines = ["t,x,y,z"]
    for ti, (x, y, z) in zip(t, pts):
        lines.append(f"{ti:.17g},{x:.17g},{y:.17g},{z:.17g}")
    return "\n".join(lines) + "\n"


def dsif_trajectory(setup: FormingSetup, forming_surface: HeightMap, path: Toolpath,
                    commanded: Optional[HeightMap] = None,
                    reference: Optional[HeightMap] = None
                    ) -> Tuple[np.ndarray, np.ndarray, Optional[float], Dict[str, Any]]:
    """(t, points, rim_pass_end, info): the support tool's trajectory knots.

    On [0, 1] the knots of the forming tool's `path` (its `t`) with
    `dsif_support_points` over `forming_surface` (the surface the path was
    made for); with `rim_pass`, the knots of `rim_pass_path` over
    `commanded` (its full command, above the plane included) and the
    outline of `reference`, on `RIM_PASS_WINDOW` by arc length, the pass
    ending (`rim_pass_end`) at its last point in contact."""
    s = setup.resolved_support()
    R2 = float(s["radius"])
    t0 = float(setup.thickness)
    pts = dsif_support_points(path, forming_surface, t0, R2, squeeze=float(s["squeeze"]),
                              thickness_law=s["thickness_law"], clearance=float(s["clearance"]))
    contact = path.level != AIR
    n = path.points - pts
    norm = np.linalg.norm(n, axis=1, keepdims=True)
    n = np.divide(n, norm, out=np.zeros_like(n), where=norm > 0)
    _check_reach(setup, pts[contact], n[contact], R2, "support tool")
    info: Dict[str, Any] = {"support_radius_m": R2, "squeeze": float(s["squeeze"]),
                            "thickness_law": s["thickness_law"],
                            "support_z_min_m": float(pts[:, 2].min()),
                            # how far the ball's body reaches in plan while in contact
                            # (a lower clamp frame starts at the free half-width)
                            "support_body_reach_m": float(np.abs(pts[contact, :2]).max() + R2)
                            if contact.any() else None,
                            "support_z_max_in_contact_m":
                                float(pts[contact, 2].max()) if contact.any() else None}
    t = path.t
    if not s["rim_pass"]:
        return t, pts, None, info
    if commanded is None or reference is None:
        raise PrecompError("the rim pass needs the commanded surface and the reference part")
    rim = rim_pass_path(commanded, reference, R2, t0, inside=float(s["rim_inside"]),
                        outside=float(s["rim_outside"]), band_spacing=float(s["rim_spacing"]),
                        spacing=float(setup.toolpath_spacing), clearance=float(s["clearance"]),
                        direction=setup.toolpath_direction)
    up = np.array([0.0, 0.0, 1.0])
    rc = rim.level != AIR
    _check_reach(setup, rim.points[rc], np.tile(up, (int(rc.sum()), 1)), R2, "rim pass")
    a, b = RIM_PASS_WINDOW
    tr = a + (b - a) * rim.t
    last = int(np.flatnonzero(rc)[-1])
    info["rim_pass"] = {"band_levels_m": rim.metadata["band_levels_m"],
                        "length_m": rim.total_length, "points": int(len(rim.points)),
                        "z_max_m": float(rim.points[rc, 2].max()),
                        "z_min_m": float(rim.points[rc, 2].min())}
    return (np.concatenate([t, tr]), np.vstack([pts, rim.points]), float(tr[last]), info)


def plan_support(setup: FormingSetup, commanded: HeightMap, forming_surface: HeightMap,
                 path: Toolpath, reference: Optional[HeightMap]) -> SupportPlan:
    """The tools, steps and files the setup's support strategy adds (see the
    module docstring). `reference` is the part the fixture is made for (the
    plate's opening and the rim pass band follow its outline); None takes
    the commanded surface."""
    plan = SupportPlan(info={"support": setup.support})
    if setup.support == "none":
        return plan
    ref = commanded if reference is None else reference
    if not ref.grid.matches(commanded.grid):
        ref = ref.resample(commanded.grid)
    plan.info["outline_from"] = "commanded" if reference is None else "target"
    if setup.support == "backing_plate":
        tool, info = plate_tool(setup, ref)
        plan.tools.append(tool)
        plan.form_tools.append(PLATE_TOOL)
        plan.info["plate"] = info
        return plan
    t, pts, rim_end, info = dsif_trajectory(setup, forming_surface, path, commanded, ref)
    plan.tools.append(support_tool(setup))
    plan.form_tools.append(SUPPORT_TOOL)
    plan.files[SUPPORT_PATH_FILE] = _csv(t, pts)
    plan.rim_pass_end = rim_end
    plan.info["dsif"] = info
    return plan


def rim_band(setup: FormingSetup, target: HeightMap) -> Tuple[np.ndarray, np.ndarray]:
    """(band, strip), (ny, nx) masks on `target`'s grid: the part within
    `rim_inside` of the outline and the flange within `rim_outside` of it -
    what the rim pass sweeps. Both empty without a DSIF rim pass."""
    band = np.zeros(target.grid.shape, dtype=bool)
    strip = band.copy()
    if setup.support != "dsif" or not setup.resolved_support()["rim_pass"]:
        return band, strip
    s = setup.resolved_support()
    sd = signed_outline_distance(target).z
    part = target.mask & (target.z < -1e-6)
    band = part & (sd >= -float(s["rim_inside"]))
    strip = target.mask & ~part & (sd <= float(s["rim_outside"]))
    return band, strip


def compensation_masks(setup: FormingSetup, target: HeightMap
                       ) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """(hold, adjust) for displacement adjustment with this setup, or
    (None, None) for the defaults (hold the flange, adjust the part). The rim
    pass pushes the flange strip it sweeps from below, so that strip is
    adjusted like the part instead of held at the sheet plane."""
    band, strip = rim_band(setup, target)
    if not strip.any():
        return None, None
    part = target.mask & (target.z < -1e-6)
    flange = target.mask & ~part
    return flange & ~strip, part | strip


def command_upper_bound(setup: FormingSetup, target: HeightMap):
    """The highest command the setup's tools can realise, as a function of
    the command: ``bound(z) -> (ny, nx)`` [m] on `target`'s grid, for
    `displacement_adjustment(upper_bound=...)` and `check_command`.

    0 (the sheet plane) everywhere for every strategy but the DSIF rim
    pass: the forming tool only pushes down, a plate at the sheet's bottom
    and the synchronised support push nothing above the plane. With the rim
    pass, on the band and strip it sweeps (`rim_band`), the command may rise
    above the plane by at most `rim_max_raise`, and only as far as the
    support ball can push its underside from below: the reach of the ball
    under the command's underside (`precomp.toolpath.reach_from_below` of
    z - t; the vertical thickness is t under the sine law) - a raised band
    narrower than the ball, or a corner seen from below tighter than it,
    cannot be pushed up by it and is cut to what the ball reaches (never
    below the plane). The cut changes the underside, so it is repeated
    until nothing changes (at most 8 times).
    """
    band, strip = rim_band(setup, target)
    cap = np.zeros(target.grid.shape)
    if not (band.any() or strip.any()):
        return lambda z: cap
    s = setup.resolved_support()
    zone = band | strip
    cap[zone] = float(s["rim_max_raise"])
    R = float(s["radius"])
    t = float(setup.thickness)

    def bound(z: np.ndarray) -> np.ndarray:
        z = np.asarray(z, dtype=float)
        out = cap.copy()
        cur = np.minimum(z, cap)
        for _ in range(8):
            if not (cur > BOUND_TOLERANCE).any():
                break
            reach = reach_from_below(target.with_z(cur - t), R).z + t
            out = np.where(zone, np.minimum(cap, np.maximum(reach, 0.0)), 0.0)
            new = np.minimum(z, out)
            if np.abs(new - cur).max() <= BOUND_TOLERANCE:
                break
            cur = new
        return out

    return bound


def check_command(setup: FormingSetup, commanded: HeightMap,
                  reference: Optional[HeightMap]) -> None:
    """PrecompError if `commanded` rises above the sheet plane where no tool
    of the setup can push it (`command_upper_bound`); a command above the
    plane needs the target the fixture is made for."""
    top = float(commanded.z[commanded.mask].max()) if commanded.mask.any() else 0.0
    if top <= BOUND_TOLERANCE:
        return
    if reference is None:
        raise PrecompError(
            f"the commanded surface rises {top * 1e3:.3f} mm above the sheet plane; only a "
            "rim pass can realise that, in a band around the target's outline - pass the "
            "target to build the deck")
    ref = reference if reference.grid.matches(commanded.grid) else \
        reference.resample(commanded.grid)
    bound = command_upper_bound(setup, ref)(commanded.z)
    over = commanded.mask & (commanded.z > bound + 1e-6)
    if over.any():
        worst = float((commanded.z - bound)[over].max())
        raise PrecompError(
            f"the commanded surface rises above what the tools can realise at {int(over.sum())} "
            f"nodes (up to {worst * 1e3:.3f} mm too high): the forming tool only pushes down, "
            f"and with support {setup.support!r} nothing pushes the sheet up there "
            "(precomp.fea.support.command_upper_bound)")


__all__ = ["PLATE_TOOL", "SUPPORT_TOOL", "SUPPORT_PATH_FILE", "RIM_PASS_WINDOW",
           "SupportPlan", "plan_support", "plate_nodes", "plate_tool", "support_tool",
           "dsif_trajectory", "command_upper_bound", "check_command", "compensation_masks",
           "rim_band", "bottom_node_grid", "part_distance"]
