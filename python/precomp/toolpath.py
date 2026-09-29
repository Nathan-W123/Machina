"""Tool paths for single-point incremental forming with a ball tool.

Everything here is the path of the tool *centre* [m], which is what SparLab's
forming trajectory and a robot's tool-centre point take.

Drop cutter. `tool_center_surface` computes, on the grid of a height map, the
height at which a ball of radius R centred above each node first touches the
surface when lowered:

    c_z(x) = max over |d| <= R of [ z(x + d) + sqrt(R^2 - |d|^2) ],

as a grey-scale dilation with a non-flat spherical structuring element
(`scipy.ndimage.grey_dilation`). It is exact over the grid nodes - the ball
touches at least one node and no node lies inside it. Between nodes the
surface is only known through interpolation; a bilinearly interpolated surface
can rise above the node samples by up to about ``h^2 / (8 R cos^3 theta)``
at a wall of angle theta, which is the gouge the grid can hide (0.5 um for
h = 0.25 mm, R = 5 mm, theta = 45 deg). Outside the grid the surface is
taken as continuing flat at the border height (border mode 'nearest').

Paths. `contour_toolpath` follows z-level loops: the tool tip descends by
`step_down` per level, i.e. level k has the tool centre at the constant height
``R - k step_down``, and the loop is the level set of the drop-cutter surface
at that height, so every point touches the target without gouging.
`spiral_toolpath` descends continuously: along each revolution the height
falls linearly with arc length from one level to the next, and the x-y
position is found on the segment between corresponding points of the two
loops where the drop-cutter surface has exactly that height (bisection), so
the spiral is gouge-free as well and its height never increases.
Loops come from `contourpy`, are closed, oriented consistently, start near
the previous level's start and are resampled by arc length.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import ndimage

from ._util import PathLike, PrecompError, require_positive
from .geometry.heightmap import HeightMap

#: Level value of points where the tool moves in the air (approach, retract,
#: traverse between separate loops).
AIR = -1


# ---------------------------------------------------------------------------
# Drop cutter
# ---------------------------------------------------------------------------
def spherical_structure(tool_radius: float, spacing: float) -> Tuple[np.ndarray, np.ndarray]:
    """(footprint, structure) of a ball of radius `tool_radius` on a grid of
    spacing `spacing` [m]: footprint is |d| <= R, structure sqrt(R^2 - |d|^2)."""
    R = require_positive("tool_radius", tool_radius)
    m = int(math.floor(R / spacing + 1e-12))
    k = np.arange(-m, m + 1) * spacing
    DX, DY = np.meshgrid(k, k, indexing="xy")
    r2 = DX * DX + DY * DY
    footprint = r2 <= R * R * (1 + 1e-12)
    structure = np.where(footprint, np.sqrt(np.maximum(R * R - r2, 0.0)), 0.0)
    return footprint, structure


def tool_center_surface(hm: HeightMap, tool_radius: float) -> HeightMap:
    """Drop-cutter tool-centre height c_z(x, y) [m] of a ball tool over `hm`.

    The tool centre is (x, y, c_z); the ball touches the surface (exact over
    the grid nodes, see the module docstring). On a flat region c_z = z + R.
    The result carries metadata {"quantity": "tool_centre_height",
    "tool_radius": R}. Invalid nodes of `hm` are treated as surface at their
    fill height.
    """
    footprint, structure = spherical_structure(tool_radius, hm.grid.h)
    cz = ndimage.grey_dilation(hm.z, footprint=footprint, structure=structure, mode="nearest")
    return HeightMap(hm.grid, cz, hm.mask, {"quantity": "tool_centre_height",
                                            "tool_radius": float(tool_radius)})


# ---------------------------------------------------------------------------
# Toolpath container
# ---------------------------------------------------------------------------
@dataclass
class Toolpath:
    """An ordered tool-centre path.

    points : (n, 3) tool-centre positions [m].
    level : (n,) int, the z level a point belongs to (1 = first level below
        the sheet), `AIR` (-1) where the tool moves in the air.
    tool_radius : ball radius [m].
    metadata : style, step-down, spacing, level heights, ...

    Consecutive duplicate points are removed on construction, so the arc
    length `s` and the pseudo-time `t = s / total_length` are strictly
    increasing.
    """

    points: np.ndarray
    level: np.ndarray
    tool_radius: float
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        p = np.asarray(self.points, dtype=float)
        lv = np.asarray(self.level, dtype=int)
        if p.ndim != 2 or p.shape[1] != 3 or len(p) < 2:
            raise ValueError("a tool path needs an (n >= 2, 3) array of points")
        if lv.shape != (len(p),):
            raise ValueError("level must hold one entry per point")
        if not np.all(np.isfinite(p)):
            raise ValueError("tool path points must be finite")
        step = np.linalg.norm(np.diff(p, axis=0), axis=1)
        keep = np.concatenate([[True], step > 1e-12])
        self.points = p[keep]
        self.level = lv[keep]
        if len(self.points) < 2:
            raise ValueError("a tool path needs at least two distinct points")
        self.tool_radius = float(self.tool_radius)
        self.metadata = dict(self.metadata)

    @property
    def s(self) -> np.ndarray:
        """(n,) cumulative 3-D arc length of the tool centre [m], s[0] = 0."""
        return np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(self.points, axis=0),
                                                               axis=1))])

    @property
    def total_length(self) -> float:
        """Path length [m]."""
        return float(self.s[-1])

    @property
    def t(self) -> np.ndarray:
        """(n,) pseudo-time s / total_length in [0, 1], strictly increasing."""
        s = self.s
        return s / s[-1]

    @property
    def num_levels(self) -> int:
        return int(len(np.unique(self.level[self.level != AIR])))

    def summary(self, feed: Optional[float] = None) -> Dict[str, Any]:
        """Length [m], points, levels, height range [m], longest segment [m],
        in-contact length [m] and, with `feed` [m/s], the time at that feed [s]."""
        p = self.points
        seg = np.linalg.norm(np.diff(p, axis=0), axis=1)
        contact = (self.level[1:] != AIR) & (self.level[:-1] != AIR)
        out = {
            "length_m": self.total_length,
            "contact_length_m": float(seg[contact].sum()),
            "points": int(len(p)),
            "levels": self.num_levels,
            "z_min_m": float(p[:, 2].min()),
            "z_max_m": float(p[:, 2].max()),
            "max_segment_m": float(seg.max()),
            "tool_radius_m": self.tool_radius,
        }
        if feed is not None:
            out["feed_m_per_s"] = require_positive("feed", feed)
            out["time_s"] = self.total_length / feed
        return out

    def to_sparlab_csv(self, path: PathLike) -> Path:
        """Write the SparLab trajectory table `t,x,y,z` (t in [0, 1], metres).

        Values are written with 17 significant digits, so the file round-trips
        exactly and its content hash depends only on the path.
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        t = self.t
        lines = ["t,x,y,z"]
        for ti, (x, y, z) in zip(t, self.points):
            lines.append(f"{ti:.17g},{x:.17g},{y:.17g},{z:.17g}")
        path.write_text("\n".join(lines) + "\n", encoding="ascii")
        return path

    def to_robot_csv(self, path: PathLike, feed: float,
                     axis: Sequence[float] = (0.0, 0.0, 1.0)) -> Path:
        """Write robot waypoints `x,y,z,i,j,k,feed`.

        x, y, z : tool-centre point [m]; i, j, k : unit tool-axis direction
        (default +z, a vertical tool; the axis points from the ball towards
        the spindle); feed : path speed [m/s]. Strict SI like the rest of the
        package: a controller that expects mm and mm/min needs the conversion
        on its side.
        """
        feed = require_positive("feed", feed)
        a = np.asarray(axis, dtype=float)
        if a.shape != (3,) or not np.linalg.norm(a) > 0:
            raise ValueError("axis must be a non-zero 3-vector")
        a = a / np.linalg.norm(a)
        frame = pd.DataFrame({"x": self.points[:, 0], "y": self.points[:, 1],
                              "z": self.points[:, 2], "i": a[0], "j": a[1], "k": a[2],
                              "feed": feed})
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(path, index=False, float_format="%.12g")
        return path

    @classmethod
    def from_sparlab_csv(cls, path: PathLike, tool_radius: float) -> "Toolpath":
        """Read a `t,x,y,z` trajectory (levels are unknown and set to 0)."""
        path = Path(path)
        if not path.is_file():
            raise PrecompError(f"{path} does not exist")
        frame = pd.read_csv(path, float_precision="round_trip")
        if list(frame.columns) != ["t", "x", "y", "z"]:
            raise PrecompError(f"{path}: expected the columns t,x,y,z, found "
                               f"{','.join(frame.columns)}")
        t = frame["t"].to_numpy()
        if np.any(np.diff(t) <= 0) or t[0] < 0 or t[-1] > 1 + 1e-12:
            raise PrecompError(f"{path}: t must increase strictly within [0, 1]")
        return cls(frame[["x", "y", "z"]].to_numpy(), np.zeros(len(frame), dtype=int),
                   tool_radius, {"source": str(path)})

    def at(self, t: np.ndarray) -> np.ndarray:
        """Tool-centre positions (k, 3) [m] at pseudo-times `t` (linear)."""
        tt = self.t
        t = np.asarray(t, dtype=float)
        return np.stack([np.interp(t, tt, self.points[:, i]) for i in range(3)], axis=-1)


# ---------------------------------------------------------------------------
# Loops
# ---------------------------------------------------------------------------
def _signed_area(loop: np.ndarray) -> float:
    x, y = loop[:, 0], loop[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def _loop_length(loop: np.ndarray) -> float:
    closed = np.vstack([loop, loop[:1]])
    return float(np.linalg.norm(np.diff(closed, axis=0), axis=1).sum())


def _nearest_on_loop(loop: np.ndarray, point: np.ndarray) -> Tuple[int, float]:
    """(segment index, fraction) of the point of the closed polyline nearest `point`."""
    a = loop
    b = np.roll(loop, -1, axis=0)
    ab = b - a
    denom = np.maximum((ab * ab).sum(axis=1), 1e-300)
    u = np.clip(((point - a) * ab).sum(axis=1) / denom, 0.0, 1.0)
    proj = a + u[:, None] * ab
    k = int(np.argmin(((proj - point) ** 2).sum(axis=1)))
    return k, float(u[k])


def _rotate_to(loop: np.ndarray, seg: int, frac: float) -> np.ndarray:
    """The closed polyline started at the point `frac` along segment `seg`."""
    a = loop[seg]
    b = loop[(seg + 1) % len(loop)]
    start = a + frac * (b - a)
    rest = np.roll(loop, -(seg + 1), axis=0)
    out = np.vstack([start, rest])
    step = np.linalg.norm(np.diff(np.vstack([out, out[:1]]), axis=0), axis=1)
    return out[np.concatenate([[True], step[:-1] > 1e-15])]


def resample_loop(loop: np.ndarray, spacing: float, n: Optional[int] = None) -> np.ndarray:
    """Resample a closed polyline (m, 2) [m] at equal arc-length steps.

    Returns `n` points (default ceil(length / spacing), at least 3) starting
    at the loop's first point; the loop stays implicitly closed (the last
    point is not repeated).
    """
    closed = np.vstack([loop, loop[:1]])
    seg = np.linalg.norm(np.diff(closed, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    total = s[-1]
    if n is None:
        n = max(3, int(math.ceil(total / spacing - 1e-9)))
    snew = np.arange(n) * (total / n)
    return np.stack([np.interp(snew, s, closed[:, 0]), np.interp(snew, s, closed[:, 1])], -1)


def extract_loops(cz: HeightMap, level: float, min_loop_length: float) -> List[np.ndarray]:
    """Closed contour loops (m, 2) [m] of the height map `cz` at height `level`.

    Uses `contourpy` (marching squares with linear interpolation along cell
    edges). A contour that ends on the grid border raises PrecompError: the
    part runs off the grid and the loop cannot be closed. Loops shorter than
    `min_loop_length` [m] are dropped (the tool cannot follow them).
    """
    from contourpy import LineType, contour_generator

    g = cz.grid
    gen = contour_generator(x=g.x, y=g.y, z=cz.z, line_type=LineType.Separate)
    loops = []
    for line in gen.lines(level):
        line = np.asarray(line, dtype=float)
        if len(line) < 3:
            continue
        if not np.allclose(line[0], line[-1], rtol=0.0, atol=1e-9 * g.h):
            raise PrecompError(
                f"the contour at height {level:.6g} m reaches the grid border; the grid "
                "must contain the part with a margin of at least the tool radius")
        loop = line[:-1]
        if _loop_length(loop) >= min_loop_length:
            loops.append(loop)
    return loops


def _centre_levels(cz: HeightMap, tool_radius: float, step_down: float) -> List[float]:
    """Tool-centre heights of the levels: R - k step_down above the deepest
    point, then one last level step_down / 1000 above it."""
    cmin = float(cz.z.min())
    delta = 1e-3 * step_down
    levels = []
    k = 1
    while tool_radius - k * step_down > cmin + delta:
        levels.append(tool_radius - k * step_down)
        k += 1
    if not levels or levels[-1] - (cmin + delta) > delta:
        levels.append(cmin + delta)
    if tool_radius - levels[0] < 0.5 * delta:
        raise PrecompError("the target is flat: there is nothing to form")
    return levels


def _ordered_loops(cz: HeightMap, levels: Sequence[float], spacing: float, direction: str,
                   alternate: bool, start_angle_deg: float,
                   min_loop_length: float) -> List[Tuple[float, List[np.ndarray]]]:
    """(height, loops) of every level that has a loop: the loops oriented,
    ordered by nearest neighbour and started at the point nearest the current
    tool position; not yet resampled."""
    if direction not in ("ccw", "cw"):
        raise ValueError("direction must be 'ccw' or 'cw'")
    out: List[Tuple[float, List[np.ndarray]]] = []
    position: Optional[np.ndarray] = None
    for k, lev in enumerate(levels):
        loops = extract_loops(cz, lev, min_loop_length)
        if not loops:
            continue
        ccw = (direction == "ccw") != (alternate and k % 2 == 1)
        oriented = []
        for loop in loops:
            if (_signed_area(loop) > 0) != ccw:
                loop = loop[::-1]
            oriented.append(loop)
        if position is None:
            a = math.radians(start_angle_deg)
            d = np.array([math.cos(a), math.sin(a)])
            first = max(oriented, key=lambda lp: float((lp @ d).max()))
            position = first[int(np.argmax(first @ d))]
        ordered = []
        remaining = list(oriented)
        while remaining:
            dists = []
            for lp in remaining:
                seg, frac = _nearest_on_loop(lp, position)
                a0 = lp[seg]
                b0 = lp[(seg + 1) % len(lp)]
                dists.append(float(np.linalg.norm(a0 + frac * (b0 - a0) - position)))
            j = int(np.argmin(dists))
            lp = remaining.pop(j)
            seg, frac = _nearest_on_loop(lp, position)
            lp = _rotate_to(lp, seg, frac)
            ordered.append(lp)
            position = lp[0]
        out.append((float(lev), ordered))
    if not out:
        raise PrecompError("no tool path loop found; check the step-down and the target depth")
    return out


def _link(cz: HeightMap, a: np.ndarray, b: np.ndarray, spacing: float,
          safe_z: float) -> Tuple[np.ndarray, bool]:
    """Points strictly between tool positions a and b (3-vectors) [m].

    A straight move whose height never has to exceed the higher end to stay
    above the drop-cutter surface is used directly (its height is lifted to
    that surface where needed); otherwise the tool retracts to `safe_z`,
    traverses and plunges. Returns (points, is_air_move).
    """
    dist = float(np.linalg.norm(b[:2] - a[:2]))
    n = max(1, int(math.ceil(dist / spacing)))
    u = (np.arange(1, n) / n)[:, None]
    xy = a[:2] + u * (b[:2] - a[:2])
    zl = a[2] + u[:, 0] * (b[2] - a[2])
    if len(xy):
        cs = cz.interpolate(xy[:, 0], xy[:, 1], masked=False)
        z = np.maximum(zl, cs)
        if np.all(z <= max(a[2], b[2]) + 1e-12):
            return np.column_stack([xy, z]), False
    else:
        return np.empty((0, 3)), False
    up = np.array([a[0], a[1], safe_z])
    over = np.array([b[0], b[1], safe_z])
    return np.vstack([up, over]), True


def _lift_to_surface(cz: HeightMap, pts: np.ndarray) -> np.ndarray:
    """Raise points below the (bilinear) drop-cutter surface onto it."""
    cs = cz.interpolate(pts[:, 0], pts[:, 1], masked=False)
    out = pts.copy()
    out[:, 2] = np.maximum(out[:, 2], cs)
    return out


def _min_loop(target: HeightMap, spacing: float, given: Optional[float]) -> float:
    if given is not None:
        return float(given)
    return max(2.0 * spacing, 4.0 * math.pi * target.grid.h)


def contour_toolpath(target: HeightMap, tool_radius: float, step_down: float,
                     spacing: float = 1e-3, *, direction: str = "ccw",
                     alternate: bool = False, start_angle_deg: float = 0.0,
                     clearance: float = 2e-3,
                     min_loop_length: Optional[float] = None) -> Toolpath:
    """Z-level contour tool path over the target surface.

    Parameters
    ----------
    target : the tool-side target surface [m].
    tool_radius : ball radius R [m].
    step_down : vertical increment of the tool tip between levels [m].
    spacing : arc-length spacing of the points along a loop [m] (the actual
        spacing is the loop length over ceil(length / spacing)).
    direction : 'ccw' or 'cw' seen from +z (from the tool side).
    alternate : reverse the direction on every other level.
    start_angle_deg : the first loop starts at its extreme point in this
        direction [deg from +x]; later loops start at their point nearest to
        the previous loop's start.
    clearance : height of the tool tip above the sheet plane for the
        approach, the final retract and traverses between separate loops [m].
    min_loop_length : loops shorter than this are skipped [m] (default the
        larger of 2 spacings and the circumference of a circle of radius two
        grid spacings, below which a loop is not resolved by the grid - e.g.
        the last loop around the pole of a dome).

    Level k (k = 1, 2, ...) holds the tool centre at height R - k step_down
    (the tip at -k step_down); the last level lies step_down / 1000 above the
    deepest point, so the deepest loop traces the boundary of a flat floor
    (the floor itself is not swept, as is usual in SPIF). Moves between loops
    descend directly when unobstructed, otherwise through a retract.
    """
    R = require_positive("tool_radius", tool_radius)
    step_down = require_positive("step_down", step_down)
    spacing = require_positive("spacing", spacing)
    cz = tool_center_surface(target, R)
    levels = _centre_levels(cz, R, step_down)
    per_level = _ordered_loops(cz, levels, spacing, direction, alternate, start_angle_deg,
                               _min_loop(target, spacing, min_loop_length))
    safe_z = R + clearance
    pts: List[np.ndarray] = []
    lvl: List[np.ndarray] = []
    current: Optional[np.ndarray] = None
    for k, (height, loops) in enumerate(per_level, start=1):
        for loop in loops:
            xy = resample_loop(loop, spacing)
            ring = np.column_stack([np.vstack([xy, xy[:1]]),
                                    np.full(len(xy) + 1, height)])
            ring = _lift_to_surface(cz, ring)
            if current is None:
                start = np.array([ring[0, 0], ring[0, 1], safe_z])
                pts.append(start[None])
                lvl.append(np.array([AIR]))
                current = start
            link, air = _link(cz, current, ring[0], spacing, safe_z)
            if len(link):
                pts.append(link)
                lvl.append(np.full(len(link), AIR if air else k))
            pts.append(ring)
            lvl.append(np.full(len(ring), k))
            current = ring[-1]
    end = np.array([current[0], current[1], safe_z])
    pts.append(end[None])
    lvl.append(np.array([AIR]))
    meta = {"style": "contour", "step_down": step_down, "spacing": spacing,
            "direction": direction, "alternate": bool(alternate),
            "level_heights": [h for h, _ in per_level],
            "loops_per_level": [len(l) for _, l in per_level], "clearance": clearance}
    path = Toolpath(np.vstack(pts), np.concatenate(lvl), R, meta)
    path.metadata["max_gouge_m"] = max_gouge(path, cz)
    return path


def spiral_toolpath(target: HeightMap, tool_radius: float, step_down: float,
                    spacing: float = 1e-3, *, direction: str = "ccw",
                    start_angle_deg: float = 0.0, clearance: float = 2e-3,
                    final_loop: bool = True,
                    min_loop_length: Optional[float] = None) -> Toolpath:
    """Continuous (helical) tool path: the height falls linearly with arc
    length through each revolution, one step_down per revolution.

    Revolution k runs from loop k - 1 to loop k (loop 0 is loop 1 at the
    sheet plane, where the tool tip touches z = 0). At normalised arc length
    lambda the height is ``c_(k-1) + lambda (c_k - c_(k-1))`` and the x-y
    position is found on the line from the point of loop k - 1 at lambda
    towards its nearest point on loop k, where the drop-cutter surface has
    exactly that height (bisection on its bilinear interpolant), so the tool
    touches the target all along the revolution and the height never
    increases. With `final_loop` a last
    revolution at the deepest level closes the path. Parameters as for
    `contour_toolpath`; every level must have exactly one loop (one pocket
    without islands), otherwise ValueError.
    """
    R = require_positive("tool_radius", tool_radius)
    step_down = require_positive("step_down", step_down)
    spacing = require_positive("spacing", spacing)
    cz = tool_center_surface(target, R)
    levels = _centre_levels(cz, R, step_down)
    per_level = _ordered_loops(cz, levels, spacing, direction, False, start_angle_deg,
                               _min_loop(target, spacing, min_loop_length))
    for k, (_, loops) in enumerate(per_level, start=1):
        if len(loops) != 1:
            raise ValueError(f"spiral_toolpath needs one loop per level; level {k} has "
                             f"{len(loops)}. Use contour_toolpath for this part.")
    levels = [h for h, _ in per_level]
    loops = [l[0] for _, l in per_level]
    heights = [R] + list(levels)
    loops = [loops[0]] + loops
    pts = []
    lvl = []
    for k in range(1, len(loops)):
        la, lb = loops[k - 1], loops[k]
        n = max(3, int(math.ceil(max(_loop_length(la), _loop_length(lb)) / spacing)))
        A = resample_loop(la, spacing, n)
        B = A.copy() if la is lb else _nearest_on_polyline(lb, A, spacing)
        lam = np.arange(n) / n
        z = heights[k - 1] + lam * (heights[k] - heights[k - 1])
        xy = _blend_to_height(cz, A, B, z)
        pts.append(np.column_stack([xy, z]))
        lvl.append(np.full(n, k))
    if final_loop:
        la, lb = loops[-2], loops[-1]
        n = max(3, int(math.ceil(max(_loop_length(la), _loop_length(lb)) / spacing)))
        A = resample_loop(la, spacing, n)
        B = A.copy() if la is lb else _nearest_on_polyline(lb, A, spacing)
        z = np.full(n + 1, heights[-1])
        xy = _blend_to_height(cz, np.vstack([A, A[:1]]), np.vstack([B, B[:1]]), z)
        pts.append(np.column_stack([xy, z]))
        lvl.append(np.full(n + 1, len(loops) - 1))
    body = np.vstack(pts)
    safe_z = R + clearance
    start = np.array([[body[0, 0], body[0, 1], safe_z]])
    end = np.array([[body[-1, 0], body[-1, 1], safe_z]])
    meta = {"style": "spiral", "step_down": step_down, "spacing": spacing,
            "direction": direction, "level_heights": [float(v) for v in levels],
            "clearance": clearance, "final_loop": bool(final_loop)}
    path = Toolpath(np.vstack([start, body, end]),
                    np.concatenate([[AIR], np.concatenate(lvl), [AIR]]), R, meta)
    path.metadata["max_gouge_m"] = max_gouge(path, cz)
    return path


def max_gouge(path: Toolpath, cz: HeightMap) -> float:
    """Largest depth [m] of a path point below the bilinear drop-cutter
    surface `cz` (0 when no point gouges); points off the grid are ignored."""
    cs = cz.interpolate(path.points[:, 0], path.points[:, 1], masked=False)
    d = cs - path.points[:, 2]
    d = d[np.isfinite(d)]
    return float(max(0.0, d.max())) if len(d) else 0.0


def _nearest_on_polyline(loop: np.ndarray, pts: np.ndarray, spacing: float) -> np.ndarray:
    """For each point (k, 2), the nearest point of the closed polyline `loop`
    (to within spacing / 20: the loop is densified and searched with a KD-tree)."""
    from scipy.spatial import cKDTree

    dense = resample_loop(loop, spacing / 20.0)
    _, idx = cKDTree(dense).query(pts)
    return dense[idx]


def _first_crossing(cz: HeightMap, A: np.ndarray, D: np.ndarray, z: np.ndarray,
                    mu: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """(found, mu_root): per row, the first sign change of c_z(A + mu D) - z
    from >= 0 to < 0 over the sampled `mu` (k,), bisected 50 times, keeping
    the side where c_z < z."""
    xy = A[:, None, :] + mu[None, :, None] * D[:, None, :]
    f = cz.interpolate(xy[..., 0], xy[..., 1], masked=False) - z[:, None]
    f = np.where(np.isfinite(f), f, -np.inf)
    change = (f[:, :-1] >= 0) & (f[:, 1:] < 0)
    found = change.any(axis=1) & (np.linalg.norm(D, axis=1) > 0)
    idx = np.argmax(change, axis=1)
    lo = mu[idx].astype(float)
    hi = mu[idx + 1].astype(float)
    rows = np.flatnonzero(found)
    a, d, zz = A[rows], D[rows], z[rows]
    l, h = lo[rows], hi[rows]
    for _ in range(50):
        mid = 0.5 * (l + h)
        p = a + mid[:, None] * d
        up = cz.interpolate(p[:, 0], p[:, 1], masked=False) - zz >= 0
        l = np.where(up, mid, l)
        h = np.where(up, h, mid)
    root = np.zeros(len(A))
    root[rows] = h
    return found, root


def _blend_to_height(cz: HeightMap, A: np.ndarray, B: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Per row, an x-y point where the bilinear drop-cutter surface equals z.

    The point is searched along the line from A through B: first from a
    quarter of |AB| behind A to six times |AB| beyond it on a grid of |AB| / 8,
    then, for rows without a crossing, within three grid spacings of B on a
    grid of h / 8 (a narrow crossing, such as the last small loop around the
    pole of a dome). The first sign change of c_z - z from >= 0 to < 0 is
    bisected, keeping the side where c_z < z, so the tool never ends inside
    the surface. Where A and B coincide (the first revolution, or a single
    level) the direction is the downhill direction of c_z at A (smoothed over
    two grid spacings), over four grid spacings. A row with no crossing keeps
    A: there the surface stays below z (the tool is in the air), or, if not,
    the point shows up in `max_gouge`.
    """
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float).copy()
    h = cz.grid.h
    same = np.linalg.norm(B - A, axis=1) < 1e-12
    if np.any(same):
        gx, gy = cz.gradient(sigma=2.0 * h)
        grad = np.column_stack([cz.sample(gx, A[same, 0], A[same, 1]),
                                cz.sample(gy, A[same, 0], A[same, 1])])
        norm = np.linalg.norm(grad, axis=1, keepdims=True)
        down = np.divide(-grad, norm, out=np.zeros_like(grad), where=norm > 0)
        B[same] = A[same] + 4.0 * h * down
    D = B - A
    found, root = _first_crossing(cz, A, D, z, np.arange(-2, 49) / 8.0)
    out = A.copy()
    out[found] = A[found] + root[found, None] * D[found]
    missing = np.flatnonzero(~found & (np.linalg.norm(D, axis=1) > 0))
    if len(missing):
        length = np.linalg.norm(D[missing], axis=1)
        unit = D[missing] / length[:, None]
        start = B[missing] - 3.0 * h * unit
        f2, r2 = _first_crossing(cz, start, unit, z[missing], np.arange(49) * (h / 8.0))
        rows = missing[f2]
        out[rows] = start[f2] + r2[f2, None] * unit[f2]
    return out


def dsif_support_path(primary: Toolpath, target: HeightMap, thickness: float,
                      support_radius: Optional[float] = None) -> Toolpath:
    """EXPERIMENTAL: the supporting (bottom) tool of double-sided incremental forming.

    Each point of the primary path is offset against the local surface normal
    by ``R1 + thickness + R2`` (R1 the primary, R2 the support radius; equal
    radii give the ``thickness + 2R`` of the design), so the support ball
    touches the underside of a sheet of that thickness opposite the primary
    tool. The normal is that of the drop-cutter surface at the primary
    centre's x-y, which equals the target's normal at the contact point for a
    smooth target. Thinning is ignored (the sheet is taken at its initial
    thickness); in the air the support simply stays below the primary tool.
    """
    t = require_positive("thickness", thickness)
    R1 = primary.tool_radius
    R2 = R1 if support_radius is None else require_positive("support_radius", support_radius)
    cz = tool_center_surface(target, R1)
    n = cz.sample(cz.normals(), primary.points[:, 0], primary.points[:, 1])
    n = np.where(np.isfinite(n), n, np.array([0.0, 0.0, 1.0]))
    n[primary.level == AIR] = np.array([0.0, 0.0, 1.0])
    n /= np.linalg.norm(n, axis=1, keepdims=True)
    pts = primary.points - (R1 + t + R2) * n
    meta = dict(primary.metadata)
    meta.update({"style": f"{primary.metadata.get('style', 'path')}_dsif_support",
                 "experimental": True, "thickness": t, "primary_radius": R1})
    return Toolpath(pts, primary.level.copy(), R2, meta)
