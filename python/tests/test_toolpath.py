"""Drop cutter and tool paths: no gouging, contact, closed loops, monotone spiral."""

import numpy as np
import pandas as pd
import pytest

from precomp.geometry import Freeform, Grid, HeightMap, TruncatedCone, zeros
from precomp.toolpath import (AIR, Toolpath, contour_toolpath, dsif_support_path, max_gouge,
                              spherical_structure, spiral_toolpath, tool_center_surface)

R = 0.005


def test_drop_cutter_on_a_plane_is_exactly_R_above_it():
    g = Grid.centered(0.05, 3e-4)
    for level in (0.0, -0.0123):
        cz = tool_center_surface(zeros(g).with_z(np.full(g.shape, level)), R)
        assert np.array_equal(cz.z, np.full(g.shape, level + R))


def test_spherical_structure_is_the_ball_cap():
    fp, st = spherical_structure(R, 1e-3)
    assert fp.shape == (11, 11) and fp[5, 5] and fp[5, 0] and not fp[0, 0]
    assert st[5, 5] == R and st[5, 0] == pytest.approx(0.0, abs=1e-12)


def test_drop_cutter_never_gouges_a_cone_and_touches_it_everywhere():
    """For tool centres over the walls, fillets and floor of a cone, the
    distance from the centre to the exact (analytic) surface is never below
    R - tol and equals R at the grid-node contact: the ball touches without
    gouging. tol bounds what a grid of spacing h can hide between nodes,
    about h^2 / (8 R cos^3 theta) = 6 um here (h = 0.25 mm, theta = 50 deg)."""
    cone = TruncatedCone(0.03, 50.0, 0.015, 0.004, 0.004)
    g = Grid.centered(0.1, 2.5e-4)
    cz = tool_center_surface(cone.heightmap(g), R)
    X, Y = g.mesh()
    rng = np.random.default_rng(3)
    cand = np.flatnonzero((np.hypot(X, Y) < cone.footprint_radius() + R).ravel())
    pick = rng.choice(cand, 120, replace=False)
    hf = g.h / 4.0                                  # dense offsets include the nodes
    k = np.arange(-int(1.05 * R / hf), int(1.05 * R / hf) + 1) * hf
    DX, DY = np.meshgrid(k, k, indexing="xy")
    tol = 1e-5
    dmin = []
    for idx in pick:
        c = np.array([X.ravel()[idx], Y.ravel()[idx], cz.z.ravel()[idx]])
        sx, sy = c[0] + DX, c[1] + DY
        sz = cone.height(sx, sy)
        dmin.append(np.sqrt((sx - c[0]) ** 2 + (sy - c[1]) ** 2 + (sz - c[2]) ** 2).min())
    dmin = np.array(dmin)
    assert dmin.min() >= R - tol                     # never penetrates
    assert np.abs(dmin - R).max() <= tol              # touches at every centre
    assert np.sum(np.abs(dmin - R) < 1e-12) > 10      # exactly R where the contact is a node


def _loops(path: Toolpath):
    """The closed ring of each level, in order: the in-contact run of a level
    is its link from the previous level followed by the ring, which starts at
    the first occurrence of its closing point."""
    out = []
    lv = path.level
    start = 0
    for i in range(1, len(lv) + 1):
        if i == len(lv) or lv[i] != lv[start]:
            if lv[start] != AIR:
                run = path.points[start:i]
                j = int(np.flatnonzero(np.all(run == run[-1], axis=1))[0])
                out.append((lv[start], run[j:]))
            start = i
    return out


def test_contour_toolpath_levels_are_closed_equally_spaced_ccw_loops(small_cone):
    g = Grid.centered(0.1, 5e-4)
    target = small_cone.heightmap(g)
    step, spacing = 1e-3, 1e-3
    path = contour_toolpath(target, R, step, spacing)
    heights = path.metadata["level_heights"]
    assert heights[:5] == pytest.approx([R - (k + 1) * step for k in range(5)])
    cz = tool_center_surface(target, R)
    assert heights[-1] == pytest.approx(cz.z.min() + 1e-3 * step)
    assert max_gouge(path, cz) == 0.0 and path.metadata["max_gouge_m"] == 0.0
    loops = [pts for lv, pts in _loops(path) if len(pts) > 10]
    assert len(loops) == len(heights)
    for ring in loops:
        closed = ring[-1] - ring[0]
        assert np.linalg.norm(closed) < 1e-12                          # closed loop
        seg = np.linalg.norm(np.diff(ring[:, :2], axis=0), axis=1)
        assert seg.max() <= spacing * (1 + 1e-9) and seg.min() > 0.5 * spacing
        x, y = ring[:-1, 0], ring[:-1, 1]
        area = 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)
        assert area > 0                                               # counter-clockwise
    # consecutive loops start close to each other (aligned starts)
    starts = np.array([pts[0, :2] for pts in loops])
    assert np.linalg.norm(np.diff(starts, axis=0), axis=1).max() < 3 * step / np.tan(0.5) + spacing
    # alternate direction flips every other level
    alt = contour_toolpath(target, R, step, spacing, alternate=True)
    signs = []
    for pts in [p for lv, p in _loops(alt) if len(p) > 10][:4]:
        x, y = pts[:-1, 0], pts[:-1, 1]
        signs.append(np.sign(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)))
    assert signs == [1, -1, 1, -1]


def test_spiral_descends_monotonically_and_touches_the_surface(small_cone):
    g = Grid.centered(0.1, 5e-4)
    target = small_cone.heightmap(g)
    path = spiral_toolpath(target, R, 1e-3, 1e-3)
    contact = path.level != AIR
    z = path.points[contact, 2]
    assert np.all(np.diff(z) <= 0.0)                                   # never climbs
    assert z[0] == R                                                   # tip on the sheet
    assert z[-1] == pytest.approx(path.metadata["level_heights"][-1])
    cz = tool_center_surface(target, R)
    gap = cz.interpolate(path.points[contact, 0], path.points[contact, 1]) - z
    assert gap.max() <= 1e-12                                          # never gouges
    first_rev = path.level[contact] == 1
    assert np.abs(gap[~first_rev]).max() < 1e-9                        # touches after rev. 1
    assert path.points[0, 2] == pytest.approx(R + 2e-3)                # approach from above
    assert path.points[-1, 2] == pytest.approx(R + 2e-3)               # and retract


def test_two_pockets_give_two_loops_per_level_and_air_moves():
    part = Freeform(((-0.03, 0.0, 0.012, 0.008), (0.03, 0.0, 0.012, 0.008)), 0.05, 0.06)
    g = Grid.centered(0.14, 5e-4)
    target = part.heightmap(g)
    path = contour_toolpath(target, R, 1e-3, 1e-3)
    assert max(path.metadata["loops_per_level"]) == 2
    air = np.flatnonzero(path.level == AIR)
    assert len(air) > 2                                                # traverses between pockets
    assert np.all(path.points[air, 2] == pytest.approx(R + 2e-3))
    with pytest.raises(ValueError, match="one loop per level"):
        spiral_toolpath(target, R, 1e-3, 1e-3)


def test_csv_exports_round_trip(tmp_path, small_cone):
    target = small_cone.heightmap(Grid.centered(0.1, 5e-4))
    path = spiral_toolpath(target, R, 1e-3, 1e-3)
    csv = path.to_sparlab_csv(tmp_path / "toolpath.csv")
    frame = pd.read_csv(csv)
    assert list(frame.columns) == ["t", "x", "y", "z"]
    t = frame["t"].to_numpy()
    assert t[0] == 0.0 and t[-1] == 1.0 and np.all(np.diff(t) > 0)
    back = Toolpath.from_sparlab_csv(csv, R)
    assert np.array_equal(back.points, path.points)
    robot = pd.read_csv(path.to_robot_csv(tmp_path / "robot.csv", feed=0.02))
    assert list(robot.columns) == ["x", "y", "z", "i", "j", "k", "feed"]
    assert np.all(robot["k"] == 1.0) and np.all(robot["feed"] == 0.02)
    summ = path.summary(feed=0.02)
    assert summ["time_s"] == pytest.approx(path.total_length / 0.02)
    assert np.allclose(path.at(np.array([0.0, 1.0])), path.points[[0, -1]])


def test_dsif_support_tool_sits_below_the_sheet(small_cone):
    target = small_cone.heightmap(Grid.centered(0.1, 5e-4))
    path = contour_toolpath(target, R, 1e-3, 2e-3)
    t = 1e-3
    support = dsif_support_path(path, target, thickness=t, support_radius=0.004)
    d = np.linalg.norm(support.points - path.points, axis=1)
    assert np.allclose(d, R + t + 0.004)
    assert support.metadata["experimental"] is True
    assert np.all(support.points[:, 2] < path.points[:, 2])


def test_flat_target_is_refused():
    with pytest.raises(Exception, match="flat"):
        contour_toolpath(zeros(Grid.centered(0.05, 1e-3)), R, 1e-3)
    with pytest.raises(ValueError):
        Toolpath(np.zeros((1, 3)), np.zeros(1), R)
    with pytest.raises(ValueError):
        HeightMap(Grid.centered(0.05, 1e-3), np.zeros((2, 2)))
