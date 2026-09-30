"""Drop cutter and tool paths: no gouging, contact, closed loops, monotone spiral."""

import numpy as np
import pandas as pd
import pytest

from precomp.geometry import Freeform, Grid, HeightMap, TruncatedCone, zeros
from precomp.toolpath import (AIR, Toolpath, contour_toolpath, dsif_support_path,
                              dsif_support_points, lift_heights, max_gouge, rim_pass_path,
                              signed_outline_distance,
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


def test_dsif_support_moves_with_the_tool_opposite_it_and_never_inside_the_sheet():
    """The support ball of DSIF: one centre per point of the forming path
    (the same pseudo-time); opposite the tool along the surface normal with
    the sine-law gap where the sheet is already formed; never inside the
    in-process sheet (the first revolution starts on a flat sheet); below
    the whole part in the air."""
    t, R2 = 1e-3, 0.004
    target = TruncatedCone(0.012, 40.0, 0.006, 0.002, 0.002).heightmap(Grid.centered(0.06, 2.5e-4))
    path = spiral_toolpath(target, R, 1e-3, 1e-3)
    pts = dsif_support_points(path, target, t, R2)
    assert pts.shape == path.points.shape
    air = path.level == AIR
    contact = ~air
    # in the air: below the deepest point, clear by the clearance
    assert np.all(pts[air, 2] + R2 <= -target.depth - t - 2e-3 + 1e-12)
    # never inside the in-process sheet (the surface above the tip's height,
    # flat at it below; the underside t lower), nor the finished one (to the
    # 1 um the last level keeps above the floor) ...
    xy = pts[contact, :2].T
    now = lift_heights(target, R2, *xy, offset=t, floor=path.points[contact, 2] - R)
    assert np.all(pts[contact, 2] <= now + 1e-12)
    assert np.all(pts[contact, 2] <= lift_heights(target, R2, *xy, offset=t) + 2e-6)
    # ... nor inside the flat sheet at the first contact
    first = np.flatnonzero(contact)[0]
    assert pts[first, 2] + R2 <= -t + 1e-9
    # where the sheet is formed (the last revolutions) the balls hold the
    # sine-law thickness t cos(theta) between them, along the normal
    d = np.linalg.norm(path.points - pts, axis=1)
    cz = tool_center_surface(target, R)
    nz = cz.sample(cz.normals(), path.points[:, 0], path.points[:, 1])[:, 2]
    late = contact & (path.points[:, 2] - R < -0.004)
    assert late.sum() > 50
    gap = d[late] - R - R2
    assert np.median(np.abs(gap - t * nz[late])) < 1e-6
    assert np.all(gap >= t * nz[late] - 1e-9)          # never squeezed more than asked
    # squeeze moves it towards the tool by squeeze t_n (where it is opposite)
    sq = dsif_support_points(path, target, t, R2, squeeze=0.2)
    d2 = np.linalg.norm(path.points - sq, axis=1)
    opposite = np.abs(gap - t * nz[late]) < 1e-6
    assert opposite.sum() > 0.5 * late.sum()
    assert np.allclose((d[late] - d2[late])[opposite], 0.2 * t * nz[late][opposite], atol=2e-6)
    # the Toolpath wrapper keeps the points (its own time is not the tool's)
    wrapped = dsif_support_path(path, target, thickness=t, support_radius=R2)
    assert wrapped.metadata["synchronised_with_primary"] is False
    assert np.all(wrapped.points[:, 2] < path.points[:len(wrapped.points), 2].max())
    with pytest.raises(ValueError):
        dsif_support_points(path, target, t, R2, thickness_law="thin")


def test_rim_pass_traces_the_band_from_below_outside_in():
    t, R2 = 1e-3, 0.004
    target = TruncatedCone(0.012, 40.0, 0.004, 0.002, 0.002).heightmap(Grid.centered(0.06, 2.5e-4))
    rim = rim_pass_path(target, target, R2, t, inside=2e-3, outside=1e-3, band_spacing=1e-3,
                        spacing=1e-3)
    assert rim.metadata["band_levels_m"] == pytest.approx([1e-3, 0.0, -1e-3, -2e-3])
    assert rim.level[0] == AIR and rim.level[-1] == AIR
    on = rim.level != AIR
    # every loop point touches the commanded underside from below and never enters it
    b = lift_heights(target, R2, rim.points[on, 0], rim.points[on, 1], offset=t)
    assert np.all(rim.points[on, 2] <= b + 1e-5)
    loops = rim.level[on]
    first = rim.points[on][loops == 1]
    assert np.median(np.abs(first[:, 2] - b[loops == 1])) < 1e-6
    # outside in: the loops shrink, around the part's axis
    radii = [np.median(np.hypot(*rim.points[on][loops == k][:, :2].T)) for k in (1, 2, 3, 4)]
    assert all(a > b_ for a, b_ in zip(radii, radii[1:]))
    sd = signed_outline_distance(target)
    edge = np.median(np.abs(sd.interpolate(first[:, 0], first[:, 1]) - 1e-3))
    assert edge < 2.5e-4
    # a raised command lifts the ball with it
    raised = target.with_z(np.where(sd.z < 0, target.z + 5e-4, target.z))
    up = rim_pass_path(raised, target, R2, t, inside=2e-3, outside=1e-3, band_spacing=1e-3,
                       spacing=1e-3)
    assert up.points[up.level == 3, 2].mean() > rim.points[rim.level == 3, 2].mean() + 1e-4


def test_flat_target_is_refused():
    with pytest.raises(Exception, match="flat"):
        contour_toolpath(zeros(Grid.centered(0.05, 1e-3)), R, 1e-3)
    with pytest.raises(ValueError):
        Toolpath(np.zeros((1, 3)), np.zeros(1), R)
    with pytest.raises(ValueError):
        HeightMap(Grid.centered(0.05, 1e-3), np.zeros((2, 2)))
