"""Displacement adjustment on synthetic linear springback operators with known inverses."""

import numpy as np
import pytest

from precomp.compensation import (CompositePredictor, SurrogatePredictor,
                                  displacement_adjustment, limit_wall_angle,
                                  linear_springback_predictor, smooth_update, update_from_scan)
from precomp.geometry import Grid, Pyramid, TruncatedCone
from precomp.metrology import flange_mask, part_mask, rotation_matrix


@pytest.fixture
def target():
    return TruncatedCone(0.035, 45.0, 0.015, 0.005, 0.005).heightmap(Grid.centered(0.1, 1e-3))


def test_da_converges_geometrically_at_the_spectral_radius(target):
    """formed = c + G c with G = -0.2 S (S a Gaussian smoothing, ||S|| <= 1):
    on the adjusted region e_(k+1) = -G_PP e_k, so the error falls by at
    most rho(G) = 0.2 per iteration with alpha = 1."""
    pred = linear_springback_predictor(0.8, smoothing=0.004)
    res = displacement_adjustment(target, pred, iterations=7, alpha=1.0,
                                  max_wall_angle_deg=None)
    rms = np.array(res.rms_history)
    ratios = rms[1:] / rms[:-1]
    assert np.all(ratios <= 0.2 + 1e-9)
    assert np.all(ratios > 0.15)                       # and not faster: G is that strong
    assert rms[-1] < 1e-4 * rms[0]
    assert res.best_iteration == len(rms) - 1
    assert res.commanded.depth > target.depth          # over-bent to beat springback
    # a half step converges more slowly: e_(k+1) = (I - alpha (I + G)) e_k
    half = displacement_adjustment(target, pred, iterations=4, alpha=0.5,
                                   max_wall_angle_deg=None)
    r = np.array(half.rms_history)
    assert np.all(r[1:] / r[:-1] < 0.61) and np.all(r[1:] / r[:-1] > 0.5)


def test_da_holds_the_flange_and_stays_below_the_sheet(target):
    seen = []

    def spy(k, c, f, err):
        seen.append(c.z.copy())

    pred = linear_springback_predictor(0.7, smoothing=0.01)   # spills onto the flange
    res = displacement_adjustment(target, pred, iterations=4, smoothing=0.002, callback=spy)
    fl = flange_mask(target)
    for z in seen + [res.proposed.z]:
        assert np.array_equal(z[fl], target.z[fl])
        assert z.max() <= 0.0
    assert res.rms_history[-1] < 0.2 * res.rms_history[0]


def test_formability_projection_only_raises_and_bounds_the_wall_angle():
    part = Pyramid(0.045, 0.04, 60.0, 0.02, 0.04, 0.003, 0.003)
    hm = part.heightmap(Grid.centered(0.12, 1e-3))
    steep = hm.with_z(1.8 * hm.z)
    assert np.degrees(steep.wall_angle().max()) > 70
    out = limit_wall_angle(steep, 65.0)
    assert np.all(out.z >= steep.z)
    assert np.degrees(out.wall_angle().max()) <= 65.0 + 1e-9
    assert out.metadata["wall_angle_limit"]["nodes_raised"] > 0
    same = limit_wall_angle(hm, 65.0)
    assert np.array_equal(same.z, hm.z)                 # formable input is untouched
    res = displacement_adjustment(hm, linear_springback_predictor(0.6), iterations=3)
    assert np.degrees(res.proposed.wall_angle().max()) <= 65.0 + 1e-9


def test_smoothing_does_not_leak_across_the_region_boundary():
    g = Grid.centered(0.02, 1e-3)
    region = np.zeros(g.shape, bool)
    region[:, :10] = True
    u = np.where(region, 1.0, 50.0)
    s = smooth_update(u, region, 3e-3, g.h)
    assert np.allclose(s[region], 1.0) and np.all(s[~region] == 0.0)


def test_normal_direction_da_converges_too(target):
    pred = linear_springback_predictor(0.85)
    res = displacement_adjustment(target, pred, iterations=4, direction="normal",
                                  max_wall_angle_deg=None)
    r = res.rms_history
    assert r[-1] < 0.01 * r[0]


def test_update_from_scan_is_one_da_step(target):
    """A scan of the formed part taken at the grid nodes (where linear
    gridding is exact whatever the triangulation), misplaced in the fixture
    by a small tilt and lift: after alignment on the flange the update is
    exactly c - (formed - target) on the part."""
    commanded = target.copy()
    formed = linear_springback_predictor(0.8)(commanded)
    pts = formed.to_points()
    R = rotation_matrix(np.radians([0.2, -0.1, 0.0]))
    scan = pts @ R.T + np.array([0.0, 0.0, 4e-4])
    new, report = update_from_scan(commanded, scan, target, alpha=1.0, fixture="flange",
                                   max_wall_angle_deg=None)
    assert report["alignment"]["rank"] == 3
    assert report["alignment"]["rms_m"] < 1e-10               # the misplacement is undone
    assert report["coverage"] == 1.0
    part = part_mask(target) & ~flange_mask(target)
    exact = commanded.z - (formed.z - target.z)
    # A flat flange cannot fix the in-plane position (rank 3): it stays where
    # the pivot of the tilt puts it, ~1e-7 m off here, which on a 45 deg wall
    # reads as the same in z. Everything else is exact.
    assert np.abs(new.z - exact)[part].max() < 3e-7
    assert np.abs(new.z - exact)[part].mean() < 3e-8
    half, _ = update_from_scan(commanded, scan, target, alpha=0.5, fixture="flange",
                               max_wall_angle_deg=None)
    assert np.abs(half.z - (commanded.z - 0.5 * (formed.z - target.z)))[part].max() < 3e-7


class _ConstantModel:
    """A stand-in learned model: dz = -0.1 z_commanded, std 1e-5 m."""

    def predict_deviation(self, commanded, setup):
        return -0.1 * commanded.z, np.full(commanded.grid.shape, 1e-5)


def test_surrogate_and_composite_predictors(target):
    sp = SurrogatePredictor(_ConstantModel(), setup=None)
    formed, std = sp.predict_with_uncertainty(target)
    assert np.allclose(formed.z, 0.9 * target.z) and np.all(std == 1e-5)
    res = displacement_adjustment(target, sp, iterations=4, max_wall_angle_deg=None)
    assert res.rms_history[-1] < 1e-3 * res.rms_history[0]
    comp = CompositePredictor(linear_springback_predictor(1.0), _ConstantModel(), None)
    assert np.allclose(comp(target).z, 0.9 * target.z)
    with pytest.raises(TypeError):
        SurrogatePredictor(object(), None)
