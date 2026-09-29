"""The high-level API and the `precomp` command, end to end with the test-double solver."""

import json

import numpy as np
import pytest

from conftest import run_count
from precomp.api import compare_scan, compensate, predict
from precomp.cli import main
from precomp.geometry import EllipticCone, Grid, HeightMap
from precomp.metrology import part_mask, rotation_matrix


class _Model:
    """A stand-in surrogate: dz = -0.1 z, with an interval and an OOD report."""

    def predict_deviation(self, commanded, setup):
        return -0.1 * commanded.z, np.full(commanded.grid.shape, 2e-5)

    def predict_interval(self, commanded, setup, level):
        m = -0.1 * commanded.z
        return m - 4e-5, m + 4e-5

    def assess(self, commanded, setup):
        return {"in_envelope": True, "reasons": []}


def test_predict_fea_and_surrogate(tmp_path, small_setup, small_cone, node_grid, counter):
    cmd = small_cone.heightmap(node_grid)
    p = predict(cmd, small_setup, method="fea", work_dir=tmp_path / "w")
    assert np.abs(p.formed.z - 0.9 * cmd.z).max() < 1e-15
    assert p.details["fea"]["cache_hit"] is False and p.std is None
    s = predict(cmd, small_setup, _Model(), method="surrogate")
    assert np.allclose(s.formed.z, 0.9 * cmd.z) and np.all(s.std == 2e-5)
    with pytest.raises(ValueError, match="work_dir"):
        predict(cmd, small_setup, method="fea")
    with pytest.raises(ValueError, match="needs a model"):
        predict(cmd, small_setup, method="surrogate")


def test_compensate_with_fea_converges_to_the_springback_inverse(tmp_path, small_setup,
                                                                 small_cone, node_grid,
                                                                 counter):
    """The double's formed shape is 0.9 x commanded at the nodes, so DA with
    alpha = 1 cuts the error tenfold per iteration and the compensated shape
    tends to target / 0.9."""
    target = small_cone.heightmap(node_grid)
    res = compensate(target, small_setup, method="fea", iterations=4, work_dir=tmp_path / "w")
    rms = [h["error"]["rms"] for h in res.history]
    assert all(r1 / r0 == pytest.approx(0.1, rel=1e-6) for r0, r1 in zip(rms, rms[1:]))
    part = part_mask(target)
    assert np.abs(res.compensated.z - target.z / 0.9)[part].max() < 1e-3 * target.depth
    assert run_count(counter) == 4
    v = res.verification
    assert v["source"].startswith("fea") and v["vertical_deviation"]["part"]["rms"] == rms[-1]
    out = res.save(tmp_path / "out")
    assert (out / "toolpath.csv").is_file() and (out / "compensated.npz").is_file()
    doc = json.loads((out / "compensation.json").read_text())
    assert doc["method"] == "fea" and len(doc["history"]) == 4


def test_compensate_with_a_surrogate_reports_interval_and_ood(tmp_path, small_setup, small_cone,
                                                               node_grid, counter):
    target = small_cone.heightmap(node_grid)
    res = compensate(target, small_setup, _Model(), method="surrogate", iterations=3,
                     verify=True, work_dir=tmp_path / "w")
    lo, hi = res.interval
    assert np.allclose(hi.z - lo.z, 8e-5)
    assert np.allclose(0.5 * (lo.z + hi.z), res.predicted.z)
    assert res.ood == {"in_envelope": True, "reasons": []}
    assert run_count(counter) == 1                     # one verification run
    assert res.verification["vertical_deviation"]["part"]["rms"] < 1e-5


def test_compare_scan_measures_a_misplaced_scan(tmp_path, rng):
    part = EllipticCone(0.05, 0.04, 45.0, 0.012, 0.004, 0.004)
    target = part.heightmap(Grid.centered(0.12, 5e-4))
    formed = target.with_z(0.95 * target.z)               # 5 % shallower everywhere
    X, Y = Grid.centered(0.11, 5e-4).mesh()
    pts = np.column_stack([X.ravel(), Y.ravel(), formed.interpolate(X.ravel(), Y.ravel())])
    # held by the fixture: tilted and lifted, not slid (a flat flange could
    # not see an in-plane shift - rank 3)
    R = rotation_matrix(np.radians([0.3, -0.2, 0.0]))
    scan_mm = (pts @ R.T + np.array([0.0, 0.0, 2e-3])) * 1e3
    np.savetxt(tmp_path / "scan.xyz", scan_mm, fmt="%.9f")
    cmp = compare_scan(tmp_path / "scan.xyz", target, scale=1e-3, fixture="flange",
                       tolerance=1e-4)
    assert cmp.alignment.rank == 3
    assert cmp.metrics["flange"]["max_abs"] < 1e-6
    vert = formed.z - target.z                               # expected: ~ +0.05 depth
    sel = part_mask(target) & cmp.deviation.mask
    n = target.normals()[..., 2]
    assert np.abs(cmp.deviation.z - vert * n)[sel].max() < 0.1 * np.abs(vert).max()
    assert cmp.metrics["part"]["bias"] > 0                   # shallower: towards the tool


def test_cli_end_to_end(tmp_path, small_setup, counter, capsys):
    d = tmp_path
    assert main(["part", "--family", "truncated_cone", "--param", "top_radius=0.03",
                 "wall_angle_deg=45", "depth=0.012", "top_fillet=0.004", "bottom_fillet=0.004",
                 "--size", "0.12", "--spacing", "2.5e-3", "--out", str(d / "t.npz"),
                 "--stl", str(d / "t.stl"), "--json-out", str(d / "t.json")]) == 0
    assert HeightMap.load(d / "t.npz").depth == pytest.approx(0.012)
    (d / "setup.json").write_text(json.dumps(small_setup.to_dict()))
    assert main(["toolpath", "--surface", str(d / "t.npz"), "--setup", str(d / "setup.json"),
                 "--out", str(d / "path.csv"), "--robot", str(d / "robot.csv"),
                 "--feed", "0.02"]) == 0
    assert (d / "path.csv").read_text().startswith("t,x,y,z\n")
    assert main(["simulate", "--setup", str(d / "setup.json"), "--commanded", str(d / "t.npz"),
                 "--work-dir", str(d / "w"), "--out", str(d / "formed.npz")]) == 0
    assert main(["compensate", "--target", str(d / "t.npz"), "--setup", str(d / "setup.json"),
                 "--iterations", "2", "--work-dir", str(d / "w"), "--out",
                 str(d / "comp")]) == 0
    assert (d / "comp" / "report" / "report.md").is_file()
    assert run_count(counter) == 2              # iteration 0 is the cached simulate run
    target = HeightMap.load(d / "t.npz")
    np.save(d / "scan.npy", HeightMap.load(d / "formed.npz").to_points())
    assert main(["scan", "compare", "--scan", str(d / "scan.npy"), "--target", str(d / "t.npz"),
                 "--align", "none", "--tolerance", "1e-4", "--out", str(d / "cmp")]) == 0
    assert (d / "cmp" / "report.md").is_file() and (d / "cmp" / "deviation.npz").is_file()
    assert main(["scan", "update", "--scan", str(d / "scan.npy"), "--commanded", str(d / "t.npz"),
                 "--target", str(d / "t.npz"), "--align", "none", "--out",
                 str(d / "next.npz")]) == 0
    nxt = HeightMap.load(d / "next.npz")
    assert nxt.depth > target.depth
    assert main(["report", "--target", str(d / "t.npz"), "--formed", str(d / "formed.npz"),
                 "--out", str(d / "rep")]) == 0
    assert main(["material", "DC04", "--sparlab"]) == 0
    assert '"youngs_modulus"' in capsys.readouterr().out
    assert main(["setup", "--material", "DC04", "--set", "thickness=8e-4",
                 "--out", str(d / "s2.json")]) == 0
    assert json.loads((d / "s2.json").read_text())["thickness"] == 8e-4


def test_cli_errors_are_reported_not_raised(tmp_path, capsys):
    assert main(["part", "--family", "nonsense"]) == 2
    assert "unknown family" in capsys.readouterr().err
    assert main(["train", "--anything"]) == 2
    assert "precomp.ml" in capsys.readouterr().err
    assert main(["simulate", "--setup", str(tmp_path / "none.json"), "--commanded", "x.npz",
                 "--work-dir", str(tmp_path)]) == 2
