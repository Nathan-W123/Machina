"""The result loader against small synthetic result directories with known answers."""

import json

import numpy as np
import pandas as pd
import pytest

from precomp import PrecompError
from precomp.fea import FORCE_COLUMNS, load_result
from precomp.geometry import Grid


def write_result(d, *, n=8, nz=2, L=0.04, t=1e-3, stretch=0.01, first_k=0,
                 steps=("form", "release"), thinning=0.0, summary_steps=None, forces=True):
    """A structured n x n x nz sheet on [-L/2, L/2]^2 x [-t, 0].

    Displacement: u = (stretch X, stretch Y, w(X, Y)) with w = -0.01 (1 - r^2 / L^2)
    at the top; bottom nodes additionally move up by `thinning` w, so the
    thickness is t (1 + thinning w / t ...) - see `thickness_expected`.
    """
    d.mkdir(parents=True, exist_ok=True)
    xs = np.linspace(-L / 2, L / 2, n + 1)
    zs = np.linspace(-t, 0.0, nz + 1)
    Z, Y, X = np.meshgrid(zs, xs, xs, indexing="ij")
    X, Y, Z = X.ravel(), Y.ravel(), Z.ravel()
    w = -0.01 * (1 - (X ** 2 + Y ** 2) / L ** 2)
    uz = w + np.where(np.isclose(Z, -t), -thinning * w, 0.0)
    for i, name in enumerate(steps):
        factor = 1.0 if i == 0 else 0.9
        frame = pd.DataFrame({"node": np.arange(len(X)), "X": X, "Y": Y, "Z": Z,
                              "ux": stretch * X, "uy": stretch * Y, "uz": factor * uz})
        frame.to_csv(d / f"step_{i + first_k}_{name}_nodes.csv", index=False)
        ne = n * n * nz
        pd.DataFrame({"element": np.arange(ne), "eq_plastic_strain": 0.2,
                      "von_mises": 2e8}).to_csv(d / f"step_{i + first_k}_{name}_elements.csv",
                                                index=False)
    names = steps if summary_steps is None else summary_steps
    (d / "summary.json").write_text(json.dumps({"completed": True,
                                                "steps": [{"name": s} for s in names]}))
    if forces:
        tt = np.linspace(0, 1, 11)
        pd.DataFrame({"step": 0, "increment": np.arange(11), "t": tt, "tool": "tool",
                      "cx": 0.0, "cy": 0.0, "cz": -tt * 0.01, "fx": 0.0, "fy": 10 * tt,
                      "fz": -100 - 400 * tt, "active_nodes": 4})[FORCE_COLUMNS].to_csv(
            d / "tool_forces.csv", index=False)
    return d


@pytest.mark.parametrize("first_k", [0, 1])
def test_formed_surface_is_the_deformed_top_surface(tmp_path, first_k):
    d = write_result(tmp_path / "r", first_k=first_k)
    res = load_result(d)
    assert res.step_names == ["form", "release"]
    assert res.completed
    # grid = the deformed node positions (stretched by 1 %): exact at the nodes
    s = 1.01
    g = Grid(-0.02 * s, -0.02 * s, 9, 9, 0.005 * s)
    f = res.formed_surface("release", grid=g)
    X, Y = g.mesh()
    Xr, Yr = X / s, Y / s
    exact = 0.9 * -0.01 * (1 - (Xr ** 2 + Yr ** 2) / 0.04 ** 2)
    assert f.mask.all()
    assert np.abs(f.z - exact).max() < 1e-15
    assert f.metadata["step"] == "release"
    # a finer grid interpolates linearly inside the deformed triangulation;
    # nodes beyond the deformed plan are invalid
    fine = res.formed_surface(-1, grid=Grid.centered(0.05, 1e-3))
    Xf, Yf = fine.grid.mesh()
    inside = (np.abs(Xf) <= 0.02 * s + 1e-12) & (np.abs(Yf) <= 0.02 * s + 1e-12)
    assert np.array_equal(fine.mask, inside)
    ex = 0.9 * -0.01 * (1 - ((Xf / s) ** 2 + (Yf / s) ** 2) / 0.04 ** 2)
    assert np.abs(fine.z - ex)[inside].max() < 0.9 * 0.01 * (0.005 * s) ** 2 / 0.04 ** 2
    assert res.step(0).name == "form" and res.step("form").index == first_k


def test_thickness_map_measures_node_columns(tmp_path):
    d = write_result(tmp_path / "r", thinning=0.02, stretch=0.0)
    res = load_result(d)
    th = res.thickness_map("form")
    X, Y = th.grid.mesh()
    w = -0.01 * (1 - (X ** 2 + Y ** 2) / 0.04 ** 2)
    assert np.abs(th.z - (1e-3 + 0.02 * w)).max() < 1e-15   # bottom moved by -0.02 w
    assert th.metadata["quantity"] == "thickness"


def test_forces_and_element_tables(tmp_path):
    res = load_result(write_result(tmp_path / "r"))
    f = res.forming_forces("tool")
    assert len(f) == 11 and f["f"].iloc[-1] == pytest.approx(np.hypot(10, 500))
    assert res.plastic_strain(-1)["eq_plastic_strain"].iloc[0] == 0.2
    with pytest.raises(PrecompError, match="no forces"):
        res.forming_forces("other")
    res2 = load_result(write_result(tmp_path / "s", forces=False))
    with pytest.raises(PrecompError, match="tool_forces.csv"):
        res2.forming_forces()


def test_malformed_results_are_refused(tmp_path):
    with pytest.raises(PrecompError, match="not a directory"):
        load_result(write_result(tmp_path / "a") / "missing")
    d = write_result(tmp_path / "b")
    (d / "summary.json").unlink()
    with pytest.raises(PrecompError, match="summary.json"):
        load_result(d)
    d = write_result(tmp_path / "c", summary_steps=["form", "other"])
    with pytest.raises(PrecompError, match="not listed"):
        load_result(d)
    d = write_result(tmp_path / "d")
    frame = pd.read_csv(d / "step_0_form_nodes.csv").drop(columns=["uz"])
    frame.to_csv(d / "step_0_form_nodes.csv", index=False)
    with pytest.raises(PrecompError, match="missing columns"):
        load_result(d)
    d = write_result(tmp_path / "e")
    frame = pd.read_csv(d / "step_1_release_nodes.csv")
    frame.loc[3, "ux"] = np.nan
    frame.to_csv(d / "step_1_release_nodes.csv", index=False)
    with pytest.raises(PrecompError, match="missing values"):
        load_result(d)
    d = write_result(tmp_path / "f")
    frame = pd.read_csv(d / "step_1_release_nodes.csv")
    frame.loc[3, "X"] += 1e-3
    frame.to_csv(d / "step_1_release_nodes.csv", index=False)
    with pytest.raises(PrecompError, match="reference coordinates"):
        load_result(d)
    d = tmp_path / "g"
    d.mkdir()
    (d / "summary.json").write_text("{}")
    with pytest.raises(PrecompError, match="no step_"):
        load_result(d)


def test_a_folded_surface_is_refused(tmp_path):
    d = write_result(tmp_path / "r")
    frame = pd.read_csv(d / "step_1_release_nodes.csv")
    frame["ux"] = -2.0 * frame["X"]                     # mirror x only: every triangle flips
    frame.to_csv(d / "step_1_release_nodes.csv", index=False)
    res = load_result(d)
    res.formed_surface("form")
    with pytest.raises(PrecompError, match="fold"):
        res.formed_surface("release")
