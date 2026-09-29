"""The result loader against small synthetic result directories with known answers."""

import json

import numpy as np
import pandas as pd
import pytest

from precomp import PrecompError
from precomp.fea import FORCE_COLUMNS, load_result
from precomp.fea.results import files_stem
from precomp.geometry import Grid


def write_result(d, *, n=8, nz=2, L=0.04, t=1e-3, stretch=0.01,
                 steps=("form", "release"), thinning=0.0, completed=None, forces=True):
    """A structured n x n x nz sheet on [-L/2, L/2]^2 x [-t, 0], written as
    sparlab_form writes a result (docs/forming.md, section 3).

    Displacement: u = (stretch X, stretch Y, w(X, Y)) with w = -0.01 (1 - r^2 / L^2)
    at the top, times 0.9 after the first step; bottom nodes additionally move
    up by `thinning` w. `completed` (default all True) marks the steps that
    completed; only those get files.
    """
    d.mkdir(parents=True, exist_ok=True)
    xs = np.linspace(-L / 2, L / 2, n + 1)
    zs = np.linspace(-t, 0.0, nz + 1)
    Z, Y, X = np.meshgrid(zs, xs, xs, indexing="ij")
    X, Y, Z = X.ravel(), Y.ravel(), Z.ravel()
    w = -0.01 * (1 - (X ** 2 + Y ** 2) / L ** 2)
    uz = w + np.where(np.isclose(Z, -t), -thinning * w, 0.0)
    completed = [True] * len(steps) if completed is None else list(completed)
    entries = []
    for i, name in enumerate(steps):
        stem = f"step_{i + 1}_{name}"
        entries.append({"name": name, "type": "form" if i == 0 else "release",
                        "completed": completed[i], "files_stem": stem})
        if not completed[i]:
            continue
        factor = 1.0 if i == 0 else 0.9
        frame = pd.DataFrame({"node": np.arange(len(X)), "X": X, "Y": Y, "Z": Z,
                              "ux": stretch * X, "uy": stretch * Y, "uz": factor * uz})
        frame.to_csv(d / f"{stem}_nodes.csv", index=False)
        ne = n * n * nz
        pd.DataFrame({"element": np.arange(ne), "eq_plastic_strain": 0.2,
                      "von_mises_Pa": 2e8}).to_csv(d / f"{stem}_elements.csv", index=False)
    (d / "summary.json").write_text(json.dumps({
        "completed": all(completed), "steps": entries,
        "termination": "completed every step" if all(completed) else "stopped"}))
    if forces:
        tt = np.linspace(0, 1, 11)
        pd.DataFrame({"step": 1, "increment": np.arange(1, 12), "t": tt, "tool": "tool",
                      "cx": 0.0, "cy": 0.0, "cz": -tt * 0.01, "fx": 0.0, "fy": -10 * tt,
                      "fz": 100 + 400 * tt, "active_nodes": 4,
                      "max_penetration_m": 1e-6})[FORCE_COLUMNS].to_csv(
            d / "tool_forces.csv", index=False)
    return d


def test_formed_surface_is_the_deformed_top_surface(tmp_path):
    d = write_result(tmp_path / "r")
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
    assert res.step(0).name == "form" and res.step("form").index == 1
    assert res.step("release").info["files_stem"] == "step_2_release"


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
    assert (f["fz"] > 0).all()                      # the force on the tool
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
    d = write_result(tmp_path / "c")
    (d / "step_2_release_nodes.csv").unlink()
    with pytest.raises(PrecompError, match="node table of completed step 'release'"):
        load_result(d)
    d = write_result(tmp_path / "d")
    frame = pd.read_csv(d / "step_1_form_nodes.csv").drop(columns=["uz"])
    frame.to_csv(d / "step_1_form_nodes.csv", index=False)
    with pytest.raises(PrecompError, match="missing columns"):
        load_result(d)
    d = write_result(tmp_path / "e")
    frame = pd.read_csv(d / "step_2_release_nodes.csv")
    frame.loc[3, "ux"] = np.nan
    frame.to_csv(d / "step_2_release_nodes.csv", index=False)
    with pytest.raises(PrecompError, match="missing values"):
        load_result(d)
    d = write_result(tmp_path / "f")
    frame = pd.read_csv(d / "step_2_release_nodes.csv")
    frame.loc[3, "X"] += 1e-3
    frame.to_csv(d / "step_2_release_nodes.csv", index=False)
    with pytest.raises(PrecompError, match="reference coordinates"):
        load_result(d)
    d = write_result(tmp_path / "g")
    frame = pd.read_csv(d / "step_1_form_elements.csv").rename(
        columns={"von_mises_Pa": "von_mises"})
    frame.to_csv(d / "step_1_form_elements.csv", index=False)
    with pytest.raises(PrecompError, match="von_mises_Pa"):
        load_result(d)
    d = write_result(tmp_path / "h")
    doc = json.loads((d / "summary.json").read_text())
    doc["steps"][1]["files_stem"] = "step_1_release"          # not step_<k>_<name>
    (d / "summary.json").write_text(json.dumps(doc))
    with pytest.raises(PrecompError, match="files_stem"):
        load_result(d)
    d = tmp_path / "i"
    d.mkdir()
    (d / "summary.json").write_text("{}")
    with pytest.raises(PrecompError, match="no 'steps'"):
        load_result(d)


def test_a_stopped_run_loads_its_completed_steps_only(tmp_path):
    d = write_result(tmp_path / "r", steps=("form", "unload", "release"),
                     completed=(True, True, False))
    res = load_result(d)
    assert res.step_names == ["form", "unload"] and not res.completed
    assert res.termination == "stopped"
    d = write_result(tmp_path / "s", completed=(False, False))
    with pytest.raises(PrecompError, match="no step completed"):
        load_result(d)


def test_step_file_stems_follow_sparlab_form():
    # forming_step_stem (src/io/FormingWriter.cpp): bytes of the name other than
    # ASCII letters, digits, '-' and '_' become '_'; an empty name is "step"
    assert files_stem(1, "form") == "step_1_form"
    assert files_stem(3, "re lease/2") == "step_3_re_lease_2"
    assert files_stem(2, "\u00e9") == "step_2___"                  # two UTF-8 bytes
    assert files_stem(4, "") == "step_4_step"


def test_snapshot_files_are_not_steps(tmp_path):
    d = write_result(tmp_path / "r")
    frame = pd.read_csv(d / "step_1_form_nodes.csv")
    frame.to_csv(d / "step_1_form_inc_5_nodes.csv", index=False)
    assert load_result(d).step_names == ["form", "release"]


def test_a_folded_surface_is_refused(tmp_path):
    d = write_result(tmp_path / "r")
    frame = pd.read_csv(d / "step_2_release_nodes.csv")
    frame["ux"] = -2.0 * frame["X"]                     # mirror x only: every triangle flips
    frame.to_csv(d / "step_2_release_nodes.csv", index=False)
    res = load_result(d)
    res.formed_surface("form")
    with pytest.raises(PrecompError, match="fold"):
        res.formed_surface("release")
