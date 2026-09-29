"""The sparlab_form deck: structure, supports, trajectory file and content hash."""

import json

import numpy as np
import pandas as pd
import pytest

from precomp import PrecompError
from precomp.fea import (EXECUTABLE_ENV, FormingSetup, build_deck, deck_document, deck_hash,
                         support_321)
from precomp.fea.deck import check_toolpath, make_toolpath
from precomp.geometry import Grid, TruncatedCone
from precomp.materials import get_material


def test_deck_document_matches_the_contract(small_setup):
    doc = deck_document(small_setup)
    m = doc["mesh"]
    assert m["type"] == "structured_hex"
    assert (m["nx"], m["ny"], m["nz"]) == (48, 48, 1)
    assert m["lx"] == pytest.approx(0.12) and m["x0"] == pytest.approx(-0.06)
    assert m["z0"] == -small_setup.thickness and m["lz"] == small_setup.thickness
    assert doc["material"] == small_setup.material.to_sparlab()
    assert doc["nonlinear"]["enabled"] is True
    tool = doc["forming"]["tools"][0]
    assert tool["shape"] == "sphere" and tool["radius"] == small_setup.tool_radius
    assert tool["trajectory"] == {"file": "toolpath.csv"}
    assert tool["surface"]["box"]["zmin"] == 0.0
    steps = doc["forming"]["steps"]
    assert [s["name"] for s in steps] == ["form", "unload", "release"]
    assert [s["type"] for s in steps] == ["form", "release", "release"]
    assert steps[0]["tools"] == ["tool"] and steps[1]["tools"] == []
    assert steps[0]["boundary_conditions"][0]["fix"] == ["x", "y", "z"]
    json.dumps(doc, allow_nan=False)


def test_321_support_removes_exactly_the_six_rigid_modes(small_setup):
    bcs = support_321(small_setup)
    assert sum(len(b["fix"]) for b in bcs) == 6
    p = [np.array(b["region"]["nearest_node"]) for b in bcs]
    assert bcs[0]["fix"] == ["x", "y", "z"] and bcs[1]["fix"] == ["y", "z"] and \
        bcs[2]["fix"] == ["z"]
    # the y-fixed pair lies along x (stops the rotation about z), the three z-fixed
    # points are not collinear (stop the rotations about x and y)
    assert p[1][1] == p[0][1] and p[1][0] > p[0][0]
    assert abs(np.cross(p[1] - p[0], p[2] - p[0])[2]) > 0
    clamped = FormingSetup(small_setup.material, blank_size=0.12, clamp_margin=0.01,
                           element_size=2.5e-3, release="clamped_only")
    assert [s["name"] for s in deck_document(clamped)["forming"]["steps"]] == ["form", "unload"]


def test_build_deck_writes_the_files_and_the_hash_tracks_the_physics(tmp_path, small_setup,
                                                                       small_cone):
    cmd = small_cone.heightmap(Grid.centered(0.12, 1e-3))
    d1 = build_deck(small_setup, cmd, tmp_path / "a")
    for name in ("deck.json", "toolpath.csv", "commanded.npz", "precomp_deck.json"):
        assert (d1 / name).is_file()
    traj = pd.read_csv(d1 / "toolpath.csv")
    assert list(traj.columns) == ["t", "x", "y", "z"]
    assert np.all(np.diff(traj["t"]) > 0) and traj["t"].iloc[-1] == 1.0
    h1 = deck_hash(d1, "v1")
    d2 = build_deck(small_setup.replace(timeout=5.0, threads=4), cmd, tmp_path / "b")
    assert deck_hash(d2, "v1") == h1                   # execution settings do not count
    assert deck_hash(d1, "v2") != h1                   # the solver build does
    d3 = build_deck(small_setup.replace(friction=0.2), cmd, tmp_path / "c")
    assert deck_hash(d3, "v1") != h1                   # the physics does
    d4 = build_deck(small_setup, cmd.with_z(cmd.z * 1.01), tmp_path / "d")
    assert deck_hash(d4, "v1") != h1                   # and so does the commanded shape


def test_a_path_reaching_the_clamp_is_refused(small_setup):
    big = TruncatedCone(0.05, 45.0, 0.01, 0.003, 0.003)
    cmd = big.heightmap(Grid.centered(0.14, 1e-3))
    path = make_toolpath(small_setup, cmd)
    with pytest.raises(PrecompError, match="unclamped"):
        check_toolpath(small_setup, path)


def test_setup_round_trip_validation_and_executable_resolution(tmp_path, monkeypatch):
    s = FormingSetup("DC04", thickness=0.8e-3, contact={"penalty": 1e12})
    back = FormingSetup.from_dict(json.loads(json.dumps(s.to_dict())))
    assert back == s and back.material == get_material("DC04")
    assert "executable" not in s.physics_dict()
    with pytest.raises(ValueError):
        FormingSetup("DC04", element="tet10")
    with pytest.raises(ValueError):
        FormingSetup("DC04", clamp_margin=1e-3, element_size=2.5e-3)
    with pytest.raises(ValueError, match="unknown keys"):
        FormingSetup.from_dict({"material": "DC04", "colour": 1})
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv(EXECUTABLE_ENV, str(tmp_path / "nope"))
    with pytest.raises(PrecompError, match="not found"):
        s.resolved_executable()
    exe = tmp_path / "sf"
    exe.write_text("#!/bin/sh\n")
    exe.chmod(0o755)
    monkeypatch.setenv(EXECUTABLE_ENV, str(exe))
    assert s.resolved_executable() == exe.resolve()
    assert s.replace(executable=str(exe)).resolved_executable() == exe.resolve()
