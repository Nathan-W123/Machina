"""python/scripts/compare_forming_runs.py: the comparison of an explicit
forming run with its implicit reference (docs/forming.md, "Explicit forming"),
on two synthetic result directories whose differences are known exactly."""

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "compare_forming_runs.py"


def load_script():
    spec = importlib.util.spec_from_file_location("compare_forming_runs", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_run(d: Path, depth: float, springback: float, force: float, ripple: float) -> None:
    """A 3 x 3 grid of top nodes (z = 0) and bottom nodes (z = -1 mm): the
    formed shape sinks the centre node by `depth`, the release lifts it by
    `springback`; the forming force is `force` plus a ripple of `ripple`."""
    d.mkdir(parents=True)
    xs = np.array([-1.0, 0.0, 1.0]) * 1e-3
    rows = []
    node = 0
    for z in (0.0, -1e-3):
        for y in xs:
            for x in xs:
                rows.append((node, x, y, z))
                node += 1
    centre = [r[0] for r in rows if r[1] == 0.0 and r[2] == 0.0]

    def write_nodes(name: str, uz_centre: float) -> None:
        with open(d / name, "w") as f:
            f.write("node,X,Y,Z,ux,uy,uz\n")
            for n, x, y, z in rows:
                uz = uz_centre if n in centre else 0.0
                f.write(f"{n},{x},{y},{z},0,0,{uz}\n")

    write_nodes("step_1_form_nodes.csv", -depth)
    write_nodes("step_2_retract_nodes.csv", -depth)
    write_nodes("step_3_unclamp_nodes.csv", -depth + springback)
    steps = [{"name": "form"}, {"name": "retract"}, {"name": "unclamp"}]
    (d / "summary.json").write_text(json.dumps({"steps": steps}))
    with open(d / "tool_forces.csv", "w") as f:
        f.write("step,increment,t,tool,cx,cy,cz,fx,fy,fz,active_nodes,max_penetration_m\n")
        for i, t in enumerate(np.linspace(0.0, 1.0, 401)):
            fz = force + ripple * (1.0 if i % 2 else -1.0)
            f.write(f"1,{i},{t},tool,0,0,0,0,0,{fz},3,0\n")


def test_the_comparison_reports_the_known_shape_springback_and_force_differences(tmp_path):
    compare = load_script().compare
    write_run(tmp_path / "ref", depth=2.0e-3, springback=0.2e-3, force=400.0, ripple=0.0)
    write_run(tmp_path / "exp", depth=2.1e-3, springback=0.15e-3, force=420.0, ripple=50.0)
    out = compare(tmp_path / "ref", tmp_path / "exp", "form", "unclamp", 0.1)
    assert out["top_nodes"] == 9
    assert out["reference_springback"]["max_mm"] == pytest.approx(0.2, rel=1e-12)
    assert out["explicit_springback"]["max_mm"] == pytest.approx(0.15, rel=1e-12)
    # Only the centre node differs: max = its difference, rms = that / 3.
    assert out["formed_shape_difference"]["max_mm"] == pytest.approx(0.1, rel=1e-9)
    assert out["formed_shape_difference"]["rms_mm"] == pytest.approx(0.1 / 3.0, rel=1e-9)
    assert out["springback_difference"]["max_mm"] == pytest.approx(0.05, rel=1e-9)
    assert out["final_shape_difference"]["max_mm"] == pytest.approx(0.15, rel=1e-9)
    assert out["reference_depth_mm"] == pytest.approx(2.0, rel=1e-12)
    # The ripple averages out of the moving average; the mean is 5 % high.
    tf = out["tool_force"]
    assert tf["mean_difference"] == pytest.approx(0.05, abs=1e-3)
    assert tf["moving_average_max_difference"] == pytest.approx(0.05, abs=2e-3)
    assert tf["explicit_raw_rms_scatter"] > 0.05


def test_runs_on_different_meshes_are_refused(tmp_path):
    compare = load_script().compare
    write_run(tmp_path / "ref", depth=2.0e-3, springback=0.2e-3, force=400.0, ripple=0.0)
    write_run(tmp_path / "exp", depth=2.0e-3, springback=0.2e-3, force=400.0, ripple=0.0)
    nodes = tmp_path / "exp" / "step_3_unclamp_nodes.csv"
    nodes.write_text(nodes.read_text().replace("\n0,-0.001", "\n0,-0.002", 1))
    with pytest.raises(SystemExit, match="do not share the mesh"):
        compare(tmp_path / "ref", tmp_path / "exp", "form", "unclamp", 0.1)
