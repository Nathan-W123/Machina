"""The ML command line: dataset generate / info, train, evaluate, active, and
compensate --method surrogate --verify-fea, end to end on proxy data (NOT
physics) with the test-double solver for the verification run."""

import json

import pandas as pd

from conftest import ML_CREATED_AT, run_count
from precomp.cli import main
from precomp.fea import FormingSetup
from precomp.geometry import Grid, TruncatedCone
from precomp.materials import get_material


def test_the_ml_commands_end_to_end(tmp_path, capsys, ml_threads, fake_solver, counter):
    d = tmp_path
    assert main(["dataset", "generate", "--out", str(d / "data"), "--simulator", "proxy",
                 "--families", "truncated_cone", "dome", "pyramid", "--n-per-family", "3",
                 "--materials", "AA5754-O", "--process", "--grid-spacing", "4e-3",
                 "--seed", "2", "--created-at", ML_CREATED_AT]) == 0
    out = capsys.readouterr()
    assert json.loads(out.out)["created"] == 27 and "NOT physics" in out.err
    assert main(["dataset", "info", "--data", str(d / "data")]) == 0
    info = json.loads(capsys.readouterr().out)
    assert info["samples"] == 27 and info["parts"] == 9
    assert info["data_source"] == "proxy - not physics"
    assert main(["train", "--data", str(d / "data"), "--model", "gbm", "--out",
                 str(d / "model"), "--param", "n_members=2", "max_iter=40",
                 "--time-source", "depth", "--points-per-sample", "200",
                 "--test-fraction", "0.3", "--calibration-fraction", "0.4",
                 "--created-at", ML_CREATED_AT]) == 0
    tr = json.loads(capsys.readouterr().out)
    assert tr["data_source"] == "proxy - not physics" and tr["n_test"] == 9
    assert tr["metrics"]["held_out"]["data_source"] == "proxy - not physics"
    man = json.loads((d / "model" / "manifest.json").read_text())
    assert man["created_at"] == ML_CREATED_AT and man["model_class"] == "GBMEnsemble"
    split = json.loads((d / "model" / "split.json").read_text())
    assert set(split["train"]).isdisjoint(split["test"])
    assert main(["evaluate", "--model", str(d / "model"), "--data", str(d / "data"),
                 "--out", str(d / "eval")]) == 0
    capsys.readouterr()
    pp = pd.read_csv(d / "eval" / "per_part.csv")
    assert len(pp) == 9 and set(pp["data_source"]) == {"proxy - not physics"}
    assert set(pp["sample_id"]) == set(split["test"])
    assert (d / "eval" / "calibration.csv").is_file()
    assert main(["active", "--model", str(d / "model"), "--n-candidates", "6", "--select", "3",
                 "--families", "truncated_cone", "dome", "--materials", "AA5754-O",
                 "--grid-spacing", "4e-3", "--out", str(d / "rank.csv")]) == 0
    capsys.readouterr()
    rank = pd.read_csv(d / "rank.csv")
    assert len(rank) == 6 and sorted(rank["rank"].dropna()) == [1, 2, 3]
    assert set(rank["model_data_source"]) == {"proxy - not physics"}
    # compensation on the surrogate, verified by the (test-double) solver
    setup = FormingSetup(get_material("AA5754-O"), executable=str(fake_solver), layers=1,
                         element_size=4e-3, step_down=1.5e-3, toolpath_spacing=3e-3)
    (d / "setup.json").write_text(json.dumps(setup.to_dict()))
    TruncatedCone(0.045, 45.0, 0.015, 0.005, 0.005).heightmap(
        Grid.centered(setup.meshed_blank_size, 4e-3)).save(d / "t.npz")
    base = ["compensate", "--target", str(d / "t.npz"), "--setup", str(d / "setup.json"),
            "--method", "surrogate", "--model", str(d / "model"), "--iterations", "2",
            "--out", str(d / "comp")]
    assert main(base + ["--verify-fea", "--no-verify"]) == 2
    assert "contradict" in capsys.readouterr().err
    assert main(base + ["--verify-fea"]) == 2                  # no --work-dir
    capsys.readouterr()
    assert main(base + ["--verify-fea", "--work-dir", str(d / "runs")]) == 0
    capsys.readouterr()
    assert run_count(counter) == 1
    doc = json.loads((d / "comp" / "compensation.json").read_text())
    assert doc["method"] == "surrogate" and doc["verification"]["source"] == "fea"
    assert doc["ood"]["model_data_source"] == "proxy - not physics"
    assert (d / "comp" / "predicted_lower.npz").is_file()
    report = (d / "comp" / "report" / "report.md").read_text()
    assert "proxy - not physics" in report


def test_ml_usage_errors_exit_with_status_2(tmp_path, capsys):
    assert main(["train", "--data", str(tmp_path / "missing"), "--out", str(tmp_path / "m")]) == 2
    assert "not a precomp.ml data set" in capsys.readouterr().err
    assert main(["active"]) == 2
    assert "precomp.ml" in capsys.readouterr().err
    assert main(["dataset", "generate", "--out", str(tmp_path / "d"), "--simulator", "sparlab",
                 "--n-per-family", "1", "--families", "dome"]) == 2
    assert "--work-dir" in capsys.readouterr().err
