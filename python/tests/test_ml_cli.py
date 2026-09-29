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
    (d / "base.json").write_text(json.dumps({"toolpath_style": "contour"}))
    assert main(["dataset", "generate", "--out", str(d / "data"), "--simulator", "proxy",
                 "--families", "truncated_cone", "dome", "pyramid", "--n-per-family", "3",
                 "--materials", "AA5754-O", "--process", "--grid-spacing", "4e-3",
                 "--setup", str(d / "base.json"), "--seed", "2",
                 "--created-at", ML_CREATED_AT]) == 0
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
    # no calibration: coverage cannot be computed, and says so (null, not NaN)
    assert main(["train", "--data", str(d / "data"), "--model", "gbm", "--out",
                 str(d / "nocal"), "--param", "n_members=2", "max_iter=20",
                 "--time-source", "depth", "--points-per-sample", "100",
                 "--test-fraction", "0.3", "--calibration-fraction", "0",
                 "--created-at", ML_CREATED_AT]) == 0
    assert json.loads(capsys.readouterr().out)["metrics"]["held_out"]["coverage_mean"] is None
    assert json.loads((d / "nocal" / "manifest.json").read_text())["conformal"] is None
    assert main(["evaluate", "--model", str(d / "nocal"), "--data", str(d / "data"),
                 "--out", str(d / "eval_nocal")]) == 0
    capsys.readouterr()
    assert json.loads((d / "eval_nocal" / "summary.json").read_text())["held_out"][
        "coverage_mean"] is None
    assert main(["evaluate", "--model", str(d / "model"), "--data", str(d / "data"),
                 "--out", str(d / "eval")]) == 0
    capsys.readouterr()
    pp = pd.read_csv(d / "eval" / "per_part.csv")
    assert len(pp) == 9 and set(pp["data_source"]) == {"proxy - not physics"}
    assert set(pp["sample_id"]) == set(split["test"])
    assert (d / "eval" / "calibration.csv").is_file()
    # candidates in the training setup by default: contour paths, AA5754-O,
    # the fixed process, the 4 mm grid - they differ in geometry only
    assert main(["active", "--model", str(d / "model"), "--n-candidates", "6", "--select", "3",
                 "--families", "truncated_cone", "dome", "--out", str(d / "rank.csv")]) == 0
    capsys.readouterr()
    rank = pd.read_csv(d / "rank.csv")
    assert len(rank) == 6 and sorted(rank["rank"].dropna()) == [1, 2, 3]
    assert set(rank["model_data_source"]) == {"proxy - not physics"}
    assert set(rank["material"]) == {"AA5754-O"}
    reasons = ";".join(rank["ood_reasons"].fillna(""))
    for key in ("setup.", "step_down", "tool_radius", "thickness", "friction"):
        assert key not in reasons, reasons
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
    # the double's setup is not the training setup: refused without the override
    assert main(base + ["--verify-fea", "--work-dir", str(d / "runs")]) == 2
    assert "outside the model's training envelope" in capsys.readouterr().err
    assert main(base + ["--verify-fea", "--work-dir", str(d / "runs"),
                        "--allow-out-of-envelope"]) == 0
    err = capsys.readouterr().err
    assert "training envelope: target and compensated shape OUTSIDE" in err
    # two calibration parts cannot support a 90 % interval: none, and why
    assert "no prediction interval" in err and "needs at least 9 held-out parts" in err
    assert run_count(counter) == 1
    doc = json.loads((d / "comp" / "compensation.json").read_text())
    assert doc["method"] == "surrogate" and doc["verification"]["source"] == "fea"
    assert doc["ood"]["model_data_source"] == "proxy - not physics"
    assert doc["model_data_source"] == "proxy - not physics" and doc["in_envelope"] is False
    assert "not simulations" in doc["history_quantity"]
    assert not (d / "comp" / "predicted_lower.npz").exists()
    assert "needs at least 9" in doc["interval_note"]
    assert (d / "comp" / "verified.npz").is_file()
    report = (d / "comp" / "report" / "report.md").read_text()
    assert "proxy - not physics" in report and "OUTSIDE" in report
    assert "sparlab_form simulation of the compensated shape" in report


def test_ml_usage_errors_exit_with_status_2(tmp_path, capsys):
    assert main(["train", "--data", str(tmp_path / "missing"), "--out", str(tmp_path / "m")]) == 2
    assert "not a precomp.ml data set" in capsys.readouterr().err
    assert main(["active"]) == 2
    assert "precomp.ml" in capsys.readouterr().err
    assert main(["dataset", "generate", "--out", str(tmp_path / "d"), "--simulator", "sparlab",
                 "--n-per-family", "1", "--families", "dome"]) == 2
    assert "--work-dir" in capsys.readouterr().err
