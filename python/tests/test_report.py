"""Reports: figures are written and every length in the text is labelled."""


from precomp.compensation import displacement_adjustment, linear_springback_predictor
from precomp.geometry import Grid
from precomp.metrology import signed_deviation
from precomp.report import markdown_metrics, metrics_table, write_report


def test_report_writes_figures_and_labelled_tables(tmp_path, small_cone):
    target = small_cone.heightmap(Grid.centered(0.08, 1e-3))
    formed = target.with_z(0.9 * target.z)
    dev = signed_deviation(formed, target)
    da = displacement_adjustment(target, linear_springback_predictor(0.9), iterations=3)
    path = write_report(tmp_path / "rep", target, title="Test part", deviation=dev,
                        formed=formed, commanded=da.commanded, tolerance=2e-4,
                        history=da.history, notes=["synthetic"], provenance={"seed": 1})
    text = path.read_text()
    for name in ("deviation_map.png", "deviation_histogram.png", "section_x.png",
                 "history.png"):
        assert (path.parent / name).stat().st_size > 1000
        assert name in text
    assert "RMS [mm]" in text and "Within tolerance [%]" in text and '"seed": 1' in text
    rows = metrics_table(dev, target, 2e-4)
    assert [r["region"] for r in rows] == ["all", "part", "wall", "flange"]
    part = next(r for r in rows if r["region"] == "part")
    assert f"{part['rms'] * 1e3:.4f}" in markdown_metrics(rows)
    # the flange holds the nodes within 1 um of the sheet plane: 0.9 z differs by < 0.1 um
    assert next(r for r in rows if r["region"] == "flange")["max_abs"] < 1e-7
