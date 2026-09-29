"""Running decks through the content-addressed cache, with the test-double solver."""

import numpy as np
import pytest

from conftest import run_count
from precomp import PrecompError
from precomp.fea import (FormingError, cache_entry, run_deck, simulate, simulate_many,
                         sparlab_version)
from precomp.fea.deck import build_deck


def test_simulate_runs_once_and_then_hits_the_cache(tmp_path, small_setup, small_cone,
                                                    node_grid, counter):
    cmd = small_cone.heightmap(node_grid)
    res = simulate(small_setup, cmd, tmp_path / "work")
    assert run_count(counter) == 1
    assert res.provenance["cache_hit"] is False
    key = res.provenance["key"]
    entry = cache_entry(tmp_path / "work", key)
    assert (entry / "COMPLETE").is_file() and (entry / "run.json").is_file()
    assert res.directory == entry / "output"
    formed = res.formed_surface(grid=node_grid)
    assert np.abs(formed.z - 0.9 * cmd.z).max() < 1e-15      # the double's springback
    again = simulate(small_setup, cmd, tmp_path / "work")
    assert run_count(counter) == 1 and again.provenance["cache_hit"] is True
    assert again.provenance["key"] == key
    simulate(small_setup, cmd, tmp_path / "work", cache=False)
    assert run_count(counter) == 2
    assert not list((tmp_path / "work" / "staging").iterdir())


def test_a_failure_is_recorded_and_not_rerun_unless_asked(tmp_path, small_setup, small_cone,
                                                          node_grid, counter, monkeypatch):
    cmd = small_cone.heightmap(node_grid)
    monkeypatch.setenv("FAKE_SPARLAB_FAIL", "1")
    with pytest.raises(FormingError, match="exited with 3") as err:
        simulate(small_setup, cmd, tmp_path / "w")
    assert "failure requested" in err.value.record["log_tail"]
    assert run_count(counter) == 1
    with pytest.raises(FormingError, match="failed before"):
        simulate(small_setup, cmd, tmp_path / "w")
    assert run_count(counter) == 1
    monkeypatch.delenv("FAKE_SPARLAB_FAIL")
    res = simulate(small_setup, cmd, tmp_path / "w", retry_failed=True)
    assert run_count(counter) == 2
    entry = cache_entry(tmp_path / "w", res.provenance["key"])
    assert (entry / "COMPLETE").is_file() and not (entry / "FAILED").exists()


@pytest.mark.parametrize("executor", ["thread", "process"])
def test_simulate_many_keeps_order_records_failures_and_shares_the_cache(
        tmp_path, small_setup, small_cone, node_grid, counter, executor):
    a = small_cone.heightmap(node_grid)
    b = a.with_z(0.8 * a.z)
    bad = small_setup.replace(executable=str(tmp_path / "missing"))
    jobs = [(small_setup, a), (small_setup, b), (bad, a), (small_setup, a)]
    out = simulate_many(jobs, tmp_path / "w", max_workers=2, executor=executor)
    assert [o.index for o in out] == [0, 1, 2, 3]
    assert [o.ok for o in out] == [True, True, False, True]
    assert "not found" in out[2].error
    assert out[0].key == out[3].key != out[1].key
    assert run_count(counter) <= 3                   # the duplicate may race but never corrupts
    formed = out[1].load().formed_surface(grid=node_grid)
    assert np.abs(formed.z - 0.9 * b.z).max() < 1e-15
    with pytest.raises(PrecompError):
        out[2].load()


def test_run_deck_and_version_errors(tmp_path, small_setup, small_cone, node_grid, fake_solver,
                                     counter):
    d = build_deck(small_setup, small_cone.heightmap(node_grid), tmp_path / "deck")
    res = run_deck(d, executable=fake_solver, timeout=60)
    assert res.step_names == ["form", "unload", "release"]
    assert sparlab_version(fake_solver).startswith("sparlab_form")
    broken = tmp_path / "broken"
    broken.write_text("#!/bin/sh\nexit 1\n")
    broken.chmod(0o755)
    with pytest.raises(PrecompError, match="--version exited"):
        sparlab_version(broken)
    with pytest.raises(PrecompError, match="does not exist"):
        sparlab_version(tmp_path / "none")
    slow = tmp_path / "slow"
    slow.write_text("#!/bin/sh\nsleep 5\n")
    slow.chmod(0o755)
    with pytest.raises(FormingError, match="timed out"):
        run_deck(d, executable=slow, timeout=0.2)
