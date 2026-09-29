"""Surrogate compensation, verified by the ProxySimulator (NOT physics) and,
for the plumbing of the SparLab path, by the test-double solver."""

import json

import numpy as np

from conftest import ML_CREATED_AT, run_count
from precomp import api
from precomp.fea import FormingSetup
from precomp.geometry import Grid, TruncatedCone
from precomp.materials import get_material
from precomp.ml import (ProxySimulator, save_model, surrogate_compensate, verify_with_fea,
                        verify_with_simulator)


def test_surrogate_compensation_cuts_the_verified_deviation_by_a_large_factor(
        gbm_surrogate, proxy_split):
    """On held-out parts: DA on the surrogate, then the compensated shape
    formed by the proxy. The deviation from the target falls several-fold
    against forming the target itself (both measured by the proxy)."""
    parts = [s for s in proxy_split["test"] if s.kind == "uncompensated"][:3]
    factors = []
    for s in parts:
        res = surrogate_compensate(s.target, s.setup, gbm_surrogate, iterations=6)
        assert res.stopped in ("stagnation", "iterations", "tolerance")
        assert res.data_source == "proxy - not physics"
        assert res.compensated.depth > s.target.depth         # over-bent against springback
        lo, hi = res.residual_interval
        assert np.all(lo.z <= res.residual.z + 1e-15) and np.all(res.residual.z <= hi.z + 1e-15)
        assert res.ood is not None and "in_envelope" in res.ood
        v = verify_with_simulator(res, s.setup, ProxySimulator())
        assert v["data_source"] == "proxy - not physics" and res.verification is v
        # the prediction is optimistic, the achieved deviation is what counts
        assert v["achieved_part_rms_m"] >= 0.5 * v["predicted_part_rms_m"]
        factors.append(v["improvement_factor"])
    assert min(factors) > 2.5 and np.median(factors) > 4.0, factors


def test_stopping_rules(tmp_path, gbm_surrogate, proxy_split):
    s = proxy_split["test"][0]
    res = surrogate_compensate(s.target, s.setup, gbm_surrogate, iterations=5, tolerance=1.0)
    assert res.stopped == "tolerance" and len(res.history) == 1
    res = surrogate_compensate(s.target, s.setup, gbm_surrogate, iterations=3, stagnation=0.0)
    assert len(res.history) <= 3 and res.best_iteration < len(res.history)
    out = res.save(tmp_path / "c")
    assert (out / "compensated.npz").is_file() and (out / "residual_upper.npz").is_file()
    doc = json.loads((out / "compensation.json").read_text())
    assert doc["model_data_source"] == "proxy - not physics"
    assert doc["predicted_residual"]["quantity"] == "predicted residual (surrogate)"


def test_verification_through_sparlab_form_and_the_api(tmp_path, gbm_surrogate, fake_solver,
                                                       counter):
    """The plumbing of the SparLab check, with the test double (formed = 0.9 x
    commanded): the record names the solver build and the deck, and a bundle
    directory is accepted by precomp.api and the command line."""
    setup = FormingSetup(get_material("AA5754-O"), executable=str(fake_solver), layers=1,
                         element_size=4e-3, clamp_margin=0.02, step_down=1.5e-3,
                         toolpath_spacing=3e-3)
    target = TruncatedCone(0.045, 45.0, 0.015, 0.005, 0.005).heightmap(
        Grid.centered(setup.meshed_blank_size, 4e-3))
    res = surrogate_compensate(target, setup, gbm_surrogate, iterations=3)
    v = verify_with_fea(res, setup, tmp_path / "runs", max_workers=1)
    assert v["data_source"] == "SparLab simulation" and "test double" in v["sparlab_version"]
    assert len(v["deck_hash"]) == 64 and run_count(counter) == 2
    # the double springs back 10 % of the depth, the proxy-trained model
    # predicted otherwise: the check exposes the mismatch instead of echoing
    # the prediction
    assert v["achieved_part_rms_m"] > 5 * v["predicted_part_rms_m"]
    assert v["uncompensated"]["vertical_deviation"]["part"]["rms"] > 0
    d = save_model(gbm_surrogate, tmp_path / "model", {"created_at": ML_CREATED_AT})
    r = api.compensate(target, setup, str(d), method="surrogate", iterations=2,
                       work_dir=tmp_path / "runs")
    assert r.ood is not None and r.interval is not None and r.verification["source"] == "fea"
    p = api.predict(target, setup, d, method="surrogate")
    assert p.details["model_data_source"] == "proxy - not physics" and p.std is not None
