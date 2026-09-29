"""precomp.ml models, conformal intervals and the envelope, on ProxySimulator data.

Every number here is measured on proxy data - an analytic stand-in, NOT
physics. The shared surrogate (conftest `gbm_surrogate`) is trained on the
in-distribution families without the held-out one.
"""

import numpy as np
import pytest

from conftest import ML_CREATED_AT
from precomp.ml import (GBMEnsemble, ProxySimulator, ResidualModel, build_table,
                        evaluate_surrogate, simulate_samples, train_surrogate, transfer_surrogate)
from precomp.ml.generate import DesignSpace, design_points
from precomp.ml.uncertainty import ConformalCalibrator, partition_coverage


def _held_out_arrays(sur, samples):
    mus, sds, ys, gs = [], [], [], []
    for s in samples:
        mu, sd = sur.predict_deviation(s.commanded, s.setup)
        sel = sur.feature_maps(s.commanded, s.setup).part & s.valid
        mus.append(mu[sel])
        sds.append(sd[sel])
        ys.append(s.dz[sel])
        gs.append(np.full(sel.sum(), s.part_id, dtype=object))
    return [np.concatenate(v) for v in (mus, sds, ys, gs)]


def test_gbm_ensemble_learns_the_proxy_on_held_out_parts(gbm_surrogate, proxy_split):
    rep = evaluate_surrogate(gbm_surrogate, proxy_split["test"], level=0.9, levels=(0.9,),
                             assess=False)
    pp = rep.per_part
    assert set(pp["data_source"]) == {"proxy - not physics"}
    assert pp["part_id"].nunique() >= 10
    # dz is ~0.5-1 mm; the held-out error a small fraction of it
    assert pp["dz_rms"].mean() > 4e-4
    assert pp["rel_rms"].mean() < 0.2, pp[["sample_id", "rel_rms"]]
    assert pp["err_rms"].mean() < 1.2e-4
    assert np.all(pp["std_mean_m"] > 0)


def test_conformal_coverage_is_nominal_on_held_out_parts(gbm_surrogate, proxy_split):
    """Expected coverage of grouped split-conformal intervals, over random
    calibration/test partitions of held-out parts, within 3 % of nominal.
    (A single draw of ~12 calibration parts lands several per cent away: the
    errors of one part are correlated, so parts, not nodes, count.)"""
    mu, sd, y, g = _held_out_arrays(gbm_surrogate, proxy_split["pool"])
    for level in (0.8, 0.9):
        pc = partition_coverage(mu, sd, y, g, level=level, n_partitions=60, seed=3)
        assert abs(pc["mean"] - level) <= 0.03, pc
        assert pc["n_parts"] >= 20
    # the finite-sample quantile: the ceil((n + 1) level)-th smallest score
    cal = ConformalCalibrator(eps=1.0).fit(np.zeros(9), np.zeros(9), np.arange(1.0, 10.0),
                                           np.arange(9))
    assert cal.quantile(0.8) == pytest.approx(8.0 / 1.0)     # k = ceil(10 * 0.8) = 8
    with pytest.raises(Exception, match="cannot support"):
        cal.quantile(0.95)                                    # k = 10 > n = 9


def test_the_envelope_flags_an_unseen_family_more_than_known_ones(gbm_surrogate, proxy_split):
    known = [gbm_surrogate.assess(s.commanded, s.setup) for s in proxy_split["test"]]
    unseen = [gbm_surrogate.assess(s.commanded, s.setup) for s in proxy_split["family"]]
    flag_known = np.mean([not a["in_envelope"] for a in known])
    flag_unseen = np.mean([not a["in_envelope"] for a in unseen])
    assert flag_known <= 0.25 and flag_unseen >= flag_known + 0.5, (flag_known, flag_unseen)
    assert np.median([a["part_score"] for a in unseen]) > np.median(
        [a["part_score"] for a in known])
    out = [a for a in unseen if not a["in_envelope"]]
    assert all(a["reasons"] for a in out)                     # it says why
    assert known[0]["model_data_source"] == "proxy - not physics"
    s = proxy_split["test"][0]
    assert gbm_surrogate.assess(s)["part_score"] == known[0]["part_score"]   # a Sample too


@pytest.fixture(scope="module")
def shifted_samples(ml_threads):
    """A 'real' process that differs from the simulated one: rim under-forming
    +40 %, unclamping curvature +30 % (still proxy data)."""
    real = ProxySimulator().shifted(rim_amplitude=1.12e-3, curvature=0.13)
    space = DesignSpace(grid_spacing=4e-3, materials=("AA5754-O",), process={},
                        families=("truncated_cone", "pyramid", "dome", "two_level"))
    pts = design_points(space, 3, seed=29)
    return simulate_samples(pts, real, created_at=ML_CREATED_AT, kinds=("uncompensated",))


def test_transfer_model_equals_its_base_without_real_data_and_improves_with_some(
        gbm_surrogate, shifted_samples, ml_features):
    base = gbm_surrogate
    t0 = transfer_surrogate(base, [])
    s = shifted_samples[0]
    a, b = base.predict_deviation(s.commanded, s.setup), t0.predict_deviation(s.commanded,
                                                                              s.setup)
    assert np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])
    assert t0.model.n_real == 0
    real, test = shifted_samples[:4], shifted_samples[4:]
    t4 = transfer_surrogate(base, real, points_per_sample=400, seed=0)
    assert t4.model.n_real == 4 and t4.training["transfer"]["n_real"] == 4
    tb = build_table(test, ml_features, points_per_sample=None)
    err_base = np.sqrt(np.mean((base.model.predict(tb.X)[0] - tb.y) ** 2))
    mu, sd = t4.model.predict(tb.X)
    err_tr = np.sqrt(np.mean((mu - tb.y) ** 2))
    assert err_tr < 0.7 * err_base, (err_tr, err_base)
    assert np.all(sd >= base.model.predict(tb.X)[1] - 1e-15)  # uncertainty only grows
    with pytest.raises(Exception, match="no conformal"):
        t4.predict_interval(s.commanded, s.setup, 0.9)        # not calibrated on real data


def test_residual_model_learns_what_its_prior_misses(proxy_split, shifted_samples, ml_features):
    """A proxy prior of the wrong strength plus a GBM on the residual: the
    hybrid path. On the shifted data it beats the prior alone by far."""
    prior = ProxySimulator()
    model = ResidualModel(prior, GBMEnsemble(2, max_iter=60, target_scale=None))
    sur = train_surrogate(model, shifted_samples[:8], features=ml_features,
                          points_per_sample=300, envelope=False)
    assert sur.feature_names[-1] == "prior_dz"
    err_prior, err_hyb = [], []
    for s in shifted_samples[8:]:
        sel = s.commanded.z < -1e-6
        mu, _ = sur.predict_deviation(s.commanded, s.setup)
        err_hyb.append(np.sqrt(np.mean((mu - s.dz)[sel] ** 2)))
        err_prior.append(np.sqrt(np.mean((prior.deviation(s.commanded, s.setup) - s.dz)[sel]
                                         ** 2)))
    assert np.mean(err_hyb) < 0.5 * np.mean(err_prior)


def test_gbm_is_deterministic_and_validates_its_input():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(600, 3))
    y = X[:, 0] + 0.1 * rng.normal(size=600)
    g = np.repeat(np.arange(6), 100)
    a = GBMEnsemble(3, max_iter=30, target_scale=None, seed=5).fit(X, y, g)
    b = GBMEnsemble(3, max_iter=30, target_scale=None, seed=5).fit(X, y, g)
    assert np.array_equal(a.predict(X)[0], b.predict(X)[0])
    with pytest.raises(ValueError, match="target_scale"):
        GBMEnsemble(2, max_iter=5).fit(X, y, g)               # no elastic_curvature column
    with pytest.raises(ValueError, match="NaN"):
        a.predict(np.full((2, 3), np.nan))
    with pytest.raises(ValueError, match="3 features"):
        a.predict(X[:, :2])
    q = GBMEnsemble(2, max_iter=30, target_scale=None, quantiles=(0.05, 0.95)).fit(X, y, g)
    band = q.predict_quantiles(X)
    assert np.mean((y >= band[0.05]) & (y <= band[0.95])) > 0.8


def test_mlp_ensemble_fits_and_is_reproducible(proxy_split, ml_features):
    pytest.importorskip("torch")
    from precomp.ml import MLPEnsemble

    tr = build_table(proxy_split["train"][:24], ml_features, points_per_sample=200, rng=0)
    te = build_table(proxy_split["test"][:6], ml_features, points_per_sample=200, rng=1)
    kw = dict(hidden=(32, 32), epochs=8, batch_size=512, seed=3)
    m = MLPEnsemble(2, **kw).fit(tr.X, tr.y, tr.groups, feature_names=tr.feature_names)
    mu, sd = m.predict(te.X)
    assert np.all(np.isfinite(mu)) and np.all(sd > 0)
    err = np.sqrt(np.mean((mu - te.y) ** 2))
    assert err < 0.5 * np.sqrt(np.mean((te.y - tr.y.mean()) ** 2))
    m2 = MLPEnsemble(2, **kw).fit(tr.X, tr.y, tr.groups, feature_names=tr.feature_names)
    assert np.array_equal(m2.predict(te.X)[0], mu)


def test_field_unet_overfits_one_sample(proxy_split):
    pytest.importorskip("torch")
    from precomp.ml import FieldUNet
    from precomp.ml.features import FeatureConfig

    s = proxy_split["train"][0]
    net = FieldUNet(resolution=24, channels=(8, 16), dropout=0.0, epochs=150, mse_warmup=150,
                    lr=5e-3, mc_samples=1, features=FeatureConfig(time_source="depth"), seed=0)
    net.fit_samples([s])
    mu, sd = net.predict_field(s.commanded, s.setup)
    sel = s.commanded.z < -1e-6
    err = np.sqrt(np.mean((mu - s.dz)[sel] ** 2))
    assert mu.shape == s.grid.shape and np.all(np.isfinite(sd))
    assert err < 0.2 * np.sqrt(np.mean(s.dz[sel] ** 2)), err
    assert net.history[-1]["loss"] < 0.05 * net.history[0]["loss"]
    with pytest.raises(TypeError, match="fit_samples"):
        net.fit(np.zeros((2, 2)), np.zeros(2), np.zeros(2))
