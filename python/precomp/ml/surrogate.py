"""The deployable surrogate: features + model + conformal intervals + envelope.

`DeviationSurrogate` wraps a trained point or field model with the feature
configuration it was trained under, its conformal calibrator and its OOD
envelope, and states the data source it learned from. It implements what
`precomp.compensation` and `precomp.api` expect of a model:

* ``predict_deviation(commanded, setup) -> (mean, std)`` - dz [m] at every
  node of the commanded grid (the `FieldModel` protocol);
* ``predict_interval(commanded, setup, level) -> (lower, upper)`` - conformal
  bounds of dz [m];
* ``assess(commanded, setup) -> dict`` - the envelope report.

A point model is trained on the nodes of `domain` ("part" by default); it is
evaluated at every node so the formed surface is defined everywhere, but only
the domain is in-distribution - the flange outside it is an extrapolation used
for display (displacement adjustment holds the flange and measures the error
on the part only).

`train_surrogate` does the whole fit: table, model, envelope, calibration on
held-out whole parts.
"""

from __future__ import annotations

import hashlib
from collections import OrderedDict
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from .._util import PrecompError, canonical_json
from ..geometry.heightmap import HeightMap
from .dataset import Sample, build_table, source_label
from .features import (DEFAULT_CONFIG, PRIOR_FEATURE, FeatureConfig, FeatureMaps, as_setup,
                       feature_maps, global_features)
from .models import ResidualModel, TransferModel, model_from_state
from .ood import OODEnvelope
from .uncertainty import ConformalCalibrator

SURROGATE_FORMAT = "precomp.ml.DeviationSurrogate/1"


def _key(commanded: HeightMap, setup: Any) -> str:
    h = hashlib.sha256()
    h.update(np.ascontiguousarray(commanded.z).tobytes())
    h.update(np.ascontiguousarray(commanded.mask).tobytes())
    h.update(canonical_json(commanded.grid.to_dict()).encode())
    h.update(canonical_json(as_setup(setup).physics_dict()).encode())
    return h.hexdigest()


def model_prior(model: Any) -> Any:
    """The physics prior a model needs as its last column, or None."""
    if isinstance(model, ResidualModel):
        return model.prior
    if isinstance(model, TransferModel):
        return model_prior(model.base)
    return None


class DeviationSurrogate:
    """A trained deviation model ready for prediction (see the module docstring).

    model : a fitted point model (`models.DeviationModel`) or field model.
    features : the FeatureConfig it was trained with.
    data_source : the label of its training data ("SparLab simulation",
        "proxy - not physics", "scan" or "mixed: ...").
    training : provenance of the fit (sample ids, counts, families, seed ...).
    calibrator, ood : optional fitted ConformalCalibrator / OODEnvelope.
    domain : the training mask ("part" or "all").
    """

    def __init__(self, model: Any, features: FeatureConfig = DEFAULT_CONFIG, *,
                 data_source: str, training: Optional[Dict[str, Any]] = None,
                 calibrator: Optional[ConformalCalibrator] = None,
                 ood: Optional[OODEnvelope] = None, domain: str = "part"):
        if getattr(model, "kind", None) not in ("point", "field"):
            raise TypeError("model must be a point or field model of precomp.ml.models")
        self.model = model
        self.features = features
        self.data_source = data_source
        self.training = dict(training or {})
        self.calibrator = calibrator
        self.ood = ood
        self.domain = domain
        #: the manifest when loaded from a bundle (`registry.load_model`)
        self.manifest: Optional[Dict[str, Any]] = None
        self._cache: "OrderedDict[str, Any]" = OrderedDict()

    # -- state -----------------------------------------------------------------
    @property
    def kind(self) -> str:
        return self.model.kind

    @property
    def feature_names(self) -> List[str]:
        names = self.features.names
        if model_prior(self.model) is not None:
            names = names + [PRIOR_FEATURE]
        return names

    @property
    def model_class(self) -> str:
        return type(self.model).__name__

    def to_state(self) -> Dict[str, Any]:
        return {"format": SURROGATE_FORMAT, "model": self.model.to_state(),
                "features": self.features.to_dict(), "data_source": self.data_source,
                "training": self.training, "domain": self.domain,
                "calibrator": None if self.calibrator is None else self.calibrator.to_state(),
                "ood": None if self.ood is None else self.ood.to_state()}

    @classmethod
    def from_state(cls, state: Mapping[str, Any]) -> "DeviationSurrogate":
        if state.get("format") != SURROGATE_FORMAT:
            raise PrecompError(f"not a surrogate state (format {state.get('format')!r})")
        return cls(model_from_state(state["model"]), FeatureConfig.from_dict(state["features"]),
                   data_source=state["data_source"], training=state.get("training"),
                   calibrator=None if state.get("calibrator") is None
                   else ConformalCalibrator.from_state(state["calibrator"]),
                   ood=None if state.get("ood") is None
                   else OODEnvelope.from_state(state["ood"]),
                   domain=state.get("domain", "part"))

    # -- evaluation ------------------------------------------------------------
    #: Cached entries per kind: feature maps are large (ny x nx x F), the
    #: predictions and envelope reports small. The key is a hash of the
    #: commanded surface and the setup's physics.
    CACHE_SIZES = {"maps": 4, "pred": 128, "assess": 128}

    def _cached(self, key: str, make):
        kind = key.split(":", 1)[0]
        if key in self._cache:
            self._cache.move_to_end(key)
            return self._cache[key]
        val = make()
        self._cache[key] = val
        keys = [k for k in self._cache if k.startswith(kind + ":")]
        for k in keys[:max(0, len(keys) - self.CACHE_SIZES.get(kind, 8))]:
            del self._cache[k]
        return val

    def feature_maps(self, commanded: HeightMap, setup: Any) -> FeatureMaps:
        return self._cached("maps:" + _key(commanded, setup),
                            lambda: feature_maps(commanded, setup, None, self.features))

    def _point_matrix(self, commanded: HeightMap, setup: Any) -> np.ndarray:
        X = self.feature_maps(commanded, setup).gather(None)
        prior = model_prior(self.model)
        if prior is not None:
            p = np.asarray(prior.prior_deviation(commanded, as_setup(setup)), float)
            if p.shape != commanded.grid.shape or not np.all(np.isfinite(p)):
                raise PrecompError("the model's prior returned a field of the wrong shape or "
                                   "non-finite values")
            X = np.column_stack([X, p.ravel()])
        return X

    def predict_deviation(self, commanded: HeightMap, setup: Any
                          ) -> Tuple[np.ndarray, np.ndarray]:
        """(mean, std) of dz [m], (ny, nx), at every node of the commanded grid."""
        def make():
            if self.kind == "field":
                mu, sd = self.model.predict_field(commanded, setup)
                return np.asarray(mu, float), np.asarray(sd, float)
            X = self._point_matrix(commanded, setup)
            mu, sd = self.model.predict(X)
            shape = commanded.grid.shape
            return np.asarray(mu, float).reshape(shape), np.asarray(sd, float).reshape(shape)
        mu, sd = self._cached("pred:" + _key(commanded, setup), make)
        return mu.copy(), sd.copy()

    def predict_interval(self, commanded: HeightMap, setup: Any, level: float = 0.9
                         ) -> Tuple[np.ndarray, np.ndarray]:
        """Conformal (lower, upper) bounds of dz [m] at coverage `level`."""
        if self.calibrator is None or not self.calibrator.fitted:
            raise PrecompError("this surrogate has no conformal calibration; train it with "
                               "calibration samples")
        mu, sd = self.predict_deviation(commanded, setup)
        region = self.feature_maps(commanded, setup).region if self.calibrator.mondrian \
            else None
        return self.calibrator.intervals(mu, sd, level, region)

    def assess(self, commanded: HeightMap, setup: Any, *, max_points: int = 1500
               ) -> Dict[str, Any]:
        """The envelope report of a commanded part (see `ood.OODEnvelope`) on
        at most `max_points` part nodes (evenly spread), with the data source
        the model learned from."""
        if self.ood is None:
            raise PrecompError("this surrogate has no OOD envelope")

        def make() -> Dict[str, Any]:
            fm = self.feature_maps(commanded, setup)
            nodes = np.flatnonzero(fm.part.ravel())
            if nodes.size > max_points:
                nodes = nodes[np.linspace(0, nodes.size - 1, max_points).round().astype(int)]
            d, _ = global_features(commanded, setup, self.features)
            rep = self.ood.assess_features(fm.gather(nodes), d)
            rep["model_data_source"] = self.data_source
            if self.kind == "field" and hasattr(self.model, "fraction_outside_window"):
                rep["fraction_outside_window"] = self.model.fraction_outside_window(
                    commanded.grid)
            return rep
        import copy
        return copy.deepcopy(self._cached(f"assess:{max_points}:" + _key(commanded, setup),
                                          make))

    def describe(self) -> Dict[str, Any]:
        return {"model_class": self.model_class, "kind": self.kind,
                "data_source": self.data_source, "domain": self.domain,
                "features": self.features.to_dict(), "training": self.training,
                "calibration": None if self.calibrator is None else self.calibrator.summary(),
                "ood": None if self.ood is None else self.ood.summary()}


def fit_envelope(samples: Sequence[Sample], features: FeatureConfig, *,
                 points_per_sample: int = 400, seed: int = 0, n_jobs: int = 1,
                 **kwargs: Any) -> OODEnvelope:
    """An OOD envelope from the training samples (their part nodes and
    descriptors, grouped by part id)."""
    table = build_table(samples, features, points_per_sample=points_per_sample, mask="part",
                        rng=seed + 101, stratify=True, n_jobs=n_jobs)
    D, names = [], None
    for s in samples:
        d, names = global_features(s.commanded, s.setup, features)
        D.append(d)
    return OODEnvelope(seed=seed, **kwargs).fit(
        table.X, table.part_id, np.stack(D), np.array([s.part_id for s in samples]),
        point_names=table.feature_names, part_names=names)


def calibrate(surrogate: DeviationSurrogate, samples: Sequence[Sample], *,
              points_per_sample: Optional[int] = None, mondrian: bool = False,
              seed: int = 0) -> ConformalCalibrator:
    """Fit a ConformalCalibrator on held-out samples: every part node (or a
    uniform subsample of `points_per_sample` per part) of each sample."""
    mus, sds, ys, gs, rs = [], [], [], [], []
    rng = np.random.default_rng(seed)
    for s in samples:
        mu, sd = surrogate.predict_deviation(s.commanded, s.setup)
        fm = surrogate.feature_maps(s.commanded, s.setup)
        sel = np.flatnonzero((fm.part & s.valid).ravel())
        if points_per_sample is not None and sel.size > points_per_sample:
            sel = np.sort(rng.choice(sel, size=points_per_sample, replace=False))
        mus.append(mu.ravel()[sel])
        sds.append(sd.ravel()[sel])
        ys.append(s.dz.ravel()[sel])
        gs.append(np.full(sel.size, s.sample_id, dtype=object))
        rs.append(fm.region.ravel()[sel])
    cal = ConformalCalibrator(mondrian=mondrian, seed=seed)
    return cal.fit(np.concatenate(mus), np.concatenate(sds), np.concatenate(ys),
                   np.concatenate(gs), region=np.concatenate(rs),
                   data_source=source_label(s.source for s in samples))


def _ids_hash(ids: Sequence[str]) -> str:
    return hashlib.sha256("\n".join(sorted(ids)).encode()).hexdigest()


def train_surrogate(model: Any, train: Sequence[Sample], *,
                    calibration: Optional[Sequence[Sample]] = None,
                    features: FeatureConfig = DEFAULT_CONFIG,
                    points_per_sample: Optional[int] = 2000, mask: str = "part",
                    seed: int = 0, envelope: bool = True, mondrian: bool = False,
                    calibration_points: Optional[int] = None, n_jobs: int = 1
                    ) -> DeviationSurrogate:
    """Fit `model` on `train`, then the envelope, then calibrate on `calibration`.

    Point models get a stratified table (`points_per_sample` per sample, the
    physics prior's column appended for residual models); field models are
    fitted on the whole samples. `calibration` must hold parts that are not
    in `train` (checked by part id); without it the surrogate has no
    intervals. The training provenance - sample ids (and their hash), part
    ids, families, sources, the data-source label, seed and sizes - is kept
    in `surrogate.training`.
    """
    train = list(train)
    if not train:
        raise ValueError("no training samples")
    cal = list(calibration or [])
    overlap = {s.part_id for s in train} & {s.part_id for s in cal}
    if overlap:
        raise PrecompError(f"{len(overlap)} part(s) are in both the training and the "
                           f"calibration set, e.g. {sorted(overlap)[:3]}; calibrate on held-out "
                           "parts")
    label = source_label(s.source for s in train)
    info: Dict[str, Any] = {
        "sample_ids": sorted(s.sample_id for s in train),
        "sample_ids_sha256": _ids_hash([s.sample_id for s in train]),
        "n_samples": len(train), "n_parts": len({s.part_id for s in train}),
        "families": sorted({s.family for s in train}),
        "sources": sorted({s.source for s in train}),
        "fidelities": sorted({s.fidelity for s in train}), "data_source": label,
        "sparlab_versions": sorted({str(s.provenance["sparlab_version"]) for s in train
                                    if s.provenance.get("sparlab_version")}),
        "seed": int(seed), "mask": mask, "points_per_sample": points_per_sample,
        "calibration_sample_ids": sorted(s.sample_id for s in cal)}
    if getattr(model, "kind", None) == "field":
        model.fit_samples(train)
    else:
        table = build_table(train, features, points_per_sample=points_per_sample, mask=mask,
                            rng=seed, stratify=True, prior=model_prior(model), n_jobs=n_jobs)
        model.fit(table.X, table.y, table.groups, feature_names=table.feature_names)
        info["n_rows"] = len(table)
    sur = DeviationSurrogate(model, features, data_source=label, training=info, domain=mask)
    if envelope:
        sur.ood = fit_envelope(train, features, seed=seed, n_jobs=n_jobs)
    if cal:
        sur.calibrator = calibrate(sur, cal, points_per_sample=calibration_points,
                                   mondrian=mondrian, seed=seed)
    return sur


def transfer_surrogate(base: DeviationSurrogate, real: Sequence[Sample], *,
                       calibration: Optional[Sequence[Sample]] = None,
                       correction: str = "bayesian_ridge",
                       points_per_sample: Optional[int] = 2000, seed: int = 0,
                       mondrian: bool = False) -> DeviationSurrogate:
    """A `TransferModel` surrogate: `base`'s model frozen, corrected on the few
    `real` samples (scans, or a higher-fidelity simulation).

    The envelope stays the base's (the correction does not widen what the
    base has seen). Intervals are calibrated only on `calibration` - held-out
    real parts; without them the transfer surrogate has no intervals rather
    than intervals calibrated on another data source. With no real samples
    the model equals the base.
    """
    if base.kind != "point":
        raise TypeError("transfer learning is implemented for point models")
    real = list(real)
    model = TransferModel(base.model, correction, seed=seed)
    info = dict(base.training)
    info.update(transfer={"n_real": 0, "real_sample_ids": [], "correction": correction})
    if real:
        table = build_table(real, base.features, points_per_sample=points_per_sample,
                            mask=base.domain, rng=seed, prior=model_prior(base.model))
        model.fit(table.X, table.y, table.groups, feature_names=table.feature_names)
        info["transfer"].update(n_real=model.n_real,
                                real_sample_ids=sorted(s.sample_id for s in real),
                                real_sources=sorted({s.source for s in real}))
    sources = set(base.training.get("sources", [])) | {s.source for s in real}
    sur = DeviationSurrogate(model, base.features, data_source=source_label(sources),
                             training=info, ood=base.ood, domain=base.domain)
    cal = list(calibration or [])
    if cal:
        overlap = {s.part_id for s in real} & {s.part_id for s in cal}
        if overlap:
            raise PrecompError(f"calibration parts {sorted(overlap)[:3]} were used for the "
                               "correction")
        sur.calibrator = calibrate(sur, cal, mondrian=mondrian, seed=seed)
    return sur


__all__ = ["DeviationSurrogate", "train_surrogate", "transfer_surrogate", "calibrate",
           "fit_envelope", "model_prior", "SURROGATE_FORMAT"]
