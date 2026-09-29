"""High-level entry points: predict, compensate, compare a scan.

These compose the modules below them and add nothing numerical of their own.
A learned model enters through the `FieldModel` protocol of
`precomp.compensation` (``predict_deviation(commanded, setup) -> (mean,
std)``), optionally extended by

* ``predict_interval(commanded, setup, level) -> (lower, upper)``: bounds of
  the vertical deviation dz [m] at the given coverage level (e.g. from
  conformal calibration);
* ``assess(commanded, setup) -> dict``: an out-of-distribution report.

When the model has them, `compensate` reports the interval and the OOD
assessment; otherwise those fields are None. `precomp.ml` provides such
models (`DeviationSurrogate`); it is imported only when `model` is given as
the directory of a saved bundle, which is then loaded with
`precomp.ml.registry.load_model`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

from ._util import PathLike, PrecompError, to_jsonable, write_json
from .compensation import (CompositePredictor, DAResult, FEAPredictor, FieldModel,
                           SurrogatePredictor, displacement_adjustment)
from .fea.deck import make_toolpath
from .fea.runner import simulate
from .fea.setup import FormingSetup
from .geometry.heightmap import HeightMap
from .metrology import (Alignment, align, flange_mask, metrics, read_point_cloud,
                        region_masks, signed_deviation, vertical_deviation)
from .toolpath import Toolpath

METHODS = ("fea", "surrogate", "hybrid")


@dataclass
class PredictionResult:
    """A predicted formed surface.

    formed : the formed tool-side surface [m] on the commanded grid.
    std : (ny, nx) standard deviation of the prediction [m], when the model
        gives one. method : "fea", "surrogate" or "hybrid".
    details : provenance (the run's hash and cache directory for FEA).
    """

    formed: HeightMap
    std: Optional[np.ndarray]
    method: str
    details: Dict[str, Any] = field(default_factory=dict)


def resolve_model(model: Union[None, FieldModel, PathLike]) -> Optional[FieldModel]:
    """A model object, or a `precomp.ml` bundle directory loaded with
    `precomp.ml.registry.load_model` (imported only then)."""
    if model is None or hasattr(model, "predict_deviation"):
        return model
    if isinstance(model, (str, Path)) or hasattr(model, "__fspath__"):
        from .ml.registry import load_model
        return load_model(model)
    raise TypeError("model must implement predict_deviation(commanded, setup) or be a "
                    "precomp.ml model directory")


def _predictor(setup: FormingSetup, model: Optional[FieldModel], method: str,
               work_dir: Optional[PathLike]):
    if method not in METHODS:
        raise ValueError(f"method must be one of {METHODS}, got {method!r}")
    if method in ("fea", "hybrid") and work_dir is None:
        raise ValueError(f"method {method!r} runs sparlab_form and needs a work_dir")
    if method in ("surrogate", "hybrid") and model is None:
        raise ValueError(f"method {method!r} needs a model")
    if method == "fea":
        return FEAPredictor(setup, work_dir)
    if method == "surrogate":
        return SurrogatePredictor(model, setup)
    return CompositePredictor(FEAPredictor(setup, work_dir), model, setup)


def predict(commanded: HeightMap, setup: FormingSetup, model: Optional[FieldModel] = None, *,
            method: str = "fea", work_dir: Optional[PathLike] = None) -> PredictionResult:
    """Predict the surface formed when `commanded` is followed with `setup`.

    method "fea" simulates with sparlab_form (cached in `work_dir`);
    "surrogate" evaluates `model`; "hybrid" simulates and adds the model's
    learned residual. `model` may be a model object or a `precomp.ml` bundle
    directory; its training data source is reported in `details`.
    """
    model = resolve_model(model)
    pred = _predictor(setup, model, method, work_dir)
    std = None
    if isinstance(pred, SurrogatePredictor):
        formed, std = pred.predict_with_uncertainty(commanded)
    else:
        formed = pred(commanded)
    details: Dict[str, Any] = {}
    if model is not None and method != "fea":
        details["model_data_source"] = getattr(model, "data_source", "unknown")
    fea = pred if isinstance(pred, FEAPredictor) else getattr(pred, "base", None)
    if isinstance(fea, FEAPredictor) and fea.results:
        details["fea"] = {k: v for k, v in fea.results[-1].provenance.items()
                          if isinstance(v, (str, int, float, bool)) or v is None}
    return PredictionResult(formed, std, method, details)


@dataclass
class CompensationResult:
    """What `compensate` produces.

    compensated : the commanded surface to form [m].
    toolpath : its tool path (setup's style and parameters).
    predicted : the predicted formed surface of `compensated` [m].
    interval : (lower, upper) formed surfaces [m] at the model's interval
        level, when the model provides `predict_interval`; else None.
    ood : the model's out-of-distribution assessment, or None.
    verification : when verified by simulation, the metrics [m] of the
        simulated formed surface against the target per region ("all",
        "part", "wall", "flange"), plus the run's content hash; else None.
    history : the displacement-adjustment history (see `DAResult`).
    method : "fea", "surrogate" or "hybrid".
    """

    target: HeightMap
    compensated: HeightMap
    toolpath: Toolpath
    predicted: HeightMap
    interval: Optional[Tuple[HeightMap, HeightMap]]
    ood: Optional[Dict[str, Any]]
    verification: Optional[Dict[str, Any]]
    history: List[Dict[str, Any]]
    method: str

    def save(self, out_dir: PathLike) -> Path:
        """Write compensated.npz, predicted.npz, toolpath.csv (SparLab
        trajectory) and compensation.json (history, verification, OOD)."""
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        self.compensated.save(out / "compensated.npz")
        self.predicted.save(out / "predicted.npz")
        self.toolpath.to_sparlab_csv(out / "toolpath.csv")
        if self.interval is not None:
            self.interval[0].save(out / "predicted_lower.npz")
            self.interval[1].save(out / "predicted_upper.npz")
        write_json(out / "compensation.json", to_jsonable({
            "method": self.method, "history": self.history,
            "verification": self.verification, "ood": self.ood,
            "toolpath": self.toolpath.summary()}))
        return out


def _region_metrics(dev: HeightMap, target: HeightMap,
                    tolerance: Optional[float] = None) -> Dict[str, Dict[str, float]]:
    out = {}
    for name, m in region_masks(target).items():
        sel = m & dev.mask
        if sel.any():
            out[name] = metrics(dev, sel, tolerance)
    return out


def compensate(target: HeightMap, setup: FormingSetup, model: Optional[FieldModel] = None, *,
               method: str = "fea", iterations: int = 3, alpha: float = 1.0,
               verify: bool = True, work_dir: Optional[PathLike] = None,
               smoothing: Optional[float] = None, direction: str = "vertical",
               tolerance: Optional[float] = None, interval_level: float = 0.9) -> CompensationResult:
    """Compensate `target` for springback by displacement adjustment.

    method "fea": every DA iteration is a sparlab_form simulation (cached in
    `work_dir`); the result is the best simulated iterate, which is already
    verified. "surrogate": DA on the model's predictions; with `verify` the
    compensated surface is then simulated once (needs `work_dir` and the
    executable). "hybrid": DA on simulation plus learned residual; `verify`
    simulates as for "surrogate".

    iterations, alpha, smoothing [m], direction, tolerance [m]: as for
    `displacement_adjustment`. interval_level: coverage passed to the model's
    `predict_interval`. `model` may be a model object or a `precomp.ml`
    bundle directory.
    """
    model = resolve_model(model)
    pred = _predictor(setup, model, method, work_dir)
    da: DAResult = displacement_adjustment(target, pred, iterations=iterations, alpha=alpha,
                                           direction=direction, smoothing=smoothing,
                                           tolerance=tolerance)
    comp = da.commanded
    path = make_toolpath(setup, comp)
    interval = None
    ood = None
    if model is not None and method in ("surrogate", "hybrid"):
        if hasattr(model, "predict_interval"):
            # The interval of the learned deviation, placed around the
            # prediction: predicted + (bound - mean). None when the model is
            # not calibrated or has too few calibration parts for the level.
            try:
                lo, hi = model.predict_interval(comp, setup, interval_level)
            except PrecompError:
                lo = hi = None
            if lo is not None:
                mean = np.asarray(model.predict_deviation(comp, setup)[0], float)
                interval = (da.formed.with_z(da.formed.z + np.asarray(lo, float) - mean),
                            da.formed.with_z(da.formed.z + np.asarray(hi, float) - mean))
        if hasattr(model, "assess"):
            ood = dict(model.assess(comp, setup))
    verification = None
    if method == "fea":
        fea_res = pred.results[da.best_iteration]
        dev = signed_deviation(da.formed, target)
        verification = {"source": "fea (best DA iterate)",
                        "key": fea_res.provenance.get("key"),
                        "normal_deviation": _region_metrics(dev, target),
                        "vertical_deviation": _region_metrics(vertical_deviation(da.formed, target),
                                                              target)}
    elif verify:
        if work_dir is None:
            raise ValueError("verify=True needs a work_dir to simulate the compensated part")
        res = simulate(setup, comp, work_dir)
        formed = res.formed_surface(grid=target.grid)
        dev = signed_deviation(formed, target)
        verification = {"source": "fea", "key": res.provenance.get("key"),
                        "normal_deviation": _region_metrics(dev, target),
                        "vertical_deviation": _region_metrics(vertical_deviation(formed, target),
                                                              target)}
    return CompensationResult(target, comp, path, da.formed, interval, ood, verification,
                              da.history, method)


@dataclass
class ScanComparison:
    """A scan compared with its target.

    alignment : the `Alignment` (scan -> target frame), None for mode "none"
        on pre-aligned data. aligned : the aligned points (n, 3) [m].
    formed : the aligned scan gridded on the target grid [m] (invalid where
        the scan has no data). deviation : signed normal deviation [m] on the
        target grid. metrics : per region ("all", "part", "wall", "flange"),
        in metres.
    """

    alignment: Optional[Alignment]
    aligned: np.ndarray
    formed: HeightMap
    deviation: HeightMap
    metrics: Dict[str, Dict[str, float]]
    tolerance: Optional[float] = None

    def summary(self) -> Dict[str, Any]:
        return to_jsonable({"alignment": None if self.alignment is None
                            else self.alignment.summary(),
                            "metrics": self.metrics, "tolerance": self.tolerance,
                            "points": len(self.aligned)})


def compare_scan(scan: Union[PathLike, np.ndarray], target: HeightMap, *, scale: float = 1.0,
                 align_mode: str = "rigid", fixture: str = "all",
                 tolerance: Optional[float] = None, max_gap: Optional[float] = None,
                 method: str = "grid") -> ScanComparison:
    """Align a scan to `target` and measure its signed normal deviation.

    scan : a point-cloud file read with `scale` (e.g. 0.001 for mm), or an
        (n, 3) array, also multiplied by `scale`.
    align_mode : "rigid", "translation_z" or "none". fixture : "all" aligns
        on every point, "flange" on the flange only (as the cell fixtures the
        part; the in-plane position is then unobservable and kept).
    method : "grid" grids the aligned scan on the target grid (nodes farther
        than `max_gap` [m], default 3 spacings, from the scan are left
        without data) and solves the normal intersection on it; "cloud" fits
        local planes to the raw points (`signed_deviation`).
    tolerance [m] adds the share of nodes within +-tolerance to the metrics.
    """
    pts = read_point_cloud(scan, scale=scale) if not isinstance(scan, np.ndarray) \
        else np.asarray(scan, dtype=float) * scale
    if fixture not in ("all", "flange"):
        raise ValueError("fixture must be 'all' or 'flange'")
    al = align(pts, target, align_mode, fixture_mask=flange_mask(target)
               if fixture == "flange" else None)
    aligned = al.apply(pts)
    gap = 3.0 * target.grid.h if max_gap is None else max_gap
    formed = HeightMap.from_points(aligned, target.grid, max_gap=gap,
                                   metadata={"source": "scan"})
    if method == "grid":
        dev = signed_deviation(formed, target)
    elif method == "cloud":
        dev = signed_deviation(aligned, target)
    else:
        raise ValueError("method must be 'grid' or 'cloud'")
    return ScanComparison(al, aligned, formed, dev, _region_metrics(dev, target, tolerance),
                          tolerance)
