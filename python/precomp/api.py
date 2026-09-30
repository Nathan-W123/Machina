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
assessment; otherwise those fields are None. A model with `assess` is held
to its envelope: `compensate` refuses a target, or a compensated shape,
outside it unless `allow_out_of_envelope`. `precomp.ml` provides such models
(`DeviationSurrogate`); it is imported only when `model` is given as the
directory of a saved bundle, which is then loaded with
`precomp.ml.registry.load_model`.

A model may declare what it predicts in `target`: "dz", the total vertical
deviation of the formed surface from the commanded one (every precomp.ml
bundle), or "residual", a correction of a simulation. Method "hybrid" adds
the model to a simulation, so it refuses a "dz" model - that would count
the springback twice; models without the attribute are taken as residuals.
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
from .fea.support import command_upper_bound, compensation_masks
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
               work_dir: Optional[PathLike], target: Optional[HeightMap] = None):
    if method not in METHODS:
        raise ValueError(f"method must be one of {METHODS}, got {method!r}")
    if method in ("fea", "hybrid") and work_dir is None:
        raise ValueError(f"method {method!r} runs sparlab_form and needs a work_dir")
    if method in ("surrogate", "hybrid") and model is None:
        raise ValueError(f"method {method!r} needs a model")
    if method == "fea":
        return FEAPredictor(setup, work_dir, target=target)
    if method == "surrogate":
        return SurrogatePredictor(model, setup)
    if getattr(model, "target", "residual") == "dz":
        raise ValueError(
            "method 'hybrid' adds the model's prediction to a simulation, but this model "
            "predicts the total deviation dz (formed - commanded), not a residual over a "
            "simulation: the springback would be counted twice. Use method 'surrogate' (a "
            "precomp.ml ResidualModel runs its FE prior itself)")
    return CompositePredictor(FEAPredictor(setup, work_dir, target=target), model, setup)


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
    ood : the model's out-of-distribution assessment of `compensated`, or None.
    verification : when verified by simulation, the metrics [m] of the
        simulated formed surface against the target per region ("all",
        "part", "wall", "flange"), plus the run's content hash; else None.
    history : the displacement-adjustment history (see `DAResult`) - for
        "surrogate" and "hybrid", predictions of the model.
    method : "fea", "surrogate" or "hybrid".
    ood_target : the assessment of the target itself, or None.
    verified : the simulated formed surface of `compensated` when verified.
    model_data_source : what the model learned from (None for "fea").
    allowed_out_of_envelope : whether the envelope was overridden.
    interval_note : why there is no interval although the model has
        `predict_interval` (e.g. too few calibration parts for the level).
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
    ood_target: Optional[Dict[str, Any]] = None
    verified: Optional[HeightMap] = None
    model_data_source: Optional[str] = None
    allowed_out_of_envelope: bool = False
    interval_note: Optional[str] = None

    @property
    def in_envelope(self) -> Optional[bool]:
        """Target and compensated shape both inside the envelope (None without one)."""
        if self.ood is None:
            return None
        return bool(self.ood["in_envelope"] and (self.ood_target is None
                                                 or self.ood_target["in_envelope"]))

    def save(self, out_dir: PathLike) -> Path:
        """Write compensated.npz, predicted.npz (and verified.npz when
        verified), toolpath.csv (SparLab trajectory) and compensation.json
        (history, verification, envelope, with what each number is)."""
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        self.compensated.save(out / "compensated.npz")
        self.predicted.save(out / "predicted.npz")
        if self.verified is not None:
            self.verified.save(out / "verified.npz")
        self.toolpath.to_sparlab_csv(out / "toolpath.csv")
        if self.interval is not None:
            self.interval[0].save(out / "predicted_lower.npz")
            self.interval[1].save(out / "predicted_upper.npz")
        what = ("simulations (sparlab_form)" if self.method == "fea" else
                f"predictions of the model (data source: {self.model_data_source}), "
                "not simulations")
        write_json(out / "compensation.json", to_jsonable({
            "method": self.method, "model_data_source": self.model_data_source,
            "history_quantity": what, "history": self.history,
            "verification": self.verification, "in_envelope": self.in_envelope,
            "allowed_out_of_envelope": self.allowed_out_of_envelope,
            "ood": self.ood, "ood_target": self.ood_target,
            "interval_note": self.interval_note, "toolpath": self.toolpath.summary()}))
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
               tolerance: Optional[float] = None, interval_level: float = 0.9,
               allow_out_of_envelope: bool = False) -> CompensationResult:
    """Compensate `target` for springback by displacement adjustment.

    method "fea": every DA iteration is a sparlab_form simulation (cached in
    `work_dir`); the result is the best simulated iterate, which is already
    verified. "surrogate": DA on the model's predictions; with `verify` the
    compensated surface is then simulated once (needs `work_dir` and the
    executable). "hybrid": DA on simulation plus learned residual; `verify`
    simulates as for "surrogate".

    iterations, alpha, smoothing [m], direction, tolerance [m]: as for
    `displacement_adjustment`. A setup with a support strategy
    (`FormingSetup.support`) forms every command with its fixture made for
    `target`, and lets the command rise above the sheet plane where the
    support can push it up (`precomp.fea.support.command_upper_bound`).
    interval_level: coverage passed to the model's
    `predict_interval`. `model` may be a model object or a `precomp.ml`
    bundle directory.

    allow_out_of_envelope : with a model that has `assess` (an envelope), a
        target outside it, or a compensated shape outside it, raises
        PrecompError unless this is True; the verdicts are in `ood_target`,
        `ood` and each history entry either way. (`precomp.ml.compensate`
        stops at the first iterate outside instead.)
    """
    model = resolve_model(model)
    pred = _predictor(setup, model, method, work_dir, target)
    upper = command_upper_bound(setup, target) if setup.support != "none" else None
    hold, adjust = compensation_masks(setup, target)
    learned = model is not None and method in ("surrogate", "hybrid")
    assess = learned and hasattr(model, "assess") and getattr(model, "ood", True) is not None
    ood_target = None
    if assess:
        ood_target = dict(model.assess(target, setup))
        if not ood_target["in_envelope"] and not allow_out_of_envelope:
            raise PrecompError(
                "the target is outside the model's training envelope (reasons: "
                f"{', '.join(ood_target.get('reasons', [])[:5]) or 'not given'}): the "
                "model would extrapolate. Pass allow_out_of_envelope=True "
                "(--allow-out-of-envelope) to compensate anyway")
    verdicts: List[Dict[str, Any]] = []

    def record(k: int, c: HeightMap, f: HeightMap, e: HeightMap) -> None:
        if assess:
            verdicts.append(dict(model.assess(c, setup)))

    da: DAResult = displacement_adjustment(target, pred, iterations=iterations, alpha=alpha,
                                           direction=direction, smoothing=smoothing,
                                           tolerance=tolerance, callback=record,
                                           upper_bound=upper, hold_mask=hold,
                                           adjust_mask=adjust)
    comp = da.commanded
    history = [dict(h) for h in da.history]
    for h, v in zip(history, verdicts):
        h["in_envelope"] = bool(v["in_envelope"])
        h["ood_part_score"] = v.get("part_score")
    if learned:
        for h in history:
            h["predicted"] = True
    path = make_toolpath(setup, comp)
    interval = None
    interval_note = None
    ood = verdicts[da.best_iteration] if assess else None
    if ood is not None and not ood["in_envelope"] and not allow_out_of_envelope:
        raise PrecompError(
            f"the compensated shape (iteration {da.best_iteration}) is outside the model's "
            f"training envelope (reasons: {', '.join(ood.get('reasons', [])[:5]) or 'not given'}"
            "): its prediction is an extrapolation. Use fewer iterations, or pass "
            "allow_out_of_envelope=True (--allow-out-of-envelope)")
    if learned and hasattr(model, "predict_interval"):
        # The interval of the learned deviation, placed around the
        # prediction: predicted + (bound - mean).
        try:
            lo, hi = model.predict_interval(comp, setup, interval_level)
        except PrecompError as exc:     # not calibrated, or too few parts for the level
            interval_note = str(exc)
        else:
            mean = np.asarray(model.predict_deviation(comp, setup)[0], float)
            interval = (da.formed.with_z(da.formed.z + np.asarray(lo, float) - mean),
                        da.formed.with_z(da.formed.z + np.asarray(hi, float) - mean))
    verification = None
    verified = None
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
        res = simulate(setup, comp, work_dir, target=target)
        verified = res.formed_surface(grid=target.grid)
        dev = signed_deviation(verified, target)
        verification = {"source": "fea", "key": res.provenance.get("key"),
                        "normal_deviation": _region_metrics(dev, target),
                        "vertical_deviation": _region_metrics(vertical_deviation(verified,
                                                                                 target),
                                                              target)}
    return CompensationResult(target, comp, path, da.formed, interval, ood, verification,
                              history, method, ood_target=ood_target, verified=verified,
                              model_data_source=getattr(model, "data_source", "unknown")
                              if learned else None,
                              allowed_out_of_envelope=bool(allow_out_of_envelope),
                              interval_note=interval_note)


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
