"""Compensation on a learned surrogate, and its honest end-to-end check.

`surrogate_compensate` runs displacement adjustment (`precomp.compensation`)
with the surrogate as the predictor of the formed surface,

    formed(c) = c + dz_hat(c),

one iteration at a time so that it can stop when the PREDICTED deviation from
the target is within tolerance or stops improving (stagnation). It returns the
compensated commanded surface (the best predicted iterate), the predicted
formed surface, the predicted residual deviation from the target with its
conformal interval, the out-of-distribution report of the compensated shape,
and the history - all predictions, labelled with the data source the model
learned from.

A prediction is not a result. `verify_with_fea` simulates the compensated
shape with sparlab_form and reports the achieved deviation from the target
(and, optionally, that of the uncompensated target for the improvement
factor): the check to run before trusting a surrogate compensation.
`verify_with_simulator` does the same through any `Simulator`, e.g. the
ProxySimulator for CI - labelled "proxy - not physics".
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from .._util import PathLike, PrecompError, to_jsonable, write_json
from ..compensation import SurrogatePredictor as _BaseSurrogatePredictor
from ..compensation import displacement_adjustment
from ..geometry.heightmap import HeightMap
from ..geometry.parts import MAX_WALL_ANGLE_DEG
from ..metrology import metrics, part_mask, region_masks, signed_deviation, vertical_deviation
from .dataset import SOURCE_LABELS
from .features import as_setup


class SurrogatePredictor(_BaseSurrogatePredictor):
    """`precomp.compensation.SurrogatePredictor` that keeps the std map of its
    last prediction (`last_std`, [m]) and counts its calls."""

    def __init__(self, model: Any, setup: Any):
        super().__init__(model, setup)
        self.last_std: Optional[np.ndarray] = None
        self.calls = 0

    def __call__(self, commanded: HeightMap) -> HeightMap:
        formed, std = self.predict_with_uncertainty(commanded)
        self.last_std = std
        self.calls += 1
        return formed


@dataclass
class SurrogateCompensation:
    """The outcome of `surrogate_compensate` (all [m]; predictions, not results).

    compensated : the commanded surface to form. predicted : its predicted
    formed surface. residual : predicted - target (vertical) on the target
    grid. residual_interval : (lower, upper) of the residual at `level`, or
    None without calibration. std : the model's std map at `compensated`.
    ood : the envelope report of `compensated` (and `ood_target` of the
    uncompensated target). history : per iteration, the predicted error over
    the part and the update size. stopped : "tolerance", "stagnation" or
    "iterations". data_source : what the model learned from.
    verification : filled by `verify_with_fea` / `verify_with_simulator`.
    """

    target: HeightMap
    compensated: HeightMap
    predicted: HeightMap
    residual: HeightMap
    residual_interval: Optional[Tuple[HeightMap, HeightMap]]
    std: Optional[np.ndarray]
    level: float
    ood: Optional[Dict[str, Any]]
    ood_target: Optional[Dict[str, Any]]
    history: List[Dict[str, Any]]
    stopped: str
    best_iteration: int
    data_source: str
    verification: Optional[Dict[str, Any]] = None

    def predicted_metrics(self, tolerance: Optional[float] = None) -> Dict[str, Any]:
        """Metrics of the predicted residual over the part, labelled as a prediction."""
        m = metrics(self.residual, part_mask(self.target), tolerance)
        return {"quantity": "predicted residual (surrogate)",
                "model_data_source": self.data_source, **m}

    def save(self, out_dir: PathLike) -> Path:
        """compensated.npz, predicted.npz, residual.npz (+ bounds), compensation.json."""
        d = Path(out_dir)
        d.mkdir(parents=True, exist_ok=True)
        self.compensated.save(d / "compensated.npz")
        self.predicted.save(d / "predicted.npz")
        self.residual.save(d / "residual.npz")
        if self.residual_interval is not None:
            self.residual_interval[0].save(d / "residual_lower.npz")
            self.residual_interval[1].save(d / "residual_upper.npz")
        write_json(d / "compensation.json", to_jsonable({
            "method": "surrogate displacement adjustment", "stopped": self.stopped,
            "best_iteration": self.best_iteration, "level": self.level,
            "model_data_source": self.data_source, "history": self.history,
            "predicted_residual": self.predicted_metrics(), "ood": self.ood,
            "ood_target": self.ood_target, "verification": self.verification}))
        return d


def surrogate_compensate(target: HeightMap, setup: Any, surrogate: Any, *,
                         iterations: int = 8, alpha: float = 1.0,
                         tolerance: Optional[float] = None, stagnation: float = 0.02,
                         smoothing: Optional[float] = None, direction: str = "vertical",
                         level: float = 0.9,
                         max_wall_angle_deg: Optional[float] = MAX_WALL_ANGLE_DEG
                         ) -> SurrogateCompensation:
    """Displacement adjustment of `target` on `surrogate` (see the module docstring).

    iterations : at most this many surrogate predictions.
    tolerance : stop when the predicted RMS error over the part is <= this [m].
    stagnation : stop when an iteration improves the predicted RMS error by
        less than this fraction (0.02 = 2 %), or makes it worse.
    alpha, smoothing [m], direction, max_wall_angle_deg : as for
        `precomp.compensation.displacement_adjustment`.
    level : coverage of the residual interval (needs a calibrated surrogate;
        without one the interval is None and the result says so).
    """
    if iterations < 1:
        raise ValueError("iterations must be >= 1")
    if not 0 <= stagnation < 1:
        raise ValueError("stagnation must lie in [0, 1)")
    setup = as_setup(setup)
    pred = SurrogatePredictor(surrogate, setup)
    c = target.copy()
    history: List[Dict[str, Any]] = []
    best: Optional[Tuple[float, int, HeightMap, HeightMap, Optional[np.ndarray]]] = None
    stopped = "iterations"
    prev = None
    for k in range(iterations):
        da = displacement_adjustment(target, pred, iterations=1, alpha=alpha,
                                     direction=direction, smoothing=smoothing,
                                     max_wall_angle_deg=max_wall_angle_deg, initial=c)
        h = dict(da.history[0])
        h["iteration"] = k
        h["predicted"] = True
        rms = h["error"]["rms"]
        history.append(h)
        if best is None or rms < best[0]:
            best = (rms, k, da.commanded, da.formed, pred.last_std)
        if tolerance is not None and rms <= tolerance:
            stopped = "tolerance"
            break
        if prev is not None and rms > (1.0 - stagnation) * prev:
            stopped = "stagnation"
            break
        prev = rms
        c = da.proposed
    assert best is not None
    _, kbest, comp, formed, std = best
    comp = comp.copy()
    comp.metadata["compensation"] = {"method": "surrogate displacement adjustment",
                                     "iteration": kbest, "alpha": alpha,
                                     "model_data_source": getattr(surrogate, "data_source",
                                                                  None)}
    residual = vertical_deviation(formed, target)
    residual.metadata.update(quantity="predicted residual (surrogate)",
                             model_data_source=getattr(surrogate, "data_source", None))
    interval = None
    cal = getattr(surrogate, "calibrator", None)
    if hasattr(surrogate, "predict_interval") and cal is not None and cal.fitted:
        lo, hi = surrogate.predict_interval(comp, setup, level)
        interval = (residual.with_z(comp.z + lo - target.z),
                    residual.with_z(comp.z + hi - target.z))
    ood = dict(surrogate.assess(comp, setup)) if getattr(surrogate, "ood", None) else None
    ood_t = dict(surrogate.assess(target, setup)) if getattr(surrogate, "ood", None) else None
    return SurrogateCompensation(target, comp, formed, residual, interval, std, level, ood,
                                 ood_t, history, stopped, kbest,
                                 str(getattr(surrogate, "data_source", "unknown")))


def _deviation_report(formed: HeightMap, target: HeightMap,
                      tolerance: Optional[float]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    vert = vertical_deviation(formed, target)
    norm = signed_deviation(formed, target)
    for key, dev in (("vertical_deviation", vert), ("normal_deviation", norm)):
        out[key] = {}
        for name, m in region_masks(target).items():
            sel = m & dev.mask
            if sel.any():
                out[key][name] = metrics(dev, sel, tolerance)
    return out


def verify_with_simulator(result: SurrogateCompensation, setup: Any, simulator: Any, *,
                          uncompensated: bool = True, tolerance: Optional[float] = None
                          ) -> Dict[str, Any]:
    """Form the compensated surface with `simulator` (a `generate.Simulator`)
    and measure the achieved deviation from the target, per region, vertical
    and normal [m]; with `uncompensated` also form the target itself and
    report the improvement factor (uncompensated / compensated RMS vertical
    deviation over the part). The record states its data source and is stored
    in `result.verification`."""
    setup = as_setup(setup)
    jobs = [(setup, result.compensated)]
    if uncompensated:
        jobs.append((setup, result.target))
    outcomes = simulator.run(jobs)
    for oc in outcomes:
        if not oc.ok:
            raise PrecompError(f"verification run failed: {oc.error}")
    source = outcomes[0].provenance.get("source", getattr(simulator, "source", "sim"))
    rec: Dict[str, Any] = {"data_source": SOURCE_LABELS.get(source, source),
                           "simulator": outcomes[0].provenance.get("fidelity"),
                           "compensated": _deviation_report(outcomes[0].formed, result.target,
                                                            tolerance)}
    for key in ("sparlab_version", "deck_hash", "runtime_s", "cache_hit"):
        if key in outcomes[0].provenance:
            rec[key] = outcomes[0].provenance[key]
    if uncompensated:
        rec["uncompensated"] = _deviation_report(outcomes[1].formed, result.target, tolerance)
        a = rec["uncompensated"]["vertical_deviation"]["part"]["rms"]
        b = rec["compensated"]["vertical_deviation"]["part"]["rms"]
        rec["improvement_factor"] = a / b if b > 0 else float("inf")
        if "deck_hash" in outcomes[1].provenance:
            rec["uncompensated_deck_hash"] = outcomes[1].provenance["deck_hash"]
    pred = result.predicted_metrics(tolerance)["rms"]
    rec["predicted_part_rms_m"] = pred
    rec["achieved_part_rms_m"] = rec["compensated"]["vertical_deviation"]["part"]["rms"]
    result.verification = rec
    return rec


def verify_with_fea(result: SurrogateCompensation, setup: Any, work_dir: PathLike, *,
                    uncompensated: bool = True, tolerance: Optional[float] = None,
                    max_workers: int = 2) -> Dict[str, Any]:
    """`verify_with_simulator` with sparlab_form through the run cache in
    `work_dir` (data source "SparLab simulation", with the solver version
    and deck hash)."""
    from .generate import SparlabSimulator

    return verify_with_simulator(result, setup,
                                 SparlabSimulator(work_dir, max_workers=max_workers,
                                                  executor="serial" if max_workers == 1
                                                  else "process"),
                                 uncompensated=uncompensated, tolerance=tolerance)


__all__ = ["SurrogatePredictor", "SurrogateCompensation", "surrogate_compensate",
           "verify_with_fea", "verify_with_simulator"]
