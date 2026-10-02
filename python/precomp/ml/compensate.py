"""Compensation on a learned surrogate, and its honest end-to-end check.

`surrogate_compensate` runs displacement adjustment (`precomp.compensation`)
with the surrogate as the predictor of the formed surface,

    formed(c) = c + dz_hat(c),

one iteration at a time so that it can stop when the PREDICTED deviation from
the target is within tolerance or stops improving (stagnation) - or, as an
option, within the model's own interval half-width. It returns the
compensated commanded surface (the best predicted iterate), the predicted
formed surface, the predicted residual deviation from the target with its
conformal interval, the out-of-distribution report of the compensated shape,
and the history - all predictions, labelled with the data source the model
learned from.

The envelope is enforced, not only reported: a target outside the model's
training envelope is refused, and every iterate is assessed before the model
is evaluated on it - the first one outside ends the loop ("left_envelope"),
keeping the best iterate inside - unless `allow_out_of_envelope`, in which
case the verdicts are recorded in the history and the result.

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

from .._util import PathLike, PrecompError, call_with_target, to_jsonable, write_json
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

    def __init__(self, model: Any, setup: Any, target: Optional[HeightMap] = None):
        super().__init__(model, setup, target)
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
    the part, the update size, the envelope verdict and the interval
    half-width. stopped : "tolerance", "within_interval", "stagnation",
    "left_envelope" or "iterations". data_source : what the model learned
    from.
    verification : filled by `verify_with_fea` / `verify_with_simulator`.
    allowed_out_of_envelope : whether the envelope was overridden.
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
    allowed_out_of_envelope: bool = False

    @property
    def in_envelope(self) -> Optional[bool]:
        """Whether the target and the compensated shape are both inside the
        envelope (None without one)."""
        if self.ood is None or self.ood_target is None:
            return None
        return bool(self.ood["in_envelope"] and self.ood_target["in_envelope"])

    def predicted_metrics(self, tolerance: Optional[float] = None) -> Dict[str, Any]:
        """Metrics of the predicted residual over the part, labelled as a
        prediction, with the mean half-width of its interval (the model's
        own uncertainty, which a predicted residual near 0 does not remove)."""
        pm = part_mask(self.target)
        m = metrics(self.residual, pm, tolerance)
        out = {"quantity": "predicted residual (surrogate)",
               "model_data_source": self.data_source, **m}
        if self.residual_interval is not None:
            lo, hi = self.residual_interval
            out["interval_level"] = self.level
            out["interval_halfwidth_mean_m"] = float(np.mean(0.5 * (hi.z - lo.z)[pm]))
        return out

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
            "model_data_source": self.data_source,
            "history_quantity": "predictions of the surrogate (model data source: "
                                f"{self.data_source}), not simulations",
            "history": self.history,
            "predicted_residual": self.predicted_metrics(), "in_envelope": self.in_envelope,
            "allowed_out_of_envelope": self.allowed_out_of_envelope, "ood": self.ood,
            "ood_target": self.ood_target, "verification": self.verification}))
        return d


def surrogate_compensate(target: HeightMap, setup: Any, surrogate: Any, *,
                         iterations: int = 8, alpha: float = 1.0,
                         tolerance: Optional[float] = None, stagnation: float = 0.02,
                         smoothing: Optional[float] = None, direction: str = "vertical",
                         level: float = 0.9,
                         max_wall_angle_deg: Optional[float] = MAX_WALL_ANGLE_DEG,
                         allow_out_of_envelope: bool = False,
                         interval_stop: bool = False,
                         upper_bound: Any = None,
                         hold_mask: Optional[np.ndarray] = None,
                         adjust_mask: Optional[np.ndarray] = None,
                         tool_radius: Optional[float] = None) -> SurrogateCompensation:
    """Displacement adjustment of `target` on `surrogate` (see the module docstring).

    iterations : at most this many surrogate predictions.
    tolerance : stop when the predicted RMS error over the part is <= this [m].
    stagnation : stop when an iteration improves the predicted RMS error by
        less than this fraction (0.02 = 2 %), or makes it worse.
    alpha, smoothing [m], direction, max_wall_angle_deg : as for
        `precomp.compensation.displacement_adjustment`.
    level : coverage of the residual interval (needs a calibrated surrogate;
        without one the interval is None and the result says so).
    allow_out_of_envelope : compensate a target outside the training
        envelope, and keep iterating when an iterate leaves it (PrecompError
        for such a target, and a stop at such an iterate, otherwise). Without
        an envelope there is nothing to enforce.
    interval_stop : with a calibrated surrogate, stop ("within_interval")
        once the predicted RMS error over the part is within the RMS
        half-width of the model's interval at `level` at that iterate - the
        model cannot resolve less. It saves about two thirds of the
        predictions (2 instead of 6), but whether it helps is unclear: on
        proxy data the verified deviation was 4 % smaller on the 23 held-out
        targets of the tests (0.054 against 0.056 mm), 9 % larger on the 11
        of the demonstration (0.049 against 0.045 mm) and the same on the 10
        of its freeform run (0.040 mm). Off by default.
    upper_bound, hold_mask, adjust_mask, tool_radius : as for
        `displacement_adjustment` (a support's command bound and masks,
        `precomp.fea.support`; the tool's reach). None: the defaults.
    """
    if iterations < 1:
        raise ValueError("iterations must be >= 1")
    if not 0 <= stagnation < 1:
        raise ValueError("stagnation must lie in [0, 1)")
    setup = as_setup(setup)
    has_env = getattr(surrogate, "ood", None) is not None and hasattr(surrogate, "assess")
    ood_t = dict(surrogate.assess(target, setup)) if has_env else None
    if ood_t is not None and not ood_t["in_envelope"] and not allow_out_of_envelope:
        raise PrecompError(
            "the target is outside the model's training envelope (reasons: "
            f"{', '.join(ood_t['reasons'][:5]) or 'fraction of nodes out'}; part score "
            f"{ood_t['part_score']:.3g}): the surrogate would extrapolate. Train on such parts, "
            "or pass allow_out_of_envelope=True (--allow-out-of-envelope) to compensate "
            "anyway")
    pred = SurrogatePredictor(surrogate, setup, target)
    cal = getattr(surrogate, "calibrator", None)
    can_interval = (hasattr(surrogate, "predict_interval") and cal is not None
                    and cal.fitted and cal.supports(level))
    pm = part_mask(target)
    c = target.copy()
    history: List[Dict[str, Any]] = []
    best: Optional[Tuple[float, int, HeightMap, HeightMap, Optional[np.ndarray]]] = None
    stopped = "iterations"
    prev = None
    verdicts: Dict[int, Dict[str, Any]] = {}
    for k in range(iterations):
        a = ood_t if k == 0 else (dict(surrogate.assess(c, setup)) if has_env else None)
        if a is not None:
            verdicts[k] = a
            if not a["in_envelope"] and k > 0 and not allow_out_of_envelope:
                history.append({"iteration": k, "predicted": False, "in_envelope": False,
                                "ood_part_score": a["part_score"],
                                "ood_reasons": a["reasons"][:5],
                                "note": "iterate outside the training envelope: not evaluated"})
                stopped = "left_envelope"
                break
        da = displacement_adjustment(target, pred, iterations=1, alpha=alpha,
                                     direction=direction, smoothing=smoothing,
                                     max_wall_angle_deg=max_wall_angle_deg, initial=c,
                                     upper_bound=upper_bound, hold_mask=hold_mask,
                                     adjust_mask=adjust_mask, tool_radius=tool_radius)
        h = dict(da.history[0])
        h["iteration"] = k
        h["predicted"] = True
        if a is not None:
            h["in_envelope"] = bool(a["in_envelope"])
            h["ood_part_score"] = a["part_score"]
        rms = h["error"]["rms"]
        history.append(h)
        if best is None or rms < best[0]:
            best = (rms, k, da.commanded, da.formed, pred.last_std)
        if tolerance is not None and rms <= tolerance:
            stopped = "tolerance"
            break
        if interval_stop and can_interval:
            lo, hi = call_with_target(surrogate.predict_interval, da.commanded, setup, level,
                                      target=target)
            hw = float(np.sqrt(np.mean((0.5 * (hi - lo))[pm] ** 2)))
            h["interval_halfwidth_rms_m"] = hw
            if rms <= hw:
                stopped = "within_interval"
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
    if can_interval:
        lo, hi = call_with_target(surrogate.predict_interval, comp, setup, level,
                                  target=target)
        interval = (residual.with_z(comp.z + lo - target.z),
                    residual.with_z(comp.z + hi - target.z))
    ood = verdicts.get(kbest) if has_env else None
    return SurrogateCompensation(target, comp, formed, residual, interval, std, level, ood,
                                 ood_t, history, stopped, kbest,
                                 str(getattr(surrogate, "data_source", "unknown")),
                                 allowed_out_of_envelope=bool(allow_out_of_envelope))


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
    from .generate import sim_job

    setup = as_setup(setup)
    jobs = [sim_job(setup, result.compensated, result.target)]
    if uncompensated:
        jobs.append(sim_job(setup, result.target, result.target))
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
