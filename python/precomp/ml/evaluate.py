"""Evaluation on held-out parts: errors, calibration and the envelope, as tidy frames.

Protocol. A model is only ever scored on WHOLE parts it has not seen: grouped
K-fold cross-validation (`cross_validate`), a held-out split, or a held-out
family (`family_holdout`), where every sample of one family is removed from
training and calibration. Within each training portion a further grouped split
holds out the calibration parts of the conformal intervals. The error is that
of dz on the part nodes of the commanded surface (depth > part_eps):

* per part - rms, mae, max_abs, p95_abs, bias and the share within a
  tolerance of the dz error [m] (`precomp.metrology.metrics`), the rms of dz
  itself and their ratio, the interval coverage and width at the report level,
  the mean predicted std, and the envelope verdict;
* per region - the same pooled over every node of a region (rim, wall, base);
* calibration - empirical coverage against the nominal level, pooled and per
  part;
* envelope - the share of parts flagged out of the envelope per split: a
  held-out family must be flagged more often than in-distribution parts.

Every frame carries `data_source`, the label of the evaluation samples
("SparLab simulation", "proxy - not physics", "scan").
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from .._util import PathLike, PrecompError
from ..metrology import metrics
from .dataset import (Sample, family_split, grouped_kfold, grouped_split, index_frame,
                      source_label)
from .features import DEFAULT_CONFIG, REGION_NAMES, FeatureConfig
from .surrogate import DeviationSurrogate, train_surrogate
from .uncertainty import REPORT_LEVELS

#: Default tolerance [m] of the "share within tolerance" metric (0.2 mm).
DEFAULT_TOLERANCE = 2e-4


@dataclass
class EvaluationReport:
    """Tidy frames of one evaluation (see the module docstring)."""

    per_part: pd.DataFrame
    per_region: pd.DataFrame
    calibration: pd.DataFrame
    ood: pd.DataFrame
    data_source: str
    info: Dict[str, Any] = field(default_factory=dict)

    def summary(self) -> Dict[str, Any]:
        """Headline numbers per split: mean and median per-part rms error [m],
        relative error, coverage at the report level, share flagged."""
        out: Dict[str, Any] = {"data_source": self.data_source}
        for split, g in self.per_part.groupby("split"):
            level = g["level"].iloc[0]
            out[split] = {"n_parts": int(g["part_id"].nunique()), "n_samples": int(len(g)),
                          "err_rms_mean_m": float(g["err_rms"].mean()),
                          "err_rms_median_m": float(g["err_rms"].median()),
                          "err_max_abs_mean_m": float(g["err_max_abs"].mean()),
                          "rel_rms_mean": float(g["rel_rms"].mean()),
                          "dz_rms_mean_m": float(g["dz_rms"].mean()),
                          "coverage_level": float(level),
                          "coverage_mean": float(g["coverage"].mean()),
                          "ood_flagged_share": float((~g["ood_in_envelope"].astype(bool)).mean())
                          if g["ood_in_envelope"].notna().all() else None}
        return out

    def to_csv(self, out_dir: PathLike, prefix: str = "") -> Dict[str, Path]:
        d = Path(out_dir)
        d.mkdir(parents=True, exist_ok=True)
        paths = {}
        for name in ("per_part", "per_region", "calibration", "ood"):
            p = d / f"{prefix}{name}.csv"
            getattr(self, name).to_csv(p, index=False, float_format="%.6g")
            paths[name] = p
        return paths

    @staticmethod
    def concatenate(reports: Sequence["EvaluationReport"]) -> "EvaluationReport":
        if not reports:
            raise ValueError("no reports")
        cat = lambda a: pd.concat([getattr(r, a) for r in reports], ignore_index=True)  # noqa
        labels = sorted({r.data_source for r in reports})
        return EvaluationReport(cat("per_part"), cat("per_region"), cat("calibration"),
                                cat("ood"), labels[0] if len(labels) == 1 else
                                "mixed: " + " + ".join(labels),
                                {"parts": [r.info for r in reports]})


def evaluate_surrogate(surrogate: DeviationSurrogate, samples: Sequence[Sample], *,
                       tolerance: float = DEFAULT_TOLERANCE, level: float = 0.9,
                       levels: Sequence[float] = REPORT_LEVELS, split: str = "held_out",
                       assess: bool = True, check_unseen: bool = True) -> EvaluationReport:
    """Score `surrogate` on `samples` (see the module docstring).

    check_unseen : refuse samples whose part the surrogate was trained or
        calibrated on (by sample id and part id recorded in its training).
    """
    samples = list(samples)
    if not samples:
        raise ValueError("no samples to evaluate")
    if check_unseen:
        seen = set(surrogate.training.get("sample_ids", [])) | set(
            surrogate.training.get("calibration_sample_ids", []))
        leaked = sorted(s.sample_id for s in samples if s.sample_id in seen)
        if leaked:
            raise PrecompError(f"{len(leaked)} evaluation sample(s) were used to train or "
                               f"calibrate the model, e.g. {leaked[:3]}")
    label = source_label(s.source for s in samples)
    has_cal = surrogate.calibrator is not None and surrogate.calibrator.fitted
    rows, pooled = [], []
    for s in samples:
        try:
            mu, sd = surrogate.predict_deviation(s.commanded, s.setup)
            fm = surrogate.feature_maps(s.commanded, s.setup)
        except (ValueError, PrecompError) as exc:
            raise PrecompError(f"sample {s.sample_id}: {exc}") from exc
        sel = fm.part & s.valid
        if not sel.any():
            raise PrecompError(f"sample {s.sample_id} has no valid part node")
        y = s.dz[sel]
        m = mu[sel]
        err = m - y
        em = metrics(err, tolerance=tolerance)
        dz_rms = float(np.sqrt(np.mean(y * y)))
        row: Dict[str, Any] = {
            "split": split, "sample_id": s.sample_id, "part_id": s.part_id,
            "family": s.family, "kind": s.kind, "source": s.source, "data_source": label,
            "n_nodes": int(sel.sum()), **{f"err_{k}": v for k, v in em.items()
                                          if k not in ("n", "tolerance")},
            "tolerance_m": tolerance, "dz_rms": dz_rms,
            "rel_rms": em["rms"] / dz_rms if dz_rms > 0 else np.nan,
            "std_mean_m": float(sd[sel].mean()), "level": level}
        block = {"sample_id": s.sample_id, "family": s.family, "y": y, "mu": m,
                 "sd": sd[sel], "region": fm.region[sel]}
        if has_cal:
            for lv in levels:
                lo, hi = surrogate.calibrator.intervals(mu, sd, lv, fm.region)
                inside = (y >= lo[sel]) & (y <= hi[sel])
                block[f"in_{lv:g}"] = inside
                block[f"width_{lv:g}"] = (hi - lo)[sel]
            row["coverage"] = float(block[f"in_{level:g}"].mean()) if f"in_{level:g}" in block \
                else np.nan
            row["width_mean_m"] = float(block[f"width_{level:g}"].mean()) \
                if f"width_{level:g}" in block else np.nan
        else:
            row["coverage"] = np.nan
            row["width_mean_m"] = np.nan
        if assess and surrogate.ood is not None:
            a = surrogate.assess(s.commanded, s.setup)
            row.update(ood_in_envelope=bool(a["in_envelope"]), ood_part_score=a["part_score"],
                       ood_fraction_points_out=a["fraction_points_out"],
                       ood_reasons=";".join(a["reasons"]))
        else:
            row.update(ood_in_envelope=None, ood_part_score=np.nan,
                       ood_fraction_points_out=np.nan, ood_reasons="")
        rows.append(row)
        pooled.append(block)
    per_part = pd.DataFrame(rows)
    # pooled per region
    reg_rows = []
    y_all = np.concatenate([b["y"] for b in pooled])
    e_all = np.concatenate([b["mu"] for b in pooled]) - y_all
    r_all = np.concatenate([b["region"] for b in pooled])
    for name, sel in [("all", np.ones(y_all.size, bool))] + [
            (REGION_NAMES[c], r_all == c) for c in np.unique(r_all)]:
        em = metrics(e_all[sel], tolerance=tolerance)
        r = {"split": split, "region": name, "n_nodes": int(sel.sum()), "data_source": label,
             **{f"err_{k}": v for k, v in em.items() if k not in ("n", "tolerance")},
             "tolerance_m": tolerance,
             "dz_rms": float(np.sqrt(np.mean(y_all[sel] ** 2)))}
        if has_cal:
            ins = np.concatenate([b[f"in_{level:g}"] for b in pooled])
            r["coverage"] = float(ins[sel].mean())
            r["level"] = level
        reg_rows.append(r)
    per_region = pd.DataFrame(reg_rows)
    # calibration table
    cal_rows = []
    if has_cal:
        for lv in levels:
            ins = np.concatenate([b[f"in_{lv:g}"] for b in pooled])
            per = np.array([b[f"in_{lv:g}"].mean() for b in pooled])
            width = np.concatenate([b[f"width_{lv:g}"] for b in pooled])
            cal_rows.append({"split": split, "level": lv, "coverage": float(ins.mean()),
                             "coverage_minus_level": float(ins.mean() - lv),
                             "part_coverage_min": float(per.min()),
                             "part_coverage_median": float(np.median(per)),
                             "mean_width_m": float(width.mean()), "n_nodes": int(ins.size),
                             "n_parts": len(pooled),
                             "calibrated_on": surrogate.calibrator.data_source,
                             "data_source": label})
    calibration = pd.DataFrame(cal_rows)
    # envelope
    ood_rows = []
    if assess and surrogate.ood is not None:
        for fam, g in [("all", per_part)] + list(per_part.groupby("family")):
            ood_rows.append({"split": split, "family": fam, "n_samples": int(len(g)),
                             "flagged_share": float((~g["ood_in_envelope"].astype(bool)).mean()),
                             "part_score_median": float(g["ood_part_score"].median()),
                             "fraction_points_out_median":
                                 float(g["ood_fraction_points_out"].median()),
                             "data_source": label})
    ood = pd.DataFrame(ood_rows)
    info = {"model_class": surrogate.model_class, "model_data_source": surrogate.data_source,
            "n_samples": len(samples), "split": split}
    return EvaluationReport(per_part, per_region, calibration, ood, label, info)


ModelFactory = Callable[[], Any]


def _fit(make_model: ModelFactory, train: List[Sample], cal: List[Sample],
         features: FeatureConfig, points_per_sample: Optional[int], seed: int,
         n_jobs: int) -> DeviationSurrogate:
    return train_surrogate(make_model(), train, calibration=cal, features=features,
                           points_per_sample=points_per_sample, seed=seed, n_jobs=n_jobs)


def _cal_split(train: List[Sample], fraction: float, seed: int
               ) -> tuple:
    byid = {s.sample_id: s for s in train}
    tr, cal = grouped_split(index_frame(train), test_fraction=fraction, seed=seed)
    return [byid[i] for i in tr], [byid[i] for i in cal]


def cross_validate(samples: Sequence[Sample], make_model: ModelFactory, *, n_splits: int = 5,
                   seed: int = 0, features: FeatureConfig = DEFAULT_CONFIG,
                   points_per_sample: Optional[int] = 2000, calibration_fraction: float = 0.25,
                   tolerance: float = DEFAULT_TOLERANCE, level: float = 0.9,
                   n_jobs: int = 1) -> EvaluationReport:
    """Grouped K-fold cross-validation by part id: each fold is scored by a
    model trained (and calibrated, on a grouped `calibration_fraction` of the
    training parts) without it. The frames gain a `fold` column."""
    samples = list(samples)
    byid = {s.sample_id: s for s in samples}
    reports = []
    for k, (tr_ids, te_ids) in enumerate(grouped_kfold(samples, n_splits=n_splits,
                                                       seed=seed)):
        train, cal = _cal_split([byid[i] for i in tr_ids], calibration_fraction, seed + k)
        sur = _fit(make_model, train, cal, features, points_per_sample, seed + k, n_jobs)
        rep = evaluate_surrogate(sur, [byid[i] for i in te_ids], tolerance=tolerance,
                                 level=level, split="cross_validation")
        for f in (rep.per_part, rep.per_region, rep.calibration, rep.ood):
            f.insert(0, "fold", k)
        reports.append(rep)
    return EvaluationReport.concatenate(reports)


def family_holdout(samples: Sequence[Sample], family: str, make_model: ModelFactory, *,
                   seed: int = 0, features: FeatureConfig = DEFAULT_CONFIG,
                   points_per_sample: Optional[int] = 2000, test_fraction: float = 0.25,
                   calibration_fraction: float = 0.3, tolerance: float = DEFAULT_TOLERANCE,
                   level: float = 0.9, n_jobs: int = 1,
                   model: Optional[DeviationSurrogate] = None) -> EvaluationReport:
    """Leave-family-out: train on the other families (minus a grouped
    `test_fraction` of their parts and a calibration split), then score the
    held-out in-distribution parts (split "in_distribution") and the held-out
    family (split "held_out_family"). Pass `model` to score an already
    trained surrogate (it must not have seen the family)."""
    samples = list(samples)
    byid = {s.sample_id: s for s in samples}
    rest, fam = family_split(samples, family)
    tr_ids, te_ids = grouped_split(index_frame([byid[i] for i in rest]),
                                   test_fraction=test_fraction, seed=seed)
    if model is None:
        train, cal = _cal_split([byid[i] for i in tr_ids], calibration_fraction, seed + 1)
        model = _fit(make_model, train, cal, features, points_per_sample, seed, n_jobs)
    elif family in model.training.get("families", []):
        raise PrecompError(f"the model was trained on family {family!r}")
    a = evaluate_surrogate(model, [byid[i] for i in te_ids], tolerance=tolerance, level=level,
                           split="in_distribution")
    b = evaluate_surrogate(model, [byid[i] for i in fam], tolerance=tolerance, level=level,
                           split="held_out_family")
    rep = EvaluationReport.concatenate([a, b])
    rep.info = {"family": family, "model_class": model.model_class,
                "model_data_source": model.data_source,
                "in_distribution": a.info, "held_out_family": b.info}
    return rep


__all__ = ["EvaluationReport", "evaluate_surrogate", "cross_validate", "family_holdout",
           "DEFAULT_TOLERANCE"]
