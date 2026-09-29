"""Calibrated prediction intervals by grouped split-conformal prediction.

A model gives a mean mu and a standard deviation s at every node; s is an
ensemble spread or a predicted variance and is not, by itself, calibrated.
`ConformalCalibrator` is fitted on held-out WHOLE parts (never on parts the
model was trained on): with the normalised non-conformity score

    r = |y - mu| / (s + eps),

eps a small stabiliser (by default a tenth of the median s on the calibration
set, recorded), the interval at level 1 - alpha is mu +- q (s + eps), q the
ceil((n + 1)(1 - alpha))-th smallest of the n calibration scores (the
finite-sample correction of split conformal prediction). Scores of the same
part are correlated, so the guarantee - marginal coverage >= 1 - alpha for a
new node of a new part exchangeable with the calibration parts - holds
approximately, with the number of parts rather than of nodes as the effective
sample size; `coverage_report` measures it on held-out parts, overall, per
region and per family. With `mondrian=True` the quantile is taken per region
(rim / wall / base / flange), which equalises coverage across regions.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .._util import PrecompError
from .features import REGION_NAMES

#: Levels reported in manifests and calibration tables.
REPORT_LEVELS = (0.5, 0.8, 0.9, 0.95, 0.99)


class ConformalCalibrator:
    """Grouped split-conformal calibration of (mu, std) (see the module docstring).

    eps : the stabiliser [m]; None sets it to `eps_fraction` x the median std
        of the calibration set (or of |y - mu| when every std is zero).
    mondrian : one quantile per region code (needs `region` at fit and use).
    max_scores : calibration scores kept (a seeded random subset beyond it);
        the finite-sample correction uses the full count.
    """

    def __init__(self, *, eps: Optional[float] = None, eps_fraction: float = 0.1,
                 mondrian: bool = False, max_scores: int = 200_000, seed: int = 0):
        if eps is not None and not eps > 0:
            raise ValueError("eps must be > 0")
        self.eps = eps
        self.eps_fraction = float(eps_fraction)
        self.mondrian = bool(mondrian)
        self.max_scores = int(max_scores)
        self.seed = int(seed)
        self.scores: Dict[str, np.ndarray] = {}
        self.counts: Dict[str, int] = {}
        self.n_groups = 0
        self.data_source: Optional[str] = None

    @property
    def fitted(self) -> bool:
        return bool(self.scores)

    def fit(self, mu: np.ndarray, std: np.ndarray, y: np.ndarray, groups: np.ndarray, *,
            region: Optional[np.ndarray] = None,
            data_source: Optional[str] = None) -> "ConformalCalibrator":
        mu, std, y = (np.asarray(a, dtype=float).ravel() for a in (mu, std, y))
        groups = np.asarray(groups).ravel()
        if not (mu.shape == std.shape == y.shape == groups.shape):
            raise ValueError("mu, std, y and groups must have one entry per node")
        if mu.size == 0:
            raise ValueError("no calibration nodes")
        if not (np.all(np.isfinite(mu)) and np.all(np.isfinite(std)) and np.all(np.isfinite(y))):
            raise ValueError("calibration data contain NaN or infinity")
        if np.any(std < 0):
            raise ValueError("std must be >= 0")
        self.n_groups = int(np.unique(groups).size)
        if self.n_groups < 2:
            raise ValueError("conformal calibration needs at least two held-out parts")
        if self.eps is None:
            med = float(np.median(std))
            base = med if med > 0 else float(np.median(np.abs(y - mu)))
            self.eps = self.eps_fraction * base if base > 0 else 1e-12
        s = np.abs(y - mu) / (std + self.eps)
        rng = np.random.default_rng(self.seed)
        keys: List[Tuple[str, np.ndarray]] = [("all", np.ones(s.size, bool))]
        if self.mondrian:
            if region is None:
                raise ValueError("mondrian calibration needs the region of every node")
            region = np.asarray(region).ravel()
            keys += [(REGION_NAMES[r], region == r) for r in np.unique(region)]
        self.scores, self.counts = {}, {}
        for name, sel in keys:
            v = s[sel]
            self.counts[name] = int(v.size)
            if v.size > self.max_scores:
                v = rng.choice(v, size=self.max_scores, replace=False)
            self.scores[name] = np.sort(v)
        self.data_source = data_source
        return self

    def quantile(self, level: float, region: Optional[str] = None) -> float:
        """The score quantile q at coverage `level` (finite-sample corrected).

        Raises PrecompError when the calibration set is too small for the
        level (ceil((n + 1) level) > n): the interval would be unbounded.
        """
        if not self.fitted:
            raise PrecompError("the calibrator is not fitted")
        if not 0 < level < 1:
            raise ValueError("level must lie in (0, 1)")
        key = region if (self.mondrian and region is not None) else "all"
        if key not in self.scores:
            key = "all"
        v = self.scores[key]
        n = self.counts[key]
        k = math.ceil((n + 1) * level)
        if k > n:
            raise PrecompError(f"{n} calibration nodes cannot support level {level}")
        # position of the k-th of n in the (possibly subsampled) sorted scores
        pos = (k - 1) / max(n - 1, 1) * (v.size - 1)
        return float(np.interp(pos, np.arange(v.size), v))

    def intervals(self, mu: np.ndarray, std: np.ndarray, level: float,
                  region: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray]:
        """(lower, upper) of dz [m] at coverage `level`."""
        mu = np.asarray(mu, dtype=float)
        std = np.asarray(std, dtype=float)
        if self.mondrian and region is not None:
            half = np.empty_like(mu)
            region = np.asarray(region)
            for r in np.unique(region):
                sel = region == r
                half[sel] = self.quantile(level, REGION_NAMES[int(r)]) * (std[sel] + self.eps)
        else:
            half = self.quantile(level) * (std + self.eps)
        return mu - half, mu + half

    def summary(self, levels: Sequence[float] = REPORT_LEVELS) -> Dict[str, Any]:
        out: Dict[str, Any] = {"eps_m": self.eps, "n_nodes": self.counts.get("all", 0),
                               "n_parts": self.n_groups, "mondrian": self.mondrian,
                               "data_source": self.data_source, "quantiles": {}}
        for lv in levels:
            try:
                out["quantiles"][f"{lv:g}"] = self.quantile(lv)
            except PrecompError:
                out["quantiles"][f"{lv:g}"] = None
        return out

    def to_state(self) -> Dict[str, Any]:
        return {"eps": self.eps, "eps_fraction": self.eps_fraction, "mondrian": self.mondrian,
                "max_scores": self.max_scores, "seed": self.seed, "scores": dict(self.scores),
                "counts": dict(self.counts), "n_groups": self.n_groups,
                "data_source": self.data_source}

    @classmethod
    def from_state(cls, state: Mapping[str, Any]) -> "ConformalCalibrator":
        c = cls(eps=state["eps"], eps_fraction=state["eps_fraction"],
                mondrian=state["mondrian"], max_scores=state["max_scores"], seed=state["seed"])
        c.scores = {k: np.asarray(v) for k, v in state["scores"].items()}
        c.counts = {k: int(v) for k, v in state["counts"].items()}
        c.n_groups = int(state["n_groups"])
        c.data_source = state.get("data_source")
        return c


def intervals(mu: np.ndarray, std: np.ndarray, level: float, calibrator: ConformalCalibrator,
              region: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray]:
    """(lower, upper) of dz [m] at `level` from a fitted calibrator."""
    return calibrator.intervals(mu, std, level, region)


def partition_coverage(mu: np.ndarray, std: np.ndarray, y: np.ndarray, groups: np.ndarray, *,
                       level: float, n_partitions: int = 50, calibration_fraction: float = 0.5,
                       seed: int = 0, **calibrator_kwargs: Any) -> Dict[str, float]:
    """Coverage of repeated grouped split-conformal calibration.

    The held-out parts (`groups`, e.g. part ids) are split at random
    `n_partitions` times into calibration and test parts; each time a
    calibrator is fitted on the calibration parts and its coverage measured
    on the test parts. The mean estimates the expected coverage - what the
    conformal guarantee is about; the spread shows how far one calibration
    draw with this many parts can land from it (with a few dozen parts and
    part-correlated errors, several per cent).
    """
    mu, std, y = (np.asarray(a, dtype=float).ravel() for a in (mu, std, y))
    groups = np.asarray(groups).ravel()
    uniq = np.unique(groups)
    if uniq.size < 4:
        raise ValueError("partition_coverage needs at least four parts")
    rng = np.random.default_rng(seed)
    n_cal = int(min(max(2, round(calibration_fraction * uniq.size)), uniq.size - 2))
    covs = []
    for _ in range(int(n_partitions)):
        cal = np.isin(groups, rng.permutation(uniq)[:n_cal])
        c = ConformalCalibrator(**calibrator_kwargs).fit(mu[cal], std[cal], y[cal],
                                                         groups[cal])
        lo, hi = c.intervals(mu[~cal], std[~cal], level)
        covs.append(float(np.mean((y[~cal] >= lo) & (y[~cal] <= hi))))
    c = np.array(covs)
    return {"level": float(level), "mean": float(c.mean()), "std": float(c.std()),
            "min": float(c.min()), "max": float(c.max()), "n_partitions": int(n_partitions),
            "n_parts": int(uniq.size), "n_calibration_parts": n_cal}


def coverage_report(y: np.ndarray, lower: np.ndarray, upper: np.ndarray, *, level: float,
                    data_source: str, groups: Optional[np.ndarray] = None,
                    region: Optional[np.ndarray] = None,
                    family: Optional[np.ndarray] = None) -> pd.DataFrame:
    """Empirical coverage of intervals: overall, per region and per family.

    Columns: scope ("all", "region", "family"), name, level, n_nodes, n_parts,
    coverage (share of nodes inside), part_coverage_min / _median (over parts,
    when `groups` is given), mean_width_m, data_source.
    """
    y, lower, upper = (np.asarray(a, dtype=float).ravel() for a in (y, lower, upper))
    inside = (y >= lower) & (y <= upper)
    width = upper - lower
    rows = []

    def row(scope: str, name: str, sel: np.ndarray) -> None:
        if not sel.any():
            return
        r = {"scope": scope, "name": name, "level": level, "n_nodes": int(sel.sum()),
             "coverage": float(inside[sel].mean()), "mean_width_m": float(width[sel].mean()),
             "data_source": data_source}
        if groups is not None:
            g = np.asarray(groups).ravel()[sel]
            per = pd.Series(inside[sel]).groupby(g).mean()
            r.update(n_parts=int(per.size), part_coverage_min=float(per.min()),
                     part_coverage_median=float(per.median()))
        rows.append(r)

    row("all", "all", np.ones(y.size, bool))
    if region is not None:
        reg = np.asarray(region).ravel()
        for code in np.unique(reg):
            row("region", REGION_NAMES[int(code)], reg == code)
    if family is not None:
        fam = np.asarray(family).ravel()
        for f in np.unique(fam):
            row("family", str(f), fam == f)
    return pd.DataFrame(rows)


__all__ = ["ConformalCalibrator", "intervals", "coverage_report", "partition_coverage",
           "REPORT_LEVELS"]
