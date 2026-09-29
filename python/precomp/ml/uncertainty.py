"""Calibrated prediction intervals by split-conformal prediction over whole parts.

A model gives a mean mu and a standard deviation s at every node; s is an
ensemble spread or a predicted variance and is not, by itself, calibrated.
`ConformalCalibrator` is fitted on held-out WHOLE parts (never on parts the
model was trained on). With the normalised non-conformity score

    r = |y - mu| / (s + eps),

eps a small stabiliser (by default a tenth of the median s on the calibration
set, recorded), the interval at level L is mu +- q (s + eps) with

    q = the smallest t such that  sum_j F_j(t) >= L (k + 1),

F_j the empirical distribution of the scores of calibration part j (every
part weighs the same, whatever its number of nodes) and k the number of
calibration PARTS. The errors of the nodes of one part are strongly
correlated - most of a part's error is one smooth field - so the part, not
the node, is the exchangeable unit. The guarantee: for a new part
exchangeable with the calibration parts, the EXPECTED share of its nodes
inside the interval is at least L. (Proof: with the new part's F included
and weighted like the others, the threshold would be symmetric in all k + 1
parts, so the new part's expected coverage would equal the average, >= L;
leaving the new part out, i.e. counting its F as 0, can only raise q.) A
level needs k >= L / (1 - L) parts - 9 at 90 %, 19 at 95 %, 99 at 99 % - and
asking more is an error. The price of validity with few parts: when parts
behave alike the expected coverage approaches L (k + 1) / k (97.5 % at a
nominal 90 % with 12 parts). What it is NOT: a guarantee for any one part - a
part whose error is a large offset can be covered far less (the evaluation
reports the per-part coverage and the share of parts reaching the level).

With `mondrian=True` the quantile is taken per region (rim / wall / base),
each over the parts that have that region. The calibrator also keeps the
scores of a constant-width interval mu +- q0 (q0 [m] from |y - mu| the same
way): comparing widths at the same level shows whether the model's s carries
information (`baseline_intervals`, and the width ratio in the evaluation).
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


def min_parts(level: float) -> int:
    """Calibration parts needed for `level`: the smallest k with k >= level (k + 1)."""
    if not 0 < level < 1:
        raise ValueError("level must lie in (0, 1)")
    return int(math.ceil(level / (1.0 - level) - 1e-9))


def _part_sketch(values: np.ndarray, groups: np.ndarray, max_per_part: int
                 ) -> Tuple[np.ndarray, np.ndarray, int]:
    """(sorted values, weights, number of parts): every part's scores reduced
    to at most `max_per_part` evenly spaced order statistics, each part
    weighing 1 in total."""
    vals, weights = [], []
    uniq = np.unique(groups)
    for g in uniq:
        v = np.sort(values[groups == g])
        if v.size > max_per_part:
            v = v[np.linspace(0, v.size - 1, max_per_part).round().astype(int)]
        vals.append(v)
        weights.append(np.full(v.size, 1.0 / v.size))
    v = np.concatenate(vals)
    w = np.concatenate(weights)
    order = np.argsort(v, kind="stable")
    return v[order], w[order], int(uniq.size)


def _part_quantile(values: np.ndarray, weights: np.ndarray, k: int, level: float,
                   what: str) -> float:
    """The smallest t with sum_j F_j(t) >= level (k + 1) (see the module docstring)."""
    need = level * (k + 1)
    if need > k * (1 + 1e-12):
        raise PrecompError(f"{k} calibration part(s) cannot support level {level:g} "
                           f"({what}): it needs at least {min_parts(level)} held-out parts")
    cw = np.cumsum(weights)
    i = int(np.searchsorted(cw, need - 1e-9 * k))
    return float(values[min(i, values.size - 1)])


class ConformalCalibrator:
    """Split-conformal calibration of (mu, std) over whole parts (see the
    module docstring).

    eps : the stabiliser [m]; None sets it to `eps_fraction` x the median std
        of the calibration set (or of |y - mu| when every std is zero).
    mondrian : one quantile per region code (needs `region` at fit and use).
    max_per_part : scores kept per part (evenly spaced order statistics).
    """

    def __init__(self, *, eps: Optional[float] = None, eps_fraction: float = 0.1,
                 mondrian: bool = False, max_per_part: int = 2000):
        if eps is not None and not eps > 0:
            raise ValueError("eps must be > 0")
        if int(max_per_part) < 10:
            raise ValueError("max_per_part must be >= 10")
        self.eps = eps
        self.eps_fraction = float(eps_fraction)
        self.mondrian = bool(mondrian)
        self.max_per_part = int(max_per_part)
        #: per key ("all" or a region name): (sorted scores, weights, n parts)
        self.scores: Dict[str, Tuple[np.ndarray, np.ndarray, int]] = {}
        #: the same for |y - mu| [m], the constant-width baseline ("all" only)
        self.baseline: Optional[Tuple[np.ndarray, np.ndarray, int]] = None
        self.n_nodes = 0
        self.n_parts = 0
        self.data_source: Optional[str] = None

    @property
    def fitted(self) -> bool:
        return bool(self.scores)

    def fit(self, mu: np.ndarray, std: np.ndarray, y: np.ndarray, groups: np.ndarray, *,
            region: Optional[np.ndarray] = None,
            data_source: Optional[str] = None) -> "ConformalCalibrator":
        """`groups`: the PART id of every node (the exchangeable unit)."""
        mu, std, y = (np.asarray(a, dtype=float).ravel() for a in (mu, std, y))
        groups = np.asarray(groups).astype(str).ravel()
        if not (mu.shape == std.shape == y.shape == groups.shape):
            raise ValueError("mu, std, y and groups must have one entry per node")
        if mu.size == 0:
            raise ValueError("no calibration nodes")
        if not (np.all(np.isfinite(mu)) and np.all(np.isfinite(std)) and np.all(np.isfinite(y))):
            raise ValueError("calibration data contain NaN or infinity")
        if np.any(std < 0):
            raise ValueError("std must be >= 0")
        self.n_parts = int(np.unique(groups).size)
        if self.n_parts < 2:
            raise ValueError("conformal calibration needs at least two held-out parts "
                             f"(got {self.n_parts}; the variants of one part count once)")
        self.n_nodes = int(mu.size)
        if self.eps is None:
            med = float(np.median(std))
            base = med if med > 0 else float(np.median(np.abs(y - mu)))
            self.eps = self.eps_fraction * base if base > 0 else 1e-12
        a = np.abs(y - mu)
        s = a / (std + self.eps)
        keys: List[Tuple[str, np.ndarray]] = [("all", np.ones(s.size, bool))]
        if self.mondrian:
            if region is None:
                raise ValueError("mondrian calibration needs the region of every node")
            region = np.asarray(region).ravel()
            keys += [(REGION_NAMES[r], region == r) for r in np.unique(region)]
        self.scores = {name: _part_sketch(s[sel], groups[sel], self.max_per_part)
                       for name, sel in keys}
        self.baseline = _part_sketch(a, groups, self.max_per_part)
        self.data_source = data_source
        return self

    def _require(self) -> None:
        if not self.fitted:
            raise PrecompError("the calibrator is not fitted")

    def quantile(self, level: float, region: Optional[str] = None) -> float:
        """The score quantile q at coverage `level` (see the module docstring).

        Raises PrecompError when there are too few calibration parts for the
        level (k < level / (1 - level)): the interval would be unbounded.
        """
        self._require()
        if not 0 < level < 1:
            raise ValueError("level must lie in (0, 1)")
        key = region if (self.mondrian and region is not None) else "all"
        if key not in self.scores:
            key = "all"
        v, w, k = self.scores[key]
        return _part_quantile(v, w, k, level, key)

    def supports(self, level: float) -> bool:
        """Whether the calibration parts suffice for `level` (every key)."""
        self._require()
        return all(k >= min_parts(level) for _, _, k in self.scores.values())

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

    def baseline_halfwidth(self, level: float) -> float:
        """Half-width q0 [m] of the constant-width interval mu +- q0 at `level`,
        calibrated on the same parts: the reference a useful std must beat."""
        self._require()
        if self.baseline is None:
            raise PrecompError("this calibrator has no constant-width baseline")
        v, w, k = self.baseline
        return _part_quantile(v, w, k, level, "constant-width baseline")

    def baseline_intervals(self, mu: np.ndarray, level: float) -> Tuple[np.ndarray, np.ndarray]:
        q0 = self.baseline_halfwidth(level)
        mu = np.asarray(mu, dtype=float)
        return mu - q0, mu + q0

    def summary(self, levels: Sequence[float] = REPORT_LEVELS) -> Dict[str, Any]:
        out: Dict[str, Any] = {"method": "split conformal over whole parts, (k + 1) correction",
                               "eps_m": self.eps, "n_nodes": self.n_nodes,
                               "n_parts": self.n_parts, "mondrian": self.mondrian,
                               "data_source": self.data_source, "quantiles": {},
                               "baseline_halfwidth_m": {},
                               "min_parts": {f"{lv:g}": min_parts(lv) for lv in levels}}
        for lv in levels:
            try:
                out["quantiles"][f"{lv:g}"] = self.quantile(lv)
                out["baseline_halfwidth_m"][f"{lv:g}"] = self.baseline_halfwidth(lv)
            except PrecompError:
                out["quantiles"][f"{lv:g}"] = None
                out["baseline_halfwidth_m"][f"{lv:g}"] = None
        return out

    def to_state(self) -> Dict[str, Any]:
        return {"eps": self.eps, "eps_fraction": self.eps_fraction, "mondrian": self.mondrian,
                "max_per_part": self.max_per_part,
                "scores": {k: {"values": v, "weights": w, "n_parts": n}
                           for k, (v, w, n) in self.scores.items()},
                "baseline": None if self.baseline is None else
                {"values": self.baseline[0], "weights": self.baseline[1],
                 "n_parts": self.baseline[2]},
                "n_nodes": self.n_nodes, "n_parts": self.n_parts,
                "data_source": self.data_source}

    @classmethod
    def from_state(cls, state: Mapping[str, Any]) -> "ConformalCalibrator":
        c = cls(eps=state["eps"], eps_fraction=state["eps_fraction"],
                mondrian=state["mondrian"], max_per_part=state["max_per_part"])
        c.scores = {k: (np.asarray(d["values"]), np.asarray(d["weights"]), int(d["n_parts"]))
                    for k, d in state["scores"].items()}
        b = state.get("baseline")
        c.baseline = None if b is None else (np.asarray(b["values"]), np.asarray(b["weights"]),
                                             int(b["n_parts"]))
        c.n_nodes = int(state["n_nodes"])
        c.n_parts = int(state["n_parts"])
        c.data_source = state.get("data_source")
        return c


def intervals(mu: np.ndarray, std: np.ndarray, level: float, calibrator: ConformalCalibrator,
              region: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray]:
    """(lower, upper) of dz [m] at `level` from a fitted calibrator."""
    return calibrator.intervals(mu, std, level, region)


def partition_coverage(mu: np.ndarray, std: np.ndarray, y: np.ndarray, groups: np.ndarray, *,
                       level: float, n_partitions: int = 50, calibration_fraction: float = 0.5,
                       seed: int = 0, **calibrator_kwargs: Any) -> Dict[str, float]:
    """Coverage of repeated split-conformal calibration over whole parts.

    The held-out parts (`groups`, part ids) are split at random
    `n_partitions` times into calibration and test parts; each time a
    calibrator is fitted on the calibration parts and measured on the test
    parts: the mean over test parts of each part's coverage (what the
    guarantee is about: its expectation is >= level), the share of test
    parts covered at >= level, the worst part, and the mean width - also for
    the constant-width baseline calibrated on the same parts, so a std that
    carries no information shows as a width ratio near or above 1.
    """
    mu, std, y = (np.asarray(a, dtype=float).ravel() for a in (mu, std, y))
    groups = np.asarray(groups).astype(str).ravel()
    uniq = np.unique(groups)
    if uniq.size < 4:
        raise ValueError("partition_coverage needs at least four parts")
    rng = np.random.default_rng(seed)
    n_cal = int(min(max(2, round(calibration_fraction * uniq.size)), uniq.size - 2))
    if n_cal < min_parts(level):
        raise PrecompError(f"{n_cal} calibration parts per partition cannot support level "
                           f"{level:g} (needs {min_parts(level)}); raise calibration_fraction "
                           "or hold out more parts")
    rows = []
    for _ in range(int(n_partitions)):
        cal = np.isin(groups, rng.permutation(uniq)[:n_cal])
        c = ConformalCalibrator(**calibrator_kwargs).fit(mu[cal], std[cal], y[cal],
                                                         groups[cal])
        lo, hi = c.intervals(mu[~cal], std[~cal], level)
        q0 = c.baseline_halfwidth(level)
        yt, mt, gt = y[~cal], mu[~cal], groups[~cal]
        per = pd.Series((yt >= lo) & (yt <= hi)).groupby(gt).mean()
        per0 = pd.Series(np.abs(yt - mt) <= q0).groupby(gt).mean()
        rows.append((per.mean(), (per >= level).mean(), per.min(), float(np.mean(hi - lo)),
                     per0.mean(), 2.0 * q0))
    r = np.array(rows)
    return {"level": float(level), "mean": float(r[:, 0].mean()), "std": float(r[:, 0].std()),
            "min": float(r[:, 0].min()), "max": float(r[:, 0].max()),
            "share_parts_at_level": float(r[:, 1].mean()),
            "worst_part_median": float(np.median(r[:, 2])),
            "mean_width_m": float(r[:, 3].mean()),
            "baseline_mean": float(r[:, 4].mean()), "baseline_width_m": float(r[:, 5].mean()),
            "width_ratio": float(r[:, 3].mean() / r[:, 5].mean()),
            "n_partitions": int(n_partitions), "n_parts": int(uniq.size),
            "n_calibration_parts": n_cal}


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
           "min_parts", "REPORT_LEVELS"]
