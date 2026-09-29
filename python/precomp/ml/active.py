"""Which part or process to simulate (or try on the machine) next.

`rank_candidates` scores candidate (commanded surface, setup) pairs by the
information a new sample would add:

* the model's mean predictive std over the part [m] - where it is unsure;
* the envelope's part score - how far outside the training data it lies;

each normalised to [0, 1] over the candidates and combined with weights,
then selects greedily with max-min diversity in the standardised part-descriptor
space (global, process and material features): every pick maximises

    acquisition_i + diversity * d_i / median(d),

d_i the distance to the nearest part already chosen or in the training set,
so the batch does not spend several runs on near-identical parts.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from ..fea.setup import FormingSetup
from ..geometry.heightmap import HeightMap
from .features import global_features
from .surrogate import DeviationSurrogate


@dataclass
class Candidate:
    """A candidate run: the commanded surface, its setup and a label (e.g. the
    design point id); `info` is carried into the ranking table."""

    commanded: HeightMap
    setup: FormingSetup
    label: str
    info: Dict[str, Any] = field(default_factory=dict)


def _normalise(v: np.ndarray) -> np.ndarray:
    lo, hi = float(np.min(v)), float(np.max(v))
    return np.zeros_like(v) if hi <= lo else (v - lo) / (hi - lo)


def rank_candidates(surrogate: DeviationSurrogate, candidates: Sequence[Candidate], *,
                    n_select: Optional[int] = None, weight_std: float = 1.0,
                    weight_ood: float = 1.0, diversity: float = 0.5) -> pd.DataFrame:
    """Rank candidates for the next simulations (see the module docstring).

    Returns one row per candidate: label, mean_std_m (over the part),
    ood_part_score, in_envelope, acquisition, min_distance (to the training
    parts and the earlier picks, standardised units), rank (1 = first pick;
    candidates beyond `n_select` keep rank NaN), plus the candidates' info.
    """
    cands = list(candidates)
    if not cands:
        raise ValueError("no candidates")
    if surrogate.ood is None:
        raise ValueError("ranking needs the surrogate's OOD envelope (train it with envelope)")
    n_select = len(cands) if n_select is None else int(min(n_select, len(cands)))
    rows, D = [], []
    for c in cands:
        mu, sd = surrogate.predict_deviation(c.commanded, c.setup)
        fm = surrogate.feature_maps(c.commanded, c.setup)
        a = surrogate.assess(c.commanded, c.setup)
        d, _ = global_features(c.commanded, c.setup, surrogate.features)
        D.append(d)
        rows.append({"label": c.label, "mean_std_m": float(sd[fm.part].mean()),
                     "ood_part_score": float(a["part_score"]),
                     "in_envelope": bool(a["in_envelope"]),
                     "ood_reasons": ";".join(a["reasons"]), **c.info})
    frame = pd.DataFrame(rows)
    acq = weight_std * _normalise(frame["mean_std_m"].to_numpy()) + \
        weight_ood * _normalise(frame["ood_part_score"].to_numpy())
    frame["acquisition"] = acq
    env = surrogate.ood
    Z = np.clip((np.stack(D) - env.d_mean) / env.d_std, -1e6, 1e6)
    ref = env.parts.reference
    dmin = np.sqrt(((Z[:, None, :] - ref[None, :, :]) ** 2).sum(-1)).min(axis=1) \
        if len(ref) else np.full(len(Z), np.inf)
    scale = float(np.median(dmin)) if np.isfinite(dmin).all() and np.median(dmin) > 0 else 1.0
    chosen: List[int] = []
    ranks = np.full(len(cands), np.nan)
    at_pick = np.full(len(cands), np.nan)
    remaining = set(range(len(cands)))
    for r in range(n_select):
        idx = sorted(remaining)
        score = acq[idx] + diversity * dmin[idx] / scale
        i = idx[int(np.argmax(score))]
        chosen.append(i)
        ranks[i] = r + 1
        at_pick[i] = dmin[i]
        remaining.discard(i)
        dmin = np.minimum(dmin, np.sqrt(((Z - Z[i]) ** 2).sum(-1)))
    frame["min_distance"] = at_pick
    frame["rank"] = ranks
    return frame.sort_values(["rank", "acquisition"], ascending=[True, False],
                             na_position="last").reset_index(drop=True)


__all__ = ["Candidate", "rank_candidates"]
