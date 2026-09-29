"""Deviation models: dz [m] from features, with an uncertainty.

Two protocols:

* `DeviationModel` (a point model) - ``fit(X, y, groups)`` on a feature table
  (rows = grid nodes, `groups` = sample ids), ``predict(X) -> (mean, std)`` in
  metres, ``feature_names``, ``to_state()`` / ``from_state(state)``. The
  implementations are `GBMEnsemble`, `MLPEnsemble`, `ResidualModel` and
  `TransferModel`.
* `FieldDeviationModel` - trained on whole parts, ``fit_samples(samples)``
  and ``predict_field(commanded, setup) -> (mean, std)`` on the commanded grid:
  `FieldUNet`, which sees the whole part at once and so can represent the
  non-local unclamping springback that a point model only reaches through its
  global descriptors.

`precomp.ml.surrogate.DeviationSurrogate` turns either into what the
compensation code consumes (``predict_deviation(commanded, setup)``).

Determinism: every model takes a `seed`; bootstrap draws, initialisations,
shuffles and dropout masks derive from it, so a fit on the same table with the
same seed gives the same model (on the same library versions and hardware).
Threads: `n_threads` caps OpenMP (scikit-learn) and torch threads during fit and
predict; None uses $PRECOMP_ML_THREADS, else the libraries' defaults.

The torch models import torch lazily and raise ImportError with the install
command when it is missing; they run on the CPU.
"""

from __future__ import annotations

import contextlib
import copy
import math
import os
from dataclasses import dataclass, field
from typing import (Any, Dict, Iterator, List, Mapping, Optional, Protocol, Sequence, Tuple,
                    runtime_checkable)

import numpy as np
# Imported here, not lazily: threadpoolctl caps only the OpenMP runtimes that
# are already loaded when a limit is entered, and scikit-learn loads its own
# with this module. A lazy import inside a limited block ran unlimited and,
# on a loaded machine, twenty times slower.
from sklearn.ensemble import HistGradientBoostingRegressor

from .._util import PrecompError
from ..geometry.heightmap import Grid, HeightMap
from .features import (DEFAULT_CONFIG, PRIOR_FEATURE, FeatureConfig, as_setup,
                       global_features)

THREADS_ENV = "PRECOMP_ML_THREADS"


# ---------------------------------------------------------------------------
# Protocols and helpers
# ---------------------------------------------------------------------------
@runtime_checkable
class DeviationModel(Protocol):
    """A point model: rows of features -> dz [m] with a standard deviation [m]."""

    kind: str
    feature_names: List[str]

    def fit(self, X: np.ndarray, y: np.ndarray, groups: np.ndarray, *,
            feature_names: Optional[Sequence[str]] = None) -> "DeviationModel": ...

    def predict(self, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]: ...

    def to_state(self) -> Dict[str, Any]: ...


@runtime_checkable
class FieldDeviationModel(Protocol):
    """A field model: a whole commanded surface -> dz map [m] and its std [m]."""

    kind: str

    def fit_samples(self, samples: Sequence[Any]) -> "FieldDeviationModel": ...

    def predict_field(self, commanded: HeightMap, setup: Any
                      ) -> Tuple[np.ndarray, np.ndarray]: ...

    def to_state(self) -> Dict[str, Any]: ...


def resolve_threads(n_threads: Optional[int]) -> Optional[int]:
    if n_threads is not None:
        if int(n_threads) < 1:
            raise ValueError("n_threads must be >= 1")
        return int(n_threads)
    env = os.environ.get(THREADS_ENV)
    if env:
        try:
            return max(1, int(env))
        except ValueError as exc:
            raise ValueError(f"${THREADS_ENV} must be an integer, got {env!r}") from exc
    return None


@contextlib.contextmanager
def thread_limit(n_threads: Optional[int]) -> Iterator[None]:
    """Cap OpenMP/BLAS threads (threadpoolctl) and torch's intra-op threads."""
    n = resolve_threads(n_threads)
    if n is None:
        yield
        return
    from threadpoolctl import threadpool_limits

    torch_prev = None
    try:
        import sys
        if "torch" in sys.modules:
            import torch
            torch_prev = torch.get_num_threads()
            torch.set_num_threads(n)
    except Exception:  # torch is optional; a failure to cap it is not fatal
        torch_prev = None
    try:
        with threadpool_limits(limits=n):
            yield
    finally:
        if torch_prev is not None:
            import torch
            torch.set_num_threads(torch_prev)


def _check_table(X: np.ndarray, y: Optional[np.ndarray] = None,
                 groups: Optional[np.ndarray] = None) -> Tuple[np.ndarray, ...]:
    X = np.asarray(X, dtype=float)
    if X.ndim != 2:
        raise ValueError(f"X must be 2-D (rows, features), got shape {X.shape}")
    if not np.all(np.isfinite(X)):
        raise ValueError("X contains NaN or infinity")
    out: List[np.ndarray] = [X]
    if y is not None:
        y = np.asarray(y, dtype=float).ravel()
        if y.shape[0] != X.shape[0]:
            raise ValueError(f"y has {y.shape[0]} rows, X has {X.shape[0]}")
        if not np.all(np.isfinite(y)):
            raise ValueError("y contains NaN or infinity")
        out.append(y)
    if groups is not None:
        groups = np.asarray(groups)
        if groups.shape[0] != X.shape[0]:
            raise ValueError(f"groups has {groups.shape[0]} rows, X has {X.shape[0]}")
        out.append(groups)
    return tuple(out)


def _check_names(model: Any, X: np.ndarray) -> None:
    names = getattr(model, "feature_names", None)
    if names and X.shape[1] != len(names):
        raise ValueError(f"{type(model).__name__} was fitted on {len(names)} features, "
                         f"X has {X.shape[1]} columns")


def _require_torch(what: str):
    try:
        import torch  # noqa: F401
    except ImportError as exc:
        raise ImportError(f"{what} needs torch: pip install -e '.[torch]'") from exc
    return torch


#: Default target scale: dz is learned divided by the elastic unloading
#: curvature 3 sigma_0.2 / (E t) of the sheet, the classical springback scale.
DEFAULT_TARGET_SCALE = "elastic_curvature"


def _target_scale(model: Any, X: np.ndarray) -> np.ndarray:
    """Per-row divisor of the target: the feature column `target_scale`
    (ones when it is None). ValueError when the column is absent or not > 0."""
    name = model.params.get("target_scale")
    if name is None:
        return np.ones(len(X))
    if name not in model.feature_names:
        raise ValueError(f"target_scale {name!r} is not among the feature names; pass "
                         "feature_names to fit, or target_scale=None")
    s = X[:, model.feature_names.index(name)]
    if not np.all(s > 0):
        raise ValueError(f"target_scale column {name!r} must be > 0 on every row")
    return s


def _group_bootstrap_rows(groups: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Row indices of a bootstrap resample of the groups (whole samples, drawn
    with replacement; a group drawn twice contributes its rows twice)."""
    uniq, inv = np.unique(groups, return_inverse=True)
    counts = np.bincount(rng.integers(0, uniq.size, size=uniq.size), minlength=uniq.size)
    return np.repeat(np.arange(len(groups)), counts[inv])


def _grouped_holdout(groups: np.ndarray, fraction: float, rng: np.random.Generator
                     ) -> Tuple[np.ndarray, np.ndarray]:
    uniq = np.unique(groups)
    if uniq.size < 2 or fraction <= 0:
        rows = np.arange(len(groups))
        return rows, rows[:0]
    n_val = int(min(max(1, round(fraction * uniq.size)), uniq.size - 1))
    val = set(rng.permutation(uniq)[:n_val])
    is_val = np.array([g in val for g in groups])
    return np.flatnonzero(~is_val), np.flatnonzero(is_val)


# ---------------------------------------------------------------------------
# Gradient-boosting ensemble
# ---------------------------------------------------------------------------
class GBMEnsemble:
    """M HistGradientBoostingRegressor members on bootstraps of SAMPLES.

    Each member is fitted on a resample of the groups (whole samples, with
    replacement) with its own seed; the prediction is the members' mean and
    the std their spread - an epistemic uncertainty that conformal
    calibration (`uncertainty.ConformalCalibrator`) scales to a coverage.
    `quantiles` (e.g. (0.05, 0.95)) adds quantile-loss members on the whole
    table for a direct 5/95 % band (`predict_quantiles`).

    `target_scale` names a feature column the target is divided by before
    fitting (and the prediction multiplied by): by default the elastic
    unloading curvature 3 sigma_0.2 / (E t), which takes the leading
    material and thickness dependence of springback out of what the trees
    must learn from a few dozen parts. None fits dz as it is.
    """

    kind = "point"

    def __init__(self, n_members: int = 8, *, max_iter: int = 300,
                 learning_rate: float = 0.08, max_leaf_nodes: int = 31,
                 min_samples_leaf: int = 40, l2_regularization: float = 1e-3,
                 max_features: float = 0.8, quantiles: Optional[Sequence[float]] = None,
                 target_scale: Optional[str] = DEFAULT_TARGET_SCALE,
                 seed: int = 0, n_threads: Optional[int] = None):
        if n_members < 1:
            raise ValueError("n_members must be >= 1")
        self.params = dict(n_members=int(n_members), max_iter=int(max_iter),
                           learning_rate=float(learning_rate),
                           max_leaf_nodes=int(max_leaf_nodes),
                           min_samples_leaf=int(min_samples_leaf),
                           l2_regularization=float(l2_regularization),
                           max_features=float(max_features),
                           quantiles=None if quantiles is None else [float(q) for q in quantiles],
                           target_scale=target_scale, seed=int(seed))
        for q in self.params["quantiles"] or []:
            if not 0 < q < 1:
                raise ValueError(f"quantiles must lie in (0, 1), got {q}")
        self.n_threads = n_threads
        self.feature_names: List[str] = []
        self.members: List[Any] = []
        self.quantile_members: Dict[float, Any] = {}
        self.y_mean = 0.0
        self.y_scale = 1.0
        self.n_train_groups = 0

    def _estimator(self, seed: int, **extra: Any):
        p = self.params
        kw = dict(max_iter=p["max_iter"], learning_rate=p["learning_rate"],
                  max_leaf_nodes=p["max_leaf_nodes"], min_samples_leaf=p["min_samples_leaf"],
                  l2_regularization=p["l2_regularization"], early_stopping=False,
                  random_state=seed)
        if p["max_features"] < 1.0:
            kw["max_features"] = p["max_features"]
        kw.update(extra)
        return HistGradientBoostingRegressor(**kw)

    def fit(self, X: np.ndarray, y: np.ndarray, groups: np.ndarray, *,
            feature_names: Optional[Sequence[str]] = None) -> "GBMEnsemble":
        X, y, groups = _check_table(X, y, groups)
        if len(y) == 0:
            raise ValueError("GBMEnsemble.fit: empty table")
        self.feature_names = list(feature_names) if feature_names is not None else \
            [f"x{i}" for i in range(X.shape[1])]
        yn = y / _target_scale(self, X)
        self.y_mean = float(yn.mean())
        self.y_scale = float(yn.std()) or 1.0
        t = (yn - self.y_mean) / self.y_scale
        rng = np.random.default_rng(self.params["seed"])
        self.n_train_groups = int(np.unique(groups).size)
        self.members = []
        with thread_limit(self.n_threads):
            for m in range(self.params["n_members"]):
                rows = _group_bootstrap_rows(groups, rng) if self.n_train_groups > 1 else \
                    np.arange(len(y))
                est = self._estimator(self.params["seed"] * 1009 + m)
                est.fit(X[rows], t[rows])
                self.members.append(est)
            self.quantile_members = {}
            for q in self.params["quantiles"] or []:
                est = self._estimator(self.params["seed"] * 1009 + 997, loss="quantile",
                                      quantile=q)
                est.fit(X, t)
                self.quantile_members[q] = est
        return self

    def _require_fitted(self) -> None:
        if not self.members:
            raise PrecompError("GBMEnsemble is not fitted")

    def predict_members(self, X: np.ndarray) -> np.ndarray:
        """(M, n) member predictions [m]."""
        self._require_fitted()
        (X,) = _check_table(X)
        _check_names(self, X)
        with thread_limit(self.n_threads):
            P = np.stack([est.predict(X) for est in self.members])
        return (self.y_mean + self.y_scale * P) * _target_scale(self, X)

    def predict(self, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        P = self.predict_members(X)
        std = P.std(axis=0, ddof=1) if len(P) > 1 else np.zeros(P.shape[1])
        return P.mean(axis=0), std

    def predict_quantiles(self, X: np.ndarray) -> Dict[float, np.ndarray]:
        """{q: dz [m]} of the quantile members (empty without `quantiles`)."""
        self._require_fitted()
        (X,) = _check_table(X)
        s = _target_scale(self, X)
        with thread_limit(self.n_threads):
            return {q: (self.y_mean + self.y_scale * est.predict(X)) * s
                    for q, est in self.quantile_members.items()}

    def to_state(self) -> Dict[str, Any]:
        self._require_fitted()
        return {"class": "GBMEnsemble", "params": dict(self.params),
                "feature_names": list(self.feature_names), "members": list(self.members),
                "quantile_members": {str(q): e for q, e in self.quantile_members.items()},
                "y_mean": self.y_mean, "y_scale": self.y_scale,
                "n_train_groups": self.n_train_groups}

    @classmethod
    def from_state(cls, state: Mapping[str, Any]) -> "GBMEnsemble":
        m = cls(**state["params"])
        m.feature_names = list(state["feature_names"])
        m.members = list(state["members"])
        m.quantile_members = {float(q): e for q, e in state["quantile_members"].items()}
        m.y_mean, m.y_scale = float(state["y_mean"]), float(state["y_scale"])
        m.n_train_groups = int(state["n_train_groups"])
        return m


# ---------------------------------------------------------------------------
# Deep ensemble of MLPs (torch)
# ---------------------------------------------------------------------------
LOGVAR_RANGE = (-14.0, 6.0)


def _standardiser(X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    tiny = 1e-12 * np.maximum(np.abs(mu), 1.0)
    sd = np.where(sd > tiny, sd, 1.0)                  # constant columns pass through centred
    return mu, sd


def _state_to_numpy(sd: Mapping[str, Any]) -> Dict[str, np.ndarray]:
    return {k: v.detach().cpu().numpy().copy() for k, v in sd.items()}


def _state_from_numpy(sd: Mapping[str, Any]):
    import torch
    return {k: torch.as_tensor(np.asarray(v)) for k, v in sd.items()}


class MLPEnsemble:
    """A deep ensemble of small MLPs with a Gaussian negative-log-likelihood head.

    Each member maps standardised features to (mean, log-variance) of the
    standardised target and is trained with Adam on the Gaussian NLL (after
    `mse_warmup` epochs of squared error on the mean), with early stopping on a
    grouped validation split (whole samples) and its own seed. The prediction
    is the mixture: mean of the means, variance = mean of (variance + mean^2)
    - mean^2. Feature standardisation, the target scale (`target_scale`, as
    for `GBMEnsemble`) and the target standardisation are part of the state.
    CPU only.
    """

    kind = "point"

    def __init__(self, n_members: int = 5, *, hidden: Sequence[int] = (128, 128),
                 epochs: int = 40, batch_size: int = 1024, lr: float = 2e-3,
                 weight_decay: float = 1e-5, patience: int = 6, val_fraction: float = 0.2,
                 mse_warmup: int = 3, bootstrap: bool = False,
                 target_scale: Optional[str] = DEFAULT_TARGET_SCALE, seed: int = 0,
                 n_threads: Optional[int] = None):
        if n_members < 1:
            raise ValueError("n_members must be >= 1")
        self.params = dict(n_members=int(n_members), hidden=[int(h) for h in hidden],
                           epochs=int(epochs), batch_size=int(batch_size), lr=float(lr),
                           weight_decay=float(weight_decay), patience=int(patience),
                           val_fraction=float(val_fraction), mse_warmup=int(mse_warmup),
                           bootstrap=bool(bootstrap), target_scale=target_scale,
                           seed=int(seed))
        self.n_threads = n_threads
        self.feature_names: List[str] = []
        self.x_mean = self.x_std = None
        self.y_mean = 0.0
        self.y_scale = 1.0
        self.states: List[Dict[str, np.ndarray]] = []
        self.history: List[Dict[str, Any]] = []
        self._nets: Optional[list] = None

    def _net(self, n_in: int):
        import torch.nn as nn
        layers: List[Any] = []
        prev = n_in
        for h in self.params["hidden"]:
            layers += [nn.Linear(prev, h), nn.SiLU()]
            prev = h
        layers.append(nn.Linear(prev, 2))
        return nn.Sequential(*layers)

    @staticmethod
    def _nll(out, t, mse: bool):
        import torch
        mu = out[:, 0]
        if mse:
            return torch.mean((mu - t) ** 2)
        lv = out[:, 1].clamp(*LOGVAR_RANGE)
        return torch.mean(0.5 * (lv + (t - mu) ** 2 * torch.exp(-lv)))

    def fit(self, X: np.ndarray, y: np.ndarray, groups: np.ndarray, *,
            feature_names: Optional[Sequence[str]] = None) -> "MLPEnsemble":
        torch = _require_torch("MLPEnsemble")
        X, y, groups = _check_table(X, y, groups)
        if len(y) < 2:
            raise ValueError("MLPEnsemble.fit needs at least two rows")
        self.feature_names = list(feature_names) if feature_names is not None else \
            [f"x{i}" for i in range(X.shape[1])]
        self.x_mean, self.x_std = _standardiser(X)
        yn = y / _target_scale(self, X)
        self.y_mean = float(yn.mean())
        self.y_scale = float(yn.std()) or 1.0
        Xs = ((X - self.x_mean) / self.x_std).astype(np.float32)
        ts = ((yn - self.y_mean) / self.y_scale).astype(np.float32)
        p = self.params
        self.states, self.history = [], []
        with thread_limit(self.n_threads):
            for m in range(p["n_members"]):
                seed = p["seed"] * 1009 + m
                rng = np.random.default_rng(seed)
                tr, va = _grouped_holdout(groups, p["val_fraction"], rng)
                if p["bootstrap"]:
                    tr = tr[_group_bootstrap_rows(groups[tr], rng)]
                torch.manual_seed(seed)
                net = self._net(X.shape[1])
                opt = torch.optim.Adam(net.parameters(), lr=p["lr"],
                                       weight_decay=p["weight_decay"])
                gen = torch.Generator().manual_seed(seed)
                Xt, Tt = torch.from_numpy(Xs[tr]), torch.from_numpy(ts[tr])
                Xv, Tv = torch.from_numpy(Xs[va]), torch.from_numpy(ts[va])
                best, best_state, bad = math.inf, None, 0
                epochs_run = 0
                for epoch in range(p["epochs"]):
                    mse = epoch < p["mse_warmup"]
                    net.train()
                    perm = torch.randperm(len(Tt), generator=gen)
                    for i in range(0, len(Tt), p["batch_size"]):
                        b = perm[i:i + p["batch_size"]]
                        opt.zero_grad()
                        loss = self._nll(net(Xt[b]), Tt[b], mse)
                        loss.backward()
                        opt.step()
                    epochs_run = epoch + 1
                    if mse:
                        continue
                    net.eval()
                    with torch.no_grad():
                        vl = float(self._nll(net(Xv), Tv, False)) if len(Tv) else float(
                            self._nll(net(Xt), Tt, False))
                    if vl < best - 1e-4:
                        best, best_state, bad = vl, copy.deepcopy(net.state_dict()), 0
                    else:
                        bad += 1
                        if bad >= p["patience"]:
                            break
                if best_state is None:
                    best_state = net.state_dict()
                self.states.append(_state_to_numpy(best_state))
                self.history.append({"member": m, "epochs": epochs_run, "best_val_nll": best,
                                     "n_train_rows": int(len(tr)), "n_val_rows": int(len(va))})
        self._nets = None
        return self

    def _networks(self):
        if not self.states:
            raise PrecompError("MLPEnsemble is not fitted")
        if self._nets is None:
            _require_torch("MLPEnsemble")
            nets = []
            for sd in self.states:
                net = self._net(len(self.feature_names))
                net.load_state_dict(_state_from_numpy(sd))
                net.eval()
                nets.append(net)
            self._nets = nets
        return self._nets

    def predict(self, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        import torch
        (X,) = _check_table(X)
        _check_names(self, X)
        nets = self._networks()
        Xs = torch.from_numpy(((X - self.x_mean) / self.x_std).astype(np.float32))
        mus, vars_ = [], []
        with thread_limit(self.n_threads), torch.no_grad():
            for net in nets:
                out = torch.cat([net(Xs[i:i + 65536]) for i in range(0, len(Xs), 65536)]) \
                    if len(Xs) else torch.zeros((0, 2))
                mus.append(out[:, 0].double().numpy())
                vars_.append(np.exp(out[:, 1].clamp(*LOGVAR_RANGE).double().numpy()))
        M = np.stack(mus)
        V = np.stack(vars_)
        mean = M.mean(axis=0)
        var = np.maximum((V + M * M).mean(axis=0) - mean * mean, 0.0)
        s = _target_scale(self, X)
        return (self.y_mean + self.y_scale * mean) * s, self.y_scale * np.sqrt(var) * s

    def to_state(self) -> Dict[str, Any]:
        if not self.states:
            raise PrecompError("MLPEnsemble is not fitted")
        return {"class": "MLPEnsemble", "params": dict(self.params),
                "feature_names": list(self.feature_names), "x_mean": self.x_mean,
                "x_std": self.x_std, "y_mean": self.y_mean, "y_scale": self.y_scale,
                "history": list(self.history),
                "torch_state_dicts": {f"member{m}": sd for m, sd in enumerate(self.states)}}

    @classmethod
    def from_state(cls, state: Mapping[str, Any]) -> "MLPEnsemble":
        m = cls(**state["params"])
        m.feature_names = list(state["feature_names"])
        m.x_mean, m.x_std = np.asarray(state["x_mean"]), np.asarray(state["x_std"])
        m.y_mean, m.y_scale = float(state["y_mean"]), float(state["y_scale"])
        m.history = list(state.get("history", []))
        sds = state["torch_state_dicts"]
        m.states = [dict(sds[f"member{i}"]) for i in range(len(sds))]
        return m


# ---------------------------------------------------------------------------
# Field U-Net (torch)
# ---------------------------------------------------------------------------
#: Height scale [m] of the commanded-surface channel of the U-Net.
UNET_Z_SCALE = 0.03


def _build_unet(c_in: int, chans: Sequence[int], dropout: float):
    import torch
    import torch.nn as nn
    import torch.nn.functional as F

    def block(a: int, b: int) -> nn.Module:
        g = 4 if b % 4 == 0 else 1
        return nn.Sequential(
            nn.Conv2d(a, b, 3, padding=1, padding_mode="replicate"), nn.GroupNorm(g, b),
            nn.SiLU(),
            nn.Conv2d(b, b, 3, padding=1, padding_mode="replicate"), nn.GroupNorm(g, b),
            nn.SiLU())

    class UNet(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            self.enc = nn.ModuleList()
            prev = c_in
            for c in chans:
                self.enc.append(block(prev, c))
                prev = c
            self.dec = nn.ModuleList()
            for k in range(len(chans) - 1, 0, -1):
                self.dec.append(block(chans[k] + chans[k - 1], chans[k - 1]))
            self.drop = nn.Dropout2d(dropout)
            self.head = nn.Conv2d(chans[0], 2, 1)

        def forward(self, x):
            skips = []
            for i, e in enumerate(self.enc):
                if i > 0:
                    x = F.avg_pool2d(x, 2, ceil_mode=True)
                x = e(x)
                skips.append(x)
            x = self.drop(x)
            for j, d in enumerate(self.dec):
                skip = skips[-2 - j]
                x = F.interpolate(x, size=skip.shape[-2:], mode="bilinear", align_corners=False)
                x = d(torch.cat([x, skip], dim=1))
                x = self.drop(x)
            return self.head(x)

    return UNet()


class FieldUNet:
    """A small U-Net from the commanded surface to the dz map (+ log-variance).

    Inputs, resampled to a fixed training grid of `resolution`^2 nodes over
    the square window [-W, W]^2 (W = `window`, default the largest half-width
    covered by every training grid): the commanded height / 0.03 m, the depth
    fraction, the wall angle at one tool radius, the distance to the clamped
    frame / W, the part indicator, and the part-level descriptors (global,
    process and material features, standardised) broadcast as constant
    channels. Output: dz / scale and its log-variance, trained with a Gaussian
    NLL weighted 1 on the part and `flange_weight` elsewhere; the scale is
    the part's `target_scale` descriptor (as for `GBMEnsemble`) times a
    constant fixed at fit. Uncertainty by
    Monte-Carlo dropout (`mc_samples` passes, seeded): the variance is the mean
    predicted variance plus the variance of the passes. The prediction is
    resampled bilinearly back to the commanded grid; nodes outside the window
    take the nearest window value.
    """

    kind = "field"

    def __init__(self, *, resolution: int = 64, channels: Sequence[int] = (16, 32, 64),
                 dropout: float = 0.1, epochs: int = 200, batch_size: int = 8,
                 lr: float = 2e-3, weight_decay: float = 1e-5, mse_warmup: int = 20,
                 mc_samples: int = 16, flange_weight: float = 0.2,
                 window: Optional[float] = None, features: FeatureConfig = DEFAULT_CONFIG,
                 target_scale: Optional[str] = DEFAULT_TARGET_SCALE,
                 seed: int = 0, n_threads: Optional[int] = None):
        if resolution < 8:
            raise ValueError("resolution must be >= 8")
        self.params = dict(resolution=int(resolution), channels=[int(c) for c in channels],
                           dropout=float(dropout), epochs=int(epochs),
                           batch_size=int(batch_size), lr=float(lr),
                           weight_decay=float(weight_decay), mse_warmup=int(mse_warmup),
                           mc_samples=int(mc_samples), flange_weight=float(flange_weight),
                           window=None if window is None else float(window),
                           target_scale=target_scale, seed=int(seed))
        self.features = features
        self.n_threads = n_threads
        self.scalar_names: List[str] = []
        self.s_mean = self.s_std = None
        self.y_scale = 1.0
        self.state: Optional[Dict[str, np.ndarray]] = None
        self.history: List[Dict[str, float]] = []
        self._net = None

    @property
    def feature_names(self) -> List[str]:
        return ["z_norm", "depth_frac", "wall_angle_R", "clamp_dist_norm", "part"] + \
            list(self.scalar_names)

    # protocol guards
    def fit(self, *args: Any, **kwargs: Any):
        raise TypeError("FieldUNet trains on whole parts: use fit_samples(samples)")

    def predict(self, *args: Any, **kwargs: Any):
        raise TypeError("FieldUNet predicts whole fields: use predict_field(commanded, setup)")

    @property
    def grid(self) -> Grid:
        W = self.params["window"]
        if W is None:
            raise PrecompError("FieldUNet is not fitted (no window)")
        n = self.params["resolution"]
        return Grid(-W, -W, n, n, 2.0 * W / (n - 1))

    def _maps(self, commanded: HeightMap, setup: Any) -> Tuple[np.ndarray, np.ndarray,
                                                               HeightMap]:
        s = as_setup(setup)
        grid = self.grid
        cmd = commanded.resample(grid) if not commanded.grid.matches(grid) else commanded
        eps = self.features.part_eps
        depth = np.maximum(-cmd.z, 0.0)
        part = depth > eps
        dmax = float(depth[part].max()) if part.any() else 1.0
        X, Y = grid.mesh()
        W = self.params["window"]
        maps = np.stack([cmd.z / UNET_Z_SCALE, depth / dmax, cmd.wall_angle(s.tool_radius),
                         (s.free_half_width - np.maximum(np.abs(X), np.abs(Y))) / W,
                         part.astype(float)])
        scal, names = global_features(commanded, s, self.features)
        return maps, scal, cmd

    def _part_scale(self, scal: np.ndarray, names: Sequence[str]) -> float:
        name = self.params.get("target_scale")
        if name is None:
            return 1.0
        if name not in names:
            raise ValueError(f"target_scale {name!r} is not a part descriptor; known: "
                             f"{list(names)}")
        v = float(scal[list(names).index(name)])
        if not v > 0:
            raise ValueError(f"target_scale {name!r} must be > 0, got {v}")
        return v

    def _input(self, maps: np.ndarray, scal: np.ndarray) -> np.ndarray:
        s = (scal - self.s_mean) / self.s_std
        n = self.params["resolution"]
        return np.concatenate([maps, np.broadcast_to(s[:, None, None], (len(s), n, n))]
                              ).astype(np.float32)

    def fit_samples(self, samples: Sequence[Any]) -> "FieldUNet":
        torch = _require_torch("FieldUNet")
        samples = list(samples)
        if not samples:
            raise ValueError("FieldUNet.fit_samples: no samples")
        p = self.params
        if p["window"] is None:
            p["window"] = float(min(min(-s.grid.x0, s.grid.xmax, -s.grid.y0, s.grid.ymax)
                                    for s in samples))
        grid = self.grid
        maps, scals, targets, weights = [], [], [], []
        for smp in samples:
            m, sc, cmd = self._maps(smp.commanded, smp.setup)
            dz = HeightMap(smp.grid, smp.dz, smp.valid).resample(grid)
            maps.append(m)
            scals.append(sc)
            targets.append(dz.z)
            w = np.where(m[4] > 0.5, 1.0, p["flange_weight"]) * dz.mask
            weights.append(w)
        S = np.stack(scals)
        _, self.scalar_names = global_features(samples[0].commanded, samples[0].setup,
                                               self.features)
        self.s_mean, self.s_std = _standardiser(S)
        Y = np.stack([t / self._part_scale(sc, self.scalar_names)
                      for t, sc in zip(targets, scals)])
        self.y_scale = float(np.sqrt(np.mean(Y ** 2))) or 1.0
        Xin = np.stack([self._input(m, sc) for m, sc in zip(maps, scals)])
        T = (Y / self.y_scale).astype(np.float32)
        Wt = np.stack(weights).astype(np.float32)
        seed = p["seed"]
        with thread_limit(self.n_threads):
            torch.manual_seed(seed)
            net = _build_unet(Xin.shape[1], p["channels"], p["dropout"])
            opt = torch.optim.Adam(net.parameters(), lr=p["lr"], weight_decay=p["weight_decay"])
            sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, max(1, p["epochs"]))
            gen = torch.Generator().manual_seed(seed)
            Xt, Tt, Wtt = torch.from_numpy(Xin), torch.from_numpy(T), torch.from_numpy(Wt)
            self.history = []
            for epoch in range(p["epochs"]):
                net.train()
                perm = torch.randperm(len(Xt), generator=gen)
                total = 0.0
                for i in range(0, len(Xt), p["batch_size"]):
                    b = perm[i:i + p["batch_size"]]
                    out = net(Xt[b])
                    mu, lv = out[:, 0], out[:, 1].clamp(*LOGVAR_RANGE)
                    w = Wtt[b]
                    if epoch < p["mse_warmup"]:
                        per = (mu - Tt[b]) ** 2
                    else:
                        per = 0.5 * (lv + (Tt[b] - mu) ** 2 * torch.exp(-lv))
                    loss = (per * w).sum() / w.sum().clamp_min(1.0)
                    opt.zero_grad()
                    loss.backward()
                    opt.step()
                    total += float(loss.detach()) * len(b)
                sched.step()
                self.history.append({"epoch": epoch, "loss": total / len(Xt)})
        self.state = _state_to_numpy(net.state_dict())
        self._net = None
        return self

    def _network(self):
        if self.state is None:
            raise PrecompError("FieldUNet is not fitted")
        if self._net is None:
            _require_torch("FieldUNet")
            net = _build_unet(5 + len(self.scalar_names), self.params["channels"],
                              self.params["dropout"])
            net.load_state_dict(_state_from_numpy(self.state))
            self._net = net
        return self._net

    def predict_window(self, commanded: HeightMap, setup: Any) -> Tuple[np.ndarray, np.ndarray]:
        """(mean, std) [m] on the training grid (`grid`)."""
        import torch
        net = self._network()
        maps, scal, _ = self._maps(commanded, setup)
        x = torch.from_numpy(self._input(maps, scal)[None])
        T = max(1, self.params["mc_samples"])
        mus, vs = [], []
        with thread_limit(self.n_threads), torch.random.fork_rng(), torch.no_grad():
            torch.manual_seed(self.params["seed"] + 7919)
            net.train(T > 1 and self.params["dropout"] > 0)   # MC dropout (no batch norm)
            for _ in range(T):
                out = net(x)[0]
                mus.append(out[0].double().numpy())
                vs.append(np.exp(out[1].clamp(*LOGVAR_RANGE).double().numpy()))
            net.eval()
        M, V = np.stack(mus), np.stack(vs)
        mean = M.mean(axis=0)
        var = V.mean(axis=0) + M.var(axis=0)
        k = self.y_scale * self._part_scale(scal, self.scalar_names)
        return k * mean, k * np.sqrt(var)

    def predict_field(self, commanded: HeightMap, setup: Any) -> Tuple[np.ndarray, np.ndarray]:
        mean, std = self.predict_window(commanded, setup)
        grid = self.grid
        if commanded.grid.matches(grid):
            return mean, std
        mu = HeightMap(grid, mean).resample(commanded.grid)
        sd = HeightMap(grid, std).resample(commanded.grid)
        return mu.z, sd.z

    def fraction_outside_window(self, grid: Grid) -> float:
        """Share of `grid`'s nodes outside the training window (extrapolated)."""
        W = self.params["window"]
        X, Y = grid.mesh()
        return float(np.mean((np.abs(X) > W * (1 + 1e-9)) | (np.abs(Y) > W * (1 + 1e-9))))

    def to_state(self) -> Dict[str, Any]:
        if self.state is None:
            raise PrecompError("FieldUNet is not fitted")
        return {"class": "FieldUNet", "params": dict(self.params),
                "features": self.features.to_dict(), "scalar_names": list(self.scalar_names),
                "s_mean": self.s_mean, "s_std": self.s_std, "y_scale": self.y_scale,
                "history": list(self.history), "torch_state_dicts": {"unet": self.state}}

    @classmethod
    def from_state(cls, state: Mapping[str, Any]) -> "FieldUNet":
        m = cls(features=FeatureConfig.from_dict(state["features"]), **state["params"])
        m.scalar_names = list(state["scalar_names"])
        m.s_mean, m.s_std = np.asarray(state["s_mean"]), np.asarray(state["s_std"])
        m.y_scale = float(state["y_scale"])
        m.history = list(state.get("history", []))
        m.state = dict(state["torch_state_dicts"]["unet"])
        return m


# ---------------------------------------------------------------------------
# Hybrid physics + ML
# ---------------------------------------------------------------------------
@dataclass
class FEAPrior:
    """A coarse finite-element prior: dz of a sparlab_form run of the commanded
    surface with `overrides` applied to the setup (e.g. {"element_size": 5e-3,
    "layers": 1}), through the run cache in `work_dir`. Needs the executable
    whenever the prior is evaluated (training and prediction)."""

    work_dir: str
    overrides: Dict[str, Any] = field(default_factory=dict)
    step: Any = -1

    #: what the prior's values are
    data_source = "SparLab simulation"

    def prior_deviation(self, commanded: HeightMap, setup: Any) -> np.ndarray:
        from ..fea.runner import simulate
        s = as_setup(setup).replace(**self.overrides)
        res = simulate(s, commanded, self.work_dir)
        formed = res.formed_surface(self.step, grid=commanded.grid)
        return formed.z - commanded.z

    def describe(self) -> Dict[str, Any]:
        return {"name": "FEAPrior", "overrides": dict(self.overrides), "step": self.step,
                "data_source": self.data_source}


class ResidualModel:
    """prior + learned residual: the learner fits y - prior on the table.

    `prior` is an object with ``prior_deviation(commanded, setup) -> (ny, nx)``
    dz [m] (a coarse FE run, `FEAPrior`, or a closed-form estimate); the
    surrogate evaluates it on each commanded surface and appends its value
    as the last feature column, `prior_dz`. The learner (a point model) sees
    that column as a feature too unless `prior_as_feature` is False. mean =
    prior + learner mean; std = learner std (the prior counts as exact).
    """

    kind = "point"

    def __init__(self, prior: Any, learner: DeviationModel, *, prior_as_feature: bool = True):
        if not hasattr(prior, "prior_deviation"):
            raise TypeError("the prior must implement prior_deviation(commanded, setup)")
        self.prior = prior
        self.learner = learner
        self.prior_as_feature = bool(prior_as_feature)
        self.feature_names: List[str] = []

    def _split(self, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        if self.feature_names and self.feature_names[-1] != PRIOR_FEATURE:
            raise ValueError(f"ResidualModel expects '{PRIOR_FEATURE}' as the last feature "
                             "column (build the table with prior=...)")
        p = X[:, -1]
        return (X if self.prior_as_feature else X[:, :-1]), p

    def fit(self, X: np.ndarray, y: np.ndarray, groups: np.ndarray, *,
            feature_names: Optional[Sequence[str]] = None) -> "ResidualModel":
        X, y, groups = _check_table(X, y, groups)
        self.feature_names = list(feature_names) if feature_names is not None else \
            [f"x{i}" for i in range(X.shape[1] - 1)] + [PRIOR_FEATURE]
        Xl, p = self._split(X)
        names = self.feature_names if self.prior_as_feature else self.feature_names[:-1]
        self.learner.fit(Xl, y - p, groups, feature_names=names)
        return self

    def predict(self, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        (X,) = _check_table(X)
        _check_names(self, X)
        Xl, p = self._split(X)
        mu, sd = self.learner.predict(Xl)
        return p + mu, sd

    def to_state(self) -> Dict[str, Any]:
        return {"class": "ResidualModel", "prior": self.prior,
                "learner": self.learner.to_state(), "prior_as_feature": self.prior_as_feature,
                "feature_names": list(self.feature_names)}

    @classmethod
    def from_state(cls, state: Mapping[str, Any]) -> "ResidualModel":
        m = cls(state["prior"], model_from_state(state["learner"]),
                prior_as_feature=state["prior_as_feature"])
        m.feature_names = list(state["feature_names"])
        return m


class TransferModel:
    """A frozen simulation-trained base plus a correction fitted on a few real parts.

    A handful of real parts supports only a low-dimensional correction. (An
    earlier default regressed the residual on the ~60 standardised features:
    from one to three real parts it made the base two to five times worse on
    proxy data, also when the "real" data equalled the simulation.)

    correction :
      "linear" (default) - y = a + (1 + b) mu, mu the base mean: a relative
          change of the springback magnitude, b ~ N(0, prior_gain^2), and an
          offset a ~ N(0, (prior_offset s)^2), s the RMS base prediction over
          the real parts;
      "gain" - y = (1 + b) mu, the gain alone;
      "gbm" - a small, heavily regularised gradient-boosting correction of
          y - mu on [features, base mean, base std], whose std is its
          in-sample residual std (optimistic).
    The linear corrections are Bayesian regressions with a FIXED prior in
    which every real part weighs one observation (its nodes are not
    independent: they share weight 1), the noise variance estimated from the
    part-weighted residuals; the predictive std combines the base std and
    the parameter uncertainty in quadrature.

    Validation: a correction is kept only if its leave-one-real-part-out
    error - fitted without a part, scored on it, for every part - is below
    the base's on the same parts. Otherwise, and always with fewer than two
    real parts (nothing to validate on), the model predicts exactly as the
    base; `accepted` and `validation` say which and why. `n_real` is the
    number of real parts (groups) seen; with none the model IS the base.

    On proxy data (a GBM base, ten held-out parts, 1-4 "real" parts in three
    draws each): with no shift, or 50 um noise, the validated "linear"
    correction was never accepted (error = the base's), while "gain" was
    accepted - wrongly - in 2 of 12 cases (+21-38 % error); with the rim and
    curvature 30-40 % stronger (base 0.18 mm) it cut the error to
    0.09-0.15 mm from two parts on; with one part nothing changes.
    """

    kind = "point"
    CORRECTIONS = ("gain", "linear", "gbm")

    def __init__(self, base: DeviationModel, correction: str = "linear", *, seed: int = 0,
                 prior_gain: float = 0.5, prior_offset: float = 0.5, min_parts: int = 2):
        if correction not in self.CORRECTIONS:
            raise ValueError(f"correction must be one of {self.CORRECTIONS}, got {correction!r} "
                             "(the former 'bayesian_ridge' on every feature extrapolated "
                             "wildly from a few parts and was removed)")
        if not (prior_gain > 0 and prior_offset > 0):
            raise ValueError("prior_gain and prior_offset must be > 0")
        if int(min_parts) < 2:
            raise ValueError("min_parts must be >= 2 (one part cannot be validated)")
        self.base = base
        self.correction = correction
        self.seed = int(seed)
        self.prior_gain = float(prior_gain)
        self.prior_offset = float(prior_offset)
        self.min_parts = int(min_parts)
        self.n_real = 0
        self.real_groups: List[str] = []
        self.accepted = False
        self.validation: Dict[str, Any] = {"reason": "no real parts"}
        self.theta: Optional[np.ndarray] = None
        self.cov: Optional[np.ndarray] = None
        self.estimator: Any = None
        self.resid_std = 0.0

    @property
    def feature_names(self) -> List[str]:
        return list(getattr(self.base, "feature_names", []))

    # -- the corrections ------------------------------------------------------
    def _design(self, mu: np.ndarray) -> np.ndarray:
        return mu[:, None] if self.correction == "gain" else np.column_stack(
            [np.ones_like(mu), mu])

    def _fit_linear(self, mu: np.ndarray, y: np.ndarray, groups: np.ndarray
                    ) -> Tuple[np.ndarray, np.ndarray]:
        uniq, inv = np.unique(groups, return_inverse=True)
        k = uniq.size
        w = 1.0 / np.bincount(inv)[inv]                      # each part weighs 1
        r = y - mu
        s = float(np.sqrt(np.sum(w * mu * mu) / k)) or 1.0
        Phi = self._design(mu)
        prec = np.array([1.0 / self.prior_gain ** 2]) if self.correction == "gain" else \
            np.array([1.0 / (self.prior_offset * s) ** 2, 1.0 / self.prior_gain ** 2])
        Lam = np.diag(prec)
        PtW = (Phi * w[:, None]).T
        s2 = max(float(np.sum(w * r * r)) / k, (1e-3 * s) ** 2)
        for _ in range(5):                                   # noise variance and theta
            A = PtW @ Phi / s2 + Lam
            theta = np.linalg.solve(A, PtW @ r / s2)
            s2 = max(float(np.sum(w * (r - Phi @ theta) ** 2)) / k, (1e-3 * s) ** 2)
        A = PtW @ Phi / s2 + Lam
        theta = np.linalg.solve(A, PtW @ r / s2)
        return theta, np.linalg.inv(A)

    def _fit_gbm(self, X: np.ndarray, mu: np.ndarray, sd: np.ndarray, y: np.ndarray):
        est = HistGradientBoostingRegressor(max_iter=60, learning_rate=0.05, max_depth=3,
                                            min_samples_leaf=200, l2_regularization=10.0,
                                            early_stopping=False, random_state=self.seed)
        Z = np.column_stack([X, mu, sd])
        est.fit(Z, y - mu)
        return est, float(np.std(y - mu - est.predict(Z)))

    def _fit_correction(self, X: np.ndarray, mu: np.ndarray, sd: np.ndarray, y: np.ndarray,
                        groups: np.ndarray) -> Any:
        if self.correction == "gbm":
            return self._fit_gbm(X, mu, sd, y)
        return self._fit_linear(mu, y, groups)

    def _apply(self, fitted: Any, X: np.ndarray, mu: np.ndarray, sd: np.ndarray
               ) -> Tuple[np.ndarray, np.ndarray]:
        if self.correction == "gbm":
            est, rs = fitted
            c = est.predict(np.column_stack([X, mu, sd]))
            return mu + c, np.sqrt(sd * sd + rs * rs)
        theta, cov = fitted
        Phi = self._design(mu)
        var = np.einsum("ij,jk,ik->i", Phi, cov, Phi)
        return mu + Phi @ theta, np.sqrt(sd * sd + np.maximum(var, 0.0))

    # -- protocol -------------------------------------------------------------
    def fit(self, X: np.ndarray, y: np.ndarray, groups: np.ndarray, *,
            feature_names: Optional[Sequence[str]] = None) -> "TransferModel":
        """`groups`: the real PART of every row (the unit of validation)."""
        y = np.asarray(y, dtype=float).ravel()
        self.theta = self.cov = self.estimator = None
        self.accepted = False
        if len(y) == 0:
            self.n_real, self.real_groups = 0, []
            self.validation = {"reason": "no real parts"}
            return self
        X, y, groups = _check_table(X, y, groups)
        _check_names(self.base, X)
        groups = np.asarray(groups).astype(str)
        self.real_groups = sorted(np.unique(groups).tolist())
        self.n_real = len(self.real_groups)
        mu, sd = self.base.predict(X)
        rms = lambda a, b: float(np.sqrt(np.mean((a - b) ** 2)))  # noqa: E731
        if self.n_real < self.min_parts:
            self.validation = {"reason": f"{self.n_real} real part(s): at least "
                               f"{self.min_parts} are needed to validate a correction, so "
                               "the model is the base", "n_parts": self.n_real}
            return self
        per = []
        for g in self.real_groups:
            out = groups == g
            fitted = self._fit_correction(X[~out], mu[~out], sd[~out], y[~out], groups[~out])
            c, _ = self._apply(fitted, X[out], mu[out], sd[out])
            per.append({"part": g, "base_rms_m": rms(mu[out], y[out]),
                        "corrected_rms_m": rms(c, y[out])})
        b = float(np.mean([p["base_rms_m"] for p in per]))
        c = float(np.mean([p["corrected_rms_m"] for p in per]))
        self.accepted = c < b
        self.validation = {"reason": "leave-one-real-part-out error "
                           + ("below" if self.accepted else "not below") + " the base's",
                           "n_parts": self.n_real, "lopo_base_rms_m": b,
                           "lopo_corrected_rms_m": c, "per_part": per}
        if self.accepted:
            fitted = self._fit_correction(X, mu, sd, y, groups)
            if self.correction == "gbm":
                self.estimator, self.resid_std = fitted
            else:
                self.theta, self.cov = fitted
                self.validation["theta"] = [float(t) for t in self.theta]
        return self

    def predict(self, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        if not self.accepted:
            return self.base.predict(X)
        (X,) = _check_table(X)
        mu, sd = self.base.predict(X)
        fitted = (self.estimator, self.resid_std) if self.correction == "gbm" else \
            (self.theta, self.cov)
        return self._apply(fitted, X, mu, sd)

    def to_state(self) -> Dict[str, Any]:
        return {"class": "TransferModel", "base": self.base.to_state(),
                "correction": self.correction, "seed": self.seed,
                "prior_gain": self.prior_gain, "prior_offset": self.prior_offset,
                "min_parts": self.min_parts, "n_real": self.n_real,
                "real_groups": list(self.real_groups), "accepted": self.accepted,
                "validation": self.validation, "theta": self.theta, "cov": self.cov,
                "estimator": self.estimator, "resid_std": self.resid_std}

    @classmethod
    def from_state(cls, state: Mapping[str, Any]) -> "TransferModel":
        m = cls(model_from_state(state["base"]), state["correction"], seed=state["seed"],
                prior_gain=state["prior_gain"], prior_offset=state["prior_offset"],
                min_parts=state["min_parts"])
        m.n_real = int(state["n_real"])
        m.real_groups = list(state["real_groups"])
        m.accepted = bool(state["accepted"])
        m.validation = dict(state["validation"])
        m.theta = None if state["theta"] is None else np.asarray(state["theta"])
        m.cov = None if state["cov"] is None else np.asarray(state["cov"])
        m.estimator = state["estimator"]
        m.resid_std = float(state["resid_std"])
        return m


MODEL_CLASSES = {"GBMEnsemble": GBMEnsemble, "MLPEnsemble": MLPEnsemble,
                 "FieldUNet": FieldUNet, "ResidualModel": ResidualModel,
                 "TransferModel": TransferModel}


def model_from_state(state: Mapping[str, Any]):
    """Rebuild any model of this module from its `to_state()`."""
    cls = MODEL_CLASSES.get(state.get("class"))
    if cls is None:
        raise PrecompError(f"unknown model class {state.get('class')!r}; known: "
                           f"{sorted(MODEL_CLASSES)}")
    return cls.from_state(state)


__all__ = ["DeviationModel", "FieldDeviationModel", "GBMEnsemble", "MLPEnsemble", "FieldUNet",
           "ResidualModel", "TransferModel", "FEAPrior", "model_from_state", "thread_limit",
           "MODEL_CLASSES", "THREADS_ENV"]
