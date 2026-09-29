"""precomp.ml - learned springback prediction and surrogate-driven pre-compensation.

The target is the vertical deviation dz(x, y) = z_formed - z_commanded [m] of
the tool-side surface after springback, as a function of the COMMANDED surface,
the process and the material; the models are queried at compensated shapes
during displacement adjustment, so the data sets contain such shapes too.

Modules
-------
features     the versioned per-node feature schema (FEATURE_SCHEMA_VERSION,
             FeatureSpec, FeatureConfig, point_features, feature_maps)
dataset      Sample, the on-disk Dataset, feature tables, grouped splits
generate     design of experiments, commanded-shape variants, SparlabSimulator,
             ProxySimulator (analytic, NOT physics), data generation
models       GBMEnsemble, MLPEnsemble, FieldUNet, ResidualModel, TransferModel
uncertainty  grouped split-conformal calibration, intervals, coverage reports
ood          the training envelope (Mahalanobis + kNN, point and part level)
surrogate    DeviationSurrogate (model + features + calibration + envelope),
             train_surrogate, transfer_surrogate
registry     save_model / load_model: bundles with manifests
evaluate     held-out, cross-validated and leave-family-out evaluation
compensate   surrogate displacement adjustment, verify_with_fea
active       ranking of candidate parts for the next simulations
cli          `precomp dataset|train|evaluate|active`

Every metric states its data source: "SparLab simulation", "proxy - not
physics" or "scan". torch is needed only for MLPEnsemble and FieldUNet.
"""

from .dataset import (SOURCE_LABELS, Dataset, Sample, Table, build_table, family_split,
                      grouped_kfold, grouped_split, source_label)
from .features import (FEATURE_SCHEMA_VERSION, FeatureConfig, FeatureSpec, feature_maps,
                       global_features, point_features, region_labels)
from .generate import (DesignPoint, DesignSpace, PerturbationSpec, ProxyParams, ProxySimulator,
                       SparlabSimulator, design_points, generate, simulate_samples)

__all__ = [
    "FEATURE_SCHEMA_VERSION", "FeatureConfig", "FeatureSpec", "point_features", "feature_maps",
    "global_features", "region_labels",
    "Sample", "Dataset", "Table", "build_table", "grouped_split", "grouped_kfold",
    "family_split", "source_label", "SOURCE_LABELS",
    "DesignSpace", "DesignPoint", "design_points", "PerturbationSpec", "ProxyParams",
    "ProxySimulator", "SparlabSimulator", "generate", "simulate_samples",
]
