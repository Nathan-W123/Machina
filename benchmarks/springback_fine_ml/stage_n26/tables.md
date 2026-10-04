# stage_n26: SparLab simulations (not experiments)

Backing plate. 72 train / 24 calibration parts, 192 samples at penalty 3 (springback_fine_bulk); 8 test parts verified at penalty 10 (springback_fine).

| method | FE runs / part | parts | vertical RMS mean | median | max mean | rim sag | interior RMS (bias) | RMS / uncomp. | better than uncomp. | better than FE-DA-1 | minus FE-DA-1 (sd) |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| uncompensated | 1 | 8 | **0.300** | 0.298 | 0.583 | -0.436 | 0.149 (-0.05) | 1.00 | 0 | 1 | 0.037 (0.043) |
| FE-DA-1 | 2 | 8 | **0.263** | 0.252 | 0.535 | -0.334 | 0.162 (0.08) | 0.87 | 7 | 0 | 0.000 (0.000) |
| FE-DA-2 | 3 | 8 | **0.272** | 0.259 | 0.595 | -0.318 | 0.188 (0.12) | 0.90 | 6 | 3 | 0.009 (0.020) |
| FE-DA-1 fixed | 2 | 8 | **0.264** | 0.256 | 0.538 | -0.337 | 0.163 (0.09) | 0.88 | 7 | 1 | 0.001 (0.004) |
| ML-GBM | 1 | 8 | **0.256** | 0.249 | 0.523 | -0.349 | 0.130 (0.06) | 0.85 | 8 | 4 | -0.007 (0.016) |
| ML-MLP | 1 | 8 | **0.248** | 0.246 | 0.517 | -0.346 | 0.115 (0.04) | 0.82 | 8 | 7 | -0.015 (0.016) |

Per part, vertical RMS [mm] (max |dev| in brackets):

| part | uncompensated | FE-DA-1 | FE-DA-2 | FE-DA-1 fixed | ML-GBM | ML-MLP |
|---|--:|--:|--:|--:|--:|--:|
| dome-s2026-0000 | 0.307 (0.50) | 0.213 (0.41) | 0.203 (0.40) | 0.213 (0.41) | 0.215 (0.42) | 0.211 (0.41) |
| dome-s2026-0001 | 0.290 (0.59) | 0.223 (0.51) | 0.216 (0.50) | 0.223 (0.51) | 0.224 (0.52) | 0.217 (0.52) |
| elliptic_cone-s2026-0000 | 0.277 (0.57) | 0.252 (0.50) | 0.294 (0.69) | 0.261 (0.53) | 0.261 (0.53) | 0.258 (0.53) |
| elliptic_cone-s2026-0001 | 0.305 (0.61) | 0.266 (0.56) | 0.260 (0.55) | 0.270 (0.55) | 0.259 (0.55) | 0.258 (0.55) |
| pyramid-s2026-0000 | 0.290 (0.57) | 0.213 (0.48) | 0.217 (0.48) | 0.210 (0.48) | 0.201 (0.48) | 0.196 (0.48) |
| pyramid-s2026-0001 | 0.308 (0.69) | 0.299 (0.61) | 0.304 (0.63) | 0.299 (0.61) | 0.302 (0.62) | 0.268 (0.56) |
| truncated_cone-s2026-0000 | 0.280 (0.55) | 0.251 (0.46) | 0.257 (0.61) | 0.251 (0.46) | 0.239 (0.48) | 0.234 (0.50) |
| truncated_cone-s2026-0001 | 0.344 (0.60) | 0.386 (0.74) | 0.426 (0.89) | 0.386 (0.74) | 0.344 (0.60) = | 0.344 (0.60) = |

`*` target outside the training envelope (override); `=` the surrogate kept the target (no compensation).

dz prediction error on the test samples (penalty-10 runs, uncompensated and fixed FE-DA-1):

| model | samples | dz error RMS [mm] | dz RMS [mm] | relative | 90 % coverage | flagged by the envelope |
|---|--:|--:|--:|--:|--:|--:|
| gbm | 16 | 0.046 | 0.289 | 0.16 | 0.94 | 0.00 |
| mlp | 16 | 0.047 | 0.289 | 0.16 | 0.88 | 0.00 |

ML-GBM compensation per test part:

| part | stopped | predictions | best iterate | predicted RMS | verified RMS | commanded depth | reach raise | envelope (target / shot) |
|---|---|--:|--:|--:|--:|--:|--:|---|
| truncated_cone-s2026-0000 | stagnation | 3 | 1 | 0.272 | 0.239 | 3.10 | 0.42 | True / True |
| truncated_cone-s2026-0001 | stagnation | 2 | 0 | 0.334 | 0.344 | 2.63 | 0.00 | True / True |
| pyramid-s2026-0000 | stagnation | 4 | 3 | 0.209 | 0.201 | 3.39 | 0.35 | True / True |
| pyramid-s2026-0001 | stagnation | 3 | 1 | 0.284 | 0.302 | 2.72 | 0.66 | True / True |
| dome-s2026-0000 | stagnation | 4 | 3 | 0.196 | 0.215 | 2.18 | 0.06 | True / True |
| dome-s2026-0001 | stagnation | 3 | 1 | 0.214 | 0.224 | 3.76 | 0.02 | True / True |
| elliptic_cone-s2026-0000 | left_envelope | 2 | 1 | 0.235 | 0.261 | 2.24 | 0.49 | True / True |
| elliptic_cone-s2026-0001 | left_envelope | 3 | 2 | 0.323 | 0.259 | 3.28 | 0.42 | True / True |

ML-MLP compensation per test part:

| part | stopped | predictions | best iterate | predicted RMS | verified RMS | commanded depth | reach raise | envelope (target / shot) |
|---|---|--:|--:|--:|--:|--:|--:|---|
| truncated_cone-s2026-0000 | stagnation | 6 | 5 | 0.210 | 0.234 | 3.12 | 0.31 | True / True |
| truncated_cone-s2026-0001 | stagnation | 2 | 0 | 0.320 | 0.344 | 2.63 | 0.00 | True / True |
| pyramid-s2026-0000 | stagnation | 4 | 3 | 0.206 | 0.196 | 3.43 | 0.34 | True / True |
| pyramid-s2026-0001 | stagnation | 5 | 3 | 0.268 | 0.268 | 2.60 | 0.58 | True / True |
| dome-s2026-0000 | stagnation | 6 | 4 | 0.156 | 0.211 | 2.19 | 0.08 | True / True |
| dome-s2026-0001 | stagnation | 7 | 6 | 0.199 | 0.217 | 3.71 | 0.13 | True / True |
| elliptic_cone-s2026-0000 | stagnation | 3 | 1 | 0.247 | 0.258 | 2.27 | 0.49 | True / True |
| elliptic_cone-s2026-0001 | left_envelope | 3 | 2 | 0.308 | 0.258 | 3.33 | 0.47 | True / True |

Failures:

none
