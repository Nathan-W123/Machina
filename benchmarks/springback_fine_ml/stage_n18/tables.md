# stage_n18: SparLab simulations (not experiments)

Backing plate. 48 train / 16 calibration parts, 128 samples at penalty 3 (springback_fine_bulk); 8 test parts verified at penalty 10 (springback_fine).

| method | FE runs / part | parts | vertical RMS mean | median | max mean | rim sag | interior RMS (bias) | RMS / uncomp. | better than uncomp. | better than FE-DA-1 | minus FE-DA-1 (sd) |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| uncompensated | 1 | 8 | **0.300** | 0.298 | 0.583 | -0.436 | 0.149 (-0.05) | 1.00 | 0 | 1 | 0.037 (0.043) |
| FE-DA-1 | 2 | 8 | **0.263** | 0.252 | 0.535 | -0.334 | 0.162 (0.08) | 0.87 | 7 | 0 | 0.000 (0.000) |
| FE-DA-2 | 3 | 8 | **0.272** | 0.259 | 0.595 | -0.318 | 0.188 (0.12) | 0.90 | 6 | 3 | 0.009 (0.020) |
| FE-DA-1 fixed | 2 | 8 | **0.264** | 0.256 | 0.538 | -0.337 | 0.163 (0.09) | 0.88 | 7 | 1 | 0.001 (0.004) |
| ML-GBM | 1 | 8 | **0.253** | 0.242 | 0.521 | -0.348 | 0.126 (0.06) | 0.84 | 8 | 5 | -0.010 (0.018) |
| ML-MLP | 1 | 8 | **0.250** | 0.246 | 0.519 | -0.350 | 0.116 (0.04) | 0.83 | 8 | 6 | -0.012 (0.017) |

Per part, vertical RMS [mm] (max |dev| in brackets):

| part | uncompensated | FE-DA-1 | FE-DA-2 | FE-DA-1 fixed | ML-GBM | ML-MLP |
|---|--:|--:|--:|--:|--:|--:|
| dome-s2026-0000 | 0.307 (0.50) | 0.213 (0.41) | 0.203 (0.40) | 0.213 (0.41) | 0.214 (0.41) | 0.217 (0.41) |
| dome-s2026-0001 | 0.290 (0.59) | 0.223 (0.51) | 0.216 (0.50) | 0.223 (0.51) | 0.222 (0.52) | 0.215 (0.52) |
| elliptic_cone-s2026-0000 | 0.277 (0.57) | 0.252 (0.50) | 0.294 (0.69) | 0.261 (0.53) | 0.263 (0.53) * | 0.261 (0.53) * |
| elliptic_cone-s2026-0001 | 0.305 (0.61) | 0.266 (0.56) | 0.260 (0.55) | 0.270 (0.55) | 0.262 (0.54) | 0.265 (0.56) |
| pyramid-s2026-0000 | 0.290 (0.57) | 0.213 (0.48) | 0.217 (0.48) | 0.210 (0.48) | 0.197 (0.48) | 0.202 (0.49) |
| pyramid-s2026-0001 | 0.308 (0.69) | 0.299 (0.61) | 0.304 (0.63) | 0.299 (0.61) | 0.301 (0.62) | 0.269 (0.56) |
| truncated_cone-s2026-0000 | 0.280 (0.55) | 0.251 (0.46) | 0.257 (0.61) | 0.251 (0.46) | 0.221 (0.47) | 0.231 (0.49) |
| truncated_cone-s2026-0001 | 0.344 (0.60) | 0.386 (0.74) | 0.426 (0.89) | 0.386 (0.74) | 0.344 (0.60) = | 0.344 (0.60) = |

`*` target outside the training envelope (override); `=` the surrogate kept the target (no compensation).

dz prediction error on the test samples (penalty-10 runs, uncompensated and fixed FE-DA-1):

| model | samples | dz error RMS [mm] | dz RMS [mm] | relative | 90 % coverage | flagged by the envelope |
|---|--:|--:|--:|--:|--:|--:|
| gbm | 16 | 0.049 | 0.289 | 0.17 | 0.95 | 0.06 |
| mlp | 16 | 0.052 | 0.289 | 0.18 | 0.88 | 0.06 |

ML-GBM compensation per test part:

| part | stopped | predictions | best iterate | predicted RMS | verified RMS | commanded depth | reach raise | envelope (target / shot) |
|---|---|--:|--:|--:|--:|--:|--:|---|
| truncated_cone-s2026-0000 | stagnation | 3 | 2 | 0.273 | 0.221 | 3.08 | 0.35 | True / True |
| truncated_cone-s2026-0001 | stagnation | 2 | 0 | 0.341 | 0.344 | 2.63 | 0.00 | True / True |
| pyramid-s2026-0000 | stagnation | 4 | 3 | 0.218 | 0.197 | 3.37 | 0.34 | True / True |
| pyramid-s2026-0001 | stagnation | 3 | 1 | 0.277 | 0.301 | 2.72 | 0.66 | True / True |
| dome-s2026-0000 | left_envelope | 2 | 1 | 0.242 | 0.214 | 2.20 | 0.04 | True / True |
| dome-s2026-0001 | stagnation | 3 | 1 | 0.210 | 0.222 | 3.76 | 0.02 | True / True |
| elliptic_cone-s2026-0000 | stagnation | 3 | 1 | 0.240 | 0.263 | 2.24 | 0.51 | False / True (override) |
| elliptic_cone-s2026-0001 | left_envelope | 2 | 1 | 0.345 | 0.262 | 3.24 | 0.54 | True / True |

ML-MLP compensation per test part:

| part | stopped | predictions | best iterate | predicted RMS | verified RMS | commanded depth | reach raise | envelope (target / shot) |
|---|---|--:|--:|--:|--:|--:|--:|---|
| truncated_cone-s2026-0000 | stagnation | 5 | 4 | 0.209 | 0.231 | 3.12 | 0.30 | True / True |
| truncated_cone-s2026-0001 | stagnation | 2 | 0 | 0.330 | 0.344 | 2.63 | 0.00 | True / True |
| pyramid-s2026-0000 | stagnation | 5 | 3 | 0.167 | 0.202 | 3.44 | 0.31 | True / True |
| pyramid-s2026-0001 | stagnation | 4 | 3 | 0.271 | 0.269 | 2.61 | 0.59 | True / True |
| dome-s2026-0000 | left_envelope | 3 | 2 | 0.184 | 0.217 | 2.24 | 0.03 | True / True |
| dome-s2026-0001 | stagnation | 4 | 3 | 0.211 | 0.215 | 3.70 | 0.05 | True / True |
| elliptic_cone-s2026-0000 | stagnation | 3 | 1 | 0.236 | 0.261 | 2.31 | 0.48 | False / True (override) |
| elliptic_cone-s2026-0001 | left_envelope | 2 | 1 | 0.339 | 0.265 | 3.29 | 0.58 | True / True |

Failures:

none
