# Springback benchmark test parts on the fine mesh: SparLab simulations, not experiments

Setup `springback_fine` (0.833 mm incompatible-mode Hex8, 1 layer, 5 thickness points, penalty 10); the 8 test parts of benchmarks/springback. Deviation of the released part from the target over the part [mm]; rim sag = mean vertical deviation of the part less than 1 mm deep (negative = too deep); interior = the part deeper than that.

70 simulations complete, 1 failed records.

| strategy | method | FE runs / part | parts | vertical RMS mean | median | max mean | normal RMS | rim sag | upper band RMS | interior RMS (bias) | flange bias | runtime / run [s] | CPU h / part | Newton its | failed |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| none | uncompensated | 1 | 8 | **0.667** | 0.660 | 1.218 | 0.662 | -0.994 | 1.007 | 0.291 (-0.15) | -0.223 | 1652 | 0.46 | 3586 | 0 |
| none | FE-DA-1 | 2 | 8 | **0.569** | 0.565 | 1.044 | 0.573 | -0.793 | 0.816 | 0.295 (0.11) | -0.202 | 1273 | 0.81 | 2710 | 0 |
| none | FE-DA-2 | 3 | 8 | **0.556** | 0.558 | 1.079 | 0.560 | -0.809 | 0.833 | 0.252 (0.01) | -0.200 | 1288 | 1.17 | 2677 | 0 |
| backing_plate | uncompensated | 1 | 8 | **0.307** | 0.301 | 0.628 | 0.287 | -0.436 | 0.445 | 0.166 (-0.03) | -0.046 | 3526 | 0.98 | 5927 | 0 |
| backing_plate | FE-DA-1 | 2 | 8 | **0.305** | 0.253 | 0.566 | 0.295 | -0.321 | 0.343 | 0.244 (0.17) | -0.044 | 3000 | 1.81 | 4974 | 0 |
| backing_plate | FE-DA-2 | 3 | 8 | **0.272** | 0.256 | 0.537 | 0.260 | -0.332 | 0.356 | 0.191 (0.04) | -0.045 | 3276 | 2.72 | 5390 | 0 |
| dsif | uncompensated | 1 | 8 | **0.295** | 0.294 | 0.666 | 0.276 | -0.148 | 0.284 | 0.279 (0.20) | 0.039 | 3074 | 0.85 | 6567 | 0 |
| dsif | FE-DA-1 | 2 | 7 | **0.355** | 0.290 | 0.696 | 0.346 | -0.222 | 0.297 | 0.341 (0.25) | 0.049 | 2248 | 1.42 | 4682 | 1 |
| dsif | FE-DA-2 | 3 | 7 | **0.274** | 0.305 | 0.629 | 0.262 | -0.233 | 0.307 | 0.238 (0.07) | 0.063 | 2737 | 2.18 | 5670 | 0 |

The same over the 7 parts where every strategy and method has a result (dome-s2026-0000, dome-s2026-0001, elliptic_cone-s2026-0000, elliptic_cone-s2026-0001, pyramid-s2026-0000, truncated_cone-s2026-0000, truncated_cone-s2026-0001):

| strategy | method | vertical RMS mean | max mean | rim sag | interior RMS (bias) |
|---|---|--:|--:|--:|--:|
| none | uncompensated | **0.670** | 1.205 | -0.992 | 0.298 (-0.16) |
| none | FE-DA-1 | **0.568** | 1.047 | -0.805 | 0.272 (0.07) |
| none | FE-DA-2 | **0.550** | 1.056 | -0.801 | 0.242 (0.04) |
| backing_plate | uncompensated | **0.307** | 0.619 | -0.432 | 0.171 (-0.03) |
| backing_plate | FE-DA-1 | **0.288** | 0.544 | -0.323 | 0.212 (0.13) |
| backing_plate | FE-DA-2 | **0.261** | 0.519 | -0.322 | 0.176 (0.07) |
| dsif | uncompensated | **0.298** | 0.672 | -0.170 | 0.277 (0.19) |
| dsif | FE-DA-1 | **0.355** | 0.696 | -0.222 | 0.341 (0.25) |
| dsif | FE-DA-2 | **0.274** | 0.629 | -0.233 | 0.238 (0.07) |

The 2 mm benchmark (benchmarks/springback, stage n18, single-point forming) on the same parts:

| method | vertical RMS mean | max mean | rim sag | upper band RMS | interior RMS (bias) |
|---|--:|--:|--:|--:|--:|
| uncompensated | 0.731 | 1.341 | -1.086 | 1.099 | 0.328 (-0.14) |
| FE-DA-1 | 0.623 | 1.129 | -0.856 | 0.879 | 0.346 (0.13) |
| FE-DA-2 | 0.602 | 1.156 | -0.884 | 0.905 | 0.272 (0.02) |

2 mm runs: 243 s mean, 227 s median per run (stage n18, two at a time on a shared machine).

Per part, vertical RMS [mm] (max |dev| in brackets; rim sag after the slash):

| part | 2 mm uncomp. | backing_plate uncompensated | backing_plate FE-DA-1 | backing_plate FE-DA-2 | dsif uncompensated | dsif FE-DA-1 | dsif FE-DA-2 | none uncompensated | none FE-DA-1 | none FE-DA-2 |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| dome-s2026-0000 | 0.654 (1.07) / -0.83 | 0.298 (0.49) / -0.37 | 0.211 (0.41) / -0.26 | 0.202 (0.39) / -0.24 | 0.341 (0.61) / -0.42 | 0.252 (0.57) / -0.29 | 0.183 (0.39) / -0.09 | 0.638 (0.98) / -0.81 | 0.487 (0.83) / -0.63 | 0.465 (0.81) / -0.60 |
| dome-s2026-0001 | 0.776 (1.40) / -1.16 | 0.358 (0.95) / -0.46 | 0.262 (0.52) / -0.36 | 0.226 (0.50) / -0.34 | 0.336 (1.08) / -0.29 | 0.290 (0.67) / -0.42 | 0.310 (0.85) / -0.45 | 0.662 (1.22) / -0.99 | 0.542 (1.12) / -0.85 | 0.498 (1.08) / -0.81 |
| elliptic_cone-s2026-0000 | 0.621 (1.21) / -0.96 | 0.277 (0.57) / -0.42 | 0.239 (0.50) / -0.32 | 0.266 (0.61) / -0.27 | 0.247 (0.54) / -0.27 | 0.385 (0.67) / +0.08 | 0.363 (0.78) / -0.34 | 0.602 (1.17) / -0.93 | 0.481 (1.00) / -0.71 | 0.465 (0.95) / -0.64 |
| elliptic_cone-s2026-0001 | 0.841 (1.51) / -1.24 | 0.305 (0.61) / -0.44 | 0.285 (0.55) / -0.37 | 0.279 (0.55) / -0.36 | 0.319 (0.83) / +0.17 | 0.359 (0.79) / -0.38 | 0.328 (0.69) / -0.13 | 0.751 (1.36) / -1.12 | 0.666 (1.23) / -0.98 | 0.657 (1.21) / -0.95 |
| pyramid-s2026-0000 | 0.791 (1.45) / -1.14 | 0.290 (0.57) / -0.41 | 0.245 (0.48) / -0.28 | 0.236 (0.47) / -0.27 | 0.221 (0.49) / +0.04 | 0.268 (0.69) / -0.12 | 0.187 (0.56) / -0.13 | 0.684 (1.27) / -0.99 | 0.553 (1.10) / -0.81 | 0.544 (1.08) / -0.79 |
| pyramid-s2026-0001 | 0.699 (1.38) / -1.08 | 0.308 (0.69) / -0.46 | 0.426 (0.73) / -0.31 | 0.349 (0.67) / -0.40 | 0.269 (0.62) / +0.01 | FAILED | - | 0.648 (1.31) / -1.00 | 0.577 (1.03) / -0.71 | 0.599 (1.24) / -0.87 |
| truncated_cone-s2026-0000 | 0.766 (1.43) / -1.15 | 0.280 (0.55) / -0.42 | 0.226 (0.48) / -0.32 | 0.245 (0.53) / -0.31 | 0.259 (0.50) / -0.00 | 0.250 (0.44) / -0.10 | 0.244 (0.52) / -0.17 | 0.698 (1.25) / -1.05 | 0.586 (1.13) / -0.89 | 0.573 (1.11) / -0.85 |
| truncated_cone-s2026-0001 | 0.703 (1.29) / -1.12 | 0.344 (0.60) / -0.51 | 0.550 (0.86) / -0.34 | 0.371 (0.59) / -0.47 | 0.364 (0.65) / -0.41 | 0.679 (1.03) / -0.31 | 0.305 (0.62) / -0.32 | 0.657 (1.19) / -1.05 | 0.659 (0.93) / -0.76 | 0.648 (1.14) / -0.97 |

Single-point forming, fine minus 2 mm, mean over the parts [mm]:

| method | parts | vertical RMS | max | rim sag | interior RMS |
|---|--:|--:|--:|--:|--:|
| uncompensated | 8 | -0.064 | -0.123 | +0.092 | -0.037 |
| FE-DA-1 | 8 | -0.054 | -0.085 | +0.064 | -0.052 |
| FE-DA-2 | 8 | -0.046 | -0.077 | +0.075 | -0.019 |

Failures:

* pyramid-s2026-0001 dsif FE-DA-1: FormingError: this deck failed before (exit code -15 (failed)); pass retry_failed=True to run it again. Entry: /home/user/wt/fine/benchmarks/springback_fine/work/runs/runs/52/52457c00c71d5f677af99695a1584d36fbb76ae888a89a5cc1d1b25e0b645111
