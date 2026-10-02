# Springback benchmark test parts on the fine mesh: SparLab simulations, not experiments

Setup `springback_fine` (0.833 mm incompatible-mode Hex8, 1 layer, 5 thickness points, penalty 10); the 8 test parts of benchmarks/springback. Deviation of the released part from the target over the part [mm]; rim sag = mean vertical deviation of the part less than 1 mm deep (negative = too deep); interior = the part deeper than that.

69 simulations complete, 3 failed records.

| strategy | method | FE runs / part | parts | vertical RMS mean | median | max mean | normal RMS | rim sag | upper band RMS | interior RMS (bias) | flange bias | runtime / run [s] | CPU h / part | Newton its | failed |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| none | uncompensated | 1 | 8 | **0.665** | 0.652 | 1.218 | 0.659 | -0.994 | 1.008 | 0.284 (-0.18) | -0.223 | 1657 | 0.46 | 3626 | 0 |
| none | FE-DA-1 | 2 | 8 | **0.555** | 0.555 | 1.093 | 0.559 | -0.819 | 0.842 | 0.232 (0.02) | -0.208 | 1345 | 0.83 | 2868 | 0 |
| none | FE-DA-2 | 3 | 5 | **0.535** | 0.489 | 1.036 | 0.540 | -0.760 | 0.788 | 0.232 (0.08) | -0.195 | 1282 | 1.12 | 2548 | 3 |
| backing_plate | uncompensated | 1 | 8 | **0.300** | 0.298 | 0.583 | 0.279 | -0.436 | 0.445 | 0.149 (-0.05) | -0.047 | 3568 | 0.99 | 5959 | 0 |
| backing_plate | FE-DA-1 | 2 | 8 | **0.263** | 0.252 | 0.535 | 0.250 | -0.334 | 0.356 | 0.162 (0.08) | -0.048 | 3352 | 1.92 | 5608 | 0 |
| backing_plate | FE-DA-2 | 3 | 8 | **0.272** | 0.259 | 0.595 | 0.260 | -0.318 | 0.344 | 0.188 (0.12) | -0.056 | 3336 | 2.85 | 5366 | 0 |
| dsif | uncompensated | 1 | 8 | **0.273** | 0.250 | 0.576 | 0.258 | -0.215 | 0.275 | 0.240 (0.17) | 0.029 | 2576 | 0.72 | 5867 | 0 |
| dsif | FE-DA-1 | 2 | 8 | **0.265** | 0.263 | 0.624 | 0.254 | -0.292 | 0.343 | 0.169 (0.08) | 0.017 | 2706 | 1.47 | 6209 | 0 |
| dsif | FE-DA-2 | 3 | 8 | **0.277** | 0.278 | 0.679 | 0.268 | -0.326 | 0.372 | 0.164 (0.09) | 0.006 | 2838 | 2.26 | 6580 | 0 |

The same over the 5 parts where every strategy and method has a result (dome-s2026-0000, dome-s2026-0001, elliptic_cone-s2026-0000, elliptic_cone-s2026-0001, truncated_cone-s2026-0000):

| strategy | method | vertical RMS mean | max mean | rim sag | interior RMS (bias) |
|---|---|--:|--:|--:|--:|
| none | uncompensated | **0.666** | 1.196 | -0.980 | 0.293 (-0.21) |
| none | FE-DA-1 | **0.549** | 1.084 | -0.810 | 0.205 (-0.00) |
| none | FE-DA-2 | **0.535** | 1.036 | -0.760 | 0.232 (0.08) |
| backing_plate | uncompensated | **0.292** | 0.562 | -0.422 | 0.138 (-0.08) |
| backing_plate | FE-DA-1 | **0.241** | 0.489 | -0.319 | 0.120 (0.06) |
| backing_plate | FE-DA-2 | **0.246** | 0.552 | -0.306 | 0.136 (0.08) |
| dsif | uncompensated | **0.276** | 0.598 | -0.225 | 0.221 (0.12) |
| dsif | FE-DA-1 | **0.253** | 0.608 | -0.301 | 0.121 (0.03) |
| dsif | FE-DA-2 | **0.265** | 0.670 | -0.346 | 0.101 (0.04) |

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
| dome-s2026-0000 | 0.654 (1.07) / -0.83 | 0.307 (0.50) / -0.37 | 0.213 (0.41) / -0.26 | 0.203 (0.40) / -0.24 | 0.347 (0.62) / -0.43 | 0.281 (0.62) / -0.34 | 0.298 (0.60) / -0.37 | 0.641 (0.98) / -0.81 | 0.477 (0.82) / -0.61 | 0.468 (0.81) / -0.60 |
| dome-s2026-0001 | 0.776 (1.40) / -1.16 | 0.290 (0.59) / -0.46 | 0.223 (0.51) / -0.36 | 0.216 (0.50) / -0.35 | 0.225 (0.63) / -0.33 | 0.310 (0.77) / -0.49 | 0.258 (0.76) / -0.36 | 0.638 (1.22) / -1.00 | 0.543 (1.24) / -0.87 | 0.489 (1.08) / -0.80 |
| elliptic_cone-s2026-0000 | 0.621 (1.21) / -0.96 | 0.277 (0.57) / -0.42 | 0.252 (0.50) / -0.30 | 0.294 (0.69) / -0.26 | 0.254 (0.58) / -0.31 | 0.289 (0.66) / -0.39 | 0.327 (0.73) / -0.43 | 0.602 (1.17) / -0.93 | 0.490 (1.02) / -0.73 | 0.483 (0.93) / -0.58 |
| elliptic_cone-s2026-0001 | 0.841 (1.51) / -1.24 | 0.305 (0.61) / -0.44 | 0.266 (0.56) / -0.36 | 0.260 (0.55) / -0.37 | 0.307 (0.67) / +0.00 | 0.200 (0.54) / -0.15 | 0.228 (0.62) / -0.28 | 0.751 (1.36) / -1.12 | 0.664 (1.25) / -0.98 | 0.659 (1.23) / -0.97 |
| pyramid-s2026-0000 | 0.791 (1.45) / -1.14 | 0.290 (0.57) / -0.41 | 0.213 (0.48) / -0.29 | 0.217 (0.48) / -0.27 | 0.195 (0.48) / -0.04 | 0.190 (0.50) / -0.24 | 0.148 (0.50) / -0.12 | 0.684 (1.27) / -0.99 | 0.566 (1.15) / -0.85 | FAILED |
| pyramid-s2026-0001 | 0.699 (1.38) / -1.08 | 0.308 (0.69) / -0.46 | 0.299 (0.61) / -0.37 | 0.304 (0.63) / -0.37 | 0.246 (0.50) / -0.16 | 0.245 (0.61) / -0.17 | 0.309 (0.78) / -0.22 | 0.648 (1.31) / -1.00 | 0.538 (1.12) / -0.79 | FAILED |
| truncated_cone-s2026-0000 | 0.766 (1.43) / -1.15 | 0.280 (0.55) / -0.42 | 0.251 (0.46) / -0.31 | 0.257 (0.61) / -0.30 | 0.246 (0.49) / -0.05 | 0.187 (0.46) / -0.13 | 0.212 (0.64) / -0.29 | 0.698 (1.25) / -1.05 | 0.572 (1.09) / -0.85 | 0.576 (1.12) / -0.85 |
| truncated_cone-s2026-0001 | 0.703 (1.29) / -1.12 | 0.344 (0.60) / -0.51 | 0.386 (0.74) / -0.42 | 0.426 (0.89) / -0.38 | 0.364 (0.65) / -0.40 | 0.415 (0.84) / -0.41 | 0.438 (0.80) / -0.53 | 0.657 (1.19) / -1.05 | 0.592 (1.05) / -0.86 | FAILED |

Single-point forming, fine minus 2 mm, mean over the parts [mm]:

| method | parts | vertical RMS | max | rim sag | interior RMS |
|---|--:|--:|--:|--:|--:|
| uncompensated | 8 | -0.066 | -0.123 | +0.092 | -0.044 |
| FE-DA-1 | 8 | -0.068 | -0.036 | +0.038 | -0.115 |
| FE-DA-2 | 5 | -0.041 | -0.065 | +0.077 | -0.014 |

Failures:

* pyramid-s2026-0000 none FE-DA-2: PrecompError: the tool path reaches 0.024 m from the centre (tool radius included), beyond the unclamped half-width 0.015 m; enlarge the blank or reduce the clamp margin
* pyramid-s2026-0001 none FE-DA-2: PrecompError: the tool path reaches 0.024 m from the centre (tool radius included), beyond the unclamped half-width 0.015 m; enlarge the blank or reduce the clamp margin
* truncated_cone-s2026-0001 none FE-DA-2: PrecompError: the tool path reaches 0.024 m from the centre (tool radius included), beyond the unclamped half-width 0.015 m; enlarge the blank or reduce the clamp margin
