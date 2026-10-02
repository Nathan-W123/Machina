# Springback benchmark test parts on the fine mesh: SparLab simulations, not experiments

Setup `springback_fine` (0.833 mm incompatible-mode Hex8, 1 layer, 5 thickness points, penalty 10); the 8 test parts of benchmarks/springback. Deviation of the released part from the target over the part [mm]; rim sag = mean vertical deviation of the part less than 1 mm deep (negative = too deep); interior = the part deeper than that.

42 simulations complete, 0 failed records.

| strategy | method | FE runs / part | parts | vertical RMS mean | median | max mean | normal RMS | rim sag | upper band RMS | interior RMS (bias) | flange bias | runtime / run [s] | CPU h / part | Newton its | failed |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| backing_plate | uncompensated | 1 | 8 | **0.300** | 0.298 | 0.583 | 0.279 | -0.436 | 0.445 | 0.149 (-0.05) | -0.047 | 3568 | 0.99 | 5959 | 0 |
| backing_plate | FE-DA-1 | 2 | 8 | **0.263** | 0.252 | 0.535 | 0.250 | -0.334 | 0.356 | 0.162 (0.08) | -0.048 | 3352 | 1.92 | 5608 | 0 |
| backing_plate | FE-DA-2 | 3 | 8 | **0.272** | 0.259 | 0.595 | 0.260 | -0.318 | 0.344 | 0.188 (0.12) | -0.056 | 3336 | 2.85 | 5366 | 0 |
| dsif | uncompensated | 1 | 6 | **0.262** | 0.250 | 0.577 | 0.247 | -0.195 | 0.260 | 0.222 (0.14) | 0.028 | 2465 | 0.68 | 5621 | 0 |
| dsif | FE-DA-1 | 2 | 6 | **0.243** | 0.241 | 0.591 | 0.234 | -0.292 | 0.345 | 0.116 (0.03) | 0.005 | 2639 | 1.42 | 6062 | 0 |
| dsif | FE-DA-2 | 3 | 6 | **0.245** | 0.243 | 0.642 | 0.238 | -0.309 | 0.359 | 0.103 (0.05) | 0.000 | 2789 | 2.19 | 6481 | 0 |

The same over the 6 parts where every strategy and method has a result (dome-s2026-0000, dome-s2026-0001, elliptic_cone-s2026-0000, elliptic_cone-s2026-0001, pyramid-s2026-0000, truncated_cone-s2026-0000):

| strategy | method | vertical RMS mean | max mean | rim sag | interior RMS (bias) |
|---|---|--:|--:|--:|--:|
| backing_plate | uncompensated | **0.291** | 0.564 | -0.420 | 0.144 (-0.08) |
| backing_plate | FE-DA-1 | **0.236** | 0.488 | -0.315 | 0.118 (0.05) |
| backing_plate | FE-DA-2 | **0.241** | 0.540 | -0.300 | 0.138 (0.08) |
| dsif | uncompensated | **0.262** | 0.577 | -0.195 | 0.222 (0.14) |
| dsif | FE-DA-1 | **0.243** | 0.591 | -0.292 | 0.116 (0.03) |
| dsif | FE-DA-2 | **0.245** | 0.642 | -0.309 | 0.103 (0.05) |

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
| dome-s2026-0000 | 0.654 (1.07) / -0.83 | 0.307 (0.50) / -0.37 | 0.213 (0.41) / -0.26 | 0.203 (0.40) / -0.24 | 0.347 (0.62) / -0.43 | 0.281 (0.62) / -0.34 | 0.298 (0.60) / -0.37 | - | - | - |
| dome-s2026-0001 | 0.776 (1.40) / -1.16 | 0.290 (0.59) / -0.46 | 0.223 (0.51) / -0.36 | 0.216 (0.50) / -0.35 | 0.225 (0.63) / -0.33 | 0.310 (0.77) / -0.49 | 0.258 (0.76) / -0.36 | - | - | - |
| elliptic_cone-s2026-0000 | 0.621 (1.21) / -0.96 | 0.277 (0.57) / -0.42 | 0.252 (0.50) / -0.30 | 0.294 (0.69) / -0.26 | 0.254 (0.58) / -0.31 | 0.289 (0.66) / -0.39 | 0.327 (0.73) / -0.43 | - | - | - |
| elliptic_cone-s2026-0001 | 0.841 (1.51) / -1.24 | 0.305 (0.61) / -0.44 | 0.266 (0.56) / -0.36 | 0.260 (0.55) / -0.37 | 0.307 (0.67) / +0.00 | 0.200 (0.54) / -0.15 | 0.228 (0.62) / -0.28 | - | - | - |
| pyramid-s2026-0000 | 0.791 (1.45) / -1.14 | 0.290 (0.57) / -0.41 | 0.213 (0.48) / -0.29 | 0.217 (0.48) / -0.27 | 0.195 (0.48) / -0.04 | 0.190 (0.50) / -0.24 | 0.148 (0.50) / -0.12 | - | - | - |
| pyramid-s2026-0001 | 0.699 (1.38) / -1.08 | 0.308 (0.69) / -0.46 | 0.299 (0.61) / -0.37 | 0.304 (0.63) / -0.37 | - | - | - | - | - | - |
| truncated_cone-s2026-0000 | 0.766 (1.43) / -1.15 | 0.280 (0.55) / -0.42 | 0.251 (0.46) / -0.31 | 0.257 (0.61) / -0.30 | 0.246 (0.49) / -0.05 | 0.187 (0.46) / -0.13 | 0.212 (0.64) / -0.29 | - | - | - |
| truncated_cone-s2026-0001 | 0.703 (1.29) / -1.12 | 0.344 (0.60) / -0.51 | 0.386 (0.74) / -0.42 | 0.426 (0.89) / -0.38 | - | - | - | - | - | - |

Failures:

none
