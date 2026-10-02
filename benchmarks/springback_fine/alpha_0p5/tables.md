# Springback benchmark test parts on the fine mesh: SparLab simulations, not experiments

Setup `springback_fine` (0.833 mm incompatible-mode Hex8, 1 layer, 5 thickness points, penalty 10); the 8 test parts of benchmarks/springback. Deviation of the released part from the target over the part [mm]; rim sag = mean vertical deviation of the part less than 1 mm deep (negative = too deep); interior = the part deeper than that.

9 simulations complete, 0 failed records.

| strategy | method | FE runs / part | parts | vertical RMS mean | median | max mean | normal RMS | rim sag | upper band RMS | interior RMS (bias) | flange bias | runtime / run [s] | CPU h / part | Newton its | failed |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| none | uncompensated | 1 | 1 | **0.657** | 0.657 | 1.187 | 0.652 | -1.052 | 1.060 | 0.205 (-0.00) | -0.221 | 1841 | 0.51 | 4082 | 0 |
| none | FE-DA-1 | 2 | 1 | **0.601** | 0.601 | 1.099 | 0.600 | -0.912 | 0.925 | 0.288 (0.13) | -0.204 | 1529 | 0.94 | 2847 | 0 |
| none | FE-DA-2 | 3 | 1 | **0.593** | 0.593 | 1.050 | 0.593 | -0.884 | 0.897 | 0.311 (0.16) | -0.213 | 1429 | 1.33 | 2487 | 0 |
| backing_plate | uncompensated | 1 | 1 | **0.344** | 0.344 | 0.595 | 0.319 | -0.507 | 0.512 | 0.194 (0.10) | -0.048 | 5063 | 1.41 | 8107 | 0 |
| backing_plate | FE-DA-1 | 2 | 1 | **0.366** | 0.366 | 0.650 | 0.345 | -0.446 | 0.454 | 0.306 (0.19) | -0.057 | 2055 | 1.98 | 3022 | 0 |
| backing_plate | FE-DA-2 | 3 | 1 | **0.368** | 0.368 | 0.702 | 0.343 | -0.479 | 0.491 | 0.276 (0.12) | -0.070 | 2607 | 2.70 | 3668 | 0 |
| dsif | uncompensated | 1 | 1 | **0.364** | 0.364 | 0.648 | 0.342 | -0.399 | 0.433 | 0.319 (0.27) | 0.010 | 2967 | 0.82 | 6845 | 0 |
| dsif | FE-DA-1 | 2 | 1 | **0.401** | 0.401 | 0.786 | 0.380 | -0.390 | 0.411 | 0.395 (0.31) | 0.040 | 2315 | 1.47 | 4604 | 0 |
| dsif | FE-DA-2 | 3 | 1 | **0.396** | 0.396 | 0.734 | 0.374 | -0.429 | 0.454 | 0.359 (0.26) | 0.016 | 2654 | 2.20 | 5062 | 0 |

The same over the 1 parts where every strategy and method has a result (truncated_cone-s2026-0001):

| strategy | method | vertical RMS mean | max mean | rim sag | interior RMS (bias) |
|---|---|--:|--:|--:|--:|
| none | uncompensated | **0.657** | 1.187 | -1.052 | 0.205 (-0.00) |
| none | FE-DA-1 | **0.601** | 1.099 | -0.912 | 0.288 (0.13) |
| none | FE-DA-2 | **0.593** | 1.050 | -0.884 | 0.311 (0.16) |
| backing_plate | uncompensated | **0.344** | 0.595 | -0.507 | 0.194 (0.10) |
| backing_plate | FE-DA-1 | **0.366** | 0.650 | -0.446 | 0.306 (0.19) |
| backing_plate | FE-DA-2 | **0.368** | 0.702 | -0.479 | 0.276 (0.12) |
| dsif | uncompensated | **0.364** | 0.648 | -0.399 | 0.319 (0.27) |
| dsif | FE-DA-1 | **0.401** | 0.786 | -0.390 | 0.395 (0.31) |
| dsif | FE-DA-2 | **0.396** | 0.734 | -0.429 | 0.359 (0.26) |

The 2 mm benchmark (benchmarks/springback, stage n18, single-point forming) on the same parts:

| method | vertical RMS mean | max mean | rim sag | upper band RMS | interior RMS (bias) |
|---|--:|--:|--:|--:|--:|
| uncompensated | 0.731 | 1.341 | -1.086 | 1.099 | 0.328 (-0.14) |
| FE-DA-1 | 0.623 | 1.129 | -0.856 | 0.879 | 0.346 (0.13) |
| FE-DA-2 | 0.602 | 1.156 | -0.884 | 0.905 | 0.272 (0.02) |

2 mm runs: 243 s mean, 227 s median per run (stage n18, two at a time on a shared machine).

Per part, vertical RMS [mm] (max |dev| in brackets; rim sag after the slash):

| part | 2 mm uncomp. | none uncompensated | none FE-DA-1 | none FE-DA-2 | backing_plate uncompensated | backing_plate FE-DA-1 | backing_plate FE-DA-2 | dsif uncompensated | dsif FE-DA-1 | dsif FE-DA-2 |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| truncated_cone-s2026-0001 | 0.703 (1.29) / -1.12 | 0.657 (1.19) / -1.05 | 0.601 (1.10) / -0.91 | 0.593 (1.05) / -0.88 | 0.344 (0.60) / -0.51 | 0.366 (0.65) / -0.45 | 0.368 (0.70) / -0.48 | 0.364 (0.65) / -0.40 | 0.401 (0.79) / -0.39 | 0.396 (0.73) / -0.43 |

Single-point forming, fine minus 2 mm, mean over the parts [mm]:

| method | parts | vertical RMS | max | rim sag | interior RMS |
|---|--:|--:|--:|--:|--:|
| uncompensated | 1 | -0.046 | -0.100 | +0.070 | -0.027 |
| FE-DA-1 | 1 | -0.099 | +0.115 | -0.142 | -0.355 |
| FE-DA-2 | 1 | -0.093 | -0.197 | +0.180 | +0.015 |

Failures:

none
