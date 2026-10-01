# Springback benchmark test parts on the fine mesh: SparLab simulations, not experiments

Setup `springback_fine` (0.833 mm incompatible-mode Hex8, 1 layer, 5 thickness points, penalty 10); the 8 test parts of benchmarks/springback. Deviation of the released part from the target over the part [mm]; rim sag = mean vertical deviation of the part less than 1 mm deep (negative = too deep); interior = the part deeper than that.

30 simulations complete, 0 failed records.

| strategy | method | FE runs / part | parts | vertical RMS mean | median | max mean | normal RMS | rim sag | upper band RMS | interior RMS (bias) | flange bias | runtime / run [s] | CPU h / part | Newton its | failed |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| backing_plate | uncompensated | 1 | 8 | **0.307** | 0.301 | 0.628 | 0.287 | -0.436 | 0.445 | 0.166 (-0.03) | -0.046 | 3526 | 0.98 | 5927 | 0 |
| backing_plate | FE-DA-1 | 2 | 8 | **0.305** | 0.253 | 0.566 | 0.295 | -0.321 | 0.343 | 0.244 (0.17) | -0.044 | 3000 | 1.81 | 4974 | 0 |
| backing_plate | FE-DA-2 | 3 | 8 | **0.272** | 0.256 | 0.537 | 0.260 | -0.332 | 0.356 | 0.191 (0.04) | -0.045 | 3276 | 2.72 | 5390 | 0 |
| dsif | uncompensated | 1 | 2 | **0.339** | 0.339 | 0.846 | 0.332 | -0.357 | 0.383 | 0.259 (0.03) | -0.052 | 1816 | 0.50 | 3932 | 0 |
| dsif | FE-DA-1 | 2 | 2 | **0.271** | 0.271 | 0.621 | 0.265 | -0.355 | 0.393 | 0.110 (-0.04) | -0.065 | 1885 | 1.03 | 3990 | 0 |
| dsif | FE-DA-2 | 3 | 2 | **0.246** | 0.246 | 0.619 | 0.238 | -0.272 | 0.357 | 0.134 (0.02) | -0.039 | 1867 | 1.55 | 3738 | 0 |

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
| dome-s2026-0000 | 0.654 (1.07) / -0.83 | - | - | - | 0.298 (0.49) / -0.37 | 0.211 (0.41) / -0.26 | 0.202 (0.39) / -0.24 | 0.341 (0.61) / -0.42 | 0.252 (0.57) / -0.29 | 0.183 (0.39) / -0.09 |
| dome-s2026-0001 | 0.776 (1.40) / -1.16 | - | - | - | 0.358 (0.95) / -0.46 | 0.262 (0.52) / -0.36 | 0.226 (0.50) / -0.34 | 0.336 (1.08) / -0.29 | 0.290 (0.67) / -0.42 | 0.310 (0.85) / -0.45 |
| elliptic_cone-s2026-0000 | 0.621 (1.21) / -0.96 | - | - | - | 0.277 (0.57) / -0.42 | 0.239 (0.50) / -0.32 | 0.266 (0.61) / -0.27 | - | - | - |
| elliptic_cone-s2026-0001 | 0.841 (1.51) / -1.24 | - | - | - | 0.305 (0.61) / -0.44 | 0.285 (0.55) / -0.37 | 0.279 (0.55) / -0.36 | - | - | - |
| pyramid-s2026-0000 | 0.791 (1.45) / -1.14 | - | - | - | 0.290 (0.57) / -0.41 | 0.245 (0.48) / -0.28 | 0.236 (0.47) / -0.27 | - | - | - |
| pyramid-s2026-0001 | 0.699 (1.38) / -1.08 | - | - | - | 0.308 (0.69) / -0.46 | 0.426 (0.73) / -0.31 | 0.349 (0.67) / -0.40 | - | - | - |
| truncated_cone-s2026-0000 | 0.766 (1.43) / -1.15 | - | - | - | 0.280 (0.55) / -0.42 | 0.226 (0.48) / -0.32 | 0.245 (0.53) / -0.31 | - | - | - |
| truncated_cone-s2026-0001 | 0.703 (1.29) / -1.12 | - | - | - | 0.344 (0.60) / -0.51 | 0.550 (0.86) / -0.34 | 0.371 (0.59) / -0.47 | - | - | - |

Failures:

none
