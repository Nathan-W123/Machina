# Springback benchmark test parts on the fine mesh: SparLab simulations, not experiments

Setup `springback_fine_bulk` (0.833 mm incompatible-mode Hex8, 1 layer, 5 thickness points, penalty 3); the 8 test parts of benchmarks/springback. Deviation of the released part from the target over the part [mm]; rim sag = mean vertical deviation of the part less than 1 mm deep (negative = too deep); interior = the part deeper than that.

8 simulations complete, 0 failed records.

| strategy | method | FE runs / part | parts | vertical RMS mean | median | max mean | normal RMS | rim sag | upper band RMS | interior RMS (bias) | flange bias | runtime / run [s] | CPU h / part | Newton its | failed |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| backing_plate | uncompensated | 1 | 4 | **0.291** | 0.289 | 0.553 | 0.272 | -0.410 | 0.419 | 0.145 (-0.09) | -0.050 | 1164 | 0.32 | 1932 | 0 |
| dsif | uncompensated | 1 | 4 | **0.261** | 0.250 | 0.541 | 0.248 | -0.215 | 0.264 | 0.216 (0.11) | 0.016 | 1992 | 0.55 | 3906 | 0 |

The same over the 4 parts where every strategy and method has a result (dome-s2026-0000, elliptic_cone-s2026-0000, pyramid-s2026-0000, truncated_cone-s2026-0000):

| strategy | method | vertical RMS mean | max mean | rim sag | interior RMS (bias) |
|---|---|--:|--:|--:|--:|
| backing_plate | uncompensated | **0.291** | 0.553 | -0.410 | 0.145 (-0.09) |
| dsif | uncompensated | **0.261** | 0.541 | -0.215 | 0.216 (0.11) |

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
| dome-s2026-0000 | 0.654 (1.07) / -0.83 | 0.307 (0.50) / -0.37 | - | - | 0.349 (0.62) / -0.43 | - | - | - | - | - |
| elliptic_cone-s2026-0000 | 0.621 (1.21) / -0.96 | 0.279 (0.58) / -0.43 | - | - | 0.256 (0.59) / -0.32 | - | - | - | - | - |
| pyramid-s2026-0000 | 0.791 (1.45) / -1.14 | 0.290 (0.58) / -0.41 | - | - | 0.192 (0.47) / -0.05 | - | - | - | - | - |
| truncated_cone-s2026-0000 | 0.766 (1.43) / -1.15 | 0.289 (0.56) / -0.43 | - | - | 0.245 (0.49) / -0.06 | - | - | - | - | - |

Failures:

none
