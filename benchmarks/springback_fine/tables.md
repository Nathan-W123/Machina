# Springback benchmark test parts on the fine mesh: SparLab simulations, not experiments

Setup `springback_fine` (0.833 mm incompatible-mode Hex8, 1 layer, 5 thickness points, penalty 10); the 8 test parts of benchmarks/springback. Deviation of the released part from the target over the part [mm]; rim sag = mean vertical deviation of the part less than 1 mm deep (negative = too deep); interior = the part deeper than that.

11 simulations complete, 1 failed records.

| strategy | method | FE runs / part | parts | vertical RMS mean | median | max mean | normal RMS | rim sag | upper band RMS | interior RMS (bias) | flange bias | runtime / run [s] | CPU h / part | Newton its | failed |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| backing_plate | uncompensated | 1 | 4 | **0.295** | 0.297 | 0.562 | 0.275 | -0.411 | 0.420 | 0.150 (-0.08) | -0.048 | 2028 | 0.56 | 3564 | 0 |
| backing_plate | FE-DA-1 | 2 | 4 | **0.236** | 0.233 | 0.487 | 0.226 | -0.304 | 0.327 | 0.122 (0.04) | -0.039 | 2972 | 1.39 | 5051 | 0 |
| backing_plate | FE-DA-2 | 3 | 3 | **0.227** | 0.217 | 0.478 | 0.217 | -0.294 | 0.319 | 0.100 (0.06) | -0.050 | 2125 | 1.84 | 3785 | 1 |

The same over the 3 parts where every strategy and method has a result (dome-s2026-0000, elliptic_cone-s2026-0001, pyramid-s2026-0000):

| strategy | method | vertical RMS mean | max mean | rim sag | interior RMS (bias) |
|---|---|--:|--:|--:|--:|
| backing_plate | uncompensated | **0.300** | 0.560 | -0.408 | 0.167 (-0.10) |
| backing_plate | FE-DA-1 | **0.230** | 0.484 | -0.305 | 0.099 (0.03) |
| backing_plate | FE-DA-2 | **0.227** | 0.478 | -0.294 | 0.100 (0.06) |

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
| dome-s2026-0000 | 0.654 (1.07) / -0.83 | - | - | - | 0.307 (0.50) / -0.37 | 0.213 (0.41) / -0.26 | 0.203 (0.40) / -0.24 | - | - | - |
| elliptic_cone-s2026-0000 | 0.621 (1.21) / -0.96 | - | - | - | 0.277 (0.57) / -0.42 | 0.252 (0.50) / -0.30 | FAILED | - | - | - |
| elliptic_cone-s2026-0001 | 0.841 (1.51) / -1.24 | - | - | - | 0.305 (0.61) / -0.44 | 0.266 (0.56) / -0.36 | 0.260 (0.55) / -0.37 | - | - | - |
| pyramid-s2026-0000 | 0.791 (1.45) / -1.14 | - | - | - | 0.290 (0.57) / -0.41 | 0.213 (0.48) / -0.29 | 0.217 (0.48) / -0.27 | - | - | - |

Failures:

* elliptic_cone-s2026-0000 backing_plate FE-DA-2: PrecompError: the tool path reaches 0.024 m from the centre (tool radius included), beyond the unclamped half-width 0.015 m; enlarge the blank or reduce the clamp margin
