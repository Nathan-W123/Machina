# Springback benchmark test parts on the fine mesh: SparLab simulations, not experiments

Setup `springback_fine` (0.833 mm incompatible-mode Hex8, 1 layer, 5 thickness points, penalty 10); the 8 test parts of benchmarks/springback. Deviation of the released part from the target over the part [mm]; rim sag = mean vertical deviation of the part less than 1 mm deep (negative = too deep); interior = the part deeper than that.

8 simulations complete, 1 failed records.

| strategy | method | FE runs / part | parts | vertical RMS mean | median | max mean | normal RMS | rim sag | upper band RMS | interior RMS (bias) | flange bias | runtime / run [s] | CPU h / part | Newton its | failed |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| backing_plate | uncompensated | 1 | 3 | **0.296** | 0.305 | 0.559 | 0.278 | -0.413 | 0.423 | 0.142 (-0.07) | -0.048 | 1727 | 0.48 | 3058 | 0 |
| backing_plate | FE-DA-1 | 2 | 3 | **0.244** | 0.252 | 0.489 | 0.234 | -0.307 | 0.331 | 0.126 (0.05) | -0.041 | 2308 | 1.12 | 3933 | 0 |
| backing_plate | FE-DA-2 | 3 | 2 | **0.231** | 0.231 | 0.475 | 0.222 | -0.307 | 0.331 | 0.076 (0.04) | -0.050 | 2365 | 1.44 | 4200 | 1 |

The same over the 2 parts where every strategy and method has a result (dome-s2026-0000, elliptic_cone-s2026-0001):

| strategy | method | vertical RMS mean | max mean | rim sag | interior RMS (bias) |
|---|---|--:|--:|--:|--:|
| backing_plate | uncompensated | **0.306** | 0.553 | -0.408 | 0.165 (-0.10) |
| backing_plate | FE-DA-1 | **0.239** | 0.486 | -0.310 | 0.093 (0.05) |
| backing_plate | FE-DA-2 | **0.231** | 0.475 | -0.307 | 0.076 (0.04) |

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

Failures:

* elliptic_cone-s2026-0000 backing_plate FE-DA-2: PrecompError: the tool path reaches 0.024 m from the centre (tool radius included), beyond the unclamped half-width 0.015 m; enlarge the blank or reduce the clamp margin
