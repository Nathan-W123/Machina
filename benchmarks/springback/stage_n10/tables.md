# stage_n10: SparLab simulations (not experiments)

24 train / 8 calibration / 8 test parts; 80 samples; 104 simulations complete, 0 failed.

| Method | FE runs per part | n parts | vertical RMS mean [mm] | median [mm] | vertical max mean [mm] | normal RMS mean [mm] | RMS / uncompensated | better than uncompensated | better than DA-1 |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| uncompensated | 1 | 8 | 0.731 | 0.735 | 1.341 | 0.729 | 1.00 | 0 | 0 |
| FE-DA-1 | 2 | 8 | 0.623 | 0.620 | 1.129 | 0.629 | 0.85 | 8 | 0 |
| FE-DA-2 | 3 | 8 | 0.602 | 0.618 | 1.156 | 0.608 | 0.82 | 8 | 7 |
| ML-GBM | 1 | 8 | 0.615 | 0.606 | 1.115 | 0.622 | 0.84 | 8 | 4 |
| ML-MLP | 1 | 8 | 0.613 | 0.608 | 1.118 | 0.620 | 0.84 | 8 | 5 |

Per part, vertical RMS [mm] (max |dev| in brackets):

| part | uncompensated | FE-DA-1 | FE-DA-2 | ML-GBM | ML-MLP |
|---|--:|--:|--:|--:|--:|
| dome-s2026-0000 | 0.654 (1.07) | 0.532 (0.90) | 0.485 (0.83) | 0.539 (0.90) | 0.493 (0.85) |
| dome-s2026-0001 | 0.776 (1.40) | 0.625 (1.25) | 0.555 (1.17) | 0.581 (1.20) | 0.588 (1.21) |
| elliptic_cone-s2026-0000 | 0.621 (1.21) | 0.512 (1.06) | 0.493 (1.01) | 0.493 (1.00) | 0.512 (1.06) |
| elliptic_cone-s2026-0001 | 0.841 (1.51) | 0.726 (1.30) | 0.714 (1.29) | 0.731 (1.32) | 0.730 (1.31) |
| pyramid-s2026-0000 | 0.791 (1.45) | 0.615 (1.19) | 0.604 (1.17) | 0.597 (1.16) | 0.600 (1.17) |
| pyramid-s2026-0001 | 0.699 (1.38) | 0.615 (1.09) | 0.648 (1.33) | 0.615 (1.09) | 0.616 (1.10) |
| truncated_cone-s2026-0000 | 0.767 (1.43) | 0.658 (1.25) | 0.632 (1.20) | 0.665 (1.26) | 0.663 (1.26) |
| truncated_cone-s2026-0001 | 0.703 (1.29) | 0.700 (0.98) | 0.686 (1.25) | 0.699 (0.98) | 0.699 (0.99) |

`*` target outside the model's training envelope (compensated with the override).

Mean over the test parts of the RMS vertical deviation per region [mm] (upper band: part nodes less than 1 mm deep; share: of the squared error over the part):

| Method | upper band RMS | bias | deep part RMS | bias | flange RMS | share of error in upper band |
|---|--:|--:|--:|--:|--:|--:|
| uncompensated | 1.099 | -1.086 | 0.328 | -0.137 | 0.376 | 0.87 |
| FE-DA-1 | 0.879 | -0.856 | 0.346 | 0.126 | 0.324 | 0.79 |
| FE-DA-2 | 0.905 | -0.884 | 0.272 | 0.021 | 0.329 | 0.87 |
| ML-GBM | 0.869 | -0.845 | 0.342 | 0.142 | 0.322 | 0.79 |
| ML-MLP | 0.869 | -0.845 | 0.335 | 0.140 | 0.321 | 0.79 |
