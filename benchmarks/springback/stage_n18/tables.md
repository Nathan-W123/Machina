# stage_n18: SparLab simulations (not experiments)

48 train / 16 calibration / 8 test parts; 144 samples; 166 simulations complete, 0 failed.

| Method | FE runs per part | n parts | vertical RMS mean [mm] | median [mm] | vertical max mean [mm] | normal RMS mean [mm] | RMS / uncompensated | better than uncompensated | better than DA-1 |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| uncompensated | 1 | 8 | 0.731 | 0.735 | 1.341 | 0.729 | 1.00 | 0 | 0 |
| FE-DA-1 | 2 | 8 | 0.623 | 0.620 | 1.129 | 0.629 | 0.85 | 8 | 0 |
| FE-DA-2 | 3 | 8 | 0.602 | 0.618 | 1.156 | 0.608 | 0.82 | 8 | 7 |
| ML-GBM | 1 | 8 | 0.634 | 0.634 | 1.147 | 0.640 | 0.87 | 7 | 4 |
| ML-MLP | 1 | 8 | 0.626 | 0.631 | 1.165 | 0.631 | 0.86 | 7 | 5 |

Per part, vertical RMS [mm] (max |dev| in brackets):

| part | uncompensated | FE-DA-1 | FE-DA-2 | ML-GBM | ML-MLP |
|---|--:|--:|--:|--:|--:|
| dome-s2026-0000 | 0.654 (1.07) | 0.532 (0.90) | 0.485 (0.83) | 0.654 (1.07) | 0.654 (1.07) |
| dome-s2026-0001 | 0.776 (1.40) | 0.625 (1.25) | 0.555 (1.17) | 0.592 (1.22) | 0.616 (1.25) |
| elliptic_cone-s2026-0000 | 0.621 (1.21) | 0.512 (1.06) | 0.493 (1.01) | 0.513 (1.06) * | 0.495 (1.03) * |
| elliptic_cone-s2026-0001 | 0.841 (1.51) | 0.726 (1.30) | 0.714 (1.29) | 0.733 (1.32) | 0.728 (1.31) |
| pyramid-s2026-0000 | 0.791 (1.45) | 0.615 (1.19) | 0.604 (1.17) | 0.601 (1.16) | 0.593 (1.15) |
| pyramid-s2026-0001 | 0.699 (1.38) | 0.615 (1.09) | 0.648 (1.33) | 0.614 (1.09) | 0.614 (1.09) |
| truncated_cone-s2026-0000 | 0.767 (1.43) | 0.658 (1.25) | 0.632 (1.20) | 0.670 (1.26) | 0.659 (1.25) |
| truncated_cone-s2026-0001 | 0.703 (1.29) | 0.700 (0.98) | 0.686 (1.25) | 0.699 (0.98) | 0.647 (1.19) |

`*` target outside the model's training envelope (compensated with the override).

Mean over the test parts of the RMS vertical deviation per region [mm] (upper band: part nodes less than 1 mm deep; share: of the squared error over the part):

| Method | upper band RMS | bias | deep part RMS | bias | flange RMS | share of error in upper band |
|---|--:|--:|--:|--:|--:|--:|
| uncompensated | 1.099 | -1.086 | 0.328 | -0.137 | 0.376 | 0.87 |
| FE-DA-1 | 0.879 | -0.856 | 0.346 | 0.126 | 0.324 | 0.79 |
| FE-DA-2 | 0.905 | -0.884 | 0.272 | 0.021 | 0.329 | 0.87 |
| ML-GBM | 0.896 | -0.874 | 0.353 | 0.115 | 0.328 | 0.79 |
| ML-MLP | 0.908 | -0.886 | 0.330 | 0.054 | 0.329 | 0.82 |
