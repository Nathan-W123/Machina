# stage_n14: SparLab simulations (not experiments)

36 train / 12 calibration / 8 test parts; 112 samples; 136 simulations complete, 0 failed.

| Method | FE runs per part | n parts | vertical RMS mean [mm] | median [mm] | vertical max mean [mm] | normal RMS mean [mm] | RMS / uncompensated | better than uncompensated | better than DA-1 |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| uncompensated | 1 | 8 | 0.731 | 0.735 | 1.341 | 0.729 | 1.00 | 0 | 0 |
| FE-DA-1 | 2 | 8 | 0.623 | 0.620 | 1.129 | 0.629 | 0.85 | 8 | 0 |
| FE-DA-2 | 3 | 8 | 0.602 | 0.618 | 1.156 | 0.608 | 0.82 | 8 | 7 |
| ML-GBM | 1 | 8 | 0.621 | 0.616 | 1.127 | 0.627 | 0.85 | 8 | 3 |
| ML-MLP | 1 | 8 | 0.605 | 0.606 | 1.136 | 0.611 | 0.83 | 8 | 7 |

Per part, vertical RMS [mm] (max |dev| in brackets):

| part | uncompensated | FE-DA-1 | FE-DA-2 | ML-GBM | ML-MLP |
|---|--:|--:|--:|--:|--:|
| dome-s2026-0000 | 0.654 (1.07) | 0.532 (0.90) | 0.485 (0.83) | 0.544 (0.91) | 0.503 (0.86) |
| dome-s2026-0001 | 0.776 (1.40) | 0.625 (1.25) | 0.555 (1.17) | 0.583 (1.20) | 0.598 (1.22) |
| elliptic_cone-s2026-0000 | 0.621 (1.21) | 0.512 (1.06) | 0.493 (1.01) | 0.514 (1.06) * | 0.492 (1.00) * |
| elliptic_cone-s2026-0001 | 0.841 (1.51) | 0.726 (1.30) | 0.714 (1.29) | 0.731 (1.32) | 0.731 (1.32) |
| pyramid-s2026-0000 | 0.791 (1.45) | 0.615 (1.19) | 0.604 (1.17) | 0.617 (1.19) | 0.599 (1.16) |
| pyramid-s2026-0001 | 0.699 (1.38) | 0.615 (1.09) | 0.648 (1.33) | 0.615 (1.09) | 0.614 (1.09) |
| truncated_cone-s2026-0000 | 0.767 (1.43) | 0.658 (1.25) | 0.632 (1.20) | 0.662 (1.25) | 0.655 (1.24) |
| truncated_cone-s2026-0001 | 0.703 (1.29) | 0.700 (0.98) | 0.686 (1.25) | 0.699 (0.99) | 0.651 (1.20) |

`*` target outside the model's training envelope (compensated with the override).

Mean over the test parts of the RMS vertical deviation per region [mm] (upper band: part nodes less than 1 mm deep; share: of the squared error over the part):

| Method | upper band RMS | bias | deep part RMS | bias | flange RMS | share of error in upper band |
|---|--:|--:|--:|--:|--:|--:|
| uncompensated | 1.099 | -1.086 | 0.328 | -0.137 | 0.376 | 0.87 |
| FE-DA-1 | 0.879 | -0.856 | 0.346 | 0.126 | 0.324 | 0.79 |
| FE-DA-2 | 0.905 | -0.884 | 0.272 | 0.021 | 0.329 | 0.87 |
| ML-GBM | 0.881 | -0.858 | 0.337 | 0.127 | 0.325 | 0.80 |
| ML-MLP | 0.888 | -0.865 | 0.302 | 0.069 | 0.324 | 0.83 |
