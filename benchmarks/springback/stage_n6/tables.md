# stage_n6: SparLab simulations (not experiments)

12 train / 4 calibration / 8 test parts; 48 samples; 55 simulations complete, 0 failed.

| Method | FE runs per part | n parts | vertical RMS mean [mm] | median [mm] | vertical max mean [mm] | normal RMS mean [mm] | RMS / uncompensated | better than uncompensated | better than DA-1 |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| uncompensated | 1 | 8 | 0.731 | 0.735 | 1.341 | 0.729 | 1.00 | 0 | 0 |
| FE-DA-1 | 2 | 8 | 0.623 | 0.620 | 1.129 | 0.629 | 0.85 | 8 | 0 |
| ML-GBM | 1 | 8 | 0.627 | 0.628 | 1.136 | 0.634 | 0.86 | 7 | 5 |

Per part, vertical RMS [mm] (max |dev| in brackets):

| part | uncompensated | FE-DA-1 | ML-GBM |
|---|--:|--:|--:|
| dome-s2026-0000 | 0.654 (1.07) | 0.532 (0.90) | 0.654 (1.07) |
| dome-s2026-0001 | 0.776 (1.40) | 0.625 (1.25) | 0.578 (1.21) |
| elliptic_cone-s2026-0000 | 0.621 (1.21) | 0.512 (1.06) | 0.495 (1.01) |
| elliptic_cone-s2026-0001 | 0.841 (1.51) | 0.726 (1.30) | 0.735 (1.33) |
| pyramid-s2026-0000 | 0.791 (1.45) | 0.615 (1.19) | 0.600 (1.17) |
| pyramid-s2026-0001 | 0.699 (1.38) | 0.615 (1.09) | 0.615 (1.10) |
| truncated_cone-s2026-0000 | 0.767 (1.43) | 0.658 (1.25) | 0.641 (1.21) |
| truncated_cone-s2026-0001 | 0.703 (1.29) | 0.700 (0.98) | 0.699 (0.99) |

`*` target outside the model's training envelope (compensated with the override).

Mean over the test parts of the RMS vertical deviation per region [mm] (upper band: part nodes less than 1 mm deep; share: of the squared error over the part):

| Method | upper band RMS | bias | deep part RMS | bias | flange RMS | share of error in upper band |
|---|--:|--:|--:|--:|--:|--:|
| uncompensated | 1.099 | -1.086 | 0.328 | -0.137 | 0.376 | 0.87 |
| FE-DA-1 | 0.879 | -0.856 | 0.346 | 0.126 | 0.324 | 0.79 |
| ML-GBM | 0.883 | -0.860 | 0.356 | 0.127 | 0.325 | 0.78 |
