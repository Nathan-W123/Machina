### Table 5 - second check part truncated_cone-s2026-0001 (wall 54 deg) [mm]

Released-shape difference against `c1_h0.625_L1_im_tp5`; errors against the band of Richardson limits of `richardson_fit_c1.csv` (empty until three IM meshes exist).

| Case | vert. RMS | vert. max | rim sag | springback RMS | released vs finest RMS / max | err. RMS | err. max | err. rim sag | wall [s] (load) |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 2 mm, 2 std layers (**benchmark**) | 0.703 | 1.287 | -1.122 | 0.150 | 0.080 / 0.256 | +0.057..+0.061 | +0.111..+0.180 | -0.086..-0.096 | 178 (4.0-4.0) |
| 2 mm, 1 IM layer 2x2x5 | 0.748 | 1.365 | -1.198 | 0.167 | 0.114 / 0.208 | +0.103..+0.106 | +0.189..+0.258 | -0.162..-0.172 | 190 (4.0-4.0) |
| 1 mm, 1 IM layer 2x2x5 | 0.664 | 1.209 | -1.065 | 0.171 | 0.019 / 0.060 | +0.019..+0.023 | +0.033..+0.102 | -0.029..-0.040 | 1147 (4.2-4.0) |
| 0.833 mm, 1 IM layer 2x2x5 | 0.657 | 1.187 | -1.052 | 0.164 | 0.012 / 0.044 | +0.011..+0.015 | +0.011..+0.080 | -0.016..-0.027 | 1965 (4.0-4.2) |
| 0.833 mm, 1 IM layer 2x2x5, penalty 3 | 0.659 | 1.192 | -1.055 | 0.165 | 0.015 / 0.055 | +0.013..+0.017 | +0.016..+0.085 | -0.019..-0.030 | 1172 (4.2-4.2) |
| 0.625 mm, 1 IM layer 2x2x5 (finest) | 0.651 | 1.178 | -1.043 | 0.165 | 0 (finest) | +0.005..+0.009 | +0.002..+0.071 | -0.006..-0.017 | 3924 (3.9-3.0) |
