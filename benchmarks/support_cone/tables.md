# Rim support on one small cone (SparLab simulations)

Deviation of the released part from the target [mm]; rim sag = the upper band's bias (part < 1 mm deep; negative = too deep).

| case | FE runs | vertical RMS | max | bias | normal RMS | upper band RMS | rim sag (upper band bias) | deep RMS | flange RMS | runtime [s] | increments | iterations | cuts | warnings |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| dsif | 1 | 0.275 | 0.537 | 0.070 | 0.263 | 0.258 | -0.221 | 0.284 | 0.069 | 407 | 516 | 3282 | 6 | 1 |
| dsif+DA1 | 2 | 0.298 | 0.646 | 0.194 | 0.278 | 0.161 | -0.053 | 0.355 | 0.142 | 702 | 571 | 3680 | 11 | 1 |
| dsif_rim | 1 | 0.361 | 0.692 | 0.250 | 0.340 | 0.142 | -0.037 | 0.442 | 0.093 | 477 | 729 | 3848 | 7 | 2 |
| dsif_rim+DA1 | 2 | 0.258 | 0.597 | 0.106 | 0.248 | 0.188 | -0.142 | 0.291 | 0.131 | 793 | 750 | 3891 | 10 | 2 |
| dsif_t0 | 1 | 0.357 | 0.745 | -0.129 | 0.340 | 0.511 | -0.494 | 0.222 | 0.099 | 397 | 525 | 3255 | 11 | 1 |
| dsif_t0+DA1 | 2 | 0.318 | 0.669 | 0.148 | 0.300 | 0.229 | -0.163 | 0.360 | 0.092 | 550 | 421 | 2928 | 6 | 1 |
| plate | 1 | 0.370 | 0.758 | -0.202 | 0.344 | 0.553 | -0.541 | 0.193 | 0.098 | 773 | 523 | 5007 | 4 | 1 |
| plate+DA1 | 2 | 0.297 | 0.626 | -0.071 | 0.285 | 0.423 | -0.399 | 0.187 | 0.084 | 335 | 370 | 3605 | 2 | 1 |
| spif | 1 | 0.767 | 1.426 | -0.545 | 0.766 | 1.164 | -1.151 | 0.367 | 0.399 | 302 | 288 | 1906 | 9 | 1 |
| spif+DA1 | 2 | 0.658 | 1.248 | -0.363 | 0.667 | 1.011 | -0.992 | 0.295 | 0.364 | 143 | 136 | 1037 | 2 | 1 |

Tools (force of the sheet on the tool; fz < 0: pushed down):

| case | tool | max nodes in contact | max penetration [um] | peak force [N] | fz min [N] | fz max [N] | steps |
|---|---|--:|--:|--:|--:|--:|---|
| dsif | tool | 4 | 0.83 | 1911 | 0 | 1700 | 1 |
| dsif | support | 3 | 1.12 | 1897 | -1572 | 0 | 1 |
| dsif+DA1 | tool | 6 | 0.93 | 2044 | 0 | 1768 | 1 |
| dsif+DA1 | support | 3 | 1.22 | 2114 | -1654 | 0 | 1 |
| dsif_rim | tool | 4 | 0.83 | 1911 | 0 | 1700 | 1 |
| dsif_rim | support | 3 | 1.12 | 1897 | -1572 | 0 | 1,3 |
| dsif_rim+DA1 | tool | 6 | 1.02 | 2133 | 0 | 1777 | 1 |
| dsif_rim+DA1 | support | 3 | 1.27 | 1973 | -1622 | 0 | 1,3 |
| dsif_t0 | tool | 3 | 0.65 | 1233 | 0 | 1191 | 1 |
| dsif_t0 | support | 2 | 0.52 | 886 | -809 | 0 | 1 |
| dsif_t0+DA1 | tool | 4 | 0.64 | 1246 | 0 | 1184 | 1 |
| dsif_t0+DA1 | support | 3 | 0.58 | 993 | -910 | 0 | 1 |
| plate | tool | 3 | 0.71 | 1362 | 0 | 1318 | 1 |
| plate | plate | 27 | 0.49 | 2113 | -2111 | 0 | 1 |
| plate+DA1 | tool | 4 | 0.70 | 1532 | 0 | 1499 | 1 |
| plate+DA1 | plate | 26 | 0.45 | 2221 | -2221 | 0 | 1 |
| spif | tool | 3 | 0.61 | 973 | 0 | 968 | 1 |
| spif+DA1 | tool | 3 | 0.59 | 1003 | 0 | 988 | 1 |
