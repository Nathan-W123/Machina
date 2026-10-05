# ML first shot + ML second round (8 test parts)

**SparLab simulation**, DSIF + rim pass, penalty 10. Tolerance 0.15 mm (vertical RMS). A part above it after the ML first shot gets a second round: the surrogate is corrected by the first shot's simulated model error (smoothed at 1 mm) and the compensation is re-run on it; kept only if better. 'Engine round' = the stopped full-step displacement-adjustment round (combined_fallback), for comparison. 8 parts: a small sample.

| part | ML first shot | engine round (full DA step) | ML round predicted | ML round | final | simulations | decision |
|---|--:|--:|--:|--:|--:|--:|---|
| dome-s2026-0000 | 0.071 | - | - | - | 0.071 | 1 | within tolerance after the ML first shot |
| dome-s2026-0001 | 0.145 | - | - | - | 0.145 | 1 | within tolerance after the ML first shot |
| elliptic_cone-s2026-0000 | 0.202 | 0.256 | 0.178 | 0.201 | 0.201 | 2 | ML round kept (better) |
| elliptic_cone-s2026-0001 | 0.205 | 0.251 | 0.194 | 0.193 | 0.193 | 2 | ML round kept (better) |
| pyramid-s2026-0000 | 0.103 | - | - | - | 0.103 | 1 | within tolerance after the ML first shot |
| pyramid-s2026-0001 | 0.216 | - | 0.213 | 0.226 | 0.216 | 2 | ML round rejected (not better: keep the first shot) |
| truncated_cone-s2026-0000 | 0.205 | - | 0.173 | 0.226 | 0.205 | 2 | ML round rejected (not better: keep the first shot) |
| truncated_cone-s2026-0001 | 0.365 | - | 0.333 | 0.351 | 0.351 | 2 | ML round kept (better) |
| **mean** | **0.189** | | | | **0.186** | **1.62** | |

Made by `python/scripts/combined_ml2.py`.
