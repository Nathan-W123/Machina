# ML first shot + automatic engine fallback (8 test parts)

**SparLab simulation**, DSIF + rim pass, penalty 10. Tolerance 0.15 mm (vertical RMS): a part above it after the ML first shot (combined v1) gets ONE displacement-adjustment step from its simulated error, simulated again; the result is kept only if better. 8 parts: a small sample.

| part | ML first shot | engine round | final | simulations | decision |
|---|--:|--:|--:|--:|---|
| dome-s2026-0000 | 0.071 | - | 0.071 | 1 | within tolerance after the ML first shot |
| dome-s2026-0001 | 0.145 | - | 0.145 | 1 | within tolerance after the ML first shot |
| elliptic_cone-s2026-0000 | 0.202 | 0.256 | 0.202 | 2 | engine round rejected (not better: keep the first shot) |
| elliptic_cone-s2026-0001 | 0.205 | 0.251 | 0.205 | 2 | engine round rejected (not better: keep the first shot) |
| pyramid-s2026-0000 | 0.103 | - | 0.103 | 1 | within tolerance after the ML first shot |
| pyramid-s2026-0001 | - | - | - | nan | not run |
| truncated_cone-s2026-0000 | - | - | - | nan | not run |
| truncated_cone-s2026-0001 | - | - | - | nan | not run |
| **mean** | **0.145** | | **0.145** | **1.40** | |

Made by `python/scripts/combined_fallback.py`.
