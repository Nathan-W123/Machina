# Smoke test on pyramid-s2026-0001: optimiser and DSIF rim pass

**SparLab simulation** (sparlab_form; predictions are the stage_n26 MLP surrogate's). **One part: not statistically meaningful.** The part was chosen in advance because the ML first shot barely helped it.

All simulated rows at penalty 10 (`springback_fine` preset). Deviation formed - target, vertical, in mm; upper band = the top 1 mm of depth (next to the clamp), deep = the rest of the part.

| case | kind | FE runs | vertical RMS | upper-band RMS | upper-band bias | deep RMS | deep bias | note |
|---|---|--:|--:|--:|--:|--:|--:|---|
| plate uncompensated | simulated | 1 | 0.308 | 0.478 | -0.464 | 0.130 | -0.012 |  |
| plate FE-DA-1 | simulated | 2 | 0.299 | 0.394 | -0.366 | 0.226 | 0.128 |  |
| dsif uncompensated | simulated | 1 | 0.246 | 0.209 | -0.157 | 0.266 | 0.240 |  |
| ML-MLP n26 (surrogate DA) | predicted (MLP) | 0 | 0.268 | 0.264 | -0.218 | 0.271 | 0.212 |  |
| ML-MLP n26 (surrogate DA) | simulated | 1 | 0.268 | 0.359 | -0.329 | 0.196 | 0.143 |  |
| optimizer (from ML-MLP DA) | predicted (MLP) | 0 | 0.228 | 0.285 | -0.245 | 0.186 | 0.086 |  |
| optimizer (from ML-MLP DA) | simulated | 1 | 0.259 | 0.390 | -0.366 | 0.131 | 0.041 |  |
| dsif + rim pass uncompensated | simulated | 1 | 0.362 | 0.136 | 0.069 | 0.444 | 0.432 |  |

Optimiser: 5 evaluation batches, 263 s, stopped 'optimizer', in envelope True, largest change of the command from the DA command 0.253 mm.

Made by `python/scripts/smoke_optimizer_edge.py` (optimize, verify, edge, report). The optimiser's command was verified on the ML first shot's path (backing plate, 1 mm clearance); the rim-pass run is benchmarks/springback_fine's dsif setup with `rim_pass: true` (benchmarks/support_cone's dsif_rim).

## Reading (one part only)

* **Optimiser**: predicted 0.268 -> 0.228 mm. Simulated 0.268 -> 0.259 mm, about a quarter of the
  predicted gain. All of the gain is in the deep region (deep bias +0.14 -> +0.04 mm). The upper
  band got slightly worse (bias -0.33 -> -0.37 mm), as the MLP predicted (-0.22 -> -0.24). The
  surrogate sees no lever on the rim with the backing plate either, and it predicts the rim sag
  about 0.1 mm too small.
* **DSIF rim pass** (uncompensated): the rim pass removes the upper-band error (bias +0.07 mm,
  RMS 0.14 mm against 0.48 with the plate and 0.21 with DSIF alone). But the deep region comes
  out 0.43 mm too shallow (bias +0.43, plain DSIF +0.24), so the vertical RMS is the worst of the
  table at 0.362 mm. Only 5 % of the remaining squared error is in the upper band. This is the
  opposite of the plate's error. A compensated rim-pass part (DA on the deep region, rim held
  by the pass) is the natural next test.
