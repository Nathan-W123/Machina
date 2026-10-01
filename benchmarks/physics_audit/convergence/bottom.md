* **No, the benchmark discretisation is not converged.** On the
  representative part (truncated cone t-0000, wall 38.6 deg) its released
  shape is 0.10-0.12 mm RMS and up to 0.24 mm (worst point) away from the
  converged answer of the *same* model. It overstates the RMS deviation from
  target by 0.085-0.107 mm (0.767 against a Richardson limit of 0.66-0.68 mm,
  +14 %), the largest deviation by 0.19-0.24 mm and the rim sag by
  0.12-0.15 mm (-1.151 against -1.00 to -1.03 mm), and puts the untouched
  floor 0.20 mm too high. A **second test part** (truncated cone t-0001,
  wall 54 deg, added on resume) shows the same picture: 0.080 mm RMS /
  0.256 mm max from its finest mesh, RMS deviation +0.06 mm, max +0.11-0.18,
  rim sag 0.09 mm too much sag (Table 5). The springback itself
  (0.12-0.17 mm) is converged to 0.015-0.021 mm RMS; the error is in the
  formed shape.
* **The in-plane element size is the parameter that matters** (observed
  order ~2 between 2 and 0.625 mm, but 1.3-6 depending on the mesh triple:
  the convergence is not cleanly asymptotic, so every error below is quoted
  as a band over 6-12 Richardson estimates, `richardson_fit*.csv`). Layers
  (2 -> 4 standard, 1 -> 2 IM), points through the thickness (5 -> 7), tool
  travel (1 -> 0.5 mm), penalty (10 -> 30), Newton tolerance (1e-6 -> 1e-8)
  and release increments (10 -> 40) each change the released shape by
  <= 0.008 mm RMS, <= 0.027 mm max; penalty 3 instead of 10 by 0.005-0.009 mm
  RMS / 0.014-0.021 mm max. The 2 mm standard Hex8 is close only through a
  cancellation between bending locking and the coarse mesh (the locking-free
  element at 2 mm is *worse*, 0.11-0.15 mm).
* **The rim sag is physical in this model**: refinement reduces it by
  8-12 % (to ~1.02-1.03 mm on both parts), so the benchmark's conclusion
  that the rim band dominates survives.
* **Cheapest setup below 0.05 mm on the final shape:** 0.833 mm in plane
  (48 x 48), one incompatible-mode (IM) layer, 2 x 2 x 5 points, penalty 10:
  released-shape field error 0.021-0.033 mm RMS, error of the RMS deviation
  +0.011..+0.038, of the rim sag -0.016..-0.055, of the worst-point deviation
  +0.011..+0.080 mm (both parts; central estimates +0.015-0.023 / -0.027 to
  -0.031 / +0.031-0.036). 1 586-1 965 s per run on one thread
  (7-11 times the benchmark). It meets 0.05 mm on the RMS and the field
  with margin, on the rim sag and the worst point only on the central
  estimates; 0.625 mm (3 924-4 441 s, 20-22 times) is needed to bound the
  worst point near 0.05 mm with the pessimistic estimates. Penalty 3 halves
  the cost (798-1 172 s) for +0.005-0.009 mm.
* **Explicit is not the route here:** 0.017-0.019 mm from implicit on the
  final shape but 13-20 % less springback (the deficit grows with
  refinement), 1.5-2.8 times the implicit cost on these sheets, and no fast
  kernel for the incompatible-mode element.
