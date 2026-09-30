# Verification and validation

Two different questions, kept separate throughout:

**Verification** - *does the code solve the equations it claims to solve?*
Judged against exact answers for the discrete problem. These must pass to tight
tolerances; a failure is a bug.

**Validation** - *are the modelling assumptions appropriate?* Judged against an
independent theory whose assumptions differ from the model's. A finite gap is
**expected**, and the useful result is its size and sign, not agreement. Where
the gap has a known cause, the test asserts the *sign* of the discrepancy, which
is a stronger statement than asserting agreement.

Reproduce everything below with:

```bash
make test              # the Catch2 suite: 293 cases, 21 454 assertions
make verify            # the studies, which exit non-zero if any tolerance is missed
make cross-validation  # the same problems in CalculiX and scikit-fem, node by node
```

All numbers in this document come from `results/verification/summary.json`,
`results/cross_validation/summary.json` and the test run;
`docs/results/README.md` carries the machine-generated tables.

## Study results

| Study | Kind | Metric | Value | Tolerance | Result |
|-------|------|--------|-------|-----------|--------|
| Patch test, distorted mesh | verification | max relative error in `u`, strain and stress | `4.07e-15` | `1e-10` | PASS |
| Solver agreement (7 solvers, multigrid and `auto` included) | verification | max relative displacement difference vs dense LU | `8.97e-12` | `1e-8` | PASS |
| Topology sensitivity | verification | min over steps of the max relative gradient error | `2.18e-08` | `1e-5` | PASS |
| Mesh convergence | verification + validation | relative tip-deflection error vs Timoshenko, finest mesh, `nu = 0` | `4.78e-04` | `0.02` | PASS |
| Modal frequencies | validation | `f1` relative error vs Euler-Bernoulli, finest mesh | `5.52e-03` | `0.02` | PASS |
| Patch test 3-D, distorted Hex8 mesh | verification | max relative error in `u`, strain and stress | `1.63e-15` | `1e-10` | PASS |
| Topology sensitivity 3-D, Hex8 | verification | min over steps of the max relative gradient error | `1.56e-08` | `1e-5` | PASS |
| Mesh convergence 3-D, Hex8 cantilever | verification + validation | relative tip-deflection error vs Timoshenko, finest mesh | `5.42e-03` | `0.03` | PASS |
| Modal frequencies 3-D, Hex8 cantilever | validation | `f1` (weak axis) relative error vs Euler-Bernoulli, finest mesh | `1.06e-02` | `0.03` | PASS |
| Patch test, distorted Tri3 and Tet4 meshes | verification | max relative error in `u`, strain and stress | `9.17e-15` | `1e-10` | PASS |
| Mesh convergence, Tri3 and Tet4 cantilevers | verification + validation | larger of the two finest-mesh errors vs Timoshenko | `1.66e-02` | `0.03` | PASS |
| Multigrid CG vs Cholesky, Hex8 and Tet4 | verification | max relative displacement difference vs LDL^T; iteration growth under refinement at most `1.6x` | `4.37e-12` | `1e-8` | PASS |
| Sensitivity through the Heaviside projection, Q4 and Tet4 | verification | worst over `beta = 2, 8, 32` of the best scaled-entry or directional error | `9.65e-08` | `1e-5` | PASS |
| Quadratic patch test (pure bending), distorted Tet10 meshes | verification | max relative error in `u` and element stress | `1.22e-14` | `1e-9` | PASS |
| Mesh convergence, Tet10 cantilever | verification + validation | relative tip-deflection error vs Timoshenko, finest mesh | `6.36e-06` | `0.01` | PASS |
| Linear buckling of a clamped column, Q4 and Tet10 | verification + validation | larger finest-mesh error of `lambda_1` vs Euler-Engesser | `6.64e-03` | `0.01` | PASS |
| Buckling-constraint sensitivity, Q4 and Hex8 | verification | worst case of the best-step max scaled error (KS aggregate and `lambda_1`) | `1.72e-06` | `1e-5` | PASS |
| Sensitivity through the overhang filter, Q4 and Hex8 | verification | worst over `beta = 0, 4, 16` of the best scaled-entry or directional error | `2.81e-06` | `1e-5` | PASS |
| Thick cylinder under internal pressure vs Lame (plane strain; Q4, Tri3, Hex8, Tet10) | verification | largest shortfall of the RMS displacement order below the element's | `-0.001` | `0.3` | PASS |
| Rotating disk (plane stress) and cylinder (plane strain) vs exact | verification | largest shortfall of the RMS displacement order | `-0.006` | `0.3` | PASS |
| Conduction and thermal stress in a thick cylinder vs exact | verification | largest shortfall of the RMS temperature and displacement orders | `0.005` | `0.3` | PASS |
| Bimetallic strip under uniform heating vs Timoshenko | verification | finest-mesh relative curvature error (order `2.00` also required) | `3.01e-04` | `1e-3` | PASS |
| Bar hanging under its own weight vs exact | verification | Tet10 displacement error (its space holds the exact field); linear-element RMS orders also required | `3.78e-13` | `1e-10` | PASS |
| Large-deflection cantilever vs Euler's elastica (Tet10, `k <= 10`) | verification + validation | largest tip-displacement error on the finest mesh, in `L` (observed order `>= 2.7` also required) | `1.36e-04` | `1e-3` | PASS |
| Thick tube at finite strain vs exact: follower pressure, spin, heating (Q4, Tri3, Hex8, Tet10) | verification | largest shortfall of the RMS displacement order below the element's | `9.0e-04` | `0.3` | PASS |
| Shallow-arch snap-through: arc length vs displacement control | verification | largest crown-force difference at equal deflection / limit force (inertia and the load-control bracket also required) | `2.82e-11` | `1e-6` | PASS |
| Thick tube to plastic collapse vs the exact limit load (Q4, Hex8, Tet10; mean dilatation and locking) | verification | largest collapse-pressure error of the default elements, finest mesh (order `>= 1.8` and the fully plastic stress field also required) | `4.82e-04` | `1e-3` | PASS |
| Pure bending: moment-curvature and residual stress vs exact (plane-stress Q4) | verification | largest moment error over `M_p` on the loading branch, finest mesh (order, residual moment and residual stress also required) | `6.17e-04` | `1e-3` | PASS |
| Uniaxial cycle with combined hardening vs exact (distorted Hex8) | verification | largest stress error over `sigma_y` along the cycle | `1.98e-14` | `1e-9` | PASS |
| Uniaxial cycle with Chaboche backstresses vs exact (distorted Hex8, Q4) | verification | largest stress error over `sigma_y` along the cycle | `3.86e-15` | `1e-9` | PASS |
| Hill48 directional yield stress and r-value vs exact (distorted Q4, Hex8) | verification | largest relative error of the stress, the r-value and the homogeneous field, 0 to 90 deg | `5.63e-15` | `1e-9` | PASS |
| Large-strain uniaxial tension and cycle in the logarithmic strain vs exact (distorted Hex8) | verification | largest Kirchhoff-stress error over `sigma_y` and relative r-value error, to a stretch of 2 | `1.15e-13` | `1e-12` | PASS |
| Thick tube to collapse at finite strain, logarithmic strain (Q4, Hex8 with mean dilatation of `ln J`; Tet10) | verification | peak-pressure difference of the averaged Q4 and Hex8 on the finest mesh to the finest Tet10 (order `>= 1.8` also required) | `4.73e-04` | `1e-3` | PASS |
| HHT-alpha transient vs the exact discrete modal solution (Q4, Hex8; consistent and lumped mass) | verification | largest relative displacement difference over the steps, models, masses and cases (trapezoidal energy balance `<= 1e-10` and positive numerical dissipation also required) | `3.28e-10` | `1e-9` | PASS |
| Harmonic response of a rod vs the exact discrete and continuum solutions (Q4, Hex8) | verification | largest relative difference to the exact discrete solution (continuum order `>= 1.9` on the finest pair also required) | `2.94e-10` | `1e-9` | PASS |
| Transient of a rod under a ramped end force vs the exact continuum solution (Q4, Hex8) | verification | smallest observed convergence order, `h` and `dt` halved together | `2.004` | `>= 1.9` | PASS |
| Non-linear oscillators, finite-strain elastic and elastoplastic, vs exact motion (Q4, Hex8) | verification | largest relative difference to the scalar HHT-alpha recursion (order `>= 1.8` to the exact motion also required) | `1.74e-11` | `1e-9` | PASS |
| Slender cantilever, one incompatible-mode Hex8 through the depth, vs Timoshenko | validation | `|tip / Timoshenko - 1|` on 20 x 1 x 1, `nu = 0` (every `nx >= 20`, `nu` and `L/h` up to 1000 within 5 % also required) | `1.62e-03` | `0.02` | PASS |
| MacNeal-Harder straight beam, incompatible-mode Hex8 | validation | worst `|tip / reference - 1|`: rectangular cells, both loads; parallelogram cells, in plane | `2.71e-02` | `0.03` | PASS |
| Thin square plate, `t/a = 1/50`, two incompatible-mode Hex8 layers, vs Kirchhoff | validation | worst `|centre / Kirchhoff - 1|` on 16 x 16 x 2, clamped and simply supported | `9.25e-03` | `0.03` | PASS |
| Elastoplastic bending and springback, one incompatible-mode layer, 2 x 2 x 5 rule, vs a converged fine standard mesh | verification | relative springback error (2 x 2 x 7 within 2 %, the reference within 1 % and half the standard one-layer error also required) | `3.20e-02` | `0.05` | PASS |

Supporting measurements from the same runs:

| Quantity | Value |
|----------|-------|
| Observed convergence order, tip deflection | `2.31` at `nu = 0`, `2.28` at `nu = 0.3` (Q4); `2.63` (Hex8); `2.27` (Tri3); `3.22` (Tet4, still pre-asymptotic - section 15) |
| Multigrid CG iterations, coarsest to finest multi-level mesh | `14 -> 16` (Hex8), `16 -> 20` (Tet4); Jacobi CG `238 -> 699` and `481 -> 729` on the same meshes |
| Mass conservation, `sum(M)/dim` vs `rho V` | `<= 4.84e-14` relative (Q4), `<= 4.89e-14` (Hex8) |
| Axial mode vs fixed-free rod theory | `4.02e-06` relative |
| Elements excluded from the FD check at active bounds | 6 of 72 (Q4), 2 of 36 (Hex8) |
| Tet10 cantilever error vs Timoshenko, coarsest to finest grid | `6.1e-04 -> 4.7e-05 -> 6.4e-06`; Hex8 3.0 % and Tet4 11 % on the finest grid |
| Column `lambda_1` error vs Euler-Engesser, finest mesh | Q4 0.66 %, Tet10 0.29 % (Richardson limit of the Tet10 errors 0.21 %), Hex8 11.1 %, Tet4 39.7 % |
| Elastica: observed Tet10 order, continuum-elastica gap | deflection `3.01 -> 3.49`, shortening `3.06 -> 3.65` from `k = 1` to 10; gap `2.2e-05 L -> 1.46e-04 L` |
| Finite-strain tube: exact bore hoop stretch, reference checks | `1.48077` (inflation), `1.17256` (spin), `1.02502` (heating); equilibrium `<= 3.7e-7`, Lame limit `3.2e-11` |
| Arch: limit force, load-control bracket | `1588.2464 N`; load control stops at `1588.135 N` and rejects `1588.257 N` |
| Plastic tube: collapse-pressure error, finest mesh (order) | Q4 mean dilatation `1.20e-04` (2.00), Hex8 `4.82e-04` (2.00), Tet10 `2.99e-06` (3.02); fully integrated Q4 `1.20e-02` with a rising plateau, Tet10 with mean dilatation `-7.12e-04` (2.05, from below) |
| Plastic bending: moment error order, residual moment, residual stress | `1.94`; `1.63e-04 M_p`; `2.33e-03 sigma_y` (order `1.65`) |
| Transient: trapezoidal energy balance; its fall from double to 80-bit arithmetic; numerical dissipation over the work | `<= 2.98e-11`; `1.71e3` times for `2.05e3` in eps; `0.36 %` to `0.40 %` (`alpha = -0.1`, harmonic load), `36 %` to `49 %` (`alpha = -0.3`, sudden load) |
| Rod: continuum error on 160 elements; transient order and error on 160 elements | harmonic `<= 3.61e-03` (third resonance), order `2.00`; transient order `2.01` / `2.00`, error `2.59e-05` / `8.00e-05` (consistent / lumped) |
| Oscillators: order to the exact motion; plastic dissipation at 320 steps a period | `1.98` to `2.00` (elastic), `1.94` to `2.04` (plastic); `174.794 J` against the exact `174.802 J` |

And from the cross-validation against two independent codes (section 14):

| Problem | Reference | Max relative nodal-displacement difference | Tolerance | Result |
|---------|-----------|-------------------------------------------:|----------:|--------|
| Plane cantilever, 1 440 Q4 | scikit-fem 12.0.2, `ElementQuad1` | `1.50e-10` | `1e-7` | PASS |
| Plane cantilever, 1 440 Q4 | CalculiX 2.21, `CPS4` | `8.75e-07` | `1e-5` | PASS |
| Solid block, tip load, 1 280 Hex8 | scikit-fem, `ElementHex1` | `5.96e-12` | `1e-7` | PASS |
| Solid block, tip load, 1 280 Hex8 | CalculiX, `C3D8` | `3.39e-06` | `1e-5` | PASS |
| Solid block, top pressure | scikit-fem, `ElementHex1` | `1.23e-11` | `1e-7` | PASS |
| Solid block, top pressure | CalculiX, `C3D8` | `2.43e-06` | `1e-5` | PASS |
| Plane cantilever, 2 880 Tri3, `nu = 0` | scikit-fem, `ElementTriP1` | `6.96e-12` | `1e-7` | PASS |
| Plane cantilever, 2 880 Tri3, `nu = 0` | CalculiX, `CPS3` | `8.86e-07` | `1e-5` | PASS |
| Solid block, 7 680 Tet4, tip load / top pressure | scikit-fem, `ElementTetP1` | `1.34e-12` / `7.75e-12` | `1e-7` | PASS |
| Solid block, 7 680 Tet4, tip load / top pressure | CalculiX, `C3D4` | `3.64e-06` / `2.97e-06` | `1e-5` | PASS |
| Gmsh lug bracket, 20 336 Tri3, `nu = 0`, two load cases | scikit-fem, `ElementTriP1` | `5.26e-13` / `8.66e-14` | `1e-7` | PASS |
| Gmsh lug bracket, 20 336 Tri3, `nu = 0`, two load cases | CalculiX, `CPS3` | `2.76e-06` / `4.29e-06` | `1e-5` | PASS |
| Gmsh lug bracket, `nu = 0.33` (the benchmark material) | scikit-fem, `ElementTriP1` | `3.90e-13` / `2.86e-13` | `1e-7` | PASS |
| Gmsh lug bracket, `nu = 0.33` | CalculiX, `CPS3` | `5.63e-04` / `1.14e-03` | - | INFO |
| Gmsh engine mount, 39 936 Tet4, two load cases | scikit-fem, `ElementTetP1` | `1.52e-12` / `5.62e-13` | `1e-7` | PASS |
| Gmsh engine mount, 39 936 Tet4, two load cases | CalculiX, `C3D4` | `1.71e-06` / `1.98e-06` | `1e-5` | PASS |
| Gmsh engine mount, 13 918 curved Tet10, two load cases | scikit-fem, `ElementTetP2` on `MeshTet2` | `1.16e-12` / `7.95e-13` | `1e-7` | PASS |
| Gmsh engine mount, 13 918 curved Tet10, two load cases | CalculiX, `C3D10` | `1.48e-06` / `1.65e-06` | `1e-5` | PASS |
| Axial column: 640 Hex8 / 480 Tet4 / 480 Tet10 | scikit-fem | `2.99e-11` / `6.63e-12` / `5.85e-11` | `1e-7` | PASS |
| Axial column: 640 Hex8 / 480 Tet4 / 480 Tet10 | CalculiX, `C3D8` / `C3D4` / `C3D10` | `1.05e-06` / `1.95e-06` / `2.46e-06` | `1e-5` | PASS |
| Block, 500 Hex8 / 648 Tet10: self-weight, pressure, rotation, body force, combined | scikit-fem integrating the loads itself | `<= 1.32e-12` / `<= 1.79e-12` | `1e-7` | PASS |
| Block, 500 Hex8 / 648 Tet10, the same five cases | CalculiX, `C3D8` / `C3D10`, its own load cards | `<= 3.25e-06` / `<= 3.29e-06` | `1e-5` | PASS |
| Two-material plate, 400 Hex8 / 300 Tet10 / 800 Q4 plane strain: conducted, uniform and regional temperatures | scikit-fem integrating the thermal load itself | `<= 7.64e-13` | `1e-7` | PASS |
| The same plates | CalculiX, `C3D8` / `C3D10` / `CPE4`, `*EXPANSION` and `*TEMPERATURE` | `<= 3.57e-06` | `1e-5` | PASS |
| The same plates, conducted temperature | CalculiX `*HEAT TRANSFER` (temperature, relative to its range) | `3.17e-05` / `3.47e-05` / `1.12e-05` | `1e-4` | PASS |
| Gmsh engine mount, 13 918 curved Tet10: bore pressure, self-weight, rotation, conduction | scikit-fem integrating the loads itself | `<= 8.92e-12` | `1e-7` | PASS |
| The same | CalculiX `C3D10`; conduction `*HEAT TRANSFER` | `<= 2.45e-06`; `5.21e-06` | `1e-5`; `1e-4` | PASS |
| Elastoplastic, small strain: a Hex8 cantilever loaded past yield and unloaded / a plane-strain Q4 strip / a Tet10 punch partly relieved / a Tet10 plate heated past yield and cooled | scikit-fem, an independent J2 solve | `5.00e-13` / `3.90e-15` / `9.94e-13` / `6.85e-12` | `1e-7` | PASS |
| The same four | CalculiX `*PLASTIC`, `C3D8` / `CPE4` / `C3D10` / `C3D10` | `1.76e-06` / `1.94e-06` / `3.15e-06` / `3.25e-06` | `1e-5` | PASS |
| Elastoplastic cycles: Hex8 with mean dilatation and combined hardening / plane-stress Q4 with kinematic and Voce hardening | scikit-fem, an independent J2 solve | `3.24e-14` / `6.45e-15` | `1e-7` | PASS |
| Elastoplastic, finite kinematics: a 960-Tet10 cantilever deflected a tenth of its span / a beam clamped at both ends, driven into membrane action and back with combined hardening, 320 Hex8 with E-bar / 160 Q4 plane strain with E-bar / 160 Q4 plane stress | scikit-fem, an independent finite-kinematics J2 solve | `1.58e-14` / `3.42e-13` / `3.37e-13` / `2.28e-13` | `1e-7` | PASS |
| Elastoplastic large deflection, 960 Tet10 | CalculiX `*PLASTIC` under `NLGEOM` (a different finite-strain model) | `1.81e-04` | - | INFO |
| Transient: a Hex8 cantilever (HHT `alpha = -0.05`, Rayleigh damping) / a Tet10 cantilever driven harmonically (trapezoidal rule) / a plane-strain Q4 strip shaken at its root (lumped mass) | scikit-fem, an independent HHT-alpha integration | `3.27e-10` / `3.46e-09` / `1.50e-09` | `1e-7` | PASS |
| The Hex8 and Tet10 transients | CalculiX `*DYNAMIC`, `C3D8` / `C3D10` | `1.64e-06` / `2.08e-06` | `1e-5` | PASS |
| Non-linear transient: a Hex8 cantilever loaded past yield (small-strain J2) / a slender Hex8 strip swinging through large deflection (Saint Venant-Kirchhoff) | scikit-fem, an independent non-linear HHT-alpha integration | `2.17e-08` / `4.90e-10` | `1e-7` | PASS |
| The same two | CalculiX `*DYNAMIC`, `C3D8`, with `*PLASTIC` / under `NLGEOM` | `6.21e-07` / `1.46e-06` | `1e-5` | PASS |
| Harmonic response: a Q4 cantilever plate through four resonances / a Hex8 block on a shaken base (lumped mass) | scikit-fem, a direct complex solve | `1.02e-10` / `8.15e-13` | `1e-7` | PASS |
| The Q4 plate's static solution (plane stress, `nu = 0.3`) | CalculiX, `CPS4` | `1.52e-03` | - | INFO |

And the linear buckling load factors of the same three columns, four modes
each (section 20):

| Column | Reference | Max relative load-factor difference | Tolerance | Result |
|--------|-----------|------------------------------------:|----------:|--------|
| 640 Hex8 | scikit-fem: `K_G` of its own static solution, dense eigensolve | `7.22e-10` | `1e-7` | PASS |
| 640 Hex8 | CalculiX `*BUCKLE`, `C3D8` | `8.26e-05` | `1e-4` | PASS |
| 480 Tet4 | scikit-fem | `8.54e-11` | `1e-7` | PASS |
| 480 Tet4 | CalculiX `*BUCKLE`, `C3D4` | `6.43e-06` | `1e-4` | PASS |
| 480 Tet10 | scikit-fem | `8.41e-10` | `1e-7` | PASS |
| 480 Tet10 | CalculiX `*BUCKLE`, `C3D10` | `6.70e-05` | `1e-4` | PASS |

The `INFO` rows are not a disagreement between codes but between
idealisations or models: CalculiX expands its plane elements into a layer of solid
elements, which reproduces plane stress only at `nu = 0` - the same mesh at
`nu = 0` agrees to the `.frd` rounding floor (section 14); the elastoplastic
large-deflection row sets SparLab's J2 return in the Green-Lagrange strain
against CalculiX's finite-strain plasticity (section 24). Where CalculiX's own
formulation of a load differs from SparLab's, its row is judged against
scikit-fem solving CalculiX's problem (section 22); SparLab's own loads are
judged by scikit-fem integrating them independently, to `1e-12`.

## 1. Element-level verification

Every item is a unit test in `tests/test_element.cpp`.

| Property | How it is checked | Result |
|----------|-------------------|--------|
| `K_e` symmetry | `max\|K - K^T\| / max\|K\|` on a uniform and a distorted element | `< 1e-15` |
| Rigid-body modes | `K_e r_i` and `r_i^T K_e r_i` for the two translations and the rotation | `< 1e-14` relative |
| Null-space dimension | eigenvalues of `K_e`: exactly three below `1e-9 max\|K\|`, the rest positive | exact |
| Quadrature sufficiency | 2x2, 3x3 and 4x4 stiffness on a rectangle | agree to `1e-13` |
| Linear-field reproduction | `B u_e` for `u = a + G x` at five parametric points including the corners, on a distorted element | `< 1e-17` absolute |
| Jacobian integration | `sum_g w_g detJ` vs the shoelace area of a distorted quadrilateral | `1e-12` relative |
| Thickness / modulus scaling | `K_e(3E) = 3 K_e(E)`, `K_e(3t) = 3 K_e(t)` | `< 1e-14` |
| `M_e` symmetry and definiteness | symmetry plus the smallest eigenvalue | symmetric to `1e-15`, all eigenvalues `> 0` |
| Mass conservation | `sum(M_e) = 2 rho t A`, and `t_x^T M_e t_x = rho t A` | `1e-12` relative |
| Lumping conserves mass | row sums of `M_e` add to `sum(M_e)` | `1e-14` |
| Consistent edge load | total equals `traction x length x thickness`, split evenly on a straight edge | exact |
| Inverted / collapsed elements | clockwise and zero-height elements | `MeshError` raised |

## 2. Assembly, partitioning and solvers

`tests/test_assembly_solver.cpp`.

* **Global symmetry and rigid-body modes.** The assembled `K` is symmetric to
  `1e-14` relative and annihilates all three rigid-body vectors to `1e-13`. An
  unconstrained `K` is correctly reported as *not* positive definite.
* **Element-cache equivalence.** The same uniform mesh assembled through the
  cached path and through the generic per-element path gives matrices agreeing
  to `1e-15` relative, for both `K` and `M`. The optimisation loop therefore
  cannot silently diverge from the general case.
* **Linear scaling of the SIMP factor.** `assemble_stiffness(0.25)` equals
  `0.25 K` to `1e-15`.
* **Positive definiteness after valid constraints.** `K_ff` on a properly
  clamped plate: every LDL^T pivot positive, and the smallest eigenvalue of the
  dense equivalent strictly positive.
* **DOF round-trip.** `expand(restrict_to_free(u))` reproduces `u` at free DOFs
  and the prescribed values at constrained ones.
* **Solver agreement.** The direct backends - `SimplicialLDLT`,
  `SimplicialLLT`, `SparseLU`, dense `PartialPivLU` - and the iterative ones -
  diagonally preconditioned CG, multigrid CG, and `auto` - agree on one 6x4
  plate within `1e-8` relative on the displacement field and `1e-9` on the
  compliance. The dedicated study on a 12x4 plate (with the multigrid coarse
  size lowered to 16 unknowns, so it builds a real hierarchy on 130 DOFs)
  measures at most `8.97e-12`, from the multigrid solver in 41 iterations.
* **Cached-pattern assembly.** The assembly that scatters into the recorded
  sparsity pattern is bitwise identical to the triplet assembly for `K` and
  `M`, on Q4, Tri3, Hex8 and Tet4 meshes, with and without SIMP scale
  factors (`tests/test_multigrid.cpp`).

### Reaction-force equilibrium

On a 24x6 cantilever with a 1 kN tip resultant:

| Check | Measured |
|-------|----------|
| Vertical reaction vs applied load | `1e-10` relative |
| Horizontal reaction (should vanish) | `< 1e-8` N |
| Global force-balance residual | `< 1e-10` relative |
| Root moment vs `P L` | `1e-10` relative |
| Moment-balance residual | `< 1e-10` relative |
| Reactions at free DOFs (should vanish) | `0` |
| `C = 2U` | `1e-10` relative |
| Scaled solve residual | `< 1e-10` |

With a **non-zero prescribed displacement** (a `1e-4` m enforced stretch on a
`nu = 0` bar), the reaction resultant matches the closed form `E A eps` to `1e-9`
relative, the two edge reactions cancel to `1e-6` N, the compliance is exactly
zero (the applied load vector is zero) and the strain energy matches
`1/2 F delta` to `1e-9`. That last pair is the reason compliance and strain
energy are reported separately rather than assuming `C = 2U`.

## 3. Detection of ill-posed models

`ModelDiagnostics` runs before any factorisation. `tests/test_assembly_solver.cpp`
covers each case:

| Model | Detected | Reported null dimension |
|-------|----------|-------------------------|
| No constraints at all | yes | 3 |
| One node pinned in x and y | yes | 1 (rotation) |
| Bottom edge fixed in y only | yes | 1 (x translation) |
| Clamped edge | well posed | 0 |
| Pin plus a roller at another node | well posed | 0 |
| Two blocks, only one constrained | yes | the floating group is named |

`StaticAnalysis` refuses to construct on an ill-posed model, so the failure
appears before any solve. Where the analytic check cannot see the problem - an
internal mechanism, two blocks joined at a single node - it surfaces as a
non-positive Cholesky pivot, which `LinearSolver` reports as a singular system
with the three usual modelling causes named.

Explicitly exercised solver failures: a zero matrix, an indefinite matrix, a
matrix containing NaN, an empty system, and the dense solver's size refusal.
Each raises `SolverError` rather than returning a plausible-looking answer.

## 4. Patch test

The strongest single verification here. A constant-strain field

```
  u = (1e-4, -2e-4) + [[3e-4, 1e-4], [1e-4, -2e-4]] x
```

is prescribed on every boundary node of a 4x3 mesh, and the interior solution,
the recovered strain and the recovered stress are compared with the exact field.
Run at four interior-node perturbations, so the isoparametric mapping is
genuinely exercised: on a *uniform* grid a Q4 reproduces linear fields for
trivial reasons.

| Perturbation (fraction of cell size) | Max relative error in `u`, strain and stress |
|--------------------------------------|----------------------------------------------|
| 0.00 | `< 1e-11` |
| 0.15 | `< 1e-11` |
| 0.30 | `< 1e-11` |
| 0.40 | `< 1e-11` |

The dedicated study measures `4.07e-15`, i.e. round-off. This case is also why
`prescribed_displacement_only` exists on a load case: the test has no applied
force at all, and declaring that explicitly keeps an accidentally empty load
case an error.

## 5. Stress recovery against a known field

A linear displacement field is imposed directly on a distorted 3x2 mesh,
bypassing the solve, and every recovered quantity is compared with its closed
form:

| Quantity | Agreement |
|----------|-----------|
| Element strain | `< 1e-12` relative |
| Element stress | `< 1e-12` relative |
| von Mises | `1e-11` relative |
| Principal-stress sum vs `sxx + syy` | `1e-11` relative |
| Nodal averaged stress (constant field) | `< 1e-11` relative |
| Total strain energy vs `1/2 eps^T sigma V` | `1e-10` relative |
| SIMP scaling: macroscopic stress scales, solid stress and strain do not | `1e-10` |

Separately, the element field is reproduced from independent pointwise
`element_stress_at` calls to `1e-12`, confirming the documented averaging rule.

von Mises is checked against four closed forms: uniaxial tension gives the axial
stress, pure shear gives `sqrt(3) tau`, equal biaxial tension gives the same
value in plane stress, and plane strain gives the smaller value that including
`szz = nu(sxx + syy)` implies. The principal stresses are checked against both
invariants of the 2-D stress tensor.

## 6. Mesh convergence and beam theory

Cantilever: `L = 1 m`, `h = 0.1 m` (`L/h = 10`), `t = 10 mm`, `E = 70 GPa`, 1 kN
tip resultant spread over the tip edge. Six meshes from 10x2 to 320x64, at two
Poisson ratios.

References:

```
  Euler-Bernoulli:  delta = P L^3 / (3 E I)                        (bending only)
  Timoshenko:       delta = P L^3 / (3 E I) + P L / (k G A),  k = 5/6
```

**Verification part** - self-convergence against the finest mesh isolates the
discretisation error from the beam-theory modelling gap. The observed order is
`2.31` at `nu = 0` and `2.28` at `nu = 0.3`, consistent with the second-order
convergence expected of the Q4 displacement, and the error shrinks monotonically
under refinement (asserted in the test suite).

**Validation part** - on the finest mesh at `nu = 0`:

| Reference | Relative difference |
|-----------|--------------------|
| Timoshenko (bending + shear) | `4.78e-04` |
| Euler-Bernoulli (bending only) | about `6e-03` |

The FEM deflection is larger than Euler-Bernoulli in magnitude, which is the
right sign: beam theory omits the shear deformation the 2-D model includes. At
`nu = 0.3` a further gap appears, from the anticlastic (Poisson) coupling that
1-D beam theory also omits - a **modelling** difference, not a code error, which
is why the studies run both Poisson ratios side by side.

Two more properties are asserted rather than merely measured:

* compliance converges **from below** under refinement, because the
  displacement-based Q4 is a lower bound on the true compliance;
* the distributed-traction resultant is mesh independent to `1e-12` at every
  resolution, and the compliance changes by only a few percent between 8x2 and
  32x8.

## 7. Modal verification and validation

**Mass matrix** (`tests/test_modal.cpp`): symmetric to `1e-15`, positive
definite on the free DOFs, and `sum(M) = 2 rho V` to `1e-11` for both the
consistent and the lumped form, on a *distorted* mesh. The dedicated study
measures `4.84e-14`. A rigid translation carries exactly the total mass. A zero
density raises `ModelError` rather than producing a singular pair.

**Eigen solver.** On a 40x8 cantilever (698 free DOFs, above the 400-DOF
threshold, so the iterative path runs) the subspace iteration matches a dense
`GeneralizedSelfAdjointEigenSolver` to `1e-8` relative on all five requested
eigenvalues. Eigenvalues come out ascending, all eigenpair residuals are below
`1e-8`, the modes satisfy `phi^T M phi = 1` to `1e-8`, and every mode vanishes
at prescribed DOFs.

**Thickness invariance.** Scaling the thickness by 0.37 leaves all three
frequencies unchanged to `1e-9` while the mass scales by exactly 0.37. This is
the property that makes "equal-mass uniform plate" a well-defined baseline in
the topology study, and it is measured rather than assumed - the topology runs
repeat the check and record the maximum difference.

**Lumped vs consistent.** Row-sum lumping gives strictly lower frequencies (mass
moved toward the diagonal) and agrees with the consistent matrix within 5% on a
converged mesh. Both the sign and the size are asserted.

**Bending validation.** Cantilever bending frequencies at `nu = 0` against
Euler-Bernoulli:

| Mode | FEM (120x12) | Theory | Difference |
|------|--------------|--------|-----------|
| 1 | within 1% | `beta_1 L = 1.87510407` | below theory |
| 2 | within 6% | `beta_2 L = 4.69409113` | below theory |

Both are asserted to be *below* theory, since the 2-D model includes shear
deformation and rotary inertia that Euler-Bernoulli omits, and the gap must grow
with mode number. The dedicated study measures `5.52e-03` for `f1` on the finest
mesh.

**Axial validation.** A plane-stress cantilever also carries axial modes,
interleaved with the bending series. Against fixed-free rod theory
`f_n = (2n-1)/(4L) sqrt(E/rho)` - which is *exact* for a bar in uniaxial stress,
and plane stress at `nu = 0` is exactly that - the measured agreement is
`4.02e-06` relative. The mode is identified by its motion being more than 90%
axial, not by its index, since where it falls in the spectrum depends on the
aspect ratio. In the `cantilever_analysis` deck it appears as mode 3 at
1272.9 Hz against a theoretical 1272.9 Hz.

**Rigid-body screening.** An unconstrained plate either produces warnings
naming the near-zero eigenvalues, with `lambda_0 < 1e-6 lambda_3`, or raises
`SolverError`; it never returns rigid-body modes silently as physical
frequencies. Requesting more modes than there are free DOFs warns and truncates.

## 8. Topology sensitivity against finite differences

The required regression test, run both as a study
(`sparlab_verify --study sensitivity`) and as a test case
(`tests/test_topopt.cpp`).

**Setup.** A 12x6 cantilever (72 elements), density filter of radius 1.5 cells,
`p = 3`, `emin_ratio = 1e-9`, direct solver with a `1e-9` residual tolerance. The
design point is deliberately non-uniform,
`x = 0.35 + 0.30 sin(7x) cos(5y)` on free elements: a uniform field makes every
sensitivity nearly identical and would hide an indexing mistake. One passive
solid patch and one passive void disc are included so the exclusion logic is
exercised. **Every** element is tested, at eight step sizes.

**Comparison.** `dc/dx_e` from the analytical chain rule against
`[c(x + h e_e) - c(x - h e_e)] / (2h)`, where `c` is the *filtered* SIMP
compliance - i.e. the whole objective, filter included.

| Step `h` | Max relative error | Behaviour |
|----------|--------------------|-----------|
| `1e-1` ... `1e-3` | falls as `h^2` | truncation-error dominated |
| `1e-5` | `2.18e-08` (best) | balance point |
| `1e-7`, `1e-8` | rises again | round-off dominated: differencing two compliances that agree to ~10 digits |

The `h^2` fall and the subsequent rise is the signature of a correct gradient;
a wrong gradient shows a plateau at its own error level instead. Three error
measures are reported at every step: the maximum and RMS over elements, and the
relative error of the directional derivative along the gradient, which is the
most sensitive single scalar.

**Exclusions.** 6 of 72 elements are excluded, each with a recorded reason:

* passive variables (lower bound equals upper bound) - no admissible
  perturbation exists;
* free variables within `h` of 0 or 1 - a central difference would leave the
  feasible box, so the one-sided value is not a valid central difference.

**Pass criterion.** Maximum relative error `<= 1e-5` at the best step. Measured
`2.18e-08`, i.e. three orders of magnitude of margin. The test-suite version
additionally requires the RMS error below `1e-6`, the directional error below
`1e-5`, and asserts that passive elements *were* excluded - so a regression that
quietly stopped excluding them would fail.

The check is also run **without** a filter, where `dc/dx` must equal
`dc/drho` exactly (verified to `1e-14`), and on a **two-load-case** problem with
unequal weights, where the weighted objective's gradient is verified the same
way.

## 9. Topology optimiser verification

`tests/test_topopt.cpp`.

| Property | Check | Result |
|----------|-------|--------|
| SIMP derivative | central differences at five densities | `1e-6` relative |
| SIMP clamping | `rho` outside `[0,1]` | clamped, no NaN from `pow` |
| Mass interpolation | `E/m` at `rho = 1e-3` under the matched law | order 1, and above the linear law's value |
| Filter row sums | every row of `Hhat` | `1` to `1e-12` |
| Filter preserves a constant | `Hhat c = c` | `1e-14` |
| Filter adjoint | `<H x, v> == <x, H^T v>` | `1e-12` |
| Filter range | a 0/1 checkerboard stays in `[0,1]` and its range shrinks | yes |
| Volume constraint per OC step | achieved vs target volume over five steps | `1e-8` relative, `volume_converged` true |
| Bounds and move limit | every variable after each step | respected to `1e-12` |
| Passive regions | solid pinned at 1, void at 0, after a full run | exact |
| Gradient signs | `dc/dx <= 0`, `dv/dx > 0` | yes |
| Adding material lowers compliance | `c(x + 0.1) < c(x)` | yes |
| Strain-energy identity | `sum_e U_e = c/2` for one load case | `1e-10` |
| Multi-load weighting | `c = sum_l w_l c_l` with normalised weights | `1e-12` |
| Continuation | penalty non-decreasing, ends at the target | exact |
| Non-convergence reporting | 3-iteration cap with an unreachable tolerance | `converged = false`, warning names the iteration cap |
| Filter radius sets the length scale | total variation of the density field at 1.2 vs 3.0 cells | larger radius gives a smoother design |

**End-to-end** on a 30x15 cantilever at 40% volume: converges in 123
iterations, the volume fraction equals the target to `1e-6` at *every* recorded
iteration, compliance improves by a factor of 5.9 over the uniform start, the
grey level settles below 0.40, and the compliance history is non-increasing over
the final third of the run to within `1e-3` relative.

**KKT consistency** (reported in `docs/benchmarks.md` rather than asserted): at
the converged MBB design, `B_e = -(dc/dx_e)/(dv/dx_e)` across elements strictly
inside their bounds has a 5-95 percentile spread of about 9%, which is the
residual expected of OC with a move limit and a finite change tolerance.

## 10. Configuration and I/O

`tests/test_io.cpp` - 30+ cases. Worth naming:

* the JSON parser handles the full value grammar, comments, trailing commas and
  `\u` escapes, and reports the **line and column** of a syntax error (checked
  on a multi-line document);
* duplicate keys, unterminated strings and comments, trailing content and
  malformed numbers are all errors;
* serialisation round-trips, including a non-finite number, which becomes
  `null` rather than invalid JSON;
* `ConfigNode` type errors name the dotted path (`a.c must be a number, got
  string`);
* **unknown keys are reported.** A deck with `"thicknes"` instead of
  `"thickness"` is an error under `--strict-config` and a warning otherwise -
  and the test asserts the misspelling really did cause the default to be taken,
  which is exactly the failure mode the check exists to catch;
* a region selecting no nodes, a region naming two primitives at once, a
  negative load-case weight, modal analysis without a density, and an empty
  load case are all rejected with specific messages;
* the CSV writer enforces its declared column count; the VTK writer's header,
  cell block and data blocks are checked by parsing the file back.

## 11. The solid (Hex8) path

Every study above has a 3-D counterpart, run by `sparlab_verify --study
patch-test-3d | sensitivity-3d | mesh-convergence-3d | modal-3d` and covered
by `tests/test_solid.cpp` (27 cases).

**Element level.** Partition of unity and the nodal delta property of the
trilinear shape functions; the cube Gauss rules as tensor products with total
weight 8; `K_e` symmetric with exactly six zero eigenvalues (three
translations, three rotations `e_j x r`) on a unit cube and on a distorted
cell; the 2x2x2 rule exact on a box (agrees with 3x3x3 to round-off); the
strain operator reproducing a linear field with six independent constant
strains on a distorted cell; `sum_g w_g detJ` against the closed-form volume
of a box and of a sheared box; inverted, degenerate and mis-sized cells
rejected; the consistent mass symmetric, positive definite and summing to
`3 rho V`; a constant face traction split evenly over the four corners with
the exact resultant.

**Patch test 3-D.** `u = (1e-4, -2e-4, 0.5e-4) + G x` with a symmetric `G`
carrying all six constant strain components, prescribed on the boundary of a
3x3x3 mesh whose interior nodes are perturbed; the interior displacement,
strain and stress agree with the exact field to `1.63e-15` - round-off, as in
2-D.

**Diagnostics 3-D.** An unconstrained solid reports a null dimension of six
(and the assembled `K` annihilates all six rigid-body vectors to `1e-13`), a
single fixed node leaves three rotations, a face fixed in one component
leaves the in-face translations and the rotation about the normal, a clamped
face is well posed, and so is a 3-2-1 set of three non-collinear point
supports; a 2-D model rejects any out-of-plane input with a message naming
the key.

**Mesh convergence 3-D.** The solid cantilever `L = 1 m, h = 0.1 m, b = 0.05 m`
at `nu = 0` on four meshes from 16x2x1 to 96x12x6 (306 to 26 481 DOFs), tip
resultant 1 kN over the tip face. The fully integrated Hex8 is stiff in
bending on coarse meshes - the coarsest is 16 % below Timoshenko - and
approaches the reference from below as the section is refined: the finest
mesh is within `5.42e-03` of Timoshenko, and the self-convergence order
against it is `2.63`. The deflection approaching from *below* is the sign
shear locking must have, and the test suite asserts it.

**Modal 3-D.** The same bar at `rho = 2700 kg/m^3` has two first bending
modes, about the weak (0.05 m) and the strong (0.1 m) axis. Both are compared
with Euler-Bernoulli on three meshes: the strong-axis mode is within 2.1 %
on the coarsest mesh, the weak-axis one needs the finest mesh to reach
`1.06e-02` because it bends through only a few elements. Mass is conserved to
`4.89e-14` (`sum(M) = 3 rho V`).

**Sensitivity 3-D.** The compliance gradient of a 6x3x2 Hex8 cantilever with
a density filter, one passive solid and one passive void cell, at a
non-uniform design point, against central differences at eight steps: best
step `1e-5`, maximum relative error `1.56e-08`, the same `h^2` fall and
round-off rise as in 2-D.

**Solver, faces and I/O.** Face tractions give a mesh-independent resultant
and the reactions balance it to `1e-10`; 3-D selectors (box with z limits,
sphere, cylinder about any axis), sub-meshes and face connectivity are
exercised; a `structured_hex` deck parses into a 3-D model; dimension
mismatches in a deck are rejected with a reason; the writers produce z
columns, six stress components and VTK hexahedra.

**End to end.** A short 3-D optimisation (OC) runs, converges on its
criteria, thresholds to a single connected structure and re-solves it.

## 12. MMA and the stress constraint

`tests/test_mma.cpp` (12 cases):

| Property | Check | Result |
|----------|-------|--------|
| Separable problem with a known optimum | `min sum x_j^2 s.t. sum x_j >= 1` converges to `x_j = 1/n`, multiplier `2/n` | `1e-6` |
| Two active constraints | optimum `(0.3, 0.7)` with multipliers `(0.6, 0.8)` | `1e-5` |
| Badly scaled constraint | a constraint violated by `1e5` is scaled to the cap and its multiplier scaled back | subproblem converges; `lambda ~ 1e-5` |
| Malformed input | wrong sizes, non-finite entries, bounds crossed | rejected |
| von Mises as a quadratic form | `sigma^T V sigma` against `von_mises()` in plane stress, plane strain and 3-D | `1e-12` |
| p-norm aggregate | bounded by the maximum from below and by `n^(1/P)` times it from above; centre stress on a uniform mesh equals the recovered element average; the adaptive scale maps it to the maximum | exact / `1e-12` |
| Stress-constraint gradient, 2-D | adjoint gradient of the aggregate, two load cases with unequal weights and a passive pad, vs central differences | relative error `< 1e-5` per entry (scale-aware) |
| Stress-constraint gradient, 3-D | the same on a Hex8 cantilever | `< 1e-5` |
| No filter | `dg/dx == dg/drho` | exact |
| Sensitivity filter | refused with a `ConfigError` | yes |
| MMA vs OC | the same 24x12 compliance problem | compliances within 5 % |
| Stress-limited L-bracket | 32x32, limit at 70 % of the unconstrained design's relaxed peak, which the test locates at the re-entrant corner and not at the load pad | feasible; relaxed ratio within 5 % of 1; the peak at least a fifth below the unconstrained one; compliance higher |
| Infeasible limit | limit `1e3` Pa, far below anything reachable | `feasible = false`, warning; the run is reported, not hidden |
| Deck parsing | every `optimizer.mma.*` and `topology.stress.*` key, plus the two incompatibilities | round-trip and rejection |

## 13. Geometry export

`tests/test_geometry.cpp`: a plane mesh extrudes to a closed surface of
volume `A t` and area `2A + P t` to `1e-12`; a hex mesh's boundary is closed,
outward (checked against the centroid) and encloses exactly its volume; a
thresholded L of cells on a distorted mesh is watertight and agrees with its
cell volume to the flat-triangle error of its bilinear cut faces (`1.5e-4`,
and *not* to round-off - the test asserts both); on the uniform grid the same
cut is exact to `1e-12`; dropping one triangle is reported as three unmatched
edges; two cells touching along an edge are reported as closed with one
non-manifold edge; a binary STL round-trips through `read_stl` with the
standard 84-byte header layout.

## 14. Cross-validation against independent codes

The one item earlier versions of this document listed as absent. The same
discrete problem - nodes, connectivity, supports, consistent nodal loads - is
solved by two independent codes and the nodal displacements compared node by
node (`python/scripts/cross_validate.py`, `make cross-validation`):

* **scikit-fem 12.0.2** rebuilds the problem with `ElementQuad1`,
  `ElementTriP1`, `ElementHex1`, `ElementTetP1` or, for quadratic
  tetrahedra, `ElementTetP2` on an isoparametric `MeshTet2` built from all
  ten nodes of every cell (so curved cells stay curved), with the same Lame
  constants (`lambda* = 2 lambda G/(lambda+2G)` for plane stress) and the
  same integration order - the 4-point rule for the Tet10. It is the *same
  element formulation* in an independent implementation, so the only
  expected difference is linear-solver round-off - and that is what is
  measured: from `2.4e-14` to `1.5e-10` relative on 60 of the 63 load
  cases, and `1.4e-9`, `8.0e-9` and `2.2e-8` on the three worst-conditioned
  systems, a plane-strain strip and the two slender Tet10 cantilevers
  (formulation, section 5). Each comparison records the round-off scale of
  its system, `kappa_1(K) eps` - the 1-norm condition number of the free
  stiffness matrix (Hager and Higham's estimate, equal to the exact value on
  the four decks also checked densely) times machine epsilon, about the most
  two backward-stable solutions of the system can differ by - and every
  difference lies below it: at most `0.47` of it, and the three largest at
  `0.12`, `0.008` and `0.015` of theirs (`1.1e-8`, `1.1e-6`, `1.5e-6`).
  scikit-fem's Tet10 node order is checked against SparLab's cell by cell
  before anything is solved;
* **CalculiX 2.21** (`ccx`) runs the exported `.inp` decks. `C3D8`, `C3D4`
  and `C3D10` are the same trilinear hexahedron and linear and quadratic
  tetrahedra as SparLab's; the differences, `5.2e-07` to `3.6e-06` (round-off
  where the static case is a rigid translation of a shaken base), are
  within the six-significant-digit rounding of its `.frd` result file (floor
  `5e-6`), i.e. as close as the file format allows one to see - including
  the 39 936-tetrahedron engine mount read from a Gmsh file, which SparLab
  solves with the multigrid solver, and the same part on 13 918 curved
  quadratic tetrahedra. Exporting the quadratic part found a writer bug:
  CalculiX reads at most 20 characters per field, and a round-off corner
  load written with 17 digits (`1.4408030432795649e-18`) is 22; the writer
  now trims digits until a field fits (a test checks every field of a Tet10
  deck).
* **CalculiX's plane elements are a different idealisation.** `CPS4` and
  `CPS3` are expanded internally into a layer of solid elements with the
  plane-stress condition imposed on it. That reproduces plane stress only
  for `nu = 0`, and the measurement says so: the Gmsh lug bracket at `nu = 0`
  agrees to `2.8e-06` and `4.3e-06`, within the `.frd` floor, while the same
  mesh at its real `nu = 0.33` differs by `5.6e-04` and `1.1e-03`. The plane
  decks judged against CalculiX therefore run at `nu = 0` (the cantilevers
  and `lug_bracket_nu0_analysis`); the `nu = 0.33` comparison is recorded as
  informational (`passed: null`, `INFO` in the tables), and scikit-fem - the
  same plane element - is the verification for it.

**Buckling load factors.** Where the run computed them (the three column
decks, `sparlab_solve` with a `buckling` section), the lowest four load
factors are compared as well:

* scikit-fem assembles the geometric stiffness `int sigma_ij phi_k,i
  phi_k,j` of *its own* static solution at the same quadrature points and
  solves the pencil `(-K_G) phi = mu K phi` densely - an independent
  implementation of the same discrete problem. The load factors agree to
  `7e-11` to `8e-10`;
* CalculiX runs the exported deck as a `*BUCKLE` step (same mesh, supports
  and nodal loads, which define the load pattern) and prints the factors to
  seven digits. For `C3D4`, whose stress is constant in an element, it
  agrees to `5e-7` to `6.4e-6`; for `C3D8` and `C3D10`, whose stress varies
  within an element, its factors lie `3.6e-5` to `8.3e-5` *above* SparLab's
  and scikit-fem's, consistently across modes, while its static solution of
  the same deck agrees to `1e-6`. Averaging the stress over each element or
  smoothing it through the nodes moves SparLab's value down, not up, so
  neither explains the gap; which detail of CalculiX's stress stiffness does
  has not been identified. The tolerance, `1e-4`, records the measured size;
  it does not explain it, and the verification of the geometric stiffness
  rests on the scikit-fem comparison and the Euler-Engesser study
  (section 20).

**Large-deflection states.** The four non-linear decks are compared at the
full load with scikit-fem's own total Lagrangian solver and CalculiX's
`*STEP, NLGEOM` (section 23).

**Transient and harmonic responses.** The seven dynamic decks are compared
at every step or frequency with an HHT-alpha integration and a direct
complex solve on scikit-fem's own matrices, and the solid-element transients
with CalculiX's `*DYNAMIC` (section 25).

Tolerances are `1e-7` for scikit-fem (displacements, load factors,
non-linear states, transient and harmonic responses) and `1e-5` for CalculiX
displacements, linear, non-linear and transient, `1e-4` for its load
factors, all recorded in the summary with the `.frd` floor. The comparison
exits non-zero if any judged pair exceeds its tolerance, or if a non-linear
run stopped short of its load. Of the 198 comparisons on 38 decks and 63
load cases, 192 pass and 6 are informational. CI runs the scikit-fem half -
displacements, load factors, the dead-load non-linear states, the
elastoplastic states and the transient and harmonic responses - on every
push.

## 15. The linear simplices (Tri3, Tet4)

`tests/test_unstructured.cpp` and `sparlab_verify --study patch-test-simplex |
mesh-convergence-simplex`.

**Element level.** `K_e` symmetric, exactly three (Tri3) or six (Tet4) zero
eigenvalues and the rest positive, a linear field reproduced exactly, the
closed-form consistent mass summing to `dim rho V`, a constant edge or face
traction split equally with the exact resultant, and inverted or flat cells
rejected with a `MeshError`.

**Patch test.** The constant-strain fields of sections 4 and 11 on triangle and tetrahedron
meshes split from perturbed Q4 and Hex8 grids (interior nodes moved by up to
40 % and 30 % of the cell size, and undistorted for comparison):
interior displacement, strain and stress agree with the exact field to
`9.17e-15` - round-off, as for the bilinear and trilinear elements.

**Mesh convergence.** The same two cantilevers as sections 6 and 11, meshed
by splitting the structured cells, so every simplex mesh has exactly the
nodes of a Q4 or Hex8 mesh:

| Element | DOFs | Error vs Timoshenko, coarsest to finest | Q4 / Hex8 on the same beam |
|---------|------|-----------------------------------------|-----------------------------|
| Tri3 (`nu = 0.3`) | 126 to 21 186 | `37.7 %, 13.4 %, 3.87 %, 1.14 %, 0.42 %` | `0.85 %` at 2 754 DOFs, `0.22 %` at 41 730 |
| Tet4 (`nu = 0`) | 306 to 59 211 | `50.2 %, 20.7 %, 6.25 %, 2.90 %, 1.66 %` | `1.21 %` at 8 775 DOFs, `0.54 %` at 26 481 |

The constant-strain elements lock in bending exactly as expected: they
approach the beam value from below, and on the *same nodes* the Tet4 error is
about five times the Hex8 error (`6.25 %` against `1.21 %` at 8 775 DOFs,
`2.90 %` against `0.54 %` at 26 481). The Tri3 needs roughly three times the
unknowns of the Q4 for the same error (interpolating the Q4 curve: 2.6x to
2.8x over the measured range). The order of the error against
Timoshenko rises towards 2 under refinement - `1.28, 1.72, 1.89, 1.95` for
the Tet4 - and the Tri3's last step (`1.44`) is already limited by the
modelling gap between the beam and plane elasticity, the same gap the Q4
curve levels off at. The self-convergence order the study also records
(`2.27` Tri3, `3.22` Tet4) is measured against the finest mesh, which for
the Tet4 is itself still `1.7 %` from the beam value, so it overstates the
rate; the error against the independent reference is the number to read.
The study passes at `1.66e-02` against a tolerance of `0.03`.

## 16. Meshes read from files

`tests/test_unstructured.cpp` (the reader cases) and the two Gmsh parts:

* Gmsh 2.2 and 4.1 files, written by the test itself, round-trip triangle
  and tetrahedron meshes node for node, with their physical groups as node
  and element sets;
* a clockwise triangle and a mirrored tetrahedron are re-ordered and
  counted; unused nodes are dropped and counted; a mesh in millimetres read
  with `scale = 0.001` has the right extent; coincident nodes are reported,
  and merged on request;
* malformed and unsupported files fail with a message naming the cause and,
  where there is one, the Gmsh option that fixes it: mixed triangles and
  quadrilaterals (`Mesh.RecombineAll`), second-order cells
  (`Mesh.ElementOrder = 1`), a binary file or an unsupported version, a
  plane mesh off `z = 0`, a cell referencing a missing node, a file holding
  only boundary elements (`Mesh.SaveAll`), a truncated section, a folded
  cell, a missing file and an unknown extension;
* SparLab's own CalculiX decks read back node for node through the Abaqus
  reader, which also handles `*INCLUDE`, `*NSET` / `*ELSET` with `GENERATE`,
  continuation lines, and refuses `*INSTANCE` translations and multi-part
  assemblies it cannot place;
* a deck addresses a group by name for supports, loads and passive regions,
  resolves a relative mesh path against its own directory, and a group that
  does not exist, or a node set used as an element region, is an error
  listing the sets the mesh has.

The two committed parts (`configs/meshes/`) are then solved and
cross-validated like the structured meshes (section 14): the 20 336-triangle
lug bracket and the 39 936-tetrahedron engine mount agree with scikit-fem to
`1e-12` and better, and with CalculiX to its `.frd` rounding wherever the two
codes solve the same problem.

## 17. The multigrid solver

`tests/test_multigrid.cpp` and `sparlab_verify --study multigrid`.

**The hierarchy's defining identities**, on a distorted Q4 and a Tet4
mesh: a near-null space of dimension 3 and 6; the tentative prolongator
reproduces it exactly (`P_hat B_c = B` to `1e-12` relative) with
orthonormal columns; every coarse operator equals `P^T A P`; the levels
shrink; the operator complexity lies between 1 and 3. On a distorted Hex8
mesh the V-cycle is symmetric (`<M^-1 x, y> = <x, M^-1 y>` to `1e-11`) and
positive definite, and reduces the error, for both smoothers - which is
what plain CG needs of a preconditioner.

**Agreement with the direct solver** on every element type, including a
SIMP-contrast model (densities down to the `1e-9` floor): the displacement
differs from LDL^T by less than the requested tolerance allows. A reused
hierarchy - aggregates kept, numbers refreshed - is bitwise equal to one
built from scratch for the new matrix. Warm starts reach the same answer in
fewer iterations (none, from the exact solution), and one thread and many
give bitwise identical results. An unsupported model is reported as
under-constrained by the coarse factorisation, and malformed input is
rejected.

**Iterations under refinement.** The study solves a cantilever block
(`nx x nx/2 x nx/4` cells, tip load) to a relative residual of `1e-10` with
multigrid CG, Jacobi CG and LDL^T:

| Element | DOFs | Levels | MG iterations | Jacobi iterations | Difference MG vs LDL^T |
|---------|------|-------:|--------------:|------------------:|-----------------------:|
| Hex8 | 2 295 | 2 | 14 | 238 | `7.2e-13` |
| Hex8 | 15 147 | 2 | 14 | 469 | `2.2e-12` |
| Hex8 | 47 775 | 3 | 16 | 699 | `3.7e-12` |
| Tet4 | 6 825 | 2 | 16 | 481 | `2.8e-12` |
| Tet4 | 21 090 | 3 | 20 | 729 | `4.4e-12` |

Jacobi's count grows with `1/h`, as `sqrt(kappa)` must; the multigrid
count barely moves, which is the property the solver exists for. The
smallest meshes of the study (up to 1 092 DOFs) sit below the coarse-grid
size, where the "hierarchy" is the direct coarse solve and converges in one
iteration; the growth limit is judged over the multi-level meshes only
(worst `1.25x`, against a limit of `1.6x`). Timings are in
`results/verification/multigrid_scaling.csv` and the larger comparison in
`docs/benchmarks.md`: for one solve at these sizes Jacobi CG is about as fast
as multigrid, because its cheap iterations cost about what the multigrid
setup does; the direct solver is about 40 times slower at 47 775 Hex8 DOFs
(31.4 s against 0.72 s in the run behind `docs/results`; wall-clock times
vary from run to run on a shared machine, the iteration counts do not).

## 18. The Heaviside projection

`tests/test_projection.cpp` and `sparlab_verify --study
sensitivity-projection`.

* The map keeps `rho_bar(0) = 0` and `rho_bar(1) = 1` for every `beta`
  from 0.5 to 64, is monotone, is within `1e-6` of the identity at
  `beta = 1e-3` and a step at `beta = 256`, and its derivative matches a
  central difference of the map itself.
* The compliance and volume gradients through filter and projection match
  central differences on a distorted Q4, a Hex8 and a Tet4 mesh, and the
  stress-constraint gradient on a Q4 and a Hex8 mesh.
* The `beta` schedule steps on its interval, early on convergence when asked,
  caps at `beta_max`, and a deck switches the projection on with its
  schedule; the sensitivity filter is refused.
* End to end, a projected run finishes at `beta_max` with a much lower grey
  level than the same run unprojected.

The study checks the gradient at `beta = 2, 8, 32` on a 12x6 Q4 and a
5x2x2-cell Tet4 mesh (72 and 120 design variables), at three steps each:

| Element | `beta = 2` | `beta = 8` | `beta = 32` |
|---------|-----------:|-----------:|------------:|
| Q4 | `6.3e-08` | `8.4e-09` | `9.7e-08` |
| Tet4 | `9.4e-09` | `1.2e-08` | `5.0e-08` |

Each entry is the best step's worse of the scaled entry error and the
directional-derivative error. At `beta = 32` the plain relative error of the
worst entry is between `1.37` and `1.76` on the Q4 mesh, depending on the
step: far from the threshold the projection's derivative is tiny, so the
analytical entry is near zero and its central difference is round-off.
Judging each entry against
`max(|analytical|, |FD|, 1e-3 ||gradient||_inf)` separates that from a real
error, and the directional derivative along the full gradient, which no
entry can hide in, agrees to `1e-9`.

## 19. The quadratic tetrahedron (Tet10)

`tests/test_tet10.cpp` and `sparlab_verify --study patch-test-quadratic`,
`--study mesh-convergence-tet10`.

* **Element identities.** The ten shape functions partition unity, are 1 at
  their own node and 0 at the other nine, and their natural gradients match
  central differences; the collapsed Gauss rules integrate every monomial up
  to their stated degree exactly on the reference tetrahedron; the element
  volume from `det J` equals the corner tetrahedron's for a straight-sided
  cell and the exact volume of a cell with a curved edge whose volume is
  known in closed form. The stiffness matrix is symmetric, positive
  semi-definite with exactly six zero eigenvalues (the rigid-body modes),
  and the consistent mass sums to `rho V` per direction; the HRZ lumped mass
  is positive with the corner / edge split `1/36`, `4/27`.
* **Quadratic patch test** (the study). Pure bending `u_x = -k x z, u_y = nu
  k y z, u_z = k/2 (x^2 + nu (z^2 - y^2))` prescribed on the boundary nodes of
  a 3 x 3 x 3 box, interior nodes solved: on undistorted and randomly
  distorted Tet10 meshes (edge nodes at the midpoints of the distorted
  edges) the displacements and the linear element stresses are reproduced to
  `1.22e-14` relative. The linear elements cannot pass it: their smallest
  error on the distorted meshes is `1.9e-3` (Hex8) and `6.2e-3` (Tet4). On
  the undistorted grid the linear elements do hit the nodal values exactly -
  nodal superconvergence of a symmetric grid, not a passed patch test, as
  the Tet4's element stresses (25 % off) show - so the study judges them on
  the distorted meshes only.
* **Mesh convergence** (the study). The solid cantilever of section 11 under
  a consistent tip traction, on the same grids for all three elements: the
  Tet10 error against Timoshenko falls from `6.1e-4` on the coarsest grid
  (10 x 2 x 1 cells) to `6.4e-6`, where the Hex8 is at 3.0 % and the Tet4 at
  11 % - the locking the quadratic element removes.
* **Meshes and files.** Tet4 meshes elevate to Tet10 with shared edge nodes
  and their named sets; boundary faces are keyed by their corners; C3D10
  Abaqus/CalculiX decks and Gmsh type-11 cells (whose last two edge nodes
  are swapped on reading) round-trip node for node; a cell folded by a
  badly curved edge is refused with its minimum Jacobian named; the VTK
  writer emits quadratic tetrahedra (type 24) and the STL writer splits
  every 6-node face into four triangles.

The engine mount on curved Tet10 cells is cross-validated in section 14,
and the Tet4-against-Tet10 comparison on it is `docs/benchmarks.md`,
section 11.

## 20. Linear buckling

`tests/test_buckling.cpp` and `sparlab_verify --study buckling-euler`,
`--study sensitivity-buckling`.

* **Geometric stiffness.** On distorted Q4, Tri3, Hex8, Tet4 and Tet10 cells the element
  `K_G` is symmetric and annihilates rigid translations; a uniform uniaxial
  stress `s0` contracted with a linear transverse mode gives `s0` times the
  element measure; and `phi^T K_G(u) phi = g(phi)^T u`, the identity the
  constraint's adjoint rests on.
* **Eigensolver.** On a column with more than 400 free DOFs the subspace
  iteration agrees with a dense generalised eigensolve of the same matrices
  to `1e-8`; its modes satisfy `phi^T K phi = 1`, residuals below `1e-6` and
  zeros at the supports; halving the load doubles the load factors; a warm
  start from the converged subspace takes fewer iterations. A pure tension
  is reported as having no positive load factor. Under a large tension plus
  a small transverse load the pencil is crowded with negative
  (reversed-load) eigenvalues: the solver switches to the spectral
  transformation with a shift in `(0, lambda_1)` and still matches the
  dense solve to `1e-7`. The Q4 column converges on Engesser's value from
  above under refinement, its second mode near 9 times the first; on one
  20 x 2 x 2 solid mesh the Tet10 is within 1 % while the Hex8 and the Tet4
  lock (43 % and 155 % high), and the square section's two lowest modes are
  equal on the Hex8 and Tet10 meshes.
* **Euler-Engesser** (the study). A 1 m steel column, clamped, under a unit
  axial tip traction, so `lambda_1` is the critical load in newtons, against
  `P = pi^2 E I / (4 L^2)` with Engesser's shear correction:

  | Element | Finest mesh | DOFs | Error of `lambda_1` |
  |---------|-------------|-----:|--------------------:|
  | Q4 (plane stress, 50 x 20 mm section) | 160 x 16 | 5 474 | 0.66 % |
  | Hex8 (50 mm square section) | 40 x 4 x 4 | 3 075 | 11.1 % |
  | Tet4 | 40 x 4 x 4 cells | 3 075 | 39.7 % |
  | Tet10 | 40 x 4 x 4 cells | 19 683 | 0.29 % |

  Every value converges from above, as a displacement-based element must.
  The Tet10 meshes halve `h` exactly, and Richardson extrapolation of their
  three errors (order 1.96) puts the limit at 0.21 %: the gap between the
  solid model and the beam formula, not a discretisation error.
* **Constraint gradient** (the study). The KS-aggregated constraint and
  `lambda_1` against central differences through the density filter and the
  projection, on a 16 x 8 Q4 mesh at `beta = 0` and 4 and a 5 x 2 x 2 Hex8
  mesh: worst best-step scaled error `1.72e-6`. The test suite repeats the
  check on about 25 variables per case, at load factors separated enough
  for the simple-eigenvalue formula to hold.
* **Cross-validation**: load factors against scikit-fem and CalculiX,
  section 14.

## 21. Robust formulation, overhang filter and the manufacturing checks

`tests/test_manufacturing.cpp` and `sparlab_verify --study
sensitivity-overhang`.

* **Overhang filter.** The supports of an element are the three cells below
  it (two at the domain edge) and none on the first layer; solid material
  passes through up to the smooth minimum's offset of `sqrt(epsilon)/2`; a
  bar floating over void is removed, with everything it would have carried;
  a 45-degree staircase is self-supporting built +y and -y; passive elements
  keep their density and support what rests on them. The filter's adjoint
  recursion matches central differences on Q4 and Hex8 meshes in four build
  directions (worst `1.04e-6`), and the study repeats it through the
  projection at `beta = 0, 4, 16` (worst `2.81e-6`). Near-void densities (`1e-8`, `1e-300`, 0) keep the output and
  the gradient finite: an MBB run met exactly this - a direct sum of
  `1e-320` terms whose derivative overflowed - and the smooth maximum is now
  evaluated relative to the largest support. Where nothing underflows, the
  result equals the textbook formula to `1e-13`.
* **Overhang check.** On the floating bar it counts the bar's six elements,
  not the row resting on it, and reports the lowest unsupported layer; on
  the staircase it counts none. The figure script's independent count must
  equal the run's.
* **Robust projection.** At every density the eroded design is below the
  blueprint and the blueprint below the dilated design; the compliance
  gradient of the eroded design and the volume gradient of the dilated one
  match central differences through the overhang filter; a robust run
  reports the eroded, blueprint and dilated designs with thinner parts
  softer; the erosion check of a non-robust run evaluates the same three
  thresholds without touching the returned design.
* **Optimality criteria with a moving target.** A target outside the
  move-limited box is an error when it is fixed and the nearest reachable
  volume when it moves (the robust formulation's), reported as such.
* **Length-scale check.** The opening removes a one-cell bar and keeps a
  six-cell one; the closing fills a one-cell gap and leaves the open field;
  the scan measures the bar and the gap at one cell. Opening also rounds
  convex corners - three cells per corner of a square block at a two-cell
  probe - which is why a probe passes while it removes at most a tolerance
  of the volume (the test's square blocks measure one cell at a 1 %
  tolerance and three cells at 5 %).

## 22. Pressure, volume and thermal loads

**Unit tests against exact answers** (`tests/test_loads_thermal.cpp`):

* the self-weight of every element type sums to the model's mass times `g`
  to round-off, on straight and curved cells;
* a bar hanging under its own weight is nodally exact on Q4 (the problem is
  one-dimensional at `nu = 0`) and exact everywhere on Tet10, whose space holds
  the quadratic field;
* a rotating bar carries the exact centrifugal stretch; a plane model refuses
  an in-plane rotation axis, and gravity on a massless model is refused;
* the pressure on the curved outer face of a Tet10 quarter cylinder has the
  resultant `-p` times the projected area in `x` and `y` and zero in `z`, to
  `1e-12`; the follower-pressure stiffness `d f / d x` matches central
  differences of the load to `1e-7` on distorted faces of all four shapes;
* free thermal expansion is stress-free on every element and idealisation; a
  fully restrained block carries the hydrostatic thermal stress
  `-E alpha dT / (1 - 2 nu)` (Hex8) and `-E alpha dT / (1 - nu)` (Q4 plane
  stress); a linear temperature gradient bends a free Tet10 bar without
  stress; the thermal prestress of a restrained bar gives the same buckling
  load as the equivalent mechanical compression;
* steady conduction reproduces exact one-dimensional profiles with a heat
  source (Tet10, exact everywhere), convection (Q4) and a surface flux (Hex8),
  and refuses a problem with fluxes alone; a conducted field drives the
  thermal strain of its load case;
* two materials in series give the exact tip displacement.

**Studies against exact continuum solutions** (`apps/verify_loads.cpp`).
Each reference solves the same continuum problem the model discretises, so
the error must vanish at the element's rate: `O(h^2)` in the nodal
displacements and temperatures of the linear elements, `O(h^3)` for the
Tet10, whose edge nodes lie on the curved surfaces. The measured error is the
RMS nodal error relative to the RMS exact field, a discrete L2 norm, and each
study checks the order measured between its two finest meshes against the
element's order less 0.3 (the pre-asymptotic margin). The largest nodal error
is reported too. The curved models are quarter sections with symmetry
supports; the Hex8 and Tet10 sections are one cell deep with `u_z = 0` at
every node, which is exactly plane strain - and on the same mesh the Hex8
errors equal the Q4 errors to every printed digit in both plane-strain
studies, the 3-D pressure and thermal loads reproducing the 2-D ones.

| Study | Element | Finest mesh (`n_r x n_theta`), DOFs | RMS displacement error | Order | Stress error | Order |
|-------|---------|------|-----:|-----:|-----:|-----:|
| Lame cylinder, internal pressure (plane strain) | Q4 | 64 x 128, 16 770 | `4.61e-05` | `2.00` | `7.22e-05` | `1.99` |
| | Tri3 | 64 x 128, 16 770 | `1.02e-04` | `2.00` | `1.30e-02` | `0.98` |
| | Hex8 | 32 x 64, 12 870 | `1.84e-04` | `2.00` | `2.87e-04` | `1.98` |
| | Tet10 | 32 x 64, 75 465 | `4.53e-07` | `3.11` | `8.15e-05` | `1.91` |
| Rotating disk (plane stress) / cylinder (plane strain) | Q4 | 64 x 128, 16 770 | `5.77e-05` | `2.01` | `1.35e-04` | `1.90` |
| | Tri3 | 64 x 128, 16 770 | `1.54e-04` | `2.01` | `1.53e-02` | `0.97` |
| | Hex8 | 32 x 64, 12 870 | `3.17e-04` | `2.01` | `5.90e-04` | `1.80` |
| | Tet10 | 32 x 64, 75 465 | `1.72e-06` | `3.10` | `3.62e-04` | `1.88` |
| Conduction + thermal stress (plane strain) | Q4 | 64 x 128, 16 770 | `6.23e-05` | `2.00` | `1.09e-04` | `1.98` |
| | Tri3 | 64 x 128, 16 770 | `5.46e-05` | `2.00` | `2.20e-03` | `0.97` |
| | Hex8 | 32 x 64, 12 870 | `2.49e-04` | `2.00` | `4.30e-04` | `1.97` |
| | Tet10 | 32 x 64, 75 465 | `3.04e-07` | `3.10` | `1.38e-04` | `1.98` |

* **Lame's thick cylinder**: `a = 0.1 m`, `b = 0.2 m`, `100 MPa` internal
  pressure on the curved bore (`pressures`), steel. The mean radial
  displacement of the bore nodes converges to Lame's `9.53333e-05 m`
  (Tet10: to all six printed digits).
* **Rotating disk and cylinder**: `a = 0.1 m`, `b = 0.3 m`, `600 rad/s`
  about `z` (`centrifugal`), free surfaces; the Q4 and Tri3 models are a
  10 mm disk in plane stress, the solid ones a long cylinder in plane strain,
  each against its exact solution.
* **Conduction and thermal stress**: `k = 45 W/(m K)`, `1 MW/m^3` generation,
  the bore held at 400 K, convection `h = 500 W/(m^2 K)` to 300 K from the
  outer surface, stress free at 300 K. The RMS error of the temperature change
  converges at `2.00` (Q4, Tri3, Hex8) and `3.01` (Tet10, `1.83e-07`), and the
  heat leaving through the bore - from the reactions of the prescribed
  temperatures, against the exact `8 635 W/m` - at `2.00` (`1.53e-05` Q4) and
  to `1.3e-11` on the Tet10.
* **Stresses** are compared at element centroids with the element's
  quadrature-point average, which converges at `O(h^2)` for the Q4, Hex8 and
  Tet10 and at `O(h)` for the constant-strain Tri3, as expected.

The reference itself is checked in each run: one exact axisymmetric solution
(`u = [c_T I(r) - c_b rho omega^2 (r^4 - a^4)/8] / r + C1 r + C2 / r`, with
`I(r)` the weighted integral of the temperature change and `C1`, `C2` from the
surface pressures) matches Lame's and the rotating-disk closed forms of
Timoshenko and Goodier to `8.9e-16` and `6.4e-16`; its equilibrium residual,
by central differences, is below `1.3e-7` of the largest stress; and the
temperature satisfies its convection condition and heat balance to `3e-16`.

* **Bimetallic strip**: 20 mm of two layers, 0.4 mm (`E = 140 GPa`,
  `alpha = 1.5e-6 /K`) under 0.6 mm (`E = 100 GPa`, `alpha = 1.9e-5 /K`),
  heated by 50 K in plane stress with supports that react nothing
  (`max 5.9e-8 N`). Away from the free ends, uniform curvature is an exact
  plane-stress state, and its curvature is Timoshenko's (1925),
  `1.29468 1/m` (the force and moment balance and Timoshenko's formula agree
  to `1.7e-16`). The curvature fitted to the middle half of the bottom edge
  converges at order `1.98 -> 2.00` to an error of `3.0e-04` on 800 x 40
  square cells (tolerance `1e-3`).
* **Hanging bar**: a 0.5 m steel bar of 0.1 m square section (a plate in
  plane stress for Q4 and Tri3) under its own weight, carried by the traction
  `rho g L` on its top face, with three-two-one supports where the exact field
  vanishes - they react `1e-12` of the weight. The exact field is quadratic:
  the Tet10 reproduces it to `4.3e-13`, and the RMS errors of the linear
  elements converge at `2.19` (Q4), `2.00` (Tri3), `2.15` (Hex8) and `1.90`
  (Tet4, 70 227 DOFs, solved by multigrid CG). The largest nodal error
  converges more slowly (orders `1.36 -> 1.64` on the Q4): it sits at the four
  corners where the ends meet the free sides, as the maximum-norm estimate for
  bilinear elements, `O(h^2 |log h|)`, allows - the largest error between
  `0.2 L` and `0.8 L` converges at `1.97, 2.02, 2.00` (Q4) and
  `2.02, 2.01, 2.00` (Tri3). This is why every study judges the RMS error.

Writing the self-weight study exposed a defect in the static solver's force
balance: it divided the residual resultant by the *resultant* of the applied
loads, so any self-equilibrated load - a thermal strain, or a self-weight
carried by a traction - failed on round-off (a relative error of `20` on a
resultant of `5e-13 N`). It now divides by the gross size `sum |f_n|` of the
applied nodal forces, the scale of the round-off in the sum
(`docs/formulation.md`, section 4); for loads acting in one direction the two
are equal, and every earlier study reproduced its value exactly.

**Cross-validation of the loads.** Six decks
(`configs/verification/block_loads_*`, `plate_thermal_*`,
`engine_mount_tet10_loads_analysis.json`) export every load in CalculiX's own
form - `P` faces, `GRAV`, `BX/BY/BZ` on element sets, `CENTRIF`, `*EXPANSION`
with `*TEMPERATURE`, one `*MATERIAL` per material - so CalculiX integrates
the loads with its own code, and a conducted temperature as a steady
`*HEAT TRANSFER` job CalculiX solves itself. scikit-fem integrates the same
loads independently from the same deck, with order-6 rules (the thermal load
with SparLab's stiffness rule). Every comparison of the cross-validation
passes - 137 with the large-deflection decks of section 23, 2 of them
informational, as before:

* scikit-fem's independent integration agrees with SparLab to `1.4e-13` -
  `1.8e-12` on the structured blocks and plates (self-weight, pressure,
  rotation, regional body force, two materials under conducted, uniform and
  regional temperatures) and to `4e-12` - `9e-12` on the curved Tet10 engine
  mount, whose bore pressure acts on curved six-node faces;
* CalculiX agrees to `6.7e-07` - `3.6e-06` wherever its formulation is
  SparLab's, and its conduction solution to `5e-06` - `3.5e-05` of the
  temperature range, at its `.frd` rounding (`5e-4 K` on temperatures near
  300 K);
* three of CalculiX's formulation choices differ from SparLab's. Each was
  identified by reproducing it in scikit-fem, which brings CalculiX to within
  `5.8e-7` - `3.4e-6` of scikit-fem, and each is now reproduced in the
  comparison, so CalculiX is judged against its own problem:
  * for a first-order hexahedron (`C3D8`, and the hexahedra it expands a
    `CPE4` into) CalculiX evaluates the thermal strain at the element-average
    temperature, where SparLab integrates the interpolated temperature - the
    consistent load, which scikit-fem reproduces to `3e-13`. The two agree for
    a uniform temperature (`2.7e-06`) and differ by `6.0e-03` and `4.2e-02` on
    the conducted and the regional field of the Hex8 plate, the difference
    vanishing with refinement for both;
  * for the `C3D10` CalculiX integrates the centrifugal load with its
    four-point rule, not exact for the cubic integrand `N_a rho omega^2 r`
    (SparLab's consistent-mass integration is exact): `2.5e-05` on the block;
  * for the `C3D10` it integrates a face pressure with a three-point rule,
    exact on a flat face but not for the degree-4 integrand of a curved one:
    `2.8e-04` on the engine mount's bore, where SparLab's three-point-per-
    direction rule agrees with scikit-fem's order-4 to order-8 rules to
    `4e-12`;
* CalculiX's plane-strain expansion (`CPE4`) is exact - the two-material Q4
  plate agrees to `2.1e-06` at `nu = 0.3 / 0.33` - so plane-strain rows are
  judged; plane stress at `nu != 0` stays informational;
* CalculiX 2.21 ignores `DC2D4` / `DC2D3` heat-transfer cards (it reads no
  integration point for them), so a plane conduction deck uses the plane
  element itself in the `*HEAT TRANSFER` step.

scikit-fem's check of the stiffness now assembles each element with its own
material's Lame constants, so the two-material decks are the same discrete
problem there too.

## 23. Geometrically non-linear statics

**Unit tests against exact answers** (`tests/test_nonlinear.cpp`, 17 cases):

* *the laws*: for both laws in plane stress (SVK), plane strain and 3-D, at
  a displacement gradient of order 0.3 with rotation, the energy derivative
  is the stress (`1e-7`), the stress derivative is the tangent (`1e-6`), and
  the tangent is symmetric (`1e-12`); at a strain of `1e-7` both reduce to
  linear elasticity, the stress to `1e-6` and the tangent to `1e-5` (the
  stress differs at second order in the strain, relatively `1e-7`); the
  neo-Hookean law refuses plane stress, a temperature and an inverted point;
* *the element*: on all five element types and both laws, the tangent is the
  central difference of the internal force to `1e-6` at a state of up to
  10 % strain and a 0.4 rad rotation, and the internal force is the
  derivative of the energy; on Q4, Hex8 and Tet10 a rigid 90-degree rotation
  leaves no internal force (`1e-12` of `E`) and no energy, and the tangent at
  zero displacement is the linear stiffness to `1e-10`;
* *the assembled system*: on a distorted Hex8 block carrying a follower
  pressure, the deformed-position rotation, self-weight and a temperature
  (`alpha dT = 0.05` at `lambda = 1`), the tangent is the central difference
  of the residual to `1e-6` and the load rate is `-dR/dlambda` to `1e-6` - and
  the tangent is measurably non-symmetric, as a follower pressure over free
  edges makes it;
* *a homogeneous large deformation* (30 % stretch, 20 % shear, a rotation)
  prescribed on the boundary of distorted Q4, Hex8 and Tet10 meshes is
  reproduced at every interior node to `1e-9` m, with the law's Cauchy stress
  in every element (`1e-7`) - the non-linear patch test, both laws;
* *a free cube under pressure on every face* takes the exact homogeneous
  stretch: `s = (-p + sqrt(p^2 + K3^2)) / K3` for a follower pressure (a Cauchy
  stress `-p`) at `p = 0.1 K3`, and the root of `K3 s (s^2 - 1)/2 = -p` for a
  dead one (a first Piola-Kirchhoff stress `-p`), to `1e-9` m, with zero
  reactions, a symmetric assembled tangent over the closed surface and no
  negative pivot;
* *the dead pressure destabilises the cube*: the smallest eigenvalue of the
  tangent at the homogeneous state changes sign at `lambda_c = 0.0456314`
  (found by bisection, independently of the path-following); load control
  must stop below it, with the bracket enclosing it and no wider than 1 %,
  every recorded state homogeneous and stable - which it does, saying that
  the tangent loses positive definiteness;
* *load control and arc length* reach the same state of a Tet10 cantilever
  at `k = 3` to `1e-7`, both with the force and moment balance to `1e-9`;
* *the small-load limit*: on a cantilever, the relative difference from the
  linear solution falls in proportion to `P` (the shortening of the bent
  beam, absent from linear theory) and that of the tip deflection to `P^2`
  (the elastica's `1 - c k^2`), both over two decades of load, down to
  `2e-8` of the deflection;
* *thermal*: a free body heated to `alpha dT = 0.05` takes `F = 1.05 I`
  stress-free on Q4, Hex8 and Tet10; a cube between two walls takes the
  split's closed-form stretch and stress (`1e-10`, formulation section 7c);
* *rotation*: translating a spinning body by `d` normal to the axis changes
  its total centrifugal force by exactly `rho omega^2 V d` on Q4, Hex8 and
  Tet10 (`1e-12`), and the tangent's spin-softening term carries the same
  `-rho omega^2 V d` (`1e-10`) - the test written after CalculiX exposed
  `rho` counted twice in that term (below);
* *the deck*: the `nonlinear` block parses, passes through its `load_factors`
  exactly, its reaction monitor balances the tip force at every step
  (`1e-8`), and eleven malformed blocks are refused; the arc-length method
  refuses a prescribed displacement.

**Studies** (`apps/verify_nonlinear.cpp`; `docs/results/README.md` has the
full tables).

*Euler's elastica* (`--study elastica`). A Tet10 cantilever `L = 1 m` of
square section `h = 0.01 m` (`E = 210 GPa`, `nu = 0`), clamped, carries a dead
shear traction on its tip face whose resultant `P` gives
`k = P L^2 / EI = 1 ... 10` at the load factors `0.1 ... 1`, which load
control passes through exactly. It is meshed with 25, 50 and 100 cells along
the span and 2 x 2 in the section (3 825 to 15 075 DOFs). The reference is
the elastica, `theta'' = -k cos(theta)`, solved by shooting on `theta'(0)`
with fourth-order Runge-Kutta; 20 000 and 40 000 steps agree to `1.2e-14`.
The tip deflection, the shortening and the rotation of the tip face (its
chord from bottom to top edge) are compared at every `k`. A continuum is not
a beam: shear adds `0.6 (h/L)^2` of the deflection at `nu = 0`, and the Saint
Venant-Kirchhoff moment `M = EI c (1 - 0.3 (c h)^2)` softens with the
curvature `c` (`c h <= 0.045` here). The study therefore separates two
errors. The discretisation error must vanish at the Tet10's order, which
Richardson's method measures from the three meshes alone; its limit, the
gap between the converged continuum and the elastica, must be small.

| `k` | Tip deflection / `L`, elastica | Error, 100 cells | Order | Continuum gap | Shortening order |
|----:|------:|------:|-----:|------:|-----:|
| 1 | `0.3017208` | `1.95e-05` | `3.01` | `2.21e-05` | `3.06` |
| 2 | `0.4934575` | `4.00e-05` | `3.29` | `4.43e-05` | `3.36` |
| 5 | `0.7137915` | `8.14e-05` | `3.46` | `8.82e-05` | `3.58` |
| 10 | `0.8106090` | `1.36e-04` | `3.49` | `1.46e-04` | `3.65` |

Every order lies between `3.01` and `3.65` (required: at least 2.7), the
finest-mesh error stays below `1.4e-04 L` in deflection and shortening and
`2.1e-05 rad` in rotation (tolerance `1e-3`), and the error approaches the
gap from below. What remains on 100 cells is therefore the model's gap, not
discretisation. At `k = 1` shear alone accounts for about `1.8e-05 L` of the
`2.2e-05 L`. The linear analysis of the same load would put the tip at
`3.3 L`. Each mesh takes 85 Newton iterations over the ten steps, with no
halving; the largest Green-Lagrange strain is `0.019`.

![Euler's elastica: the tip against the elastica, and the error against the continuum gap](figures/verify_elastica.png)

*A thick tube at finite strain* (`--study hyperelastic-cylinder`). A
quarter section of a tube `a = 0.1 m`, `b = 0.2 m` in plane strain, with
symmetry supports, is loaded three ways in ten steps. The Hex8 and Tet10
sections are one cell deep with `u_z = 0`.

* *inflation*: neo-Hookean (`E = 10 MPa`, `nu = 0.3`) under a bore pressure
  of 1.5 MPa that follows the bore, to a bore hoop stretch of `1.48077`
  (largest Green strain 0.59);
* *spin*: the same material (`rho = 1100 kg/m^3`) at 200 rad/s, the
  centrifugal load at the deformed radius, to a bore hoop stretch of
  `1.17256`;
* *heating*: Saint Venant-Kirchhoff (`E = 1 GPa`, `nu = 0.3`,
  `alpha = 5e-4 /K`) under the conducted temperature of a bore at +100 K and
  an outer surface at 0 K, the multiplicative split, to a bore hoop stretch
  of `1.02502`.

The exact solution is the radial equilibrium of the finite deformation, with
body force and temperature gradient. It is a two-point boundary-value
problem for the deformed radius, solved by shooting with fourth-order
Runge-Kutta (4 000 steps; halving them changes the displacement by less
than `1e-14`). It checks itself in each run: it meets its outer boundary
condition and its own equilibrium equation, by central differences, to
`3.7e-7` of the stress scale, and without load it reduces to Lame's
solution (`2 d(p) - d(2p) = 3.2e-11` at a small pressure `p`). RMS
displacement errors on the finest meshes (64 x 128 for Q4 and Tri3, 32 x 64
for Hex8, 16 x 32 for Tet10), with the order between the two finest:

| Load | Q4 | Tri3 | Hex8 | Tet10 |
|------|----|------|------|-------|
| Inflation | `5.31e-05` (`2.000`) | `1.33e-04` (`2.000`) | `2.12e-04` (`1.999`) | `1.57e-06` (`3.098`) |
| Spin | `5.35e-05` (`2.001`) | `1.07e-04` (`2.001`) | `2.14e-04` (`2.001`) | `2.08e-06` (`3.140`) |
| Heating | `1.10e-04` (`2.004`) | `7.92e-05` (`2.006`) | `4.42e-04` (`2.007`) | `6.04e-06` (`3.237`) |

Every element converges at its order, with a largest shortfall of `9.0e-4`
(tolerance 0.3). The Tet10's bore hoop stretch is `1.4807644` against the
exact `1.4807651`. Every run factorises its tangent by `LDL^T`: the
non-symmetric part of a follower-pressure stiffness is a skew block at the
end nodes of the loaded surface (it integrates `(N_a N_b)_s`), and there the
symmetry supports fix one of the two components, which removes it.

![The tube at finite strain: RMS displacement error under refinement](figures/verify_finite_strain_tube.png)

*Snap-through of a shallow arch* (`--study arch-snap-through`). A clamped
circular arch (half span 1 m, rise 0.1 m, section 0.02 x 0.02 m, plane
stress, `E = 70 GPa`, `nu = 0.3`) is modelled as a half with a symmetry
support at the crown: 60 x 4 Q4 cells, Saint Venant-Kirchhoff, a crown force
of 20 kN at `lambda = 1`. The arc-length method follows its path in 22 steps
(89 iterations, no halving): up to a limit point, down through the
snap-through, whose lowest sample is 752 N at a crown deflection of 0.124 m,
and up again to 20 kN at 0.224 m. Eight of the steps are unstable (their
tangent has a negative pivot). Four checks:

1. *Two methods, one path.* A displacement-controlled run drives the crown
   to exactly the arc-length deflections, and its reaction is the crown
   force: the two agree to `2.8e-11` of the limit force at all 22 points
   (tolerance `1e-6`).
2. *Inertia.* Where the force falls as the crown goes down, the
   force-controlled tangent has exactly one more negative pivot than the
   displacement-controlled one, and elsewhere the same number: 18 of 18
   points (Haynsworth's inertia additivity, section 7c of the formulation).
3. *The limit load.* Displacement control through 60 stations around the
   highest arc-length sample puts the limit at `1588.2464 +- 0.0005 N`, at a
   crown deflection of 0.0374 m. A parabola through the three arc-length
   samples around it, 7 mm apart, would say 1589.12 N (0.055 % high).
4. *Load control stops and brackets it.* Load control of the same arch stops
   at 1588.135 N and reports that the path turns at a limit point, every
   step beyond it landing on a distant branch. Its last converged load and
   its lowest rejected one, 1588.135 N and 1588.257 N, enclose the limit
   load. It does not jump silently to the far branch at 20 kN, as it did
   before the checks of section 7c existed.

![Snap-through of a shallow arch](figures/verify_arch_snap_through.png)

**Cross-validation of the large-deflection decks**
(`configs/verification/*_nonlinear.json`, `python/scripts/cross_validate.py`).
Each deck is solved by `sparlab_solve` and its state at the full load is
compared node by node with two independent solutions:

* an independent total Lagrangian solver written on scikit-fem for this
  purpose. It forms the first Piola-Kirchhoff stress `P = F S` and its
  consistent tangent with scikit-fem's tensor helpers, uses SparLab's
  quadrature, and runs Newton with a tangent predictor and an energy line
  search through SparLab's own load factors. It covers every case whose
  loads are dead; a follower pressure, a rotation or a temperature is
  outside it;
* CalculiX's `*STEP, NLGEOM`, for Saint Venant-Kirchhoff. Its `*ELASTIC` is
  Saint Venant-Kirchhoff under `NLGEOM`, and it uses its own follower
  `*DLOAD` pressure and deformed-position `CENTRIF`, with increments of at
  most 1/50 of the load and residual and correction controls of `1e-6`.

A CalculiX run that did not complete its step is refused (CalculiX writes
its last state to the `.frd` anyway, so the script reads the `.sta` file),
and so is a SparLab run that stopped short of `lambda = 1`.

| Deck | Load case | State at the full load | scikit-fem non-linear | CalculiX `NLGEOM` |
|------|-----------|------------------------|----------------------:|------------------:|
| Elastica, 576 Tet10, `L/h = 50`, `k = 10` | dead tip traction | tip at `0.811 L` down, `0.555 L` back | `6.10e-15` | `6.27e-07` (`C3D10`) |
| Soft block, 256 Hex8 | follower pressure and self-weight | tip down 0.156 m, Green strain 0.16 | - | `3.72e-06` (`C3D8`) |
| | rotation, 300 rad/s | 10 % stretch | - | `1.63e-06` |
| Plane-strain strip, 160 Q4 | dead tip load | tip down 0.075 m | `1.36e-15` | `6.80e-07` (`CPE4`) |
| | follower pressure | tip down 0.057 m | - | `8.74e-07` |
| Rubber block, 192 Tet10, neo-Hookean | tip driven 0.08 m, self-weight | Green strain 0.19 | `8.75e-15` | - |

The independent solver agrees to `1e-14`, the level of the Newton tolerance
and round-off. CalculiX agrees to `6.3e-07` - `3.7e-06`, the size of its
`.frd` rounding. The follower pressures' non-symmetric tangents, factorised
by LU, meet CalculiX as closely as the symmetric cases do. The neo-Hookean
deck has no CalculiX counterpart, since its `NEO HOOKE` is a different strain
energy. The elastica deck runs at `L/h = 50`, where the linear solve of the
study's `L/h = 100` beam would leave the *linear* scikit-fem comparison at
`6.6e-8`, too near its `1e-7` tolerance for a check run on other machines
(formulation, section 5). At 50 it is `8.0e-9`.

What CalculiX was found to do, measured on probe decks while building the
comparison:

* on a cube under a dead pressure of `0.1 K3` it reproduces the exact Saint
  Venant-Kirchhoff stretch to `6.6e-8`, but it follows that homogeneous
  branch past the bifurcation at 4.6 % of the load, where SparLab stops -
  CalculiX assesses no stability;
* its centrifugal load acts at the deformed position but lags within an
  increment. On the spinning block its difference from SparLab falls with
  the increment: `1.03e-5` at 10 increments, `3.5e-6` at 20, `1.8e-6` at
  40 and `1.4e-6` at 80, which is why the export uses at most 1/50;
* its default convergence controls leave errors of order `1e-4` of the
  displacement: `5.6e-5` on the elastica deck at `L/h = 100`, against
  `6.4e-7` with the `1e-6` controls;
* a free thermal expansion comes out exactly `1.05` (`C3D8`, `C3D10`), but a
  cube held between two walls reaches `-43.018 MPa`, where the
  multiplicative split gives `-43.070 MPa` and the version without `1/theta`
  gives `-45.224 MPa`. Its finite-strain thermal model is a third one, and a
  thermal `NLGEOM` comparison is informational.

**What writing these checks found**, each fixed before the numbers above
were taken:

* *the linear solve refused a slender beam.* A Tet10 cantilever at
  `L/h = 100` left a scaled residual of `6.45e-7`, above the `1e-8`
  tolerance, although its backward error was `2.4e-16`: the residual is the
  rounding of sums `2.7e9` times larger than the load. Solves are now
  accepted on either measure (formulation, section 5);
* *Newton stalled on the elastica.* Backtracking on the residual norm cut
  good steps to a quarter, iteration after iteration, because a large
  rotation stirs up axial forces. The energy line search took the run from
  158 iterations and a halving to 87 and none;
* *load control jumped over the arch's limit point* and converged on the far
  branch at 20 kN, a snap-through presented as a static solution. The two
  rejection tests and the bracketing of section 7c stop it now;
* *the spin-softening stiffness counted the density twice*, `rho omega^2 M P`
  with an `M` that already carries `rho`. CalculiX exposed it: on the
  spinning block SparLab reported an instability at `lambda = 0.00017` where
  CalculiX found none. It is fixed, and a unit test measures the term;
* *the thermal strain at finite strain was wrong twice.* Adding `alpha dT`
  to the Green strain gave a free body the stretch `sqrt(1 + 2 alpha dT)`.
  Subtracting the Green strain of the free stretch fixed that but gave a
  restrained body a stress `theta` times too high. The multiplicative split
  is right on both, and is what the study and the unit tests now check;
* *a freely heated body could not converge*: with no applied load and no
  reaction the residual scale was the `1e-300` floor. The thermal forces
  now enter the scale and the round-off floor;
* *the elastica deck stalled at a residual of `3.4e-10`* of the load
  against a tolerance of `1e-10`. At that level the residual is the rounding
  of the displacement itself, about `0.3 eps || |K||u| ||`, so the round-off
  floor now includes `64 eps || |K||u| ||`;
* *a planned test case was unstable.* A cube under a dead compressive
  pressure of `0.1 K3` was meant to test the homogeneous stretch, but the
  cube loses stability at 4.6 % of that pressure. This is physics, not a
  defect, and it became the test of load control's stop.

## 24. Plasticity

**Unit tests against exact answers** (`tests/test_plasticity.cpp`, 20 cases):

* *the hardening laws*: `sigma_y'` is the derivative of `sigma_y` and the
  stored hardening energy its integral (`1e-6`); parameters are validated
  (a negative or non-finite one, a saturation stress without a rate,
  hardening without a yield stress);
* *uniaxial stress*, the lateral strains found by Newton at the point: with
  linear isotropic, kinematic or combined hardening the stress at 8 times
  the yield strain is `E (sigma_y0 + (H + H_kin) eps) / (E + H + H_kin)` to
  `1e-12`, in one step or forty, with the plastic strain, its lateral
  components `-eps_p / 2` and the stored energy exact; reversed to `-eps` it
  follows the isotropic `E (H eps - sigma_y0 - 2 H eps_p1) / (E + H)` and the
  kinematic `E (H_kin eps - sigma_y0) / (E + H_kin)` (the Bauschinger
  effect), the surface keeping its size (kinematic) or growing (isotropic);
  with Voce saturation the stress solves `sigma = sigma_y(eps - sigma/E)` at
  three strains to `1e-11`;
* *far outside the surface* (deviatoric trial increments of order one, as
  the Newton iterates of a large-strain run reach, where one ulp of
  `||xi_tr||` exceeds the multiplier's `1e-13`-of-the-radius tolerance and
  the Voce return used to throw on up to 6 of 401 such strains): three Voce
  laws, with and without Prager's modulus, from the virgin state at 401
  strains `(ln 2.2, ln 1.3, e33)`, `e33` in `[-1, 1]`, and 41 isochoric ones
  with shear, and a plane-stress step to `(ln 2.2, ln 1.3)` - every return
  converges onto its yield surface, to `2.5e-13 sigma_y` (tolerance
  `1e-12`), and the plane-stress one to its `sigma_33` tolerance;
* *the consistent tangent* is the central difference of the returned stress
  to `1e-6` - in 3-D, plane strain and plane stress, with every hardening
  mechanism at once, on a non-proportional second step from a state with
  plastic strain, back stress and accumulated strain, and on an elastic step
  back - and is symmetric; after a plastic step a zero increment gives the
  continuum elastoplastic tangent, which a further increment of `1e-9`
  along the flow follows to `1e-4`;
* *plane stress* leaves `sigma_33` below `1e-11` of the stress and returns
  the 3-D state at the `eps_33` it found; below yield its condensed tangent
  is the plane-stress elasticity matrix;
* *temperature*: free expansion is stress-free at any temperature; a bar
  held axially carries `-E alpha dT` and then `-sigma_y`, cools to a
  residual `E alpha dT - sigma_y`, or yields back to `+sigma_y` when
  `E alpha dT > 2 sigma_y`;
* *the element*, small strain: on all five element types (plane strain,
  plane stress, 3-D), with and without mean dilatation, the tangent is the
  central difference of the internal force at a plastic state to `1e-6`,
  and the thermal load rate the derivative with respect to the load factor;
  an elastic element is the linear one (stiffness and energy to `1e-12`),
  and mean dilatation leaves a constant dilatation alone;
* *the element*, finite kinematics: the same derivative checks at a state
  turned by 0.5 rad, with the geometric stiffness and the Green-strain mean
  dilatation, and the thermal load rate with the Green strain of the thermal
  stretch; a rigidly rotated deformed state keeps its plastic state and
  energy and turns its forces; at zero displacement the tangent is the
  linear stiffness;
* *a bar pulled and pushed back* on distorted Q4, Hex8 and Tet10 meshes
  along the load path `1 -> -1` matches the exact cyclic curve at every step
  - the end force and every element's stress to `1e-9`, the accumulated
  plastic strain to `1e-8`;
* *a homogeneous finite elastoplastic deformation* (3 % stretch, shear and a
  rigid turn) prescribed on the boundary of distorted Q4, Hex8 and Tet10
  meshes is reproduced at every node to `1e-10` m, and every element holds
  the point history of the same Green strains in the solver's steps;
* *heating with finite kinematics*: a body on statically determinate
  supports heated by `alpha dT = 6e-3` expands to `u = alpha dT x` at every
  node to `1e-12` m, without stress or yielding (a linear thermal strain in
  the Green-strain return would stretch it `1.8e-5` short), and the
  homogeneous deformation above, heated as it is loaded, holds at every
  element the return of the Green strain less the thermal stretch's
  `alpha dT (1 + alpha dT / 2)`;
* *plastic collapse*: a perfectly plastic bar under `1.2 sigma_y` stops
  within 1 % below `lambda = 1/1.2` with the collapse bracketed; with
  hardening, arc length and load control end in the same state (`1e-9`) and
  the exact plastic strain;
* *small-strain kinematics with elastic materials* is the linear analysis
  (`1e-10`, one iteration per step), reports its neglected rotation, and
  refuses the neo-Hookean law;
* *the deck* parses the `plasticity` block and the `kinematics`,
  `mean_dilatation` and `load_path` keys, a plastic strip loaded and unloaded
  keeps a permanent set, and thirteen malformed blocks are refused; *the
  CalculiX export* writes one fixed-increment step per leg and a hardening
  table whose chords stay within `1e-4 Q` of the Voce curve.

**Unit tests of Hill48 and Chaboche** (`tests/test_anisotropy.cpp`, 15
cases; the measured worst case of each check in brackets):

* *the coefficients and the frame*: F, G, H, N from the r-values
  (`1.9, 1.5, 2.3`) and `r0 = H/G`, `r90 = H/F`, `r45` (exact to `1e-15`);
  `r = 1` gives von Mises and its matrix (exact); P annihilates the pressure
  and is symmetric in a tilted frame (`2.2e-16`), whose axes are orthonormal
  and right-handed; a quarter turn of RD leaves the criterion unchanged when
  `r0 = r90` (`1.5e-17`); the stress-ratio and explicit calibrations of the
  same sheet reproduce it (`3.9e-16`); thirteen malformed parameter sets are
  refused (non-positive r-values, a non-convex quadratic form, RD along ND,
  five backstresses or four and Prager's modulus, a backstress without
  modulus or with negative recovery, a criterion or backstresses without a
  yield stress);
* *uniaxial stress at 0 to 90 deg from RD* (3-D, five stress components
  found by Newton; plane stress, two), perfectly plastic past yield in one
  step: the stress is `sigma_y phi(theta)^(-1/2)` and the plastic width over
  thickness strain Hill's `r_theta` (`1.1e-13` and `3.3e-14`, tolerance
  `1e-12`), the flow is isochoric with no out-of-plane shear, and the
  accumulated strain is work-conjugate; equibiaxial plane stress yields at
  `sigma_y sqrt(r90 (1 + r0) / (r0 + r90))` with `eps_22 / eps_11 = r0 / r90`
  (`2.3e-15`);
* *Hill48 with linear hardening* along 0 to 90 deg follows
  `E (k sigma_y0 + k^2 H eps) / (E + k^2 H)` in one step and in forty
  (`6.2e-14`, tolerance `1e-12`), with `alpha = k eps_p` and the stored energy;
* *frame indifference*: a sheet and a strain history turned together (in
  the sheet plane, and about a general axis) return the turned stress,
  plastic strain and back stress with Voce hardening and recovering
  backstresses on a non-proportional path (`9.2e-16`);
* *the von Mises limit*: the general return with `r = 1`, and with Prager's
  rule as one backstress without recovery, reproduces the radial return -
  stress, tangent, plastic strain, back stress, energy - along the random
  non-proportional path of the J2 tangent test with every hardening
  mechanism and a temperature, in 3-D, plane strain and plane stress
  (`1.2e-15`, tolerance `1e-12`);
* *Chaboche with Voce* (three backstresses, `gamma = 600, 60, 0`): the
  monotonic curve `sigma_y(p) + sum C_i/gamma_i (1 - exp(-gamma_i p))` at
  strains of `0.2 %`, `2 %` and `20 %`, in one step and 25 (`3.7e-14`,
  tolerance `1e-11`), with every backstress on its own curve; backward Euler
  at first order - errors `6.6e-3, 3.4e-3, 1.7e-3, 8.6e-4` at 10, 20, 40, 80
  steps, ratios within 2.3 % of 2;
* *the reversed branch* `X_i = nu C_i/gamma_i + (X_i0 - nu C_i/gamma_i)
  exp(-gamma_i D)`: loading, the early reverse yield at `sum X_i - sigma_y`
  (the Bauschinger effect), and the reversed branch in one step and twenty
  (`8.3e-16`); *strain cycles* of `+-0.8 %`: each half-cycle's backstress
  follows the recursion (`1.5e-15`), and a stabilised loop started at its
  tip peaks at `sigma_y + (C/gamma) tanh(gamma D / 2)` (`3.5e-15`);
* *the consistent tangent* against central differences (`1.3e-10`,
  tolerance `1e-8`) for Hill48 with Voce, Chaboche with Voce, Hill48 with
  recovering and linear backstresses, Prager's modulus and Voce (in a
  tilted frame and in the sheet plane), exponential and backward Euler, and
  Hill48 with linear backstresses - in 3-D, plane strain and plane stress
  (the tilted frame in 3-D only; the plane states refuse it), with a
  temperature: symmetric to round-off without recovery (`7.2e-16`),
  non-symmetric (above `1e-4`) with it; the converged point is on its
  surface to `1e-13 sigma_y`, so the return from it at the same strain is
  elastic; plane stress returns the 3-D state at its `eps_33`; after a
  plastic step a zero increment gives the continuum tangent, which a
  further increment of `1e-9` along the flow follows (tolerance `1e-6`);
* *the returned state on its surface*: 100 random Hill48 sheets with Voce, a
  recovering and a linear backstress, each along six random increments of
  up to `1e-3`, `1e-2` and `1e-1` per strain component, in 3-D and plane
  stress (3,055 plastic returns): the yield function of the returned state
  is at most `9.1e-15 sigma_y` up to `1e-2` (tolerance `1e-13`) and
  `9.3e-14 sigma_y` at `1e-1`, where the trial stress is some hundred
  `sigma_y` (tolerance `3e-13`); not one re-return from it is plastic, and
  each reproduces the stress to `1e-12`;
* *an independent reference on a non-proportional path*: Hill48 (r-values
  1.9, 1.5, 2.3, L = 1.3, M = 1.7) in a frame tilted out of every
  coordinate plane, Voce and linear hardening, two recovering backstresses,
  a linear one and Prager's modulus, along three straight 3-D strain
  segments with a reversal, against the rate equations integrated by the
  classical Runge-Kutta method in 40,000 substeps a segment (its own Hill
  matrix, polarised from the quadratic form in its own frame; it leaves the
  surface by `2.2e-14 sigma_y`): the return converges at first order,
  errors `1.27e-3, 6.32e-4, 3.16e-4, 1.58e-4` (exponential) and
  `1.36e-3, 6.81e-4, 3.41e-4, 1.70e-4` (backward Euler) at 100 to 800 steps
  a segment, ratios `2.000 +- 0.003` (tolerance `2 +- 0.05`);
* *the plane models*: a Hill48 frame whose axes miss z is refused by the
  plane-strain and plane-stress return (`ConfigError`) and by the
  non-linear analysis before its first step, as by the deck;
* *the element*: Hill48 with recovering backstresses on Q4 (plane strain,
  plane stress), Hex8 and Tet10, small strain and finite kinematics turned
  by 0.5 rad, with and without mean dilatation - the unsymmetrised tangent
  is the central difference of the internal force (`5.8e-10`, tolerance
  `1e-8`) and is non-symmetric; without recovery it stays symmetric; a
  rigidly rotated finite state of a sheet in a tilted frame keeps its
  plastic state and energy and turns its forces;
* *a homogeneous deformation* of Hill48 (a tilted frame in 3-D, a rolling
  angle of 35 deg in plane strain) with recovering backstresses and Voce,
  loaded and reversed past reverse yield on distorted Q4, Hex8 and Tet10
  meshes with small strain and finite kinematics: every node to `1e-12` m
  (`9.5e-18` m) and every element's stress - every component, `sigma_33` of
  plane strain too - and plastic strain the point history in the solver's
  steps (`1.8e-13`, tolerance `1e-11`);
* *the deck*: the Hill48, Chaboche and integration keys parse in their three
  calibrations and with a rolling angle, a Hill48-Chaboche strip is loaded
  and unloaded through the non-symmetric tangent (at most 5 Newton
  iterations a step; 7 with the tangent symmetrised), which its plastic
  steps factorise by LU (no inertia reported) where the same strip without
  recovery keeps LDL^T, `summary.json`
  reports the law, the CalculiX export refuses it with the reason, and
  twenty-three malformed blocks are refused.

**Studies** (`apps/verify_plasticity.cpp`; `docs/results/README.md` has the
full tables).

*A thick tube to plastic collapse* (`--study plastic-cylinder`). A quarter
section of a tube `a = 0.1 m`, `b = 0.2 m` in plane strain (`E = 200 GPa`,
`nu = 0.3`, `sigma_y = 250 MPa`, no hardening), small strain, is loaded by its
bore pressure along the arc-length path past its collapse. The exact
collapse pressure of a von Mises tube in plane strain,
`p_L = (2/sqrt 3) sigma_y ln(b/a) = 200.094 MPa`, holds for any Poisson
ratio - at collapse the elastic strain rates vanish and the flow is
isochoric, so `sigma_z = (sigma_r + sigma_theta)/2` - and on the plateau the
whole wall carries the exact fully plastic field,
`sigma_r = (2/sqrt 3) sigma_y ln(r/b)`,
`sigma_theta = sigma_r + 2 sigma_y / sqrt 3`. First yield is at
`0.540 p_L` (Lame's field with `sigma_z = 2 nu A`). The largest load factor
of the path is a lower bound on the discrete collapse load (every converged
state is an equilibrium within yield) which the plateau attains; the study
also records how much the load factor rose over the last ten steps, zero on
a true collapse plateau.

| Element | Collapse error, finest mesh | Order | Plateau rise | Plateau stress error / `sigma_y` (order) |
|---------|---------------------------:|------:|-------------:|------------------------------------------|
| Q4, mean dilatation (the default), `n_r = 32` | `1.20e-04` | 2.00 | `9e-10` | `1.49e-04` (2.01) |
| Hex8, mean dilatation (the default), `n_r = 16` | `4.82e-04` | 2.00 | `9e-10` | `6.00e-04` (2.00) |
| Tet10 (the default, no mean dilatation), `n_r = 8` | `2.99e-06` | 3.02 | `2e-07` | `2.33e-03` (0.96) |
| Q4 fully integrated, `n_r = 16` | `1.20e-02` | 1.99 | `2.0e-03` | `1.21e-02` (1.98) |
| Tet10 with mean dilatation, `n_r = 8` | `-7.12e-04` (below) | 2.05 | `4e-08` | `2.37e-02` (0.85) |
| Tri3 (checkerboard split), `n_r = 16` | `3.93e-04` | 2.00 | `1.5e-09` | `9.03e-03` (1.01) |

The one-cell-deep Hex8 section held at `u_z = 0` is exactly plane strain and
reproduces the Q4 to the last digit. With mean dilatation the Q4 converges
to the exact collapse pressure at second order and its plateau is flat;
fully integrated it locks - 18 % high on the coarsest mesh, still 1.2 % at
16 cells, its plateau rising - as Nagtegaal, Parks and Rice (1974) predict.
The measurement also settled the default for Tet10: its four-point element
converges at third order without mean dilatation, while with it the
collapse load converges at second order from below and the constant element
pressure oscillates (the stress error converging at order 0.85). The
checkerboard split of the structured Tri3 mesh - the crossed pattern whose
triangles, unlike a single diagonal's, leave divergence-free fields enough
freedom - does not lock here. The pass criteria judge the defaults: the
collapse error below `1e-3` at order `>= 1.8` and the plateau stress within
`0.01 sigma_y`.

*Pure bending* (`--study plastic-bending`). A beam `0.2 m x 0.05 m`, `0.01 m`
thick, in plane stress (`sigma_y = 250 MPa`, no hardening) is bent by end
displacements `u_x = -k (x - L/2)(y - h/2)` on square Q4 cells (4 to 32
through the depth). The section is in uniaxial stress, which plane stress
makes compatible with the free transverse strain, so the exact moment is
`E I k` up to `k_y = 2 sigma_y / (E h)` and `M_p (1 - (k_y/k)^2 / 3)` beyond.
The moment from the end reactions at `k = 0.5, 1, 1.5, 2, 3 k_y` converges at
order 1.94 to `6.17e-04 M_p`. Loaded to `3 k_y` and unloaded along the load
path to `k_res = 1.5556 k_y`, where the exact moment vanishes, the beam keeps
`1.63e-04 M_p` and the residual stress - the loaded profile less the elastic
unloading, which stays short of reverse yield - to `2.33e-03 sigma_y`
(order 1.65).

*A uniaxial cycle* (`--study plastic-cycle`). A bar on a distorted
`4 x 2 x 2` Hex8 mesh (mean dilatation) is strained through
`0 -> 1 % -> -1 % -> 1 %` in 60 steps, with linear isotropic hardening
(1 GPa), a Voce saturation (100 MPa at rate 30) and Prager kinematic
hardening (4 GPa). The end force over the area matches the exact uniaxial
response - the algebraic equations of the uniaxial state, solved at every
step - to `1.98e-14 sigma_y`: the reversal yields early (the Bauschinger
effect) and the loop grows as the isotropic hardening accumulates.

*A uniaxial cycle with Chaboche backstresses* (`--study chaboche-cycle`).
The same cycle, `0 -> 1 % -> -1 % -> 1 %` in 20 steps a leg, on a
distorted `4 x 2 x 2` Hex8 bar (mean dilatation) and a distorted `4 x 2`
plane-stress Q4 bar, of a steel with `sigma_y = 170 MPa`, Voce 90 MPa at
rate 12 and three Armstrong-Frederick backstresses
`(C, gamma) = (60 GPa, 600), (8 GPa, 60), (1 GPa, 0)`, integrated
exponentially. The exact response solves, step by step, the branch solution
`X_i = nu C_i/gamma_i + (X_i0 - nu C_i/gamma_i) exp(-gamma_i D)` with the
consistency condition; the end force over the area matches it to
`3.9e-15 sigma_y` (tolerance `1e-9`) over the 120 steps. The reverse yield
stress at the first turn is `-18.6 MPa` against a forward flow stress of
about 400 MPa - the Bauschinger effect of the saturating backstresses - and
Newton needs at most 4 iterations a step on the non-symmetric tangent.

*Hill48 directional properties* (`--study hill-directional`). A unit
square (distorted `4 x 4` plane-stress Q4) and a `1 x 1 x 0.5` block
(distorted `3 x 3 x 2` Hex8) of a sheet with `r0 = 1.9, r45 = 1.5,
r90 = 2.3`, `sigma_y = 170 MPa` and linear hardening 1 GPa, rolled at
`-theta` from x for `theta = 0, 15, ..., 90 deg`, are pulled along x to a
strain of 1 % in 10 steps, with supports that leave the lateral strains and
the in-plane shear of off-axis flow free so that uniaxial stress is exact.
The end force over the area matches the directional hardening curve
`E (k sigma_y + k^2 H eps) / (E + k^2 H)`, `k = phi(theta)^(-1/2)`, to
`1.2e-15`; the plastic width over thickness strain, from the corner
displacements (and incompressibility in plane stress), matches Hill's
`r_theta` to `5.6e-15`; every node lies on the homogeneous field to
`1.0e-15` of the applied displacement (tolerance `1e-9` on all three).

**Cross-validation** (`python/scripts/cross_validate.py`, the ten
`plastic_*` decks). scikit-fem solves every deck with an independent J2
implementation - the radial return in 3 x 3 tensor form (Newton on the
multiplier for Voce, on the thickness strain in plane stress), a material
tangent by central differences of that return so that no tangent formula is
shared, the mean dilatation where SparLab applied it, Newton converged at
each of SparLab's load factors with the internal variables committed there;
with finite kinematics the return in the Green-Lagrange strain
`E = (H + H^T + H^T H) / 2`, the internal force `int dE : S dV_0` with the
displacement-dependent strain operator `sym(F^T grad du)`, the geometric
stiffness and, where SparLab applied it, E-bar with its geometric term. On
the six small-strain decks it agrees to `3.9e-15` - `6.8e-12`; with finite
kinematics to `1.6e-14` on a slender Tet10 cantilever deflected a tenth of
its span, and to `3.4e-13`, `3.4e-13` and `2.3e-13` on a steel beam clamped
at both ends whose mid-span section is driven 20 mm down - ten times its
first-yield deflection, into membrane action (a 150 kN end tension on
Hex8) - and back to zero, which pushes the plastically stretched beam into
compression and reversed yielding, with isotropic, Voce and kinematic
hardening: on Hex8 and plane-strain Q4 with E-bar, and on plane-stress Q4
(its thickness strain from the Green-strain return). These comparisons
discriminate: switched off in the reference, E-bar moves the result by
`1.2` (Hex8) and `1.0` (Q4) of the largest residual displacement, and
small-strain kinematics by `0.09` to `0.12`. CalculiX's `*PLASTIC` (an
isotropic hardening table, without `NLGEOM`, the same fixed increments, one
`*STEP` per leg of the load path) agrees to `1.8e-06` - `3.3e-06`, within
its `.frd` rounding, on a Hex8 cantilever loaded past yield and unloaded, a
plane-strain strip, a Tet10 punch partly relieved and a Tet10 plate heated
past yield and cooled (its ramped `*TEMPERATURE` and the additive thermal
strain). Against CalculiX's finite-strain plasticity under `NLGEOM` the
Tet10 cantilever differs by `1.8e-04` - two different models, which agree
to the order of the plastic strain (up to `2.7e-3` here); informational.

What CalculiX was found to do: its `HARDENING=KINEMATIC`, given the table of
a linear rule (250 MPa at zero plastic strain, 2 250 MPa at 0.1), softens a
single `C3D8` in uniaxial tension at `-E H_kin / (E - H_kin)` - the rate at
which Prager's rule hardens, with the opposite sign - and `COMBINED`
saturates; its users report similar anomalies
(https://calculix.discourse.group/t/isotropic-hardening-kinematic-hardening/790).
Its isotropic hardening reproduces the exact reversed curve to its printed
digits (409.0909, -516.5289 and -698.3471 MPa). Kinematic hardening is
therefore not exported, and the two decks with it are compared with
scikit-fem only.

**What writing these checks found**, each fixed before the numbers above
were taken:

* *load control "converged" past the collapse load.* Driven to 1.05 `p_L`,
  the perfectly plastic tube reported completion at strains of `1e11`: as
  the displacement ran away, the round-off floor `64 eps || |K||u| ||` grew
  with it past the load itself, and every residual passed. A floor now
  counts only below `1e-6` of the residual's scale (formulation, section
  7c), and the run stops with the collapse bracketed. The same run exposed
  that successes short of a collapse reset the halving count, so load
  control crept towards it until the step budget ran out; the nearest load
  factor no step reached is now tracked like a critical point
  (`unreached_load_factor`);
* *the default mean dilatation was wrong for Tet10*: it was on for every
  multi-point element until the tube showed the Tet10 better without it;
* *a linear static summary reported a force-balance error of `6.8e21`* for
  a case driven by prescribed displacements alone, dividing by an applied
  force of zero; such a case is now measured against its reactions (formulation,
  section 4);
* *a homogeneous finite-deformation test with a rigid turn of 0.5 rad
  failed*, correctly: prescribing `lambda (F - I) X` runs along the chord of
  the rotation and compresses the body by 3 % at mid-path, a pressure of
  5 GPa whose geometric stiffness outweighs the plastic tangent, and the
  homogeneous state loses stability. The test uses 0.1 rad;
* *the independent J2 solver diverged* on a prescribed tip displacement,
  whose jump strains the next row of elements far past yield, until it was
  given the tangent predictor SparLab uses; a comparison against a
  `summary.json` without the materials' plasticity parameters (written
  before they were recorded) compared a plastic run with an elastic one,
  which the script now refuses;
* *the independent finite-kinematics solver stalled* on the slender Tet10
  cantilever at `3e-11` of the load - below the rounding of `u` itself,
  `eps || |K||u| || = 7.2e-7 N`, and above its `1e-11` tolerance, the same
  floor SparLab's run sat on (`3.4e-11`, within its `1e-10`). It now accepts,
  as SparLab does, a residual at its round-off floor, once Newton has
  stopped reducing it.

### Logarithmic-strain finite plasticity

`kinematics: finite_logarithmic` (formulation, section 7d). **Unit tests**
(`tests/test_logarithmic.cpp`, 14 cases):

* *the strain and its derivatives* at C = Q diag(lambda) Q^T in a general
  frame: distinct eigenvalues, two equal, `C = I`, `C = s^2 I`, two and three
  eigenvalues within `1e-9`, `1e-7` and `1e-4`, three spread over `4e-3`
  (the Taylor series) and `6e-3` (the difference quotient), and a stretch of
  2 with two lateral eigenvalues `1e-8` apart. `E_log` matches an independent
  matrix logarithm (Eigen's Schur-Pade) to `4.8e-15` (tolerance `1e-14`) and
  `ln J` its trace; `P` and `T : L`, for a stress not coaxial with C, match
  fourth-order central differences of `E_log(E)` and of `S(E) = P^T T` to
  `4.1e-12` and `4.3e-12` (tolerance `1e-11`: the differences' own
  accuracy, their round-off `eps / h`); `T : L` is symmetric and `S = P^T T` to `1e-15`; two eigenvalues
  `1e-9` apart give the `P` and `T : L` of the equal pair to `2e-9`; `C = I`
  gives `P = I`; at strains of `1e-6` the strain is `E - E^2 + (4/3) E^3` to
  `1e-15` (the `log1p` of the eigenvalues of E); an E that no deformation
  has (`C` not positive definite) or a non-finite one throws `SolverError`;
* *the element tangent* of Hex8, Q4 in plane strain and in plane stress and
  Tet10, with and without mean dilatation, elastic (Hencky), J2 with Voce,
  linear and Prager hardening, Hill48 with Voce, and Chaboche (non-symmetric),
  from a plastic state at a stretch of 1.35 with shear and a turn of 0.6 rad
  to a second one: the central difference of the internal force to `1.1e-9`
  (tolerance `1e-8`), symmetric exactly when the law is; the thermal load
  rate the derivative of the force with respect to the temperature scale
  (`2.1e-11` in 3-D, `9.3e-10` in plane stress; tolerance `1e-8`);
* *objectivity*: Hill48 in a tilted frame with a recovering backstress, Hex8
  and Tet10, with and without mean dilatation - a superposed rotation of
  0.9 rad leaves the energy (`1e-12`), the plastic strain, back stress and
  accumulated strain (`1e-11`), `S` and the log strain unchanged, and turns
  the forces, the Kirchhoff and the Cauchy stress (`1e-11`);
* *inverted points*: Hex8, Q4 in plane strain and plane stress and Tet10,
  with and without mean dilatation - a reflection `diag(-1/2, 1, 1)` (whose
  C, log strain and energy are those of a valid compression) and a general
  inverted F are refused by the element and by the stress recovery with
  `SolverError`, while a compression to `det F = 0.05` evaluates and a turn
  by pi (`det F = +1`, negative entries) gives the energy of the unturned
  state (`1e-12`); a Hex8 block, elastic and with J2, whose end is pushed to
  `u_x = -1.5 L` in 10 steps stops short of zero volume (`lambda < 2/3`)
  with every converged state valid - without the check the elastic run
  converged through it to `lambda = 1` at `det F = -0.76`;
* *the small-strain limit*: at strains of `1e-6` the forces and tangents of
  the logarithmic and the Green-Lagrange Hex8, Tet10 and Q4 differ by
  `3.3` to `11.9` times the strain (tolerance 20 times), elastic or yielding,
  and ten times less at `1e-7` (ratio `10.000`, within 5 %) - first order;
* *uniaxial tension to a stretch of 2* (log strain 0.69) of one Hex8 and of
  a distorted `2 x 2 x 2` patch on symmetry planes, J2 with linear and with
  Voce hardening, 20 steps: the Kirchhoff stress (end force times stretch
  over the reference area) against the 1-D law at `ln(l1)` at every step to
  `4.3e-14` (tolerance `1e-12`), the two lateral log strains equal and adding
  with the axial one to the elastic volume change `(1 - 2 nu) tau / E`
  (`7.2e-16`), the accumulated plastic strain exact (`6.4e-16`), and every
  element's Kirchhoff stress, log strain and Cauchy stress `tau / J` those of
  the exact state (`4.8e-14`; tolerance `1e-12` throughout); no small-strain
  warning;
* *a tension-compression cycle to `+-0.3` log strain* with Voce hardening
  and two Armstrong-Frederick backstresses (exponential), through
  `0 -> e^0.3 -> e^-0.3 -> e^0.3`: the Kirchhoff stress follows the exact
  branch solutions to `7.0e-14 sigma_y` and the accumulated plastic strain to
  `3.6e-15` (tolerance `1e-12`), with reverse yielding;
* *Hill48 along RD and along TD* to a stretch of 1.6: the Kirchhoff stress
  follows `E (k sigma_y + k^2 H eps) / (E + k^2 H)` to `2.4e-14` and the
  plastic lateral log strains keep the ratio `r0 = 1.9` and `r90 = 2.3` to
  `6.3e-15` (tolerance `1e-12`) at every plastic step;
* *plane stress*: a distorted Q4 strip pulled to a stretch of 1.6 follows the
  1-D law (`9.8e-15`), and its thickness log strain, found by the return,
  equals the in-plane lateral one (`2.1e-15`), with `J = l1 l2
  exp(E_log,33)` (`1.2e-14`; tolerance `1e-12` throughout);
* *free heating*: Q4 (plane stress), Hex8 and Tet10 blocks heated by 500 K on
  determinate supports take `u = alpha dT x` (to `1e-13` m), stress-free
  (`2.8e-6` Pa against `E alpha dT = 1.2 GPa`; tolerance `1e-12 E alpha dT`),
  with the log strain `ln(1 + alpha dT)` in every direction, the thickness's
  included (`2.2e-15`; tolerance `1e-12`);
* *the verification decks*: every deck of `configs/verification` with a
  non-linear static or transient analysis, switched to
  `finite_logarithmic` - small-strain decks included - completes (17 runs);
  the neo-Hookean one is refused. The plastic clamped beam, driven 20 mm
  down (strains to 2.1 %, 2 % plastic) and back, agrees with its
  Green-Lagrange run to the order of that strain: mid-span force within
  `1.42e-2` of its largest value, end tension within `1.57e-2 sigma_y A`,
  accumulated plastic strain within `0.5 %`;
* *the deck*: `"kinematics": "finite_logarithmic"` parses, the run reports
  `hencky` and the largest log strain, the stress CSV and VTK carry the
  Kirchhoff stress and the log strain, CalculiX export says why it cannot,
  and `neo_hookean` with it is refused.

**Studies.** *Large-strain uniaxial tension and cycle*
(`--study logarithmic-uniaxial`): the tests' five uniaxial runs on the
distorted `2 x 2 x 2` Hex8 patch (J2 linear and Voce to a stretch of 2,
Chaboche to `+-0.3` log strain, Hill48 along RD and TD to 1.6; 122 steps,
accumulated plastic strain up to 1.48): the Kirchhoff stress against the 1-D
law in the log strain to `1.15e-13 sigma_y` and the r-values to `6.3e-15`
(tolerance `1e-12`).

*The thick tube at finite strain* (`--study logarithmic-tube`): the tube of
`plastic-cylinder` with logarithmic kinematics and a follower bore pressure,
along the arc-length path past its peak (strains to 4 - 5 %). At finite
strain the peak falls below the small-strain limit as the wall thins:
`0.99160 p_L` and `0.99156 p_L` for Tet10 with 4 and 8 cells through the
wall, the reference. Q4 and Hex8 with the mean dilatation of `ln J` reach it
at second order - `7.6e-3`, `1.9e-3`, `4.7e-4` above it on 4, 8, 16 cells
(order `2.01`; tolerance `1e-3` on the finest) - while the fully integrated
Q4 locks: `0.117`, `8.8e-3`, `1.9e-3` above it.

What remains: no cross-validation of the logarithmic kinematics against
another code (CalculiX's finite-strain plasticity is multiplicative, and its
elasticity under `NLGEOM` Saint Venant-Kirchhoff), no exact solution with
rotating principal axes at large strain, and the elastica deck, which
completes, needs 74 steps instead of 10: the first Newton iterate from the
linear predictor lands where the Hencky tangent has negative pivots, and
load control halves the first step three times.

## 25. Dynamics

**Unit tests against exact answers** (`tests/test_dynamics.cpp`, 12 cases;
three more in `tests/test_io.cpp` for the deck blocks, the output files and
the CalculiX `*DYNAMIC` decks):

* *amplitudes*: the step, table and harmonic values and rates (the table's
  slope to the right of a point), the HHT-alpha parameters and their range;
* *the scalar trapezoidal recursion* - the reference of the next tests - has
  the exact discrete period, a phase advance `tan(theta / 2) = omega dt / 2`
  per step at constant amplitude, to `1e-12` over 200 steps;
* *a transient run is the HHT-alpha recursion of every mode*: on a Q4 plate
  and a Hex8 block, consistent and lumped mass, `alpha = 0` and `-0.1`,
  Rayleigh damping, a harmonic load from rest and a preload released, the
  global solution equals the sum of every mode integrated by the scalar
  recursion in acceleration form (the solver uses the displacement form) to
  `1e-10` of the largest displacement at every step;
* *energy*: the trapezoidal rule balances `E_0 + W - T - U - D` to `1e-12`,
  undamped and damped; HHT-alpha dissipates;
* *prescribed motion*: a free plate shaken at one edge balances its energy to
  `1e-12`, and the sum of its reactions equals the rate of change of its
  momentum `sum (M a + C v)` to `1e-10`;
* *the harmonic response is the sum over all modes*, with structural and
  Rayleigh damping, at and between resonances, consistent and lumped, Q4 and
  Hex8, to `1e-10`; velocity and acceleration monitors are `i omega U` and
  `-omega^2 U` to `1e-14`; the largest displacement of a node over a cycle
  is the semi-major axis of its orbit (in phase, circular and elliptical
  orbits, and against a 3 600-point sampling of the cycle);
* *the non-linear transient* of a linear model (small strain, elastic) is the
  linear transient to `1e-10` in at most two Newton iterations a step; with
  finite kinematics it converges at second order in the step (observed
  orders above 1.8 for the motion and the energy balance); an elastoplastic
  plate loaded slowly and damped comes to rest at the non-linear static
  state, the gap closing as the ramp lengthens, with positive plastic
  dissipation;
* *an undamped response at a natural frequency* is flagged (or the LU fails
  outright), and with damping it is finite and unflagged;
* the input checks: a duration that is not a whole number of steps, `alpha`
  outside `[-1/3, 0]`, negative damping, no frequencies, a massless model;
* a prescribed rigid motion of a body passes the static equilibrium check
  (`tests/test_static_verification.cpp`): its reactions are round-off, now
  judged against the gross force `|K| |u|` they are formed from - a shaken
  base, whose static analysis precedes its dynamic one, used to fail that
  check with a relative error of 0.13 of round-off over round-off.

**The transient against the exact discrete solution** (`transient-modal`).
Q4 (20 x 4, plane stress) and Hex8 (10 x 2 x 2) steel cantilevers,
consistent and lumped mass, five cases each - a sudden tip load with the
trapezoidal rule and with `alpha = -0.3`, a harmonic load at `1.3 f1` with
2 % Rayleigh damping at `alpha = 0` and `-0.1`, and a static preload released
with `alpha = -0.05` - 120 steps of `T1 / 40`, against every mode of a dense
generalised eigensolve, each integrated by the scalar recursion. The
reference is computed in 80-bit extended precision: in double precision the
modal sum is itself only good to about `omega_max^2 / omega_1^2 eps`, and the
differences then measured (up to `7.9e-10`) were mostly the reference's. Against the
extended-precision reference the largest difference over the 20 runs is
`3.28e-10` of the largest displacement (`1.83e-10` over the first half of the
runs: the solver's double-precision round-off accumulates step by step - run
once with 60, 120 and 240 steps, the undamped Q4 case measured `6.5e-12`,
`1.75e-11` and `3.4e-11`). The trapezoidal energy balance is at most
`2.98e-11`; the same recursion carried out by the study in double and in
80-bit arithmetic gives `1.08e-11` and `6.33e-15`, a ratio of `1.71e3` for a
ratio of `2.05e3` in eps - round-off. `alpha < 0` dissipates, as it should:
0.36 % to 0.40 % of the work of the harmonic load at `alpha = -0.1`, 36 % to 49 % of the
sudden load's at `alpha = -0.3` - a sudden point load excites every mode, and
the method damps those the step does not resolve.

**A fixed-free rod** (`rod-harmonic`, `rod-transient`). A steel rod 1 m long,
0.05 m square, `nu = 0`, its lateral displacements held, one element across
- an exactly one-dimensional model on Q4 and on Hex8, which give the same
numbers. Harmonic: driven by an end force or by its base, undamped between
resonances, with structural damping `eta = 0.02` and with 1 % Rayleigh
damping at the first and third resonances; on 10 to 160 elements with
consistent and lumped mass. Every node equals the exact solution of the
discrete equations - `u_j = g cos(j theta) + D sin(j theta)` from the
dispersion relation of the stencil - to `2.94e-10`, and the reaction as
well. The end amplitude converges to the exact damped continuum solution,
`u = g cos(k (L - x)) / cos(kL) + F sin(kx) / (E* A k cos(kL))` with the
complex modulus and density of the damping, at order `2.00` on the finest
pair in every case; the largest error on 160 elements is `3.6e-3`, at the
third resonance, where a frequency error of `O(h^2)` is amplified by the
sharpness of the peak. Transient: an end force ramped as `sin^2` over
`0.6 T1`, 2.5 periods, the trapezoidal rule at Courant number 1 (`dt = h /
c`): the end displacement converges at order `2.01` (consistent mass) and
`2.00` (lumped) to the continuum's modal series, whose static part is summed
exactly and whose 4 000-mode truncation changes it by `6.4e-13` of `F L / (E
A)`; on 160 elements the error is `2.6e-5` (consistent) and `8.0e-5`
(lumped) of the largest displacement, and the energy balance at most
`3.9e-12`. `docs/figures/verify_rod_dynamics.png` shows the harmonic sweep
through three resonances against the exact response, both convergence
studies and the end displacement against the series.

**Non-linear oscillators** (`nonlinear-oscillator`). One element (Q4 in plane
strain, Hex8) with its lateral displacements held: uniaxial strain, a single
degree of freedom `m u'' + N(u) = F`, under a sudden load, 20 to 320 steps
per elastic period, `alpha = 0` and `-0.1`, lumped and consistent mass. The
Saint Venant-Kirchhoff law at finite strain (`N = A (lambda + 2 mu) F11
(F11^2 - 1) / 2`, a peak strain of 9 %) and J2 small-strain plasticity with
linear hardening (0.8 of the yield force, doubled by the sudden
application: the element yields on its first swing and then oscillates
elastically about its permanent set). Against the scalar HHT-alpha recursion
of the same equation - Newton to machine precision, the plastic law's radial
return in closed form - every run agrees to `1.74e-11`. Against the exact
motion - classical Runge-Kutta at 1/64 of the finest step for the elastic
law (the Runge-Kutta reference agrees with itself at half its step to
`1.3e-14`), the piecewise closed form for the plastic one (harmonic with the
elastic slope, then with `K + (4/3) mu H / (3 mu + H)`, then elastic again
about the shifted equilibrium) - the error falls at second order: `1.98` to
`2.00` for the elastic law, `1.94` to `2.04` for the plastic one, whose
order wanders as the yield point falls at a different place within a step.
The run's final energy balance - the energy not found as kinetic, stored or
external work - holds the plastic dissipation: `174.794 J` at 320 steps per
period against the exact `sigma_y alpha_p V = 174.802 J` (`4.5e-5` at
`alpha = 0`, `4.6e-5` at `-0.1`). The gap shrinks with the step, but
unevenly - the yield point moves within a step - and at `-0.1` the method's
numerical dissipation enters the balance too (its coarsest run, 20 steps a
period, overshoots to `174.91 J`).
`docs/figures/verify_nonlinear_oscillator.png` shows both motions and the
convergence.

**Cross-validation** (section 14's two codes, seven decks, `configs/
verification/transient_*.json` and `frequency_response_*.json`). scikit-fem
integrates every transient again with its own `K` and `M` - consistent, or
lumped as SparLab lumps (row sums; HRZ for the Tet10) - and HHT-alpha in the
acceleration form, and solves every harmonic deck directly; CalculiX runs the
exported `*DYNAMIC, DIRECT, ALPHA` decks of the solid-element transients:

| Deck | Reference | Max relative difference | Tolerance |
|------|-----------|------------------------:|----------:|
| Hex8 cantilever, 3-D tip load ramped, HHT `alpha = -0.05`, Rayleigh damping, 300 steps | scikit-fem HHT-alpha / CalculiX `C3D8 *DYNAMIC` | `3.27e-10` / `1.64e-06` | `1e-7` / `1e-5` |
| Tet10 cantilever driven harmonically from rest, trapezoidal rule, 400 steps | scikit-fem / CalculiX `C3D10 *DYNAMIC` | `3.46e-09` / `2.08e-06` | `1e-7` / `1e-5` |
| Plane-strain Q4 strip shaken at its root, lumped mass, Rayleigh damping | scikit-fem | `1.50e-09` | `1e-7` |
| Q4 cantilever plate, harmonic tip force through four resonances, structural and Rayleigh damping | scikit-fem direct complex solve | `1.02e-10` | `1e-7` |
| Hex8 block on a shaken base, lumped mass, structural damping, 80 frequencies to 8 kHz | scikit-fem direct complex solve | `8.15e-13` | `1e-7` |
| Hex8 cantilever loaded dynamically past yield, small-strain J2 (plastic strain `1.0e-3`) | scikit-fem's own J2 in Newton / CalculiX `C3D8 *DYNAMIC *PLASTIC` | `2.17e-08` / `6.21e-07` | `1e-7` / `1e-5` |
| Slender Hex8 strip swinging through large deflection (tip down 0.36 L), Saint Venant-Kirchhoff | scikit-fem's own finite-kinematics system / CalculiX `C3D8 *DYNAMIC NLGEOM` | `4.90e-10` / `1.46e-06` | `1e-7` / `1e-5` |

The scikit-fem comparisons are judged on the full-precision data - the
monitors at every step (or frequency) and the final displacement, velocity
and acceleration (the complex fields at the snapshot frequencies) - with the
snapshot fields read from the VTK series (nine significant digits, agreeing
to `5e-9` or better) checked beside them. The largest of them is the final
acceleration of the non-linear runs: `a = (u1 - u~) / (beta dt^2)` magnifies
the Newton-tolerance difference of `u` by about `1 / (beta (omega dt)^2)`,
120 for the plastic cantilever, whose displacements agree to `1.2e-12`.
The CalculiX comparisons are of the displacement fields at every snapshot,
read from its `.frd` file (six significant digits, a rounding of `5e-6`).

What the cross-validation found, and how it was resolved:

* **CalculiX's expanded plane elements in `*DYNAMIC`.** A CPS4 cantilever
  whose first period is 12.1 ms (CalculiX's own `*FREQUENCY` of the same
  deck agrees with SparLab's 82.3 Hz to `8e-4`) peaks under a sudden tip
  load after 4.1 ms rather than half a period; CPE4 alike, the two 86 % and
  99 % from SparLab, while a cantilever of C3D8 elements agrees to the
  output rounding.
  CalculiX's dynamic plane elements are therefore not used; the plane decks
  are compared with scikit-fem, the solid ones with both codes.
* **A load acting at `t = 0`.** SparLab starts in equilibrium with it (the
  initial acceleration `M^-1 A(0) f`); CalculiX does not (measured: 2.6 %
  apart in the first step of a sudden load). The export refuses such a run,
  with the reason, and the decks ramp their loads from zero.
* **`ALPHA`.** CalculiX's HHT-alpha has SparLab's sign convention and Newmark
  parameters: `alpha = 0` and `-0.1`, with and without Rayleigh `*DAMPING`,
  agree to the `.dat` rounding (`1.6e-7`) while the two `alpha` differ by
  `2.1e-4`; without `ALPHA` it takes `-0.05` (reproduced digit for digit).
* **Mean dilatation.** CalculiX's C3D8 averages no dilatation, so the plastic
  deck integrates in full, as CalculiX does, to be judged against it: the
  fully integrated Hex8 is stiffer (the tip load had to rise from 13.5 kN to
  20 kN to yield the root), which is the element, not the dynamics.

## 26. Incremental forming

The forming analysis (`docs/forming.md`) is checked in `tests/test_forming.cpp`
(`[forming]`), each claim against an exact answer:

* **Trajectories**: interpolation, clamping outside the span, the right
  derivative, monotone and random lookups through the cached segment, and
  the CSV reader's errors, which name the line.
* **The contact tangent** against central differences (step 1e-9 m) of the
  contact residual at fixed friction history, for a sphere, a tilted plane
  and a cylinder, frictionless, sticking, slipping and newly in contact
  (13 to 215 nodes in contact), and a 2-D circle in slip: the largest error
  relative to the largest tangent entry is `7.2e-10` (sphere `7.0e-10`, plane
  `5.4e-10`, cylinder `7.2e-10`, 2-D `1.7e-10`); every node of the slip cases
  slips and every node of the stick cases sticks.
* **A flat punch on an elastic column** (one Hex8 in section, so every top
  node has the same penalty stiffness; uniaxial stress): the punch force
  against `E A delta / H / (1 + E / (kappa H))`, the penalty in series,
  `4.5e-12`, `1.9e-11` and `3.6e-10` at `s` = 10, 100, 1000; against rigid
  contact `1.96e-2`, `2.0e-3`, `2.0e-4` - exactly `1/s`. A 3 x 3 section
  (edge nodes stiffer) lies within the compliance bound.
* **The symmetric friction tangent** (`friction_tangent: "symmetric"`): at
  the tangent test's states it gives the same forces and symmetric node
  blocks, positive semi-definite for the plane in stick and slip, and the
  flat-punch ironing below completes with it without a cut (48 iterations,
  the exact tangent 31) at the same tool force to `1e-8`.
* **Ironing** with `mu = 0.2`: a flat punch pressed into an elastic block and
  dragged 10 mm at 30 degrees to x, and a sphere (R = 20 mm) dragged along x.
  In steady sliding every node in contact slips (231 and 9 nodes). The
  punch's friction resultant lies along the drag to `7e-16` with
  `|F_T| = mu F_z` to `2e-16`, and every node's friction force lies along the
  drag to `3e-15` (all tested to `1e-6`) - x and y in the drag's proportion,
  which a friction force that did not follow the slip would miss. The
  sphere's `-F_x / F_z = 0.19986`, `mu` to `6.9e-4` (tested to `3e-3`), its
  contact normals tilted.
* **Double-sided forming**: two spheres on the two faces of a clamped elastic
  sheet, their paths mirror images about its mid-plane, pressed in and moved
  together with friction, take equal and opposite normal forces and equal
  friction forces to `1e-6`, each face touching only its own tool; the
  resultant of the clamp's reactions equals the total force on the tools to
  `1e-6`.
* **Rigid-body invariance of a release**: a stress-free block released onto
  3-2-1 supports stays exactly where it is, and follows a translation of the
  supports to `1e-12`, in one iteration per increment; a plastically bent
  strip (finite kinematics) released twice from the same restarted state,
  the second time with the supports' values moved by a translation and by a
  0.03 rad rotation, ends moved rigidly: largest deviation `6.3e-18` and
  `1.1e-17 m` (tested to `1e-14 m`), its plastic history unchanged, in the
  held release's 41 iterations without a cut (the rigid predictor).
* **Tools between steps**: a tool active in two steps whose windows leave a
  gap over which its path moves is refused before the run (and on a restart
  from the first step's state), naming the steps and the 2.4 mm jump; held
  over the gap, or taken away by a release in between, it is accepted; a
  tool that becomes active in contact is warned about.
* **A start state that cannot be evaluated** (a sphere whose path starts on
  the surface): the step is recorded, not completed, with the state it
  started from, and the run's files are written.
* **Springback against the exact unloading**: an elastic-perfectly plastic
  plane-stress beam (L = 8 h) bent by end displacements to `3 k_y`, where
  `M = M_p (1 - (k_y/k)^2/3)` (the moment-curvature relation of section 24),
  then released onto 3-2-1 supports; the curvature of the central half
  springs back by `M(3 k_y)/(E I) = 1.4444 k_y`, measured `1.4565`, `1.4482`,
  `1.4453 k_y` with 4, 8 and 16 Q4 through the depth: errors `8.4e-3`,
  `2.6e-3`, `5.6e-4` (order 1.7, then 2.2); the supports' reactions end below
  `1e-6` of the reference force.
* **Single-point incremental forming** (the smoke case of `docs/forming.md`
  on a coarse 12 x 12 x 2 mesh with 12 chords a contour, about 17 s on one
  thread): completes; the blank's centre sinks 1.70 mm, the release moves it
  by up to 0.16 mm; the tool force is positive in z at every increment in
  contact. (The 20 x 20 x 2 deck, `configs/forming/spif_smoke.json`, sinks
  1.65 mm and springs back 0.18 mm.)
* **The deck and the output contract**: every key of the `forming` block read
  strictly, a trajectory file beside the deck, each refusal, and a run
  writing every file of the contract with its exact header and row counts.

Every one of them also passes in a build without SuiteSparse
(`-DSPARLAB_WITH_CHOLMOD=OFF`).

## 27. The Hex8 with incompatible modes and the rule through a sheet

The element of `docs/formulation.md` ("The Hex8 with incompatible modes",
"A sheet's rule through its thickness") is checked in
`tests/test_incompatible_modes.cpp` (24 cases, 2 337 assertions; the whole
suite now has 350 cases and 29 090 assertions) and by four studies
(`apps/verify_elements.cpp`) and a variant of a fifth.

**Unit tests**, each against an exact answer, a central difference or an
invariant:

* *quadrature*: the 1- to 7-point line rules integrate polynomials of
  degree `2n - 1` exactly; the box rule with equal orders is the cube rule
  bit for bit, and puts its thickness points along the axis asked for; a
  Hex8 with 5 thickness points has 20 stiffness and stress points, and on
  a box (where 2 x 2 x 2 is exact) the same stiffness; a 3-D deck's
  thickness rule gives 20 plastic points per element;
* *the condensed stiffness* is symmetric with exactly six zero eigenvalues
  on a box and a distorted cell; the modes take no part in a linear field on
  a distorted cell (`alpha = 0`, exact strains); one element in pure bending
  (`nu = 0.3`, `L/h = 1, 10, 100`) reproduces the exact 3-D field at its
  nodes to round-off where the standard Hex8 locks; the stiffness is frame
  invariant (`K(RX) = R K R^T`) and independent of the node numbering; 2-,
  3- and 4-point rules agree on a box;
* *patch tests*: the linear constant-strain patch on meshes perturbed by 0,
  0.15 and 0.30, a free thermal expansion (uniform and linear) stress-free,
  and a homogeneous finite deformation (elastic, J2, Hill48, all three
  kinematics) with `alpha = 0` and the standard element's forces;
* *consistent tangents* by central differences (`h = 1e-8`, relative
  `1e-6`) after a committed step: small strain, finite (Green-Lagrange) and
  logarithmic, each elastic, J2, Hill48 (through a 3- or 5-point thickness
  rule) and Chaboche (non-symmetric, asserted so); the elastic forces are the
  derivative of the stored energy; Saint Venant-Kirchhoff and neo-Hookean
  through the total Lagrangian path; the thermal load rate is the derivative
  of the force in the load factor;
* *objectivity*: a rigid turn of 0.7 rad gives no force and no modes, and a
  turn superposed on a plastic state turns the force and the parameters;
  the tangent at rest is the linear `K*` under every kinematics;
* *the local iteration*: one iteration for elastic small strain; at most
  three for every converged step of a strip bent past yield and unloaded,
  whose springback equals the linear elastic solution (`1e-8`); a reversal
  to six times the strain of a committed plastic step in one increment
  converges in 3 (small strain, J2 and Chaboche), 9 (logarithmic, J2) and
  12 (logarithmic, Chaboche) iterations; a reversal of 2 % logarithmic
  strain to 12 to 30 % (perfectly plastic J2, hardening J2, Chaboche,
  Hill48; 5 and 7 thickness points), where `K_aa` is indefinite, in 9 to
  18;
* *the non-linear small-strain solution* of an elastic cantilever equals
  the linear one (`1e-9`), within 2 % of Timoshenko on 8 x 1 x 1 cells;
* *the history*: `evaluate` is bitwise identical on one and three threads;
  a forming restart from the state after a bend reproduces the one-run
  result to `1e-12`, and without the committed parameters to `1e-8`; a
  state of the other formulation is refused;
* *stability*: a neo-Hookean cube in uniaxial compression - the tangent of
  the compression test loses definiteness between the stretches 0.66 and
  0.65 (hourglassing), the standard Hex8's stays positive to 0.3; equibiaxial
  thinning to a thickness stretch below 0.5 stays stable;
* *configuration*: the deck keys, the refusals (unknown formulation, 8 or -1
  thickness points, an unknown direction, a Tet or plane mesh, mean
  dilatation `all`).

**Slender cantilever** (`--study im-cantilever`). `L = 1 m`, `h = 0.1 m` (one
element), `b = 0.05 m` (one element), `E = 70 GPa`, clamped root face, a
`-1000 N` tip resultant; mean tip deflection over Timoshenko's:

| Mesh | Incompatible modes, `nu = 0` | `nu = 0.3` | Standard, `nu = 0` | `nu = 0.3` |
|------|---:|---:|---:|---:|
| 10 x 1 x 1 | `0.9965` | `0.9870` | `0.666` | `0.658` |
| 20 x 1 x 1 | `0.9984` | `0.9930` | `0.888` | `0.818` |
| 40 x 1 x 1 | `0.9989` | `0.9953` | `0.969` | `0.872` |

Against Euler-Bernoulli the incompatible 20 x 1 x 1 mesh gives `1.0044`. At
`L/h = 100` and `1000` (20 x 1 x 1, `nu = 0`) the incompatible modes give
`0.99937` and `0.99957` of Timoshenko, the standard Hex8 `0.0741` and
`0.00080`: its error grows as `(L/h)^2`, the modes' does not. (At
`L/h = 1000` the relative force balance of the direct solve is `2e-4`, from
the conditioning; the solve accepts `1e-3`.)

**MacNeal-Harder straight beam** (`--study im-macneal-harder`; MacNeal and
Harder 1985). 6 x 0.2 x 0.1, `E = 1e7`, `nu = 0.3`, 6 x 1 x 1 cells, unit tip
loads; tip deflection over the references `0.1081` (in plane) and
`0.4321` (out of plane):

| Cells | Incompatible, in plane | out of plane | Standard, in plane | out of plane |
|-------|---:|---:|---:|---:|
| rectangular | `0.978` | `0.973` | `0.093` | `0.025` |
| parallelogram (45 deg) | `0.978` | `0.836` | `0.021` | `0.013` |
| trapezoidal (+-45 deg) | `0.047` | `0.030` | `0.026` | `0.010` |

Judged: rectangular cells and the parallelogram in plane (affine cells, on
which Taylor's modes are exact in pure bending). Reported: the parallelogram
out of plane, whose inclined clamped root couples bending with twist (the
same body on 48 x 8 x 8 cells gives `0.945` out of plane and `0.950` in
plane: it is not the reference beam), and the trapezoidal cells, which are
not affine and on which the modes lose most of their effect (MacNeal 1987).

**Thin plate** (`--study im-plate`). A square plate `a = 1 m`, `t = 0.02 m`
(`t/a = 1/50`), steel, `1e4 Pa` on the top face, a quarter on `n x n x 2`
cells; the centre deflection over Kirchhoff's `0.00126 q a^4/D` (clamped side
faces) and `0.00406 q a^4/D` (`u_z = 0` on the edges' mid-surface nodes):

| `n` | Incompatible, clamped | simply supported | Standard, clamped | simply supported |
|----:|---:|---:|---:|---:|
| 4 | `0.831` | `0.919` | `0.091` | `0.128` |
| 8 | `0.977` | `1.000` | `0.276` | `0.366` |
| 16 | `0.996` | `1.009` | `0.587` | `0.682` |

The 3-D plate differs from Kirchhoff's by its transverse shear and its 3-D
supports, of order `(t/a)^2`; the simply supported plate converges from
above for that reason.

**Springback of a strip** (`--study im-springback`). A steel strip 20 x 2 x
2 mm (`E = 200 GPa`, `nu = 0.3`, `sigma_y = 250 MPa`, linear hardening
`1 GPa`), small strain, bent by end rotations to `4 k_y` and released onto
statically determinate supports by the forming driver; springback = change
of the end-face relative rotation over `L`. Reference: the standard Hex8 on
4 and 8 layers of cubic cells (3 075 and 19 683 unknowns),
Richardson-extrapolated to `1.8468 1/m` (the 8-layer value is `3.6e-3` from
it):

| Model | Unknowns | Springback [1/m] | Error |
|-------|---:|---:|---:|
| incompatible, 1 layer, 2 x 2 x 5 | 132 | `1.7877` | `3.2e-2` |
| incompatible, 1 layer, 2 x 2 x 7 | 132 | `1.8356` | `6.0e-3` |
| incompatible, 1 layer, 2 x 2 x 2 | 132 | `2.1792` | `0.18` |
| incompatible, 2 layers, 2 x 2 x 3 | 198 | `1.8261` | `1.1e-2` |
| standard, 1 layer, 2 x 2 x 5 | 132 | `2.1788` | `0.18` |

One incompatible-mode layer represents the linear strain through the
thickness and unloads exactly (`M/(E I)`); what is left is the Gauss
quadrature of the kinked stress: the moment of a rigid-plastic shell outside
the elastic core (`|z| < t/8`) is short by `3.5e-2` with 5 points and
`8.8e-3` with 7. Two points a layer are far from enough. The 8-layer
reference takes 125 s of the study's 131 s (2 threads, load average 2 to 4).

**Plastic collapse** (`--study plastic-cylinder`, reported variant). The
incompatible-mode Hex8 without mean dilatation does not lock under the
isochoric flow: collapse pressure within `3.52e-4` of the exact limit load
on the finest mesh (order `2.00`), plateau rise `9e-10`, fully plastic
stress `5.4e-4` of `sigma_y` (the averaged Hex8: `4.82e-4`, `6.0e-4`).

## What is not covered

Stated plainly, since the absence matters as much as the presence:

* **no comparison against experiment**;
* the forming analysis is verified on exact solutions only; no forming run
  is compared with another code or with a measured part, and the springback
  of the SPIF cases is not mesh-converged (their meshes are coarse for a
  5 mm tool, and Hex8 through a sheet's thickness is stiff in bending);
* the cross-validation covers linear static displacements on thirty-eight
  problems, five of them on the three meshes read from files - six of them
  under pressure, body and thermal loads, with conducted temperatures on four
  - linear buckling load factors on three, the final large-deflection
  states of four, eight comparisons in all, the final elastoplastic states
  of ten, fifteen comparisons, the transient histories of five, nine
  comparisons (two of the decks non-linear), and the harmonic responses of
  two. Stresses, natural frequencies, the non-linear static load paths (only
  the final states are compared) and the optimised designs are not compared
  with another code, and CalculiX's
  `*BUCKLE` factors for
  `C3D8` and `C3D10` differ from SparLab's by up to `8.3e-5` for a reason
  not identified (section 14);
* linear buckling is bifurcation of the perfect geometry, and the plane
  models only buckle in their plane. The non-linear analysis has no branch
  switching, so a post-buckling path is computed only from an imperfect mesh,
  and none is verified;
* plasticity is small-strain J2 or Hill48 with rate-independent linear,
  Voce, Prager and Chaboche (Armstrong-Frederick) hardening: no finite-strain
  plasticity (with `finite` kinematics
  the return works in the Green-Lagrange strain, sound for large rotation
  with small strain; it is cross-validated against an independent
  implementation of that model, and against CalculiX's finite-strain model
  only as an informational comparison), no rate dependence, creep, damage,
  fracture or non-associative flow, no non-quadratic anisotropic yield. The
  exact solutions cover the
  collapse load and fully plastic field of a tube, pure bending with
  unloading and uniaxial cycles (Prager and Chaboche), and for Hill48 the
  directional yield stress and r-value of homogeneous states; locking was
  measured on the tube only, so
  the Tri3 result there (no locking with the checkerboard split) does not
  carry to other meshes. Kinematic hardening is not cross-validated against
  CalculiX (its implementation does not reproduce Prager's rule), only
  against scikit-fem and the exact uniaxial solutions; Hill48 and the
  Chaboche backstresses are not cross-validated against another code at
  all; on a non-proportional path the general return is measured against
  an independent Runge-Kutta integration of the rate equations (first
  order), not against an exact solution;
* the non-linear analysis is verified against exact solutions of a beam
  theory (the elastica, small strain), of plane-strain finite elasticity
  (the tube: inflation, spin, heating) and by the consistency of three
  solution methods on one snap-through. There is no exact 3-D finite-strain
  solution with shear-dominated deformation among them, and no experiment.
  The neo-Hookean law is cross-validated against scikit-fem only
  (CalculiX's `NEO HOOKE` is a different strain energy). A follower pressure
  over free edges, whose tangent is non-symmetric, is checked against
  CalculiX `NLGEOM` on two decks and in the unit tests. The finite-strain
  thermal split is checked against the exact tube and the unit tests only:
  CalculiX's thermal model at finite strain is a third one (its stress on a
  restrained cube lies 0.12 % from the split's), so the cross-validation
  records a thermal `NLGEOM` comparison as informational, and none of the
  committed non-linear decks carries a temperature. The bracketing of a
  critical point is verified on
  one limit point (the arch) and one bifurcation (the cube);
* dynamics is verified against the exact solution of the discrete equations
  (modal superposition, the rod's dispersion relation, the scalar recursion
  of a one-element oscillator), the rod's continuum solutions and the exact
  motion of the one-element oscillators. A non-linear transient with many
  degrees of freedom - large deflection, plasticity - is cross-validated
  against scikit-fem and CalculiX, not verified against an exact answer, and
  its accuracy in time is only as good as the step: the energy balance
  bounds nothing about the error of the motion. The harmonic response is
  cross-validated against scikit-fem only (CalculiX's steady-state dynamics
  is modal); CalculiX's plane elements are not used in dynamics (their
  `*DYNAMIC` response contradicts CalculiX's own `*FREQUENCY`); there is no
  experiment, no damping identified from a test, and no validation of
  Rayleigh or structural damping as a model of a real structure's damping;
* the overhang filter and check are verified for the 3- and 5-element
  stencils of structured square and cubic grids; the robust formulation for
  uniform erosion and dilation only. Neither is a process simulation;
* on a curved Tet10 cell the 4-point stiffness rule is not exact (the
  integrand is rational); the same rule in scikit-fem and CalculiX
  reproduces SparLab's answer, which verifies the implementation, not the
  rule;
* the thermal and pressure studies run in plane strain and the rotating disk
  in plane stress; the conduction studies cover prescribed temperatures,
  convection and generation on curved boundaries, and surface fluxes only in
  the unit tests and the cross-validation decks;
* the sensitivity checks run on 72- and 36-element meshes (they need two
  extra solves per element per step); the gradients are not FD-verified at
  benchmark resolution, although they are the same code path;
* the stress-constraint gradient is FD-verified; the *constraint's* claim
  - that the relaxed aggregate bounds the stress of a part - is checked only
  by the re-solve of the thresholded structure, which is a consistency check,
  not a proof;
* the Hex8 mesh-convergence study still stops at 26 481 DOFs and runs on
  the direct solver, so its numbers stay comparable with earlier runs; the
  Tet4 study reaches 59 211 DOFs, where `auto` has switched to multigrid;
* the multigrid solver is verified against the direct solver up to 47 775
  DOFs in the study and at the benchmark sizes in `docs/benchmarks.md`;
  beyond the direct solver's reach (the 356 475- and 830 115-DOF runs) there
  is no second solution to compare with, only the residual checks;
* the mesh readers are tested on files written by the tests and by Gmsh
  4.15.2; files from other generators (Abaqus/CAE, HyperMesh, Salome) use the
  same keywords but have not been tried;
* the incompatible-mode Hex8 is verified on patch tests, exact pure bending,
  tangents and invariants, and validated in bending against beam and plate
  theory and one converged springback case; the export to CalculiX
  (`C3D8I`) is not cross-validated (whether CalculiX applies Taylor's
  correction on distorted cells is not known here), and no forming run with
  it is compared with another code or with a measured part;
* no convergence study of the *optimised topology* against mesh size in the
  verification suite - that lives in the design study, where the
  `mesh_fixed_r` and `mesh_fixed_cells` arms address it directly.
