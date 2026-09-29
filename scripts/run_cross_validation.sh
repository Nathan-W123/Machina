#!/usr/bin/env bash
# Cross-validate the static solver against CalculiX and scikit-fem.
#
# usage: scripts/run_cross_validation.sh [results-dir]
#
# Solves the reference analysis decks - the 2-D cantilever and the 3-D block
# on Q4 / Hex8 and on Tri3 / Tet4 elements, the two Gmsh parts (the lug
# bracket at its real Poisson ratio and at nu = 0, the engine mount on Tet4
# and on curved Tet10 cells), and an axially loaded column on Hex8, Tet4 and
# Tet10 with its linear buckling check - and the decks of the volume,
# pressure and thermal loads (a block under self-weight, pressure, rotation
# and a body force on Hex8 and Tet10; a two-material plate with conducted,
# uniform and regional temperatures on Hex8, Tet10 and plane-strain Q4; the
# curved Tet10 engine mount under a bore pressure, self-weight, rotation and
# conduction) - and the geometrically non-linear decks (the elastica on
# Tet10, a soft Hex8 block under a follower pressure and a rotation, a
# plane-strain Q4 strip under a dead load and a follower pressure, a
# neo-Hookean Tet10 block with a driven tip) - with sparlab_solve (exporting
# CalculiX decks), then compares the nodal displacements node by node, the
# conducted temperatures, the buckling load factors mode by mode and the
# final non-linear states, with CalculiX (ccx: static, *BUCKLE, *HEAT
# TRANSFER and *STEP, NLGEOM, each load in CalculiX's own form) and
# scikit-fem (with SparLab's load vector, with the loads integrated by
# scikit-fem, and an independent total Lagrangian solve). Exits non-zero if
# any comparison exceeds its documented tolerance.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

RESULTS="${1:-$SPARLAB_ROOT/results}"
PYTHON="${PYTHON:-python3}"
require_binaries sparlab_solve

# The decks of the volume, pressure and thermal loads.
LOAD_CASES="block_loads_hex_analysis block_loads_tet10_analysis plate_thermal_hex_analysis \
plate_thermal_tet10_analysis plate_thermal_q4_analysis engine_mount_tet10_loads_analysis"
# The geometrically non-linear decks.
NONLINEAR_CASES="elastica_tet10_nonlinear block_hex_nonlinear strip_q4_nonlinear \
block_tet10_neohookean_nonlinear"

banner "cross-validation: solving the reference decks"
for case in cantilever_analysis block_3d_analysis cantilever_tri_analysis \
            block_tet_analysis lug_bracket_nu0_analysis engine_mount_tet10_analysis \
            column_hex_buckling_analysis column_tet4_buckling_analysis \
            column_tet10_buckling_analysis $LOAD_CASES $NONLINEAR_CASES; do
  "$BIN_DIR/sparlab_solve" --config "$SPARLAB_ROOT/configs/verification/$case.json" \
                           --output "$RESULTS/$case" --export-calculix
done
# The Gmsh parts as static analyses (their topology sections are ignored).
for case in lug_bracket_2d engine_mount_3d; do
  "$BIN_DIR/sparlab_solve" --config "$SPARLAB_ROOT/configs/benchmarks/$case.json" \
                           --output "$RESULTS/${case}_analysis" --export-calculix
done

banner "cross-validation: CalculiX and scikit-fem"
cd "$SPARLAB_ROOT"
CASES=()
for case in cantilever_analysis block_3d_analysis cantilever_tri_analysis \
            block_tet_analysis lug_bracket_nu0_analysis lug_bracket_2d_analysis \
            engine_mount_3d_analysis engine_mount_tet10_analysis \
            column_hex_buckling_analysis column_tet4_buckling_analysis \
            column_tet10_buckling_analysis $LOAD_CASES $NONLINEAR_CASES; do
  CASES+=(--case "$RESULTS/$case")
done
"$PYTHON" python/scripts/cross_validate.py "${CASES[@]}" --output "$RESULTS/cross_validation"
