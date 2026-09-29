#!/usr/bin/env python3
"""Write the contour toolpath of a single-point incremental forming run.

The tool's centre (the reference point of a sphere tool of radius ``R``)
over a sheet whose top face is at ``z = top``: approach above the first
contour, plunge to the first level, then one full circle per level, each
level ``dz`` deeper and ``dr`` tighter than the one before (a cone with a
wall angle of ``atan(dz / dr)``), a straight step-down between levels, and
a retract. Pseudo-time: the plunge ends at t = 1, each circle takes 1, each
step-down and the retract 0.1. The circles are polygons of ``--segments``
chords (the analysis reaches every knot exactly).

The CSV has the header ``t,x,y,z`` (s, m) that ``forming.tools[].trajectory
.file`` reads, preceded by ``#`` comment lines with the parameters; the
times at which the forming and the retract end are printed, for the steps'
``time`` windows. Run from the repository root, e.g. the smoke case::

    python3 python/scripts/make_forming_toolpath.py --tool-radius 0.005 \\
        --radius 0.008 --dr 0.001 --dz 0.001 --levels 2 --segments 24 \\
        --output configs/forming/spif_smoke_toolpath.csv
"""

from __future__ import annotations

import argparse
import math


def contour_path(tool_radius: float, radius: float, dr: float, dz: float, levels: int,
                 segments: int, top: float, clearance: float, retract: float):
    """Knots (t, x, y, z) of the path, and the times the forming and the
    retract end."""
    knots = [(0.0, radius, 0.0, top + tool_radius + clearance)]
    t = 1.0
    for level in range(levels):
        r = radius - level * dr
        z = top + tool_radius - (level + 1) * dz
        if level == 0:
            knots.append((t, r, 0.0, z))  # the plunge
        else:
            t += 0.1
            knots.append((t, r, 0.0, z))  # the step-down
        for i in range(1, segments + 1):
            a = 2.0 * math.pi * i / segments
            knots.append((t + i / segments, r * math.cos(a), r * math.sin(a), z))
        t += 1.0
    formed = t
    last = knots[-1]
    knots.append((t + 0.1, last[1], last[2], top + tool_radius + retract))
    return knots, formed, t + 0.1


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--tool-radius", type=float, required=True, help="sphere radius [m]")
    parser.add_argument("--radius", type=float, required=True,
                        help="radius of the first contour of the tool's centre [m]")
    parser.add_argument("--dr", type=float, required=True, help="radius step per level [m]")
    parser.add_argument("--dz", type=float, required=True, help="step-down per level [m]")
    parser.add_argument("--levels", type=int, required=True, help="number of contour levels")
    parser.add_argument("--segments", type=int, default=32, help="chords per circle")
    parser.add_argument("--top", type=float, default=0.0, help="z of the sheet's top face [m]")
    parser.add_argument("--clearance", type=float, default=5.0e-4,
                        help="height of the start above the sheet [m]")
    parser.add_argument("--retract", type=float, default=1.0e-3,
                        help="height of the retract above the sheet [m]")
    parser.add_argument("--output", required=True, help="CSV file to write")
    args = parser.parse_args()
    if args.levels < 1 or args.segments < 3:
        parser.error("--levels must be >= 1 and --segments >= 3")
    if args.radius - (args.levels - 1) * args.dr <= 0.0:
        parser.error("the last contour's radius is not positive")
    knots, formed, retracted = contour_path(args.tool_radius, args.radius, args.dr, args.dz,
                                            args.levels, args.segments, args.top,
                                            args.clearance, args.retract)
    with open(args.output, "w", encoding="utf-8") as out:
        out.write("# SPIF contour toolpath of the tool's centre (python/scripts/"
                  "make_forming_toolpath.py)\n")
        out.write(f"# tool radius {args.tool_radius} m, first contour {args.radius} m, "
                  f"dr {args.dr} m, dz {args.dz} m, {args.levels} level(s), "
                  f"{args.segments} chords per circle, sheet top at z = {args.top} m\n")
        out.write(f"# forming ends at t = {formed:.6g} s, the retract at t = {retracted:.6g} s\n")
        out.write("t,x,y,z\n")
        for t, x, y, z in knots:
            out.write(f"{t:.10g},{x:.10g},{y:.10g},{z:.10g}\n")
    length = sum(math.dist(a[1:], b[1:]) for a, b in zip(knots, knots[1:]))
    print(f"wrote {args.output}: {len(knots)} knots, path length {length * 1e3:.1f} mm; "
          f"forming ends at t = {formed:.6g} s, the retract at t = {retracted:.6g} s")


if __name__ == "__main__":
    main()
