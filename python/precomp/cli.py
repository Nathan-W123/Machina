"""The `precomp` command.

    precomp part        make a target part (a family and its parameters, or a sample)
    precomp material    list the material library or print one material
    precomp setup       write a FormingSetup JSON with defaults and overrides
    precomp toolpath    tool path of a surface (SparLab trajectory / robot CSV)
    precomp simulate    simulate forming a commanded surface with sparlab_form
    precomp compensate  displacement-adjustment compensation of a target
    precomp scan compare   align a scan to its target and measure the deviation
    precomp scan update    one DA step of the commanded surface from a scan
    precomp report      Markdown report and figures of target / formed / deviation
    precomp dataset generate, precomp train, precomp evaluate
                        machine learning (delegated to precomp.ml, imported only
                        when one of these runs)

Every number on the command line is SI (metres, pascals, seconds); a point
cloud in millimetres is read with --scale 0.001. Configuration files are
JSON. Exit status: 0 on success, 2 on a usage or input error, 3 when a
simulation fails.
"""

from __future__ import annotations

import argparse
import importlib
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from . import __version__
from ._util import PrecompError, canonical_json, read_json, to_jsonable, write_json

ML_COMMANDS = ("dataset", "train", "evaluate")


def _parse_value(text: str) -> Any:
    """A command-line value: JSON if it parses (numbers, booleans, lists,
    objects), else the plain string."""
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return text


def _parse_assignments(items: Optional[Sequence[str]]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for item in items or []:
        if "=" not in item:
            raise ValueError(f"expected key=value, got {item!r}")
        key, value = item.split("=", 1)
        out[key.strip()] = _parse_value(value.strip())
    return out


def _load_setup(path: str):
    from .fea.setup import FormingSetup
    return FormingSetup.from_dict(read_json(path))


def _load_map(path: str):
    from .geometry.heightmap import HeightMap
    return HeightMap.load(path)


# -- part ---------------------------------------------------------------------
def cmd_part(args: argparse.Namespace) -> int:
    import numpy as np

    from .geometry.heightmap import Grid
    from .geometry.parts import families, part_from_dict

    fams = families()
    if args.list:
        for name, cls in sorted(fams.items()):
            print(f"{name}: {cls.__doc__.strip().splitlines()[0]}")
            for key, (lo, hi) in cls.bounds.items():
                print(f"    {key}: [{lo:g}, {hi:g}]")
        return 0
    if args.json_in:
        part = part_from_dict(read_json(args.json_in))
    else:
        if not args.family:
            raise ValueError("give --family (or --json-in, or --list)")
        if args.family not in fams:
            raise ValueError(f"unknown family {args.family!r}; known: {sorted(fams)}")
        cls = fams[args.family]
        if args.sample:
            part = cls.sample(np.random.default_rng(args.seed))
        else:
            part = cls(**_parse_assignments(args.param))
    if not args.out and not args.json_out and not args.stl:
        print(canonical_json(part.to_dict()), end="")
        return 0
    grid = Grid.centered(args.size, args.spacing)
    hm = part.heightmap(grid)
    if args.out:
        hm.save(args.out)
    if args.stl:
        hm.to_stl(args.stl)
    if args.json_out:
        write_json(args.json_out, part.to_dict())
    print(f"{part.family}: depth {hm.depth:.6g} m, steepest wall "
          f"{part.max_wall_angle_deg():.3f} deg, grid {grid.nx} x {grid.ny} at {grid.h:g} m")
    return 0


# -- material / setup ----------------------------------------------------------
def cmd_material(args: argparse.Namespace) -> int:
    from .materials import get_material, library

    if args.name is None:
        for name, mat in library().items():
            print(f"{name}: E = {mat.youngs_modulus:.4g} Pa, yield = {mat.yield_stress:.4g} Pa"
                  f" ({mat.source})")
        return 0
    mat = get_material(args.name)
    doc = mat.to_sparlab() if args.sparlab else mat.to_dict()
    print(canonical_json(doc), end="")
    return 0


def cmd_setup(args: argparse.Namespace) -> int:
    from .fea.setup import FormingSetup
    from .materials import Material, get_material

    material = (Material.from_dict(read_json(args.material_json)) if args.material_json
                else get_material(args.material))
    setup = FormingSetup(material, **_parse_assignments(args.set))
    write_json(args.out, setup.to_dict())
    print(f"wrote {args.out}")
    return 0


# -- toolpath ------------------------------------------------------------------
def cmd_toolpath(args: argparse.Namespace) -> int:
    from .toolpath import contour_toolpath, spiral_toolpath

    surface = _load_map(args.surface)
    radius, step, spacing, style = args.tool_radius, args.step_down, args.spacing, args.style
    if args.setup:
        s = _load_setup(args.setup)
        radius = radius or s.tool_radius
        step = step or s.step_down
        spacing = spacing or s.toolpath_spacing
        style = style or s.toolpath_style
    if not (radius and step):
        raise ValueError("give --tool-radius and --step-down, or --setup")
    spacing = spacing or 1e-3
    style = style or "spiral"
    fn = spiral_toolpath if style == "spiral" else contour_toolpath
    path = fn(surface, radius, step, spacing)
    if args.out:
        path.to_sparlab_csv(args.out)
    if args.robot:
        if not args.feed:
            raise ValueError("--robot needs --feed [m/s]")
        path.to_robot_csv(args.robot, args.feed)
    print(canonical_json(path.summary(args.feed)), end="")
    return 0


# -- simulate ------------------------------------------------------------------
def cmd_simulate(args: argparse.Namespace) -> int:
    from .fea.runner import simulate

    setup = _load_setup(args.setup)
    commanded = _load_map(args.commanded)
    res = simulate(setup, commanded, args.work_dir, cache=not args.no_cache,
                   retry_failed=args.retry_failed)
    step = int(args.step) if args.step.lstrip("-").isdigit() else args.step
    formed = res.formed_surface(step, grid=commanded.grid)
    if args.out:
        formed.save(args.out)
    print(canonical_json({"result": str(res.directory), "steps": res.step_names,
                          "provenance": to_jsonable({k: v for k, v in res.provenance.items()
                                                     if not isinstance(v, dict)}),
                          "formed_depth_m": formed.depth}), end="")
    return 0


# -- compensate ----------------------------------------------------------------
def cmd_compensate(args: argparse.Namespace) -> int:
    from .api import compensate
    from .report import write_report

    setup = _load_setup(args.setup)
    target = _load_map(args.target)
    model = None
    if args.method in ("surrogate", "hybrid"):
        if not args.model:
            raise ValueError(f"--method {args.method} needs --model DIR")
        registry = _import_ml("precomp.ml.registry")
        model = registry.load_model(args.model)
    result = compensate(target, setup, model, method=args.method, iterations=args.iterations,
                        alpha=args.alpha, verify=not args.no_verify, work_dir=args.work_dir,
                        smoothing=args.smoothing, direction=args.direction)
    out = result.save(args.out)
    from .metrology import signed_deviation
    write_report(out / "report", target, title="Compensation report",
                 deviation=signed_deviation(result.predicted, target),
                 formed=result.predicted, commanded=result.compensated,
                 tolerance=args.tolerance, history=result.history,
                 provenance={"method": result.method, "verification": result.verification,
                             "setup": setup.physics_dict()})
    print(f"wrote {out}")
    return 0


# -- scan ----------------------------------------------------------------------
def cmd_scan_compare(args: argparse.Namespace) -> int:
    from .api import compare_scan
    from .report import write_report

    target = _load_map(args.target)
    cmp = compare_scan(args.scan, target, scale=args.scale, align_mode=args.align,
                       fixture=args.fixture, tolerance=args.tolerance, method=args.method)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cmp.deviation.save(out / "deviation.npz")
    cmp.formed.save(out / "scan_gridded.npz")
    write_json(out / "comparison.json", cmp.summary())
    write_report(out, target, title=f"Scan comparison: {Path(args.scan).name}",
                 deviation=cmp.deviation, formed=cmp.formed, tolerance=args.tolerance,
                 provenance={"scan": str(args.scan), "scale": args.scale,
                             "alignment": cmp.summary()["alignment"]})
    print(canonical_json(cmp.summary()["metrics"]), end="")
    return 0


def cmd_scan_update(args: argparse.Namespace) -> int:
    from .compensation import update_from_scan

    target = _load_map(args.target)
    commanded = _load_map(args.commanded)
    new, report = update_from_scan(commanded, args.scan, target, alpha=args.alpha,
                                   scale=args.scale, align_mode=args.align,
                                   fixture=args.fixture, smoothing=args.smoothing)
    new.save(args.out)
    print(canonical_json(to_jsonable(report)), end="")
    return 0


# -- report --------------------------------------------------------------------
def cmd_report(args: argparse.Namespace) -> int:
    from .metrology import signed_deviation
    from .report import write_report

    target = _load_map(args.target)
    formed = _load_map(args.formed) if args.formed else None
    commanded = _load_map(args.commanded) if args.commanded else None
    dev = _load_map(args.deviation) if args.deviation else (
        signed_deviation(formed, target) if formed is not None else None)
    path = write_report(args.out, target, title=args.title, deviation=dev, formed=formed,
                        commanded=commanded, tolerance=args.tolerance)
    print(f"wrote {path}")
    return 0


# -- machine learning (delegated) -----------------------------------------------
def _import_ml(module: str = "precomp.ml"):
    try:
        return importlib.import_module(module)
    except ImportError as exc:
        raise PrecompError(f"the machine-learning commands need precomp.ml, which could not be "
                           f"imported ({exc}); install the package with its ML modules") from exc


def run_ml(argv: List[str]) -> int:
    """Hand an ML command line to `precomp.ml.cli.main(argv) -> int`."""
    try:
        ml_cli = importlib.import_module("precomp.ml.cli")
    except ImportError as exc:
        raise PrecompError(f"'precomp {argv[0]}' needs precomp.ml.cli, which could not be "
                           f"imported ({exc}); install the package with its ML modules") from exc
    if not hasattr(ml_cli, "main"):
        raise PrecompError("precomp.ml.cli has no main(argv) entry point")
    return int(ml_cli.main(list(argv)) or 0)


# -- parser ---------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="precomp", description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", action="version", version=f"precomp {__version__}")
    p.add_argument("--debug", action="store_true", help="show tracebacks")
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("part", help="make a target part")
    s.add_argument("--family", help="part family (see --list)")
    s.add_argument("--param", nargs="*", metavar="KEY=VALUE", help="family parameters [SI]")
    s.add_argument("--sample", action="store_true", help="draw the parameters at random")
    s.add_argument("--seed", type=int, default=0, help="seed for --sample")
    s.add_argument("--json-in", help="read the part from a JSON file (Part.to_dict)")
    s.add_argument("--list", action="store_true", help="list families and parameter bounds")
    s.add_argument("--size", type=float, default=0.2, help="grid side [m] (default 0.2)")
    s.add_argument("--spacing", type=float, default=5e-4, help="grid spacing [m] (default 5e-4)")
    s.add_argument("--out", help="height map .npz")
    s.add_argument("--stl", help="also write an STL surface")
    s.add_argument("--json-out", help="also write the part description as JSON")
    s.set_defaults(func=cmd_part)

    s = sub.add_parser("material", help="list or print library materials")
    s.add_argument("name", nargs="?", help="material name (omit to list)")
    s.add_argument("--sparlab", action="store_true", help="print the SparLab material block")
    s.set_defaults(func=cmd_material)

    s = sub.add_parser("setup", help="write a FormingSetup JSON")
    s.add_argument("--material", default="AA5754-O", help="library material name")
    s.add_argument("--material-json", help="material from a JSON file (Material.to_dict)")
    s.add_argument("--set", nargs="*", metavar="KEY=VALUE", help="override setup fields [SI]")
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_setup)

    s = sub.add_parser("toolpath", help="tool path of a surface")
    s.add_argument("--surface", "--target", dest="surface", required=True,
                   help="height map .npz (the commanded surface)")
    s.add_argument("--setup", help="take tool radius, step-down, spacing, style from a setup")
    s.add_argument("--tool-radius", type=float, help="[m]")
    s.add_argument("--step-down", type=float, help="[m]")
    s.add_argument("--spacing", type=float, help="point spacing [m] (default 1e-3)")
    s.add_argument("--style", choices=["spiral", "contour"])
    s.add_argument("--out", help="SparLab trajectory CSV (t,x,y,z)")
    s.add_argument("--robot", help="robot waypoint CSV (x,y,z,i,j,k,feed)")
    s.add_argument("--feed", type=float, help="feed rate [m/s]")
    s.set_defaults(func=cmd_toolpath)

    s = sub.add_parser("simulate", help="simulate forming with sparlab_form")
    s.add_argument("--setup", required=True)
    s.add_argument("--commanded", required=True, help="height map .npz")
    s.add_argument("--work-dir", required=True, help="run cache directory")
    s.add_argument("--step", default="-1", help="step whose surface to write (index or name)")
    s.add_argument("--out", help="formed surface .npz")
    s.add_argument("--no-cache", action="store_true")
    s.add_argument("--retry-failed", action="store_true")
    s.set_defaults(func=cmd_simulate)

    s = sub.add_parser("compensate", help="displacement-adjustment compensation")
    s.add_argument("--target", required=True)
    s.add_argument("--setup", required=True)
    s.add_argument("--method", choices=["fea", "surrogate", "hybrid"], default="fea")
    s.add_argument("--model", help="model directory (surrogate / hybrid; needs precomp.ml)")
    s.add_argument("--iterations", type=int, default=3)
    s.add_argument("--alpha", type=float, default=1.0)
    s.add_argument("--smoothing", type=float, help="update smoothing length [m]")
    s.add_argument("--direction", choices=["vertical", "normal"], default="vertical")
    s.add_argument("--tolerance", type=float, help="report tolerance [m]")
    s.add_argument("--work-dir", help="run cache directory")
    s.add_argument("--no-verify", action="store_true")
    s.add_argument("--out", required=True, help="output directory")
    s.set_defaults(func=cmd_compensate)

    scan = sub.add_parser("scan", help="scanned parts")
    ssub = scan.add_subparsers(dest="scan_command", required=True)
    s = ssub.add_parser("compare", help="align a scan and measure its deviation")
    s.add_argument("--scan", required=True, help=".xyz/.csv/.txt/.ply/.stl/.npy")
    s.add_argument("--scale", type=float, default=1.0, help="0.001 for a scan in mm")
    s.add_argument("--target", required=True)
    s.add_argument("--align", choices=["rigid", "translation_z", "none"], default="rigid")
    s.add_argument("--fixture", choices=["all", "flange"], default="all")
    s.add_argument("--method", choices=["grid", "cloud"], default="grid")
    s.add_argument("--tolerance", type=float, help="[m]")
    s.add_argument("--out", required=True, help="output directory")
    s.set_defaults(func=cmd_scan_compare)
    s = ssub.add_parser("update", help="one DA step of the commanded surface from a scan")
    s.add_argument("--scan", required=True)
    s.add_argument("--scale", type=float, default=1.0)
    s.add_argument("--commanded", required=True)
    s.add_argument("--target", required=True)
    s.add_argument("--alpha", type=float, default=1.0)
    s.add_argument("--smoothing", type=float)
    s.add_argument("--align", choices=["rigid", "translation_z", "none"], default="rigid")
    s.add_argument("--fixture", choices=["all", "flange"], default="flange")
    s.add_argument("--out", required=True, help="new commanded surface .npz")
    s.set_defaults(func=cmd_scan_update)

    s = sub.add_parser("report", help="Markdown report and figures")
    s.add_argument("--target", required=True)
    s.add_argument("--formed")
    s.add_argument("--commanded")
    s.add_argument("--deviation", help="a deviation .npz (default: computed from --formed)")
    s.add_argument("--tolerance", type=float, help="[m]")
    s.add_argument("--title", default="Part report")
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_report)

    for name, text in (("dataset", "simulated data sets (precomp.ml)"),
                       ("train", "train a surrogate (precomp.ml)"),
                       ("evaluate", "evaluate a surrogate (precomp.ml)")):
        s = sub.add_parser(name, help=text, add_help=False)
        s.add_argument("ml_args", nargs=argparse.REMAINDER)
        s.set_defaults(func=None)
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    debug = "--debug" in argv
    if debug:
        argv.remove("--debug")
    try:
        if argv and argv[0] in ML_COMMANDS:
            return run_ml(argv)
        args = build_parser().parse_args(argv)
        return int(args.func(args) or 0)
    except KeyboardInterrupt:
        return 130
    except Exception as exc:  # reported, not swallowed: message + exit status
        from .fea.runner import FormingError
        if debug:
            raise
        print(f"precomp: error: {exc}", file=sys.stderr)
        if isinstance(exc, FormingError):
            return 3
        if isinstance(exc, (PrecompError, ValueError, KeyError, FileNotFoundError, TypeError)):
            return 2
        raise


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
