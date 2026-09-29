"""Figures and a Markdown report of a formed, simulated or scanned part.

Figures are human-facing, so their axes and colour bars are in millimetres,
always labelled "[mm]"; everything computed stays in metres. The colour rules
are those of `sparlab_viz.style` (validated there):

* signed deviation - `coolwarm`, limits symmetric about zero so the neutral
  midpoint is zero; nodes without data are left blank, never coloured;
* heights - `viridis` (sequential, monotone in lightness);
* several series - the categorical slots in fixed order, with a legend;
* one y-axis per chart: quantities of different scale go in separate figures.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

import matplotlib
matplotlib.use("Agg")  # figures are written to files, never shown

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402

from ._util import PathLike, canonical_json, to_jsonable  # noqa: E402
from .geometry.heightmap import HeightMap  # noqa: E402
from .metrology import metrics, region_masks  # noqa: E402

SURFACE = "#fcfcfb"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#8a8982"
GRID = "#e3e2dc"
#: Categorical slots in fixed order (the palette of `sparlab_viz.style`).
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7",
          "#e34948"]
SIGNED_CMAP = "coolwarm"
FIELD_CMAP = "viridis"
MM = 1e3  # metres -> millimetres, for labels only

_RC = {
    "figure.facecolor": SURFACE, "figure.dpi": 130, "savefig.dpi": 160,
    "savefig.facecolor": SURFACE, "savefig.bbox": "tight", "axes.facecolor": SURFACE,
    "axes.edgecolor": INK_MUTED, "axes.labelcolor": INK_SECONDARY,
    "axes.titlecolor": INK_PRIMARY, "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "axes.labelsize": 9.5, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "axes.axisbelow": True,
    "grid.color": GRID, "grid.linewidth": 0.7, "xtick.color": INK_SECONDARY,
    "ytick.color": INK_SECONDARY, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "legend.frameon": False, "legend.fontsize": 8.5, "lines.linewidth": 2.0,
}


def series_color(index: int) -> str:
    """Categorical colour of series `index` (fixed order, never cycled)."""
    if index >= len(SERIES):
        raise ValueError(f"series index {index} exceeds the {len(SERIES)}-slot palette; "
                         "use small multiples instead")
    return SERIES[index]


def _save(fig, path: Optional[PathLike]) -> None:
    if path is not None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path)


def _extent_mm(hm: HeightMap):
    x0, x1, y0, y1 = hm.grid.extent
    h = hm.grid.h / 2
    return [(x0 - h) * MM, (x1 + h) * MM, (y0 - h) * MM, (y1 + h) * MM]


def plot_deviation_map(dev: HeightMap, path: Optional[PathLike] = None, *,
                       title: str = "Signed normal deviation", limit: Optional[float] = None,
                       mask: Optional[np.ndarray] = None, label: Optional[str] = None):
    """Map of a deviation field [m], shown in mm with a zero-centred diverging scale.

    `limit` [m] sets the colour range +-limit (default: the 99th percentile
    of |d|, so a few outliers do not wash the map out; values beyond are
    clipped to the end colours). Nodes invalid in `dev` or outside `mask`
    are blank.
    """
    sel = dev.mask.copy() if mask is None else dev.mask & np.asarray(mask, bool)
    values = np.ma.masked_array(dev.z * MM, mask=~sel)
    if limit is None:
        a = np.abs(dev.z[sel])
        limit = float(np.percentile(a, 99)) if a.size else 1e-6
    lim = max(limit * MM, 1e-9)
    with plt.rc_context(_RC):
        fig, ax = plt.subplots(figsize=(5.6, 4.6))
        im = ax.imshow(values, origin="lower", extent=_extent_mm(dev), cmap=SIGNED_CMAP,
                       norm=TwoSlopeNorm(vcenter=0.0, vmin=-lim, vmax=lim),
                       interpolation="nearest")
        ax.grid(False)
        ax.set_aspect("equal")
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
        ax.set_title(title)
        cb = fig.colorbar(im, ax=ax, shrink=0.85, extend="both")
        cb.set_label(label or "deviation [mm] (+ = towards the tool)")
        cb.outline.set_edgecolor(INK_MUTED)
        _save(fig, path)
    return fig


def plot_height_map(hm: HeightMap, path: Optional[PathLike] = None, *,
                    title: str = "Tool-side surface"):
    """Map of a height field [m] in mm (sequential scale; invalid nodes blank)."""
    values = np.ma.masked_array(hm.z * MM, mask=~hm.mask)
    with plt.rc_context(_RC):
        fig, ax = plt.subplots(figsize=(5.6, 4.6))
        im = ax.imshow(values, origin="lower", extent=_extent_mm(hm), cmap=FIELD_CMAP,
                       interpolation="nearest")
        ax.grid(False)
        ax.set_aspect("equal")
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
        ax.set_title(title)
        cb = fig.colorbar(im, ax=ax, shrink=0.85)
        cb.set_label("z [mm]")
        cb.outline.set_edgecolor(INK_MUTED)
        _save(fig, path)
    return fig


def plot_histogram(dev: HeightMap, path: Optional[PathLike] = None, *,
                   mask: Optional[np.ndarray] = None, tolerance: Optional[float] = None,
                   title: str = "Deviation distribution", bins: int = 60):
    """Histogram of the deviations [m] (in mm) with the +-tolerance band marked."""
    sel = dev.mask.copy() if mask is None else dev.mask & np.asarray(mask, bool)
    v = dev.z[sel] * MM
    with plt.rc_context(_RC):
        fig, ax = plt.subplots(figsize=(5.6, 3.4))
        ax.hist(v, bins=bins, color=series_color(0), edgecolor=SURFACE, linewidth=0.6)
        if tolerance is not None:
            for s in (-1, 1):
                ax.axvline(s * tolerance * MM, color=INK_MUTED, linestyle="--", linewidth=1.2)
            ax.text(tolerance * MM, ax.get_ylim()[1] * 0.95, f" +/-{tolerance * MM:.3g} mm",
                    color=INK_SECONDARY, fontsize=8.5, va="top")
        ax.set_xlabel("deviation [mm]")
        ax.set_ylabel("grid nodes")
        ax.set_title(title)
        _save(fig, path)
    return fig


def plot_profile(surfaces: Mapping[str, HeightMap], path: Optional[PathLike] = None, *,
                 axis: str = "x", at: float = 0.0, title: Optional[str] = None):
    """Sections of several surfaces [m] along x (at y = `at`) or y (at x = `at`),
    in mm, one colour per surface in fixed order, with a legend."""
    if axis not in ("x", "y"):
        raise ValueError("axis must be 'x' or 'y'")
    with plt.rc_context(_RC):
        fig, ax = plt.subplots(figsize=(6.4, 3.4))
        for i, (name, hm) in enumerate(surfaces.items()):
            g = hm.grid
            if axis == "x":
                s = g.x
                z = hm.interpolate(s, np.full_like(s, at))
            else:
                s = g.y
                z = hm.interpolate(np.full_like(s, at), s)
            ax.plot(s * MM, z * MM, color=series_color(i), label=name)
        ax.set_xlabel(f"{axis} [mm]")
        ax.set_ylabel("z [mm]")
        other = "y" if axis == "x" else "x"
        ax.set_title(title or f"Section at {other} = {at * MM:.3g} mm")
        if len(surfaces) >= 2:
            ax.legend(loc="lower right")
        _save(fig, path)
    return fig


def plot_history(history: Sequence[Mapping[str, Any]], path: Optional[PathLike] = None, *,
                 title: str = "Compensation convergence"):
    """RMS error [m] (in mm, log scale) against the compensation iteration."""
    it = [h["iteration"] for h in history]
    rms = [h["error"]["rms"] * MM for h in history]
    with plt.rc_context(_RC):
        fig, ax = plt.subplots(figsize=(5.6, 3.4))
        ax.semilogy(it, rms, color=series_color(0), marker="o", markersize=5)
        ax.set_xlabel("iteration")
        ax.set_ylabel("RMS deviation over the part [mm]")
        ax.set_title(title)
        ax.set_xticks(it)
        _save(fig, path)
    return fig


def _fmt_mm(v: float) -> str:
    return f"{v * MM:.4f}"


def metrics_table(dev: HeightMap, target: HeightMap, tolerance: Optional[float] = None
                  ) -> List[Dict[str, Any]]:
    """Metrics [m] of `dev` over the target's regions (all, part, wall, flange);
    regions without valid nodes are skipped."""
    rows = []
    for name, m in region_masks(target).items():
        sel = m & dev.mask
        if not sel.any():
            continue
        row = {"region": name}
        row.update(metrics(dev, sel, tolerance))
        rows.append(row)
    return rows


def markdown_metrics(rows: Sequence[Mapping[str, Any]]) -> str:
    """A Markdown table of `metrics_table` rows, lengths in mm (labelled)."""
    has_tol = any("pct_within_tol" in r for r in rows)
    head = ("| Region | Nodes | RMS [mm] | Mean abs [mm] | Max abs [mm] | P95 abs [mm] "
            "| Bias [mm] |" + (" Within tolerance [%] |" if has_tol else ""))
    sep = "|---|---:|---:|---:|---:|---:|---:|" + ("---:|" if has_tol else "")
    lines = [head, sep]
    for r in rows:
        line = (f"| {r['region']} | {r['n']} | {_fmt_mm(r['rms'])} | {_fmt_mm(r['mae'])} | "
                f"{_fmt_mm(r['max_abs'])} | {_fmt_mm(r['p95_abs'])} | {_fmt_mm(r['bias'])} |")
        if has_tol:
            line += f" {r.get('pct_within_tol', float('nan')):.1f} |"
        lines.append(line)
    return "\n".join(lines)


def write_report(out_dir: PathLike, target: HeightMap, *, title: str = "Part report",
                 deviation: Optional[HeightMap] = None, formed: Optional[HeightMap] = None,
                 commanded: Optional[HeightMap] = None, tolerance: Optional[float] = None,
                 history: Optional[Sequence[Mapping[str, Any]]] = None,
                 notes: Optional[Sequence[str]] = None,
                 provenance: Optional[Mapping[str, Any]] = None) -> Path:
    """Write `report.md` and its figures into `out_dir`; return the report path.

    Contents: the metrics of `deviation` [m] per region (tables in mm), the
    deviation map and histogram, a section through the centre of the target,
    formed and commanded surfaces, the compensation history, free-text notes
    and the provenance as JSON. Every length in the report is labelled with
    its unit.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    lines = [f"# {title}", ""]
    lines.append(f"Target: {target.grid.nx} x {target.grid.ny} nodes at "
                 f"{target.grid.h * MM:.4g} mm spacing, depth {target.depth * MM:.4g} mm.")
    lines.append("")
    if deviation is not None:
        kind = deviation.metadata.get("quantity", "deviation")
        lines += ["## Deviation", "",
                  f"Quantity: `{kind}`; positive means the formed surface lies on the tool "
                  "side of the target." + (f" Tolerance: +/-{tolerance * MM:.4g} mm."
                                           if tolerance is not None else ""), ""]
        rows = metrics_table(deviation, target, tolerance)
        lines += [markdown_metrics(rows), ""]
        plt.close(plot_deviation_map(deviation, out / "deviation_map.png",
                                     mask=deviation.mask))
        plt.close(plot_histogram(deviation, out / "deviation_histogram.png",
                                 mask=region_masks(target)["part"] if rows else None,
                                 tolerance=tolerance,
                                 title="Deviation distribution over the part"))
        lines += ["![Deviation map](deviation_map.png)", "",
                  "![Deviation histogram over the part](deviation_histogram.png)", ""]
    surfaces = {"target": target}
    if formed is not None:
        surfaces["formed"] = formed
    if commanded is not None:
        surfaces["commanded"] = commanded
    plt.close(plot_profile(surfaces, out / "section_x.png", axis="x", at=0.0))
    lines += ["## Section through the centre", "", "![Section at y = 0](section_x.png)", ""]
    if history:
        plt.close(plot_history(history, out / "history.png"))
        lines += ["## Compensation history", "",
                  "| Iteration | RMS [mm] | Max abs [mm] | Update RMS [mm] |",
                  "|---:|---:|---:|---:|"]
        for h in history:
            lines.append(f"| {h['iteration']} | {_fmt_mm(h['error']['rms'])} | "
                         f"{_fmt_mm(h['error']['max_abs'])} | "
                         f"{_fmt_mm(h.get('update_rms_m', float('nan')))} |")
        lines += ["", "![Convergence](history.png)", ""]
    if notes:
        lines += ["## Notes", ""] + [f"- {n}" for n in notes] + [""]
    if provenance:
        lines += ["## Provenance", "", "```json", canonical_json(to_jsonable(dict(provenance)))
                  .rstrip(), "```", ""]
    path = out / "report.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path
