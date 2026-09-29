"""STL input and output for height fields.

`read_stl` reads binary and ASCII STL into an (m, 3, 3) array of triangle
corners; `write_stl` writes a binary (default) or ASCII file with facet normals
from the right-hand winding. `raycast_top` turns any triangulated surface into
a height field on a grid by casting a ray along -z through each node and
keeping the highest hit - the tool-side surface of a part, whatever the
triangulation. Lengths are whatever the file holds; the callers apply a unit
`scale` explicitly (`HeightMap.from_stl`).
"""

from __future__ import annotations

import struct
from pathlib import Path
from typing import TYPE_CHECKING, Tuple

import numpy as np

from .._util import PathLike, PrecompError

if TYPE_CHECKING:  # pragma: no cover
    from .heightmap import Grid, HeightMap

#: Triangle pairs per chunk of the ray cast: bounds its memory (~100 MB).
_RAYCAST_CHUNK_PAIRS = 2_000_000


def heightmap_triangles(hm: "HeightMap") -> Tuple[np.ndarray, np.ndarray]:
    """Vertices (n, 3) [m] and triangles (m, 3) of the valid cells of `hm`.

    Cell (i, j) with all four nodes valid is split along its (i, j)-(i+1, j+1)
    diagonal into two triangles wound counter-clockwise seen from +z, so the
    right-hand normals point towards the tool.
    """
    g = hm.grid
    X, Y = g.mesh()
    vertices = np.stack([X.ravel(), Y.ravel(), hm.z.ravel()], axis=-1)
    idx = np.arange(g.nx * g.ny).reshape(g.shape)
    n00 = idx[:-1, :-1].ravel()
    n10 = idx[:-1, 1:].ravel()
    n01 = idx[1:, :-1].ravel()
    n11 = idx[1:, 1:].ravel()
    m = hm.mask
    cell_ok = (m[:-1, :-1] & m[:-1, 1:] & m[1:, :-1] & m[1:, 1:]).ravel()
    t1 = np.stack([n00, n10, n11], axis=-1)[cell_ok]
    t2 = np.stack([n00, n11, n01], axis=-1)[cell_ok]
    triangles = np.concatenate([t1, t2], axis=0)
    if len(triangles) == 0:
        raise PrecompError("the height map has no cell with four valid nodes")
    return vertices, triangles


def _facet_normals(tri: np.ndarray) -> np.ndarray:
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    length = np.linalg.norm(n, axis=1, keepdims=True)
    return np.divide(n, length, out=np.zeros_like(n), where=length > 0)


def write_stl(path: PathLike, vertices: np.ndarray, triangles: np.ndarray,
              binary: bool = True, header: str = "precomp") -> Path:
    """Write triangles (m, 3) indexing `vertices` (n, 3) as an STL file.

    Binary STL stores 32-bit floats: coordinates in metres keep about seven
    significant digits (0.1 um on a 1 m part). ASCII is written with 17
    significant digits and round-trips exactly.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tri = np.asarray(vertices, dtype=float)[np.asarray(triangles, dtype=int)]
    normals = _facet_normals(tri)
    if binary:
        rec = np.zeros(len(tri), dtype=[("n", "<f4", (3,)), ("v", "<f4", (3, 3)),
                                        ("attr", "<u2")])
        rec["n"] = normals
        rec["v"] = tri
        head = header.encode("ascii", "replace")[:80].ljust(80, b" ")
        with open(path, "wb") as handle:
            handle.write(head)
            handle.write(struct.pack("<I", len(tri)))
            handle.write(rec.tobytes())
    else:
        lines = [f"solid {header}"]
        for n, t in zip(normals, tri):
            lines.append(f"  facet normal {n[0]:.9e} {n[1]:.9e} {n[2]:.9e}")
            lines.append("    outer loop")
            for v in t:
                lines.append(f"      vertex {v[0]:.17g} {v[1]:.17g} {v[2]:.17g}")
            lines.append("    endloop")
            lines.append("  endfacet")
        lines.append(f"endsolid {header}")
        path.write_text("\n".join(lines) + "\n", encoding="ascii")
    return path


def read_stl(path: PathLike) -> np.ndarray:
    """Read a binary or ASCII STL file into an (m, 3, 3) float array of corners.

    A file is binary when its size is exactly 84 + 50 m bytes for the facet
    count m in its header; otherwise it must be ASCII ('solid ... facet ...
    vertex ...'). Facet normals in the file are ignored.
    """
    path = Path(path)
    if not path.is_file():
        raise PrecompError(f"{path} does not exist")
    data = path.read_bytes()
    if len(data) >= 84:
        (count,) = struct.unpack("<I", data[80:84])
        if len(data) == 84 + 50 * count:
            rec = np.frombuffer(data, dtype=[("n", "<f4", (3,)), ("v", "<f4", (3, 3)),
                                             ("attr", "<u2")], count=count, offset=84)
            return rec["v"].astype(float)
    text = data.decode("ascii", errors="replace")
    if not text.lstrip().lower().startswith("solid"):
        raise PrecompError(f"{path} is neither a binary STL (size mismatch) nor an ASCII STL")
    values = []
    for line in text.splitlines():
        parts = line.split()
        if parts and parts[0].lower() == "vertex":
            if len(parts) != 4:
                raise PrecompError(f"{path}: malformed vertex line '{line.strip()}'")
            values.append([float(p) for p in parts[1:]])
    if not values or len(values) % 3:
        raise PrecompError(f"{path}: expected a multiple of three vertices, found {len(values)}")
    return np.asarray(values, dtype=float).reshape(-1, 3, 3)


def raycast_top(triangles: np.ndarray, grid: "Grid") -> Tuple[np.ndarray, np.ndarray]:
    """Top-most intersection of a vertical ray through every grid node.

    Parameters
    ----------
    triangles : (m, 3, 3) corners [m].
    grid : the nodes to cast through.

    Returns
    -------
    z : (ny, nx) height of the highest triangle above each node [m] (-inf
        where nothing is hit).
    hit : (ny, nx) bool.

    Each triangle is tested against the nodes inside its x-y bounding box
    with barycentric coordinates (a node on a shared edge or vertex is inside
    both neighbours, to a relative 1e-6 - enough for the float32 corners of a
    binary STL); triangles seen edge-on from above
    (zero x-y area) are skipped - their edges belong to neighbouring facets.
    The work is vectorised over (triangle, node) pairs in bounded chunks.
    """
    tri = np.asarray(triangles, dtype=float)
    if tri.ndim != 3 or tri.shape[1:] != (3, 3):
        raise ValueError("triangles must be an (m, 3, 3) array")
    h = grid.h
    zbuf = np.full(grid.nx * grid.ny, -np.inf)
    a, b, c = tri[:, 0], tri[:, 1], tri[:, 2]
    v0 = b[:, :2] - a[:, :2]
    v1 = c[:, :2] - a[:, :2]
    den = v0[:, 0] * v1[:, 1] - v1[:, 0] * v0[:, 1]
    scale = np.maximum(np.abs(v0).max(axis=1), np.abs(v1).max(axis=1))
    keep = np.abs(den) > 1e-12 * np.maximum(scale, 1e-300) ** 2
    xmin = tri[:, :, 0].min(axis=1)
    xmax = tri[:, :, 0].max(axis=1)
    ymin = tri[:, :, 1].min(axis=1)
    ymax = tri[:, :, 1].max(axis=1)
    eps = 1e-6
    i0 = np.clip(np.ceil((xmin - grid.x0) / h - eps), 0, grid.nx - 1).astype(np.int64)
    i1 = np.clip(np.floor((xmax - grid.x0) / h + eps), -1, grid.nx - 1).astype(np.int64)
    j0 = np.clip(np.ceil((ymin - grid.y0) / h - eps), 0, grid.ny - 1).astype(np.int64)
    j1 = np.clip(np.floor((ymax - grid.y0) / h + eps), -1, grid.ny - 1).astype(np.int64)
    wi = np.maximum(i1 - i0 + 1, 0)
    wj = np.maximum(j1 - j0 + 1, 0)
    wi = np.where((xmax < grid.x0) | (xmin > grid.xmax), 0, wi)
    wj = np.where((ymax < grid.y0) | (ymin > grid.ymax), 0, wj)
    counts = np.where(keep, wi * wj, 0)
    ids = np.flatnonzero(counts)
    if len(ids) == 0:
        return zbuf.reshape(grid.shape), np.zeros(grid.shape, dtype=bool)
    csum = np.cumsum(counts[ids])
    start = 0
    while start < len(ids):
        base = csum[start - 1] if start else 0
        stop = int(np.searchsorted(csum, base + _RAYCAST_CHUNK_PAIRS, side="right"))
        stop = max(stop, start + 1)
        sel = ids[start:stop]
        cnt = counts[sel]
        t = np.repeat(sel, cnt)
        offs = np.repeat(np.cumsum(cnt) - cnt, cnt)
        k = np.arange(int(cnt.sum())) - offs
        ii = i0[t] + k % wi[t]
        jj = j0[t] + k // wi[t]
        px = grid.x0 + ii * h - a[t, 0]
        py = grid.y0 + jj * h - a[t, 1]
        d = den[t]
        u = (px * v1[t, 1] - v1[t, 0] * py) / d
        v = (v0[t, 0] * py - px * v0[t, 1]) / d
        tol = 1e-6
        inside = (u >= -tol) & (v >= -tol) & (u + v <= 1.0 + tol)
        zz = a[t, 2] + u * (b[t, 2] - a[t, 2]) + v * (c[t, 2] - a[t, 2])
        flat = (jj * grid.nx + ii)[inside]
        np.maximum.at(zbuf, flat, zz[inside])
        start = stop
    z = zbuf.reshape(grid.shape)
    return z, np.isfinite(z)
