"""Readers for the artefacts a SparLab run writes.

Nothing here recomputes physics. Every loader is a thin, validated reader, so a
missing or malformed file produces a clear error instead of a silently empty
plot. Meshes are two- or three-dimensional; the dimension is read from
`mesh.json` and every accessor that only makes sense in one of the two says so.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd


class ResultError(RuntimeError):
    """Raised when a result directory is missing or inconsistent."""


def load_json(path: str) -> Dict[str, Any]:
    """Read a JSON document written by SparLab."""
    if not os.path.isfile(path):
        raise ResultError(f"{path} does not exist; run the corresponding case first")
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def load_csv(path: str, required: bool = True) -> Optional[pd.DataFrame]:
    """Read one of the CSV tables. Returns None when optional and absent."""
    if not os.path.isfile(path):
        if required:
            raise ResultError(f"{path} does not exist; run the corresponding case first")
        return None
    frame = pd.read_csv(path)
    if frame.empty:
        raise ResultError(f"{path} contains no rows")
    return frame


#: Local node lists of the six faces of a Hex8, in the VTK convention used by
#: the C++ mesh layer (each face wound so its normal points outward).
HEX_FACES = np.array([[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4],
                      [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]], dtype=int)

#: The four outward-wound faces of a positively oriented Tet4 (same table as
#: the C++ mesh layer).
TET_FACES = np.array([[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]], dtype=int)

#: The four 6-node faces of a Tet10: the Tet4 corners, then the edge nodes of
#: the face's corner pairs (0, 1), (1, 2), (2, 0) (same table as the C++ mesh
#: layer).
TET10_FACES = np.array([[0, 2, 1, 6, 5, 4], [0, 1, 3, 4, 8, 7],
                        [1, 2, 3, 5, 9, 8], [0, 3, 2, 7, 9, 6]], dtype=int)

#: A 6-node face as four triangles through its edge nodes, with its winding
#: (the split the STL writer uses), so a curved Tet10 face is drawn at the
#: resolution of its nodes.
TRI6_SPLIT = np.array([[0, 3, 5], [3, 1, 4], [5, 4, 2], [3, 4, 5]], dtype=int)

#: Face table of each solid element type.
SOLID_FACES = {"Hex8": HEX_FACES, "Tet4": TET_FACES, "Tet10": TET10_FACES}

#: Corner nodes per face, which alone identify a face shared by two cells.
FACE_CORNERS = {"Hex8": 4, "Tet4": 3, "Tet10": 3}


@dataclass
class Mesh:
    """Nodes, connectivity, prescribed DOFs and applied loads of one case."""

    nodes: np.ndarray          # (num_nodes, dim) [m]
    elements: np.ndarray       # (num_elements, nodes_per_element)
    element_type: str
    prescribed: List[Dict[str, Any]] = field(default_factory=list)
    load_cases: List[Dict[str, Any]] = field(default_factory=list)
    # Index into the summary's "materials" per element; None for one material.
    element_materials: Optional[np.ndarray] = None

    @property
    def dim(self) -> int:
        return int(self.nodes.shape[1])

    @property
    def num_nodes(self) -> int:
        return int(self.nodes.shape[0])

    @property
    def num_elements(self) -> int:
        return int(self.elements.shape[0])

    @property
    def extent(self) -> tuple:
        """(xmin, xmax, ymin, ymax[, zmin, zmax]) of the nodal coordinates [m]."""
        out = []
        for k in range(self.dim):
            out.extend((float(self.nodes[:, k].min()), float(self.nodes[:, k].max())))
        return tuple(out)

    @property
    def polygons(self) -> np.ndarray:
        """(num_elements, nodes_per_element, 2) array of element corner points.

        Only a plane mesh has element polygons; a solid mesh exposes its
        boundary faces through `boundary_faces` instead.
        """
        if self.dim != 2:
            raise ResultError("element polygons exist for 2-D meshes only; use "
                              "boundary_faces() for a solid mesh")
        return self.nodes[self.elements]

    @property
    def element_centroids(self) -> np.ndarray:
        return self.nodes[self.elements].mean(axis=1)

    @property
    def is_simplex(self) -> bool:
        """True for triangles and tetrahedra (linear or quadratic)."""
        return self.element_type in ("Tri3", "Tet4", "Tet10")

    def boundary_faces(self, mask: Optional[np.ndarray] = None,
                       return_owners: bool = False):
        """Outward-wound boundary faces of a solid mesh (or of a subset of it).

        Returns an (n_faces, k) array of node indices - quads (k = 4) for a
        Hex8 mesh, triangles (k = 3) for a Tet4 mesh, and for a Tet10 mesh
        each 6-node face split into four triangles through its edge nodes:
        the faces owned by exactly one element of the subset, wound so the
        right-hand normal points out of the material. `mask` selects the elements (all when
        None). With `return_owners`, also returns the index (into the full
        element list) of the element each face belongs to, which is what maps
        an element field onto the surface.
        """
        if self.dim != 3:
            raise ResultError("boundary faces are defined for 3-D meshes only")
        table = SOLID_FACES.get(self.element_type)
        if table is None:
            raise ResultError(f"no face table for element type {self.element_type}")
        k = table.shape[1]
        ids = (np.arange(self.num_elements) if mask is None
               else np.flatnonzero(np.asarray(mask, bool)))
        elements = self.elements[ids]
        if elements.size == 0:
            empty = np.empty((0, 3 if k == 6 else k), dtype=int)
            return (empty, np.empty(0, dtype=int)) if return_owners else empty
        faces = elements[:, table].reshape(-1, k)               # (f n, k)
        owners = np.repeat(ids, table.shape[0])
        keys = np.sort(faces[:, :FACE_CORNERS[self.element_type]], axis=1)
        _unique, first, counts = np.unique(keys, axis=0, return_index=True,
                                           return_counts=True)
        keep = first[counts == 1]
        faces, owners = faces[keep], owners[keep]
        if k == 6:
            faces = faces[:, TRI6_SPLIT].reshape(-1, 3)
            owners = np.repeat(owners, TRI6_SPLIT.shape[0])
        if return_owners:
            return faces, owners
        return faces

    def constrained_nodes(self, component: Optional[str] = None) -> np.ndarray:
        """Indices of nodes with a prescribed DOF, optionally for one component."""
        ids = [
            entry["node"]
            for entry in self.prescribed
            if component is None or entry["component"] == component
        ]
        return np.unique(np.asarray(ids, dtype=int)) if ids else np.empty(0, dtype=int)

    def nodal_forces(self, load_case: str) -> np.ndarray:
        """(k, 1 + dim) array of [node, fx, fy(, fz)] for the named load case [N]."""
        keys = ["fx_N", "fy_N", "fz_N"][: self.dim]
        for case in self.load_cases:
            if case["name"] != load_case:
                continue
            rows = [
                (entry["node"], *[entry[k] for k in keys])
                for entry in case["nodal_forces"]
            ]
            return (np.asarray(rows, dtype=float) if rows
                    else np.empty((0, 1 + self.dim)))
        raise ResultError(f"load case '{load_case}' is not present in mesh.json")

    @property
    def load_case_names(self) -> List[str]:
        return [case["name"] for case in self.load_cases]


def load_mesh(directory: str) -> Mesh:
    """Read `mesh.json` from a result directory."""
    doc = load_json(os.path.join(directory, "mesh.json"))
    nodes = np.asarray(doc["nodes_m"], dtype=float)
    elements = np.asarray(doc["elements"], dtype=int)
    if nodes.ndim != 2 or nodes.shape[1] not in (2, 3):
        raise ResultError("mesh.json: nodes_m must be an array of [x, y] or [x, y, z] rows")
    declared = int(doc.get("dim", nodes.shape[1]))
    if declared != nodes.shape[1]:
        raise ResultError(f"mesh.json: dim is {declared} but nodes carry "
                          f"{nodes.shape[1]} coordinates")
    if elements.ndim != 2:
        raise ResultError("mesh.json: elements must be a rectangular array")
    if elements.size and (elements.min() < 0 or elements.max() >= nodes.shape[0]):
        raise ResultError("mesh.json: connectivity references a node outside the mesh")
    materials = doc.get("element_materials")
    if materials is not None:
        materials = np.asarray(materials, dtype=int)
        if materials.shape != (elements.shape[0],):
            raise ResultError("mesh.json: element_materials must hold one index per element")
    return Mesh(
        nodes=nodes,
        elements=elements,
        element_type=doc.get("element_type", "Quad4"),
        prescribed=doc.get("prescribed_dofs", []),
        load_cases=doc.get("load_cases", []),
        element_materials=materials,
    )


def load_summary(directory: str) -> Dict[str, Any]:
    """Read `summary.json` from a result directory."""
    return load_json(os.path.join(directory, "summary.json"))


@dataclass
class CaseResults:
    """Everything one result directory holds, loaded lazily where it is large."""

    directory: str
    mesh: Mesh
    summary: Dict[str, Any]

    def path(self, name: str) -> str:
        return os.path.join(self.directory, name)

    def has(self, name: str) -> bool:
        return os.path.isfile(self.path(name))

    @property
    def dim(self) -> int:
        return self.mesh.dim

    # -- static fields -----------------------------------------------------
    def displacement(self, load_case: str) -> pd.DataFrame:
        return load_csv(self.path(f"displacement_{_safe(load_case)}.csv"))

    def temperature(self, load_case: str) -> pd.DataFrame:
        return load_csv(self.path(f"temperature_{_safe(load_case)}.csv"))

    def stress(self, load_case: str) -> pd.DataFrame:
        return load_csv(self.path(f"stress_{_safe(load_case)}.csv"))

    def reactions(self, load_case: str) -> pd.DataFrame:
        return load_csv(self.path(f"reactions_{_safe(load_case)}.csv"))

    # -- modal -------------------------------------------------------------
    def modes(self, tag: str = "") -> Optional[pd.DataFrame]:
        suffix = f"_{_safe(tag)}" if tag else ""
        return load_csv(self.path(f"modes{suffix}.csv"), required=False)

    def mode_shapes(self, tag: str = "") -> Optional[pd.DataFrame]:
        suffix = f"_{_safe(tag)}" if tag else ""
        return load_csv(self.path(f"mode_shapes{suffix}.csv"), required=False)

    # -- buckling ----------------------------------------------------------
    def buckling(self, tag: str) -> Optional[pd.DataFrame]:
        """Load factors of the buckling check: tag "solid" or "topology"."""
        return load_csv(self.path(f"buckling_{_safe(tag)}.csv"), required=False)

    def buckling_mode(self, tag: str, load_case: str, mode: int) -> Optional["VtkGrid"]:
        """The mesh and shape of one buckling mode (1-based `mode`), from the
        VTK file the run wrote (with output.vtk and output.mode_shapes)."""
        path = self.path(f"buckling_{_safe(tag)}_{_safe(load_case)}_{mode}.vtk")
        return read_legacy_vtk(path) if os.path.isfile(path) else None

    # -- topology ----------------------------------------------------------
    def history(self) -> Optional[pd.DataFrame]:
        return load_csv(self.path("history.csv"), required=False)

    def density(self) -> Optional[pd.DataFrame]:
        return load_csv(self.path("density_final.csv"), required=False)

    def density_history(self) -> Optional[pd.DataFrame]:
        return load_csv(self.path("density_history.csv"), required=False)

    @property
    def load_case_names(self) -> List[str]:
        return self.mesh.load_case_names

    @property
    def name(self) -> str:
        return self.summary.get("provenance", {}).get("case", os.path.basename(self.directory))

    @property
    def is_topology_run(self) -> bool:
        return "optimization_result" in self.summary


@dataclass
class VtkGrid:
    """An unstructured grid read from a legacy ASCII VTK file."""

    title: str
    points: np.ndarray                 # (num_points, 3)
    cells: List[np.ndarray]            # node lists, one array per cell
    cell_types: np.ndarray             # VTK cell type ids
    point_data: Dict[str, np.ndarray] = field(default_factory=dict)
    cell_data: Dict[str, np.ndarray] = field(default_factory=dict)


def read_legacy_vtk(path: str) -> VtkGrid:
    """Read the ASCII UNSTRUCTURED_GRID files SparLab writes: POINTS, CELLS,
    CELL_TYPES, and SCALARS / VECTORS blocks under POINT_DATA and CELL_DATA."""
    if not os.path.isfile(path):
        raise ResultError(f"{path} does not exist; run the corresponding case first")
    with open(path, "r", encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    if len(lines) < 4 or not lines[0].startswith("# vtk DataFile"):
        raise ResultError(f"{path} is not a legacy VTK file")
    if lines[2].strip().upper() != "ASCII":
        raise ResultError(f"{path} is not an ASCII VTK file")
    title = lines[1].strip()
    tokens: List[str] = []
    for line in lines[3:]:
        tokens.extend(line.split())
    pos = 0

    def take(count: int) -> List[str]:
        nonlocal pos
        out = tokens[pos:pos + count]
        if len(out) != count:
            raise ResultError(f"{path} ends inside a data block")
        pos += count
        return out

    points = np.zeros((0, 3))
    cells: List[np.ndarray] = []
    cell_types = np.zeros(0, dtype=int)
    point_data: Dict[str, np.ndarray] = {}
    cell_data: Dict[str, np.ndarray] = {}
    section: Optional[Dict[str, np.ndarray]] = None
    section_size = 0
    while pos < len(tokens):
        key = tokens[pos].upper()
        pos += 1
        if key == "DATASET":
            kind = take(1)[0].upper()
            if kind != "UNSTRUCTURED_GRID":
                raise ResultError(f"{path}: dataset {kind} is not supported")
        elif key == "POINTS":
            count = int(take(2)[0])
            points = np.array(take(3 * count), dtype=float).reshape(count, 3)
        elif key == "CELLS":
            count, size = (int(v) for v in take(2))
            flat = np.array(take(size), dtype=int)
            at = 0
            for _ in range(count):
                n = flat[at]
                cells.append(flat[at + 1:at + 1 + n])
                at += n + 1
        elif key == "CELL_TYPES":
            count = int(take(1)[0])
            cell_types = np.array(take(count), dtype=int)
        elif key in ("POINT_DATA", "CELL_DATA"):
            section_size = int(take(1)[0])
            section = point_data if key == "POINT_DATA" else cell_data
        elif key in ("SCALARS", "VECTORS"):
            if section is None:
                raise ResultError(f"{path}: {key} outside POINT_DATA / CELL_DATA")
            name = take(1)[0]
            take(1)  # data type
            components = 3 if key == "VECTORS" else 1
            if key == "SCALARS":
                # Optional component count, then the LOOKUP_TABLE line.
                if tokens[pos].upper() != "LOOKUP_TABLE":
                    components = int(take(1)[0])
                take(2)
            values = np.array(take(section_size * components), dtype=float)
            section[name] = values.reshape(section_size, components) if components > 1 \
                else values
        else:
            raise ResultError(f"{path}: unexpected keyword {key}")
    return VtkGrid(title, points, cells, cell_types, point_data, cell_data)


def load_case(directory: str) -> CaseResults:
    """Load the mesh and summary of a result directory."""
    if not os.path.isdir(directory):
        raise ResultError(f"{directory} is not a directory")
    return CaseResults(
        directory=directory,
        mesh=load_mesh(directory),
        summary=load_summary(directory),
    )


def _safe(name: str) -> str:
    """Mirror the file-name sanitisation done by the C++ ResultWriter."""
    return "".join(
        c if (c.isalnum() or c in "-_") else "_" for c in name
    ) or "unnamed"
