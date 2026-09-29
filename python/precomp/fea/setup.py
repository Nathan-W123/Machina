"""FormingSetup: everything about a forming run except the commanded shape.

A setup is JSON-serialisable (`to_dict` / `from_dict`) so that a data set can
record exactly how every sample was simulated. Fields that change the physics
go into the deck and hence into the content hash of a run; the fields that
only change how the run is executed (`executable`, `timeout`, `threads`) do
not.
"""

from __future__ import annotations

import dataclasses
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional

from .._util import PrecompError, require_nonnegative, require_positive
from ..materials import Material, get_material

#: Environment variable naming the sparlab_form executable.
EXECUTABLE_ENV = "PRECOMP_SPARLAB_FORM"

#: Default executable, relative to the working directory or the repository root.
DEFAULT_EXECUTABLE = "build/bin/sparlab_form"

#: Fields that only affect how a run is executed, not its result.
EXECUTION_FIELDS = ("executable", "timeout", "threads")


def repository_root() -> Optional[Path]:
    """The SparLab checkout this package lives in (an editable install), or None."""
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "CMakeLists.txt").is_file() and (parent / "python" / "precomp").is_dir():
            return parent
    return None


@dataclass
class FormingSetup:
    """Process, discretisation and solver settings of one SPIF simulation (SI).

    Blank and mesh
      blank_size : side of the square blank [m], centred on the origin.
      clamp_margin : width of the clamped frame along the blank edges [m]
          (at least one element).
      thickness : initial sheet thickness t [m]; the blank occupies
          z in [-t, 0], its tool side at z = 0.
      element_size : in-plane element size [m] (the blank gets
          round(blank_size / element_size) elements per side).
      layers : elements through the thickness.
      element : "hex8" (structured_hex) or "tet4" (structured_tet).
    Material
      material : a `Material`.
    Tool and path
      tool_radius [m], friction (Coulomb coefficient) [-], step_down [m],
      toolpath_style ("spiral" | "contour"), toolpath_spacing [m],
      toolpath_direction ("ccw" | "cw").
      contact : parameters passed verbatim into each tool's `contact` object
          (penalty / complementarity settings of the C++ contact; SI).
    Steps
      max_tool_travel : increment control, the largest tool-centre travel per
          increment [m].
      release : "321" - after the tool is withdrawn (step "unload", still
          clamped) the clamp is replaced by a statically determinate 3-2-1
          support (step "release"), so the part springs back freely;
          "clamped_only" - springback within the fixture only.
      kinematics : "finite" or "small_strain", the deck's nonlinear kinematics.
      solver : extra keys for the deck's `nonlinear` block (tolerances, ...).
    Execution (not part of the physics hash)
      executable : path of sparlab_form; None means $PRECOMP_SPARLAB_FORM,
          else build/bin/sparlab_form relative to the working directory or
          the repository root.
      timeout : wall-clock limit per run [s]; threads : OMP_NUM_THREADS.
    """

    material: Material
    blank_size: float = 0.20
    clamp_margin: float = 0.02
    thickness: float = 1.0e-3
    element_size: float = 2.5e-3
    layers: int = 2
    element: str = "hex8"
    tool_radius: float = 5.0e-3
    friction: float = 0.1
    step_down: float = 0.5e-3
    toolpath_style: str = "spiral"
    toolpath_spacing: float = 1.0e-3
    toolpath_direction: str = "ccw"
    contact: Dict[str, Any] = field(default_factory=dict)
    max_tool_travel: float = 1.0e-3
    release: str = "321"
    kinematics: str = "finite"
    solver: Dict[str, Any] = field(default_factory=dict)
    name: str = "spif"
    executable: Optional[str] = None
    timeout: float = 24 * 3600.0
    threads: int = 1

    def __post_init__(self) -> None:
        if isinstance(self.material, dict):
            self.material = Material.from_dict(self.material)
        elif isinstance(self.material, str):
            self.material = get_material(self.material)
        if not isinstance(self.material, Material):
            raise ValueError("material must be a Material, a Material dict or a library name")
        for key in ("blank_size", "thickness", "element_size", "tool_radius", "step_down",
                    "toolpath_spacing", "max_tool_travel", "timeout"):
            require_positive(key, getattr(self, key))
        require_positive("clamp_margin", self.clamp_margin)
        require_nonnegative("friction", self.friction)
        if not (isinstance(self.layers, int) and self.layers >= 1):
            raise ValueError(f"layers must be an integer >= 1, got {self.layers!r}")
        if not (isinstance(self.threads, int) and self.threads >= 1):
            raise ValueError(f"threads must be an integer >= 1, got {self.threads!r}")
        if self.element not in ("hex8", "tet4"):
            raise ValueError(f"element must be 'hex8' or 'tet4', got {self.element!r}")
        if self.toolpath_style not in ("spiral", "contour"):
            raise ValueError("toolpath_style must be 'spiral' or 'contour'")
        if self.toolpath_direction not in ("ccw", "cw"):
            raise ValueError("toolpath_direction must be 'ccw' or 'cw'")
        if self.release not in ("321", "clamped_only"):
            raise ValueError("release must be '321' or 'clamped_only'")
        if self.kinematics not in ("finite", "small_strain"):
            raise ValueError("kinematics must be 'finite' or 'small_strain'")
        if self.clamp_margin < self.element_size * (1 - 1e-9):
            raise ValueError("clamp_margin must be at least one element_size")
        if 2 * self.clamp_margin >= self.blank_size:
            raise ValueError("the clamped frame covers the whole blank")

    # -- derived -----------------------------------------------------------
    @property
    def elements_per_side(self) -> int:
        """In-plane elements per blank side, round(blank_size / element_size)."""
        return max(1, int(round(self.blank_size / self.element_size)))

    @property
    def meshed_blank_size(self) -> float:
        """Side of the meshed blank [m]: `blank_size` rounded to a whole number
        of elements (elements_per_side * element_size)."""
        return self.elements_per_side * self.element_size

    @property
    def free_half_width(self) -> float:
        """Half-width of the unclamped window [m]: meshed_blank_size / 2 - clamp_margin."""
        return 0.5 * self.meshed_blank_size - self.clamp_margin

    def resolved_executable(self) -> Path:
        """The sparlab_form executable this setup runs (see the class docstring).

        Raises PrecompError naming every place looked at when none exists.
        """
        tried = []
        candidates = []
        if self.executable:
            candidates.append(Path(self.executable))
        elif os.environ.get(EXECUTABLE_ENV):
            candidates.append(Path(os.environ[EXECUTABLE_ENV]))
        else:
            candidates.append(Path.cwd() / DEFAULT_EXECUTABLE)
            root = repository_root()
            if root is not None:
                candidates.append(root / DEFAULT_EXECUTABLE)
        for c in candidates:
            tried.append(str(c))
            if c.is_file() and os.access(c, os.X_OK):
                return c.resolve()
        raise PrecompError(
            "the SparLab forming executable was not found (looked at: " + ", ".join(tried)
            + f"); build it with `make build`, or set {EXECUTABLE_ENV} or "
            "FormingSetup.executable")

    # -- serialisation -----------------------------------------------------
    def to_dict(self) -> Dict[str, Any]:
        doc = {f.name: getattr(self, f.name) for f in dataclasses.fields(self)}
        doc["material"] = self.material.to_dict()
        doc["contact"] = dict(self.contact)
        doc["solver"] = dict(self.solver)
        return doc

    def physics_dict(self) -> Dict[str, Any]:
        """`to_dict` without the execution-only fields."""
        doc = self.to_dict()
        for key in EXECUTION_FIELDS:
            doc.pop(key)
        return doc

    @classmethod
    def from_dict(cls, doc: Dict[str, Any]) -> "FormingSetup":
        names = {f.name for f in dataclasses.fields(cls)}
        unknown = set(doc) - names
        if unknown:
            raise ValueError(f"FormingSetup: unknown keys {sorted(unknown)}")
        if "material" not in doc:
            raise ValueError("FormingSetup: 'material' is required")
        return cls(**doc)

    def replace(self, **changes: Any) -> "FormingSetup":
        return dataclasses.replace(self, **changes)

    def check_part_fits(self, footprint_radius: float) -> None:
        """Raise ValueError unless a part of that footprint radius [m], plus
        the tool radius, fits inside the unclamped window."""
        need = footprint_radius + self.tool_radius
        if need > self.free_half_width:
            raise ValueError(
                f"the part (footprint radius {footprint_radius:.4g} m plus tool radius "
                f"{self.tool_radius:.4g} m) does not fit in the unclamped window of "
                f"half-width {self.free_half_width:.4g} m")
