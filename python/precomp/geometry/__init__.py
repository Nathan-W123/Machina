"""Geometry: the tool-side surface as a height field, STL input/output, part families.

See `heightmap` for the surface convention (flat blank at z = 0, formed part
z <= 0, arrays indexed [y, x] with x fastest) and `parts` for the parametric
target families.
"""

from .heightmap import Grid, HeightMap, zeros
from .parts import (
    MAX_WALL_ANGLE_DEG,
    Dome,
    EllipticCone,
    FilletedProfile,
    Freeform,
    Part,
    Pyramid,
    Saddle,
    TruncatedCone,
    TwoLevel,
    ellipse_distance,
    families,
    part_from_dict,
    rounded_rectangle_distance,
)
from .stl import raycast_top, read_stl, write_stl

__all__ = [
    "Grid", "HeightMap", "zeros",
    "MAX_WALL_ANGLE_DEG", "Part", "TruncatedCone", "Pyramid", "Dome", "EllipticCone",
    "TwoLevel", "Freeform", "Saddle", "FilletedProfile", "ellipse_distance",
    "rounded_rectangle_distance", "families", "part_from_dict",
    "read_stl", "write_stl", "raycast_top",
]
