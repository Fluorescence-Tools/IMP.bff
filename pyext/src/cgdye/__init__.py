"""Coarse-Grained Dye Simulations for IMP."""

from .rotamer import RotamerFRET
from .system import System
from .labeling import attach_dyes, resolve_site, align_hierarchies, backbone_frame, place_dye_from_coords
from .sampling.rotamer import apply_rotamer_coords, sample_rotamer_index
from .analysis.fret import compute_exact_efficiency, calculate_fret_exact
from .io import read_ff_system, write_ff_system, read_cgdye_template

__all__ = [
    "RotamerFRET",
    "System",
    "attach_dyes",
    "resolve_site",
    "align_hierarchies",
    "backbone_frame",
    "place_dye_from_coords",
    "apply_rotamer_coords",
    "sample_rotamer_index",
    "compute_exact_efficiency",
    "calculate_fret_exact",
    "read_ff_system",
    "write_ff_system",
    "read_cgdye_template",
]
