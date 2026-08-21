"""``IMP.bff.representation`` -- how a dye's configuration space is described.

A dye can be represented as an accessible volume (a region plus a mobility
field), a rotamer library (listed states with weights), a coarse-grained or
all-atom molecular model (coordinates plus a force field), or a Gaussian. Each
supplies the same thing to everything downstream -- **states**: positions,
orientations and weights -- so distances, kappa^2 and the interaction terms are
written once and work for all of them.

Only the *representation* parameters live here. ``linker_length``,
``allowed_sphere_radius`` and ``simulation_grid_resolution`` belong to the
accessible volume; a rotamer library has none of them, and none of them are
properties of the dye (``IMP.bff.Dye``) or of where it is attached
(``IMP.bff.label``).

Port status (PRD-117): 15 of 17 names are C++-covered (the rotamer ensemble,
``RotamerDistance``/``RotamerPosition``, ``distances_from_ensembles``,
``load_rotamer_library``/``*_fps``, ``compute_av``, ``compute_av_from_structure``,
``compute_avs_for_structure``, etc.); the remaining 2 (``av_pair_statistics``,
``histogram_rda``) are %pythoncode in ``avmodel.i`` and still to port.
Pure re-export.
"""

from __future__ import annotations

from IMP.bff import (  # noqa: F401
    ACV, AccessibleVolume, DyeDistributionNormal, LabelDistribution,
    LabelDistributionAV, States, av_pair_statistics, compute_av,
    compute_av_from_structure, compute_avs_for_structure,
    distance_between_mean_positions, histogram_rda, resample_av,
    standard_deviation_of_distances,
    states_average_distance as average_distance,
    states_mean_fret_distance as mean_fret_distance,
)

__all__ = [
    "States", "AccessibleVolume", "ACV", "resample_av",
    "compute_av", "compute_av_from_structure", "compute_avs_for_structure",
    "LabelDistribution", "LabelDistributionAV", "DyeDistributionNormal",
    "av_pair_statistics", "average_distance",
    "distance_between_mean_positions", "histogram_rda",
    "mean_fret_distance", "standard_deviation_of_distances",
]
