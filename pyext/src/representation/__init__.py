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
properties of the dye (``IMP.bff.dye``) or of where it is attached
(``IMP.bff.label``).
"""

from __future__ import annotations

from .pathmap import PathMapReading, resample_av

__all__ = ["PathMapReading", "resample_av"]
