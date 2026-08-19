"""``IMP.bff.dye`` -- the dye as a species, independent of how it is modelled.

A dye's spectra, quantum yield, extinction coefficient, transition dipole and
formal charges are what it *is*; an accessible volume, a rotamer library, a
coarse-grained conformer set and an MD trajectory are ways of *representing* it.
Only the former lives here -- see ``IMP.bff.representation`` for the latter.

The practical consequence is that the Förster radius is derived rather than
supplied: :func:`forster_radius` takes two dyes, kappa^2 and the medium's
refractive index, so no call site needs to carry ``forster_radius=52.0``.
"""

from __future__ import annotations

from .species import Dye, Spectrum, FLRCIF_ITEMS
from IMP.bff.dye.library import find_dye, available_dyes
from .cif import read_dye_library, DYE_LIBRARY_CIF
from .spectra import (
    DEFAULT_REFRACTIVE_INDEX,
    forster_radius,
    forster_radius_from_spectra,
    spectral_overlap,
    normalize_dye_name,
)

__all__ = [
    "Dye", "Spectrum", "FLRCIF_ITEMS",
    "find_dye", "available_dyes", "read_dye_library", "DYE_LIBRARY_CIF",
    "forster_radius", "spectral_overlap", "forster_radius_from_spectra",
    "normalize_dye_name", "DEFAULT_REFRACTIVE_INDEX",
]
