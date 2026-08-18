"""Moved to :mod:`IMP.bff.photophysics.kappa2` by PRD-113 stage 4.

kappa^2 is a photophysical process parameter, not a property of the FRET engine.
It now sits beside the orientation *distributions* it feeds
(:mod:`IMP.bff.photophysics.orientation`), which lived in an unexported
``spectroscopy`` package and had no consumers at all.

Re-exported here so the call sites inside this package keep resolving while the
restructure lands; they move in the same stage as their own modules.
"""

from IMP.bff.photophysics.kappa2 import (  # noqa: F401
    kappa2_from_dipoles,
    kappa2_isotropic,
)

__all__ = ["kappa2_from_dipoles", "kappa2_isotropic"]
