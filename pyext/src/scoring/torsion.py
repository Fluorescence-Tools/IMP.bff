"""The torsion potential, and the sign convention that made it wrong once.

CHARMM writes a dihedral term as ``k(1 + cos(n*phi - delta))``; IMP's
``IMP.core.Cosine`` writes it as ``k(1 + cos(n*phi + delta))``. The two differ
by the sign of the phase, and getting it backwards makes a planar torsion
*minimal* at 90 degrees instead of 0 -- which is what PRD-108's Langevin sampler
found by way of a thermodynamic pin, not by inspection.

Moved out of ``cgdye/topology/dye.py`` in the four-stage restructure: a torsion
term is an energy, so it is scoring, and leaving it in the topology builder made
``scoring`` depend on a representation.
"""

from __future__ import annotations

import math

import IMP.core

__all__ = ["torsion_cosine"]


def torsion_cosine(torsion_type):
    """The ``IMP.core.Cosine`` for a torsion type stored in the CHARMM convention.

    cgdye's ``torsion_types`` (``k``, ``periodicity`` n, ``phase_rad`` δ) mean
    the CHARMM/AMBER form ``V = k (1 + cos(n φ − δ))``: ``T_PI`` (n = 2,
    δ = π) is minimal at the planar 0°/180°, ``T_LINK`` (n = 3, δ = 0) at
    the staggered ±60°/180°. ``IMP.core.Cosine(k, n, δ')`` scores
    ``k (1 − cos(n φ − δ'))`` — the *opposite* sign — so δ' = δ + π. Passing
    δ straight through (as the code did) put every conjugated torsion's
    minimum at 90° and every linker torsion's at the eclipsed 0°/±120°.
    """
    return IMP.core.Cosine(float(torsion_type["k"]), int(torsion_type["periodicity"]),
                           float(torsion_type["phase_rad"]) + math.pi)
