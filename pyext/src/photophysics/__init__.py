"""``IMP.bff.photophysics`` -- the processes that deactivate or depolarise a dye.

Interaction terms with an arity, in the shape a force field has them: a
functional form plus parameters looked up by type.

* **1-body** -- radiative decay, internal conversion, rotational relaxation.
  These depend on the dye and the solvent, not on anything else in the system.
* **2-body** -- PET (a dye and a quencher), FRET (a dye and a dye). Their
  parameters are *pair* properties: ``kQ`` depends on the redox potentials of
  both partners, and R0 on both spectra plus the medium.
* **N-body** -- homo-FRET and multi-chromophore transfer, where the rate is not
  a sum of pairs.

Rate constants **add**, because the channels are parallel -- which is why they
are objects that compose rather than numbers computed inside an observable.

Everything here is **representation-agnostic**: a term consumes
:class:`~IMP.bff.representation.States` -- positions, weights, orientations -- so
one implementation serves an accessible volume, a rotamer library, a
coarse-grained model and an MD trajectory alike.
"""

from __future__ import annotations

from .kappa2 import kappa2_from_dipoles, kappa2_isotropic
from .terms import (
    FRETTerm, InteractionTerm, PETTerm, RadiativeTerm, total_rate,
)
from . import orientation

__all__ = [
    "kappa2_from_dipoles", "kappa2_isotropic", "orientation",
    "InteractionTerm", "RadiativeTerm", "PETTerm", "FRETTerm", "total_rate",
]
