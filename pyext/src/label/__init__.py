"""``IMP.bff.label`` -- what is attached where, and what quenches it.

The *system* layer: a :class:`Label` is a dye at a position on a structure, a
:class:`Quencher` is a moiety that deactivates it. Neither says how the dye's
positions are enumerated -- that is ``IMP.bff.representation`` -- nor what the
rates are, which is ``IMP.bff.photophysics``.

Names follow the FLR dictionaries where an item exists, checked against
``../mmfdb/src/mmfdb/data/*.dic`` rather than recalled: ``_flr_poly_probe_position``
for the position, ``_flr_sample_probe_details`` for what the probe is doing.
Quenching has no item anywhere in the stack, so those names are bff-native.

.. note::
   ``LabelDistribution`` and friends used to live here and now do not. A label
   *distribution* is one way of representing where the dye can be -- a
   representation -- while a :class:`Label` says which dye is attached where.
   The same word meant both; they moved to ``IMP.bff.representation`` in
   PRD-113 stage 3d and this package keeps only the system meaning.
"""

from __future__ import annotations

from .site import Label, FLUOROPHORE_TYPES
from .quencher import (
    Quencher, PETParameters, reference_quenchers,
    reference_pet_parameters, REFERENCE_DYE,
)

__all__ = [
    "Label",
    "Quencher",
    "PETParameters",
    "reference_quenchers",
    "reference_pet_parameters",
    "REFERENCE_DYE",
    "FLUOROPHORE_TYPES",
]
