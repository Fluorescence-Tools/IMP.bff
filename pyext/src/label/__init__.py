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
   ``LabelDistribution`` and friends are re-exported here for now but they are a
   **different sense of the word**: a label *distribution* is one way of
   representing where the dye can be, i.e. a representation, while a
   :class:`Label` is a system entity saying which dye is attached where. They
   move to ``IMP.bff.representation`` in PRD-113 stage 3, at which point this
   package keeps only the system meaning.
"""

from __future__ import annotations

from .site import Label, FLUOROPHORE_TYPES
from .quencher import (
    Quencher, PETParameters, reference_quenchers,
    reference_pet_parameters, REFERENCE_DYE,
)
from IMP.bff.label.distribution import (
    LabelDistribution,
    LabelDistributionAV,
    DyeDistributionNormal,
)

__all__ = [
    "Label",
    "Quencher",
    "PETParameters",
    "reference_quenchers",
    "reference_pet_parameters",
    "REFERENCE_DYE",
    "FLUOROPHORE_TYPES",
    # moving to IMP.bff.representation in stage 3
    "LabelDistribution",
    "LabelDistributionAV",
    "DyeDistributionNormal",
]
