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

Attachment came the other way. :mod:`~IMP.bff.label.attachment` and
:mod:`~IMP.bff.label.backbone_frame` -- putting a dye on a residue, and the
local frame that orients it -- were ``cgdye/labeling/`` until the 2026-08-18
cleanup. Attaching a label *is* the system layer; it was filed under the
explicit-dye package because that is what first needed it.
"""

from __future__ import annotations

from .attachment import (
    attach_dyes, resolve_dye_site, align_hierarchies,
    place_dye_from_coords, place_dye_from_rotamer_cb,
    strip_sidechain_at_site, SITE_KEEP_ATOM_NAMES,
)
from .backbone_frame import (
    backbone_frame, backbone_frame_from_coords,
    backbone_transformation, backbone_transformation_from_coords,
)
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
    "attach_dyes", "resolve_dye_site", "align_hierarchies",
    "place_dye_from_coords", "place_dye_from_rotamer_cb",
    "strip_sidechain_at_site", "SITE_KEEP_ATOM_NAMES",
    "backbone_frame", "backbone_frame_from_coords",
    "backbone_transformation", "backbone_transformation_from_coords",
]
