"""IMP.bff.label — Dye label distribution models.

Provides :class:`LabelDistributionAV` (accessible-volume based) and
:class:`DyeDistributionNormal` (Gaussian) for modelling the 3-D
positional distribution of fluorescent labels.

These classes build on the AV infrastructure in ``IMP.bff.av`` and
are designed as non-breaking additions to ``IMP.bff``.
"""

from __future__ import annotations

from IMP.bff.label.distribution import (
    LabelDistribution,
    LabelDistributionAV,
    DyeDistributionNormal,
)

__all__ = [
    "LabelDistribution",
    "LabelDistributionAV",
    "DyeDistributionNormal",
]
