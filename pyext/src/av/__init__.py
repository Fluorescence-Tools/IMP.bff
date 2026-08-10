"""IMP.bff.av — Accessible Volume computation and analysis.

This module provides Accessible Volume (AV) functionality as an extension
to ``IMP.bff``.  It offers the following:

* :class:`BasicAV` — point-cloud-based AV with distance calculations
* :class:`ACV` — Accessible Contact Volume with trapped fraction
* :func:`compute_av` — standalone AV computation from raw coordinates
* :mod:`IMP.bff.av._kernels` — Numba-accelerated distance kernels

Usage
-----
>>> from IMP.bff.av import BasicAV, ACV, compute_av
"""

from __future__ import annotations

from IMP.bff.av.basic import BasicAV
from IMP.bff.av.acv import ACV
from IMP.bff.av.compute import compute_av, AccessibleVolume

__all__ = [
    "BasicAV",
    "ACV",
    "compute_av",
    "AccessibleVolume",
]
