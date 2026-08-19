"""The accessible volume: where a tethered dye can reach.

The implicit representation — a region the dye may occupy, uniform on that
region, with a mobility field over it — against the rotamer library's discrete
states and the coarse-grained dye's explicit coordinates. All three are
representations, which is why they are siblings here rather than three
top-level packages.

Two front doors onto the same C++ core, and **both are needed**:

* :func:`~IMP.bff.representation.av.compute.compute_av` takes **arrays** —
  coordinates and radii — and knows nothing about files. This is the entry point
  for a caller that already has a structure in memory.
* :mod:`~IMP.bff.representation.av.structure` takes a **structure and an fps
  position**, which is what a labelling file describes.

They grew independently (imp-tricks and ChiSurf), and each carries fixes the
other lacks: the array door has the AV-on-its-own-particle handling and the
build/resample lock split, the structure door has the ``disc_step`` trap and
source clearance. PRD-113 stage 3 set out to merge them and did not; what it did
do is put them side by side, which is the precondition.

:class:`~IMP.bff.representation.av.basic.BasicAV` is the point cloud with its
distance methods, and :class:`~IMP.bff.representation.av.acv.ACV` splits it into
contact and free volumes for a dye that sticks.
"""

from __future__ import annotations

from IMP.bff.representation.av.acv import ACV
from IMP.bff.representation.av.basic import BasicAV
from IMP.bff.representation.av.compute import AccessibleVolume, compute_av

__all__ = ["ACV", "AccessibleVolume", "BasicAV", "compute_av"]
