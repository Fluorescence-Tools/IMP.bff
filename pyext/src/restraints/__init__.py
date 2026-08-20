"""IMP.bff.restraints — scoring restraints for FRET and labeling data.

Three restraints, all scoring a *structure* against derived fluorescence
observables, which is what IMP.bff is for. One is Python, because what it is
*for* is IMP.pmi: :class:`AVNetworkRestraintWrapper`, the fps.json-driven AV
network.

The other two are C++ and reachable as ``IMP.bff.SimpleAVNetworkRestraint`` and
``IMP.bff.DirectLabelingRestraint`` -- chi-squared over a network of AV-pair
distances without an fps.json, and chi-squared over attachment-atom distances
with no volume at all. They are re-exported here because this is still where
the name says they are.

They arrived from imp-tricks in 2026-08, shipped there as additions to this
package: legal under the add-never-replace rule, but it put half of one
namespace in another repository -- and a subpackage resolves to the first
matching directory and stops, so for five days ``AVNetworkRestraintWrapper``
was unreachable and recorded as a deleted upstream API. Everything under
``IMP.bff`` now ships from imp.bff.
"""

from __future__ import annotations

from IMP.bff import (
    AVMeasurement, DirectLabelingRestraint, LabelingSite,
    SimpleAVNetworkRestraint,
)

from .network import AVNetworkRestraintWrapper

__all__ = [
    "AVNetworkRestraintWrapper",
    "SimpleAVNetworkRestraint",
    "AVMeasurement",
    "DirectLabelingRestraint",
    "LabelingSite",
]
