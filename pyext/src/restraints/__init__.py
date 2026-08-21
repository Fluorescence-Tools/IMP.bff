"""IMP.bff.restraints — scoring restraints for FRET and labeling data.

**Nothing is defined here.** Every restraint in this package is C++ and lives on
the flat surface; this module re-exports the four names under the path they have
always had, so `IMP.bff.restraints.AVNetworkRestraintWrapper` keeps resolving.

- :class:`AVNetworkRestraintWrapper` — the fps.json-driven AV network, wrapped
  for IMP.pmi. The wrapper is `%pythoncode` in `pyext/IMP_bff.avmeandistance.i`
  and is built **lazily**, because it subclasses `IMP.pmi.restraints.RestraintBase`
  and IMP.pmi is not one of this module's `required_modules`.
- :class:`SimpleAVNetworkRestraint` — chi-squared over a network of AV-pair
  distances, with no fps.json.
- :class:`DirectLabelingRestraint` — chi-squared over attachment-atom distances,
  with no volume at all.
- :class:`AVMeanDistanceRestraint` — the separation of two mean dye positions,
  with derivatives.

They arrived from imp-tricks in 2026-08, shipped there as additions to this
package: legal under the add-never-replace rule, but it put half of one
namespace in another repository -- and a subpackage resolves to the first
matching directory and stops, so for five days ``AVNetworkRestraintWrapper``
was unreachable and recorded as a deleted upstream API. Everything under
``IMP.bff`` now ships from imp.bff.

Port status (PRD-117): 6 of 7 names are C++-covered (the restraints and
``LabelingSite``/``AVMeasurement``); only ``AVNetworkRestraintWrapper`` remains
%pythoncode (lazily) in ``avmeandistance.i`` because it subclasses an IMP.pmi
class. Pure re-export -- holds no logic of its own to move.
"""

from __future__ import annotations

from IMP.bff import (
    AVMeanDistanceRestraint, AVMeasurement, AVNetworkRestraintWrapper,
    DirectLabelingRestraint, LabelingSite, SimpleAVNetworkRestraint,
)

__all__ = [
    "AVMeanDistanceRestraint",
    "AVMeasurement",
    "AVNetworkRestraintWrapper",
    "DirectLabelingRestraint",
    "LabelingSite",
    "SimpleAVNetworkRestraint",
]
