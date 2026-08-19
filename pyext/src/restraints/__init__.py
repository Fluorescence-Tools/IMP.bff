"""IMP.bff.restraints — scoring restraints for FRET and labeling data.

Three restraints, all scoring a *structure* against derived fluorescence
observables, which is what IMP.bff is for:

* :class:`AVNetworkRestraintWrapper` — the fps.json-driven AV network, wrapped
  for IMP.pmi.
* :class:`SimpleAVNetworkRestraint` — chi-squared scoring of a network of
  AV-pair distance measurements, without an fps.json.
* :class:`DirectLabelingRestraint` — fast scoring from attachment-atom
  positions, suitable for MCMC inner loops.

The last two arrived here from imp-tricks in 2026-08. They had been shipped
there as additions to this package, which is legal under the add-never-replace
rule but put half of one namespace in another repository -- and a subpackage
resolves to the first matching directory and stops, so for five days
``AVNetworkRestraintWrapper`` was unreachable and recorded as a deleted upstream
API. Everything under ``IMP.bff`` now ships from imp.bff.
"""

from __future__ import annotations

from .network import AVNetworkRestraintWrapper
from .simple_av_network import SimpleAVNetworkRestraint, AVMeasurement
from .direct_labeling import DirectLabelingRestraint, LabelingSite

__all__ = [
    "AVNetworkRestraintWrapper",
    "SimpleAVNetworkRestraint",
    "AVMeasurement",
    "DirectLabelingRestraint",
    "LabelingSite",
]
