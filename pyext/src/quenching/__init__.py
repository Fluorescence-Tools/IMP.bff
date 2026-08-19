"""PET quenching of a dye diffusing in its accessible volume.

The model QuEst was built on, moved here by PRD-109 so that ``IMP.bff`` owns it:
a dye sphere diffuses on an accessible-volume grid, is slowed near sticky
residues, and is quenched by photo-induced electron transfer whenever it comes
within contact distance of a redox-active side chain. The observable is the
donor's fluorescence decay, and comparing it with a measured lifetime is what
calibrates the accessible *contact* volume.

Three layers, each usable on its own:

``pet``, ``asa``
    The chemistry: which residues quench, how fast, where on the residue, and
    how buried it is.
``grids``, ``maps``, ``fret_trace``
    The fields: stamp the quenching rate, the mobility and the FRET rate onto
    the accessible-volume grid, and read them along a trajectory.
``model``
    :class:`DyeDiffusionSimulation` and :class:`QuenchedDonorDecay`, which put
    the layers together for one labelling site.

The **integrators moved out** in PRD-113 stage 7. The Brownian walk, the
Smoluchowski field solver and the excited-state Monte Carlo are now
:mod:`IMP.bff.sampling`: none of them is specific to PET quenching -- a Brownian
walk in a volume is a Brownian walk in a volume -- and filing a general
integrator under the first physics that used it is how it comes to look like a
detail of one model.

Not to be confused with :class:`IMP.bff.LangevinDyeSampler` (PRD-108), which
integrates an **explicit all-atom dye** under a force field. That is a different
model at a different scale; this one is a sphere on a grid, cheap enough to scan
every position in a protein.
"""

from __future__ import annotations
