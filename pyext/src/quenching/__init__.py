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
``grids``, ``diffusion``, ``photon``, ``fret_trace``
    The kernels: stamp the rate and stickiness maps, run the Brownian walk,
    turn a per-frame rate into a decay curve or a photon histogram.
``model``
    :class:`DyeDiffusionSimulation` and :class:`QuenchedDonorDecay`, which put
    the three together for one labelling site.

Not to be confused with :class:`IMP.bff.LangevinDyeSampler` (PRD-108), which
integrates an **explicit all-atom dye** under a force field. That is a different
model at a different scale; this one is a sphere on a grid, cheap enough to scan
every position in a protein.
"""

from __future__ import annotations
