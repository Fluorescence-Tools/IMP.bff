"""Stage 2 — is this dye configuration allowed, and how heavily does it count?

The score of a *configuration*, from the structure alone. A rotamer that clashes
with the backbone is not a rotamer the dye adopts; a conformer at high internal
energy is one it adopts rarely. This is what turns a raw set of candidate
positions into a weighted ensemble, and it is the same question whichever
representation produced the candidates:

======================  ==============================================
representation          what scoring it means
======================  ==============================================
accessible volume       occupancy: a voxel is reachable or it is not
rotamer library         clash energy of each conformer against the site
coarse-grained dye      the simplified force field's internal energy
======================  ==============================================

**Not the same thing as** :mod:`IMP.bff.restraints`. That answers *does this
model agree with a measurement* -- it takes experimental distances and returns
a chi-squared. This package never sees an experiment: it takes coordinates and
returns an energy or a weight. Two different questions that both get called
"scoring", which is why they are two packages and why this paragraph exists.

* :mod:`~IMP.bff.scoring.lennard_jones` -- sterics, and the CHARMM36 table.
  Shared: the rotamer scorer and the explicit dye's force field read the same
  parameters, so they cannot drift apart.
* :mod:`~IMP.bff.scoring.rotamer` -- a conformer against its site.
* :mod:`~IMP.bff.scoring.dye_lj` -- LJ pair scoring and exclusion lists for an
  explicit dye.
* :mod:`~IMP.bff.scoring.dye_restraints` -- the same, as IMP restraints.
* :mod:`~IMP.bff.scoring.boltzmann` -- weights from energies.
* :mod:`~IMP.bff.scoring.mean_field` -- weights solved self-consistently, for
  when two dyes see each other.
"""

from __future__ import annotations
