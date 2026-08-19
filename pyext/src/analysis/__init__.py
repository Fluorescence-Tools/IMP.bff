"""Stage 4 — what the sampled ensemble looks like, as a structural quantity.

The dye's occupancy in space, and quantities read off it. Distinct from
:mod:`IMP.bff.observables`, deliberately: that package owns the projection onto
what an *experiment* would see, with a contract that keeps convolution, pileup,
binning and counting noise outside. This one answers structural questions — how
the dye density is distributed, how a kinetic ensemble transfers — without
committing to a measurement.

Both are stage 4. Keeping them apart means the experiment-neutral contract has a
package to live in, rather than being a rule about part of a larger one.

* :mod:`~IMP.bff.analysis.density` — occupancy grids for a mobile component,
  region-resolved, and their projections.
* :mod:`~IMP.bff.analysis.kinetic_fret` — FRET over a kinetic rotamer ensemble,
  where the transfer is resolved per conformer rather than averaged first.
"""

from __future__ import annotations
