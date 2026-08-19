"""What this package emits: experiment-neutral quantities, and nothing else.

`IMP.bff` is a **forward model**. It takes a structure, a solvent and labels and
answers questions about the *sample*. It does not answer questions about an
*instrument*, and the boundary between the two is this module's subject.

The contract
------------

Every observable here is one of:

* a :class:`LifetimeSpectrum` -- ``(amplitude, rate constant)`` pairs, rates in
  1/ns;
* **rate constants** -- FRET, PET, radiative, in 1/ns, per state or per frame;
* a **distance** or a distribution of them, in Angstrom;
* a **kappa^2** distribution.

Deliberately absent, and to stay absent:

* convolution with an instrument response -- an IRF belongs to a measurement,
  and folding one in here makes the model's output un-reusable by any other
  measurement;
* pileup, dead time, afterpulsing, scatter, background;
* counting noise -- a spectrum is the expectation, and Poisson sampling it is
  the experiment's job;
* binning -- a TAC's bin width and time range are instrument settings. That is
  why the primary output is a spectrum and not a histogram: a curve on a grid
  has already chosen both.

Unconvolved curves *are* allowed, through
:meth:`LifetimeSpectrum.decay` -- evaluating the model on a time axis the caller
chose is not the same as folding in an instrument. It is a convenience over the
spectrum, never a substitute for it.

When a spectrum is exact, and when it is an approximation
---------------------------------------------------------

A sum of exponentials is **exact** whenever each species keeps its rate constant
for the whole excited-state lifetime -- one species per accessible-volume point,
per rotamer, per conformer. This is the static (quasi-static) limit, and it is
where most of the structure-based modelling in this package lives.

It is **wrong** when the dye reorganises fast enough to average over rate
constants during the lifetime. Then the decay is not multi-exponential at all
and no set of amplitudes reproduces it; the population has to be propagated,
which is what :class:`~IMP.bff.sampling.smoluchowski.GridDiffusionSolver` and the
Brownian walk are for. Their output is a curve, and
:func:`lifetime_spectrum_from_states` is not the right reduction for it.

Saying which of the two a number came from is part of the answer, so the
builders below record it rather than leaving the caller to guess.
"""

from IMP.bff.observables.spectrum import (
    LifetimeSpectrum,
    fret_efficiency_from_lifetimes,
)
from IMP.bff.observables.reduce import (
    lifetime_spectrum_from_rates,
    lifetime_spectrum_from_states,
    rate_constants,
)

__all__ = [
    "LifetimeSpectrum",
    "fret_efficiency_from_lifetimes",
    "lifetime_spectrum_from_rates",
    "lifetime_spectrum_from_states",
    "rate_constants",
]
