"""Building a spectrum: from rate constants, and from states plus interaction terms.

The reductions that turn a forward model's internals into the contract stated in
:mod:`IMP.bff.observables`.
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

from IMP.bff.observables.spectrum import LifetimeSpectrum

__all__ = ["lifetime_spectrum_from_rates", "lifetime_spectrum_from_states",
           "rate_constants"]


def lifetime_spectrum_from_rates(
    rates: np.ndarray,
    weights: Optional[np.ndarray] = None,
    exact: bool = True,
) -> LifetimeSpectrum:
    """One species per rate constant.

    The static limit taken literally: if a population of weight ``w_i`` decays
    at ``k_i`` and keeps that rate, then ``F(t) = sum w_i exp(-k_i t)`` is the
    decay, with no fitting and no approximation. Whether that premise holds is
    the caller's to know, and *exact* records the answer.

    :param rates: total deactivation rate per species, 1/ns.
    :param weights: population per species; uniform if omitted.
    :param exact: see :class:`LifetimeSpectrum`.
    """
    k = np.asarray(rates, dtype=np.float64).ravel()
    w = (np.ones_like(k) if weights is None
         else np.asarray(weights, dtype=np.float64).ravel())
    if w.shape != k.shape:
        raise ValueError(f"one weight per rate: {w.shape} against {k.shape}")
    return LifetimeSpectrum(w, k, exact=exact)


def rate_constants(terms: Sequence, *participants, **kwargs) -> np.ndarray:
    """Total deactivation rate per state, in 1/ns.

    A thin, named pass-through to
    :func:`IMP.bff.photophysics.terms.total_rate`, kept here because *the rate
    constants are themselves an observable* -- one of the four this package
    emits. A caller that wants FRET rates rather than a decay should be able to
    ask for them without going through a spectrum.

    The channels add because they are parallel.
    """
    from IMP.bff.photophysics.terms import total_rate
    return total_rate(terms, *participants, **kwargs)


def lifetime_spectrum_from_states(
    terms: Sequence,
    *participants,
    weights: Optional[np.ndarray] = None,
    exact: bool = True,
    **kwargs,
) -> LifetimeSpectrum:
    """The spectrum of a state ensemble under a set of interaction terms.

    One species per state -- per accessible-volume point, per rotamer, per
    conformer -- carrying the summed rate of every channel acting on it. This is
    the reduction that connects :mod:`IMP.bff.representation` and
    :mod:`IMP.bff.photophysics` to an experiment-neutral answer, and it is
    representation-agnostic for the same reason the terms are: it consumes
    states.

    .. warning::
       This is the **static** limit. It is exact when each state holds its rate
       for the whole excited-state lifetime, and wrong when the dye reorganises
       fast enough to average over rates -- then the decay is not a sum of
       exponentials at all and the population has to be propagated instead
       (:class:`~IMP.bff.sampling.smoluchowski.GridDiffusionSolver`, or the Brownian
       walk). Pass ``exact=False`` when using it outside that limit, so the
       spectrum says what it is.

    :param terms: the interaction terms to sum.
    :param participants: states, in the terms' arity order -- the donor's
        states first, then an acceptor's or the quencher atoms.
    :param weights: population per state. Taken from the first participant's
        ``weights`` when omitted, which is what an accessible volume's
        occupancy already is.
    """
    k = np.asarray(rate_constants(terms, *participants, **kwargs),
                   dtype=np.float64).ravel()
    if weights is None and participants:
        w = getattr(participants[0], "weights", None)
        weights = None if w is None else np.asarray(w, dtype=np.float64).ravel()
    return lifetime_spectrum_from_rates(k, weights, exact=exact)
