"""The lifetime spectrum: ``(amplitude, rate constant)`` pairs, and its reductions.

The output contract is stated in :mod:`IMP.bff.observables`. This module is the
type that carries it.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

import IMP.bff

__all__ = ["LifetimeSpectrum", "fret_efficiency_from_lifetimes"]


@dataclass(frozen=True)
class LifetimeSpectrum:
    """A set of deactivation channels and the population that carries each.

    :param amplitudes: population fraction per species. Not required to be
        normalised -- an unnormalised spectrum still carries the *relative*
        answer, and normalising silently would discard a quantum yield.
    :param rate_constants: total deactivation rate per species, **1/ns**.
        Rates, not lifetimes: rate constants are what the physics adds
        (``k = k_radiative + k_PET + k_FRET + ...``), and a species with zero
        rate is a legitimate non-decaying population while a species with
        infinite lifetime is an awkward one.
    :param exact: ``True`` when each species genuinely holds its rate for the
        whole excited-state lifetime, so the sum of exponentials is the decay
        rather than an approximation to it. Recorded rather than assumed --
        see :mod:`IMP.bff.observables`.
    """

    amplitudes: np.ndarray
    rate_constants: np.ndarray
    exact: bool = True

    def __post_init__(self):
        a = np.ascontiguousarray(np.asarray(self.amplitudes, dtype=np.float64).ravel())
        k = np.ascontiguousarray(np.asarray(self.rate_constants, dtype=np.float64).ravel())
        if a.shape != k.shape:
            raise ValueError(
                f"a spectrum needs one rate per amplitude: {a.shape} against {k.shape}")
        if np.any(k < 0.0):
            raise ValueError("a negative rate constant is a growing population, "
                             "not a decaying one")
        object.__setattr__(self, "amplitudes", a)
        object.__setattr__(self, "rate_constants", k)

    # -- basic properties ----------------------------------------------------

    @property
    def n_species(self) -> int:
        return int(self.amplitudes.size)

    @property
    def lifetimes(self) -> np.ndarray:
        """``1 / k`` in ns; ``inf`` for a non-decaying species."""
        with np.errstate(divide="ignore"):
            return np.where(self.rate_constants > 0.0,
                            1.0 / np.where(self.rate_constants > 0.0,
                                           self.rate_constants, 1.0),
                            np.inf)

    @property
    def total_amplitude(self) -> float:
        return float(self.amplitudes.sum())

    def normalized(self) -> "LifetimeSpectrum":
        """The same spectrum with amplitudes summing to 1."""
        total = self.total_amplitude
        if total == 0.0:
            return self
        return LifetimeSpectrum(self.amplitudes / total, self.rate_constants, self.exact)

    # -- the two averages ----------------------------------------------------

    @property
    def species_averaged_lifetime(self) -> float:
        r"""\f$\langle\tau\rangle_x = \sum a_i \tau_i / \sum a_i\f$.

        The average over *molecules*. This is the one that enters a FRET
        efficiency, because efficiency is defined per molecule.
        """
        total = self.total_amplitude
        if total == 0.0:
            return 0.0
        finite = np.isfinite(self.lifetimes)
        if not finite.all():
            return np.inf
        return float(np.dot(self.amplitudes, self.lifetimes) / total)

    @property
    def intensity_averaged_lifetime(self) -> float:
        r"""\f$\langle\tau\rangle_f = \sum a_i \tau_i^2 / \sum a_i \tau_i\f$.

        The average over *photons*, which is what a long-lived species
        dominates. Equal to the species average only for a single exponential;
        the gap between the two is a measure of how heterogeneous the sample is,
        and reporting one under the other's name is a common way to be wrong by
        tens of percent.
        """
        tau = self.lifetimes
        if not np.isfinite(tau).all():
            return np.inf
        denom = float(np.dot(self.amplitudes, tau))
        if denom == 0.0:
            return 0.0
        return float(np.dot(self.amplitudes, tau ** 2) / denom)

    # -- derived curves and reductions ---------------------------------------

    def decay(self, time: np.ndarray) -> np.ndarray:
        """``F(t) = sum_i a_i exp(-k_i t)`` on the caller's axis. **Unconvolved.**

        No instrument response, no pileup, no counting noise. Folding an IRF in
        belongs to whatever owns the instrument -- see :mod:`IMP.bff.observables`.
        """
        t = np.asarray(time, dtype=np.float64)
        out = IMP.bff.lifetime_spectrum_decay(
            self.amplitudes, self.rate_constants,
            np.ascontiguousarray(t.ravel()))
        return np.asarray(out, dtype=np.float64).reshape(t.shape)

    def coarse_grain(self, n_bins: int = 128) -> "LifetimeSpectrum":
        """Reduce many species to few, preserving the first two moments.

        An accessible volume gives one species per point -- tens of thousands,
        almost none of them distinguishable. Species are binned on a uniform
        grid in the *rate constant* and each bin keeps its summed amplitude and
        its amplitude-weighted mean rate, so ``sum(a)`` and ``sum(a*k)`` are
        exact for any *n_bins*: the population and the initial slope ``F'(0)``
        survive. Curvature does not -- it is under-stated by an amount falling
        as the square of the bin width, which is why this returns a spectrum
        marked inexact unless nothing was actually merged.
        """
        rates = IMP.bff.VectorDouble()
        amps = IMP.bff.lifetime_spectrum_coarse_grain(
            self.amplitudes, self.rate_constants, int(n_bins), rates)
        a = np.asarray(amps, dtype=np.float64)
        return LifetimeSpectrum(a, np.asarray(rates, dtype=np.float64),
                                exact=self.exact and a.size == self.n_species)

    def __len__(self) -> int:
        return self.n_species

    def __repr__(self) -> str:
        tau = self.species_averaged_lifetime
        kind = "exact" if self.exact else "approximate"
        return (f"LifetimeSpectrum({self.n_species} species, "
                f"<tau>_x = {tau:.4g} ns, {kind})")


def fret_efficiency_from_lifetimes(
    donor_only: "LifetimeSpectrum",
    donor_acceptor: "LifetimeSpectrum",
) -> float:
    r"""\f$E = 1 - \langle\tau\rangle_{x,DA} / \langle\tau\rangle_{x,D}\f$.

    The **species** average, not the intensity average. Efficiency is defined
    per molecule -- the fraction of excitations that transfer -- and the
    intensity average weights each molecule by how many photons it emitted,
    which is exactly the thing transfer suppresses. Using it gives an efficiency
    that is systematically too low, by tens of percent on a heterogeneous
    sample.

    :param donor_only: the donor without an acceptor.
    :param donor_acceptor: the same donor with one.
    :returns: efficiency in [0, 1]; 0 if the donor-only lifetime is zero.
    """
    tau_d = donor_only.species_averaged_lifetime
    if tau_d <= 0.0 or not np.isfinite(tau_d):
        return 0.0
    return float(1.0 - donor_acceptor.species_averaged_lifetime / tau_d)
