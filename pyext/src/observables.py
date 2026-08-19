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

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Sequence
import csv

import numpy as np

import IMP.bff.representation.distance as _dist
import IMP.bff.representation.av as _av
import IMP.bff

__all__ = [
    'LifetimeSpectrum',
    'compute_distance_distributions',
    'fret_efficiency_from_lifetimes',
    'lifetime_spectrum_from_rates',
    'lifetime_spectrum_from_states',
    'rate_constants',
]

# --------------------------------------------------------------------------
# pair_distribution
# --------------------------------------------------------------------------
"""Full FRET distance distributions P(R_DA) for a structure (FPS-style).

For each experimental pair, compute the accessible volumes of the two dyes on a
(docked) structure and the distribution of donor-acceptor distances between the
two AV point clouds — the P(R_DA) that FPS reports, richer than the single mean
distance used for fast docking.
"""

def compute_distance_distributions(
    pdb_path: str,
    positions: Dict,
    distances: Dict,
    out_csv: Optional[str] = None,
    *,
    rda_min: float = 1.0,
    rda_max: float = 200.0,
    n_bins: int = 100,
) -> Dict:
    """Return per-pair ``P(R_DA)`` distributions for the structure at ``pdb_path``.

    Parameters
    ----------
    pdb_path : str
        Structure (e.g. a docked PDB) on which to compute the AVs.
    positions, distances : dict
        The ``Positions`` / ``Distances`` sections of the fps.json.
    out_csv : str, optional
        Write a table with an ``R_DA`` column plus one probability column per
        pair.

    Returns
    -------
    dict
        ``{"rda_axis": [...], "pairs": {name: {"p": [...], "mean": float}},
        "distributions_csv": path}``.
    """
    atoms = _av.load_structure_with_vdw(pdb_path)
    avs = _av.compute_avs_for_structure(atoms, positions, pdb_path=pdb_path)

    rda_axis = None
    pairs: Dict[str, Dict] = {}
    for name, d in distances.items():
        av1 = avs.get(d.get("position1_name"))
        av2 = avs.get(d.get("position2_name"))
        if av1 is None or av2 is None or not av1.has_volume or not av2.has_volume:
            continue
        # histogram_rda returns (histogram, bin_edges)
        p, edges = _dist.histogram_rda(
            av1, av2, rda_min=rda_min, rda_max=rda_max, n_rda_bins=n_bins,
            normalize=True)
        centers = 0.5 * (np.asarray(edges[:-1]) + np.asarray(edges[1:]))
        rda_axis = centers
        p = np.asarray(p, dtype=float)
        total = float(p.sum())
        mean = float(np.sum(centers * p) / total) if total > 0 else float("nan")
        pairs[name] = {"p": p.tolist(), "mean": mean}

    if out_csv and pairs and rda_axis is not None:
        names = list(pairs)
        with open(out_csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["R_DA"] + names)
            for i, r in enumerate(rda_axis):
                w.writerow([round(float(r), 2)]
                           + [round(pairs[n]["p"][i], 6) for n in names])

    return {
        "rda_axis": rda_axis.tolist() if rda_axis is not None else [],
        "pairs": pairs,
        "distributions_csv": out_csv if (out_csv and pairs) else None,
    }


# --------------------------------------------------------------------------
# spectrum
# --------------------------------------------------------------------------
"""The lifetime spectrum: ``(amplitude, rate constant)`` pairs, and its reductions.

The output contract is stated in :mod:`IMP.bff.observables`. This module is the
type that carries it.
"""

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


# --------------------------------------------------------------------------
# reduce
# --------------------------------------------------------------------------
"""Building a spectrum: from rate constants, and from states plus interaction terms.

The reductions that turn a forward model's internals into the contract stated in
:mod:`IMP.bff.observables`.
"""

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
    from IMP.bff.photophysics import total_rate
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
