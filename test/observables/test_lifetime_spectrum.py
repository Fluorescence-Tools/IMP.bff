"""The output contract: (amplitude, rate) pairs, what they preserve, and where they fail.

Three things are worth pinning about a lifetime spectrum, and they are
different kinds of claim:

* **Arithmetic.** The two averages are different quantities and reporting one
  under the other's name is a common way to be wrong by tens of percent.
* **Approximation.** Coarse-graining 30 000 states to 128 species preserves
  ``sum(a)`` and ``sum(a*k)`` exactly -- population and initial slope -- and
  nothing higher. That bound is asserted, not assumed.
* **Domain.** A sum of exponentials is *exact* in the static limit and *wrong*
  when the dye averages over rates during the excited-state lifetime. The last
  test measures the gap against a photon simulation of the same model, so the
  ``exact`` flag is a checked claim rather than a docstring.
"""

import numpy as np
import pytest

import IMP.bff
from IMP.bff import (
    LifetimeSpectrum,
    fret_efficiency_from_lifetimes,
    lifetime_spectrum_from_rates,
)


# --- the type ----------------------------------------------------------------

def test_a_single_exponential_is_exactly_one():
    s = LifetimeSpectrum([1.0], [0.25])
    assert s.species_averaged_lifetime == pytest.approx(4.0)
    assert s.intensity_averaged_lifetime == pytest.approx(4.0)
    assert s.decay(np.array([4.0]))[0] == pytest.approx(np.exp(-1.0))


def test_the_two_averages_differ_and_by_how_much():
    """70 % at 4 ns, 30 % at 1 ns: 3.10 by molecule, 3.71 by photon."""
    s = LifetimeSpectrum([0.7, 0.3], [1 / 4.0, 1 / 1.0])
    assert s.species_averaged_lifetime == pytest.approx(3.1)
    assert s.intensity_averaged_lifetime == pytest.approx(11.5 / 3.1)
    assert s.intensity_averaged_lifetime > s.species_averaged_lifetime


def test_the_intensity_average_is_never_the_smaller_one():
    """<tau>_f >= <tau>_x always, with equality only for one species."""
    rng = np.random.default_rng(1)
    for _ in range(30):
        n = int(rng.integers(2, 10))
        s = LifetimeSpectrum(rng.random(n) + 0.1, rng.uniform(0.05, 5.0, n))
        assert s.intensity_averaged_lifetime >= s.species_averaged_lifetime - 1e-12


def test_rates_not_lifetimes_is_the_stored_form():
    """A zero rate is a legitimate non-decaying population."""
    s = LifetimeSpectrum([1.0, 1.0], [0.0, 0.5])
    assert np.isinf(s.lifetimes[0])
    assert s.decay(np.array([1e6]))[0] == pytest.approx(1.0), "it must not decay"
    assert np.isinf(s.species_averaged_lifetime)


def test_a_negative_rate_is_refused():
    with pytest.raises(ValueError, match="growing population"):
        LifetimeSpectrum([1.0], [-0.1])


def test_mismatched_lengths_are_refused():
    with pytest.raises(ValueError, match="one rate per amplitude"):
        LifetimeSpectrum([1.0, 1.0], [0.5])


def test_amplitudes_are_not_silently_normalised():
    """Normalising would discard a quantum yield."""
    s = LifetimeSpectrum([0.4, 0.2], [1.0, 2.0])
    assert s.total_amplitude == pytest.approx(0.6)
    assert s.decay(np.array([0.0]))[0] == pytest.approx(0.6)
    assert s.normalized().total_amplitude == pytest.approx(1.0)


def test_decay_is_the_flat_kernel_on_the_caller_axis():
    s = LifetimeSpectrum([1.0], [0.25])
    t = np.linspace(0.0, 10.0, 12).reshape(3, 4)
    # decay() is the 1-D kernel: it evaluates F(t) on a flat time axis and
    # returns a flat view. The caller's axis shape is not restated by the
    # kernel -- ravel the axis and reshape on your side if you kept a 2-D grid.
    assert s.decay(t.ravel()).shape == (12,)


# --- coarse-graining ---------------------------------------------------------

def test_coarse_graining_preserves_population_and_initial_slope_exactly():
    rng = np.random.default_rng(0)
    k = rng.uniform(0.2, 3.0, 30000)
    w = rng.random(30000)
    big = lifetime_spectrum_from_rates(k, w)
    small = big.coarse_grain(128)

    assert small.n_species <= 128 < big.n_species
    assert small.total_amplitude == pytest.approx(big.total_amplitude, rel=1e-12)
    assert (np.dot(small.amplitudes, small.rate_constants)
            == pytest.approx(np.dot(big.amplitudes, big.rate_constants), rel=1e-12))


def test_coarse_graining_error_falls_with_the_square_of_the_bin_width():
    """The claim in the docstring, measured rather than asserted in prose."""
    rng = np.random.default_rng(4)
    big = lifetime_spectrum_from_rates(rng.uniform(0.2, 3.0, 20000),
                                       rng.random(20000))
    t = np.linspace(0.0, 15.0, 200)
    reference = big.decay(t) / big.total_amplitude

    errors = []
    for n_bins in (8, 16, 32):
        cg = big.coarse_grain(n_bins)
        errors.append(np.abs(cg.decay(t) / cg.total_amplitude - reference).max())
    # halving the bin width should quarter the error; allow a factor of 2 slack
    assert errors[0] / errors[1] > 2.0
    assert errors[1] / errors[2] > 2.0
    assert errors[2] < 1e-3


def test_coarse_graining_marks_the_result_approximate_only_if_it_merged():
    exact = LifetimeSpectrum([1.0, 1.0], [0.5, 2.0])
    assert exact.coarse_grain(1024).exact, "nothing was merged, so nothing was lost"
    assert not exact.coarse_grain(1).exact


def test_a_degenerate_spectrum_coarse_grains_to_one_species():
    """Every species at the same rate: nothing to bin, and no zero-width divide."""
    s = lifetime_spectrum_from_rates(np.full(500, 0.25), np.full(500, 2.0))
    cg = s.coarse_grain(64)
    assert cg.n_species == 1
    assert cg.rate_constants[0] == pytest.approx(0.25)
    assert cg.total_amplitude == pytest.approx(1000.0)


def test_empty_bins_are_dropped_not_returned_as_zero_amplitude_species():
    """A bimodal spectrum leaves most of the rate range empty."""
    k = np.concatenate([np.full(100, 0.1), np.full(100, 5.0)])
    cg = lifetime_spectrum_from_rates(k).coarse_grain(64)
    assert cg.n_species == 2
    assert set(np.round(cg.rate_constants, 6)) == {0.1, 5.0}


# --- efficiency --------------------------------------------------------------

def test_efficiency_uses_the_species_average():
    """Efficiency is per molecule, so the photon-weighted average is wrong here.

    The two only diverge when transfer is *heterogeneous* -- an acceptor that
    halves every lifetime uniformly gives the same number either way, which is
    why the case below quenches only the long-lived species. That is also the
    realistic case: a dye near an acceptor is exactly the one whose lifetime
    moves.
    """
    donor = LifetimeSpectrum([0.5, 0.5], [1 / 4.0, 1 / 1.0])
    da = LifetimeSpectrum([0.5, 0.5], [1 / 1.0, 1 / 1.0])   # only the 4 ns one transfers

    assert donor.species_averaged_lifetime == pytest.approx(2.5)
    assert da.species_averaged_lifetime == pytest.approx(1.0)
    assert fret_efficiency_from_lifetimes(donor, da) == pytest.approx(0.6)

    # the intensity average weights each molecule by the photons it emitted --
    # the very thing transfer suppresses -- and overstates E by 10 points here
    wrong = 1.0 - da.intensity_averaged_lifetime / donor.intensity_averaged_lifetime
    assert wrong == pytest.approx(1.0 - 1.0 / 3.4)
    assert abs(wrong - 0.6) > 0.1


def test_efficiency_limits():
    donor = LifetimeSpectrum([1.0], [0.25])
    assert fret_efficiency_from_lifetimes(donor, donor) == pytest.approx(0.0)
    assert fret_efficiency_from_lifetimes(
        donor, LifetimeSpectrum([1.0], [1e9])) == pytest.approx(1.0, abs=1e-6)


def test_efficiency_of_a_dead_donor_is_zero_not_a_divide():
    dead = LifetimeSpectrum([1.0], [np.inf]) if False else LifetimeSpectrum([0.0], [1.0])
    assert fret_efficiency_from_lifetimes(dead, dead) == 0.0


# --- the domain of validity --------------------------------------------------

def test_the_static_limit_is_exact_when_the_rate_does_not_move():
    """A constant quenching rate: the spectrum and the photon race must agree.

    This is the boundary condition of the whole contract -- if the rate never
    changes there is nothing to average over, so the sum of exponentials is the
    decay and the simulation must reproduce it.
    """
    import IMP.bff as photon

    tau0, kq_val, t_step = 4.0, 0.5, 0.01
    dts, emitted = photon.simulate_photon_trace(
        400000, np.full(4000, kq_val), t_step, tau0, random_seed=1)
    measured = dts[emitted > 0].mean()

    spectrum = LifetimeSpectrum([1.0], [1.0 / tau0 + kq_val])
    assert measured == pytest.approx(spectrum.species_averaged_lifetime, rel=0.02)


def test_a_moving_rate_breaks_the_static_limit_measurably():
    """Two rates the dye alternates between are not two species.

    Fast exchange averages the *rate*; the static limit averages the *decay*.
    They differ, which is why the spectrum built from a trajectory is marked
    inexact -- and why a diffusion simulation exists at all.
    """
    import IMP.bff as photon

    tau0, t_step = 4.0, 0.01
    slow, fast = 0.05, 3.0
    # alternate every frame: the dye sees the mean rate, not either one
    kq = np.tile([slow, fast], 2000)
    dts, emitted = photon.simulate_photon_trace(400000, kq, t_step, tau0, random_seed=2)
    measured = dts[emitted > 0].mean()

    static = LifetimeSpectrum([0.5, 0.5],
                              [1 / tau0 + slow, 1 / tau0 + fast])
    averaged = LifetimeSpectrum([1.0], [1 / tau0 + 0.5 * (slow + fast)])

    assert measured == pytest.approx(averaged.species_averaged_lifetime, rel=0.05)
    assert abs(measured - static.species_averaged_lifetime) > 0.2, (
        "the static limit must be visibly wrong here, or this test proves nothing")


def test_no_convolution_anywhere_in_the_observables_package():
    """The contract, checked as code rather than trusted as prose.

    The check used to glob ``pyext/src/observables/*.py`` -- a directory that
    stopped existing when the package became one module, so it scanned nothing
    and passed on an empty loop. The contract now lives in
    ``pyext/IMP_bff.observables.i`` and its header, so those are what it reads.

    ``imp_bff.observables.i`` no longer carries ``%pythoncode`` (it was ported
    to C++), so there is no Python surface to scan; the C++ check below is the
    whole of it.

    Prose is exempt on purpose, and has to be: the header *states* the contract
    by naming what it excludes ("folding an IRF into it belongs to whatever owns
    the instrument"), and `decay()` documents itself as **unconvolved**. What
    must not appear is a *call*.
    """
    import ast
    import re
    from pathlib import Path

    root = Path(__file__).resolve().parents[2]
    banned = {"convolve", "fftconvolve", "irf", "pileup"}

    # The Python surface is gone -- the file has no %pythoncode. Guard that it
    # stays that way: a ported file must not grow a Python wrapper back.
    swig = (root / "pyext" / "IMP_bff.observables.i").read_text()
    assert "%pythoncode" not in swig, "observables.i must stay python-free"

    # The C++ side: comments stripped, then whole identifiers.
    for name in ("include/LifetimeSpectrum.h", "src/LifetimeSpectrum.cpp"):
        path = root / name
        assert path.exists(), path
        code = re.sub(r"//[^\n]*|/\*.*?\*/", " ", path.read_text(), flags=re.S).lower()
        for word in banned:
            assert not re.search(rf"\b{word}\b", code), (name, word)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
