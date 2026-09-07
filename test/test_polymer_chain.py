"""Polymer chain distributions, and the linker convolution that used to raise.

``worm_like_chain_linker`` had **no test**, which is why it went on raising
``TypingError`` for a whole commit: it was numba-jitted and called
``normal_distribution``, which had become a C++ delegation numba cannot type.
The C++ port fixed it, and this file is what stops it recurring silently.
"""

import math

import numpy as np
import pytest

import IMP.bff
from IMP.bff import (
    distance_between_gaussian, gaussian_chain, gaussian_chain_ree, i0,
    i0_array, ising_chain, saw_nu, worm_like_chain, worm_like_chain_linker,
)

R = np.linspace(0.5, 120.0, 241)


def test_ideal_chain_scales_as_sqrt_n():
    assert gaussian_chain_ree(3.8, 20) == pytest.approx(3.8 * np.sqrt(20))
    assert gaussian_chain_ree(3.8, 80) == pytest.approx(2.0 * gaussian_chain_ree(3.8, 20))


def test_ideal_chain_peaks_near_the_rms_extension():
    ree = gaussian_chain_ree(3.8, 20)
    pdf = gaussian_chain(R, 3.8, 20)
    assert pdf.min() >= 0.0
    peak = R[int(np.argmax(pdf))]
    assert 0.5 * ree < peak < 1.5 * ree


def test_a_worm_like_chain_cannot_exceed_its_contour():
    pdf = worm_like_chain(R, 0.4, chain_length=60.0)
    assert np.all(pdf[R >= 60.0] == 0.0)
    assert pdf[R < 60.0].sum() == pytest.approx(1.0)


@pytest.mark.parametrize("kappa", [0.05, 0.125, 0.4, 2.0])
def test_it_is_a_distribution_on_both_sides_of_the_branch(kappa):
    """The closed form branches at kappa = 0.125; both branches must behave."""
    pdf = worm_like_chain(R, kappa)
    assert np.all(np.isfinite(pdf)) and pdf.min() >= 0.0
    assert pdf.sum() == pytest.approx(1.0)


def test_the_linker_convolution_runs_and_broadens():
    """The regression: this raised TypingError and nothing noticed."""
    bare = worm_like_chain(R, 0.4)
    broad = worm_like_chain_linker(R, 0.4, sigma=6.0)
    assert np.all(np.isfinite(broad))
    assert broad.sum() == pytest.approx(1.0)
    sd = lambda p: np.sqrt(np.sum(p * (R - np.sum(p * R)) ** 2))
    assert sd(broad) > sd(bare), "linker broadening must widen the distribution"


def test_the_linker_convolution_matches_an_independent_one():
    """Checked against numpy rather than against the version it replaced.

    The numba original could not run at all by then, so there was nothing to
    compare with -- an equality gate needs a reference that works.

    **Corrected 2026-09-02, and worth reading as a cautionary note.** This test
    passed for the whole time `worm_like_chain_linker` was wrong, because its
    "independent" reference made the same two mistakes as the code: it
    convolved with a plain Gaussian (`g` below used to be the normal density)
    where the kernel is the *distance between two Gaussian clouds*, and it fed
    in a chain distribution carrying the r^2 factor that the convolution must
    not see. An independent reference is only independent if it is derived from
    the definition, not written next to the implementation.
    """
    def two_gaussian_distance(x, separation, sigma):
        """P(distance) between two points each Gaussian about their centres."""
        g = (lambda v, loc: np.exp(-((v - loc) ** 2) / (2 * sigma * sigma))
             / (np.sqrt(2 * np.pi) * sigma))
        if separation > 0.0:
            return x / separation * (g(x, separation) - g(x, -separation))
        return 2.0 * x ** 2 / sigma ** 2 * g(x, 0.0)

    pr = worm_like_chain(R, 0.4, 0.0, True, False)     # bare, not r^2-weighted
    ref = sum(pr[i] * two_gaussian_distance(R, R[i], 6.0)
              for i in range(len(R)) if pr[i] != 0.0)
    ref = ref / ref.sum()
    np.testing.assert_allclose(worm_like_chain_linker(R, 0.4, sigma=6.0), ref,
                               rtol=1e-12, atol=0)


# --- the Ising two-state chain -------------------------------------------
#
# Ported from `chisurf.core.math.functions.rdf.ising_chain` because it was the
# single largest piece of model compute ChiSurf still ran in the interpreter:
# 55 ms a curve, 0.6 s per LM iteration at its ten free parameters, which is
# more than half of everything `test/minimizer/bench_models.py` measures as
# movable. The transfer-matrix product is the whole cost and it is a Python
# loop numpy cannot vectorise -- ~80k iterations of 2x2 arithmetic.

ISING_R = np.linspace(1.0, 200.0, 96)


def test_the_ising_chain_is_a_normalised_distribution():
    pr = ising_chain(ISING_R, 40, 4.0, 8.0)
    assert np.all(np.isfinite(pr)) and pr.min() >= 0.0
    assert np.trapezoid(pr, ISING_R) == pytest.approx(1.0)


def test_one_bond_length_reduces_to_the_ideal_chain():
    """b_S = b_U removes the two states, so the Ising chain *is* Gaussian.

    The reduction is the cheapest check that the transfer matrix is being
    contracted correctly: with one bond variance the Ising weights cancel out
    of the ratio phi(k)/phi(0) entirely, whatever J and h are.
    """
    pr = ising_chain(ISING_R, 40, 6.0, 6.0, coupling=1.5, field=0.3)
    ideal = gaussian_chain(ISING_R, 6.0, 40)
    ideal = ideal / np.trapezoid(ideal, ISING_R)
    assert np.max(np.abs(pr - ideal)) < 2e-3


def test_the_ising_weights_drop_out_when_the_bonds_are_equal():
    """Same reduction, stated as an invariance: J and h cannot matter."""
    a = ising_chain(ISING_R, 30, 5.0, 5.0, coupling=0.2, field=-2.0)
    b = ising_chain(ISING_R, 30, 5.0, 5.0, coupling=4.0, field=3.0)
    assert np.max(np.abs(a - b)) < 1e-12


@pytest.mark.parametrize("field,expected", [(-6.0, 9.0), (6.0, 3.0)])
def test_the_field_drives_the_chain_to_one_state(field, expected):
    """h biases towards structured, so a large |h| picks out one bond length.

    Written as the *mean* extension rather than as a curve comparison: the
    two limits are the two ideal chains, and which one it lands on is the
    sign convention this kernel has to get right.
    """
    pr = ising_chain(ISING_R, 40, 3.0, 9.0, coupling=0.0, field=field)
    mean = np.trapezoid(pr * ISING_R, ISING_R)
    reference = gaussian_chain(ISING_R, expected, 40)
    reference = reference / np.trapezoid(reference, ISING_R)
    assert mean == pytest.approx(
        np.trapezoid(reference * ISING_R, ISING_R), rel=0.02)


def test_a_degenerate_chain_is_zero_rather_than_a_nan():
    assert np.all(ising_chain(ISING_R, 0, 4.0, 8.0) == 0.0)
    assert np.all(np.isfinite(ising_chain(ISING_R, 1, 4.0, 8.0)))


def test_a_long_strongly_coupled_chain_stays_finite():
    """The transfer product is formed unscaled, so overflow is possible.

    `include/PolymerChain.h` says why it is not rescaled -- a per-k rescaling
    would change the ratio it is used in. What this pins is the consequence:
    an overflowed phi must leave zeros behind, never a NaN that propagates
    into the FRET rate spectrum downstream.
    """
    pr = ising_chain(ISING_R, 400, 4.0, 8.0, coupling=8.0, field=2.0)
    assert np.all(np.isfinite(pr)) and pr.min() >= 0.0


# --- I0, and the worm-like chain that was computing exp instead ------------
#
# `worm_like_chain` multiplied by `exp(x)` where the Becker-Rosa-Everaers form
# has `I0(x)`, from the original transcription until 2026-09-02. Nothing but
# this file called it -- ChiSurf ran its own Python copy -- so the error was
# latent, and it would have gone live the moment that copy was forwarded here,
# which is exactly the change that found it.


def test_i0_is_even_which_is_the_whole_point():
    """I0(-x) = I0(x). exp(-x) does not, and that was the bug."""
    for x in (0.5, 3.0, 3.75, 8.0, 25.0):
        assert i0(-x) == pytest.approx(i0(x), rel=1e-15)


def test_i0_matches_its_series_on_both_sides_of_the_branch():
    """The polynomial branches at |x| = 3.75; check against the I0 series.

    The tolerance is 3e-4, not machine precision, and that is the honest
    number: this is the Abramowitz & Stegun polynomial, kept deliberately
    inexact so that fitted distributions do not move (see `SpecialFunctions.h`).
    What the test is for is a transcription slip, which moves it by far more.
    """
    for x in (0.1, 1.0, 3.7, 3.8, 6.0):
        series = sum((x / 2.0) ** (2 * k) / float(math.factorial(k)) ** 2
                     for k in range(80))
        assert i0(x) == pytest.approx(series, rel=3e-4)


def test_i0_array_agrees_with_the_scalar():
    x = np.linspace(-20.0, 20.0, 101)
    assert np.allclose(i0_array(x), [i0(v) for v in x], rtol=0, atol=0)


def test_the_worm_like_chain_extends_as_it_stiffens():
    """The regression, stated as physics rather than as a number.

    <R^2>/L^2 must RISE with kappa = lp/L -- a stiffer chain is a more extended
    one, and Kratky-Porod gives the exact value. With `exp` in place of `I0`
    this sequence fell instead: 0.070, 0.026, 0.007 against an exact
    0.095, 0.736, 0.852. Any future transcription slip in that factor moves the
    trend, so the trend is what is pinned.
    """
    r = np.linspace(1e-6, 1.0 - 1e-9, 20001)
    moments = []
    for kappa in (0.05, 0.1, 0.25, 0.5, 1.0, 2.0):
        pdf = worm_like_chain(r, kappa, 1.0, True, True)
        moments.append(float(np.sum(pdf * r * r) / np.sum(pdf)))
    assert np.all(np.diff(moments) > 0.0), moments
    # and it must stay in the right neighbourhood of the exact result
    for kappa, got in zip((0.05, 0.1, 0.25, 0.5, 1.0, 2.0), moments):
        exact = 2 * kappa - 2 * kappa ** 2 * (1 - np.exp(-1.0 / kappa))
        assert 0.5 * exact < got < 1.2 * exact


def test_the_contour_cut_is_a_prefix_not_a_mask():
    """The reference `break`s at the first r >= L, so the tail stays zero.

    On an unsorted axis a prefix and a mask are different functions, and
    nothing requires the axis to be sorted. A short distance placed *after* a
    long one must therefore still be zero.
    """
    axis = np.array([10.0, 20.0, 70.0, 30.0, 40.0])
    pdf = worm_like_chain(axis, 0.4, 60.0, False, False)
    assert pdf[0] > 0.0 and pdf[1] > 0.0
    assert pdf[2] == 0.0                      # reaches the contour length
    assert pdf[3] == 0.0 and pdf[4] == 0.0    # after it, though both are short


# --- the self-avoiding walk ------------------------------------------------


def test_saw_nu_reduces_to_the_ideal_chain():
    """nu = 0.5 with gamma_exp = 1 is the Gaussian chain, by construction."""
    axis = np.linspace(0.5, 200.0, 400)
    r_rms = gaussian_chain_ree(3.8, 40)
    saw = saw_nu(axis, r_rms, 0.5, 1.0)
    ideal = gaussian_chain(axis, 3.8, 40)
    saw = saw / np.trapezoid(saw, axis)
    ideal = ideal / np.trapezoid(ideal, axis)
    assert np.max(np.abs(saw - ideal)) < 1e-3


@pytest.mark.parametrize("nu", [0.33, 0.5, 0.588, 0.75])
def test_saw_nu_hits_the_rms_distance_it_was_given(nu):
    """r0 is fixed so that sqrt(<r^2>) = r_rms; that is the whole scaling."""
    axis = np.linspace(1e-3, 600.0, 60001)
    pdf = saw_nu(axis, 50.0, nu)
    assert np.all(np.isfinite(pdf)) and pdf.min() >= 0.0
    rms = np.sqrt(np.trapezoid(pdf * axis ** 2, axis)
                  / np.trapezoid(pdf, axis))
    assert rms == pytest.approx(50.0, rel=1e-3)


def test_saw_nu_rejects_an_impossible_exponent():
    axis = np.linspace(1.0, 100.0, 50)
    assert np.all(saw_nu(axis, 50.0, 0.0) == 0.0)
    assert np.all(saw_nu(axis, 50.0, 1.0) == 0.0)
    assert np.all(saw_nu(axis, -1.0, 0.588) == 0.0)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
