"""Polymer chain distributions, and the linker convolution that used to raise.

``worm_like_chain_linker`` had **no test**, which is why it went on raising
``TypingError`` for a whole commit: it was numba-jitted and called
``normal_distribution``, which had become a C++ delegation numba cannot type.
The C++ port fixed it, and this file is what stops it recurring silently.
"""

import numpy as np
import pytest

import IMP.bff
from IMP.bff import (
    gaussian_chain, gaussian_chain_ree, worm_like_chain, worm_like_chain_linker,
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
    """
    pr = worm_like_chain(R, 0.4)
    g = lambda x, loc, s: np.exp(-((x - loc) ** 2) / (2 * s * s)) / (np.sqrt(2 * np.pi) * s)
    ref = sum(pr[i] * g(R, R[i], 6.0) for i in range(len(R)) if pr[i] != 0.0)
    ref = ref / ref.sum()
    np.testing.assert_allclose(worm_like_chain_linker(R, 0.4, sigma=6.0), ref, atol=1e-15)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
