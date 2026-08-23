"""Tests for ``IMP.bff.representation.distance`` — polymer chain models."""

from __future__ import annotations

import numpy as np
import pytest

from IMP.bff import (
    gaussian_chain_ree,
    gaussian_chain,
    worm_like_chain,
    worm_like_chain_linker,
)


class TestGaussianChain:
    """Gaussian (ideal) chain end-to-end distributions."""

    def test_ree_formula(self):
        ree = gaussian_chain_ree(segment_length=10.0, number_of_segments=9)
        assert ree == pytest.approx(30.0, abs=1e-9)

    def test_pdf_non_negative(self):
        r = np.linspace(0, 100, 101)
        pdf = gaussian_chain(r, segment_length=10.0, number_of_segments=9)
        assert np.all(pdf >= 0.0)

    def test_pdf_normalises_approximately(self):
        r = np.linspace(0, 200, 201)
        pdf = gaussian_chain(r, segment_length=10.0, number_of_segments=9)
        # Trapezoidal integration over r
        dr = r[1] - r[0]
        # numpy 2 removed np.trapz in favour of np.trapezoid
        trapezoid = getattr(np, "trapezoid", None) or np.trapz
        integral = trapezoid(pdf, r)
        # Not perfectly normalised because the grid is finite; just check ballpark
        assert integral > 0.8


class TestWormLikeChain:
    """Becker–Rosa–Everaers worm-like chain solution."""

    def test_pdf_non_negative(self):
        r = np.linspace(0, 0.99, 50)
        pdf = worm_like_chain(r, kappa=1.0, chain_length=1.0, normalize=False)
        assert np.all(pdf >= 0.0)

    def test_normalisation(self):
        r = np.linspace(0, 0.99, 100)
        pdf = worm_like_chain(r, kappa=1.0, chain_length=1.0, normalize=True)
        assert pdf.sum() == pytest.approx(1.0, abs=1e-3)

    def test_linker_broadening_non_negative(self):
        r = np.linspace(0, 0.99, 50)
        pdf = worm_like_chain_linker(
            r, kappa=1.0, chain_length=1.0, sigma=0.1, normalize=True
        )
        assert np.all(pdf >= 0.0)

    def test_linker_broadening_normalises(self):
        r = np.linspace(0, 0.99, 100)
        pdf = worm_like_chain_linker(
            r, kappa=1.0, chain_length=1.0, sigma=0.1, normalize=True
        )
        assert pdf.sum() == pytest.approx(1.0, abs=1e-2)


# IMP runs every .py under test/ as a standalone script; a file of bare pytest
# functions would import cleanly and exit 0, reporting success without running
# a single assertion. Hand it to pytest so a failure here fails ctest.
if __name__ == "__main__":
    import sys
    try:
        import pytest
    except ImportError:
        print("pytest not installed; skipping", __file__)
        sys.exit(0)
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
