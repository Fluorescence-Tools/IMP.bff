"""Greedy Olga: the chi-squared tail, and the candidate scoring that never allocates.

Experiment planning after Olga -- repeatedly add the FRET pair that most reduces
the expected RMSD between the true structure and the one the measurements would
pick out. The cost is entirely in scoring every remaining candidate against every
ordered pair of ensemble frames.

The vectorised Python built a ``(chunk, n_frames, n_frames)`` array to do it: at
700 frames that is 251 MB per 64-candidate chunk, allocated only to be reduced
away. ``expected_rmsd_after_adding`` runs the reduction inside the loop, so its
working set is ``O(n_frames)`` per candidate.

Measured, single-threaded (see ``test_the_speed_figures_are_single_threaded``):

    n=700 m=120  ndof even   4111 ms -> 1311 ms   3.14x
    n=700 m=120  ndof odd    7346 ms -> 4444 ms   1.65x

The odd branch is the slower one because it needs ``erfc`` per element, which
numpy vectorises well and a scalar loop does not. Recorded rather than hidden:
the win there is the 251 MB, not the clock.
"""

import numpy as np
import pytest

import IMP.bff
from IMP.bff.fret import greedy_olga as go

scipy_special = pytest.importorskip("scipy.special")


# --- the chi-squared right tail ---------------------------------------------

@pytest.mark.parametrize("ndof", [1, 2, 3, 4, 5, 7, 12, 33, 100, 199, 201, 400, 801])
def test_the_tail_matches_the_reference_series_and_scipy(ndof):
    """Three independent routes to the same function must agree.

    The closed forms here, the Python they were ported from, and scipy's
    general ``gammaincc``. Agreement of all three is what makes this a check on
    the port rather than on one implementation's habits.
    """
    rng = np.random.default_rng(5)
    x = np.concatenate([[0.0], rng.uniform(0.0, 900.0, 2000)])
    got = go._chisq_rt_cdf(x, ndof)
    np.testing.assert_allclose(got, go._chisq_rt_cdf_python(x, ndof), atol=1e-12)
    np.testing.assert_allclose(got, scipy_special.gammaincc(0.5 * ndof, 0.5 * x),
                               atol=1e-12)


def test_the_tail_is_one_at_zero_and_falls_monotonically():
    for ndof in (1, 2, 5, 40):
        v = go._chisq_rt_cdf(np.array([0.0, 1.0, 5.0, 20.0, 100.0]), ndof)
        assert v[0] == 1.0
        assert np.all(np.diff(v) <= 0.0)
        assert np.all((v >= 0.0) & (v <= 1.0))


def test_the_large_ndof_branch_needs_both_series_and_continued_fraction():
    """The continued fraction alone returns 1.0 for every x below a+1.

    At ndof = 201 that is every chi-squared the algorithm ever sees, so taking
    only that branch -- as a first cut of the kernel did -- is not a small error
    at the edges but a wrong answer everywhere.
    """
    ndof = 201
    x = np.array([0.5, 5.0, 50.0, 150.0])          # all well below a + 1 = 101.5
    got = go._chisq_rt_cdf(x, ndof)
    assert np.all(got > 0.99), "these are far into the tail; Q must be near 1"
    assert not np.all(got == 1.0), "but not exactly 1 -- that was the bug"
    np.testing.assert_allclose(got, scipy_special.gammaincc(0.5 * ndof, 0.5 * x),
                               atol=1e-12)


# --- the candidate scorer ----------------------------------------------------

def _reference(effs, rmsds, chi2, inv, ndof, dw):
    """The vectorised Python, kept verbatim."""
    e = np.atleast_2d(effs.T.astype(np.float64))
    c = e[:, :, None] - e[:, None, :]
    np.square(c, out=c)
    c *= inv
    c += chi2
    w = go._chisq_rt_cdf_python(c, ndof)
    sw = w.sum(axis=-2)
    return (np.einsum('kij,ij->kj', w, rmsds) / (sw - 1.0 + dw)).mean(axis=-1)


def _symmetric(rng, n, hi):
    m = rng.uniform(0, hi, (n, n))
    m = 0.5 * (m + m.T)
    np.fill_diagonal(m, 0.0)
    return m


@pytest.mark.parametrize("ndof", [1, 3, 6, 11])
def test_scoring_matches_the_python_it_replaced(ndof):
    rng = np.random.default_rng(11)
    n, m = 60, 25
    effs = rng.random((n, m))
    rmsds = _symmetric(rng, n, 12.0)
    chi2 = _symmetric(rng, n, 5.0)
    got = np.asarray(IMP.bff.expected_rmsd_after_adding(
        np.ascontiguousarray(rmsds).ravel(), np.ascontiguousarray(chi2).ravel(),
        np.ascontiguousarray(effs.T).ravel(), 400.0, ndof, 0.99, n, m))
    np.testing.assert_allclose(got, _reference(effs, rmsds, chi2, 400.0, ndof, 0.99),
                               rtol=1e-11, atol=1e-13)


def test_a_candidate_that_separates_nothing_scores_worst():
    """A pair with the same efficiency in every frame adds no information, so
    the expected RMSD it leaves must be the largest."""
    rng = np.random.default_rng(3)
    n = 40
    rmsds = _symmetric(rng, n, 10.0)
    chi2 = np.zeros((n, n))
    flat = rng.random(n)
    effs = np.column_stack([np.full(n, 0.5), flat, np.linspace(0.0, 1.0, n)])
    got = np.asarray(IMP.bff.expected_rmsd_after_adding(
        np.ascontiguousarray(rmsds).ravel(), chi2.ravel(),
        np.ascontiguousarray(effs.T).ravel(), 400.0, 1, 0.99, n, 3))
    assert got[0] == max(got), "the constant-efficiency candidate must be worst"


def test_empty_inputs_do_not_crash():
    e = np.empty(0)
    assert list(IMP.bff.expected_rmsd_after_adding(e, e, e, 1.0, 1, 0.99, 0, 0)) == []


def test_the_greedy_selection_still_runs_end_to_end():
    rng = np.random.default_rng(8)
    n, m = 30, 12
    effs = rng.random((n, m))
    rmsds = _symmetric(rng, n, 8.0)
    out = go.select_informative_pairs(effs, rmsds, err=0.05, max_pairs=4)
    idx = out[0] if isinstance(out, tuple) else out
    idx = np.asarray(idx).ravel()
    assert idx.size == 4
    assert len(set(idx.tolist())) == 4, "greedy with unique_only must not repeat"
    assert np.all((idx >= 0) & (idx < m))


def test_the_speed_figures_are_single_threaded():
    """Every measurement recorded for this module was taken with OpenMP off.

    The kernels carry ``#pragma omp parallel for``, but IMP's CMake leaves
    ``OpenMP_CXX_FLAGS`` empty in this build, so the pragmas are inert. If this
    ever starts failing, the build gained OpenMP -- and every speed figure in
    this repository's commit messages became a lower bound.
    """
    assert IMP.bff.parallel_threads() >= 1
    if not IMP.bff.built_with_openmp():
        assert IMP.bff.parallel_threads() == 1


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
