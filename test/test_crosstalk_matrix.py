"""The excitation / emission crosstalk matrices and their algebra.

Moved here from ChiSurf's ``core/fluorescence/crosstalk.py`` (owner,
2026-09-04: the definition lives in bff, not in chisurf). What is pinned:

* the *label* semantics — ``select`` re-orders, sub-samples, and zero-fills a
  requested label the matrix does not carry, because a missing element of a
  light path is a dark element, not a broken one;
* the algebra against the references a Python caller would reach for —
  ``numpy.linalg.lstsq`` for the plain inverse, the closed Tikhonov form for
  the ridge, ``scipy.optimize.nnls`` for the non-negative one;
* the photon shuffle's contract — integer, non-negative, exactly
  count-preserving, seed-reproducible, and unbiased against the NNLS
  estimate on average;
* the Eigen-5 trap this module already fell into once: a ``Block`` wrapped in
  ``asDiagonal()`` inside a matrix product silently mis-evaluates under
  Eigen 5.0.1, so the posterior broadcast must stay an array expression.
"""

import numpy as np
import pytest
import scipy.optimize as so

import IMP.bff as bff


# --------------------------------------------------------- the labelled matrix

PAYLOAD_ROWS = ["D", "A"]
PAYLOAD_COLUMNS = ["green", "red"]
PAYLOAD_VALUES = [0.9, 0.1, 0.05, 0.95]


def test_select_reorders_and_zero_fills():
    m = bff.CrosstalkMatrix(PAYLOAD_ROWS, PAYLOAD_COLUMNS, PAYLOAD_VALUES)
    sel = m.select(["A", "D", "X"], ["red", "green"])
    assert list(sel.get_rows()) == ["A", "D", "X"]
    assert list(sel.get_columns()) == ["red", "green"]
    assert np.allclose(
        np.asarray(sel.get_values()).reshape(3, 2),
        [[0.95, 0.05], [0.10, 0.90], [0.00, 0.00]],
    )


def test_get_by_label_zero_for_missing():
    m = bff.CrosstalkMatrix(PAYLOAD_ROWS, PAYLOAD_COLUMNS, PAYLOAD_VALUES)
    assert m.get("D", "green") == pytest.approx(0.9)
    assert m.get("D", "violet") == 0.0
    assert m.get("X", "green") == 0.0


def test_scale_row_and_column():
    m = bff.CrosstalkMatrix(PAYLOAD_ROWS, PAYLOAD_COLUMNS, PAYLOAD_VALUES)
    m.scale_row("A", 0.5)          # a quantum yield
    m.scale_column("red", 0.8)     # a detection efficiency
    assert m.get("A", "red") == pytest.approx(0.95 * 0.5 * 0.8)
    assert m.get("A", "green") == pytest.approx(0.05 * 0.5)
    assert m.get("D", "red") == pytest.approx(0.1 * 0.8)
    assert m.get("D", "green") == pytest.approx(0.9)


def test_constructor_checks_the_shape():
    with pytest.raises(Exception):
        bff.CrosstalkMatrix(["D"], ["green", "red"], [1.0])


# ------------------------------------------------------------------- the algebra

M = np.array([[0.9, 0.1], [0.2, 0.8]])


def test_forward_mixing_matches_the_convention():
    sources = np.array([3.0, 5.0])
    measured = np.asarray(bff.crosstalk_apply_mixing(M.ravel(), 2, 2, sources))
    assert np.allclose(measured, M.T @ sources)


def test_plain_inverse_is_the_least_squares_solution():
    rng = np.random.default_rng(7)
    for _ in range(10):
        n_src = int(rng.integers(2, 5))
        n_det = n_src + int(rng.integers(0, 3))
        mx = rng.uniform(0.0, 1.0, size=(n_src, n_det))
        k = int(rng.integers(1, 6))
        src = rng.uniform(0.1, 2.0, size=(n_src, k))
        meas = mx.T @ src
        rec = np.asarray(
            bff.crosstalk_invert_mixing(mx.ravel(), n_src, n_det,
                                        meas.ravel(), False, 0.0))
        ref = np.linalg.lstsq(mx.T, meas, rcond=None)[0]
        assert np.allclose(rec.reshape(n_src, k), ref, atol=1e-8)


def test_ridge_is_the_tikhonov_closed_form():
    rng = np.random.default_rng(11)
    lam = 0.3
    mx = rng.uniform(0.1, 1.0, size=(3, 4))
    meas = rng.uniform(0.0, 10.0, size=(4, 5))
    rec = np.asarray(
        bff.crosstalk_invert_mixing(mx.ravel(), 3, 4, meas.ravel(),
                                    False, lam)).reshape(3, 5)
    a = mx.T
    ref = np.linalg.solve(a.T @ a + lam * np.eye(3), a.T @ meas)
    assert np.allclose(rec, ref, atol=1e-8)


def test_nonneg_matches_scipy_and_never_returns_negative():
    rng = np.random.default_rng(13)
    mx = rng.uniform(0.1, 1.0, size=(2, 3))
    a = mx.T
    for _ in range(10):
        meas = np.abs(a @ rng.uniform(0.0, 2.0, size=2)) * \
            rng.choice([-1.0, 1.0], size=3)
        rec = np.asarray(
            bff.crosstalk_invert_mixing(mx.ravel(), 2, 3, meas, True, 0.0))
        assert np.all(rec >= 0.0)
        assert np.allclose(rec, so.nnls(a, meas)[0], atol=1e-6)


def test_nonneg_survives_ridge_augmentation():
    # the augmentation rows carry zero targets; zeroing them inside the loop
    # would wipe the measurement when there is no augmentation at all
    mx = np.array([[0.9, 0.1], [0.2, 0.8]])
    meas = np.array([1000.0, 500.0])
    rec = np.asarray(
        bff.crosstalk_invert_mixing(mx.ravel(), 2, 2, meas, True, 0.5))
    a = mx.T
    ref = so.nnls(np.vstack([a, np.sqrt(0.5) * np.eye(2)]),
                  np.concatenate([meas, np.zeros(2)]))[0]
    assert np.allclose(rec, ref, atol=1e-6)


# --------------------------------------------------------------- photon shuffle

def test_shuffle_is_integer_nonneg_and_count_preserving():
    counts = np.array([1000.0, 500.0])
    out = np.asarray(bff.crosstalk_shuffle_unmix(M.ravel(), 2, 2, counts, [], 1))
    assert np.all(out >= 0) and np.all(out == np.round(out))
    assert out.sum() == counts.sum()


def test_shuffle_identity_is_lossless():
    eye = np.eye(2)
    counts = np.array([1000.0, 500.0])
    out = np.asarray(bff.crosstalk_shuffle_unmix(eye.ravel(), 2, 2, counts, [], 1))
    assert np.array_equal(out, counts)


def test_shuffle_is_seed_reproducible_and_unbiased():
    big = np.array([100000.0, 100000.0])
    reps = [np.asarray(bff.crosstalk_shuffle_unmix(M.ravel(), 2, 2, big, [], i))
            for i in range(60)]
    for i, r in enumerate(reps):
        assert np.array_equal(
            r, np.asarray(bff.crosstalk_shuffle_unmix(M.ravel(), 2, 2, big, [], i)))
    mean = np.mean(reps, axis=0)
    ref = np.asarray(bff.crosstalk_invert_mixing(M.ravel(), 2, 2, big, True, 0.0))
    assert np.all(np.abs(mean - ref) / ref < 0.02)


def test_shuffle_honours_supplied_abundances():
    eye = np.eye(2)
    counts = np.array([1000.0, 500.0])
    abundances = np.array([1000.0, 500.0])
    out = np.asarray(
        bff.crosstalk_shuffle_unmix(eye.ravel(), 2, 2, counts, abundances, 3))
    assert np.array_equal(out, counts)


def test_shuffle_abundance_size_is_checked():
    with pytest.raises(Exception):
        bff.crosstalk_shuffle_unmix(M.ravel(), 2, 2, np.array([1.0, 2.0]),
                                    np.array([1.0, 2.0, 3.0]), 0)


# --------------------------------------------------- the Eigen-5 diagonal trap

def test_posterior_broadcast_stays_an_array_expression():
    """``a * b.col(d).asDiagonal()`` mis-evaluates under Eigen 5.0.1.

    The diagonal wraps a Block, and the product silently degrades to
    ``a * b(0, d)``: with the identity matrix and abundances ``[1000, 500]``
    the detector-0 posterior reads [2/3, 1/3] -- the abundance ratio -- and
    the shuffle leaks 1/3 of a donor photon into the acceptor. The module
    uses an array broadcast; this pins the arithmetic it replaced.
    """
    eye = np.eye(2)
    counts = np.array([1000.0, 500.0])
    abundances = np.array([1000.0, 500.0])
    for seed in range(20):
        out = np.asarray(
            bff.crosstalk_shuffle_unmix(eye.ravel(), 2, 2, counts, abundances,
                                        seed))
        # exact under the identity: no detector-0 photon may leak into
        # source 1 (the broken product leaks about a third of them)
        assert np.array_equal(out, counts), (seed, out)
