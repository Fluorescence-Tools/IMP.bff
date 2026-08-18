"""Mean-field rotamer weights: the pair energies, computed once and in the right order.

Two things changed when this moved to C++, and only one of them is a port.

**The port.** Every conformer of one dye against every conformer of another is
an all-pairs problem with a bounding-box pre-filter, and it is now
``rotamer_pair_energy_matrix``. Single dye, 60 clusters x 14 atoms against 900
protein atoms: 27 ms -> 7 ms.

**The algorithm.** Those energies do not depend on the weights, so they are
computed once instead of inside the iteration. The Python rebuilt the entire
dye-dye energy matrix on all ten iterations for the same answer: two dyes,
20 x 18 clusters, 267 ms -> 6 ms, a factor of 47.

**And a bug.** The old code stored one parameter matrix per *unordered* dye pair
and used it for both orders. With ``d1 > d2`` that multiplied an
``(n_1, n_2)``-shaped distance matrix by ``(n_2, n_1)``-shaped parameters. The
flattened lengths match, so numpy broadcast it without complaint and paired
every atom with the wrong partner's parameters. On two dyes whose elements
differ sharply in radius the cross energy came out **98.9 % wrong** -- 6.35
where the answer is 577.75.
"""

import numpy as np
import pytest

import IMP.bff
from IMP.bff.scoring.mean_field import (
    _aabb_overlap,
    _build_cross_lj_params,
    _lj_energy_pairs,
    _pair_energy_matrix,
    rotamer_mean_field_weights,
    rotamer_mean_field_weights_multi_dye,
)


def _reference_single(rot, w0, prot, de, pe, K=1.0, n_iter=10, pad=3.5, rc=12.0):
    """The pre-port implementation, verbatim."""
    n = rot.shape[0]
    q = w0.copy().astype(float)
    rmin, eps = _build_cross_lj_params(de, pe)
    pbox = np.concatenate([prot.min(0) - pad, prot.max(0) + pad])
    boxes = np.concatenate([rot.min(1) - pad, rot.max(1) + pad], axis=1)
    E = np.zeros(n)
    for i in range(n):
        if _aabb_overlap(boxes[i], pbox):
            E[i] = _lj_energy_pairs(rot[i], prot, rmin, eps, r_cutoff=rc)
    for _ in range(n_iter):
        lq = np.log(np.maximum(q, 1e-300)) - K * E
        lq -= lq.max()
        qn = np.exp(lq)
        s = qn.sum()
        q = qn / s if s > 0 else np.ones(n) / n
    return q


def test_single_dye_matches_the_python_it_replaced():
    rng = np.random.default_rng(4)
    n_c, n_a, n_p = 40, 12, 500
    rot = rng.normal(0, 4, (n_c, n_a, 3))
    prot = rng.normal(0, 7, (n_p, 3))
    de = list(rng.choice(["C", "N", "O", "S", "H"], n_a))
    pe = list(rng.choice(["C", "N", "O", "S"], n_p))
    w0 = np.full(n_c, 1.0 / n_c)

    np.testing.assert_allclose(
        rotamer_mean_field_weights(rot, w0, prot, de, pe),
        _reference_single(rot, w0, prot, de, pe), rtol=1e-12, atol=1e-15)


def test_the_pair_matrix_is_symmetric_under_exchange():
    """E(a_i, b_j) computed from a's side must equal it from b's side.

    This is the property the old code broke: it reused one parameter matrix for
    both orders, which is only the transpose if the two atom orderings happen to
    agree.
    """
    rng = np.random.default_rng(7)
    ca = rng.normal(0, 3, (6, 4, 3))
    cb = rng.normal(0, 3, (5, 7, 3)) + 2.0
    ea = list(rng.choice(["C", "N", "O", "S"], 4))
    eb = list(rng.choice(["H", "C", "S"], 7))

    ab = _pair_energy_matrix(ca, cb, ea, eb, 12.0, 3.5)
    ba = _pair_energy_matrix(cb, ca, eb, ea, 12.0, 3.5)
    assert ab.shape == (6, 5) and ba.shape == (5, 6)
    np.testing.assert_allclose(ab, ba.T, rtol=1e-12, atol=1e-15)


def test_the_transposed_parameter_bug_is_measurable():
    """Kept as a test because it is the reason the helper takes ordered elements.

    If this ever stops differing, the parameter build has become
    order-insensitive and the ordered helper is no longer needed -- which would
    itself be worth knowing.
    """
    rng = np.random.default_rng(2)
    e1, e2 = ["H", "H", "S"], ["S", "H"]
    c1 = rng.normal(0, 1.0, (3, 3))
    c2 = rng.normal(0, 1.0, (2, 3)) + 1.5

    r12, p12 = _build_cross_lj_params(e1, e2)
    r21, p21 = _build_cross_lj_params(e2, e1)
    wrong = _lj_energy_pairs(c2, c1, r12, p12, r_cutoff=12.0)
    right = _lj_energy_pairs(c2, c1, r21, p21, r_cutoff=12.0)
    assert abs(wrong - right) / abs(right) > 0.5

    # and the matrix helper agrees with the *right* one
    m = _pair_energy_matrix(c2[None, :, :], c1[None, :, :], e2, e1, 12.0, 3.5)
    assert m[0, 0] == pytest.approx(right)


def test_multi_dye_weights_are_normalised_per_dye():
    rng = np.random.default_rng(9)
    prot = rng.normal(0, 7, (300, 3))
    pe = list(rng.choice(["C", "N", "O"], 300))
    rots = [rng.normal(0, 3, (14, 5, 3)), rng.normal(0, 3, (11, 8, 3)) + 5.0]
    els = [list(rng.choice(["C", "N", "O"], 5)), list(rng.choice(["S", "C", "H"], 8))]
    w0s = [np.full(14, 1 / 14), np.full(11, 1 / 11)]

    out = rotamer_mean_field_weights_multi_dye(rots, w0s, prot, els, pe)
    assert len(out) == 2
    for w, n in zip(out, (14, 11)):
        assert w.shape == (n,)
        assert w.sum() == pytest.approx(1.0)
        assert np.all(w >= 0.0) and np.isfinite(w).all()


def test_dyes_far_apart_do_not_influence_each_other():
    """The bounding-box filter should skip every cross pair, leaving the
    single-dye answer."""
    rng = np.random.default_rng(12)
    prot = rng.normal(0, 6, (200, 3))
    pe = list(rng.choice(["C", "N"], 200))
    a = rng.normal(0, 2, (10, 4, 3))
    b = rng.normal(0, 2, (9, 4, 3)) + 500.0        # far away
    els = [["C", "N", "O", "S"], ["C", "N", "O", "S"]]
    w0 = [np.full(10, 0.1), np.full(9, 1 / 9)]

    joint = rotamer_mean_field_weights_multi_dye([a, b], w0, prot, els, pe)
    alone = rotamer_mean_field_weights(a, w0[0], prot, els[0], pe)
    np.testing.assert_allclose(joint[0], alone, rtol=1e-12, atol=1e-15)


def test_a_conformer_pair_beyond_the_boxes_contributes_nothing():
    ca = np.zeros((1, 2, 3))
    cb = np.full((1, 2, 3), 1000.0)
    m = _pair_energy_matrix(ca, cb, ["C", "C"], ["C", "C"], 12.0, 3.5)
    assert m.shape == (1, 1) and m[0, 0] == 0.0


def test_the_energy_is_repulsive_only():
    """Beyond Rmin the LJ tail is attractive; the mean-field term drops it,
    because it asks what is blocked rather than what is bound."""
    close = _pair_energy_matrix(np.zeros((1, 1, 3)), np.array([[[1.0, 0.0, 0.0]]]),
                                ["C"], ["C"], 12.0, 3.5)
    far = _pair_energy_matrix(np.zeros((1, 1, 3)), np.array([[[6.0, 0.0, 0.0]]]),
                              ["C"], ["C"], 12.0, 3.5)
    assert close[0, 0] > 0.0, "inside Rmin the wall repels"
    assert far[0, 0] == 0.0, "outside Rmin nothing is added"


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
