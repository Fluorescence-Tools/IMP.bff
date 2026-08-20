"""FRET over every donor/acceptor state pair — and the hazards of a numpy view.

The pair matrices are an ``N1 x N2`` problem, and the vectorised Python paid for
it twice: an ``(N1, N2, 3)`` array of separation vectors (which it returned as
``r_vectors``, and **nothing ever read**), and five full-size temporaries in the
efficiency expression.

The kernel forms each separation vector and discards it. At 900 x 800 with
dipoles: 44.3 ms -> 4.68 ms.

.. warning::
   Most of that win is **not** the arithmetic. It is the return path. Handing a
   ``std::vector`` back makes SWIG build one Python float per element and numpy
   walk them back — **35-40 ns each**, which on a 400 x 350 matrix was ~20 ms
   against about 1 ms of actual work. Returning a numpy *view* over the kernel's
   own buffer costs nothing.

   Views are free and dangerous, and every hazard below is silent:

   * **Ownership.** ``ARGOUTVIEWM`` releases the buffer with ``free``, so it
     must be ``malloc``/``calloc``, never ``new[]``. Mixing them is undefined
     behaviour and invisible from Python.
   * **Which typemap binds.** ``(double **output, int *n_output)`` is claimed by
     *both* the managed and the unmanaged typemap in ``IMP_bff.types.i``.
     Whichever is declared last wins, so reordering two lines turns every call
     into a leak of the whole array. The kernel uses ``out_view``, a name only
     the managed typemap claims.
   * **Aliasing.** Inputs arrive as ``double*`` into numpy's own buffer — the
     typemap cannot express ``const``. A kernel that writes to one silently
     corrupts the caller's array.

   The tests below check all three, because none of them announces itself.
"""

import gc
import resource

import numpy as np
import pytest

import IMP.bff
from IMP.bff.photophysics import kappa2_from_dipoles
from IMP.bff import fret_pair_efficiencies, fret_pair_geometry


def _reference(p1, w1, p2, w2, mu1=None, mu2=None):
    """The pre-port implementation, verbatim."""
    p1 = np.asarray(p1, float)[:, :3]
    p2 = np.asarray(p2, float)[:, :3]
    weight = np.outer(np.asarray(w1, float), np.asarray(w2, float))
    total = weight.sum()
    if total > 0:
        weight = weight / total
    r_vectors = p1[:, None, :] - p2[None, :, :]
    r = np.linalg.norm(r_vectors, axis=2)
    if mu1 is None or mu2 is None:
        kappa2 = np.full(r.shape, 2.0 / 3.0)
    else:
        m1 = np.asarray(mu1, float)
        m2 = np.asarray(mu2, float)
        m1 = m1 / np.linalg.norm(m1, axis=1, keepdims=True)
        m2 = m2 / np.linalg.norm(m2, axis=1, keepdims=True)
        kappa2 = kappa2_from_dipoles(m1, m2, r_vectors)
    return r, kappa2, weight


@pytest.mark.parametrize("oriented", [False, True])
@pytest.mark.parametrize("n1,n2", [(1, 1), (7, 3), (60, 45), (200, 180)])
def test_matches_the_python_it_replaced(oriented, n1, n2):
    rng = np.random.default_rng(3)
    p1 = rng.normal(0, 8, (n1, 3))
    p2 = rng.normal(50, 8, (n2, 3))
    w1, w2 = rng.random(n1), rng.random(n2)
    mu1 = rng.normal(0, 1, (n1, 3)) if oriented else None
    mu2 = rng.normal(0, 1, (n2, 3)) if oriented else None

    want_r, want_k2, want_w = _reference(p1, w1, p2, w2, mu1, mu2)
    got = fret_pair_geometry(p1, w1, p2, w2, mu1, mu2)

    np.testing.assert_allclose(got["R"], want_r, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(got["kappa2"], want_k2, rtol=1e-11, atol=1e-13)
    np.testing.assert_allclose(got["weight"], want_w, rtol=1e-14)
    assert got["kappa2_avg"] == pytest.approx(float(np.sum(want_k2 * want_w)))


def test_efficiencies_match_the_python_they_replaced():
    rng = np.random.default_rng(5)
    n1, n2 = 40, 30
    p1, p2 = rng.normal(0, 9, (n1, 3)), rng.normal(45, 9, (n2, 3))
    w1, w2 = rng.random(n1), rng.random(n2)
    mu1, mu2 = rng.normal(0, 1, (n1, 3)), rng.normal(0, 1, (n2, 3))
    g = fret_pair_geometry(p1, w1, p2, w2, mu1, mu2)
    r, kappa2, weight, ka = g["R"], g["kappa2"], g["weight"], g["kappa2_avg"]

    with np.errstate(divide="ignore", invalid="ignore"):
        ratio6 = np.power(r / 52.0, 6)
        rate_ratio = 1.5 * kappa2 / ratio6
        e_pair = 1.0 / (1.0 + 2.0 / 3.0 * ratio6 / kappa2)
        e_dyn1 = 1.0 / (1.0 + 2.0 / 3.0 / ka * ratio6)
    e_pair = np.nan_to_num(e_pair, nan=1.0, posinf=1.0)
    e_dyn1 = np.nan_to_num(e_dyn1, nan=1.0, posinf=1.0)
    rate_avg = float(np.sum(np.nan_to_num(rate_ratio, posinf=0.0) * weight))

    got = fret_pair_efficiencies(g, 52.0)
    np.testing.assert_allclose(got["E"], e_pair, rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(got["rate_ratio"], rate_ratio, rtol=1e-12)
    assert got["static"] == pytest.approx(float(np.sum(e_pair * weight)))
    assert got["dynamic1"] == pytest.approx(float(np.sum(e_dyn1 * weight)))
    assert got["dynamic2"] == pytest.approx(rate_avg / (rate_avg + 1.0))


def test_coincident_states_mean_complete_transfer():
    """R = 0 makes the efficiency 0/0 or 1/(1+0). Both mean E = 1."""
    p = np.zeros((1, 3))
    g = fret_pair_geometry(p, [1.0], p, [1.0])
    assert g["R"][0, 0] == 0.0
    e = fret_pair_efficiencies(g, 52.0)
    assert e["E"][0, 0] == pytest.approx(1.0)
    assert np.isinf(e["rate_ratio"][0, 0]), "an infinite rate at zero separation is true"


def test_r_vectors_is_gone():
    """It was the largest allocation in the call and had no reader."""
    g = fret_pair_geometry(np.zeros((2, 3)), [1, 1], np.ones((2, 3)), [1, 1])
    assert not hasattr(g, "r_vectors")
    with pytest.raises(KeyError):
        g["r_vectors"]
    for key in ("R", "kappa2", "weight", "kappa2_avg"):
        assert g[key] is not None


# --- the hazards of a numpy view --------------------------------------------

def test_the_view_is_managed_not_leaked():
    """`base` must be a capsule that frees the buffer.

    If the unmanaged typemap bound instead, this would leak the whole matrix on
    every call and nothing else would notice.
    """
    rng = np.random.default_rng(0)
    n1 = n2 = 300
    p1 = np.ascontiguousarray(rng.normal(0, 8, (n1, 3))).ravel()
    p2 = np.ascontiguousarray(rng.normal(50, 8, (n2, 3))).ravel()
    empty = np.empty(0)

    first = IMP.bff.fret_pair_matrices(p1, p2, empty, empty, n1, n2)
    assert isinstance(first, np.ndarray)
    assert not first.flags.owndata
    assert first.base is not None, "an unmanaged view has no owner and leaks"

    # A leak is 1.4 MB *per call*, so it is linear and unmissable; what is not
    # unmissable is the allocator's first few hundred blocks, which raise the
    # high-water mark once and never again. Measuring the **second** batch
    # separates the two: a leak still shows ~290 MB there, a one-off step
    # shows nothing. Measuring the first batch made this fail only when the
    # whole suite ran, and only because the heap was already fragmented.
    per_call_mb = n1 * n2 * 2 * 8 / 1e6

    def one_batch():
        for _ in range(200):
            del_me = IMP.bff.fret_pair_matrices(p1, p2, empty, empty, n1, n2)
            del del_me
        gc.collect()
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)

    one_batch()                                   # reach steady state
    before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)
    grew = one_batch() - before
    assert grew < per_call_mb * 10, (
        f"200 further calls of {per_call_mb:.1f} MB grew RSS by {grew:.1f} MB "
        "-- the view is not being freed")


def test_the_kernel_does_not_write_to_its_inputs():
    """Inputs arrive as `double*` into numpy's buffer; the typemap cannot say
    `const`. A kernel that wrote to one would corrupt the caller's array."""
    rng = np.random.default_rng(1)
    n1, n2 = 50, 40
    p1 = np.ascontiguousarray(rng.normal(0, 8, (n1, 3))).ravel()
    p2 = np.ascontiguousarray(rng.normal(50, 8, (n2, 3))).ravel()
    mu1 = np.ascontiguousarray(rng.normal(0, 1, (n1, 3))).ravel()
    mu2 = np.ascontiguousarray(rng.normal(0, 1, (n2, 3))).ravel()
    keep = [a.copy() for a in (p1, p2, mu1, mu2)]

    IMP.bff.fret_pair_matrices(p1, p2, mu1, mu2, n1, n2)

    for got, want, name in zip((p1, p2, mu1, mu2), keep,
                               ("points1", "points2", "mu1", "mu2")):
        np.testing.assert_array_equal(got, want, err_msg=f"{name} was modified")


def test_dipoles_need_not_be_normalised():
    """They are normalised inside, so scaling one must not change kappa^2."""
    rng = np.random.default_rng(4)
    p1, p2 = rng.normal(0, 8, (5, 3)), rng.normal(40, 8, (4, 3))
    w = np.ones(5), np.ones(4)
    mu1, mu2 = rng.normal(0, 1, (5, 3)), rng.normal(0, 1, (4, 3))
    a = fret_pair_geometry(p1, w[0], p2, w[1], mu1, mu2)
    b = fret_pair_geometry(p1, w[0], p2, w[1], mu1 * 17.0, mu2 * 0.03)
    np.testing.assert_allclose(a["kappa2"], b["kappa2"], rtol=1e-12)


def test_an_empty_ensemble_returns_empty_matrices():
    g = fret_pair_geometry(np.zeros((0, 3)), [], np.ones((3, 3)), [1, 1, 1])
    assert g["R"].shape == (0, 3)
    assert g["kappa2"].shape == (0, 3)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
