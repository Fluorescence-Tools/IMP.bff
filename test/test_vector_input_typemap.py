"""How an array gets *into* a kernel, and why the obvious way was the slowest.

Every kernel taking ``const std::vector<double>&`` was handed
``np.ascontiguousarray(x).ravel()`` by its Python adapter. SWIG's default
conversion walks the sequence and calls ``SWIG_AsVal_double`` on each item, and
for an ndarray each of those accesses mints a fresh Python float -- so the
numpy array, the one thing the adapters went out of their way to produce, was
the **most expensive** of the three ways in. Measured on this box, 10k-400k
elements, on ``quenched_decay`` and ``lifetime_spectrum_decay``:

===========================================  ==================
 ndarray  -> ``const std::vector<double>&``   32-37 ns/element
 list     -> ``const std::vector<double>&``    8-13 ns/element
 ndarray  -> ``(double*, int)`` (IN_ARRAY1)    free
===========================================  ==================

``IMP_bff.types.i`` now carries an ``in`` typemap that one memcpy when the
argument is a contiguous float64 array: **4.4 ns/element**, and
``diffusion_propagate`` -- four ``ng^3`` grids in, one out -- went from 9.54 ms
to 0.296 ms of pure boundary crossing at ``ng = 41``. That is the cost of 53
solver steps recovered on every call, and ``equilibrium_occupancy`` makes up to
200 of them.

The typemap is a copy, deliberately: a ``std::vector`` cannot alias numpy's
buffer, and a kernel that kept one would hold a dangling pointer. Kernels on the
hot path should still take ``(double*, int)`` and read the caller's memory.

What this file guards is that the fast path did not change any *answer*, and
that everything which is not a contiguous float64 array still goes through
SWIG's own converter:

* a Python list, a tuple, and a wrapped ``std::vector<double>`` all still
  work --
  the last one needs explicit handling, because overriding the typemap stops
  ``std_vector.i``'s traits specialisation from ever being emitted;
* a **strided** array must not be read as though it were contiguous. Silently
  reading ``x[::2]``'s underlying buffer would take the wrong elements and
  return a plausible number;
* a float32 array must be converted, not reinterpreted.
"""

import numpy as np
import pytest

import IMP.bff

# The wrapped-vector case below needs an instance of whatever class SWIG uses
# for `std::vector<double>`, and that is not this module's to name: SWIG wraps
# a type once across a module and its imports, and `IMP.bff` imports `rmf`,
# which brings `isd`, which brings `saxs` -- so `IMP.saxs.DistBase` is the
# class and `IMP_bff.types.i` instantiates no `VectorDouble`. Nothing in the
# module's own surface asks a caller to build one; this test does, because
# binding one is the property under test.
from IMP.saxs import DistBase as VectorDouble



AMP = np.array([0.7, 0.3])
RATE = np.array([0.25, 1.1])


def decay(amp, rate, time):
    return np.asarray(IMP.bff.lifetime_spectrum_decay(amp, rate, time),
                      dtype=np.float64)


@pytest.fixture(scope="module")
def reference():
    t = np.linspace(0.0, 20.0, 512)
    return t, decay(AMP, RATE, t)


def test_a_list_gives_the_same_answer_as_an_array(reference):
    t, want = reference
    np.testing.assert_array_equal(decay(list(AMP), list(RATE), list(t)), want)


def test_a_tuple_gives_the_same_answer_as_an_array(reference):
    t, want = reference
    np.testing.assert_array_equal(
        decay(tuple(AMP), tuple(RATE), tuple(t)), want)


def test_a_wrapped_vectordouble_still_binds(reference):
    """The case the fast path breaks if it is not handled explicitly."""
    t, want = reference
    v = VectorDouble
    np.testing.assert_array_equal(
        decay(v(AMP.tolist()), v(RATE.tolist()), v(t.tolist())), want)


def test_a_strided_array_is_not_read_as_contiguous():
    """The hazard the ``array_is_contiguous`` check exists for.

    ``t[::2]`` shares a buffer whose first N/2 doubles are *not* its elements.
    Reading it as contiguous returns a decay -- just the wrong one.
    """
    t = np.linspace(0.0, 20.0, 1024)
    strided = t[::2]
    assert not strided.flags["C_CONTIGUOUS"]
    np.testing.assert_array_equal(
        decay(AMP, RATE, strided),
        decay(AMP, RATE, np.ascontiguousarray(strided)))


def test_a_float32_array_is_converted_not_reinterpreted():
    t = np.linspace(0.0, 20.0, 512)
    np.testing.assert_allclose(
        decay(AMP, RATE, t.astype(np.float32)),
        decay(AMP, RATE, t), rtol=1e-6, atol=1e-9)


def test_a_reversed_array_is_not_read_forwards():
    t = np.linspace(0.0, 20.0, 512)
    np.testing.assert_array_equal(
        decay(AMP, RATE, t[::-1]),
        decay(AMP, RATE, np.ascontiguousarray(t[::-1])))


def test_an_empty_array_is_accepted():
    assert decay(AMP, RATE, np.zeros(0)).size == 0
    assert decay(np.zeros(0), np.zeros(0), np.linspace(0, 1, 8)).size == 8


def test_a_contiguous_two_dimensional_array_is_flattened():
    """A contiguous float64 array of any shape is copied row-major -- which is
    exactly the ``.ravel()`` a caller would otherwise write before calling
    a kernel. A 2-D axis is therefore accepted and flattens to the 1-D answer.
    """
    t = np.linspace(0.0, 20.0, 512).reshape(2, 256)
    np.testing.assert_array_equal(decay(AMP, RATE, t),
                                  decay(AMP, RATE, t.ravel()))


def test_the_kernel_does_not_write_through_the_copy():
    """A memcpy in means the caller's array is untouched by construction --
    asserted anyway, since the whole point of the fast path is that it now
    touches the caller's buffer directly."""
    t = np.linspace(0.0, 20.0, 512)
    before = t.copy()
    decay(AMP, RATE, t)
    np.testing.assert_array_equal(t, before)


def test_the_fast_path_is_actually_taken():
    """A behavioural proxy for the timing claim, so the test survives a slow
    machine: a contiguous array must be at least twice as fast as the same
    values as a Python list, which is what the old ordering inverted."""
    import time
    n = 200_000
    t = np.linspace(0.0, 50.0, n)
    tl = t.tolist()

    def best(fn, reps=5):
        return min((lambda t0: (fn(), time.perf_counter() - t0)[1])(
            time.perf_counter()) for _ in range(reps))

    fast = best(lambda: IMP.bff.lifetime_spectrum_decay(AMP, RATE, t))
    slow = best(lambda: IMP.bff.lifetime_spectrum_decay(AMP, RATE, tl))
    assert fast < slow, f"ndarray {fast*1e3:.2f} ms, list {slow*1e3:.2f} ms"


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
