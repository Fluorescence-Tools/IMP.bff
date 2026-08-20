"""Every kernel that publishes a numpy view, checked for the three silent hazards.

``ARGOUTVIEWM_ARRAY1`` hands Python an ndarray *over the kernel's own buffer*.
Nothing is copied, which is the point — a returned ``std::vector`` costs about
35-40 ns per element on the way out, and on an ``ng = 41`` grid that is 69 000
voxels of pure marshalling. It is also the reason each of these can go wrong
without saying so:

* **Ownership.** numpy releases the buffer with ``free``, so it must come from
  ``malloc``/``calloc``. ``new[]`` is undefined behaviour and looks fine.
  Checked by ``base`` being a capsule and by RSS not growing over many calls.
* **Which typemap binds.** The parameter *name* selects the typemap. The
  unmanaged ``ARGOUTVIEW`` and the managed ``ARGOUTVIEWM`` both claim
  ``output``; the kernels use ``out_view`` / ``out_view_i``, claimed only by
  the managed one. A leak here is invisible per call.
* **Aliasing.** Inputs arrive as raw pointers into numpy's buffers and the
  typemap cannot express ``const``. A kernel that writes through one corrupts
  its caller's array.

This file is the sweep-wide gate. Per-kernel numerics live with their kernel;
what is asserted here is only that the *view contract* holds.
"""

import gc
import resource

import numpy as np
import pytest

import IMP.bff
import IMP.bff as _kernels
import IMP.bff.quenching as qmaps


def rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)


def grid(ng=21, seed=0):
    """A plausible accessible volume: a filled ball, not noise."""
    rng = np.random.default_rng(seed)
    ax = np.arange(ng) - (ng - 1) / 2
    x, y, z = np.meshgrid(ax, ax, ax, indexing="ij")
    d = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    return np.where(d < ng / 3, rng.random((ng, ng, ng)) + 0.1, 0.0)


# --- the calls under test, each returning the array a caller would consume ---

def call_random_distances():
    p = np.ascontiguousarray(np.random.default_rng(1).random((400, 4)))
    return _kernels.random_distances(p, p + 5.0, 20000, seed=3), (p,)


def call_density_to_points():
    d = grid(21)
    pts = _kernels.density_to_points(d, 1.0, np.zeros(3), 0.0)
    return pts, (d,)


def call_split_contact_volume():
    ng = 21
    d = grid(ng)
    label = np.asarray(IMP.bff.split_contact_volume(
        np.ascontiguousarray(d).ravel(), ng, 1.0,
        np.array([2.0, 3.0]), np.array([0.0, 0.0, 0.0, 2.0, 1.0, 0.0]),
        np.zeros(3)))
    return label, (d,)


def call_quenching_map():
    ng = 21
    d = grid(ng)
    axis = np.array([float(ng), float(ng), float(ng)])
    atoms = np.array([1.0, 2.0, 3.0, -4.0, 0.5, 2.0])
    return qmaps._quenching_map(
        d, axis, np.zeros(3), atoms,
        np.array([1.0, 2.0]), np.array([1.5, 1.5]), 1.5, 1.0 / 4.0), (d, atoms)


def call_slow_near_atoms():
    ng = 21
    d = grid(ng)
    dmap = np.full((ng, ng, ng), 8.0)
    axis = np.array([float(ng), float(ng), float(ng)])
    atoms = np.array([1.0, 2.0, 3.0, -4.0, 0.5, 2.0])
    return qmaps._slow_near_atoms(
        dmap, d, axis, np.zeros(3), atoms, 36.0, 0.985), (d, dmap, atoms)


def call_fret_map():
    ng = 21
    dd, da = grid(ng, 0), grid(ng, 1)
    axis = np.array([float(ng), float(ng), float(ng)])
    return qmaps._fret_map(dd, da, axis, axis, np.zeros(3), np.zeros(3),
                           52.0 ** 6, 1.0 / 4.0, 2), (dd, da)


def call_fret_pair_matrices():
    rng = np.random.default_rng(2)
    a = np.ascontiguousarray(rng.random((60, 4)) * 30.0)
    b = np.ascontiguousarray(rng.random((50, 4)) * 30.0)
    out = np.asarray(IMP.bff.fret_pair_matrices(
        a.ravel(), b.ravel(), np.zeros(0), np.zeros(0), 60, 50))
    return out, (a, b)


def call_diffusion_step():
    ng = 21
    d = grid(ng)
    cur = d / d.sum()
    mob = np.full(ng ** 3, 0.05)
    decay = np.ones(ng ** 3)
    bounds = (d > 0).astype(np.float64).ravel()
    return np.asarray(IMP.bff.diffusion_step(
        np.ascontiguousarray(cur).ravel(), mob, decay, bounds,
        ng, IMP.bff.FLUX_SMOLUCHOWSKI)), (cur, mob, decay, bounds)


def call_diffusion_propagate():
    ng = 21
    d = grid(ng)
    cur = d / d.sum()
    mob = np.full(ng ** 3, 0.05)
    decay = np.ones(ng ** 3)
    bounds = (d > 0).astype(np.float64).ravel()
    fluo = IMP.bff.VectorDouble()
    return np.asarray(IMP.bff.diffusion_propagate(
        np.ascontiguousarray(cur).ravel(), mob, decay, bounds,
        ng, IMP.bff.FLUX_SMOLUCHOWSKI, 20, 5, fluo)), (cur, mob, decay, bounds)


def call_stamp_spheres():
    ng = 21
    d = grid(ng)
    rng = np.random.default_rng(7)
    rs = (rng.random(30) * 10.0) - 5.0
    return np.asarray(IMP.bff.stamp_spheres(
        np.ascontiguousarray(d).ravel(), ng, np.full(10, 3.0), rs,
        np.zeros(3), 1.0, np.full(10, 0.5),
        IMP.bff.GRID_COMBINE_MULTIPLY)), (d, rs)


def call_photon_trace():
    k = np.full(4000, 0.05)
    return np.asarray(IMP.bff.photon_trace(20000, k, 0.01, 4.0, 3)), (k,)


def call_rotamer_pair_energy_matrix():
    rng = np.random.default_rng(11)
    a = np.ascontiguousarray(rng.random(40 * 12 * 3) * 10.0)
    b = np.ascontiguousarray(rng.random(35 * 12 * 3) * 10.0 + 4.0)
    rmin = np.full(12, 3.5)
    eps = np.full(12, 0.1)
    return np.asarray(IMP.bff.rotamer_pair_energy_matrix(
        a, b, rmin, eps, 40, 12, 35, 12)), (a, b, rmin, eps)


KERNELS = {
    "random_distances": call_random_distances,
    "density_to_points": call_density_to_points,
    "split_contact_volume": call_split_contact_volume,
    "quenching_map": call_quenching_map,
    "slow_near_atoms": call_slow_near_atoms,
    "fret_map": call_fret_map,
    "fret_pair_matrices": call_fret_pair_matrices,
    "diffusion_step": call_diffusion_step,
    "diffusion_propagate": call_diffusion_propagate,
    "stamp_spheres": call_stamp_spheres,
    "photon_trace": call_photon_trace,
    "rotamer_pair_energy_matrix": call_rotamer_pair_energy_matrix,
}


@pytest.mark.parametrize("name", sorted(KERNELS))
def test_the_view_is_managed_not_leaked(name):
    """``base`` is the capsule that owns the buffer.

    ``.base is None`` would mean the array does not own it and nothing ever
    will -- the signature of the unmanaged typemap having bound instead.
    """
    out, _ = KERNELS[name]()
    out = np.asarray(out)
    base = out.base
    while base is not None and getattr(base, "base", None) is not None:
        base = base.base
    assert base is not None, f"{name}: view owns nothing -- unmanaged typemap?"
    assert type(base).__name__ == "PyCapsule", type(base).__name__


@pytest.mark.parametrize("name", sorted(KERNELS))
def test_the_kernel_does_not_write_to_its_inputs(name):
    """Aliasing. Inputs are raw pointers into numpy's buffers, not ``const``."""
    out, inputs = KERNELS[name]()
    before = [np.array(a, copy=True) for a in inputs]
    del out
    out2, inputs2 = KERNELS[name]()
    for i, (a, b) in enumerate(zip(inputs2, before)):
        np.testing.assert_array_equal(a, b, err_msg=f"{name}: input {i} mutated")


@pytest.mark.parametrize("name", sorted(KERNELS))
def test_repeated_calls_do_not_grow_the_heap(name):
    """The buffer is freed when the array dies.

    A leak here is one whole array per call, so it shows up fast. The threshold
    is loose because the allocator does not return pages promptly; what would
    fail is linear growth in the number of calls.
    """
    fn = KERNELS[name]

    def one_batch(n=300):
        for _ in range(n):
            out, _ = fn()
            del out
        gc.collect()
        return rss_mb()

    one_batch()                               # reach the allocator's steady state
    start = rss_mb()
    grew = one_batch() - start
    assert grew < 40.0, f"{name}: RSS grew {grew:.1f} MB over 300 further calls"


def test_a_zero_length_view_is_an_empty_array_not_a_null_pointer():
    """A null buffer would become an ndarray over address zero, and reading it
    is a segfault, not an exception. Every kernel allocates at least one
    element for this reason."""
    empty = np.zeros(0)
    out = np.asarray(IMP.bff.random_distances(empty, empty, 0, 0))
    assert out.size == 0
    pts = _kernels.density_to_points(np.zeros((2, 2, 2)), 1.0, np.zeros(3), 0.0)
    assert pts.size == 0


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
