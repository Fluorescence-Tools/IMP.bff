"""The GPU backend, where there is one.

`libimp_bff_wgpu` is a plugin: found at import, used if it loads, and silently
absent otherwise. So these tests skip rather than fail where no GPU is
reachable, which is most of CI -- but where one *is* reachable they are the
only thing standing between a wrong kernel and a plausible-looking answer.

The tolerance is not roundoff. The GPU runs the stencil in f32 with f16
weights where that is sound, against the CPU's f64, so the two agree to about
1e-6 on the trace and 3e-6 on the density. See
okf/validation/gpu_diffusion_is_worth_it.md for where those numbers come from
and `use_f16_for` in the plugin for why one flux form gets f16 and the other
does not.
"""

import numpy as np
import pytest

import IMP.bff

pytestmark = pytest.mark.skipif(
    IMP.bff.get_compute_backend_name() == "cpu",
    reason="no GPU backend loaded (%s)" % (IMP.bff.get_compute_backend_error() or "none found"))

TRACE_TOL = 5e-6
DENSITY_TOL = 2e-5


def _fields(ng, seed=3):
    rng = np.random.default_rng(seed)
    i = np.arange(ng) - ng // 2
    X, Y, Z = np.meshgrid(i, i, i, indexing="ij")
    inside = (X**2 + Y**2 + Z**2) <= (ng // 2 - 2) ** 2
    cur = (inside / inside.sum()).astype(float)
    d = (0.10 + 0.05 * rng.random((ng, ng, ng))) * inside
    decay = np.exp(-(0.02 + 0.20 * (X**2 + Y**2 + Z**2) / ng**2) * 0.005) * inside + (~inside)
    return [a.ravel().tolist() for a in (cur, d, decay, inside.astype(float))]


def _both(fields, ng, flux, n_steps, n_out):
    """The same call on the GPU and then on the CPU."""
    gpu = IMP.bff.diffusion_propagate(*fields, ng, flux, n_steps, n_out)
    IMP.bff.reset_compute_backend()
    try:
        cpu = IMP.bff.diffusion_propagate(*fields, ng, flux, n_steps, n_out)
    finally:
        IMP.bff.enable_gpu(quiet=True)
    return [np.asarray(x) for x in gpu], [np.asarray(x) for x in cpu]


def test_the_backend_names_the_adapter():
    """A backend that answered must say so, or somebody measures a CPU and
    reports a GPU."""
    name = IMP.bff.get_compute_backend_name()
    assert name.startswith("wgpu:"), name
    assert name.count(":") >= 2, "the name carries backend and adapter: %s" % name


@pytest.mark.parametrize("flux", [IMP.bff.FLUX_SMOLUCHOWSKI, IMP.bff.FLUX_ITO])
def test_it_agrees_with_the_cpu(flux):
    ng = 25
    fields = _fields(ng)
    (fl_g, dn_g), (fl_c, dn_c) = _both(fields, ng, flux, 1000, 100)
    assert len(fl_g) == 11
    assert np.abs(fl_g - fl_c).max() / np.abs(fl_c).max() < TRACE_TOL
    assert np.abs(dn_g - dn_c).max() / np.abs(dn_c).max() < DENSITY_TOL


def test_it_declines_a_grid_too_small_to_be_worth_the_transfer():
    """Below the size where the transfer pays, the backend returns non-zero
    and the CPU runs -- so the answer is the CPU's, bit for bit."""
    ng = 9
    fields = _fields(ng)
    (fl_g, dn_g), (fl_c, dn_c) = _both(fields, ng, IMP.bff.FLUX_SMOLUCHOWSKI, 200, 50)
    np.testing.assert_array_equal(fl_g, fl_c)
    np.testing.assert_array_equal(dn_g, dn_c)


def test_the_cpu_can_be_taken_back_and_returned():
    IMP.bff.reset_compute_backend()
    try:
        assert IMP.bff.get_compute_backend_name() == "cpu"
    finally:
        IMP.bff.enable_gpu(quiet=True)
    assert IMP.bff.get_compute_backend_name().startswith("wgpu:")


def test_the_density_outside_the_domain_stays_zero():
    """The GPU writes only the active voxels, so anything left outside would
    live for ever -- the CPU sweep clears the whole grid every step."""
    ng = 25
    cur, d, decay, bounds = _fields(ng)
    b = np.asarray(bounds)
    _, dn = IMP.bff.diffusion_propagate(cur, d, decay, bounds, ng,
                                        IMP.bff.FLUX_SMOLUCHOWSKI, 500, 100)
    assert np.abs(np.asarray(dn)[b == 0.0]).max() == 0.0
