"""PRD-140: the propagation as one matrix exponential instead of n_steps sweeps.

``diffusion_propagate`` takes four thousand explicit Euler sweeps because
forward Euler on a diffusion operator is stable only for ``dt <= dg^2/(6D)``.
``diffusion_propagate_krylov`` takes none: this operator is constant over a
run, so the whole propagation is ``exp(-Lt)p0``, and it is symmetrisable, so
one Lanczos projection answers every reported time at once.

**The two do not agree exactly, and must not be tested as if they did.** The
stepping route carries a forward-Euler and a first-order Lie-splitting error
of its own -- about 2e-6 on these fixtures, measured against
``scipy.sparse.linalg.expm_multiply`` in ``benchmark/large_time_steps.py``.
The Krylov route does not carry it. So the tolerance here is the *stepping*
scheme's error, and the Krylov answer is the better of the two.

What is checked:

* both flux forms (the Ito one only after the change of variables
  ``q = sqrt(d) p``, which is what makes it symmetric -- get that wrong and
  the trace is silently off in the third digit),
* the density as well as the trace, since the second pass that rebuilds it
  reruns the recursion from stored coefficients and could drift from the first,
* that ``diffusion_trace_krylov`` gives the first pass's trace exactly,
* that raising the Krylov dimension makes it converge rather than wander,
* that the Ito form refuses rather than divides by zero where ``d`` vanishes
  inside the domain.
"""

import numpy as np
import pytest

from IMP.bff import (
    FLUX_ITO,
    FLUX_SMOLUCHOWSKI,
    diffusion_propagate,
    diffusion_propagate_krylov,
    diffusion_trace_krylov,
)

# What the stepping scheme is itself wrong by on these fields; the Krylov
# answer is nearer the truth, so this is a distance, not an error budget.
STEPPING_ERROR = 1e-5


def _fields(ng, seed=3):
    rng = np.random.default_rng(seed)
    i = np.arange(ng) - ng // 2
    X, Y, Z = np.meshgrid(i, i, i, indexing="ij")
    inside = (X**2 + Y**2 + Z**2) <= (ng // 2 - 2) ** 2
    cur = (inside / inside.sum()).astype(float)
    d = (0.10 + 0.05 * rng.random((ng, ng, ng))) * inside
    decay = np.exp(-(0.02 + 0.20 * (X**2 + Y**2 + Z**2) / ng**2) * 0.005) * inside + (~inside)
    return [a.ravel().tolist() for a in (cur, d, decay, inside.astype(float))]


class TestKrylovPropagate:
    def test_matches_the_stepping_route_both_flux_forms(self):
        ng, n_steps, n_out = 21, 2000, 100
        fields = _fields(ng)
        for flux in (FLUX_SMOLUCHOWSKI, FLUX_ITO):
            fl_step, dn_step = diffusion_propagate(*fields, ng, flux, n_steps, n_out)
            fl_kry, dn_kry = diffusion_propagate_krylov(
                *fields, ng, flux, n_steps, n_out, 128)
            fl_step, dn_step = np.asarray(fl_step), np.asarray(dn_step)
            fl_kry, dn_kry = np.asarray(fl_kry), np.asarray(dn_kry)
            assert len(fl_kry) == n_steps // n_out + 1
            assert np.abs(fl_kry - fl_step).max() / np.abs(fl_step).max() < STEPPING_ERROR, "trace, flux_form=%d" % flux
            assert np.abs(dn_kry - dn_step).max() / np.abs(dn_step).max() < STEPPING_ERROR, "density, flux_form=%d" % flux

    def test_the_trace_only_entry_is_the_full_one_without_the_density(self):
        ng = 17
        fields = _fields(ng)
        full, _ = diffusion_propagate_krylov(*fields, ng, FLUX_SMOLUCHOWSKI, 1000, 100, 64)
        trace = diffusion_trace_krylov(*fields, ng, FLUX_SMOLUCHOWSKI, 1000, 100, 64)
        np.testing.assert_array_equal(np.asarray(full), np.asarray(trace))

    def test_more_krylov_vectors_converge(self):
        """The sequence has to settle, not merely change. A basis that is too
        small is the one failure mode with no symptom other than a wrong
        answer, so the check is that the step from 96 to 128 is far smaller
        than the step from 16 to 32."""
        ng = 21
        fields = _fields(ng)
        traces = {m: np.asarray(diffusion_trace_krylov(
            *fields, ng, FLUX_SMOLUCHOWSKI, 2000, 100, m)) for m in (16, 32, 96, 128)}
        coarse = np.abs(traces[32] - traces[16]).max()
        fine = np.abs(traces[128] - traces[96]).max()
        assert fine < coarse * 1e-3
        assert fine < 1e-9

    def test_zero_density_propagates_to_zero(self):
        ng = 13
        cur, d, decay, bounds = _fields(ng)
        fl = diffusion_trace_krylov([0.0] * len(cur), d, decay, bounds,
                                    ng, FLUX_SMOLUCHOWSKI, 500, 100, 32)
        np.testing.assert_allclose(np.asarray(fl), 0.0)

    def test_ito_refuses_where_the_mobility_vanishes(self):
        """`q = sqrt(d) p` is what makes the Ito form symmetric, so a voxel
        that is inside the domain with zero mobility has no image. Refusing is
        the honest answer; dividing by zero would return a NaN field."""
        ng = 13
        cur, d, decay, bounds = _fields(ng)
        d = list(d)
        centre = ((ng // 2) * ng + ng // 2) * ng + ng // 2
        assert bounds[centre] != 0.0
        d[centre] = 0.0
        with pytest.raises(ValueError):
            diffusion_trace_krylov(cur, d, decay, bounds, ng, FLUX_ITO, 500, 100, 32)
        # the Smoluchowski form has no such requirement
        fl = diffusion_trace_krylov(cur, d, decay, bounds, ng,
                                    FLUX_SMOLUCHOWSKI, 500, 100, 32)
        assert np.all(np.isfinite(np.asarray(fl)))
