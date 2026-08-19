/*
 * The field picture of a tethered dye, as a C++ object.
 *
 * The Python held four ng^3 grids and folded `dt` into two of them on every
 * call -- `d = D dt/dg^2` and `decay = exp(-k dt)` -- then handed all four
 * across the boundary. At ng = 41 that is 276 000 doubles per call, and
 * `equilibrium()` makes up to 200 of them.
 *
 * `run()` and `gradient()` each return one buffer rather than a tuple of
 * arrays, because a numpy view is one array; the shim below splits it, which
 * costs two slices of a view and no copy.
 */

IMP_SWIG_VALUE(IMP::bff, GridDiffusionSolver, GridDiffusionSolvers);

// The two result records live here rather than in `IMP.bff.sampling`, because
// the methods that build them are now here and importing the Python module
// from the generated one would be circular. `sampling` re-exports both, so
// `from IMP.bff.sampling import GridDiffusionResult` keeps working.
%pythoncode %{
from collections import namedtuple as _namedtuple

GridDiffusionResult = _namedtuple(
    "GridDiffusionResult", ["time", "fluorescence", "density"])
GridDiffusionResult.__doc__ = (
    "Time axis, surviving excited-state fraction, and the final density.")

GridDiffusionGradient = _namedtuple(
    "GridDiffusionGradient", ["d_diffusion", "d_rate", "d_density"])
GridDiffusionGradient.__doc__ = (
    "What GridDiffusionSolver.gradient() returns -- dL/dD, dL/dk and dL/dp0 "
    "per voxel, each (ng, ng, ng), in the units of the solver's inputs "
    "(A^2/ns, 1/ns, and the unnormalised initial density).")
%}

%feature("shadow") IMP::bff::GridDiffusionSolver::GridDiffusionSolver %{
def __init__(self, diffusion_map, bounds, density, rate_map=None,
             t_step=1.0, dg=1.0, check_stability=True,
             flux_form="smoluchowski"):
    """Propagate an excited-state density on a masked 3-D grid.

    The parameter order is the one this class has always had -- diffusion map
    first, then the domain, then the population -- and it is kept because
    callers pass positionally. It is *not* the order the C++ constructor takes,
    which groups the four grids in the order they are used.
    """
    import numpy as _np
    def _flat(a):
        return (_np.zeros(0) if a is None
                else _np.ascontiguousarray(_np.asarray(a, dtype=_np.float64)).ravel())
    if isinstance(flux_form, str):
        if flux_form not in ("smoluchowski", "ito"):
            raise ValueError(
                "flux_form must be 'smoluchowski' or 'ito', not %r" % (flux_form,))
        flux = (_IMP_bff.FLUX_SMOLUCHOWSKI if flux_form == "smoluchowski"
                else _IMP_bff.FLUX_ITO)
        name = flux_form
    else:
        flux = int(flux_form)
        name = "smoluchowski" if flux == _IMP_bff.FLUX_SMOLUCHOWSKI else "ito"
    d = _flat(diffusion_map)
    rates = _flat(rate_map)
    if rates.size == 0:
        rates = _np.zeros(d.size)
    _IMP_bff.GridDiffusionSolver_swiginit(self, _IMP_bff.new_GridDiffusionSolver(
        _flat(density), _flat(bounds), d, rates,
        float(dg), float(t_step), flux, bool(check_stability)))
    self._flux_form = name
%}

%include "IMP/bff/GridDiffusionSolver.h"

%extend IMP::bff::GridDiffusionSolver {
    %pythoncode %{
        def _cube(self, flat):
            ng = self.get_ng()
            return flat.reshape(ng, ng, ng)

        @property
        def density(self):
            return self._cube(_IMP_bff.GridDiffusionSolver_get_density(self))

        @density.setter
        def density(self, value):
            self.set_density(
                np.ascontiguousarray(np.asarray(value, dtype=np.float64)).ravel())

        @property
        def bounds(self):
            return self._cube(_IMP_bff.GridDiffusionSolver_get_bounds(self))

        @property
        def diffusion_map(self):
            return self._cube(_IMP_bff.GridDiffusionSolver_get_diffusion_map(self))

        @diffusion_map.setter
        def diffusion_map(self, value):
            self.set_diffusion_map(
                np.ascontiguousarray(np.asarray(value, dtype=np.float64)).ravel())

        @property
        def rate_map(self):
            return self._cube(_IMP_bff.GridDiffusionSolver_get_rate_map(self))

        @rate_map.setter
        def rate_map(self, value):
            self.set_rate_map(
                np.ascontiguousarray(np.asarray(value, dtype=np.float64)).ravel())

        @property
        def dg(self):
            return self.get_dg()

        @property
        def t_step(self):
            return self.get_t_step()

        @property
        def n_iterations(self):
            return self.get_n_iterations()

        @property
        def flux_form(self):
            return getattr(self, "_flux_form", "smoluchowski")

        def run(self, n_steps, n_out=10):
            """Integrate *n_steps* steps, reporting every *n_out*."""
            n_steps = int(n_steps)
            n_out = max(1, int(n_out))
            n_reports = n_steps // n_out + 1
            ng = self.get_ng()
            out = _IMP_bff.GridDiffusionSolver_run(self, n_steps, n_out)
            fluorescence = out[:n_reports]
            final = out[n_reports:].reshape(ng, ng, ng)
            time = np.arange(n_reports, dtype=np.float64) * self.t_step * n_out
            return GridDiffusionResult(time, fluorescence, final)

        def gradient(self, dL_dF, n_steps, n_out=10, density=None):
            """Gradient of a loss on the decay, every voxel at once."""
            ng = self.get_ng()
            flat = _IMP_bff.GridDiffusionSolver_gradient(
                self,
                np.ascontiguousarray(np.asarray(dL_dF, dtype=np.float64)).ravel(),
                int(n_steps), max(1, int(n_out)),
                np.zeros(0) if density is None else
                np.ascontiguousarray(np.asarray(density, dtype=np.float64)).ravel())
            g = flat.reshape(3, ng, ng, ng)
            return GridDiffusionGradient(g[0], g[1], g[2])

        def equilibrium(self, n_steps=20000, tolerance=1e-8, n_check=100):
            """Propagate with no decay until the occupancy stops moving."""
            ng = self.get_ng()
            return _IMP_bff.GridDiffusionSolver_equilibrium(
                self, int(n_steps), float(tolerance), int(n_check)
            ).reshape(ng, ng, ng)
    %}
}
