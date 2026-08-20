/*
 * The particle picture of a tethered dye, as a C++ object.
 *
 * The Python held the occupancy grid, the mobility field, the rate map and the
 * trajectory as numpy arrays and drove the walk kernel from a
 * `ThreadPoolExecutor`. Every array crossed the boundary on every access, and
 * `sample_grid` -- tens of millions of indexed reads along the trajectory --
 * was 0.18 s of a profiled suite slice.
 *
 * The threading moved with it, to `std::thread` rather than OpenMP: this
 * build's `OpenMP_CXX_FLAGS` is empty, so every `#pragma omp` in the module is
 * inert, and the Python got real parallelism only because the kernel releases
 * the GIL. Losing it in the move would have been a regression that no test
 * would have reported.
 */

IMP_SWIG_VALUE(IMP::bff, DyeDiffusionSimulation, DyeDiffusionSimulations);

%pythoncode %{
def _as_cube(flat):
    """A flat grid as (ng, ng, ng), taking ng from its *own* length.

    Not from the occupancy grid's edge: the fields are stored separately and
    are only conventionally the same size, so swapping the volume for a
    smaller one must not make the rate map unreadable.
    """
    n = int(flat.size)
    ng = int(round(n ** (1.0 / 3.0)))
    for cand in (ng, ng - 1, ng + 1):
        if cand > 0 and cand ** 3 == n:
            return flat.reshape(cand, cand, cand)
    raise ValueError("a grid of %d values is not a cube" % n)
%}

// Every argument optional and by keyword, and `None` for an absent grid --
// which SWIG cannot express for a function with defaults.
%feature("shadow") IMP::bff::DyeDiffusionSimulation::DyeDiffusionSimulation %{
def __init__(self, density=None, dg=0.5, x0=None, slow_density=None,
             slow_factor_map=None, quenching_rate_map=None):
    import numpy as _np
    def _i(a):
        return (_np.zeros(0, dtype=_np.int32) if a is None
                else _np.ascontiguousarray(_np.asarray(a, dtype=_np.int32)).ravel())
    def _d(a):
        return (_np.zeros(0) if a is None
                else _np.ascontiguousarray(_np.asarray(a, dtype=_np.float64)).ravel())
    _IMP_bff.DyeDiffusionSimulation_swiginit(
        self, _IMP_bff.new_DyeDiffusionSimulation(
            _i(density), float(dg),
            _np.zeros(3) if x0 is None else _d(x0),
            _i(slow_density), _d(slow_factor_map), _d(quenching_rate_map)))
%}

%include "IMP/bff/DyeDiffusion.h"

%extend IMP::bff::DyeDiffusionSimulation {
    %pythoncode %{
        @property
        def trajectory(self):
            """``(n_frames, 3)`` in the structure's frame, or ``None``.

            ``None`` rather than an empty array before ``run()``: every caller
            tests it, and an empty ``(0, 3)`` would silently produce empty
            answers downstream instead of saying the walk has not happened.
            """
            flat = _IMP_bff.DyeDiffusionSimulation_get_trajectory(self)
            return None if flat.size == 0 else flat.reshape(-1, 3)

        @property
        def x0(self):
            return _IMP_bff.DyeDiffusionSimulation_get_x0(self)

        @property
        def n_frames(self):
            return self.get_n_frames()

        @property
        def n_accepted(self):
            return self.get_n_accepted()

        @property
        def n_rejected(self):
            return self.get_n_rejected()

        @property
        def dg(self):
            return self.get_dg()

        @property
        def t_step(self):
            return self.get_t_step() or None

        @t_step.setter
        def t_step(self, value):
            self.set_t_step(float(value))

        @property
        def density(self):
            """The occupancy grid, ``(ng, ng, ng)`` uint8."""
            m = _IMP_bff.DyeDiffusionSimulation_get_density(self)
            if m.size == 0:
                return None
            return _as_cube(m).astype(np.uint8)

        @density.setter
        def density(self, grid):
            self.set_density(
                np.ascontiguousarray(np.asarray(grid, dtype=np.int32)).ravel())

        @property
        def slow_factor_map(self):
            m = _IMP_bff.DyeDiffusionSimulation_get_slow_factor_map(self)
            if m.size == 0:
                return None
            return _as_cube(m)

        @property
        def slow_density(self):
            m = _IMP_bff.DyeDiffusionSimulation_get_slow_density(self)
            if m.size == 0:
                return None
            return _as_cube(m).astype(np.uint8)

        @property
        def quenching_rate_map(self):
            m = _IMP_bff.DyeDiffusionSimulation_get_quenching_rate_map(self)
            if m.size == 0:
                return None
            return _as_cube(m)

        @property
        def mean_position(self):
            return _IMP_bff.DyeDiffusionSimulation_get_mean_position(self)

        @property
        def k_quench(self):
            """The quenching rate the dye sees, frame by frame, 1/ns.

            ``float32``, deliberately. It halves the memory the photon race
            walks, and that race is the bottleneck once a trajectory is long --
            it makes tens of millions of random reads into this array, so
            80 MB against 160 MB decides whether it fits in cache. The
            precision costs nothing real: a PET rate constant is a transferable
            starting value known to perhaps two significant figures, and
            float32 carries seven. The fused kernel holds its trace at the same
            width, which is what lets the fused and three-call paths be
            compared for *equality* rather than to a tolerance.

            The **rounding** is the kernel's, not this cast's -- `get_k_quench`
            already returns float-rounded values, so that a C++ caller and a
            Python one race against the same rates. This narrows the dtype and
            changes no value.
            """
            return _IMP_bff.DyeDiffusionSimulation_get_k_quench(self).astype(np.float32)

        @property
        def quenched(self):
            """Which frames the dye spent in contact with a quencher."""
            return self.k_quench > 0.0

        @property
        def collision_fraction(self):
            return self.get_collision_fraction()

        def sample_grid(self, grid):
            """Read a per-voxel field along the trajectory."""
            values = np.ascontiguousarray(np.asarray(grid, dtype=np.float64))
            ng = int(values.shape[0])
            out = _IMP_bff.DyeDiffusionSimulation_sample_grid(
                self, values.ravel(), ng)
            return out.astype(np.float32)

        def run(self, D=40.0, slow_fact=0.01, t_step=0.002, t_max=10000.0,
                n_trajectories=-1, random_seed=None):
            """Simulate the walk; returns the trajectory or ``None``."""
            n = _IMP_bff.DyeDiffusionSimulation_run(
                self, float(D), float(slow_fact), float(t_step), float(t_max),
                -1 if n_trajectories is None else int(n_trajectories),
                -1 if random_seed is None else int(random_seed))
            return None if n == 0 else self.trajectory
    %}
}
