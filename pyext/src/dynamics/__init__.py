"""How a model is realised in time: the integrators.

A representation says *where* the dye can be and with what weight. Dynamics says
how it gets there, and it is a separate choice: the same accessible volume can
be sampled by a Brownian walk, propagated as a density, or held still. Keeping
the two apart is what lets a new integrator land without touching the model --
the split MD engines have made for decades, and the reason a propagator does not
belong inside the thing it propagates.

Two kinds of dynamics live here, and they are about different variables:

**The dye's position.**

* :mod:`IMP.bff.dynamics.brownian` -- one trajectory at a time, rejection-sampled
  inside an occupancy grid. Resolves the dye's *history*, which is what a
  correlation function needs.
* :mod:`IMP.bff.dynamics.smoluchowski` -- the same dynamics as a density,
  ``dp/dt = div(D grad p) - k p`` on the grid. Deterministic and free of shot
  noise, at the cost of a stability limit on the time step, and it gives the
  equilibrium occupancy in closed form.

  These are a Langevin trajectory and its Fokker-Planck density: **one model,
  sampled or integrated**, not two rival ones. They must agree on ``D``, and
  ``test/quenching/test_dynamics_cpp.py`` asserts they do -- they did not until
  2026-08-18, when the walk's step width was found to be a factor of three too
  large (``okf/validation/particle_vs_field_diffusion.md``).

**The excited state.**

* :mod:`IMP.bff.dynamics.excited_state` -- kinetic Monte Carlo over an
  excitation's fate: race the intrinsic decay against quenching frame by frame,
  photon by photon, or propagate the population for the noise-free curve.

Moved out of ``quenching/`` by PRD-113 stage 7. None of this is specific to PET
quenching -- a Brownian walk in a volume is a Brownian walk in a volume -- and
filing it under the first physics that used it is how a general integrator comes
to look like a detail of one model.

Not here yet: the Langevin/MD sampler for explicit dyes, which is still
:func:`IMP.bff.make_langevin_simulator` in ``cgdye``. It belongs here.
"""

from IMP.bff.dynamics.brownian import (  # noqa: F401
    DyeDiffusionTrajectory,
    simulate_dye_diffusion,
)
from IMP.bff.dynamics.excited_state import (  # noqa: F401
    simulate_photon_trace,
    simulate_quenched_decay,
)
from IMP.bff.dynamics.smoluchowski import (  # noqa: F401
    GridDiffusionResult,
    GridDiffusionSolver,
    diffusion_stability_limit,
    equilibrium_occupancy,
)

__all__ = [
    "DyeDiffusionTrajectory", "simulate_dye_diffusion",
    "simulate_photon_trace", "simulate_quenched_decay",
    "GridDiffusionResult", "GridDiffusionSolver",
    "diffusion_stability_limit", "equilibrium_occupancy",
]
