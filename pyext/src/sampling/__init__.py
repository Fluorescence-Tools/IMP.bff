"""Stage 3 — how configurations are drawn, enumerated, or propagated.

A representation says *where* the dye can be. Sampling says how you get at it,
and the method follows from the representation:

======================  ==============================================
representation          how it is sampled
======================  ==============================================
accessible volume       pathfinding on the grid, then a walk in it
rotamer library         screening the listed conformers
coarse-grained dye      molecular dynamics under the force field
======================  ==============================================

* :mod:`~IMP.bff.sampling.brownian` — one rejection-sampled trajectory at a
  time, which resolves the dye's *history*.
* :mod:`~IMP.bff.sampling.smoluchowski` — the same dynamics as a density,
  ``dp/dt = div(D grad p) - k p``. A Langevin trajectory and its Fokker-Planck
  density are **one model** sampled or integrated, not two rival ones; they must
  agree on ``D``, and ``test/dynamics`` asserts they do.
* :mod:`~IMP.bff.sampling.excited_state` — kinetic Monte Carlo over an
  excitation's fate, which is a different variable being sampled: not where the
  dye is, but whether the photon escapes.
* :mod:`~IMP.bff.sampling.rotamer_library` — loading a library and drawing from
  it by weight.

This package was ``dynamics/`` until the layout was reconciled with the four
stages. "Dynamics" named the integrators and left library screening and
pathfinding without a home; sampling covers all three representations, which is
the point of having a stage rather than a package per technique.

The explicit dye's MD samplers are still in ``cgdye/sampling`` — that package is
off the domain layout, and moving its samplers here without the model they
sample would split it in half.
"""

from IMP.bff.sampling.brownian import (  # noqa: F401
    DyeDiffusionTrajectory,
    simulate_dye_diffusion,
)
from IMP.bff.sampling.excited_state import (  # noqa: F401
    simulate_photon_trace,
    simulate_quenched_decay,
)
from IMP.bff.sampling.rotamer_library import (  # noqa: F401
    apply_rotamer_coordinates,
    load_rotamer_library_dcd,
    sample_rotamer_index,
)
from IMP.bff.sampling.smoluchowski import (  # noqa: F401
    GridDiffusionGradient,
    GridDiffusionResult,
    GridDiffusionSolver,
    diffusion_stability_limit,
    equilibrium_occupancy,
)

__all__ = [
    "DyeDiffusionTrajectory", "simulate_dye_diffusion",
    "simulate_photon_trace", "simulate_quenched_decay",
    "apply_rotamer_coordinates", "load_rotamer_library_dcd",
    "sample_rotamer_index",
    "GridDiffusionGradient", "GridDiffusionResult", "GridDiffusionSolver",
    "diffusion_stability_limit", "equilibrium_occupancy",
]
