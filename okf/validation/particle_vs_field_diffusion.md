---
type: validation
title: "The particle walk and the field solver disagree on D by exactly 3x"
description: The Brownian walk in quenching/diffusion.py puts 6*D*dt into the per-component step variance where the standard convention is 2*D*dt, so at the same nominal D it diffuses three times faster than GridDiffusionSolver -- which reproduces <x^2> = 2Dt exactly. Found during the numba-to-C++ port (PRD-113 stage 5), which preserved the convention rather than changing it. Open for the owner.
resource: /Users/tpeulen/dev/imp.bff
tags: [validation, imp.bff, quenching, diffusion, brownian, convention, open]
timestamp: '2026-08-18T00:00:00Z'
---

# The particle walk and the field solver disagree on D by exactly 3×

**Status: open.** The C++ port preserved the existing behaviour. Nothing has
been changed; this page records what was measured and what the decision is.

## What was measured

The step width in `_simulate_scalar` / `_simulate_grid`, now
`brownian_walk_in_volume`:

```python
sigma = np.sqrt(2 * D * 3 * t_step) / dg      # per Cartesian component
```

so the per-component step variance is `6 D dt`. The Brownian convention is
`2 D dt` per component — `6 D dt` is the *total* three-dimensional mean squared
displacement, used here as if it were the width of one component.

Measured 2026-08-18, `D = 8 Å²/ns`, `dt = 0.005 ns`, `dg = 1 Å`, 600 walks of
60 steps in a 101³ box (large enough that no walk touched a wall):

| | per-component ⟨Δx²⟩ | expected |
|---|---|---|
| numba (pre-port) | 13.98 | `2Dt` = 4.80 |
| C++ (post-port) | 13.36 | `2Dt` = 4.80 |
| | | `6Dt` = 14.40 |

Both give `6Dt`. **The port is faithful**; the convention predates it.

The field solver, on the same nominal `D`, a point source on an 81³ grid after
300 steps:

```
grid solver: <x^2> = 40.0000    2Dt = 40.0000    ratio 1.0000
```

Exact. So the two models of the same dye, given the same `D`, diffuse at
`D` and `3D` respectively.

## Why it matters

PRD-109 pinned `GridDiffusionSolver` to `⟨x²⟩ = 2Dt` to 0.1 % and treats the
two as the same dynamics sampled two ways — a Langevin trajectory and its
Fokker-Planck density. They are not, at the same parameter value. Anything that
compares them, or that carries a `D` from one to the other, is off by three.

## Why it was not simply fixed

Two reasons, neither of which is "it might be right":

1. **`simulate_dye_diffusion` defaults to `D = 40 Å²/ns`**, which is already
   high for a tethered dye — free Alexa488 in water is around 4 Å²/ns. If that
   default was chosen by matching simulated to measured decays *through this
   code*, it absorbed the factor of three, and correcting σ without revisiting
   the default would break agreement that currently holds.
2. **Every quenching number computed through the particle path moves**, by a
   factor of three in the mobility. PRD-113's rule for a port is that nothing
   moves; a change of this size is its own piece of work with its own
   re-validation.

## The decision

Three options:

* **Correct σ to `sqrt(2 D dt)` and re-derive the default `D`.** Makes the two
  models agree on what `D` means, at the cost of re-validating the particle
  path against whatever it was last calibrated on.
* **Keep σ and rename the parameter** — if the quantity really is a
  three-dimensional MSD rate, call it that, and give the walk a `D` that means
  the same thing the solver's does.
* **Keep both and document the ratio.** The weakest option: it leaves a factor
  of three between two functions in the same module.

Recommendation: the first, but it needs whatever data the `D = 40` default came
from, which is not in this repository.

## Where the check lives

`test/quenching/test_dynamics_cpp.py`:

* `test_free_walk_diffuses_at_three_times_D` — asserts the discrepancy, so it
  cannot be changed silently.
* `test_grid_solver_uses_the_standard_convention` — asserts the other half.

Both must be updated together if the convention is settled.
