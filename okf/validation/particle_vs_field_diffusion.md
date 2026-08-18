---
type: validation
title: "The particle walk's step width was 3x too large — a transcription error from a teaching wiki, fixed"
description: The Brownian walk in quenching/diffusion.py put 6*D*dt into the per-component step variance where the convention is 2*D*dt, so at the same nominal D it diffused three times faster than GridDiffusionSolver. Found during the numba-to-C++ port (PRD-113 stage 5). Traced to a 2019 QuEst docstring that transcribes a Berkeley teaching page's step *magnitude* as a per-component sigma. Fixed 2026-08-18; the D defaults were left alone, because D = 40 A^2/ns is free Alexa488 entered correctly and no calibration of D exists anywhere in the stack.
resource: /Users/tpeulen/dev/imp.bff
tags: [validation, imp.bff, quenching, diffusion, brownian, convention, units, resolved]
timestamp: '2026-08-18T00:00:00Z'
---

# The particle walk diffused at 3D — fixed, and the defaults were right all along

**Status: resolved 2026-08-18.** The step width is corrected. `D` is unchanged.

## What was measured

The step width in `_simulate_scalar` / `_simulate_grid`, now
`brownian_walk_in_volume`:

```python
sigma = np.sqrt(2 * D * 3 * t_step) / dg      # per Cartesian component
```

Per-component step variance `6 D dt`, where the Brownian convention is `2 D dt`
— `6 D dt` is the *total* three-dimensional MSD, used here as the width of one
component.

Measured 2026-08-18, `D = 8 Å²/ns`, `dt = 0.005 ns`, `dg = 1 Å`:

| | per-component ⟨Δx²⟩ | expected |
|---|---|---|
| numba (pre-port) | 13.98 | `2Dt` = 4.80 |
| C++ (faithful port) | 13.36 | `2Dt` = 4.80 |
| | | `6Dt` = 14.40 |

The field solver on the same nominal `D`, a point source on an 81³ grid:

```
grid solver: <x^2> = 40.0000    2Dt = 40.0000    ratio 1.0000
```

Exact. So the two models of the same dye, given the same `D`, diffused at `D`
and `3D`.

## Where the factor of three came from

QuEst's initial commit, `d02f41a` (2019-07-31),
`quenching/mfm/fluorescence/fps/_fps.pyx`, `simulate_traj_point` — the ancestor
of every version of this walk. Its docstring *is* the derivation:

```
dimensions = 2;         % two dimensional simulation
tau = .1;               % time interval in seconds
k = sqrt(D * dimensions * tau);

http://labs.physics.berkeley.edu/mediawiki/index.php/Simulating_Brownian_Motion
```

followed immediately by

```python
sigma = np.sqrt(2 * diffusion_coefficient * 3 * t_step) / dg
```

The cited page gives a 2-D recipe in which `k` is the **magnitude** of a step
whose *direction* is then drawn at random. It was transcribed as a
**per-Cartesian-component Gaussian σ**, with `dimensions` raised to 3 and an
extra factor of 2. That is the whole provenance: a units/convention misreading
of a teaching page, not a physical choice.

`git log -S "2 * D * 3" --all` in `quest` returns five commits, all of them pure
moves — Cython → numba (`dc3f88f`) → C++. The expression was never *changed*,
only relocated.

## Why `D` was left at 40

The first draft of this page argued the correction might be unsafe because
`D = 40 Å²/ns` looked "already high for a tethered dye — free Alexa488 in water
is around 4 Å²/ns". **That comparison was wrong by a factor of ten**:

```
1 Å²/ns = 1e-20 m² / 1e-9 s = 1e-11 m²/s = 10 µm²/s
free Alexa488 ≈ 400 µm²/s   = 40 Å²/ns
```

`D = 40` is exactly the free-solution value, entered in the parameter's
documented units. And this parameter has always meant free solution — QuEst's
`settings/parameter_catalog.json` defines it as *"Diffusion coefficient of the
donor dye **in solution**, when not interacting with the protein surface"*;
surface stickiness enters separately, through the mobility field.

The second reason to leave it: **no calibration of `D` exists anywhere in the
stack.** Searched `quest`, `chisurf`, `imp-tricks` and `imp.bff` — no fit
script, notebook, dataset, changelog entry or commit message ties any `D` to a
measurement. Every document that discusses these parameters says the opposite:
`quest/okf/references/pet-quenching-theory.md` calls them *"transferable
starting values, not constants"*, and PRD-110 scopes fitting `kQ` and `D` to a
measured decay as work **not yet done**. `40 / 3 = 13.3` matches nothing in any
repository.

So there was nothing for the factor of three to have been absorbed into. The
width was corrected and the defaults stand.

## What is now inconsistent, and was hidden before

Correcting σ exposes a real disagreement between the two models' *defaults*:

| | parameter | default |
|---|---|---|
| particle model (`simulate_dye_diffusion`) | `D` | **40** Å²/ns — free solution |
| field model (`quenching/dynamic.py`) | `free_diffusion` | **8.0** Å²/ns — tethered |
| QuEst's shipping project default | `"D"` | **7.5** Å²/ns |
| QuEst's per-dye repository (Alexa488) | | **7** Å²/ns |

Same name, same units, 5× apart. The `6 D dt` width used to mask this: it made
the particle model's effective mobility 120 against the field model's 8, so the
two were *already* 15× apart and nobody was comparing them.

This is a defaults question, not a physics one, and it is **not** decided here:
40 is a defensible free-solution value for a parameter documented as free
solution, and 7–8 is a defensible tethered value for a dye on a linker. What is
not defensible is calling both `free_diffusion` and shipping them 5× apart.
Recorded for whoever settles it.

## Where the checks live

`test/quenching/test_dynamics_cpp.py`:

* `test_free_walk_diffuses_at_the_stated_D` — per-component step variance is
  `2 D dt` and the free MSD is `2 D t`.
* `test_grid_solver_uses_the_standard_convention` — the same for the field
  model, so the two cannot drift apart again.
* `test_the_two_models_agree_on_D` — measures both and asserts they match.
