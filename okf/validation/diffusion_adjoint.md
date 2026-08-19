# The adjoint of the lattice field solver: exact to roundoff, 4.4× a forward

Recorded 2026-08-19 (PRD-115 stage 0, T-20260819-02). Kernel in
`src/DiffusionSolver.cpp` (`diffusion_propagate_adjoint`, `adjoint_sweep`),
Python `GridDiffusionSolver.gradient` in `pyext/src/sampling/smoluchowski.py`;
pinned by `test/quenching/test_diffusion_adjoint.py`.

## What it computes

For `L = Σ_k dL_dF[k] · F(t_k)` — any loss on the decay once its sensitivity to
each reported point is known — the gradient with respect to **every voxel** of
the mobility `D`, the rate `k` and the initial density, from one backward pass.
The forward sweep is linear in the density, `p_{n+1} = A p_n`, with `A` built
from `d = D dt/dg²`, `decay = e^{−k dt}` and `bounds`; the adjoint density obeys
`p̄_n = Aᵀ p̄_{n+1}` plus `dL_dF[k]` on every voxel at a reported step (the
observable is a plain sum), and the parameter gradients are local products of
forward state and adjoint at each step. Written as a gather over voxels — the
diagonal of `A` for the voxel itself, the off-diagonal from each active
neighbour that read it — so each output is owned by one iteration.

The forward keeps no history; the adjoint checkpoints every ⌈√n_steps⌉ steps
and re-runs each segment forward on the way back (two segments' worth of
memory instead of `n_steps · ng³` doubles, which is 2.6 GB at PRD-111 size).

## The check: dot-product identity against the forward

Forward and reverse share no code, so `⟨dL/dθ, v⟩` from the adjoint is compared
with a central difference of the *forward* solver along a random direction `v`,
for `θ ∈ {d, decay, cur}`:

| flux form | n_steps / n_out | domain | rel. agreement (d, decay, cur) |
|---|---|---|---|
| Smoluchowski | 47 / 5 (checkpoint length 7 does not divide) | sphere with a bite, off the shell | 3.5e-11, 2.8e-9, 1.2e-12 |
| Smoluchowski | 100 / 10 | same, touching the shell (kernel-level) | 1.9e-9, 2.1e-9, 1.3e-11 |
| Itô | 47 / 5 | off the shell | 1.4e-8, 3.1e-9, 1.3e-12 |
| Itô | 100 / 10 | touching the shell | 2.5e-8, 2.6e-8, 7.1e-13 |

(The `cur` direction is exactly linear, hence 1e-12; the others are limited by
the finite difference's own truncation at `h = 1e-5`.) The Python test repeats
this through `GridDiffusionSolver.gradient` — i.e. through the chain rule for
`dt/dg²`, `−dt·e^{−k dt}` and the `p b / Σ(p b)` normalisation `run()` applies —
at 1e-7, and reproduces the finite-difference Jacobian on the five PRD-111
parameters (`free_diffusion, slow_factor, contact_distance, kQ_scale, rC`)
through smooth field formulas to 1e-4 relative.

Two things the identity caught while writing it: a neighbour-interior test
that assumed the *centre* voxel was interior (wrong for a domain on the shell —
`0.5 %` disagreement, only in the shell-touching case), and nothing else. The
Smoluchowski stencil is self-adjoint up to the `bounds` mask; the Itô one is
not, and both are written out rather than derived from a symmetry.

## Cost

41³ grid, 14 000 active voxels, Python-level wall clock (loaded machine, ratio
is what matters):

| n_steps | forward | adjoint | ratio |
|---|---|---|---|
| 1000 | 0.09 s | 0.37 s | 4.3× |
| 4000 | 0.31 s | 1.38 s | 4.4× |

Where it goes: one forward with checkpoints (1×), one forward re-run inside the
segments (1×), and the reverse sweep at ~2.5× a forward sweep — it reads five
fields per neighbour (`p, d, bounds, decay, p̄`) where the forward reads three,
and it is memory-bound: a first version that tabulated the stencil
coefficients per voxel was *slower* (more memory), and templating the flux form
out of the inner loop gained ~10 %. So the ≤ 3× written into the ticket was a
guess and 4.4× is the floor of this formulation; a finite-difference gradient
costs one forward **per parameter**, so on a 14 000-voxel field the adjoint is
~3 000× cheaper, which is the number that matters.

Not done, on purpose: OpenMP is inert in this build
([`openmp_is_not_enabled.md`](openmp_is_not_enabled.md)), so no threading claim;
storing all forward states when memory allows (would save the 1× re-run) — a
`memory_budget` argument if a consumer ever needs it.
