---
type: validation
title: "The diffusion solver without time steps: a Krylov exponential, implemented, 14-24x on one thread"
description: The explicit solver's 4000 steps come from the forward-Euler stability bound, and the shipped scheme sits at 77% of it. A Krylov (Lanczos) exponential does the whole propagation in 64-128 operator applications instead, because the operator is constant over a run and symmetrisable. Now implemented in C++ as diffusion_propagate_krylov / diffusion_trace_krylov and measured: 24x on one thread for the trace, 14x for the safe default basis, 8x with the density as well, and it agrees with the exact matrix exponential to 1e-13 while the stepping route it replaces is 2e-6 away from it. RKL2 super-time-stepping was the alternative and is worth only 2.1x.
resource: /Users/tpeulen/dev/imp.bff
tags: [validation, imp.bff, numerics, diffusion, quenching, performance, krylov]
timestamp: '2026-09-08T00:00:00Z'
---
# The step count is not a law of nature

`diffusion_propagate` takes four thousand explicit Euler sweeps to cover 50 ns
because `dt <= dg^2/(6D)` says it must. That bound is why the solver costs
`dg^-5`, why calibration needs a surrogate, and why the GPU work in
`gpu_diffusion_is_worth_it.md` was worth doing at all. This note asks what can
be done about the bound rather than about the cost of a step.

Measured with `benchmark/large_time_steps.py`, at ng = 41, against
`scipy.sparse.linalg.expm_multiply` -- that is, against the exact solution of
the semi-discrete equation, not against the shipped scheme.

## The budget

**The shipped scheme's own error is 2.04e-6.** Forward Euler on the diffusion
plus a first-order Lie splitting against the decay: that is the accuracy that
ships today, and it is the number any faster method has to fit inside. It is
not a target to reproduce. A method that is *more* accurate is a behaviour
change even though it is an improvement, and pinned test values would move.

Two other facts, both measured:

| | |
|---|---|
| `rho(A+K)` | 1.530, against a forward-Euler limit of 2 -- the shipped step is at **77% of the limit**, so there is no slack to take by simply raising `dt` |
| stiffness `\|\|L\|\| T` | 6 121 |

## Why leapfrog is the wrong family, and what is the right one

Leapfrog is for time-reversible second-order dynamics. Applied to a parabolic
equation, centred-in-time-and-space is the textbook *unstable* scheme --
unstable for every step size, not merely for large ones. DuFort-Frankel is the
historical repair and it is unconditionally stable, but only conditionally
*consistent*: the step still has to stay small or the scheme solves a
different equation. It has been superseded.

What works on a parabolic problem divides into three families: stabilised
explicit methods (larger stability region), implicit methods (no bound, at the
price of a solve), and exponential methods (no time stepping at all). Two of
the three were measured here. Implicit methods were not, and the reason is in
the last section.

## Two properties of this operator that a general PDE does not have

Both are what make the third family available, and both are worth stating
because they are the reason the answer is unusually good:

* **The operator does not change over a run.** `d`, `decay` and `bounds` are
  fixed for the whole propagation, so it is not a sequence of different
  problems -- it is one matrix exponential, `p(t) = exp(-L t) p0`.
* **It is symmetrisable.** The Smoluchowski flux form
  `0.5(d_c + d_m)(p_c - p_m) b_m` is symmetric as it stands (measured
  asymmetry 0.0). The Ito form `(d_c p_c - d_m p_m) b_m` is not (5.5e-2), but
  `D^(1/2) A D^(-1/2)` is symmetric to 6e-17 and positive semidefinite -- the
  change of variables `q = sqrt(d) p` is a diagonal multiply. So the symmetric
  Krylov machinery serves **both** flux forms.

## Super-time-stepping (RKL2): a drop-in worth 2.1x

An explicit method whose stability polynomial is a shifted Legendre
polynomial: `s` stencil applications advance the solution by up to
`s(s+1)/4` explicit steps, because the real-axis stability interval grows as
`s^2` instead of staying at 2. Nothing else changes -- the same kernel, the
same buffers, a different outer loop with closed-form coefficients.

| s | steps per cycle | matvecs | vs 4 000 | error |
|---|---|---|---|---|
| 4 | 5.6 | 2 864 | 1.4× | 1.7e-7 |
| **6** | **12.4** | **1 938** | **2.1×** | **8.1e-7** |
| 8 | 21.7 | 1 480 | 2.7× | 2.7e-6 — over budget |
| 16 | 83.8 | 768 | 5.2× | 3.9e-5 — over budget |
| 32 | 327.2 | 416 | 9.6× | 5.9e-4 — over budget |

Stability is not what limits it -- every row above is stable. **Accuracy is.**
RKL2 is second order, so the error grows as the square of the superstep, and
past s = 6 it exceeds what ships today. 2.1× is the honest figure, for
perhaps a day of work and no new concepts.

## A Krylov exponential: 31-62x, and it gets better as the grid grows

Because `L` is constant and symmetric, the whole propagation is one Lanczos
run. Project onto an `m`-dimensional Krylov space, take the exponential of the
resulting `m x m` tridiagonal, and **every** output time comes out of the same
projection:

    F(t) = beta * u^T exp(-T_m t) e_1 = beta * sum_i c_i exp(-lam_i t)

| m | matvecs | vs 4 000 | error f64 | error f32 |
|---|---|---|---|---|
| 32 | 32 | 125× | 1.2e-10 | 3.8e-7 |
| **64** | **64** | **62×** | 2.3e-12 | **3.7e-7** |
| **128** | **128** | **31×** | 5.1e-14 | **3.5e-7** |
| 256 | 256 | 16× | 5.0e-14 | 3.4e-7 |

The f32 column plateaus at 3.4e-7 -- that is the recursion's single-precision
floor, and it is still six times better than what ships today.

**It survives the cases that were meant to break it.** The first result came
from a smooth start on a mild field, which is the easy case, so the sweep was
repeated over both flux forms, a point-source start, and a quencher raising
`k` by 300× in a four-voxel ball. The smallest `m` that stays inside the
2.04e-6 budget, in f32:

| flux form | start | quenching | m |
|---|---|---|---|
| Smoluchowski | smooth | mild | 64 |
| Smoluchowski | smooth | strong | 64 |
| Smoluchowski | point | mild | 128 |
| Smoluchowski | point | strong | 128 |
| Ito | smooth | mild | 64 |
| Ito | smooth | strong | 128 |
| Ito | point | mild | 128 |
| Ito | point | strong | 128 |

A point source is the hard case and it doubles `m`; nothing doubles it twice.

**And `m` scales as the square root of the stiffness**, not linearly. Doubling
the number of steps costs about 1.4× more Krylov vectors, where the explicit
scheme costs 2×. The advantage therefore *grows* exactly where the solver
hurts -- which is the same thing the GPU work found, for a different reason.

### Three practicalities, all measured

* **Reorthogonalisation is unnecessary.** It would cost `m^2 n` flops -- at
  m = 256, n = 92 k that is 6e9, about nine thousand stencil applications, more
  than the four thousand steps being replaced. Without it the results are
  identical to fourteen digits at every `m` tested. (Known: for `exp(-tL)v`
  the ghost copies belong to already-converged eigenvalues.)
* **The trace needs no stored basis.** The observable is a plain sum, so the
  only thing each Krylov vector has to leave behind is its column sum
  `u_j = 1^T v_j`, one scalar. Three vectors of memory for the whole run.
* **The final density does need the basis**, so it costs either `m` stored
  vectors (245 MB at ng = 81, m = 256, f32) or a second pass over the
  recursion with the coefficients already known -- `2m` matvecs, still 15-30×.

### Restarting makes it worse, which is worth knowing

The natural instinct is to chop the interval into pieces with a small `m`
each. Measured, it is a mistake: 40 restarts of m = 16 (640 matvecs) gave
2.7e-5 where one basis of m = 128 gave 4.7e-12. Lanczos converges
superlinearly once `m` passes a threshold set by the stiffness of the
*interval*, and restarting resets you below the threshold every time. One long
basis is both cheaper and far more accurate.

## What this would cost on the GPU, and the one open risk

Each matvec is exactly the stencil kernel that already exists, so the 13-14×
from `gpu_diffusion_is_worth_it.md` applies to each of the 128 rather than to
each of the 4 000.

**But the arithmetic of the run changes shape.** Lanczos needs `v^T w` and
`||w||` each step -- both reductions over the same vector, so one fused
reduction per matvec. At ng = 81 the two-stage reduction measures ~90 us
including its pass break, so 128 of them is ~11.5 ms against ~2.7 ms of
matvecs. **The method would turn a bandwidth-bound stencil loop into something
dominated by global reductions**, which is a different optimisation problem
from the one just solved. Rough estimate for ng = 81: ~14 ms against 98 ms
today, so ~7× on the GPU and ~100× against the threaded CPU -- but that is an
estimate, not a measurement, and it is the first thing to measure before
committing to the design.

Implicit methods were not measured for the same reason in reverse: Crank-Nicolson
needs a linear solve per step, and unpreconditioned CG converges in about
`sqrt(kappa)` iterations, which here is the same order as RKL2's stage count --
so it would have to beat RKL2's 2.1×, not the 4 000 steps. ADI would turn the
3-D solve into tridiagonal sweeps along lines, which suits a GPU badly and
suits the irregular domain worse. Neither is competitive with 31× once the
exponential route is on the table.

## It is implemented

`src/KrylovDiffusion.cpp`, reached as `diffusion_propagate_krylov()` (trace and
final density) and `diffusion_trace_krylov()` (trace only, half the work,
which is what calibration wants). Wall clock on the same fixture, ng = 41/61/81,
4 000 steps, reported every 100 -- one thread, which is what the library
actually shipped until the OpenMP fix:

| grid | stepping | trace, m=64 | trace, m=128 | trace + density, m=128 |
|---|---|---|---|---|
| 41³ | 410 ms | 17.8 ms (**23×**) | 30.9 ms (13×) | 50.9 ms (8.1×) |
| 61³ | 1 595 ms | 64.9 ms (**25×**) | 112.5 ms (14×) | 189.0 ms (8.4×) |
| 81³ | 4 139 ms | 170.7 ms (**24×**) | 292.5 ms (14×) | 490.1 ms (8.4×) |

and on eight threads, where the stepping route has more to gain from the
threads than this does:

| grid | stepping | trace, m=64 | trace, m=128 |
|---|---|---|---|
| 41³ | 335 ms | 28.1 ms (12×) | 47.0 ms (7.1×) |
| 61³ | 796 ms | 51.3 ms (16×) | 83.1 ms (9.6×) |
| 81³ | 1 680 ms | 99.2 ms (17×) | 158.5 ms (11×) |

Note the 41³ row on eight threads: the Krylov route is *slower* threaded than
serial (28.1 against 17.8 ms), the same fork-and-join tax the stepping route
pays at that size.

The measured factor is below the 31-62× the operator count promises, and the
reason is worth recording: **the vector operations, not the stencil.** Written
the obvious way -- serial loops over the full cube for each dot product, axpy
and norm -- they made eight threads worth nothing at all (4.3× at ng = 81
where one thread got 12×), because the generator scaled with the threads and
they did not. Walking a list of the active voxels instead, in parallel, and
fusing the two reductions into one pass took ng = 81 from 12.0× to 15.5× on
one thread and from 4.3× to 11.1× on eight. What is left of the gap is the
same thing: three passes over the active voxels per Lanczos step against the
generator's one pass with seven reads.

**Correctness.** Against `expm_multiply` the C++ agrees to 1.2e-13 at m = 128,
and against the stepping route it differs by exactly 2.04e-6 -- the number the
numpy study predicted for the stepping route's own error, arrived at through
entirely different code. `test/quenching/test_krylov_propagate.py` gates both
flux forms, the density as well as the trace, convergence in `m`, and the
refusal when Ito meets a zero mobility.

## What would change for callers

The Krylov route does not reproduce today's numbers -- it is more accurate
than they are, by about 2e-6. Anything pinning `diffusion_propagate` output
below that would move. It also produces the decay as **a sum of m
exponentials with explicit rates and amplitudes**, which is the form a
fluorescence decay is fitted with, so the natural interface is arguably not
"a trace sampled every 100 steps" at all.

## Sources

- Meyer, Balsara & Aslam, *A stabilized Runge-Kutta-Legendre method for
  explicit super-time-stepping of parabolic and mixed equations*,
  J. Comput. Phys. 257 (2014) 594.
- Vaidya et al., *Scalable explicit implementation of anisotropic diffusion
  with Runge-Kutta-Legendre super-time-stepping*, MNRAS 472 (2017) 3147.
- Niesen & Wright, `phipm` -- adaptive Krylov `exp`/`phi` evaluation, and
  Expokit before it, for the adaptive-`m` variant not implemented here.
