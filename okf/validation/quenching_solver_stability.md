# What constrains the explicit time step (and what does not)

Recorded 2026-08-18. Kernel and criterion in `pyext/src/quenching/solver.py`;
pinned by `test/quenching/test_quenching_field.py::RateStabilityTests`.

## Two thresholds, not one

The explicit update leaves a voxel with the coefficient

```
c = 1 − 6 D dt / dg²  −  k dt
```

and the most oscillatory mode is amplified by `|c|` per step. So:

| sum `6D dt/dg² + k dt` | what happens |
|---|---|
| ≤ 1 | `c ≥ 0`; the density stays non-negative |
| 1 – 2 | `c < 0`: bounded, but the density oscillates negative |
| **> 2** | **`|c| > 1`: diverges** |

`diffusion_stability_limit` returns the **positivity** bound, `dt ≤ 1/(6D/dg² +
k_max)` — the stricter of the two. A negative probability density is not an
acceptable answer either, so the factor of two is deliberately not spent. The
historical `dg²/(6D)` was this same bound for pure diffusion.

## The defect: the rate term was omitted

The criterion validated only the diffusion half, and on a quenched site the rate
half is far larger. Measured at T4L site 19, 2.5 Å, with a 25 Å²/ns ceiling on
`D` setting the step:

| | contribution |
|---|---|
| diffusion | 0.16 |
| **quenching** | **2.01** |

Sum 2.17 — past *divergence*, not merely past positivity — and the decay reached
**7 × 10³⁶** while the solver reported the step as safe.

**The divergence need not look like one.** Site 124, at a sum of 2.10, returned a
smooth, finite, monotone, entirely plausible decay: 8.841 × 10⁻⁵ at 25 ns against
9.067 × 10⁻⁵ correct. **2.6 % wrong, with nothing to see.** That is the reason
this is pinned by a test rather than left to inspection — an unstable run that
announces itself is a nuisance; one that does not is a wrong result.

## The fix: integrate the rate exactly

The kernels now apply `exp(-k dt)` as a factor instead of subtracting `k dt`.
That is the exact solution of `dp/dt = -k p` over the step, so the rate
contributes **no stability constraint at all** and only the diffusion CFL binds.
A strongly quenched site now costs the same as a weak one, where before it either
needed a smaller step or silently lied.

Pinned: a uniform rate reproduces `exp(-k t)` to 1e-10 *whatever* the step, and a
rate at `k dt = 1` — fatal when subtracted — stays finite, non-negative and
monotone.

## A false alarm, recorded because it was published

On finding this I claimed PRD-111 stages 0–1 were contaminated, with the contact
sites unstable above `kQ_scale` = 1.02–1.10 and the Fisher Jacobian's ±5 % step
reaching through a diverging integrator. **That was wrong twice.**

* I used the **positivity** bound as though it were the divergence bound, halving
  every threshold.
* I computed the thresholds with the **new ChiSurf harness's** step, 1.6× larger
  than the one stages 0–1 actually used, because that harness caps
  `free_diffusion` at 25 rather than 40 Å²/ns.

Stages 0–1 ran at sums of **0.92–0.98**, inside the positivity bound. Re-measuring
after the fix reproduces them to about 1 %:

| | pre-fix | post-fix |
|---|---|---|
| joint eigenvalues | 8.209e5, 2.994e4, 735, 128, 2.659 | 8.223e5, 2.992e4, 736, 129, 2.673 |
| `contact_distance` | ±61.1 % | ±60.9 % |
| rotation gain | 396× | 401× |

So the divergence was real but belongs to the **new harness**: tightening the `D`
ceiling to save compute is what raised the step and pushed the contact sites past
2. A performance change created the instability, and the first diagnosis
attributed it backwards to work that predated it.

## One hazard the fix does remove

At 2.0 Å the fits were free to search `kQ_scale` up to 10, where the sum reaches
**8.6** and the old scheme would have diverged outright. They converged near 1.0
and never went there — but nothing was stopping them, and a fit that wandered
would have returned a number rather than an error.
