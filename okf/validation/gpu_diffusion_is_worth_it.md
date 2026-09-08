---
type: validation
title: "The diffusion stencil on the GPU: 13-14x the eight-threaded CPU, and the last 16% was a reduction running on one core"
description: The explicit Smoluchowski sweep written in WGSL and run through wgpu on an M1 Pro against the same solve on the CPU. The naive kernel is 3.8-5.3x the eight-threaded CPU; six optimisations take it to 13.3x/13.2x/14.0x at ng=41/61/81. Most of them were about bytes rather than arithmetic: the neighbour index table is redundant, f16 weights are exact if the self term is built from the rounded ones, and the population sum was running in a single workgroup for 16% of the whole run. Fixing that also took the deviation from 1.1e-5 to 1.3e-7, so most of what looked like f32 error was a serial f32 accumulator. Remaining time: at ng=41 41% is the dispatch floor, at ng=81 85% is the kernel.
resource: /Users/tpeulen/dev/imp.bff
tags: [validation, imp.bff, performance, gpu, webgpu, diffusion, openmp]
timestamp: '2026-09-08T00:00:00Z'
---
# The GPU is worth it, and most worth it where it hurts

`benchmark/gpu_diffusion_wgsl.py` runs the same explicit Smoluchowski
propagation two ways: `IMP::bff::diffusion_propagate` on the CPU in f64, and
the identical stencil in WGSL on the GPU in f32. 4 000 steps, a report every
100, a position-dependent mobility and decay. Apple M1 Pro, Metal backend,
wgpu-native 29.0.1.1 through wgpu-py.

The first kernel written — a thread per voxel, a branch for what is outside:

| grid | CPU, 1 thread | CPU, 8 threads | GPU, naive | vs 1 | vs 8 | deviation |
|---|---|---|---|---|---|---|
| 41³ | 418 ms | 378 ms | 66.5 ms | 6.3× | 5.7× | 1.2 × 10⁻⁶ |
| 61³ | 1 578 ms | 788 ms | 169.9 ms | 9.3× | 4.6× | 2.9 × 10⁻⁶ |
| 81³ | 4 181 ms | 1 550 ms | 364.2 ms | 11.5× | 4.3× | 1.1 × 10⁻⁵ |

**Read the last two columns together.** Against one thread the factor *grows*,
6.3 → 11.5, because OpenMP itself only starts paying at the larger grids
(1.1× at 41³, 2.0× at 61³, 2.7× at 81³ — a 41-voxel x-slab is too little work
to pay for a fork and join four thousand times, and that is the same shape
recorded in `openmp_is_not_enabled.md`). Against the eight threads it
*shrinks*, 5.7 → 4.3: this kernel scaled worse than the CPU it was beating,
which is the signal that it was not compute-bound at all. The next section is
what that turned out to mean.

## Five optimisations, measured one at a time

The first kernel was the obvious one: a thread per voxel of the cube, a
branch to skip what is outside. Five changes, each measured against the one
before:

| | worth | why |
|---|---|---|
| dispatch only the active voxels | 1.3–1.5× | a third to two fifths of the cube is inside the volume, and the rest was being dispatched and then discarded by a branch |
| workgroup 256 instead of 64 | ~1.04× | — |
| precompute the stencil weights | ~1.6× | `d`, `decay` and `bounds` do not change over a propagation, so neither do the coefficients |
| drop the neighbour index table | 1.34–1.59× | the six neighbours are `c ± 1`, `c ± ng`, `c ± ng²`; the table was redundant |
| f16 weights | 1.09–1.23× | with the table gone the weights are most of what is read |
| two-stage population sum | 1.16–1.19× | the report ran in one workgroup over every active voxel |

**Most of them are about bytes and parallelism, not arithmetic**, which is
the shape of this kernel: seven or eight numbers are read per voxel and six
multiply-adds are done with them.

The last one is the interesting one. Writing the step as

    nxt = [decay·(1 − Σ a_q)]·p0 + Σ [decay·a_q]·cur[m],   a_q = ½(d0+d_m)·b_m

and computing the six `a_q` once turns **eighteen scattered reads per voxel
into six** — the mobilities and the mask are gone from the inner loop — and
the six weights that replace them are read in a coalesced layout. It is worth
1.6× on its own, and it is exact: the deviation against the f64 CPU is
unchanged, and against the naive f32 kernel it is 1e-7.

With all six, against the eight-threaded CPU measured in the same session
(351 / 692 / 1 378 ms; the single thread is 418 / 1 578 / 4 181 ms):

| grid | naive | + weights | + no table | + f16 | deviation |
|---|---|---|---|---|---|
| 41³ | 5.3× | 8.8× | 12.2× | **13.3×** | 2.6 × 10⁻⁶ |
| 61³ | 4.1× | 7.0× | 11.1× | **13.2×** | 6.1 × 10⁻⁷ |
| 81³ | 3.8× | 6.7× | 10.9× | **14.0×** | 2.9 × 10⁻⁷ |

Read the first column against the last. The naive kernel *loses* ground as
the grid grows, 5.6 → 4.1; the tuned one holds at about 13. A kernel that
reads less scales better, and `dg⁻⁵` means the large grids are the ones that
matter.

## The neighbour table was twenty-four bytes of nothing

Compacting the dispatch onto the active voxels brought an index list with it,
and with the index list came a table of each voxel's six neighbours — six
`u32`, the largest single read in the kernel. It should never have been
written. Compaction changes *which* voxels are visited; it does not change how
the grid is numbered, so the neighbours of `c` are still `c ± 1`, `c ± ng` and
`c ± ng²`. Three additions replace twenty-four bytes, and the result is
bit-identical: the deviation does not move at any grid size.

## f16 for the weights: worth 1.2×, and free if the row sum is protected

Halving the weights is worth 1.09× at ng = 41 and 1.23× at ng = 81 — the
larger the grid, the more of the time is spent reading them. Done the obvious
way it also costs three orders of magnitude of accuracy, 7 × 10⁻⁷ → 1 × 10⁻³,
which is where this stopped looking worthwhile.

The reason is not the weights themselves; f16 carries about 5 × 10⁻⁴ of
relative precision and the trace comes out at 1 × 10⁻³. It is that **the
quantity the propagation conserves is the row sum**: `self + Σw = decay`, one
identity per voxel. Round the seven numbers independently and the sum is wrong
by 5 × 10⁻⁴, every one of four thousand steps applies the same wrong sum, and
the error accumulates instead of averaging out.

So do not round them independently. Round the six weights to f16, then build
the self term **from the rounded weights** — `self = decay − Σ f16(w)` — and
keep that one number in f32:

```python
wh = w.astype(np.float16)
self_ = (k0 - wh.astype(np.float64).sum(0)).astype(np.float32)
```

The row sum is then exact to f32, and the trace is back at the f32 deviation:
2.6 × 10⁻⁶ / 3.0 × 10⁻⁶ / 1.1 × 10⁻⁵ against 7.2 × 10⁻⁷ / 3.0 × 10⁻⁶ /
1.1 × 10⁻⁵. Six of the seven numbers are still half width, so nearly all the
bandwidth is still saved.

Two variants that do **not** work, both measured:

| | deviation at 41³ |
|---|---|
| everything f16, no compensation | 1.0 × 10⁻³ |
| everything f16, compensated | 1.8 × 10⁻³ |
| f16 weights, f32 self, no compensation | 5.5 × 10⁻⁵ |
| **f16 weights, f32 self, compensated** | **1.0 × 10⁻⁶** |

The second row is the instructive one: compensating and then rounding the
self term is *worse than not compensating at all*, because compensation
gathers the whole row error into precisely the number that is then thrown
away. The compensation and the f32 self term are one measure, not two.

## The population sum was running on one core

The report -- the fluorescence value written every hundred steps -- summed the
density over every active voxel in **one workgroup**. At ng = 81 that is 897
serial iterations in each of 256 threads while the other thirty-one cores of
the GPU sit idle. Forty such reports cost 23 ms of a 119 ms run: **16% of the
whole propagation for forty numbers.**

Split into NPART partial sums and a second tiny pass over those, it costs
3.7 ms. NPART anywhere from 16 to 256 measures the same (1.6–4.4 ms, which is
the noise on this machine); what mattered was being more than one.

It also bought two orders of magnitude of accuracy, which nobody was looking
for. The deviation at ng = 81 fell from 1.1 × 10⁻⁵ to 1.3 × 10⁻⁷, and at
ng = 61 from 3.0 × 10⁻⁶ to 1.6 × 10⁻⁷. The old kernel was accumulating a
quarter of a million f32 values in series per thread; the new one accumulates
fourteen and then folds in a tree. **Most of what had been recorded as the
cost of f32 was the reduction, not the stencil** -- which strengthens
`f32_holds_the_diffusion_stencil.md` rather than weakening it, and is a
reminder that a parity number measures the whole path, not the part being
questioned.

## Where the remaining time goes

Measured with a no-op sweep (same dispatch count, empty body) and with the
reports switched off:

| | ng = 41 | ng = 81 |
|---|---|---|
| total | 26.5 ms | 98.3 ms |
| dispatch floor | 10.8 ms (41%) | 11.5 ms (12%) |
| reports | 2.9 ms | 3.7 ms |
| kernel | ~13 ms | ~83 ms |

The floor is 2.6–3.0 µs per dispatch and barely depends on the grid, so the
two ends of the range want opposite things. **At ng = 41 the kernel is already
cheaper than the command stream that launches it**, and no amount of kernel
work will show; only fewer dispatches can. At ng = 81 the kernel is 85% of the
run, and roughly half of that is the six neighbour gathers.

Of the floor, `set_bind_group` before each dispatch is 2.8 ms of 24.8 at
ng = 41 and nothing measurable at ng = 81; the rest is the dispatch itself.

## The backend exists now, and f16 turned out to be flux-form-specific

`gpu/imp_bff_wgpu.c` is the plugin the door in `Compute.h` was built for: a C
library that links nothing, resolves every wgpu entry point with dlopen from a
path the Python side finds, and installs itself beside the extension module.
`IMP.bff.get_compute_backend_name()` answers `wgpu:Metal:Apple M1 Pro:f16`.
Through the real API, against the eight-threaded CPU:

| grid | Smoluchowski | Ito |
|---|---|---|
| 41³ | 10.2× | 10.0× |
| 61³ | 9.2× | 7.8× |
| 81³ | 9.2× | 7.6× |

Lower than the 13× the harness measured, because the plugin does per call what
the harness did once: build the active-voxel list and the weights, create and
upload the buffers, and read back. The shader and the pipelines are cached
across calls, which is what matters for calibration's 100-1 000 solves.

**f16 weights are right for one flux form and wrong for the other**, and it
took a measurement to see why. The compensation that makes f16 exact keeps the
*row sum* exact, so the step is `(self + sum w) * p0` plus a term in the
differences between neighbouring voxels: the rounding cancels to first order
exactly when neighbours hold nearly the same density.

| flux form | equilibrium | neighbour-to-neighbour variation of the density | deviation with f16 |
|---|---|---|---|
| Smoluchowski | uniform, whatever D | 0.06% (median) | 4.5e-7 after 32 steps |
| Ito | `p ~ 1/D` | **12%** (median) | 2.3e-4 after 32 steps |

The Ito density inherits the mobility's own voxel-to-voxel variation — 12%
against the mobility's 11.9% on this field — so there is nothing for the
compensation to cancel against and the raw f16 error survives. It is a
per-step error, not an accumulating one: 32 steps and 4 000 steps are within a
factor of seven of each other, where an accumulating error would be 125 times
apart. The plugin therefore uses f16 for Smoluchowski and f32 for Ito, which
costs Ito about 1.2× and buys back three orders of magnitude
(1.2e-7 trace, 5.0e-7 density). `IMP_BFF_GPU_F16=off` turns it off everywhere.

## What has not been tried

In rough order of what looks most promising, with the headroom each one is
actually competing for:

- **Tiling with workgroup memory** — competing for the ~40 ms of neighbour
  gathers at ng = 81, so at best about 1.4× there and nothing at ng = 41. The
  six neighbour reads of `cur` are still scattered over the full cube; a tile
  with a halo in workgroup storage is the standard answer for a 3-D stencil.
  It fights the compaction, though: tiles are dense and the compaction is what
  made the kernel fast, so it needs tile-level compaction to keep both.
- **Temporal blocking — two or more steps per dispatch**, with a halo wide
  enough to cover them. The only idea that attacks *both* remaining costs: it
  divides the dispatch floor and the global traffic by the number of fused
  steps. It is also the hardest, and it is the one that helps at ng = 41,
  where 41% of the time is command overhead.
- **Renumbering the density into the compacted layout**, so that neighbours
  are near each other in memory rather than a plane apart. Measured headroom
  is smaller than it looks: the base access alone (index, self weight, own
  density, write) is ~45 ms of the 83 ms kernel at ng = 81, and within a
  z-run the active voxels are already contiguous.
- **On the CPU side**, hoisting the OpenMP parallel region out of the step
  loop: that is what makes threading worth nothing at ng = 41 today.
- **Algorithmically**, the explicit scheme's stability bound `dt ≤ dg²/(6D)`
  is what forces four thousand steps. An exponential or implicit integrator
  changes the *number* of steps rather than their cost, which is where orders
  of magnitude live rather than factors — and it is a numerics project, not an
  optimisation.

## What is trustworthy here and what is not

- **The deviations are trustworthy.** They match what
  `f32_holds_the_diffusion_stencil.md` predicted from numpy alone, which is
  two independent routes to the same answer: single precision does not cost
  this scheme anything a fitted decay curve could notice.
- **The ratios are reasonably trustworthy**: consistent across three grid
  sizes and across repetitions, taken as minima.
- **The absolute times are not.** The machine was at load average 15–30
  throughout. A CPU time swung by 12× between runs of the same configuration;
  the GPU times moved by a few per cent, which is itself worth noticing -- the
  GPU path is insulated from whatever else the box is doing. The CPU rows of
  the table above were therefore taken in the same session as the GPU rows,
  with the OpenMP build on `PYTHONPATH`; an earlier table in this note
  compared against a CPU measured on a busier machine and read 9.6× where
  this one reads 9.2× for the same kernel.
- **Which CPU is in the denominator has to be checked, not assumed.** The
  wheel in the test venv is built without OpenMP, so running the benchmark
  against it silently compares the GPU to *one* thread and inflates every
  ratio by two to three. `main()` now prints `parallel_threads()` for that
  reason.
- **This is not the backend.** It is wgpu driven from Python; the C plugin has
  its own dispatch overhead, likely lower. And the adjoint is not measured at
  all.
