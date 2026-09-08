---
type: validation
title: "The diffusion stencil on the GPU: 13x the eight-threaded CPU, and four of the five wins were about bytes"
description: The explicit Smoluchowski sweep written in WGSL and run through wgpu on an M1 Pro against the same solve on the CPU. The naive kernel is 4.3-5.7x the eight-threaded CPU; five optimisations take it to 12.5-13.1x, and against the single thread that shipped until today it is 14-35x. The two that paid best were both about bytes, not arithmetic: the neighbour index table is redundant (the grid numbering does not change), and f16 weights are exact if the self term is built from the rounded ones. The fluorescence trace agrees to 2.6e-6 (ng=41) and 1.1e-5 (ng=81).
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

**Four of the five are about bytes, not arithmetic**, which is the shape of
this kernel: seven or eight numbers are read per voxel and six multiply-adds
are done with them.

The last one is the interesting one. Writing the step as

    nxt = [decay·(1 − Σ a_q)]·p0 + Σ [decay·a_q]·cur[m],   a_q = ½(d0+d_m)·b_m

and computing the six `a_q` once turns **eighteen scattered reads per voxel
into six** — the mobilities and the mask are gone from the inner loop — and
the six weights that replace them are read in a coalesced layout. It is worth
1.6× on its own, and it is exact: the deviation against the f64 CPU is
unchanged, and against the naive f32 kernel it is 1e-7.

With all five, against the eight-threaded CPU measured in the same session
(378 / 788 / 1 550 ms; the single thread is 418 / 1 578 / 4 181 ms):

| grid | naive | + weights | + no table | + f16 | deviation |
|---|---|---|---|---|---|
| 41³ | 5.6× | 9.2× | 11.9× | **13.1×** | 2.6 × 10⁻⁶ |
| 61³ | 4.8× | 7.4× | 11.4× | **13.3×** | 3.0 × 10⁻⁶ |
| 81³ | 4.1× | 6.5× | 10.2× | **12.5×** | 1.1 × 10⁻⁵ |

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

## What has not been tried

In rough order of what looks most promising:

- **Tiling with workgroup memory.** The six neighbour reads of `cur` are still
  scattered over the full cube. A tile with a halo in workgroup storage is the
  standard answer for a 3-D stencil and the largest remaining structural win.
- **Renumbering the density into the compacted layout**, so that neighbours
  are near each other in memory rather than a plane apart.
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
