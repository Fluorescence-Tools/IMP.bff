---
type: validation
title: "The diffusion stencil on the GPU: 5-12x, and the lead grows with the grid"
description: The explicit Smoluchowski sweep written in WGSL and run through wgpu on an M1 Pro against the same solve on the CPU. Against the eight-threaded CPU it is 5.3-5.6x at every grid size measured; against the single thread that shipped until today it is 5.9x at ng=41 and 12.3x at ng=81. The fluorescence trace agrees to 1.2e-6 (ng=41) and 1.1e-5 (ng=81). Machine at load average 30, so minima over repetitions; the ratios are consistent across sizes, which the absolute times are not.
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

| grid | CPU, 1 thread | CPU, 8 threads | GPU | vs 1 | vs 8 | deviation |
|---|---|---|---|---|---|---|
| 41³ | 448 ms | 440 ms | 76 ms | 5.9× | 5.6× | 1.2 × 10⁻⁶ |
| 61³ | 1 725 ms | 951 ms | 178 ms | 9.6× | 5.3× | 2.9 × 10⁻⁶ |
| 81³ | 4 641 ms | 1 992 ms | 372 ms | 12.3× | 5.4× | 1.1 × 10⁻⁵ |

**Read the last two columns together.** Against the eight threads the factor
is flat at about 5.4 — both sides scale with the grid, so the GPU is simply
the faster machine. Against one thread it *grows*, 5.9 → 12.3, because
OpenMP itself only starts paying at the larger grids (1.0× at 41³, 1.8× at
61³, 2.3× at 81³ — a 41-voxel x-slab is too little work to pay for a fork and
join four thousand times, and that is the same shape recorded in
`openmp_is_not_enabled.md`).

That growth is the point. The solver's cost goes as `dg⁻⁵`, so the resolutions
that matter scientifically are exactly the ones where the CPU is worst and the
GPU's lead is largest.

## What is trustworthy here and what is not

- **The deviations are trustworthy.** They match what
  `f32_holds_the_diffusion_stencil.md` predicted from numpy alone, which is
  two independent routes to the same answer: single precision does not cost
  this scheme anything a fitted decay curve could notice.
- **The ratios are reasonably trustworthy**: consistent across three grid
  sizes and across repetitions, taken as minima.
- **The absolute times are not.** The machine was at load average 30. A CPU
  time swung by 12× between runs of the same configuration; the GPU times
  moved by a few per cent, which is itself worth noticing -- the GPU path is
  insulated from whatever else the box is doing.
- **This is not the backend.** It is wgpu driven from Python; the C plugin has
  its own dispatch overhead, likely lower. And the adjoint is not measured at
  all.
