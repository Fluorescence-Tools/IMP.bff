---
type: validation
title: "The diffusion stencil holds in f32: 1e-6 after four thousand steps"
description: WGSL has no f64, so a WebGPU backend for diffusion_propagate would run the explicit Smoluchowski sweep in single precision. Measured against the f64 kernel on the grids the tests use, the fluorescence trace deviates by 1.2e-6 relative at ng=41 over 4000 steps and the final density by 1.5e-6 -- six significant digits, far below the Poisson noise the curves are fitted against. The forward solve is therefore not blocked by precision. The adjoint is not covered by this and needs its own measurement.
resource: /Users/tpeulen/dev/imp.bff
tags: [validation, imp.bff, performance, gpu, webgpu, precision, diffusion]
timestamp: '2026-09-08T00:00:00Z'
---
# f32 is enough for the forward solve

**Status: measured, and it decides a design.** A GPU backend for
`diffusion_propagate` was proposed on WebGPU, and WGSL has **no f64** — only
f32, with f16 behind an extension. The solver is `std::vector<double>`
throughout. So before writing any wgpu plumbing, the question had to be
answered: does the scheme survive single precision?

## What was measured

The same iteration as `internal::lattice_sweep` — explicit Smoluchowski flux
over the six face neighbours, masked by `bounds`, multiplied by the per-voxel
`decay` — run in `float64` and in `float32` from identical initial
conditions, with a position-dependent mobility and decay rather than uniform
fields, and compared on the two quantities a caller reads back: the
fluorescence trace and the final density.

| grid | steps | fluorescence, max relative | final density, max relative | surviving population |
|---|---|---|---|---|
| ng = 21 | 2 000 | 5.2 × 10⁻⁶ | 5.5 × 10⁻⁶ | 0.689 |
| ng = 41 | 4 000 | 1.2 × 10⁻⁶ | 1.5 × 10⁻⁶ | 0.423 |

## Why it does not accumulate

The scheme is contractive: diffusion averages neighbouring voxels and the
decay multiplies by a number below one, so a rounding error introduced at
step *k* is damped by every later step rather than amplified. That is why
four thousand steps are *better* than two thousand in the table above rather
than worse — the longer run has spread and decayed more of its own noise.

## What this does not say

- **The adjoint is not covered.** `lattice_propagate_adjoint` runs the
  recursion backwards and gradients can amplify what the forward pass damps.
  Any GPU adjoint needs this same comparison before it is trusted, and it may
  well come out differently.
- **A GPU is not numpy.** These numbers come from IEEE binary32 in numpy. A
  GPU may fuse multiply-adds and may flush denormals, both of which move the
  last bits. The comparison to run before shipping is the real backend
  against the real CPU kernel; this measurement says the *scheme* is not the
  obstacle, which is what was in doubt.
- Nothing here argues that a GPU is faster. That is a separate measurement,
  and on a machine at load average 18.9 it cannot honestly be taken at all.
