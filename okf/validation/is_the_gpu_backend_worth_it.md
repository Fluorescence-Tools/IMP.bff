---
type: validation
title: "A critical review of the WGSL backend: the Krylov solver on one CPU beats it"
description: Asked whether the GPU work is needed and whether it makes sense. Measured head to head, the Krylov solver added the same day is 1.28-1.53x faster than the WGSL stepping backend at every grid size, on the CPU, with no plugin and no wgpu. The GPU path keeps two things Krylov does not have -- the exact stepping semantics and the adjoint -- and costs ~700 lines of C against a moving Rust ABI, a vendored 6766-line header, and a Vulkan/D3D12 path that has never run. The network path has no caller at all. What is unambiguously worth keeping is the door, the Krylov solver, and the measurements.
resource: /Users/tpeulen/dev/imp.bff
tags: [validation, imp.bff, gpu, webgpu, review, prd-139]
timestamp: '2026-09-08T00:00:00Z'
---
# Is the GPU backend needed?

Written after building it, and it does not flatter what was built.

## The head-to-head that matters

Same fixture, same session, 4 000 steps reported every 100, eight CPU threads:

| grid | WGSL stepping (GPU) | stepping (CPU) | **Krylov trace (CPU)** | Krylov vs GPU |
|---|---|---|---|---|
| 41³ | 32 ms | 344 ms | **25 ms** | 1.28× faster |
| 61³ | 74 ms | 716 ms | **50 ms** | 1.49× faster |
| 81³ | 149 ms | 1 437 ms | **97 ms** | 1.53× faster |

**The Krylov solver on the CPU beats the GPU backend at every size**, and it
needs no plugin, no wgpu-native, no vendored header, no plugin ABI, no f16
subtlety and no GPU. Both landed the same day; had the order been reversed,
the WGSL diffusion backend would probably not have been written.

## What the GPU path still has that Krylov does not

Two things, and they are not nothing:

- **The exact stepping semantics.** Krylov solves the semi-discrete equation
  and is therefore about 2e-6 away from what `diffusion_propagate` returns --
  more accurate, and a different answer. A caller pinned to the stepping
  scheme's own numbers is served only by the stepping path.
- **The adjoint.** `diffusion_propagate_adjoint` is checkpointed reverse-mode
  over the sweeps. Krylov has no steps to accumulate parameter gradients over,
  and nothing has been worked out. Every gradient in the quenching stack goes
  through it.

So the GPU backend is not useless. It is, on the forward solve that
calibration actually hammers, **slower than a CPU alternative that costs
nothing to ship.**

## What it costs to keep

| | |
|---|---|
| the plugin | ~700 lines of C against wgpu-native's Rust ABI, which moves between major versions |
| the vendored header | 6 766 lines, kept in step with whatever wgpu-py ships |
| the contract | restated in C because `Compute.h` is C++, kept in sync by a text-matching test |
| the traps | three abort-class validation errors found so far (workgroup count, binding size, a default that is not zero); wgpu validation failures panic in a callback that cannot unwind, so they kill the process rather than declining |
| **platforms** | **it has only ever been compiled and run on macOS/Metal.** Vulkan and D3D12 have never executed a line of it |

**What it does *not* cost, corrected 2026-09-08.** An earlier version of this
note, and of PRD-139, said the plugin could not exist in the conda `imp.bff`
module package because IMP's tooling cannot build a second shared library.
That confuses building with shipping. The plugin is one C file, 67 KB
compiled, that links nothing but libc and is located **by path at run time**.
IMP's tooling never has to know about it: the recipe compiles it in one
command and drops it beside the module, and `enable_gpu()` finds it there
exactly as it does in the wheel -- and still only uses it if a wgpu library is
present. Both recipes now do this
(`conda-recipe/build.sh`, `conda-recipe/bld.bat`), verified by building into a
prefix and loading the result. So the capability is available to **both**
packages, and the only question left about the backend is whether it is worth
having at all.

That last row is the sharpest. `gpu/` was wired into the standalone build
unconditionally, and CI builds wheels with cibuildwheel on Linux and Windows,
so **CI would have been the first compiler ever to see it on two of three
platforms.** It now has an off switch (`IMPBFF_WITH_GPU`, default ON) so an
optional run-time plugin cannot break the build that does not want it. That is
a fix for something this work introduced, not a feature.

## The network path has no caller

`MlpCore.h` is vendored and **nothing in the library calls it**. PRD-115 is a
plan, not code. `IMP::bff::NeuralNet` and the WGSL forward pass were built for
a consumer that does not exist yet, and 12× on a network nobody evaluates is
worth zero today.

The `NeuralNet` class itself is defensible on its own terms -- the vendored
header says in as many words that the bindings should go through such a face,
and it is tested against an independent numpy forward. **The GPU half is
speculative generality**, and it also cost an ABI bump (1 → 2) on a published
plugin contract for a feature with no consumer. Harmless today because no
third-party plugin exists; still churn.

## What is worth keeping regardless of the above

- **The door.** `Compute.h` is small, and taking the whole loop rather than
  one step is the right granularity -- the same judgement saved the network
  path from a per-GEMM design that would have moved a hundred megabytes per
  layer.
- **The Krylov solver.** 14–24× on one thread, no dependency, more accurate
  than what it replaces.
- **The measurements**, which are durable whatever happens to the code: that a
  single-workgroup reduction was 16% of the run *and* two orders of magnitude
  of the error attributed to f32; that f16 weights are exact for one flux form
  and not the other, for a reason about the physics; that a GPU crossing costs
  a flat 1.2 ms of synchronisation on this machine; that residency, not
  kernels, is where the remaining wins are.

## Recommendation

1. **Do not extend the WGSL work.** Freeze the network path until PRD-115 has
   a caller.
2. **Decide the diffusion backend rather than letting it drift.** Either
   retire it in favour of Krylov and keep the door for the next candidate, or
   measure Krylov-on-GPU -- the only version that would beat both, and the one
   PRD-140 already flags as reduction-dominated and unmeasured.
3. **Before it ships anywhere but macOS**, build and run it on Linux. Until
   then it is a macOS feature with a cross-platform build.
4. The next candidate that needs none of this is `get_av` (2.8 ms at 1.5 Å,
   50.5 ms at 0.5 Å): one call, one answer, comfortably above the 1.2 ms
   synchronisation floor.
