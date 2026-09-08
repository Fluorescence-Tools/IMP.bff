---
type: validation
title: "torch on the same propagation: uncompetitive eager, competitive compiled, and it wins on a different axis"
description: The diffusion propagation implemented in PyTorch on the same fixture as the C++ core and the WGSL prototype. Eager torch is 6-10x slower than the C++ on the CPU and its MPS backend is slower than eight CPU threads, because 4000 steps of seven elementwise ops is 28000 kernel launches. torch.compile is worth 5-13x and makes it competitive. On the Krylov route torch is nearly level with the C++, because 128 matvecs amortise framework overhead 30x less often -- which is the interesting result. Three measurement hazards recorded, one of which silently returns a trace 2x too small.
resource: /Users/tpeulen/dev/imp.bff
tags: [validation, imp.bff, performance, torch, gpu, diffusion, krylov]
timestamp: '2026-09-08T00:00:00Z'
---
# What torch would be competing with

torch is the obvious thing to reach for: it is already installed wherever
`quench_pinn` is worked on, it has MPS here and CUDA on the compute box, and it
needs no plugin, no shader and no C ABI. Before writing a wgpu backend it is
worth knowing what that backend would be beating.

`benchmark/torch_diffusion.py` implements both routes on the same fields as
`benchmark/gpu_diffusion_wgsl.py`: the stepping route as seven shifted
multiply-adds with the stencil weights precomputed -- the same arithmetic the
tuned WGSL kernel does -- and the Krylov route the C++ core now implements.
4 000 steps, reported every 100, Apple M1 Pro, all in one session.

## The stepping route: torch is not close, until it is compiled

| grid | C++ 1 thr | C++ 8 thr | WGSL | torch cpu f32 | torch cpu f32 compiled | torch MPS f32 | torch MPS compiled |
|---|---|---|---|---|---|---|---|
| 41³ | 410 ms | 356 ms | **25 ms** | 2 411 ms | 487 ms | 1 008 ms | 170 ms |
| 61³ | 1 594 ms | 701 ms | **51 ms** | 4 553 ms | 754 ms | 2 235 ms | 247 ms |
| 81³ | 4 133 ms | 1 315 ms | **99 ms** | 10 544 ms | 1 322 ms | 4 548 ms | 345 ms |

**Eager torch is not in the running.** It is 6-8× slower than the C++ on the
same CPU, and its *GPU* backend is 2.8-3.5× slower than eight CPU threads of
the C++. The reason is not the arithmetic, which is seven multiply-adds per
voxel: it is that eager torch launches each of the seven operations separately
with a temporary in between, four thousand times over -- of the order of
28 000 kernel launches for a job whose whole content is one fused loop.

**`torch.compile` is worth 5.3-13.2×** and changes the verdict. Fused, torch
on MPS lands at 170-345 ms: faster than the C++ on eight CPU threads, and
3.5-6.7× behind the WGSL prototype. That is a real gap but a defensible one,
and it is the number the wgpu plan should be justified against -- not the
eager one.

f64 on the CPU costs another 1.7-1.8× over f32 at the larger grids, and MPS
has no f64 at all, which is the same constraint the WGSL work ran into.

## The Krylov route: torch is nearly level, and that is the interesting part

| grid | C++ 1 thr | C++ 8 thr | torch cpu f32 | torch cpu f32 compiled | torch MPS f32 |
|---|---|---|---|---|---|
| 41³ | 17.9 ms | 25.2 ms | 74.8 ms | 36.0 ms | 69.3 ms |
| 61³ | 65.2 ms | **47.6 ms** | 125.0 ms | 64.5 ms | 194.5 ms |
| 81³ | 166.7 ms | **94.0 ms** | 324.3 ms | 153.6 ms | 252.6 ms |

Compiled on the CPU, torch is level with the C++ on one thread and 1.4-1.6×
behind it on eight. Against the stepping route's 6-8× that is a different
league, and the reason is structural: **the Krylov route does 128 matvecs
where the stepping route does 4 000, so per-call framework overhead is
amortised thirty times less often.** The method that removes the step count
also removes most of torch's disadvantage.

Read the two tables together and the ranking is not what one would guess.
The best number on this machine for the whole propagation is **94 ms**, from
the C++ Krylov solver on eight CPU threads at ng = 81 -- level with the WGSL
stepping kernel on the GPU (99 ms) and 44× the 4 133 ms the library actually
shipped. A better method on a CPU matched a tuned kernel on a GPU.

## Why torch might still be worth a path, on an axis that is not speed

- **Autodiff.** PRD-140's open item is the adjoint: the current one accumulates
  parameter gradients as local products at each step, and a method with no
  steps has none. torch would give it by differentiating through 128 Lanczos
  iterations, which is 128 stored grids -- large but bounded, and free to
  write.
- **CUDA.** The compute box has it; wgpu would reach it too, but torch reaches
  it today with no plugin.
- **It is already installed** wherever the surrogate is worked on.

Against that: PRD-110's refusal stands -- torch cannot enter the conda
dependency list -- so this would be an optional path like the wgpu one, and it
is 3.5-6.7× off the pace on the route it would accelerate.

## Three measurement hazards, all of which produce plausible wrong numbers

1. **`torch.as_tensor` shares memory with the numpy array.** The stepping loop
   writes into its buffers, so the first run overwrote the caller's initial
   density with the propagated one, and every later run on the same fields
   started from a decayed field. It reported a trace uniformly **2× too small**
   -- right shape, right decay, wrong amplitude -- rather than an error.
   `torch.tensor` copies; `benchmark/torch_diffusion.py` uses it and says why.
2. **`torch.compile` silently falls back to eager** when a process compiles the
   same function for too many shape/device combinations: dynamo's
   recompilation limit is hit, a warning goes to stderr, and the timing is
   simply the eager timing. It is detectable only by noticing that the
   "compiled" number equals the eager one to the tenth of a millisecond, which
   is what happened to the ng = 81 MPS row here; the value in the table comes
   from a separate process that compiled only two variants.
3. **torch's bundled libomp collides with the one `libimp_bff` links**, and the
   process aborts on import with OMP error #15. `KMP_DUPLICATE_LIB_OK=TRUE`
   silences it and is documented as unsafe; it was also unstable here.
   `benchmark/torch_diffusion.py` therefore imports no `IMP.bff`, and the C++
   figures above come from a separate run. **This is a packaging fact, not
   only a benchmarking one**: a user who imports both in one process meets it.
