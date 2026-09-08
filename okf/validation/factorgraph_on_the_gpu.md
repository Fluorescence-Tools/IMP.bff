---
type: validation
title: "The curve engine on the GPU: hopeless per call, 9-78x if it stays resident"
description: Measured before writing any shader, then measured again when asked what happens if the graph stays resident. Called per curve it is hopeless: the curve costs 24 us and a crossing costs 1.3 ms. Batched but still crossing with the data, under 2x. Resident, with only parameters going in and chi2 coming back, 9.1x at 200 datasets and 58.6x when a Levenberg-Marquardt step's eight independent evaluations share one synchronisation. The wall is not bandwidth or arithmetic but a flat ~1.2 ms synchronisation, measured three ways; recording the work costs 46-72 us whatever its size. FactorGraph proper is structural and has no arithmetic at all.
resource: /Users/tpeulen/dev/imp.bff
tags: [validation, imp.bff, performance, gpu, webgpu, expression, factorgraph]
timestamp: '2026-09-08T00:00:00Z'
---
# What the factor graph evaluates, and why a GPU cannot help

Asked whether the factor-graph evaluation is the next thing to put in WGSL.
Measured first, because the last two candidates were decided by measurement
and this one turns out the other way.

## Two different things wear that name

`include/FactorGraph.h` is **structural**: moralisation, a greedy elimination
order, cliques, treewidth, connected components. It manipulates sets of
variable keys. There is no arithmetic in it to accelerate and no array
anywhere -- a GPU has nothing to do here at all.

The thing that *evaluates* is the dataflow side: `Node` and `Port`, whose
callbacks are elementwise add and multiply over `std::vector<double>` ports,
and `Expression`, which evaluates a model equation over a curve. That is the
candidate, and it is the one measured below.

## The numbers, and they are not close

A whole model curve, `benchmark/expression_curves.py`, at the lengths that
actually occur:

| equation | 128 | 512 | 1024 | 2048 | 4096 |
|---|---|---|---|---|---|
| two exponentials | 1.6 µs | 3.8 µs | 7.0 µs | 13.3 µs | 24.4 µs |
| mixed transcendental | 1.0 | 1.5 | 2.2 | 3.5 | 6.3 µs |
| a line | 0.8 | 1.0 | 1.2 | 1.7 | 2.6 µs |
| FCS | 1.0 | 1.7 | 2.6 | 4.4 | 8.2 µs |

Against that, **what one crossing costs on this machine**, measured two ways:

| | |
|---|---|
| upload *n* floats, one trivial kernel, read *n* back (wgpu-py) | 1 393 µs at n=128, 1 397 µs at n=4096 |
| the smallest job the C plugin accepts (ng=12, 32 steps -- 8 224 voxel-steps, so the arithmetic is nothing) | **1 304 µs** |

Both agree: **a crossing costs about 1.3 milliseconds and hardly depends on
how much data goes with it.** The most expensive real curve costs 24
microseconds. The GPU is fifty-seven times the wrong way round, and no shader
changes that, because none of the 1.3 ms is arithmetic.

Note also what the curve lengths mean for parallelism: 4 096 points is sixteen
workgroups. The device is idle either way.

## Batching, and why it still does not pay

The engine takes columns, so a global fit could concatenate every dataset's
curve into one call with the per-dataset parameters passed as columns rather
than scalars. Measured on the two-exponential model at 4 096 points a curve:

| curves | points | CPU | ceiling against a 1.3 ms crossing |
|---|---|---|---|
| 1 | 4 096 | 27 µs | 0.02× |
| 10 | 41 k | 256 µs | 0.17× |
| **50** | **205 k** | **1 391 µs** | **1.0× — break-even** |
| 100 | 410 k | 2 868 µs | 1.9× |
| 500 | 2.0 M | 15 533 µs | 10× |
| 2 000 | 8.2 M | 62 925 µs | 42× |

**Those are ceilings, not speedups**: they assume the arithmetic and the
transfer are free. The transfer is not. At 500 curves the four input columns
are 32 MB and the result 8 MB, and this machine moved 1.6 MB for about 370 µs
of the round trip above, so roughly 4.4 GB/s -- about 9 ms for 40 MB, against
15.5 ms of CPU. **Under 2×**, for a fit-loop rewrite.

## Residency: measured, and it changes the verdict

Asked what happens if the graph *stays* on the GPU rather than crossing per
call. It is a different problem, and the answer is much better than the
paragraph that stood here before it was measured.

Across a fit's iterations the curve axis does not change and neither does the
data or the weights; only a handful of parameters per dataset do. Put those on
the device once, send the parameters in and reduce chi² there, and the only
traffic per iteration is a few kilobytes in and four bytes out.
Two-exponential model, 4 096 points a curve, model and residual and reduction
fused into one kernel:

| datasets | CPU | GPU, one evaluation per sync | vs CPU | recording only, no wait |
|---|---|---|---|---|
| 10 | 376 µs | 1 411 µs | 0.27× | 45 µs |
| 50 | 2 612 µs | 1 406 µs | 1.9× | 52 µs |
| 100 | 5 305 µs | 1 340 µs | 4.0× | 50 µs |
| 500 | 24 740 µs | 1 413 µs | 17.5× | 49 µs |
| 2 000 | 110 898 µs | 1 372 µs | 41× | 72 µs |

Read the last column. **Recording the work costs 46–72 µs whether it is ten
curves or two thousand** — 8.2 million points of transcendental arithmetic
lands in the noise. Everything else in the middle column is one fixed cost.

## The wall is a synchronisation, not bandwidth and not arithmetic

That fixed cost is ~1.2 ms and it is neither the transfer nor the map.
Measured three ways that agree:

| | |
|---|---|
| upload, trivial kernel, download (wgpu-py) | 1 393 µs at n=128, 1 397 µs at n=4 096 |
| submit and wait for completion, **mapping nothing** | 1 102–1 335 µs, flat in problem size |
| the smallest job the C plugin accepts (ng=12, 32 steps) | 1 304 µs |

So it is the price of asking the device "are you finished" once. Whether that
is intrinsic or is wgpu's polling granularity is **not settled here** — 1.2 ms
is very long for a Metal submit, and if it were reducible to tens of
microseconds every number above improves by an order of magnitude. Worth
finding out before building anything.

## The design that follows: more evaluations per synchronisation

If one synchronisation costs 1.2 ms, put more independent work inside it. A
Levenberg–Marquardt step needs the value and one finite difference per
parameter: `n_params + 1` evaluations that do not depend on each other and can
all be recorded before the single wait. Measured at 200 datasets × 4 096
points, K independent parameter sets per crossing:

| K | total | per evaluation | vs CPU |
|---|---|---|---|
| 1 | 1 513 µs | 1 513 µs | 9.1× |
| 2 | 1 536 µs | 768 µs | 17.9× |
| 4 | 1 626 µs | 406 µs | 33.8× |
| **8** | **1 874 µs** | **234 µs** | **58.6×** |
| 16 | 2 927 µs | 183 µs | 75.1× |
| 32 | 5 626 µs | 176 µs | 78.1× |

K = 8 is a seven-parameter LM step. It is worth **58.6×**.

## So what would actually have to be built

Not a shader for one equation — the measurement above hand-wrote the
two-exponential model. The engine would have to **compile an expression to
WGSL**, which is a code generator over the existing AST, and own device
buffers with a lifetime across fit iterations. Both are real work and neither
is the kind of thing to start before the 1.2 ms synchronisation is understood,
because that single number decides whether the payoff is 9× or 90×.

It is also the same shape as the open question in PRD-140: the Krylov solver
would likewise be dominated by what crosses rather than by what computes.

## A candidate that needs none of that

`get_av`, the accessible-volume search, on a 3 000-atom obstacle field:

| grid | time |
|---|---|
| 1.5 Å | 2.8 ms |
| 1.0 Å | 7.5 ms |
| 0.7 Å | 19.4 ms |
| 0.5 Å | 50.5 ms |

At the resolutions that matter this clears the 1.3 ms synchronisation by two
to thirty times on its own, so it needs no residency and no compiler -- it is
a single call with a single answer, the shape the diffusion backend already
has. Harder kernel (a masked path search rather than a stencil, and the CPU
side already has a NEON inner loop), but the easier project.

## Housekeeping

`src/standalone/Expression.cpp` and `include/Expression.h` were being modified
by another agent when this was measured, so the absolute expression timings
are against that working tree rather than against HEAD. The conclusion does
not depend on them: it would take a hundredfold regression in the engine
before a GPU crossing became competitive.
