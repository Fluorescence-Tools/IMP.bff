---
type: validation
title: "The curve engine is not a WGSL candidate: the crossing costs 1.3 ms and the curve costs 24 us"
description: Measured before writing any shader. A whole model curve at the lengths that occur (128-4096 points) evaluates in 0.8-24 us on the CPU, and the smallest job the GPU plugin will accept costs 1304 us of which almost all is the crossing. That is 57x the wrong way. Batching does not rescue it below ~50 curves per crossing, and once the transfer is counted even 500 curves leaves under 2x. The win, if there is one, is residency rather than a kernel, and it is a redesign of the fit loop. FactorGraph proper is structural and has no arithmetic at all. The accessible-volume search, at 2.8-50 ms a site, is the better next candidate.
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

## Where a win would actually come from, if anyone wants one

Not from a kernel. From **residency**. Across a fit's iterations the curve
axis does not change and neither does the data; only a handful of parameters
per dataset do. If those lived on the device and the only traffic per
iteration were the parameters in (a few kilobytes) and chi² back (one number,
reduced on the device), the transfer term disappears and the ceiling becomes
the arithmetic ratio -- of order 50× at 500 datasets.

That is a redesign of the fit loop and it makes the expression engine own
device memory and a lifetime. It is not a shader, and it should not be started
as one. It is also the same shape as the open question in PRD-140: the Krylov
solver would likewise be dominated by what crosses rather than by what
computes.

## A better next candidate

`get_av`, the accessible-volume search, on a 3 000-atom obstacle field:

| grid | time |
|---|---|
| 1.5 Å | 2.8 ms |
| 1.0 Å | 7.5 ms |
| 0.7 Å | 19.4 ms |
| 0.5 Å | 50.5 ms |

At the resolutions that matter this clears the 1.3 ms crossing by two to
thirty times, which is the first condition the curve engine fails. It is a
harder kernel -- a masked path search rather than a stencil, and the CPU side
already has a NEON inner loop -- but it is at least the right size, and it is
called once per site rather than once per fit iteration.

## Housekeeping

`src/standalone/Expression.cpp` and `include/Expression.h` were being modified
by another agent when this was measured, so the absolute expression timings
are against that working tree rather than against HEAD. The conclusion does
not depend on them: it would take a hundredfold regression in the engine
before a GPU crossing became competitive.
