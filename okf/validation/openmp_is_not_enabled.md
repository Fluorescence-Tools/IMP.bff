---
type: validation
title: "Every `#pragma omp` in IMP.bff is inert: the build has no -fopenmp"
description: The kernels carry OpenMP directives over their outer loops, and IMP's CMake leaves OpenMP_CXX_FLAGS empty in the local arm64 build, so none of them thread. Every speed figure recorded for this module is therefore single-threaded against single-threaded numpy, and roughly 8x remains available on this machine. Found 2026-08-19 while a benchmark's arithmetic did not add up.
resource: /Users/tpeulen/dev/imp.bff
tags: [validation, imp.bff, performance, openmp, build, open]
timestamp: '2026-08-19T00:00:00Z'
---

# The OpenMP pragmas do nothing in this build

**Status: recorded, not fixed.** Fixing it means changing IMP's build
configuration, which is outside this repository.

## How it surfaced

Not by inspection. `expected_rmsd_after_adding` at 700 frames and 120
candidates took 1311 ms for 5.9 × 10⁷ tail evaluations. At the ~20 ns each one
costs that is about 1.2 s of work — which is the whole wall clock, on a machine
with 8 cores and a `#pragma omp parallel for` over the candidate loop. The
arithmetic only closes if exactly one thread ran.

## The cause

```
$ grep OpenMP ~/dev/imp/cmake-build-arm64/CMakeCache.txt
OpenMP_CXX_FLAGS:STRING=
OpenMP_CXX_INCLUDE_DIR:PATH=/Users/tpeulen/mambaforge/envs/arm64/include
OpenMP_CXX_LIB_NAMES:STRING=libomp
```

`libomp` was found and the include path resolved, but the **flags string is
empty**, so no `-fopenmp` reaches the compiler. This is the usual AppleClang
case: it needs `-Xpreprocessor -fopenmp`, which CMake's `FindOpenMP` does not
always fill in.

A `#pragma` a compiler does not recognise is a comment. Nothing warns. The code
reads as parallel and runs as serial, and the two are indistinguishable from
the source.

## What it means for every number in this repository

**Every speed figure recorded for `IMP.bff` is single-threaded.** They are
honest comparisons — single-threaded C++ against single-threaded numpy — and
they are lower bounds:

| kernel | recorded | threads used |
|---|---|---|
| rotamer interaction energies | 7.9× | 1 |
| mean-field pair energies | 47× | 1 |
| fused walk → rate → photons | 10.5× | 1 |
| greedy Olga candidate scoring | 3.1× | 1 |

One earlier claim needs correcting with this in hand: the quenching suite's
recovery from 558 s to 18 s was attributed partly to "OpenMP parallelisation
over x-slabs". It was not. The whole of that win was moving the per-step loop
out of Python; the threading contributed nothing, because there was none.

## Making it visible

`IMP.bff.built_with_openmp()` and `IMP.bff.parallel_threads()` report the truth
at runtime, and `test/fret/test_greedy_olga.py::test_the_speed_figures_are_single_threaded`
asserts the two agree. If the build ever gains OpenMP that test still passes —
but the assertion documents that the recorded figures became lower bounds.

## What to do

Enabling it is a change to `~/dev/imp`'s CMake configuration, which this
repository does not own and which affects every IMP module. The pragmas are
already correct and already placed over the right loops — the outer loop in
every case, with per-thread accumulators where a reduction is involved — so the
day the flag appears, they work.

Worth checking first that the reductions really are safe under threading. Two
were written with that in mind and say so: `quenched_decay`'s blocked
accumulation (the bin index is data-dependent) and `expected_rmsd_after_adding`'s
per-candidate column accumulators. The rest write to disjoint output slots and
are safe by construction.
