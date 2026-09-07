---
okf_version: "0.2"
---

# AV solver: the timing baseline, and why it is provisional

`benchmark/benchmark_av.py` exists because step 2 of the independent-core work
rewrites the AV hot path — `PathMap` loses its `IMP::em::SampledDensityMap`
base for a bff-owned grid — and the standing rule is that performance only
improves. Until this benchmark there was nothing to measure against:
`benchmark_av_screening.py` times `ProbeNetworkRestraint`, not the solver.

## What it measures

Two doors, because they cost differently and can regress independently:

* **structure** — `compute_av_from_structure`: parse, strip, resample. What a
  caller with a file pays.
* **array** — `compute_av` on the obstacle array. Isolates the path search from
  the parsing, so a regression can be attributed to one or the other.

Three grid steps (1.5 / 1.0 / 0.5 Å) across two structure sizes (T4L, 1293
atoms; mGBP2A, 4606). Cost goes as the voxel count, so 0.5 Å is roughly 27x the
work of 1.5 Å — the sweep exists so a change that is neutral on a coarse grid
and quadratic on a fine one cannot hide.

## The numbers are PROVISIONAL — do not gate on them

Taken 2026-09-07 at `09269d5`, on a machine carrying **load average 14.3** with
two other agent sessions at ~120 % CPU each.

| case | structure (ms) | array (ms) |
|---|---:|---:|
| T4L-132 @1.5 | 7.8 | 4.1 |
| T4L-132 @1.0 | 15.8 | 10.8 |
| T4L-132 @0.5 | 84.1 | 78.1 |
| T4L-150 @1.5 | 9.1 | 5.2 |
| T4L-150 @1.0 | 15.3 | 10.7 |
| T4L-150 @0.5 | 82.8 | 75.4 |
| GBP-400 @1.5 | 22.3 | 10.2 |
| GBP-400 @1.0 | 32.0 | 19.4 |
| GBP-400 @0.5 | 108.5 | 88.4 |

**`--noise` reads 1.12x on this machine**, against a 1.08x gate — that is the
same binary compared with itself, so any verdict `--compare` gives right now is
contention, not code. A self-comparison swung as far as 1.19x.

## The rule this establishes

Run `--noise` **before** taking a baseline or gating a change. If it does not
read close to 1.00x, the box is too busy and the measurement is worthless. This
is why the benchmark reports the floor instead of leaving the caller to assume
one: the failure mode it prevents is a real regression hidden under noise, or —
just as bad — a clean change reverted because the machine was loaded.

    python benchmark/benchmark_av.py --noise                    # is it quiet?
    python benchmark/benchmark_av.py --repeats 7 --json base.json
    python benchmark/benchmark_av.py --repeats 7 --compare base.json

Re-take the baseline on a quiet machine before step 2 lands, and replace the
table above with it, keeping the commit.


## Addendum, 2026-09-07: paired measurement replaces the quiet machine

The machine did not go quiet, so the gate stopped requiring it.
`benchmark/paired_ab.py` measures two builds **alternately within each round**,
so a burst of load lands on both, and reports the median of per-round *ratios*.
Absolute times still wander; the ratio does not.

Each build is a self-contained directory (`IMP/`, `_IMP_bff.so`,
`libimp_bff.0.dylib` with install names rewritten to point inside itself via
`install_name_tool` + `/usr/bin/codesign` -- **Apple's codesign, not conda's
cctools-port shim**, which fails). Neither run touches the shared build tree,
so measuring needs no build lock.

Two things this cost, both worth keeping:

* **Judge on the minimum, not just the median.** A case counts as a regression
  only when its median exceeds the threshold *and its best round was still
  slower*. A real slowdown moves the whole distribution. Without the second
  half, a self-comparison at load 17 reported a 1.07x "regression" on one case
  of eighteen.
* **Eleven rounds, not seven.** At load 13-18 a seven-round run called
  `GBP-400@1.0/array` a 1.081x regression, slower in every round; eleven rounds
  put it at 0.988. A gate whose verdict depends on sample count teaches people
  to re-roll until it passes.

Result for the PathMap change (`v_before` 12:10 vs `v_after` 13:08, 11 rounds,
load 13-18): **worst median 1.038, no case slower in every round.** Medians
span 0.923-1.038 and centre on 1.00.

Caveat on attribution: the two builds are an hour apart in a source tree
several agents share, so the delta is this change *plus* whatever else landed
in that window. The stronger evidence that the change is behaviour-preserving
is the AV oracle -- 18/18 bit-identical density and point-cloud hashes.
