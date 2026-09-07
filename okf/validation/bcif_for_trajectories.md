---
type: validation
title: "BinaryCIF for the rotamer libraries — what it actually buys, after two wrong answers"
description: BinaryCIF replaced the 95 DCD rotamer libraries. It is lossless float32 and 7 % smaller, not the 4x a first pass claimed — that figure came from a 0.1 A grid that shifts transition-dipole directions by 1.6 degrees and breaks the FRETpredict pins. The XTC left data/ because nothing reads it, not because BinaryCIF beat it.
resource: /Users/tpeulen/dev/imp.bff
tags: [validation, imp.bff, io, bcif, rotamer-library, trajectory, dcd, xtc]
timestamp: '2026-08-19T00:00:00Z'
---

# BinaryCIF for the rotamer libraries

> **Update (2026-08-24).** The shipped *rotamer libraries* moved on again, to
> `.drot` (PRD-118) -- an internal-coordinate store that carries its own
> template, at 0.65x the BinaryCIF. BinaryCIF remains this package's
> **trajectory** format and the `.bcif` libraries still ship and still read.
> The lesson below is what decided `.drot`'s default too: the lossless rung is
> the default there for exactly the dipole-direction reason recorded here, and
> the torsion-only draft of `.drot` failed the same pins in the same way.

**BinaryCIF is the trajectory format as of 2026-08-19.** The 95 DCD libraries
in `data/rotamer_library/` are now `.bcif`, read by
`IMP::bff::read_bcif_trajectory` through the C parser IMP already vendors.

The size result is modest and the road to it went through two wrong answers,
both recorded here because each was wrong in a way that looked right.

## What it buys

| | | |
|---|---|---|
| 95 DCD libraries | 32.33 MB | 4.31 bytes/coordinate |
| 95 BinaryCIF, **lossless float32** | **30.05 MB** | 4.01 bytes/coordinate |

7 % smaller, **bit-exact**, and one format instead of two. That is the whole
size story for the shipped data. What BinaryCIF really buys here is not bytes:
it is that the C++ side reads the libraries directly, through a parser that is
already in the build, with no new dependency.

## The first wrong answer: 0.1 A, and 4x

A first pass quantised to a 0.1 A grid and reported **8.12 MB, 4.0x smaller**.
The grid was chosen by measuring how far a mean FRET efficiency moved when the
*centre of mass* of a 92-atom dye was perturbed: 1.7e-5, comfortably below the
~1e-2 an experiment resolves.

**That was the wrong quantity.** kappa^2 comes from transition-dipole
*directions*, which are differences between atoms about 1.7 A apart. A 0.05 A
coordinate error is then ~1.6 degrees of orientation error:

| grid | max coordinate error | median dipole angle error |
|---|---|---|
| 0.001 A | 0.0005 A | 0.016 deg |
| 0.01 A | 0.005 A | 0.158 deg |
| 0.1 A | 0.050 A | **1.582 deg** |

The FRETpredict reference values this package is pinned against moved by
**2.3e-3 in E, against a 2e-5 tolerance**. Averaging over 92 atoms hid it; a
dipole does not average. The lesson is not "0.1 A is too coarse" — it is that
a proxy chosen for convenience measured something the code does not compute.

## The second wrong answer: a finer grid

Refining to 0.001 A brought the pins to within 2.3e-5 to 5.2e-5 — still outside
their 1e-5 tolerance — and the *file grew to 256 MB*, eight times the DCDs.

Delta encoding was the cause. It assumes the next value resembles the last, and
in a rotamer **library** that is false: consecutive frames are independent
conformers, not a time series. The deltas are as large as the coordinates,
IntegerPacking spends its output on escape runs, and the file explodes.

The encoder now measures both chains and keeps the smaller
(`scripts/trajectory_to_bcif.py`).

## The thing that made the grid question moot

| grid | max int | type | total | bytes/coordinate |
|---|---|---|---|---|
| 0.001 A | 28 089 | int16 | 15.05 MB | 2.01 |
| 0.005 A | 5 618 | int16 | **15.05 MB** | 2.01 |
| 0.01 A | 2 809 | int16 | **15.05 MB** | 2.01 |

**The size is set by the integer type, not by the grid.** The largest
coordinate in the corpus is 28.1 A, so the quantised values fit `int16` at
every grid from 0.001 A down. Coarsening from 0.001 A to 0.005 A buys **zero
bytes** and costs five times the error. There is no reason to quantise more
coarsely than 0.001 A, and no reason to quantise at all unless 15 MB against 30
matters more than the pins do.

Quantisation is therefore an *option* in the encoder, not the default. The
shipped libraries are lossless.

## The XTC

`data/rotamer_library/A48_C1R/traj.xtc` (12.44 MB) has left `data/` — **but not
because BinaryCIF beat it.** It does not: XTC's own coordinate codec gives 1.60
bytes/coordinate at 0.01 A, and lossless BinaryCIF would be 31 MB, two and a
half times larger.

It left because **nothing reads it**. It is source material — the MD trajectory
the library was derived from — and the installed `data/` directory is for what
the library reads at runtime (see [`../data-layout.md`](../data-layout.md)).
Recoverable from commit `8fac573`.

## Two traps in the encoding, for whoever touches this next

* **The `encoding` array is stored in *encode* order.** The C reader prepends
  each entry as it parses, so its list comes out reversed and it decodes
  correctly. Writing decode order fails with `FixedPoint not given integers as
  input`.
* **IntegerPacking\'s int8 sentinels are legitimate values.** A delta of exactly
  +-127 must be emitted as an escape plus a remainder, or the three columns
  decode to different lengths — `Column size mismatch`.
