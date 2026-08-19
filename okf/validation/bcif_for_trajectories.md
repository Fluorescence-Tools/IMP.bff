---
type: validation
title: "BinaryCIF replaces DCD and XTC for the rotamer libraries — measured"
description: Whether BinaryCIF can carry rotamer-library coordinates in place of the 95 DCD files and the one XTC. It can, at 1.27 bytes per coordinate against DCD's 4.31 and XTC's 1.60, decoded exactly by the C reader IMP already vendors.
resource: /Users/tpeulen/dev/imp.bff
tags: [validation, imp.bff, io, bcif, rotamer-library, trajectory, dcd, xtc]
timestamp: '2026-08-19T00:00:00Z'
---

# BinaryCIF for the rotamer-library trajectories

**Question.** The rotamer libraries ship as 95 `*.dcd` files (32.33 MB) plus one
`traj.xtc` (12.44 MB). Can BinaryCIF carry the same coordinates, so the package
reads one format instead of three?

**Answer: yes, and it is smaller than both.**

| format | bytes/coordinate | this corpus |
|---|---|---|
| DCD, raw float32 + per-frame headers | 4.31 | 32.33 MB, 95 files |
| XTC, its own 3-D compression | 1.60 | 12.44 MB, 1 file |
| **BinaryCIF**, FixedPoint→Delta→IntegerPacking→ByteArray | **1.27** | **9.87 MB** measured |

Measured on `data/rotamer_library/A48_C1R/traj.xtc`: 28 110 frames × 92 atoms
× 3 = 2 586 120 coordinates.

## It decodes with what IMP already has

`IMP` vendors the C implementation of `ihm` at
`modules/core/dependency/python-ihm/src/` (`ihm_format.c` + `cmp.c`), and it
implements **all seven** BinaryCIF encodings — ByteArray, Delta, FixedPoint,
IntegerPacking, IntervalQuantization, RunLength, StringArray. `bff`'s
`src/CMakeLists.txt` already carries `${PYTHON-IHM_INCLUDE_PATH}`, and
`libimp_atom` already exports the 22 `ihm_*` symbols, so a `bff` translation
unit can read BinaryCIF with **no build change and no new dependency**.

Verified end to end: the file above was written from Python and read back by a
C++ program using that reader — **2 586 120 rows in 68.3 ms, exact round trip**
(`max |decoded - expected| = 0.00 Å`). For comparison, mdtraj loads the XTC in
115 ms.

## The precision it costs, and why it does not matter

The win comes from a **0.1 Å grid** (`FixedPoint` factor 10 per Å). That is
lossy where XTC's 0.01 Å is effectively not. It moves nothing that is used:

| grid | max coordinate error | Δ⟨R_DA⟩ | Δ⟨E⟩ |
|---|---|---|---|
| 0.01 Å | 0.000 Å | 4.8e-7 Å | 1.1e-8 |
| 0.05 Å | 0.020 Å | 1.5e-7 Å | 1.1e-7 |
| **0.10 Å** | **0.050 Å** | **8.1e-4 Å** | **1.7e-5** |
| 0.50 Å | 0.250 Å | 7.7e-4 Å | 1.7e-5 |

⟨E⟩ moves by 2e-5 at the chosen grid. A FRET efficiency is measured to about
1e-2, and the accessible-volume grids these libraries feed are built at 0.5–2 Å.

At 0.01 Å — XTC's own precision — BinaryCIF is 15.52 MB, **25 % larger** than
XTC. The format only wins because the precision can be relaxed, and it can be
relaxed because nothing downstream reads that many digits.

## Two things that will bite whoever writes the encoder

* **The `encoding` array is stored in *encode* order.** The C reader prepends
  each entry as it parses (`ihm_format.c:2075`), so its linked list comes out
  reversed and `decode_bcif_data` walks it forward. Writing the list in decode
  order fails with `FixedPoint not given integers as input`.
* **`IntegerPacking`'s sentinel is a real value.** For int8 the decoder treats
  `127` and `-128` as escape markers, so a *legitimate* delta of exactly ±127
  must still be emitted as an escape plus a remainder. Using `>` instead of
  `>=` in the packing loop produced a file whose three columns decoded to three
  different lengths — `Column size mismatch 2579635 != 2580332`.

`python-ihm`'s own `BinaryCifWriter` implements only ByteArray, Delta,
RunLength and the string/mask encoders — **no FixedPoint and no
IntegerPacking** — so it cannot produce this file. The encoder has to be
written, in Python or C++; the *reader* is free.

## The encoder, and the corpus

`scripts/trajectory_to_bcif.py` writes it. The whole shipped corpus has been
run through it, every file verified by decoding it back and comparing against
the quantised input:

| | before | after | |
|---|---|---|---|
| 95 `*.dcd` | 32.33 MB | **8.12 MB** | 4.0x |
| 1 `traj.xtc` | 12.44 MB | 9.87 MB | 1.26x |
| **total** | **44.78 MB** | **17.99 MB** | **2.5x** |

The DCD files compress far better than the XTC because DCD stores raw float32
and does nothing else; the XTC is already compressed, so BinaryCIF is competing
with a real coordinate codec rather than with a bare array.

> **Correction.** An earlier version of this page estimated the total at "about
> 12 MB". That was arithmetic done from a single file's ratio rather than the
> corpus, and it is wrong: the measured figure is **17.99 MB**.

## Reading it

`IMP::bff::read_bcif_trajectory` (`include/IMP/bff/TrajectoryIO.h`) is the
reader, in C++, over the vendored parser. It takes only the header — the
symbols come from `libimp_atom`, which `bff` already links and which exports
them because `IMP::atom::read_mmcif` compiles the same reader in. Compiling
`ihm_format.c` into `bff` as well would also work and would duplicate it.

Storage is atom-major so the deltas run along an atom's own frame series; the
reader returns frame-major, because every consumer wants `(frame, atom, 3)` and
the storage order exists to make the deltas small rather than to match anyone's
indexing.

Gated against the DCD reader on `A64_C2R_cutoff10` (798 frames × 122 atoms):
identical shape, agreement with the quantised DCD to **3.6e-15 Å**, and a
maximum deviation from the raw DCD of **0.050 Å** — exactly half the grid,
which is the most a 0.1 Å rounding can be wrong by. Read time is **2.9 ms
either way**, from a file a quarter the size.

## Status

**BinaryCIF is the trajectory format from 2026-08-19.** The encoder
(`scripts/trajectory_to_bcif.py`), the C++ reader and their tests
(`test/io/test_bcif_trajectory.py`) are in place.

`traj.xtc` and the 95 DCD files remain in the repository and remain what the
shipped code paths read; they are marked for removal and nothing has been
deleted. What is still open is converting the shipped corpus and repointing
`read_rotamer_library` at the converted files.
