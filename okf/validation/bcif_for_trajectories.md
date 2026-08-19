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

## Status

Nothing has been converted or deleted. `traj.xtc` and the 95 DCD files are
**marked for removal**, pending the encoder and a conversion of the shipped
data. Estimated result: 44.8 MB of trajectory data becomes about 12 MB in one
format.
