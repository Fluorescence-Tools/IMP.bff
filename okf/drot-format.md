---
type: reference
title: "`.drot` v10 — the rotamer-library container, byte for byte"
description: The format the shipped dye rotamer libraries use since 2026-08-24 - a PTO (EBML) container of a template, a Z-matrix and per-conformer internal coordinates. What each object holds, how the numbers are laid out, which versions read, and why the default rung is lossless.
resource: /Users/tpeulen/dev/imp.bff
tags: [reference, imp.bff, drot, rotamer-library, format, pto, ebml, brotli, zmatrix, prd-118]
timestamp: '2026-08-24T00:00:00Z'
---

# `.drot` v10

A `.drot` file is one dye+linker rotamer library: an ensemble of conformers
with a weight each, stored as internal coordinates. It is written by
`IMP::bff::write_drot` (`include/DrotWriter.h`, program
`bin/imp_bff_traj2drot`) and read by `IMP::bff::read_drot`
(`include/DrotReader.h`). Both sides go through ptolib's built-in
brotli codec (`compress_bytes`/`decompress_bytes` under the name `"brotli"` —
ptolib 0.4.0 embeds the codec in the vendored PTO core
(`include/internal/ptolib.h`, a verbatim copy of
[ptolib](https://github.com/tpeulen/ptolib), the header tttrlib carries too;
`include/Pto.h` is the thin `PtoWriter`/`PtoReader` face over it), so neither
adds a dependency — the container is shared with tttrlib by one header and with
chimol by specification, not by linkage
([the convergence record](../../chimol/okf/references/pto-chm-convergence.md)). The design ledger — what was measured and what was rejected — is
[PRD-118](prds/prd-118.md); the prototype that arrived at the format is
`prototypes/drot_rotlib/` and its `CONVERTING.md` describes the older v8.

## Three families, three containers

The rotamer libraries ship as one container per family, because a family is
what changes together:

| file | holds | size |
|---|---|---|
| `dyes.drot.pto` | 95 dye+linker ensembles (FRETpredict; Montepietra et al., GPLv3) | 19.7 MB |
| `spinlabels.drot.pto` | 10 spin-label ensembles (DEER-PREdict; same group, GPLv3) | 0.24 MB |
| `sidechains.drot.pto` | the Dunbrack-2010 backbone-dependent table, 18 residues | 3.3 MB |

Inside, each library is a set of objects named `<library>/<member>`, and a
`drot.catalog` object at the head lists them. One library is addressed by a
**locator** — `dyes.drot.pto::A48_C1R_cutoff10`, the container, `::`, the
library — which is what `resolve_rotamer_library_path` answers with and what
`read_drot` takes. A plain path still means a container holding one library,
which is what a user's own file is.

Bundling is a copy: `write_drot_bundle` moves payload bytes across unread, so
a family container holds exactly what its parts held (measured: 19,675,419
bytes against the 19,683,305 of the 95 separate files — the difference is
their EBML headers, now shared). Reading one library out of it costs that
library, because the reader walks the framing by seeking and pulls only the
objects it was asked for; the old objection to one archive — that reading one
member means inflating all of them — was a property of the tar, not of
bundling.

**The side-chain family is a different payload in the same envelope.** A dye
library is an *ensemble*: conformers with weights. A side-chain library is a
*distribution over backbone conformation*: for each residue and each (phi,
psi) bin, the rotamers with their chi angles, their spreads and a probability.
So it carries kinds of its own — `rot.bbdep.header` and, per residue,
`rot.bbdep.records` — and no element ID was invented for it either. Its
records are FASPR's own 20-byte layout, byte for byte, which makes
`write_dunbrack_bin` an exact inverse and lets the vendored FASPR engine (it
seeks in that file) run against the shipped container with its parity pin
intact. See [`include/DunbrackLibrary.h`](../include/DunbrackLibrary.h).

## The container

A **PTO** document — EBML (RFC 8794), `DocType "pto"` — written through
[`include/Pto.h`](../include/Pto.h), the face over the vendored
[ptolib](https://github.com/tpeulen/ptolib) header (`include/internal/ptolib.h`).
Libraries written before 2026-09-07 carry one SeekHead and every object in one
`Attachments` element; ptolib reads those as it reads its own, read-only. One file per library, named
`<stem>.drot.pto`: `.pto` because that is the container the whole stack shares
(`.mmfbd.pto` for photons in tttrlib, `.chm.pto` for structures in chimol), and
`.drot` because a profile suffix names the primary payload kind so a human, a
shell and a file dialog can all see which profile a container carries. The
grammar is untouched by this profile — the members are attached objects whose
`PtoKind` says what they are, which is the extension point PTO already
provides, so no element ID was invented here.

Each object's bytes are brotli (quality 11, window 24) on their own rather than
one stream over all of them. Measured over the 95 shipped libraries that costs
nothing — 19.60 MB against 19.64 MB — because tar's 512-byte headers go away and
pay for the joint context that is given up; and it buys a reader the header
without inflating the grids. Payloads sit on 8-byte boundaries (`Void`
padding), so a mapping reader hands out slices instead of copies.

In a family container every name below is prefixed `<library>/`.

| object (`FileName`) | `PtoKind` | holds |
|---|---|---|
| `drot.catalog` | `drot.catalog` | the libraries this container holds, as JSON — the only object a listing reads |
| `drot.json` | `drot.header` | the header: version, shapes, the Z-matrix base, one entry per grid, the weight encoding |
| `template.cif` | `drot.template` | a minimal mmCIF `_atom_site` loop — atom names, residue names (optional), elements, and the first conformer's coordinates |
| `rows.json` | `drot.rows` | the Z-matrix as a flat array, four indices per row: `atom, parent, grandparent, great-grandparent` |
| `base.f32` | `drot.grid` | base-atom coordinates, `n_base * 3` per conformer |
| `lens.f32` | `drot.grid` | bond length per row per conformer |
| `theta.f32` | `drot.grid` | bond angle per row per conformer, degrees |
| `phi.f32` | `drot.grid` | dihedral per row per conformer, degrees |
| `weights.bin` | `drot.weights` | one weight per conformer, varint or float32 |

Every object also carries a `PtoEncoding` naming how its bytes are coded —
`f32.col+brotli`, `json+brotli`, `varint+brotli` — and the grammar never looks
inside a payload, which is what keeps one container serving three domains.

A v10 file begins with the EBML magic `1A 45 DF A3`, so a reader tells it from
the older envelope by looking. v5–v9 were a ustar tar compressed as one brotli
stream (`brotli -dc x.drot | tar -t`), and brotli streams carry no magic, so a
reader that must sniff those looks for the xz magic first (v2–v7 used xz) and
falls back to brotli.

Nothing in the file carries a timestamp, so the same ensemble always produces
the same bytes — re-encoding a corpus is a no-op in a diff.

**Conformance is not checked against our own reader.** A consistent misreading
of RFC 8794 would satisfy writer and reader together, so a `.drot.pto` is
validated by `../tttrlib/test/tools/pto_ebml_check.cpp`, which walks it with
**libebml** — the reference implementation Matroska is built on — and by
chimol's independent Python walker (`../chimol/chimol/render/pto.py`). Every
shipped container passes both, with `--aligned` enforced.

That checker wants **libebml 2.0**, and says so nowhere: it uses
`EbmlId::FromBuffer` and expects `EMaxSizeLength`/`EDocType` from
`EbmlHead.h`, neither of which is in 1.4.7 (what Homebrew installs today), so
against 1.4.7 it fails to compile rather than to validate. There is no
packaged 2.0 to install either — build it from the clone:

```sh
cmake -S ../tttrlib/junk/libebml -B /tmp/ebml-build -DCMAKE_INSTALL_PREFIX=/tmp/ebml \
      -DBUILD_SHARED_LIBS=OFF && cmake --build /tmp/ebml-build -j && cmake --install /tmp/ebml-build
c++ -std=c++17 -I/tmp/ebml/include ../tttrlib/test/tools/pto_ebml_check.cpp \
    -L/tmp/ebml/lib -lebml -o /tmp/pto_ebml_check
```

## `drot.json`

```json
{"format": "drot", "version": 10, "container": "pto",
 "producer": "IMP.bff write_drot",
 "n_atoms": 83, "n_rotamers": 711, "n_rows": 73, "n_base": 10,
 "base_offset": [46, 47, 48, 50, 79, 49, 51, 52, 53, 80],
 "grids": {
   "base":  {"file": "base.f32",  "dtype": "f32", "order": "col", "shuffle": 1},
   "lens":  {"file": "lens.f32",  "dtype": "f32", "order": "col", "shuffle": 1},
   "theta": {"file": "theta.f32", "dtype": "f32", "order": "col", "shuffle": 1},
   "phi":   {"file": "phi.f32",   "dtype": "f32", "order": "col", "shuffle": 1}},
 "weights": "varint:711"}
```

`base_offset` gives the atom index of each base atom, in the order `base.f32`
stores them. A grid entry with `"dtype": "i16"` carries `"scale"` as well, the
step one count is worth. `"weights"` is `varint:<n>` (LEB128, what integral
cluster populations encode as) or `f32:<n>`.

## How the numbers are laid out

**Column-major**: the value of degree of freedom `d` for conformer `k` sits at
`d * n_rotamers + k`. All conformers of one dihedral are contiguous. This was
measured, not assumed: same-row values cluster around their own mean, and
contiguity is what the entropy coder's context model needs to see — −7.4 %
corpus-wide over the rotamer-major order.

**Byte-plane shuffled** (float rungs): a `n * 4`-byte array is written as all
first bytes, then all second bytes, and so on. The exponent bytes of a column
are nearly constant once separated from the mantissa — −8 %.

**The degrees of freedom are complete.** Per conformer the file holds its own
base-atom coordinates *and* its own bond length, bond angle and dihedral per
row. `ZMatrix::decode_internals` rebuilds from exactly that and reproduces the
source to ~1e-13 Å; the only error left is what the number format rounds off.

## Rungs

| rung | grids | size vs `.bcif` | reconstruction |
|---|---|---|---|
| lossless (default) | float32, shuffled | 0.65x corpus-wide (0.47x–0.92x per library; 1.8x for the one single-rotamer library, where the template outweighs the frames) | ~1e-6 Å mean, ~1e-5 Å max |
| compact (`--grid`) | int16 at `grid_a` / `grid_deg` | ~0.30x | ~3e-3 Å |

**The default is lossless and should stay that way.** κ² comes from
transition-dipole directions — differences between atoms about 1.7 Å apart —
so a coordinate error that looks harmless as a displacement is not harmless as
an angle, and it does not average away over an ensemble. The same reasoning
already decided `imp_bff_traj2bcif`'s default
([`validation/bcif_for_trajectories.md`](validation/bcif_for_trajectories.md)),
and the torsion-only draft of this format proved it again: rebuilding on
*template* bond lengths cost 0.02 Å and moved the FRETpredict parity pins by
1.3e-3 in E against a 2e-5 tolerance.

## Versions

`read_drot` reads v5, v7, v8, v9 and v10. **v10 is v9's payload in the PTO
envelope** — every number is laid out identically, and the version bumped
because a v9 reader cannot find `drot.json` in a file it cannot walk. v9 is the
same members in a brotli-compressed tar. v5–v8 are the prototype's: int16 grids
with no `lens` member, so they rebuild on template bond lengths (v5 stores
grids rotamer-major, v7/v8 column-major, v8 is v7 in a brotli container rather
than xz). Only v10 is written. v2–v4 and v6 are not read here; the prototype's
`validate_drot.py` handles those.

## Reading it without IMP.bff

EBML is a handful of variable-width integers, and the specification
(`../tttrlib/okf/specs/pto-binary-decoding.md`) ships a ~200-line C99 walker for
exactly this. In Python, any PTO reader does it — here with chimol's, which was
written from that spec and shares no code with ours:

```python
import brotli
from chimol.render.pto import PtoReader          # or tttrlib's, or your own

with PtoReader("A48_C1R_cutoff10.drot.pto") as r:
    for o in r.objects():
        print(o.name, o.kind, o.encoding, o.offset, o.size)
    header = brotli.decompress(bytes(r.data(r.find("drot.json"))))
```

The grids are plain little-endian arrays once un-shuffled; `rows.json` and
`template.cif` are text. Nothing in the format needs the writer to read it
back, which is the point of storing the Z-matrix rather than re-deriving it.
