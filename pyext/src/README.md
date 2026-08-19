# `pyext/src` — the Python side, and where it is going

## The shape this is converging on

Every other IMP module has one of two shapes:

| module | py lines | py files | subdirs | cpp lines |
|---|---|---|---|---|
| `atom` | 314 | 1 | 0 | 13,137 |
| `core` | 495 | 2 | 0 | 10,305 |
| `npctransport` | 168 | 1 | 0 | 8,835 |
| `isd` | 5,407 | 13 | 0 | 7,387 |
| `pmi` (the most Python-heavy in IMP) | 23,394 | 26 | 5 | 530 |
| **`bff` (here)** | **~26,000** | **~113** | **~20** | **5,818** |

**C++-carried** — `atom`, `core`, `em`, `isd`: Python is a 1–13 file shim over
7–13k lines of C++. **Python-carried** — `pmi`: 26 files and five subdirectories,
each a genuine family (`dof io plotting restraints topology`).

`bff` is neither. It has four times `pmi`'s file count at similar total Python,
so its files average ~230 lines against `pmi`'s ~900. **The tree is deep because
the files are fragmented**, and no choice of directory scheme fixes that.

**The target is the C++-carried shape.** Numerics keep moving into C++ until the
Python is a thin shim of roughly a dozen flat files; the four stages then become
file names rather than directories. Restructuring happens *last*, once the
Python has stopped shrinking — reorganising a tree that is still losing half its
contents is work that gets redone.

## The four stages

The conceptual model, and the axis the eventual file names follow. Each stage
has a variant per representation, which is why it is a matrix and not a tree:

|  | accessible volume | rotamer library | coarse-grained dye |
|---|---|---|---|
| **1 representation** | the allowed region | listed conformers | beads under a simplified force field |
| **2 scoring** | occupancy | clash energy | the force field's internal energy |
| **3 sampling** | pathfinding, then a walk | library screening | MD / Langevin / RRT |
| **4 analysis** | dye density, and the projection onto an experiment | | |

Full version: [`okf/architecture.md`](../../okf/architecture.md).

Three things cut across every model and are none of the four stages: **`dye`**
(the species), **`label`** (where it is attached, and what quenches it), and
**`photophysics`** (the rate constants). **`observables`** owns the
experiment-neutral output contract — the projection half of stage 4, kept
separate so that contract has somewhere to live. **`io`** speaks the formats and
**`restraints`** scores against measured structural data, which is a different
question from stage 2's "is this configuration allowed".

## What is here now

Directories in transition. `av/`, `fret/` and `dynamics/` are gone — the first
two decomposed into the stages, the third folded into `sampling/`.

* **Stages**: `representation/` (with `av/` and `rotamer/`), `scoring/`,
  `sampling/`, `analysis/`
* **Across the stages**: `dye/`, `label/`, `photophysics/`, `observables/`
* **Formats and data**: `io/`, `restraints/`
* **Assembled models**, applications of the stages rather than stages
  themselves: `quenching/` (the PET model for one site), `cgdye/` (explicit
  all-atom dye MD, off the domain layout — nothing in the package imports it at
  module scope)
* **Utility**: `tools/` (paths to shipped data), `cli/` (what `bin/imp_bff`
  imports)
* `api.py` — the public surface. `BY_DOMAIN` is authored; `EXPORTS` is derived.

## What holds it together

Two tests, and they earn their keep — between them they have caught a cyclic
import, a duplicated κ² module, exports filed under the wrong domain, a
test-module basename collision that silently skipped a file, and four stale
imports that 808 passing tests never touched:

* `test/test_import_discipline.py` — the subpackage graph is acyclic, every
  subpackage imports standalone in a fresh interpreter, **every module imports
  at all**, no domain reaches another's privates, test basenames are unique.
* `test/test_public_api_names.py` — every export resolves, is documented, and
  lives in the domain it is filed under.
