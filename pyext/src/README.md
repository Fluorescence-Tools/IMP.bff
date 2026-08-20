# `pyext/src` — the Python side

## The shape

Every other IMP module has one of two shapes:

| module | py lines | py files | subdirs | cpp lines |
|---|---|---|---|---|
| `atom` | 314 | 1 | 0 | 13,137 |
| `core` | 495 | 2 | 0 | 10,305 |
| `npctransport` | 168 | 1 | 0 | 8,835 |
| `isd` | 5,407 | 13 | 0 | 7,387 |
| `pmi` (the most Python-heavy in IMP) | 23,394 | 26 | 5 | 530 |
| **`bff` (here)** | **~26,900** | **50** | **6** | **~6,000** |

**C++-carried** — `atom`, `core`, `em`, `isd`: Python is a 1–13 file shim over
7–13k lines of C++. **Python-carried** — `pmi`: 26 files and five
subdirectories, each a genuine family (`dof io plotting restraints topology`).

`bff` used to be neither: **113 files averaging 236 lines**, against `pmi`'s 26
averaging 900, at a similar total. The tree was deep because the files were
fragmented, and no choice of directory scheme fixes that — so the files were
merged rather than re-filed. It is now 50 files averaging ~510, which is the
`pmi` shape; the target remains the C++-carried one, and it gets closer every
time a kernel moves.

### What decided each case

A domain became **one flat module** unless it had a reason not to. Two reasons
counted:

* **Something binds its submodules by path.** Of the 99 dotted `IMP.bff.*`
  names referenced across `chisurf`, `imp-tricks`, `quest` and `ucfret`, only
  19 still resolve — and they are concentrated in `quenching`, `restraints` and
  `cgdye`. Those stayed packages.
* **Size.** `representation` is 4,900 lines and `io` is 3,800; one module each
  would be larger than anything in IMP (`pmi/macros.py`, at 2,803, is the
  biggest). They stayed packages holding a handful of substantial modules.

Everything else — `scoring`, `sampling`, `analysis`, `dye`, `label`,
`photophysics`, `observables`, `tools`, `cli` — is one file.

**A CLI is not a stage, and it does not live beside the code it drives.** A
`click` command is a decorated function, so `import click` runs at module scope;
merging `representation/rotamer/cli.py` into `representation/rotamer.py` made
`import IMP.bff.representation.rotamer` require `click`. Every command-line
entry point outside `cgdye` is in `cli.py`, and
`test_import_is_lazy_and_click_free` is what keeps it there.

### What merging found

Names that only ever differed by which module they sat in, and would have
shadowed each other silently:

* **four different `FLRCIF_ITEMS`** — one per dataclass, mapping its fields to
  flrCIF items, in `dye/species`, `dye/spectra`, `label/site` and
  `label/quencher`. Now `DYE_`, `FORSTER_RADIUS_`, `LABEL_` and `QUENCHER_`
  prefixed.
* **two `compute_av`** — one taking arrays, one taking a PDB plus an fps
  position. `IMP.bff.compute_av` resolved to the second while
  `IMP.bff.representation.av.compute_av` resolved to the first, so the *route*
  decided which function you got. The structure one is now
  `compute_av_from_structure`; both flat names keep their meaning.
* **`validate` against a `validate` flag** — `io/fps_schema.validate(payload)`
  merged into a module whose `read_fps_json` takes `validate: bool`, so the
  call became `validate(payload)` on a boolean. Now `validate_fps`.

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

* **Stages**: `representation/` (`av.py`, `rotamer.py`, `distance.py`,
  `pathmap.py`), `scoring.py`, `sampling.py`
* **Across the stages**: `label.py`, `photophysics.py`
* **Formats and data**: `io/` (`cif.py`, `fps.py`, `structure.py`),
  `restraints/` (`docking.py`, `network.py`,
  `simple_av_network.py`, `direct_labeling.py`)
* **Assembled models**, applications of the stages rather than stages
  themselves: `quenching/` (the PET model for one site), `cgdye/` (explicit
  all-atom dye MD, off the domain layout — nothing in the package imports it at
  module scope)
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
