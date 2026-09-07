---
type: note
title: "Where data belongs — IMP's three places, and where bff puts things"
description: IMP modules keep runtime parameters in data/, test fixtures in test/input/, and example inputs in examples/. bff has 54 MB in data/ against every other module's 1-14 flat files, no test/input/ at all, and a 17.3 MB unreferenced archive in examples/. What each of bff's trees is, and where it should go.
resource: /Users/tpeulen/dev/imp.bff
tags: [note, imp.bff, data, layout, conventions, packaging]
timestamp: '2026-08-19T00:00:00Z'
---

# Where data belongs

## IMP's three places

Surveyed across `../imp/modules`, and they are not interchangeable:

| directory | what it holds | installed? | typical size |
|---|---|---|---|
| `<module>/data/` | **runtime parameters** the installed library reads — `atom/data/par.lib`, `top.lib`, `radii.lib` are the CHARMM force field | **yes**, via `get_data_path()` | 1–14 files, **flat** |
| `<module>/test/input/` | test fixtures | **no** | `em` 68 MB, `pmi` 24 MB, `atom` 10 MB |
| `<module>/examples/` | example scripts and the inputs they need | **yes**, via `get_example_path()` | 224 KB – 4.5 MB |

Two things follow that are easy to miss:

* **`data/` is flat everywhere.** Of the fifteen IMP modules that ship data,
  eleven have **zero** subdirectories and the rest have one or two. `atom`
  keeps the whole CHARMM force field as six files at the top level.
* **`test/input/` is where bulk belongs.** It is not installed, which is why
  `em` can keep 68 MB there without it reaching a single user. It is the
  pressure valve that keeps `data/` small.

## Where bff stands

| | files | subdirs | size |
|---|---|---|---|
| every other IMP module's `data/` | 1–14 | 0–2 | ≤ 12 MB |
| **`bff/data/`** | **290** | **15** | **54 MB** |
| `bff/test/input/` | — | — | **does not exist** |
| `bff/examples/` | 82 | | **24 MB** |

`data/` is symlinked into the build tree entry by entry and installed, so all
54 MB reaches every conda-forge user of IMP.

## What is actually in it

| tree | size | reads it | belongs |
|---|---|---|---|
| `data/rotamer_library/` | 45 MB | the package, by dye name | `data/` — it *is* runtime reference data |
| `data/cgprobe/templates/` | 8.2 MB | the package, via `get_template_dir()` | `data/` — dye templates are parameters |
| `data/cgprobe/inputs/structures/` | 1.0 MB | the **CLI**, as default inputs | `examples/` |
| `data/cgprobe/inputs/dyes`, `test_systems` | 32 KB | nothing, by path | fixtures |
| `data/cgprobe/inputs/restraints/` | 32 KB | **nothing** — orphaned when the NMR restraints were retired | moved to `junk/nmr_restraints/data/` |
| `examples/structure/GBP/mGBP2_AF_dimer.result.zip` | **17.3 MB** | **nothing** | it is a *result*, not an example input |

That last row is 72 % of `examples/`, is installed, and no file in the
repository references it.

## The thing that actually fixes `data/`

The 45 MB rotamer library is the bulk, and it is legitimately runtime data — so
it cannot be moved out. It can be made **small**:
[`bcif_for_trajectories.md`](validation/bcif_for_trajectories.md) measured the
95 DCD files at 32.33 MB → **8.12 MB** and the XTC at 12.44 MB → 9.87 MB, all
verified exact. Converting the corpus takes `data/` from 54 MB to about 27 MB
without moving anything or changing what the package can read.

Reorganising is the smaller half of the problem. Re-encoding is the larger one.

## Open, and needing a decision

* **The 17.3 MB archive.** Unreferenced and installed. It may be a deliberate
  reference dataset for a publication, which is why it has not been touched.
* **`data/cgprobe/inputs/structures/`.** The CLI defaults to reading `1DG3.pdb`
  and `alexa488_r48.mol2` from there through `IMP.bff.get_structure_dir`.
  Under IMP's convention those are example inputs, but moving them changes what
  `get_structure_dir()` means and touches both CLIs, so it is a code change
  rather than a file move.
