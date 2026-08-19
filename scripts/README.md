# `scripts/`

One-shot tools that operate on the repository's **data**, not on its code.

Nothing here is imported, built, installed or exercised by the test suite.
These are run by hand, occasionally, usually once — a conversion, a migration,
a regeneration. That is what separates this directory from the others:

* `pyext/src/` is the package;
* `bin/` is what an end user gets on their `PATH`;
* `utility/` is IMP's own convention for programs used *during the build*;
* `scripts/` is for changing what is in `data/`.

A script here may depend on things the package does not — `mdtraj` for reading
an XTC, for instance. That is the point of keeping them apart: the package's
dependency list is a public contract through conda-forge, and a converter run
once must not enter it.

| script | what it does |
|---|---|
| `trajectory_to_bcif.py` | DCD or XTC → BinaryCIF, the trajectory format from 2026-08-19. Verifies the round trip by default. |

## `trajectory_to_bcif.py`

```
scripts/trajectory_to_bcif.py traj.xtc traj.bcif --top conf_ed.gro
scripts/trajectory_to_bcif.py lib.dcd lib.bcif
scripts/trajectory_to_bcif.py --all data/rotamer_library
```

Measured over the shipped corpus, every file verified exact on its grid:

| | before | after | |
|---|---|---|---|
| 95 `*.dcd` | 32.33 MB | **8.12 MB** | 4.0× |
| 1 `traj.xtc` | 12.44 MB | 9.87 MB | 1.26× |
| total | 44.78 MB | **17.99 MB** | 2.5× |

It exists because `python-ihm`'s `BinaryCifWriter` implements only ByteArray,
Delta, RunLength and the string/mask encoders. The compression here needs
**FixedPoint** and **IntegerPacking**, which it does not have. The reader side
needs nothing written: `ihm_format.c`, which IMP vendors and `libimp_atom`
exports, implements all seven encodings.

The reasoning, the precision it costs, and the two encoding traps are in
[`okf/validation/bcif_for_trajectories.md`](../okf/validation/bcif_for_trajectories.md).
