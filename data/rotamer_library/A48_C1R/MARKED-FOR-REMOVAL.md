# `traj.xtc` and the sibling `*.dcd` libraries are marked for removal

Superseded by BinaryCIF, which is smaller than either and readable by the C
parser IMP already vendors, with no new dependency:

| format | bytes/coordinate |
|---|---|
| DCD (95 files, 32.33 MB) | 4.31 |
| XTC (this file, 12.44 MB) | 1.60 |
| BinaryCIF at a 0.1 Å grid | **1.27** |

Measured, decoded exactly, and written up in
[`okf/validation/bcif_for_trajectories.md`](../../../okf/validation/bcif_for_trajectories.md).

**Nothing is converted or deleted yet.** These files are still the shipped data
and are still read by `IMP.bff.io.structure.read_dcd`. Removal waits on a
BinaryCIF encoder — `python-ihm`'s writer cannot produce the needed
FixedPoint + IntegerPacking chain — and on converting the corpus.
