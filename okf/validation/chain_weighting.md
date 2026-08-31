# Chain weighting: what the shipped table actually does to a dye AV

**Date: 2026-08-31.** `chain_weighting` was a field in this module's fps.json
schema, documented as *"Weight AV grid points by linker-chain statistics instead
of uniformly"*, that **nothing read**. There was no AV attribute for it and one
grep hit in the whole tree. A file that asked for it got the unweighted volume
and no complaint. It is implemented now, and this is what it does.

## What it is

An accessible volume treats every reachable voxel as equally likely. A real
linker does not: a random coil rarely sits fully collapsed against its
attachment. The correction is a tabulated \(P(\ell)\) over the **path length**
\(\ell\) -- the geodesic through free space that `PathMap` already computes per
voxel -- multiplied into the density of every reached voxel, so every quantity
read from the volume sees it: mean position, distances, efficiencies, the
exported grid.

`data/linker/chain_weighting.csv` is the table, from Olga
(`src/weighting_function.csv`, same author). Its header row carries 58 axis
keys from 79.2 to 293.8; a column is chosen by the **first key not less than the
linker length**, which is the reference implementation's rule.

## The finding: the shipped table does not cover a dye linker

A dye linker is 10-25 A. Every one of them selects the *first* column, key 79.2,
whose distribution peaks at a path length of **53 A**. Measured on the shipped
table for a 20.5 A linker:

| path length | weight |
|---|---|
| 0-10 A | 0.00000 |
| 20 A | 0.00022 |
| 40 A | 0.25527 |
| 53.5 A | 0.99943 (the peak) |

**0.02 % of the column's weight lies within a 20.5 A linker's reach.** So on a
T4L site the correction does not reweight the volume, it selects the volume's
outer shell and discards the rest: the mean position of the 132 AV moves
**13.7 A**.

That is not a refinement, and it should not pass for one. `AV::set_chain_weighting`
therefore measures `LinkerWeighting::get_supported_fraction(linker_length)` and
**warns** when it is below 1 %, naming the linker length and the fraction. The
capability is there, the table for a dye linker is not.

Two readings, and this note does not choose between them: either the column keys
are not linker lengths in angstrom and the reference's `lower_bound(linkerLength)`
is a unit mismatch for short linkers, or the table is for long tethers (PEG,
protein fusions) and a dye linker is out of its domain. Either way, a caller who
wants chain weighting on a dye AV should supply their own table through
`read_linker_weighting(path, linker_length)`.

## Where the weight is applied, and why not where it was first put

The weight multiplies the density **where the density is read** -- in
`get_xyz_density`, `get_xyz_density_soa` and `sync_tiles_from_soa` -- and is not
folded into `density_soa_` once per raster.

The first implementation did fold it, guarded by an "already applied" flag, and
it was wrong in a way that only one test caught: the **point cloud carried the
weighting and the exported grid did not**. `sync_tiles_from_soa` copies
densities into tiles exactly once per raster and returns early afterwards, and a
resample reaches that copy at a different point than it reaches the fold. A
second flag to track whether the tiles had seen the fold did not fix it either;
the interleaving across a multi-stage resample has more orders than the flags
had states.

Multiplying at the point of use has no state to get out of step and costs one
multiply per voxel per read. `PathMap::path_weight(i)` is the one place.

The test is `test_chain_weighting_reaches_the_exported_grid`: write the same
site with and without weighting, parse both OpenDX files, and assert the grids
differ and the weighted maximum is the smaller. It also found a second defect --
`write_av`'s `.dx` branch was reading `PathMap::get_value`, which is the
**obstacle raster** `SampledDensityMap` was built from, not the accessible
volume. Every `.dx` this module wrote for an AV was the wrong map.

## Two deliberate differences from the reference

* **Out-of-range path lengths are clamped**, not NaN. The reference's
  `TabulatedFunction::value` returns NaN outside `[xMin, xMax]`; one NaN weight
  makes the mean position, every distance and every efficiency computed from
  that volume NaN too.
* **A linker longer than every key gets the uniform weighting**, not the last
  column. Extrapolating a fitted distribution past its axis is a guess, and a
  silent one.

## Tests

`test/label/test_av_kernels_and_io.py` -- the uniform weighting is the
unweighted volume; a table interpolates linearly and clamps; a table needs two
points; the shipped table loads and peaks at 1.0; a linker longer than the table
gets no weighting; and the fps.json field reaches the volume, which is the bug
this closes.
