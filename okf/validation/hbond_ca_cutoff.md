---
type: validation
title: "The hydrogen-bond potential's C-alpha cutoff never engaged"
description: The ported H-bond kernel compared a plain C-alpha distance against the square of its cutoff, so its 8 A prefilter passed every residue pair. Applied properly in imp.bff's HydrogenBondRestraint; the counts and energies move for extended geometries, and the difference is pinned.
resource: /Users/tpeulen/dev/imp.bff
tags: [validation, imp.bff, potentials, hydrogen-bond, chisurf, imp-tricks, coarse-grained]
timestamp: '2026-08-26T00:00:00Z'
---

# The hydrogen-bond potential's C-alpha cutoff never engaged

`IMP.cgmol.statpot._kernels._hbond_kernel` (imp-tricks), which chisurf's
`HPotential` calls, takes a C-alpha distance matrix and a **squared** cutoff:

```python
nHbond, Ehbond = _hbond_kernel(s1.l_res, s1.dist_ca, s1.xyz, self._hPot,
                               cca2, ch2)          # cca2 = cutoff_ca ** 2
...
    for rj in range(ri + 1, n_res):
        if ca_dist[ri, rj] > cutoff_ca2:
            continue
```

`dist_ca` does not hold squared distances. `atom_dist`
(`chisurf/core/structure/protein.py:81`) writes

```python
d12 = math.sqrt(b1*b1 + b2*b2 + b3*b3)
aDist[i, j] = d12
```

so the comparison is a plain distance against 64 A^2. A protein whose
C-alphas are all within 64 A of each other -- which is every protein under
~120 kDa -- passes the prefilter entirely, and the 8 A cutoff the parameter
advertises is never applied.

## Where it matters and where it does not

The same mismatch is in `_mj_kernel` (`ca_dist[i, j] > cut_ca2`) and
`_asa_kernel` (`ca_dist[i, j] < cutoff2`), and **there it costs only time**: a
pair the prefilter should have dropped still fails the real test -- a C-beta
contact needs the C-alphas within `cutoff + 3.2 A`, and a sphere at 2.5 A
cannot be occluded by a neighbour more than 6 A away. The prefilters are
optimisations in those two.

In the hydrogen-bond kernel the C-alpha cutoff is **not** an optimisation. A
carbonyl O and an amide H can be 3 A apart with their C-alphas much further,
in an extended or distorted geometry, so the cutoff decides which of those
count. Passing everything counts them all.

## What the port does

`IMP::bff::HydrogenBondRestraint` applies it: a plain distance against a plain
cutoff. On a compact structure nothing changes -- every bonded pair's C-alphas
are inside 8 A -- and on the test structure the two agree:

| cutoff_ca | bonds | energy |
|---|---|---|
| 8 A (the default, applied) | 3 | -0.1324432499910165 |
| 8 A (the Python, not applied) | 3 | -0.1324432499910165 |
| 4 A (applied) | 2 | 0.32062226315756825 |

The third row is the proof that it engages;
`test/potentials/test_potentials.py::test_the_c_alpha_cutoff_is_applied` pins
both.

## The UNRES table is filled in one triangle, and was read in both

`unres.npy` is `(20, 20, 400)` -- two residue types and a distance -- and it is
**not symmetric**. Past the short-range repulsion the entries with
`type_i < type_j` average -0.129 and their transposes -0.003: one half carries
the potential and the other is all but empty.

`centroid2` indexes it with `potential[residue_type_i, residue_type_j, bin]`
where `i` and `j` are the two residues' positions **in the chain**, not their
type order. So which half a contact reads depends on which of the two residues
comes first in the sequence, and about one contact in two reads the empty half
and scores zero.

Measured on T4 lysozyme, over the same 4 983 scored pairs:

| how the pair is indexed | energy |
|---|---|
| upper triangle (one value per unordered pair) | **-371.031** |
| the residues' order in the chain (the Python) | -207.198 |

A pair potential has one value per *unordered* pair. `data/potentials.pto`
stores the filled half, and `IMP.core.StatisticalPairScore` orders every pair
before it looks -- so a contact scores the same whichever residue the file
happens to reach first. The manifest records `"triangle": "upper"`, and
`test_the_unres_table_is_one_value_per_unordered_pair` holds it there.

## Two more, while reading the same kernels

- **`potentials.py:go()` carried its accumulator across pairs.** `tmp -= sr`
  runs on whatever the previous pair left in `tmp`, so every pair after the
  first is wrong by the sum before it. `_go_kernel` in imp-tricks fixed this
  with a local; imp.bff ports that one, and `GoRestraint` on a native fold
  returns exactly minus the sum of the well depths, which is what a Go model
  at its own contacts must give.
- **`_asa_kernel` reads its representative atom from column 6** of the residue
  lookup table, documented as "representative atom (C-alpha / centroid)".
  Column 6 is `CG` in the layout it is fed
  (`chisurf/core/structure/protein.py:internal_atom_numbers`), which only
  proline and a few other residues carry -- so the area it returned was of
  those residues alone. In imp.bff the atom is an argument
  (`IMP::atom::AT_CA` by default).
