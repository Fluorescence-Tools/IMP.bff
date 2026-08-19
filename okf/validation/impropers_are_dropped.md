# The combined-system builder drops every improper

**Status:** closed 2026-08-19 — fixed in `2e08e06`. Kept because the reasoning is the record of a physics change.

## What

`cgdye.topology.build_dye_protein_system` builds a combined dye+protein
force-field system from two MOL2 files and two templates. It derives bonds,
angles and dihedrals from the MOL2 connectivity, and then:

```python
impropers = []
```

The list is never filled. Both shipped templates declare improper centres, and
`_build_impropers` — which exists, is correct, and is used by the *other*
builder, `build_dye_topology` — would expand them:

| component | improper centres declared | impropers that would be built |
|---|---|---|
| atto655 | 36 | 33 |
| CX4 | 48 | 48 |

## Why it matters

atto655's template marks its ring system as feature `ring_core` with
`rb NO, md_fixed NO` — the ring atoms are ordinary movable particles, not a
rigid body and not frozen. The declared ring impropers are therefore the *only*
term keeping the conjugated system planar. Without them a Langevin or MC run
has nothing opposing pyramidalisation of the chromophore, which is the part of
the dye whose geometry the transition dipole is defined by.

The consuming machinery is complete and has always been: both
`scoring.build_dye_restraints` and `cgdye.sim._build_restraints` turn impropers
into `IMP.core.DihedralRestraint`s with a harmonic about the current dihedral.
Only the production side is missing.

## What this hid

Because `impropers` is always empty:

* `scoring.build_dye_restraints`'s improper branch had never executed, and it
  read `t["k"]` on an `FFTorsionType` — a `TypeError` waiting for its first
  iteration. Fixed and gated in 192a763.
* The two exclusion derivations (`scoring.compute_exclusions`, which adds
  improper pairs, and `cgdye.sim._derive_exclusions`, which does not) agree on
  every real system — but only because the term they differ on is empty. They
  are not equivalent; they are untested against each other.
* `test/cgdye/test_dye_topology.py` exercises `build_dye_topology`'s improper
  path only with `"impropers": []`, so the expansion is unmeasured on both
  paths.

## What it turned out to be

Not a missing capability. `cgdye/topology.py` held **three** builders of a
`DyeForceFieldSystem`, and the largest — `build_system_from_specs`, the body of
the `build-system` command — *did* expand the templates' impropers, producing
exactly the 81 predicted above. It could not be used: it never returned a
system, it wrote a CIF and printed, so its output was reachable only by writing
a file and reading it back. Every caller used one of the two wrappers, and both
hard-coded `impropers = []`.

Verified before the change: on CX4+atto655 the two builders agreed on sites
(138), bonds (146), angles (262), dihedrals (207), components, all four group
maps, `fixed_groups` and all five type tables, and differed on nothing but
impropers.

Fixed in `2e08e06` by collapsing the three into one implementation that returns
the system. Systems now carry the impropers their templates declare — 81 for
CX4+atto655, 33 for atto655 alone. This moves simulation results, and was done
on instruction after the gap was raised.

Two consequences already visible: `scoring.build_dye_restraints`'s improper
branch now executes (it had a `TypeError` in it, fixed in `192a763` while it
was still unreachable), and `DyeForceFieldSystem::get_exclusions`'s
`include_impropers` flag is live rather than inert.
