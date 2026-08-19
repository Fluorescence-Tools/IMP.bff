# The combined-system builder drops every improper

**Status:** open — a physics gap, not yet changed. Found 2026-08-19.

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

## The decision that is open

Filling `impropers` changes what the force field restrains, so it moves
simulation results. That is a physics change and is deliberately not being made
as part of the structural work. What is needed is a judgement on whether the
omission was intentional — a decision that the ring core is stiff enough
without it — or an oversight when the combined builder was written alongside
the per-dye one.
