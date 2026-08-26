"""
Coarse-grained protein potentials
=================================
The knowledge-based and steric terms a coarse-grained protein model is scored
with, on a real structure. They were ``numba`` kernels in ``IMP.cgmol`` and
classes in ChiSurf; here they are IMP scores and restraints, so they compose
with IMP's containers, scoring functions and optimizers like any other term.

What is what:

* ``IMP.bff.build_clash_restraint`` -- soft-sphere overlap, which is
  ``IMP.core.SoftSpherePairScore`` over an ``IMP.container.ClosePairContainer``
  with the bonded pairs filtered out;
* ``IMP.bff.GoRestraint`` -- a Gō model: a truncated Lennard-Jones well per
  residue pair, centred on the distance *this* structure has;
* ``IMP.bff.LennardJonesBeadPairScore`` -- one equilibrium distance for every
  bead pair, no types and no force field;
* ``IMP.bff.GeneralizedBornRestraint`` -- implicit solvation, reading charges
  from ``IMP.atom.Charged`` and radii from ``IMP.core.XYZR``;
* ``IMP.bff.residue_solvent_accessible_surface`` -- Shrake-Rupley with one
  sphere per residue;
* ``IMP.bff.MiyazawaJerniganPairScore`` and ``IMP.bff.UNRESCentroidPairScore``
  -- residue-typed contact potentials, which are
  ``IMP.core.StatisticalPairScore`` over a table.
"""

# %%
# First the imports, and a C-alpha representation of T4 lysozyme.
import json

import numpy as np
import pylab as plt

import IMP
import IMP.atom
import IMP.container
import IMP.core
import IMP.bff

model = IMP.Model()
hierarchy = IMP.atom.read_pdb(
    IMP.bff.get_example_path("structure/T4L/3GUN.pdb"), model,
    IMP.atom.CAlphaPDBSelector())
cas = [a.get_particle_index() for a in IMP.atom.get_leaves(hierarchy)]
print(f"{len(cas)} residues")

# %%
# A Gō model. The native contacts are taken from the structure as it stands, so
# the energy is at its floor here and rises as the fold is disturbed -- which
# is the point of the term.
go = IMP.bff.GoRestraint(model, cas, epsilon=1.0, cutoff=6.5,
                         nn_e_factor=0.1)
print(f"native contacts: {go.get_n_native()} of "
      f"{go.get_contacts().get_n_contacts()}")
print(f"Go energy of the native fold: {go.unprotected_evaluate(None):.2f}")

# %%
# Push one residue a long way and the term notices.
xyz = IMP.core.XYZ(model, cas[len(cas) // 2])
home = xyz.get_coordinates()
displaced = []
offsets = np.linspace(0.0, 8.0, 17)
for d in offsets:
    xyz.set_coordinates(home + IMP.algebra.Vector3D(d, 0.0, 0.0))
    displaced.append(go.unprotected_evaluate(None))
xyz.set_coordinates(home)

# %%
# The steric term over the same beads. ``clash_tolerance`` is the softness: the
# ported potential is ``((sigma - r)/t)^2`` and IMP's soft sphere is
# ``k(sigma - r)^2/2``, so the restraint uses ``k = 2/t^2``.
for a in IMP.atom.get_leaves(hierarchy):
    IMP.core.XYZR(a).set_radius(2.5)
clash = IMP.bff.build_clash_restraint(hierarchy, clash_tolerance=2.0)
print(f"clash energy: "
      f"{IMP.core.RestraintsScoringFunction([clash]).evaluate(False):.2f}")

# %%
# The solvent-accessible surface, one sphere per residue. A buried residue
# contributes nothing and an exposed one contributes most of its sphere.
area = IMP.bff.residue_solvent_accessible_surface(
    hierarchy, IMP.atom.AT_CA, n_sphere=590, probe=1.0, radius=2.5)
print(f"accessible surface: {area:.0f} A^2 over {len(cas)} residues")

# %%
# Per residue, so the exposed stretches are visible. ``residue_asa`` takes the
# representative coordinates directly, which is what the restraint does
# internally.
coords = np.array([IMP.core.XYZ(model, p).get_coordinates() for p in cas])
per_residue = np.array([
    IMP.bff.residue_asa(coords.ravel(), 590, 1.0, 2.5)
    - IMP.bff.residue_asa(np.delete(coords, i, axis=0).ravel(), 590, 1.0, 2.5)
    for i in range(len(cas))])

# %%
# And the two plots: what the Gō term says about a displacement, and where the
# surface is.
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].plot(offsets, displaced, "o-")
ax[0].set_xlabel("displacement of one residue / $\\AA$")
ax[0].set_ylabel("Gō energy")
ax[0].set_title("a native contact map resists being broken")

ax[1].plot(per_residue, lw=1)
ax[1].set_xlabel("residue")
ax[1].set_ylabel("accessible surface / $\\AA^2$")
ax[1].set_title("solvent exposure along the sequence")
plt.tight_layout()
plt.show()

# %%
# The residue-typed contact potentials. Every parameter table this module ships
# is in one container, ``data/potentials.pto``, and the default constructors
# read it -- so a caller says what it wants scored, not where the numbers live.
print(IMP.bff.potential_table_names())

# %%
# ``add_residue_type_score_data`` puts the residue type on one atom per residue
# -- the C-beta for Miyazawa-Jernigan -- and the score reads it through
# ``IMP.bff.get_residue_type_key()``, exactly as ``IMP.atom.DopePairScore``
# reads the type ``IMP.atom.add_dope_score_data`` writes. A residue without
# that atom (a glycine has no C-beta) takes no part.
m2 = IMP.Model()
all_atom = IMP.atom.read_pdb(
    IMP.bff.get_example_path("structure/T4L/3GUN.pdb"), m2,
    IMP.atom.NonWaterNonHydrogenPDBSelector())
typed = IMP.bff.add_residue_type_score_data(all_atom, IMP.atom.AT_CB)
lsc = IMP.container.ListSingletonContainer(
    m2, [p.get_index() for p in typed])
print(f"{len(typed)} residues carry a type")

mj = IMP.container.PairsRestraint(
    IMP.bff.MiyazawaJerniganPairScore(),
    IMP.container.ClosePairContainer(lsc, 6.5, 0.0))
unres = IMP.container.PairsRestraint(
    IMP.bff.UNRESCentroidPairScore(),
    IMP.container.ClosePairContainer(lsc, 15.0, 0.0))
print(f"Miyazawa-Jernigan: "
      f"{IMP.core.RestraintsScoringFunction([mj]).evaluate(False):.1f}")
print(f"UNRES centroid:    "
      f"{IMP.core.RestraintsScoringFunction([unres]).evaluate(False):.1f} "
      f"kcal/mol")

# %%
# The Ramachandran map, three channels of it -- general, proline, glycine --
# read as :math:`-\log(P/P_{max})`, so a residue at the peak of the
# distribution costs nothing.
table = IMP.bff.read_potential_table("ramachandran")
rama = IMP.bff.RamachandranRestraint(m2, all_atom, table.table,
                                     table.shape[0], table.shape[1])
print(f"Ramachandran: {rama.unprotected_evaluate(None):.1f} over "
      f"{len(typed)} residues")

# %%
# The hydrogen-bond term reads four distances per bond -- O-H, O-N, C-H, C-N --
# out of a four-channel lookup, so it needs the amide **hydrogens**. A crystal
# structure read without them scores zero, which is not a defect: there is
# nothing to donate.
hb_table = IMP.bff.read_potential_table("hbond")
hbond = IMP.bff.HydrogenBondRestraint(m2, all_atom, hb_table.table,
                                      hb_table.shape[1])
print(f"hydrogen bonds found: {hbond.get_n_hbonds()}")

# %%
# What the container carries, and what was done to each table on the way in --
# the part four loose ``.npy`` files could not say.
manifest = json.loads(IMP.bff.read_potential_manifest())
for name, entry in manifest["tables"].items():
    print(f"{name:14s} {entry['kind']:9s} {entry['source']}")
    print(f"               {entry['note']}")
