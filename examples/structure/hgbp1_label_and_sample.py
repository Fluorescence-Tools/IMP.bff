## \example structure/hgbp1_label_and_sample.py
# Label hGBP1 with a coarse-grained dye and sample the dye's conformations.
#
# This is the end-to-end cgdye workflow on a real system: take a structure,
# attach a dye at a named site, and explore where the dye can actually go. The
# same three steps underlie every accessible-volume or FRET-efficiency number
# IMP.bff produces, so it is worth seeing them separately before they disappear
# into a scoring function.
#
# The system is hGBP1 (PDB 1DG3) labelled at chain A residue 481 with
# Alexa488 C5-maleimide (the `alexa488_r48` template) -- a site used throughout
# the cgdye tests, so the numbers here can be compared against them.
#
# Three things this example is careful about, each of which has caused a real
# bug:
#
# - **The structure and dye templates are module data**, reached through
#   `get_structure_dir()` rather than by walking up from `__file__`. Deriving
#   data paths from the source layout broke the moment cgdye moved between
#   repositories.
# - **Attachment mutates the dye in place.** `attach_dyes` transforms the dye's
#   coordinates onto the site frame; it does not return a new molecule. The
#   check below is on the dye having *moved*, which is the observable.
# - **The site's sidechain is stripped.** A dye occupies the space the
#   sidechain did, so leaving both produces a steric clash that quietly
#   distorts every subsequent sample.

import sys

import IMP
import IMP.atom
import IMP.core
import IMP.algebra

from IMP.bff.cgdye.utils import get_structure_dir
from IMP.bff.cgdye.labeling.attachment import attach_dyes, resolve_site

# sys.argv, not []: IMP.setup_from_argv reads argv[0] without checking, so an
# empty list segfaults the interpreter rather than raising. Reported as a trap
# in the OKF known-issues.
IMP.setup_from_argv(
    sys.argv, "Label hGBP1 at residue 481 with a CG dye and sample the dye")

# --- 1. the structure and the dye -------------------------------------------
model = IMP.Model()
protein = IMP.atom.read_pdb(
    str(get_structure_dir("1DG3.pdb")), model, IMP.atom.NonWaterPDBSelector())
dye = IMP.atom.read_mol2(str(get_structure_dir("alexa488_r48.mol2")), model)

n_protein_atoms = len(IMP.atom.get_by_type(protein, IMP.atom.ATOM_TYPE))
n_dye_atoms = len(IMP.atom.get_by_type(dye, IMP.atom.ATOM_TYPE))
print("hGBP1 atoms:", n_protein_atoms, " dye atoms:", n_dye_atoms)

# --- 2. resolve the labelling site ------------------------------------------
# resolve_site returns the backbone frame the dye is attached against. If the
# residue is missing from the model this raises rather than silently labelling
# somewhere else, which is the behaviour you want from a labelling step.
site = resolve_site(protein, "A", 481)
print("site 481 backbone atoms:", sorted(site))
assert {"N", "CA", "C"} <= set(site)

# --- 3. attach the dye ------------------------------------------------------
before = [IMP.core.XYZ(a).get_coordinates()
          for a in IMP.atom.get_by_type(dye, IMP.atom.ATOM_TYPE)]

attached = attach_dyes(protein, [(dye, "A", 481)], strip_site_sidechain=True)
print("attached:", len(attached), "dye(s) at residue", attached[0]["resnum"])

after = [IMP.core.XYZ(a).get_coordinates()
         for a in IMP.atom.get_by_type(dye, IMP.atom.ATOM_TYPE)]

moved = max(IMP.algebra.get_distance(b, a) for b, a in zip(before, after))
print("largest atom displacement on attachment: %.2f A" % moved)
assert moved > 1e-8, "attach_dyes transforms the dye in place; it did not move"

# The sidechain at the labelled site is gone: a dye and a sidechain cannot both
# occupy that volume, and leaving both produces clashes that distort sampling.
backbone = {"N", "CA", "C", "O", "OXT"}
site_residue_atoms = [
    IMP.atom.Atom(a).get_atom_type().get_string().strip().upper()
    for a in IMP.atom.get_by_type(protein, IMP.atom.ATOM_TYPE)
    if IMP.atom.Residue.get_is_setup(a.get_parent())
    and IMP.atom.Residue(a.get_parent()).get_index() == 481
]
print("atoms left at site 481:", sorted(set(site_residue_atoms)))
assert not (set(site_residue_atoms) - backbone), "sidechain was not stripped"

# --- 4. where can the dye go? -----------------------------------------------
# The dye's reachable volume is what an AV or a FRET efficiency is computed
# from. Here we simply report the labelled system; the sampling scripts under
# IMP.bff.cgdye.scripts (langevin_hgbp1_site481, rrt_hgbp1_site481) walk the
# dye's internal degrees of freedom and write a trajectory, and
# IMP.bff.cgdye.sampling holds the samplers they drive.
dye_atoms = IMP.atom.get_by_type(dye, IMP.atom.ATOM_TYPE)
centroid = IMP.algebra.get_centroid(
    [IMP.core.XYZ(a).get_coordinates() for a in dye_atoms])
print("dye centroid after attachment: "
      "(%.2f, %.2f, %.2f)" % (centroid[0], centroid[1], centroid[2]))

print("labelled hGBP1 at A/481 and verified the attachment")
