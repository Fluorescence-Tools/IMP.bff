## \example structure/drot_rotamer_library.py
# Building a rotamer library and using it for FRET, end to end.
#
# A dye rotamer library is an ensemble of conformers of one dye+linker with a
# weight each. `.drot` (PRD-118) stores that ensemble the way a chemist thinks
# about it: one template, a Z-matrix, and per conformer its own base
# coordinates, bond lengths, bond angles and dihedrals. Two consequences --
# the file is self-contained (atom names, residue names and elements ride
# inside, so no `.pdb` has to sit beside it), and a conformer *is* a dihedral
# vector, which is what continuous side-chain sampling needs.
#
# This example runs the whole path: trajectory in, library out, library used.
# On the command line the first half is one program:
#
#     imp_bff_traj2drot lib.bcif lib.drot.pto --top lib.pdb --weights lib_w.txt
#     imp_bff_traj2drot raw.dcd  lib.drot.pto --top lib.pdb --cluster 1.0
#     imp_bff_traj2drot --all data/rotamer_library
#
# with `--cluster` taking the raw-MD path: leader-clustering (the algorithm
# FRETpredict's own libraries were built with), leaders become the rotamers
# and the cluster populations the weights.

from pathlib import Path

import numpy as np

import IMP.bff
from IMP.bff import get_structure_dir

DATA = Path(IMP.bff.get_data_path("rotamer_library"))
STEM, CUTOFF = "A48_C1R", 30          # Alexa488 + C1R linker, 3 A cutoff set

# --- 1. the conformers: a trajectory and one weight per frame -------------
# Any trajectory IMP.bff reads works here (.bcif, .dcd, .xtc through mdtraj);
# the shipped BinaryCIF is used because it is what the package already has.
template = DATA / f"{STEM}.pdb"
names, elements, resnames = [], [], []
for line in template.read_text().splitlines():
    if line.startswith(("ATOM", "HETATM")):
        names.append(line[12:16].strip())
        resnames.append(line[17:20].strip())
        elements.append(line[76:78].strip() or names[-1][:1])

frames = np.asarray(IMP.bff.read_bcif_trajectory(
    str(DATA / f"{STEM}_cutoff{CUTOFF}.bcif"), len(names), "_rotamer_coord"),
    dtype=float).reshape(-1, len(names), 3)
weights = np.loadtxt(DATA / f"{STEM}_cutoff{CUTOFF}_weights.txt").reshape(-1)
print(f"{frames.shape[0]} conformers of {frames.shape[1]} atoms")

# --- 2. write the library --------------------------------------------------
# The default encoding is lossless (float32 internal coordinates). The compact
# rung -- `encoding.lossless = False` -- is a third the size at ~3e-3 A, which
# is fine for a picture and not fine for a fit: kappa^2 comes from transition
# dipole directions between atoms ~1.7 A apart, and those errors do not
# average away.
out = Path("A48_C1R_cutoff30.drot.pto")
IMP.bff.write_probe_rotamer_drot(str(out), np.ascontiguousarray(frames).ravel(),
                   names, elements, resnames, np.ascontiguousarray(weights),
                   IMP.bff.ProbeRotamerDrotEncoding())
print(f"wrote {out} ({out.stat().st_size} bytes)")

# --- 3. read it back and check it is the same ensemble ---------------------
library = IMP.bff.read_probe_rotamer_drot(str(out))
back = library.coords
print(f"round trip: {library.n_rotamers} rotamers, "
      f"max deviation {np.abs(back - frames).max():.2e} A, "
      f"names carried in the file: {list(library.atom_names)[:3]} ...")

# --- 4. use it: two dyes on hGBP1 and the FRET between them ----------------
# A path is a library name as far as the ensemble API is concerned; the dye's
# transition-dipole and attachment selectors come from the registry entry the
# stem names. A library of a dye that is *not* in the registry is passed as a
# loaded dict with its own `metadata` instead.
pdb = str(get_structure_dir("1DG3.pdb"))
donor = IMP.bff.ProbeRotamerEnsemble.from_site(pdb, "A", 481, str(out),
                                          position_name="A481")
acceptor = IMP.bff.ProbeRotamerEnsemble.from_site(
    pdb, "A", 496, "AlexaFluor 594 C1R cutoff30", position_name="A496")
print(f"donor {donor.n_rotamers} rotamers, Z={donor.partition:.3f}, "
      f"mean position {donor.mean_position.round(1)}")

eff = donor.pair_distribution_from_probes(acceptor, "AlexaFluor 488",
                                        "AlexaFluor 594")
print(f"R0 = {eff.forster_radius / 10:.2f} nm at <kappa2> = "
      f"{eff.kappa2_avg:.3f}; E_static {eff.static_efficiency:.3f}, "
      f"E_dynamic1 {eff.dynamic1:.3f}, E_dynamic2 {eff.dynamic2:.3f}")

# The same library through the registry name gives the same numbers -- the
# store changed, the ensemble did not.
shipped = IMP.bff.ProbeRotamerEnsemble.from_site(
    pdb, "A", 481, "AlexaFluor 488 C1R cutoff30", position_name="A481")
print(f"same as the shipped library: "
      f"{abs(shipped.partition - donor.partition) < 1e-9}")

out.unlink()
