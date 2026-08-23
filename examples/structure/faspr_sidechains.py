## \example bff/structure/faspr_sidechains.py
"""
Repacking protein side chains with the FASPR port
=================================================

``IMP.bff.faspr_pack`` is a 1:1 C++ port of FASPR (Huang, Pearce & Zhang,
Bioinformatics 2020; `source <https://github.com/tommyhuangthu/FASPR>`__,
MIT): it rebuilds every rotatable side chain of a protein backbone on the
Dunbrack-2010 backbone-dependent rotamer library and searches the global
minimum-energy assignment (DEE elimination + tree decomposition) under
FASPR's energy model -- CHARMM19 vdW, hydrogen bonds, disulfides and the
rotamer prior.

Why a modelling package wants this: an explicit dye (``cgdye``) is attached
where a side chain used to be, and the *neighbouring* side chains shape the
steric and quenching environment the dye feels. Packing them properly --
instead of trusting crystal coordinates that mutation or modelling may have
invalidated -- is part of PRD-118's accurate dye-quencher geometry.

The Dunbrack library is not redistributed with IMP.bff; this example finds
FASPR's ``dun2010bbdep.bin`` through ``$IMP_BFF_FASPR_ROT_LIB`` or the
``junk/FASPR`` clone and exits with a pointer otherwise.
"""

import os
import sys
from pathlib import Path

import numpy as np

import IMP.bff


def find_rotamer_library():
    env = os.environ.get("IMP_BFF_FASPR_ROT_LIB")
    if env and Path(env).exists():
        return Path(env)
    clone = Path(__file__).resolve().parents[2] / "junk" / "FASPR" / \
        "dun2010bbdep.bin"
    if clone.exists():
        return clone
    return None


def read_atoms(path):
    atoms = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith("ATOM"):
            key = (line[21], int(line[22:26]), line[12:16].strip())
            atoms[key] = (float(line[30:38]), float(line[38:46]),
                          float(line[46:54]))
    return atoms


def dihedral(p0, p1, p2, p3):
    b0 = np.asarray(p0) - np.asarray(p1)
    b1 = np.asarray(p2) - np.asarray(p1)
    b2 = np.asarray(p3) - np.asarray(p2)
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def main():
    rotlib = find_rotamer_library()
    if rotlib is None:
        sys.exit("dun2010bbdep.bin not found -- set IMP_BFF_FASPR_ROT_LIB "
                 "or clone FASPR into junk/FASPR")

    here = Path(__file__).resolve().parent
    backbone = here / "T4L" / "3GUN.pdb"
    repacked = here / "T4L" / "3GUN_faspr_port.pdb"

    IMP.bff.faspr_pack(str(backbone), str(repacked), str(rotlib))
    print(f"repacked {backbone.name} -> {repacked.name}")

    # chi1 before/after for a few buried and surface residues
    before, after = read_atoms(backbone), read_atoms(repacked)
    for chain, resid in [("A", 3), ("A", 22), ("A", 88), ("A", 119)]:
        def chi1(atoms, cg):
            n, ca, cb = (atoms[(chain, resid, a)] for a in ("N", "CA", "CB"))
            return dihedral(n, ca, cb, atoms[(chain, resid, cg)])
        for cg in ("CG", "CG1", "OG", "OG1", "SG"):
            if (chain, resid, cg) in before and (chain, resid, cg) in after:
                print(f"  {chain}:{resid:4d} chi1  crystal "
                      f"{chi1(before, cg):7.1f}   repacked {chi1(after, cg):7.1f}"
                      " deg")
                break


if __name__ == "__main__":
    main()