"""Reading a PDB into mol2, and writing mol2 out.

The library could **parse** mol2 (:func:`IMP.bff.parse_dye_mol2`,
``io/cif.py``) but not write one, and mol2 is the format a dye topology has to
be in before a force field can be built for it. The writer lived in a CLI
driver, ``cgdye/scripts/pdb_to_mol2.py``, which is why it was invisible: a
capability reachable only by running a script is a capability nobody finds.

Lifted into the library in the PRD-113 cleanup (2026-08-18) when that script
moved to ``junk/``. The driver is gone; this is what it was for.

It landed in ``cgdye/io/`` first and moved here when ``cgdye`` left the domain
layout -- a capability rescued from a script should not then be buried in a
legacy package.

Bonds come from ``CONECT`` records when the PDB has them and from geometry when
it does not -- :func:`infer_bonds` uses a covalent-radius cutoff, which is a
guess and is why it is not the default.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

__all__ = [
    "parse_pdb_atoms",
    "parse_conect_bonds",
    "infer_bonds",
    "write_mol2",
]

_ELEM_TO_TRIPOS = {
    "C": "C.3",
    "N": "N.3",
    "O": "O.3",
    "S": "S.3",
    "H": "H",
    "P": "P.3",
    "F": "F",
    "Cl": "Cl",
    "Br": "Br",
    "I": "I",
}


def parse_pdb_atoms(path: Path) -> dict:
    atoms = {}
    with open(path) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            serial = int(line[6:11])
            name = line[12:16].strip()
            resname = line[17:20].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            elem_raw = line[76:78].strip() if len(line) > 76 else ""
            if not elem_raw:
                m = re.match(r"[A-Za-z]", name)
                elem_raw = m.group(0) if m else "C"
            atoms[serial] = dict(
                serial=serial,
                name=name,
                resname=resname,
                x=x,
                y=y,
                z=z,
                elem=elem_raw[0].upper(),
            )
    return atoms


def parse_conect_bonds(path: Path) -> set:
    bonds = set()
    with open(path) as fh:
        for line in fh:
            if line[:6].strip() != "CONECT":
                continue
            fields = line.split()
            if len(fields) < 3:
                continue
            a = int(fields[1])
            for tok in fields[2:]:
                b = int(tok)
                if a != b:
                    bonds.add(tuple(sorted((a, b))))
    return bonds


def infer_bonds(atoms: dict, scale: float = 1.22) -> set:
    """Geometry-based bond inference using covalent radii."""
    cov_r = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "P": 1.07}
    serials = sorted(atoms)
    bonds = set()
    for i, a in enumerate(serials):
        ea = atoms[a]["elem"]
        ra = cov_r.get(ea, 0.77)
        for b in serials[i + 1 :]:
            eb = atoms[b]["elem"]
            rb = cov_r.get(eb, 0.77)
            dx = atoms[a]["x"] - atoms[b]["x"]
            dy = atoms[a]["y"] - atoms[b]["y"]
            dz = atoms[a]["z"] - atoms[b]["z"]
            d = math.sqrt(dx * dx + dy * dy + dz * dz)
            if d <= scale * (ra + rb):
                bonds.add(tuple(sorted((a, b))))
    return bonds


def write_mol2(path: Path, atoms: dict, bonds: set, mol_name: str = "MOL") -> None:
    serials = sorted(atoms)
    idx = {s: i + 1 for i, s in enumerate(serials)}  # 1-based atom indices

    lines = []
    lines.append("@<TRIPOS>MOLECULE")
    lines.append(mol_name)
    lines.append(f"{len(atoms)} {len(bonds)} 0 0 0")
    lines.append("SMALL")
    lines.append("NO_CHARGES")
    lines.append("")

    lines.append("@<TRIPOS>ATOM")
    for s in serials:
        a = atoms[s]
        el = a["elem"]
        atype = _ELEM_TO_TRIPOS.get(el, el)
        lines.append(
            f"{idx[s]:6d} {a['name']:<8s}"
            f" {a['x']:10.4f} {a['y']:10.4f} {a['z']:10.4f}"
            f" {atype:<8s} 1 {a['resname']} 0.0000"
        )

    lines.append("@<TRIPOS>BOND")
    for bi, (a, b) in enumerate(sorted(bonds), 1):
        lines.append(f"{bi:6d} {idx[a]:6d} {idx[b]:6d} 1")

    lines.append("")
    path.write_text("\n".join(lines))
