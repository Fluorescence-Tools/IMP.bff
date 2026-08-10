#!/usr/bin/env python
"""Convert PDB files to TRIPOS MOL2 format for IMP bond-topology ingestion.

For molecules WITH CONECT records (atto655, atto655amin): bonds come from CONECT.
For molecules WITHOUT CONECT (cx4): bonds are inferred from geometry (covalent radii).

Usage:
    python scripts/pdb_to_mol2.py inputs/structures/atto655.pdb
    python scripts/pdb_to_mol2.py inputs/structures/cx4.pdb --infer-bonds
    python scripts/pdb_to_mol2.py --all          # converts all PDBs under inputs/structures/
"""

import math
import re
import sys
from pathlib import Path

import click

# ---------------------------------------------------------------------------
# PDB parsing helpers
# ---------------------------------------------------------------------------


def _parse_atoms(path: Path) -> dict:
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


def _parse_conect(path: Path) -> set:
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


def _infer_bonds(atoms: dict, scale: float = 1.22) -> set:
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


# ---------------------------------------------------------------------------
# MOL2 writer
# ---------------------------------------------------------------------------

# Minimal TRIPOS atom-type mapping by element
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


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("pdb_files", nargs=-1, type=click.Path(exists=True, path_type=Path))
@click.option(
    "--all",
    "convert_all",
    is_flag=True,
    default=False,
    help="Convert all *.pdb files under inputs/structures/.",
)
@click.option(
    "--infer-bonds",
    is_flag=True,
    default=False,
    help="Use geometry-based bond inference (for PDBs without CONECT records).",
)
@click.option(
    "--out-dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Output directory (default: same as input PDB).",
)
def main(pdb_files, convert_all, infer_bonds, out_dir):
    root = Path(__file__).resolve().parents[1]

    paths = list(pdb_files)
    if convert_all:
        paths = sorted((root / "inputs" / "structures").glob("*.pdb"))

    if not paths:
        click.echo("No PDB files specified. Use --all or pass filenames.", err=True)
        sys.exit(1)

    for pdb in paths:
        pdb = Path(pdb)
        atoms = _parse_atoms(pdb)
        bonds = _parse_conect(pdb)

        if not bonds or infer_bonds:
            if not bonds:
                print(f"{pdb.name}: no CONECT records — inferring bonds from geometry")
            else:
                print(
                    f"{pdb.name}: --infer-bonds requested — overriding CONECT with geometry"
                )
            bonds = _infer_bonds(atoms)

        mol_name = pdb.stem
        dest_dir = out_dir if out_dir else pdb.parent
        dest_dir = Path(dest_dir)
        dest_dir.mkdir(parents=True, exist_ok=True)
        mol2_path = dest_dir / (mol_name + ".mol2")

        write_mol2(mol2_path, atoms, bonds, mol_name=mol_name)
        print(f"  {pdb.name} -> {mol2_path}  ({len(atoms)} atoms, {len(bonds)} bonds)")


if __name__ == "__main__":
    main()
