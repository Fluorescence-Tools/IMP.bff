#!/usr/bin/env python3
"""
Convert dye PDB files to mmCIF structure files.

Task 0.2: Use IMP.atom.read_pdb() to get coordinates/topology, 
then write out standard _atom_site records.
"""

import argparse
import os
import re
import sys
from pathlib import Path

import IMP
import IMP.atom
import ihm.format


def _element_from_atom_name(atom_name):
    # Simple heuristic for PDB atom names
    m = re.match(r"([A-Za-z]+)", atom_name)
    if not m:
        return "C"
    # Take first char of name part, e.g. 'CA' -> 'C', 'NZ' -> 'N'
    # NOTE: This is sufficient for the ~100 atom small molecules in dyes
    return m.group(1)[0].upper()


def convert_pdb_to_cif(pdb_path, cif_path, dye_id=None):
    """Convert a dye PDB to mmCIF with _atom_site records."""
    model = IMP.Model()
    # Read PDB into IMP
    # Use AllPDBSelector because dyes often have HETATM or unusual resnames
    hier = IMP.atom.read_pdb(pdb_path, model, IMP.atom.AllPDBSelector())

    if dye_id is None:
        dye_id = Path(pdb_path).stem

    atoms = IMP.atom.get_by_type(hier, IMP.atom.ATOM_TYPE)

    with open(cif_path, "w") as fh:
        writer = ihm.format.CifWriter(fh)
        writer.start_block(dye_id)

        # Basic atom site loop
        cols = [
            "group_PDB",
            "id",
            "type_symbol",
            "label_atom_id",
            "label_comp_id",
            "label_asym_id",
            "label_entity_id",
            "label_seq_id",
            "Cartn_x",
            "Cartn_y",
            "Cartn_z",
            "occupancy",
            "B_iso_or_equiv",
        ]

        with writer.loop("_atom_site", cols) as l:
            for i, a in enumerate(atoms, start=1):
                pos = IMP.core.XYZ(a).get_coordinates()
                name = IMP.atom.Atom(a).get_name()
                # Find parent residue for comp_id
                res = IMP.atom.Residue(a.get_parent())
                comp_id = res.get_name() if res.get_is_setup(res) else "DYE"

                l.write(
                    group_PDB="HETATM",
                    id=i,
                    type_symbol=_element_from_atom_name(name),
                    label_atom_id=name,
                    label_comp_id=comp_id,
                    label_asym_id="A",
                    label_entity_id=1,
                    label_seq_id=1,
                    Cartn_x=pos[0],
                    Cartn_y=pos[1],
                    Cartn_z=pos[2],
                    occupancy=1.0,
                    B_iso_or_equiv=0.0,
                )


def main():
    parser = argparse.ArgumentParser(description="Convert dye PDB to mmCIF")
    parser.add_argument("pdb", help="Input PDB path")
    parser.add_argument("cif", help="Output CIF path")
    parser.add_argument("--id", help="Dye ID for CIF block")

    args = parser.parse_args()

    if not os.path.exists(args.pdb):
        print(f"Error: {args.pdb} not found")
        sys.exit(1)

    convert_pdb_to_cif(args.pdb, args.cif, args.id)
    print(f"Converted {args.pdb} -> {args.cif}")


if __name__ == "__main__":
    main()
