#!/usr/bin/env python

import os

import click
import IMP
import IMP.atom

from IMP.bff.cgdye.labeling.attachment import attach_dyes


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--protein-pdb",
    default="inputs/structures/1DG3.pdb",
    show_default=True,
    help="Protein PDB path",
)
@click.option(
    "--dye-mol2",
    default="inputs/structures/alexa488_r48.mol2",
    show_default=True,
    help="Dye MOL2 path",
)
@click.option("--chain", default="A", show_default=True, help="Attachment chain ID")
@click.option(
    "--residue",
    default=481,
    show_default=True,
    type=int,
    help="Attachment residue number",
)
@click.option(
    "--output-pdb",
    default="output/test_systems/hGBP1_1DG3_alexa488_r48_site481.pdb",
    show_default=True,
    help="Output attached-system PDB",
)
def main(protein_pdb, dye_mol2, chain, residue, output_pdb):
    model = IMP.Model()
    protein = IMP.atom.read_pdb(protein_pdb, model, IMP.atom.NonWaterPDBSelector())
    dye = IMP.atom.read_mol2(dye_mol2, model)

    attach_dyes(protein, [(dye, chain, residue)], strip_site_sidechain=True)

    root = IMP.atom.Hierarchy.setup_particle(IMP.Particle(model, "root"))
    root.add_child(protein)
    root.add_child(dye)

    out_dir = os.path.dirname(os.path.abspath(output_pdb))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    IMP.atom.write_pdb(root, output_pdb)
    click.echo(f"Wrote attached test system: {output_pdb}")


if __name__ == "__main__":
    main()
