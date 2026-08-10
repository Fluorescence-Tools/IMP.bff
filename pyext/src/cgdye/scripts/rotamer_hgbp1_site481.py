#!/usr/bin/env python

import os
import random
import sys
from pathlib import Path

import click
import IMP
import IMP.atom
import IMP.core
import IMP.rmf
import RMF



from IMP.bff.cgdye.io.rotamer_rmf import read_rotamer_library_rmf
from IMP.bff.cgdye.labeling.attachment import attach_dyes, place_dye_from_coords
from IMP.bff.cgdye.sampling.rotamer import (
from IMP.bff.cgdye.utils import get_structure_dir, get_template_dir, _data_root

def _structure(name):
    """A bundled input structure, wherever IMP.bff is installed.

    These paths used to be relative to the working directory, so a script only
    ran from one place -- and stopped running at all once the data became IMP
    module data under data/cgdye.
    """
    from IMP.bff.cgdye.utils import get_structure_dir
    return str(get_structure_dir(name))

    apply_rotamer_coords,
    sample_rotamer_index,
)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--protein-pdb", default=_structure("1DG3.pdb"), show_default=True)
@click.option("--chain", default="A", show_default=True)
@click.option("--residue", default=481, show_default=True, type=int)
@click.option("--dye-name", default="A48_C1R", show_default=True)
@click.option("--rotamer-lib-dir", default="cgdye/templates/rotamer", show_default=True)
@click.option("--n-samples", default=200, show_default=True, type=int)
@click.option("--seed", default=481, show_default=True, type=int)
@click.option(
    "--output-rmf",
    default="output/test_runs/hgbp1_1dg3_site481/rotamer/trajectory_rmf_lib.rmf3",
    show_default=True,
)
def main(
    protein_pdb,
    chain,
    residue,
    dye_name,
    rotamer_lib_dir,
    n_samples,
    seed,
    output_rmf,
):
    """Rotamer library sampling with clash detection (RMF library)."""
    rng = random.Random(seed)
    
    lib_path = os.path.join(rotamer_lib_dir, f"{dye_name}.rmf3")
    if not os.path.exists(lib_path):
        click.echo(f"Error: Rotamer library not found at {lib_path}")
        return
        
    lib = read_rotamer_library_rmf(lib_path)

    model = IMP.Model()
    protein = IMP.atom.read_pdb(protein_pdb, model, IMP.atom.NonWaterPDBSelector())
    
    # Use one of the structure files as a template for attachment.
    dye_dir = get_structure_dir()
    dye_mol2 = dye_dir / f"{dye_name}.mol2"
    if not dye_mol2.exists():
        # fallback for common naming
        if "A48" in dye_name or "alexa488" in dye_name:
            dye_mol2 = dye_dir / "alexa488_r48.mol2"
        elif "cx4" in dye_name:
            dye_mol2 = dye_dir / "cx4.mol2"
        elif "atto655" in dye_name:
            dye_mol2 = dye_dir / "atto655.mol2"
            
    dye = IMP.atom.read_mol2(str(dye_mol2), model)

    attached = attach_dyes(protein, [(dye, chain, residue)], strip_site_sidechain=True)
    site = attached[0]["site"]
    ca = IMP.core.XYZ(site["CA"]).get_coordinates()
    n = IMP.core.XYZ(site["N"]).get_coordinates()
    c = IMP.core.XYZ(site["C"]).get_coordinates()

    # Ensure bonds are present for visualization
    IMP.atom.add_bonds(protein)

    out_dir = os.path.dirname(os.path.abspath(output_rmf))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    root = IMP.atom.Hierarchy.setup_particle(IMP.Particle(model, "root"))
    root.add_child(protein)
    root.add_child(dye)

    fh = RMF.create_rmf_file(output_rmf)
    IMP.rmf.add_hierarchies(fh, [root])

    # Initial frame
    idx0 = sample_rotamer_index(lib["weight"], rng=rng)
    apply_rotamer_coords(dye, lib["coords"][idx0 + 1])
    place_dye_from_coords(dye, ca, n, c)
    IMP.rmf.save_frame(fh, f"init_rotamer_{idx0 + 1}")

    for i in range(n_samples):
        ridx = sample_rotamer_index(lib["weight"], rng=rng)
        apply_rotamer_coords(dye, lib["coords"][ridx + 1])
        place_dye_from_coords(dye, ca, n, c)
        IMP.rmf.save_frame(fh, f"rotamer_{ridx + 1}")

    click.echo(
        f"Rotamer sampling finished: library={dye_name} frames={n_samples + 1}"
    )
    click.echo(f"Wrote RMF: {output_rmf}")


if __name__ == "__main__":
    main()
