#!/usr/bin/env python
"""Generate rotamer libraries for all available dyes."""

import os
import sys
from pathlib import Path

import click



from IMP.bff.cgdye.sampling.library_gen import generate_rotamers
from IMP.bff.cgdye.io.rotamer_rmf import write_rotamer_library_rmf


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--n-steps", default=1000, show_default=True, type=int)
@click.option("--cluster-threshold", default=0.5, show_default=True, type=float)
@click.option("--output-dir", default="cgdye/templates/rotamer", show_default=True)
def main(n_steps, cluster_threshold, output_dir):
    dye_dir = Path("inputs/structures")
    # Find all MOL2 files that look like dyes (not the protein 1DG3)
    mol2_files = [f for f in dye_dir.glob("*.mol2") if "1DG3" not in f.name]
    
    if not mol2_files:
        click.echo("No dye MOL2 files found in inputs/structures/")
        return

    os.makedirs(output_dir, exist_ok=True)

    for mol2 in mol2_files:
        dye_name = mol2.stem.replace("_r48", "")
        click.echo(f"Generating library for {dye_name}...")
        
        try:
            library = generate_rotamers(
                str(mol2),
                n_steps=n_steps,
                cluster_threshold=cluster_threshold
            )
            
            output_path = os.path.join(output_dir, f"{dye_name}.rmf3")
            write_rotamer_library_rmf(output_path, library)
            
            n_rot = len(library["weight"])
            click.echo(f"  Done: {n_rot} rotamers saved to {output_path}")
        except Exception as e:
            click.echo(f"  Error generating {dye_name}: {e}")


if __name__ == "__main__":
    main()
