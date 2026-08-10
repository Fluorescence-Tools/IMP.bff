#!/usr/bin/env python

import os
import sys
from pathlib import Path

import click



from IMP.bff.cgdye.io.rotamer_cif import write_rotamer_library
from IMP.bff.cgdye.sampling.rotamer import (
    find_reference_rotamer_files,
    load_reference_rotamers,
)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--dye-name", required=True, help="Dye+linker name, e.g. A48_C1R"
)
@click.option(
    "--ref-lib", default="FRETpredict/FRETpredict/lib", show_default=True
)
@click.option("--cutoff", default=30, show_default=True, type=int)
@click.option("--max-frames", default=None, type=int)
@click.option(
    "--output-base", required=True, help="Output base path (without suffixes)"
)
def main(dye_name, ref_lib, cutoff, max_frames, output_base):
    pdb, dcd, weights = find_reference_rotamer_files(
        ref_lib, dye_name, cutoff=cutoff
    )
    rot = load_reference_rotamers(pdb, dcd, weights, max_frames=max_frames)

    lib = {
        "weight": [float(x) for x in rot["weights"]],
        "atom_names": list(rot["atom_names"]),
        "coords": {i + 1: rot["coords"][i] for i in range(rot["coords"].shape[0])},
    }
    out_dir = os.path.dirname(os.path.abspath(output_base))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    write_rotamer_library(output_base, lib)
    click.echo(
        f"Wrote rotamer library: base={output_base} n_rotamers={rot['coords'].shape[0]}"
    )


if __name__ == "__main__":
    main()
