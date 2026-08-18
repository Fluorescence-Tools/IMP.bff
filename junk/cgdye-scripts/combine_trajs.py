#!/usr/bin/env python
"""Concatenate multiple RMF3 trajectories into a single combined RMF3 file.

Usage:
    python scripts/combine_trajs.py \\
        --input-rmf run01/CX4_atto655_imp/rmfs/0.rmf3 \\
        --input-rmf run02/CX4_atto655_imp/rmfs/0.rmf3 \\
        --output-rmf combined/CX4_atto655_imp/rmfs/0.rmf3
"""
import os
import sys
from pathlib import Path

import click



import IMP
import IMP.rmf
import RMF


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--input-rmf",
    "input_rmfs",
    multiple=True,
    required=True,
    type=click.Path(exists=True),
    help="Input RMF3 file. Repeatable; frames are concatenated in order.",
)
@click.option(
    "--output-rmf",
    required=True,
    type=click.Path(),
    help="Output combined RMF3 file.",
)
def main(input_rmfs, output_rmf):
    """Merge all input RMF3 trajectories into a single output RMF3."""
    os.makedirs(os.path.dirname(os.path.abspath(output_rmf)), exist_ok=True)

    out_fh = None
    out_root = None
    total_frames = 0

    for i, rmf_path in enumerate(input_rmfs):
        print(f"  Reading {rmf_path} ...", flush=True)
        model = IMP.Model()
        IMP.set_check_level(IMP.NONE)
        in_fh = RMF.open_rmf_file_read_only(str(rmf_path))
        roots = IMP.rmf.create_hierarchies(in_fh, model)
        if not roots:
            print(f"  WARNING: no hierarchies in {rmf_path}, skipping.")
            continue

        n_frames = in_fh.get_number_of_frames()

        if out_fh is None:
            # Create the output file using the first input as the template
            out_fh = RMF.create_rmf_file(str(output_rmf))
            IMP.rmf.add_hierarchies(out_fh, roots)
            out_root = roots[0]

        for fi in range(n_frames):
            IMP.rmf.load_frame(in_fh, RMF.FrameID(fi))
            IMP.rmf.save_frame(out_fh)
            total_frames += 1

        print(f"    -> {n_frames} frames added (running total: {total_frames})", flush=True)

    if out_fh is None:
        print("ERROR: no valid input RMF files found.")
        sys.exit(1)

    print(f"\nDone. Combined trajectory written to: {output_rmf}")
    print(f"Total frames: {total_frames}")


if __name__ == "__main__":
    main()
