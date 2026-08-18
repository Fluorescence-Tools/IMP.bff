#!/usr/bin/env python
"""Calculate rotational correlation times from a kinetic rotamer library (RMF)."""

import click
import numpy as np
from IMP.bff.cgdye.io.rotamer_rmf import read_rotamer_library_rmf
from IMP.bff.cgdye.sampling.kinetic import rotamer_correlation_times, rotamer_rotational_correlation_time


@click.command()
@click.option("--input-rmf", type=click.Path(exists=True), required=True, help="Path to .rmf3 library")
@click.option("--dt", type=float, default=1.0, help="Sampling time interval (ps) used during generation")
@click.option("--show-all", is_flag=True, help="Show all relaxation times")
def main(input_rmf, dt, show_all):
    lib = read_rotamer_library_rmf(input_rmf)
    
    if "transitions" not in lib or lib["transitions"] is None:
        print(f"Error: {input_rmf} does not contain a transition matrix.")
        return

    tc = rotamer_rotational_correlation_time(lib["transitions"], dt)
    print(f"Rotational Correlation Time (slowest): {tc:.2f} ps")
    
    if show_all:
        times = rotamer_correlation_times(lib["transitions"], dt)
        print("\nRelaxation Times (ps):")
        for i, t in enumerate(times):
            if t > 0:
                print(f"  tau_{i+2}: {t:.2f}")


if __name__ == "__main__":
    main()
