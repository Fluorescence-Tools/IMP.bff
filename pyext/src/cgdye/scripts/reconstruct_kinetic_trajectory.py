#!/usr/bin/env python
"""Reconstruct a kinetic trajectory from a rotamer library."""

import os
import sys
from pathlib import Path

import click
import IMP
import IMP.atom
import IMP.rmf
import RMF



from IMP.bff.cgdye.io.rotamer_rmf import read_rotamer_library_rmf
from IMP.bff.cgdye.sampling.rotamer import apply_rotamer_coordinates
from IMP.bff.cgdye.sampling.kinetic import reconstruct_rotamer_trajectory


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--lib-rmf", required=True, help="Path to kinetic rotamer library (.rmf3)")
@click.option("--n-frames", default=100, show_default=True, type=int)
@click.option("--output-rmf", required=True, help="Path to output reconstructed trajectory")
@click.option("--seed", default=42, show_default=True, type=int)
def main(lib_rmf, n_frames, output_rmf, seed):
    click.echo(f"Loading library: {lib_rmf}")
    lib = read_rotamer_library_rmf(lib_rmf)
    
    if lib["transitions"] is None:
        click.echo("Error: Library is not kinetic (no transition matrix found)")
        return

    click.echo(f"Reconstructing {n_frames} frame trajectory...")
    indices = reconstruct_rotamer_trajectory(lib, n_frames, seed=seed)
    
    # Setup IMP model for output
    model = IMP.Model()
    root = IMP.atom.Hierarchy.setup_particle(IMP.Particle(model, "reconstructed"))
    
    # Create a dummy chain and residue to hold the atoms
    chain = IMP.atom.Chain.setup_particle(IMP.Particle(model, "A"), "A")
    root.add_child(chain)
    res = IMP.atom.Residue.setup_particle(IMP.Particle(model, "DYE"), IMP.atom.ResidueType("DYE"), 1)
    chain.add_child(res)

    # Create dye atoms based on library names
    particles = []
    for name in lib["atom_names"]:
        p = IMP.Particle(model, name)
        IMP.atom.Atom.setup_particle(p, IMP.atom.AtomType("C"))
        IMP.core.XYZR.setup_particle(p)
        IMP.core.XYZR(p).set_radius(1.0)
        res.add_child(IMP.atom.Hierarchy.setup_particle(p))
        particles.append(p)
        
    IMP.atom.add_bonds(root)
    
    fh = RMF.create_rmf_file(output_rmf)
    IMP.rmf.add_hierarchies(fh, [root])
    
    for i, ridx in enumerate(indices):
        # apply_rotamer_coordinates expect 0-based list/array for coords
        # but lib["coords"] is 1-based dict in my reader
        coords = lib["coords"][ridx + 1]
        for p, c in zip(particles, coords):
            IMP.core.XYZ(p).set_coordinates(IMP.algebra.Vector3D(c[0], c[1], c[2]))
            
        IMP.rmf.save_frame(fh, f"step_{i}_rotamer_{ridx + 1}")
        
    click.echo(f"Wrote reconstructed trajectory: {output_rmf}")


if __name__ == "__main__":
    main()
