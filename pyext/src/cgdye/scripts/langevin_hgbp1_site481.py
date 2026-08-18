#!/usr/bin/env python
"""Langevin / Brownian dynamics of Alexa488 (C5-maleimide, R48 linker) on hGBP1 site 481.

The reference explicit-dye workflow (1DG3 chain A, residue 481) run through
``IMP.bff.LangevinDyeSampler``: real stochastic dynamics on the dye force
field with the protein as soft-sphere obstacles, the anchor fixed on the
residue. ``--integrator md`` (Langevin-thermostat MD, 2 fs) or ``bd``
(Brownian dynamics, 0.5 fs). The collision-gated random walk lives in
``dof_walk_hgbp1_site481.py``.
"""

from __future__ import annotations

import os

import click
import IMP
import IMP.atom
import numpy as np

from IMP.bff.label.attachment import attach_dyes
from IMP.bff.cgdye.sampling.langevin import LangevinDyeSampler
from IMP.bff.tools.paths import get_structure_dir


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--protein-pdb", default=str(get_structure_dir("1DG3.pdb")), show_default=True)
@click.option("--dye-mol2", default=str(get_structure_dir("alexa488_r48.mol2")), show_default=True)
@click.option("--chain", default="A", show_default=True)
@click.option("--residue", default=481, show_default=True, type=int)
@click.option("--integrator", type=click.Choice(["md", "bd"]), default="md", show_default=True)
@click.option("--temperature", default=300.0, show_default=True, type=float)
@click.option("--timestep-fs", default=None, type=float)
@click.option("--friction-ps", default=10.0, show_default=True, type=float)
@click.option("--n-steps", default=20000, show_default=True, type=int)
@click.option("--write-every", default=100, show_default=True, type=int)
@click.option("--minimize-steps", default=200, show_default=True, type=int)
@click.option("--interaction-sphere", default=25.0, show_default=True, type=float)
@click.option("--seed", default=481, show_default=True, type=int)
@click.option("--output-rmf", default="output/test_runs/hgbp1_1dg3_site481/langevin/trajectory.rmf3", show_default=True)
def main(protein_pdb, dye_mol2, chain, residue, integrator, temperature, timestep_fs, friction_ps,
         n_steps, write_every, minimize_steps, interaction_sphere, seed, output_rmf):
    """Langevin (md) / Brownian (bd) dynamics of the dye at hGBP1 481."""
    model = IMP.Model()
    protein = IMP.atom.read_pdb(protein_pdb, model, IMP.atom.NonWaterPDBSelector())
    dye = IMP.atom.read_mol2(dye_mol2, model)
    attach_dyes(protein, [(dye, chain, residue)], strip_site_sidechain=True)
    IMP.atom.add_bonds(protein)
    sampler = LangevinDyeSampler(protein, dye, dye_mol2, chain, residue, integrator=integrator,
                                 temperature=temperature, timestep_fs=timestep_fs, friction_ps=friction_ps,
                                 interaction_sphere=interaction_sphere, seed=seed)
    if minimize_steps > 0:
        sampler.minimize(minimize_steps)
    out_dir = os.path.dirname(os.path.abspath(output_rmf))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    traj = sampler.run(n_steps, write_every=write_every, out_rmf=output_rmf)
    line = (f"{integrator} dynamics finished: frames={traj.n_frames} steps={n_steps} dt={sampler.timestep_fs} fs "
            f"mobile={len(sampler.mobile)} obstacles={len(sampler.obstacles)} <E_pot>={np.mean(traj.potential_energy):.1f}")
    if integrator == "md":
        line += f" <T_kin>={np.nanmean([sampler.kinetic_temperature(k) for k in traj.kinetic_energy]):.0f} K"
    click.echo(line)
    click.echo(f"Wrote RMF: {output_rmf}")


if __name__ == "__main__":
    main()
