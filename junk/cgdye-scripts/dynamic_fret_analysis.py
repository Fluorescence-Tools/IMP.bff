#!/usr/bin/env python
"""Dynamic FRET analysis using site-specific clustering and exact Master Equation.

System: hGBP1
- Donor: Alexa488 at residue 488
- Acceptor: Alexa647 at residue 344

This script performs site-specific sampling, clusters the conformations,
and computes the PHYSICALLY EXACT FRET efficiency by solving the master equation
for conformational transitions and fluorescence decay.
"""

import os
import sys
from pathlib import Path
import numpy as np
import click

import IMP
import IMP.atom
import IMP.core
import IMP.algebra



from IMP.bff.cgdye.sampling.library_gen import LinkerSampler
from IMP.bff.cgdye.sampling.clustering import cluster_frames_leader, assign_frames_to_clusters
from IMP.bff.cgdye.sampling.boltzmann import rotamer_cluster_weights
from IMP.bff.cgdye.sampling.kinetic import rotamer_transition_matrix
from IMP.bff.cgdye.analysis.fret import fret_efficiency_exact_kinetic_pair
from IMP.bff.cgdye.utils import get_structure_dir, get_template_dir, _data_root

def _structure(name):
    """A bundled input structure, wherever IMP.bff is installed.

    These paths used to be relative to the working directory, so a script only
    ran from one place -- and stopped running at all once the data became IMP
    module data under data/cgdye.
    """
    from IMP.bff.cgdye.utils import get_structure_dir
    return str(get_structure_dir(name))



def compute_kappa2(mu_d, mu_a, r_vec):
    """Compute kappa^2 orientation factor."""
    r_norm = r_vec / np.linalg.norm(r_vec)
    cos_da = np.dot(mu_d, mu_a)
    cos_dr = np.dot(mu_d, r_norm)
    cos_ar = np.dot(mu_a, r_norm)
    return (cos_da - 3 * cos_dr * cos_ar)**2


@click.command()
@click.option("--n-steps", default=1000, show_default=True)
@click.option("--threshold", default=2.0, show_default=True, help="Clustering threshold (A)")
@click.option("--r0", default=52.0, show_default=True, help="Förster radius (A)")
@click.option("--tau0", default=4.0, show_default=True, help="Donor lifetime (ns)")
@click.option("--dt", default=0.1, show_default=True, help="Simulation step timescale (ns)")
def main(n_steps, threshold, r0, tau0, dt):
    dye_dir = get_structure_dir()
    pdb_path = _structure("1DG3.pdb")
    
    def sample_at_site(mol2_name, resnum):
        click.echo(f"Sampling site {resnum}...")
        sampler = LinkerSampler(str(dye_dir / mol2_name))
        model = sampler.model
        protein = IMP.atom.read_pdb(pdb_path, model, IMP.atom.NonWaterPDBSelector())
        
        from IMP.bff.cgdye.labeling.attachment import attach_dyes, place_dye_from_coords
        attached = attach_dyes(protein, [(sampler.hier, "A", resnum)], strip_site_sidechain=True)
        site = attached[0]["site"]
        ca, n, c = [IMP.core.XYZ(site[k]).get_coordinates() for k in ["CA", "N", "C"]]
        
        place_dye_from_coords(sampler.hier, ca, n, c)
        for i, p in sampler.idx_to_particle.items():
            pos = IMP.core.XYZ(p).get_coordinates()
            sampler.serial_to_pos0[i] = np.array([pos[0], pos[1], pos[2]])
            
        all_coords = sampler.sample(n_steps=n_steps, write_every=1) # High res for transitions
        return all_coords

    # 1. Generate trajectories
    coords_d = sample_at_site("alexa488_r48.mol2", 488)
    coords_a = sample_at_site("atto655.mol2", 344)
    
    # 2. Site-Specific Clustering
    def cluster_and_track(coords):
        centers = cluster_frames_leader(coords, threshold)
        assignments = assign_frames_to_clusters(coords, centers)
        
        # Build transition matrix
        n_c = len(centers)
        trans_counts = np.zeros((n_c, n_c), dtype=int)
        for i in range(len(assignments) - 1):
            trans_counts[assignments[i], assignments[i+1]] += 1
        
        p_matrix = rotamer_transition_matrix(trans_counts.tolist())
        
        # Simple uniform prior weight for stochastic walk segments
        weights = np.ones(len(coords)) / len(coords)
        c_weights = rotamer_cluster_weights(assignments, weights, n_c)
        
        return coords[centers], p_matrix, c_weights

    click.echo("Clustering donor ensemble...")
    clust_d, p_d, w_d = cluster_and_track(coords_d)
    click.echo("Clustering acceptor ensemble...")
    clust_a, p_a, w_a = cluster_and_track(coords_a)
    
    n_d, n_a = len(w_d), len(w_a)
    click.echo(f"States: {n_d} donor, {n_a} acceptor. Solving system ({n_d*n_a} joint states)...")

    # 3. Compute Distance and Kappa2 Matrices
    idx_r, idx_mu = 7, [0, 12]
    dist_mat = np.zeros((n_d, n_a))
    k2_mat = np.zeros((n_d, n_a))
    
    for i in range(n_d):
        rd = clust_d[i, idx_r]
        mu_d = (clust_d[i, idx_mu[1]] - clust_d[i, idx_mu[0]])
        mu_d /= np.linalg.norm(mu_d)
        for j in range(n_a):
            ra = clust_a[j, idx_r]
            mu_a = (clust_a[j, idx_mu[1]] - clust_a[j, idx_mu[0]])
            mu_a /= np.linalg.norm(mu_a)
            r_vec = ra - rd
            dist_mat[i, j] = np.linalg.norm(r_vec)
            k2_mat[i, j] = compute_kappa2(mu_d, mu_a, r_vec)

    # 4. Solve FRET using all approaches
    from IMP.bff.cgdye.analysis.fret import fret_efficiency_regimes
    
    e_exact = fret_efficiency_exact_kinetic_pair(dist_mat, k2_mat, p_d, p_a, w_d, w_a, R0=r0, tau0=tau0, dt=dt)
    regimes = fret_efficiency_regimes(dist_mat, k2_mat, w_d, w_a, R0=r0)
    
    click.echo("\n--- FRET Efficiency Comparison ---")
    click.echo(f"Physically Exact (Master Eq): {e_exact:.4f}")
    click.echo(f"Static Approximation:         {regimes['static']:.4f}")
    click.echo(f"Dynamic Approximation:        {regimes['dynamic']:.4f}")
    click.echo(f"Dynamic+ Approximation:       {regimes['dynamic_plus']:.4f}")
    
    click.echo("\n--- Ensemble Statistics ---")
    click.echo(f"Avg Distance (RDA):             {np.sum(np.outer(w_d, w_a) * dist_mat):.2f} A")
    click.echo(f"Mean Kappa^2:                   {results['kappa2_avg'] if 'results' in locals() else regimes['kappa2_avg']:.4f}")


if __name__ == "__main__":
    main()
