#!/usr/bin/env python
"""Example of FRET efficiency calculation regimes (Static, Dynamic, Dynamic+).

System: hGBP1
- Donor: Alexa488 at residue 488
- Acceptor: Alexa647 at residue 344
"""

import os
import sys
from pathlib import Path
import numpy as np

import IMP
import IMP.atom
import IMP.core
import IMP.algebra



from IMP.bff.cgdye.io.rotamer_rmf import read_rotamer_library_rmf
from IMP.bff.cgdye.labeling.attachment import backbone_frame, place_dye_from_coords


def compute_kappa2(mu_d, mu_a, r_vec):
    """Compute kappa^2 orientation factor."""
    r_norm = r_vec / np.linalg.norm(r_vec)
    cos_da = np.dot(mu_d, mu_a)
    cos_dr = np.dot(mu_d, r_norm)
    cos_ar = np.dot(mu_a, r_norm)
    return (cos_da - 3 * cos_dr * cos_ar)**2


def main():
    # 1. Load system
    model = IMP.Model()
    protein = IMP.atom.read_pdb("inputs/structures/1DG3.pdb", model, IMP.atom.NonWaterPDBSelector())
    
    # 2. Setup backbone frames for labeling sites
    # NOTE: user requested 488 and 344
    site_d = backbone_frame(protein, "A", 488)
    site_a = backbone_frame(protein, "A", 344)
    
    # 3. Load rotamer libraries
    lib_d = read_rotamer_library_rmf("cgdye/templates/rotamer/A48_C1R.rmf3")
    lib_a = read_rotamer_library_rmf("cgdye/templates/rotamer/A64_C2R.rmf3")
    
    # Dipole atoms (from libraries.yml / template.cif)
    # A48: C2, C13. center: C7
    # A64: C4, C15. center: C29
    def get_indices(atom_names, target_names):
        return [atom_names.index(n) for n in target_names]
    
    idx_d_mu = get_indices(lib_d["atom_names"], ["C2", "C13"])
    idx_d_r = lib_d["atom_names"].index("C7")
    
    idx_a_mu = get_indices(lib_a["atom_names"], ["C4", "C15"])
    idx_a_r = lib_a["atom_names"].index("C29")
    
    # 4. Transform all rotamers to protein frame
    def transform_coords(lib, frame):
        trans = frame.get_transformation_to()
        out = []
        for rid in sorted(lib["coords"].keys()):
            coords = lib["coords"][rid]
            # Convert to IMP transformation for consistency or use numpy
            # We'll use simple numpy math for speed here
            # R = trans.get_rotation(), T = trans.get_translation()
            rot = trans.get_rotation()
            # IMP rotation to 3x3 matrix
            R = np.array([[rot.get_rotated(IMP.algebra.Vector3D(1,0,0))[j] for j in range(3)],
                          [rot.get_rotated(IMP.algebra.Vector3D(0,1,0))[j] for j in range(3)],
                          [rot.get_rotated(IMP.algebra.Vector3D(0,0,1))[j] for j in range(3)]]).T
            T = np.array([trans.get_translation()[j] for j in range(3)])
            
            # transform all atoms
            t_coords = np.dot(coords, R.T) + T
            out.append(t_coords)
        return np.array(out)

    coords_d = transform_coords(lib_d, site_d)
    coords_a = transform_coords(lib_a, site_a)
    
    weights_d = np.array(lib_d["weight"])
    weights_a = np.array(lib_a["weight"])
    
    # 5. Compute pairwise distances and kappa2
    n_d = len(weights_d)
    n_a = len(weights_a)
    
    # assumed R0 for A48-A64 (kappa2=2/3)
    R0 = 52.0 
    
    eff_matrix = np.zeros((n_d, n_a))
    kappa2_matrix = np.zeros((n_d, n_a))
    dist_matrix = np.zeros((n_d, n_a))
    
    print(f"Computing FRET for {n_d} donor x {n_a} acceptor rotamers...")
    
    for i in range(n_d):
        rd = coords_d[i, idx_d_r]
        mu_d = coords_d[i, idx_d_mu[1]] - coords_d[i, idx_d_mu[0]]
        mu_d /= np.linalg.norm(mu_d)
        
        for j in range(n_a):
            ra = coords_a[j, idx_a_r]
            mu_a = coords_a[j, idx_a_mu[1]] - coords_a[j, idx_a_mu[0]]
            mu_a /= np.linalg.norm(mu_a)
            
            r_vec = ra - rd
            r = np.linalg.norm(r_vec)
            k2 = compute_kappa2(mu_d, mu_a, r_vec)
            
            dist_matrix[i, j] = r
            kappa2_matrix[i, j] = k2
            
            # Rate constant k_fret / k_rad = (R0/r)^6 * (k2 / (2/3))
            rate_ratio = (R0 / r)**6 * (k2 / (2/3.0))
            eff_matrix[i, j] = rate_ratio / (1 + rate_ratio)

    # Combined weights
    w_ij = np.outer(weights_d, weights_a)
    
    # 6. Calculate Regimes
    
    # Static: Average of efficiencies
    e_static = np.sum(w_ij * eff_matrix)
    
    # Dynamic (Dynamic1 in FRETpredict): 
    # Efficiency computed with <k2> for each rotamer pair
    k2_avg = np.sum(w_ij * kappa2_matrix)
    rate_ratio_dyn = (R0 / dist_matrix)**6 * (k2_avg / (2/3.0))
    eff_dyn = rate_ratio_dyn / (1 + rate_ratio_dyn)
    e_dynamic = np.sum(w_ij * eff_dyn)
    
    # Dynamic+ (Dynamic2 in FRETpredict):
    # Efficiency computed from average rate
    rate_avg = np.sum(w_ij * (R0 / dist_matrix)**6 * (kappa2_matrix / (2/3.0)))
    e_dynamic_plus = rate_avg / (1 + rate_avg)
    
    print("\n--- FRET Efficiency Regimes ---")
    print(f"Static:       {e_static:.4f}")
    print(f"Dynamic:      {e_dynamic:.4f}")
    print(f"Dynamic+:     {e_dynamic_plus:.4f}")
    print(f"<kappa^2>:    {k2_avg:.4f}")
    print(f"<distance>:   {np.sum(w_ij * dist_matrix):.2f} A")


if __name__ == "__main__":
    main()
