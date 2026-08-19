/**
 *  \file IMP/bff/RotamerEnergy.h
 *  \brief The interaction energy of each rotamer with its protein.
 *
 * A rotamer library is a list of candidate dye conformers; scoring turns it
 * into a weighted ensemble by asking how hard each one collides with the
 * structure it is attached to. That is an all-pairs problem — every dye atom
 * against every protein atom, for every conformer — and it is the inner loop
 * of `RotamerFRET`.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_ROTAMERENERGY_H
#define IMPBFF_ROTAMERENERGY_H

#include <IMP/bff/bff_config.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Which functional form the steric term takes.
enum RotamerPotential {
    ROTAMER_POTENTIAL_LJ = 0,     //!< 12-6 Lennard-Jones
    ROTAMER_POTENTIAL_GAUSS = 1   //!< a soft Gaussian well, for a coarse dye
};

//! Steric and electrostatic energy of every conformer against the structure.
/*!
    \param[in] rotamer_coords,n_rotamer_coords flat,
               `n_rotamers * n_dye_atoms * 3`. A raw buffer, so numpy's array
               passes straight through -- see #brownian_walk_in_volume on why.
    \param[in] protein_coords,n_protein_coords flat, `n_protein_atoms * 3`
    \param[in] rmin_ij combined \f$R_{min}\f$ per (dye atom, protein atom),
               flat `n_dye_atoms * n_protein_atoms`
    \param[in] eps_ij combined well depth, same shape
    \param[in] q_dye per dye atom, or **empty** to skip electrostatics
    \param[in] q_protein per protein atom, or empty
    \param[in] n_rotamers,n_dye_atoms,n_protein_atoms shapes
    \param[in] potential a #RotamerPotential
    \param[in] steric_cutoff pairs beyond this contribute nothing, Angstrom
    \param[in] coulomb_cutoff the same for the screened Coulomb term
    \param[in] debye_length screening length of the Debye-Huckel term, Angstrom
    \param[in] coulomb_prefactor the `q_i q_j / r` coefficient
    \return two values per rotamer: steric energy, then electrostatic energy.
            Separate rather than summed because the Boltzmann factor treats them
            differently -- the steric term is divided by RT and the screened
            Coulomb term is already in units of it.
*/
IMPBFFEXPORT std::vector<double> rotamer_interaction_energies(
        double* rotamer_coords, int n_rotamer_coords,
        double* protein_coords, int n_protein_coords,
        double* rmin_ij, int n_rmin_ij,
        double* eps_ij, int n_eps_ij,
        const std::vector<double>& q_dye,
        const std::vector<double>& q_protein,
        int n_rotamers, int n_dye_atoms, int n_protein_atoms,
        int potential,
        double steric_cutoff = 10.0,
        double coulomb_cutoff = 20.0,
        double debye_length = 10.0,
        double coulomb_prefactor = 7.0);

//! Repulsive-only LJ energy between every pair of conformers of two sets.
/*!
    The mean-field weight update needs the interaction energy of conformer *i*
    of one dye against conformer *j* of another (or against a protein, which is
    one "conformer"). Those energies **do not depend on the weights**, so they
    are computed once here and the iteration becomes a matrix-vector product.
    The Python this replaces recomputed the whole matrix inside the iteration
    loop, doing ten times the work for the same answer.

    Repulsive-only: the attractive tail (\f$r \ge R_{min}\f$) is dropped, which
    is what the mean-field treatment wants -- it is asking what is *blocked*,
    not what is bound.

    Conformer pairs whose padded axis-aligned bounding boxes do not overlap are
    skipped entirely, which is most of them once two dyes are more than a linker
    apart.

    \param[in] coords_a flat, `n_a_conf * n_a_atoms * 3`
    \param[in] coords_b flat, `n_b_conf * n_b_atoms * 3`
    \param[in] rmin,eps combined parameters, flat `n_a_atoms * n_b_atoms`
    \param[in] r_cutoff pairs beyond this contribute nothing, Angstrom
    \param[in] aabb_pad padding on each bounding box, Angstrom
    \param[in] r_floor distance clamp keeping the energy finite
    \param[out] out_view,n_out_view flat `n_a_conf * n_b_conf`. Placed before
                the defaulted parameters, not after: a parameter without a
                default may not follow one that has it, and numpy's typemap
                binds on the pair's *names*, not its position.
*/
IMPBFFEXPORT void rotamer_pair_energy_matrix(
        const std::vector<double>& coords_a,
        const std::vector<double>& coords_b,
        const std::vector<double>& rmin,
        const std::vector<double>& eps,
        int n_a_conf, int n_a_atoms, int n_b_conf, int n_b_atoms,
        double** out_view, int* n_out_view,
        double r_cutoff = 12.0,
        double aabb_pad = 3.5,
        double r_floor = 0.01);

//! Internal LJ energy of each frame, over an explicit list of atom pairs.
/*!
    The bonded exclusions have already decided which atom pairs interact, so
    this takes the list rather than rediscovering it. One energy per frame.

    The vectorised form it replaces gathered `(n_frames, n_pairs, 3)` twice --
    once for each end of every pair -- before taking a single difference. At ten
    thousand frames and twelve hundred pairs that is 288 MB per gather.

    \param[in] coords,n_coords flat, `n_frames * n_atoms * 3`, zero-copy
    \param[in] index_a,index_b atom indices of each pair
    \param[in] rmin,eps combined parameters per pair
    \param[in] n_frames,n_atoms,n_pairs shapes
    \param[in] repulsive_only drop the attractive tail beyond \f$R_{min}\f$
    \param[in] r_floor distance clamp keeping the energy finite
    \return one energy per frame
*/
IMPBFFEXPORT std::vector<double> lj_pair_energies(
        double* coords, int n_coords,
        const std::vector<int>& index_a,
        const std::vector<int>& index_b,
        const std::vector<double>& rmin,
        const std::vector<double>& eps,
        int n_frames, int n_atoms, int n_pairs,
        bool repulsive_only = true,
        double r_floor = 0.01);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_ROTAMERENERGY_H
