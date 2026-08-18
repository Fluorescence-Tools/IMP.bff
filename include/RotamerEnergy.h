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
    \param[in] rotamer_coords flat, `n_rotamers * n_dye_atoms * 3`
    \param[in] protein_coords flat, `n_protein_atoms * 3`
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
        const std::vector<double>& rotamer_coords,
        const std::vector<double>& protein_coords,
        const std::vector<double>& rmin_ij,
        const std::vector<double>& eps_ij,
        const std::vector<double>& q_dye,
        const std::vector<double>& q_protein,
        int n_rotamers, int n_dye_atoms, int n_protein_atoms,
        int potential,
        double steric_cutoff = 10.0,
        double coulomb_cutoff = 20.0,
        double debye_length = 10.0,
        double coulomb_prefactor = 7.0);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_ROTAMERENERGY_H
