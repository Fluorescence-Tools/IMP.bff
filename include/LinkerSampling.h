/**
 *  \file IMP/bff/LinkerSampling.h
 *  \brief Metropolis sampling of a dye's linker, and the library it makes.
 *
 * A dye on a protein is a chromophore on a flexible tether, and what the
 * tether can do is a handful of torsions and bond angles. Sampling those --
 * rather than every Cartesian degree of freedom -- is what makes a rotamer
 * library cheap enough to build per dye.
 *
 * \note There is **no IMP model here**. Sampling reads and writes plain
 * coordinate arrays: routing each trial configuration through
 * `IMP::core::XYZ` decorators would cost thousands of boundary crossings per
 * step for coordinates #IMP::bff::LinkerGeometry::apply has just returned.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_LINKERSAMPLING_H
#define IMPBFF_LINKERSAMPLING_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/LinkerGeometry.h>
#include <IMP/bff/Mol2IO.h>
#include <IMP/bff/RotamerLibrary.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! What a Metropolis run over a linker produced.
struct IMPBFFEXPORT LinkerSamplingResult {
    //! Flat `n_frames * n_atoms * 3`, one kept configuration per frame.
    std::vector<double> coordinates;
    //! The internal energy of each kept frame, kcal/mol.
    std::vector<double> energies;
    int n_frames, n_atoms;
    //! How many proposals were accepted, of how many made.
    int n_accepted, n_steps;

    LinkerSamplingResult()
        : n_frames(0), n_atoms(0), n_accepted(0), n_steps(0) {}

    void get_coordinates(double** out_view, int* n_out_view) const;
    void get_energies(double** out_view, int* n_out_view) const;
    //! The fraction of proposals accepted -- a run far from ~0.3 is a run
    //! whose step sizes are wrong, and it says so here rather than in a plot.
    double get_acceptance() const;

    IMP_SHOWABLE_INLINE(LinkerSamplingResult,
                        out << "LinkerSamplingResult(" << n_frames
                            << " frames, acceptance " << get_acceptance()
                            << ")");
};
IMP_VALUES(LinkerSamplingResult, LinkerSamplingResults);

//! Which torsions and bond angles of a dye can turn, and what turns with them.
/*!
    The rules, in one place:

    - a bond is rotatable unless it touches the backbone anchor (`N`, `CA`,
      `C`), lies inside a ring, or ends in a hydrogen;
    - an angle is rotatable unless its centre is in a ring or either arm is a
      hydrogen;
    - what *moves* is the side of the bond away from the anchor, which is what
      makes the anchor the fixed frame the whole library is expressed in.

    \param[in] mol2_path the dye's MOL2
    \param[in] anchor_atom the atom held fixed; the first atom when the
               structure has none by that name
    \return the geometry, ready to `apply` a configuration to
    \throw IOException when the MOL2 cannot be read
*/
IMPBFFEXPORT LinkerGeometry linker_geometry_from_mol2(
        const std::string& mol2_path, const std::string& anchor_atom = "CA");

//! Metropolis sampling of a linker's torsions and bond angles.
/*!
    Each step perturbs every degree of freedom by a Gaussian -- torsions and
    angles have their own widths, because a radian of torsion and a radian of
    bond angle are not the same size of move -- and accepts by the Metropolis
    criterion on the dye's internal energy
    (#IMP::bff::IntramolecularEnergy, which excludes the 1-2, 1-3 and
    1-4 pairs the bonded terms already hold).

    \param[in] mol2_path the dye's MOL2
    \param[in] n_steps proposals to make
    \param[in] write_every keep a frame this often
    \param[in] step_size_dih,step_size_ang the proposal widths, radians
    \param[in] temperature K
    \param[in] seed the random seed, so a run repeats
    \param[in] anchor_atom the atom held fixed
*/
IMPBFFEXPORT LinkerSamplingResult sample_linker(
        const std::string& mol2_path, int n_steps = 10000,
        int write_every = 10, double step_size_dih = 0.1,
        double step_size_ang = 0.02, double temperature = 298.15,
        int seed = 42, const std::string& anchor_atom = "CA");

//! Sample a linker, cluster the frames and weight the clusters.
/*!
    The library a dye gets when nobody has run molecular dynamics for it:
    sample the torsions, cluster the frames by RMSD without superposition
    (#IMP::bff::cluster_frames_leader), and weight each cluster by the
    Boltzmann factor of its representative's internal energy.

    \param[in] mol2_path the dye's MOL2
    \param[in] n_steps,write_every,step_size_dih,step_size_ang,temperature,seed
               as for #sample_linker
    \param[in] cluster_threshold the RMSD a cluster spans, A
    \param[in] anchor_atom the atom held fixed
    \param[in] protein_coords flat `(n, 3)` static protein atoms; non-empty
               turns on the mean-field reweighting, which lets the site the
               dye sits on decide which conformers survive
    \param[in] protein_elements one per protein atom; empty calls them carbon
    \param[in] mean_field_k,mean_field_n_iter the reweighting's strength and
               how many times it is iterated
    \return the library: one conformer per cluster, weights summing to one,
            and the jump counts between clusters along the trajectory

    A cluster's weight is the **sum of its members'** Boltzmann factors, not
    its representative's -- a broad shallow basin holds more of the ensemble
    than a narrow deep one, and only the sum says so.
*/
IMPBFFEXPORT RotamerLibrary generate_linker_rotamers(
        const std::string& mol2_path, int n_steps = 10000,
        int write_every = 10, double step_size_dih = 0.1,
        double step_size_ang = 0.02, double cluster_threshold = 0.5,
        double temperature = 298.15, int seed = 42,
        const std::string& anchor_atom = "CA",
        const std::vector<double>& protein_coords = std::vector<double>(),
        const std::vector<std::string>& protein_elements =
                std::vector<std::string>(),
        double mean_field_k = 1.0, int mean_field_n_iter = 10);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_LINKERSAMPLING_H
