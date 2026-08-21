/**
 *  \file IMP/bff/Scoring.h
 *  \brief Stage-2 scoring: the CHARMM36 table, LJ, Boltzmann, and the
 *         rotamer-protein score that turns conformers into a weighted ensemble.
 *
 * The inner kernels (all-pairs steric/electrostatic energy, the pair energy
 * matrix, the per-frame LJ over an explicit pair list) are in
 * \ref RotamerEnergy.h. What lives here is the orchestration around them: the
 * CHARMM36 parameter table and the Lorentz-Berthelot combining rules, the
 * Boltzmann weight, the axis-aligned bounding-box pre-filter, and the
 * end-to-end rotamer score that builds masks, assembles parameters, calls the
 * kernel, and returns normalised weights.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_SCORING_H
#define IMPBFF_SCORING_H

#include <IMP/bff/bff_config.h>

#include <vector>
#include <string>
#include <array>

IMPBFF_BEGIN_NAMESPACE

// --------------------------------------------------------------------------
// CHARMM36 Lennard-Jones parameters
// --------------------------------------------------------------------------

//! The CHARMM36 LJ parameters for the five elements bff dyes carry.
/*!
    \return {rmin_half, epsilon} for C, N, O, S, H. Unknown elements fall back
    to carbon, as \ref lj_params does.
*/
IMPBFFEXPORT std::array<double, 2> charmm36_lj(const std::string& element);

//! Lorentz-Berthelot cross parameters (rmin, eps) for two elements.
/*!
    \c rmin = rmin_half_i + rmin_half_j (Angstrom),
    \c eps = sqrt(eps_i * eps_j) (kcal/mol, positive well depth).
*/
IMPBFFEXPORT std::array<double, 2> lj_cross(const std::string& elem_i,
                                             const std::string& elem_j);

//! (rmin_half, epsilon) arrays for a sequence of element symbols.
/*!
    Unknown elements fall back to carbon. The two output vectors are resized
    to \c n.
*/
IMPBFFEXPORT void lj_parameter_arrays(const std::vector<std::string>& elements,
                                       std::vector<double>& rmin_half,
                                       std::vector<double>& epsilon);

//! 12-6 Lennard-Jones energy eps * ((rmin/r)^12 - 2 (rmin/r)^6).
/*!
    Vectorised over flat arrays. \c r is clamped to \c r_floor to keep the
    energy finite. \c repulsive_only zeroes the attractive tail (r >= rmin),
    \c cutoff zeroes pairs beyond that distance.
*/
IMPBFFEXPORT std::vector<double> lj_energy(const std::vector<double>& r,
                                            const std::vector<double>& rmin,
                                            const std::vector<double>& eps,
                                            bool repulsive_only = false,
                                            double cutoff = 0.0,
                                            double r_floor = 0.01);

// --------------------------------------------------------------------------
// Boltzmann weights
// --------------------------------------------------------------------------

//! Normalised Boltzmann weights from energies (log-sum-exp stabilised).
/*!
    \param energies in kcal/mol (or consistent units with kT)
    \param temperature in Kelvin
    \return weights summing to 1.0
*/
IMPBFFEXPORT std::vector<double> boltzmann_weights(
        const std::vector<double>& energies, double temperature = 298.15);

//! Aggregate frame weights into cluster weights.
/*!
    \param assignments cluster index per frame (n_frames,)
    \param frame_weights normalised weight per frame (n_frames,)
    \param n_clusters total number of clusters
    \return normalised cluster weights (n_clusters,)
*/
IMPBFFEXPORT std::vector<double> rotamer_cluster_weights(
        const std::vector<int>& assignments,
        const std::vector<double>& frame_weights,
        int n_clusters);

// --------------------------------------------------------------------------
// Axis-aligned bounding-box pre-filter
// --------------------------------------------------------------------------

//! Build padded AABBs for a batch of coordinate frames.
/*!
    \param coords flat, n_frames * n_atoms * 3
    \param n_frames, n_atoms, pad shapes and padding in Angstrom
    \return flat n_frames * 6: [xmin, ymin, zmin, xmax, ymax, zmax] per frame
*/
IMPBFFEXPORT std::vector<double> aabb_build(
        const std::vector<double>& coords,
        int n_frames, int n_atoms, double pad);

//! Build a single padded AABB for a static coordinate set.
/*!
    \param coords flat, n_atoms * 3
    \param n_atoms, pad
    \return 6 values: [xmin, ymin, zmin, xmax, ymax, zmax]
*/
IMPBFFEXPORT std::vector<double> aabb_build_single(
        const std::vector<double>& coords, int n_atoms, double pad);

//! Vectorised AABB intersection test of many frames vs one reference.
/*!
    \param rot_boxes n_frames * 6
    \param ref_box 6
    \return bool mask, n_frames: true if frame could clash
*/
IMPBFFEXPORT std::vector<int> aabb_intersects_reference(
        const std::vector<double>& rot_boxes,
        const std::vector<double>& ref_box);

// --------------------------------------------------------------------------
// Rotamer-protein score (the end-to-end orchestration)
// --------------------------------------------------------------------------

//! The result of rotamer-protein scoring.
struct RotamerScoreResult {
    std::vector<double> weights;   //!< normalised Boltzmann weights
    double partition;              //!< partition function before normalisation
    std::vector<double> energies;   //!< raw potential + electrostatic energies

    const std::vector<double>& get_weights() const { return weights; }
    const std::vector<double>& get_energies() const { return energies; }
    double get_partition() const { return partition; }
    void show(std::ostream& out) const {
        out << "RotamerScoreResult(partition=" << partition
            << ", n_weights=" << weights.size() << ")";
    }
};

//! A list of rotamer score results (the SWIG plural type).
typedef std::vector<RotamerScoreResult> RotamerScoreResults;

//! Compute Boltzmann weights for a rotamer ensemble against a protein.
/*!
    Builds the hydrogen/backbone/site masks, assembles the per-pair LJ
    parameters, calls the all-pairs steric/electrostatic kernel, and returns
    normalised Boltzmann weights. This is the C++ end of
    \c compute_rotamer_score; the mask-building helpers that operate on
    Python lists stay in the \c .i file as \c %pythoncode.

    \param rotamer_coords flat, n_rotamers * n_dye_atoms * 3
    \param protein_coords flat, n_protein_atoms * 3
    \param rmin_ij, eps_ij combined parameters, flat n_dye_atoms * n_protein_atoms
    \param q_dye per dye atom, or empty to skip electrostatics
    \param q_protein per protein atom, or empty
    \param library_weights external rotamer population weights, or empty
    \param n_rotamers, n_dye_atoms, n_protein_atoms shapes
    \param potential 0 = LJ, 1 = Gaussian
    \param temperature in K
    \return RotamerScoreResult
*/
IMPBFFEXPORT RotamerScoreResult compute_rotamer_score(
        const std::vector<double>& rotamer_coords,
        const std::vector<double>& protein_coords,
        const std::vector<double>& rmin_ij,
        const std::vector<double>& eps_ij,
        const std::vector<double>& q_dye,
        const std::vector<double>& q_protein,
        const std::vector<double>& library_weights,
        int n_rotamers, int n_dye_atoms, int n_protein_atoms,
        int potential, double temperature);

IMPBFF_END_NAMESPACE

#endif // IMPBFF_SCORING_H
