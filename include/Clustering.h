/**
 * \file IMP/bff/Clustering.h
 * \brief Leader clustering of conformers by unaligned RMSD.
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_CLUSTERING_H
#define IMPBFF_CLUSTERING_H

#include <IMP/bff/bff_config.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! RMSD between two conformers of `n_atoms`, without superposition.
/** The atoms are taken in the order given -- this is the metric a rotamer
    library is clustered under, where the conformers already share a frame. */
IMPBFFEXPORT double rmsd_no_align(
        double* coords_a, int n_coords_a,
        double* coords_b, int n_coords_b);

//! RMSD between every pair of frames, as an `(n_frames, n_frames)` matrix.
/*! The second input a greedy pair selection needs
    (#IMP::bff::select_informative_pairs): how far apart the ensemble places
    its conformers, against which a candidate measurement's power to tell them
    apart is weighed. Symmetric and zero on the diagonal, so only the
    `n_frames * (n_frames - 1) / 2` distinct pairs are computed.

    \param[in] cluster_coords,n_frames,n_atoms,n_dim the ensemble,
               `(n_frames, n_atoms, 3)`
    \param[in] superpose Kabsch-align each pair before measuring. False gives
               the metric #IMP::bff::cluster_frames_leader uses, for conformers
               that already share a frame.
    \param[out] output,n_output1,n_output2 the matrix
    \throw ValueException when \p n_dim is not 3
*/
IMPBFFEXPORT void pairwise_rmsd(
        double* cluster_coords, int n_frames, int n_atoms, int n_dim,
        bool superpose,
        double** output, int* n_output1, int* n_output2);

//! Greedy leader clustering: the representative frames of `coords`.
/** `coords` is `(n_frames, n_atoms, 3)` flattened. Frame 0 is the first
    leader; each later frame joins the first leader within `threshold` and
    becomes a leader itself if there is none. The order of the sweep is part of
    the answer -- this is not a k-medoids that would find better centres.

    \note This is NOT the algorithm FRETpredict's shipped libraries were built
    with, despite what this comment used to claim. Theirs (kept verbatim under
    `junk/fretpredict_protocol/rotamer_libraries.py`) clusters on *linker
    dihedrals*: peaks per dihedral via `find_peaks`, cluster centres from the
    combinations of those peaks, two rounds of assign-and-average, then a
    representative frame per centre. What the two share is only the last step,
    the filter on minimum cluster population that `cutoff10/20/30` names.

    Measured on FRETpredict's own A48_C1R trajectory, leader clustering at a
    5 A threshold followed by a minimum-population-10 filter yields 1064
    conformers against their 711, and reproduces the r(CA->dye) distribution to
    an overlap of 0.825 -- at the 0.815 ceiling that their own trajectory
    scores against their own library. Equivalent in effect, different in
    construction; conformer counts will not match theirs and are not expected
    to. */
IMPBFFEXPORT void cluster_frames_leader(
        double* cluster_coords, int n_frames, int n_atoms, int n_dim,
        double threshold,
        int** out_view_i, int* n_out_view_i);

//! The nearest leader for every frame, by the same metric.
IMPBFFEXPORT void assign_frames_to_clusters(
        double* cluster_coords, int n_frames, int n_atoms, int n_dim,
        int* cluster_centers, int n_cluster_centers,
        int** out_view_i, int* n_out_view_i);

//! Raw transition counts to jump probabilities, rows summing to 1.
//! \note The three functions below are Markov arithmetic over a state
//!       index. They were named `rotamer_*` for their first caller; the
//!       states can be rotamers, conformers, MSM microstates or anything
//!       else counted into an `n x n` table.
/*!
    \param[in] counts flat `n * n`, row `i` the destinations from state `i`
    \param[in] n the state count
    \param[out] out_view,n_out_view flat `n * n`, row-major
    A state with no outgoing counts becomes absorbing (jumps to itself).
*/
IMPBFFEXPORT void transition_matrix_from_counts(
        const std::vector<double>& counts, int n,
        double** out_view, int* n_out_view);

//! Relaxation times of the transition matrix, `-timestep / ln(|lambda_i|)`.
/*!
    The eigenvalues of the transpose, sorted by magnitude descending; the
    leading 1.0 (the stationary state) is dropped, so \p n - 1 times come
    back in that order. A mode with |lambda| outside (1e-10, 0.99999999)
    contributes 0.0 -- infinitely slow or fully relaxed within a step.
*/
IMPBFFEXPORT void relaxation_times(
        const std::vector<double>& counts, int n, double timestep,
        double** out_view, int* n_out_view);

//! A trajectory of state indices drawn from a transition matrix.
/*!
    The Markov walk a set of jump counts implies: start somewhere, then keep
    drawing the next state from the current one's row. What it is for is
    turning a *library* -- conformers and how often they followed one another
    -- back into something that varies in time, which is what a photon
    simulation or a correlation function needs.

    \param[in] counts flat `n * n` jump counts, as
               #IMP::bff::transition_matrix_from_counts takes them
    \param[in] n the state count
    \param[in] weights one per state, for choosing where to start; empty or
               summing to zero starts at \p start_index
    \param[in] n_frames how long a trajectory to draw
    \param[in] start_index the first state; negative draws it from \p weights
    \param[in] seed the random seed, so a run repeats
    \return \p n_frames state indices
    \throw ValueException when \p counts is not `n * n`
*/
IMPBFFEXPORT std::vector<int> markov_state_trajectory(
        const std::vector<double>& counts, int n,
        const std::vector<double>& weights, int n_frames,
        int start_index = -1, int seed = 0);

//! The slowest relaxation time of the transition matrix, or 0.0.
/*! It was `rotamer_rotational_correlation_time`, which claimed more than the
    arithmetic does: this is the slowest mode of a jump process, and it is a
    *rotational* correlation time only when the states are orientations. */
IMPBFFEXPORT double slowest_relaxation_time(
        const std::vector<double>& counts, int n, double timestep);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_CLUSTERING_H
