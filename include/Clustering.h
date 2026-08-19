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

//! Greedy leader clustering: the representative frames of `coords`.
/** `coords` is `(n_frames, n_atoms, 3)` flattened. Frame 0 is the first
    leader; each later frame joins the first leader within `threshold` and
    becomes a leader itself if there is none. The order of the sweep is part of
    the answer -- this is the algorithm FRETpredict's libraries are built with,
    not a k-medoids that would find better centres. */
IMPBFFEXPORT void cluster_frames_leader(
        double* cluster_coords, int n_frames, int n_atoms, int n_dim,
        double threshold,
        int** out_view_i, int* n_out_view_i);

//! The nearest leader for every frame, by the same metric.
IMPBFFEXPORT void assign_frames_to_clusters(
        double* cluster_coords, int n_frames, int n_atoms, int n_dim,
        int* cluster_centers, int n_cluster_centers,
        int** out_view_i, int* n_out_view_i);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_CLUSTERING_H
