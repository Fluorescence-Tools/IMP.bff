/**
 * \file Clustering.cpp
 * \brief Leader clustering of conformers by unaligned RMSD.
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/Clustering.h>
#include <IMP/bff/internal/OutputView.h>

#include <cmath>
#include <stdexcept>

IMPBFF_BEGIN_NAMESPACE

namespace {
//! Squared unaligned RMSD between two conformers laid out `(n_atoms, 3)`.
inline double msd(const double* a, const double* b, int n) {
    double acc = 0.0;
    for (int i = 0; i < n * 3; ++i) {
        const double d = a[i] - b[i];
        acc += d * d;
    }
    return acc / (double) n;
}
}

double rmsd_no_align(double* coords_a, int n_coords_a, double* coords_b, int n_coords_b) {
    if (n_coords_a != n_coords_b || n_coords_a % 3 != 0) {
        throw std::invalid_argument(
            "rmsd_no_align: both conformers must be (n_atoms, 3) and the same size");
    }
    return std::sqrt(msd(coords_a, coords_b, n_coords_a / 3));
}

void cluster_frames_leader(double* cluster_coords, int n_frames, int n_atoms, int n_dim,
                           double threshold, int** out_view_i, int* n_out_view_i) {
    if (n_dim != 3) throw std::invalid_argument("coords must be (n_frames, n_atoms, 3)");
    const int stride = n_atoms * 3;
    const double t2 = threshold * threshold;

    std::vector<int> centers;
    if (n_frames > 0) centers.push_back(0);
    for (int i = 1; i < n_frames; ++i) {
        bool is_new = true;
        for (size_t c = 0; c < centers.size(); ++c) {
            // compare squared: the threshold test is the whole inner loop, and
            // a square root per pair is the only thing in it
            if (msd(cluster_coords + (size_t) i * stride,
                    cluster_coords + (size_t) centers[c] * stride, n_atoms) < t2) {
                is_new = false;
                break;
            }
        }
        if (is_new) centers.push_back(i);
    }

    int* view = internal::new_int_view(centers.size(), out_view_i, n_out_view_i);
    if (view == nullptr) return;
    for (size_t i = 0; i < centers.size(); ++i) view[i] = centers[i];
}

void assign_frames_to_clusters(double* cluster_coords, int n_frames, int n_atoms, int n_dim,
                               int* cluster_centers, int n_cluster_centers,
                               int** out_view_i, int* n_out_view_i) {
    if (n_dim != 3) throw std::invalid_argument("coords must be (n_frames, n_atoms, 3)");
    const int stride = n_atoms * 3;

    int* view = internal::new_int_view(n_frames < 0 ? 0 : n_frames, out_view_i, n_out_view_i);
    if (view == nullptr) return;
    for (int i = 0; i < n_frames; ++i) {
        int best = 0;
        double best_d = -1.0;
        for (int c = 0; c < n_cluster_centers; ++c) {
            const double d = msd(cluster_coords + (size_t) i * stride,
                                 cluster_coords + (size_t) cluster_centers[c] * stride, n_atoms);
            if (best_d < 0.0 || d < best_d) { best_d = d; best = c; }
        }
        // the index *into centers*, which is what the Python returned
        view[i] = best;
    }
}

IMPBFF_END_NAMESPACE
