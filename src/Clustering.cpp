/**
 * \file Clustering.cpp
 * \brief Leader clustering of conformers by unaligned RMSD.
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/Clustering.h>

#include <IMP/bff/StructureIO.h>

#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/IMPCompatibility.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <IMP/bff/internal/OutputView.h>

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

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

void pairwise_rmsd(double* cluster_coords, int n_frames, int n_atoms, int n_dim,
                   bool superpose,
                   double** output, int* n_output1, int* n_output2) {
    if (n_dim != 3) {
        IMP_THROW("coords must be (n_frames, n_atoms, 3)", ValueException);
    }
    // The module build's no-output wrapper reaches here with all three
    // pointers NULL (the header's defaults); publishing through them is the
    // crash histogram_rda and density_to_points already taught us about.
    if (output == nullptr || n_output1 == nullptr || n_output2 == nullptr) {
        int scratch = 0;
        double* out = internal::new_double_view(
                (std::size_t) n_frames * (std::size_t) n_frames, nullptr,
                &scratch);
        std::free(out);
        return;
    }
    int published = 0;
    double* out = internal::new_double_view(
            (std::size_t) n_frames * (std::size_t) n_frames, output, &published);
    *n_output1 = out == nullptr ? 0 : n_frames;
    *n_output2 = out == nullptr ? 0 : n_frames;
    if (out == nullptr || n_frames == 0) return;

    // One copy of each frame, not one per pair: `get_rmsd` takes vectors,
    // and at n_frames^2 pairs the slicing would cost more than the SVD.
    const std::size_t stride = (std::size_t) n_atoms * 3;
    std::vector<std::vector<double> > frames(n_frames);
    for (int i = 0; i < n_frames; ++i) {
        frames[i].assign(cluster_coords + (std::size_t) i * stride,
                         cluster_coords + (std::size_t) (i + 1) * stride);
    }
    const std::vector<int> every_atom;
    for (int i = 0; i < n_frames; ++i) {
        for (int j = i + 1; j < n_frames; ++j) {
            const double r = get_rmsd(frames[i], frames[j], every_atom, superpose);
            out[(std::size_t) i * n_frames + j] = r;
            out[(std::size_t) j * n_frames + i] = r;
        }
    }
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
    if (out_view_i == nullptr || n_out_view_i == nullptr) {
        std::free(view);  // C++ scratch: the header's NULL defaults
        return;
    }
    if (view == nullptr) return;
    for (int i = 0; i < n_frames; ++i) {
        int best = 0;
        double best_d = -1.0;
        for (int c = 0; c < n_cluster_centers; ++c) {
            const double d = msd(cluster_coords + (size_t) i * stride,
                                 cluster_coords + (size_t) cluster_centers[c] * stride, n_atoms);
            if (best_d < 0.0 || d < best_d) { best_d = d; best = c; }
        }
        // the index *into centers*
        view[i] = best;
    }
}

namespace {
//! The row-normalised matrix of raw counts, filled into `out` (`n * n`).
void transition_into(const std::vector<double>& counts, int n, double* out) {
    for (int i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < n; j++) {
            row_sum += counts[static_cast<std::size_t>(i) * n + j];
        }
        for (int j = 0; j < n; j++) {
            out[static_cast<std::size_t>(i) * n + j] =
                    row_sum == 0.0 ? (i == j ? 1.0 : 0.0)
                                   : counts[static_cast<std::size_t>(i) * n + j] / row_sum;
        }
    }
}

std::vector<double> correlation_times_impl(const std::vector<double>& counts,
                                           int n, double timestep) {
    // P, row-normalised; eigenvalues of the transpose (right-eigenvectors of
    // the row-stochastic matrix's adjoint -- what numpy eig(p.T) gives)
    Eigen::MatrixXd p(n, n);
    std::vector<double> flat(static_cast<std::size_t>(n) * n);
    transition_into(counts, n, flat.data());
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) p(i, j) = flat[j * n + i];  // transpose
    }
    const Eigen::EigenSolver<Eigen::MatrixXd> solver(p);
    std::vector<double> mags(n);
    for (int i = 0; i < n; i++) mags[i] = std::abs(solver.eigenvalues()[i]);
    std::sort(mags.begin(), mags.end(), std::greater<double>());

    std::vector<double> times;
    times.reserve(n > 0 ? static_cast<std::size_t>(n) - 1 : 0);
    for (int i = 1; i < n; i++) {
        const double v = mags[i];
        times.push_back((v > 1e-10 && v < 0.99999999)
                                ? -timestep / std::log(v) : 0.0);
    }
    return times;
}
}  // namespace

void transition_matrix_from_counts(const std::vector<double>& counts, int n,
                               double** out_view, int* n_out_view) {
    double* out = internal::new_double_view(
            static_cast<std::size_t>(n) * n, out_view, n_out_view);
    if (out == NULL) return;
    transition_into(counts, n, out);
}

void relaxation_times(const std::vector<double>& counts, int n,
                               double timestep, double** out_view,
                               int* n_out_view) {
    const std::vector<double> times = correlation_times_impl(counts, n,
                                                             timestep);
    double* out = internal::new_double_view(times.size(), out_view,
                                            n_out_view);
    if (out == NULL) return;
    for (std::size_t i = 0; i < times.size(); i++) out[i] = times[i];
}

double slowest_relaxation_time(const std::vector<double>& counts,
                                           int n, double timestep) {
    const std::vector<double> times = correlation_times_impl(counts, n,
                                                             timestep);
    double best = 0.0;
    for (double t : times) best = std::max(best, t);
    return best;
}

namespace {

//! Draw an index from `probabilities`, which need not be normalised.
int draw_index(const std::vector<double>& probabilities, double u) {
    double total = 0.0;
    for (std::size_t i = 0; i < probabilities.size(); ++i) {
        total += probabilities[i];
    }
    if (total <= 0.0) return 0;
    double cumulative = 0.0;
    for (std::size_t i = 0; i < probabilities.size(); ++i) {
        cumulative += probabilities[i] / total;
        if (u <= cumulative) return static_cast<int>(i);
    }
    return static_cast<int>(probabilities.size()) - 1;
}

}  // namespace

std::vector<int> markov_state_trajectory(const std::vector<double>& counts,
                                         int n,
                                         const std::vector<double>& weights,
                                         int n_frames, int start_index,
                                         int seed) {
    std::vector<int> out;
    if (n <= 0 || n_frames <= 0) return out;
    if (static_cast<int>(counts.size()) != n * n) {
        IMP_THROW("markov_state_trajectory: " << counts.size()
                  << " counts for " << n << " states", ValueException);
    }

    internal::OwnedView matrix;
    transition_matrix_from_counts(counts, n, &matrix.data, &matrix.size);

    boost::mt19937 rng(static_cast<boost::uint32_t>(seed));
    boost::uniform_real<double> unit(0.0, 1.0);

    int current = start_index;
    if (current < 0 || current >= n) {
        current = weights.empty() ? 0 : draw_index(weights, unit(rng));
    }
    out.push_back(current);
    for (int frame = 1; frame < n_frames; ++frame) {
        const std::vector<double> row(
                matrix.data + static_cast<std::size_t>(current) * n,
                matrix.data + static_cast<std::size_t>(current) * n + n);
        current = draw_index(row, unit(rng));
        out.push_back(current);
    }
    return out;
}

IMPBFF_END_NAMESPACE
