/**
 * \file States.cpp
 * \brief A label's states: the kernels, the cloud, the distances between two.
 *
 * Sections in the order of IMP/bff/States.h; each is marked with the file it
 * came from.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

// -------- from AVDistance.cpp --------
/**
 * (formerly AVDistance.cpp, now a section of this file)
 * \brief Distances and reductions over accessible-volume point clouds.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/States.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/IMPCompatibility.h>

#include <cmath>
#include <limits>
#include <algorithm>
#include <random>

IMPBFF_BEGIN_NAMESPACE

namespace {
//! The one fps.json distance vocabulary, read both ways.
const std::pair<const char*, int>* distance_type_table(std::size_t& n) {
    static const std::pair<const char*, int> table[] = {
            std::make_pair("RDAMean", PROBE_PAIR_DISTANCE_MEAN),
            std::make_pair("RDAMeanE", PROBE_PAIR_DISTANCE_E),
            std::make_pair("Rmp", PROBE_PAIR_DISTANCE_MP),
            std::make_pair("Efficiency", PROBE_PAIR_EFFICIENCY),
            std::make_pair("pRDA", PROBE_PAIR_DISTANCE_DISTRIBUTION),
            std::make_pair("Rmin", PROBE_PAIR_DISTANCE_MIN)};
    n = sizeof(table) / sizeof(table[0]);
    return table;
}
}  // namespace

std::string probe_pair_distance_type_name(int distance_type) {
    std::size_t n = 0;
    const std::pair<const char*, int>* table = distance_type_table(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (table[i].second == distance_type) return table[i].first;
    }
    return std::string();
}

int probe_pair_distance_type(const std::string& name) {
    std::size_t n = 0;
    const std::pair<const char*, int>* table = distance_type_table(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (name == table[i].first) return table[i].second;
    }
    return -1;
}

std::vector<double> points_weighted_mean(const std::vector<double>& points) {
    std::vector<double> mean(3, 0.0);
    const std::size_t n = points.size() / 4;
    if (n == 0) return mean;
    double w_sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double w = points[4 * i + 3];
        w_sum += w;
        mean[0] += points[4 * i + 0] * w;
        mean[1] += points[4 * i + 1] * w;
        mean[2] += points[4 * i + 2] * w;
    }
    if (w_sum > 0.0) {
        mean[0] /= w_sum; mean[1] /= w_sum; mean[2] /= w_sum;
    }
    return mean;
}

namespace {
//! Draw `n_samples` (distance, weight-product) pairs into `out`.
/*! Shared by the public kernel, which publishes a numpy view, and by the two
    reductions below, which consume the samples and never hand them back --
    those would otherwise allocate a view only to free it. */
void sample_pairs(const std::vector<double>& p1, const std::vector<double>& p2,
                  int n_samples, int seed, double* out) {
    const std::size_t n1 = p1.size() / 4;
    const std::size_t n2 = p2.size() / 4;
    if (n_samples <= 0 || n1 == 0 || n2 == 0) return;
    std::mt19937_64 rng(static_cast<std::uint64_t>(seed));
    std::uniform_int_distribution<std::size_t> d1(0, n1 - 1);
    std::uniform_int_distribution<std::size_t> d2(0, n2 - 1);
    for (int i = 0; i < n_samples; ++i) {
        const std::size_t i1 = d1(rng);
        const std::size_t i2 = d2(rng);
        const double dx = p1[4 * i1 + 0] - p2[4 * i2 + 0];
        const double dy = p1[4 * i1 + 1] - p2[4 * i2 + 1];
        const double dz = p1[4 * i1 + 2] - p2[4 * i2 + 2];
        out[2 * i + 0] = std::sqrt(dx * dx + dy * dy + dz * dz);
        out[2 * i + 1] = p1[4 * i1 + 3] * p2[4 * i2 + 3];
    }
}
}  // namespace

void random_distances(
        double* p1, int n_p1, int n_p1c,
        double* p2, int n_p2, int n_p2c,
        int n_samples, int seed, double** output, int* n_output1,
        int* n_output2) {
    // **Not bit-comparable with the numba it replaces.** numba draws from its
    // own Mersenne stream, which no other generator reproduces, so a sampled
    // estimator can only be checked distributionally -- two estimates of the
    // same quantity, agreeing to the sampling error of ~1/sqrt(n_samples).
    const std::vector<double> v1(p1, p1 + static_cast<std::size_t>(n_p1) * n_p1c);
    const std::vector<double> v2(p2, p2 + static_cast<std::size_t>(n_p2) * n_p2c);
    const int ns = std::max(0, n_samples);
    // `new_double_view` sizes the flat buffer; the published shape is (ns, 2),
    // so dim1 is the sample count, not the element count.
    int flat = 0;
    double* out = internal::new_double_view(
            static_cast<std::size_t>(ns) * 2, output, &flat);
    if (out == nullptr) { *n_output1 = 0; *n_output2 = 2; return; }
    *n_output1 = ns;
    *n_output2 = 2;
    sample_pairs(v1, v2, ns, seed, out);
}

double average_distance(const std::vector<double>& p1,
                        const std::vector<double>& p2, int n_samples, int seed) {
    std::vector<double> d(static_cast<std::size_t>(std::max(0, n_samples)) * 2, 0.0);
    sample_pairs(p1, p2, n_samples, seed, d.data());
    double w_sum = 0.0, rda = 0.0;
    for (std::size_t i = 0; i < d.size() / 2; ++i) {
        w_sum += d[2 * i + 1];
        rda += d[2 * i + 0] * d[2 * i + 1];
    }
    return w_sum > 0.0 ? rda / w_sum : 0.0;
}

double mean_fret_distance(const std::vector<double>& p1,
                          const std::vector<double>& p2, double forster_radius,
                          int n_samples, int seed) {
    std::vector<double> d(static_cast<std::size_t>(std::max(0, n_samples)) * 2, 0.0);
    sample_pairs(p1, p2, n_samples, seed, d.data());
    double w_sum = 0.0, mean_e = 0.0;
    for (std::size_t i = 0; i < d.size() / 2; ++i) {
        const double r = d[2 * i + 0];
        const double w = d[2 * i + 1];
        // Average the *efficiency*, not the distance: 1/r^6 weights close pairs
        // far more heavily, so this is always the shorter of the two averages.
        const double e = 1.0 / (1.0 + std::pow(r / forster_radius, 6.0));
        w_sum += w;
        mean_e += e * w;
    }
    if (w_sum > 0.0) mean_e /= w_sum;
    // Both limits are degenerate, and they are not the same degeneracy:
    // E = 1 is zero separation (0), E = 0 is no transfer at all (infinity).
    if (mean_e >= 1.0) return 0.0;
    if (mean_e <= 0.0) return std::numeric_limits<double>::infinity();
    return forster_radius * std::pow(1.0 / mean_e - 1.0, 1.0 / 6.0);
}

std::vector<double> distance_sample_statistics(
        const std::vector<double>& distances, const std::vector<double>& weights,
        double forster_radius) {
    std::vector<double> out(4, 0.0);
    const std::size_t n = std::min(distances.size(), weights.size());
    double w_sum = 0.0, sum_r = 0.0, sum_r2 = 0.0, sum_e = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double r = distances[i], w = weights[i];
        w_sum += w;
        sum_r += r * w;
        sum_r2 += r * r * w;
        sum_e += w / (1.0 + std::pow(r / forster_radius, 6.0));
    }
    if (w_sum <= 0.0) return out;
    const double mean_r = sum_r / w_sum;
    const double mean_e = sum_e / w_sum;
    out[0] = mean_r;
    if (mean_e >= 1.0) out[1] = 0.0;
    else if (mean_e <= 0.0) out[1] = std::numeric_limits<double>::infinity();
    else out[1] = forster_radius * std::pow(1.0 / mean_e - 1.0, 1.0 / 6.0);
    out[2] = mean_e;
    // Clamped: the two-pass variance is exact but this one-pass form can go
    // slightly negative on a narrow distribution, and a NaN width would
    // propagate silently through every caller.
    out[3] = std::sqrt(std::max(sum_r2 / w_sum - mean_r * mean_r, 0.0));
    return out;
}

void density_to_points(
        double* density, int nx, int ny, int nz, double dg,
        const std::vector<double>& r0, double threshold,
        double** output, int* n_output1, int* n_output2) {
    const std::size_t n_total =
            static_cast<std::size_t>(nx) * ny * nz;
    const std::vector<double> dens(density, density + n_total);
    // Two passes. A view has to be sized before it is filled, and the count of
    // occupied voxels is not known in advance -- so count, allocate exactly,
    // then fill. The counting pass is a scan with no writes and is far cheaper
    // than the growth it replaces.
    std::size_t kept = 0;
    for (std::size_t k = 0; k < dens.size(); ++k) {
        if (dens[k] > threshold) ++kept;
    }
    double* points = internal::new_double_view(kept * 4, output, n_output1);
    if (output == nullptr || n_output1 == nullptr || n_output2 == nullptr) {
        // A C++ caller taking the header's defaults: the buffer is scratch,
        // publishing through NULL pointers is the crash the module build's
        // 3-argument wrapper found (new_double_view already declined to
        // publish; this function wrote through them a second time).
        std::free(points);
        return;
    }
    if (points == nullptr) { *n_output1 = 0; *n_output2 = 4; return; }
    // Shape (kept, 4): dim1 is the point count, not the element count.
    *n_output1 = static_cast<int>(kept);
    *n_output2 = 4;
    std::size_t w = 0;
    for (int ix = 0; ix < nx; ++ix) {
        const double x = dg * ix + r0[0];
        for (int iy = 0; iy < ny; ++iy) {
            const double y = dg * iy + r0[1];
            for (int iz = 0; iz < nz; ++iz) {
                const double v = dens[(static_cast<std::size_t>(ix) * ny + iy) * nz + iz];
                if (v > threshold) {
                    points[w++] = x;
                    points[w++] = y;
                    points[w++] = dg * iz + r0[2];
                    points[w++] = v;
                }
            }
        }
    }
}

void split_contact_volume(
        const std::vector<double>& density, int ng, double dg,
        const std::vector<double>& rad, const std::vector<double>& rs,
        const std::vector<double>& r0, int** output_i, int* dim1, int* dim2,
        int* dim3) {
    const std::size_t n = static_cast<std::size_t>(ng);
    int* label = internal::new_int_view(n * n * n, output_i, dim1);
    if (label == nullptr) { *dim1 = *dim2 = *dim3 = ng; return; }
    // Shape (ng, ng, ng): dim1 is the cube side, not the element count.
    *dim1 = ng;
    *dim2 = ng;
    *dim3 = ng;
    const std::size_t n_centre = rad.size();
    // Integer offset and `floor`, both deliberate: see the header. Hoisted out
    // of the voxel loop -- they depend only on the centre.
    const int half = (ng - 1) / 2;
    std::vector<int> ix0(n_centre), iy0(n_centre), iz0(n_centre), r2(n_centre);
    for (std::size_t c = 0; c < n_centre; ++c) {
        ix0[c] = static_cast<int>(std::floor((rs[3 * c + 0] - r0[0]) / dg)) + half;
        iy0[c] = static_cast<int>(std::floor((rs[3 * c + 1] - r0[1]) / dg)) + half;
        iz0[c] = static_cast<int>(std::floor((rs[3 * c + 2] - r0[2]) / dg)) + half;
        const int ri = static_cast<int>(rad[c] / dg);
        r2[c] = ri * ri;
    }
    for (int ix = 0; ix < ng; ++ix) {
        for (int iy = 0; iy < ng; ++iy) {
            for (int iz = 0; iz < ng; ++iz) {
                const std::size_t k =
                        (static_cast<std::size_t>(ix) * n + iy) * n + iz;
                if (density[k] <= 0.0) continue;
                label[k] = AV_VOXEL_FREE;
                for (std::size_t c = 0; c < n_centre; ++c) {
                    const int dx = ix - ix0[c], dy = iy - iy0[c], dz = iz - iz0[c];
                    if (dx * dx + dy * dy + dz * dz < r2[c]) {
                        label[k] = AV_VOXEL_CONTACT;
                        break;
                    }
                }
            }
        }
    }

}

void split_contact_volume_masks(
        double* density, int ng, int ng2, int ng3, double dg,
        const std::vector<double>& radius, double* rs, int n_rs, int n_rsc,
        const std::vector<double>& r0, unsigned char** contact,
        int* contact_dim1, int* contact_dim2, int* contact_dim3,
        unsigned char** free, int* free_dim1, int* free_dim2, int* free_dim3) {
    const std::size_t n_centre = static_cast<std::size_t>(n_rs) * n_rsc / 3;
    std::vector<double> rad = radius;
    if (rad.size() != n_centre) {
        // A scalar radius broadcasts over every centre when the counts disagree.
        const double r0_val = rad.empty() ? 0.0 : rad[0];
        rad.assign(n_centre, r0_val);
    }
    const std::size_t n_voxels =
            static_cast<std::size_t>(ng) * ng2 * ng3;
    const std::vector<double> centres(rs, rs + n_centre * 3);

    int* label = nullptr;
    int n = 0, d2 = 0, d3 = 0;
    split_contact_volume(std::vector<double>(density, density + n_voxels), ng,
                         dg, rad, centres, r0, &label, &n, &d2, &d3);
    unsigned char* c = static_cast<unsigned char*>(
            std::calloc(n_voxels ? n_voxels : 1, sizeof(unsigned char)));
    unsigned char* f = static_cast<unsigned char*>(
            std::calloc(n_voxels ? n_voxels : 1, sizeof(unsigned char)));
    *contact = c;
    *free = f;
    *contact_dim1 = *free_dim1 = ng;
    *contact_dim2 = *free_dim2 = ng2;
    *contact_dim3 = *free_dim3 = ng3;
    if (c != nullptr && f != nullptr && label != nullptr) {
        for (std::size_t i = 0; i < n_voxels; ++i) {
            if (label[i] == AV_VOXEL_CONTACT) c[i] = 1;
            else if (label[i] == AV_VOXEL_FREE) f[i] = 1;
        }
    }
    std::free(label);
}

double minimum_distance(const std::vector<double>& p1,
                        const std::vector<double>& p2) {
    const std::size_t n1 = p1.size() / 4, n2 = p2.size() / 4;
    double best = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < n1; ++i) {
        if (p1[4 * i + 3] <= 0.0) continue;     // a zero-weight point is not
        for (std::size_t j = 0; j < n2; ++j) {  // in the volume
            if (p2[4 * j + 3] <= 0.0) continue;
            const double dx = p1[4 * i + 0] - p2[4 * j + 0];
            const double dy = p1[4 * i + 1] - p2[4 * j + 1];
            const double dz = p1[4 * i + 2] - p2[4 * j + 2];
            const double d2 = dx * dx + dy * dy + dz * dz;
            if (d2 < best) best = d2;
        }
    }
    return std::isinf(best) ? best : std::sqrt(best);
}

double cloud_overlap(const std::vector<double>& points,
                     const std::vector<double>& reference_coords,
                     double radius) {
    if (reference_coords.size() % 3 != 0) {
        IMP_THROW("reference_coords must be three per point, got "
                          << reference_coords.size(), ValueException);
    }
    const std::size_t n = points.size() / 4, m = reference_coords.size() / 3;
    if (n == 0 || m == 0) return 0.0;
    const double r2 = radius * radius;
    double total = 0.0, touching = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double w = points[4 * i + 3];
        if (w <= 0.0) continue;
        total += w;
        for (std::size_t j = 0; j < m; ++j) {
            const double dx = points[4 * i + 0] - reference_coords[3 * j + 0];
            const double dy = points[4 * i + 1] - reference_coords[3 * j + 1];
            const double dz = points[4 * i + 2] - reference_coords[3 * j + 2];
            if (dx * dx + dy * dy + dz * dz <= r2) { touching += w; break; }
        }
    }
    return total > 0.0 ? touching / total : 0.0;
}

namespace {
//! `\\langle R_{DA}\\rangle`, `R_E` or `\\langle E\\rangle` of a fixed pair sample.
/*! The sample is a set of separation *vectors*, so the same sample answers for
    any rigid shift of the second cloud -- which is what makes the root of
    rmp_from_model_distance() smooth. */
double reduce_shifted(const std::vector<double>& dvec,
                      const std::vector<double>& weight,
                      const double* shift, int distance_type,
                      double forster_radius) {
    double w_sum = 0.0, acc = 0.0;
    const std::size_t n = weight.size();
    for (std::size_t i = 0; i < n; ++i) {
        const double dx = dvec[3 * i + 0] - shift[0];
        const double dy = dvec[3 * i + 1] - shift[1];
        const double dz = dvec[3 * i + 2] - shift[2];
        const double d2 = dx * dx + dy * dy + dz * dz;
        const double w = weight[i];
        w_sum += w;
        switch (distance_type) {
            case PROBE_PAIR_DISTANCE_MEAN:
                acc += std::sqrt(d2) * w; break;
            default: {
                const double r = std::sqrt(d2);
                acc += w / (1.0 + std::pow(r / forster_radius, 6.0));
                break;
            }
        }
    }
    if (w_sum <= 0.0) return 0.0;
    if (distance_type == PROBE_PAIR_DISTANCE_MEAN) return acc / w_sum;
    const double e = acc / w_sum;
    if (distance_type == PROBE_PAIR_EFFICIENCY) return e;
    if (e >= 1.0) return 0.0;
    if (e <= 0.0) return std::numeric_limits<double>::infinity();
    return forster_radius * std::pow(1.0 / e - 1.0, 1.0 / 6.0);
}

//! One fixed sample of separation vectors between two clouds.
void sample_separation_vectors(const std::vector<double>& p1,
                               const std::vector<double>& p2, int n_samples,
                               int seed, std::vector<double>& dvec,
                               std::vector<double>& weight) {
    const std::size_t n1 = p1.size() / 4, n2 = p2.size() / 4;
    if (n_samples <= 0 || n1 == 0 || n2 == 0) return;
    dvec.reserve(static_cast<std::size_t>(n_samples) * 3);
    weight.reserve(static_cast<std::size_t>(n_samples));
    std::mt19937_64 rng(static_cast<std::uint64_t>(seed));
    std::uniform_int_distribution<std::size_t> d1(0, n1 - 1), d2(0, n2 - 1);
    for (int i = 0; i < n_samples; ++i) {
        const std::size_t a = d1(rng), b = d2(rng);
        dvec.push_back(p1[4 * a + 0] - p2[4 * b + 0]);
        dvec.push_back(p1[4 * a + 1] - p2[4 * b + 1]);
        dvec.push_back(p1[4 * a + 2] - p2[4 * b + 2]);
        weight.push_back(p1[4 * a + 3] * p2[4 * b + 3]);
    }
}
}  // namespace

double cloud_model_distance(const std::vector<double>& p1,
                            const std::vector<double>& p2, int distance_type,
                            double forster_radius, int n_samples, int seed) {
    if (distance_type == PROBE_PAIR_DISTANCE_DISTRIBUTION) {
        // pRDA is a distribution. Returning its R_E would be a number that
        // looks like an answer and is a different quantity.
        IMP_THROW("pRDA is a distribution, not a distance: use "
                  "IMP::bff::cloud_distance_distribution() for the histogram, "
                  "or IMP::bff::random_distances() for the samples",
                  ValueException);
    }
    if (distance_type == PROBE_PAIR_DISTANCE_MIN) {
        return minimum_distance(p1, p2);
    }
    if (distance_type == PROBE_PAIR_DISTANCE_MP ||
        distance_type == PROBE_PAIR_XYZ_DISTANCE) {
        const std::vector<double> m1 = points_weighted_mean(p1);
        const std::vector<double> m2 = points_weighted_mean(p2);
        const double dx = m1[0] - m2[0], dy = m1[1] - m2[1], dz = m1[2] - m2[2];
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }
    std::vector<double> dvec, weight;
    sample_separation_vectors(p1, p2, n_samples, seed, dvec, weight);
    const double no_shift[3] = {0.0, 0.0, 0.0};
    return reduce_shifted(dvec, weight, no_shift, distance_type, forster_radius);
}

std::vector<double> cloud_distance_distribution(
        const std::vector<double>& p1, const std::vector<double>& p2,
        const std::vector<double>& axis, int n_samples, int seed,
        bool normalize) {
    if (axis.size() < 2) {
        IMP_THROW("a distribution needs at least two bin edges, got "
                          << axis.size(), ValueException);
    }
    for (std::size_t i = 1; i < axis.size(); ++i) {
        if (!(axis[i] > axis[i - 1])) {
            IMP_THROW("the bin edges must ascend; edge " << i << " is "
                                                          << axis[i] << " after "
                                                          << axis[i - 1],
                      ValueException);
        }
    }
    std::vector<double> hist(axis.size() - 1, 0.0);
    std::vector<double> d(static_cast<std::size_t>(std::max(0, n_samples)) * 2, 0.0);
    sample_pairs(p1, p2, n_samples, seed, d.data());

    double total = 0.0;
    for (std::size_t i = 0; i < d.size() / 2; ++i) {
        const double r = d[2 * i + 0], w = d[2 * i + 1];
        if (w <= 0.0 || r < axis.front() || r >= axis.back()) continue;
        // upper_bound gives the first edge above r; the bin is the one before
        const std::size_t bin = static_cast<std::size_t>(
                std::upper_bound(axis.begin(), axis.end(), r) - axis.begin() - 1);
        if (bin < hist.size()) { hist[bin] += w; total += w; }
    }
    if (normalize && total > 0.0) {
        for (double& v : hist) v /= total;
    }
    return hist;
}

double rmp_from_model_distance(const std::vector<double>& p1,
                               const std::vector<double>& p2,
                               double target_distance, int distance_type,
                               double forster_radius, double accuracy,
                               int n_samples, int seed) {
    if (p1.empty() || p2.empty()) return 0.0;
    if (distance_type == PROBE_PAIR_DISTANCE_MP ||
        distance_type == PROBE_PAIR_XYZ_DISTANCE) {
        return target_distance;   // already the quantity asked for
    }
    if (distance_type == PROBE_PAIR_DISTANCE_MIN) {
        // The root would have to be found on the *exact* minimum, which is the
        // full O(n1 n2) double sum at every bisection step. The sampled
        // minimum this function could afford is biased high, and a restraint
        // built on a biased bound is worse than no restraint.
        IMP_THROW("Rmin cannot be inverted from a sample: its minimum is a "
                  "bound, and a sampled one is biased high",
                  ValueException);
    }
    if (distance_type == PROBE_PAIR_DISTANCE_DISTRIBUTION) {
        IMP_THROW("pRDA is a distribution, not a distance", ValueException);
    }

    const std::vector<double> m1 = points_weighted_mean(p1);
    const std::vector<double> m2 = points_weighted_mean(p2);
    double ux = m2[0] - m1[0], uy = m2[1] - m1[1], uz = m2[2] - m1[2];
    const double r_mp0 = std::sqrt(ux * ux + uy * uy + uz * uz);
    if (r_mp0 > 1e-9) { ux /= r_mp0; uy /= r_mp0; uz /= r_mp0; }
    else { ux = 1.0; uy = uz = 0.0; }   // coincident means: any axis will do

    std::vector<double> dvec, weight;
    sample_separation_vectors(p1, p2, n_samples, seed, dvec, weight);
    if (weight.empty()) return 0.0;

    // Separations are p1 - p2, so moving the second cloud to r_mp along u
    // subtracts (r_mp - r_mp0) * u from every one of them.
    auto model_at = [&](double r_mp) {
        const double t = r_mp - r_mp0;
        const double shift[3] = {t * ux, t * uy, t * uz};
        return reduce_shifted(dvec, weight, shift, distance_type,
                              forster_radius);
    };

    // Efficiency falls with separation; every other convention rises with it.
    const bool rising = distance_type != PROBE_PAIR_EFFICIENCY;
    double lo = 0.0, hi = std::max(2.0 * r_mp0, 10.0);
    for (int i = 0; i < 60; ++i) {   // expand until the target is bracketed
        const double f_hi = model_at(hi);
        if (rising ? (f_hi >= target_distance) : (f_hi <= target_distance)) break;
        hi *= 1.5;
    }
    const double f_lo = model_at(lo), f_hi = model_at(hi);
    const bool bracketed = rising ? (f_lo <= target_distance && target_distance <= f_hi)
                                  : (f_hi <= target_distance && target_distance <= f_lo);
    if (!bracketed) {
        IMP_THROW("no separation reproduces "
                          << target_distance << " ("
                          << probe_pair_distance_type_name(distance_type)
                          << "): the two volumes span [" << std::min(f_lo, f_hi)
                          << ", " << std::max(f_lo, f_hi) << "]",
                  ValueException);
    }
    for (int i = 0; i < 200 && (hi - lo) > 1e-9; ++i) {
        const double mid = 0.5 * (lo + hi);
        const double f = model_at(mid);
        if (std::fabs(f - target_distance) <= accuracy) return mid;
        if ((f < target_distance) == rising) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

double chi2_score(double model_distance, double experimental_distance,
                  double error_neg, double error_pos) {
    // Residual is model minus data: a model distance that is too large is a
    // positive deviation, and a positive deviation is judged against the
    // error on that side, error_pos. Squaring hides the sign, so what this
    // convention decides is which error bar divides -- and with asymmetric
    // errors that is the whole answer.
    const double residual = model_distance - experimental_distance;
    const double err = residual < 0.0 ? error_neg : error_pos;
    if (err <= 0.0) return 0.0;
    const double z = residual / err;
    return z * z;
}

double chi2_score_capped(double model_distance, double experimental_distance,
                         double error_neg, double error_pos,
                         double max_force) {
    if (!(max_force > 0.0)) {
        return chi2_score(model_distance, experimental_distance, error_neg,
                          error_pos);
    }
    const double residual = model_distance - experimental_distance;
    const double err = residual < 0.0 ? error_neg : error_pos;
    if (err <= 0.0) return 0.0;
    // FPS's own spelling, kept so the two can be compared line by line:
    // k = 2/err^2 (SpringEngine.cs:137), drmax = MaxForce/k (:141).
    const double k = 2.0 / (err * err);
    const double drmax = max_force / k;
    const double a = std::fabs(residual);
    // Inside the knee 0.5*k*a^2 is exactly (a/err)^2, i.e. chi2_score; outside,
    // the two pieces meet at 0.5*max_force*drmax with slope max_force, so the
    // join is C^1 and the tail is straight.
    if (a <= drmax) return 0.5 * k * a * a;
    return max_force * (a - 0.5 * drmax);
}

double fret_efficiency(double distance, double forster_radius) {
    const double x = distance / forster_radius;
    const double x3 = x * x * x;
    return 1.0 / (1.0 + x3 * x3);
}

double distance_from_fret_efficiency(double efficiency, double forster_radius) {
    return forster_radius * std::pow(1.0 / efficiency - 1.0, 1.0 / 6.0);
}

IMPBFF_END_NAMESPACE

// -------- from ProbeAccessibleVolume.cpp (the States base class) --------
#include <IMP/bff/StructureIO.h>
#include <IMP/bff/internal/GridShape.h>

#include <cstring>

IMPBFF_BEGIN_NAMESPACE

// --------------------------------------------------------------------------
// States
// --------------------------------------------------------------------------

namespace {

void check_cloud(const std::vector<double>& points) {
    if (!points.empty() && points.size() % 4 != 0) {
        IMP_THROW("a point cloud is four values per point (x, y, z, weight); "
                          << points.size() << " is not a multiple of four",
                  IMP::ValueException);
    }
}

}  // namespace

States::States(const std::vector<double>& points,
               const std::vector<double>& attachment_point,
               const std::vector<double>& orientations,
               const std::string& position_name,
               const std::map<std::string, std::string>& params)
        : points_(points), attachment_point_(attachment_point),
          orientations_(orientations), position_name_(position_name),
          params_(params) {
    check_cloud(points_);
    if (!attachment_point_.empty() && attachment_point_.size() != 3) {
        IMP_THROW("an attachment point is three coordinates, not "
                          << attachment_point_.size(),
                  IMP::ValueException);
    }
    if (!orientations_.empty() &&
        orientations_.size() != static_cast<std::size_t>(get_n_points()) * 3) {
        IMP_THROW("one transition dipole per state: " << orientations_.size() / 3
                          << " against " << get_n_points() << " states",
                  IMP::ValueException);
    }
}

void States::get_points(double** out_view, int* n_out_view) const {
    internal::copy_to_view(points_, out_view, n_out_view);
}

void States::get_attachment_point(double** out_view, int* n_out_view) const {
    internal::copy_to_view(attachment_point_, out_view, n_out_view);
}

void States::get_orientations(double** out_view, int* n_out_view) const {
    internal::copy_to_view(orientations_, out_view, n_out_view);
}

void States::get_mean_position(double** out_view, int* n_out_view) const {
    // `points_weighted_mean` returns the origin for an empty or zero-weight
    // cloud. A label that has an attachment point has a better answer than the
    // origin, and every caller of this wants that one -- a buried site whose
    // volume came back empty is *at its attachment atom*, not at (0, 0, 0).
    double total_weight = 0.0;
    for (std::size_t i = 3; i < points_.size(); i += 4) total_weight += points_[i];
    if ((points_.empty() || total_weight == 0.0) &&
        attachment_point_.size() == 3) {
        internal::copy_to_view(attachment_point_, out_view, n_out_view);
        return;
    }
    internal::copy_to_view(points_weighted_mean(points_), out_view, n_out_view);
}

void States::set_points(const std::vector<double>& points) {
    check_cloud(points);
    points_ = points;
    // The dipoles were one per state; a new cloud invalidates that pairing, and
    // a stale one would silently misalign every kappa^2 read afterwards.
    if (!orientations_.empty() &&
        orientations_.size() != points_.size() / 4 * 3) {
        orientations_.clear();
    }
}

void States::set_orientations(const std::vector<double>& orientations) {
    if (!orientations.empty() &&
        orientations.size() != static_cast<std::size_t>(get_n_points()) * 3) {
        IMP_THROW("one transition dipole per state: " << orientations.size() / 3
                          << " against " << get_n_points() << " states",
                  IMP::ValueException);
    }
    orientations_ = orientations;
}

void States::set_attachment_point(const std::vector<double>& xyz) {
    if (!xyz.empty() && xyz.size() != 3) {
        IMP_THROW("an attachment point is three coordinates, not " << xyz.size(),
                  IMP::ValueException);
    }
    attachment_point_ = xyz;
}

double States::dRmp(const States& other) const {
    double *a = nullptr, *b = nullptr;
    int na = 0, nb = 0;
    get_mean_position(&a, &na);
    other.get_mean_position(&b, &nb);
    double d = 0.0;
    if (na == 3 && nb == 3) {
        const double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
        d = std::sqrt(dx * dx + dy * dy + dz * dz);
    }
    std::free(a);
    std::free(b);
    return d;
}

double States::dRDA(const States& other, int n_samples) const {
    return average_distance(points_, other.points_, n_samples, 0);
}

double States::dRDAE(const States& other, double forster_radius,
                     int n_samples) const {
    return mean_fret_distance(points_, other.points_, forster_radius,
                              n_samples, 0);
}

void States::pRDA(const States& other, const std::vector<double>& axis,
                  int n_samples, double** out_view, int* n_out_view) const {
    // `axis` is bin *edges*, so the histogram has one fewer bin than it has
    // edges, the way a histogram reads a bin-edge array.
    const std::size_t n_bins = axis.size() > 1 ? axis.size() - 1 : 0;
    double* out = internal::new_double_view(n_bins, out_view, n_out_view);
    if (out == nullptr || n_bins == 0) return;

    double* d = nullptr;
    int n = 0, nc = 0;
    random_distances(const_cast<double*>(points_.data()),
                     static_cast<int>(points_.size()) / 4, 4,
                     const_cast<double*>(other.points_.data()),
                     static_cast<int>(other.points_.size()) / 4, 4,
                     n_samples, 0, &d, &n, &nc);
    if (d == nullptr) return;

    double total = 0.0;
    for (int i = 0; i < n; ++i) {
        const double r = d[2 * i + 0], w = d[2 * i + 1];
        if (r < axis.front() || r > axis.back()) continue;
        // Upper edge closed, as numpy's histogram has it: a sample exactly at
        // the top lands in the last bin rather than nowhere.
        std::size_t b = 0;
        while (b + 1 < n_bins && r >= axis[b + 1]) ++b;
        out[b] += w;
        total += w;
    }
    std::free(d);
    if (total > 0.0) {
        for (std::size_t b = 0; b < n_bins; ++b) out[b] /= total;
    }
}

IMPBFF_END_NAMESPACE

// -------- from StatesDistance.cpp --------
/**
 * (formerly StatesDistance.cpp, now a section of this file)
 * \brief Distances between two labels, whatever represents them.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */



#include <cstdlib>

IMPBFF_BEGIN_NAMESPACE

namespace sd {

//! A cloud as a flat (x, y, z, w) buffer this side owns.
std::vector<double> cloud(const States& s) {
    double* points = NULL;
    int n = 0;
    s.get_points(&points, &n);
    std::vector<double> out(points, points + n);
    std::free(points);
    return out;
}

std::vector<double> mean_position(const States& s) {
    double* p = NULL;
    int n = 0;
    s.get_mean_position(&p, &n);
    std::vector<double> out(p, p + n);
    std::free(p);
    out.resize(3, 0.0);
    return out;
}

//! `(distance, weight)` pairs, from the one sampler.
std::vector<double> sample(const States& s1, const States& s2, int n_samples) {
    if (s1.get_n_points() == 0 || s2.get_n_points() == 0) {
        IMP_THROW("cannot sample a distance: one or both labels have no states",
                  ValueException);
    }
    const std::vector<double> c1 = cloud(s1);
    const std::vector<double> c2 = cloud(s2);
    double* out = NULL;
    int n = 0, nc = 0;
    random_distances(const_cast<double*>(c1.data()),
                     static_cast<int>(c1.size()) / 4, 4,
                     const_cast<double*>(c2.data()),
                     static_cast<int>(c2.size()) / 4, 4,
                     n_samples, 0, &out, &n, &nc);
    std::vector<double> pairs(out, out + static_cast<std::size_t>(n) * nc);
    std::free(out);
    return pairs;
}

void split(const std::vector<double>& pairs, std::vector<double>* d,
           std::vector<double>* w) {
    d->reserve(pairs.size() / 2);
    w->reserve(pairs.size() / 2);
    for (std::size_t i = 0; i + 1 < pairs.size(); i += 2) {
        d->push_back(pairs[i]);
        w->push_back(pairs[i + 1]);
    }
}

}  // namespace sd

double states_average_distance(const States& s1, const States& s2,
                               int n_samples) {
    std::vector<double> d, w;
    sd::split(sd::sample(s1, s2, n_samples), &d, &w);
    double num = 0.0, den = 0.0;
    for (std::size_t i = 0; i < d.size(); ++i) {
        num += d[i] * w[i];
        den += w[i];
    }
    return den > 0.0 ? num / den : 0.0;
}

double states_mean_fret_distance(const States& s1, const States& s2,
                                 double forster_radius, int n_samples) {
    std::vector<double> d, w;
    sd::split(sd::sample(s1, s2, n_samples), &d, &w);
    return distance_sample_statistics(d, w, forster_radius)[1];
}

double distance_between_mean_positions(const States& s1, const States& s2) {
    const std::vector<double> a = sd::mean_position(s1);
    const std::vector<double> b = sd::mean_position(s2);
    const double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double standard_deviation_of_distances(const States& s1, const States& s2,
                                       int n_samples) {
    std::vector<double> d, w;
    sd::split(sd::sample(s1, s2, n_samples), &d, &w);
    return distance_sample_statistics(d, w)[3];
}

std::vector<double> av_pair_statistics(const States& s1, const States& s2,
                                       double forster_radius, int n_samples) {
    const double rmp = distance_between_mean_positions(s1, s2);
    std::vector<double> out(4, 0.0);
    out[0] = out[1] = out[2] = rmp;
    if (s1.get_n_points() == 0 || s2.get_n_points() == 0) return out;

    std::vector<double> d, w;
    sd::split(sd::sample(s1, s2, n_samples), &d, &w);
    double total = 0.0;
    for (std::size_t i = 0; i < w.size(); ++i) total += w[i];
    // A zero total weight is a cloud that exists but carries no occupancy; the
    // mean position is still the one thing it has.
    if (total <= 0.0) return out;

    const std::vector<double> stats =
            distance_sample_statistics(d, w, forster_radius);
    out[1] = stats[0];
    out[2] = stats[1];
    out[3] = stats[3];
    return out;
}

double model_distance(const States& s1, const States& s2,
                      const std::string& distance_type, double forster_radius,
                      int n_samples) {
    if (distance_type == "Rmp") return distance_between_mean_positions(s1, s2);
    if (distance_type == "RDAMean") {
        return states_average_distance(s1, s2, n_samples);
    }
    if (distance_type == "RDAMeanE") {
        return states_mean_fret_distance(s1, s2, forster_radius, n_samples);
    }
    IMP_THROW("Unknown distance type: " << distance_type, ValueException);
}

void histogram_rda(const States& s1, const States& s2,
                   const std::vector<double>& axis, int n_samples,
                   bool normalize, double** out_view, int* n_out_view) {
    std::vector<double> d, w;
    sd::split(sd::sample(s1, s2, n_samples), &d, &w);

    std::vector<double> used_axis = axis;
    if (used_axis.size() < 2) {
        // No axis given: span the sample. One edge per sample distance, closed
        // at the upper end -- `n_samples + 1` edges, `n_samples` bins.
        const int n_b = n_samples > 0 ? n_samples : 1;
        used_axis.resize(static_cast<std::size_t>(n_b) + 1);
        for (int b = 0; b <= n_b; ++b)
            used_axis[static_cast<std::size_t>(b)] =
                    static_cast<double>(b) / n_b;
        double r_max = 0.0;
        for (std::size_t i = 0; i < d.size(); ++i) r_max = std::max(r_max, d[i]);
        if (r_max <= 0.0) r_max = 1.0;
        for (std::size_t b = 0; b < used_axis.size(); ++b)
            used_axis[b] *= r_max;
    }

    const std::size_t n_bins = used_axis.size() > 1 ? used_axis.size() - 1 : 0;
    double* out = internal::new_double_view(n_bins, out_view, n_out_view);
    if (out == NULL || n_bins == 0) {
        if (out_view == NULL) std::free(out);
        return;
    }

    double total = 0.0;
    for (std::size_t i = 0; i < d.size(); ++i) {
        const double r = d[i];
        if (r < used_axis.front() || r > used_axis.back()) continue;
        // Upper edge closed, as numpy's histogram has it.
        std::size_t b = 0;
        while (b + 1 < n_bins && r >= used_axis[b + 1]) ++b;
        out[b] += w[i];
        total += w[i];
    }
    if (normalize && total > 0.0) {
        for (std::size_t b = 0; b < n_bins; ++b) out[b] /= total;
    }
    if (out_view == NULL) {
        // No output pointers: the buffer was scratch (new_double_view did
        // not publish it), so the histogram is discarded here, not leaked.
        std::free(out);
    }
}

std::vector<double> fit_transfer_polynomial(const States& s1, const States& s2,
                                            const std::string& distance_type,
                                            double forster_radius, int degree,
                                            int n_samples) {
    // y = x: there is nothing to correct.
    std::vector<double> identity(degree + 1, 0.0);
    if (degree >= 1) identity[degree - 1] = 1.0;

    const double rmp = distance_between_mean_positions(s1, s2);
    if (s1.get_n_points() == 0 || s2.get_n_points() == 0 || rmp <= 1e-6) {
        return identity;
    }

    // Separation vectors rather than distances: the fit translates one cloud
    // along the join, and `|v + t u|` needs the vector.
    const std::vector<double> c1 = sd::cloud(s1), c2 = sd::cloud(s2);
    const std::size_t n1 = c1.size() / 4, n2 = c2.size() / 4;
    std::vector<double> v_dot_u(n_samples), v_sq(n_samples), w(n_samples);

    const std::vector<double> m1 = sd::mean_position(s1);
    const std::vector<double> m2 = sd::mean_position(s2);
    const double ux = (m2[0] - m1[0]) / rmp;
    const double uy = (m2[1] - m1[1]) / rmp;
    const double uz = (m2[2] - m1[2]) / rmp;

    // The same linear congruential draw the sampler uses, so a fit and a
    // distance drawn at the same seed see the same pairs.
    unsigned int state = 1u;
    double w_sum = 0.0;
    for (int i = 0; i < n_samples; ++i) {
        state = state * 1103515245u + 12345u;
        const std::size_t i1 = (state >> 16) % n1;
        state = state * 1103515245u + 12345u;
        const std::size_t i2 = (state >> 16) % n2;
        const double vx = c2[i2 * 4 + 0] - c1[i1 * 4 + 0];
        const double vy = c2[i2 * 4 + 1] - c1[i1 * 4 + 1];
        const double vz = c2[i2 * 4 + 2] - c1[i1 * 4 + 2];
        v_dot_u[i] = vx * ux + vy * uy + vz * uz;
        v_sq[i] = vx * vx + vy * vy + vz * vz;
        w[i] = c1[i1 * 4 + 3] * c2[i2 * 4 + 3];
        w_sum += w[i];
    }
    if (w_sum <= 0.0) return identity;

    const int n_points = 7;
    const double t_min = -std::min(rmp - 5.0, 15.0);
    const double t_max = 20.0;
    std::vector<double> xs(n_points), ys(n_points);
    for (int k = 0; k < n_points; ++k) {
        const double t = t_min + (t_max - t_min) * k / (n_points - 1);
        xs[k] = rmp + t;
        double num = 0.0, e_sum = 0.0;
        for (int i = 0; i < n_samples; ++i) {
            const double d = std::sqrt(std::max(
                    v_sq[i] + 2.0 * t * v_dot_u[i] + t * t, 1e-10));
            num += d * w[i];
            if (distance_type == "RDAMeanE") {
                const double x = d / forster_radius;
                const double x3 = x * x * x;
                e_sum += w[i] / (1.0 + x3 * x3);
            }
        }
        if (distance_type == "RDAMeanE") {
            const double mean_e = e_sum / w_sum;
            if (mean_e <= 0.0) ys[k] = num / w_sum;
            else if (mean_e >= 1.0) ys[k] = 0.0;
            else ys[k] = forster_radius * std::pow(1.0 / mean_e - 1.0, 1.0 / 6.0);
        } else {
            ys[k] = num / w_sum;
        }
    }

    // Least squares on the Vandermonde normal equations, highest power first --
    // the order `numpy.polyfit` returns and `polynomial_transfer` expects.
    const int m = degree + 1;
    std::vector<double> ata(m * m, 0.0), atb(m, 0.0);
    for (int k = 0; k < n_points; ++k) {
        std::vector<double> row(m);
        double p = 1.0;
        for (int j = m - 1; j >= 0; --j) { row[j] = p; p *= xs[k]; }
        for (int a = 0; a < m; ++a) {
            atb[a] += row[a] * ys[k];
            for (int b = 0; b < m; ++b) ata[a * m + b] += row[a] * row[b];
        }
    }
    // Gaussian elimination with partial pivoting; m is 4 in every caller.
    for (int col = 0; col < m; ++col) {
        int pivot = col;
        for (int r = col + 1; r < m; ++r) {
            if (std::fabs(ata[r * m + col]) > std::fabs(ata[pivot * m + col])) {
                pivot = r;
            }
        }
        if (std::fabs(ata[pivot * m + col]) < 1e-12) return identity;
        if (pivot != col) {
            for (int c = 0; c < m; ++c) {
                std::swap(ata[col * m + c], ata[pivot * m + c]);
            }
            std::swap(atb[col], atb[pivot]);
        }
        for (int r = col + 1; r < m; ++r) {
            const double f = ata[r * m + col] / ata[col * m + col];
            for (int c = col; c < m; ++c) ata[r * m + c] -= f * ata[col * m + c];
            atb[r] -= f * atb[col];
        }
    }
    std::vector<double> coeffs(m, 0.0);
    for (int r = m - 1; r >= 0; --r) {
        double acc = atb[r];
        for (int c = r + 1; c < m; ++c) acc -= ata[r * m + c] * coeffs[c];
        coeffs[r] = acc / ata[r * m + r];
    }
    return coeffs;
}

double gaussian_rmp_to_rda_mean(double rmp, double sigma) {
    if (!(rmp > 0.0)) return 0.0;
    return rmp + sigma * sigma / rmp;
}

double polynomial_transfer_ascending(double rmp,
                                     const std::vector<double>& coeffs) {
    // One evaluator underneath: reverse and delegate.
    std::vector<double> descending(coeffs.rbegin(), coeffs.rend());
    return polynomial_transfer(rmp, descending);
}

double mean_position_distance(const std::vector<double>& points_a,
                              const std::vector<double>& points_b,
                              const std::vector<double>& weights_a,
                              const std::vector<double>& weights_b) {
    const std::size_t na = points_a.size() / 3, nb = points_b.size() / 3;
    if (na == 0 || nb == 0) {
        IMP_THROW("both point clouds must be non-empty to define R_mp",
                  ValueException);
    }
    double ma[3] = {0, 0, 0}, mb[3] = {0, 0, 0}, wa = 0.0, wb = 0.0;
    for (std::size_t i = 0; i < na; ++i) {
        const double w = weights_a.empty() ? 1.0 : weights_a[i];
        for (int c = 0; c < 3; ++c) ma[c] += w * points_a[i * 3 + c];
        wa += w;
    }
    for (std::size_t i = 0; i < nb; ++i) {
        const double w = weights_b.empty() ? 1.0 : weights_b[i];
        for (int c = 0; c < 3; ++c) mb[c] += w * points_b[i * 3 + c];
        wb += w;
    }
    if (wa == 0.0 || wb == 0.0) {
        IMP_THROW("a cloud with zero total weight has no mean position",
                  ValueException);
    }
    double d2 = 0.0;
    for (int c = 0; c < 3; ++c) {
        const double d = ma[c] / wa - mb[c] / wb;
        d2 += d * d;
    }
    return std::sqrt(d2);
}

// --------------------------------------------------------------------------
// the Gaussian pair converter
// --------------------------------------------------------------------------

FRETDistanceConverter::FRETDistanceConverter(double forster_radius,
                                             double sigma, double distance_min,
                                             double distance_max,
                                             int n_distances)
    : forster_radius_(forster_radius), sigma_(sigma) {
    distances_.resize(n_distances);
    for (int i = 0; i < n_distances; ++i) {
        distances_[i] = distance_min + (distance_max - distance_min) * i /
                                               std::max(1, n_distances - 1);
    }
    update_efficiencies();
    update_lookup();
}

void FRETDistanceConverter::update_efficiencies() {
    efficiencies_.resize(distances_.size());
    for (std::size_t i = 0; i < distances_.size(); ++i) {
        const double x = distances_[i] / forster_radius_;
        const double x3 = x * x * x;
        efficiencies_[i] = 1.0 / (1.0 + x3 * x3);
    }
}

void FRETDistanceConverter::update_lookup() {
    const std::size_t n = distances_.size();
    d_mean_.assign(n, 0.0);
    d_mean_fret_.assign(n, 0.0);
    e_mean_.assign(n, 0.0);

    const double norm = 1.0 / (sigma_ * std::sqrt(2.0 * M_PI));
    for (std::size_t k = 0; k < n; ++k) {
        const double dcc = distances_[k];
        // Two mirrored normals: the separation of two isotropic clouds is
        // symmetric about zero, and only its magnitude is observable.
        std::vector<double> p(n);
        double total = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double a = (distances_[i] + dcc) / sigma_;
            const double b = (distances_[i] - dcc) / sigma_;
            p[i] = norm * (std::exp(-0.5 * a * a) + std::exp(-0.5 * b * b));
            total += p[i];
        }
        if (total > 0.0) {
            for (std::size_t i = 0; i < n; ++i) p[i] /= total;
        }
        double dm = 0.0, em = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            dm += p[i] * distances_[i];
            em += p[i] * efficiencies_[i];
        }
        d_mean_[k] = dm;
        e_mean_[k] = em;
        d_mean_fret_[k] =
                em > 0.0 ? forster_radius_ * std::pow(1.0 / em - 1.0, 1.0 / 6.0)
                         : dm;
    }
}

void FRETDistanceConverter::set_forster_radius(double v) {
    forster_radius_ = v;
    update_efficiencies();
    update_lookup();
}

void FRETDistanceConverter::set_sigma(double v) {
    sigma_ = v;
    update_lookup();
}

namespace sd {

//! Linear interpolation on an ascending axis, clamped at both ends.
double interp(double x, const std::vector<double>& xs,
              const std::vector<double>& ys) {
    if (xs.empty()) return 0.0;
    if (x <= xs.front()) return ys.front();
    if (x >= xs.back()) return ys.back();
    std::size_t i = 0;
    while (i + 1 < xs.size() && xs[i + 1] < x) ++i;
    const double t = (x - xs[i]) / (xs[i + 1] - xs[i]);
    return ys[i] + t * (ys[i + 1] - ys[i]);
}

}  // namespace sd

double FRETDistanceConverter::get_distance_mean(double dcc) const {
    return sd::interp(dcc, distances_, d_mean_);
}

double FRETDistanceConverter::get_distance_mean_fret(double dcc) const {
    return sd::interp(dcc, distances_, d_mean_fret_);
}

double FRETDistanceConverter::get_fret_efficiency_mean(double dcc) const {
    return sd::interp(dcc, distances_, e_mean_);
}

double FRETDistanceConverter::get_effective_distance(double value,
                                                     int distance_type) const {
    if (distance_type == PROBE_PAIR_DISTANCE_MEAN) return get_distance_mean(value);
    if (distance_type == PROBE_PAIR_DISTANCE_E) return get_distance_mean_fret(value);
    if (distance_type == PROBE_PAIR_DISTANCE_MP) return value;
    IMP_THROW("unknown probe-pair distance type " << distance_type,
              ValueException);
}

double effective_distance(double rmp,
                          const std::string& transfer_function_type,
                          double sigma_rda,
                          const std::vector<double>& coeffs) {
    if (transfer_function_type == "None") return rmp;
    if (transfer_function_type == "Gaussian") {
        return sigma_rda > 0.0 ? gaussian_rmp_to_rda_mean(rmp, sigma_rda) : rmp;
    }
    if (transfer_function_type == "Polynomial") {
        if (!coeffs.empty()) return polynomial_transfer_ascending(rmp, coeffs);
        // A calibration that names a polynomial and carries no coefficients has
        // not been fitted; sigma is the parametric stand-in.
        return sigma_rda > 0.0 ? gaussian_rmp_to_rda_mean(rmp, sigma_rda) : rmp;
    }
    return rmp;
}

// --------------------------------------------------------------------------
// Empirical corrections from a computed distance to a measured one
// --------------------------------------------------------------------------

double polynomial_transfer(double x, const std::vector<double>& coefficients) {
    double y = 0.0;
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        y = y * x + coefficients[i];
    }
    return y;
}

void polynomial_transfer_vector(
        const std::vector<double>& x, const std::vector<double>& coefficients,
        double** out_view, int* n_out_view) {
    std::vector<double> y(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        y[i] = polynomial_transfer(x[i], coefficients);
    }
    internal::copy_to_view(y, out_view, n_out_view);
}

IMPBFF_END_NAMESPACE
