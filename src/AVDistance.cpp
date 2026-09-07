/**
 * \file AVDistance.cpp
 * \brief Distances and reductions over accessible-volume point clouds.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/AVDistance.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/Base.h>

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
