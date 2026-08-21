/**
 * \file QuenchingGrid.cpp
 * \brief Stamping spheres of influence onto an accessible-volume grid.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/QuenchingGrid.h>

#include <IMP/bff/AVDistance.h>
#include <IMP/bff/internal/OutputView.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>

IMPBFF_BEGIN_NAMESPACE

void center_grid_indices(
        const std::vector<double>& rs, const std::vector<double>& r0, double dg,
        int ng, const std::vector<double>& radius, std::vector<int>& ix0,
        std::vector<int>& iy0, std::vector<int>& iz0,
        std::vector<int>& radius_idx) {
    const std::size_t n = rs.size() / 3;
    ix0.assign(n, 0); iy0.assign(n, 0); iz0.assign(n, 0); radius_idx.assign(n, 0);
    const int offset = (ng - 1) / 2;
    for (std::size_t i = 0; i < n; ++i) {
        ix0[i] = static_cast<int>(std::floor((rs[3 * i + 0] - r0[0]) / dg)) + offset;
        iy0[i] = static_cast<int>(std::floor((rs[3 * i + 1] - r0[1]) / dg)) + offset;
        iz0[i] = static_cast<int>(std::floor((rs[3 * i + 2] - r0[2]) / dg)) + offset;
        radius_idx[i] = static_cast<int>(radius[i] / dg);
    }
}

namespace {
std::vector<double> stamp_spheres_impl(
        const std::vector<double>& density, int ng,
        const std::vector<double>& radius, const std::vector<double>& rs,
        const std::vector<double>& r0, double dg,
        const std::vector<double>& values, int combine) {
    const bool multiply = (combine == GRID_COMBINE_MULTIPLY);
    const double identity = multiply ? 1.0 : 0.0;
    const std::size_t n_vox = static_cast<std::size_t>(ng) * ng * ng;
    std::vector<double> factors(n_vox, identity);

    std::vector<int> ix0, iy0, iz0, r_idx;
    center_grid_indices(rs, r0, dg, ng, radius, ix0, iy0, iz0, r_idx);
    const std::size_t n_centre = ix0.size();
    if (n_centre == 0) {
        // Nothing stamped: every accessible voxel keeps the identity, and the
        // inaccessible ones do too, which is what the Python returns.
        return factors;
    }

    std::vector<double> slab(static_cast<std::size_t>(ng) * ng, identity);
    for (int ix = 0; ix < ng; ++ix) {
        std::fill(slab.begin(), slab.end(), identity);
        bool touched = false;
        for (std::size_t c = 0; c < n_centre; ++c) {
            const double value = values[c];
            if (!multiply && value == 0.0) continue;   // an added zero is a no-op
            const int ri = r_idx[c];
            const int dx = ix - ix0[c];
            const int remaining = ri * ri - dx * dx;
            if (remaining <= 0) continue;
            const int y_lo = std::max(0, iy0[c] - ri);
            const int y_hi = std::min(ng - 1, iy0[c] + ri);
            for (int iy = y_lo; iy <= y_hi; ++iy) {
                const int dy = iy - iy0[c];
                const int span2 = remaining - dy * dy;
                if (span2 <= 0) continue;
                int span = static_cast<int>(std::sqrt(static_cast<double>(span2)));
                // Membership is strict (d^2 < r^2), so drop the boundary voxel
                // when span2 is a perfect square.
                if (span * span >= span2) --span;
                if (span < 0) continue;
                const int z_lo = std::max(0, iz0[c] - span);
                const int z_hi = std::min(ng - 1, iz0[c] + span);
                for (int iz = z_lo; iz <= z_hi; ++iz) {
                    double& s = slab[static_cast<std::size_t>(iy) * ng + iz];
                    if (multiply) s *= value; else s += value;
                    touched = true;
                }
            }
        }
        if (!touched) continue;
        for (int iy = 0; iy < ng; ++iy) {
            for (int iz = 0; iz < ng; ++iz) {
                const std::size_t k =
                        (static_cast<std::size_t>(ix) * ng + iy) * ng + iz;
                if (density[k] != 0.0) {
                    factors[k] = slab[static_cast<std::size_t>(iy) * ng + iz];
                }
            }
        }
    }
    return factors;
}
}  // namespace

void stamp_spheres(const std::vector<double>& density, int ng,
        const std::vector<double>& radius, const std::vector<double>& rs,
        const std::vector<double>& r0, double dg,
        const std::vector<double>& values, int combine, double** out_view, int* n_out_view) {
    internal::copy_to_view(stamp_spheres_impl(density, ng, radius, rs, r0, dg, values, combine),
                           out_view, n_out_view);
}

void slow_factor_grid(const std::vector<double>& density, int ng, double dg,
                      const std::vector<double>& slow_radius,
                      const std::vector<double>& rs,
                      const std::vector<double>& r0,
                      const std::vector<double>& slow_fact,
                      double** out_view, int* n_out_view) {
    stamp_spheres(density, ng, slow_radius, rs, r0, dg, slow_fact,
                  GRID_COMBINE_MULTIPLY, out_view, n_out_view);
}

void quenching_rate_grid(const std::vector<double>& density, int ng, double dg,
                         const std::vector<double>& radius,
                         const std::vector<double>& rs,
                         const std::vector<double>& r0,
                         const std::vector<double>& values,
                         double** out_view, int* n_out_view) {
    stamp_spheres(density, ng, radius, rs, r0, dg, values, GRID_COMBINE_ADD,
                  out_view, n_out_view);
}

void av_contact_mask(const std::vector<double>& density, int ng, double dg,
                     const std::vector<double>& slow_radius,
                     const std::vector<double>& rs,
                     const std::vector<double>& r0,
                     int** out_view_i, int* n_out_view_i) {
    int* labels = nullptr;
    int n_labels = 0, d2 = 0, d3 = 0;
    split_contact_volume(density, ng, dg, slow_radius, rs, r0, &labels,
                         &n_labels, &d2, &d3);
    const std::size_t total = static_cast<std::size_t>(n_labels) * d2 * d3;
    int* out = internal::new_int_view(total, out_view_i, n_out_view_i);
    if (out == nullptr) { std::free(labels); return; }
    for (std::size_t i = 0; i < total; ++i) {
        out[i] = labels[i] == AV_VOXEL_CONTACT ? 1 : 0;
    }
    std::free(labels);
}

IMPBFF_END_NAMESPACE
