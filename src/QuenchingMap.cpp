/**
 * \file QuenchingMap.cpp
 * \brief Mobility, quenching-rate and FRET-rate fields on an AV grid.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/QuenchingMap.h>
#include <IMP/bff/internal/OutputView.h>

#include <IMP/exception.h>

#include <cmath>

IMPBFF_BEGIN_NAMESPACE

// Named, not anonymous: IMP compiles this module as one translation unit.
namespace qmap {

//! The side of a cube of `n` voxels, rounded to the nearest integer.
/*! `cbrt` of a perfect cube can land a hair below it in double arithmetic, and
    truncating then gives `ng - 1` -- one axis short of the grid. */
int grid_side(std::size_t n) {
    return static_cast<int>(std::floor(std::cbrt(static_cast<double>(n)) + 0.5));
}

}  // namespace qmap

void slow_near_atoms(
        const std::vector<double>& d_map, const std::vector<double>& density,
        const std::vector<double>& axis, const std::vector<double>& r0,
        const std::vector<double>& atoms_xyz, double min_distance_sq,
        double factor, double** out_view, int* n_out_view) {
    const std::size_t ng = axis.size();
    const std::size_t n_atoms = atoms_xyz.size() / 3;
    double* out = internal::new_double_view(density.size(), out_view, n_out_view);
    if (out == nullptr) return;
    for (std::size_t ix = 0; ix < ng; ++ix) {
        const double x = axis[ix] + r0[0];
        for (std::size_t iy = 0; iy < ng; ++iy) {
            const double y = axis[iy] + r0[1];
            for (std::size_t iz = 0; iz < ng; ++iz) {
                const std::size_t k = (ix * ng + iy) * ng + iz;
                if (density[k] <= 0.0) continue;
                const double z = axis[iz] + r0[2];
                double slow = 1.0;
                for (std::size_t a = 0; a < n_atoms; ++a) {
                    const double dx = atoms_xyz[3 * a + 0] - x;
                    const double dy = atoms_xyz[3 * a + 1] - y;
                    const double dz = atoms_xyz[3 * a + 2] - z;
                    // Applied once per contacting atom, so it compounds.
                    if (dx * dx + dy * dy + dz * dz < min_distance_sq) slow *= factor;
                }
                out[k] = d_map[k] * slow;
            }
        }
    }
}

void quenching_map(
        const std::vector<double>& density, const std::vector<double>& axis,
        const std::vector<double>& r0, const std::vector<double>& atoms_xyz,
        const std::vector<double>& kQ, const std::vector<double>& rC,
        double dye_radius, double inv_tau0,
        double** out_view, int* n_out_view) {
    const std::size_t ng = axis.size();
    const std::size_t n_atoms = atoms_xyz.size() / 3;
    double* out = internal::new_double_view(density.size(), out_view, n_out_view);
    if (out == nullptr) return;
    for (std::size_t ix = 0; ix < ng; ++ix) {
        const double x = axis[ix] + r0[0];
        for (std::size_t iy = 0; iy < ng; ++iy) {
            const double y = axis[iy] + r0[1];
            for (std::size_t iz = 0; iz < ng; ++iz) {
                const std::size_t k = (ix * ng + iy) * ng + iz;
                // Outside the volume the rate stays zero, not 1/tau0.
                if (density[k] <= 0.0) continue;
                const double z = axis[iz] + r0[2];
                double v = inv_tau0;
                for (std::size_t a = 0; a < n_atoms; ++a) {
                    if (kQ[a] == 0.0 || rC[a] == 0.0) continue;
                    const double dx = atoms_xyz[3 * a + 0] - x;
                    const double dy = atoms_xyz[3 * a + 1] - y;
                    const double dz = atoms_xyz[3 * a + 2] - z;
                    const double d = std::sqrt(dx * dx + dy * dy + dz * dz) - dye_radius;
                    v += kQ[a] * std::exp(-d / rC[a]);
                }
                out[k] = v;
            }
        }
    }
}

void fret_map(
        const std::vector<double>& density_d, const std::vector<double>& density_a,
        const std::vector<double>& axis_d, const std::vector<double>& axis_a,
        const std::vector<double>& r0_d, const std::vector<double>& r0_a,
        double r0_6, double kf, int step,
        double** out_view, int* n_out_view) {
    const std::size_t ng_d = axis_d.size();
    const std::size_t ng_a = axis_a.size();
    double* out = internal::new_double_view(density_d.size(), out_view, n_out_view);
    if (out == nullptr) return;
    if (step < 1) step = 1;
    for (std::size_t ixd = 0; ixd < ng_d; ++ixd) {
        const double x = axis_d[ixd] + r0_d[0];
        for (std::size_t iyd = 0; iyd < ng_d; ++iyd) {
            const double y = axis_d[iyd] + r0_d[1];
            for (std::size_t izd = 0; izd < ng_d; ++izd) {
                const std::size_t kd = (ixd * ng_d + iyd) * ng_d + izd;
                if (density_d[kd] <= 0.0) continue;
                const double z = axis_d[izd] + r0_d[2];
                double t_ret = 0.0;   // accumulate transfer *time*, not rate
                double weight = 0.0;
                for (std::size_t ixa = 0; ixa < ng_a; ixa += step) {
                    const double ddx = axis_a[ixa] + r0_a[0] - x;
                    const double sx = ddx * ddx;
                    for (std::size_t iya = 0; iya < ng_a; iya += step) {
                        const double ddy = axis_a[iya] + r0_a[1] - y;
                        const double sy = ddy * ddy;
                        for (std::size_t iza = 0; iza < ng_a; iza += step) {
                            const std::size_t ka = (ixa * ng_a + iya) * ng_a + iza;
                            const double da = density_a[ka];
                            if (da <= 0.0) continue;
                            const double ddz = axis_a[iza] + r0_a[2] - z;
                            const double rda2 = sx + sy + ddz * ddz;
                            const double r2 = r0_6 / (rda2 * rda2 * rda2);
                            t_ret += da / (r2 * kf);
                            weight += da;
                        }
                    }
                }
                // Harmonic mean: the weighted mean transfer time, inverted.
                if (weight > 0.0 && t_ret > 0.0) out[kd] = weight / t_ret;
            }
        }
    }
}

std::vector<double> grid_axis(int ng, double dg) {
    std::vector<double> axis(ng > 0 ? ng : 0);
    const int centre = (ng - 1) / 2;
    for (int i = 0; i < ng; ++i) axis[i] = (i - centre) * dg;
    return axis;
}

void diffusion_coefficient_map(const std::vector<double>& density,
                               const std::vector<double>& r0, double dg,
                               const std::vector<double>& atoms_xyz,
                               double free_diffusion, double min_distance,
                               double slow_factor,
                               const std::vector<double>& base,
                               double** out_view, int* n_out_view) {
    const int ng = qmap::grid_side(density.size());
    std::vector<double> d_map;
    if (base.empty()) {
        d_map.assign(density.size(), free_diffusion);
    } else {
        if (base.size() != density.size()) {
            IMP_THROW("the base coefficient map must have as many voxels as the "
                      "density: " << base.size() << " against " << density.size(),
                      ValueException);
        }
        d_map = base;
    }
    slow_near_atoms(d_map, density, grid_axis(ng, dg), r0, atoms_xyz,
                    min_distance * min_distance, slow_factor, out_view,
                    n_out_view);
}

void quenching_rate_map(const std::vector<double>& density,
                        const std::vector<double>& r0, double dg,
                        const std::vector<double>& atoms_xyz,
                        const std::vector<double>& kQ,
                        const std::vector<double>& rC, double tau0,
                        double dye_radius, double** out_view,
                        int* n_out_view) {
    const int ng = qmap::grid_side(density.size());
    quenching_map(density, grid_axis(ng, dg), r0, atoms_xyz, kQ, rC, dye_radius,
                  tau0 > 0.0 ? 1.0 / tau0 : 0.0, out_view, n_out_view);
}

void fret_rate_map(const std::vector<double>& density_donor,
                   const std::vector<double>& density_acceptor,
                   const std::vector<double>& r0_donor,
                   const std::vector<double>& r0_acceptor, double dg_donor,
                   double dg_acceptor, double forster_radius, double kf,
                   int acceptor_step, double** out_view, int* n_out_view) {
    if (!(kf > 0.0)) {
        IMP_THROW("kf (the donor's radiative rate) must be positive, not " << kf,
                  ValueException);
    }
    const double r0_6 = std::pow(forster_radius, 6);
    fret_map(density_donor, density_acceptor,
             grid_axis(qmap::grid_side(density_donor.size()), dg_donor),
             grid_axis(qmap::grid_side(density_acceptor.size()), dg_acceptor),
             r0_donor, r0_acceptor, r0_6, kf,
             acceptor_step > 1 ? acceptor_step : 1, out_view, n_out_view);
}

IMPBFF_END_NAMESPACE
