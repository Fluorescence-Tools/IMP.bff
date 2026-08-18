/**
 * \file QuenchingMap.cpp
 * \brief Mobility, quenching-rate and FRET-rate fields on an AV grid.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/QuenchingMap.h>

#include <cmath>

IMPBFF_BEGIN_NAMESPACE

std::vector<double> slow_near_atoms(
        const std::vector<double>& d_map, const std::vector<double>& density,
        const std::vector<double>& axis, const std::vector<double>& r0,
        const std::vector<double>& atoms_xyz, double min_distance_sq,
        double factor) {
    const std::size_t ng = axis.size();
    const std::size_t n_atoms = atoms_xyz.size() / 3;
    std::vector<double> out(density.size(), 0.0);
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
    return out;
}

std::vector<double> quenching_map(
        const std::vector<double>& density, const std::vector<double>& axis,
        const std::vector<double>& r0, const std::vector<double>& atoms_xyz,
        const std::vector<double>& kQ, const std::vector<double>& rC,
        double dye_radius, double inv_tau0) {
    const std::size_t ng = axis.size();
    const std::size_t n_atoms = atoms_xyz.size() / 3;
    std::vector<double> out(density.size(), 0.0);
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
    return out;
}

std::vector<double> fret_map(
        const std::vector<double>& density_d, const std::vector<double>& density_a,
        const std::vector<double>& axis_d, const std::vector<double>& axis_a,
        const std::vector<double>& r0_d, const std::vector<double>& r0_a,
        double r0_6, double kf, int step) {
    const std::size_t ng_d = axis_d.size();
    const std::size_t ng_a = axis_a.size();
    std::vector<double> out(density_d.size(), 0.0);
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
    return out;
}

IMPBFF_END_NAMESPACE
