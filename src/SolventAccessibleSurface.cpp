/**
 * \file SolventAccessibleSurface.cpp
 * \brief Shrake-Rupley solvent-accessible surface, and per-frame quenching rates.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/SolventAccessibleSurface.h>

#include <cmath>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

std::vector<double> sphere_points(int n) {
    std::vector<double> out;
    if (n <= 0) return out;
    out.resize(static_cast<std::size_t>(n) * 3);
    const double inc = M_PI * (3.0 - std::sqrt(5.0));
    const double offset = 2.0 / static_cast<double>(n);
    for (int k = 0; k < n; ++k) {
        const double y = k * offset - 1.0 + offset / 2.0;
        const double rd = std::sqrt(std::max(0.0, 1.0 - y * y));
        const double phi = k * inc;
        out[3 * k + 0] = std::cos(phi) * rd;
        out[3 * k + 1] = y;
        out[3 * k + 2] = std::sin(phi) * rd;
    }
    return out;
}

std::vector<double> solvent_accessible_surface_area(
        const std::vector<double>& xyz, const std::vector<double>& vdw,
        const std::vector<int>& probe_atom_indices,
        const std::vector<double>& points, double probe, double radius) {
    const std::size_t n_probe = probe_atom_indices.size();
    const std::size_t n_sphere = points.size() / 3;
    const std::size_t n_atoms = xyz.size() / 3;
    std::vector<double> asa(n_probe, 0.0);
    if (n_probe == 0 || n_sphere == 0 || n_atoms == 0) return asa;

    const double c = 4.0 * M_PI / static_cast<double>(n_sphere);
    const double neighbour_cutoff2 = (2.0 * radius + probe) * (2.0 * radius + probe);
    const double occlusion2 = (radius + probe) * (radius + probe);
    std::vector<std::size_t> neighbours(n_atoms);

    for (std::size_t i = 0; i < n_probe; ++i) {
        const std::size_t p = static_cast<std::size_t>(probe_atom_indices[i]);
        const double ax = xyz[3 * p + 0], ay = xyz[3 * p + 1], az = xyz[3 * p + 2];
        std::size_t n_neighbour = 0;
        for (std::size_t a = 0; a < n_atoms; ++a) {
            if (a == p) continue;
            const double dx = ax - xyz[3 * a + 0];
            const double dy = ay - xyz[3 * a + 1];
            const double dz = az - xyz[3 * a + 2];
            if (dx * dx + dy * dy + dz * dz < neighbour_cutoff2) {
                neighbours[n_neighbour++] = a;
            }
        }
        std::size_t n_accessible = 0;
        for (std::size_t j = 0; j < n_sphere; ++j) {
            const double sx = points[3 * j + 0] * radius + ax;
            const double sy = points[3 * j + 1] * radius + ay;
            const double sz = points[3 * j + 2] * radius + az;
            bool accessible = true;
            for (std::size_t k = 0; k < n_neighbour; ++k) {
                const std::size_t b = neighbours[k];
                const double dx = sx - xyz[3 * b + 0];
                const double dy = sy - xyz[3 * b + 1];
                const double dz = sz - xyz[3 * b + 2];
                if (dx * dx + dy * dy + dz * dz < occlusion2) {
                    accessible = false;
                    break;
                }
            }
            if (accessible) ++n_accessible;
        }
        asa[i] = c * static_cast<double>(n_accessible) * vdw[p] * vdw[p];
    }
    return asa;
}

std::vector<double> quenching_rate_per_frame(
        const std::vector<int>& collided, int n_frames,
        const std::vector<double>& k_quench) {
    std::vector<double> out(static_cast<std::size_t>(std::max(0, n_frames)), 0.0);
    if (n_frames <= 0 || k_quench.empty()) return out;
    const std::size_t n_atoms = k_quench.size();
    for (std::size_t f = 0; f < out.size(); ++f) {
        double total = 0.0;
        for (std::size_t a = 0; a < n_atoms; ++a) {
            if (collided[f * n_atoms + a]) total += k_quench[a];
        }
        out[f] = total;
    }
    return out;
}

IMPBFF_END_NAMESPACE
