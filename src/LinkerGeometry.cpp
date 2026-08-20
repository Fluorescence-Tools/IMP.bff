/**
 * \file LinkerGeometry.cpp
 * \brief Applying a torsion/angle configuration to a linker.
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/LinkerGeometry.h>

#include <cmath>
#include <stdexcept>

IMPBFF_BEGIN_NAMESPACE

namespace {

struct V3 {
    double x, y, z;
    V3() : x(0), y(0), z(0) {}
    V3(double a, double b, double c) : x(a), y(b), z(c) {}
    V3 operator-(const V3& o) const { return V3(x - o.x, y - o.y, z - o.z); }
    double norm() const { return std::sqrt(x * x + y * y + z * z); }
};

V3 cross(const V3& a, const V3& b) {
    return V3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

V3 at(const std::vector<double>& c, int i) {
    return V3(c[3 * i], c[3 * i + 1], c[3 * i + 2]);
}

//! Rodrigues rotation of `p` about the axis through `origin` along `axis`.
V3 rotate_about(const V3& p, const V3& origin, const V3& axis, double angle) {
    const double len = axis.norm();
    const V3 u(axis.x / len, axis.y / len, axis.z / len);
    const V3 d = p - origin;
    const double c = std::cos(angle), s = std::sin(angle);
    const V3 k = cross(u, d);
    const double dot = u.x * d.x + u.y * d.y + u.z * d.z;
    return V3(origin.x + d.x * c + k.x * s + u.x * dot * (1.0 - c),
              origin.y + d.y * c + k.y * s + u.y * dot * (1.0 - c),
              origin.z + d.z * c + k.z * s + u.z * dot * (1.0 - c));
}

}  // namespace

void LinkerGeometry::add_torsion(int fixed, int moving, const std::vector<int>& moves) {
    torsion_fixed_.push_back(fixed);
    torsion_moving_.push_back(moving);
    torsion_sets_.push_back(moves);
}

void LinkerGeometry::add_angle(int b, int c, int a, const std::vector<int>& moves) {
    angle_b_.push_back(b);
    angle_c_.push_back(c);
    angle_a_.push_back(a);
    angle_sets_.push_back(moves);
}

std::vector<double> LinkerGeometry::apply(const std::vector<double>& config) const {
    const size_t n_torsions = torsion_fixed_.size();
    const size_t n_angles = angle_b_.size();
    if (config.size() != n_torsions + n_angles) {
        throw std::invalid_argument(
            "config must hold one value per torsion followed by one per angle");
    }
    std::vector<double> c(base_);

    for (size_t t = 0; t < n_torsions; ++t) {
        const V3 cf = at(c, torsion_fixed_[t]);
        const V3 cm = at(c, torsion_moving_[t]);
        const V3 axis = cm - cf;
        // two atoms at the same point define no rotation
        if (axis.norm() < 1e-8) continue;
        const double angle = config[t];
        const std::vector<int>& moves = torsion_sets_[t];
        for (size_t i = 0; i < moves.size(); ++i) {
            const V3 p = rotate_about(at(c, moves[i]), cf, axis, angle);
            c[3 * moves[i]] = p.x; c[3 * moves[i] + 1] = p.y; c[3 * moves[i] + 2] = p.z;
        }
    }

    for (size_t a = 0; a < n_angles; ++a) {
        const V3 cb = at(c, angle_b_[a]);
        const V3 cc = at(c, angle_c_[a]);
        const V3 ca = at(c, angle_a_[a]);
        const V3 axis = cross(ca - cb, cc - cb);      // normal of the a-b-c plane
        if (axis.norm() < 1e-8) continue;
        const double angle = config[n_torsions + a];
        const std::vector<int>& moves = angle_sets_[a];
        for (size_t i = 0; i < moves.size(); ++i) {
            const V3 p = rotate_about(at(c, moves[i]), cb, axis, angle);
            c[3 * moves[i]] = p.x; c[3 * moves[i] + 1] = p.y; c[3 * moves[i] + 2] = p.z;
        }
    }
    return c;
}

IMPBFF_END_NAMESPACE
