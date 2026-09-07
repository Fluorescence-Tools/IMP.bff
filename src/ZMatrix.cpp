/**
 * \file ZMatrix.cpp
 * \brief Internal coordinates: measurement, placement, and the Z-matrix.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/ZMatrix.h>

#include <IMP/algebra/Vector3D.h>

#include <algorithm>
#include <cmath>
#include <deque>
#include <map>
#include <queue>
#include <stdexcept>

IMPBFF_BEGIN_NAMESPACE

//! Covalent radii (A) by element symbol; unknown -> carbon's.
double covalent_radius(const std::string& e) {
    if (e == "H") return 0.37;
    if (e == "C") return 0.77;
    if (e == "N") return 0.75;
    if (e == "O") return 0.73;
    if (e == "P") return 1.10;
    if (e == "S") return 1.02;
    return 0.77;
}

double bond_tolerance() { return 0.35; }

namespace {
constexpr double kBondTol = 0.35;
}  // namespace

double dihedral_deg(const algebra::Vector3D& a, const algebra::Vector3D& b,
                    const algebra::Vector3D& c, const algebra::Vector3D& d) {
    // praxeolitic form: project the outer bonds onto the plane through the
    // b-c axis and take the signed angle between the projections.
    const algebra::Vector3D b0 = a - b;
    algebra::Vector3D b1 = c - b;
    const algebra::Vector3D b2 = d - c;
    const double nb1 = b1.get_magnitude();
    if (nb1 < 1e-12) return 0.0;
    b1 /= nb1;
    const algebra::Vector3D v = b0 - b1 * (b0 * b1);
    const algebra::Vector3D w = b2 - b1 * (b2 * b1);
    const algebra::Vector3D cr = algebra::get_vector_product(b1, v);
    return std::atan2(cr * w, v * w) * 180.0 / M_PI;
}

double bond_angle_rad(const algebra::Vector3D& a, const algebra::Vector3D& b,
                      const algebra::Vector3D& c) {
    const algebra::Vector3D u = a - b, w = c - b;
    const double nu = u.get_magnitude(), nw = w.get_magnitude();
    // Tetrahedral for a collapsed arm; see the header for why not zero.
    if (nu < 1e-12 || nw < 1e-12) return 1.9106332362490186;
    double cosv = (u * w) / (nu * nw);
    cosv = std::max(-1.0, std::min(1.0, cosv));
    return std::acos(cosv);
}

double bond_angle_deg(const algebra::Vector3D& a, const algebra::Vector3D& b,
                      const algebra::Vector3D& c) {
    return bond_angle_rad(a, b, c) * 180.0 / M_PI;
}

algebra::Vector3D internal2cartesian(
        const algebra::Vector3D& ggp, const algebra::Vector3D& gp,
        const algebra::Vector3D& parent,
        double bond_length, double bond_angle_deg_, double phi_deg) {
    // FASPR's basis: e along gp->parent; v from gp toward ggp made
    // perpendicular to e (the azimuth origin); q = e x v. The new atom sits
    // at parent + r * (-cos(theta) e + sin(theta) (cos(phi) v + sin(phi) q)),
    // which is the exact inverse of dihedral_deg on (ggp, gp, parent, new).
    algebra::Vector3D e = parent - gp;
    const double ne = e.get_magnitude();
    if (ne < 1e-12) return parent;
    e /= ne;
    algebra::Vector3D v = ggp - gp;  // gp toward ggp (the azimuth origin)
    v = v - e * (v * e);
    double nv = v.get_magnitude();
    if (nv < 1e-10) {
        // degenerate plane (ggp on the gp->parent axis): any perpendicular
        v = algebra::Vector3D(1, 0, 0) - e * e[0];
        nv = v.get_magnitude();
        if (nv < 1e-10) {
            v = algebra::Vector3D(0, 1, 0) - e * e[1];
            nv = v.get_magnitude();
        }
    }
    v /= nv;
    const algebra::Vector3D q = algebra::get_vector_product(e, v);
    const double th = bond_angle_deg_ * M_PI / 180.0;
    const double ph = phi_deg * M_PI / 180.0;
    const algebra::Vector3D dir =
            e * (-std::cos(th)) +
            (v * std::cos(ph) + q * std::sin(ph)) * std::sin(th);
    return parent + dir * bond_length;
}

std::vector<std::pair<int, int> > perceive_bonds(
        const std::vector<algebra::Vector3D>& xyz,
        const std::vector<std::string>& elements) {
    std::vector<std::pair<int, int> > out;
    const int n = static_cast<int>(xyz.size());
    if (static_cast<int>(elements.size()) != n) return out;
    std::vector<double> radii(n);
    for (int i = 0; i < n; ++i) radii[i] = covalent_radius(elements[i]);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const double d = (xyz[i] - xyz[j]).get_magnitude();
            if (d > 0.1 && d < radii[i] + radii[j] + kBondTol) {
                out.push_back(std::make_pair(i, j));
            }
        }
    }
    return out;
}

void ZMatrix::set_template(const std::vector<algebra::Vector3D>& xyz,
                           const std::vector<std::string>& elements,
                           int root) {
    set_template_bonds(xyz, perceive_bonds(xyz, elements), root);
}

void ZMatrix::set_template_bonds(
        const std::vector<algebra::Vector3D>& xyz,
        const std::vector<std::pair<int, int> >& bonds, int root) {
    n_ = static_cast<unsigned int>(xyz.size());
    template_ = xyz;
    bonds_ = bonds;
    adj_.assign(n_, {});
    for (unsigned int k = 0; k < bonds_.size(); ++k) {
        const int i = bonds_[k].first, j = bonds_[k].second;
        adj_[i].push_back(j);
        adj_[j].push_back(i);
    }
    for (unsigned int i = 0; i < n_; ++i) {
        std::sort(adj_[i].begin(), adj_[i].end());
    }
    if (root < 0 || root >= static_cast<int>(n_)) root = 0;
    build_(root);
}

void ZMatrix::build_(int root) {
    // breadth-first spanning tree; order = (depth, index) so a atom's three
    // references are always placed before it.
    std::vector<int> parent(n_, -1), depth(n_, -1), order;
    order.reserve(n_);
    std::deque<int> queue;
    parent[root] = root;  // sentinel: itself, so "no parent" stays -1-free
    depth[root] = 0;
    queue.push_back(root);
    while (!queue.empty()) {
        const int u = queue.front();
        queue.pop_front();
        order.push_back(u);
        for (unsigned int k = 0; k < adj_[u].size(); ++k) {
            const int v = adj_[u][k];
            if (depth[v] == -1) {
                depth[v] = depth[u] + 1;
                parent[v] = u;
                queue.push_back(v);
            }
        }
    }
    if (static_cast<int>(order.size()) != static_cast<int>(n_)) {
        throw std::invalid_argument(
                "ZMatrix: template is not a connected molecule "
                "(a spanning tree covers only " +
                std::to_string(order.size()) + " of " +
                std::to_string(n_) + " atoms)");
    }
    // deterministic ordering: (depth, atom index)
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        if (depth[a] != depth[b]) return depth[a] < depth[b];
        return a < b;
    });
    // base = first three atoms in order; rows for the rest whose ancestor
    // chain reaches depth >= 3; atoms with a shorter chain join the base.
    std::vector<bool> is_base(n_, false);
    for (int k = 0; k < 3 && k < static_cast<int>(order.size()); ++k) {
        is_base[order[k]] = true;
    }
    rows_flat_.clear();
    base_.clear();
    std::vector<int> row_of(n_, -1);
    for (unsigned int k = 0; k < order.size(); ++k) {
        const int i = order[k];
        if (is_base[i]) {
            base_.push_back(i);
            continue;
        }
        const int p1 = parent[i];
        const int p2 = (p1 == root) ? -1 : parent[p1];
        const int p3 = (p2 == -1 || p2 == root) ? -1 : parent[p2];
        if (p1 == -1 || p2 == -1 || p3 == -1) {
            is_base[i] = true;
            base_.push_back(i);
            continue;
        }
        row_of[i] = static_cast<int>(rows_flat_.size() / 4);
        rows_flat_.push_back(i);
        rows_flat_.push_back(p1);
        rows_flat_.push_back(p2);
        rows_flat_.push_back(p3);
    }
    // template internals per row
    const unsigned int nrow = rows_flat_.size() / 4;
    lens_.resize(nrow);
    angs_.resize(nrow);
    chi_template_.resize(nrow);
    for (unsigned int r = 0; r < nrow; ++r) {
        const int i = rows_flat_[4 * r], p1 = rows_flat_[4 * r + 1],
                  p2 = rows_flat_[4 * r + 2], p3 = rows_flat_[4 * r + 3];
        lens_[r] = (template_[i] - template_[p1]).get_magnitude();
        angs_[r] = bond_angle_deg(template_[p2], template_[p1], template_[i]);
        chi_template_[r] = dihedral_deg(template_[p3], template_[p2],
                                        template_[p1], template_[i]);
    }
}

std::vector<double> ZMatrix::encode(
        const std::vector<algebra::Vector3D>& xyz) const {
    const unsigned int nrow = rows_flat_.size() / 4;
    std::vector<double> chi(nrow);
    if (xyz.size() != n_) return chi;
    for (unsigned int r = 0; r < nrow; ++r) {
        const int i = rows_flat_[4 * r], p1 = rows_flat_[4 * r + 1],
                  p2 = rows_flat_[4 * r + 2], p3 = rows_flat_[4 * r + 3];
        chi[r] = dihedral_deg(xyz[p3], xyz[p2], xyz[p1], xyz[i]);
    }
    return chi;
}

std::vector<algebra::Vector3D> ZMatrix::decode(
        const std::vector<double>& chi,
        const std::vector<algebra::Vector3D>* base) const {
    const unsigned int nrow = rows_flat_.size() / 4;
    std::vector<algebra::Vector3D> out = template_;
    if (base != nullptr && base->size() == n_) out = *base;
    if (chi.size() != nrow) return out;
    for (unsigned int r = 0; r < nrow; ++r) {
        const int i = rows_flat_[4 * r], p1 = rows_flat_[4 * r + 1],
                  p2 = rows_flat_[4 * r + 2], p3 = rows_flat_[4 * r + 3];
        out[i] = internal2cartesian(out[p3], out[p2], out[p1],
                                    lens_[r], angs_[r], chi[r]);
    }
    return out;
}

void ZMatrix::frame_internals(const std::vector<algebra::Vector3D>& xyz,
                              std::vector<double>& lengths,
                              std::vector<double>& angles,
                              std::vector<double>& chi) const {
    const unsigned int nrow = rows_flat_.size() / 4;
    lengths.assign(nrow, 0.0);
    angles.assign(nrow, 0.0);
    chi.assign(nrow, 0.0);
    if (xyz.size() != n_) return;
    for (unsigned int r = 0; r < nrow; ++r) {
        const int i = rows_flat_[4 * r], p1 = rows_flat_[4 * r + 1],
                  p2 = rows_flat_[4 * r + 2], p3 = rows_flat_[4 * r + 3];
        lengths[r] = (xyz[i] - xyz[p1]).get_magnitude();
        angles[r] = bond_angle_deg(xyz[p2], xyz[p1], xyz[i]);
        chi[r] = dihedral_deg(xyz[p3], xyz[p2], xyz[p1], xyz[i]);
    }
}

std::vector<algebra::Vector3D> ZMatrix::decode_internals(
        const std::vector<double>& lengths, const std::vector<double>& angles,
        const std::vector<double>& chi,
        const std::vector<algebra::Vector3D>* base) const {
    const unsigned int nrow = rows_flat_.size() / 4;
    std::vector<algebra::Vector3D> out = template_;
    if (base != nullptr && base->size() == n_) out = *base;
    if (chi.size() != nrow || lengths.size() != nrow || angles.size() != nrow) {
        return out;
    }
    for (unsigned int r = 0; r < nrow; ++r) {
        const int i = rows_flat_[4 * r], p1 = rows_flat_[4 * r + 1],
                  p2 = rows_flat_[4 * r + 2], p3 = rows_flat_[4 * r + 3];
        out[i] = internal2cartesian(out[p3], out[p2], out[p1],
                                    lengths[r], angles[r], chi[r]);
    }
    return out;
}

IMPBFF_END_NAMESPACE
