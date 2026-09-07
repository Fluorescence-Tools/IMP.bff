/**
 * \file RRT.cpp
 * \brief Rapidly-exploring random trees over a dye's degrees of freedom.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/RRT.h>

#include <IMP/bff/internal/OutputView.h>

#include <IMP/algebra/Rotation3D.h>
#include <IMP/constants.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>

#include <algorithm>
#include <cmath>

IMPBFF_BEGIN_NAMESPACE

void RRTTree::get_configurations(double** out_view, int* n_out_view) const {
    internal::copy_to_view(configurations, out_view, n_out_view);
}

std::vector<double> RRTTree::get_configuration(int node) const {
    if (node < 0 || node >= n_nodes || n_dof == 0) return std::vector<double>();
    const std::size_t offset = static_cast<std::size_t>(node) * n_dof;
    return std::vector<double>(configurations.begin() + offset,
                               configurations.begin() + offset + n_dof);
}

std::vector<int> RRTTree::get_path_to_root(int node) const {
    std::vector<int> path;
    int current = node;
    while (current >= 0 && current < static_cast<int>(parents.size())) {
        path.push_back(current);
        current = parents[current];
    }
    std::reverse(path.begin(), path.end());
    return path;
}

namespace {

//! An angle in [-pi, pi): the torus has no other representative.
double wrap_angle(double a) {
    const double two_pi = 2.0 * IMP::PI;
    return std::fmod(std::fmod(a + IMP::PI, two_pi) + two_pi, two_pi) - IMP::PI;
}

double torsion_distance(const std::vector<double>& a,
                        const std::vector<double>& b) {
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size() && i < b.size(); ++i) {
        const double d = wrap_angle(a[i] - b[i]);
        sum += d * d;
    }
    return std::sqrt(sum);
}

//! The rigid metric: a translation distance plus a weighted rotation distance.
double rigid_distance(const std::vector<double>& a,
                      const std::vector<double>& b, double rot_weight) {
    double dt = 0.0, dr = 0.0;
    for (int i = 0; i < 3; ++i) dt += (a[i] - b[i]) * (a[i] - b[i]);
    for (int i = 3; i < 6; ++i) dr += (a[i] - b[i]) * (a[i] - b[i]);
    return std::sqrt(dt) + rot_weight * std::sqrt(dr);
}

//! The node of \p tree closest to \p target under \p distance.
template <typename Distance>
int nearest_node(const RRTTree& tree, const std::vector<double>& target,
                 Distance distance) {
    int best = 0;
    double best_distance = std::numeric_limits<double>::infinity();
    for (int i = 0; i < tree.n_nodes; ++i) {
        const double d = distance(tree.get_configuration(i), target);
        if (d < best_distance) {
            best_distance = d;
            best = i;
        }
    }
    return best;
}

void append_node(RRTTree& tree, const std::vector<double>& configuration,
                 int parent) {
    tree.configurations.insert(tree.configurations.end(),
                               configuration.begin(), configuration.end());
    tree.parents.push_back(parent);
    ++tree.n_nodes;
}

}  // namespace

RRTTree grow_torsion_rrt(int n_dof, int n_iter, double step_size,
                         RRTCollision* collision,
                         const std::vector<double>& start,
                         const std::vector<double>& goal, double goal_bias,
                         double goal_tolerance, int seed) {
    RRTTree tree;
    tree.n_dof = std::max(0, n_dof);
    if (tree.n_dof == 0) return tree;

    std::vector<double> root = start;
    root.resize(tree.n_dof, 0.0);
    append_node(tree, root, -1);

    boost::mt19937 rng(static_cast<boost::uint32_t>(seed));
    boost::uniform_real<double> unit(0.0, 1.0);
    const bool has_goal = static_cast<int>(goal.size()) == tree.n_dof;

    for (int step = 0; step < n_iter; ++step) {
        std::vector<double> target(tree.n_dof);
        if (has_goal && unit(rng) < goal_bias) {
            target = goal;
        } else {
            for (int i = 0; i < tree.n_dof; ++i) {
                target[i] = -IMP::PI + 2.0 * IMP::PI * unit(rng);
            }
        }

        const int near = nearest_node(tree, target, torsion_distance);
        const std::vector<double> from = tree.get_configuration(near);
        const double d = torsion_distance(from, target);
        std::vector<double> candidate(tree.n_dof);
        const double alpha = d > 0.0 ? std::min(1.0, step_size / d) : 0.0;
        for (int i = 0; i < tree.n_dof; ++i) {
            candidate[i] =
                    wrap_angle(from[i] + alpha * wrap_angle(target[i] - from[i]));
        }

        if (collision != NULL && collision->is_collision(candidate)) continue;
        append_node(tree, candidate, near);
        if (has_goal &&
            torsion_distance(candidate, goal) <= goal_tolerance) {
            tree.goal_node = tree.n_nodes - 1;
            break;
        }
    }
    return tree;
}

std::vector<double> transformation_to_configuration(
        const IMP::algebra::Transformation3D& transformation) {
    const IMP::algebra::Vector3D t = transformation.get_translation();
    const IMP::algebra::FixedXYZ e =
            IMP::algebra::get_fixed_xyz_from_rotation(
                    transformation.get_rotation());
    std::vector<double> out;
    out.push_back(t[0]);
    out.push_back(t[1]);
    out.push_back(t[2]);
    out.push_back(e.get_x());
    out.push_back(e.get_y());
    out.push_back(e.get_z());
    return out;
}

IMP::algebra::Transformation3D configuration_to_transformation(
        const std::vector<double>& c) {
    if (c.size() < 6) return IMP::algebra::Transformation3D();
    return IMP::algebra::Transformation3D(
            IMP::algebra::get_rotation_from_fixed_xyz(c[3], c[4], c[5]),
            IMP::algebra::Vector3D(c[0], c[1], c[2]));
}

RRTTree grow_rigid_body_rrt(const std::vector<double>& start,
                            const std::vector<double>& bounds_low,
                            const std::vector<double>& bounds_high, int n_iter,
                            double step_size, RRTCollision* collision,
                            const std::vector<double>& goal, double goal_bias,
                            double goal_tolerance, double rot_weight,
                            int seed) {
    RRTTree tree;
    tree.n_dof = 6;
    if (bounds_low.size() < 6 || bounds_high.size() < 6) {
        IMP_THROW("grow_rigid_body_rrt needs six bounds, got "
                          << bounds_low.size() << " and " << bounds_high.size(),
                  ValueException);
    }
    std::vector<double> root = start;
    root.resize(6, 0.0);
    append_node(tree, root, -1);

    boost::mt19937 rng(static_cast<boost::uint32_t>(seed));
    boost::uniform_real<double> unit(0.0, 1.0);
    const std::vector<double> goal_configuration =
            goal.size() >= 6 ? goal : std::vector<double>();

    for (int step = 0; step < n_iter; ++step) {
        std::vector<double> target(6);
        if (!goal_configuration.empty() && unit(rng) < goal_bias) {
            target = goal_configuration;
        } else {
            for (int i = 0; i < 6; ++i) {
                target[i] = bounds_low[i] +
                            (bounds_high[i] - bounds_low[i]) * unit(rng);
            }
        }

        const int near = nearest_node(
                tree, target, [rot_weight](const std::vector<double>& a,
                                           const std::vector<double>& b) {
                    return rigid_distance(a, b, rot_weight);
                });
        const std::vector<double> from = tree.get_configuration(near);
        const double d = rigid_distance(from, target, rot_weight);
        std::vector<double> candidate(6);
        if (d <= step_size) {
            candidate = target;
        } else {
            const double alpha = d > 0.0 ? step_size / d : 0.0;
            for (int i = 0; i < 6; ++i) {
                candidate[i] = from[i] + alpha * (target[i] - from[i]);
            }
        }

        if (collision != NULL && collision->is_collision(candidate)) continue;
        append_node(tree, candidate, near);
        if (!goal_configuration.empty() &&
            rigid_distance(candidate, goal_configuration, rot_weight) <=
                    goal_tolerance) {
            tree.goal_node = tree.n_nodes - 1;
            break;
        }
    }
    return tree;
}

IMPBFF_END_NAMESPACE
