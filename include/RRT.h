/**
 *  \file IMP/bff/RRT.h
 *  \brief Rapidly-exploring random trees over a dye's degrees of freedom.
 *
 * An RRT grows a tree from a starting configuration toward random targets,
 * one bounded step at a time, keeping only the steps that do not clash. What
 * it is for here is finding a *path* rather than a minimum: how a dye gets
 * from one conformation to another without passing through the protein, which
 * a Metropolis walk answers only by accident.
 *
 * Two configuration spaces, two functions, one tree: rigid-body placements
 * (a translation and three Euler angles) and vectors of torsions (periodic,
 * so the distance is measured on the torus). They are kept apart rather than
 * generalised because their *metrics* differ in kind -- the rigid one adds a
 * rotation distance to a translation distance, which no single weighted norm
 * expresses.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_RRT_H
#define IMPBFF_RRT_H

#include <IMP/bff/bff_config.h>
// IMP_OBJECT_METHODS, used below; Base.h resolves the IMP macros in either
// configuration.
#include <IMP/bff/Base.h>

#include <IMP/Object.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Does this configuration clash?
/*! Subclass it -- in C++ or, through the director, in Python -- and answer
    true for a configuration the tree may not grow into. A planner asks it
    once per candidate step, so what it costs is what the planning costs.

    The callback may apply the configuration to a molecule as a side effect;
    a caller that wants the coordinates of a node builds them the same way. */
class IMPBFFEXPORT RRTCollision : public IMP::Object {
public:
    RRTCollision(std::string name = "RRTCollision%1%") : IMP::Object(name) {}
    virtual bool is_collision(const std::vector<double>& configuration) {
        return false;
    }
    IMP_OBJECT_METHODS(RRTCollision);
};

//! A tree of configurations, each reached from its parent by one step.
struct IMPBFFEXPORT RRTTree {
    //! Flat `n_nodes * n_dof`; node 0 is the start.
    std::vector<double> configurations;
    //! The node each was grown from; -1 for the root.
    std::vector<int> parents;
    int n_nodes, n_dof;
    //! The node that reached the goal, or -1 when none did.
    /*! -1 is not a failure to report as an error: a planner that did not
        reach the goal in the iterations it was given has still explored, and
        the tree is what it explored. */
    int goal_node;

    RRTTree() : n_nodes(0), n_dof(0), goal_node(-1) {}

    void get_configurations(double** out_view, int* n_out_view) const;
    //! The configuration of one node.
    std::vector<double> get_configuration(int node) const;
    //! The nodes from the root to \p node, in that order.
    std::vector<int> get_path_to_root(int node) const;

    IMP_SHOWABLE_INLINE(RRTTree,
                        out << "RRTTree(" << n_nodes << " nodes, "
                            << (goal_node >= 0 ? "goal reached" : "no goal")
                            << ")");
};
IMP_VALUES(RRTTree, RRTTrees);

//! Grow an RRT over a vector of periodic torsions.
/*!
    Distances are measured on the torus: an angle of \f$-\pi\f$ and one of
    \f$+\pi\f$ are the same place, and a planner that does not know that will
    refuse to step across the wrap.

    \param[in] n_dof how many torsions
    \param[in] n_iter how many growth attempts
    \param[in] step_size the most one step may move, radians
    \param[in] collision asked about every candidate; null never collides
    \param[in] start the root; empty starts at all zeros
    \param[in] goal what to aim at; empty explores without a goal
    \param[in] goal_bias how often to aim at the goal rather than at random
    \param[in] goal_tolerance how close counts as reaching it
    \param[in] seed the random seed, so a run repeats
*/
IMPBFFEXPORT RRTTree grow_torsion_rrt(
        int n_dof, int n_iter, double step_size, RRTCollision* collision,
        const std::vector<double>& start = std::vector<double>(),
        const std::vector<double>& goal = std::vector<double>(),
        double goal_bias = 0.1, double goal_tolerance = 0.1, int seed = 0);

//! Grow an RRT over rigid-body placements.
/*!
    A configuration is six numbers -- a translation and three fixed-xyz Euler
    angles -- and the distance between two is the translation distance plus
    \p rot_weight times the rotation distance. Two distances added rather than
    one norm over six numbers, because Ångström and radian are not the same
    thing and pretending otherwise makes the weight meaningless.

    It plans in the six numbers, not in `Transformation3D`: a placement is the
    caller's way of holding one, and #transformation_to_configuration converts.
    (IMP's SWIG layer will not pass a value type like `Transformation3D`
    through a pointer either, so an optional goal has to be a vector anyway --
    and an empty vector says "no goal" more plainly than a null pointer.)

    \param[in] start the root placement, six numbers
    \param[in] bounds_low,bounds_high six each: the box random targets come from
    \param[in] n_iter,step_size,collision,goal_bias,goal_tolerance,seed as above
    \param[in] goal the placement to reach; empty explores without one
    \param[in] rot_weight what a radian is worth in Ångström
*/
IMPBFFEXPORT RRTTree grow_rigid_body_rrt(
        const std::vector<double>& start,
        const std::vector<double>& bounds_low,
        const std::vector<double>& bounds_high, int n_iter, double step_size,
        RRTCollision* collision,
        const std::vector<double>& goal = std::vector<double>(),
        double goal_bias = 0.1, double goal_tolerance = 1.0,
        double rot_weight = 0.25, int seed = 0);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_RRT_H
