/**
 * \file IMP/bff/MolecularGraph.h
 * \brief Bond connectivity, and everything derived from it.
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_MOLECULARGRAPH_H
#define IMPBFF_MOLECULARGRAPH_H

#include <IMP/bff/bff_config.h>

#include <map>
#include <utility>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Bond connectivity over node ids, and the terms derived from it.
/** Angles, torsions and rings are all functions of which node is bonded to
    which, and were computed in four places in Python -- `cgdye.topology` twice
    over MOL2 serials, `cgdye.sim` and `scoring` over site ids. One
    implementation, in C++, over whatever integer ids the caller has.

    Node ids need not be dense or start anywhere in particular: MOL2 serials
    start at 1 and can be sparse, which is why the adjacency is a map. */
class IMPBFFEXPORT MolecularGraph {
    std::map<int, std::vector<int> > adj_;

public:
    MolecularGraph() {}

    //! From bonded pairs. Order within a pair does not matter.
    explicit MolecularGraph(const std::vector<std::pair<int, int> >& bonds);

    //! The nodes, ascending.
    std::vector<int> get_nodes() const;

    //! Neighbours of `node`, ascending. Empty for an unknown node.
    std::vector<int> get_neighbors(int node) const;

    //! Every `(a, b, c)` with `b` bonded to both, `a < c`.
    std::vector<std::vector<int> > get_angles() const;

    //! Every proper torsion `(a, b, c, d)` over four distinct nodes.
    /** Each is emitted once, in whichever of `(a,b,c,d)` and `(d,c,b,a)` is
        smaller, so a torsion and its reverse are one entry. */
    std::vector<std::vector<int> > get_dihedrals() const;

    //! Simple cycles of at most `max_len` nodes, each in canonical form.
    std::vector<std::vector<int> > get_rings(int max_len = 8) const;

    //! Whether `a` reaches `b` in at most `max_depth` bonds.
    bool is_within_bonds(int a, int b, int max_depth) const;
};

IMPBFF_END_NAMESPACE

#endif //IMPBFF_MOLECULARGRAPH_H
