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
#include <string>
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

    //! The nodes that lie on any ring of at most `max_len`.
    std::vector<int> get_ring_nodes(int max_len = 8) const;

    //! Improper quadruples about `centers`, by kind.
    /** The four kinds a cgdye template can declare, expanded against the bond
        graph. Each returns `(n1, centre, n2, n3)` quadruples.

        - `ring`   — a centre on a ring: two ring neighbours and a third
                     substituent, preferring one off the ring.
        - `pi`     — a three-coordinate C or N: its three neighbours.
        - `flat`   — a sulfur with at least three oxygens: three of them.
        - `orient` — a sulfur with exactly one carbon and at least two
                     oxygens: that carbon and two oxygens.

        `elements` and `atom_names` are parallel to `nodes`: the element and
        MOL2 atom name of each. `flat` and `orient` key off the *atom name*
        starting with S rather than the element, which is what the Python did
        and is not the same test for a MOL2 whose types are unreliable. */
    std::vector<std::vector<int> > expand_impropers(
            const std::string& kind,
            const std::vector<int>& centers,
            const std::vector<int>& nodes,
            const std::vector<std::string>& elements,
            const std::vector<std::string>& atom_names,
            int max_ring_len = 8) const;
};

IMPBFF_END_NAMESPACE

#endif //IMPBFF_MOLECULARGRAPH_H
