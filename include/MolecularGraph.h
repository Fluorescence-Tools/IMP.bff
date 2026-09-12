/**
 * \file IMP/bff/MolecularGraph.h
 * \brief Bond connectivity, and everything derived from it.
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_MOLECULARGRAPH_H
#define IMPBFF_MOLECULARGRAPH_H

#include <IMP/bff/bff_config.h>

#include <IMP/bff/IMPCompatibility.h>

#include <map>
#include <string>
#include <utility>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Bond connectivity over node ids, and the terms derived from it.
/** Angles, torsions and rings are all functions of which node is bonded to
    which, and were computed in four places in Python -- `cgprobe.topology` twice
    over MOL2 serials, `cgprobe.sim` and `scoring` over site ids. One
    implementation, in C++, over whatever integer ids the caller has.

    Node ids need not be dense or start anywhere in particular: MOL2 serials
    start at 1 and can be sparse, which is why the adjacency is a map. */
//! Which side of a bond or angle turns, and what goes with it.
/** The anchor is the atom a sampler holds still -- the attachment point of a
    dye, say. `is_rotatable` is false when both sides reach it, which is what
    a ring bond does. */
struct IMPBFFEXPORT GraphRotor {
    //! The end that stays put, and the end whose side turns.
    int fixed;
    int moving;
    //! The nodes that travel with `moving`, ascending.
    std::vector<int> moving_nodes;
    bool is_rotatable;

    GraphRotor() : fixed(-1), moving(-1), is_rotatable(false) {}

    IMP_SHOWABLE_INLINE(GraphRotor,
                        out << "GraphRotor(" << fixed << "->" << moving << ", "
                            << moving_nodes.size() << " moving)");
};
IMP_VALUES(GraphRotor, GraphRotors);

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

    //! The nodes reachable from `start` when the edge `(u, v)` is cut.
    /** Ascending. The edge is cut in both directions; every other edge stays.
        This is what says which half of a molecule a rotatable bond moves. */
    std::vector<int> get_component_without_edge(int start, int u, int v) const;

    //! Which side of the bond `(a, b)` turns when `anchor` may not move.
    /** \return `(fixed, moving)` and the moving side's nodes; `is_rotatable`
        is false when the anchor is reachable from both ends -- a bond inside
        a ring, which no torsion can turn. */
    GraphRotor get_bond_rotor(int a, int b, int anchor) const;

    //! The same for the angle `a-b-c`: the axis sits at `b`.
    GraphRotor get_angle_rotor(int a, int b, int c, int anchor) const;

    //! Heavy atoms lying on a ring whose size is one of `ring_sizes`.
    /** `elements` is parallel to `nodes`, as in #expand_impropers. Hydrogens
        are dropped: a ring's hydrogens ride with it and are not what a
        rotamer generator means by "ring atom". */
    std::vector<int> get_ring_atoms(
            const std::vector<int>& nodes,
            const std::vector<std::string>& elements,
            const std::vector<int>& ring_sizes) const;

    //! Improper quadruples about `centers`, by kind.
    /** The four kinds a cgprobe template can declare, expanded against the bond
        graph. Each returns `(n1, centre, n2, n3)` quadruples.

        - `ring`   — a centre on a ring: two ring neighbours and a third
                     substituent, preferring one off the ring.
        - `pi`     — a three-coordinate C or N: its three neighbours.
        - `flat`   — a sulfur with at least three oxygens: three of them.
        - `orient` — a sulfur with exactly one carbon and at least two
                     oxygens: that carbon and two oxygens.

        `elements` and `atom_names` are parallel to `nodes`: the element and
        MOL2 atom name of each. `flat` and `orient` key off the *atom name*
        starting with S rather than the element -- deliberately, and not the
        same test for a MOL2 whose types are unreliable. */
    std::vector<std::vector<int> > expand_impropers(
            const std::string& kind,
            const std::vector<int>& centers,
            const std::vector<int>& nodes,
            const std::vector<std::string>& elements,
            const std::vector<std::string>& atom_names,
            int max_ring_len = 8) const;
};

//! The same graph over a caller's labels rather than over integers.
/*!
    #MolecularGraph keys its nodes by `int`, which is what a MOL2 serial is;
    site ids are strings. Bridging that in a caller means numbering the
    labels, calling the graph and mapping the answers back -- and doing it
    twice, since a first-appearance numbering and a sorted one give different
    orders out. The numbering
    belongs with the graph.

    Nodes are numbered in **sorted label order**, so an angle comes back with
    its ends ascending by label and a torsion in whichever direction is
    smaller, which is what #MolecularGraph promises about integers.
*/
class IMPBFFEXPORT LabelledGraph {
    MolecularGraph graph_;
    std::vector<std::string> labels_;
    std::map<std::string, int> index_;

    std::vector<std::string> labelled(const std::vector<int>& nodes) const;

public:
    LabelledGraph() {}

    //! From bonded label pairs. Order within a pair does not matter.
    explicit LabelledGraph(
            const std::vector<std::pair<std::string, std::string> >& bonds);

    //! The nodes, ascending by label.
    std::vector<std::string> get_nodes() const { return labels_; }

    //! Neighbours of `node`, ascending. Empty for a label the graph lacks.
    std::vector<std::string> get_neighbors(const std::string& node) const;

    //! Every `(a, b, c)` with `b` bonded to both, `a < c` by label.
    std::vector<std::vector<std::string> > get_angles() const;

    //! Every proper torsion, once, in whichever direction is smaller.
    std::vector<std::vector<std::string> > get_dihedrals() const;

    //! Simple cycles of at most `max_len` nodes.
    std::vector<std::vector<std::string> > get_rings(int max_len = 8) const;

    //! The integer graph underneath, for callers that want its other methods.
    const MolecularGraph& get_graph() const { return graph_; }
};

IMPBFF_END_NAMESPACE

#endif //IMPBFF_MOLECULARGRAPH_H
