/**
 * \file MolecularGraph.cpp
 * \brief Bond connectivity, and everything derived from it.
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/MolecularGraph.h>

#include <algorithm>
#include <set>

IMPBFF_BEGIN_NAMESPACE

MolecularGraph::MolecularGraph(const std::vector<std::pair<int, int> >& bonds) {
    std::map<int, std::set<int> > uniq;
    for (size_t i = 0; i < bonds.size(); ++i) {
        uniq[bonds[i].first].insert(bonds[i].second);
        uniq[bonds[i].second].insert(bonds[i].first);
    }
    for (std::map<int, std::set<int> >::const_iterator it = uniq.begin();
         it != uniq.end(); ++it) {
        adj_[it->first] = std::vector<int>(it->second.begin(), it->second.end());
    }
}

std::vector<int> MolecularGraph::get_nodes() const {
    std::vector<int> out;
    out.reserve(adj_.size());
    for (std::map<int, std::vector<int> >::const_iterator it = adj_.begin();
         it != adj_.end(); ++it) out.push_back(it->first);
    return out;
}

std::vector<int> MolecularGraph::get_neighbors(int node) const {
    std::map<int, std::vector<int> >::const_iterator it = adj_.find(node);
    return it == adj_.end() ? std::vector<int>() : it->second;
}

std::vector<std::vector<int> > MolecularGraph::get_angles() const {
    std::set<std::vector<int> > out;
    for (std::map<int, std::vector<int> >::const_iterator it = adj_.begin();
         it != adj_.end(); ++it) {
        const std::vector<int>& nn = it->second;          // already ascending
        for (size_t i = 0; i < nn.size(); ++i) {
            for (size_t j = i + 1; j < nn.size(); ++j) {
                std::vector<int> a(3);
                a[0] = nn[i]; a[1] = it->first; a[2] = nn[j];
                out.insert(a);
            }
        }
    }
    return std::vector<std::vector<int> >(out.begin(), out.end());
}

std::vector<std::vector<int> > MolecularGraph::get_dihedrals() const {
    std::set<std::vector<int> > out;
    for (std::map<int, std::vector<int> >::const_iterator it = adj_.begin();
         it != adj_.end(); ++it) {
        const int b = it->first;
        for (size_t k = 0; k < it->second.size(); ++k) {
            const int c = it->second[k];
            if (b > c) continue;                          // each central bond once
            const std::vector<int> left = get_neighbors(b);
            const std::vector<int> right = get_neighbors(c);
            for (size_t i = 0; i < left.size(); ++i) {
                const int a = left[i];
                if (a == c) continue;
                for (size_t j = 0; j < right.size(); ++j) {
                    const int d = right[j];
                    if (d == b) continue;
                    if (a == b || a == d || b == d || c == a || c == d) continue;
                    std::vector<int> f(4), r(4);
                    f[0] = a; f[1] = b; f[2] = c; f[3] = d;
                    r[0] = d; r[1] = c; r[2] = b; r[3] = a;
                    out.insert(std::min(f, r));           // a torsion and its reverse are one
                }
            }
        }
    }
    return std::vector<std::vector<int> >(out.begin(), out.end());
}

namespace {
//! Smallest of a cycle's rotations and reflections: one ring, one entry,
//! however the walk entered it.
std::vector<int> canonical(const std::vector<int>& cyc) {
    std::vector<int> best;
    const size_t n = cyc.size();
    for (int rev = 0; rev < 2; ++rev) {
        std::vector<int> base(cyc);
        if (rev) std::reverse(base.begin(), base.end());
        for (size_t i = 0; i < n; ++i) {
            std::vector<int> rot(base.begin() + i, base.end());
            rot.insert(rot.end(), base.begin(), base.begin() + i);
            if (best.empty() || rot < best) best = rot;
        }
    }
    return best;
}

void ring_dfs(const MolecularGraph& g, int start, int cur,
              std::set<int>& visited, std::vector<int>& path, int max_len,
              std::set<std::vector<int> >& out) {
    const std::vector<int> nbrs = g.get_neighbors(cur);
    for (size_t i = 0; i < nbrs.size(); ++i) {
        const int nb = nbrs[i];
        if (nb == start && path.size() >= 3) { out.insert(canonical(path)); continue; }
        // `nb < start` keeps each cycle to the walk beginning at its smallest
        // node; the cap is on the path length, not the recursion depth
        if (visited.count(nb) || nb < start || (int)path.size() >= max_len) continue;
        visited.insert(nb);
        path.push_back(nb);
        ring_dfs(g, start, nb, visited, path, max_len, out);
        path.pop_back();
        visited.erase(nb);
    }
}
}

std::vector<std::vector<int> > MolecularGraph::get_rings(int max_len) const {
    std::set<std::vector<int> > found;
    const std::vector<int> nodes = get_nodes();
    for (size_t i = 0; i < nodes.size(); ++i) {
        std::set<int> visited;
        visited.insert(nodes[i]);
        std::vector<int> path(1, nodes[i]);
        ring_dfs(*this, nodes[i], nodes[i], visited, path, max_len, found);
    }
    return std::vector<std::vector<int> >(found.begin(), found.end());
}

bool MolecularGraph::is_within_bonds(int a, int b, int max_depth) const {
    if (a == b) return true;
    std::set<int> seen;
    seen.insert(a);
    std::vector<int> frontier(1, a);
    for (int depth = 0; depth < max_depth && !frontier.empty(); ++depth) {
        std::vector<int> next;
        for (size_t i = 0; i < frontier.size(); ++i) {
            const std::vector<int> nbrs = get_neighbors(frontier[i]);
            for (size_t j = 0; j < nbrs.size(); ++j) {
                if (nbrs[j] == b) return true;
                if (seen.insert(nbrs[j]).second) next.push_back(nbrs[j]);
            }
        }
        frontier.swap(next);
    }
    return false;
}

IMPBFF_END_NAMESPACE
