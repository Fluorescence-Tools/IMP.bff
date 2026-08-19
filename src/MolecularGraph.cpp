/**
 * \file MolecularGraph.cpp
 * \brief Bond connectivity, and everything derived from it.
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/MolecularGraph.h>

#include <algorithm>
#include <cctype>
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

std::vector<int> MolecularGraph::get_ring_nodes(int max_len) const {
    const std::vector<std::vector<int> > rings = get_rings(max_len);
    std::set<int> on_ring;
    for (size_t i = 0; i < rings.size(); ++i)
        on_ring.insert(rings[i].begin(), rings[i].end());
    return std::vector<int>(on_ring.begin(), on_ring.end());
}

std::vector<std::vector<int> > MolecularGraph::expand_impropers(
        const std::string& kind,
        const std::vector<int>& centers,
        const std::vector<int>& nodes,
        const std::vector<std::string>& elements,
        const std::vector<std::string>& atom_names,
        int max_ring_len) const {
    std::map<int, std::string> element_of, name_of;
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (i < elements.size()) element_of[nodes[i]] = elements[i];
        if (i < atom_names.size()) name_of[nodes[i]] = atom_names[i];
    }
    const std::set<int> known(nodes.begin(), nodes.end());

    std::set<int> on_ring;
    if (kind == "ring") {
        const std::vector<int> r = get_ring_nodes(max_ring_len);
        on_ring.insert(r.begin(), r.end());
    }

    std::vector<int> ordered(centers);
    std::sort(ordered.begin(), ordered.end());
    ordered.erase(std::unique(ordered.begin(), ordered.end()), ordered.end());

    // `seen` dedups on (centre, {n1,n2,n3}) for the two kinds that can produce
    // the same quadruple twice; `flat` and `orient` cannot, and the Python did
    // not dedup them either.
    std::set<std::pair<int, std::vector<int> > > seen;
    std::vector<std::vector<int> > out;

    for (size_t ci = 0; ci < ordered.size(); ++ci) {
        const int centre = ordered[ci];
        std::vector<int> nbrs;
        const std::vector<int> all = get_neighbors(centre);
        for (size_t i = 0; i < all.size(); ++i)
            if (known.count(all[i])) nbrs.push_back(all[i]);   // already ascending

        std::vector<int> quad;
        if (kind == "ring") {
            if (!on_ring.count(centre) || nbrs.size() < 3) continue;
            std::vector<int> ring_nbrs;
            for (size_t i = 0; i < nbrs.size(); ++i)
                if (on_ring.count(nbrs[i])) ring_nbrs.push_back(nbrs[i]);
            if (ring_nbrs.size() < 2) continue;
            const int n1 = ring_nbrs[0], n2 = ring_nbrs[1];
            std::vector<int> rest;
            for (size_t i = 0; i < nbrs.size(); ++i)
                if (nbrs[i] != n1 && nbrs[i] != n2) rest.push_back(nbrs[i]);
            if (rest.empty()) continue;
            int n3 = rest[0];
            for (size_t i = 0; i < rest.size(); ++i)      // prefer one off the ring
                if (!on_ring.count(rest[i])) { n3 = rest[i]; break; }
            quad.push_back(n1); quad.push_back(centre);
            quad.push_back(n2); quad.push_back(n3);
        } else if (kind == "pi") {
            std::map<int, std::string>::const_iterator e = element_of.find(centre);
            if (e == element_of.end() || (e->second != "C" && e->second != "N")) continue;
            if (nbrs.size() != 3) continue;
            quad.push_back(nbrs[0]); quad.push_back(centre);
            quad.push_back(nbrs[1]); quad.push_back(nbrs[2]);
        } else if (kind == "flat" || kind == "orient") {
            std::map<int, std::string>::const_iterator nm = name_of.find(centre);
            if (nm == name_of.end()) continue;
            std::string trimmed = nm->second;
            const size_t first = trimmed.find_first_not_of(" \t");
            if (first == std::string::npos) continue;
            trimmed = trimmed.substr(first);
            if (std::toupper((unsigned char) trimmed[0]) != 'S') continue;

            std::vector<int> c_nbrs, o_nbrs;
            for (size_t i = 0; i < nbrs.size(); ++i) {
                std::map<int, std::string>::const_iterator e = element_of.find(nbrs[i]);
                if (e == element_of.end()) continue;
                if (e->second == "C") c_nbrs.push_back(nbrs[i]);
                else if (e->second == "O") o_nbrs.push_back(nbrs[i]);
            }
            if (kind == "flat") {
                if (o_nbrs.size() < 3) continue;
                quad.push_back(o_nbrs[0]); quad.push_back(centre);
                quad.push_back(o_nbrs[1]); quad.push_back(o_nbrs[2]);
            } else {
                if (c_nbrs.size() != 1 || o_nbrs.size() < 2) continue;
                quad.push_back(c_nbrs[0]); quad.push_back(centre);
                quad.push_back(o_nbrs[0]); quad.push_back(o_nbrs[1]);
            }
        } else {
            continue;
        }

        if (kind == "ring" || kind == "pi") {
            std::vector<int> key;
            key.push_back(quad[0]); key.push_back(quad[2]); key.push_back(quad[3]);
            std::sort(key.begin(), key.end());
            if (!seen.insert(std::make_pair(centre, key)).second) continue;
        }
        out.push_back(quad);
    }
    return out;
}

IMPBFF_END_NAMESPACE
