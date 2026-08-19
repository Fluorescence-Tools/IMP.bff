/**
 * \file DyeForceField.cpp
 * \brief A coarse-grained dye/protein system: sites, bonds, and their terms.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/DyeForceField.h>
#include <IMP/bff/MolecularGraph.h>

#include <algorithm>
#include <set>

IMPBFF_BEGIN_NAMESPACE

std::vector<int> DyeForceFieldSystem::get_component_sites(
        const std::string& component) const {
    std::vector<int> out;
    for (std::size_t i = 0; i < sites_.size(); ++i) {
        if (sites_[i].component == component) {
            out.push_back(static_cast<int>(i));
        }
    }
    return out;
}

namespace {
std::vector<std::string> with_role(
        const std::map<std::string, FFComponent>& components,
        const std::string& role) {
    std::vector<std::string> out;
    for (std::map<std::string, FFComponent>::const_iterator it =
                 components.begin(); it != components.end(); ++it) {
        if (it->second.role == role) out.push_back(it->first);
    }
    return out;
}
}  // namespace

std::vector<std::string> DyeForceFieldSystem::get_fixed_components() const {
    return with_role(components_, "fixed");
}

std::vector<std::string> DyeForceFieldSystem::get_mobile_components() const {
    return with_role(components_, "mobile");
}


// --------------------------------------------------------------------------
// The molecular graph
// --------------------------------------------------------------------------

std::map<std::string, std::vector<std::string> >
DyeForceFieldSystem::get_bonded_neighbors() const {
    std::map<std::string, std::set<std::string> > adj;
    for (size_t i = 0; i < bonds_.size(); ++i) {
        adj[bonds_[i].site_a].insert(bonds_[i].site_b);
        adj[bonds_[i].site_b].insert(bonds_[i].site_a);
    }
    std::map<std::string, std::vector<std::string> > out;
    for (std::map<std::string, std::set<std::string> >::const_iterator it = adj.begin();
         it != adj.end(); ++it) {
        out[it->first] = std::vector<std::string>(it->second.begin(), it->second.end());
    }
    return out;
}

namespace {
std::pair<std::string, std::string> ordered_pair(const std::string& a,
                                                 const std::string& b) {
    return a <= b ? std::make_pair(a, b) : std::make_pair(b, a);
}
}

std::vector<std::pair<std::string, std::string> >
DyeForceFieldSystem::get_exclusions(bool include_impropers) const {
    std::set<std::pair<std::string, std::string> > excl;
    for (size_t i = 0; i < bonds_.size(); ++i)
        excl.insert(ordered_pair(bonds_[i].site_a, bonds_[i].site_b));
    for (size_t i = 0; i < angles_.size(); ++i)
        excl.insert(ordered_pair(angles_[i].site_a, angles_[i].site_c));
    for (size_t i = 0; i < dihedrals_.size(); ++i)
        excl.insert(ordered_pair(dihedrals_[i].site_a, dihedrals_[i].site_d));
    if (include_impropers) {
        for (size_t i = 0; i < impropers_.size(); ++i) {
            const FFTorsion& t = impropers_[i];
            excl.insert(ordered_pair(t.site_a, t.site_c));
            excl.insert(ordered_pair(t.site_a, t.site_d));
            excl.insert(ordered_pair(t.site_b, t.site_d));
        }
    }
    return std::vector<std::pair<std::string, std::string> >(excl.begin(), excl.end());
}

std::vector<std::vector<std::string> >
DyeForceFieldSystem::find_rings(int max_len) const {
    // through `MolecularGraph`, so the ring walk exists once: site ids are
    // numbered in sorted order and mapped back, which is what makes the
    // canonical form come out ascending by site id rather than by number.
    std::vector<std::string> ids;
    std::map<std::string, int> index;
    const std::map<std::string, std::vector<std::string> > adj = get_bonded_neighbors();
    for (std::map<std::string, std::vector<std::string> >::const_iterator it = adj.begin();
         it != adj.end(); ++it) {
        index[it->first] = (int)ids.size();
        ids.push_back(it->first);
    }
    std::vector<std::pair<int, int> > edges;
    for (size_t i = 0; i < bonds_.size(); ++i) {
        std::map<std::string, int>::const_iterator a = index.find(bonds_[i].site_a);
        std::map<std::string, int>::const_iterator b = index.find(bonds_[i].site_b);
        if (a != index.end() && b != index.end())
            edges.push_back(std::make_pair(a->second, b->second));
    }
    const MolecularGraph g(edges);
    const std::vector<std::vector<int> > rings = g.get_rings(max_len);
    std::vector<std::vector<std::string> > out;
    out.reserve(rings.size());
    for (size_t i = 0; i < rings.size(); ++i) {
        std::vector<std::string> r;
        r.reserve(rings[i].size());
        for (size_t j = 0; j < rings[i].size(); ++j) r.push_back(ids[rings[i][j]]);
        out.push_back(r);
    }
    return out;
}

bool DyeForceFieldSystem::is_within_bonds(const std::string& a,
                                          const std::string& b,
                                          int max_depth) const {
    if (a == b) return true;
    const std::map<std::string, std::vector<std::string> > adj = get_bonded_neighbors();
    std::set<std::string> seen;
    seen.insert(a);
    std::vector<std::string> frontier(1, a);
    for (int depth = 0; depth < max_depth && !frontier.empty(); ++depth) {
        std::vector<std::string> next;
        for (size_t i = 0; i < frontier.size(); ++i) {
            std::map<std::string, std::vector<std::string> >::const_iterator it =
                adj.find(frontier[i]);
            if (it == adj.end()) continue;
            for (size_t j = 0; j < it->second.size(); ++j) {
                const std::string& n = it->second[j];
                if (n == b) return true;
                if (seen.insert(n).second) next.push_back(n);
            }
        }
        frontier.swap(next);
    }
    return false;
}

IMPBFF_END_NAMESPACE
