/**
 * \file ProbeForceField.cpp
 * \brief A coarse-grained dye/protein system: sites, bonds, and their terms.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/ProbeForceField.h>
#include <IMP/bff/MolecularGraph.h>
#include <IMP/bff/internal/json.h>

#include <IMP/bff/Base.h>

#include <algorithm>
#include <set>

IMPBFF_BEGIN_NAMESPACE

namespace {

double ffjson_num(const nlohmann::json& j, const char* k, double def) {
    if (j.contains(k) && j[k].is_number()) return j[k].get<double>();
    return def;
}

int ffjson_int(const nlohmann::json& j, const char* k, int def) {
    if (j.contains(k) && j[k].is_number_integer()) return j[k].get<int>();
    return def;
}

//! The first non-empty string of two keys.
std::string ffjson_pick(const nlohmann::json& j, const char* a, const char* b) {
    for (const char* k : {a, b}) {
        if (j.contains(k) && j[k].is_string() && !j[k].get<std::string>().empty())
            return j[k].get<std::string>();
    }
    return std::string();
}

//! A value as a string: strings verbatim, numbers by their JSON spelling,
//! null as empty.
std::string ffjson_val(const nlohmann::json& v) {
    if (v.is_string()) return v.get<std::string>();
    if (v.is_null()) return std::string();
    return v.dump();
}

//! A string field, whatever the JSON type: a numeric id (`"id": 0`) is the
//! string "0", not an error.
std::string ffjson_str(const nlohmann::json& j, const char* k) {
    return j.contains(k) ? ffjson_val(j[k]) : std::string();
}

}  // namespace

ProbeForceFieldSystem forcefield_system_from_json(const std::string& json) {
    nlohmann::json d = nlohmann::json::parse(json);
    if (!d.is_object()) {
        IMP_THROW("forcefield_system_from_json needs a JSON object",
                  ValueException);
    }
    ProbeForceFieldSystem out(ffjson_str(d, "name"));

    std::map<std::string, FFComponent> comps;
    if (d.contains("components") && d["components"].is_object()) {
        for (auto it = d["components"].begin(); it != d["components"].end();
             ++it) {
            FFComponent c;
            c.mol2_path = ffjson_pick(it.value(), "mol2_path", "mol2");
            c.pdb_path = ffjson_pick(it.value(), "pdb_path", "pdb");
            c.role = ffjson_pick(it.value(), "role", "");
            comps[it.key()] = c;
        }
    }
    out.set_components(comps);

    const nlohmann::json sites_in =
            d.contains("sites") && d["sites"].is_array() ? d["sites"]
                                                         : nlohmann::json::array();
    // A group member or bonded term may name a site by its number rather than
    // its id; resolve the ones the sites section declares.
    std::map<int, std::string> by_no;
    for (int i = 0; i < static_cast<int>(sites_in.size()); ++i) {
        const std::string id = ffjson_str(sites_in[i], "id");
        by_no[i] = id;
        if (sites_in[i].contains("site_no") &&
            sites_in[i]["site_no"].is_number_integer()) {
            by_no[sites_in[i]["site_no"].get<int>()] = id;
        }
    }
    auto token = [&](const nlohmann::json& x) -> std::string {
        if (x.is_string()) return x.get<std::string>();
        if (x.is_number_integer()) {
            auto f = by_no.find(x.get<int>());
            return f != by_no.end() ? f->second : x.dump();
        }
        return x.dump();
    };

    std::vector<FFSite> sites;
    for (int i = 0; i < static_cast<int>(sites_in.size()); ++i) {
        const nlohmann::json& row = sites_in[i];
        FFSite st;
        st.id = ffjson_str(row, "id");
        st.component = ffjson_str(row, "component");
        st.atom_name = ffjson_str(row, "atom_name");
        st.element = ffjson_str(row, "element");
        st.site_no = row.contains("site_no") &&
                             row["site_no"].is_number_integer()
                             ? row["site_no"].get<int>()
                             : i;
        st.site_serial = ffjson_int(row, "site_serial", 0);
        st.radius = ffjson_num(row, "radius", 1.7);
        st.mass = ffjson_num(row, "mass", 12.0);
        sites.push_back(st);
    }
    out.set_sites(sites);

    auto groups_from = [&](const char* key) {
        std::map<std::string, std::vector<std::string> > m;
        if (d.contains(key) && d[key].is_object()) {
            for (auto it = d[key].begin(); it != d[key].end(); ++it) {
                std::vector<std::string> v;
                if (it.value().is_array()) {
                    for (const auto& x : it.value()) v.push_back(token(x));
                }
                m[it.key()] = v;
            }
        }
        return m;
    };
    out.set_groups(groups_from("groups"));
    out.set_rb_groups(groups_from("rb_groups"));
    out.set_md_fixed_groups(groups_from("md_fixed_groups"));

    std::vector<std::string> fixed;
    if (d.contains("fixed_groups") && d["fixed_groups"].is_array()) {
        for (const auto& g : d["fixed_groups"]) {
            fixed.push_back(g.is_string() ? g.get<std::string>() : g.dump());
        }
    }
    out.set_fixed_groups(fixed);

    auto scalar_types = [&](const char* key) {
        std::map<std::string, double> m;
        if (d.contains(key) && d[key].is_object()) {
            for (auto it = d[key].begin(); it != d[key].end(); ++it) {
                m[it.key()] = ffjson_num(it.value(), "k", 0.0);
            }
        }
        return m;
    };
    out.set_bond_types(scalar_types("bond_types"));
    out.set_angle_types(scalar_types("angle_types"));

    auto torsion_types = [&](const char* key) {
        std::map<std::string, FFTorsionType> m;
        if (d.contains(key) && d[key].is_object()) {
            for (auto it = d[key].begin(); it != d[key].end(); ++it) {
                FFTorsionType t;
                t.periodicity = ffjson_int(it.value(), "periodicity", 1);
                t.phase = ffjson_num(it.value(), "phase_rad", 0.0);
                t.k = ffjson_num(it.value(), "k", 0.0);
                m[it.key()] = t;
            }
        }
        return m;
    };
    out.set_torsion_types(torsion_types("torsion_types"));
    out.set_improper_types(torsion_types("improper_types"));

    std::map<std::string, FFLJType> lj;
    if (d.contains("lj_types") && d["lj_types"].is_object()) {
        for (auto it = d["lj_types"].begin(); it != d["lj_types"].end(); ++it) {
            FFLJType t;
            t.element = ffjson_str(it.value(), "element");
            t.rmin_half = ffjson_num(it.value(), "rmin_half", 0.0);
            t.epsilon = ffjson_num(it.value(), "epsilon", 0.0);
            lj[it.key()] = t;
        }
    }
    out.set_lj_types(lj);

    if (d.contains("bonds") && d["bonds"].is_array()) {
        std::vector<FFBond> rows;
        for (const auto& row : d["bonds"]) {
            if (!row.is_array() || row.size() != 4) {
                IMP_THROW("a bond is [site_a, site_b, length, type_id]",
                          ValueException);
            }
            FFBond r;
            r.site_a = token(row[0]);
            r.site_b = token(row[1]);
            r.length = row[2].is_number() ? row[2].get<double>() : 0.0;
            r.type_id = ffjson_val(row[3]);
            rows.push_back(r);
        }
        out.set_bonds(rows);
    }
    if (d.contains("angles") && d["angles"].is_array()) {
        std::vector<FFAngle> rows;
        for (const auto& row : d["angles"]) {
            if (!row.is_array() || row.size() != 5) {
                IMP_THROW("an angle is [site_a, site_b, site_c, theta, type_id]",
                          ValueException);
            }
            FFAngle r;
            r.site_a = token(row[0]);
            r.site_b = token(row[1]);
            r.site_c = token(row[2]);
            r.theta = row[3].is_number() ? row[3].get<double>() : 0.0;
            r.type_id = ffjson_val(row[4]);
            rows.push_back(r);
        }
        out.set_angles(rows);
    }
    auto torsions_from = [&](const char* key) {
        std::vector<FFTorsion> rows;
        if (d.contains(key) && d[key].is_array()) {
            for (const auto& row : d[key]) {
                if (!row.is_array() || row.size() != 5) {
                    IMP_THROW("a torsion is [site_a, site_b, site_c, site_d,"
                              " type_id]",
                              ValueException);
                }
                FFTorsion r;
                r.site_a = token(row[0]);
                r.site_b = token(row[1]);
                r.site_c = token(row[2]);
                r.site_d = token(row[3]);
                r.type_id = ffjson_val(row[4]);
                rows.push_back(r);
            }
        }
        return rows;
    };
    out.set_dihedrals(torsions_from("dihedrals"));
    out.set_impropers(torsions_from("impropers"));

    if (d.contains("probes") && d["probes"].is_array()) {
        std::vector<FFProbe> probes;
        for (const auto& row : d["probes"]) {
            FFProbe pr;
            pr.id = ffjson_int(row, "id", 0);
            pr.name = ffjson_str(row, "name");
            pr.origin = ffjson_str(row, "origin");
            pr.link_type = ffjson_str(row, "link_type");
            if (pr.origin.empty()) pr.origin = "extrinsic";
            if (pr.link_type.empty()) pr.link_type = "covalent";
            probes.push_back(pr);
        }
        out.set_probes(probes);
    }

    FFNonbonded nb;
    const nlohmann::json empty = nlohmann::json::object();
    const nlohmann::json& nbs =
            d.contains("nonbonded") && d["nonbonded"].is_object()
                    ? d["nonbonded"] : empty;
    nb.enabled = !(nbs.contains("enabled") && nbs["enabled"].is_boolean() &&
                   !nbs["enabled"].get<bool>());
    nb.k = ffjson_num(nbs, "k", 5.0);
    nb.cutoff = ffjson_num(nbs, "cutoff_A", 6.0);
    out.set_nonbonded(nb);

    const nlohmann::json& sp =
            d.contains("sampling") && d["sampling"].is_object()
                    ? d["sampling"] : empty;
    FFSampling s;
    s.temperature_K = ffjson_num(sp, "temperature_K", 300.0);
    s.friction_ps = ffjson_num(sp, "friction_ps", 10.0);
    s.timestep_fs = ffjson_num(sp, "timestep_fs", 0.25);
    s.n_steps = ffjson_int(sp, "n_steps", 500000);
    s.write_every = ffjson_int(sp, "write_every", 1000);
    s.minimize_steps = ffjson_int(sp, "minimize_steps", 200);
    out.set_sampling(s);
    return out;
}

std::vector<int> ProbeForceFieldSystem::get_component_sites(
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

std::vector<std::string> ProbeForceFieldSystem::get_fixed_components() const {
    return with_role(components_, "fixed");
}

std::vector<std::string> ProbeForceFieldSystem::get_mobile_components() const {
    return with_role(components_, "mobile");
}


// --------------------------------------------------------------------------
// The molecular graph
// --------------------------------------------------------------------------

std::map<std::string, std::vector<std::string> >
ProbeForceFieldSystem::get_bonded_neighbors() const {
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

std::set<std::pair<std::string, std::string> >
ProbeForceFieldSystem::get_exclusions(bool include_impropers) const {
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
    return excl;
}

std::vector<std::vector<std::string> >
ProbeForceFieldSystem::find_rings(int max_len) const {
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

bool ProbeForceFieldSystem::is_within_bonds(const std::string& a,
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

std::string ProbeForceFieldSystem::get_inconsistency() const {
    if (components_.empty()) return "system requires components";
    if (sites_.empty()) return "system requires sites";

    std::set<std::string> ids;
    for (size_t i = 0; i < sites_.size(); ++i) {
        if (!ids.insert(sites_[i].id).second)
            return "duplicate site ids";
    }
    for (size_t i = 0; i < sites_.size(); ++i) {
        if (components_.find(sites_[i].component) == components_.end())
            return "site " + sites_[i].id + " unknown component " + sites_[i].component;
    }
    for (size_t i = 0; i < bonds_.size(); ++i) {
        if (!ids.count(bonds_[i].site_a) || !ids.count(bonds_[i].site_b))
            return "bond references unknown site";
    }
    for (size_t i = 0; i < angles_.size(); ++i) {
        if (!ids.count(angles_[i].site_a) || !ids.count(angles_[i].site_b) ||
            !ids.count(angles_[i].site_c))
            return "angle references unknown site";
    }
    const std::vector<FFTorsion>* blocks[2] = {&dihedrals_, &impropers_};
    const char* names[2] = {"dihedrals", "impropers"};
    for (int b = 0; b < 2; ++b) {
        for (size_t i = 0; i < blocks[b]->size(); ++i) {
            const FFTorsion& t = (*blocks[b])[i];
            if (!ids.count(t.site_a) || !ids.count(t.site_b) ||
                !ids.count(t.site_c) || !ids.count(t.site_d))
                return std::string(names[b]) + " references unknown site";
        }
    }
    return "";
}

IMPBFF_END_NAMESPACE
