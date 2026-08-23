/**
 * \file TopologyBuild.cpp
 * \brief The one force-field system builder, from components to the typed
 *        value.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/TopologyBuild.h>

#include <IMP/bff/CifIO.h>
#include <IMP/bff/Mol2IO.h>
#include <IMP/bff/MolecularGraph.h>
#include <IMP/bff/Scoring.h>
#include <IMP/bff/internal/json.h>

#include <IMP/exception.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

namespace {

std::string tb_upper(std::string s) {
    for (auto& c : s) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    return s;
}

//! `0 -> "A" ... 25 -> "Z", 26 -> "AA"`: the spreadsheet column spelling.
std::string alpha_suffix(int i) {
    std::string chars;
    long n = i;
    while (true) {
        chars.push_back(static_cast<char>('A' + (n % 26)));
        n = n / 26 - 1;
        if (n < 0) break;
    }
    std::reverse(chars.begin(), chars.end());
    return chars;
}

//! One component's read, with the graph and the deduplicated site names.
struct Component {
    std::string name, mol2, role, template_path;
    std::map<int, Mol2Atom> atoms;                  // serial -> atom
    std::vector<std::pair<int, int>> bonds;         // sorted serial pairs
    std::map<int, std::string> site_names;          // serial -> deduped name
    bool has_template = false;
    ComponentTemplate tmpl;

    std::string sid(int serial) const {
        return name + "/" + site_names.at(serial);
    }
};

//! The graph of a component's bonds, numbered by sorted serial.
MolecularGraph graph_of(const Component& c) {
    std::vector<std::pair<int, int>> edges;
    for (const auto& b : c.bonds) edges.push_back({b.first, b.second});
    return MolecularGraph(edges);
}

double tb_distance(const Mol2Atom& a, const Mol2Atom& b) {
    const double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double angle_value(const Mol2Atom& a, const Mol2Atom& b, const Mol2Atom& c) {
    const double bax = a.x - b.x, bay = a.y - b.y, baz = a.z - b.z;
    const double bcx = c.x - b.x, bcy = c.y - b.y, bcz = c.z - b.z;
    const double ba2 = bax * bax + bay * bay + baz * baz;
    const double bc2 = bcx * bcx + bcy * bcy + bcz * bcz;
    // the tetrahedral angle when an arm has zero length: a harmonic minimum
    // at a collapsed angle is not a default anyone wants
    if (ba2 == 0.0 || bc2 == 0.0) return 1.910633;
    const double dot = bax * bcx + bay * bcy + baz * bcz;
    const double cosang = std::max(-1.0, std::min(
            1.0, dot / std::sqrt(ba2 * bc2)));
    return std::acos(cosang);
}

//! The serials a template's improper kind centres on.
std::set<int> center_serials(const ComponentTemplate& tmpl,
                             const std::string& kind,
                             const std::map<int, Mol2Atom>& atoms) {
    std::set<std::string> center_names;
    for (const auto& imp : tmpl.impropers) {
        if (imp.type == kind) center_names.insert(imp.center_atom);
    }
    std::set<int> out;
    for (const auto& [serial, atom] : atoms) {
        if (center_names.count(atom.atom_name)) out.insert(serial);
    }
    return out;
}

//! The site ids a template feature selects, by (atom name, occurrence).
std::vector<std::string> resolve_feature_ids(
        const ComponentTemplate& tmpl, const std::string& feature_id,
        const std::string& comp_name, const std::map<int, Mol2Atom>& atoms,
        const std::map<int, std::string>& serial_to_site) {
    std::map<std::string, std::vector<int>> by_name;
    for (const auto& [serial, atom] : atoms) {
        by_name[atom.atom_name].push_back(serial);
    }
    std::set<std::string> out;
    const auto fit = tmpl.features.find(feature_id);
    if (fit != tmpl.features.end()) {
        for (const auto& entry : fit->second.atoms) {
            const auto found = by_name.find(entry.name);
            if (found == by_name.end()) continue;
            const int idx = entry.occurrence - 1;
            if (idx >= 0 && idx < static_cast<int>(found->second.size())) {
                const int serial = found->second[idx];
                out.insert(comp_name + "/" + serial_to_site.at(serial));
            }
        }
    }
    return std::vector<std::string>(out.begin(), out.end());
}

}  // namespace

std::map<int, std::string> serial_to_site_atom_names(
        const std::map<int, std::string>& serial_to_name) {
    std::map<std::string, std::vector<int>> by_name;
    for (const auto& [serial, name] : serial_to_name) {
        by_name[name].push_back(serial);
    }
    std::map<int, std::string> out;
    for (const auto& [name, serials] : by_name) {
        if (serials.size() == 1) {
            out[serials[0]] = name;
            continue;
        }
        for (int i = 0; i < static_cast<int>(serials.size()); i++) {
            out[serials[i]] = name + alpha_suffix(i);
        }
    }
    return out;
}

DyeForceFieldSystem build_forcefield_system(
        const std::string& components_json, double bond_k, double angle_k,
        double pi_dihedral_k, double linker_dihedral_k, double ring_improper_k,
        double pi_improper_k, double flat_improper_k, double orient_improper_k,
        int n_steps, int write_every, double default_radius, double default_mass,
        double nonbonded_k, double nonbonded_cutoff, int minimize_steps,
        const std::string& relative_to) {
    const nlohmann::json specs = nlohmann::json::parse(components_json);
    if (!specs.is_array()) {
        IMP_THROW("components_json must be a JSON array of component specs",
                  ValueException);
    }

    // -- read every component once ------------------------------------------
    std::vector<Component> components;
    int n_fixed = 0;
    for (const auto& spec : specs) {
        Component c;
        c.name = spec.value("name", std::string());
        c.mol2 = spec.value("mol2", std::string());
        c.role = spec.value("role", std::string());
        if (spec.contains("template") && spec["template"].is_string()) {
            c.template_path = spec["template"].get<std::string>();
        }
        if (c.name.empty() || c.mol2.empty() || c.role.empty()) {
            IMP_THROW("Component spec must include name, mol2, role",
                      ValueException);
        }
        if (c.role == "fixed") n_fixed++;
        if (n_fixed > 1) {
            IMP_THROW("At most one fixed component is supported",
                      ValueException);
        }
        const Mol2Component read = read_mol2_component(c.mol2, c.name);
        for (const auto& atom : read.atoms) c.atoms[atom.serial] = atom;
        c.bonds = read.bonds;
        std::map<int, std::string> serial_to_name;
        for (const auto& [serial, atom] : c.atoms) {
            serial_to_name[serial] = atom.atom_name;
        }
        c.site_names = serial_to_site_atom_names(serial_to_name);
        if (!c.template_path.empty()) {
            c.tmpl = read_component_template_cif(c.template_path);
            c.has_template = true;
        }
        std::printf("%s (%s): %zu atoms, %zu bonds\n",
                    c.role == "fixed" ? "Fixed" : "Mobile", c.name.c_str(),
                    c.atoms.size(), c.bonds.size());
        components.push_back(c);
    }
    // fixed first, mobile after -- the site order the Python emitted
    std::stable_sort(components.begin(), components.end(),
                     [](const Component& a, const Component& b) {
                         return (a.role == "fixed") && (b.role != "fixed");
                     });

    // -- the JSON the Python dict held, then the one conversion path --------
    nlohmann::json system;
    std::string system_name;
    for (const auto& c : components) {
        if (!system_name.empty()) system_name += "_";
        system_name += c.name;
    }
    system["name"] = system_name;

    nlohmann::json comp_out = nlohmann::json::object();
    for (const auto& c : components) {
        nlohmann::json entry;
        entry["role"] = c.role;
        if (!relative_to.empty()) {
            // os.path.relpath semantics: both absolute, relative to out_dir
            std::string abs = c.mol2, base = relative_to;
            const std::size_t common = [&] {
                std::size_t i = 0;
                while (i < abs.size() && i < base.size() &&
                       abs[i] == base[i]) {
                    if (abs[i] == '/') ++i;
                    else {
                        // extend to the end of the matching path segment
                        const std::size_t na = abs.find('/', i);
                        const std::size_t nb = base.find('/', i);
                        if (na == nb && na != std::string::npos) i = na + 1;
                        else break;
                    }
                }
                return i;
            }();
            std::string rel_a = abs.substr(common);
            std::string rel_b = base.substr(common);
            int ups = 0;
            for (char ch : rel_b) ups += (ch == '/');
            std::string rel;
            for (int u = 0; u < ups; u++) rel += "../";
            rel += rel_a;
            entry["mol2"] = rel.empty() ? "." : rel;
        } else {
            entry["mol2"] = c.mol2;
        }
        comp_out[c.name] = entry;
    }
    system["components"] = comp_out;

    // probes: the mobile components, one each
    nlohmann::json probes = nlohmann::json::array();
    for (const auto& c : components) {
        if (c.role != "mobile") continue;
        probes.push_back({{"id", static_cast<int>(probes.size()) + 1},
                          {"name", c.name},
                          {"origin", "extrinsic"},
                          {"link_type", "covalent"}});
    }
    system["probes"] = probes;

    // sites
    nlohmann::json sites = nlohmann::json::array();
    for (const auto& c : components) {
        for (const auto& [serial, _] : c.atoms) {
            sites.push_back({{"id", c.sid(serial)},
                             {"component", c.name},
                             {"atom_name", c.site_names.at(serial)},
                             {"site_serial", serial},
                             {"radius", default_radius},
                             {"mass", default_mass}});
        }
    }
    system["sites"] = sites;

    // bonds
    nlohmann::json bonds = nlohmann::json::array();
    for (const auto& c : components) {
        for (const auto& [a, b] : c.bonds) {
            bonds.push_back({c.sid(a), c.sid(b),
                             tb_distance(c.atoms.at(a), c.atoms.at(b)), "B1"});
        }
    }
    system["bonds"] = bonds;

    // angles
    nlohmann::json angles = nlohmann::json::array();
    for (const auto& c : components) {
        const MolecularGraph g = graph_of(c);
        for (const auto& a : g.get_angles()) {
            angles.push_back({c.sid(a[0]), c.sid(a[1]), c.sid(a[2]),
                              angle_value(c.atoms.at(a[0]), c.atoms.at(a[1]),
                                          c.atoms.at(a[2])),
                              "A1"});
        }
    }
    system["angles"] = angles;

    // dihedrals (mobile only): pi when both central atoms are C/N, else link
    nlohmann::json dihedrals = nlohmann::json::array();
    for (const auto& c : components) {
        if (c.role != "mobile") continue;
        const MolecularGraph g = graph_of(c);
        for (const auto& d : g.get_dihedrals()) {
            const std::string eb = c.atoms.at(d[1]).element;
            const std::string ec = c.atoms.at(d[2]).element;
            const bool pi = (eb == "C" || eb == "N") && (ec == "C" || ec == "N");
            dihedrals.push_back({c.sid(d[0]), c.sid(d[1]), c.sid(d[2]),
                                 c.sid(d[3]), pi ? "T_PI" : "T_LINK"});
        }
    }
    system["dihedrals"] = dihedrals;

    // impropers: each template kind expanded against the bond graph.
    // MolecularGraph keys nodes by the raw pair values -- the serials -- so
    // centers, nodes, elements and names are all in serial spelling.
    nlohmann::json impropers = nlohmann::json::array();
    const char* kinds[] = {"ring", "pi", "flat", "orient"};
    const char* itypes[] = {"I_RING", "I_PI", "I_FLAT", "I_ORIENT"};
    for (const auto& c : components) {
        if (!c.has_template) continue;
        const MolecularGraph g = graph_of(c);
        std::vector<int> nodes;
        std::vector<std::string> elements, names;
        for (const auto& [serial, atom] : c.atoms) {
            nodes.push_back(serial);
            elements.push_back(atom.element);
            names.push_back(atom.atom_name);
        }
        for (int k = 0; k < 4; k++) {
            const std::set<int> centers = center_serials(c.tmpl, kinds[k],
                                                         c.atoms);
            if (centers.empty()) continue;
            const auto quads = g.expand_impropers(
                    kinds[k],
                    std::vector<int>(centers.begin(), centers.end()),
                    nodes, elements, names, 8);
            for (const auto& q : quads) {
                impropers.push_back({c.sid(q[0]), c.sid(q[1]), c.sid(q[2]),
                                     c.sid(q[3]), itypes[k]});
            }
        }
    }
    system["impropers"] = impropers;

    // groups: comp_all, plus one per template dof feature
    nlohmann::json groups = nlohmann::json::object();
    nlohmann::json rb_groups = nlohmann::json::object();
    nlohmann::json md_fixed_groups = nlohmann::json::object();
    for (const auto& c : components) {
        nlohmann::json all = nlohmann::json::array();
        for (const auto& [serial, _] : c.atoms) all.push_back(c.sid(serial));
        groups[c.name + "_all"] = all;

        if (!c.has_template) continue;
        for (const auto& [fid, spec] : c.tmpl.features) {
            if (spec.feature_type != "dof") continue;
            const std::vector<std::string> ids = resolve_feature_ids(
                    c.tmpl, fid, c.name, c.atoms, c.site_names);
            if (ids.empty()) continue;
            if (spec.rb) {
                rb_groups[c.name + "_" + fid + "_rb"] = ids;
            } else if (spec.md_fixed) {
                md_fixed_groups[c.name + "_" + fid + "_md_fixed"] = ids;
            } else {
                groups[c.name + "_" + fid] = ids;
            }
        }
    }
    system["groups"] = groups;
    system["rb_groups"] = rb_groups;
    system["md_fixed_groups"] = md_fixed_groups;

    nlohmann::json fixed_groups = nlohmann::json::array();
    for (const auto& c : components) {
        if (c.role == "fixed") fixed_groups.push_back(c.name + "_all");
    }
    system["fixed_groups"] = fixed_groups;

    // type tables
    system["bond_types"] = {{"B1", {{"k", bond_k}}}};
    system["angle_types"] = {{"A1", {{"k", angle_k}}}};
    system["torsion_types"] = {
            {"T_PI", {{"periodicity", 2}, {"phase_rad", M_PI},
                      {"k", pi_dihedral_k}}},
            {"T_LINK", {{"periodicity", 3}, {"phase_rad", 0.0},
                        {"k", linker_dihedral_k}}}};
    system["improper_types"] = {
            {"I_RING", {{"periodicity", 2}, {"phase_rad", 0.0},
                        {"k", ring_improper_k}}},
            {"I_PI", {{"periodicity", 2}, {"phase_rad", 0.0},
                      {"k", pi_improper_k}}},
            {"I_FLAT", {{"periodicity", 2}, {"phase_rad", 0.0},
                        {"k", flat_improper_k}}},
            {"I_ORIENT", {{"periodicity", 2}, {"phase_rad", 0.0},
                          {"k", orient_improper_k}}}};

    // LJ types: one per element the site names imply, from the single C++
    // table (build_lj_type_table -> charmm36), so there is no second source
    std::set<std::string> all_elements;
    for (const auto& c : components) {
        for (const auto& [serial, atom] : c.atoms) {
            all_elements.insert(element_from_atom_name(
                    c.site_names.at(serial)));
        }
    }
    nlohmann::json lj_types = nlohmann::json::object();
    for (const auto& [key, t] :
             build_lj_type_table(std::vector<std::string>(all_elements.begin(),
                                                          all_elements.end()))) {
        lj_types[key] = {{"element", t.element},
                         {"rmin_half", t.rmin_half},
                         {"epsilon", t.epsilon}};
    }
    system["lj_types"] = lj_types;

    system["nonbonded"] = {{"enabled", true},
                           {"k", nonbonded_k},
                           {"cutoff_A", nonbonded_cutoff}};
    system["sampling"] = {{"n_steps", n_steps},
                          {"write_every", write_every},
                          {"minimize_steps", minimize_steps}};

    return forcefield_system_from_json(system.dump());
}

IMPBFF_END_NAMESPACE
