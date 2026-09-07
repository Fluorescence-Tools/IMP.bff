/**
 * \file HierarchyBridge.cpp
 * \brief The Hierarchy and Particle overloads of the core's functions.
 *
 * Sections in the order of IMP/bff/HierarchyBridge.h; each is marked with
 * the file it came from.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/HierarchyBridge.h>
#include <IMP/bff/AVBuilder.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/internal/Text.h>
#include <IMP/bff/internal/json.h>

#include <IMP/atom/Chain.h>
#include <IMP/atom/Residue.h>
#include <IMP/atom/Mass.h>
#include <IMP/atom/pdb.h>
#include <IMP/core/XYZ.h>
#include <IMP/core/XYZR.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

// -------- from HierarchyFrame.cpp --------
namespace {
std::vector<double> hierarchy_atom_coordinates_impl(IMP::atom::Hierarchy hierarchy) {
    const IMP::ParticlesTemp particles = xyz_leaves(hierarchy);
    std::vector<double> out(particles.size() * 4, 0.0);
    for (std::size_t i = 0; i < particles.size(); ++i) {
        const IMP::algebra::Vector3D v =
                IMP::core::XYZ(particles[i]).get_coordinates();
        out[4 * i + 0] = v[0];
        out[4 * i + 1] = v[1];
        out[4 * i + 2] = v[2];
        double residue_index = -1.0;
        if (IMP::atom::Atom::get_is_setup(particles[i])) {
            IMP::Particle* parent =
                    IMP::atom::Hierarchy(particles[i]).get_parent().get_particle();
            if (parent && IMP::atom::Residue::get_is_setup(parent)) {
                residue_index = IMP::atom::Residue(parent).get_index();
            }
        }
        out[4 * i + 3] = residue_index;
    }
    return out;
}
}  // namespace

void hierarchy_atom_coordinates(IMP::atom::Hierarchy hierarchy, double** out_view, int* n_out_view) {
    internal::copy_to_view(hierarchy_atom_coordinates_impl(hierarchy),
                           out_view, n_out_view);
}

std::string atom_name(IMP::atom::Atom atom) {
    std::string name = atom.get_atom_type().get_string();
    const std::size_t het = name.find("HET:");
    if (het != std::string::npos) name.erase(het, 4);
    return internal::trimmed(name);
}

std::vector<std::string> hierarchy_atom_metadata(IMP::atom::Hierarchy hierarchy) {
    const IMP::ParticlesTemp particles = xyz_leaves(hierarchy);
    std::vector<std::string> out;
    out.reserve(particles.size() * 4);
    for (std::size_t i = 0; i < particles.size(); ++i) {
        IMP::Particle* p = particles[i];
        std::string name, type, resname, chain;
        if (IMP::atom::Atom::get_is_setup(p)) {
            const std::string full = p->get_name();
            // IMP names an atom "Atom CB of residue 1" -- five fields, not two.
            // The atom name is the SECOND field; the last is the residue
            // number, which is easy to return by accident for every atom.
            const std::size_t first = full.find_first_not_of(" \t");
            const std::size_t gap =
                    first == std::string::npos ? std::string::npos
                                               : full.find_first_of(" \t", first);
            if (gap == std::string::npos) {
                name = first == std::string::npos ? full : full.substr(first);
            } else {
                const std::size_t start = full.find_first_not_of(" \t", gap);
                const std::size_t end = start == std::string::npos
                                                ? std::string::npos
                                                : full.find_first_of(" \t", start);
                name = start == std::string::npos
                               ? full.substr(first, gap - first)
                               : full.substr(start, end == std::string::npos
                                                            ? std::string::npos
                                                            : end - start);
            }
            type = IMP::atom::Atom(p).get_atom_type().get_string();
            IMP::Particle* res =
                    IMP::atom::Hierarchy(p).get_parent().get_particle();
            if (res && IMP::atom::Residue::get_is_setup(res)) {
                resname = IMP::atom::Residue(res).get_residue_type().get_string();
                IMP::Particle* ch =
                        IMP::atom::Hierarchy(res).get_parent().get_particle();
                if (ch && IMP::atom::Chain::get_is_setup(ch)) {
                    chain = IMP::atom::Chain(ch).get_id();
                }
            }
        } else {
            name = p->get_name();
            // A bare particle has no atom type; its first character is what the
            // Python used, and an empty name falls back to carbon.
            type = name.empty() ? std::string("C") : name.substr(0, 1);
        }
        out.push_back(name);
        out.push_back(type);
        out.push_back(resname);
        out.push_back(chain);
    }
    return out;
}

ProteinFrame protein_frame_from_hierarchy(IMP::atom::Hierarchy hierarchy) {
    ProteinFrame frame;
    const std::vector<double> packed = hierarchy_atom_coordinates_impl(hierarchy);
    const int n_atoms = static_cast<int>(packed.size() / 4);
    frame.coords.reserve(static_cast<std::size_t>(n_atoms) * 3);
    frame.residue_indices.reserve(n_atoms);
    for (int i = 0; i < n_atoms; ++i) {
        frame.coords.push_back(packed[i * 4]);
        frame.coords.push_back(packed[i * 4 + 1]);
        frame.coords.push_back(packed[i * 4 + 2]);
        frame.residue_indices.push_back(static_cast<int>(packed[i * 4 + 3]));
    }
    const std::vector<std::string> meta = hierarchy_atom_metadata(hierarchy);
    const std::size_t n = meta.size() / 4;
    frame.atom_names.reserve(n);
    frame.atom_types.reserve(n);
    frame.resnames.reserve(n);
    frame.chain_ids.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        frame.atom_names.push_back(meta[i * 4]);
        frame.atom_types.push_back(meta[i * 4 + 1]);
        frame.resnames.push_back(meta[i * 4 + 2]);
        frame.chain_ids.push_back(meta[i * 4 + 3]);
    }
    return frame;
}


// -------- from StructureIO.cpp --------
void LoadedStructure::get_coords(double** output, int* n_output1,
                                 int* n_output2) const {
    double* out = internal::new_double_view(coords_.size(), output, n_output1);
    if (out == NULL) return;
    if (!coords_.empty()) {
        std::memcpy(out, coords_.data(), coords_.size() * sizeof(double));
    }
    *n_output1 = static_cast<int>(coords_.size() / 3);
    *n_output2 = 3;
}

IMP::atom::Hierarchy read_pdb_hierarchy(const std::string& path,
                                        IMP::Model* m) {
    return IMP::atom::read_pdb(path, m, new IMP::atom::NonWaterPDBSelector());
}

void structure_coordinates(IMP::atom::Hierarchy hierarchy, double** out_view,
                           int* n_out_view) {
    const IMP::atom::Hierarchies leaves = IMP::atom::get_leaves(hierarchy);
    double* out = internal::new_double_view(leaves.size() * 3, out_view,
                                           n_out_view);
    if (out == nullptr) return;
    for (unsigned int i = 0; i < leaves.size(); ++i) {
        const IMP::algebra::Vector3D xyz =
                IMP::core::XYZ(leaves[i]).get_coordinates();
        for (int k = 0; k < 3; ++k) out[3 * i + k] = xyz[k];
    }
}

LoadedStructure load_structure_with_particles(const std::string& path,
                                              IMP::Model* model) {
    const IMP::atom::Hierarchy hier = read_pdb_hierarchy(path, model);
    const IMP::ParticlesTemp leaves = IMP::atom::get_leaves(hier);
    std::vector<double> coords;
    coords.reserve(leaves.size() * 3);
    for (unsigned int i = 0; i < leaves.size(); ++i) {
        const IMP::algebra::Vector3D v =
                IMP::core::XYZ(leaves[i]).get_coordinates();
        coords.push_back(v[0]);
        coords.push_back(v[1]);
        coords.push_back(v[2]);
    }
    return LoadedStructure(hier, leaves, coords);
}

FlexFitSelection read_angle_file(IMP::atom::Hierarchy hier,
                                 const std::string& flexfit_json) {
    nlohmann::json block;
    try {
        block = nlohmann::json::parse(flexfit_json);
    } catch (const std::exception& e) {
        IMP_THROW("read_angle_file: the FlexFit block does not parse: "
                  << e.what(), IMP::ValueException);
    }
    if (!block.contains("Flexible residues") || !block.contains("Bonds")) {
        IMP_THROW("read_angle_file: a FlexFit block needs 'Flexible residues' "
                  "and 'Bonds'", IMP::ValueException);
    }
    IMP::Model* model = hier.get_model();

    std::vector<AtomReference> residues;
    for (const nlohmann::json& fr : block["Flexible residues"]) {
        residues.push_back(AtomReference(
                fr.value("chain_identifier", std::string()),
                fr.value("residue_seq_number", 0)));
    }
    std::vector<AtomReference> ends;
    for (const nlohmann::json& bond : block["Bonds"]) {
        for (int e = 0; e < 2 && e < static_cast<int>(bond.size()); ++e) {
            const nlohmann::json& end = bond[e];
            ends.push_back(AtomReference(
                    end.value("chain_identifier", std::string()),
                    end.value("residue_seq_number", 0),
                    end.value("atom_name", std::string())));
        }
    }

    IMP::ParticlesTemp residue_particles;
    const IMP::ParticleIndexes flexible =
            select_flexible_residues(hier, residues);
    for (unsigned int i = 0; i < flexible.size(); ++i) {
        residue_particles.push_back(model->get_particle(flexible[i]));
    }
    IMP::atom::Bonds bond_decorators;
    const IMP::ParticleIndexes bonds = create_named_bonds(hier, ends);
    for (unsigned int i = 0; i < bonds.size(); ++i) {
        bond_decorators.push_back(IMP::atom::Bond(model, bonds[i]));
    }
    return FlexFitSelection(residue_particles, bond_decorators);
}

IMP::ParticleIndexes select_flexible_residues(
        IMP::atom::Hierarchy hierarchy,
        const std::vector<AtomReference>& residues) {
    IMP::ParticleIndexes out;
    for (std::size_t i = 0; i < residues.size(); ++i) {
        IMP::atom::Selection sel(hierarchy);
        sel.set_chain_ids(
                IMP::Strings(1, residues[i].chain_identifier));
        sel.set_residue_indexes(IMP::Ints(1, residues[i].residue_seq_number));
        const IMP::ParticleIndexes found = sel.get_selected_particle_indexes(false);
        if (found.empty()) {
            IMP_THROW("no residue " << residues[i].chain_identifier << ":"
                                    << residues[i].residue_seq_number
                                    << " in the hierarchy",
                      ValueException);
        }
        out.push_back(found[0]);
    }
    return out;
}

IMP::ParticleIndexes create_named_bonds(
        IMP::atom::Hierarchy hierarchy,
        const std::vector<AtomReference>& atom_pairs) {
    if (atom_pairs.size() % 2 != 0) {
        IMP_THROW("a bond takes two atoms; got " << atom_pairs.size()
                                                 << " references",
                  ValueException);
    }
    IMP::ParticleIndexes out;
    for (std::size_t i = 0; i < atom_pairs.size(); i += 2) {
        IMP::ParticleIndex ends[2];
        for (int k = 0; k < 2; ++k) {
            const AtomReference& ref = atom_pairs[i + k];
            IMP::atom::Selection sel(hierarchy);
            sel.set_chain_ids(IMP::Strings(1, ref.chain_identifier));
            sel.set_residue_indexes(IMP::Ints(1, ref.residue_seq_number));
            sel.set_atom_type(IMP::atom::AtomType(ref.atom_name));
            const IMP::ParticleIndexes found = sel.get_selected_particle_indexes();
            if (found.empty()) {
                IMP_THROW("no atom " << ref.chain_identifier << ":"
                                     << ref.residue_seq_number << ":"
                                     << ref.atom_name << " in the hierarchy",
                          ValueException);
            }
            ends[k] = found[0];
        }
        IMP::atom::Bonded b1(hierarchy.get_model(), ends[0]);
        IMP::atom::Bonded b2(hierarchy.get_model(), ends[1]);
        out.push_back(IMP::atom::create_bond(b1, b2, IMP::atom::Bond::SINGLE)
                              .get_particle_index());
    }
    return out;
}


// -------- from SelectionExpression.cpp --------

std::vector<SelectionAtom> selection_atoms(IMP::atom::Hierarchy hierarchy) {
    std::vector<SelectionAtom> out;
    const IMP::atom::Hierarchies leaves = IMP::atom::get_leaves(hierarchy);
    out.reserve(leaves.size());
    int index = 0;
    for (IMP::atom::Hierarchy leaf : leaves) {
        SelectionAtom a;
        a.index = ++index;
        if (IMP::core::XYZ::get_is_setup(leaf)) {
            const IMP::algebra::Vector3D xyz =
                    IMP::core::XYZ(leaf).get_coordinates();
            a.x = xyz[0]; a.y = xyz[1]; a.z = xyz[2];
        }
        if (IMP::atom::Atom::get_is_setup(leaf)) {
            IMP::atom::Atom atom(leaf);
            a.name = internal::trimmed(atom.get_atom_type().get_string());
            a.id = atom.get_input_index();
            a.elem = IMP::atom::get_element_table().get_name(atom.get_element());
            a.hetatm = a.name.compare(0, 4, "HET:") == 0;
            if (a.hetatm) a.name = a.name.substr(4);
        } else {
            a.name = internal::trimmed(leaf->get_name());
        }
        IMP::atom::Hierarchy parent = leaf.get_parent();
        if (parent && IMP::atom::Residue::get_is_setup(parent)) {
            IMP::atom::Residue residue(parent);
            a.resn = internal::trimmed(residue.get_residue_type().get_string());
            a.resi = residue.get_index();
            const char icode = residue.get_insertion_code();
            if (icode != ' ' && icode != '\0') {
                a.resi_text = std::to_string(a.resi) + std::string(1, icode);
            }
            IMP::atom::Hierarchy chain_h = parent.get_parent();
            if (chain_h && IMP::atom::Chain::get_is_setup(chain_h)) {
                a.chain = IMP::atom::Chain(chain_h).get_id();
            }
        }
        out.push_back(a);
    }
    return out;
}

std::vector<int> select_atom_indices(IMP::atom::Hierarchy hierarchy,
                                     const std::string& expression) {
    const std::vector<SelectionAtom> atoms = selection_atoms(hierarchy);
    const SelectionExpression selection(expression);
    const std::vector<int> mask = selection.evaluate(atoms);
    std::vector<int> out;
    for (std::size_t i = 0; i < mask.size(); ++i) {
        if (mask[i]) out.push_back(static_cast<int>(i));
    }
    return out;
}


// -------- from StripMask.cpp --------
std::vector<double> strip_obstacles(IMP::atom::Hierarchy hierarchy,
                                    const std::string& mask,
                                    const std::string& keep) {
    const std::vector<SelectionAtom> atoms = selection_atoms(hierarchy);
    const IMP::atom::Hierarchies leaves = IMP::atom::get_leaves(hierarchy);

    std::vector<int> drop(atoms.size(), 0);
    const std::string text = internal::trimmed(mask);
    if (!text.empty()) drop = SelectionExpression(text).evaluate(atoms);

    std::vector<double> out;
    out.reserve(atoms.size() * 4);
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        double radius = 0.0;
        bool stripped = drop[i] != 0;
        if (stripped && !keep.empty()) {
            const std::string id = atoms[i].chain + "/" +
                    std::to_string(atoms[i].resi) + "/" + atoms[i].name;
            if (id == keep) stripped = false;
        }
        if (!stripped && IMP::core::XYZR::get_is_setup(leaves[i])) {
            radius = IMP::core::XYZR(leaves[i]).get_radius();
        }
        out.push_back(atoms[i].x);
        out.push_back(atoms[i].y);
        out.push_back(atoms[i].z);
        out.push_back(radius);
    }
    return out;
}

StripReport strip_report(IMP::atom::Hierarchy hierarchy,
                         const std::string& mask) {
    StripReport report;
    report.mask = mask;
    const std::vector<SelectionAtom> atoms = selection_atoms(hierarchy);
    report.n_atoms = static_cast<int>(atoms.size());
    if (internal::trimmed(mask).empty()) return report;

    const SelectionExpression sel = parse_strip_mask(mask);
    const std::vector<int> selected = sel.evaluate(atoms);
    std::set<std::string> seen;
    for (std::size_t i = 0; i < selected.size(); ++i) {
        if (!selected[i]) continue;
        const SelectionAtom& a = atoms[i];
        const std::string residue =
                a.chain + "/" + (a.resi_text.empty() ? std::to_string(a.resi)
                                                     : a.resi_text);
        report.atoms.push_back(residue + "/" + a.name);
        if (seen.insert(residue).second) report.residues.push_back(residue);
    }
    report.n_selected = static_cast<int>(report.atoms.size());
    return report;
}


// -------- from ProbeSampling.cpp --------
void apply_coordinates(const IMP::atom::Hierarchy hierarchy,
                               const std::vector<double>& coords) {
    IMP::atom::Hierarchies leaves = IMP::atom::get_leaves(hierarchy);
    const std::size_t n_atoms = leaves.size();
    if (coords.size() != n_atoms * 3) {
        IMP_THROW("Atom count mismatch: coords=" << coords.size() / 3
                          << " hierarchy=" << n_atoms,
                  ValueException);
    }
    for (std::size_t i = 0; i < n_atoms; ++i) {
        IMP::core::XYZ atom(leaves[i]);
        atom.set_coordinates(IMP::algebra::Vector3D(
                coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]));
    }
}


// -------- from VdwRadii.cpp --------
double olga_vdw_particle_radius(IMP::Particle *p) {
    if (p == nullptr) return get_olga_vdw_fallback_radius();
    if (!IMP::atom::Atom::get_is_setup(p)) {
        return get_olga_vdw_fallback_radius();
    }
    return olga_vdw_radius(
            IMP::atom::Atom(p).get_atom_type().get_string());
}

std::vector<std::string> olga_vdw_unknown_atom_names(
        const IMP::ParticlesTemp &ps) {
    const std::map<std::string, double> &t = get_olga_vdw_radii();
    std::set<std::string> missing;
    for (std::size_t i = 0; i < ps.size(); ++i) {
        if (!IMP::atom::Atom::get_is_setup(ps[i])) continue;
        const std::string n =
                IMP::atom::Atom(ps[i]).get_atom_type().get_string();
        if (t.find(n) == t.end()) missing.insert(n);
    }
    return std::vector<std::string>(missing.begin(), missing.end());
}

std::vector<double> olga_vdw_radii(const IMP::ParticlesTemp &ps) {
    std::vector<double> out;
    out.reserve(ps.size());
    for (std::size_t i = 0; i < ps.size(); ++i) {
        out.push_back(olga_vdw_particle_radius(ps[i]));
    }
    return out;
}

// -------- from PathMap.cpp --------
void set_path_map_particles(PathMap* map, const IMP::ParticlesTemp& ps) {
    IMP_USAGE_CHECK(map, "set_path_map_particles: no map");
    IMP::ParticlesTemp held(ps);
    map->set_sphere_source([held](GridSpheres& spheres) {
        spheres.clear();
        spheres.reserve(held.size());
        for (std::size_t i = 0; i < held.size(); ++i) {
            IMP::core::XYZR xyzr(held[i]);
            spheres.push_back(GridSphere(xyzr.get_coordinates(), xyzr.get_radius()));
        }
    });
}

IMPBFF_END_NAMESPACE
