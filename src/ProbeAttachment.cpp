/**
 * \file ProbeAttachment.cpp
 * \brief Putting a probe on a residue.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/ProbeAttachment.h>

#include <IMP/bff/HierarchyFrame.h>
#include <IMP/bff/internal/Text.h>
#include <IMP/bff/internal/json.h>

#include <IMP/algebra/Rotation3D.h>
#include <IMP/atom/Chain.h>
#include <IMP/atom/Residue.h>
#include <IMP/atom/Selection.h>
#include <IMP/atom/hierarchy_tools.h>
#include <IMP/core/XYZ.h>
#include <IMP/exception.h>
#include <IMP/rotamer/RotamerCalculator.h>
#include <IMP/rotamer/RotamerLibrary.h>

#include <algorithm>
#include <cmath>
#include <sstream>

IMPBFF_BEGIN_NAMESPACE

namespace {

//! `(chain, residue, NAME)` of an atom; an empty chain when it has none.
void atom_fields(IMP::atom::Atom atom, std::string& chain, int& residue,
                 std::string& name) {
    name = internal::trimmed(atom_name(atom));
    for (std::size_t i = 0; i < name.size(); ++i) {
        name[i] = static_cast<char>(
                std::toupper(static_cast<unsigned char>(name[i])));
    }
    chain.clear();
    residue = -1;
    const IMP::atom::Hierarchy parent = IMP::atom::Hierarchy(atom).get_parent();
    if (parent && IMP::atom::Residue::get_is_setup(parent)) {
        residue = IMP::atom::Residue(parent).get_index();
        const IMP::atom::Hierarchy grandparent = parent.get_parent();
        if (grandparent && IMP::atom::Chain::get_is_setup(grandparent)) {
            chain = IMP::atom::Chain(grandparent).get_id();
        }
    }
}

//! The first atom of \p name in `(chain, residue)`, or null.
IMP::Particle* find_atom(IMP::atom::Hierarchy hierarchy,
                         const std::string& chain_id, int residue,
                         const std::string& name) {
    const IMP::atom::Hierarchies atoms =
            IMP::atom::get_by_type(hierarchy, IMP::atom::ATOM_TYPE);
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        std::string chain, atom_name;
        int index = -1;
        atom_fields(IMP::atom::Atom(atoms[i]), chain, index, atom_name);
        if (index != residue) continue;
        if (!chain_id.empty() && chain != chain_id) continue;
        if (atom_name != name) continue;
        return atoms[i].get_particle();
    }
    return NULL;
}

}  // namespace

IMP::algebra::Transformation3D backbone_frame(
        const IMP::algebra::Vector3D& ca, const IMP::algebra::Vector3D& n,
        const IMP::algebra::Vector3D& c) {
    // x toward N, z out of the peptide plane, y completing the set. The same
    // three cross products `backbone_rotation` publishes as a flat matrix.
    const IMP::algebra::Vector3D x = (n - ca).get_unit_vector();
    const IMP::algebra::Vector3D toward_c = (c - ca).get_unit_vector();
    const IMP::algebra::Vector3D z =
            IMP::algebra::get_vector_product(x, toward_c).get_unit_vector();
    const IMP::algebra::Vector3D y = IMP::algebra::get_vector_product(z, x);
    return IMP::algebra::Transformation3D(
            IMP::algebra::get_rotation_from_x_y_axes(x, y), ca);
}

IMP::algebra::Transformation3D backbone_frame(IMP::atom::Hierarchy hierarchy,
                                              const std::string& chain_id,
                                              int residue) {
    const IMP::ParticlesTemp site = resolve_probe_site(hierarchy, chain_id,
                                                     residue);
    return backbone_frame(IMP::core::XYZ(site[0]).get_coordinates(),
                          IMP::core::XYZ(site[1]).get_coordinates(),
                          IMP::core::XYZ(site[2]).get_coordinates());
}

IMP::ParticlesTemp resolve_probe_site(IMP::atom::Hierarchy hierarchy,
                                    const std::string& chain_id, int residue) {
    static const char* kNames[3] = {"CA", "N", "C"};
    IMP::ParticlesTemp out;
    for (int i = 0; i < 3; ++i) {
        IMP::Particle* p = find_atom(hierarchy, chain_id, residue, kNames[i]);
        if (p == NULL) {
            IMP_THROW("atom " << kNames[i] << " not found for chain "
                              << (chain_id.empty() ? "*" : chain_id)
                              << " residue " << residue,
                      ValueException);
        }
        out.push_back(p);
    }
    return out;
}

IMP::ParticlesTemp select_atoms(IMP::atom::Hierarchy hierarchy,
                                const std::string& mask_text) {
    const StripSelection mask = parse_strip_mask(mask_text);
    IMP::ParticlesTemp out;
    const IMP::atom::Hierarchies atoms =
            IMP::atom::get_by_type(hierarchy, IMP::atom::ATOM_TYPE);
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        std::string chain, name;
        int residue = -1;
        atom_fields(IMP::atom::Atom(atoms[i]), chain, residue, name);
        if (mask.matches(chain, residue, name)) {
            out.push_back(atoms[i].get_particle());
        }
    }
    return out;
}

std::vector<int> strip_keep_mask(const std::vector<std::string>& chains,
                                 const std::vector<int>& resseqs,
                                 const std::vector<std::string>& names,
                                 const std::string& mask_text) {
    if (chains.size() != resseqs.size() || chains.size() != names.size()) {
        IMP_THROW("strip_keep_mask: " << chains.size() << " chains, "
                                      << resseqs.size() << " residues and "
                                      << names.size() << " names",
                  ValueException);
    }
    const StripSelection mask = parse_strip_mask(mask_text);
    std::vector<int> keep;
    keep.reserve(chains.size());
    for (std::size_t i = 0; i < chains.size(); ++i) {
        keep.push_back(mask.matches(chains[i], resseqs[i], names[i]) ? 0 : 1);
    }
    return keep;
}

int strip_hierarchy(IMP::atom::Hierarchy hierarchy,
                    const std::string& mask) {
    const IMP::ParticlesTemp victims = select_atoms(hierarchy, mask);
    for (std::size_t i = 0; i < victims.size(); ++i) {
        IMP::atom::Hierarchy atom(victims[i]);
        IMP::atom::Hierarchy parent = atom.get_parent();
        if (parent) parent.remove_child(atom);
        IMP::atom::destroy(atom);
    }
    return static_cast<int>(victims.size());
}

std::vector<std::string> site_keep_atom_names() {
    std::vector<std::string> out;
    out.push_back("N");
    out.push_back("CA");
    out.push_back("C");
    out.push_back("O");
    out.push_back("OXT");
    return out;
}

int strip_sidechain_at_site(IMP::atom::Hierarchy hierarchy,
                            const std::string& chain_id, int residue,
                            const std::vector<std::string>& keep_atom_names) {
    const std::vector<std::string> keep =
            keep_atom_names.empty() ? site_keep_atom_names() : keep_atom_names;
    return strip_hierarchy(hierarchy,
                           site_strip_mask(chain_id, residue, keep));
}

void place_probe(IMP::atom::Hierarchy probe,
                 const IMP::algebra::Transformation3D& frame) {
    IMP::atom::transform(probe, frame);
}

void place_probe_from_coords(IMP::atom::Hierarchy probe,
                             const IMP::algebra::Vector3D& ca,
                             const IMP::algebra::Vector3D& n,
                             const IMP::algebra::Vector3D& c) {
    place_probe(probe, backbone_frame(ca, n, c));
}

void align_hierarchies(IMP::atom::Hierarchy source, const std::string& chain_id,
                       int residue, const IMP::algebra::Vector3D& ca,
                       const IMP::algebra::Vector3D& n,
                       const IMP::algebra::Vector3D& c) {
    // Out of the source residue's frame and into the target's; what is left
    // is the move that takes one onto the other.
    const IMP::algebra::Transformation3D from =
            backbone_frame(source, chain_id, residue);
    const IMP::algebra::Transformation3D to = backbone_frame(ca, n, c);
    IMP::atom::transform(source, to * from.get_inverse());
}

std::vector<ProbeAttachment> attach_probes(
        IMP::atom::Hierarchy protein,
        const std::vector<ProbeAttachment>& attachments,
        bool strip_site_sidechain) {
    std::vector<ProbeAttachment> out = attachments;
    for (std::size_t i = 0; i < out.size(); ++i) {
        // The frame is read before anything is stripped. The side chain goes
        // and the backbone that defines the frame stays, so it would not
        // matter -- reading first is what keeps that an invariant rather than
        // an assumption about which atoms a mask happens to name.
        out[i].set_frame(backbone_frame(protein, out[i].get_chain(),
                                        out[i].get_residue()));
        if (strip_site_sidechain) {
            out[i].set_n_stripped(strip_sidechain_at_site(
                    protein, out[i].get_chain(), out[i].get_residue(),
                    std::vector<std::string>()));
        }
        place_probe(out[i].get_probe(), out[i].get_frame());
    }
    return out;
}

std::vector<std::string> fluorophore_types() {
    std::vector<std::string> out;
    out.push_back("donor");
    out.push_back("acceptor");
    out.push_back("unspecified");
    return out;
}

std::string probe_position_flrcif_items() {
    nlohmann::json out;
    out["asym_id"] = "_flr_poly_probe_position.asym_id";
    out["seq_id"] = "_flr_poly_probe_position.seq_id";
    out["comp_id"] = "_flr_poly_probe_position.comp_id";
    out["atom_id"] = "_flr_poly_probe_position.atom_id";
    out["auth_name"] = "_flr_poly_probe_position.auth_name";
    out["mutation_flag"] = "_flr_poly_probe_position.mutation_flag";
    out["modification_flag"] = "_flr_poly_probe_position.modification_flag";
    out["fluorophore_type"] = "_flr_sample_probe_details.fluorophore_type";
    out["description"] = "_flr_sample_probe_details.description";
    out["dye"] = "_flr_sample_probe_details.probe_id";
    return out.dump();
}

std::string quencher_flrcif_items() {
    nlohmann::json out;
    out["comp_id"] = "_flr_poly_probe_position.comp_id";
    out["asym_id"] = "_flr_poly_probe_position.asym_id";
    out["seq_id"] = "_flr_poly_probe_position.seq_id";
    out["atom_ids"] = "_flr_poly_probe_position.atom_id";
    // bff-native: no quenching category anywhere in the stack's dictionaries.
    out["rate_constant"] = nullptr;
    out["attenuation_length"] = nullptr;
    out["contact_distance"] = nullptr;
    out["dye"] = nullptr;
    return out.dump();
}

ProbePosition::ProbePosition(const std::string& asym_id, int seq_id, const std::string& atom_id,
             const std::string& probe, const std::string& fluorophore_type)
    : asym_id(asym_id), seq_id(seq_id), atom_id(atom_id), probe(probe),
      fluorophore_type(fluorophore_type), mutation_flag(false),
      modification_flag(false) {
    const std::vector<std::string> roles = fluorophore_types();
    if (std::find(roles.begin(), roles.end(), fluorophore_type) ==
        roles.end()) {
        IMP_THROW("fluorophore_type must be donor, acceptor or unspecified, "
                  "not '" << fluorophore_type << "'",
                  ValueException);
    }
}

void ProbePosition::set_dye(const Dye& d) {
    dye = d;
    probe = d.name;
}

std::string ProbePosition::get_key() const {
    std::ostringstream out;
    out << asym_id << ":" << seq_id << ":" << atom_id;
    return out.str();
}

std::string ProbePosition::get_source_info() const {
    nlohmann::json out;
    out["chain_identifier"] = asym_id;
    out["residue_seq_number"] = seq_id;
    out["atom_name"] = atom_id;
    return out.dump();
}

bool ProbePosition::operator==(const ProbePosition& other) const {
    return asym_id == other.asym_id && seq_id == other.seq_id &&
           atom_id == other.atom_id && probe == other.probe &&
           fluorophore_type == other.fluorophore_type;
}

ProbePosition probe_position_from_source_info(const std::string& source_info_json,
                             const Dye& dye) {
    nlohmann::json info = nlohmann::json::parse(source_info_json, NULL, false);
    if (info.is_discarded() || !info.is_object()) info = nlohmann::json::object();
    ProbePosition out;
    out.asym_id = info.value("chain_identifier", std::string());
    out.seq_id = info.value("residue_seq_number", 0);
    out.atom_id = info.value("atom_name", std::string("CB"));
    out.set_dye(dye);
    return out;
}

namespace {

//! The residue decorator at (chain, resnum), or an unset one.
IMP::atom::Residue residue_at(IMP::atom::Hierarchy hierarchy,
                              const std::string& chain_id, int resnum) {
    IMP::atom::Selection sel(hierarchy);
    sel.set_chain_id(chain_id);
    sel.set_residue_index(resnum);
    sel.set_resolution(IMP::atom::ALL_RESOLUTIONS);
    const IMP::ParticlesTemp selected = sel.get_selected_particles();
    for (unsigned int i = 0; i < selected.size(); ++i) {
        if (IMP::atom::Residue::get_is_setup(selected[i])) {
            return IMP::atom::Residue(selected[i]);
        }
    }
    // A selection resolves to atoms when the hierarchy has no residue-level
    // representation to select; the residue is then the atom's parent.
    for (unsigned int i = 0; i < selected.size(); ++i) {
        IMP::atom::Hierarchy parent =
                IMP::atom::Hierarchy(selected[i]).get_parent();
        if (parent && IMP::atom::Residue::get_is_setup(parent)) {
            return IMP::atom::Residue(parent);
        }
    }
    return IMP::atom::Residue();
}

}  // namespace

IMP::algebra::Vector3D get_anchor_cb_position(IMP::atom::Hierarchy hierarchy,
                                              std::string chain_id, int resnum,
                                              std::string rl_path,
                                              double prob_threshold) {
    IMP::atom::Residue residue = residue_at(hierarchy, chain_id, resnum);
    IMP_USAGE_CHECK(residue, "no residue " << chain_id << resnum
                                           << " in the structure");

    IMP_NEW(IMP::rotamer::RotamerLibrary, library, ());
    if (!rl_path.empty()) library->read_library_file(rl_path);
    IMP_NEW(IMP::rotamer::RotamerCalculator, calculator, (library));
    const IMP::rotamer::ResidueRotamer rotamer =
            calculator->get_rotamer(residue, prob_threshold);

    // Index 0 is the deposited geometry and index 1 the best rotamer, so a
    // size of one means the library had nothing to say about this residue --
    // glycine, or a residue it does not carry.
    if (rotamer.get_size() > 1 &&
        rotamer.get_atom_exists(IMP::atom::AT_CB)) {
        return rotamer.get_coordinates(1, IMP::atom::AT_CB);
    }

    const IMP::ParticlesTemp cb = select_atoms(
            hierarchy, "chain " + chain_id + " and resid " +
                               std::to_string(resnum) + " and name CB");
    IMP_USAGE_CHECK(!cb.empty(), "residue " << chain_id << resnum
                                            << " has no C-beta, in the "
                                            << "library or in the structure");
    return IMP::core::XYZ(cb[0]).get_coordinates();
}

IMPBFF_END_NAMESPACE
