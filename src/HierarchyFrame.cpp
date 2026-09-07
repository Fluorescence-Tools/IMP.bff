/**
 * \file HierarchyFrame.cpp
 * \brief Reading one frame out of an IMP hierarchy, in one call.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/HierarchyFrame.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/internal/Text.h>

#include <IMP/Model.h>
#include <IMP/atom/pdb.h>

#include <IMP/atom/Atom.h>
#include <IMP/atom/Chain.h>
#include <IMP/atom/Residue.h>
#include <IMP/core/XYZ.h>

IMPBFF_BEGIN_NAMESPACE

std::string atom_name(IMP::atom::Atom atom) {
    std::string name = atom.get_atom_type().get_string();
    const std::size_t het = name.find("HET:");
    if (het != std::string::npos) name.erase(het, 4);
    return internal::trimmed(name);
}

namespace {
//! The XYZ leaves, in hierarchy order.
IMP::ParticlesTemp xyz_leaves(IMP::atom::Hierarchy hierarchy) {
    IMP::ParticlesTemp out;
    for (IMP::atom::Hierarchy leaf : IMP::atom::get_leaves(hierarchy)) {
        IMP::Particle* p = leaf.get_particle();
        if (IMP::core::XYZ::get_is_setup(p)) out.push_back(p);
    }
    return out;
}
}  // namespace

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

// --------------------------------------------------------------------------
// ProteinFrame
// --------------------------------------------------------------------------

void ProteinFrame::get_coords(double** out_view, int* n_out_view) const {
    internal::copy_to_view(coords, out_view, n_out_view);
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

std::vector<ProteinFrame> load_protein_frames(const std::string& path,
                                              int max_frames) {
    if (!internal::ends_with(path, ".pdb") && !internal::ends_with(path, ".ent")) {
        IMP_THROW("not a PDB file: " << path
                                     << " (RMF is read through IMP.rmf and "
                                        "protein_frame_from_hierarchy)",
                  IMP::ValueException);
    }
    IMP_NEW(IMP::Model, model, ());
    IMP::atom::Hierarchies hierarchies = IMP::atom::read_multimodel_pdb(
            path, model, new IMP::atom::NonWaterPDBSelector());
    std::vector<ProteinFrame> frames;
    for (std::size_t i = 0; i < hierarchies.size(); ++i) {
        if (max_frames >= 0 && static_cast<int>(frames.size()) >= max_frames) {
            break;
        }
        frames.push_back(protein_frame_from_hierarchy(hierarchies[i]));
    }
    return frames;
}

IMPBFF_END_NAMESPACE
