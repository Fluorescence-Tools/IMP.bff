/**
 * \file HierarchyFrame.cpp
 * \brief Reading one frame out of an IMP hierarchy, in one call.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/HierarchyFrame.h>

#include <IMP/atom/Atom.h>
#include <IMP/atom/Chain.h>
#include <IMP/atom/Residue.h>
#include <IMP/core/XYZ.h>

IMPBFF_BEGIN_NAMESPACE

namespace {
//! The XYZ leaves, in hierarchy order -- the same set the Python selected.
IMP::ParticlesTemp xyz_leaves(IMP::atom::Hierarchy hierarchy) {
    IMP::ParticlesTemp out;
    for (IMP::atom::Hierarchy leaf : IMP::atom::get_leaves(hierarchy)) {
        IMP::Particle* p = leaf.get_particle();
        if (IMP::core::XYZ::get_is_setup(p)) out.push_back(p);
    }
    return out;
}
}  // namespace

std::vector<double> hierarchy_atom_coordinates(IMP::atom::Hierarchy hierarchy) {
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
            // The atom name is the SECOND, which is what the Python took; the
            // last is the residue number, which a first cut of this returned
            // for every atom in the structure.
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

IMPBFF_END_NAMESPACE
