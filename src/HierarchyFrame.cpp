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

// PRD-137 step 5c residue: load_protein_frames reads the PDB through
// IMP::atom::read_multimodel_pdb and converts each hierarchy with the
// connection layer's protein_frame_from_hierarchy (HierarchyBridge.h). The
// core reader (AVBuilder.h's read_pdb_records) is the replacement.
ProteinFrame protein_frame_from_hierarchy(IMP::atom::Hierarchy hierarchy);

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


// --------------------------------------------------------------------------
// ProteinFrame
// --------------------------------------------------------------------------

void ProteinFrame::get_coords(double** out_view, int* n_out_view) const {
    internal::copy_to_view(coords, out_view, n_out_view);
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
