/**
 * \file src/imp/DyeDynamics.cpp
 * \brief The dye roads, reached by file path rather than by IMP hierarchy.
 *
 * The implementation is the connection layer's own -- ProbeAttachment.h and
 * ProbeDynamics.h, IMP's hierarchies and integrators underneath. What this
 * file adds is the door: the model, the two hierarchies and the particle
 * indexes are made here and die here, so a caller passes paths and gets
 * arrays back and never holds an IMP object.
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/DyeDynamics.h>

#include <IMP/bff/HierarchyBridge.h>
#include <IMP/bff/ProbeAttachment.h>
#include <IMP/bff/ProbeDynamics.h>

#include <IMP/Model.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/pdb.h>

IMPBFF_BEGIN_NAMESPACE

namespace {

//! The structure, the dye placed on it, and the model that owns both.
/*! Both roads below need exactly this, and the model has to outlive the
    hierarchies, so it is one place rather than two. */
struct LabelledSite {
    IMP::Pointer<IMP::Model> model;
    IMP::atom::Hierarchy protein, dye;
    int n_stripped;
};

LabelledSite label_site(const std::string& protein_pdb, const std::string& dye_pdb,
                        const std::string& chain, int residue,
                        bool strip_site_sidechain) {
    LabelledSite s;
    s.model = new IMP::Model();
    s.protein = read_pdb_hierarchy(protein_pdb, s.model);
    s.dye = read_pdb_hierarchy(dye_pdb, s.model);
    const std::vector<ProbeAttachment> done = attach_probes(
            s.protein, std::vector<ProbeAttachment>(
                               1, ProbeAttachment(s.dye, chain, residue)),
            strip_site_sidechain);
    s.n_stripped = done.empty() ? 0 : done[0].get_n_stripped();
    return s;
}

}  // namespace

int attach_dye_to_pdb(const std::string& protein_pdb, const std::string& dye_pdb,
                      const std::string& chain, int residue,
                      const std::string& out_pdb, bool strip_site_sidechain) {
    LabelledSite s = label_site(protein_pdb, dye_pdb, chain, residue,
                                strip_site_sidechain);
    if (!out_pdb.empty()) {
        IMP::atom::write_pdb(s.protein, out_pdb);
    }
    return s.n_stripped;
}

LangevinTrajectory run_dye_langevin(
        const std::string& protein_pdb, const std::string& dye_pdb,
        const std::string& dye_mol2, const std::string& chain, int residue,
        int n_steps, int write_every, int minimize_steps,
        const std::string& integrator, double temperature, double timestep_fs,
        double friction_ps, double interaction_sphere, double repulsion_k,
        int seed, bool strip_site_sidechain) {
    LabelledSite s = label_site(protein_pdb, dye_pdb, chain, residue,
                                strip_site_sidechain);
    AttachedProbeDynamics dynamics(s.protein, s.dye, dye_mol2, chain, residue,
                                   integrator, temperature, timestep_fs,
                                   friction_ps, interaction_sphere, repulsion_k,
                                   seed);
    if (minimize_steps > 0) dynamics.minimize(minimize_steps);
    return dynamics.run(n_steps, write_every);
}

IMPBFF_END_NAMESPACE
