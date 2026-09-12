/**
 * \file src/imp/MolecularProbeSimulation.cpp
 * \brief A dye on a structure, behind arrays instead of hierarchies.
 *
 * The implementation is the connection layer's own -- ProbeAttachment.h and
 * ProbeDynamics.h, IMP's hierarchies and integrators underneath. What this
 * file adds is the shape: the model, the two hierarchies and the particle
 * indexes are made here and die here, so a caller passes arrays, drives an
 * ordinary simulation object, and never holds an IMP object.
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/MolecularProbeSimulation.h>

#include <IMP/bff/IMPHierarchyBridge.h>
#include <IMP/bff/ProbeAttachment.h>
#include <IMP/bff/ProbeDynamics.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/internal/json.h>

#include <IMP/Model.h>
#include <IMP/atom/Hierarchy.h>

#include <memory>

IMPBFF_BEGIN_NAMESPACE

//! Everything of IMP the simulation needs, out of the header's sight.
struct MolecularProbeSimulation::Impl {
    IMP::Pointer<IMP::Model> model;
    IMP::atom::Hierarchy protein, dye;
    std::unique_ptr<AttachedProbeDynamics> dynamics;
    int n_stripped;
    // what is in force, kept so that get_parameters() reports the run rather
    // than the defaults
    std::string integrator;
    double temperature, timestep_fs, friction_ps, interaction_sphere, repulsion_k;
    int seed;
    bool strip_site_sidechain;

    Impl()
        : n_stripped(0), integrator("md"), temperature(300.0), timestep_fs(-1.0),
          friction_ps(10.0), interaction_sphere(25.0), repulsion_k(10.0),
          seed(-1), strip_site_sidechain(true) {}

    //! Read a partial JSON object over the fields; an unknown key is an error.
    void read(const std::string& text) {
        if (text.empty()) return;
        nlohmann::json j;
        try {
            j = nlohmann::json::parse(text);
        } catch (const std::exception& e) {
            IMP_THROW("the parameters are not JSON: " << e.what(), ValueException);
        }
        if (!j.is_object()) {
            IMP_THROW("the parameters must be a JSON object, not "
                      << j.type_name(), ValueException);
        }
        for (nlohmann::json::const_iterator it = j.begin(); it != j.end(); ++it) {
            const std::string& k = it.key();
            if (k == "integrator") integrator = it.value().get<std::string>();
            else if (k == "temperature") temperature = it.value().get<double>();
            else if (k == "timestep_fs") timestep_fs = it.value().get<double>();
            else if (k == "friction_ps") friction_ps = it.value().get<double>();
            else if (k == "interaction_sphere") interaction_sphere = it.value().get<double>();
            else if (k == "repulsion_k") repulsion_k = it.value().get<double>();
            else if (k == "seed") seed = it.value().get<int>();
            else if (k == "strip_site_sidechain") strip_site_sidechain = it.value().get<bool>();
            else {
                IMP_THROW("a dye simulation has no parameter '" << k
                          << "'; it has integrator, temperature, timestep_fs, "
                             "friction_ps, interaction_sphere, repulsion_k, "
                             "seed and strip_site_sidechain", ValueException);
            }
        }
        if (integrator != "md" && integrator != "bd") {
            IMP_THROW("integrator is 'md' or 'bd', not '" << integrator << "'",
                      ValueException);
        }
        if (temperature <= 0.0) {
            IMP_THROW("temperature is above zero, not " << temperature,
                      ValueException);
        }
    }

    std::string write() const {
        nlohmann::json j;
        j["integrator"] = integrator;
        j["temperature"] = temperature;
        j["timestep_fs"] = timestep_fs;
        j["friction_ps"] = friction_ps;
        j["interaction_sphere"] = interaction_sphere;
        j["repulsion_k"] = repulsion_k;
        j["seed"] = seed;
        j["strip_site_sidechain"] = strip_site_sidechain;
        return j.dump();
    }
};

MolecularProbeSimulation::MolecularProbeSimulation(const ProteinFrame& protein, const ProteinFrame& dye,
                             const std::string& dye_mol2,
                             const std::string& chain, int residue,
                             const std::string& parameters)
    : impl_(new Impl()) {
    impl_->read(parameters);
    impl_->model = new IMP::Model();
    impl_->protein = hierarchy_from_protein_frame(protein, impl_->model);
    impl_->dye = hierarchy_from_protein_frame(dye, impl_->model);

    const std::vector<ProbeAttachment> done = attach_probes(
            impl_->protein,
            std::vector<ProbeAttachment>(
                    1, ProbeAttachment(impl_->dye, chain, residue)),
            impl_->strip_site_sidechain);
    impl_->n_stripped = done.empty() ? 0 : done[0].get_n_stripped();

    impl_->dynamics.reset(new AttachedProbeDynamics(
            impl_->protein, impl_->dye, dye_mol2, chain, residue,
            impl_->integrator, impl_->temperature, impl_->timestep_fs,
            impl_->friction_ps, impl_->interaction_sphere, impl_->repulsion_k,
            impl_->seed));
}

std::string MolecularProbeSimulation::get_parameters() const { return impl_->write(); }

void MolecularProbeSimulation::set_parameters(const std::string& json) {
    // The integrator, the thermostat and the wall of protein spheres are
    // built once, when the dye is placed. Changing them afterwards would
    // report one thing and integrate another, so it is refused; build
    // another simulation instead.
    IMP_THROW("a dye simulation's parameters are fixed once it is built"
              " (asked to set " << json << "); construct another one",
              ValueException);
}

int MolecularProbeSimulation::get_n_atoms() const {
    return static_cast<int>(impl_->dynamics->get_atom_names().size());
}

void MolecularProbeSimulation::get_positions(double** out_view, int* n_out_view) const {
    impl_->dynamics->get_coordinates(out_view, n_out_view);
}

void MolecularProbeSimulation::set_positions(const std::vector<double>& xyz) {
    const std::size_t expected = static_cast<std::size_t>(get_n_atoms()) * 3;
    if (xyz.size() != expected) {
        IMP_THROW("the dye has " << get_n_atoms() << " atoms, so " << expected
                  << " coordinates, not " << xyz.size(), ValueException);
    }
    // the dye's own hierarchy is what the integrator reads coordinates from
    apply_coordinates(impl_->dye, xyz);
}

double MolecularProbeSimulation::minimize(int n_steps) {
    return impl_->dynamics->minimize(n_steps);
}

void MolecularProbeSimulation::step(int n_steps) {
    // The frames are what run() is for; stepping keeps only the state.
    impl_->dynamics->run(n_steps, n_steps > 0 ? n_steps : 1);
}

ProbeSimulationTrajectory MolecularProbeSimulation::run(int n_steps, int write_every) {
    return impl_->dynamics->run(n_steps, write_every);
}

double MolecularProbeSimulation::get_potential_energy() const {
    return impl_->dynamics->get_energy();
}

std::vector<std::string> MolecularProbeSimulation::get_atom_names() const {
    return impl_->dynamics->get_atom_names();
}

int MolecularProbeSimulation::get_n_stripped() const { return impl_->n_stripped; }

ProteinFrame MolecularProbeSimulation::get_labelled_frame() const {
    // The dye keeps its own root -- attaching moves it, it does not adopt it
    // -- so the labelled structure is the two of them, the dye last.
    ProteinFrame out = protein_frame_from_hierarchy(impl_->protein);
    const ProteinFrame label = protein_frame_from_hierarchy(impl_->dye);
    out.coords.insert(out.coords.end(), label.coords.begin(), label.coords.end());
    out.atom_names.insert(out.atom_names.end(), label.atom_names.begin(),
                          label.atom_names.end());
    out.atom_types.insert(out.atom_types.end(), label.atom_types.begin(),
                          label.atom_types.end());
    out.resnames.insert(out.resnames.end(), label.resnames.begin(),
                        label.resnames.end());
    out.chain_ids.insert(out.chain_ids.end(), label.chain_ids.begin(),
                         label.chain_ids.end());
    out.residue_indices.insert(out.residue_indices.end(),
                               label.residue_indices.begin(),
                               label.residue_indices.end());
    return out;
}

double MolecularProbeSimulation::kinetic_temperature(double kinetic_energy) const {
    return impl_->dynamics->kinetic_temperature(kinetic_energy);
}

IMPBFF_END_NAMESPACE
