/**
 * \file ProbeDynamics.cpp
 * \brief Langevin and Brownian dynamics of an attached explicit dye.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/ProbeDynamics.h>
// IMP::atom::Atom, used by value below and so needed complete
#include <IMP/atom/Atom.h>
// atom_name
#include <IMP/bff/IMPHierarchyBridge.h>
// create_probe_restraints
#include <IMP/bff/ProbePotentialRestraints.h>

#include <IMP/bff/ProbeAccessibleVolumeBuilder.h>
#include <IMP/bff/Mol2IO.h>
#include <IMP/bff/RotamerScoring.h>
#include <IMP/bff/ProbeTopology.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/HierarchyFrame.h>
#include <IMP/bff/internal/Text.h>

#include <IMP/atom/BrownianDynamics.h>
#include <IMP/atom/Chain.h>
#include <IMP/atom/Diffusion.h>
#include <IMP/atom/estimates.h>
#include <IMP/atom/LangevinThermostatOptimizerState.h>
#include <IMP/atom/MolecularDynamics.h>
#include <IMP/atom/Residue.h>
#include <IMP/atom/Mass.h>
#include <IMP/container/CloseBipartitePairContainer.h>
#include <IMP/container/ListSingletonContainer.h>
#include <IMP/container/PairsRestraint.h>
#include <IMP/core/ConjugateGradients.h>
#include <IMP/core/RestraintsScoringFunction.h>
#include <IMP/core/SphereDistancePairScore.h>
#include <IMP/core/XYZR.h>
#include <IMP/bff/IMPCompatibility.h>
#include <IMP/random.h>

#include <algorithm>
#include <cmath>
#include <limits>

IMPBFF_BEGIN_NAMESPACE

namespace {
//! What the attached molecule's component is called in its force field, and
//! therefore what its anchor group is called.
const char* const kComponentName = "label";
}  // namespace


void prepare_particles(IMP::Model* model, const IMP::ParticleIndexes& particles,
                       double temperature, const std::string& integrator,
                       const std::vector<double>& radii) {
    for (std::size_t i = 0; i < particles.size(); ++i) {
        IMP::Particle* p = model->get_particle(particles[i]);
        const bool has_radius = i < radii.size();
        if (!IMP::core::XYZR::get_is_setup(p)) {
            IMP::core::XYZR::setup_particle(p, has_radius ? radii[i] : 1.7);
        } else if (has_radius) {
            IMP::core::XYZR(p).set_radius(radii[i]);
        }
        if (!IMP::atom::Mass::get_is_setup(p)) {
            IMP::atom::Mass::setup_particle(p, 12.011);
        }
        if (integrator == "md") {
            if (!IMP::atom::LinearVelocity::get_is_setup(p)) {
                IMP::atom::LinearVelocity::setup_particle(
                        p, IMP::algebra::Vector3D(0, 0, 0));
            }
        } else {
            // Einstein's coefficient for this particle's own radius, so a big
            // atom diffuses more slowly than a small one without anybody
            // saying so. A^2/fs.
            const double d = IMP::atom::get_einstein_diffusion_coefficient(
                    IMP::core::XYZR(p).get_radius(), temperature);
            if (!IMP::atom::Diffusion::get_is_setup(p)) {
                IMP::atom::Diffusion::setup_particle(p, d);
            } else {
                IMP::atom::Diffusion(p).set_diffusion_coefficient(d);
            }
        }
    }
}

IMP::atom::Simulator* make_langevin_simulator(
        IMP::Model* model, const IMP::ParticleIndexes& mobile,
        IMP::ScoringFunction* scoring_function, const std::string& integrator,
        double temperature, double timestep_fs, double friction_ps, int seed) {
    if (integrator != "md" && integrator != "bd") {
        IMP_THROW("integrator must be 'md' or 'bd', not '" << integrator << "'",
                  ValueException);
    }
    if (seed >= 0) IMP::random_number_generator.seed(seed);
    for (std::size_t i = 0; i < mobile.size(); ++i) {
        IMP::core::XYZ(model, mobile[i]).set_coordinates_are_optimized(true);
    }
    prepare_particles(model, mobile, temperature, integrator);

    IMP::ParticlesTemp particles;
    for (std::size_t i = 0; i < mobile.size(); ++i) {
        particles.push_back(model->get_particle(mobile[i]));
    }

    if (integrator == "bd") {
        IMP_NEW(IMP::atom::BrownianDynamics, bd, (model));
        bd->set_particles(particles);
        bd->set_scoring_function(scoring_function);
        bd->set_maximum_time_step(timestep_fs);
        bd->set_temperature(temperature);
        return bd.release();
    }

    IMP_NEW(IMP::atom::MolecularDynamics, md, (model));
    md->set_particles(particles);
    md->set_scoring_function(scoring_function);
    md->set_maximum_time_step(timestep_fs);
    md->set_temperature(temperature);
    IMP_NEW(IMP::atom::LangevinThermostatOptimizerState, thermostat,
            (model, particles, temperature, friction_ps));
    thermostat->set_period(1);
    md->add_optimizer_state(thermostat);
    md->assign_velocities(temperature);
    return md.release();
}

// The trajectory's array views live with the record itself, in
// src/ProbeSimulation.cpp: ProbeSimulationTrajectory is the module's one
// trajectory type, shared by every simulation.

namespace {

}  // namespace

AttachedProbeDynamics::AttachedProbeDynamics(
        IMP::atom::Hierarchy protein_hier, IMP::atom::Hierarchy label_hier,
        const std::string& label_mol2, const std::string& chain, int residue,
        const std::string& integrator, double temperature, double timestep_fs,
        double friction_ps, double interaction_sphere, double repulsion_k,
        int seed)
    : model_(protein_hier.get_model()), protein_(protein_hier), dye_(label_hier),
      integrator_(integrator), temperature_(temperature),
      friction_ps_(friction_ps) {
    if (integrator != "md" && integrator != "bd") {
        IMP_THROW("integrator must be 'md' or 'bd', not '" << integrator << "'",
                  ValueException);
    }
    // Overdamped `bd` with the same stiff bonds is stable only below about
    // 2 fs, so the two integrators do not share a default.
    timestep_fs_ = timestep_fs > 0.0 ? timestep_fs
                                     : (integrator == "md" ? 2.0 : 0.5);

    // The label's atoms in MOL2 serial order -- the order its force-field
    // system's sites are in.
    IMP::atom::Hierarchies label_atoms =
            IMP::atom::get_by_type(label_hier, IMP::atom::ATOM_TYPE);
    std::sort(label_atoms.begin(), label_atoms.end(),
              [](IMP::atom::Hierarchy a, IMP::atom::Hierarchy b) {
                  return IMP::atom::Atom(a).get_input_index() <
                         IMP::atom::Atom(b).get_input_index();
              });
    for (std::size_t i = 0; i < label_atoms.size(); ++i) {
        label_particles_.push_back(label_atoms[i].get_particle_index());
    }

    system_ = probe_forcefield_system(label_mol2, kComponentName);
    const std::vector<FFSite>& sites = system_.get_sites();
    if (sites.size() != label_particles_.size()) {
        IMP_THROW("MOL2 has " << sites.size() << " atoms, the dye hierarchy "
                              << label_particles_.size(),
                  ValueException);
    }

    std::set<std::string> anchor;
    const std::map<std::string, std::vector<std::string> >& groups =
            system_.get_groups();
    // The component is named `label` above, so its anchor group is
    // `label_anchor`: the two are one name and must be written as one.
    std::map<std::string, std::vector<std::string> >::const_iterator g =
            groups.find(std::string(kComponentName) + "_anchor");
    if (g != groups.end()) anchor.insert(g->second.begin(), g->second.end());

    std::vector<double> radii;
    for (std::size_t i = 0; i < sites.size(); ++i) {
        site_ids_.push_back(sites[i].id);
        atom_names_.push_back(sites[i].atom_name);
        radii.push_back(vdw_radius(element_from_atom_name(sites[i].atom_name)));
        IMP::Particle* p = model_->get_particle(label_particles_[i]);
        if (!IMP::core::XYZR::get_is_setup(p)) {
            IMP::core::XYZR::setup_particle(p, radii.back());
        } else {
            IMP::core::XYZR(p).set_radius(radii.back());
        }
        if (!IMP::atom::Mass::get_is_setup(p)) {
            IMP::atom::Mass::setup_particle(p, sites[i].mass);
        } else {
            IMP::atom::Mass(p).set_mass(sites[i].mass);
        }
        if (anchor.count(sites[i].id) > 0) {
            fixed_.push_back(label_particles_[i]);
            IMP::core::XYZ(p).set_coordinates_are_optimized(false);
        } else {
            mobile_.push_back(label_particles_[i]);
        }
    }

    // The protein's heavy atoms near the site are the wall. The labelled
    // residue is not an obstacle to its own dye.
    bool found_ca = false;
    IMP::ParticleIndexes candidates;
    const IMP::atom::Hierarchies protein_atoms =
            IMP::atom::get_by_type(protein_hier, IMP::atom::ATOM_TYPE);
    for (std::size_t i = 0; i < protein_atoms.size(); ++i) {
        const IMP::atom::Atom atom(protein_atoms[i]);
        const std::string name = atom_name(atom);
        // Atom is a Hierarchy: the wrap-then-ask construction is
        // ambiguous under gcc (copy vs conversion), the call is not.
        const IMP::atom::Hierarchy parent = atom.get_parent();
        const bool is_residue =
                parent && IMP::atom::Residue::get_is_setup(parent);
        std::string chain_id;
        if (is_residue && parent.get_parent() &&
            IMP::atom::Chain::get_is_setup(parent.get_parent())) {
            chain_id = IMP::atom::Chain(parent.get_parent()).get_id();
        }
        const bool same_site =
                is_residue && IMP::atom::Residue(parent).get_index() == residue &&
                (chain.empty() || chain_id == chain);
        if (same_site && name == "CA") {
            site_ca_ = IMP::core::XYZ(atom).get_coordinates();
            found_ca = true;
        }
        if (same_site || (!name.empty() && (name[0] == 'H' || name[0] == 'h'))) {
            continue;
        }
        candidates.push_back(atom.get_particle_index());
    }
    if (!found_ca) {
        IMP_THROW("site " << chain << ":" << residue << " has no CA",
                  ValueException);
    }
    for (std::size_t i = 0; i < candidates.size(); ++i) {
        IMP::Particle* p = model_->get_particle(candidates[i]);
        if ((IMP::core::XYZ(p).get_coordinates() - site_ca_).get_magnitude() >
            interaction_sphere) {
            continue;
        }
        if (!IMP::core::XYZR::get_is_setup(p)) {
            IMP::core::XYZR::setup_particle(p, 1.7);
        }
        IMP::core::XYZ(p).set_coordinates_are_optimized(false);
        obstacles_.push_back(candidates[i]);
    }

    IMP::Restraints restraints = create_probe_restraints(model_, system_,
                                                      site_ids_,
                                                      label_particles_);
    if (!obstacles_.empty()) {
        IMP_NEW(IMP::container::ListSingletonContainer, label_container,
                (model_, mobile_));
        IMP_NEW(IMP::container::ListSingletonContainer, protein_container,
                (model_, obstacles_));
        IMP_NEW(IMP::container::CloseBipartitePairContainer, close,
                (label_container, protein_container, 3.0, 1.0));
        restraints.push_back(new IMP::container::PairsRestraint(
                new IMP::core::SoftSpherePairScore(repulsion_k), close,
                "label-protein"));
    }
    IMP_NEW(IMP::core::RestraintsScoringFunction, sf, (restraints));
    scoring_function_ = sf;
    simulator_ = make_langevin_simulator(model_, mobile_, scoring_function_,
                                         integrator_, temperature_,
                                         timestep_fs_, friction_ps_, seed);
}

namespace {
IMP::ParticlesTemp particles_of(IMP::Model* model,
                                const IMP::ParticleIndexes& indexes) {
    IMP::ParticlesTemp out;
    out.reserve(indexes.size());
    for (std::size_t i = 0; i < indexes.size(); ++i) {
        out.push_back(model->get_particle(indexes[i]));
    }
    return out;
}
}  // namespace

IMP::ParticlesTemp AttachedProbeDynamics::get_probe_particles() const {
    return particles_of(model_, label_particles_);
}
IMP::ParticlesTemp AttachedProbeDynamics::get_mobile() const {
    return particles_of(model_, mobile_);
}
IMP::ParticlesTemp AttachedProbeDynamics::get_fixed() const {
    return particles_of(model_, fixed_);
}
IMP::ParticlesTemp AttachedProbeDynamics::get_obstacles() const {
    return particles_of(model_, obstacles_);
}

void AttachedProbeDynamics::get_coordinates(double** out_view,
                                         int* n_out_view) const {
    std::vector<double> out;
    out.reserve(label_particles_.size() * 3);
    for (std::size_t i = 0; i < label_particles_.size(); ++i) {
        const IMP::algebra::Vector3D v =
                IMP::core::XYZ(model_, label_particles_[i]).get_coordinates();
        out.push_back(v[0]);
        out.push_back(v[1]);
        out.push_back(v[2]);
    }
    internal::copy_to_view(out, out_view, n_out_view);
}

double AttachedProbeDynamics::get_energy() const {
    return scoring_function_->evaluate(false);
}

double AttachedProbeDynamics::minimize(int n_steps) {
    IMP_NEW(IMP::core::ConjugateGradients, cg, (model_));
    cg->set_scoring_function(scoring_function_);
    return cg->optimize(std::max(0, n_steps));
}

ProbeSimulationTrajectory AttachedProbeDynamics::run(int n_steps, int write_every) {
    const int stride = std::max(1, write_every);
    const int n_frames = std::max(1, n_steps / stride);
    ProbeSimulationTrajectory out;
    out.n_frames = n_frames;
    out.n_atoms = static_cast<int>(label_particles_.size());
    out.atom_names = atom_names_;
    out.integrator = integrator_;
    out.temperature = temperature_;
    out.timestep_fs = timestep_fs_;
    out.coordinates.reserve(static_cast<std::size_t>(n_frames) *
                            label_particles_.size() * 3);
    const double nan = std::numeric_limits<double>::quiet_NaN();

    for (int frame = 0; frame < n_frames; ++frame) {
        simulator_->optimize(stride);
        for (std::size_t i = 0; i < label_particles_.size(); ++i) {
            const IMP::algebra::Vector3D v =
                    IMP::core::XYZ(model_, label_particles_[i]).get_coordinates();
            out.coordinates.push_back(v[0]);
            out.coordinates.push_back(v[1]);
            out.coordinates.push_back(v[2]);
        }
        out.times_fs.push_back((frame + 1) * stride * timestep_fs_);
        out.potential_energy.push_back(get_energy());
        // Brownian dynamics has no velocities, so it has no kinetic energy --
        // NaN and not zero, which would read as a system at absolute zero.
        out.kinetic_energy.push_back(
                integrator_ == "md"
                        ? dynamic_cast<IMP::atom::MolecularDynamics*>(
                                  simulator_.get())
                                  ->get_kinetic_energy()
                        : nan);
    }
    return out;
}

double AttachedProbeDynamics::kinetic_temperature(double kinetic_energy) const {
    const std::size_t n = mobile_.size();
    if (n == 0) return 0.0;
    return 2.0 * kinetic_energy / (3.0 * static_cast<double>(n) * kb_kcal());
}

IMPBFF_END_NAMESPACE
