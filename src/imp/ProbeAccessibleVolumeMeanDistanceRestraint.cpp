/**
 * \file ProbeAccessibleVolumeMeanDistanceRestraint.cpp
 * \brief A FRET restraint on the distance between two mean dye positions.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/ProbeAccessibleVolumeMeanDistanceRestraint.h>
#include <IMP/bff/States.h>
#include <cmath>
#include <limits>
#include <fstream>
#include <IMP/atom/Chain.h>
#include <IMP/bff/internal/json.h>
#include <IMP/atom/Residue.h>
#include <IMP/atom/Atom.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/algebra/Vector3D.h>
#include <IMP/core/Harmonic.h>
#include <IMP/core/DistanceRestraint.h>
#include <IMP/bff/internal/Text.h>

#include <IMP/atom/Mass.h>
#include <IMP/core/XYZ.h>
#include <IMP/core/XYZR.h>
#include <IMP/core/rigid_bodies.h>

#include <algorithm>
#include <map>

IMPBFF_BEGIN_NAMESPACE

ProbeAccessibleVolumeMeanDistanceRestraint::ProbeAccessibleVolumeMeanDistanceRestraint(
        IMP::Model* m, IMP::ParticleIndexAdaptor p1,
        IMP::ParticleIndexAdaptor p2,
        const AVPairDistanceMeasurement& measurement, double sigma,
        double weight, double max_force)
    : IMP::Restraint(m, "AVMeanDistanceRestraint %1%"), p1_(p1), p2_(p2),
      measurement_(measurement),
      // The lookup runs to 2.5 R0. Past that the efficiency is 2e-4 and the
      // inverse is not a function any more, so extending the table buys
      // resolution in a region no measurement constrains.
      converter_(measurement.forster_radius, sigma, 1.0,
                 2.5 * measurement.forster_radius),
      weight_(weight), max_force_(max_force) {}

double ProbeAccessibleVolumeMeanDistanceRestraint::score_at(double d_mp) const {
    const double model =
            converter_.get_effective_distance(d_mp, measurement_.distance_type);
    if (max_force_ > 0.0) {
        // FPS's capped form, halved to stay in this module's units: a Gaussian
        // restraint is -log L = chi2/2, which is what score_model returns and
        // what every other term here is weighted against. The *knee* is at
        // MaxForce*err^2/2 either way -- the factor scales the slope, not the
        // place the parabola stops (see IMP::bff::chi2_score_capped).
        if (std::isnan(model)) return std::numeric_limits<double>::infinity();
        return 0.5 * chi2_score_capped(model, measurement_.distance,
                                       measurement_.error_neg,
                                       measurement_.error_pos, max_force_);
    }
    return measurement_.score_model(model);
}

double ProbeAccessibleVolumeMeanDistanceRestraint::unprotected_evaluate(
        IMP::DerivativeAccumulator* accum) const {
    IMP::core::XYZ d1(get_model(), p1_);
    IMP::core::XYZ d2(get_model(), p2_);
    const IMP::algebra::Vector3D r = d1.get_coordinates() - d2.get_coordinates();
    const double d_mp = r.get_magnitude();
    const double score = score_at(d_mp);

    if (accum != nullptr && d_mp > 1e-7) {
        // dS/dd_mp by central difference through the cached converter, then
        // onto the coordinates analytically: d(d_mp)/dx = +/- r_hat.
        const double eps = 1e-3;
        const double ds = (score_at(d_mp + eps) - score_at(d_mp - eps)) /
                          (2.0 * eps);
        const IMP::algebra::Vector3D grad = r * (weight_ * ds / d_mp);
        d1.add_to_derivatives(grad, *accum);
        d2.add_to_derivatives(-grad, *accum);
    }
    return weight_ * score;
}

void add_avs_to_rigid_bodies(const ProbeAccessibleVolumeDecorators& avs) {
    for (unsigned int i = 0; i < avs.size(); ++i) {
        ProbeAccessibleVolumeDecorator av = avs[i];
        av.resample();
        IMP::Particle* source = av.get_source();
        if (IMP::core::RigidBodyMember::get_is_setup(source)) {
            IMP::core::RigidBodyMember(source).get_rigid_body().add_member(
                    av.get_particle());
        }
    }
}

void set_av_xyzr_mass(const ProbeAccessibleVolumeDecorators& avs) {
    for (unsigned int i = 0; i < avs.size(); ++i) {
        ProbeAccessibleVolumeDecorator av = avs[i];
        IMP::Particle* p = av.get_particle();
        const IMP::algebra::Vector3D r = ProbeAccessibleVolumeDecorator(av).get_radii();
        const double radius = std::max(std::max(r[0], r[1]), r[2]);
        if (IMP::core::XYZR::get_is_setup(p)) {
            IMP::core::XYZR(p).set_radius(radius);
        } else {
            IMP::core::XYZR::setup_particle(p).set_radius(radius);
        }
        if (IMP::atom::Mass::get_is_setup(p)) {
            IMP::atom::Mass(p).set_mass(radius * 2.0);
        } else {
            IMP::atom::Mass::setup_particle(p, radius * 2.0);
        }
    }
}

AVFlatBottomRestraint::AVFlatBottomRestraint(
        IMP::Model* m, IMP::ParticleIndexAdaptor p1, IMP::ParticleIndexAdaptor p2,
        double r1, double r2, double r3, double r4, double k2, double k3,
        std::string name)
        : IMP::Restraint(m, name), p1_(p1), p2_(p2),
          r1_(r1), r2_(r2), r3_(r3), r4_(r4), k2_(k2), k3_(k3) {
    if (!(r1 <= r2 && r2 <= r3 && r3 <= r4)) {
        IMP_THROW("the bounds must be non-decreasing, got (" << r1 << ", " << r2
                                                             << ", " << r3
                                                             << ", " << r4
                                                             << ")",
                  ValueException);
    }
}

void AVFlatBottomRestraint::set_bounds(double r1, double r2, double r3,
                                       double r4) {
    if (!(r1 <= r2 && r2 <= r3 && r3 <= r4)) {
        IMP_THROW("the bounds must be non-decreasing, got (" << r1 << ", " << r2
                                                             << ", " << r3
                                                             << ", " << r4
                                                             << ")",
                  ValueException);
    }
    r1_ = r1; r2_ = r2; r3_ = r3; r4_ = r4;
}

std::vector<double> AVFlatBottomRestraint::get_bounds() const {
    std::vector<double> b(4);
    b[0] = r1_; b[1] = r2_; b[2] = r3_; b[3] = r4_;
    return b;
}

std::vector<double> AVFlatBottomRestraint::get_force_constants() const {
    std::vector<double> k(2);
    k[0] = k2_; k[1] = k3_;
    return k;
}

IMP::ParticleIndexes AVFlatBottomRestraint::get_particle_indexes() const {
    IMP::ParticleIndexes v;
    v.push_back(p1_);
    v.push_back(p2_);
    return v;
}

double AVFlatBottomRestraint::get_distance() const {
    return IMP::algebra::get_distance(
            IMP::core::XYZ(get_model(), p1_).get_coordinates(),
            IMP::core::XYZ(get_model(), p2_).get_coordinates());
}

double AVFlatBottomRestraint::unprotected_evaluate(
        IMP::DerivativeAccumulator* accum) const {
    IMP::core::XYZ a(get_model(), p1_), b(get_model(), p2_);
    const IMP::algebra::Vector3D v = a.get_coordinates() - b.get_coordinates();
    const double d = v.get_magnitude();

    double score = 0.0, dscore_dd = 0.0;
    if (d < r1_) {
        // linear continuation, tangent to the parabola at r1: the force is
        // capped rather than growing, which is what a first MD step needs
        const double w = r1_ - r2_;
        dscore_dd = 2.0 * k2_ * w;
        score = dscore_dd * (d - r1_) + k2_ * w * w;
    } else if (d < r2_) {
        const double w = d - r2_;
        score = k2_ * w * w;
        dscore_dd = 2.0 * k2_ * w;
    } else if (d <= r3_) {
        return 0.0;                     // inside the data: no score, no force
    } else if (d <= r4_) {
        const double w = d - r3_;
        score = k3_ * w * w;
        dscore_dd = 2.0 * k3_ * w;
    } else {
        const double w = r4_ - r3_;
        dscore_dd = 2.0 * k3_ * w;
        score = dscore_dd * (d - r4_) + k3_ * w * w;
    }

    if (accum != nullptr && d > 1e-9) {
        const IMP::algebra::Vector3D u = v / d;   // d(d)/da = u, d(d)/db = -u
        a.add_to_derivatives(u * dscore_dd, *accum);
        b.add_to_derivatives(-u * dscore_dd, *accum);
    }
    return score;
}

IMP::ModelObjectsTemp AVFlatBottomRestraint::do_get_inputs() const {
    IMP::ModelObjectsTemp ret;
    ret.push_back(get_model()->get_particle(p1_));
    ret.push_back(get_model()->get_particle(p2_));
    return ret;
}

MDRestraintSystem md_flat_bottom_restraints(
        const IMP::core::Hierarchy& hier, const std::string& fps_json,
        std::string score_set, double f_max, double tether_k,
        double probe_mass, std::string tether_atom) {
    if (!internal::file_exists(fps_json)) {
        IMP_THROW("fps.json not found: " << fps_json, IOException);
    }
    IMP::Model* m = hier.get_model();
    IMP::Pointer<ProbeNetworkRestraint> net(
            new ProbeNetworkRestraint(hier, fps_json, "md_flat_bottom",
                                      score_set));
    net->unprotected_evaluate(nullptr);          // build the volumes once
    const ProbeAccessibleVolumeDecorators avs = net->get_used_avs();

    IMP::Pointer<IMP::RestraintSet> rs(new IMP::RestraintSet(m, "AVFlatBottom"));
    AVFlatBottomRestraints wells;
    IMP::Vector<IMP::Pointer<IMP::core::Harmonic> > tethers;
    ProbeParticles probes;
    std::vector<std::string> pair_names;
    std::vector<double> bounds;

    // One probe per labelling position, tethered to its attachment atom at the
    // offset the volume itself has: the probe is where the dye is, and it
    // follows the site rather than the site following it.
    std::map<std::string, IMP::ParticleIndex> probe_of;
    for (unsigned int i = 0; i < avs.size(); ++i) {
        const ProbeAccessibleVolumeDecorator& av = avs[i];
        const std::string pname = av.get_particle()->get_name();
        IMP::Particle* p = new IMP::Particle(m, "probe_" + pname);
        const IMP::algebra::Vector3D mean =
                IMP::core::XYZ(av.get_particle()).get_coordinates();
        IMP::core::XYZR probe = IMP::core::XYZR::setup_particle(p);
        probe.set_coordinates(mean);
        probe.set_radius(std::max(1.0, av.get_radius1()));
        IMP::atom::Mass::setup_particle(p, probe_mass);
        // A probe an optimiser will not move is a probe the restraints cannot
        // act through: the well would push on nothing and the tether would
        // never transmit anything to the structure.
        probe.set_coordinates_are_optimized(true);

        IMP::ParticleIndex source = av.get_source()->get_index();
        if (!tether_atom.empty()) {
            // Tether to a named atom of the same residue instead: a probe
            // bound to an atom the optimiser does not move cannot follow.
            IMP::atom::Residue res = IMP::atom::get_residue(
                    IMP::atom::Atom(m, source));
            IMP::atom::Atom alt = IMP::atom::get_atom(
                    res, IMP::atom::AtomType(tether_atom));
            if (!alt) {
                IMP_THROW("residue " << res.get_index() << " has no atom "
                                     << tether_atom << " to tether "
                                     << pname << " to",
                          ValueException);
            }
            source = alt.get_particle_index();
        }
        const IMP::algebra::Vector3D at =
                IMP::core::XYZ(m, source).get_coordinates();
        const double offset = IMP::algebra::get_distance(mean, at);
        IMP::Pointer<IMP::core::Harmonic> tether(
                new IMP::core::Harmonic(offset, tether_k));
        rs->add_restraint(new IMP::core::DistanceRestraint(
                m, tether, m->get_particle(source), p));
        tethers.push_back(tether);

        ProbeParticle pp;
        pp.particle = p->get_index();
        pp.attachment = source;
        pp.position_name = pname;
        pp.tether_k = tether_k;
        probes.push_back(pp);
        probe_of[pname] = p->get_index();
    }

    const std::map<std::string, AVPairDistanceMeasurement> distances =
            net->get_used_distances();
    std::map<std::string, ProbeAccessibleVolumeDecorator> by_name;
    for (unsigned int i = 0; i < avs.size(); ++i) {
        by_name[avs[i].get_particle()->get_name()] = avs[i];
    }
    for (auto it = distances.begin(); it != distances.end(); ++it) {
        const AVPairDistanceMeasurement& d = it->second;
        auto a1 = by_name.find(d.position_1), a2 = by_name.find(d.position_2);
        if (a1 == by_name.end() || a2 == by_name.end()) {
            IMP_THROW("distance " << it->first << " names a position the "
                                  << "network has no volume for",
                      ValueException);
        }
        const std::vector<double> b =
                rmp_flat_bottom_bounds(a1->second, a2->second, d);
        const double lo = b[0], centre = b[1], hi = b[2];
        // k = f_max / (2 sigma): one error bar out of the well costs f_max,
        // and past r1/r4 the force stops growing.
        const double s_lo = std::max(1e-6, centre - lo);
        const double s_hi = std::max(1e-6, hi - centre);
        const double k2 = f_max / (2.0 * s_lo), k3 = f_max / (2.0 * s_hi);
        IMP::Pointer<AVFlatBottomRestraint> well(new AVFlatBottomRestraint(
                m, probe_of[d.position_1], probe_of[d.position_2],
                std::max(0.0, lo - s_lo), lo, hi, hi + s_hi, k2, k3,
                it->first));
        rs->add_restraint(well);
        wells.push_back(well);
        pair_names.push_back(it->first);
        bounds.push_back(std::max(0.0, lo - s_lo));
        bounds.push_back(lo);
        bounds.push_back(hi);
        bounds.push_back(hi + s_hi);
    }
    return MDRestraintSystem(rs, net, wells, tethers, probes, pair_names, bounds);
}

namespace {
//! A -> nm, and kcal/mol/A^2 -> kJ/mol/nm^2. OpenMM's units, stated once.
const double ANGSTROM_TO_NM = 0.1;
const double KCAL_PER_A2_TO_KJ_PER_NM2 = 4.184 * 100.0;
}  // namespace

std::string openmm_flat_bottom_energy() {
    // The five branches of AVFlatBottomRestraint, innermost last. `select(c,
    // a, b)` is a when c is non-zero; `step(x)` is 1 for x >= 0. Written so
    // that every branch is evaluated safely for any r -- OpenMM evaluates both
    // arms of a select.
    return "select(step(r1-r), k2*(r1-r2)^2 + 2*k2*(r1-r2)*(r-r1),"
           "select(step(r2-r), k2*(r-r2)^2,"
           "select(step(r-r4), k3*(r4-r3)^2 + 2*k3*(r4-r3)*(r-r4),"
           "select(step(r-r3), k3*(r-r3)^2, 0))))";
}

void write_openmm_restraints(const MDRestraintSystem& system,
                             const std::string& path,
                             const IMP::core::Hierarchy& hier) {
    IMP::Model* m = hier.get_model();
    nlohmann::json doc;
    doc["format"] = "IMP.bff OpenMM flat-bottom FRET restraints";
    doc["version"] = 1;
    doc["units"]["length"] = "nanometer";
    doc["units"]["energy"] = "kilojoule_per_mole";
    doc["units"]["mass"] = "dalton";
    doc["energy_expression"] = openmm_flat_bottom_energy();
    doc["per_bond_parameters"] = {"r1", "r2", "r3", "r4", "k2", "k3"};

    std::map<IMP::ParticleIndex, std::string> probe_name_of;
    for (unsigned int i = 0; i < system.get_probes().size(); ++i) {
        const ProbeParticle& p = system.get_probes()[i];
        probe_name_of[p.particle] = p.position_name;
        const IMP::algebra::Vector3D xyz =
                IMP::core::XYZ(m, p.particle).get_coordinates();

        nlohmann::json j;
        j["name"] = p.position_name;
        j["mass"] = IMP::atom::Mass(m, p.particle).get_mass();
        j["position"] = {xyz[0] * ANGSTROM_TO_NM, xyz[1] * ANGSTROM_TO_NM,
                         xyz[2] * ANGSTROM_TO_NM};

        // named, not indexed: an index depends on how the reader built its
        // topology, and a spec that names atoms by position in someone else's
        // file silently restrains the wrong atoms
        IMP::atom::Atom at(m, p.attachment);
        IMP::atom::Residue res = IMP::atom::get_residue(at);
        IMP::atom::Chain chain = IMP::atom::get_chain(res);
        j["attachment"]["chain"] = chain ? chain.get_id() : std::string();
        j["attachment"]["residue"] = res.get_index();
        j["attachment"]["atom"] = at.get_atom_type().get_string();
        j["tether"]["length"] = system.get_tether_length(i) * ANGSTROM_TO_NM;
        j["tether"]["k"] = p.tether_k * KCAL_PER_A2_TO_KJ_PER_NM2;
        doc["probes"].push_back(j);
    }

    const AVFlatBottomRestraints& wells = system.get_wells();
    const std::vector<double>& b = system.get_bounds();
    for (unsigned int i = 0; i < wells.size(); ++i) {
        const std::vector<double> k = wells[i]->get_force_constants();
        nlohmann::json j;
        j["name"] = system.get_pair_names()[i];
        const IMP::ParticleIndexes ends = wells[i]->get_particle_indexes();
        j["probe_1"] = probe_name_of.count(ends[0]) ? probe_name_of[ends[0]]
                                                    : std::string();
        j["probe_2"] = probe_name_of.count(ends[1]) ? probe_name_of[ends[1]]
                                                    : std::string();
        j["r1"] = b[4 * i + 0] * ANGSTROM_TO_NM;
        j["r2"] = b[4 * i + 1] * ANGSTROM_TO_NM;
        j["r3"] = b[4 * i + 2] * ANGSTROM_TO_NM;
        j["r4"] = b[4 * i + 3] * ANGSTROM_TO_NM;
        j["k2"] = k[0] * KCAL_PER_A2_TO_KJ_PER_NM2;
        j["k3"] = k[1] * KCAL_PER_A2_TO_KJ_PER_NM2;
        doc["restraints"].push_back(j);
    }

    std::ofstream out(path.c_str());
    if (!out) IMP_THROW("cannot write " << path, IOException);
    out << doc.dump(2) << "\n";
}

AVRebuildOptimizerState::AVRebuildOptimizerState(
        IMP::Model* m, const MDRestraintSystem& system, unsigned int period)
        : IMP::OptimizerState(m, "AVRebuildOptimizerState%1%"),
          network_(system.get_network()), wells_(system.get_wells()),
          probes_(system.get_probes()), tethers_(system.get_tethers()),
          pair_names_(system.get_pair_names()), n_updates_(0) {
    set_period(period);
}

void AVRebuildOptimizerState::update_now() {
    IMP::Model* m = get_model();
    // Rebuild every volume at the coordinates the trajectory has reached.
    network_->unprotected_evaluate(nullptr);
    const ProbeAccessibleVolumeDecorators avs = network_->get_used_avs();
    std::map<std::string, ProbeAccessibleVolumeDecorator> by_name;
    for (unsigned int i = 0; i < avs.size(); ++i) {
        by_name[avs[i].get_particle()->get_name()] = avs[i];
    }

    // A volume's offset from its attachment moves with its shape, so the
    // tether that holds the probe at that offset has to move with it.
    for (unsigned int i = 0; i < probes_.size() && i < tethers_.size(); ++i) {
        std::map<std::string, ProbeAccessibleVolumeDecorator>::const_iterator it =
                by_name.find(probes_[i].position_name);
        if (it == by_name.end()) continue;
        const IMP::algebra::Vector3D mean =
                IMP::core::XYZ(it->second.get_particle()).get_coordinates();
        const IMP::algebra::Vector3D at =
                IMP::core::XYZ(m, probes_[i].attachment).get_coordinates();
        const double offset = IMP::algebra::get_distance(mean, at);
        tethers_[i]->set_mean(offset);
    }

    // And the R_mp that reproduces a measured <R_DA> moves with both shapes.
    const std::map<std::string, AVPairDistanceMeasurement> distances =
            network_->get_used_distances();
    for (unsigned int i = 0; i < wells_.size() && i < pair_names_.size(); ++i) {
        std::map<std::string, AVPairDistanceMeasurement>::const_iterator d =
                distances.find(pair_names_[i]);
        if (d == distances.end()) continue;
        std::map<std::string, ProbeAccessibleVolumeDecorator>::const_iterator a1 =
                by_name.find(d->second.position_1);
        std::map<std::string, ProbeAccessibleVolumeDecorator>::const_iterator a2 =
                by_name.find(d->second.position_2);
        if (a1 == by_name.end() || a2 == by_name.end()) continue;
        std::vector<double> b;
        try {
            b = rmp_flat_bottom_bounds(a1->second, a2->second, d->second);
        } catch (const IMP::ValueException&) {
            // The structure has moved somewhere the measurement cannot be
            // reproduced at any separation. Leaving the well where it was is
            // the conservative choice: it still pulls the right way.
            continue;
        }
        const double lo = b[0], centre = b[1], hi = b[2];
        const double s_lo = std::max(1e-6, centre - lo);
        const double s_hi = std::max(1e-6, hi - centre);
        wells_[i]->set_bounds(std::max(0.0, lo - s_lo), lo, hi, hi + s_hi);
    }
    ++n_updates_;
}

void AVRebuildOptimizerState::do_update(unsigned int) { update_now(); }

void write_openmm_script(const MDRestraintSystem& system,
                         const std::string& path,
                         const IMP::core::Hierarchy& hier,
                         const std::string& pdb_path,
                         const std::string& force_field, int n_steps) {
    IMP::Model* m = hier.get_model();
    std::ofstream out(path.c_str());
    if (!out) IMP_THROW("cannot write " << path, IOException);

    out << "#!/usr/bin/env python\n"
        << "\"\"\"FRET-restrained OpenMM run, generated by IMP.bff.\n\n"
           "The restraint table below is embedded, so this file is the whole\n"
           "input: the flat-bottom wells, the probe particles that carry them\n"
           "and the tethers holding those probes to their labelling sites.\n\n"
           "Lengths are nanometres and energies kJ/mol -- OpenMM's units, not\n"
           "the angstrom and kcal/mol the restraints were derived in.\n\n"
           "What to change first: the force field, the solvent, the integrator\n"
           "and the run length, all at the top. This is a setup that runs, not\n"
           "a production protocol.\n\"\"\"\n"
        << "import openmm\n"
        << "import openmm.app as app\n"
        << "import openmm.unit as unit\n\n"
        << "PDB = \"" << pdb_path << "\"\n"
        << "FORCE_FIELD = \"" << force_field << "\".split()\n"
        << "N_STEPS = " << n_steps << "\n"
        << "TEMPERATURE = 300.0 * unit.kelvin\n"
        << "TIME_STEP = 2.0 * unit.femtosecond\n\n"
        << "# The flat-bottom well, as IMP.bff evaluates it. One string, so the\n"
        << "# restraint that runs here is the restraint that ran there.\n"
        << "ENERGY = (\"" << openmm_flat_bottom_energy() << "\")\n"
        << "PARAMETERS = [\"r1\", \"r2\", \"r3\", \"r4\", \"k2\", \"k3\"]\n\n"
        << "# name, chain, residue, atom, mass, tether length, tether k\n"
        << "PROBES = [\n";
    for (unsigned int i = 0; i < system.get_probes().size(); ++i) {
        const ProbeParticle& p = system.get_probes()[i];
        IMP::atom::Atom at(m, p.attachment);
        IMP::atom::Residue res = IMP::atom::get_residue(at);
        IMP::atom::Chain chain = IMP::atom::get_chain(res);
        const IMP::algebra::Vector3D xyz =
                IMP::core::XYZ(m, p.particle).get_coordinates();
        out << "    (\"" << p.position_name << "\", \""
            << (chain ? chain.get_id() : std::string()) << "\", "
            << res.get_index() << ", \"" << at.get_atom_type().get_string()
            << "\", " << IMP::atom::Mass(m, p.particle).get_mass() << ", "
            << system.get_tether_length(i) * 0.1 << ", "
            << p.tether_k * 418.4 << ", ("
            << xyz[0] * 0.1 << ", " << xyz[1] * 0.1 << ", " << xyz[2] * 0.1
            << ")),\n";
    }
    out << "]\n\n# name, probe 1, probe 2, r1, r2, r3, r4, k2, k3\nWELLS = [\n";
    const std::vector<double>& b = system.get_bounds();
    std::map<IMP::ParticleIndex, std::string> probe_name_of;
    for (unsigned int i = 0; i < system.get_probes().size(); ++i) {
        probe_name_of[system.get_probes()[i].particle] =
                system.get_probes()[i].position_name;
    }
    for (unsigned int i = 0; i < system.get_wells().size(); ++i) {
        const IMP::ParticleIndexes ends =
                system.get_wells()[i]->get_particle_indexes();
        const std::vector<double> k =
                system.get_wells()[i]->get_force_constants();
        out << "    (\"" << system.get_pair_names()[i] << "\", \""
            << probe_name_of[ends[0]] << "\", \"" << probe_name_of[ends[1]]
            << "\", " << b[4 * i + 0] * 0.1 << ", " << b[4 * i + 1] * 0.1
            << ", " << b[4 * i + 2] * 0.1 << ", " << b[4 * i + 3] * 0.1
            << ", " << k[0] * 418.4 << ", " << k[1] * 418.4 << "),\n";
    }
    out << "]\n\n"
        << "pdb = app.PDBFile(PDB)\n"
        << "modeller = app.Modeller(pdb.topology, pdb.positions)\n"
        << "system = app.ForceField(*FORCE_FIELD).createSystem(\n"
        << "    modeller.topology, nonbondedMethod=app.NoCutoff,\n"
        << "    constraints=app.HBonds)\n\n"
        << "# Attachment atoms are found by name, not by index: an index would\n"
        << "# depend on how this file's topology was built -- hydrogens,\n"
        << "# waters, altlocs -- and would silently restrain the wrong atoms.\n"
        << "index_of = {(c.id, r.id, a.name): a.index\n"
        << "            for c in modeller.topology.chains()\n"
        << "            for r in c.residues() for a in r.atoms()}\n\n"
        << "positions = list(modeller.positions)\n"
        << "probe_index, tether = {}, openmm.HarmonicBondForce()\n"
        << "for name, chain, residue, atom, mass, length, k, xyz in PROBES:\n"
        << "    key = (chain, str(residue), atom)\n"
        << "    if key not in index_of:\n"
        << "        raise SystemExit(f\"{key} is not in {PDB}\")\n"
        << "    probe_index[name] = system.addParticle(mass * unit.dalton)\n"
        << "    positions.append(openmm.Vec3(*xyz) * unit.nanometer)\n"
        << "    tether.addBond(index_of[key], probe_index[name], length, k)\n"
        << "system.addForce(tether)\n\n"
        << "wells = openmm.CustomBondForce(ENERGY)\n"
        << "for parameter in PARAMETERS:\n"
        << "    wells.addPerBondParameter(parameter)\n"
        << "for name, p1, p2, *values in WELLS:\n"
        << "    wells.addBond(probe_index[p1], probe_index[p2], values)\n"
        << "system.addForce(wells)\n"
        << "print(f\"{len(PROBES)} probes, {wells.getNumBonds()} FRET wells\")\n\n"
        << "integrator = openmm.LangevinMiddleIntegrator(\n"
        << "    TEMPERATURE, 1.0 / unit.picosecond, TIME_STEP)\n"
        << "simulation = app.Simulation(modeller.topology, system, integrator)\n"
        << "simulation.context.setPositions(positions)\n"
        << "simulation.minimizeEnergy()\n"
        << "simulation.reporters.append(app.StateDataReporter(\n"
        << "    \"run.log\", 1000, step=True, potentialEnergy=True,\n"
        << "    temperature=True))\n"
        << "simulation.reporters.append(app.DCDReporter(\"run.dcd\", 1000))\n"
        << "simulation.step(N_STEPS)\n\n"
        << "# What the run made of the data: a distance is satisfied when it\n"
        << "# sits inside its flat bottom.\n"
        << "state = simulation.context.getState(getPositions=True)\n"
        << "final = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)\n"
        << "satisfied = 0\n"
        << "for name, p1, p2, r1, r2, r3, r4, k2, k3 in WELLS:\n"
        << "    d = ((final[probe_index[p1]] - final[probe_index[p2]]) ** 2).sum() ** 0.5\n"
        << "    satisfied += r2 <= d <= r3\n"
        << "print(f\"{satisfied} of {len(WELLS)} distances satisfied\")\n";
}

IMP::RestraintSet* probe_network_restraint_set(
        const IMP::core::Hierarchy& hier, const std::string& fps_json,
        std::string name, std::string score_set,
        bool mean_position_restraint, double sigma_DA, bool occupy_volume,
        double weight) {
    if (!internal::file_exists(fps_json)) {
        IMP_THROW("fps.json not found: " << fps_json, IOException);
    }
    IMP::Model* m = hier.get_model();
    IMP::Pointer<IMP::RestraintSet> rs(
            new IMP::RestraintSet(m, "ProbeNetworkRestraint"));
    IMP::Pointer<ProbeNetworkRestraint> net(
            new ProbeNetworkRestraint(hier, fps_json, name, score_set));
    const ProbeAccessibleVolumeDecorators avs = net->get_used_avs();

    if (!mean_position_restraint) {
        rs->add_restraint(net);
    } else {
        // The volumes are carried by the rigid bodies of the atoms they are
        // attached to, which is what makes scoring their mean positions --
        // rather than rebuilding them -- an approximation and not a fiction.
        add_avs_to_rigid_bodies(avs);
        std::map<std::string, ProbeAccessibleVolumeDecorator> by_name;
        for (unsigned int i = 0; i < avs.size(); ++i) {
            by_name[avs[i].get_particle()->get_name()] = avs[i];
        }
        const std::map<std::string, AVPairDistanceMeasurement> distances =
                net->get_used_distances();
        for (std::map<std::string, AVPairDistanceMeasurement>::const_iterator
                     it = distances.begin();
             it != distances.end(); ++it) {
            const AVPairDistanceMeasurement& d = it->second;
            std::map<std::string, ProbeAccessibleVolumeDecorator>::const_iterator a1 =
                    by_name.find(d.position_1);
            std::map<std::string, ProbeAccessibleVolumeDecorator>::const_iterator a2 =
                    by_name.find(d.position_2);
            IMP_USAGE_CHECK(a1 != by_name.end() && a2 != by_name.end(),
                            "distance " << it->first << " names a position "
                                        << "the network has no volume for");
            rs->add_restraint(new ProbeAccessibleVolumeMeanDistanceRestraint(
                    m, a1->second.get_particle(), a2->second.get_particle(), d,
                    sigma_DA));
        }
    }
    if (occupy_volume) set_av_xyzr_mass(avs);
    rs->set_weight(weight);
    return rs.release();
}

IMP::ModelObjectsTemp ProbeAccessibleVolumeMeanDistanceRestraint::do_get_inputs() const {
    IMP::ModelObjectsTemp out;
    out.push_back(get_model()->get_particle(p1_));
    out.push_back(get_model()->get_particle(p2_));
    return out;
}

IMPBFF_END_NAMESPACE
