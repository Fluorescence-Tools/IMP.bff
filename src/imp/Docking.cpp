/**
 * \file Docking.cpp
 * \brief What a docking run is told, and what it reports.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/Docking.h>
// olga_vdw_particle_radius
#include <IMP/bff/HierarchyBridge.h>

#include <IMP/bff/States.h>
#include <IMP/bff/AVMeanDistanceRestraint.h>
#include <IMP/bff/FPS.h>
#include <IMP/bff/internal/json.h>

#include <IMP/algebra/geometric_alignment.h>
#include <IMP/atom/Atom.h>
#include <IMP/atom/Chain.h>
#include <IMP/atom/Mass.h>
#include <IMP/atom/Residue.h>
#include <IMP/atom/pdb.h>
#include <IMP/container/ListSingletonContainer.h>
#include <IMP/algebra/vector_generators.h>
#include <IMP/core/ConjugateGradients.h>
#include <IMP/core/ExcludedVolumeRestraint.h>
#include <IMP/core/XYZR.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/bff/Base.h>
#include <IMP/bff/internal/Text.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <unistd.h>
#include <set>
#include <fstream>
#include <iomanip>
#include <limits>
#include <random>
#include <sstream>

IMPBFF_BEGIN_NAMESPACE

namespace {

//! `round(x, n)` as the CSV writer spelled it.
double rounded(double v, int digits) {
    if (!std::isfinite(v)) return v;
    const double scale = std::pow(10.0, digits);
    return std::floor(v * scale + 0.5) / scale;
}

//! A CSV cell: numbers plainly, and text quoted only when it has to be.
std::string cell(const std::string& text) {
    if (text.find_first_of(",\"\n") == std::string::npos) return text;
    std::string out = "\"";
    for (std::size_t i = 0; i < text.size(); ++i) {
        if (text[i] == '"') out += '"';
        out += text[i];
    }
    return out + "\"";
}

std::string number(double v) {
    if (std::isnan(v)) return "";
    std::ostringstream out;
    out << v;
    return out.str();
}

}  // namespace

double PairDistance::get_chi2() const {
    return chi2_score(distance_model, distance_exp, error_neg, error_pos);
}

double PairDistance::get_efficiency_model() const {
    return fret_efficiency(distance_model, forster_radius);
}

double PairDistance::get_efficiency_exp() const {
    return fret_efficiency(distance_exp, forster_radius);
}

std::string PairDistance::get_json() const {
    nlohmann::json out;
    out["name"] = name;
    out["position1"] = position1;
    out["position2"] = position2;
    out["distance_exp"] = distance_exp;
    out["distance_model"] = distance_model;
    out["distance_type"] = distance_type;
    out["forster_radius"] = forster_radius;
    out["error_neg"] = error_neg;
    out["error_pos"] = error_pos;
    out["is_bond"] = is_bond;
    out["residual"] = get_residual();
    out["chi2"] = get_chi2();
    out["E_exp"] = get_efficiency_exp();
    out["E_model"] = get_efficiency_model();
    return out.dump();
}

std::string DockingResult::get_json() const {
    nlohmann::json out;
    out["score"] = score;
    out["n_avs"] = n_avs;
    out["n_distances"] = n_distances;
    out["e_bond"] = e_bond;
    out["n_bonds"] = n_bonds;
    out["e_clash"] = e_clash;
    out["converged"] = converged;
    nlohmann::json rows = nlohmann::json::array();
    for (std::size_t i = 0; i < pairs.size(); ++i) {
        rows.push_back(nlohmann::json::parse(pairs[i].get_json()));
    }
    out["pairs"] = rows;
    out["output_dir"] = output_dir;
    out["rmf_file"] = rmf_file;
    out["best_pdbs"] = best_pdbs;
    out["score_csv"] = score_csv;
    out["extra"] = nlohmann::json::parse(extra.empty() ? "{}" : extra);
    out["poses"] = nlohmann::json::parse(poses.empty() ? "[]" : poses);
    return out.dump();
}

namespace {

//! The fixed part of a `PairDistance`: everything but the model distance.
PairDistance from_measurement(const std::string& name,
                              const AVPairDistanceMeasurement& e) {
    PairDistance pair;
    pair.name = name;
    pair.position1 = e.position_1;
    pair.position2 = e.position_2;
    pair.distance_exp = e.distance;
    pair.forster_radius = e.forster_radius;
    pair.error_neg = e.error_neg;
    pair.error_pos = e.error_pos;
    pair.distance_type = probe_pair_distance_type_name(e.distance_type);
    pair.distance_model = std::numeric_limits<double>::quiet_NaN();
    return pair;
}

//! The modelled observable a mean-position separation stands for.
/*! Not #IMP::bff::effective_distance, which converts by a named transfer
    function; this one goes through the same #IMP::bff::FRETDistanceConverter
    the mean-position restraint scores with. */
double modelled_distance(double d_mp, const AVPairDistanceMeasurement& e,
                         double sigma_da) {
    const FRETDistanceConverter converter(e.forster_radius, sigma_da, 1.0,
                                          2.5 * e.forster_radius);
    return converter.get_effective_distance(d_mp, e.distance_type);
}

}  // namespace

std::vector<PairDistance> pair_distances_at_positions(
        const std::map<std::string, AVPairDistanceMeasurement>& distances,
        IMP::Model* m, const std::vector<std::string>& position_names,
        const IMP::ParticleIndexes& particles, double sigma_da) {
    std::map<std::string, IMP::ParticleIndex> by_name;
    const std::size_t n = std::min(position_names.size(), particles.size());
    for (std::size_t i = 0; i < n; ++i) {
        by_name.insert(std::make_pair(position_names[i], particles[i]));
    }

    std::vector<PairDistance> out;
    for (std::map<std::string, AVPairDistanceMeasurement>::const_iterator it =
                 distances.begin();
         it != distances.end(); ++it) {
        const AVPairDistanceMeasurement& e = it->second;
        PairDistance pair = from_measurement(it->first, e);
        std::map<std::string, IMP::ParticleIndex>::const_iterator p1 =
                by_name.find(e.position_1);
        std::map<std::string, IMP::ParticleIndex>::const_iterator p2 =
                by_name.find(e.position_2);
        if (p1 != by_name.end() && p2 != by_name.end()) {
            const IMP::algebra::Vector3D r =
                    IMP::core::XYZ(m, p1->second).get_coordinates() -
                    IMP::core::XYZ(m, p2->second).get_coordinates();
            pair.distance_model = modelled_distance(r.get_magnitude(), e,
                                                     sigma_da);
        }
        out.push_back(pair);
    }
    return out;
}

std::map<std::string, AVPairDistanceMeasurement> scored_measurements(
        ProbeNetworkRestraint* network) {
    std::map<std::string, AVPairDistanceMeasurement> out;
    if (network == NULL) return out;
    out = network->get_used_distances();
    for (std::map<std::string, AVPairDistanceMeasurement>::iterator it =
                 out.begin();
         it != out.end(); ++it) {
        if (network->get_position_is_point(it->second.position_1) ||
            network->get_position_is_point(it->second.position_2)) {
            it->second.distance_type = PROBE_PAIR_DISTANCE_MP;
        }
    }
    return out;
}

std::vector<PairDistance> collect_pair_distances(ProbeNetworkRestraint* restraint,
                                                 bool mean_position,
                                                 double sigma_da) {
    std::vector<PairDistance> out;
    if (restraint == NULL) return out;

    const std::map<std::string, AVPairDistanceMeasurement> used =
            scored_measurements(restraint);
    for (std::map<std::string, AVPairDistanceMeasurement>::const_iterator it =
                 used.begin();
         it != used.end(); ++it) {
        const AVPairDistanceMeasurement& e = it->second;
        PairDistance pair = from_measurement(it->first, e);
        pair.is_bond = restraint->get_is_bond(it->first);
        if (mean_position) {
            // One lookup for both kinds of position: an AV's particle carries
            // its mean position, an `ATOM` position's particle *is* the atom
            // and an `XYZ` position's is the fixed point. This used to key on
            // the AV particle's *name*, which meant a point position simply
            // had no entry and its distance came back NaN.
            const IMP::ParticleIndex p1 =
                    restraint->get_position_particle_index(e.position_1);
            const IMP::ParticleIndex p2 =
                    restraint->get_position_particle_index(e.position_2);
            if (p1 != IMP::ParticleIndex() && p2 != IMP::ParticleIndex()) {
                IMP::Model* m = restraint->get_model();
                const IMP::algebra::Vector3D r =
                        IMP::core::XYZ(m, p1).get_coordinates() -
                        IMP::core::XYZ(m, p2).get_coordinates();
                // For a point end `e.distance_type` is already `Rmp`, so the
                // transfer function returns the separation untouched.
                pair.distance_model = modelled_distance(r.get_magnitude(), e,
                                                         sigma_da);
            }
        } else {
            pair.distance_model = restraint->get_model_distance(
                    e.position_1, e.position_2, e.forster_radius,
                    e.distance_type);
        }
        out.push_back(pair);
    }
    return out;
}

void write_score_csv(const std::string& path, double score,
                     const std::vector<PairDistance>& pairs) {
    std::ofstream out(path.c_str());
    if (!out) IMP_THROW("cannot write " << path, IOException);
    out << "# total_score," << score << "\n";
    // `is_bond` last, so a reader keyed on the earlier columns is unaffected.
    // It is a column rather than a silence because a crosslink in a table of
    // FRET pairs is exactly what a reader would otherwise mistake for one.
    out << "name,position1,position2,distance_type,forster_radius,"
           "distance_exp,distance_model,residual,chi2,E_exp,E_model,is_bond\n";
    for (std::size_t i = 0; i < pairs.size(); ++i) {
        const PairDistance& p = pairs[i];
        out << cell(p.name) << "," << cell(p.position1) << ","
            << cell(p.position2) << "," << cell(p.distance_type) << ","
            << number(p.forster_radius) << "," << number(p.distance_exp) << ","
            << number(p.distance_model) << ","
            << number(rounded(p.get_residual(), 3)) << ","
            << number(rounded(p.get_chi2(), 3)) << ","
            << number(rounded(p.get_efficiency_exp(), 4)) << ","
            << number(rounded(p.get_efficiency_model(), 4)) << ","
            << (p.is_bond ? 1 : 0) << "\n";
    }
}

void write_screening_pairs_csv(
        const std::string& path,
        const std::vector<ScreenedStructure>& structures) {
    std::ofstream out(path.c_str());
    if (!out) IMP_THROW("cannot write " << path, IOException);
    out << "pdb,name,position1,position2,distance_type,forster_radius,"
           "distance_exp,distance_model,residual,chi2\n";
    for (std::size_t i = 0; i < structures.size(); ++i) {
        const ScreenedStructure& s = structures[i];
        for (std::size_t k = 0; k < s.pairs.size(); ++k) {
            const PairDistance& p = s.pairs[k];
            out << cell(s.path) << "," << cell(p.name) << ","
                << cell(p.position1) << "," << cell(p.position2) << ","
                << cell(p.distance_type) << "," << number(p.forster_radius)
                << "," << number(p.distance_exp) << ","
                << number(rounded(p.distance_model, 4)) << ","
                << number(rounded(p.get_residual(), 3)) << ","
                << number(rounded(p.get_chi2(), 3)) << "\n";
        }
    }
}

DockingAssembly::DockingAssembly(
        IMP::Model* model, IMP::atom::Hierarchy root,
        ProbeNetworkRestraint* network, IMP::RestraintSet* restraints,
        IMP::core::RestraintsScoringFunction* scoring_function,
        const IMP::core::RigidBodies& bodies,
        const std::vector<int>& body_of_pdb, bool mean_position,
        double sigma_da, const std::string& fps_json_path, double max_force)
    : model_(model), root_(root), network_(network), restraints_(restraints),
      scoring_function_(scoring_function), bodies_(bodies),
      body_of_pdb_(body_of_pdb), mean_position_(mean_position),
      sigma_da_(sigma_da), fps_json_path_(fps_json_path),
      max_force_(max_force) {}

namespace {

/*! \brief Links an atom to the sphere the excluded-volume term measures it by.
 *
 * Absent on every atom the clash term measures directly, which is the default
 * and today's behaviour. Present only when `clash_container` was asked for a
 * radii source or a scale that is not the particles' own, in which case it
 * points at the shadow sphere that carries the substituted radius.
 *
 * The reverse direction (shadow -> atom) is deliberately not stored: the only
 * thing that needs the link is `set_bond_anchor_radii`, which starts from the
 * anchor atom and has to reach its sphere.
 */
IMP::ParticleIndexKey docking_clash_shadow_key() {
    static const IMP::ParticleIndexKey k("bff_clash_shadow");
    return k;
}

}  // namespace

int set_bond_anchor_radii(ProbeNetworkRestraint* network, double radius) {
    if (network == NULL) return 0;
    IMP::Model* m = network->get_model();
    const std::map<std::string, AVPairDistanceMeasurement> used =
            network->get_used_distances();
    std::set<IMP::ParticleIndex> anchors;
    for (std::map<std::string, AVPairDistanceMeasurement>::const_iterator it =
                 used.begin();
         it != used.end(); ++it) {
        if (!network->get_is_bond(it->first)) continue;
        anchors.insert(
                network->get_position_particle_index(it->second.position_1));
        anchors.insert(
                network->get_position_particle_index(it->second.position_2));
    }
    int changed = 0;
    for (std::set<IMP::ParticleIndex>::const_iterator it = anchors.begin();
         it != anchors.end(); ++it) {
        if (*it == IMP::ParticleIndex()) continue;
        if (!IMP::core::XYZR::get_is_setup(m, *it)) continue;
        IMP::core::XYZR(m, *it).set_radius(radius);
        // The clash term may be measuring this atom by a shadow sphere on a
        // substituted radii set (`clash_container` with a non-default source
        // or scale). The bond rule is about the *clash*, so it has to reach
        // the sphere that clashes; leaving it out would let an anchor keep a
        // full Bondi radius and repel the atom it is bonded to.
        if (m->get_has_attribute(docking_clash_shadow_key(), *it)) {
            const IMP::ParticleIndex s =
                    m->get_attribute(docking_clash_shadow_key(), *it);
            if (IMP::core::XYZR::get_is_setup(m, s)) {
                IMP::core::XYZR(m, s).set_radius(radius);
            }
        }
        ++changed;
    }
    return changed;
}

namespace {

//! Which `body_id` the \p idx-th PDB is.
/*! When the file names as many distinct bodies as there are structures, the
    two are matched in declaration order; otherwise a structure is its own
    body. */
int body_for_pdb_index(int idx, std::size_t n_pdb,
                       const nlohmann::json& positions) {
    std::set<int> ids;
    for (nlohmann::json::const_iterator it = positions.begin();
         it != positions.end(); ++it) {
        if (it.value().is_object()) {
            ids.insert(it.value().value("body_id", 0));
        }
    }
    if (ids.size() == n_pdb) {
        std::vector<int> sorted(ids.begin(), ids.end());
        return sorted[idx];
    }
    return idx;
}

}  // namespace

DockingAssembly create_docking_assembly(
        const std::vector<std::string>& pdb_paths,
        const std::string& fps_json_path, const std::string& score_set,
        bool mean_position_restraint, double ev_weight, double sigma_da,
        double clash_tolerance, double max_force,
        const std::string& clash_radii_source, double clash_radii_scale) {
    // Validated before a file is opened: an unknown radii source is a typo in
    // a command line, and reporting it after a two-minute AV calculation is a
    // worse error message for the same mistake.
    const AVRadiiSource clash_radii =
            av_radii_source_from_string(clash_radii_source);
    if (clash_radii_scale <= 0.0) {
        IMP_THROW("clash radii scale must be positive, got "
                          << clash_radii_scale,
                  ValueException);
    }
    if (pdb_paths.empty()) {
        IMP_THROW("At least one PDB file is required.", ValueException);
    }
    for (std::size_t i = 0; i < pdb_paths.size(); ++i) {
        if (!internal::file_exists(pdb_paths[i])) {
            IMP_THROW("PDB file not found: " << pdb_paths[i], IOException);
        }
    }
    if (!internal::file_exists(fps_json_path)) {
        IMP_THROW("fps.json not found: " << fps_json_path, IOException);
    }

    // A path that is not `.json` is the original FPS / C# labelling file;
    // `read_fps_json` reads either, but `ProbeNetworkRestraint` opens the file
    // itself, so a legacy input is converted once and the network is given
    // the converted copy.
    std::string fps_path = fps_json_path;
    if (!internal::ends_with(fps_path, ".json")) {
        const FPSDocument legacy = read_fps_json(fps_json_path, pdb_paths);
        fps_path = fps_json_path + ".fps.json";
        write_fps_json(fps_path, legacy.positions, legacy.distances,
                       legacy.score_sets, legacy.extra);
    }
    const FPSDocument document = read_fps_json(fps_path);
    nlohmann::json positions = nlohmann::json::parse(document.positions, NULL,
                                                     false);
    if (positions.is_discarded()) positions = nlohmann::json::object();

    IMP_NEW(IMP::Model, model, ());
    IMP::atom::Hierarchy root = IMP::atom::Hierarchy::setup_particle(
            new IMP::Particle(model, "root"));

    IMP::core::RigidBodies bodies;
    std::vector<int> body_of_pdb;
    for (std::size_t i = 0; i < pdb_paths.size(); ++i) {
        const int body_id = body_for_pdb_index(static_cast<int>(i),
                                               pdb_paths.size(), positions);
        IMP::atom::Hierarchy h = IMP::atom::read_pdb(
                pdb_paths[i], model,
                new IMP::atom::NonWaterNonHydrogenPDBSelector());
        root.add_child(h);
        body_of_pdb.push_back(body_id);
        if (IMP::atom::get_leaves(h).empty()) continue;
        IMP::core::RigidBody rb = IMP::atom::create_rigid_body(h);
        std::ostringstream name;
        name << "body_" << body_id;
        rb->set_name(name.str());
        bodies.push_back(rb);
    }
    model->update();

    // The FRET network, built after the bodies exist so the volumes attach to
    // them.
    IMP_NEW(ProbeNetworkRestraint, network,
            (root, fps_path, "ProbeNetworkRestraint%1%", score_set));

    // A fixed (`XYZ`) position is a coordinate this restraint owns rather than
    // an atom, so nothing carries it when a body moves unless it is made a
    // member of one. Its `body_id` says which. (An `ATOM` position needs none
    // of this: it *is* an atom, and is already a member.)
    {
        const std::map<std::string, IMP::ParticleIndex> points =
                network->get_point_positions();
        const std::vector<std::string> atom_points =
                network->get_atom_position_names();
        for (std::map<std::string, IMP::ParticleIndex>::const_iterator it =
                     points.begin();
             it != points.end(); ++it) {
            if (std::find(atom_points.begin(), atom_points.end(), it->first) !=
                atom_points.end()) {
                continue;
            }
            int body_id = 0;
            if (positions.contains(it->first) &&
                positions[it->first].is_object()) {
                body_id = positions[it->first].value("body_id", 0);
            }
            for (std::size_t k = 0; k < bodies.size() &&
                                    k < body_of_pdb.size(); ++k) {
                if (body_of_pdb[k] != body_id) continue;
                IMP::core::RigidBody(bodies[k]).add_member(
                        model->get_particle(it->second));
                break;
            }
        }
        if (!points.empty()) model->update();
    }

    IMP_NEW(IMP::RestraintSet, restraints, (model, "docking"));
    const IMP::bff::AVs used_avs = network->get_used_avs();
    if (mean_position_restraint) {
        // The volumes are sampled once and carried by the rigid body of the
        // atom they hang off, so a move takes them with it. The position an
        // AV is *added at* is the one the body will carry it by, which is why
        // it is resampled first.
        for (std::size_t i = 0; i < used_avs.size(); ++i) {
            IMP::bff::AV av = used_avs[i];
            av.resample();
            IMP::Particle* source = av.get_source();
            if (IMP::core::RigidBodyMember::get_is_setup(source)) {
                IMP::core::RigidBodyMember(source).get_rigid_body().add_member(
                        av.get_particle());
            }
        }
        // `scored_measurements` rather than `get_used_distances`: a distance
        // with a point end is an `Rmp` restraint, and the converter would
        // otherwise apply an <R_DA> transfer to a separation that has no
        // distribution behind it.
        const std::map<std::string, AVPairDistanceMeasurement> used =
                scored_measurements(network);
        for (std::map<std::string,
                      AVPairDistanceMeasurement>::const_iterator it =
                     used.begin();
             it != used.end(); ++it) {
            const IMP::ParticleIndex p1 =
                    network->get_position_particle_index(it->second.position_1);
            const IMP::ParticleIndex p2 =
                    network->get_position_particle_index(it->second.position_2);
            if (p1 == IMP::ParticleIndex() || p2 == IMP::ParticleIndex()) {
                continue;
            }
            restraints->add_restraint(new AVMeanDistanceRestraint(
                    model, p1, p2, it->second, sigma_da, 1.0, max_force));
        }
    } else {
        restraints->add_restraint(network);
    }

    // A volume occupies space: give each one a radius and a mass so the clash
    // term sees it.
    for (std::size_t i = 0; i < used_avs.size(); ++i) {
        IMP::bff::AV av = used_avs[i];
        const IMP::algebra::Vector3D radii = av.get_radii();
        double r_mean = 0.0;
        for (unsigned int k = 0; k < 3; ++k) {
            if (radii[k] > r_mean) r_mean = radii[k];
        }
        IMP::core::XYZR::setup_particle(av.get_particle()).set_radius(r_mean);
        IMP::atom::Mass::setup_particle(av.get_particle(), 0.1)
                .set_mass(r_mean * 2.0);
    }

    // FPS's bond rule, before the clash term is built: a distance between two
    // plain atoms is a covalent tie, and its anchors must not repel.
    set_bond_anchor_radii(network);

    // Excluded volume over every structured leaf; IMP is rigid-body aware, so
    // pairs inside one body are skipped. `clash_container(coarse = false)` is
    // exactly that leaf list -- the scoring door has never used the coarse
    // beads -- and going through it is what gives this door the radii source
    // and scale without a second copy of the rule.
    IMP::Pointer<IMP::container::ListSingletonContainer> leaves =
            clash_container(model, root, false, 2.5, clash_radii,
                            clash_radii_scale);
    // Again, and for the reason `dock_minimize` repeats it: a substituted
    // radii set puts a full van der Waals sphere on a bond anchor, and the
    // 0.4 A written above landed on the *atom*. The second call reaches the
    // shadow sphere the clash term will actually measure.
    if (clash_radii != AV_RADII_IMP || clash_radii_scale != 1.0) {
        set_bond_anchor_radii(network);
    }
    // IMP's soft sphere scores 0.5*k*overlap^2, so FPS's k = 2/ClashTolerance^2
    // makes an overlap of one tolerance cost exactly one chi-square unit --
    // which is what `ClashTolerance` *means* (`SpringEngine.cs:147`).
    const double k_clash = clash_tolerance > 0.0
                                   ? 2.0 / (clash_tolerance * clash_tolerance)
                                   : 1.0;
    IMP_NEW(IMP::core::ExcludedVolumeRestraint, clash,
            (leaves.get(), k_clash, 10.0));
    clash->set_name("excluded_volume");
    clash->set_weight(ev_weight);
    restraints->add_restraint(clash);

    IMP::RestraintsTemp all;
    all.push_back(restraints);
    IMP_NEW(IMP::core::RestraintsScoringFunction, sf, (all));

    return DockingAssembly(model, root, network, restraints, sf, bodies,
                           body_of_pdb, mean_position_restraint, sigma_da,
                           fps_path, max_force);
}

namespace {

//! The excluded-volume term of an assembly, by name.
/*! FPS keeps `Eclash` beside `E` rather than inside it and writes the two as
    separate columns; this module's score is the objective the minimiser
    descended, clash included, so the clash term has to be recoverable from a
    scored assembly rather than derived from it. */
double clash_energy_of(const DockingAssembly& assembly) {
    IMP::RestraintSet* set = assembly.get_restraints();
    if (set == NULL) return 0.0;
    for (unsigned int i = 0; i < set->get_number_of_restraints(); ++i) {
        IMP::Restraint* r = set->get_restraint(i);
        if (r != NULL && r->get_name() == "excluded_volume") {
            return r->evaluate(false);
        }
    }
    return 0.0;
}

//! Fill in `e_bond` / `n_bonds` from a table that already knows its bonds.
/*! FPS accumulates `Ebond` inside the loop that accumulates the total
    (`SpringEngine.cs:445`), so it is a **subset** of the score. Recomputing it
    from the same pairs, with the same cap, keeps that true here. */
void set_bond_totals(DockingResult& out, double max_force) {
    out.n_bonds = 0;
    out.e_bond = 0.0;
    for (std::size_t i = 0; i < out.pairs.size(); ++i) {
        const PairDistance& p = out.pairs[i];
        if (!p.is_bond) continue;
        ++out.n_bonds;
        if (!std::isfinite(p.distance_model)) continue;
        // Halved, because that is the unit every restraint here scores in.
        out.e_bond += 0.5 * chi2_score_capped(p.distance_model, p.distance_exp,
                                              p.error_neg, p.error_pos,
                                              max_force);
    }
}

}  // namespace

DockingResult score_assembly(const DockingAssembly& assembly,
                             const std::string& output_csv) {
    const double total = assembly.evaluate();
    const std::vector<PairDistance> pairs = collect_pair_distances(
            assembly.get_network(), assembly.get_mean_position(),
            assembly.get_sigma_da());
    if (!output_csv.empty()) write_score_csv(output_csv, total, pairs);

    DockingResult out;
    out.score = total;
    out.n_avs = static_cast<int>(assembly.get_network()->get_used_avs().size());
    out.n_distances = static_cast<int>(pairs.size());
    out.pairs = PairDistances(pairs.begin(), pairs.end());
    out.score_csv = output_csv;
    // The excluded-volume term by name, so a caller can report FPS's two
    // columns (`chi2` = total - clash, `chi2_clash`) without re-deriving it.
    out.e_clash = clash_energy_of(assembly);
    set_bond_totals(out, assembly.get_max_force());
    return out;
}

std::string capture_poses(const DockingAssembly& assembly) {
    nlohmann::json out = nlohmann::json::array();
    const IMP::core::RigidBodies bodies = assembly.get_rigid_bodies();
    const std::vector<int> body_of_pdb = assembly.get_body_of_pdb();
    for (std::size_t i = 0; i < bodies.size(); ++i) {
        const IMP::algebra::Transformation3D t =
                bodies[i].get_reference_frame().get_transformation_to();
        const IMP::algebra::Vector3D v = t.get_translation();
        const IMP::algebra::VectorD<4> q = t.get_rotation().get_quaternion();
        nlohmann::json pose;
        pose["body_id"] = i < body_of_pdb.size() ? body_of_pdb[i]
                                                 : static_cast<int>(i);
        pose["t"] = {v[0], v[1], v[2]};
        pose["q"] = {q[0], q[1], q[2], q[3]};
        out.push_back(pose);
    }
    return out.dump();
}

void apply_poses(const DockingAssembly& assembly,
                 const std::string& poses_json) {
    const nlohmann::json poses = nlohmann::json::parse(poses_json, NULL, false);
    if (poses.is_discarded() || !poses.is_array()) {
        IMP_THROW("poses must be a JSON array of {body_id, t, q}",
                  ValueException);
    }
    IMP::core::RigidBodies bodies = assembly.get_rigid_bodies();
    const std::vector<int> body_of_pdb = assembly.get_body_of_pdb();
    for (nlohmann::json::const_iterator it = poses.begin(); it != poses.end();
         ++it) {
        const nlohmann::json& pose = *it;
        if (!pose.is_object() || !pose.contains("t") || !pose.contains("q")) {
            continue;
        }
        const int body_id = pose.value("body_id", -1);
        std::size_t index = bodies.size();
        for (std::size_t i = 0; i < body_of_pdb.size() && i < bodies.size();
             ++i) {
            if (body_of_pdb[i] == body_id) {
                index = i;
                break;
            }
        }
        if (index >= bodies.size()) continue;
        const std::vector<double> t = pose["t"].get<std::vector<double> >();
        const std::vector<double> q = pose["q"].get<std::vector<double> >();
        if (t.size() != 3 || q.size() != 4) continue;
        const IMP::algebra::Rotation3D rotation(q[0], q[1], q[2], q[3]);
        const IMP::algebra::Transformation3D transformation(
                rotation, IMP::algebra::Vector3D(t[0], t[1], t[2]));
        bodies[index].set_reference_frame(
                IMP::algebra::ReferenceFrame3D(transformation));
    }
    assembly.get_model()->update();
}

DockingResult score_structures(const std::vector<std::string>& pdb_paths,
                               const std::string& fps_json_path,
                               const std::string& score_set,
                               bool mean_position_restraint, double sigma_da,
                               const std::string& output_csv,
                               const std::string& clash_radii_source,
                               double clash_radii_scale) {
    const DockingAssembly assembly = create_docking_assembly(
            pdb_paths, fps_json_path, score_set, mean_position_restraint, 1.0,
            sigma_da, 0.0, 0.0, clash_radii_source, clash_radii_scale);
    return score_assembly(assembly, output_csv);
}

// --------------------------------------------------------------------------
// The minimiser
// --------------------------------------------------------------------------

ScoreTrace::ScoreTrace(IMP::Model* m, IMP::ScoringFunction* scoring_function,
                       const std::string& path)
    : IMP::OptimizerState(m, "ScoreTrace%1%"),
      scoring_function_(scoring_function), path_(path), step_(0) {
    std::ofstream out(path_.c_str());
    if (out) out << "frame,score\n";
}

void ScoreTrace::do_update(unsigned int) {
    std::ofstream out(path_.c_str(), std::ios::app);
    if (!out) return;
    out << step_ << "," << scoring_function_->evaluate(false) << "\n";
    ++step_;
}

namespace {

//! One shadow sphere per atom, carrying the radius the clash term should use.
/*! Returns \p atoms unchanged for the default choice, so the ordinary path
    allocates nothing and scores exactly what it scored before this existed.

    A shadow is a bare `XYZR` particle at the atom's coordinate, added to the
    atom's rigid body so a pose change carries it and so the excluded volume's
    gradient reaches the body. It is **not** put in the hierarchy: everything
    that walks the structure -- `assembly_atoms`, `capture_poses`, the PDB
    writers -- goes through `IMP::atom::get_leaves`, and a shadow appearing
    there would silently change an RMSD.

    An atom with no radius at all is passed through as itself; there is nothing
    to substitute and dropping it would quietly shrink the clash term. */
IMP::ParticleIndexes docking_clash_spheres(IMP::Model* m,
                                           const IMP::ParticleIndexes& atoms,
                                           AVRadiiSource radii_source,
                                           double radii_scale) {
    if (radii_source == AV_RADII_IMP && radii_scale == 1.0) return atoms;
    IMP::ParticleIndexes out;
    out.reserve(atoms.size());
    int made = 0;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        const IMP::ParticleIndex pi = atoms[i];
        if (!IMP::core::XYZR::get_is_setup(m, pi)) {
            out.push_back(pi);
            continue;
        }
        // A shadow that already exists is reused: `create_docking_assembly`
        // builds the container once, but a caller is free to ask twice and a
        // second set of spheres would double every overlap.
        if (m->get_has_attribute(docking_clash_shadow_key(), pi)) {
            out.push_back(m->get_attribute(docking_clash_shadow_key(), pi));
            continue;
        }
        IMP::core::XYZR atom(m, pi);
        // The *source* is a van der Waals table keyed by atom name, so it can
        // only answer for an atom. Anything else in the clash set carries a
        // radius that is not a van der Waals radius at all (a dye sphere, a
        // coarse bead) and keeps its own; the scale still applies, because a
        // scale is a statement about the spheres and not about the table.
        const bool is_atom = IMP::atom::Atom::get_is_setup(m, pi);
        const double r = (radii_source == AV_RADII_OLGA && is_atom)
                                 ? olga_vdw_particle_radius(m->get_particle(pi))
                                 : atom.get_radius();
        IMP::Particle* s = new IMP::Particle(
                m, "clashsphere_" + m->get_particle(pi)->get_name());
        IMP::core::XYZR::setup_particle(
                s, IMP::algebra::Sphere3D(atom.get_coordinates(),
                                          r * radii_scale));
        if (IMP::core::RigidBodyMember::get_is_setup(m, pi)) {
            IMP::core::RigidBodyMember(m, pi).get_rigid_body().add_member(s);
        }
        m->add_attribute(docking_clash_shadow_key(), pi, s->get_index());
        out.push_back(s->get_index());
        ++made;
    }
    if (made > 0) {
        /* IMP caches every `ModelObject`'s input and output lists, and a rigid
           body's position constraint reports its **members** as its outputs.
           Adding members after the model has computed dependencies once leaves
           that list stale, and a scoring function built afterwards then decides
           the constraint is *not required* -- so `update_members` never runs,
           `pull_back_members_adjoints` never runs, and no derivative reaches
           the body. `RigidBody::add_member` calls `Model::clear_caches()`,
           which does not cover this.

           Measured on HIV-RT before this loop existed: with the excluded volume
           over the shadows, `ScoringFunction::get_required_score_states()`
           dropped from 3 to 1, body 1's gradient came out exactly (0, 0, 0),
           and `dock_minimize` walked the DNA 147 A to the origin while
           reporting a *better* score (31.69) because the frozen proxies went
           with it. The score was right and the gradient was gone, which is the
           worst shape a defect can have. */
        for (unsigned int i = 0; i < m->get_number_of_score_states(); ++i) {
            m->get_score_state(i)->set_has_dependencies(false);
        }
    }
    m->update();
    return out;
}

}  // namespace

IMP::container::ListSingletonContainer* clash_container(
        IMP::Model* m, IMP::atom::Hierarchy root, bool coarse,
        double bead_radius, AVRadiiSource radii_source, double radii_scale) {
    if (radii_scale <= 0.0) {
        IMP_THROW("clash radii scale must be positive, got " << radii_scale,
                  ValueException);
    }
    const bool substituted =
            radii_source != AV_RADII_IMP || radii_scale != 1.0;
    if (coarse && substituted) {
        // A coarse bead is one sphere standing for a whole residue. Olga's
        // table is keyed by atom name and has no entry for "residue", and a
        // scale on the bead is `bead_radius` under a second name -- so both
        // would be answering a question the representation cannot be asked.
        IMP_THROW("the clash radii source and scale describe atoms, and the "
                  "coarse clash representation is one bead per residue, not "
                  "an atom. Turn the coarse representation off "
                  "(DockingParameters::coarse_clash = false, or "
                  "--no-coarse-clash) to use them, or set bead_radius to "
                  "resize the beads.",
                  ValueException);
    }
    static const char* kBackbone[] = {"CA", "P", "C1'", "C4'"};
    const IMP::atom::Hierarchies leaves = IMP::atom::get_leaves(root);
    IMP::ParticleIndexes chosen;
    IMP::ParticleIndexes all;
    for (std::size_t i = 0; i < leaves.size(); ++i) {
        IMP::Particle* p = leaves[i].get_particle();
        all.push_back(p->get_index());
        if (!coarse || !IMP::atom::Atom::get_is_setup(p)) continue;
        const std::string name = internal::trimmed(
                IMP::atom::Atom(p).get_atom_type().get_string());
        for (int k = 0; k < 4; ++k) {
            if (name != kBackbone[k]) continue;
            if (IMP::core::XYZR::get_is_setup(p)) {
                IMP::core::XYZR(p).set_radius(bead_radius);
            }
            chosen.push_back(p->get_index());
            break;
        }
    }
    // A structure with no recognised backbone atom (a ligand, a coarse model)
    // is scored atomistically rather than not at all.
    const IMP::ParticleIndexes& atoms = chosen.empty() ? all : chosen;
    return new IMP::container::ListSingletonContainer(
            m, docking_clash_spheres(m, atoms, radii_source, radii_scale));
}

namespace {

//! Move each mobile body by a random translation of at most \p amplitude.
/*! `IMP.pmi.tools.shuffle_configuration` does this too, but PMI is not one of
    this module's dependencies and a shuffle that silently does nothing when it
    is absent is worse than none. A random direction and a random distance,
    applied to the bodies that are free to move. */
void shuffle_bodies(const IMP::core::RigidBodies& bodies, double amplitude) {
    if (amplitude <= 0.0) return;
    for (std::size_t i = 0; i < bodies.size(); ++i) {
        IMP::core::RigidBody rb = bodies[i];
        if (!rb.get_coordinates_are_optimized()) continue;
        const IMP::algebra::Vector3D step =
                IMP::algebra::get_random_vector_in(
                        IMP::algebra::Sphere3D(IMP::algebra::Vector3D(0, 0, 0),
                                               amplitude));
        const IMP::algebra::ReferenceFrame3D frame = rb.get_reference_frame();
        const IMP::algebra::Transformation3D t =
                frame.get_transformation_to();
        rb.set_reference_frame(IMP::algebra::ReferenceFrame3D(
                IMP::algebra::Transformation3D(t.get_rotation(),
                                               t.get_translation() + step)));
    }
}

//! Conjugate gradients in chunks, so a stop is answered promptly.
/*! Chunking costs the descent almost nothing and is the only reliable place
    to ask: raising out of an optimiser-state callback is not.
    \return true when the run was stopped early */
bool optimize_in_chunks(IMP::core::ConjugateGradients* cg, int n,
                        DockingStop* stop, bool* converged = NULL) {
    const int chunk = std::min(100, std::max(20, n / 20));
    int done = 0;
    // FPS reports `niter < MaxIterations`, which is unconditionally true under
    // SelectedThenAll (its D5) and says nothing about the gradient. Measure it
    // instead: the descent has converged when a whole chunk moves the score by
    // less than 1e-6 of it. One extra evaluate per chunk (~20 per run) against
    // hundreds of gradient evaluations inside the chunk.
    double last = std::numeric_limits<double>::infinity();
    bool settled = false;
    while (done < n) {
        if (stop != NULL && stop->should_stop()) {
            if (converged != NULL) *converged = false;
            return true;
        }
        double score = last;
        try {
            score = cg->optimize(std::min(chunk, n - done));
        } catch (const IMP::ModelException&) {
            // A near-singular gradient -- a severe clash right after a
            // shuffle, say -- ends the descent with the best pose reached
            // rather than throwing the run away.
            break;
        }
        // Observed, not acted on: the budget is still spent in full, so no
        // score this module pins can move because the flag was added.
        settled = std::isfinite(score) && std::isfinite(last) &&
                  std::fabs(score - last) <
                          1e-6 * std::max(1.0, std::fabs(score));
        last = score;
        done += chunk;
    }
    if (converged != NULL) *converged = settled;
    return false;
}

//! `OptimizeSelected` as a lower-case token; anything unknown is "selected".
std::string optimize_mode(const std::string& text) {
    std::string t;
    for (std::size_t i = 0; i < text.size(); ++i) {
        const char c = text[i];
        if (c == ' ' || c == '_' || c == '-') continue;
        t += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    if (t == "all") return "all";
    if (t == "selectedthenall") return "selectedthenall";
    return "selected";
}

//! The distance names one named score set holds; empty means "no such set".
std::set<std::string> score_set_distance_names(const std::string& fps_json,
                                               const std::string& score_set) {
    std::set<std::string> out;
    if (score_set.empty() || fps_json.empty()) return out;
    try {
        const FPSDocument document = read_fps_json(fps_json);
        const nlohmann::json sets =
                nlohmann::json::parse(document.score_sets, NULL, false);
        if (!sets.is_object() || !sets.contains(score_set)) return out;
        const nlohmann::json& listed = sets[score_set]["distances"];
        if (!listed.is_array()) return out;
        for (nlohmann::json::const_iterator it = listed.begin();
             it != listed.end(); ++it) {
            if (it->is_string()) out.insert(it->get<std::string>());
        }
    } catch (const std::exception&) {
        // A file whose score sets cannot be read scores every distance, which
        // is what "no selection" means anyway.
    }
    return out;
}

}  // namespace

DockingResult dock_minimize(const std::vector<std::string>& pdb_paths,
                            const std::string& fps_json_path,
                            const std::string& output_dir,
                            const DockingParameters& params, DockingStop* stop,
                            const std::string& initial_poses) {
    // Before the directory, the files or the volumes: a coarse clash term
    // cannot honour a per-atom radii choice, and finding that out after the
    // AVs have been computed is the same error message half an hour later.
    if (params.coarse_clash &&
        (av_radii_source_from_string(params.clash_radii_source) !=
                 AV_RADII_IMP ||
         params.clash_radii_scale != 1.0)) {
        IMP_THROW("DockingParameters: clash_radii_source / clash_radii_scale "
                  "describe atoms, and coarse_clash replaces them with one "
                  "bead per residue. Set coarse_clash = false "
                  "(--no-coarse-clash) to use them.",
                  ValueException);
    }
    internal::make_directory(output_dir);

    // FPS's `OptimizeSelected`, over *distances* (SpringEngine.cs:428). Under
    // "All" and "SelectedThenAll" the run has to see distances the score set
    // leaves out, so the network is built without the score set and the set is
    // used to say which of them are the selected ones.
    const std::string mode = optimize_mode(params.optimize_selected);
    const bool needs_all = (mode != "selected");

    // `mean_position_restraint=false`: the assembly must not attach the dye
    // particles as body members here, because the proxies below do that job in
    // the one way that keeps the gradients (see the header).
    const DockingAssembly assembly = create_docking_assembly(
            pdb_paths, fps_json_path,
            needs_all ? std::string() : params.score_set, false,
            params.ev_weight, params.sigma_da, params.clash_tolerance,
            params.max_force, params.clash_radii_source,
            params.clash_radii_scale);
    IMP::Model* model = assembly.get_model();
    IMP::atom::Hierarchy root = assembly.get_root();

    // Resume before the proxies are built, so they are placed in the resumed
    // body frames.
    if (!initial_poses.empty()) apply_poses(assembly, initial_poses);

    const IMP::core::RigidBodies bodies = assembly.get_rigid_bodies();
    const std::vector<int> body_of_pdb = assembly.get_body_of_pdb();
    int n_mobile = 0;
    for (std::size_t i = 0; i < bodies.size(); ++i) {
        const int body_id = i < body_of_pdb.size() ? body_of_pdb[i]
                                                   : static_cast<int>(i);
        const bool mobile = body_id != params.fixed_body;
        IMP::core::RigidBody(bodies[i]).set_coordinates_are_optimized(mobile);
        if (mobile) ++n_mobile;
    }
    if (n_mobile == 0 && !bodies.empty()) {
        // One body: let it move, or there is nothing to dock.
        IMP::core::RigidBody(bodies[0]).set_coordinates_are_optimized(true);
        n_mobile = 1;
    }
    if (initial_poses.empty()) {
        shuffle_bodies(bodies, params.shuffle_max_translation);
    }

    ProbeNetworkRestraint* network = assembly.get_network();
    const IMP::bff::AVs used_avs = network->get_used_avs();
    // The point rule (an `XYZ`/`ATOM` end forces `Rmp`) applied once, here, so
    // the springs and the table below cannot disagree about it.
    const std::map<std::string, AVPairDistanceMeasurement> used =
            scored_measurements(network);

    // Point-member proxies at the mean positions: they ride rigidly with the
    // body of the atom the dye hangs off and feed gradients back to it.
    std::vector<std::string> proxy_names;
    IMP::ParticleIndexes proxy_particles;
    std::map<std::string, IMP::ParticleIndex> proxy_of;
    for (std::size_t i = 0; i < used_avs.size(); ++i) {
        IMP::bff::AV av = used_avs[i];
        av.resample();
        const std::string name = av.get_particle()->get_name();
        IMP::Particle* q = new IMP::Particle(model, "avproxy_" + name);
        IMP::core::XYZ::setup_particle(
                q, IMP::core::XYZ(av).get_coordinates());
        IMP::Particle* source = av.get_source();
        if (IMP::core::RigidBodyMember::get_is_setup(source)) {
            IMP::core::RigidBodyMember(source).get_rigid_body().add_member(q);
        }
        proxy_names.push_back(name);
        proxy_particles.push_back(q->get_index());
        proxy_of.insert(std::make_pair(name, q->get_index()));
    }
    // A point position needs no proxy: it is already a plain XYZ member of its
    // rigid body (an atom, or the fixed point `create_docking_assembly`
    // attached), so it moves with the body and propagates gradients the same
    // classic way a proxy does.
    {
        const std::map<std::string, IMP::ParticleIndex> points =
                network->get_point_positions();
        for (std::map<std::string, IMP::ParticleIndex>::const_iterator it =
                     points.begin();
             it != points.end(); ++it) {
            proxy_names.push_back(it->first);
            proxy_particles.push_back(it->second);
            proxy_of.insert(*it);
        }
    }
    model->update();

    // Two sets: FPS gates *distances* on their selection flag, and this file's
    // score set is that flag. Clashes are in neither -- they are evaluated
    // globally in every phase (SpringEngine.cs:326), which is why `evr` is
    // added to both scoring functions rather than to a set.
    const std::set<std::string> selected = score_set_distance_names(
            assembly.get_fps_json_path(), params.score_set);
    IMP_NEW(IMP::RestraintSet, springs, (model, "fret_min_selected"));
    IMP_NEW(IMP::RestraintSet, springs_rest, (model, "fret_min_deselected"));
    for (std::map<std::string, AVPairDistanceMeasurement>::const_iterator it =
                 used.begin();
         it != used.end(); ++it) {
        std::map<std::string, IMP::ParticleIndex>::const_iterator q1 =
                proxy_of.find(it->second.position_1);
        std::map<std::string, IMP::ParticleIndex>::const_iterator q2 =
                proxy_of.find(it->second.position_2);
        if (q1 == proxy_of.end() || q2 == proxy_of.end()) continue;
        IMP::Restraint* r = new AVMeanDistanceRestraint(
                model, q1->second, q2->second, it->second, params.sigma_da, 1.0,
                params.max_force);
        const bool is_selected = selected.empty() || selected.count(it->first);
        (is_selected ? springs : springs_rest)->add_restraint(r);
    }

    IMP::Pointer<IMP::container::ListSingletonContainer> clash =
            clash_container(model, root, params.coarse_clash, 2.5,
                            av_radii_source_from_string(
                                    params.clash_radii_source),
                            params.clash_radii_scale);
    // After clash_container, which writes a bead radius on every backbone atom
    // it picks -- and would otherwise undo the 0.4 A on a bond anchor that
    // happens to be a CA, and which under a substituted radii set puts a full
    // van der Waals sphere on the anchor's shadow.
    const int n_bond_anchors = set_bond_anchor_radii(network);
    const double k_clash =
            params.clash_tolerance > 0.0
                    ? 2.0 / (params.clash_tolerance * params.clash_tolerance)
                    : 1.0;
    IMP_NEW(IMP::core::ExcludedVolumeRestraint, evr,
            (clash.get(), k_clash, 10.0));
    evr->set_weight(params.ev_weight);
    IMP::RestraintsTemp phase1;
    phase1.push_back(springs);
    if (mode == "all") phase1.push_back(springs_rest);
    phase1.push_back(evr);
    IMP_NEW(IMP::core::RestraintsScoringFunction, sf, (phase1));
    IMP::RestraintsTemp everything;
    everything.push_back(springs);
    everything.push_back(springs_rest);
    everything.push_back(evr);
    IMP_NEW(IMP::core::RestraintsScoringFunction, sf_all, (everything));

    const std::string convergence_csv = output_dir + "/convergence.csv";
    IMP_NEW(IMP::core::ConjugateGradients, cg, (model));
    cg->set_scoring_function(sf);
    // Under SelectedThenAll each phase gets half the budget, as FPS caps them
    // (SpringEngine.cs:255).
    const int n_iter = std::max(1, params.n_frames);
    const int n_phase1 =
            mode == "selectedthenall" ? std::max(1, n_iter / 2) : n_iter;
    // The trace follows the objective the run is finally judged by, not the
    // phase-1 one: under SelectedThenAll those differ, and a convergence plot
    // whose last row is not the reported score is a plot that invites the
    // wrong conclusion.
    IMP::core::RestraintsScoringFunction* traced =
            mode == "selected" ? sf.get() : sf_all.get();
    IMP_NEW(ScoreTrace, trace, (model, traced, convergence_csv));
    trace->set_period(std::max(1, n_phase1 / 100));
    cg->add_optimizer_state(trace);

    bool converged = false;
    bool stopped = optimize_in_chunks(cg, n_phase1, stop, &converged);
    // Phase two: the same optimiser, everything scoring.
    if (mode == "selectedthenall" && !stopped) {
        cg->set_scoring_function(sf_all);
        stopped = optimize_in_chunks(cg, n_iter - n_phase1, stop, &converged);
    }
    // What the run is finally judged by. For "Selected" the deselected
    // distances are not part of it, exactly as FPS leaves them out.
    IMP::core::RestraintsScoringFunction* final_sf =
            mode == "selected" ? sf.get() : sf_all.get();
    cg->set_scoring_function(final_sf);

    // FPS-style refinement: with the partner docked alongside, each volume is
    // resampled so inter-body occlusion moves its mean position, the proxy is
    // put back on that mean *in the body frame*, and the descent runs again.
    // The volumes are recomputed once per cycle, not once per step.
    for (int cycle = 0; cycle < params.refine_av_cycles && !stopped; ++cycle) {
        for (std::size_t i = 0; i < used_avs.size(); ++i) {
            IMP::bff::AV av = used_avs[i];
            av.resample();
            std::map<std::string, IMP::ParticleIndex>::const_iterator q =
                    proxy_of.find(av.get_particle()->get_name());
            if (q == proxy_of.end()) continue;
            IMP::Particle* proxy = model->get_particle(q->second);
            if (!IMP::core::RigidBodyMember::get_is_setup(proxy)) continue;
            IMP::core::RigidBodyMember member(proxy);
            const IMP::algebra::Vector3D local =
                    member.get_rigid_body()
                            .get_reference_frame()
                            .get_transformation_to()
                            .get_inverse()
                            .get_transformed(
                                    IMP::core::XYZ(av).get_coordinates());
            member.set_internal_coordinates(local);
        }
        model->update();
        stopped = optimize_in_chunks(cg, n_iter, stop, &converged);
    }

    const double total = final_sf->evaluate(false);
    // FPS keeps `Eclash` beside `E` rather than inside it; here `total` is the
    // objective the descent actually followed, clash included, so the clash
    // term is recorded separately and `total - e_clash` is FPS's `E`.
    const double e_clash = evr->evaluate(false);
    // The proxies are what moved with the bodies; the volumes did not, so the
    // table is read at the proxies' coordinates.
    std::vector<PairDistance> pairs = pair_distances_at_positions(
            used, model, proxy_names, proxy_particles, params.sigma_da);
    for (std::size_t i = 0; i < pairs.size(); ++i) {
        pairs[i].is_bond = network->get_is_bond(pairs[i].name);
    }

    const std::string out_pdb = output_dir + "/docked.pdb";
    IMP::atom::write_pdb(root, out_pdb);
    const std::string score_csv = output_dir + "/scores.csv";
    write_score_csv(score_csv, total, pairs);

    nlohmann::json extra;
    extra["method"] = "minimize";
    extra["iterations"] = n_iter;
    extra["n_mobile"] = n_mobile;
    extra["convergence_csv"] = convergence_csv;
    extra["coarse_clash"] = params.coarse_clash;
    extra["refine_av_cycles"] = params.refine_av_cycles;
    extra["stopped"] = stopped;
    extra["converged"] = converged;
    extra["e_clash"] = e_clash;
    // The protocol knobs, in the result: a docked pose that does not say which
    // of these it was produced under cannot be compared with another one.
    extra["optimize_selected"] = mode;
    extra["n_selected_distances"] =
            static_cast<int>(springs->get_number_of_restraints());
    extra["n_deselected_distances"] =
            static_cast<int>(springs_rest->get_number_of_restraints());
    extra["max_force"] = params.max_force;
    extra["clash_tolerance"] = params.clash_tolerance;
    extra["k_clash"] = k_clash;
    // What the spheres were sized by. `k_clash` alone does not say it, and the
    // same tolerance over two radii sets is two different terms.
    extra["clash_radii_source"] = params.clash_radii_source;
    extra["clash_radii_scale"] = params.clash_radii_scale;
    extra["n_bond_anchors"] = n_bond_anchors;
    extra["n_points"] =
            static_cast<int>(network->get_point_positions().size());

    DockingResult out;
    out.score = total;
    out.n_avs = static_cast<int>(used_avs.size());
    out.n_distances = static_cast<int>(pairs.size());
    out.pairs = PairDistances(pairs.begin(), pairs.end());
    out.output_dir = output_dir;
    out.best_pdbs.push_back(out_pdb);
    out.score_csv = score_csv;
    out.extra = extra.dump();
    out.poses = capture_poses(assembly);
    out.e_clash = e_clash;
    out.converged = converged && !stopped;
    set_bond_totals(out, params.max_force);
    return out;
}

DockingResult refine_docking(const std::vector<std::string>& pdb_paths,
                     const std::string& fps_json_path,
                     const std::string& output_dir,
                     const std::string& score_set, int steps,
                     double ev_weight) {
    internal::make_directory(output_dir);
    const DockingAssembly assembly = create_docking_assembly(
            pdb_paths, fps_json_path, score_set, true, ev_weight);
    const IMP::core::RigidBodies bodies = assembly.get_rigid_bodies();
    for (std::size_t i = 0; i < bodies.size(); ++i) {
        IMP::core::RigidBody(bodies[i]).set_coordinates_are_optimized(true);
    }

    IMP_NEW(IMP::core::ConjugateGradients, cg, (assembly.get_model()));
    cg->set_scoring_function(assembly.get_scoring_function());
    optimize_in_chunks(cg, std::max(1, steps), NULL);

    const std::string out_pdb = output_dir + "/refined.pdb";
    IMP::atom::write_pdb(assembly.get_root(), out_pdb);
    DockingResult out = score_assembly(assembly, output_dir + "/scores.csv");
    out.output_dir = output_dir;
    out.best_pdbs.push_back(out_pdb);
    out.poses = capture_poses(assembly);
    return out;
}

namespace {

//! Every atom of a structure, keyed the two ways a reference atom names one.
/*! Built once and reused: a frame of twenty atoms times a library of four
    hundred structures is eight thousand lookups, and an `IMP::atom::Selection`
    walks the whole hierarchy each time. */
struct AtomIndex {
    std::map<std::string, IMP::ParticleIndex> by_site;
    //! Sites whose chain-agnostic key matched more than one chain.
    std::set<std::string> ambiguous;
    std::map<int, IMP::ParticleIndex> by_serial;

    static std::string key(const std::string& chain, int residue,
                           const std::string& atom) {
        std::ostringstream out;
        out << chain << ":" << residue << ":" << atom;
        return out.str();
    }

    explicit AtomIndex(IMP::atom::Hierarchy root) {
        const IMP::atom::Hierarchies leaves = IMP::atom::get_leaves(root);
        for (std::size_t i = 0; i < leaves.size(); ++i) {
            if (!IMP::atom::Atom::get_is_setup(leaves[i])) continue;
            const IMP::atom::Atom atom(leaves[i]);
            const std::string name = atom.get_atom_type().get_string();
            const IMP::atom::Residue residue =
                    IMP::atom::get_residue(atom, true);
            if (residue == IMP::atom::Residue()) continue;
            const std::string chain = IMP::atom::get_chain_id(residue);
            const int index = residue.get_index();
            by_site[key(chain, index, name)] = leaves[i].get_particle_index();
            // The chain-agnostic key is a convenience for single-chain files
            // and a trap for everything else, so a collision disables it
            // rather than picking whichever chain was read first. A structure
            // whose chains have no id already *is* that key, and marking it
            // ambiguous against itself would make every atom unfindable.
            if (!chain.empty()) {
                const std::string any = key("", index, name);
                if (by_site.count(any)) {
                    ambiguous.insert(any);
                } else {
                    by_site[any] = leaves[i].get_particle_index();
                }
            }
            const int serial = atom.get_input_index();
            if (serial > 0 && !by_serial.count(serial)) {
                by_serial[serial] = leaves[i].get_particle_index();
            }
        }
    }

    //! The atom, or a null index when this structure does not have it.
    IMP::ParticleIndex find(const ReferenceAtom& ref) const {
        if (!ref.atom_name.empty() && ref.residue_seq_number > 0) {
            const std::string k = key(ref.chain_identifier,
                                      ref.residue_seq_number, ref.atom_name);
            if (!ambiguous.count(k)) {
                std::map<std::string, IMP::ParticleIndex>::const_iterator it =
                        by_site.find(k);
                if (it != by_site.end()) return it->second;
            }
        }
        if (ref.atom_serial > 0) {
            std::map<int, IMP::ParticleIndex>::const_iterator it =
                    by_serial.find(ref.atom_serial);
            if (it != by_serial.end()) return it->second;
        }
        return IMP::ParticleIndex();
    }
};

//! The fit, against an index that has already been built.
ReferenceFit fit_against(const ReferenceAtoms& atoms, IMP::Model* m,
                         const AtomIndex& index,
                         const IMP::algebra::Vector3D& point) {
    ReferenceFit out;
    out.coordinates = point;
    IMP::algebra::Vector3Ds source, target;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        const IMP::ParticleIndex pi = index.find(atoms[i]);
        if (pi == IMP::ParticleIndex()) continue;
        source.push_back(atoms[i].get_coordinates());
        target.push_back(IMP::core::XYZ(m, pi).get_coordinates());
    }
    out.n_atoms = static_cast<int>(source.size());
    if (source.size() < 3) return out;

    const IMP::algebra::Transformation3D t =
            IMP::algebra::get_transformation_aligning_first_to_second(source,
                                                                     target);
    double sum = 0.0;
    for (std::size_t i = 0; i < source.size(); ++i) {
        sum += (t.get_transformed(source[i]) - target[i]).get_squared_magnitude();
    }
    out.squared_deviation = sum;
    out.rmsd = std::sqrt(sum / static_cast<double>(source.size()));
    out.coordinates = t.get_transformed(point);
    out.fitted = true;
    return out;
}

//! The `x`/`y`/`z` of a position object, or the origin when it has none.
IMP::algebra::Vector3D fixed_point(const nlohmann::json& position) {
    return IMP::algebra::Vector3D(position.value("x", 0.0),
                                  position.value("y", 0.0),
                                  position.value("z", 0.0));
}

ReferenceAtoms reference_atoms_of(const nlohmann::json& position) {
    ReferenceAtoms out;
    if (!position.is_object() || !position.contains("reference_atoms")) {
        return out;
    }
    const nlohmann::json& listed = position["reference_atoms"];
    if (!listed.is_array()) return out;
    for (nlohmann::json::const_iterator it = listed.begin();
         it != listed.end(); ++it) {
        if (!it->is_object()) continue;
        ReferenceAtom atom;
        atom.chain_identifier = it->value("chain_identifier", std::string());
        atom.residue_seq_number = it->value("residue_seq_number", 0);
        atom.atom_name = it->value("atom_name", std::string());
        atom.atom_serial = it->value("atom_serial", 0);
        atom.x = it->value("x", 0.0);
        atom.y = it->value("y", 0.0);
        atom.z = it->value("z", 0.0);
        out.push_back(atom);
    }
    return out;
}

}  // namespace

ReferenceAtoms fps_reference_atoms(const std::string& position_json) {
    const nlohmann::json position =
            nlohmann::json::parse(position_json, NULL, false);
    if (position.is_discarded()) {
        IMP_THROW("cannot parse the position as JSON", ValueException);
    }
    return reference_atoms_of(position);
}

ReferenceFit fit_reference_atoms(const ReferenceAtoms& atoms,
                                 IMP::atom::Hierarchy structure,
                                 const IMP::algebra::Vector3D& point) {
    const AtomIndex index(structure);
    return fit_against(atoms, structure.get_model(), index, point);
}

ReferenceFits fit_reference_positions(const std::string& positions_json,
                                      IMP::atom::Hierarchy structure) {
    const nlohmann::json positions =
            nlohmann::json::parse(positions_json, NULL, false);
    ReferenceFits out;
    if (positions.is_discarded() || !positions.is_object()) return out;
    // The index is worth building only if something is going to ask it a
    // question, and most files ask none.
    bool any = false;
    for (nlohmann::json::const_iterator it = positions.begin();
         it != positions.end() && !any; ++it) {
        any = !reference_atoms_of(it.value()).empty();
    }
    if (!any) return out;

    const AtomIndex index(structure);
    for (nlohmann::json::const_iterator it = positions.begin();
         it != positions.end(); ++it) {
        const ReferenceAtoms atoms = reference_atoms_of(it.value());
        if (atoms.empty()) continue;
        ReferenceFit fit = fit_against(atoms, structure.get_model(), index,
                                       fixed_point(it.value()));
        fit.position = it.key();
        out.push_back(fit);
    }
    return out;
}

double reference_rmsd(const ReferenceFits& fits) {
    double sum = 0.0;
    int n = 0;
    for (std::size_t i = 0; i < fits.size(); ++i) {
        if (!fits[i].fitted) continue;
        sum += fits[i].squared_deviation;
        n += fits[i].n_atoms;
    }
    if (n == 0) {
        // Nothing declared is 0 (FPS's answer); declared and unfittable is
        // not, and saying 0 there would read as a perfect fit.
        return fits.empty() ? 0.0 : std::numeric_limits<double>::quiet_NaN();
    }
    return std::sqrt(sum / static_cast<double>(n));
}

std::vector<ScreenedStructure> screen_structures(
        const std::vector<std::string>& pdb_paths,
        const std::string& fps_json_path, const std::string& score_set,
        const std::string& output_csv, bool mean_position_restraint) {
    std::vector<std::string> structures;
    for (std::size_t i = 0; i < pdb_paths.size(); ++i) {
        const std::vector<std::string> expanded =
                internal::directory_entries(pdb_paths[i], ".pdb");
        if (expanded.empty()) {
            structures.push_back(pdb_paths[i]);
        } else {
            structures.insert(structures.end(), expanded.begin(),
                              expanded.end());
        }
    }

    // The reference frames are a property of the *file*, not of a structure,
    // so the positions are read once. Only fps.json is read: a legacy C# FPS
    // labelling file carries its frames as `ATOM` lines trailing an LP line
    // (`FilterEngine.ReadRefAtoms`) and `read_old_lps_txt` does not yet keep
    // them, so converting one here would find nothing and cost a pass over
    // the library's PDBs to learn it.
    std::string positions_json = "{}";
    if (internal::ends_with(fps_json_path, ".json")) {
        try {
            const FPSDocument document = read_fps_json(fps_json_path);
            if (!document.positions.empty()) {
                positions_json = document.positions;
            }
        } catch (const std::exception&) {
            // Diagnostics only; create_docking_assembly reports the real
            // failure, per structure and with the structure's name on it.
        }
    }

    std::vector<ScreenedStructure> out;
    for (std::size_t i = 0; i < structures.size(); ++i) {
        ScreenedStructure entry;
        entry.path = structures[i];
        entry.score = std::numeric_limits<double>::quiet_NaN();
        entry.chi2_r = std::numeric_limits<double>::quiet_NaN();
        entry.ref_rmsd = std::numeric_limits<double>::quiet_NaN();
        try {
            std::vector<std::string> one(1, structures[i]);
            const DockingAssembly assembly = create_docking_assembly(
                    one, fps_json_path, score_set, mean_position_restraint);
            entry.score = assembly.evaluate();
            const std::vector<PairDistance> ep_ = collect_pair_distances(assembly.get_network(),
                                                 assembly.get_mean_position(),
                                                 assembly.get_sigma_da());
            entry.pairs = PairDistances(ep_.begin(), ep_.end());

            double chi2 = 0.0;
            int scored = 0;
            for (std::size_t k = 0; k < entry.pairs.size(); ++k) {
                const PairDistance& p = entry.pairs[k];
                if (!std::isfinite(p.distance_model)) {
                    entry.invalid_r++;
                    continue;
                }
                chi2 += p.get_chi2();
                scored++;
                // FPS counts each side against its own error bar, and the
                // counts nest: a 3 sigma outlier is 1, 2 and 3 sigma at once
                // (`FilterEngine.CalculateChi2`).
                const double dr = p.get_residual();
                const double err = dr > 0.0 ? p.error_pos : p.error_neg;
                if (!(err > 0.0)) continue;
                const double z = std::fabs(dr) / err;
                if (z > 1.0) entry.sigma1++;
                if (z > 2.0) entry.sigma2++;
                if (z > 3.0) entry.sigma3++;
            }
            entry.chi2_r = scored > 0
                                   ? chi2 / static_cast<double>(scored)
                                   : std::numeric_limits<double>::quiet_NaN();
            entry.ref_rmsd = reference_rmsd(
                    fit_reference_positions(positions_json,
                                            assembly.get_root()));
        } catch (const IMP::Exception&) {
            // One unscorable structure does not end the library.
        }
        out.push_back(entry);
    }
    // Ascending, with the unscorable last: a NaN sorts nowhere on its own.
    std::stable_sort(out.begin(), out.end(), [](const ScreenedStructure& a,
                                                const ScreenedStructure& b) {
        if (std::isnan(a.score)) return false;
        if (std::isnan(b.score)) return true;
        return a.score < b.score;
    });

    if (!output_csv.empty()) {
        std::ofstream csv(output_csv.c_str());
        if (!csv) IMP_THROW("cannot write " << output_csv, IOException);
        csv << "pdb,score,chi2_r,sigma1,sigma2,sigma3,invalid_r,ref_rmsd\n";
        for (std::size_t i = 0; i < out.size(); ++i) {
            const ScreenedStructure& s = out[i];
            csv << cell(s.path) << "," << number(s.score) << ","
                << number(rounded(s.chi2_r, 6)) << "," << s.sigma1 << ","
                << s.sigma2 << "," << s.sigma3 << "," << s.invalid_r << ","
                << number(rounded(s.ref_rmsd, 4)) << "\n";
        }
    }
    return out;
}

/* ------------------------------------------------------------------------
 * Error estimation -- the parametric bootstrap (PRD-121 G2)
 * ------------------------------------------------------------------------ */

namespace {

//! \f$1/\sqrt{2\pi}\f$, the height of a standard normal at its mode.
const double kInvSqrt2Pi = 0.3989422804014327;

//! Both error bars as non-negative scales; a missing bar is a zero-width side.
void perturbation_scales(double error_neg, double error_pos, double* s_minus,
                         double* s_plus) {
    *s_minus = error_neg > 0.0 ? error_neg : 0.0;
    *s_plus = error_pos > 0.0 ? error_pos : 0.0;
}

}  // namespace

std::vector<double> sample_distance_perturbations(double error_neg,
                                                  double error_pos, int n,
                                                  unsigned int seed,
                                                  FPSPerturbationModel model) {
    std::vector<double> out;
    if (n <= 0) return out;
    double s_minus = 0.0, s_plus = 0.0;
    perturbation_scales(error_neg, error_pos, &s_minus, &s_plus);
    out.reserve(static_cast<std::size_t>(n));

    std::mt19937 rng(seed);
    std::normal_distribution<double> normal(0.0, 1.0);
    if (model == FPS_SIGN_SPLIT_NORMAL) {
        // `ErrorEstimation.cs:56-57`: one standard normal draw, scaled by the
        // error bar the *sign of the draw* selects. Each side therefore gets
        // mass 1/2 whatever the two widths are.
        for (int i = 0; i < n; ++i) {
            const double z = normal(rng);
            out.push_back(z > 0.0 ? z * s_plus : z * s_minus);
        }
        return out;
    }
    // The two-piece normal: pick the side with probability proportional to its
    // own width, then draw a half-normal of that width. That is what makes the
    // density continuous at the mode -- the wider side is visited more often in
    // exactly the proportion that compensates its lower peak.
    const double total = s_minus + s_plus;
    const double p_up = total > 0.0 ? s_plus / total : 0.5;
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    for (int i = 0; i < n; ++i) {
        const double half = std::fabs(normal(rng));
        out.push_back(uniform(rng) < p_up ? half * s_plus : -half * s_minus);
    }
    return out;
}

double perturbation_mean(double error_neg, double error_pos,
                         FPSPerturbationModel model) {
    double s_minus = 0.0, s_plus = 0.0;
    perturbation_scales(error_neg, error_pos, &s_minus, &s_plus);
    // Sign-split: E[X] = (1/2)E[|Z|]s+ - (1/2)E[|Z|]s- with E[|Z|] = 2/sqrt(2pi).
    // Split normal: the side probabilities s+-/(s++s-) turn that into twice as
    // much, because the wide side is chosen more often as well as being wider.
    const double factor = model == FPS_SIGN_SPLIT_NORMAL ? 1.0 : 2.0;
    return factor * (s_plus - s_minus) * kInvSqrt2Pi;
}

double perturbation_sd(double error_neg, double error_pos,
                       FPSPerturbationModel model) {
    double s_minus = 0.0, s_plus = 0.0;
    perturbation_scales(error_neg, error_pos, &s_minus, &s_plus);
    const double mean = perturbation_mean(error_neg, error_pos, model);
    double second = 0.0;
    if (model == FPS_SIGN_SPLIT_NORMAL) {
        second = 0.5 * (s_plus * s_plus + s_minus * s_minus);
    } else {
        // E[X^2] = (s+^3 + s-^3)/(s+ + s-) = s+^2 - s+ s- + s-^2.
        second = s_plus * s_plus - s_plus * s_minus + s_minus * s_minus;
    }
    const double variance = second - mean * mean;
    return variance > 0.0 ? std::sqrt(variance) : 0.0;
}

double perturbation_upper_mass(double error_neg, double error_pos,
                               FPSPerturbationModel model) {
    double s_minus = 0.0, s_plus = 0.0;
    perturbation_scales(error_neg, error_pos, &s_minus, &s_plus);
    if (model == FPS_SIGN_SPLIT_NORMAL) return 0.5;
    const double total = s_minus + s_plus;
    return total > 0.0 ? s_plus / total : 0.5;
}

DockingParameters fps_error_estimation_parameters() {
    DockingParameters p;
    p.max_force = 10000.0;      // ProjectData.cs:107 -- the Huber tail is off
    p.clash_tolerance = 0.5;    // ProjectData.cs:108 -- 4x harder than docking
    p.optimize_selected = "All";  // ProjectData.cs:113
    p.shuffle_max_translation = 0.0;  // no random restart in this mode
    // FPS's clash is all-atom over Bondi radii; the coarse term here is one
    // 2.5 A bead per residue, an approximation bought for docking's step cost.
    // A mode whose pose moves by fractions of an angstrom has nothing to spend
    // it on, so it is off.
    p.coarse_clash = false;
    return p;
}

/* ------------------------------------------------------------------------
 * Poses as geometry: RMSD and superposition
 * ------------------------------------------------------------------------ */

namespace {

//! Every atom of an assembly as (body index, body-local coordinate).
/*! Atoms that are not members of a rigid body cannot move with a pose and are
    left out; there are none in an assembly `create_docking_assembly` made. */
struct AssemblyAtoms {
    std::vector<int> body;
    std::vector<IMP::algebra::Vector3D> local;
};

AssemblyAtoms assembly_atoms(const DockingAssembly& assembly) {
    AssemblyAtoms out;
    IMP::Model* m = assembly.get_model();
    const IMP::core::RigidBodies bodies = assembly.get_rigid_bodies();
    std::map<IMP::ParticleIndex, int> index_of;
    for (std::size_t i = 0; i < bodies.size(); ++i) {
        index_of.insert(std::make_pair(bodies[i].get_particle_index(),
                                       static_cast<int>(i)));
    }
    const IMP::atom::Hierarchies leaves =
            IMP::atom::get_leaves(assembly.get_root());
    for (std::size_t i = 0; i < leaves.size(); ++i) {
        const IMP::ParticleIndex p = leaves[i].get_particle_index();
        if (!IMP::core::RigidBodyMember::get_is_setup(m, p)) continue;
        IMP::core::RigidBodyMember member(m, p);
        std::map<IMP::ParticleIndex, int>::const_iterator it = index_of.find(
                member.get_rigid_body().get_particle_index());
        if (it == index_of.end()) continue;
        out.body.push_back(it->second);
        out.local.push_back(member.get_internal_coordinates());
    }
    return out;
}

//! One transformation per rigid body, in `get_rigid_bodies()` order.
/*! An empty or unparseable \p poses_json is the assembly's current pose, which
    is what a caller passing "" means. */
std::vector<IMP::algebra::Transformation3D> pose_frames(
        const DockingAssembly& assembly, const std::string& poses_json) {
    const IMP::core::RigidBodies bodies = assembly.get_rigid_bodies();
    std::vector<IMP::algebra::Transformation3D> out;
    for (std::size_t i = 0; i < bodies.size(); ++i) {
        out.push_back(bodies[i].get_reference_frame().get_transformation_to());
    }
    if (poses_json.empty()) return out;
    const nlohmann::json poses = nlohmann::json::parse(poses_json, NULL, false);
    if (poses.is_discarded() || !poses.is_array()) return out;
    const std::vector<int> body_of_pdb = assembly.get_body_of_pdb();
    for (nlohmann::json::const_iterator it = poses.begin(); it != poses.end();
         ++it) {
        const nlohmann::json& pose = *it;
        if (!pose.is_object() || !pose.contains("t") || !pose.contains("q")) {
            continue;
        }
        const int body_id = pose.value("body_id", -1);
        std::size_t index = out.size();
        for (std::size_t i = 0; i < body_of_pdb.size() && i < out.size(); ++i) {
            if (body_of_pdb[i] == body_id) {
                index = i;
                break;
            }
        }
        if (index >= out.size()) continue;
        const std::vector<double> t = pose["t"].get<std::vector<double> >();
        const std::vector<double> q = pose["q"].get<std::vector<double> >();
        if (t.size() != 3 || q.size() != 4) continue;
        out[index] = IMP::algebra::Transformation3D(
                IMP::algebra::Rotation3D(q[0], q[1], q[2], q[3]),
                IMP::algebra::Vector3D(t[0], t[1], t[2]));
    }
    return out;
}

}  // namespace

double pose_rmsd(const DockingAssembly& assembly, const std::string& poses_a,
                 const std::string& poses_b, bool fps_sign_convention) {
    const AssemblyAtoms atoms = assembly_atoms(assembly);
    if (atoms.body.empty()) return 0.0;
    const std::vector<IMP::algebra::Transformation3D> a =
            pose_frames(assembly, poses_a);
    const std::vector<IMP::algebra::Transformation3D> b =
            pose_frames(assembly, poses_b);
    double sd = 0.0;
    for (std::size_t i = 0; i < atoms.body.size(); ++i) {
        const std::size_t k = static_cast<std::size_t>(atoms.body[i]);
        if (k >= a.size() || k >= b.size()) continue;
        const IMP::algebra::Vector3D& r = atoms.local[i];
        // U*r with U = R_a - R_b, kept as the difference of the two *rotated*
        // vectors so no matrix subtraction is spelled out.
        const IMP::algebra::Vector3D ur =
                a[k].get_rotation().get_rotated(r) -
                b[k].get_rotation().get_rotated(r);
        const IMP::algebra::Vector3D t =
                a[k].get_translation() - b[k].get_translation();
        const IMP::algebra::Vector3D d = fps_sign_convention ? ur - t : ur + t;
        sd += d.get_squared_magnitude();
    }
    return std::sqrt(sd / static_cast<double>(atoms.body.size()));
}

IMP::algebra::Transformation3D pose_superposition(
        const DockingAssembly& assembly, const std::string& poses,
        const std::string& reference_poses) {
    const AssemblyAtoms atoms = assembly_atoms(assembly);
    const std::vector<IMP::algebra::Transformation3D> a =
            pose_frames(assembly, poses);
    const std::vector<IMP::algebra::Transformation3D> b =
            pose_frames(assembly, reference_poses);
    IMP::algebra::Vector3Ds source, target;
    for (std::size_t i = 0; i < atoms.body.size(); ++i) {
        const std::size_t k = static_cast<std::size_t>(atoms.body[i]);
        if (k >= a.size() || k >= b.size()) continue;
        source.push_back(a[k].get_transformed(atoms.local[i]));
        target.push_back(b[k].get_transformed(atoms.local[i]));
    }
    if (source.size() < 3) return IMP::algebra::get_identity_transformation_3d();
    return IMP::algebra::get_transformation_aligning_first_to_second(source,
                                                                    target);
}

/* ------------------------------------------------------------------------
 * The bootstrap itself
 * ------------------------------------------------------------------------ */

std::string BootstrapResult::get_json() const {
    nlohmann::json out;
    out["parent_score"] = parent_score;
    out["parent_e_clash"] = parent_e_clash;
    out["parent_poses"] =
            nlohmann::json::parse(parent_poses.empty() ? "[]" : parent_poses);
    nlohmann::json truth_rows = nlohmann::json::array();
    for (std::size_t i = 0; i < truth.size(); ++i) {
        truth_rows.push_back(nlohmann::json::parse(truth[i].get_json()));
    }
    out["truth"] = truth_rows;
    nlohmann::json rows = nlohmann::json::array();
    for (std::size_t i = 0; i < replicas.size(); ++i) {
        const BootstrapReplica& r = replicas[i];
        nlohmann::json row;
        row["replica"] = r.replica;
        row["score"] = r.score;
        row["e_bond"] = r.e_bond;
        row["e_clash"] = r.e_clash;
        row["converged"] = r.converged;
        row["chi2_r_truth"] = r.chi2_r_truth;
        row["rmsd_to_parent"] = r.rmsd_to_parent;
        row["rmsd_to_parent_fps"] = r.rmsd_to_parent_fps;
        row["output_dir"] = r.output_dir;
        rows.push_back(row);
    }
    out["replicas"] = rows;
    out["rmsd_mean"] = rmsd_mean;
    out["rmsd_sd"] = rmsd_sd;
    out["rmsd_max"] = rmsd_max;
    out["n_perturbed"] = n_perturbed;
    out["n_pinned"] = n_pinned;
    out["extra"] = nlohmann::json::parse(extra.empty() ? "{}" : extra);
    return out.dump();
}

namespace {

//! `<dir>/replica_007` and the like, zero padded so a listing sorts.
std::string numbered_directory(const std::string& dir, const std::string& stem,
                               int index) {
    std::ostringstream name;
    name << dir << "/" << stem << "_" << std::setw(3) << std::setfill('0')
         << index;
    return name.str();
}

//! Delete a directory and the files directly in it; silent on failure.
/*! Only what this module wrote goes in one, and a replica directory that
    cannot be removed is not a reason to fail a run that has already produced
    its numbers. */
void remove_directory(const std::string& dir) {
    static const char* suffixes[] = {".pdb", ".csv", ".json", ".txt", ".pml"};
    for (std::size_t i = 0; i < sizeof(suffixes) / sizeof(suffixes[0]); ++i) {
        const std::vector<std::string> files =
                internal::directory_entries(dir, suffixes[i]);
        for (std::size_t k = 0; k < files.size(); ++k) {
            std::remove(files[k].c_str());
        }
    }
    ::rmdir(dir.c_str());
}

}  // namespace

BootstrapResult fps_bootstrap(const std::vector<std::string>& pdb_paths,
                              const std::string& fps_json_path,
                              const std::string& output_dir,
                              const std::string& parent_poses,
                              const DockingParameters& params,
                              const BootstrapParameters& bootstrap) {
    internal::make_directory(output_dir);
    BootstrapResult out;

    // 1. The parent. Either the caller's pose, or one this run docks first.
    std::string poses = parent_poses;
    if (poses.empty()) {
        const DockingResult parent = dock_minimize(
                pdb_paths, fps_json_path, output_dir + "/parent", params);
        poses = parent.poses;
    }

    // 2. The truth: the parent's own model distances, residuals zeroed. The
    //    assembly is built over *every* distance -- the score set decides what
    //    is optimised, not what the truth is stated for (ErrorEstimation.cs
    //    rewrites R for all of them, `:23-29`).
    //
    //    `mean_position_restraint = false` so the volumes are **not** made
    //    rigid-body members: they must be resampled once the parent pose is in
    //    place, because that is where a replica's run will compute them, and a
    //    volume carried rigidly from the input pose would state a truth the
    //    replicas cannot reproduce.
    const DockingAssembly reference = create_docking_assembly(
            pdb_paths, fps_json_path, std::string(), false, params.ev_weight,
            params.sigma_da, params.clash_tolerance, params.max_force,
            params.clash_radii_source, params.clash_radii_scale);
    apply_poses(reference, poses);
    {
        const IMP::bff::AVs avs = reference.get_network()->get_used_avs();
        for (std::size_t i = 0; i < avs.size(); ++i) {
            IMP::bff::AV(avs[i]).resample();
        }
        reference.get_model()->update();
    }
    const std::vector<PairDistance> truth = collect_pair_distances(
            reference.get_network(), true, params.sigma_da);
    out.parent_poses = poses;
    out.truth = PairDistances(truth.begin(), truth.end());
    // Stated the same way however the parent arrived: the mean-position
    // restraint sum against the *experimental* distances, in this module's
    // half-chi-square units, plus the excluded volume -- the same objective a
    // replica reports, so the two numbers are comparable.
    {
        double parent = clash_energy_of(reference);
        for (std::size_t i = 0; i < truth.size(); ++i) {
            const PairDistance& p = truth[i];
            if (!std::isfinite(p.distance_model)) continue;
            parent += 0.5 * chi2_score_capped(p.distance_model, p.distance_exp,
                                              p.error_neg, p.error_pos,
                                              params.max_force);
        }
        out.parent_score = parent;
        out.parent_e_clash = clash_energy_of(reference);
    }

    // Which distances carry noise. FPS perturbs only the selected ones and
    // scores all of them, so a deselected distance becomes a zero-noise pin;
    // `perturb_deselected` says whether to reproduce that.
    const std::set<std::string> selected = score_set_distance_names(
            reference.get_fps_json_path(), params.score_set);
    std::map<std::string, double> truth_of;
    std::set<std::string> noisy;
    for (std::size_t i = 0; i < truth.size(); ++i) {
        if (!std::isfinite(truth[i].distance_model)) continue;
        truth_of.insert(std::make_pair(truth[i].name, truth[i].distance_model));
        const bool is_selected =
                selected.empty() || selected.count(truth[i].name) > 0;
        if (is_selected || bootstrap.perturb_deselected) {
            noisy.insert(truth[i].name);
        }
    }
    out.n_perturbed = static_cast<int>(noisy.size());
    out.n_pinned = static_cast<int>(truth_of.size() - noisy.size());

    // 3. The replicas. Each is a full docking run against a perturbed copy of
    //    the file, resumed from the parent pose -- `initial_poses` non-empty is
    //    also what suppresses the random shuffle, which this mode must not do.
    const FPSDocument document = read_fps_json(reference.get_fps_json_path());
    nlohmann::json distances =
            nlohmann::json::parse(document.distances, NULL, false);
    if (distances.is_discarded()) distances = nlohmann::json::object();

    const int n = std::max(0, bootstrap.n_replicas);
    std::vector<double> rmsds;
    for (int i = 0; i < n; ++i) {
        const std::string dir = numbered_directory(output_dir, "replica", i);
        internal::make_directory(dir);

        // One RNG stream per replica, derived from the run seed, so a replica
        // is reproducible on its own and the schedule cannot change it.
        nlohmann::json perturbed = distances;
        unsigned int stream = bootstrap.seed + 1000003u * (i + 1);
        for (nlohmann::json::iterator it = perturbed.begin();
             it != perturbed.end(); ++it) {
            std::map<std::string, double>::const_iterator t =
                    truth_of.find(it.key());
            if (t == truth_of.end()) continue;
            double value = t->second;
            if (noisy.count(it.key()) > 0) {
                const double e_neg = it.value().value("error_neg", 0.0);
                const double e_pos = it.value().value("error_pos", 0.0);
                const std::vector<double> draw = sample_distance_perturbations(
                        e_neg, e_pos, 1, stream, bootstrap.perturbation);
                if (!draw.empty()) value += draw[0];
            }
            ++stream;
            it.value()["distance"] = value;
        }
        const std::string file = dir + "/perturbed.fps.json";
        write_fps_json(file, document.positions, perturbed.dump(),
                       document.score_sets, document.extra);

        const DockingResult replica =
                dock_minimize(pdb_paths, file, dir, params, NULL, poses);

        BootstrapReplica row;
        row.replica = i;
        row.score = replica.score;
        row.e_bond = replica.e_bond;
        row.e_clash = replica.e_clash;
        row.converged = replica.converged;
        row.poses = replica.poses;
        row.pairs = replica.pairs;
        row.rmsd_to_parent = pose_rmsd(reference, replica.poses, poses, false);
        row.rmsd_to_parent_fps =
                pose_rmsd(reference, replica.poses, poses, true);
        // What the replica gave up by chasing noise: its pose against the
        // *unperturbed* truth, per distance.
        double chi2 = 0.0;
        int counted = 0;
        for (std::size_t k = 0; k < replica.pairs.size(); ++k) {
            const PairDistance& p = replica.pairs[k];
            std::map<std::string, double>::const_iterator t =
                    truth_of.find(p.name);
            if (t == truth_of.end() || !std::isfinite(p.distance_model)) {
                continue;
            }
            chi2 += chi2_score_capped(p.distance_model, t->second, p.error_neg,
                                      p.error_pos, params.max_force);
            ++counted;
        }
        row.chi2_r_truth = counted > 0
                                   ? chi2 / static_cast<double>(counted)
                                   : internal::nan_value();
        if (bootstrap.keep_replica_dirs) {
            row.output_dir = dir;
        } else {
            remove_directory(dir);
        }
        out.replicas.push_back(row);
        if (std::isfinite(row.rmsd_to_parent)) rmsds.push_back(row.rmsd_to_parent);
    }

    // 4. The spread. Sample standard deviation (n-1): the replicas are draws,
    //    not a population, and with the ten FPS ships the difference is 5 %.
    double sum = 0.0, max = 0.0;
    for (std::size_t i = 0; i < rmsds.size(); ++i) {
        sum += rmsds[i];
        if (rmsds[i] > max) max = rmsds[i];
    }
    out.rmsd_max = max;
    out.rmsd_mean = rmsds.empty() ? internal::nan_value()
                                  : sum / static_cast<double>(rmsds.size());
    if (rmsds.size() > 1) {
        double ss = 0.0;
        for (std::size_t i = 0; i < rmsds.size(); ++i) {
            const double d = rmsds[i] - out.rmsd_mean;
            ss += d * d;
        }
        out.rmsd_sd = std::sqrt(ss / static_cast<double>(rmsds.size() - 1));
    }

    nlohmann::json extra;
    extra["perturbation"] = bootstrap.perturbation == FPS_SIGN_SPLIT_NORMAL
                                    ? "fps_sign_split_normal"
                                    : "split_normal";
    extra["perturb_deselected"] = bootstrap.perturb_deselected;
    extra["seed"] = bootstrap.seed;
    extra["n_replicas"] = n;
    extra["score_set"] = params.score_set;
    extra["optimize_selected"] = params.optimize_selected;
    extra["max_force"] = params.max_force;
    extra["clash_tolerance"] = params.clash_tolerance;
    extra["clash_radii_source"] = params.clash_radii_source;
    extra["clash_radii_scale"] = params.clash_radii_scale;
    extra["output_dir"] = output_dir;
    out.extra = extra.dump();
    return out;
}

IMPBFF_END_NAMESPACE
