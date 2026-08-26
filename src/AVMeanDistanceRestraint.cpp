/**
 * \file AVMeanDistanceRestraint.cpp
 * \brief A FRET restraint on the distance between two mean dye positions.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/AVMeanDistanceRestraint.h>
#include <IMP/bff/internal/Text.h>

#include <IMP/atom/Mass.h>
#include <IMP/core/XYZ.h>
#include <IMP/core/XYZR.h>
#include <IMP/core/rigid_bodies.h>

#include <algorithm>
#include <map>

IMPBFF_BEGIN_NAMESPACE

AVMeanDistanceRestraint::AVMeanDistanceRestraint(
        IMP::Model* m, IMP::ParticleIndexAdaptor p1,
        IMP::ParticleIndexAdaptor p2,
        const AVPairDistanceMeasurement& measurement, double sigma,
        double weight)
    : IMP::Restraint(m, "AVMeanDistanceRestraint %1%"), p1_(p1), p2_(p2),
      measurement_(measurement),
      // The lookup runs to 2.5 R0. Past that the efficiency is 2e-4 and the
      // inverse is not a function any more, so extending the table buys
      // resolution in a region no measurement constrains.
      converter_(measurement.forster_radius, sigma, 1.0,
                 2.5 * measurement.forster_radius),
      weight_(weight) {}

double AVMeanDistanceRestraint::score_at(double d_mp) const {
    return measurement_.score_model(
            converter_.get_effective_distance(d_mp, measurement_.distance_type));
}

double AVMeanDistanceRestraint::unprotected_evaluate(
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

void add_avs_to_rigid_bodies(const AVs& avs) {
    for (unsigned int i = 0; i < avs.size(); ++i) {
        AV av = avs[i];
        av.resample();
        IMP::Particle* source = av.get_source();
        if (IMP::core::RigidBodyMember::get_is_setup(source)) {
            IMP::core::RigidBodyMember(source).get_rigid_body().add_member(
                    av.get_particle());
        }
    }
}

void set_av_xyzr_mass(const AVs& avs) {
    for (unsigned int i = 0; i < avs.size(); ++i) {
        AV av = avs[i];
        IMP::Particle* p = av.get_particle();
        const IMP::algebra::Vector3D r = AV(av).get_radii();
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
    const AVs avs = net->get_used_avs();

    if (!mean_position_restraint) {
        rs->add_restraint(net);
    } else {
        // The volumes are carried by the rigid bodies of the atoms they are
        // attached to, which is what makes scoring their mean positions --
        // rather than rebuilding them -- an approximation and not a fiction.
        add_avs_to_rigid_bodies(avs);
        std::map<std::string, AV> by_name;
        for (unsigned int i = 0; i < avs.size(); ++i) {
            by_name[avs[i].get_particle()->get_name()] = avs[i];
        }
        const std::map<std::string, AVPairDistanceMeasurement> distances =
                net->get_used_distances();
        for (std::map<std::string, AVPairDistanceMeasurement>::const_iterator
                     it = distances.begin();
             it != distances.end(); ++it) {
            const AVPairDistanceMeasurement& d = it->second;
            std::map<std::string, AV>::const_iterator a1 =
                    by_name.find(d.position_1);
            std::map<std::string, AV>::const_iterator a2 =
                    by_name.find(d.position_2);
            IMP_USAGE_CHECK(a1 != by_name.end() && a2 != by_name.end(),
                            "distance " << it->first << " names a position "
                                        << "the network has no volume for");
            rs->add_restraint(new AVMeanDistanceRestraint(
                    m, a1->second.get_particle(), a2->second.get_particle(), d,
                    sigma_DA));
        }
    }
    if (occupy_volume) set_av_xyzr_mass(avs);
    rs->set_weight(weight);
    return rs.release();
}

IMP::ModelObjectsTemp AVMeanDistanceRestraint::do_get_inputs() const {
    IMP::ModelObjectsTemp out;
    out.push_back(get_model()->get_particle(p1_));
    out.push_back(get_model()->get_particle(p2_));
    return out;
}

IMPBFF_END_NAMESPACE
