/**
 * \file AVMeanDistanceRestraint.cpp
 * \brief A FRET restraint on the distance between two mean dye positions.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/AVMeanDistanceRestraint.h>

#include <IMP/core/XYZ.h>

IMPBFF_BEGIN_NAMESPACE

AVMeanDistanceRestraint::AVMeanDistanceRestraint(
        IMP::Model* m, IMP::ParticleIndex p1, IMP::ParticleIndex p2,
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

IMP::ModelObjectsTemp AVMeanDistanceRestraint::do_get_inputs() const {
    IMP::ModelObjectsTemp out;
    out.push_back(get_model()->get_particle(p1_));
    out.push_back(get_model()->get_particle(p2_));
    return out;
}

IMPBFF_END_NAMESPACE
