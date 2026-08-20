/**
 *  \file IMP/bff/AVMeanDistanceRestraint.h
 *  \brief A FRET restraint on the distance between two mean dye positions.
 *
 * The cheap half of #IMP::bff::AVNetworkRestraint. That one rebuilds both
 * accessible volumes on every evaluation, which is the whole cost of an
 * AV-restrained optimisation; this one scores the separation of the two volumes'
 * *mean positions*, converts it to the modelled observable through a cached
 * #IMP::bff::FRETDistanceConverter, and scores that with the measurement's
 * asymmetric chi-squared.
 *
 * That is an approximation, and its content is that the shape of the volume does
 * not change when the structure moves — true while the volumes are carried as
 * rigid-body members and false as soon as the linker's environment changes.
 *
 * **Derivatives are computed.** The score depends on the coordinates only
 * through \f$d_{mp}\f$, so \f$dS/dd_{mp}\f$ comes from a 1-D central difference
 * through the cached converter and maps onto the two particles analytically as
 * \f$\partial d_{mp}/\partial x = \pm\hat r\f$. That is what lets an IMP
 * optimiser dock by minimisation rather than only by Monte-Carlo; the Python
 * this replaces had two copies of the class, one with the gradient and one
 * without, and the one without was the one the wrapper used.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_AVMEANDISTANCERESTRAINT_H
#define IMPBFF_AVMEANDISTANCERESTRAINT_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/AV.h>
#include <IMP/bff/StatesDistance.h>

#include <IMP/Restraint.h>
#include <IMP/particle_index.h>

IMPBFF_BEGIN_NAMESPACE

//! Score the separation of two mean dye positions against a measurement.
class IMPBFFEXPORT AVMeanDistanceRestraint : public IMP::Restraint {
    IMP::ParticleIndex p1_, p2_;
    AVPairDistanceMeasurement measurement_;
    FRETDistanceConverter converter_;
    double weight_;

    //! The chi-squared at a given mean-position separation.
    double score_at(double d_mp) const;

public:
    //! \param[in] m the model the two particles belong to
    /*! \param[in] p1,p2 the two dye particles; their `XYZ` coordinates are the
               mean positions
        \param[in] measurement the experimental distance, its asymmetric error,
               its Förster radius and which distance convention it is in
        \param[in] sigma the per-component width of the separation vector, for
               the \f$R_{mp} \to \langle R_{DA}\rangle\f$ conversion
        \param[in] weight multiplies the score */
    AVMeanDistanceRestraint(IMP::Model* m, IMP::ParticleIndex p1,
                            IMP::ParticleIndex p2,
                            const AVPairDistanceMeasurement& measurement,
                            double sigma = 6.0, double weight = 1.0);

    double get_weight_factor() const { return weight_; }
    void set_weight_factor(double w) { weight_ = w; }
    const FRETDistanceConverter& get_converter() const { return converter_; }

    virtual double unprotected_evaluate(
            IMP::DerivativeAccumulator* accum) const override;
    virtual IMP::ModelObjectsTemp do_get_inputs() const override;

    IMP_OBJECT_METHODS(AVMeanDistanceRestraint);
};

IMPBFF_END_NAMESPACE

#endif //IMPBFF_AVMEANDISTANCERESTRAINT_H
