/**
 *  \file IMP/bff/AVMeanDistanceRestraint.h
 *  \brief A FRET restraint on the distance between two mean dye positions.
 *
 * The cheap half of #IMP::bff::ProbeNetworkRestraint. That one rebuilds both
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
#include <IMP/bff/ProbeNetworkRestraint.h>
#include <IMP/bff/StatesDistance.h>

#include <IMP/Restraint.h>
#include <IMP/RestraintSet.h>
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
    /*! \note The two particles are \c ParticleIndexAdaptor, which is how
        `IMP::core`'s own restraints take theirs: a `Particle`, an `AV`
        decorator or a bare index all convert. The Python wrapper used to do
        that with a `hasattr` ladder in a `%feature("shadow")`. */
    AVMeanDistanceRestraint(IMP::Model* m, IMP::ParticleIndexAdaptor p1,
                            IMP::ParticleIndexAdaptor p2,
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

//! Make each volume a member of its attachment atom's rigid body.
/*! The coordinates of an accessible volume are the mean of its density, so a
    volume is resampled before it is added: the position it is added *at* is
    the one the rigid body will carry it by. Volumes whose source atom is not
    in a rigid body are left alone. */
IMPBFFEXPORT void add_avs_to_rigid_bodies(const AVs& avs);

//! Give each volume a radius and a mass, so it occupies space.
/*! The radius is the largest of the volume's three; the mass is twice that.
    Without this a volume is a point, and nothing keeps two of them -- or a
    volume and the structure -- from sharing the same coordinates. */
IMPBFFEXPORT void set_av_xyzr_mass(const AVs& avs);

//! Every restraint an fps.json network asks for, as one set.
/*!
    This was `AVNetworkRestraintWrapper`, a `%pythoncode` class subclassing
    `IMP.pmi.restraints.RestraintBase` -- which is why `import IMP.bff` had a
    lazy door in it: IMP.pmi is not one of this module's modules. What the
    wrapper did beyond PMI's bookkeeping is here, and a caller that wants a
    `RestraintBase` subclasses one around this set, in the program that already
    imports PMI.

    \param[in] hier the structure the labelled residues are in
    \param[in] fps_json the network: positions, distances and score sets
    \param[in] name names the restraint, and is the score set's key
    \param[in] score_set which set of distances to score; empty means all
    \param[in] mean_position_restraint score the separation of the volumes'
               *mean positions* (one #IMP::bff::AVMeanDistanceRestraint per
               distance) instead of rebuilding both volumes per evaluation.
               That is an approximation whose content is that the shape of a
               volume does not change when the structure moves, so the volumes
               are made rigid-body members first.
    \param[in] sigma_DA the per-component width for the mean-position branch
    \param[in] occupy_volume give the volumes a radius and a mass
    \param[in] weight multiplies the set's score
    
eturn the set, owned by the caller (an IMP::Object; keep a reference)
    	hrow IOException when \p fps_json cannot be read
*/
IMPBFFEXPORT IMP::RestraintSet* probe_network_restraint_set(
        const IMP::core::Hierarchy& hier, const std::string& fps_json,
        std::string name = "ProbeNetworkRestraint",
        std::string score_set = "", bool mean_position_restraint = false,
        double sigma_DA = 6.0, bool occupy_volume = true, double weight = 1.0);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_AVMEANDISTANCERESTRAINT_H
