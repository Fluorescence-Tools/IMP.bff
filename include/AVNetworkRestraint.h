/**
 *  \file IMP/bff/AVNetworkRestraint.h
 *  \brief Simple restraint for networks of accessible volumes.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPBFF_AVNETWORKRESTRAINT_H
#define IMPBFF_AVNETWORKRESTRAINT_H

#include <IMP/bff/bff_config.h>
#include <IMP/score_functor/distance_pair_score_macros.h>

#include <IMP/Model.h>
#include <IMP/Restraint.h>
#include <IMP/Object.h>
#include <IMP/Pointer.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/UnaryFunction.h>

#include <IMP/bff/AV.h>
#include <IMP/bff/internal/FPSReaderWriter.h>
#include <IMP/bff/internal/json.h>

#include <vector>
#include <algorithm>

IMPBFF_BEGIN_NAMESPACE



class IMPBFFEXPORT AVNetworkRestraint : public IMP::Restraint {

private:

    /// Map to AVs that are used to compute the score
    std::map<std::string, IMP::bff::AV*> avs_{};
    /// IMP ParticleIndexes of AVs used to compute the score
    IMP::ParticleIndexes av_pi_;

    /// ParticleIndexes particles contributing to the score
    IMP::ParticleIndexes model_ps_;

    /// Map of experimental distance measurements (incl. errors)
    std::map<std::string, AVPairDistanceMeasurement> distances_;

    /// Find and decorate labeled particles with accessible volume (AVs)
    /** This is method is automatically called by the constructor.
     *  You only need to call this if you change parameters of
     *  AVs (e.g., the linker length).
     */
    std::map<std::string, IMP::bff::AV*> create_av_decorated_particles(
            nlohmann::json used_positions,
            const IMP::core::Hierarchy &hier);

    IMP::bff::AV* get_av(std::string name) const;

public:

    /// Setup the AVNetworkRestraint
    /**
     * @param[in] hier Hierarchy used to obtain particles
     * @param[in] fps_json_fn Filename of fps.json file
     * @param[in] name Name of this restraint
     * @param[in] score_set Name of score in the fps.json file (of no score
     * is provided all distances are used for scoring)
     */
    AVNetworkRestraint(
            const IMP::core::Hierarchy &hier,
            std::string fps_json_fn,
            std::string name="AVNetworkRestraint%1%",
            std::string score_set = ""
    );

    /// Returns exp(score)
    double get_probability() const {
        return exp(-unprotected_evaluate(nullptr));
    }

    const IMP::bff::AVs get_used_avs();

    /// Get used experimental distances
    const std::map<std::string, AVPairDistanceMeasurement> get_used_distances(){
        return distances_;
    }

    /// Get the model distance (or FRET efficiency) between two dyes
    double get_model_distance(
            std::string position1_name,
            std::string position2_name,
            double forster_radius,
            int distance_type
    ) const;


    ParticleIndexes const get_indexes(){
        return av_pi_;
    }

    virtual double unprotected_evaluate(IMP::DerivativeAccumulator *accum) const override;

    virtual IMP::ModelObjectsTemp do_get_inputs() const override;

    void show(std::ostream &out) const {out << "AVNetwork restraint";}

    IMP_OBJECT_METHODS(AVNetworkRestraint)

};


IMPBFF_END_NAMESPACE


#endif //IMPBFF_AVNETWORKRESTRAINT_H
