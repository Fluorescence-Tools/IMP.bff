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

#include <IMP/atom/Hierarchy.h>
#include <IMP/UnaryFunction.h>

#include <IMP/bff/AV.h>
#include <IMP/bff/internal/FPSReaderWriter.h>
#include <IMP/bff/internal/json.h>

#include <vector>
#include <algorithm>
#include <functional> // placeholders

IMPBFF_BEGIN_NAMESPACE



class IMPBFFEXPORT AVNetworkRestraint : public IMP::Restraint {

private:

    IMP::ParticleIndexes model_ps_;
    IMP::ParticleIndexes av_pi_;
    IMP::core::Hierarchy hier_;
    ParticleIndexes pis_;
    nlohmann::json fps_json_;

    std::vector<IMP::bff::AV*> avs_;
    std::map<std::string, AVPairDistanceMeasurement> distances_;
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

    /// Find and decorate labeled particles with accessible volume (AVs)
    /** This is automatically called by the constructor.
     *  You only need to call this if you change parameters of
     *  AVs (e.g., the linker length).
     */
    void create_av_decorated_particles(nlohmann::json used_positions);

    const IMP::bff::AVs get_used_avs();

    const std::map<std::string, AVPairDistanceMeasurement> get_used_distances(){
        return distances_;
    }

    ParticleIndexes const get_indexes(){
        return pis_;
    }

    virtual double unprotected_evaluate(
            IMP::DerivativeAccumulator *accum) const override;

    virtual IMP::ModelObjectsTemp do_get_inputs() const override{
        IMP::ModelObjectsTemp ret;
        for (size_t i = 0; i < model_ps_.size(); i++) {
            ret.push_back(get_model()->get_particle(model_ps_[i]));
        }
        return ret;
    }

    void show(std::ostream &out) const {out << "FRETNetwork restraint";}

    IMP_OBJECT_METHODS(AVNetworkRestraint)

};


IMPBFF_END_NAMESPACE


#endif //IMPBFF_AVNETWORKRESTRAINT_H
