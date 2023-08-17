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



/**
 * @class AVNetworkRestraint
 * @brief A restraint that uses an annotated volumetric network to score particle distances.
 *
 * The AVNetworkRestraint class represents a restraint that utilizes an annotated volumetric network
 * to score distances between particles. It is designed to be used with the IMP library.
 *
 * The restraint is initialized with a hierarchy, a filename of a fps.json file, a name, and an optional
 * score set. The hierarchy is used to obtain the particles involved in the restraint. The fps.json file
 * contains the annotated volumetric network data. The name parameter is used to assign a name to the restraint.
 * The score set parameter specifies the name of the score in the fps.json file to be used for scoring. If no
 * score set is provided, all distances are used for scoring.
 */
class IMPBFFEXPORT AVNetworkRestraint : public IMP::Restraint {

private:

    /**
     * @brief Map of AVs used to compute the score.
     *
     * This map stores the AVs used to compute the score. The keys are the names of the AVs,
     * and the values are pointers to the AV objects.
     */
    std::map<std::string, IMP::bff::AV*> avs_{};
    
    /**
     * @brief ParticleIndexes of AVs used to compute the score.
     *
     * This list stores the ParticleIndexes of AVs used to compute the score. These ParticleIndexes
     * correspond to the AVs stored in the `avs_` map.
     */
    IMP::ParticleIndexes av_pi_;

    /**
     * @brief ParticleIndexes of particles contributing to the score.
     *
     * This list stores the ParticleIndexes of particles that contribute to the score computation.
     * These particles are not necessarily AVs.
     */
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

    /**
     * Get the accessible volume (AV) for a labeled particle.
     * @param name The name of the labeled particle.
     * @return The AV associated with the labeled particle.
     */
    IMP::bff::AV* get_av(std::string name) const;

public:

    /**
     * @brief Constructs an AVNetworkRestraint object.
     * @param[in] hier The hierarchy used to obtain particles.
     * @param[in] fps_json_fn The filename of the fps.json file.
     * @param[in] name The name of this restraint. Default is "AVNetworkRestraint%1%".
     * @param[in] score_set The name of the score in the fps.json file. If not provided, all distances are used for scoring.
     */
    AVNetworkRestraint(
        const IMP::core::Hierarchy &hier,
        std::string fps_json_fn,
        std::string name = "AVNetworkRestraint%1%",
        std::string score_set = ""
    );

    /**
     * @brief Returns exp(score).
     * @return The exponential of the score.
     */
    double get_probability() const {
        return exp(-unprotected_evaluate(nullptr));
    }

    /**
     * @brief Returns the used Atom::AVs.
     * @return The used Atom::AVs.
     */
    const IMP::bff::AVs get_used_avs();

    /**
     * @brief Returns the used experimental distances.
     * @return The used experimental distances.
     */
    const std::map<std::string, AVPairDistanceMeasurement> get_used_distances(){
        return distances_;
    }

    /**
     * @brief Returns the model distance (or FRET efficiency) between two dyes.
     * @param[in] position1_name The name of the first dye position.
     * @param[in] position2_name The name of the second dye position.
     * @param[in] forster_radius The Förster radius.
     * @param[in] distance_type The type of distance calculation.
     * @return The model distance (or FRET efficiency) between the two dyes.
     */
    double get_model_distance(
            std::string position1_name,
            std::string position2_name,
            double forster_radius,
            int distance_type
    ) const;


    /**
     * @brief Returns the particle indexes of the AVs.
     * @return The particle indexes.
     */
    ParticleIndexes const get_indexes(){
        return av_pi_;
    }

    /**
     * @brief Evaluates the restraint.
     * @param[in] accum The derivative accumulator.
     * @return The score of the restraint.
     */
    virtual double unprotected_evaluate(IMP::DerivativeAccumulator *accum) const override;

    /**
     * @brief Returns the inputs required by the restraint.
     * @return The inputs required by the restraint.
     */
    virtual IMP::ModelObjectsTemp do_get_inputs() const override;

    /**
     * @brief Prints a description of the restraint.
     * @param[in] out The output stream.
     */
    void show(std::ostream &out) const {out << "AVNetwork restraint";}

    IMP_OBJECT_METHODS(AVNetworkRestraint)

};


IMPBFF_END_NAMESPACE


#endif //IMPBFF_AVNETWORKRESTRAINT_H
