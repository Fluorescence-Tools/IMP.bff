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

#include <memory>
#include <IMP/score_functor/distance_pair_score_macros.h>

#include <IMP/Model.h>
#include <IMP/Restraint.h>
#include <IMP/Object.h>
#include <IMP/Pointer.h>
#include <IMP/atom/Hierarchy.h>

#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <IMP/UnaryFunction.h>

#include <IMP/bff/AV.h>
#include <IMP/bff/internal/ThreadPool.h>
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

    friend class cereal::access;

    /* `avs_` holds bare AV decorator handles, which are views onto particles
       that the Model owns, not state of their own. They are therefore stored
       as the particle indices they decorate and rebuilt against the restored
       model on load, rather than archived as pointers. */
    template<class Archive> void serialize(Archive &ar) {
        ar(cereal::base_class<IMP::Restraint>(this),
           n_samples, av_pi_, model_ps_, distances_,
           space_fixed_, shared_map_, distance_, quad_k_);
        // On save this is built from avs_; on load ar() overwrites it and the
        // decorators are rebuilt from it below. A single serialize() (rather
        // than a save/load pair) is required here because IMP::Restraint
        // already supplies one, and cereal rejects two candidate functions.
        std::map<std::string, IMP::ParticleIndex> av_index;
        for (const auto &kv : avs_) {
            av_index[kv.first] = kv.second->get_particle_index();
        }
        ar(av_index);
        if (std::is_base_of<cereal::detail::InputArchiveBase, Archive>::value) {
            avs_.clear();
            for (const auto &kv : av_index) {
                avs_[kv.first].reset(new IMP::bff::AV(get_model(), kv.second));
            }
            registry_ = nullptr;
            configure_avs();
        }
    }

    IMP_OBJECT_SERIALIZE_DECL(AVNetworkRestraint);

private:

    /**
     * @brief Number of random samples in distance computation
     *
     * The number of random samples used to compute a distance.
     * larger numbers increase the precision of the distance computation.
     */
    int n_samples = 50000;

    /* PRD-105 evaluation modes. */
    bool space_fixed_ = true;      //!< lattice-anchored windows (else legacy)
    bool shared_map_ = true;       //!< one occupancy raster per (spacing, extra) class
    std::string distance_ = "quad";//!< "quad" (lattice quadrature) or "mc"
    int quad_k_ = 100;             //!< representative points per cloud

    //! Shared occupancy rasters (only under space_fixed && shared_map)
    IMP::Pointer<AVOccupancyRegistry> registry_;

    //! Number of unprotected_evaluate() calls
    mutable long n_evaluations_ = 0;

    //! Threads for the AVs' compute phases (0 = hardware concurrency)
    int n_threads_ = 0;
    //! Persistent workers, created on first threaded evaluation
    mutable std::shared_ptr<internal::ThreadPool> pool_;
    internal::ThreadPool &get_pool() const;

    //! Apply the mode flags to the AV handles (anchoring, registry)
    void configure_avs();

    /**
     * @brief Map of AVs used to compute the score.
     *
     * This map stores the AVs used to compute the score. The keys are the names of the AVs,
     * and the values are pointers to the AV objects.
     */
    /* Owned. These used to be bare `new AV(...)` with no destructor on this
       class at all, so every restraint leaked one decorator per labelled
       position for the lifetime of the process. */
    std::map<std::string, std::unique_ptr<IMP::bff::AV> > avs_{};
    
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
    std::map<std::string, std::unique_ptr<IMP::bff::AV> > create_av_decorated_particles(
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
     * @brief The AV handle of a labeled position, sharing this restraint's
     * path map and lattice state (a copy of the handle, not a fresh one).
     * @param name The name of the labeled position (fps.json key)
     */
    IMP::bff::AV get_used_av(std::string name) const;

    /**
     * @brief Constructs an AVNetworkRestraint object.
     * @param[in] hier The hierarchy used to obtain particles.
     * @param[in] fps_json_fn The filename of the fps.json file.
     * @param[in] name The name of this restraint. Default is "AVNetworkRestraint%1%".
     * @param[in] score_set The name of the score in the fps.json file. If not provided, all distances are used for scoring.
     */
    /**
     * @param[in] n_samples Random samples for `distance="mc"`.
     * @param[in] space_fixed Anchor every AV grid on the global lattice
     *   (default). `false` selects the deprecated legacy anchoring.
     * @param[in] shared_map One occupancy raster per (spacing, extra-radius)
     *   class shared by all AVs (default). Requires `space_fixed`.
     * @param[in] distance `"quad"` (deterministic lattice quadrature, default)
     *   or `"mc"` (random sampling, order-dependent).
     * @param[in] quad_k Representative points per cloud for `"quad"`.
     */
    AVNetworkRestraint(
        const IMP::core::Hierarchy &hier,
        std::string fps_json_fn,
        std::string name = "AVNetworkRestraint%1%",
        std::string score_set = "",
        int n_samples = 50000,
        bool space_fixed = true,
        bool shared_map = true,
        std::string distance = "quad",
        int quad_k = 100
    );

    bool get_space_fixed() const { return space_fixed_; }
    bool get_shared_map() const { return shared_map_; }
    std::string get_distance_method() const { return distance_; }
    int get_quad_k() const { return quad_k_; }
    int get_n_samples() const { return n_samples; }

    //! The shared occupancy registry (nullptr unless shared_map)
    AVOccupancyRegistry *get_occupancy_registry() const { return registry_; }

    /**
     * @brief Threads used for the AVs' searches under `space_fixed`.
     *
     * The AVs of one evaluation are independent once their occupancy
     * sources are up to date, so their compute phases (obstacle spheres,
     * Dijkstra, carve, cloud) run on threads; results are identical to the
     * serial evaluation. 0 (default) = hardware concurrency, 1 = serial.
     */
    void set_number_of_threads(int n) { n_threads_ = n; }
    int get_number_of_threads() const;

    /**
     * @brief Diagnostics of the last run as JSON.
     *
     * Per AV: skip / local / full / roll counts; per shared occupancy map:
     * skip / local / full / grow counts, moved beads (last, total), extent;
     * plus the mode flags and the number of evaluations.
     */
    std::string get_diagnostics_json() const;

    /**
     * @brief Estimate of the quadrature error of the model distances.
     *
     * Largest absolute difference, over the used distances, between the
     * distance at `quad_k` and at `reference_k` representative points
     * (both from the current maps). Costs O(reference_k^2) per pair.
     */
    double get_quad_error_estimate(int reference_k = 1000) const;

    //! Default constructor, needed to deserialize the restraint.
    AVNetworkRestraint() {}

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
