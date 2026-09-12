/**
 *  \file IMP/bff/ProbeNetworkRestraint.h
 *  \brief Simple restraint for networks of accessible volumes.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPBFF_PROBENETWORKRESTRAINT_H
#define IMPBFF_PROBENETWORKRESTRAINT_H

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
#include <cereal/types/vector.hpp>
#include <IMP/UnaryFunction.h>

#include <IMP/bff/ProbeAccessibleVolumeDecorator.h>
#include <IMP/bff/internal/ThreadPool.h>
#include <IMP/bff/internal/FPSLegacyIO.h>
#include <IMP/bff/internal/json.h>

#include <vector>
#include <algorithm>

IMPBFF_BEGIN_NAMESPACE

namespace internal { struct AVEvalJob; }



/**
 * @class ProbeNetworkRestraint
 * @brief A restraint that uses an annotated volumetric network to score particle distances.
 *
 * The ProbeNetworkRestraint class represents a restraint that utilizes an annotated volumetric network
 * to score distances between particles. It is designed to be used with the IMP library.
 *
 * The restraint is initialized with a hierarchy, a filename of a fps.json file, a name, and an optional
 * score set. The hierarchy is used to obtain the particles involved in the restraint. The fps.json file
 * contains the annotated volumetric network data. The name parameter is used to assign a name to the restraint.
 * The score set parameter specifies the name of the score in the fps.json file to be used for scoring. If no
 * score set is provided, all distances are used for scoring.
 */
class IMPBFFEXPORT ProbeNetworkRestraint : public IMP::Restraint {

    friend class cereal::access;

    /* `avs_` holds bare AV decorator handles, which are views onto particles
       that the Model owns, not state of their own. They are therefore stored
       as the particle indices they decorate and rebuilt against the restored
       model on load, rather than archived as pointers. */
    template<class Archive> void serialize(Archive &ar) {
        ar(cereal::base_class<IMP::Restraint>(this),
           n_samples, av_pi_, model_ps_, distances_,
           space_fixed_, shared_map_, distance_, quad_k_, search_grid_factor_,
           search_stencil_, search_mode_, points_, atom_points_);
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
                avs_[kv.first].reset(new IMP::bff::ProbeAccessibleVolumeDecorator(get_model(), kv.second));
            }
            registry_ = nullptr;
            configure_avs();
        }
    }

    IMP_OBJECT_SERIALIZE_DECL(ProbeNetworkRestraint);

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
    int quad_k_ = 50;              //!< representative points per cloud
    int search_grid_factor_ = 1;   //!< ProbeAccessibleVolumeDecorator::set_search_grid_factor for every AV
    int search_stencil_ = 26;      //!< ProbeAccessibleVolumeDecorator::set_search_stencil for every AV
    std::string search_mode_ = "dijkstra";  //!< ProbeAccessibleVolumeDecorator::set_search_mode for every AV

    //! Shared occupancy rasters (only under space_fixed && shared_map)
    IMP::Pointer<ProbeAccessibleVolumeOccupancyRegistry> registry_;

    //! Number of unprotected_evaluate() calls
    mutable long n_evaluations_ = 0;

    //! Cumulative wall time of the evaluation phases (seconds), diagnostics
    mutable double t_registry_ = 0, t_prepare_ = 0, t_compute_ = 0, t_pairs_ = 0;

    //! Threads for the AVs' compute phases (0 = hardware concurrency)
    int n_threads_ = 0;

    //! The three phases of an evaluation: serial begin (Model reads, occupancy
    //! classification, AV prepare), the pool run, serial finish (occupancy
    //! end_update, AV finish, score).
    std::shared_ptr<internal::AVEvalJob> begin_evaluation() const;
    void run_evaluation(internal::AVEvalJob &job) const;
    double finish_evaluation(internal::AVEvalJob &job) const;
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
    /* Owned. A bare `new AV(...)` with no destructor on this class leaks one
       decorator per labelled position for the lifetime of the process. */
    std::map<std::string, std::unique_ptr<IMP::bff::ProbeAccessibleVolumeDecorator> > avs_{};

    /**
     * @brief Positions the network carries as a **point**, not a volume.
     *
     * FPS's `AVSimlationType.None`: a position that has no cloud to average
     * over. Two spellings reach it, and they differ only in where the point
     * comes from:
     *
     *  - `simulation_type = "XYZ"` -- a mean dye position measured once, in
     *    some other structure's frame. A particle of its own is created at
     *    `x`/`y`/`z`; if the position declares `reference_atoms` the frame is
     *    Kabsch-fitted onto this structure first and the coordinate is
     *    transported through the fit (`FilterEngine.cs:225-278`).
     *  - `simulation_type = "ATOM"` -- an atom of the structure itself
     *    (`LabelingPositions.cs:203`). No new particle: the atom *is* the
     *    point, so it moves with the body for free.
     *
     * A distance touching either is scored as **R_mp** whatever the file's
     * `distance_type` says (`FilterEngine.cs:305-307`) -- there is no
     * distribution at that end for an average to be taken over.
     */
    std::map<std::string, IMP::ParticleIndex> points_{};

    //! Which of #points_ are `ATOM` positions -- backed by a structure atom.
    /*! The `XYZ` ones are not: their particle is this restraint's own. The
        distinction is what FPS's bond rule tests (`SpringEngine.cs:111-116`),
        and what says whether a radius may be written on the particle. */
    std::vector<std::string> atom_points_{};

    //! The volume of a position, or `nullptr` when it has none. Silent.
    /*! #get_av warns, which is right for a name that should have resolved and
        wrong for a position that is deliberately a point. */
    IMP::bff::ProbeAccessibleVolumeDecorator* find_av(const std::string &name) const;

    //! Where a position is right now: an AV's mean position, or a point.
    IMP::algebra::Vector3D get_position_coordinates(const std::string &name) const;


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
     *
     *  Positions of `simulation_type` `XYZ` or `ATOM` are **not** volumes and
     *  do not appear in the returned map; they go to #points_. Sending them
     *  through the labelling-site search is what used to make a network
     *  containing one unbuildable: an `XYZ` position names no chain, residue
     *  or atom, so the selection matched the whole structure and
     *  IMP::bff::search_labeling_site threw *ambiguous labelling site*.
     */
    std::map<std::string, std::unique_ptr<IMP::bff::ProbeAccessibleVolumeDecorator> > create_av_decorated_particles(
            nlohmann::json used_positions,
            const IMP::core::Hierarchy &hier);

    /**
     * Get the accessible volume (AV) for a labeled particle.
     * @param name The name of the labeled particle.
     * @return The AV associated with the labeled particle.
     */
    IMP::bff::ProbeAccessibleVolumeDecorator* get_av(std::string name) const;

public:

    /**
     * @brief The AV handle of a labeled position, sharing this restraint's
     * path map and lattice state (a copy of the handle, not a fresh one).
     * @param name The name of the labeled position (fps.json key)
     */
    IMP::bff::ProbeAccessibleVolumeDecorator get_used_av(std::string name) const;

    /**
     * @brief Constructs an ProbeNetworkRestraint object.
     * @param[in] hier The hierarchy used to obtain particles.
     * @param[in] fps_json_fn The filename of the fps.json file.
     * @param[in] name The name of this restraint. Default is "ProbeNetworkRestraint%1%".
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
     * @param[in] search_grid_factor Coarsening of the path search
     *   (ProbeAccessibleVolumeDecorator::set_search_grid_factor); 1 = exact on the AV grid.
     * @param[in] search_stencil 26 (symmetric, default) or 30 (historical,
     *   asymmetric); see ProbeAccessibleVolumeDecorator::set_search_stencil.
     * @param[in] search_mode "dijkstra" (default, path search) or
     *   "euclidean" (straight linker, visibility only); see ProbeAccessibleVolumeDecorator::set_search_mode.
     */
    ProbeNetworkRestraint(
        const IMP::core::Hierarchy &hier,
        std::string fps_json_fn,
        std::string name = "ProbeNetworkRestraint%1%",
        std::string score_set = "",
        int n_samples = 50000,
        bool space_fixed = true,
        bool shared_map = true,
        std::string distance = "quad",
        int quad_k = 50,
        int search_grid_factor = 1,
        int search_stencil = 26,
        std::string search_mode = "dijkstra"
    );

    bool get_space_fixed() const { return space_fixed_; }
    bool get_shared_map() const { return shared_map_; }
    std::string get_distance_method() const { return distance_; }
    int get_quad_k() const { return quad_k_; }
    int get_search_grid_factor() const { return search_grid_factor_; }
    int get_search_stencil() const { return search_stencil_; }
    std::string get_search_mode() const { return search_mode_; }
    int get_n_samples() const { return n_samples; }

    //! The shared occupancy registry (nullptr unless shared_map)
    ProbeAccessibleVolumeOccupancyRegistry *get_occupancy_registry() const { return registry_; }

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
     * @brief Start an evaluation and return before it finishes.
     *
     * Everything that reads the Model (coordinates, parameters, occupancy
     * classification) is done before this returns; the searches, carves and
     * pair sums then run on the pool while the caller does something else --
     * typically loading the next frame -- and wait_score() collects the
     * score and writes the AV mean positions. Between the two calls the
     * caller may change the Model freely. Bit-identical to
     * unprotected_evaluate(). Distance sets that use PROBE_PAIR_DISTANCE_MP or
     * PROBE_PAIR_XYZ_DISTANCE read the Model at pair time and are therefore
     * evaluated synchronously inside evaluate_async().
     */
    void evaluate_async() const;
    //! Finish the evaluation started by evaluate_async() and return its score
    double wait_score() const;
    //! True between evaluate_async() and wait_score()
    bool get_has_pending_evaluation() const { return (bool) job_; }

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
    ProbeNetworkRestraint() {}

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
    const IMP::bff::ProbeAccessibleVolumeDecorators get_used_avs();

    /**
     * @brief Returns the used experimental distances.
     * @return The used experimental distances.
     */
    const std::map<std::string, AVPairDistanceMeasurement> get_used_distances(){
        return distances_;
    }

    //! The positions carried as a point rather than a volume: name -> particle.
    /*! `XYZ` and `ATOM` positions (see #points_). Empty for the ordinary
        all-AV network, which is why nothing had to know about it before. */
    std::map<std::string, IMP::ParticleIndex> get_point_positions() const {
        return points_;
    }

    //! The names of #get_point_positions, sorted. The Python-facing spelling.
    std::vector<std::string> get_point_position_names() const {
        std::vector<std::string> names;
        names.reserve(points_.size());
        for(const auto &kv : points_) names.push_back(kv.first);
        return names;
    }

    //! Names of the point positions that are an atom of the structure (`ATOM`).
    std::vector<std::string> get_atom_position_names() const {
        return atom_points_;
    }

    //! The particle standing for a position, volume or point.
    /*! An AV's particle carries its mean position, an `ATOM` position's *is*
        the atom and an `XYZ` position's is a fixed point this restraint owns.
        A name the network does not know gives a default-constructed index. */
    IMP::ParticleIndex get_position_particle_index(std::string name) const;

    //! True when \p name is a position with no volume (`XYZ` or `ATOM`).
    bool get_position_is_point(std::string name) const {
        return points_.find(name) != points_.end();
    }

    //! FPS's **bond** test, on one of the used distances.
    /*!
        `SpringEngine.cs:111-116`: a distance is a bond iff **both** ends are
        plain atoms -- no dye, a real atom, and no accessible volume. Here that
        is: both positions are `ATOM` positions. Such a restraint is a
        crosslink or a covalent tie between subunits, not a FRET measurement,
        and FPS treats it differently in two ways -- its two anchor atoms are
        excluded from clash detection (#IMP::bff::set_bond_anchor_radii) and
        its energy is reported separately as `Ebond`, **a subset of the total
        and never an addition to it**.

        An `XYZ` position does not qualify: it is a coordinate, not an atom,
        and FPS's test requires `AtomID > 0`.

        \param[in] distance_name a key of #get_used_distances
        \return false for a name the network does not use
    */
    bool get_is_bond(std::string distance_name) const;

    //! The names of the used distances that are bonds, in map order.
    std::vector<std::string> get_bond_names() const;

    //! The candidate pair names, in the order #get_pair_efficiencies reports.
    /*! The score set's pairs sorted by name. Experiment planning builds an
        `(n_frames, n_pairs)` efficiency matrix over an ensemble, and which
        column is which pair has to be stated somewhere; it is stated here,
        rather than by each caller rediscovering that a `std::map` iterates
        sorted. */
    std::vector<std::string> get_pair_names() const;

    //! Mean FRET efficiency of every candidate pair, at the current coordinates.
    /*! Each pair with its own Förster radius, and as an *efficiency* whatever
        distance type the file asked to be scored on -- planning weighs a pair
        by how far apart it places conformers in the observable, and the
        expected error #IMP::bff::select_probe_pairs takes is quoted in
        efficiency units.

        The accessible volumes are re-evaluated first, so a caller that has
        just loaded a new frame does not have to score it to read it.
    */
    std::vector<double> get_pair_efficiencies() const;

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
    void show(std::ostream &out) const {out << "ProbeNetworkRestraint";}

    IMP_OBJECT_METHODS(ProbeNetworkRestraint)

private:
    /* The evaluation in flight (evaluate_async / wait_score). Declared last
       on purpose: members are destroyed in reverse order, so an unfinished
       job -- whose runner thread uses the pool and the AVs -- is joined
       (AVEvalJob's destructor) before anything it touches goes away. */
    mutable std::shared_ptr<internal::AVEvalJob> job_;

};


IMPBFF_END_NAMESPACE


#endif //IMPBFF_PROBENETWORKRESTRAINT_H
