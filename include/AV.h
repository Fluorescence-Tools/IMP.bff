/**
 * \file IMP/bff/AV.h
 * \brief Simple Accessible Volume decorator.
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_AV_H
#define IMPBFF_AV_H

#include <IMP/bff/bff_config.h>

#include <IMP/decorator_macros.h>
#include <IMP/Decorator.h>
#include <IMP/algebra/Vector3D.h>
#include <IMP/algebra/Transformation3D.h>
#include <IMP/core/XYZ.h>
#include <IMP/log.h>


#include <IMP/bff/PathMap.h>

#include <string>
#include <cmath>
#include <vector>
#include <limits>
#include <iostream>            // std::cout, std::cout, std::flush

#include <IMP/bff/internal/json.h>
#include <IMP/bff/internal/InverseSampler.h>
// requires C++14
// #include <boost/histogram.hpp> // make_histogram, regular, weight, indexed
#include <IMP/bff/internal/Histogram.h>

IMPBFF_BEGIN_NAMESPACE


/// Different types of distances between two accessible volumes
typedef enum{
    DYE_PAIR_DISTANCE_E,            /// Mean FRET averaged distance R_E
    DYE_PAIR_DISTANCE_MEAN,         /// Mean distance <R_DA>
    DYE_PAIR_DISTANCE_MP,           /// Distance between AV mean positions
    DYE_PAIR_EFFICIENCY,            /// Mean FRET efficiency
    DYE_PAIR_DISTANCE_DISTRIBUTION, /// Distance distribution
    DYE_PAIR_XYZ_DISTANCE           /// Distance between XYZ of dye particles
} DyePairMeasures;



/// Container for experimental distance measurement
class IMPBFFEXPORT AVPairDistanceMeasurement{

public:

    double distance = -1;
    double error_neg = -1;
    double error_pos = -1;
    double forster_radius = 52.0;
    int distance_type = IMP::bff::DYE_PAIR_DISTANCE_MEAN;
    std::string position_1;
    std::string position_2;

    /// Get a JSON string
    std::string get_json();

    /**
    @brief Score a model distance against the experimental distance.
    @param model The model distance to be scored.
    @return The score of the model distance.
    */
    double score_model(double model);

};



//! A decorator for a particle with accessible volume (AV).
/** Using the decorator one can get and set AV
parameters and modify derivatives.

 AV must have IMP.Hierarchy parent with XYZ -> is labeling site
 AV coordinates = AV mean position

\ingroup helper
\ingroup decorators
\include AV_Decorator.py

*/
class IMPBFFEXPORT AV : public IMP::core::Gaussian {

private:

    IMP::bff::PathMap* av_map_ = nullptr;

protected:

    IMP::algebra::VectorD<9> get_parameter() const {
        IMP::algebra::VectorD<9> v;
        v[0] = get_linker_length();
        v[1] = get_radius1();
        v[2] = get_radius2();
        v[3] = get_radius3();
        v[4] = get_linker_width();
        v[5] = get_allowed_sphere_radius();
        v[6] = get_contact_volume_thickness();
        v[7] = get_contact_volume_trapped_fraction();
        v[8] = get_simulation_grid_resolution();
        return v;
    }

public:

    /**
    @brief Creates a path map header.
    @return The created path map header.
    */
    IMP::bff::PathMapHeader create_path_map_header();

    /**
     * @brief Initializes the path map.
     *
     * This function initializes the path map used by the AV system. The path map is a 
     * data structure that stores information about the available paths or routes in the system.
     *
     * @return void
     */
    void init_path_map();

    /**
     * @brief Get the FloatKey object for a specific AV feature.
     *
     * This function returns a FloatKey object representing a specific AV feature
     * based on the given index.
     *
     * @param i The index of the AV feature.
     * @return The FloatKey object for the AV feature at the given index.
     *
     * @note The valid range for the index is 0 to 8.
     * @note If the index is out of range, an error message will be generated.
     */
    static FloatKey get_av_key(unsigned int i){
        IMP_USAGE_CHECK(i < 9, "Out of range av feature");
        static const FloatKey k[] = {
                FloatKey("linker_length"),
                FloatKey("radius1"),
                FloatKey("radius2"),
                FloatKey("radius3"),
                FloatKey("linker_width"),
                FloatKey("allowed_sphere_radius"),
                FloatKey("contact_volume_thickness"),
                FloatKey("contact_volume_trapped_fraction"),
                FloatKey("simulation_grid_resolution")
        };
        return k[i];
    }

    /**
     * @brief Sets up the attributes for a particle in a model for AV (Anisotropic Volume) calculations.
     *
     * This function sets up the attributes for a particle in a model for AV calculations. 
     * The attributes include linker length, radii, linker width, allowed sphere radius, 
     * contact volume thickness, contact volume trapped fraction, and simulation grid resolution. 
     * It also sets up a particle attribute to store the source particle index.
     *
     * @param m The model in which the particle resides.
     * @param pi The index of the particle to set up.
     * @param pi_source The index of the source particle.
     * @param linker_length The length of the linker.
     * @param radii The radii of the particle in the x, y, and z directions.
     * @param linker_width The width of the linker.
     * @param allowed_sphere_radius The radius of the allowed sphere.
     * @param contact_volume_thickness The thickness of the contact volume.
     * @param contact_volume_trapped_fraction The fraction of the contact volume that is trapped.
     * @param simulation_grid_resolution The resolution of the simulation grid.
     */
    static void do_setup_particle(Model *m, ParticleIndex pi,
                                ParticleIndex pi_source,
                                double linker_length = 20.0,
                                const algebra::Vector3D radii = algebra::Vector3D(3.5, 0, 0),
                                double linker_width = 0.5,
                                double allowed_sphere_radius = 1.5,
                                double contact_volume_thickness = 0.0,
                                double contact_volume_trapped_fraction = -1,
                                double simulation_grid_resolution = 1.5) {
        if (!IMP::core::Gaussian::get_is_setup(m, pi)) {
            IMP::core::Gaussian::setup_particle(m, pi);
        }
        m->add_attribute(get_av_key(0), pi, linker_length);
        m->add_attribute(get_av_key(1), pi, radii[0]);
        m->add_attribute(get_av_key(2), pi, radii[1]);
        m->add_attribute(get_av_key(3), pi, radii[2]);
        m->add_attribute(get_av_key(4), pi, linker_width);
        m->add_attribute(get_av_key(5), pi, allowed_sphere_radius);
        m->add_attribute(get_av_key(6), pi, contact_volume_thickness);
        m->add_attribute(get_av_key(7), pi, contact_volume_trapped_fraction);
        m->add_attribute(get_av_key(8), pi, simulation_grid_resolution);
        m->add_attribute(get_particle_key(0), pi, pi_source);
    }

    /**
     * @brief Get the particle key for the specified index.
     * @param i The index of the particle key.
     * @return The particle key at the specified index.
     */
    static ParticleIndexKey get_particle_key(unsigned int i) {
        static const ParticleIndexKey k[1] = {
            ParticleIndexKey("Source particle")
        };
        return k[i];
    }

    /**
     * @brief Get the particle index of the AV object.
     * @return The particle index of the AV object.
     */
    ParticleIndex get_particle_index() const {
        return Decorator::get_particle_index();
    }

    /**
     * @brief Get the particle index of the AV object at the specified index.
     * @param i The index of the AV object.
     * @return The particle index of the AV object at the specified index.
     */
    ParticleIndex get_particle_index(unsigned int i) const {
        return get_model()->get_attribute(get_particle_key(i), get_particle_index());
    }

    /**
     * @brief Get the particle pointer of the AV object.
     * @return The particle pointer of the AV object.
     */
    Particle* get_particle() const {
        return Decorator::get_particle();
    }

    /**
     * @brief Get the particle pointer of the AV object at the specified index.
     * @param i The index of the AV object.
     * @return The particle pointer of the AV object at the specified index.
     */
    Particle* get_particle(unsigned int i) const {
        return get_model()->get_particle(get_particle_index(i));
    }
    IMP_DECORATOR_METHODS(AV, IMP::core::Gaussian);

    /** Setup the particle with unspecified AV - uses default values. */
    IMP_DECORATOR_SETUP_1(AV, IMP::ParticleIndex, pi_source);
    IMP_DECORATOR_GET_SET(linker_length, get_av_key(0), Float, Float);
    IMP_DECORATOR_GET_SET(radius1, get_av_key(1), Float, Float);
    IMP_DECORATOR_GET_SET(radius2, get_av_key(2), Float, Float);
    IMP_DECORATOR_GET_SET(radius3, get_av_key(3), Float, Float);
    IMP_DECORATOR_GET_SET(linker_width, get_av_key(4), Float, Float);
    IMP_DECORATOR_GET_SET(allowed_sphere_radius, get_av_key(5), Float, Float);
    IMP_DECORATOR_GET_SET(contact_volume_thickness, get_av_key(6), Float, Float);
    IMP_DECORATOR_GET_SET(contact_volume_trapped_fraction, get_av_key(7), Float, Float);
    IMP_DECORATOR_GET_SET(simulation_grid_resolution, get_av_key(8), Float, Float);


    /**
     * @brief Returns the radii of an object.
     *
     * This function returns the radii of an object as a 3D vector. The radii are obtained by calling the functions get_radius1(), get_radius2(), and get_radius3() and storing the values in a Vector3D object.
     *
     * @return A Vector3D object representing the radii of the object.
     */
    IMP::algebra::Vector3D get_radii(){
        return IMP::algebra::Vector3D({get_radius1(), get_radius3(), get_radius3()});
    }

    //! Get whether the coordinates are optimized
    /** \return true only if all of them are optimized.
      */
    bool get_parameters_are_optimized() const {
        return
            get_particle()->get_is_optimized(get_av_key(0)) &&
            get_particle()->get_is_optimized(get_av_key(1)) &&
            get_particle()->get_is_optimized(get_av_key(2)) &&
            get_particle()->get_is_optimized(get_av_key(3)) &&
            get_particle()->get_is_optimized(get_av_key(6)) &&
            get_particle()->get_is_optimized(get_av_key(7));
    }

    //! Set whether the coordinates are optimized
    void set_av_parameters_are_optimized(bool tf) const {
        get_particle()->set_is_optimized(get_av_key(0), tf);
        get_particle()->set_is_optimized(get_av_key(1), tf);
        get_particle()->set_is_optimized(get_av_key(2), tf);
        get_particle()->set_is_optimized(get_av_key(3), tf);
        get_particle()->set_is_optimized(get_av_key(6), tf);
        get_particle()->set_is_optimized(get_av_key(7), tf);
    }

    /**
     * @brief Sets the AV parameter using a JSON object.
     *
     * This function takes a JSON object as input and sets the AV parameter accordingly.
     *
     * @param j The JSON object containing the AV parameter.
     */
    void set_av_parameter(const nlohmann::json &j);


    //! Get the vector of derivatives accumulated by add_to_derivatives().
    /** Somewhat suspect based on wanting a Point/Vector differentiation
        but we don't have points */
    algebra::Vector3D get_derivatives() const {
        return get_model()->get_coordinate_derivatives(get_particle_index());
    }

    static bool get_is_setup(Model *m, ParticleIndex pi) {
        return m->get_has_attribute(get_av_key(2), pi);
    }

    /**
     * @brief Get the PathMap associated with the AV object.
     * @return A pointer to the PathMap object.
     */
    IMP::bff::PathMap* get_map() const;

    /**
     * @brief Resample the AV object.
     * @param shift_xyz Flag indicating whether to shift the XYZ coordinates.
     */
    void resample(bool shift_xyz=true);

    /**
     * @brief Get the mean position of the AV object.
     * @return The mean position as a Vector3D.
     */
    IMP::algebra::Vector3D get_mean_position(bool include_source=true) const;

    /**
     * @brief Get the source coordinates of the AV object.
     * @return The source coordinates as a Vector3D.
     */
    IMP::algebra::Vector3D get_source_coordinates() const;

    /**
     * @brief Get the source particle of the AV object.
     * @return A pointer to the source Particle object.
     */
    Particle* get_source() const;

};

IMP_DECORATORS(AV, AVs, ParticlesTemp);

/**
 * @brief Computes the FRET efficiency given the distance and Forster radius.
 * @tparam T The type of the distance.
 * @param distance The distance between the two volumes.
 * @param forster_radius The Forster radius.
 * @return The FRET efficiency.
 */
template<typename T>
T inline fret_efficiency(T distance, double forster_radius){
    double rda_r0_6 = std::pow(distance / forster_radius, 6.0);
    return 1. / (1. + rda_r0_6);
}

/**
 * @brief Computes the distance between two volumes given the FRET efficiency and Forster radius.
 * @tparam T The type of the distance.
 * @param fret_efficiency The FRET efficiency.
 * @param forster_radius The Forster radius.
 * @return The distance between the two volumes.
 */
template<typename T>
T inline distance_fret(double fret_efficiency, double forster_radius){
    return forster_radius * std::pow(1. / fret_efficiency - 1.0, 1. / 6.);
}

/**
 * @brief Computes the distance to another accessible volume.
 * @param a The first accessible volume.
 * @param b The second accessible volume.
 * @param forster_radius The Forster radius.
 * @param distance_type The type of distance to compute.
 * @param n_samples The number of samples to use for distance computation.
 * @return The distance between the two accessible volumes.
 */
IMPBFFEXPORT double av_distance(
        const AV& a,
        const AV& b,
        double forster_radius = 52.0,
        int distance_type = DYE_PAIR_DISTANCE_MEAN,
        int n_samples = 10000
);

// Draw random points in AV. Returns (x,y,z,d) vector
IMPBFFEXPORT std::vector<double> av_random_points(
        const AV& av1,
        int n_samples=10000
);

//! Random sampling over AV/AV distances
IMPBFFEXPORT std::vector<double> av_random_distances(
        const AV& av1,
        const AV& av2,
        int n_samples=10000
);


//! Compute the distance to another accessible volume
IMPBFFEXPORT std::vector<double> av_distance_distribution(
        const AV& av1,
        const AV& av2,
        std::vector<double> axis,
        //double start, double stop, int n_bins, // for boost histogram
        int n_samples=10000
);


/**
 * @brief Find the particle index of a labeling site.
 *
 * This function searches for the particle index of a labeling site in a given hierarchy.
 *
 * @param[in] hier The hierarchy in which the labeling site is searched.
 * @param[in] json_str The JSON string containing the FPS.json position (optional).
 * @param[in] json_data The JSON data containing the FPS.json position (optional).
 *
 * @return The particle index of the labeling site.
 *
 * @note If both json_str and json_data are provided, json_str will be used.
 */
IMPBFFEXPORT IMP::ParticleIndex search_labeling_site(
        const IMP::core::Hierarchy& hier,
        std::string json_str = "",
        const nlohmann::json &json_data = nlohmann::json()
);


IMPBFF_END_NAMESPACE

#endif /* IMPBFF_AV_H */