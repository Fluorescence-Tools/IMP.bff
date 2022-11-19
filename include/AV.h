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



struct AVPairDistanceMeasurement{

    double distance = -1;
    double error_neg = -1;
    double error_pos = -1;
    double forster_radius = 52.0;
    int distance_type = IMP::bff::DYE_PAIR_DISTANCE_MEAN;
    std::string position_1;
    std::string position_2;

    std::string get_json(){
        // create an empty structure (null)
        nlohmann::json j;
        j["position1_name"] = position_1;
        j["position2_name"] = position_2;
        j["distance"] = distance;
        j["error_neg"] = error_neg;
        j["error_pos"] = error_pos;
        j["Forster_radius"] = forster_radius;
        j["distance_type"] = distance_type;
        return j.dump();
    }

    double score_model(double model){
        auto ev = [](auto f, auto m, auto en, auto ep){
            auto dev = m - f;
            auto w = (dev < 0) ? 1. / en : 1. / ep;
            return .5 * algebra::get_squared(dev * w);
        };
        return 0.5 * ev(model, distance, error_neg, error_pos);
    }

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

    IMP::bff::PathMapHeader create_path_map_header();

    void init_path_map();

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

    static void do_setup_particle(Model *m, ParticleIndex pi,
            ParticleIndex pi_source,
            double linker_length = 20.0,
            const algebra::Vector3D radii = algebra::Vector3D(3.5, 0, 0),
            double linker_width = 0.5,
            double allowed_sphere_radius = 1.5,
            double contact_volume_thickness = 0.0,
            double contact_volume_trapped_fraction = -1,
            double simulation_grid_resolution = 1.5
    ) {
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

    static ParticleIndexKey get_particle_key(unsigned int i) {
        static const ParticleIndexKey k[1] = {
                ParticleIndexKey("Source particle")
        };
        return k[i];
    }

    ParticleIndex get_particle_index() const {
        return Decorator::get_particle_index();
    }

    ParticleIndex get_particle_index(unsigned int i) const {
        return get_model()->get_attribute(get_particle_key(i), get_particle_index());
    }

    Particle* get_particle() const {
        return Decorator::get_particle();
    }

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

    IMP::bff::PathMap* get_map() const;

    void resample(bool shift_xyz=true);

    IMP::algebra::Vector3D get_mean_position() const;

    IMP::algebra::Vector3D get_source_coordinates() const;

    Particle* get_source() const;

};

IMP_DECORATORS(AV, AVs, ParticlesTemp);

template<typename T>
T inline fret_efficiency(T distance, double forster_radius){
    double rda_r0_6 = std::pow(distance / forster_radius, 6.0);
    return 1. / (1. + rda_r0_6);
}

template<typename T>
T inline distance_fret(double fret_efficiency, double forster_radius){
    return forster_radius * std::pow(1. / fret_efficiency - 1.0, 1. / 6.);
}

//! Compute the distance to another accessible volume
IMPBFFEXPORT double av_distance(
        const AV& a,
        const AV& b,
        double forster_radius = 52.0,
        int distance_type=DYE_PAIR_DISTANCE_MEAN,
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


/// Find the particle index of a labeling site
/*!
 *
 * @param[in] json_data A FPS.json position
 * @param[in] hier hierarchy that is searches
 * @return Particle index of the labeling site
 */
IMPBFFEXPORT IMP::ParticleIndex search_labeling_site(
        const IMP::core::Hierarchy& hier,
        std::string json_str = "",
        const nlohmann::json &json_data = nlohmann::json()
);

IMPBFF_END_NAMESPACE

#endif /* IMPBFF_AV_H */