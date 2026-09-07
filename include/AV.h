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

#include <IMP/Pointer.h>

#include <IMP/bff/Base.h>

#include <IMP/decorator_macros.h>
#include <IMP/Decorator.h>
#include <IMP/algebra/Vector3D.h>
#include <IMP/algebra/Transformation3D.h>
#include <IMP/core/XYZ.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/log.h>


#include <IMP/bff/AVModel.h>
#include <IMP/bff/AVBuilder.h>
#include <IMP/bff/FPS.h>
#include <IMP/bff/PathMap.h>
#include <IMP/bff/AVOccupancyMap.h>
#include <IMP/bff/VdwRadii.h>
#include <IMP/bff/internal/AVLatticeState.h>

#include <algorithm>
#include <memory>
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

    /* Ref-counted, not a bare pointer. AV is an IMP *decorator* -- a value
       handle constructed, copied and discarded freely -- and it was holding a
       raw `new PathMap` with no destructor on this class at all. Every fresh
       handle over the same particle built and then abandoned an entire path
       map: measured at ~1 MB a time, 21 MB for twenty handles. A destructor
       would have been worse, not better, because two copies of a decorator
       would then free the same map twice. PathMap derives from IMP::Object and
       is therefore already ref-counted, so IMP::Pointer gives copies a shared
       map and frees it when the last handle goes. */
    IMP::Pointer<IMP::bff::PathMap> av_map_;

    /* Bookkeeping of the lattice (space-fixed) evaluation path -- window,
       occupancy sources, what the tiles were computed from, quadrature cache,
       counters. Shared between copies of a handle like the map is; created
       together with the map in init_path_map(). */
    std::shared_ptr<internal::AVLatticeState> state_;

    void resample_legacy(bool shift_xyz);
    void resample_lattice(bool shift_xyz, bool force_full);
    // the three phases of resample_lattice (see AVLatticeState::pending)
    void resample_lattice_prepare(bool shift_xyz, bool force_full);
    void resample_lattice_compute(bool split_stages);
    void resample_lattice_compute_carve();
    void resample_lattice_finish();

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
                                //! Negative derives it; see
                                //! get_effective_allowed_sphere_radius().
                                double allowed_sphere_radius = -1.0,
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
    //! The accessible **contact** volume: how deep the surface layer is, and
    //! how much of the dye's time it holds.
    /*! `contact_volume_thickness` is the depth in Angstrom of the layer above
        the surface the dye cannot enter; a cloud voxel is *in contact* when
        excluded volume lies within that distance of it.
        `contact_volume_trapped_fraction` is the share of the cloud's total
        weight those voxels then carry, the rest sharing the remainder -- a dye
        that touches the protein stays there longer than free diffusion would
        put it. Both have to ask for it: a thickness of 0 or a negative
        fraction leaves the volume unweighted, and a fraction at or above 1
        (which would leave the free part weightless) is refused with a warning.

        This is Olga's ACV, not FPS's -- **FPS has no contact volume at all**.
        The rule and the two deliberate differences from Olga's source are in
        IMP::bff::PathMap::apply_contact_weighting(); the one worth knowing
        here is that the layer is quantised to whole grid steps, so a thickness
        below `simulation_grid_resolution` is no layer at all.

        The cloud does not change -- the same voxels, with different weights.
        What moves is the mean position, towards the surface. Measured on 3GUN
        with the seventeen fitted fractions of
        `examples/structure/T4L/fret.fps.json`: the mean moves 2.3 A on average
        (0.04-4.6 A) and *closer to the nearest atom* at sixteen of the
        seventeen sites, by 1.8 A. Those fractions are all well above the
        geometric contact share (0.13-0.47 of the cloud, mean 0.20), so the
        layer is up-weighted, and the 33 <R_DA> of that file all shorten, by
        3.04 A on average.

        Inert until 2026-09-01 (PRD-121 G9): both fields were read, stored, and
        never used. Wiring them up immediately found a second defect behind the
        first -- set_av_parameter() read the fraction with an `int` default, so
        every fps.json value was truncated to 0 or 1 -- and answered the +2 A
        offset against Zenodo 3376527's published <R_DA> that
        `okf/validation/fps_screening_ab.md` had recorded as unexplained: that
        table was computed *with* the contact volume its file asks for, and
        honouring it brings this module from +2.06 A to **+0.22 A** of it. */
    IMP_DECORATOR_GET_SET(contact_volume_thickness, get_av_key(6), Float, Float);
    IMP_DECORATOR_GET_SET(contact_volume_trapped_fraction, get_av_key(7), Float, Float);
    IMP_DECORATOR_GET_SET(simulation_grid_resolution, get_av_key(8), Float, Float);

    //! The selection whose atoms are not obstacles for this volume.
    /*! The fps.json `strip_mask` field. A labelling site's own side chain
        walls in its probe, so the position's parameters -- the linker length
        and above all the `allowed_sphere_radius` -- are calibrated against a
        structure with it removed. An empty mask strips nothing.

        The source atom is kept whatever the mask says: it is the anchor the
        volume grows from, and removing it would let the cloud grow through
        its own attachment point. */
    //! Whether the volume weights its voxels by the linker's chain statistics.
    /*! An accessible volume treats every reachable voxel as equally likely;
        a real linker does not. `chain_weighting` turns on the correction --
        #IMP::bff::linker_weighting() for this volume's linker length, applied
        to every voxel by its path length. Off by default, so the number a
        volume reports does not change unless it is asked for.

        This is the fps.json `chain_weighting` field. */
    static IntKey get_chain_weighting_key() {
        static const IntKey k("chain_weighting");
        return k;
    }

    bool get_chain_weighting() const {
        return get_model()->get_has_attribute(get_chain_weighting_key(),
                                              get_particle_index()) &&
               get_model()->get_attribute(get_chain_weighting_key(),
                                          get_particle_index()) != 0;
    }

    //! Turn chain weighting on or off; the volume is stale afterwards.
    void set_chain_weighting(bool on) {
        if (get_model()->get_has_attribute(get_chain_weighting_key(),
                                           get_particle_index())) {
            get_model()->set_attribute(get_chain_weighting_key(),
                                       get_particle_index(), on ? 1 : 0);
        } else {
            get_model()->add_attribute(get_chain_weighting_key(),
                                       get_particle_index(), on ? 1 : 0);
        }
        if (!on) {
            if (av_map_) av_map_->set_linker_weighting(LinkerWeighting());
            return;
        }
        const LinkerWeighting w = linker_weighting(get_linker_length());
        const double reach = w.get_supported_fraction(get_linker_length());
        // A table column fitted for a much longer tether puts none of its mass
        // inside a dye linker's reach. Weighting by it then does not correct
        // the volume, it selects the volume's outer shell -- so say so rather
        // than let a 13 A shift in the mean position pass for a refinement.
        if (reach < 0.01) {
            IMP_WARN("chain_weighting: the tabulated weighting for a linker of "
                     << get_linker_length() << " A carries only "
                     << reach * 100.0
                     << "% of its weight within that reach, so it selects the "
                        "volume's outer shell rather than reweighting it. See "
                        "okf/validation/chain_weighting.md."
                     << std::endl);
        }
        if (av_map_) av_map_->set_linker_weighting(w);
    }

    static StringKey get_strip_mask_key() {
        static const StringKey k("strip_mask");
        return k;
    }

    std::string get_strip_mask() const {
        return get_model()->get_has_attribute(get_strip_mask_key(),
                                              get_particle_index())
                       ? get_model()->get_attribute(get_strip_mask_key(),
                                                    get_particle_index())
                       : std::string();
    }

    void set_strip_mask(const std::string& mask) {
        if (get_model()->get_has_attribute(get_strip_mask_key(),
                                           get_particle_index())) {
            get_model()->set_attribute(get_strip_mask_key(),
                                       get_particle_index(), mask);
        } else {
            get_model()->add_attribute(get_strip_mask_key(),
                                       get_particle_index(), mask);
        }
        // The obstacle set is built from it, so a new mask invalidates the map.
        av_map_ = nullptr;
    }


    /**
     * @brief Returns the radii of an object.
     *
     * This function returns the radii of an object as a 3D vector. The radii are obtained by calling the functions get_radius1(), get_radius2(), and get_radius3() and storing the values in a Vector3D object.
     *
     * @return A Vector3D object representing the radii of the object.
     */
    IMP::algebra::Vector3D get_radii(){
        return IMP::algebra::Vector3D({get_radius1(), get_radius2(), get_radius3()});
    }

    /**
     * @brief The dye radii that actually take part in the carve.
     *
     * `radius1` always; `radius2` and `radius3` only when positive. **Zero is
     * this library's sentinel for "unused", not a sphere of radius zero**, and
     * the difference matters: LabelLib sorts the radii it is given and grades by
     * all of them (`FlexLabel/src/FlexLabel.cxx:235`), so a zero radius is a real
     * probe that fits everywhere and `dyeDensityAV3(r, 0, 0)` is *not*
     * `dyeDensityAV1(r)` -- it is `(1 + 1 + [r fits])/3`, measured here as 21625
     * voxels against AV1's 16370. bff writes AV1 as `(r, 0, 0)` and the fps
     * schema requires `radius2`/`radius3` to be positive for an `AV3` position
     * (see #IMP::bff::fps_json_schema), so selecting the positive radii
     * reproduces both
     * conventions: one radius gives the AV1 carve, three give the AV3 carve.
     */
    /**
     * @brief Linker-length scale that removes a stencil's path-length bias.
     *
     * A coarse stencil overestimates path length, so the accessible volume it
     * returns is that of a **shorter** linker. The bias is a property of the
     * stencil, not of the structure: measured obstacle-free over linker lengths
     * 12-25 A and dye radii 1.0-3.5 A, the scale that makes the 26 stencil
     * reproduce the 74 stencil's volume is **1.0551 +/- 0.0021** (range
     * 1.0526-1.0592). Stencil 30 measures the same, being 26's asymmetric
     * variant.
     *
     * Applying it recovers the volume without paying for the finer search: on
     * T4L 172L site 22 the 26 stencil gives 0.835-0.844 of the reference volume
     * uncompensated and **0.995** compensated, and obstacle-free it is within
     * 0.8 % at every length tested.
     *
     * 1.0 for the reference stencil, so the default configuration is untouched.
     */
    double get_stencil_length_compensation() const {
        if(!get_compensate_stencil()) return 1.0;
        switch(get_search_stencil()){
            case 26: case 30: return 1.0551;
            default: return 1.0;          // 74 is the reference
        }
    }

    //! Linker length the search actually uses: nominal, times the stencil
    //! compensation. Equal to get_linker_length() for the default stencil.
    double get_effective_linker_length() const {
        return get_linker_length() * get_stencil_length_compensation();
    }

    //! The source clearance the search actually uses.
    /*!
        A **negative** stored `allowed_sphere_radius` means *derive it*, and
        deriving is the default because the value that works depends on two
        other parameters. The path search inflates every obstacle by half the
        linker width, so the free sphere around the attachment atom has to
        clear that inflation plus a grid step of slack, or the source voxel is
        walled in and the volume comes back **empty**.

        This lived in `get_av_from_structure()` only, so the two doors onto
        the same volume disagreed: the fps.json door derived, the decorator
        door (and `imp_bff av-export` behind it) took a flat 1.5 and returned
        nothing at FPS's standard linker width of 4.5 Å. One place now.

        The FPS mapping, for the record: FPS seeds its search from a sphere of
        `LinkerInitialSphere * W` with `LinkerInitialSphere = 0.5`
        (`av_routines.cpp:136`), i.e. exactly `W/2` — the same half-width, with
        no grid slack because FPS's seed is unconditional rather than carved
        out of an inflated obstacle map.
     */
    double get_effective_allowed_sphere_radius() const {
        double r = get_allowed_sphere_radius();
        if(r >= 0.0) return r;
        return std::max(1.5, 0.5 * get_linker_width()
                                     + 0.5 * get_simulation_grid_resolution());
    }

    //! Whether to correct a coarse stencil's volume bias. **Default false**:
    //! selecting stencil 26 or 30 gives that stencil's own answer, so a caller
    //! pinning the historical metric keeps it. Turn it on when 26 is chosen for
    //! speed and the reference volume is still wanted. A no-op on stencil 74.
    bool get_compensate_stencil() const;
    void set_compensate_stencil(bool tf);
    static IntKey get_compensate_stencil_key();

    std::vector<double> get_active_radii() const {
        std::vector<double> r;
        r.push_back(get_radius1());
        if(get_radius2() > 0.0) r.push_back(get_radius2());
        if(get_radius3() > 0.0) r.push_back(get_radius3());
        return r;
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

    //! Take this volume's parameters from one fps.json position object.
    /*! \param[in] json_text the position, as JSON text -- the file's own
                   bytes, which is what every other reader in this module
                   takes. Reads `linker_length`, `radius1..3`, `linker_width`,
                   `allowed_sphere_radius`, `contact_volume_thickness`,
                   `contact_volume_trapped_fraction`,
                   `simulation_grid_resolution` and `strip_mask`.
        \throw ValueException when the text is not JSON, or the mask cannot be
               read. */
    void set_av_parameter(const std::string &json_text);


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
     *
     * Under `space_fixed` (the default) the path map is anchored on the
     * global lattice: voxel centres at integer multiples of the grid
     * spacing, window centred on the lattice-quantised source and rolled by
     * whole voxels only; occupancy is read from a shared per-class raster
     * when a registry is set (see set_occupancy_registry) or maintained
     * privately by exact deltas; the search is skipped altogether when
     * nothing that feeds it changed. With `space_fixed` off the legacy
     * source-anchored path runs unchanged.
     *
     * @param shift_xyz Flag indicating whether to shift the XYZ coordinates.
     * @param force_full Recompute everything from scratch (full raster of
     *        every occupancy source, no skip). The result is bit-identical
     *        to the incremental path; used to prove it.
     */
    void resample(bool shift_xyz=true, bool force_full=false);

    /**
     * @brief Whether the path map is anchored on the global lattice.
     *
     * Default true. `false` selects the legacy anchoring of the grid to the
     * source sub-voxel, byte-for-byte as before PRD-105; it is deprecated
     * and kept for one transition period.
     */
    bool get_space_fixed() const;
    void set_space_fixed(bool tf);

    //! IntKey holding the space_fixed flag (absent = default, true)
    static IntKey get_space_fixed_key();

    /**
     * @brief Coarsening factor of the path search (lattice path only).
     *
     * 1 (default): the search runs on the AV's own grid. f > 1: the search
     * runs on the lattice points whose index is a multiple of f (spacing
     * f*h; the same absolute lattice, so windows still roll by whole
     * voxels), and every fine tile takes the cost of the cheapest
     * neighbouring coarse point plus the straight hop to it; the obstacle
     * carve, the cloud and the mean stay on the fine grid. About f^3 fewer
     * search nodes; the AV changes at the coarse-voxel level (narrow
     * channels below f*h are lost, path lengths are quantised coarser).
     * An approximation -- not covered by the bit-exact contract across
     * factors, but exact and deterministic for a given factor.
     */
    int get_search_grid_factor() const;
    void set_search_grid_factor(int f);
    static IntKey get_search_grid_factor_key();

    /**
     * @brief Neighbour stencil of the path search (lattice path only).
     *
     * 26 (default): face, edge and corner neighbours (radius sqrt 3), no
     * length-2 jumps, so a path cannot tunnel through a one-voxel wall.
     * 30: the historical stencil -- radius 2 built by loops running
     * -2 <= d < 2, i.e. the +2 axis face missing and the tile itself
     * included; a path could cross a one-voxel wall towards -x/-y/-z only.
     * The legacy anchoring always uses 30. On T4L @2 A the 30-stencil AVs
     * are 15-25 % larger (they leak through the inflated surface layer) and
     * their means sit up to 2.5 A away from the 26-stencil ones.
     */
    int get_search_stencil() const;
    void set_search_stencil(int stencil);
    static IntKey get_search_stencil_key();

    /**
     * @brief Search algorithm of the lattice path.
     *
     * "dijkstra" (default): the path search -- a tile is reached when a
     * path through free space of length <= linker length exists.
     * "euclidean": the linker is straight -- a tile is reached when the
     * source voxel sees it (straight voxel path free of obstacles) and it
     * lies within the linker length; no path search at all. Tiles in the
     * shadow of the protein are dropped. Faster by the search's share of
     * the frame; a different, more restrictive model.
     */
    std::string get_search_mode() const;
    void set_search_mode(std::string mode);
    static IntKey get_search_mode_key();

    /**
     * @brief Which van der Waals radii the obstacles are inflated by.
     *
     * `"imp"` (**the default**) uses whatever radius each particle carries,
     * which after `IMP::atom::read_pdb` is IMP's CHARMM-derived
     * **united-atom** set: a carbon of 1.85-2.275 A, standing in for hydrogens
     * that are not in the file.
     * `"olga"` takes every atom's radius from Olga's name-keyed table,
     * #IMP::bff::olga_vdw_radius -- carbon 1.70 A, nitrogen 1.625, oxygen
     * 1.49, sulphur 1.782, phosphorus 1.86, hydrogen 1.00, and a flat 1.50 A
     * for a name the table does not carry.
     *
     * **Why IMP's is the default** (owner, 2026-09-01, reversing a one-day
     * default of `"olga"`): the excluded-volume half of a docking score is
     * `clash_container`, which measures overlap through #IMP::core::XYZR --
     * the *particles'* radii. A volume built on Olga's table and a clash term
     * built on IMP's are two halves of one score that disagree about how big
     * an atom is, and neither half reports the disagreement. Consistency with
     * IMP docking outweighs agreement with Olga's published numbers, and the
     * cost of that choice is real and written down:
     * `okf/validation/fps_screening_ab.md` measures it (nine of 33 published
     * <R_DA> lose their model value, and the agreement with the Zenodo table
     * worsens from -0.02 A / 0.71 A rmsd to +0.22 A / 0.91 A).
     *
     * `"olga"` is therefore kept and selectable, and is what a caller
     * reproducing Olga-era numbers -- the fitted
     * `contact_volume_trapped_fraction` of a `.fps.json`, the published
     * <R_DA> of Zenodo 3376527 -- should ask for. It also *opens volumes up*:
     * its heavy atoms are smaller, so sites IMP's united-atom radii wall in
     * come back with a cloud.
     *
     * The choice is per volume, but it is not free to differ *within* one
     * shared occupancy raster: the raster is one obstacle set for every volume
     * in it, so an #IMP::bff::AVOccupancyRegistry adopts the first radii
     * vector handed to it and throws if a second volume asks for another. It
     * also invalidates the volume, like every other parameter that feeds the
     * raster.
     */
    std::string get_radii_source() const;
    void set_radii_source(std::string source);
    static IntKey get_radii_source_key();

    /**
     * @brief Read occupancy from a shared per-class raster registry.
     *
     * Requires `space_fixed`. Pass nullptr to return to a private raster.
     */
    void set_occupancy_registry(AVOccupancyRegistry *registry);
    AVOccupancyRegistry *get_occupancy_registry() const;

    //! The lattice window: {kx, ky, kz, n} (lattice index of voxel 0, edge)
    std::vector<int> get_lattice_window() const;

    /**
     * @brief Announce this AV's window to the shared occupancy maps.
     *
     * Lets a registry grow once to the union of all windows before the
     * first map is rasterised, instead of once per AV. No-op without a
     * registry or under legacy anchoring.
     */
    void prepare_lattice_window();

#ifndef SWIG
    /**
     * @brief The lattice evaluation split for a caller that runs several AVs
     * concurrently: prepare (touches the Model: coordinates, parameters,
     * occupancy updates -- serial), compute (this AV's map only -- safe on a
     * thread), finish (writes the mean position -- serial).
     * resample() is prepare + compute + finish; the split is a no-op for
     * legacy anchoring, where prepare() runs the whole legacy path.
     */
    void resample_prepare(bool shift_xyz=true, bool force_full=false);
    void resample_compute();
    void resample_finish();
    //! compute() in two halves -- the search, then the carve/cloud/mean --
    //! so a scheduler can interleave the halves of different AVs
    void resample_compute_search();
    void resample_compute_carve();
    //! True between a prepare() that decided to recompute and its compute()
    bool get_has_pending_compute() const { return state_ && state_->pending; }
    //! Build (or refresh) the cached quadrature representation for `k`
    void prepare_quadrature(int k) const;
    //! Tell prepare() that the caller updates the registry maps itself
    void set_registry_driven_externally(bool tf) { get_state().registry_driven_externally = tf; }

#endif

    //! Diagnostics of the lattice path: {skip, local, full, roll} counts
    long get_number_of_skips() const;
    long get_number_of_local_updates() const;
    long get_number_of_full_updates() const;
    long get_number_of_rolls() const;

    //! Bumped whenever resample() changed the tiles
    unsigned long get_result_generation() const;

    //! Wall time of the last compute phase in seconds (scheduling hint,
    //! diagnostics)
    double get_last_compute_seconds() const {
        return state_ ? state_->last_compute_seconds : 0.0;
    }

#ifndef SWIG
    //! Lattice bookkeeping (created on first use); C++ only
    internal::AVLatticeState &get_state();

    //! The (x, y, z, density) cloud of the current map, cached per result
    const std::vector<IMP::algebra::Vector4D> &get_cloud() const;
    //! Refresh the structure-of-arrays cloud in the state (mean, quadrature)
    void refresh_cloud_soa() const;
#endif

    /**
     * @brief Weighted representative points of the AV cloud on the lattice.
     *
     * The cloud is coarsened into cubic lattice blocks, the smallest block
     * edge that leaves at most `k` non-empty blocks; each block is
     * represented by its weighted centroid and total weight, so the mean
     * position is preserved exactly. (Internally each block also carries
     * its second central moments, which av_distance_quadrature() uses for
     * a second-order correction.) Cached per resample() result.
     * @return flat (x, y, z, w) per point
     */
    std::vector<double> get_quadrature_points(int k = 50) const;

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
    // (r/R0)^6 by three multiplications; std::pow(x, 6.0) costs ~20x more
    double x = distance / forster_radius;
    double x2 = x * x;
    double rda_r0_6 = x2 * x2 * x2;
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
        int distance_type = PROBE_PAIR_DISTANCE_MEAN,
        int n_samples = 10000
);

/**
 * @brief Deterministic distance between two accessible volumes by lattice
 * quadrature (PRD-105).
 *
 * Each cloud is coarsened to at most `quad_k` weighted representative
 * points (AV::get_quadrature_points); the pair quantity is the weighted
 * double sum over the two point sets, each term corrected to second order
 * by the blocks' second central moments. Same number every call,
 * independent of call order. Distance types as in av_distance().
 * Measured on T4L @ 2.0 A, K = 100: max error vs the exact double sum
 * < 0.005 A on the mean distance (MC with 50k samples: ~0.25 A range).
 * @param quad_k maximum number of representative points per cloud
 */
//! The mean-position separation that would reproduce a measurement.
/*!
    The one call an experiment-to-restraint conversion needs: two volumes and
    what was measured between them, in, and the \f$R_{mp}\f$ to restrain to,
    out. It is #IMP::bff::rmp_from_model_distance() over the two volumes'
    point clouds, with the measurement's own distance convention and Förster
    radius, so a caller never has to unpack either.

    \param[in] a,b the two volumes; they are resampled if they are stale
    \param[in] measurement the experimental distance to reproduce
    \param[in] accuracy stop when the model distance is this close, A
    \param[in] n_samples,seed the fixed pair sample the root is found on
    \return the \f$R_{mp}\f$ that reproduces `measurement.distance`
    \throw ValueException when no separation reproduces it

    \see IMP::bff::effective_distance() for the forward, parametric
         approximation of the same relation.
*/
IMPBFFEXPORT double rmp_from_measurement(
        const AV& a, const AV& b, const AVPairDistanceMeasurement& measurement,
        double accuracy = 0.01, int n_samples = 50000, int seed = 0);

//! The flat-bottom well a measurement asks for, in mean-position coordinates.
/*!
    The measurement's value and its asymmetric errors, each converted through
    rmp_from_measurement(): what comes back is
    \f$(R_{mp}^{lo}, R_{mp}^{target}, R_{mp}^{hi})\f$ -- the bounds of the
    interval inside which a model is not penalised, and the centre.

    The three are converted **separately** rather than the errors being carried
    across unchanged, because the \f$R_{mp} \to \langle R_{DA}\rangle\f$
    relation is not linear: an error bar that is symmetric in the measured
    quantity is not symmetric in \f$R_{mp}\f$.

    \param[in] a,b the two volumes
    \param[in] measurement the experimental distance and its errors
    \param[in] accuracy,n_samples,seed as rmp_from_measurement()
    \return three values: low bound, target, high bound. A bound that no
            separation reproduces is clamped to the reachable end rather than
            throwing -- an error bar running past what the volumes allow is a
            one-sided restraint, not an error.
*/
IMPBFFEXPORT std::vector<double> rmp_flat_bottom_bounds(
        const AV& a, const AV& b, const AVPairDistanceMeasurement& measurement,
        double accuracy = 0.01, int n_samples = 50000, int seed = 0);

//! Write a volume to \p path; the format follows the file extension.
/*!
    The one call for getting a volume into a viewer.

    | extension | what is written |
    |---|---|
    | `.xyz` | the point cloud, weight in a fifth column |
    | `.pqr` | the point cloud, weight in the charge column |
    | `.dx` | the density grid, as OpenDX |
    | `.mrc`, `.map`, `.ccp4` | the density grid, through `IMP.em` |

    \param[in] av the volume; it is resampled if it is stale
    \param[in] path the output path, whose extension picks the format
    \throw ValueException on an extension that is not one of the above
    \throw IOException when \p path cannot be opened
*/
IMPBFFEXPORT void write_av(const AV& av, const std::string& path);

//! Fraction of a volume's weight within \p radius of the selected atoms.
/*!
    How much of the volume touches something. \p selection is a
    #IMP::bff::SelectionExpression in either dialect -- `"resi 50-60"` or
    `"resid 50 to 60"` -- evaluated against \p hierarchy.

    \param[in] av the volume
    \param[in] hierarchy where the reference atoms are selected from
    \param[in] selection the selection expression; empty selects every atom
    \param[in] radius the contact radius, A. It has to clear the volume's own
               exclusion: an accessible volume already excludes the structure's
               van der Waals envelope inflated by the dye radius and half the
               linker width, so **no point of it is within about 5 A of any
               atom centre** and a smaller radius returns zero for everything.
               8 A is a reasonable "touching" for a 3.5 A dye.
    \return the weight fraction in [0, 1]
*/
IMPBFFEXPORT double av_overlap(
        const AV& av, const IMP::atom::Hierarchy& hierarchy,
        const std::string& selection = "", double radius = 3.5);

IMPBFFEXPORT double av_distance_quadrature(
        const AV& a,
        const AV& b,
        double forster_radius = 52.0,
        int distance_type = PROBE_PAIR_DISTANCE_MEAN,
        int quad_k = 50
);

// Draw random points in AV. Returns (x,y,z,d) vector
IMPBFFEXPORT void av_random_points(
        const AV& av1,
        double** out_view, int* n_out_view,
        int n_samples=10000
);

//! Random sampling over AV/AV distances
IMPBFFEXPORT void av_random_distances(
        const AV& av1,
        const AV& av2,
        double** out_view, int* n_out_view,
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


// -------- the doors over an IMP::Model (formerly AVBuilder.h) --------

//! Decorate a fresh particle as an AV, resample it, and read the map.
/*!
    The AV is decorated onto its **own** particle with the source passed
    separately. Setting it up on the source particle leaves the resampled map at
    the coordinate origin — header origin (0, 0, 0), obstacles nowhere near the
    search region, and most of the grid reported accessible.

    \param[in] model holds the obstacle hierarchy and the attachment particle
    \param[in] source_particle the attachment atom, a real `XYZR` in that model
    \param[in] linker_length,linker_width,r1,r2,r3,disc_step the AV parameters
    \param[in] allowed_sphere_radius obstacles inside this radius of the
               attachment are ignored; without it the search starts inside the
               attachment atom's own neighbourhood and returns an empty volume,
               reporting nothing
    \param[in] contact_volume_thickness,contact_volume_trapped_fraction the ACV
               split, when one is wanted; see
               #IMP::bff::PathMap::apply_contact_weighting for the rule
    \param[in] search_stencil Dijkstra neighbour stencil. `0` leaves the
               decorator's own default, which is **74** — the LabelLib reference
               metric. `26` is the speed option: ~1.9x faster for ~20 % less
               volume, because a {1,2,3} stencil's isopath surface is cubic
               rather than spherical. `30` is the historical variant.
    \return the volume, carrying its cloud, its grid and the source coordinate
*/
IMPBFFEXPORT AccessibleVolume resample_av(
        IMP::Model* model, IMP::Particle* source_particle,
        double linker_length, double linker_width, double r1, double r2,
        double r3, double disc_step, double allowed_sphere_radius,
        double contact_volume_thickness = 0.0,
        double contact_volume_trapped_fraction = -1.0,
        int search_stencil = 0);

//! An accessible volume from a structure and an fps position definition.
/*!
    The **structure** front door. It applies the FPS strip before measuring
    anything: fps.json positions are calibrated for a structure whose attachment
    residue does not wall in its own dye. A declared `strip_mask` outside the
    fps dialect raises here — loudly, because the alternative is computing
    against obstacles the document said to remove.

    Source clearance scales with the linker: the path search inflates obstacles
    by half the linker width, so the free sphere has to clear that inflation
    plus a grid step of slack or the source tile is walled in and the volume
    comes back empty. A position that declares `allowed_sphere_radius` keeps its
    own value — honoured exactly, never escalated: an empty volume at the
    declared parameters is the answer, not a signal to retry at invented ones.

    \param[in] pdb_path the structure
    \param[in] chain,resseq,atom_name the attachment site
    \param[in] linker_length,linker_width,r1,r2,r3,disc_step the AV parameters
    \param[in] strip_mask an fps `strip_mask`; empty means the default strip
    \param[in] allowed_sphere_radius negative derives it from the linker width
    \param[in] contact_volume_thickness,contact_volume_trapped_fraction the ACV
               split, when one is wanted; see
               #IMP::bff::PathMap::apply_contact_weighting for the rule. **This
               door's uniform-weight convention is suspended when one is
               asked for** -- an accessible contact volume *is* the weighting,
               and flattening it would return a volume that ignored the
               request.
    \throw ValueException when the attachment site is not in the structure
*/
IMPBFFEXPORT AccessibleVolume get_av_from_structure(
        const std::string& pdb_path, const std::string& chain, int resseq,
        const std::string& atom_name, double linker_length, double linker_width,
        double r1, double r2, double r3, double disc_step = 1.5,
        const std::string& strip_mask = "",
        double allowed_sphere_radius = -1.0,
        double contact_volume_thickness = 0.0,
        double contact_volume_trapped_fraction = -1.0);

//! The same front door, from an fps.json `Positions` entry.
/*!
    A position is data, so this takes it as data (a nlohmann JSON value), not
    as separate keywords. The entry's fields -- chain_identifier,
    residue_seq_number, atom_name, linker_length, linker_width, radius1..
    radius3, strip_mask, allowed_sphere_radius, contact_volume_*, and
    simulation_grid_resolution -- are read the way the fps dictionary states
    them, and a declared `simulation_grid_resolution` that disagrees with
    \p disc_step raises: that field is *written into* the particle from
    \p disc_step, so a caller who declares it and omits the step would silently
    build at 1.5 A.

    \param[in] position_json one entry of an fps.json `Positions` section,
               serialised as JSON text
    \param[in] disc_step the voxel spacing, A; `<= 0` derives it from the entry's
               `simulation_grid_resolution` and then refuses a declared one that
               disagrees
    \throw ValueException on a disagreeing resolution or a missing attachment site
*/
IMPBFFEXPORT AccessibleVolume get_av_from_structure(
        const std::string& pdb_path, const std::string& position_json,
        double disc_step = -1.0);

//! Every position of an fps.json `Positions` section, keyed by name.
/*!
    \p pdb_path may be one path or a JSON array of paths; a position's
    `body_id` indexes into it, which is how a docking run gives each rigid body
    its own structure. A position whose attachment atom is not in its structure
    comes back as an **empty** volume rather than one computed at a guessed
    coordinate. \p positions_json is the section serialised as JSON text.

    \return `{name: AccessibleVolume}`, one per key of \p positions_json
*/
IMPBFFEXPORT std::map<std::string, AccessibleVolume> get_avs_for_structure(
        const std::string& positions_json, const std::string& pdb_path,
        double disc_step = -1.0);

IMPBFF_END_NAMESPACE

#endif /* IMPBFF_AV_H */