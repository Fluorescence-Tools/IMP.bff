/**
 *  \file IMP/bff/internal/AVLatticeState.h
 *  \brief Per-AV bookkeeping for the space-fixed (lattice) evaluation path.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_INTERNAL_AV_LATTICE_STATE_H
#define IMPBFF_INTERNAL_AV_LATTICE_STATE_H

#include <IMP/bff/bff_config.h>
#include <IMP/Pointer.h>
#include <IMP/Particle.h>
#include <IMP/algebra/Vector3D.h>
#include <IMP/algebra/VectorD.h>
#include <IMP/bff/ProbeAccessibleVolumeOccupancyMap.h>
#include <IMP/bff/PathMap.h>

#include <array>
#include <cstdint>
#include <chrono>
#include <vector>

IMPBFF_BEGIN_INTERNAL_NAMESPACE

//! State the lattice path of ProbeAccessibleVolumeDecorator::resample() carries between evaluations.
/** Lives beside the path map on the AV handle (shared between copies of the
    handle) -- it is *not* particle state, so a fresh handle starts cold. */
struct AVLatticeState {
    IMP::ParticlesTemp particles;         //!< the leaves the map was built over (the lattice itself keeps no particle)
    // The window: lattice index of voxel 0 and the (cubic) edge in voxels
    bool have_window = false;
    int k0[3] = {0, 0, 0};
    int n = 0;

    // What the current tiles were computed from
    bool have_result = false;
    IMP::algebra::Vector3D last_source;
    IMP::algebra::VectorD<9> last_parameter;
    unsigned long generation1 = 0;      // occupancy generation, pass 1
    unsigned long generation2 = 0;      // occupancy generation, pass 2
    IMP::algebra::Vector3D last_mean;

    // One warning per handle, not one per frame: the legacy anchoring cannot
    // honour an accessible-contact-volume request (ProbeAccessibleVolumeDecorator::resample_legacy) and a
    // trajectory would otherwise say so thousands of times.
    bool warned_contact_volume = false;
    bool warned_radii_source = false;

    // Bumped whenever the tiles change; consumers cache against it
    unsigned long result_generation = 0;

    // The (x, y, z, density) cloud of the current tiles, computed once per
    // result and shared by the mean position and the quadrature -- as four
    // float arrays (structure of arrays; the values are the float voxel
    // locations and the float density, widened to double where they are
    // used, so sums are bit-identical to the Vector4D form). `cloud` is the
    // Vector4D form, built on demand for external readers.
    std::vector<float> cloud_x, cloud_y, cloud_z, cloud_w;
    std::vector<IMP::algebra::Vector4D> cloud;
    unsigned long cloud_generation = 0;
    bool cloud_valid = false;          // the SoA arrays are current
    unsigned long cloud_aos_generation = 0;
    bool cloud_aos_valid = false;      // `cloud` (Vector4D) is current

    // Set by a restraint that drives the registry maps' updates itself
    // (begin_update before, end_update after the AV prepare phase); prepare
    // then neither updates them nor assumes their counts are current yet.
    bool registry_driven_externally = false;

    // Occupancy sources. `registry` set: shared maps; else private windows.
    IMP::Pointer<ProbeAccessibleVolumeOccupancyRegistry> registry;
    IMP::Pointer<ProbeAccessibleVolumeOccupancyMap> private1;
    IMP::Pointer<ProbeAccessibleVolumeOccupancyMap> private2;

    // Quadrature representation of the cloud (PRD-105 distances):
    // weighted block centroids plus per-block second central moments
    // (xx, yy, zz, xy, xz, yz) for the second-order correction.
    std::vector<IMP::algebra::Vector4D> quad_points;
    std::vector<std::array<double, 6> > quad_moments;
    int quad_k = -1;
    unsigned long quad_generation = 0;
    bool quad_valid = false;

    // A prepared-but-not-yet-computed evaluation (ProbeAccessibleVolumeDecorator::resample_lattice is
    // split into a serial prepare phase that touches the Model and a compute
    // phase that only touches this AV's own map, so a restraint can run the
    // compute phases of its AVs on threads).
    bool pending = false;
    int pending_stage = 0;               // 0: search due, 1: carve due, 2: done
    std::chrono::steady_clock::time_point compute_t0;
    bool pending_shift_xyz = true;
    double pending_ll = 0, pending_allowed = 0;
    //! The attachment atom's own radius, read in prepare().
    /*! The compute phase must not touch the Model, and dropping the
        attachment atom from the obstacle set needs its radius. */
    double pending_source_radius = 0;
    //! The attachment atom's radius **as the raster used it**, cached when the
    //! path map is built. Under `radii_source = "olga"` that is not the radius
    //! the Model carries, and subtracting the wrong sphere would leave a ring
    //! of blocked voxels round the anchor or open ones no rule opened.
    //! Negative = not cached yet (no map, or the map was just rebuilt).
    double source_obstacle_radius = -1.0;
    //! The extra radius each occupancy source was rasterised with, so the
    //! same sphere can be subtracted again: half the linker width for the
    //! linker pass, and one per dye radius for the carve.
    double pending_extra1 = 0, pending_extra2 = 0;
    std::vector<double> pending_extra_dye;
    IMP::algebra::Vector3D pending_source;
    bool pending_set_origin = false;     // window moved: recompute voxel locations
    IMP::algebra::Vector3D pending_grid_origin;
    ProbeAccessibleVolumeOccupancyMap *pending_occ1 = nullptr;
    ProbeAccessibleVolumeOccupancyMap *pending_occ2 = nullptr;
    //! AV3: the extra dye-radius occupancy sources (`pending_occ2` is the
    //! first). Empty for AV1, so that path allocates and reads nothing new.
    std::vector<ProbeAccessibleVolumeOccupancyMap *> pending_occ_dye;
    std::vector<int32_t> window_counts_dye;   //!< AV3 scratch, radii 2..n
    std::vector<IMP::Pointer<ProbeAccessibleVolumeOccupancyMap> > private_dye;
    unsigned long pending_gen1 = 0, pending_gen2 = 0;

    // Coarse search grid (search_grid_factor > 1): a small PathMap on the
    // lattice points whose fine index is a multiple of the factor.
    IMP::Pointer<IMP::bff::PathMap> coarse_map;
    int coarse_k0[3] = {0, 0, 0};    // coarse-lattice index of its voxel 0
    int coarse_n[3] = {0, 0, 0};
    int pending_factor = 1;
    std::vector<double> coarse_data;

    // Scratch: the window's occupancy counts (integer), reused per frame
    std::vector<int32_t> window_counts;

    // Wall time of the last compute phase (for longest-first scheduling)
    double last_compute_seconds = 0.0;

    // Diagnostics
    long n_skip = 0, n_local = 0, n_full = 0, n_roll = 0;
};

IMPBFF_END_INTERNAL_NAMESPACE

#endif /* IMPBFF_INTERNAL_AV_LATTICE_STATE_H */
