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
#include <IMP/algebra/Vector3D.h>
#include <IMP/algebra/VectorD.h>
#include <IMP/bff/AVOccupancyMap.h>
#include <IMP/bff/PathMap.h>

#include <array>
#include <chrono>
#include <vector>

IMPBFF_BEGIN_INTERNAL_NAMESPACE

//! State the lattice path of AV::resample() carries between evaluations.
/** Lives beside the path map on the AV handle (shared between copies of the
    handle) -- it is *not* particle state, so a fresh handle starts cold. */
struct AVLatticeState {
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

    // Bumped whenever the tiles change; consumers cache against it
    unsigned long result_generation = 0;

    // The (x, y, z, density) cloud of the current tiles, computed once per
    // result and shared by the mean position and the quadrature.
    std::vector<IMP::algebra::Vector4D> cloud;
    unsigned long cloud_generation = 0;
    bool cloud_valid = false;

    // Occupancy sources. `registry` set: shared maps; else private windows.
    IMP::Pointer<AVOccupancyRegistry> registry;
    IMP::Pointer<AVOccupancyMap> private1;
    IMP::Pointer<AVOccupancyMap> private2;

    // Quadrature representation of the cloud (PRD-105 distances):
    // weighted block centroids plus per-block second central moments
    // (xx, yy, zz, xy, xz, yz) for the second-order correction.
    std::vector<IMP::algebra::Vector4D> quad_points;
    std::vector<std::array<double, 6> > quad_moments;
    int quad_k = -1;
    unsigned long quad_generation = 0;
    bool quad_valid = false;

    // A prepared-but-not-yet-computed evaluation (AV::resample_lattice is
    // split into a serial prepare phase that touches the Model and a compute
    // phase that only touches this AV's own map, so a restraint can run the
    // compute phases of its AVs on threads).
    bool pending = false;
    int pending_stage = 0;               // 0: search due, 1: carve due, 2: done
    std::chrono::steady_clock::time_point compute_t0;
    bool pending_shift_xyz = true;
    double pending_ll = 0, pending_allowed = 0;
    IMP::algebra::Vector3D pending_source;
    bool pending_set_origin = false;     // window moved: recompute voxel locations
    IMP::algebra::Vector3D pending_grid_origin;
    AVOccupancyMap *pending_occ1 = nullptr;
    AVOccupancyMap *pending_occ2 = nullptr;
    unsigned long pending_gen1 = 0, pending_gen2 = 0;

    // Coarse search grid (search_grid_factor > 1): a small PathMap on the
    // lattice points whose fine index is a multiple of the factor.
    IMP::Pointer<IMP::bff::PathMap> coarse_map;
    int coarse_k0[3] = {0, 0, 0};    // coarse-lattice index of its voxel 0
    int coarse_n[3] = {0, 0, 0};
    int pending_factor = 1;
    std::vector<double> coarse_data;

    // Wall time of the last compute phase (for longest-first scheduling)
    double last_compute_seconds = 0.0;

    // Diagnostics
    long n_skip = 0, n_local = 0, n_full = 0, n_roll = 0;
};

IMPBFF_END_INTERNAL_NAMESPACE

#endif /* IMPBFF_INTERNAL_AV_LATTICE_STATE_H */
