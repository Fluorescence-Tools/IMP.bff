/**
 *  \file IMP/bff/AVOccupancyMap.h
 *  \brief Integer occupancy counts of inflated atoms on the global AV lattice.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_AVOCCUPANCYMAP_H
#define IMPBFF_AVOCCUPANCYMAP_H

#include <IMP/bff/bff_config.h>

#include <IMP/Object.h>
#include <IMP/Pointer.h>
#include <IMP/Particle.h>
#include <IMP/core/XYZR.h>
#include <IMP/algebra/Vector3D.h>

#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Occupancy counts of atoms, inflated by one extra radius, on the AV lattice.
/** The AV lattice (PRD-105) is absolute: voxel centres sit at integer
    multiples of the grid spacing, anchored at the world origin. This class
    keeps, for a rectangular extent of that lattice, the number of inflated
    atoms covering every lattice point -- exactly what `PathMap::sample_obstacles`
    computes with the BINARIZED_SPHERE kernel, but expressed in global integer
    lattice indices so that two extents agree bit-for-bit on any lattice point
    they share.

    Because the values are integer covering counts, moved atoms are applied as
    subtract-old / add-new deltas and remain bit-exact against a full
    re-rasterisation. `update()` classifies each refresh as a *skip* (nothing
    moved), a *local* delta, or a *full* raster (first build, extent change,
    forced, or too many atoms moved for a delta to pay).

    Windows are read out with `read_window()`: lattice points inside the extent
    return their count, points outside return zero. This is the "one shared
    occupancy raster per (spacing, extra-radius) class" of PRD-105 when the
    extent covers several AV windows, and the private per-AV raster when the
    extent is a single window.
 */
class IMPBFFEXPORT AVOccupancyMap : public IMP::Object {

    double spacing_;
    double extra_radius_;
    IMP::core::XYZRs xyzr_;

    // (x, y, z, radius) of every particle as last rasterised; empty until the
    // first full raster.
    std::vector<IMP::algebra::Vector4D> last_;

    // Current extent: lattice index of the first voxel and voxel counts.
    int k0_[3] = {0, 0, 0};
    int n_[3] = {0, 0, 0};
    bool have_extent_ = false;
    std::vector<int32_t> counts_;

    // Extent requested since the last update (union of windows).
    int req_lo_[3] = {0, 0, 0};
    int req_hi_[3] = {0, 0, 0};
    bool have_request_ = false;
    bool include_atom_reach_ = true;
    int grow_margin_ = 3;
    bool fixed_extent_ = false;
    bool force_full_ = true;

    unsigned long generation_ = 0;

    // Diagnostics
    long n_skip_ = 0, n_local_ = 0, n_full_ = 0, n_grow_ = 0, n_roll_ = 0;
    long moved_last_ = 0, moved_total_ = 0;

    void add_sphere(const IMP::algebra::Vector3D &c, double radius, int sign,
                    const int lo[3], const int hi[3]);
    void full_raster();
    void set_extent(const int k0[3], const int n[3]);
    bool covers(const int lo[3], const int hi[3]) const;
    void atom_reach(int lo[3], int hi[3]) const;

public:

    /**
     * @param spacing lattice spacing
     * @param extra_radius radius added to every particle radius before
     *        rasterisation (half the linker width, or the dye radius)
     * @param ps particles (XYZR) that act as obstacles
     */
    AVOccupancyMap(double spacing, double extra_radius,
                   const IMP::ParticlesTemp &ps,
                   std::string name = "AVOccupancyMap%1%");

    double get_spacing() const { return spacing_; }
    double get_extra_radius() const { return extra_radius_; }

    //! Ask the map to cover the window [k0, k0 + n) before the next update
    void request_window(int kx, int ky, int kz, int nx, int ny, int nz);

    //! Fix the extent to exactly one window (private per-AV raster).
    /** A window that differs from the current one by a whole-voxel shift is a
        *roll*: the retained block is shifted in place and only the newly
        exposed slab is rasterised on the next update.
     */
    void set_window(int kx, int ky, int kz, int nx, int ny, int nz);

    //! Include the reach of all atoms in the extent (shared maps, default on)
    void set_include_atom_reach(bool tf) { include_atom_reach_ = tf; }

    //! Voxels of margin added around the required extent when growing
    void set_grow_margin(int m) { grow_margin_ = m; }

    //! Bring the counts up to date with the current particle coordinates.
    /** Returns true if anything changed (generation advanced). */
    bool update(bool force_full = false);

    //! Copy the window [k0, k0 + n) into `out` (row-major x fastest, like
    //! IMP::em::DensityMap). Lattice points outside the extent read as 0.
    void read_window(int kx, int ky, int kz, int nx, int ny, int nz,
                     double *out) const;

    //! read_window() into a fresh vector (Python-friendly)
    std::vector<double> get_window(int kx, int ky, int kz,
                                   int nx, int ny, int nz) const {
        std::vector<double> out((size_t) nx * ny * nz);
        read_window(kx, ky, kz, nx, ny, nz, out.data());
        return out;
    }

    //! Increments on every update that changed at least one count
    unsigned long get_generation() const { return generation_; }

    //! Diagnostics
    long get_number_of_skips() const { return n_skip_; }
    long get_number_of_local_updates() const { return n_local_; }
    long get_number_of_full_updates() const { return n_full_; }
    long get_number_of_grows() const { return n_grow_; }
    long get_number_of_rolls() const { return n_roll_; }
    //! Particles that had moved at the last update that changed a count
    long get_number_of_moved_last() const { return moved_last_; }
    long get_number_of_moved_total() const { return moved_total_; }
    std::vector<int> get_extent() const {
        return {k0_[0], k0_[1], k0_[2], n_[0], n_[1], n_[2]};
    }
    long get_number_of_voxels() const {
        return (long) n_[0] * n_[1] * n_[2];
    }

    IMP_OBJECT_METHODS(AVOccupancyMap);
};

IMP_OBJECTS(AVOccupancyMap, AVOccupancyMaps);


//! One AVOccupancyMap per (spacing, extra-radius) class, created on demand.
/** Shared by all AVs of an AVNetworkRestraint under `shared_map=True`.
    Occupancy is a function of (atoms, lattice, extra radius) only -- the
    linker length merely masks -- so AVs with the same spacing and the same
    inflation radius read the same raster.
 */
class IMPBFFEXPORT AVOccupancyRegistry : public IMP::Object {
    IMP::ParticlesTemp ps_;
    std::map<std::pair<double, double>, IMP::Pointer<AVOccupancyMap> > maps_;
public:
    AVOccupancyRegistry(const IMP::ParticlesTemp &ps,
                        std::string name = "AVOccupancyRegistry%1%")
        : IMP::Object(name), ps_(ps) {}

    //! The map of the (spacing, extra_radius) class, created on first use
    AVOccupancyMap *get_map(double spacing, double extra_radius);

    //! All maps created so far
    AVOccupancyMaps get_maps() const;

    //! Force a full raster of every map on its next update
    void update_all(bool force_full = false);

    IMP_OBJECT_METHODS(AVOccupancyRegistry);
};

IMPBFF_END_NAMESPACE

#endif /* IMPBFF_AVOCCUPANCYMAP_H */
