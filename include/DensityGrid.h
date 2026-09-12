/**
 *  \file IMP/bff/DensityGrid.h
 *  \brief A regular grid with an origin and a spacing, and spheres sampled
 *         into it.
 *
 * `PathMap` derived from `IMP::em::SampledDensityMap` and used four things of
 * it: a `double` array over a regular lattice, the per-voxel coordinate
 * caches, a list of spheres, and `resample()` to mark the voxels those spheres
 * cover. None of that is electron microscopy. It cost this module a dependency
 * on `IMP.em` -- and through the module graph `statistics` as well -- for a
 * path search over a lattice.
 *
 * So the lattice is written out here. The member names (`data_`, `header_`,
 * `x_loc_`) and the method names are IMP's on purpose: `PathMap.cpp` is
 * fifteen hundred lines that call them, and a base class that answers to the
 * same names turns a rewrite into a swap. The names are the seam, not the
 * design.
 *
 * **The sampling is IMP's `BINARIZED_SPHERE`, reproduced exactly** -- a voxel
 * takes 1 when its centre is strictly inside a sphere, and the value is
 * accumulated, not assigned, so two overlapping spheres leave a 2 behind.
 * Nothing here reads that number except through a threshold, but it is what
 * was written before and the AV densities are checked against the ones IMP
 * produced. The Gaussian and soft-sphere kernels IMP also offers are *not*
 * reproduced: `PathMap` never asked for them (`PathMap.h` defaults to
 * `BINARIZED_SPHERE`), and a kernel nobody calls is a kernel nobody tests.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#ifndef IMPBFF_DENSITYGRID_H
#define IMPBFF_DENSITYGRID_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/IMPCompatibility.h>

#include <IMP/Object.h>
#include <IMP/algebra/Vector3D.h>

#include <cereal/access.hpp>

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! One sphere: where it is and how big it is.
/*! What `IMP::core::XYZR` was used for here, without the particle it hung on.
    The setters exist because `PathMap::sample_obstacles` inflates every radius
    by the linker half-width, samples, and puts the radii back. */
class IMPBFFEXPORT GridSphere {
    friend class cereal::access;
    IMP::algebra::Vector3D c_;
    double r_;

 public:
    GridSphere() : c_(0, 0, 0), r_(0) {}
    GridSphere(const IMP::algebra::Vector3D& c, double r) : c_(c), r_(r) {}
    const IMP::algebra::Vector3D& get_coordinates() const { return c_; }
    void set_coordinates(const IMP::algebra::Vector3D& c) { c_ = c; }
    double get_radius() const { return r_; }
    void set_radius(double r) { r_ = r; }
    template <class Archive> void serialize(Archive &ar) { ar(c_, r_); }

    IMP_SHOWABLE_INLINE(GridSphere,
                        out << "GridSphere(r=" << r_ << ")");
};
IMP_VALUES(GridSphere, GridSpheres);

//! The shape of a grid: extent, spacing and where its corner sits.
/*! A stand-in for the handful of `IMP::em::DensityHeader` fields this module
    ever touched. Everything an EM header also carries -- the cell angles, the
    map statistics, the axis order -- is absent because nothing here read it. */
class IMPBFFEXPORT GridHeader {
    friend class cereal::access;
    int nx_, ny_, nz_;
    // float, not double, and that is not an oversight. IMP's DensityHeader
    // holds these as float and returns them as float, so every voxel location
    // downstream was computed at float precision. Holding them as double here
    // and rounding on the way out is *more* accurate and therefore wrong: it
    // moved AV mean positions in the eighth significant figure, which the
    // lattice pins catch and no amount of staring at the formula explains.
    float spacing_;
    float xorigin_, yorigin_, zorigin_;
    float xtop_, ytop_, ztop_;
    float resolution_;
    bool top_calculated_;

 public:
    GridHeader()
        : nx_(0), ny_(0), nz_(0), spacing_(1.0), xorigin_(0), yorigin_(0),
          zorigin_(0), xtop_(0), ytop_(0), ztop_(0), resolution_(1.0),
          top_calculated_(false) {}

    int get_nx() const { return nx_; }
    int get_ny() const { return ny_; }
    int get_nz() const { return nz_; }

    float get_spacing() const { return spacing_; }
    void set_spacing(float s) { spacing_ = s; top_calculated_ = false; }

    float get_xorigin() const { return xorigin_; }
    float get_yorigin() const { return yorigin_; }
    float get_zorigin() const { return zorigin_; }
    void set_xorigin(float v) { xorigin_ = v; top_calculated_ = false; }
    void set_yorigin(float v) { yorigin_ = v; top_calculated_ = false; }
    void set_zorigin(float v) { zorigin_ = v; top_calculated_ = false; }

    float get_xtop() const { return xtop_; }
    float get_ytop() const { return ytop_; }
    float get_ztop() const { return ztop_; }

    float get_resolution() const { return resolution_; }
    void set_resolution(float r) { resolution_ = r; }

    //! Recompute the far corner. \p force redoes it even if it is current.
    /*! IMP's signature, including the flag: `PathMap` passes `true` after
        moving the origin, and the caching is why that flag has to exist. */
    void compute_xyz_top(bool force = false) {
        if (top_calculated_ && !force) return;
        xtop_ = xorigin_ + spacing_ * nx_;
        ytop_ = yorigin_ + spacing_ * ny_;
        ztop_ = zorigin_ + spacing_ * nz_;
        top_calculated_ = true;
    }

    void update_map_dimensions(int nnx, int nny, int nnz) {
        nx_ = nnx;
        ny_ = nny;
        nz_ = nnz;
        top_calculated_ = false;
        compute_xyz_top();
    }

    //! Cereal, because `PathMapHeader` holds one of these and serializes it.
    /*! Eleven numbers. The `IMP::em::DensityHeader` this replaced wrote about
        forty -- cell angles, axis order, microscope voltage, defocus -- none
        of which a path-search lattice ever had a value for. That does mean a
        `PathMapHeader` pickled by an older build will not load into this one:
        the archive is a different shape, and it is not versioned. */
    template <class Archive> void serialize(Archive &ar) {
        ar(nx_, ny_, nz_, spacing_, xorigin_, yorigin_, zorigin_, xtop_, ytop_,
           ztop_, resolution_, top_calculated_);
    }

    IMP_SHOWABLE_INLINE(GridHeader,
                        out << "GridHeader(" << nx_ << "x" << ny_ << "x"
                            << nz_ << " @ " << spacing_ << ")");
};
IMP_VALUES(GridHeader, GridHeaders);

//! A regular grid of doubles, with spheres sampled into it.
class IMPBFFEXPORT DensityGrid : public IMP::Object {
 protected:
    GridHeader header_;
    std::vector<double> data_;
    //! Per-voxel coordinates, cached. `float` because IMP's were.
    std::unique_ptr<float[]> x_loc_, y_loc_, z_loc_;
    bool loc_calculated_;
    //! How many voxels the location arrays were actually allocated for.
    /*! `loc_calculated_` alone is not enough: it says the caches were built,
        not that they were built for *this* shape. A grid that grows without
        clearing the flag leaves arrays that are too short, and the reads that
        follow run off the end -- which shows up as a heap corruption
        somewhere else entirely, not as a bad index here. */
    long loc_n_ = 0;
    //! The obstacles. IMP kept these as particles; they are values here.
    /*! `GridSpheres`, the `IMP_VALUES` plural, and not `std::vector`: SWIG
        matches typemaps by the exact spelling of the type, and the converter
        `IMP_SWIG_VALUE(GridSphere, GridSpheres)` installs is for the typedef.
        Declared as `std::vector<GridSphere>` the same method wrapped, and
        refused every Python list handed to it. */
    GridSpheres xyzr_;

 public:
    explicit DensityGrid(std::string name = "DensityGrid")
        : IMP::Object(name), loc_calculated_(false) {}

    const IMP::bff::GridHeader* get_header() const { return &header_; }
    IMP::bff::GridHeader* get_header_writable() { return &header_; }

    long get_number_of_voxels() const {
        return static_cast<long>(header_.get_nx()) * header_.get_ny() *
               header_.get_nz();
    }

    float get_spacing() const { return header_.get_spacing(); }

    IMP::algebra::Vector3D get_origin() const {
        return IMP::algebra::Vector3D(header_.get_xorigin(),
                                      header_.get_yorigin(),
                                      header_.get_zorigin());
    }

    //! Move the corner, and rebuild the coordinate cache immediately.
    /*! Invalidating without rebuilding is not enough, and the difference is a
        crash rather than a wrong number: the caches are `unique_ptr`s, so a
        caller that reads a location after moving the origin dereferences
        null. IMP's `DensityMap::set_origin` recomputes here for exactly that
        reason, and the AV lattice fast path depends on it -- it moves the
        origin and reads locations without asking for a recompute in between. */
    void set_origin(const IMP::algebra::Vector3D& o) {
        header_.set_xorigin(o[0]);
        header_.set_yorigin(o[1]);
        header_.set_zorigin(o[2]);
        header_.compute_xyz_top(true);
        reset_all_voxel2loc();
        calc_all_voxel2loc();
    }

    void resize(long nvox) {
        data_.assign(static_cast<std::size_t>(nvox), 0.0);
        reset_all_voxel2loc();
    }

    void reset_data(double v = 0.0) {
        std::fill(data_.begin(), data_.end(), v);
    }

    double get_value(long i) const { return data_[i]; }
    void set_value(long i, double v) { data_[i] = v; }
    double* get_data() { return data_.data(); }
    const double* get_data() const { return data_.data(); }

    void reset_all_voxel2loc() { loc_calculated_ = false; loc_n_ = 0; }

    //! Fill the per-voxel coordinate caches. Idempotent, as IMP's was.
    void calc_all_voxel2loc() {
        if (loc_calculated_) return;
        const long nvox = get_number_of_voxels();
        x_loc_.reset(new float[nvox]);
        y_loc_.reset(new float[nvox]);
        z_loc_.reset(new float[nvox]);
        loc_n_ = nvox;
        const float sp = header_.get_spacing();
        const float ox = header_.get_xorigin(), oy = header_.get_yorigin(),
                    oz = header_.get_zorigin();
        const int nx = header_.get_nx(), ny = header_.get_ny();
        int ix = 0, iy = 0, iz = 0;
        for (long ii = 0; ii < nvox; ++ii) {
            x_loc_[ii] = ix * sp + ox;
            y_loc_[ii] = iy * sp + oy;
            z_loc_[ii] = iz * sp + oz;
            if (++ix == nx) {
                ix = 0;
                if (++iy == ny) {
                    iy = 0;
                    ++iz;
                }
            }
        }
        loc_calculated_ = true;
    }

    double get_location_in_dim_by_voxel(long i, int dim) const {
        return dim == 0 ? x_loc_[i] : (dim == 1 ? y_loc_[i] : z_loc_[i]);
    }

    IMP::algebra::Vector3D get_location_by_voxel(long i) const {
        return IMP::algebra::Vector3D(x_loc_[i], y_loc_[i], z_loc_[i]);
    }

    //! The voxel index of an (x, y, z) lattice position.
    inline long xyz_ind2voxel(int x, int y, int z) const {
        return static_cast<long>(z) * header_.get_nx() * header_.get_ny() +
               static_cast<long>(y) * header_.get_nx() + x;
    }

    //! Which lattice plane a coordinate falls in, along one axis.
    /*! The `0.5 +` before the floor is IMP's and is load-bearing: it rounds to
        the *nearest* plane rather than the one below, so a coordinate sitting
        a hair under a voxel centre still belongs to that voxel. */
    int get_dim_index_by_location(float loc, int dim) const {
        // The inputs are float, as IMP's are, and the arithmetic is double,
        // as IMP's is: `DensityMap::get_dim_index_by_location` subtracts a
        // float `loc_val` from `get_origin()[ind]`, which is a *double* built
        // from the float origin, and divides by the float spacing -- every
        // step promoted. Doing the subtraction and division in float instead
        // rounds twice more, and at a half-voxel boundary that is a different
        // voxel: the legacy lattice search seeds itself with
        // `get_voxel_by_location(source)`, so a different voxel there is a
        // different path and an AV tens of angstroms away.
        const double o = dim == 0 ? header_.get_xorigin()
                                  : (dim == 1 ? header_.get_yorigin()
                                              : header_.get_zorigin());
        const double sp = header_.get_spacing();
        return static_cast<int>(std::floor(0.5 + (static_cast<double>(loc) - o) / sp));
    }

    long get_voxel_by_location(float x, float y, float z) const {
        return xyz_ind2voxel(get_dim_index_by_location(x, 0),
                             get_dim_index_by_location(y, 1),
                             get_dim_index_by_location(z, 2));
    }

    long get_voxel_by_location(const IMP::algebra::Vector3D& v) const {
        return get_voxel_by_location(v[0], v[1], v[2]);
    }

    int get_dim_index_by_voxel(long index, int dim) const {
        const int nx = header_.get_nx(), ny = header_.get_ny();
        if (dim == 0) return static_cast<int>(index % nx);
        if (dim == 1) return static_cast<int>((index / nx) % ny);
        return static_cast<int>(index / (static_cast<long>(nx) * ny));
    }

    const GridSpheres& get_spheres() const { return xyzr_; }
    void set_spheres(const GridSpheres& s) { xyzr_ = s; }

    //! Mark every voxel whose centre falls inside a sphere.
    /*! IMP's `BINARIZED_SPHERE` path, reproduced: zero the map, then for each
        sphere walk only the voxels its own bounding box covers and **add** one
        where the centre is strictly inside. Strictly, and accumulating, both
        because IMP was -- a voxel exactly on the surface stays out, and an
        overlap counts twice. */
    void resample() {
        reset_data();
        calc_all_voxel2loc();
        const long nvox_check = get_number_of_voxels();
        IMP_USAGE_CHECK(static_cast<long>(data_.size()) == nvox_check,
                        "DensityGrid::resample: data is " << data_.size()
                        << " voxels, the header says " << nvox_check);
        IMP_USAGE_CHECK(loc_n_ == nvox_check,
                        "DensityGrid::resample: location caches are for "
                        << loc_n_ << " voxels, the header says " << nvox_check);
        const int nx = header_.get_nx(), ny = header_.get_ny(),
                  nz = header_.get_nz();
        const long nxny = static_cast<long>(nx) * ny;
        const float sp = header_.get_spacing();
        const float ox = header_.get_xorigin(), oy = header_.get_yorigin(),
                    oz = header_.get_zorigin();
        for (std::size_t ii = 0; ii < xyzr_.size(); ++ii) {
            const IMP::algebra::Vector3D c = xyzr_[ii].get_coordinates();
            const double r = xyzr_[ii].get_radius();
            if (r <= 0.0) continue;  // a zero radius blocks nothing
            const double r2 = r * r;
            int iminx, iminy, iminz, imaxx, imaxy, imaxz;
            local_bounding_box(c, r, sp, ox, oy, oz, nx, ny, nz, &iminx,
                               &iminy, &iminz, &imaxx, &imaxy, &imaxz);
            for (int iz = iminz; iz <= imaxz; ++iz) {
                const long znxny = iz * nxny;
                for (int iy = iminy; iy <= imaxy; ++iy) {
                    long ivox = znxny + static_cast<long>(iy) * nx + iminx;
                    for (int ix = iminx; ix <= imaxx; ++ix, ++ivox) {
                        const double dx = x_loc_[ivox] - c[0];
                        const double dy = y_loc_[ivox] - c[1];
                        const double dz = z_loc_[ivox] - c[2];
                        if (dx * dx + dy * dy + dz * dz < r2) data_[ivox] += 1.0;
                    }
                }
            }
        }
    }

    //! This grid's data as an MRC map (write_mrc()).
    void write_mrc(const std::string& path) const;

    IMP_OBJECT_METHODS(DensityGrid);

 private:
    //! IMP's `DensityMap::lower_voxel_shift`, verbatim in its arithmetic.
    /*! `kdist` is a **float** because IMP's is: `calc_local_bounding_box`
        declares it so and the sphere's double radius is rounded on the way
        in. Both ends clamp to `[0, ndim-1]`, and both use `floor` -- the
        upper bound is not a `ceil`. None of that changes which voxels end up
        marked (the sphere test inside the box is what decides), but it is the
        same code path as before rather than an argument that it is
        equivalent, and the argument is the thing that goes stale. */
    static int lower_shift(double loc, float kdist, float orig, float sp,
                           int ndim) {
        int i = static_cast<int>(std::floor((loc - kdist - orig) / sp));
        if (i < 0) i = 0;
        if (i > ndim - 1) i = ndim - 1;
        return i;
    }
    static int upper_shift(double loc, float kdist, float orig, float sp,
                           int ndim) {
        int i = static_cast<int>(std::floor((loc + kdist - orig) / sp));
        if (i < 0) i = 0;
        if (i > ndim - 1) i = ndim - 1;
        return i;
    }

    //! The index range a sphere can reach, clamped to the grid.
    static void local_bounding_box(const IMP::algebra::Vector3D& c, double r,
                                   float sp, float ox, float oy, float oz,
                                   int nx, int ny, int nz, int* iminx,
                                   int* iminy, int* iminz, int* imaxx,
                                   int* imaxy, int* imaxz) {
        const float kdist = static_cast<float>(r);
        *iminx = lower_shift(c[0], kdist, ox, sp, nx);
        *iminy = lower_shift(c[1], kdist, oy, sp, ny);
        *iminz = lower_shift(c[2], kdist, oz, sp, nz);
        *imaxx = upper_shift(c[0], kdist, ox, sp, nx);
        *imaxy = upper_shift(c[1], kdist, oy, sp, ny);
        *imaxz = upper_shift(c[2], kdist, oz, sp, nz);
    }
};

// -------- the MRC writer --------
//! Write a lattice's values as an MRC2014 map.
/*! A 1024-byte MRC2014 header and the values as 32-bit floats, x fastest
    (columns = x, rows = y, sections = z, so `mapc/mapr/maps` = 1/2/3). The
    header says what the map is: mode 2, the grid and cell sizes with the
    cell lengths `n * spacing` and 90-degree angles, zero starts with the
    origin in the ORIGIN fields (Angstrom), `ispg` 1 (a volume), NVERSION
    20140, the density minimum, maximum, mean and RMS computed from the
    values, the "MAP " tag, a machine stamp naming this machine's byte order
    (the bytes are written in host order), and one label. `mrcfile`,
    ChimeraX and IMP::em read it. The previous writer copied what IMP::em's
    writer produced for this module's maps -- NaN statistics, `ispg` set to
    the bit pattern of 1.0f, cell lengths that ignored the spacing -- and
    that was IMP's, not correct; this one is correct (owner, 2026-09-07).
    \param[in] values `n` values, `n` = nx * ny * nz, x fastest */
IMPBFFEXPORT void write_mrc(const std::string& path, const GridHeader& header,
                            const float* values, std::size_t n);

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_DENSITYGRID_H
