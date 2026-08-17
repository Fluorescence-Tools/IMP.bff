/**
 *  \file AVOccupancyMap.cpp
 *  \brief Integer occupancy counts of inflated atoms on the global AV lattice.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#include <IMP/bff/AVOccupancyMap.h>

#include <algorithm>
#include <cmath>
#include <limits>

IMPBFF_BEGIN_NAMESPACE

AVOccupancyMap::AVOccupancyMap(
        double spacing, double extra_radius,
        const IMP::ParticlesTemp &ps, std::string name
) : IMP::Object(name), spacing_(spacing), extra_radius_(extra_radius),
    xyzr_(ps) {
    IMP_USAGE_CHECK(spacing > 0, "AVOccupancyMap: spacing must be positive");
}

void AVOccupancyMap::atom_reach(int lo[3], int hi[3]) const {
    for (int d = 0; d < 3; d++) {
        lo[d] = std::numeric_limits<int>::max();
        hi[d] = std::numeric_limits<int>::min();
    }
    for (const auto &p : xyzr_) {
        IMP::algebra::Vector3D c = p.get_coordinates();
        double R = p.get_radius() + extra_radius_;
        for (int d = 0; d < 3; d++) {
            int a = (int) std::floor((c[d] - R) / spacing_);
            int b = (int) std::ceil((c[d] + R) / spacing_);
            lo[d] = std::min(lo[d], a);
            hi[d] = std::max(hi[d], b);
        }
    }
}

bool AVOccupancyMap::covers(const int lo[3], const int hi[3]) const {
    if (!have_extent_) return false;
    for (int d = 0; d < 3; d++) {
        if (lo[d] < k0_[d]) return false;
        if (hi[d] > k0_[d] + n_[d] - 1) return false;
    }
    return true;
}

void AVOccupancyMap::set_extent(const int k0[3], const int n[3]) {
    for (int d = 0; d < 3; d++) {
        k0_[d] = k0[d];
        n_[d] = n[d];
    }
    counts_.assign((size_t) n_[0] * n_[1] * n_[2], 0);
    have_extent_ = true;
    force_full_ = true;
}

void AVOccupancyMap::request_window(int kx, int ky, int kz,
                                    int nx, int ny, int nz) {
    IMP_USAGE_CHECK(!fixed_extent_,
                    "AVOccupancyMap: request_window on a fixed-window map");
    int lo[3] = {kx, ky, kz};
    int hi[3] = {kx + nx - 1, ky + ny - 1, kz + nz - 1};
    if (!have_request_) {
        for (int d = 0; d < 3; d++) {
            req_lo_[d] = lo[d];
            req_hi_[d] = hi[d];
        }
        have_request_ = true;
    } else {
        for (int d = 0; d < 3; d++) {
            req_lo_[d] = std::min(req_lo_[d], lo[d]);
            req_hi_[d] = std::max(req_hi_[d], hi[d]);
        }
    }
}

void AVOccupancyMap::set_window(int kx, int ky, int kz,
                                int nx, int ny, int nz) {
    fixed_extent_ = true;
    include_atom_reach_ = false;
    int k0[3] = {kx, ky, kz};
    int n[3] = {nx, ny, nz};

    bool same_shape = have_extent_ &&
        n[0] == n_[0] && n[1] == n_[1] && n[2] == n_[2];
    bool same_origin = same_shape &&
        k0[0] == k0_[0] && k0[1] == k0_[1] && k0[2] == k0_[2];
    if (same_origin) return;

    // A first window, a reshaped window, or a window whose counts were
    // never rasterised: start from scratch.
    if (!same_shape || force_full_ || last_.size() != xyzr_.size()) {
        set_extent(k0, n);
        return;
    }

    // Integer-voxel roll: keep the block the old and the new window share,
    // rasterise the exposed remainder at the coordinates the retained block
    // was rasterised at, so the whole window again equals one full raster at
    // `last_` and the ordinary moved-atom delta applies afterwards.
    int ilo[3], ihi[3];
    bool overlap = true;
    for (int d = 0; d < 3; d++) {
        ilo[d] = std::max(k0[d], k0_[d]);
        ihi[d] = std::min(k0[d] + n[d], k0_[d] + n_[d]) - 1;
        if (ihi[d] < ilo[d]) overlap = false;
    }
    if (!overlap) {
        set_extent(k0, n);
        return;
    }
    n_roll_++;

    std::vector<int32_t> fresh((size_t) n[0] * n[1] * n[2], 0);
    const long onx = n_[0], onxy = (long) n_[0] * n_[1];
    const long nnx = n[0], nnxy = (long) n[0] * n[1];
    for (int z = ilo[2]; z <= ihi[2]; z++) {
        for (int y = ilo[1]; y <= ihi[1]; y++) {
            long osrc = (z - k0_[2]) * onxy + (y - k0_[1]) * onx + (ilo[0] - k0_[0]);
            long ndst = (z - k0[2]) * nnxy + (y - k0[1]) * nnx + (ilo[0] - k0[0]);
            std::copy(counts_.begin() + osrc,
                      counts_.begin() + osrc + (ihi[0] - ilo[0] + 1),
                      fresh.begin() + ndst);
        }
    }
    counts_.swap(fresh);
    int old_k0[3] = {k0_[0], k0_[1], k0_[2]};
    int old_n[3] = {n_[0], n_[1], n_[2]};
    for (int d = 0; d < 3; d++) {
        k0_[d] = k0[d];
        n_[d] = n[d];
    }

    // The exposed region is the new window minus the retained block: up to
    // one slab per axis, each a box.
    //   X slab: x outside old, all y, all z of the new window
    //   Y slab: x inside old,  y outside old, all z
    //   Z slab: x inside old,  y inside old,  z outside old
    auto raster_box = [&](const int lo[3], const int hi[3]) {
        for (int d = 0; d < 3; d++) if (hi[d] < lo[d]) return;
        for (size_t i = 0; i < xyzr_.size(); i++) {
            const IMP::algebra::Vector4D &l = last_[i];
            IMP::algebra::Vector3D c(l[0], l[1], l[2]);
            add_sphere(c, l[3] + extra_radius_, +1, lo, hi);
        }
    };
    int new_lo[3] = {k0[0], k0[1], k0[2]};
    int new_hi[3] = {k0[0] + n[0] - 1, k0[1] + n[1] - 1, k0[2] + n[2] - 1};
    int old_lo[3] = {old_k0[0], old_k0[1], old_k0[2]};
    int old_hi[3] = {old_k0[0] + old_n[0] - 1, old_k0[1] + old_n[1] - 1,
                     old_k0[2] + old_n[2] - 1};
    for (int axis = 0; axis < 3; axis++) {
        // two candidate intervals outside the old window along `axis`:
        // below old_lo and above old_hi
        int lo[3], hi[3];
        for (int d = 0; d < 3; d++) {
            if (d < axis) {          // inside old on lower axes
                lo[d] = std::max(new_lo[d], old_lo[d]);
                hi[d] = std::min(new_hi[d], old_hi[d]);
            } else {                 // all of new on higher axes
                lo[d] = new_lo[d];
                hi[d] = new_hi[d];
            }
        }
        // below
        lo[axis] = new_lo[axis];
        hi[axis] = std::min(new_hi[axis], old_lo[axis] - 1);
        raster_box(lo, hi);
        // above
        lo[axis] = std::max(new_lo[axis], old_hi[axis] + 1);
        hi[axis] = new_hi[axis];
        raster_box(lo, hi);
    }
}

void AVOccupancyMap::add_sphere(
        const IMP::algebra::Vector3D &c, double radius, int sign,
        const int lo[3], const int hi[3]) {
    if (radius <= 0) return;
    const double R2 = radius * radius;
    int a[3], b[3];
    for (int d = 0; d < 3; d++) {
        a[d] = std::max(lo[d], (int) std::ceil((c[d] - radius) / spacing_));
        b[d] = std::min(hi[d], (int) std::floor((c[d] + radius) / spacing_));
        if (b[d] < a[d]) return;
    }
    const long nx = n_[0], nxy = (long) n_[0] * n_[1];
    for (int z = a[2]; z <= b[2]; z++) {
        double dz = z * spacing_ - c[2];
        double dz2 = dz * dz;
        for (int y = a[1]; y <= b[1]; y++) {
            double dy = y * spacing_ - c[1];
            double dz2dy2 = dz2 + dy * dy;
            long row = (z - k0_[2]) * nxy + (y - k0_[1]) * nx - k0_[0];
            for (int x = a[0]; x <= b[0]; x++) {
                double dx = x * spacing_ - c[0];
                if (dz2dy2 + dx * dx < R2) {
                    counts_[row + x] += sign;
                }
            }
        }
    }
}

void AVOccupancyMap::full_raster() {
    std::fill(counts_.begin(), counts_.end(), 0);
    int lo[3] = {k0_[0], k0_[1], k0_[2]};
    int hi[3] = {k0_[0] + n_[0] - 1, k0_[1] + n_[1] - 1, k0_[2] + n_[2] - 1};
    last_.resize(xyzr_.size());
    for (size_t i = 0; i < xyzr_.size(); i++) {
        IMP::algebra::Vector3D c = xyzr_[i].get_coordinates();
        double r = xyzr_[i].get_radius();
        last_[i] = IMP::algebra::Vector4D(c[0], c[1], c[2], r);
        add_sphere(c, r + extra_radius_, +1, lo, hi);
    }
    force_full_ = false;
}

bool AVOccupancyMap::update(bool force_full) {
    // 1. Extent: grow on demand to the union of requested windows and,
    //    for shared maps, the reach of all atoms.
    if (!fixed_extent_) {
        int lo[3], hi[3];
        bool have = false;
        if (have_request_) {
            for (int d = 0; d < 3; d++) {
                lo[d] = req_lo_[d];
                hi[d] = req_hi_[d];
            }
            have = true;
        }
        if (include_atom_reach_ && !xyzr_.empty()) {
            int alo[3], ahi[3];
            atom_reach(alo, ahi);
            if (!have) {
                for (int d = 0; d < 3; d++) { lo[d] = alo[d]; hi[d] = ahi[d]; }
                have = true;
            } else {
                for (int d = 0; d < 3; d++) {
                    lo[d] = std::min(lo[d], alo[d]);
                    hi[d] = std::max(hi[d], ahi[d]);
                }
            }
        }
        have_request_ = false;
        if (have && !covers(lo, hi)) {
            int k0[3], n[3];
            for (int d = 0; d < 3; d++) {
                // never shrink: keep what is covered already
                int cur_lo = have_extent_ ? std::min(lo[d], k0_[d]) : lo[d];
                int cur_hi = have_extent_ ? std::max(hi[d], k0_[d] + n_[d] - 1) : hi[d];
                k0[d] = cur_lo - grow_margin_;
                n[d] = cur_hi - cur_lo + 1 + 2 * grow_margin_;
            }
            set_extent(k0, n);
            n_grow_++;
        }
    }
    if (!have_extent_) return false;

    // 2. Full raster when required
    if (force_full || force_full_ || last_.size() != xyzr_.size()) {
        full_raster();
        n_full_++;
        moved_last_ = (long) xyzr_.size();
        moved_total_ += moved_last_;
        generation_++;
        return true;
    }

    // 3. Otherwise classify by the moved set
    std::vector<size_t> moved;
    for (size_t i = 0; i < xyzr_.size(); i++) {
        IMP::algebra::Vector3D c = xyzr_[i].get_coordinates();
        double r = xyzr_[i].get_radius();
        const IMP::algebra::Vector4D &l = last_[i];
        if (c[0] != l[0] || c[1] != l[1] || c[2] != l[2] || r != l[3]) {
            moved.push_back(i);
        }
    }
    if (moved.empty()) {
        n_skip_++;
        return false;
    }
    moved_last_ = (long) moved.size();
    moved_total_ += moved_last_;
    // A delta costs two spheres per moved atom, a full raster one per atom.
    if (2 * moved.size() >= xyzr_.size()) {
        full_raster();
        n_full_++;
        generation_++;
        return true;
    }
    int lo[3] = {k0_[0], k0_[1], k0_[2]};
    int hi[3] = {k0_[0] + n_[0] - 1, k0_[1] + n_[1] - 1, k0_[2] + n_[2] - 1};
    for (size_t i : moved) {
        const IMP::algebra::Vector4D &l = last_[i];
        add_sphere(IMP::algebra::Vector3D(l[0], l[1], l[2]),
                   l[3] + extra_radius_, -1, lo, hi);
        IMP::algebra::Vector3D c = xyzr_[i].get_coordinates();
        double r = xyzr_[i].get_radius();
        add_sphere(c, r + extra_radius_, +1, lo, hi);
        last_[i] = IMP::algebra::Vector4D(c[0], c[1], c[2], r);
    }
    n_local_++;
    generation_++;
    return true;
}

void AVOccupancyMap::read_window(int kx, int ky, int kz,
                                 int nx, int ny, int nz, double *out) const {
    const long wnx = nx, wnxy = (long) nx * ny;
    const long enx = n_[0], enxy = (long) n_[0] * n_[1];
    for (int z = 0; z < nz; z++) {
        int gz = kz + z;
        bool zin = have_extent_ && gz >= k0_[2] && gz < k0_[2] + n_[2];
        for (int y = 0; y < ny; y++) {
            int gy = ky + y;
            bool yin = zin && gy >= k0_[1] && gy < k0_[1] + n_[1];
            double *row = out + z * wnxy + y * wnx;
            if (!yin) {
                std::fill(row, row + nx, 0.0);
                continue;
            }
            const int32_t *src = counts_.data()
                + (gz - k0_[2]) * enxy + (gy - k0_[1]) * enx;
            for (int x = 0; x < nx; x++) {
                int gx = kx + x;
                row[x] = (gx >= k0_[0] && gx < k0_[0] + n_[0])
                         ? (double) src[gx - k0_[0]] : 0.0;
            }
        }
    }
}


IMPBFF_END_NAMESPACE
