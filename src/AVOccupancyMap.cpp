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
    for (size_t i = 0; i < xyzr_.size(); i++) {
        IMP::algebra::Vector3D c = coord(i);
        double R = radius(i) + extra_radius_;
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
    auto inside = [&](int x, double dz2dy2) {
        double dx = x * spacing_ - c[0];
        return dz2dy2 + dx * dx < R2;
    };
    for (int z = a[2]; z <= b[2]; z++) {
        double dz = z * spacing_ - c[2];
        double dz2 = dz * dz;
        for (int y = a[1]; y <= b[1]; y++) {
            double dy = y * spacing_ - c[1];
            double dz2dy2 = dz2 + dy * dy;
            double rem = R2 - dz2dy2;
            if (rem <= 0) continue;
            // The x-run inside the sphere on this row: estimated from the
            // chord half-width, then pinned down with the exact per-voxel
            // test at both ends, so the covered set is the same as testing
            // every voxel.
            double half = std::sqrt(rem);
            int xa = std::max(a[0], (int) std::ceil((c[0] - half) / spacing_));
            int xb = std::min(b[0], (int) std::floor((c[0] + half) / spacing_));
            while (xa > a[0] && inside(xa - 1, dz2dy2)) xa--;
            while (xa <= xb && !inside(xa, dz2dy2)) xa++;
            while (xb < b[0] && inside(xb + 1, dz2dy2)) xb++;
            while (xb >= xa && !inside(xb, dz2dy2)) xb--;
            if (xb < xa) continue;
            long row = (z - k0_[2]) * nxy + (y - k0_[1]) * nx - k0_[0];
            int32_t *p = counts_.data() + row + xa;
            for (int x = xa; x <= xb; x++) *p++ += sign;
        }
    }
}

void AVOccupancyMap::record_change(const int lo[3], const int hi[3]) {
    const size_t cap = 256;
    if (changes_.size() >= cap) {
        changes_.erase(changes_.begin(), changes_.begin() + cap / 2);
        oldest_tracked_ = changes_.front().generation;
    }
    ChangeBox b;
    b.generation = generation_;
    for (int d = 0; d < 3; d++) { b.lo[d] = lo[d]; b.hi[d] = hi[d]; }
    changes_.push_back(b);
}

void AVOccupancyMap::record_change_all() {
    int lo[3] = {k0_[0], k0_[1], k0_[2]};
    int hi[3] = {k0_[0] + n_[0] - 1, k0_[1] + n_[1] - 1, k0_[2] + n_[2] - 1};
    record_change(lo, hi);
}

bool AVOccupancyMap::get_changed_since(
        unsigned long generation,
        int kx, int ky, int kz, int nx, int ny, int nz) const {
    if (pending_action_ != 0 && generation <= generation_) {
        // the update classified but not yet applied
        if (pending_box_all_) return true;
        const int lo[3] = {kx, ky, kz};
        const int hi[3] = {kx + nx - 1, ky + ny - 1, kz + nz - 1};
        bool overlap = true;
        for (int d = 0; d < 3; d++) {
            if (pending_hi_[d] < lo[d] || pending_lo_[d] > hi[d]) { overlap = false; break; }
        }
        if (overlap) return true;
    }
    if (generation >= generation_) return false;
    if (generation + 1 < oldest_tracked_) return true;   // history lost
    const int lo[3] = {kx, ky, kz};
    const int hi[3] = {kx + nx - 1, ky + ny - 1, kz + nz - 1};
    for (auto it = changes_.rbegin(); it != changes_.rend(); ++it) {
        if (it->generation <= generation) break;
        bool overlap = true;
        for (int d = 0; d < 3; d++) {
            if (it->hi[d] < lo[d] || it->lo[d] > hi[d]) { overlap = false; break; }
        }
        if (overlap) return true;
    }
    return false;
}

void AVOccupancyMap::full_raster() {
    std::fill(counts_.begin(), counts_.end(), 0);
    int lo[3] = {k0_[0], k0_[1], k0_[2]};
    int hi[3] = {k0_[0] + n_[0] - 1, k0_[1] + n_[1] - 1, k0_[2] + n_[2] - 1};
    last_.resize(xyzr_.size());
    for (size_t i = 0; i < xyzr_.size(); i++) {
        IMP::algebra::Vector3D c = coord(i);
        double r = radius(i);
        last_[i] = IMP::algebra::Vector4D(c[0], c[1], c[2], r);
        add_sphere(c, r + extra_radius_, +1, lo, hi);
    }
    force_full_ = false;
}

int AVOccupancyMap::begin_update(bool force_full) {
    pending_action_ = 0;
    pending_moved_.clear();
    // 1. Extent: grow on demand to the union of requested windows and,
    //    for shared maps, the reach of all atoms.
    if (!fixed_extent_) {
        // The extent only has to cover the atoms' reach: lattice points
        // beyond it hold no count and read_window() returns zero for them,
        // which is the exact occupancy there. Requested windows are honoured
        // only when no atom reach is included (a map without particles).
        int lo[3], hi[3];
        bool have = false;
        if (include_atom_reach_ && !xyzr_.empty()) {
            atom_reach(lo, hi);
            have = true;
        } else if (have_request_) {
            for (int d = 0; d < 3; d++) {
                lo[d] = req_lo_[d];
                hi[d] = req_hi_[d];
            }
            have = true;
        }
        have_request_ = false;
        if (have && !covers(lo, hi)) {
            // Re-centre on the current reach (a full raster follows anyway),
            // rather than growing the union forever as a structure drifts.
            int k0[3], n[3];
            for (int d = 0; d < 3; d++) {
                k0[d] = lo[d] - grow_margin_;
                n[d] = hi[d] - lo[d] + 1 + 2 * grow_margin_;
            }
            set_extent(k0, n);
            n_grow_++;
        }
    }
    if (!have_extent_) return 0;

    // 2. Full raster when required
    if (force_full || force_full_ || last_.size() != xyzr_.size()) {
        pending_action_ = 2;
        pending_box_all_ = true;
        moved_last_ = (long) xyzr_.size();
        moved_total_ += moved_last_;
        return 2;
    }

    // 3. Otherwise classify by the moved set
    for (size_t i = 0; i < xyzr_.size(); i++) {
        IMP::algebra::Vector3D c = coord(i);
        double r = radius(i);
        const IMP::algebra::Vector4D &l = last_[i];
        if (c[0] != l[0] || c[1] != l[1] || c[2] != l[2] || r != l[3]) {
            pending_moved_.push_back(i);
        }
    }
    if (pending_moved_.empty()) {
        n_skip_++;
        return 0;
    }
    moved_last_ = (long) pending_moved_.size();
    moved_total_ += moved_last_;
    // The box the change is confined to: old and new footprints of every
    // moved atom, clipped to the extent. Recorded whichever way the counts
    // are refreshed -- a full raster changes nothing outside it.
    for (int d = 0; d < 3; d++) {
        pending_lo_[d] = std::numeric_limits<int>::max();
        pending_hi_[d] = std::numeric_limits<int>::min();
    }
    for (size_t i : pending_moved_) {
        const IMP::algebra::Vector4D &l = last_[i];
        IMP::algebra::Vector3D c = coord(i);
        double R_old = l[3] + extra_radius_, R_new = radius(i) + extra_radius_;
        for (int d = 0; d < 3; d++) {
            pending_lo_[d] = std::min(pending_lo_[d], (int) std::ceil((std::min(l[d] - R_old, c[d] - R_new)) / spacing_));
            pending_hi_[d] = std::max(pending_hi_[d], (int) std::floor((std::max(l[d] + R_old, c[d] + R_new)) / spacing_));
        }
    }
    for (int d = 0; d < 3; d++) {
        pending_lo_[d] = std::max(pending_lo_[d], k0_[d]);
        pending_hi_[d] = std::min(pending_hi_[d], k0_[d] + n_[d] - 1);
    }
    pending_box_all_ = false;
    // A delta costs two spheres per moved atom, a full raster one per atom.
    pending_action_ = (2 * pending_moved_.size() >= xyzr_.size()) ? 2 : 1;
    return pending_action_;
}

void AVOccupancyMap::raster_slab(int z_lo, int z_hi) {
    // clear the slab, then add every atom clipped to it: bit-identical to
    // full_raster() slab by slab (integer counts, order-free)
    if (z_hi < z_lo) return;
    const long nxy = (long) n_[0] * n_[1];
    std::fill(counts_.begin() + (z_lo - k0_[2]) * nxy,
              counts_.begin() + (z_hi - k0_[2] + 1) * nxy, 0);
    int lo[3] = {k0_[0], k0_[1], z_lo};
    int hi[3] = {k0_[0] + n_[0] - 1, k0_[1] + n_[1] - 1, z_hi};
    for (size_t i = 0; i < xyzr_.size(); i++) {
        IMP::algebra::Vector3D c = coord(i);
        double r = radius(i);
        add_sphere(c, r + extra_radius_, +1, lo, hi);
    }
}

void AVOccupancyMap::apply_local() {
    int lo[3] = {k0_[0], k0_[1], k0_[2]};
    int hi[3] = {k0_[0] + n_[0] - 1, k0_[1] + n_[1] - 1, k0_[2] + n_[2] - 1};
    for (size_t i : pending_moved_) {
        const IMP::algebra::Vector4D &l = last_[i];
        add_sphere(IMP::algebra::Vector3D(l[0], l[1], l[2]),
                   l[3] + extra_radius_, -1, lo, hi);
        IMP::algebra::Vector3D c = coord(i);
        double r = radius(i);
        add_sphere(c, r + extra_radius_, +1, lo, hi);
    }
}

void AVOccupancyMap::end_update() {
    if (pending_action_ == 0) return;
    if (pending_action_ == 2) {
        last_.resize(xyzr_.size());
        for (size_t i = 0; i < xyzr_.size(); i++) {
            IMP::algebra::Vector3D c = coord(i);
            last_[i] = IMP::algebra::Vector4D(c[0], c[1], c[2], radius(i));
        }
        force_full_ = false;
        n_full_++;
    } else {
        for (size_t i : pending_moved_) {
            IMP::algebra::Vector3D c = coord(i);
            last_[i] = IMP::algebra::Vector4D(c[0], c[1], c[2], radius(i));
        }
        n_local_++;
    }
    generation_++;
    if (pending_box_all_) record_change_all();
    else record_change(pending_lo_, pending_hi_);
    pending_action_ = 0;
    pending_moved_.clear();
}

bool AVOccupancyMap::update(bool force_full) {
    int action = begin_update(force_full);
    if (action == 0) return false;
    if (action == 2) {
        raster_slab(k0_[2], k0_[2] + n_[2] - 1);
    } else {
        apply_local();
    }
    end_update();
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


void AVOccupancyMap::read_window_counts(int kx, int ky, int kz,
                                        int nx, int ny, int nz, int32_t *out) const {
    const long wnx = nx, wnxy = (long) nx * ny;
    const long enx = n_[0], enxy = (long) n_[0] * n_[1];
    for (int z = 0; z < nz; z++) {
        int gz = kz + z;
        bool zin = have_extent_ && gz >= k0_[2] && gz < k0_[2] + n_[2];
        for (int y = 0; y < ny; y++) {
            int gy = ky + y;
            bool yin = zin && gy >= k0_[1] && gy < k0_[1] + n_[1];
            int32_t *row = out + z * wnxy + y * wnx;
            if (!yin) {
                std::fill(row, row + nx, 0);
                continue;
            }
            const int32_t *src = counts_.data() + (gz - k0_[2]) * enxy + (gy - k0_[1]) * enx;
            // the in-extent x-run, copied; the rest zero
            int x0 = std::max(0, k0_[0] - kx);
            int x1 = std::min(nx, k0_[0] + n_[0] - kx);
            if (x1 <= x0) { std::fill(row, row + nx, 0); continue; }
            std::fill(row, row + x0, 0);
            std::copy(src + (kx + x0 - k0_[0]), src + (kx + x1 - k0_[0]), row + x0);
            std::fill(row + x1, row + nx, 0);
        }
    }
}

void AVOccupancyMap::read_window_strided(int kx, int ky, int kz,
                                         int nx, int ny, int nz, int stride,
                                         double *out) const {
    const long enx = n_[0], enxy = (long) n_[0] * n_[1];
    long o = 0;
    for (int z = 0; z < nz; z++) {
        int gz = kz + z * stride;
        bool zin = have_extent_ && gz >= k0_[2] && gz < k0_[2] + n_[2];
        for (int y = 0; y < ny; y++) {
            int gy = ky + y * stride;
            bool yin = zin && gy >= k0_[1] && gy < k0_[1] + n_[1];
            for (int x = 0; x < nx; x++, o++) {
                int gx = kx + x * stride;
                if (yin && gx >= k0_[0] && gx < k0_[0] + n_[0]) {
                    out[o] = (double) counts_[(gz - k0_[2]) * enxy + (gy - k0_[1]) * enx + (gx - k0_[0])];
                } else {
                    out[o] = 0.0;
                }
            }
        }
    }
}

AVOccupancyMap *AVOccupancyRegistry::get_map(double spacing, double extra_radius) {
    auto key = std::make_pair(spacing, extra_radius);
    auto it = maps_.find(key);
    if (it == maps_.end()) {
        IMP::Pointer<AVOccupancyMap> m = new AVOccupancyMap(spacing, extra_radius, ps_);
        m->set_was_used(true);
        m->set_coordinate_snapshot(&snapshot_);
        it = maps_.emplace(key, m).first;
    }
    return it->second.get();
}

void AVOccupancyRegistry::refresh_snapshot() {
    IMP::core::XYZRs xyzr(ps_);
    snapshot_.resize(xyzr.size());
    for (size_t i = 0; i < xyzr.size(); i++) {
        IMP::algebra::Vector3D c = xyzr[i].get_coordinates();
        snapshot_[i] = IMP::algebra::Vector4D(c[0], c[1], c[2], xyzr[i].get_radius());
    }
}

AVOccupancyMaps AVOccupancyRegistry::get_maps() const {
    AVOccupancyMaps out;
    for (const auto &kv : maps_) out.push_back(kv.second.get());
    return out;
}

void AVOccupancyRegistry::update_all(bool force_full) {
    refresh_snapshot();
    for (auto &kv : maps_) kv.second->update(force_full);
}

IMPBFF_END_NAMESPACE
