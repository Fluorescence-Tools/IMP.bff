/**
 * \file AVOccupancyMap.cpp
 * \brief The occupancy raster over IMP particles, and the registry that shares one.
 *
 * The raster itself is OccupancyGrid.cpp (core); this file reads the
 * particles into its snapshot and keys shared rasters by particle set.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/AVOccupancyMap.h>
#include <IMP/core/XYZR.h>

IMPBFF_BEGIN_NAMESPACE

AVOccupancyMap::AVOccupancyMap(double spacing, double extra_radius,
                               const IMP::ParticlesTemp &ps, std::string name)
    : OccupancyGrid(spacing, extra_radius,
                    std::shared_ptr<const std::vector<IMP::algebra::Vector4D> >(),
                    name),
      xyzr_(ps),
      own_(std::make_shared<std::vector<IMP::algebra::Vector4D> >()) {
    refresh_from_particles();
    OccupancyGrid::set_coordinate_snapshot(own_);
}

void AVOccupancyMap::refresh_from_particles() {
    std::vector<IMP::algebra::Vector4D> &snap = *own_;
    snap.resize(xyzr_.size());
    for (size_t i = 0; i < xyzr_.size(); i++) {
        const IMP::algebra::Vector3D c = xyzr_[i].get_coordinates();
        snap[i] = IMP::algebra::Vector4D(c[0], c[1], c[2], xyzr_[i].get_radius());
    }
}

int AVOccupancyMap::begin_update(bool force_full) {
    if (own_) refresh_from_particles();
    return OccupancyGrid::begin_update(force_full);
}

AVOccupancyMap *AVOccupancyRegistry::get_map(double spacing, double extra_radius) {
    auto key = std::make_pair(spacing, extra_radius);
    auto it = maps_.find(key);
    if (it == maps_.end()) {
        IMP::Pointer<AVOccupancyMap> m = new AVOccupancyMap(spacing, extra_radius, ps_);
        m->set_was_used(true);
        m->set_coordinate_snapshot(snapshot_);
        if (!obstacle_radii_.empty()) m->set_obstacle_radii(obstacle_radii_);
        it = maps_.emplace(key, m).first;
    }
    return it->second.get();
}

void AVOccupancyRegistry::adopt_obstacle_radii(
        AVRadiiSource source, const std::vector<double> &radii) {
    if (adopted_) {
        if (source == adopted_source_) return;  // O(1): the hot path
        /* Two volumes of one shared raster disagreeing about how big the
           atoms are. Whichever answer were kept, one of them would be scored
           against a raster it did not ask for -- and silently, because the
           volume would still come back full of plausible voxels. See
           IMP::bff::AV::set_radii_source. */
        IMP_THROW("AVOccupancyRegistry: the shared occupancy raster was built "
                  "with the \"" << av_radii_source_to_string(adopted_source_)
                  << "\" radii; a volume sharing it cannot ask for \""
                  << av_radii_source_to_string(source)
                  << "\" (see AV::set_radii_source)", IMP::ValueException);
    }
    if (!radii.empty() && radii.size() != ps_.size()) {
        IMP_THROW("AVOccupancyRegistry: " << radii.size() << " obstacle radii "
                  "for " << ps_.size() << " particles", IMP::ValueException);
    }
    adopted_ = true;
    adopted_source_ = source;
    obstacle_radii_ = radii;
    for (auto &kv : maps_) kv.second->set_obstacle_radii(radii);
}

void AVOccupancyRegistry::refresh_snapshot() {
    IMP::core::XYZRs xyzr(ps_);
    std::vector<IMP::algebra::Vector4D> &snap = *snapshot_;
    snap.resize(xyzr.size());
    for (size_t i = 0; i < xyzr.size(); i++) {
        IMP::algebra::Vector3D c = xyzr[i].get_coordinates();
        snap[i] = IMP::algebra::Vector4D(c[0], c[1], c[2], xyzr[i].get_radius());
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
