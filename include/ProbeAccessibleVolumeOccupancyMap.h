/**
 *  \file IMP/bff/ProbeAccessibleVolumeOccupancyMap.h
 *  \brief Integer occupancy counts of inflated atoms on the global AV lattice.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_PROBEACCESSIBLEVOLUMEOCCUPANCYMAP_H
#define IMPBFF_PROBEACCESSIBLEVOLUMEOCCUPANCYMAP_H

#include <IMP/bff/bff_config.h>

#include <IMP/Object.h>
#include <IMP/Pointer.h>
#include <IMP/Particle.h>
#include <IMP/core/XYZR.h>
#include <IMP/algebra/Vector3D.h>

#include <IMP/bff/VdwRadii.h>
#include <IMP/bff/OccupancyGrid.h>

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! The occupancy raster over IMP particles.
/*! #IMP::bff::OccupancyGrid over the coordinates and radii of `ps`, read
    from the particles at every update. With an external snapshot
    (set_coordinate_snapshot(), which the registry uses) the particles are
    not consulted again: the registry refreshes the snapshot once for every
    map it owns. The name and the Python surface are those of the pre-split
    class; everything but the constructor is the grid's. */
class IMPBFFEXPORT ProbeAccessibleVolumeOccupancyMap : public OccupancyGrid {
    IMP::core::XYZRs xyzr_;
    std::shared_ptr<std::vector<IMP::algebra::Vector4D> > own_;
    void refresh_from_particles();
public:
    ProbeAccessibleVolumeOccupancyMap(double spacing, double extra_radius,
                   const IMP::ParticlesTemp &ps,
                   std::string name = "AVOccupancyMap%1%");
    //! Hand over a snapshot maintained elsewhere; the particles are then not read.
    void set_coordinate_snapshot(std::shared_ptr<const std::vector<IMP::algebra::Vector4D> > snapshot) {
        own_.reset();
        OccupancyGrid::set_coordinate_snapshot(snapshot);
    }
    //! Reads the particles first (unless an external snapshot is set), then classifies the update.
    int begin_update(bool force_full = false) override;
    IMP_OBJECT_METHODS(ProbeAccessibleVolumeOccupancyMap);
};
IMP_OBJECTS(ProbeAccessibleVolumeOccupancyMap, ProbeAccessibleVolumeOccupancyMaps);

//! One AVOccupancyMap per (spacing, extra-radius) class, created on demand.
/** Shared by all AVs of an ProbeNetworkRestraint under `shared_map=True`.
    Occupancy is a function of (atoms, lattice, extra radius) only -- the
    linker length merely masks -- so AVs with the same spacing and the same
    inflation radius read the same raster.
 */
class IMPBFFEXPORT ProbeAccessibleVolumeOccupancyRegistry : public IMP::Object {
    IMP::ParticlesTemp ps_;
    std::vector<double> obstacle_radii_;
    bool adopted_ = false;
    AVRadiiSource adopted_source_ = AV_RADII_IMP;
    std::map<std::pair<double, double>, IMP::Pointer<ProbeAccessibleVolumeOccupancyMap> > maps_;
    std::shared_ptr<std::vector<IMP::algebra::Vector4D> > snapshot_ =
        std::make_shared<std::vector<IMP::algebra::Vector4D> >();
public:
    ProbeAccessibleVolumeOccupancyRegistry(const IMP::ParticlesTemp &ps,
                        std::string name = "AVOccupancyRegistry%1%")
        : IMP::Object(name), ps_(ps) {}

    //! The map of the (spacing, extra_radius) class, created on first use
    ProbeAccessibleVolumeOccupancyMap *get_map(double spacing, double extra_radius);

    //! Adopt one per-particle obstacle-radii vector for every map here.
    /*! A shared raster **is** one obstacle set, so the radii set is a property
        of the registry, not of the volume reading it. The first volume to ask
        fixes it; a second volume asking for a different #IMP::bff::AVRadiiSource
        is a modelling error, not a preference, and throws rather than silently
        getting the first one's answer. (A strip mask already takes its volume
        out of the registry entirely, for the same reason.) The source is
        carried beside the vector so the check on the hot path is one integer
        comparison rather than a walk of every atom.

        Applying it to maps that already exist is safe: `ProbeAccessibleVolumeOccupancyMap`
        remembers the radius it last rasterised each atom with, so a changed
        radius is picked up as an ordinary subtract-old / add-new delta.

        An empty vector means "the particles' own radii", which is what
        `AV_RADII_IMP` -- the default -- with no mask hands over. */
    void adopt_obstacle_radii(AVRadiiSource source,
                              const std::vector<double> &radii);

    //! The radii adopted by adopt_obstacle_radii(); empty = the model's own.
    const std::vector<double> &get_obstacle_radii() const {
        return obstacle_radii_;
    }

    //! The radii set this registry was built with; `imp` until one is adopted.
    AVRadiiSource get_radii_source() const { return adopted_source_; }

    //! All maps created so far
    ProbeAccessibleVolumeOccupancyMaps get_maps() const;

    //! Force a full raster of every map on its next update
    void update_all(bool force_full = false);

    //! Read every particle's (x, y, z, r) from the Model once, for all maps.
    //! Call before begin_update()/update() of the maps in a frame.
    void refresh_snapshot();

    IMP_OBJECT_METHODS(ProbeAccessibleVolumeOccupancyRegistry);
};

IMPBFF_END_NAMESPACE

#endif /* IMPBFF_PROBEACCESSIBLEVOLUMEOCCUPANCYMAP_H */
