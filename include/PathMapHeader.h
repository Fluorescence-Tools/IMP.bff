/**
 *  \file IMP/bff/PathMapHeader.h
 *  \brief Header class for path search class PathMap
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_PATHMAPHEADER_H
#define IMPBFF_PATHMAPHEADER_H

#include <IMP/bff/bff_config.h>

#include <cmath> /* ceil */

#include <IMP/Particle.h>
#include <IMP/algebra/Vector3D.h>
#include <IMP/atom/Atom.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/Selection.h>
#include <IMP/em/DensityHeader.h>

#include <IMP/bff/internal/json.h>
#include <IMP/bff/AV.h>

#include <algorithm>

IMPBFF_BEGIN_NAMESPACE


class PathMap;

class IMPBFFEXPORT PathMapHeader {

friend class IMP::bff::PathMap;

private:

    double grid_spacing_;
    double max_path_length_;
    double neighbor_radius_;
    double obstacle_threshold_;

    IMP::em::DensityHeader density_header_;

protected:

    IMP::algebra::Vector3D path_origin_;

public:

    ~PathMapHeader() = default;

    /*!
     *
     * @param max_path_length maximum length of path (defines also size of grid)
     * @param grid_spacing spacing between grid tiles
     * @param neighbor_radius defines size of box around tile where other
     * voxels are considered a neighbor
     * @param obstacle_threshold voxels with density larger than this
     * threshold value are considered an obstacle.
     */
    PathMapHeader(
            double max_path_length = 10.0,
            double grid_spacing = 1.0,
            double neighbor_radius = 2,
            double obstacle_threshold = std::numeric_limits<double>::epsilon()
    );

    //! Update the dimensions of the AV to be (nnx,nny,nnz)
    //! The origin of the map does not change. If not values
    //! are provided used linker length & radius to update.
    /**
        \param[in] nnx the new number of voxels on the X axis
        \param[in] nny the new number of voxels on the Y axis
        \param[in] nnz the new number of voxels on the Z axis
     */
    void update_map_dimensions(int nx=-1, int ny=-1, int nz=-1);

    /**
     * @brief Sets the origin of the path.
     * @param v The origin vector.
     */
    void set_path_origin(const IMP::algebra::Vector3D &v);

    //! Returns position of the labeling site
    IMP::algebra::Vector3D get_path_origin() const {
        return path_origin_;
    }

    /**
     * @brief Get the maximum linker/path length from origin.
     * @return The maximum linker/path length.
     */
    double get_max_path_length(){
        return max_path_length_;
    }

    /**
     * @brief Get the simulation grid resolution.
     * @return The simulation grid resolution as a double.
     */
    double get_simulation_grid_resolution();

    /**
     * @brief Set the obstacle threshold.
     * @param obstacle_threshold The obstacle threshold value
     */
    void set_obstacle_threshold(double obstacle_threshold);

    /**
     * @brief Get the obstacle threshold.
     * @return The obstacle threshold value
     */
    double get_obstacle_threshold() const{
        return obstacle_threshold_;
    }

    /**
     * @brief Set the neighbor radius.
     * @param neighbor_radius The neighbor radius value
     */
    void set_neighbor_radius(double neighbor_radius);

    /**
     * @brief Get the neighbor radius.
     * @return The neighbor radius as a double.
     */
    double get_neighbor_radius() const{
        return neighbor_radius_;
    }

    /**
     * @brief Get the size of the neighbor box.
     * @return The size of the neighbor box
     */
    int get_neighbor_box_size() const;

    //! Returns a read-only pointer to the header of the map
    const IMP::em::DensityHeader *get_density_header() const {
        return &density_header_; }

    //! Returns a pointer to the header of the map in a writable version
    IMP::em::DensityHeader *get_density_header_writable() {
        return &density_header_; }

    //! Get origin on the PathMap (the corner of the grid)
    IMP::algebra::Vector3D get_origin() const ;

    /**
     * @brief Get the edge length of the grid.
     * @return The edge length of the grid as a double.
     */
    double get_grid_edge_length();

    //! Set origin on the PathMap (the corner of the grid)
    void set_origin(float x, float y, float z);

};

IMP_OBJECTS(PathMapHeader, PathMapHeaders);

IMPBFF_END_NAMESPACE


#endif //IMPBFF_PATHMAPHEADER_H
