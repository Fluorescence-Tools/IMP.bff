/**
 *  \file IMP/bff/PathMapTile.h
 *  \brief Tile used in path search by PathMap
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_PATHMAPTILE_H
#define IMPBFF_PATHMAPTILE_H

#include <IMP/bff/bff_config.h>

#include <cmath>  /* std::sqrt */
#include <vector>
#include <utility> /* std::pair */
#include <algorithm>

#include <IMP/bff/PathMap.h>
#include <IMP/bff/PathMapTileEdge.h>

IMPBFF_BEGIN_NAMESPACE

const bool  TILE_VISITED_DEFAULT    = false;
const float TILE_PENALTY_DEFAULT    = 100000.0f;
const float TILE_COST_DEFAULT       = 100000.0f;
const float TILE_EDGE_COST_DEFAULT  = 100000.0f;

const float TILE_PENALTY_THRESHOLD  = 100000.0f;
const float TILE_OBSTACLE_THRESHOLD = 0.000001f;
const float TILE_OBSTACLE_PENALTY   = 100000.0f;


/// Value types that can be read from a PathMapTile
typedef enum{
    PM_TILE_PENALTY,             /// Write path penalty
    PM_TILE_COST,                /// Write cost
    PM_TILE_DENSITY,             /// Density of tile
    PM_TILE_COST_DENSITY,        /// Threshold path length and write tile weights
    PM_TILE_PATH_LENGTH,         /// Write path length
    PM_TILE_PATH_LENGTH_DENSITY, /// Threshold path length and write tile weights
    PM_TILE_FEATURE,             /// Threshold path length and write tile weights
    PM_TILE_ACCESSIBLE_DENSITY,  /// Density that is accessible (Path length in bounds)
    PM_TILE_ACCESSIBLE_FEATURE   /// Feature that is accessible (Path length in bounds)
} PathMapTileOutputs;


class PathMap;

class IMPBFFEXPORT PathMapTile{

friend class PathMap;

private:

    //! Compute edges of tiles.
    /*! Fill the edges tiles that correspond to voxels
     * in an AccessibleVolume.
     *
     * A edge is a neighboring tile. The neighborhood is
     * defined by a 3D box around a Tile. A tile that has
     * a visit penalty that exceeds a threshold value is
     * not added to the neighbor / edge list.
     *
     * @param av AccessibleVolume
     * @param tiles List of tiles with empty edges
     * @param neighbor_radius neighboring tiles closer than neighbor_radius
     * are connected by edges.
     * @param tile_penalty_threshold Tiles with a visit penalty
     * larger than this threshold are not added to the edge list
     */
    void update_edges(
            IMP::bff::PathMap* av,
            std::vector<PathMapTile> &tiles,
            double neighbor_radius,
            float tile_penalty_threshold = TILE_PENALTY_THRESHOLD
    );

    /**
      * @brief Updates the edges of the tile.
      * @param nx The x-coordinate of the tile.
      * @param ny The y-coordinate of the tile.
      * @param nz The z-coordinate of the tile.
      * @param tiles The vector of PathMapTile objects.
      * @param neighbor_idxs The vector of neighbor indices.
      * @param tile_penalty_threshold The tile penalty threshold.
      */
    void update_edges_2(
            int nx, int ny, int nz,
            std::vector<PathMapTile>& tiles,
            std::vector<int> neighbor_idxs,
            float tile_penalty_threshold =TILE_PENALTY_THRESHOLD
    );


protected:

    long idx;                  // tile index: corresponds to voxel index

    // Variable for path search
    float penalty;                  // penalty for visiting tile
    float cost;                     // total cost for visiting tile
    PathMapTile* previous; // tile previously visited in path search

    // Additional tile feature (e.g. av density)
    std::map<std::string, float> features;

    // Edges leaving tile and going to neighboring tiles
    std::vector<PathMapTileEdge> edges;

public:

    /// AV density
    float density;


    //! Construct an accessible volume tile
    /*!
     * An accessible volume (AV) tile relates to a voxel in
     * an AV. A set of interconnected tiles (neighboring tiles)
     * is used to compute an optimal path from the labeling site
     * to all other positions the the AV. A path is a sequence of
     * tiles. The cost of a path the the sum of all costs (associated
     * to tiles and edges connecting tiles). Visiting a tile in
     * a path adds to the cost of a path. The visiting penalty is
     * defined when constructing a tile.
     *
     * @param index Identifier of the tile (corresponds to index of voxel)
     * @param visit_penalty Penalty for visiting (used for implementing
     * obstacles)
     * @param tile_density Additional information of tile (can be used to
     * implement weighted AVs)
     */
    explicit PathMapTile(
            long index=-1,
            float visit_penalty = 0.0,
            float tile_density = 1.0
    ) :
            idx(index),
            penalty(visit_penalty),
            cost(std::numeric_limits<float>::max()),
            previous(nullptr),
            density(tile_density)
    {}


    /**
     * @brief Computes the path from a tile to the origin.
     * @return A vector of long integers representing the path.
     */
    std::vector backtrack_to_path();

   /**
    * @brief Get the value of a tile.
    *
    * A tile in an accessible volume contains information on the penalty for visiting a tile,
    * the cost of a path from the origin of a path search to the tile, the density of the tile,
    * and other user-defined information.
    *
    * When getting information from a tile, the returned values can be cropped to a range.
    *
    * @param value_type Specifies the type of the returned information (see: PathMapTileOutputs).
    * Depending on the value type, the output can be the penalty for visiting the tile,
    * the total cost of a path to the tile, or the density of the tile.
    * Additionally, user-defined content can be accessed.
    * @param bounds Bound for cropping the output values.
    * @param feature_name Name of a feature (when accessing additional information).
    * @param grid_spacing Spacing between the tiles (important to specify when accessing path length).
    * @return Value of the tile for the specified parameters.
    */
    float get_value(
            int value_type,
            std::pair<float, float> bounds = std::pair<float, float>(
                            {std::numeric_limits<float>::min(),
                             std::numeric_limits<float>::max()}),
            const std::string &feature_name="",
            float grid_spacing = 1.0
    );

    //! Set the value of a tile
    /*!
     * Sets the value of a tile.
     * @param value_type Type of the value
     * @param value value that will be written
     * @param name name of the value (only used for user-defined tile features)
     */
    void set_value(int value_type, float value, const std::string &name="");

};


IMPBFF_END_NAMESPACE

#endif //IMPBFF_PATHMAPTILE_H