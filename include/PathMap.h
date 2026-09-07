#ifndef IMPBFF_PATHMAP_H
#define IMPBFF_PATHMAP_H

/**
 *  \file IMP/bff/PathMap.h
 *  \brief Class to search path on grids
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2023 IMP Inventors. All rights reserved.
 *
 */

// -------- from PathMapHeader.h --------
/**
 *  (formerly IMP/bff/PathMapHeader.h, now a section of this file)
 *  \brief Header class for path search class PathMap
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */
#include <IMP/bff/bff_config.h>

#include <IMP/bff/Base.h>

#include <cmath> /* ceil */

#include <IMP/algebra/Vector3D.h>
#include <functional>

#include <cereal/access.hpp>
#include <IMP/bff/DensityGrid.h>

#include <IMP/bff/internal/json.h>

#include <algorithm>

IMPBFF_BEGIN_NAMESPACE

class AV;  // the decorator (connection layer) is a friend; a name is all a friend needs


class PathMap;

class IMPBFFEXPORT PathMapHeader {

friend class IMP::bff::PathMap;
friend class cereal::access;

    template<class Archive> void serialize(Archive &ar) {
        ar(grid_spacing_, max_path_length_, neighbor_radius_,
           obstacle_threshold_, density_header_, path_origin_);
    }

private:

    // mutable: get_simulation_grid_resolution() is a const getter that
    // refreshes this cache from density_header_ before returning it.
    mutable double grid_spacing_;
    double max_path_length_;
    double neighbor_radius_;
    double obstacle_threshold_;

    //! The grid this header describes.
    /*! Was an `IMP::em::DensityHeader`. The half-dozen fields that were ever
        read of it -- extent, spacing, origin -- are what `GridHeader` holds,
        and none of them is about electron microscopy. */
    GridHeader density_header_;

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

    /**
     * @brief Sets the labeling-site position and the grid origin separately.
     *
     * The lattice (space-fixed) path of the AV anchors the grid corner on
     * the global lattice rather than at `v - edge/2`.
     * @param v The labeling site
     * @param grid_origin The location of voxel (0, 0, 0)
     */
    void set_path_origin(const IMP::algebra::Vector3D &v,
                         const IMP::algebra::Vector3D &grid_origin);

    //! Returns position of the labeling site
    IMP::algebra::Vector3D get_path_origin() const {
        return path_origin_;
    }

    /**
     * @brief Get the maximum linker/path length from origin.
     * @return The maximum linker/path length.
     */
    double get_max_path_length() const {
        return max_path_length_;
    }

    /**
     * @brief Get the simulation grid resolution.
     * @return The simulation grid resolution as a double.
     */
    double get_simulation_grid_resolution() const;

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
    const IMP::bff::GridHeader *get_density_header() const {
        return &density_header_; }

    //! Returns a pointer to the header of the map in a writable version
    IMP::bff::GridHeader *get_density_header_writable() {
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


    IMP_SHOWABLE_INLINE(PathMapHeader,
                        out << "PathMapHeader(grid_spacing=" << grid_spacing_
                            << ", max_path_length=" << max_path_length_ << ")");
};

IMP_VALUES(PathMapHeader, PathMapHeaders);

// PathMapHeader is a value, not an IMP::Object -- it derives from nothing and
// is copied by value everywhere. IMP_OBJECTS declared its plural as a vector
// of ref-counted pointers, which was never right and collided the moment the
// class was declared to SWIG as the value it is.

IMPBFF_END_NAMESPACE

// -------- from PathMapTileEdge.h --------
/**
 *  (formerly IMP/bff/PathMapTileEdge.h, now a section of this file)
 *  \brief Tile edges used in path search by PathMap
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */


#include <vector>



IMPBFF_BEGIN_NAMESPACE

class PathMap;
class PathMapTile;


class PathMapTileEdge{

friend class PathMapTile;
friend class PathMap;
friend class cereal::access;

    template<class Archive> void serialize(Archive &ar) {
        ar(tile_idx, length);
    }


protected:

    int tile_idx; /// The tile the edge is pointing to
    float length; /// the path length / cost of going to the tile

public:

    /// Length of an edge (usually cartesian distance between tiles)
    float get_length() const{
        return length;
    }

    /*!
     *
     * @param edge_target
     * @param edge_cost
     */
    PathMapTileEdge(
            int edge_target = -1,
            float edge_cost = std::numeric_limits<float>::max()
    ) :
            tile_idx(edge_target), length(edge_cost){}

    IMP_SHOWABLE_INLINE(PathMapTileEdge,
                        out << "PathMapTileEdge(tile=" << tile_idx
                            << ", length=" << length << ")");
};

IMP_VALUES(PathMapTileEdge, PathMapTileEdges);


IMPBFF_END_NAMESPACE

// -------- from PathMapTile.h --------
/**
 *  (formerly IMP/bff/PathMapTile.h, now a section of this file)
 *  \brief Tile used in path search by PathMap
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */


#include <cmath>  /* std::sqrt */
#include <utility> /* std::pair */


#include <cereal/types/map.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

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
friend class cereal::access;

    /* `previous` is deliberately not archived: it points at another tile
       inside the same map and is transient path-search state, rebuilt by
       find_path_dijkstra()/find_path_astar(). Archiving a bare intra-container
       pointer would deep-copy the chain and hand each tile its own duplicate. */
    template<class Archive> void save(Archive &ar) const {
        ar(idx, penalty, cost, features, edges, density);
    }

    template<class Archive> void load(Archive &ar) {
        ar(idx, penalty, cost, features, edges, density);
        previous = nullptr;
    }

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
            const std::vector<int> &neighbor_idxs,
            float tile_penalty_threshold =TILE_PENALTY_THRESHOLD
    );


protected:

    long idx;                  // tile index: corresponds to voxel index

    // Variable for path search
    float penalty;                  // penalty for visiting tile in a path search
    float cost;                     // total cost for visiting tile in a path search (integrated cost)
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
    std::vector<long> backtrack_to_path();

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


    IMP_SHOWABLE_INLINE(PathMapTile,
                        out << "PathMapTile(idx=" << idx << ", penalty=" << penalty
                            << ", cost=" << cost << ", density=" << density << ")");
};

IMP_VALUES(PathMapTile, PathMapTiles);


IMPBFF_END_NAMESPACE

// -------- from PathMap.h --------

#include <stdlib.h>     /* malloc, free, rand */
#include <limits>
#include <cstdint>

#include <cmath>
#include <unordered_set>
#include <queue>
#include <utility>  /* std::pair */
#include <Eigen/Dense>

#include <IMP/Object.h>



namespace IMP { namespace em { class DensityMap; } }

IMPBFF_BEGIN_NAMESPACE

class PathMapTile;


//! Weight per linker path length: how a flexible linker distributes its dye.
/*!
    An accessible volume treats every reachable voxel as equally likely. A real
    linker does not: a random coil rarely sits fully collapsed against its
    attachment point, so voxels a short path from the source are less likely
    than the count of them suggests. This is that correction -- a tabulated
    \f$P(\ell)\f$ over the **path length** \f$\ell\f$, the geodesic through
    free space that #IMP::bff::PathMap already computes per voxel.

    The uniform weighting (the default) returns 1 everywhere, which is exactly
    the unweighted volume, so nothing changes unless a table is set.

    \note Out-of-range path lengths are **clamped** to the table's end values.
          The reference implementation returns NaN there, which a weight cannot
          be: one NaN voxel makes every mean position, distance and efficiency
          computed from the volume NaN too.
*/
class IMPBFFEXPORT LinkerWeighting {
    double x_min_, x_max_;
    std::vector<double> y_;      //!< empty means uniform

public:
    //! The uniform weighting: 1 at every path length.
    LinkerWeighting() : x_min_(0.0), x_max_(0.0) {}

    //! \param[in] x_min,x_max the path-length range the table spans, A
    /*! \param[in] y the weights, evenly spaced over `[x_min, x_max]`
        \throw ValueException when \p y has fewer than two entries or the range
               is empty */
    LinkerWeighting(double x_min, double x_max, const std::vector<double>& y);

    //! True when this is the uniform weighting and costs nothing to apply.
    bool get_is_uniform() const { return y_.empty(); }

    //! The weight at \p path_length, linearly interpolated and clamped.
    double get_weight(double path_length) const;

    //! Fraction of the table's total weight at path lengths at most \p x.
    /*! A linker can only reach `linker_length`, so this is how much of the
        tabulated distribution the volume actually samples. When it is
        vanishing the weighting is not a correction to the volume -- it
        selects the volume's outer shell and discards the rest. */
    double get_supported_fraction(double x) const;

    double get_x_min() const { return x_min_; }
    double get_x_max() const { return x_max_; }
    //! The tabulated weights.
    const std::vector<double>& get_weights() const { return y_; }

    IMP_SHOWABLE_INLINE(LinkerWeighting,
                        out << "LinkerWeighting("
                            << (y_.empty() ? std::string("uniform")
                                           : std::to_string(y_.size())
                                                     + " points)"));
};
IMP_VALUES(LinkerWeighting, LinkerWeightings);

//! Read a chain-weighting table and pick the column for \p linker_length.
/*!
    The file is tab-separated with a header row: its first field is a label and
    the rest are the table's axis keys. Each later row is a path length in A
    followed by one weight per key. The column chosen is the one with the
    **first key not less than** \p linker_length; a linker longer than every
    key gets the uniform weighting rather than an extrapolation.

    Lines beginning with `#` are comments.

    \param[in] path the table file
    \param[in] linker_length the linker this volume uses, A
    \throw IOException when \p path cannot be read
    \throw ValueException when the file is not in the format above
*/
IMPBFFEXPORT LinkerWeighting read_linker_weighting(const std::string& path,
                                                   double linker_length);

//! The shipped chain-weighting table, for \p linker_length.
/*! `data/linker/chain_weighting.csv`, read once and cached. */
IMPBFFEXPORT LinkerWeighting linker_weighting(double linker_length);

class IMPBFFEXPORT PathMap : public DensityGrid {

friend class PathMapTile;
friend class AV;

private:

    // used in path search
    std::vector<char>  visited;   // byte per tile: a scan is a byte load, not a bit test
    std::vector<bool>  edge_computed;
    std::vector<float> cost;

    // True while `visited` marks the tiles reached by the last search on the
    // current tiles; only they can carry an accessible density, so
    // get_xyz_density() skips the rest (invalidated by
    // update_tiles/set_data/resize).
    bool reached_valid_ = false;

    // The search loop of find_path(), generic over the frontier comparator so
    // that plain Dijkstra does not pay for a std::function heuristic call in
    // every heap comparison.
    template<class Cmp>
    void find_path_impl(long path_begin_idx, long path_end_idx, Cmp cmp);

    //! The linker weight of voxel \p i, from its path length. 1 if uniform.
    /*! Applied where the density is *read* rather than folded into
        `density_soa_`: the raster is rebuilt and the tile copies are taken at
        different points of a multi-stage resample, and a fold that has to
        happen between them is a fold that will one day happen on the wrong
        side of one. This cannot get out of step. */
    inline float path_weight(long i) const {
        if (weighting_.get_is_uniform()) return 1.0f;
        const float c = cost[i];
        if (!(c >= 0.0f && c < TILE_COST_DEFAULT)) return 0.0f;
        return (float) weighting_.get_weight(
                c * pathMapHeader_.get_simulation_grid_resolution());
    }

    // Per-tile "interior" flags for the current shape: an interior tile is
    // at least the neighbour box away from every face and needs no bounds
    // check in the search. Rebuilt when the shape or the box changes.
    std::vector<char> interior_;
    int interior_nx_ = -1, interior_ny_ = -1, interior_nz_ = -1, interior_box_ = -1;
    const char *get_interior_flags();

    // The lattice search proper on the compact arrays: obstacles are
    // encoded in `cost` (BLOCKED_COST for blocked tiles, TILE_COST_DEFAULT
    // for open ones) so the relaxation reads one array; penalties are
    // binary there. Buckets and the queued stamps are reused across calls.
    static constexpr float BLOCKED_COST = -1.0f;
    std::vector<std::vector<int> > bucket_scratch_;
    std::vector<int16_t> queued_scratch_;
    // the stencil in the form the relaxation loop reads (built with offsets_)
    std::vector<long> nb_delta_;
    std::vector<float> nb_len_;
    std::vector<int> nb_dz_, nb_dy_, nb_dx_;
    void dijkstra_lattice(long source_idx, float max_cost);

    // Candidate tiles of the lattice search: those that can lie inside the
    // block sphere around any point of the source voxel, in index order,
    // with their (ix, iy, iz); every other tile is blocked by the cost
    // template. Cached per shape / source voxel / radius (the AV window is
    // centred on the source voxel, so this is one build per AV).
    // kind 1: the tile is inside the block sphere and outside the open sphere
    // for every source position within the source voxel -- only the
    // occupancy decides; kind 0: on one of the two shells, needs the test.
    struct BallCandidate { int idx; int16_t ix, iy, iz; int16_t kind; };
    std::vector<BallCandidate> ball_cand_;
    // Runs of consecutive candidates along x (idx0, ix0, iy, iz, length):
    // contiguous in memory, so the passes over them vectorise.
    struct BallRun { int idx0; int16_t ix0, iy, iz, len; };
    std::vector<BallRun> ball_runs_;
    std::vector<float> ball_cost_template_;
    long ballc_source_ = -1; double ballc_radius_ = -1, ballc_open_ = -1;
    int ballc_nx_ = -1, ballc_ny_ = -1, ballc_nz_ = -1;
    const std::vector<BallCandidate> &get_ball_candidates(long source_idx, double radius, double open_radius);

    // Euclidean ("visible") search support: tiles of the block ball in order
    // of increasing distance from the source voxel centre, cached per
    // shape/source/radius; visibility scratch array; mode flag.
    std::vector<int> ball_order_;
    long ball_source_ = -1; double ball_radius_ = -1;
    int ball_nx_ = -1, ball_ny_ = -1, ball_nz_ = -1;
    const std::vector<int> &get_ball_order(long source_idx, double radius);
    std::vector<char> visible_;
    bool euclidean_search_ = false;

    // Compact per-tile arrays (structure-of-arrays) the lattice path works
    // on: PathMapTile is ~100 bytes and streaming 12k of them per AV per
    // frame for three floats dominated the O(window) phases. While
    // soa_valid_, `cost` (the search's cost array), penalty_soa_ and
    // density_soa_ are the truth and `tiles` is stale; sync_tiles_from_soa()
    // brings the tiles up to date before any API that reads them.
    std::vector<float> penalty_soa_;
    std::vector<float> density_soa_;
    LinkerWeighting weighting_;      //!< uniform unless a table is set
    bool soa_valid_ = false;
    void sync_tiles_from_soa();

    // The bounded exact search proper: penalties from `penalty`, results
    // into `cost`/`visited`; returns the reached tiles. `record_previous`
    // writes PathMapTile::previous (the tiles-backed variant only).
    std::vector<int> dijkstra_bounded_core(long source_idx, long end_idx,
                                           float max_cost, const float *penalty,
                                           bool record_previous);

    // Exact label-correcting Dijkstra (lazy deletion), see find_path().
    // Stops once the cheapest open tile costs >= max_cost (voxel units);
    // tiles that stay unsettled keep a cost >= max_cost, so an accessible
    // density thresholded below max_cost is unaffected.
    void find_path_dijkstra_exact(long path_begin_idx, long path_end_idx,
                                  float max_cost = std::numeric_limits<float>::infinity(),
                                  bool keep_source_cost_default = false);
    bool exact_search_ = false;
    bool symmetric_stencil_ = false;

protected:

    std::vector<PathMapTile> tiles;
    PathMapHeader pathMapHeader_;
    //! Per-particle radii used instead of the model's; 0 = transparent.
    std::vector<double> obstacle_radii_;
    std::vector<int> offsets_;
    std::vector<PathMapTileEdge>& get_edges(int tile_idx);

public:


    /**

    @brief Updates the tiles in the path map.
    *
    This function updates the tiles in the path map based on the given parameters.
    *
    @param obstacle_threshold The threshold value for considering a cell as an obstacle. Default value is -1.0.
    @param binarize A flag indicating whether to binarize the path map. Default value is true.
    @param obstacle_penalty The penalty value for obstacle cells. Default value is TILE_PENALTY_DEFAULT.
    @param reset_tile_edges A flag indicating whether to reset the edges of the tiles. Default value is true.
    */
    void update_tiles(
        float obstacle_threshold=-1.0, 
        bool binarize=true, 
        float obstacle_penalty=TILE_PENALTY_DEFAULT,
        bool reset_tile_edges=true
    );

    /**
     * @brief Resizes the PathMap object.
     *
     * This function resizes the PathMap object to accommodate the specified number of voxels.
     *
     * @param nvox The number of voxels to resize the PathMap to.
     */
    void resize(long nvox);

    /**

    @brief Sets the data for the path map.
    *
    This function sets the input data for the path map. The input data is a 1D array of doubles representing the map.
    *
    @param input Pointer to the input data array.
    @param n_input Number of elements in the input data array.
    @param obstacle_threshold The threshold value for considering a cell as an obstacle. Default value is -1.
    @param binarize Flag indicating whether to binarize the input data. Default value is true.
    @param obstacle_penalty The penalty value for obstacle cells. Default value is TILE_PENALTY_DEFAULT.
    *
    @note The input data array should be of size n_input.
    @note If binarize is set to true, the input data will be converted to binary values based on the obstacle_threshold.
    @note The obstacle_penalty is used to assign penalty values to obstacle cells in the path map.
    */
    void set_data(double *input, int n_input, 
        float obstacle_threshold=-1, bool binarize=true, 
        float obstacle_penalty=TILE_PENALTY_DEFAULT);

    /**

    Returns a vector of neighbor index offsets within a given radius.
    @param neighbor_radius The radius within which to find neighbors. If negative, the radius is obtained from the path map header.
    @return A vector of neighbor index offsets, where each offset consists of three integers (z, y, x) representing the relative position of the neighbor,
    and one integer representing the tile offset. The vector also includes the edge cost between the current tile and the neighbor.
    */
    std::vector<int> get_neighbor_idx_offsets(double neighbor_radius = -1){
        if(neighbor_radius < 0){
            neighbor_radius = get_path_map_header().get_neighbor_radius();
        }
        const int nn = ceil(neighbor_radius);
        const double nr2 = neighbor_radius * neighbor_radius;

        const GridHeader* header = get_header();
        int nx = header->get_nx();
        int ny = header->get_ny();
        int nx_ny = nx * ny;

        // Historical stencil: the loops run -nn <= z < nn, so the +nn face is
        // missing (a length-2 jump exists towards -x/-y/-z only) and the
        // tile itself is included. The symmetric stencil (see
        // set_symmetric_stencil) runs -nn..nn and drops the self offset.
        const int hi = symmetric_stencil_ ? nn : nn - 1;
        std::vector<int> offsets;
        for(int z = -nn; z <= hi; z += 1) {
            double dz2 = z * z;
            int oz = z * nx_ny;
            for(int y = -nn; y <= hi; y++) {
                int oy = y * nx;
                double dz2_dy2 = dz2 + y * y;
                for(int x = -nn; x <= hi; x++) {
                    int ox = x;
                    int dz2_dy2_dx2 = dz2_dy2 + x * x;
                    if(symmetric_stencil_ && dz2_dy2_dx2 == 0) continue;
                    if(dz2_dy2_dx2 <= nr2){
                        int d;                  // edge_cost is a float stored in an 32bit int
                        float *p = (float*) &d; // Make a float pointer point at the integer
                        *p = sqrt((float) dz2_dy2_dx2); // Pretend that the integer is a float and store the value
                        int tile_offset = oz + oy + ox;
                        offsets.emplace_back(z);
                        offsets.emplace_back(y);
                        offsets.emplace_back(x);
                        offsets.emplace_back(tile_offset);
                        offsets.emplace_back(d);
                    }
                }
            }
        }
        return offsets;
    }

    //! Get index of voxel in an axis
    /*!
     * Returns the index of a voxel on a grid in a certain
     * dimension.
     * @param index voxel index
     * @param dim dimension
     * @return index of voxel in axis of dimension
     */
    int get_dim_index_by_voxel(long index, int dim);

    /**

    @brief Returns a read-only reference to the header of the map.
    *
    @return A read-only reference to the header of the map.
    */
    const PathMapHeader &get_path_map_header() const { return pathMapHeader_; }


    /**

    @brief Returns a writable reference to the header of the map.
    *
    @return A writable reference to the header of the map.
    */
    PathMapHeader &get_path_map_header_writable() { return pathMapHeader_; }


    /**

    @brief Set the path map header.
    *
    This function sets the path map header for the path map.
    *
    @param path_map_header The path map header to set.
    @param resolution The resolution of the path map. Default value is -1.0.
    */
    void set_path_map_header(const PathMapHeader &path_map_header, float resolution = -1.0);


    /**

    @brief Get the values of all tiles.
    *
    A tile in a path map contains information on the penalty for visiting a tile, the cost of a path
    from the origin of a path search to the tile, the density of the tile, and other user-defined information.
    *
    When getting information from a tile, the returned values can be cropped to a specified range.
    *
    @param value_type Specifies the type of the returned information (see: PathMapTileOutputs).
    Depending on the value type, the output can be the penalty for visiting the tile, the total
    cost of a path to the tile, or the density of the tile. Additional user-defined content can also be accessed.
    @param bounds Bound for cropping the output values.
    @param feature_name Name of a feature when accessing additional information.
    @return A vector of values for the specified parameters.
    @relates PathMapTile::get_value
    */
    //! Weight accessible voxels by the linker's chain statistics.
    /*! The weighting multiplies the density of every reached voxel by
        #IMP::bff::LinkerWeighting::get_weight() of that voxel's path length,
        so every quantity read from the map -- mean position, distances,
        efficiencies, the exported grid -- sees it. Setting it invalidates
        whatever was folded in before, so it can be changed and re-read.

        \param[in] w the weighting; the default-constructed one is uniform */
    void set_linker_weighting(const LinkerWeighting& w);
    const LinkerWeighting& get_linker_weighting() const { return weighting_; }

    std::vector<float> get_tile_values(
            int value_type = PM_TILE_COST,
            std::pair<float, float> bounds = std::pair<float, float>(
                    {std::numeric_limits<float>::min(),
                     std::numeric_limits<float>::max()}),
            const std::string &feature_name=""
    );

    /**

    @brief Retrieves the values of the tiles in the path map.
    *
    This function retrieves the values of the tiles in the path map and stores them in the output array.
    *
    @param output A pointer to a 2D array of floats where the tile values will be stored.
    @param nx A pointer to an integer that will store the number of tiles in the x-direction.
    @param ny A pointer to an integer that will store the number of tiles in the y-direction.
    @param nz A pointer to an integer that will store the number of tiles in the z-direction.
    @param value_type The type of value to retrieve for each tile. Defaults to PM_TILE_COST.
    @param bounds A pair of floats representing the lower and upper bounds for the tile values. 
    Defaults to the minimum and maximum float values.
    @param feature_name The name of the feature for which to retrieve the tile values. 
    Defaults to an empty string.
    */
    void get_tile_values(
            float **output, int *nx, int *ny, int *nz,
            int value_type = PM_TILE_COST,
            std::pair<float, float> bounds = std::pair<float, float>(
                    {std::numeric_limits<float>::min(),
                     std::numeric_limits<float>::max()}),
            const std::string &feature_name=""
    );

    /*!
     * Values of tiles
     * @return vector of all tiles in the accessible volume
     */
    std::vector<PathMapTile>& get_tiles();

    /// Change the value of a density inside or outside of a sphere
    /*!
     * Changes the value of the density inside or outside of a sphere.
     * The density is used in the path-search (i.e., the accessible
     * volume calculation
     * @param r0 location of the sphere
     * @param radius radius of the sphere
     * @param value value inside or outside of the sphere
     * @param inverse if set to true (default) the values outside of the sphere
     * are modified. If false the values inside of the sphere are modified.
     */
    void fill_sphere(IMP::algebra::Vector3D r0, double radius, double value, bool inverse=true);

    /**

    @brief Finds a path between two indices in the path map.
    *
    This function finds a path between the specified path begin index and path end index in the path map.
    If the path end index is not specified, the function will find a path from the path begin index to the last index in the path map.
    The heuristic mode parameter determines the heuristic function to be used for path finding.
    *
    @param path_begin_idx The index of the path begin point in the path map.
    @param path_end_idx The index of the path end point in the path map. Default value is -1, which means the last index in the path map.
    @param heuristic_mode The mode of the heuristic function to be used for path finding. Default value is 0.
    *
    @return void
    */
    void find_path(long path_begin_idx, long path_end_idx = -1, int heuristic_mode = 0);


    /**

    @brief Finds the shortest path between two nodes using Dijkstra's algorithm.
    *
    This function finds the shortest path between two nodes in the path map using Dijkstra's algorithm.
    The path is calculated from the node at index path_begin_idx to the node at index path_end_idx.
    If path_end_idx is not provided, the function will calculate the path to the last node in the map.
    *
    @param path_begin_idx The index of the starting node.
    @param path_end_idx The index of the ending node (optional).
    */    
    void find_path_dijkstra(long path_begin_idx, long path_end_idx = -1);

    /**
     * @brief Select the Dijkstra variant used by find_path_dijkstra().
     *
     * The historical search pushes every tile once, at discovery, into a
     * heap ordered by the *live* cost array; a tile whose cost improves after
     * it was pushed is not re-ordered, so it can be settled early and its
     * neighbours relaxed from a non-final cost -- path lengths can come out
     * longer than the true shortest path. The exact variant re-pushes on
     * every improvement and skips stale entries (textbook lazy Dijkstra):
     * true shortest paths, and cheaper per operation. Off by default so the
     * legacy anchoring stays byte-identical; the lattice path of AV turns
     * it on.
     */
    void set_exact_search(bool tf) { exact_search_ = tf; }
    bool get_exact_search() const { return exact_search_; }

    /**
     * @brief Use a symmetric neighbour stencil.
     *
     * The historical stencil is built by loops running -nn <= d < nn: the
     * +nn face is missing (with the default radius 2 a length-2 axis jump
     * exists towards -x/-y/-z but not +x/+y/+z -- a path could tunnel
     * through a one-voxel wall in the negative directions only) and the
     * tile itself is included as a zero-length neighbour. The symmetric
     * stencil runs -nn..nn and drops the self offset. Off by default (the
     * legacy anchoring keeps its stencil); the lattice path of AV turns it
     * on together with a neighbour radius of sqrt(3), i.e. the 26 face,
     * edge and corner neighbours -- no length-2 jumps at all.
     */
    void set_symmetric_stencil(bool tf) {
        if(tf != symmetric_stencil_){ offsets_.clear(); nb_delta_.clear(); }
        symmetric_stencil_ = tf;
    }
    bool get_symmetric_stencil() const { return symmetric_stencil_; }

    /**
     * @brief Exact Dijkstra from `path_begin_idx`, bounded.
     *
     * Textbook lazy Dijkstra that stops once the cheapest open tile costs at
     * least `max_cost` (voxel units): every tile whose final cost is below
     * `max_cost` is settled exactly, tiles beyond keep a tentative or default
     * cost >= max_cost. An accessible density thresholded at
     * `max_cost * spacing` (the AV's linker length) is therefore identical to
     * an unbounded search, at a fraction of the pops. With
     * `keep_source_cost_default` the source tile keeps TILE_COST_DEFAULT,
     * as the historical search leaves it (it never wrote the source's cost).
     */
    void find_path_dijkstra_bounded(long path_begin_idx, float max_cost,
                                    bool keep_source_cost_default = false) {
        find_path_dijkstra_exact(path_begin_idx, -1, max_cost, keep_source_cost_default);
    }

    /**
     * @brief Euclidean ("visible") search instead of the path search.
     *
     * A free tile is reached iff the source voxel sees it along a straight
     * voxel path -- a 3D Bresenham chain of 26-neighbour steps towards the
     * source in which every diagonal step also requires the face neighbours
     * it passes to be free (no corner cutting) -- and its cost is the exact
     * Euclidean distance to the source. Tiles in the shadow of the protein
     * are not reached: the linker is modelled as straight. Off by default;
     * AV::set_search_mode("euclidean") turns it on.
     */
    void set_euclidean_search(bool tf) { euclidean_search_ = tf; }
    bool get_euclidean_search() const { return euclidean_search_; }

    /**
     * @brief The lattice evaluation of one AV in compact arrays.
     *
     * Equivalent to update_tiles() + find_path_dijkstra_bounded(source,
     * max_cost, true) but on structure-of-arrays storage: penalties are
     * binarised from the current data, the bounded exact search runs, and the
     * tiles are left stale until an API reads them (get_tiles,
     * get_tile_values, the tiles-backed searches). get_xyz_density() reads
     * the arrays directly.
     */
    void search_lattice(long source_idx, float max_cost);

    /**
     * @brief search_lattice() with the two source spheres folded in.
     *
     * Same result as fill_sphere(r0, block_radius, TILE_PENALTY_THRESHOLD,
     * true); fill_sphere(r0, open_radius, 0, false); search_lattice(...),
     * but the penalties are formed in one pass straight from the data (the
     * data itself is left untouched).
     */
    void search_lattice(long source_idx, float max_cost,
                        const IMP::algebra::Vector3D &r0,
                        double block_radius, double open_radius);

    //! As above, with the obstacle occupancy taken from an integer covering
    //! count per voxel (`occupancy`, window order) instead of the map data:
    //! a count above the obstacle threshold is an obstacle, exactly as a
    //! data value would be. The map data is left untouched.
    void search_lattice(long source_idx, float max_cost,
                        const IMP::algebra::Vector3D &r0,
                        double block_radius, double open_radius,
                        const int32_t *occupancy);

    //! carve_lattice() from integer counts; the counts are also stored as
    //! the map data (like the second raster of the historical path)
    void carve_lattice(const int32_t *occupancy);

    /**
     * @brief AV3 carve: density is the fraction of dye radii that fit.
     *
     * One occupancy-count array per radius. A voxel gets `k/n` where `k` of the
     * `n` probes clear the obstacles there.
     *
     * This is LabelLib's rule, taken from its source rather than inferred:
     * `Grid3DExt::excludeConcentricSpheres` (`FlexLabel/src/FlexLabel.cxx:235`,
     * LabelLib 2af43ac) **sorts** the radii, builds
     * `rhos = LinSpaced(n + 1, 0, 1)` -- `{0, 1/3, 2/3, 1}` for AV3 -- and for
     * each atom writes `ref = min(ref, rhos[iClash])` over the shell between
     * consecutive `atom_vdW + radius[iClash]`. A voxel inside the smallest probe
     * is 0, inside the middle one 1/3, inside the largest 2/3, outside all of
     * them 1: the fraction of probes that fit. `min` means the nearest atom
     * wins, which the per-radius occupancy counts reproduce.
     *
     * The sort is why the result does not depend on the order of the radii, and
     * it is checked against LabelLib directly in `test_av3_matches_labellib`.
     *
     * `n == 1` reproduces carve_lattice(const int32_t*) exactly, so the AV1
     * path is unchanged.
     */
    void carve_lattice_fractional(const int32_t *const *occupancy, int n);

    /**
     * @brief Re-weight the carved cloud towards the surface: the accessible
     *        **contact** volume (ACV).
     *
     * A voxel of the cloud is a *contact* voxel when some voxel within
     * \p thickness of it is excluded for the dye, i.e. when the dye sitting
     * there is within `thickness` of the molecular surface it cannot enter.
     * The contact voxels are then scaled so that they carry
     * \p trapped_fraction of the cloud's total weight and the rest carry the
     * remainder -- a dye that touches the protein stays there longer than free
     * diffusion would put it.
     *
     * This is Olga's rule, taken from its source rather than inferred:
     * `path2points()` (`Olga/src/AV/fretAV.cpp:236`) collects the accessible
     * points, marks a point trapped when any neighbour within
     * `contactR / gridStep` voxels is occupied in `occupancyVdWDye` -- the
     * obstacle raster inflated by the *dye* radius, which is exactly
     * \p occupancy here -- and then sets every trapped point's weight to
     * `contactRho = volFree * trappedFrac / (volTrapped * (1 - trappedFrac))`
     * with the free points left at 1.
     *
     * Two deliberate differences from that source, both stated because they
     * are choices and not accidents:
     *
     * - Olga fixes the *ratio* of the two weights; this fixes the two
     *   *shares*. They are the same thing whenever the free points all have
     *   weight one, which is Olga's case (its AV1 without chain weighting, and
     *   its AV3 always), so the reweighted cloud differs only by a global
     *   factor that cancels in every distance. Where they differ -- AV3, where
     *   a voxel's weight is the fraction of dye radii that fit, and chain
     *   weighting, where it is the linker statistics -- shares are the
     *   statement the parameter actually makes: *the dye spends
     *   `trapped_fraction` of its time in contact*. Olga instead overwrites
     *   the chain weight of every trapped point, which drops the weighting it
     *   was asked for.
     * - The contact neighbourhood is a true sphere of radius \p thickness
     *   (`|d| * step <= thickness`), where Olga rounds the radius down to
     *   whole voxels first (`deltaIlist(contactR / step)` truncates) and
     *   guards its 1-D offsets by range alone, so its shell wraps around the
     *   grid faces. The difference is sub-voxel except at the faces, where
     *   Olga's is wrong.
     *
     * Must be called **after** carve_lattice()/carve_lattice_fractional() and
     * after the search: the normalisation is over the voxels that actually
     * reach the cloud, so it needs both the density and the costs.
     *
     * A no-op unless `thickness > 0` and `0 <= trapped_fraction < 1`, and a
     * no-op when either side of the split is empty -- there is no ratio to set
     * then, and Olga does nothing in that case either.
     *
     * \param[in] occupancy per-voxel covering counts at the dye radius, in
     *            window order; a count above TILE_OBSTACLE_THRESHOLD is
     *            excluded volume. The same array carve_lattice() was given.
     * \param[in] thickness the contact layer, in Angstrom (fps.json
     *            `contact_volume_thickness`)
     * \param[in] trapped_fraction the share of the cloud's weight the contact
     *            layer carries (fps.json `contact_volume_trapped_fraction`)
     */
    void apply_contact_weighting(const int32_t *occupancy, double thickness,
                                 double trapped_fraction);

    //! get_xyz_density() as four arrays (x, y, z, density), appended to the
    //! given vectors after clearing them; same tiles, same order, same values.
    void get_xyz_density_soa(std::vector<float> &x, std::vector<float> &y,
                             std::vector<float> &z, std::vector<float> &w);

    //! set_origin() without reallocating the location arrays (same values)
    void set_origin_fast(const IMP::algebra::Vector3D &origin);

    /**
     * @brief Set the tile density from the current data: 0 where the data
     * exceeds TILE_OBSTACLE_THRESHOLD, 1 elsewhere (the dye-radius carve).
     * Structure-of-arrays like search_lattice(); the pair belongs together.
     */
    void carve_lattice();
    
    
    /**

    @brief Finds the shortest path between two indices using the A* algorithm.
    *
    This function uses the A* algorithm to find the shortest path between two indices in the path map.
    The path is stored in the 
    path
    member variable.
    *
    @param path_begin_idx The index of the starting point of the path.
    @param path_end_idx The index of the ending point of the path. If not provided, the function will use the last index in the path map.
    */
    void find_path_astar(long path_begin_idx, long path_end_idx = -1);

    /**

    @brief Get the XYZ density of the path map.
    This function returns a vector of IMP::algebra::Vector4D objects representing
    the XYZ density of the path map.
    @return std::vector The XYZ density of the path map.
    */
    std::vector<IMP::algebra::Vector4D> get_xyz_density();
    
    /**

    @brief Get the XYZ density of the path map.
    This function returns the XYZ density of the path map as a 2D array.
    @param output A pointer to a 2D array to store the XYZ density.
    @param n_output1 A pointer to an integer to store the number of rows in the output array.
    @param n_output2 A pointer to an integer to store the number of columns in the output array.
    */
    void get_xyz_density(double** output, int* n_output1, int* n_output2);

    /**

    @brief Resamples the obstacles in the path map.
    *
    This function resamples the obstacles in the path map, updating their positions and sizes.
    *
    @param extra_radius The extra radius to add to the obstacles (optional, default is 0.0).
    */
    void sample_obstacles(double extra_radius=0.0);

    //! Radii to use instead of the particles' own; empty = the particles'.
    /*! A **zero radius is transparent**: it is not inflated by the probe
        radius, and a binarized sphere of radius zero contains no point, so it
        blocks nothing. That is how a `strip_mask` reaches the obstacle set --
        the atoms stay in the list, at their own coordinates, with no size --
        so every volume indexes the same particles whatever it strips, and
        nothing is mutated on a model another thread is reading.

        \param[in] radii one per particle, or empty to clear the override */
    void set_obstacle_radii(const std::vector<double>& radii) {
        obstacle_radii_ = radii;
    }

    //! The radius override, empty when there is none.
    const std::vector<double>& get_obstacle_radii() const {
        return obstacle_radii_;
    }

    /**

    @brief Constructs a PathMap object.
    *
    @param header The PathMapHeader object.
    @param name The name of the PathMap.
    @param kt The kernel type.
    @param resolution The resolution of the PathMap.
    */
    //! A copy of this grid as an `IMP::em::DensityMap`.
    /*! `PathMap` used to *be* one, and callers passed it straight to
        `IMP.em.write_map` and friends. It is its own lattice now
        (`DensityGrid`), so that door is explicit: a new map with the same
        extent, spacing, origin and voxel values, owned by the caller. This is
        IMP integration and moves to the connection layer with
        #write_map_feature; the type is forward-declared so that this header
        still pulls in nothing from `IMP.em`. */
    IMP::em::DensityMap* create_density_map() const;

    //! Take the obstacles from IMP particles.
    /*! Where the obstacle spheres come from at every sample.
        The lattice owes nothing to particles: it keeps spheres (set_spheres()).
        A source, when installed, is asked for fresh spheres at the start of
        every sample_obstacles(), so a caller whose obstacles move -- the
        decorator over an IMP::Model, which the connection layer installs
        through set_path_map_particles() -- samples where they are now, not
        where they were when the spheres were set. With no source the spheres
        set last are used as they are. */
    void set_sphere_source(std::function<void(GridSpheres&)> source) {
        sphere_source_ = std::move(source);
        if (sphere_source_) sphere_source_(xyzr_);
    }
    bool get_has_sphere_source() const { return static_cast<bool>(sphere_source_); }
public:
    std::function<void(GridSpheres&)> sphere_source_;

    explicit PathMap(
            const PathMapHeader &header,
            std::string name = "PathMap%1%",
            float resolution = -1.0
    );
protected:

};


/**
 * @brief Writes one voxel feature of a map to a density file.
 *
 * The map is a #PathMap because that is what carries per-voxel features
 * (`PM_TILE_*`, plus named ones); what is written is a single scalar field
 * chosen by @p value_type, which is why this is not `write_path_map` -- the
 * caller picks a field and gets a density file, and nothing about that is
 * particular to a path search. Anything that fills a lattice and wants to
 * look at one of its fields -- an occupancy, a learned rate, a diffusion
 * coefficient -- writes it the same way.
 *
 * Guesses the file type from the file name. The supported file formats are:
 * - .mrc/.map
 * - .em
 * - .vol
 * - .xplor
 *
 * @param m The PathMap object to write.
 * @param filename The name of the file to write to.
 * @param value_type The value type.
 * @param bounds The bounds of the path map.
 * @param feature_name The name of the feature.
 */
IMPBFFEXPORT
void write_map_feature(
        PathMap *m,
        std::string filename,
        int value_type,
        const std::pair<float, float> bounds = std::pair<float, float>(
                std::numeric_limits<float>::min(),
                std::numeric_limits<float>::max()
        ),
        const std::string &feature_name = ""
);

// -------- the lattice window and the attachment-atom subtraction --------
// Both were statics of the AV decorator's source; the Model-free get_av in
// AVBuilder.cpp needs them too, so they live with the lattice.

//! The cubic window of the global lattice a search from `source` can reach.
/*! `n` is odd and centred on the voxel nearest the source; `k0` is the
    window's lattice origin in global voxel indices. */
IMPBFFEXPORT void lattice_window(const IMP::algebra::Vector3D &source, double ll,
                                 double h, int k0[3], int &n);

/*!
    FPS drops the attachment atom before it rasterises anything
    (`av_routines.cpp:50`: `if (i == atom_i) continue;`), and it has to be
    dropped: an atom cannot block the linker that is tied to it. `IMP.bff` kept
    it, so the source sat inside its own inflated sphere and the volume came
    back empty unless the clearance was raised past it -- a knob standing in
    for a missing rule. Measured on FPS's own reference cloud, dropping it
    takes `p_1bp` from 0.88 to **0.99** of FPS's volume (PRD-121).

    It is done here rather than through the radii override because the
    occupancy raster is **shared** between every volume of the same (spacing,
    extra-radius) class (PRD-105), and each volume drops a *different* atom.
    Giving each its own raster would cost that sharing. The raster stores a
    per-voxel atom **count**, so one atom's contribution subtracts exactly: a
    voxel only this atom covered falls to zero and opens, a voxel any other
    atom covers stays blocked.

    \param[in,out] counts the window, as read from the occupancy
    \param[in] k0,n the window's lattice origin and edge
    \param[in] spacing the lattice spacing
    \param[in] c the attachment atom's position
    \param[in] radius its own radius plus this occupancy's extra radius --
               the same `r + extra_radius` the raster used, so the voxel set
               subtracted is exactly the one that was added
*/
template <typename T>
inline void drop_source_obstruction(T *counts, const IMP::algebra::Vector3D &origin,
                             int n, double spacing,
                             const IMP::algebra::Vector3D &c,
                             double source_radius, double extra) {
    /* Nothing to subtract when the attachment atom contributed nothing. A
       radius of zero means it is not an obstacle -- the array door builds its
       obstacles from a caller's list and the source is not in it -- and
       subtracting a sphere that was never added would open voxels no rule
       opened. */
    if (source_radius <= 0.0) return;
    const double radius = source_radius + extra;
    if (radius <= 0.0) return;
    const double r2 = radius * radius;
    const long nxy = (long) n * n;
    for (int iz = 0; iz < n; iz++) {
        const double dz = origin[2] + iz * spacing - c[2];
        const double dz2 = dz * dz;
        if (dz2 >= r2) continue;
        for (int iy = 0; iy < n; iy++) {
            const double dy = origin[1] + iy * spacing - c[1];
            const double dyz2 = dz2 + dy * dy;
            if (dyz2 >= r2) continue;
            const long row = (long) iz * nxy + (long) iy * n;
            for (int ix = 0; ix < n; ix++) {
                const double dx = origin[0] + ix * spacing - c[0];
                // strict `<`, as AVOccupancyMap::add_sphere tests
                if (dyz2 + dx * dx < r2) {
                    T &v = counts[row + ix];
                    if (v > 0) v -= 1;
                }
            }
        }
    }
}

IMPBFF_END_NAMESPACE


#endif  // IMPBFF_PATHMAP_H
