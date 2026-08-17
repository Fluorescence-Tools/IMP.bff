/**
 *  \file IMP/bff/PathMap.h
 *  \brief Class to search path on grids
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2023 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_PATHMAP_H
#define IMPBFF_PATHMAP_H

#include <IMP/bff/bff_config.h>

#include <stdlib.h>     /* malloc, free, rand */
#include <limits>
#include <cstdint>

#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <queue>
#include <vector>
#include <utility>  /* std::pair */
#include <Eigen/Dense>

#include <IMP/Object.h>
#include <IMP/Particle.h>
#include <IMP/core/XYZR.h>
#include <IMP/em/SampledDensityMap.h>

#include <IMP/em/MRCReaderWriter.h>
#include <IMP/em/XplorReaderWriter.h>
#include <IMP/em/EMReaderWriter.h>
#include <IMP/em/SpiderReaderWriter.h>

#include <IMP/bff/PathMapHeader.h>
#include <IMP/bff/PathMapTile.h>
#include <IMP/bff/PathMapTileEdge.h>

IMPBFF_BEGIN_NAMESPACE

class PathMapTile;


class IMPBFFEXPORT PathMap : public IMP::em::SampledDensityMap {

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
    void resize(unsigned int nvox);

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

        const IMP::em::DensityHeader* header = get_header();
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

    //! get_xyz_density() as four arrays (x, y, z, density), appended to the
    //! given vectors after clearing them; same tiles, same order, same values.
    void get_xyz_density_soa(std::vector<float> &x, std::vector<float> &y,
                             std::vector<float> &z, std::vector<float> &w);

    //! set_origin() without reallocating the location arrays (same values)
    void set_origin_fast(const IMP::algebra::Vector3D &origin);
    long loc_size_ = -1;   //!< voxels the location arrays were computed for

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

    /**

    @brief Constructs a PathMap object.
    *
    @param header The PathMapHeader object.
    @param name The name of the PathMap.
    @param kt The kernel type.
    @param resolution The resolution of the PathMap.
    */
    explicit PathMap(
            const PathMapHeader &header,
            std::string name = "PathMap%1%",
            IMP::em::KernelType kt = IMP::em::BINARIZED_SPHERE,
            float resolution = -1.0
    );

};


/**
 * @brief Writes a path map to a file.
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
IMPEMEXPORT
void write_path_map(
        PathMap *m,
        std::string filename,
        int value_type,
        const std::pair<float, float> bounds = std::pair<float, float>(
                std::numeric_limits<float>::min(),
                std::numeric_limits<float>::max()
        ),
        const std::string &feature_name = ""
);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_PATHMAP_H
