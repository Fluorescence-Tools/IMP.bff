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
    std::vector<bool>  visited;
    std::vector<bool>  edge_computed;
    std::vector<float> cost;

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
            auto pmh = get_path_map_header();
            neighbor_radius = pmh->get_neighbor_radius();
        }
        const int nn = ceil(neighbor_radius);
        const double nr2 = neighbor_radius * neighbor_radius;

        const IMP::em::DensityHeader* header = get_header();
        int nx = header->get_nx();
        int ny = header->get_ny();
        int nx_ny = nx * ny;

        std::vector<int> offsets;
        for(int z = -nn; z < nn; z += 1) {
            double dz2 = z * z;
            int oz = z * nx_ny;
            for(int y = -nn; y < nn; y++) {
                int oy = y * nx;
                double dz2_dy2 = dz2 + y * y;
                for(int x = -nn; x < nn; x++) {
                    int ox = x;
                    int dz2_dy2_dx2 = dz2_dy2 + x * x;
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

    @brief Returns a read-only pointer to the header of the map.
    *
    @return const PathMapHeader* A read-only pointer to the header of the map.
    */
    const PathMapHeader *get_path_map_header() const { return &pathMapHeader_; }


    /**

    @brief Returns a pointer to the header of the map in a writable version.
    *
    @return A pointer to the header of the map in a writable version.
    */
    PathMapHeader *get_path_map_header_writable() { return &pathMapHeader_; }


    /**

    @brief Set the path map header.
    *
    This function sets the path map header for the path map.
    *
    @param path_map_header The path map header to set.
    @param resolution The resolution of the path map. Default value is -1.0.
    */
    void set_path_map_header(PathMapHeader &path_map_header, float resolution = -1.0);


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
            PathMapHeader &header,
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
