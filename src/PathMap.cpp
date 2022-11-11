/**
 *  \file IMP/bff/PathMap.h
 *  \brief Class to search path on grids
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */
#include <IMP/bff/PathMap.h>

IMPBFF_BEGIN_NAMESPACE

PathMap::PathMap(
        PathMapHeader &av_header,
        std::string name,
        IMP::em::KernelType kt,
        float resolution
) : SampledDensityMap(kt), pathMapHeader_(av_header)
{
    set_name(name);
    set_path_map_header(av_header, resolution);
}

void PathMap::set_path_map_header(PathMapHeader &av_header, float resolution)
{
    header_ = *av_header.get_density_header();
    header_.compute_xyz_top(true);
    if(resolution < 0)
        resolution = av_header.get_simulation_grid_resolution();
    header_.set_resolution(resolution);
    // allocate the data
    long nvox = get_number_of_voxels();
    resize(nvox);
    kernel_params_ = IMP::em::KernelParameters(header_.get_resolution());
    calc_all_voxel2loc();
}

void PathMap::update_tiles(float obstacle_penalty){
    const float obstacle_threshold = pathMapHeader_.get_obstacle_threshold();
    double* grid_data = get_data();
    for(auto &tile: tiles){
        tile.penalty = (grid_data[tile.idx] > obstacle_threshold) ? obstacle_penalty : 0.0f;
        tile.cost = std::numeric_limits<float>::max();
    }
    const int nn = pathMapHeader_.get_neighbor_box_size();
    for(auto &tile: tiles){
        tile.update_edges(this, tiles, nn, obstacle_threshold);
    }
}

void PathMap::find_path_dijkstra(
        const long path_begin_idx, const long path_end_idx){
    long n_voxel = get_number_of_voxels();
    IMP_USAGE_CHECK(path_begin_idx >= 0 && path_begin_idx < n_voxel && path_end_idx < n_voxel,
                    "find_path_dijkstra: invalid start/stop map index");
    update_tiles();

    PathMapTile* start = &tiles[path_begin_idx];
    PathMapTile* end = nullptr;
    if (path_end_idx > 0) end = &tiles[path_end_idx];
    start->cost = 0.0;
    start->previous = nullptr;

    // Using lambda to compare elements -> use priority_queue as MinHeap
    auto cmp = [](PathMapTile* left, PathMapTile* right) { return right->cost < left->cost; };
    std::priority_queue<PathMapTile*, std::vector<PathMapTile*>, decltype(cmp)> frontier(cmp);
    frontier.push(start);

    // Keep track of visited nodes, std::unordered_set is a hash set
    std::unordered_set<PathMapTile*> visited = std::unordered_set<PathMapTile*>();
    visited.insert(start);

    while(!frontier.empty()){
        PathMapTile* current = frontier.top();
        frontier.pop();
        if (current == end){
            break;
        }
        for(auto &edge : current->edges) {
            PathMapTile* neighbor = edge.tile;
            if(neighbor == nullptr) continue;
            float new_neighbor_cost = current->cost + edge.length + neighbor->penalty;
            if (new_neighbor_cost < neighbor->cost) {
                neighbor->cost = new_neighbor_cost;
                neighbor->previous = current;
            }
            if (visited.find(neighbor) == visited.end()){
                visited.insert(neighbor);
                frontier.push(neighbor);
            }
        }
    }
}

void PathMap::find_path_astar(
        const long begin_idx, const long end_idx){

    auto distance = [&](long idx1, long idx2){
        float d = 0.0f;
        IMP::algebra::Vector3D v1 = get_location_by_voxel(idx1);
        IMP::algebra::Vector3D v2 = get_location_by_voxel(idx2);
        d = IMP::algebra::get_distance(v1, v2);
        return d;
    };

    long n_voxel = get_number_of_voxels();
    IMP_USAGE_CHECK(begin_idx >= 0 && begin_idx < n_voxel && end_idx < n_voxel,
                    "find_path_dijkstra: invalid start/stop map index");
    update_tiles();

    PathMapTile* start = &tiles[begin_idx];
    PathMapTile* end = nullptr;
    if (end_idx > 0) end = &tiles[end_idx];
    start->cost = 0.0;
    start->previous = nullptr;

    // Using lambda to compare elements -> use priority_queue as MinHeap
    // A* uses an additional distance term to bias the search
    // in right direction
    auto cmp = [&](PathMapTile* left, PathMapTile* right) {
        float lhs_cost = left->cost  + distance(left->idx,  end_idx);
        float rhs_cost = right->cost + distance(right->idx, end_idx);
        return rhs_cost < lhs_cost;
    };
    std::priority_queue<PathMapTile*, std::vector<PathMapTile*>, decltype(cmp)> frontier(cmp);
    frontier.push(start);

    // Keep track of visited nodes, std::unordered_set is a hash set
    std::unordered_set<PathMapTile*> visited = std::unordered_set<PathMapTile*>();
    visited.insert(start);

    while(!frontier.empty()){
        PathMapTile* current = frontier.top();
        frontier.pop();
        if (current == end){
            break;
        }
        for(auto &edge : current->edges) {
            PathMapTile* neighbor = edge.tile;
            if(neighbor == nullptr) continue;
            float new_neighbor_cost = current->cost + edge.length + neighbor->penalty;
            if (new_neighbor_cost < neighbor->cost) {
                neighbor->cost = new_neighbor_cost;
                neighbor->previous = current;
            }
            if (visited.find(neighbor) == visited.end()){
                visited.insert(neighbor);
                frontier.push(neighbor);
            }
        }
    }

}

void PathMap::fill_sphere(
        IMP::algebra::Vector3D r0, double radius,
        double value, bool inverse
){
    double dsq = radius * radius;
    IMP::algebra::Vector3D vox_cent;
    long n_vox = get_number_of_voxels();
    // Loop over all voxels
    for (long v=0; v<n_vox; ++v) {
        vox_cent = get_location_by_voxel(v);
        double vox_dist = IMP::algebra::get_squared_distance(vox_cent, r0);
        if(!inverse)
            if(vox_dist < dsq) data_[v] = value;
        if(inverse)
            if(vox_dist >= dsq) data_[v] = value;
    }
}

int PathMap::get_dim_index_by_voxel(long index, int dim){
    IMP_USAGE_CHECK(index >= 0 && index < get_number_of_voxels(), "invalid start map index");
    int ny = header_.get_ny();
    int nx = header_.get_nx();
    switch (dim) {
        case 0:
            return static_cast<int>(index % nx);
        case 1:
            return static_cast<int>((index / nx) % ny);
        default:
            return static_cast<int>(index / (nx * ny));
    }
}

void PathMap::sample_obstacles(double extra_radius){
    set_origin(pathMapHeader_.get_origin());
    std::vector<double> radii_original;
    // 1. Update radii = radius + dye radius
    if(extra_radius!=0.0){
        radii_original.reserve(xyzr_.size());
        for(auto &p : xyzr_){
            double r = p.get_radius();
            radii_original.emplace_back(r);
            p.set_radius(r + extra_radius);
        }
    }

    // 2. Use sampled density map to place atoms (density in map)
    SampledDensityMap::resample();

    // Restore radii
    if(extra_radius!=0.0){
        radii_original.reserve(xyzr_.size());
        for (size_t i = 0; i < radii_original.size(); i++) {
            xyzr_[i].set_radius(radii_original[i]);
        }
    }
}

std::vector<IMP::algebra::Vector4D> PathMap::get_xyz_density(){
    long n_voxel = get_number_of_voxels();
    float linker_length = pathMapHeader_.get_max_path_length();
    float grid_spacing = pathMapHeader_.get_simulation_grid_resolution();

    std::vector<IMP::algebra::Vector4D> v;
    v.reserve(n_voxel);
    for(int i = 0; i < n_voxel; i++){
        IMP::algebra::Vector3D r = get_location_by_voxel(i);
        float density = tiles[i].get_value(
                PM_TILE_ACCESSIBLE_DENSITY,
                std::pair<float, float>({0.0f, linker_length}), "",
                grid_spacing
        );
        v.emplace_back(
                IMP::algebra::Vector4D({r[0], r[1], r[2], density})
        );
    }
    return v;
}

void write_path_map(
    PathMap *d,
    std::string name,
    int value_type,
    std::pair<float, float> bounds,
    const std::string &feature_name
) {
    IMP_USAGE_CHECK(name.rfind('.') != std::string::npos, "No suffix in file name: " << name);
    std::string suf = name.substr(name.rfind('.'));
    Pointer<IMP::em::MapReaderWriter> rw;
    if (suf == ".mrc" || suf == ".mrcs" || suf == ".map") {
        rw = new IMP::em::MRCReaderWriter();
    } else if (suf == ".em") {
        rw = new IMP::em::EMReaderWriter();
    } else if (suf == ".vol") {
        rw = new IMP::em::SpiderMapReaderWriter();
    } else if (suf == ".xplor") {
        rw = new IMP::em::XplorReaderWriter();
    } else {
        IMP_THROW("Unable to determine type for file " << name << " with suffix " << suf, IOException);
    }
    rw->set_was_used(true);
    d->set_was_used(true);
    std::vector<float> f_data;
    f_data = d->get_tile_values(value_type, bounds, feature_name);
    rw->write(name.c_str(), f_data.data(), *d->get_header());
}

std::vector<float> PathMap::get_tile_values(
        const int value_type,
        std::pair<float, float> bounds,
        const std::string &feature_name
){
    size_t size = tiles.size();
    std::vector<float> data;
    data.reserve(size);
    float grid_spacing = get_spacing();

    for(auto &tile : tiles){
        float value = tile.get_value(
                value_type, bounds,
                feature_name, grid_spacing);
        data.emplace_back(value);
    }
    return data;
}

void PathMap::get_tile_values(
        float **output, int *nx, int *ny, int *nz,
        int value_type,
        std::pair<float, float> bounds,
        const std::string &feature_name){
    int n_voxel = get_number_of_voxels();
    float grid_spacing = get_spacing();
    *nx = header_.get_nx();
    *ny = header_.get_ny();
    *nz = header_.get_nz();
    auto* o = static_cast<float*>(malloc(sizeof(float) * n_voxel));
    for(int i = 0; i < n_voxel; i++){
        o[i] = tiles[i].get_value(value_type, bounds, feature_name, grid_spacing);
    }
    *output = o;
}

void PathMap::resize(unsigned int nvox){
    data_.reset(new double[nvox]);
    tiles.resize(nvox);
    for(unsigned int i = 0; i < nvox; i++) tiles[i].idx = i;
}

std::vector<PathMapTile>& PathMap::get_tiles(){
    return tiles;
}


IMPBFF_END_NAMESPACE
