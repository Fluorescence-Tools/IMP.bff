/**
 *  \file IMP/bff/PathMap.h
 *  \brief Class to search path on grids
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2023 IMP Inventors. All rights reserved.
 *
 */
#include <IMP/bff/PathMap.h>

#include <cstring>

IMPBFF_BEGIN_NAMESPACE

PathMap::PathMap(
        const PathMapHeader &av_header,
        std::string name,
        IMP::em::KernelType kt,
        float resolution
) : SampledDensityMap(kt), pathMapHeader_(av_header)
{
    set_name(name);
    set_path_map_header(av_header, resolution);
}

void PathMap::set_path_map_header(const PathMapHeader &av_header, float resolution)
{
    const IMP::em::DensityHeader *nh = av_header.get_density_header();
    if(nh->get_nx() != header_.get_nx() || nh->get_ny() != header_.get_ny()
       || nh->get_nz() != header_.get_nz()){
        // neighbour offsets are linear-index deltas of the old shape
        offsets_.clear();
    }
    pathMapHeader_ = av_header;
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

template<class Cmp>
void PathMap::find_path_impl(
        const long path_begin_idx,
        const long path_end_idx,
        Cmp cmp
) {
    long n_voxel = get_number_of_voxels();

    // Set default initial values for path search
    cost.resize(0);
    cost.resize(n_voxel, TILE_COST_DEFAULT);

    // fast lookup which idx was visited
    visited.resize(0);
    visited.resize(n_voxel, false);

    // Store which idx were visited
    std::vector<int> visited_idx;
    visited_idx.reserve(1024);

    // A compact copy of the tile penalties: the relaxation reads it per
    // neighbour, and PathMapTile is far too large to stream through cache.
    std::vector<float> penalty(n_voxel);
    for(long i = 0; i < n_voxel; i++) penalty[i] = tiles[i].penalty;

    // Get start and end tile
    PathMapTile* start = &tiles[path_begin_idx];
    PathMapTile* end = nullptr;
    if (path_end_idx > 0) end = &tiles[path_end_idx];

    // priority_queue stores the tile indices in the frontier
    std::priority_queue<long, std::vector<long>, Cmp> frontier(cmp);

    // Neighbours are enumerated straight from the offset table instead of
    // materialising a per-tile edge vector every frame (get_edges() still
    // does that for callers that want the vectors). Same order, same bounds
    // and the same penalty test as PathMapTile::update_edges_2, so the
    // relaxation sequence -- and the result -- is unchanged.
    if(offsets_.empty()){
        offsets_ = get_neighbor_idx_offsets();
    }
    const std::vector<int> &off = offsets_;
    const size_t noff = off.size();
    const float penalty_threshold = pathMapHeader_.get_obstacle_threshold();
    const int nx = header_.get_nx();
    const int ny = header_.get_ny();
    const int nz = header_.get_nz();
    const int nxy = nx * ny;
    // Cheap interior test: a tile whose position is at least the neighbour
    // box away from every face needs no bounds check.
    const int box = pathMapHeader_.get_neighbor_box_size();

    // perform the search
    cost[path_begin_idx] = 0.0;
    visited[start->idx] = true;
    start->previous = nullptr;
    frontier.push(path_begin_idx);

    while(!frontier.empty()){
        const long cidx = frontier.top();
        frontier.pop();
        PathMapTile* current = &tiles[cidx];
        if (current == end)
            break;
        const int x0 = (int) (cidx % nx);
        const int y0 = (int) ((cidx / nx) % ny);
        const int z0 = (int) (cidx / nxy);
        const bool interior = x0 >= box && x0 < nx - box &&
                              y0 >= box && y0 < ny - box &&
                              z0 >= box && z0 < nz - box;
        const float ccost = cost[cidx];
        for(size_t i = 0; i < noff; i += 5){
            if(!interior){
                int iz = z0 + off[i + 0];
                int iy = y0 + off[i + 1];
                int ix = x0 + off[i + 2];
                if(iz >= nz || iz < 0) continue;
                if(iy >= ny || iy < 0) continue;
                if(ix >= nx || ix < 0) continue;
            }
            const long nidx = cidx + off[i + 3];
            const float npen = penalty[nidx];
            if(!(npen < penalty_threshold)) continue;
            float edge_length;
            std::memcpy(&edge_length, &off[i + 4], sizeof(float));
            float new_neighbor_cost = ccost + edge_length + npen;
            if (new_neighbor_cost < cost[nidx]) {
                cost[nidx] = new_neighbor_cost;
                tiles[nidx].previous = current;
            }
            if(!visited[nidx]){
                visited[nidx] = true;
                visited_idx.emplace_back(nidx);
                frontier.push(nidx);
            }
        }
    }

    for(int &idx : visited_idx){
        tiles[idx].cost = cost[idx];
    }
    // `visited` now marks every reached tile (the source included)
    reached_valid_ = true;
}

void PathMap::find_path(
        const long path_begin_idx, 
        const long path_end_idx,
        const int heuristic_mode
) {
    // std::cout << "void PathMap::find_path(" << std::endl;
    sync_tiles_from_soa();
    long n_voxel = get_number_of_voxels();
    IMP_USAGE_CHECK(
        path_begin_idx >= 0 && 
        path_begin_idx < n_voxel && 
        path_end_idx < n_voxel,
        "PathMap::find_path: invalid start/stop index"
    );

    if (heuristic_mode == 0 || (heuristic_mode != 1 && heuristic_mode != 2)) {
        // dijkstra: the frontier is ordered by the live cost alone. This is
        // the same comparison as the generic one below with a zero heuristic
        // (float + 0.0 converted back to float is the float), without the
        // std::function call.
        if (heuristic_mode != 0) {
            std::cout << "PathMap::find_path: Invalid heuristic_mode. Defaulting to Dijkstra.\n";
        }
        auto cmp = [this](long left, long right) {
            return cost[right] < cost[left];
        };
        find_path_impl(path_begin_idx, path_end_idx, cmp);
        return;
    }

    // define distance heuristic to target for priority_queue
    std::function<double(long, long)> heuristic;
    if (heuristic_mode == 1) { // A* (euclidian)
        heuristic = [&](long left, long right){
            auto dx = (x_loc_[left] - x_loc_[right]);
            auto dy = (y_loc_[left] - y_loc_[right]);
            auto dz = (z_loc_[left] - z_loc_[right]);
            auto d2 = dx*dx + dy*dy + dz*dz;
            return sqrt(d2);
        };
    } else { // A* (manhattan)
        heuristic = [&](long left, long right){
            auto ax = abs(x_loc_[left] - x_loc_[right]);
            auto ay = abs(y_loc_[left] - y_loc_[right]);
            auto az = abs(z_loc_[left] - z_loc_[right]);
            return ax + ay + az;
        };
    }
    auto cmp = [&](long left_idx, long right_idx) {
        float lhs_cost = cost[left_idx]  + heuristic(left_idx,  path_end_idx);
        float rhs_cost = cost[right_idx] + heuristic(right_idx, path_end_idx);
        return rhs_cost < lhs_cost;
    };
    find_path_impl(path_begin_idx, path_end_idx, cmp);
}

void PathMap::update_tiles(
    float obstacle_threshold, 
    bool binarize, 
    float obstacle_penalty,
    bool reset_tile_edges
){
    // Connection between tiles are computed when needed
    // Thus, it is enough to set edge_computed to false for
    // all edges.
    long nvox = get_number_of_voxels();

    if(obstacle_threshold < 0) 
        obstacle_threshold = pathMapHeader_.get_obstacle_threshold();    

    if(reset_tile_edges){
        // Was `resize(false, nvox)` -- arguments swapped, which resized the
        // flags to zero and left get_edges() reading the previous frame's
        // (freed, still-set) bits: after the first evaluation edges were never
        // recomputed and a moved structure kept the old connectivity.
        edge_computed.assign(nvox, false);
    }

    normalized_ = false;
    rms_calculated_ = false;
    reached_valid_ = false;
    soa_valid_ = false;
    for(int idx = 0; idx < nvox; idx++){
        auto value = data_[idx];
        if(binarize){
            value = (value > obstacle_threshold) ? obstacle_penalty : 0.0f;
        }
        tiles[idx].penalty = value;
        tiles[idx].cost = TILE_COST_DEFAULT;
    }

}

std::vector<PathMapTileEdge>& PathMap::get_edges(int tile_idx){
    auto tile = &tiles[tile_idx];
    if(!edge_computed[tile_idx]){
        if(offsets_.empty()){
            offsets_ = get_neighbor_idx_offsets();
        }
        const float obstacle_threshold = pathMapHeader_.get_obstacle_threshold();
        auto header = get_header();
        int nx = header->get_nx();
        int ny = header->get_ny();
        int nz = header->get_nz();
        tile->update_edges_2(nx, ny, nz, tiles, offsets_, obstacle_threshold);
        edge_computed[tile_idx] = true;
    }
    return tiles[tile_idx].edges;
}

void PathMap::find_path_dijkstra(
        const long begin_idx, 
        const long end_idx
){
    if(exact_search_){
        find_path_dijkstra_exact(begin_idx, end_idx);
    } else {
        find_path(begin_idx, end_idx, 0);
    }
}

std::vector<int> PathMap::dijkstra_bounded_core(
        const long path_begin_idx,
        const long path_end_idx,
        const float max_cost,
        const float *penalty,
        const bool record_previous
){
    long n_voxel = get_number_of_voxels();
    cost.assign(n_voxel, TILE_COST_DEFAULT);
    visited.assign(n_voxel, false);

    if(offsets_.empty()){
        offsets_ = get_neighbor_idx_offsets();
    }
    const std::vector<int> &off = offsets_;
    const size_t nnb = off.size() / 5;
    // Compact copies of the stencil for the inner loop: linear index deltas
    // and edge lengths as floats (the offset table stores lengths bit-cast
    // into ints), plus the per-axis deltas for the boundary check.
    std::vector<long> nb_delta(nnb);
    std::vector<float> nb_len(nnb);
    std::vector<int> nb_dz(nnb), nb_dy(nnb), nb_dx(nnb);
    for(size_t j = 0; j < nnb; j++){
        nb_dz[j] = off[5*j + 0]; nb_dy[j] = off[5*j + 1]; nb_dx[j] = off[5*j + 2];
        nb_delta[j] = off[5*j + 3];
        std::memcpy(&nb_len[j], &off[5*j + 4], sizeof(float));
    }
    const float penalty_threshold = pathMapHeader_.get_obstacle_threshold();
    const int nx = header_.get_nx();
    const int ny = header_.get_ny();
    const int nz = header_.get_nz();
    const int nxy = nx * ny;
    const int box = pathMapHeader_.get_neighbor_box_size();

    // Monotone bucket queue. Every edge is at least one voxel long and
    // penalties are non-negative, so a relaxation from a tile of cost c
    // produces a key >= c + 1: it always lands in a later unit bucket than
    // the one being drained. Each bucket is therefore sorted once, by
    // (cost, idx), when it becomes active -- the exact pop order of a binary
    // heap of (cost, idx) pairs, at a fraction of the heap traffic. Stale
    // entries (improved after being queued) are skipped.
    // Consequently, once bucket b becomes active every tile in it already has
    // its final cost (an improvement would have to come from a settled tile of
    // cost < b, i.e. from a bucket drained before), so the pop order inside a
    // bucket does not affect the costs. A tile is queued at most once per
    // bucket (`queued_in`); an entry whose current cost has moved to a lower
    // bucket is stale. The bucket is sorted by (cost, idx) only when
    // predecessors are recorded, to make them deterministic.
    std::vector<std::vector<long> > buckets;
    std::vector<int> queued_in(n_voxel, -1);
    auto push = [&](float c, long idx){
        int b = (int) c;
        if(queued_in[idx] == b) return;
        queued_in[idx] = b;
        if((size_t) b >= buckets.size()) buckets.resize(b + 1);
        buckets[b].push_back(idx);
    };
    std::vector<int> reached;
    reached.reserve(1024);

    cost[path_begin_idx] = 0.0f;
    visited[path_begin_idx] = true;
    if(record_previous) tiles[path_begin_idx].previous = nullptr;
    reached.push_back((int) path_begin_idx);
    push(0.0f, path_begin_idx);

    bool done = false;
    for(size_t b = 0; b < buckets.size() && !done; b++){
        // (index access throughout: push() may reallocate `buckets`; the
        // active bucket itself never grows)
        if(record_previous){
            std::sort(buckets[b].begin(), buckets[b].end(), [&](long l, long r){
                return cost[l] < cost[r] || (cost[l] == cost[r] && l < r);
            });
        }
        const size_t nb = buckets[b].size();
        for(size_t e = 0; e < nb; e++){
            const long cidx = buckets[b][e];
            const float ccost = cost[cidx];
            if((size_t) ccost != b) continue;   // stale: moved to a lower bucket
            if(ccost >= max_cost){ done = true; break; }   // nothing cheaper is left
            if(cidx == path_end_idx){ done = true; break; }
            const int x0 = (int) (cidx % nx);
            const int y0 = (int) ((cidx / nx) % ny);
            const int z0 = (int) (cidx / nxy);
            const bool interior = x0 >= box && x0 < nx - box &&
                                  y0 >= box && y0 < ny - box &&
                                  z0 >= box && z0 < nz - box;
            const float *cost_ptr = cost.data();
            for(size_t j = 0; j < nnb; j++){
                if(!interior){
                    int iz = z0 + nb_dz[j];
                    int iy = y0 + nb_dy[j];
                    int ix = x0 + nb_dx[j];
                    if(iz >= nz || iz < 0) continue;
                    if(iy >= ny || iy < 0) continue;
                    if(ix >= nx || ix < 0) continue;
                }
                const long nidx = cidx + nb_delta[j];
                const float npen = penalty[nidx];
                if(!(npen < penalty_threshold)) continue;
                const float new_cost = ccost + nb_len[j] + npen;
                if(new_cost < cost_ptr[nidx]){
                    cost[nidx] = new_cost;
                    if(record_previous) tiles[nidx].previous = &tiles[cidx];
                    if(!visited[nidx]){
                        visited[nidx] = true;
                        reached.push_back((int) nidx);
                    }
                    push(new_cost, nidx);
                }
            }
        }
    }
    reached_valid_ = true;
    return reached;
}

void PathMap::find_path_dijkstra_exact(
        const long path_begin_idx,
        const long path_end_idx,
        const float max_cost,
        const bool keep_source_cost_default
){
    long n_voxel = get_number_of_voxels();
    IMP_USAGE_CHECK(
        path_begin_idx >= 0 &&
        path_begin_idx < n_voxel &&
        path_end_idx < n_voxel,
        "PathMap::find_path_dijkstra_exact: invalid start/stop index"
    );
    sync_tiles_from_soa();
    std::vector<float> penalty(n_voxel);
    for(long i = 0; i < n_voxel; i++) penalty[i] = tiles[i].penalty;
    std::vector<int> reached = dijkstra_bounded_core(
        path_begin_idx, path_end_idx, max_cost, penalty.data(), true);
    for(int idx : reached){
        tiles[idx].cost = cost[idx];
    }
    if(keep_source_cost_default){
        tiles[path_begin_idx].cost = TILE_COST_DEFAULT;
    }
}

void PathMap::search_lattice(long source_idx, float max_cost){
    long n_voxel = get_number_of_voxels();
    IMP_USAGE_CHECK(source_idx >= 0 && source_idx < n_voxel,
                    "PathMap::search_lattice: invalid source index");
    // penalties, binarised exactly as update_tiles() does
    const float obstacle_threshold = pathMapHeader_.get_obstacle_threshold();
    penalty_soa_.resize(n_voxel);
    for(long i = 0; i < n_voxel; i++){
        penalty_soa_[i] = (data_[i] > obstacle_threshold) ? TILE_PENALTY_DEFAULT : 0.0f;
    }
    normalized_ = false;
    rms_calculated_ = false;
    dijkstra_bounded_core(source_idx, -1, max_cost, penalty_soa_.data(), false);
    // the historical search never wrote the source tile's cost
    cost[source_idx] = TILE_COST_DEFAULT;
    if(density_soa_.size() != (size_t) n_voxel){
        density_soa_.assign(n_voxel, 1.0f);
    }
    soa_valid_ = true;
}

void PathMap::search_lattice(long source_idx, float max_cost,
                             const IMP::algebra::Vector3D &r0,
                             double block_radius, double open_radius){
    long n_voxel = get_number_of_voxels();
    IMP_USAGE_CHECK(source_idx >= 0 && source_idx < n_voxel,
                    "PathMap::search_lattice: invalid source index");
    calc_all_voxel2loc();
    const float obstacle_threshold = pathMapHeader_.get_obstacle_threshold();
    const double bsq = block_radius * block_radius;
    const double osq = open_radius * open_radius;
    const float *xl = x_loc_.get();
    const float *yl = y_loc_.get();
    const float *zl = z_loc_.get();
    const double x0 = r0[0], y0 = r0[1], z0 = r0[2];
    penalty_soa_.resize(n_voxel);
    for(long v = 0; v < n_voxel; v++){
        // fill_sphere(block, inverse) then fill_sphere(open) then binarise:
        // inside the open sphere wins, then beyond the block radius, then
        // the data
        double dx = (double) xl[v] - x0, dy = (double) yl[v] - y0, dz = (double) zl[v] - z0;
        double d2 = dx*dx + dy*dy + dz*dz;
        double value;
        if(d2 < osq) value = 0.0;
        else if(d2 >= bsq) value = TILE_PENALTY_THRESHOLD;
        else value = data_[v];
        penalty_soa_[v] = (value > obstacle_threshold) ? TILE_PENALTY_DEFAULT : 0.0f;
    }
    normalized_ = false;
    rms_calculated_ = false;
    dijkstra_bounded_core(source_idx, -1, max_cost, penalty_soa_.data(), false);
    cost[source_idx] = TILE_COST_DEFAULT;
    if(density_soa_.size() != (size_t) n_voxel){
        density_soa_.assign(n_voxel, 1.0f);
    }
    soa_valid_ = true;
}

void PathMap::carve_lattice(){
    long n_voxel = get_number_of_voxels();
    density_soa_.resize(n_voxel);
    for(long i = 0; i < n_voxel; i++){
        density_soa_[i] = (data_[i] > TILE_OBSTACLE_THRESHOLD) ? 0.0f : 1.0f;
    }
}

void PathMap::sync_tiles_from_soa(){
    if(!soa_valid_) return;
    long n_voxel = get_number_of_voxels();
    for(long i = 0; i < n_voxel; i++){
        tiles[i].penalty = penalty_soa_[i];
        tiles[i].cost = cost[i];
        tiles[i].density = density_soa_[i];
        tiles[i].previous = nullptr;
    }
    soa_valid_ = false;
}

void PathMap::find_path_astar(
    const long begin_idx, 
    const long end_idx
){
    find_path(begin_idx, end_idx, 1);
}

void PathMap::set_data(
    double *input, int n_input, 
    float obstacle_threshold, 
    bool binarize, 
    float obstacle_penalty
){
    long n_voxel = get_number_of_voxels();
    IMP_USAGE_CHECK(n_voxel == n_input, "Invalid size of input data");

    for(int idx = 0; idx < n_input; idx++){
        data_[idx] = input[idx];
    }
    update_tiles(obstacle_threshold, binarize, obstacle_penalty);

}

void PathMap::fill_sphere(
        IMP::algebra::Vector3D r0, double radius,
        double value, bool inverse
){
    double dsq = radius * radius;
    long n_vox = get_number_of_voxels();
    calc_all_voxel2loc();
    // The same arithmetic as get_squared_distance(get_location_by_voxel(v), r0)
    // -- float voxel centres widened to double -- without a Vector3D per voxel
    // and without the per-call usage checks.
    const float *xl = x_loc_.get();
    const float *yl = y_loc_.get();
    const float *zl = z_loc_.get();
    const double x0 = r0[0], y0 = r0[1], z0 = r0[2];
    // Loop over all voxels
    if(inverse){
        for (long v=0; v<n_vox; ++v) {
            double dx = (double) xl[v] - x0, dy = (double) yl[v] - y0, dz = (double) zl[v] - z0;
            double vox_dist = dx*dx + dy*dy + dz*dz;
            if(vox_dist >= dsq) data_[v] = value;
        }
    } else {
        for (long v=0; v<n_vox; ++v) {
            double dx = (double) xl[v] - x0, dy = (double) yl[v] - y0, dz = (double) zl[v] - z0;
            double vox_dist = dx*dx + dy*dy + dz*dz;
            if(vox_dist < dsq) data_[v] = value;
        }
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
    if(soa_valid_){
        // same arithmetic as PathMapTile::get_value(PM_TILE_ACCESSIBLE_DENSITY)
        calc_all_voxel2loc();
        const float *xl = x_loc_.get();
        const float *yl = y_loc_.get();
        const float *zl = z_loc_.get();
        for(long i = 0; i < n_voxel; i++){
            if(!visited[i]) continue;
            float c = cost[i] * grid_spacing;
            float density = (c >= 0.0f && c < linker_length) ? 1.0 : 0.0f;
            density *= density_soa_[i];
            if(density > 0){
                v.emplace_back(IMP::algebra::Vector4D({(double) xl[i], (double) yl[i], (double) zl[i], (double) density}));
            }
        }
        return v;
    }
    auto emit = [&](long i){
        float density = tiles[i].get_value(
                PM_TILE_ACCESSIBLE_DENSITY,
                std::pair<float, float>({0.0f, linker_length}), "",
                grid_spacing
        );
        if(density > 0){
            IMP::algebra::Vector3D r = get_location_by_voxel(i);
            auto n = IMP::algebra::Vector4D({r[0], r[1], r[2], density});
            v.emplace_back(n);
        }
    };
    if(reached_valid_){
        // unreached tiles keep TILE_COST_DEFAULT and contribute nothing
        for(long i = 0; i < n_voxel; i++) if(visited[i]) emit(i);
    } else {
        for(long i = 0; i < n_voxel; i++) emit(i);
    }

    return v;
}

void PathMap::get_xyz_density(double** output, int* n_output1, int* n_output2){
    sync_tiles_from_soa();
    long n_voxel = get_number_of_voxels();
    float linker_length = pathMapHeader_.get_max_path_length();
    float grid_spacing = pathMapHeader_.get_simulation_grid_resolution();

    int n_dim = 4;
    int n = 0;
    auto* t = (double*) calloc(n_voxel * n_dim, sizeof(double)); 
    auto emit = [&](long i){
        double density = tiles[i].get_value(
                PM_TILE_ACCESSIBLE_DENSITY,
                std::pair<float, float>({0.0f, linker_length}), "",
                grid_spacing
        );
        if(density > 0){
            IMP::algebra::Vector3D r = get_location_by_voxel(i);            
            t[n * n_dim + 0] = r[0];
            t[n * n_dim + 1] = r[1];
            t[n * n_dim + 2] = r[2];
            t[n * n_dim + 3] = density;
            n += 1;
        }
    };
    if(reached_valid_){
        for(long i = 0; i < n_voxel; i++) if(visited[i]) emit(i);
    } else {
        for(long i = 0; i < n_voxel; i++) emit(i);
    }
    *n_output1 = (int) n;
    *n_output2 = (int) n_dim;
    *output = t;
}

void write_path_map(
    PathMap *d,
    std::string name,
    int value_type,
    const std::pair<float, float> bounds,
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
    sync_tiles_from_soa();
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
    sync_tiles_from_soa();
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
    reached_valid_ = false;
    soa_valid_ = false;
    penalty_soa_.clear();
    density_soa_.clear();

    edge_computed.resize(0);
    edge_computed.resize(nvox, false);

    tiles.resize(0);
    for(int i = 0; i < nvox; i++){
        auto tile = PathMapTile(i);
        tiles.emplace_back(tile);
    }

}

std::vector<PathMapTile>& PathMap::get_tiles(){
    sync_tiles_from_soa();
    return tiles;
}


IMPBFF_END_NAMESPACE
