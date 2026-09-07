/**
 *  \file PathMap.cpp
 *  \brief Class to search path on grids
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2023 IMP Inventors. All rights reserved.
 *
 */
#include <IMP/bff/PathMap.h>
#include <IMP/core/XYZR.h>

// Only `write_map_feature` needs these, and only until it moves to the
// connection layer; `PathMap` itself no longer knows what an EM map is.
#include <IMP/em/DensityHeader.h>
#include <IMP/em/DensityMap.h>
#include <IMP/em/MRCReaderWriter.h>
#include <IMP/em/XplorReaderWriter.h>
#include <IMP/em/EMReaderWriter.h>
#include <IMP/em/SpiderReaderWriter.h>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <IMP/bff/DataPaths.h>

#include <cstring>
#ifdef __ARM_NEON
#include <arm_neon.h>
#endif

IMPBFF_BEGIN_NAMESPACE

constexpr float PathMap::BLOCKED_COST;

PathMap::PathMap(
        const PathMapHeader &av_header,
        std::string name,
        float resolution
) : DensityGrid(name), pathMapHeader_(av_header)
{
    set_name(name);
    set_path_map_header(av_header, resolution);
}

void PathMap::set_path_map_header(const PathMapHeader &av_header, float resolution)
{
    const GridHeader *nh = av_header.get_density_header();
    if(nh->get_nx() != header_.get_nx() || nh->get_ny() != header_.get_ny()
       || nh->get_nz() != header_.get_nz()){
        // neighbour offsets are linear-index deltas of the old shape
        offsets_.clear();
        nb_delta_.clear();
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
    // The location arrays must match the new shape: calc_all_voxel2loc()
    // alone is a no-op when locations were computed for the old shape (and
    // the lattice path then wrote nvox entries into arrays of the old size).
    reset_all_voxel2loc();
    calc_all_voxel2loc();
}

LinkerWeighting::LinkerWeighting(double x_min, double x_max,
                                 const std::vector<double>& y)
        : x_min_(x_min), x_max_(x_max), y_(y) {
    if (y_.size() < 2) {
        IMP_THROW("a linker weighting needs at least two tabulated points",
                  ValueException);
    }
    if (!(x_max_ > x_min_)) {
        IMP_THROW("the path-length range must be positive, got [" << x_min_
                                                                  << ", "
                                                                  << x_max_
                                                                  << "]",
                  ValueException);
    }
}

double LinkerWeighting::get_weight(double path_length) const {
    if (y_.empty()) return 1.0;
    const std::size_t n = y_.size();
    if (path_length <= x_min_) return y_.front();
    if (path_length >= x_max_) return y_.back();
    const double step = (x_max_ - x_min_) / (double) (n - 1);
    const double f = (path_length - x_min_) / step;
    std::size_t i = (std::size_t) f;
    if (i >= n - 1) return y_.back();
    const double t = f - (double) i;
    return y_[i] * (1.0 - t) + y_[i + 1] * t;
}

double LinkerWeighting::get_supported_fraction(double x) const {
    if (y_.empty()) return 1.0;
    double total = 0.0, below = 0.0;
    const double step = (x_max_ - x_min_) / (double) (y_.size() - 1);
    for (std::size_t i = 0; i < y_.size(); ++i) {
        total += y_[i];
        if (x_min_ + step * (double) i <= x) below += y_[i];
    }
    return total > 0.0 ? below / total : 0.0;
}

LinkerWeighting read_linker_weighting(const std::string& path,
                                      double linker_length) {
    std::ifstream in(path.c_str());
    if (!in) IMP_THROW("cannot read " << path, IOException);

    std::vector<double> keys;
    std::vector<double> xs;
    std::vector<std::vector<double> > rows;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ls(line);
        std::string first;
        if (!(ls >> first)) continue;
        if (keys.empty()) {
            double k;                       // the header: label then the keys
            while (ls >> k) keys.push_back(k);
            if (keys.empty()) {
                IMP_THROW(path << ": the header row carries no axis keys",
                          ValueException);
            }
            continue;
        }
        double x;
        try { x = std::stod(first); } catch (...) { continue; }
        std::vector<double> row;
        double v;
        while (ls >> v) row.push_back(v);
        if (row.size() != keys.size()) {
            IMP_THROW(path << ": row at " << x << " has " << row.size()
                           << " values for " << keys.size() << " keys",
                      ValueException);
        }
        xs.push_back(x);
        rows.push_back(row);
    }
    if (xs.size() < 2) {
        IMP_THROW(path << ": fewer than two tabulated path lengths",
                  ValueException);
    }

    // the first key not less than the linker: a longer linker than the table
    // covers gets no weighting rather than an extrapolation
    std::size_t col = keys.size();
    for (std::size_t i = 0; i < keys.size(); ++i) {
        if (keys[i] >= linker_length) { col = i; break; }
    }
    if (col == keys.size()) return LinkerWeighting();

    std::vector<double> y;
    y.reserve(rows.size());
    for (const auto& r : rows) y.push_back(r[col]);
    return LinkerWeighting(xs.front(), xs.back(), y);
}

LinkerWeighting linker_weighting(double linker_length) {
    static const std::string path =
            get_data_path("linker/chain_weighting.csv");
    return read_linker_weighting(path, linker_length);
}

void PathMap::set_linker_weighting(const LinkerWeighting& w) { weighting_ = w; }


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
    calc_all_voxel2loc();
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
    const char *interior_flags = get_interior_flags();

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
            const bool interior = interior_flags[cidx] != 0;
            int x0 = 0, y0 = 0, z0 = 0;
            if(!interior){
                x0 = (int) (cidx % nx);
                y0 = (int) ((cidx / nx) % ny);
                z0 = (int) (cidx / nxy);
            }
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

void PathMap::dijkstra_lattice(long source_idx, float max_cost){
    // `cost` holds BLOCKED_COST for blocked tiles and TILE_COST_DEFAULT for
    // open ones; the source is open. Same relaxation order and arithmetic as
    // dijkstra_bounded_core with binary penalties (open = 0), so costs are
    // identical; blocked tiles are skipped by the cost comparison itself
    // (new_cost >= 0 > BLOCKED_COST).
    long n_voxel = get_number_of_voxels();
    if(offsets_.empty() || nb_delta_.size() != offsets_.size() / 5){
        if(offsets_.empty()) offsets_ = get_neighbor_idx_offsets();
        const size_t n = offsets_.size() / 5;
        nb_delta_.resize(n); nb_len_.resize(n); nb_dz_.resize(n); nb_dy_.resize(n); nb_dx_.resize(n);
        for(size_t j = 0; j < n; j++){
            nb_dz_[j] = offsets_[5*j + 0]; nb_dy_[j] = offsets_[5*j + 1]; nb_dx_[j] = offsets_[5*j + 2];
            nb_delta_[j] = offsets_[5*j + 3];
            std::memcpy(&nb_len_[j], &offsets_[5*j + 4], sizeof(float));
        }
    }
    const size_t nnb = nb_delta_.size();
    const std::vector<long> &nb_delta = nb_delta_;
    const std::vector<float> &nb_len = nb_len_;
    const std::vector<int> &nb_dz = nb_dz_, &nb_dy = nb_dy_, &nb_dx = nb_dx_;
    const int nx = header_.get_nx();
    const int ny = header_.get_ny();
    const int nz = header_.get_nz();
    const int nxy = nx * ny;
    const char *interior_flags = get_interior_flags();
#ifdef __ARM_NEON
    // Row form of the symmetric 26 stencil for interior tiles: the three
    // neighbours dx = -1, 0, 1 of each of the 9 (dz, dy) rows are consecutive
    // in memory, so one 4-lane load covers them (lane 3 is masked with an
    // infinite edge length; the tile itself in the centre row has length 0
    // and never improves on itself). Same edge lengths and float arithmetic
    // as the offset table; relaxation order does not change costs.
    const bool rows26 = symmetric_stencil_ && nnb == 26;
    long row_base[9];
    float32x4_t row_len[9];
    if(rows26){
        int r = 0;
        for(int dz = -1; dz <= 1; dz++){
            for(int dy = -1; dy <= 1; dy++, r++){
                row_base[r] = (long) dz * nxy + (long) dy * nx - 1;
                float l[4];
                for(int k = 0; k < 3; k++){
                    int dx = k - 1;
                    l[k] = std::sqrt((float) (dx*dx + dy*dy + dz*dz));
                }
                l[3] = std::numeric_limits<float>::infinity();
                row_len[r] = vld1q_f32(l);
            }
        }
    }
    const float32x4_t vmax = vdupq_n_f32(max_cost);
#endif
    std::vector<std::vector<int> > &buckets = bucket_scratch_;
    for(auto &b : buckets) b.clear();
    if(queued_scratch_.size() != (size_t) n_voxel) queued_scratch_.assign(n_voxel, -1);
    else std::fill(queued_scratch_.begin(), queued_scratch_.end(), (int16_t) -1);
    int16_t *queued_in = queued_scratch_.data();
    auto push = [&](float c, long idx){
        int b = (int) c;
        if(queued_in[idx] == (int16_t) b) return;
        queued_in[idx] = (int16_t) b;
        if((size_t) b >= buckets.size()) buckets.resize(b + 1);
        buckets[b].push_back((int) idx);
    };
    float *cost_ptr = cost.data();
    cost_ptr[source_idx] = 0.0f;
    push(0.0f, source_idx);

    bool done = false;
    for(size_t b = 0; b < buckets.size() && !done; b++){
        const size_t nb = buckets[b].size();
        for(size_t e = 0; e < nb; e++){
            const long cidx = buckets[b][e];
            const float ccost = cost_ptr[cidx];
            if((size_t) ccost != b) continue;   // stale: moved to a lower bucket
            if(ccost >= max_cost){ done = true; break; }
            const bool interior = interior_flags[cidx] != 0;
#ifdef __ARM_NEON
            if(interior && rows26){
                const float32x4_t vc = vdupq_n_f32(ccost);
                for(int r = 0; r < 9; r++){
                    const long base = cidx + row_base[r];
                    const float32x4_t cur = vld1q_f32(cost_ptr + base);
                    const float32x4_t nw = vaddq_f32(vc, row_len[r]);
                    const uint32x4_t better = vandq_u32(vcltq_f32(nw, cur), vcltq_f32(nw, vmax));
                    if(__builtin_expect(vmaxvq_u32(better) == 0, 1)) continue;
                    float nwv[4]; vst1q_f32(nwv, nw);
                    uint32_t bv[4]; vst1q_u32(bv, better);
                    for(int k = 0; k < 3; k++){
                        if(bv[k]){
                            const long nidx = base + k;
                            cost_ptr[nidx] = nwv[k];
                            push(nwv[k], nidx);
                        }
                    }
                }
                continue;
            }
#endif
            int x0 = 0, y0 = 0, z0 = 0;
            if(!interior){
                x0 = (int) (cidx % nx);
                y0 = (int) ((cidx / nx) % ny);
                z0 = (int) (cidx / nxy);
            }
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
                const float new_cost = ccost + nb_len[j];   // open tile: penalty 0
                // Tiles that would settle at or beyond the bound never carry
                // density; they are neither written nor queued (their cost
                // stays at the default).
                if(__builtin_expect(new_cost < cost_ptr[nidx] && new_cost < max_cost, 0)){   // false for blocked (< 0)
                    cost_ptr[nidx] = new_cost;
                    push(new_cost, nidx);
                }
            }
        }
    }
    reached_valid_ = true;
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
    dijkstra_bounded_core(source_idx, -1, max_cost, penalty_soa_.data(), false);
    cost[source_idx] = TILE_COST_DEFAULT;
    if(density_soa_.size() != (size_t) n_voxel){
        density_soa_.assign(n_voxel, 1.0f);
    }
    soa_valid_ = true;
}

void PathMap::search_lattice(long source_idx, float max_cost,
                             const IMP::algebra::Vector3D &r0,
                             double block_radius, double open_radius,
                             const int32_t *occupancy){
    long n_voxel = get_number_of_voxels();
    IMP_USAGE_CHECK(source_idx >= 0 && source_idx < n_voxel,
                    "PathMap::search_lattice: invalid source index");
    const float obstacle_threshold = pathMapHeader_.get_obstacle_threshold();
    const double bsq = block_radius * block_radius;
    const double osq = open_radius * open_radius;
    const double x0 = r0[0], y0 = r0[1], z0 = r0[2];
    if(!euclidean_search_){
        // Obstacles straight into `cost`: fill_sphere(block, inverse) then
        // fill_sphere(open) then binarise -- inside the open sphere wins,
        // then beyond the block radius, then the data. Only the tiles that
        // can be inside the block sphere are tested; the template blocks the
        // rest. Locations are the same float formula as the location arrays
        // (ix * spacing + origin), computed inline. penalty_soa_ is derived
        // from `cost` when the tiles API asks (sync_tiles_from_soa).
        const std::vector<BallCandidate> &cand = get_ball_candidates(source_idx, block_radius, open_radius);
        cost = ball_cost_template_;
        const float sp = header_.get_spacing();
        const float ox = header_.get_xorigin(), oy = header_.get_yorigin(), oz = header_.get_zorigin();
        float *cost_ptr = cost.data();
        // pass 1, all candidates as if the occupancy alone decided (contiguous
        // runs, vectorisable: a count above the threshold blocks)
        for(const BallRun &run : ball_runs_){
            const int32_t *oc = occupancy + run.idx0;
            float *cc = cost_ptr + run.idx0;
            const int len = run.len;
            for(int k = 0; k < len; k++){
                cc[k] = ((double) oc[k] > obstacle_threshold) ? BLOCKED_COST : TILE_COST_DEFAULT;
            }
        }
        // pass 2, the shell tiles get the exact sphere tests
        for(const BallCandidate &c : cand){
            if(c.kind == 1) continue;
            const float xf = c.ix * sp + ox, yf = c.iy * sp + oy, zf = c.iz * sp + oz;
            double dx = (double) xf - x0, dy = (double) yf - y0, dz = (double) zf - z0;
            double d2 = dx*dx + dy*dy + dz*dz;
            double value;
            if(d2 < osq) value = 0.0;
            else if(d2 >= bsq) value = TILE_PENALTY_THRESHOLD;
            else value = (double) occupancy[c.idx];
            cost_ptr[c.idx] = (value > obstacle_threshold) ? BLOCKED_COST : TILE_COST_DEFAULT;
        }
        penalty_soa_.clear();   // derived from cost on sync
        dijkstra_lattice(source_idx, max_cost);
    } else {
        calc_all_voxel2loc();
        const float *xl = x_loc_.get();
        const float *yl = y_loc_.get();
        const float *zl = z_loc_.get();
        penalty_soa_.resize(n_voxel);
        cost.resize(n_voxel);
        for(long v = 0; v < n_voxel; v++){
            double dx = (double) xl[v] - x0, dy = (double) yl[v] - y0, dz = (double) zl[v] - z0;
            double d2 = dx*dx + dy*dy + dz*dz;
            double value;
            if(d2 < osq) value = 0.0;
            else if(d2 >= bsq) value = TILE_PENALTY_THRESHOLD;
            else value = (double) occupancy[v];
            const bool blocked = value > obstacle_threshold;
            penalty_soa_[v] = blocked ? TILE_PENALTY_DEFAULT : 0.0f;
            cost[v] = blocked ? BLOCKED_COST : TILE_COST_DEFAULT;
        }
        // Exact voxel visibility: the segment from the source point to the
        // tile centre is traversed voxel by voxel (Amanatides-Woo DDA);
        // the tile is reached iff every voxel on the way is free. Its cost
        // is the Euclidean distance to the source. Tiles are taken from
        // the ball around the source voxel; the traversal stops at the
        // first obstacle, so shadowed tiles are cheap.
        const std::vector<int> &order = get_ball_order(source_idx, block_radius);
        cost.assign(n_voxel, TILE_COST_DEFAULT);
        visited.assign(n_voxel, 0);
        const int nx = header_.get_nx(), ny = header_.get_ny(), nz = header_.get_nz();
        const long nxy = (long) nx * ny;
        const double sp = header_.get_spacing();
        const double ox = header_.get_xorigin(), oy = header_.get_yorigin(), oz = header_.get_zorigin();
        const float *pen = penalty_soa_.data();
        // source voxel indices and the source point's position inside the grid
        int si = (int) (source_idx % nx), sj = (int) ((source_idx / nx) % ny), sk = (int) (source_idx / nxy);
        visited[source_idx] = 1;
        cost[source_idx] = 0.0f;
        for(int v : order){
            if(v == source_idx) continue;
            if(pen[v] != 0.0f) continue;
            const int ti = (int) (v % nx), tj = (int) ((v / nx) % ny), tk = (int) (v / nxy);
            const double ex = (double) xl[v] - x0, ey = (double) yl[v] - y0, ez = (double) zl[v] - z0;
            const double dist = std::sqrt(ex*ex + ey*ey + ez*ez);
            const float c = (float) (dist / sp);
            if(c >= max_cost) continue;
            // DDA from the source point (x0,y0,z0) towards the tile centre.
            // Voxel boundaries lie at origin + (i +- 1/2) * spacing.
            int i = si, j = sj, k = sk;
            const int stepi = (ex > 0) - (ex < 0), stepj = (ey > 0) - (ey < 0), stepk = (ez > 0) - (ez < 0);
            const double inf = std::numeric_limits<double>::infinity();
            // parametric distance (0..1 along the segment) to the next boundary
            auto first = [&](double p0, double d, int idx, double o, int step) -> double {
                if(step == 0) return inf;
                double boundary = o + (idx + (step > 0 ? 0.5 : -0.5)) * sp;
                return (boundary - p0) / d;
            };
            double tmx = first(x0, ex, i, ox, stepi), tmy = first(y0, ey, j, oy, stepj), tmz = first(z0, ez, k, oz, stepk);
            const double tdx = stepi ? sp / std::fabs(ex) : inf, tdy = stepj ? sp / std::fabs(ey) : inf, tdz = stepk ? sp / std::fabs(ez) : inf;
            bool clear = true;
            while(i != ti || j != tj || k != tk){
                if(tmx <= tmy && tmx <= tmz){ i += stepi; tmx += tdx; }
                else if(tmy <= tmz){ j += stepj; tmy += tdy; }
                else { k += stepk; tmz += tdz; }
                if(pen[(long) k * nxy + (long) j * nx + i] != 0.0f){ clear = false; break; }
            }
            if(!clear) continue;
            cost[v] = c;
            visited[v] = 1;
        }
        reached_valid_ = true;
    }
    cost[source_idx] = TILE_COST_DEFAULT;
    if(density_soa_.size() != (size_t) n_voxel){
        density_soa_.assign(n_voxel, 1.0f);
    }
    soa_valid_ = true;
}

const std::vector<PathMap::BallCandidate> &PathMap::get_ball_candidates(long source_idx, double radius, double open_radius){
    const int nx = header_.get_nx(), ny = header_.get_ny(), nz = header_.get_nz();
    if(source_idx == ballc_source_ && radius == ballc_radius_ && open_radius == ballc_open_ &&
       nx == ballc_nx_ && ny == ballc_ny_ && nz == ballc_nz_){
        return ball_cand_;
    }
    const double h = header_.get_spacing();
    const int sx = (int) (source_idx % nx), sy = (int) ((source_idx / nx) % ny), sz = (int) (source_idx / ((long) nx * ny));
    // A tile at lattice offset o from the source voxel is between
    // (|o| - sqrt(3)/2) h and (|o| + sqrt(3)/2) h from any point inside
    // that voxel (the source itself is anywhere in the voxel).
    const double half_diag = std::sqrt(3.0) / 2.0 + 1e-9;
    const double reach = radius / h + half_diag;
    const double reach2 = reach * reach;
    const double inner_block = radius / h - half_diag;     // |o| below: always inside block sphere
    const double outer_open = open_radius / h + half_diag;  // |o| above: never inside open sphere
    ball_cand_.clear();
    ball_cost_template_.assign((size_t) nx * ny * nz, BLOCKED_COST);
    for(int z = 0; z < nz; z++){
        double dz = z - sz;
        for(int y = 0; y < ny; y++){
            double dy = y - sy;
            for(int x = 0; x < nx; x++){
                double dx = x - sx;
                double o2 = dx*dx + dy*dy + dz*dz;
                if(o2 <= reach2){
                    BallCandidate c;
                    c.idx = (int) (((long) z * ny + y) * nx + x);
                    c.ix = (int16_t) x; c.iy = (int16_t) y; c.iz = (int16_t) z;
                    double o = std::sqrt(o2);
                    c.kind = (o < inner_block && o > outer_open) ? 1 : 0;
                    ball_cand_.push_back(c);
                }
            }
        }
    }
    // runs of consecutive candidates along x
    ball_runs_.clear();
    for(size_t i = 0; i < ball_cand_.size();){
        const BallCandidate &c = ball_cand_[i];
        size_t j = i + 1;
        while(j < ball_cand_.size() && ball_cand_[j].iy == c.iy && ball_cand_[j].iz == c.iz &&
              ball_cand_[j].ix == ball_cand_[j-1].ix + 1) j++;
        BallRun r; r.idx0 = c.idx; r.ix0 = c.ix; r.iy = c.iy; r.iz = c.iz; r.len = (int16_t) (j - i);
        ball_runs_.push_back(r);
        i = j;
    }
    ballc_source_ = source_idx; ballc_radius_ = radius; ballc_open_ = open_radius;
    ballc_nx_ = nx; ballc_ny_ = ny; ballc_nz_ = nz;
    return ball_cand_;
}

const std::vector<int> &PathMap::get_ball_order(long source_idx, double radius){
    const int nx = header_.get_nx(), ny = header_.get_ny(), nz = header_.get_nz();
    if(source_idx == ball_source_ && radius == ball_radius_ &&
       nx == ball_nx_ && ny == ball_ny_ && nz == ball_nz_){
        return ball_order_;
    }
    const double h = header_.get_spacing();
    const int sx = (int) (source_idx % nx), sy = (int) ((source_idx / nx) % ny), sz = (int) (source_idx / ((long) nx * ny));
    const double reach = radius / h + std::sqrt(3.0) / 2.0 + 1e-9;
    const double reach2 = reach * reach;
    std::vector<std::pair<double, int> > cand;
    for(int z = 0; z < nz; z++){
        double dz = z - sz;
        for(int y = 0; y < ny; y++){
            double dy = y - sy;
            for(int x = 0; x < nx; x++){
                double dx = x - sx;
                double d2 = dx*dx + dy*dy + dz*dz;
                if(d2 <= reach2){
                    cand.emplace_back(d2, (int) (((long) z * ny + y) * nx + x));
                }
            }
        }
    }
    std::sort(cand.begin(), cand.end());
    ball_order_.resize(cand.size());
    for(size_t i = 0; i < cand.size(); i++) ball_order_[i] = cand[i].second;
    ball_source_ = source_idx; ball_radius_ = radius;
    ball_nx_ = nx; ball_ny_ = ny; ball_nz_ = nz;
    return ball_order_;
}

void PathMap::carve_lattice(const int32_t *occupancy){
    long n_voxel = get_number_of_voxels();
    density_soa_.resize(n_voxel);
    for(long i = 0; i < n_voxel; i++){
        double d = (double) occupancy[i];
        data_[i] = d;
        density_soa_[i] = (d > TILE_OBSTACLE_THRESHOLD) ? 0.0f : 1.0f;
    }
}

void PathMap::carve_lattice_fractional(const int32_t *const *occupancy, int n){
    long n_voxel = get_number_of_voxels();
    density_soa_.resize(n_voxel);
    if(n < 1) n = 1;
    const double inv = 1.0 / (double) n;
    for(long i = 0; i < n_voxel; i++){
        int k = 0;
        for(int j = 0; j < n; j++){
            if(!((double) occupancy[j][i] > TILE_OBSTACLE_THRESHOLD)) k++;
        }
        // data_ keeps the largest probe's counts so downstream consumers that
        // read the map data see the same thing the AV1 path gave them.
        data_[i] = (double) occupancy[0][i];
        density_soa_[i] = (float) (k * inv);
    }
}

// --------------------------------------------------------------------------
// The accessible *contact* volume (PRD-121 G9)
// --------------------------------------------------------------------------
// Distinctly named rather than file-static: src/*.cpp are compiled as one
// translation unit here, so an anonymous namespace would not isolate it.

//! Voxel offsets of the contact shell around a voxel, Olga's rule.
/*! `deltaIlist(contactR / discretizationStep, edgeL)` (`Olga/src/AV/fretAV.cpp:34`)
    takes its radius as an **`int`**, so the layer is truncated to whole voxels
    first and the shell is then `dx^2 + dy^2 + dz^2 <= delta^2` in voxel units.

    That truncation is reproduced rather than corrected, and it is worth being
    explicit about why, because it is not the more accurate rule. The contact
    layer is only half of the ACV: the other half is
    `contact_volume_trapped_fraction`, a number **fitted per site against this
    discretisation**, and a fitted share is only meaningful against the region
    it was fitted for. The seventeen sites of
    `examples/structure/T4L/fret.fps.json` carry fitted fractions of 0.33-0.72
    at `thickness = 3` and `simulation_grid_resolution = 2`, where this rule
    gives a 7-offset shell one voxel deep and a true 3 A sphere gives a
    19-offset one. Measured on 3GUN with those fitted fractions, over the 33
    pairs of `chi2_C1_33p`: this rule shortens <R_DA> by **3.04 A** on average,
    the true sphere by 2.23 A.

    The measurement that decides it is against Zenodo 3376527's published
    <R_DA> for the same structure, which this rule reproduces to a bias of
    **+0.22 A** (rmsd 0.91, r 0.996) and the true sphere to +0.49 A (rmsd 1.04)
    -- see `okf/validation/fps_screening_ab.md`. The truncation is not an
    approximation to correct; it is part of what produced the reference.

    The consequence, stated so it is not rediscovered: the contact layer is
    quantised by the grid, and `thickness` below one grid step is no layer at
    all.

 */
static void path_map_contact_offsets(double thickness, double spacing,
                                     std::vector<int> &dx,
                                     std::vector<int> &dy,
                                     std::vector<int> &dz){
    dx.clear(); dy.clear(); dz.clear();
    if(!(spacing > 0.0)) return;
    const int d = (int) (thickness / spacing);      // Olga's int conversion
    const int d2 = d * d;
    for(int k = -d; k <= d; k++){
        for(int j = -d; j <= d; j++){
            for(int i = -d; i <= d; i++){
                if(i * i + j * j + k * k > d2) continue;
                dx.push_back(i); dy.push_back(j); dz.push_back(k);
            }
        }
    }
}

void PathMap::apply_contact_weighting(const int32_t *occupancy,
                                      double thickness,
                                      double trapped_fraction){
    // Both parameters have to ask for it: `thickness <= 0` is the fps.json
    // default and a negative trapped fraction is this module's "unset".
    if(occupancy == nullptr) return;
    if(!(thickness > 0.0) || !(trapped_fraction >= 0.0)) return;
    if(!(trapped_fraction < 1.0)){
        // Olga computes volFree * f / (volTrapped * (1 - f)) and divides by
        // zero here, then writes an infinite weight into every trapped point.
        // Refusing is the only answer that leaves a usable cloud.
        IMP_WARN("contact_volume_trapped_fraction=" << trapped_fraction
                 << " is not below 1, so the free part of the volume would "
                    "carry no weight at all; the contact weighting is skipped."
                 << std::endl);
        return;
    }
    const int nx = header_.get_nx(), ny = header_.get_ny(), nz = header_.get_nz();
    const long n_voxel = (long) nx * ny * nz;
    if(n_voxel <= 0 || density_soa_.size() != (size_t) n_voxel) return;
    if(cost.size() != (size_t) n_voxel) return;

    std::vector<int> ox, oy, oz;
    path_map_contact_offsets(thickness, header_.get_spacing(), ox, oy, oz);
    if(ox.empty()) return;

    const float linker_length = pathMapHeader_.get_max_path_length();
    const float grid_spacing = pathMapHeader_.get_simulation_grid_resolution();

    // Pass one: who is in the cloud, who of them is in contact, and what the
    // two sides weigh. The membership test is get_xyz_density()'s, voxel for
    // voxel -- the normalisation has to be over the points that are actually
    // emitted, or the trapped share is a share of something else.
    std::vector<char> contact((size_t) n_voxel, 0);
    double sum_contact = 0.0, sum_free = 0.0;
    for(long i = 0; i < n_voxel; i++){
        const float c = cost[i];
        if(!(c >= 0.0f && c < TILE_COST_DEFAULT)) continue;
        if(!(c * grid_spacing < linker_length)) continue;
        const double w = (double) density_soa_[i] * (double) path_weight(i);
        if(!(w > 0.0)) continue;
        const int iz = (int) (i / ((long) nx * ny));
        const int iy = (int) ((i / nx) % ny);
        const int ix = (int) (i % nx);
        bool touched = false;
        for(size_t k = 0; k < ox.size() && !touched; k++){
            const int jx = ix + ox[k], jy = iy + oy[k], jz = iz + oz[k];
            if(jx < 0 || jx >= nx || jy < 0 || jy >= ny || jz < 0 || jz >= nz){
                continue;   // Olga's 1-D offsets wrap here; this does not
            }
            const long j = (long) jx + (long) nx * (jy + (long) ny * jz);
            if((double) occupancy[j] > TILE_OBSTACLE_THRESHOLD) touched = true;
        }
        if(touched){ contact[i] = 1; sum_contact += w; }
        else        { sum_free += w; }
    }
    // Nothing to set a ratio between: a volume entirely in contact (a buried
    // site) or entirely free (a linker that never reaches the surface) has one
    // population, and scaling it is a global factor that cancels. Olga does
    // nothing in both cases too, for the same reason.
    if(!(sum_contact > 0.0) || !(sum_free > 0.0)) return;

    const double c_scale = trapped_fraction / sum_contact;
    const double f_scale = (1.0 - trapped_fraction) / sum_free;
    for(long i = 0; i < n_voxel; i++){
        if(density_soa_[i] <= 0.0f) continue;
        density_soa_[i] = (float) (density_soa_[i] *
                                   (contact[i] ? c_scale : f_scale));
    }
}

void PathMap::carve_lattice(){
    long n_voxel = get_number_of_voxels();
    density_soa_.resize(n_voxel);
    for(long i = 0; i < n_voxel; i++){
        density_soa_[i] = (data_[i] > TILE_OBSTACLE_THRESHOLD) ? 0.0f : 1.0f;
    }
}

const char *PathMap::get_interior_flags(){
    const int nx = header_.get_nx();
    const int ny = header_.get_ny();
    const int nz = header_.get_nz();
    const int box = pathMapHeader_.get_neighbor_box_size();
    if(nx != interior_nx_ || ny != interior_ny_ || nz != interior_nz_ || box != interior_box_){
        interior_.assign((size_t) nx * ny * nz, 0);
        for(int z = box; z < nz - box; z++)
            for(int y = box; y < ny - box; y++)
                for(int x = box; x < nx - box; x++)
                    interior_[((long) z * ny + y) * nx + x] = 1;
        interior_nx_ = nx; interior_ny_ = ny; interior_nz_ = nz; interior_box_ = box;
    }
    return interior_.data();
}

void PathMap::set_origin_fast(const IMP::algebra::Vector3D &origin){
    header_.set_xorigin(origin[0]);
    header_.set_yorigin(origin[1]);
    header_.set_zorigin(origin[2]);
    header_.compute_xyz_top();
    // The location arrays are dropped and recomputed by whoever needs them
    // (calc_all_voxel2loc); the lattice path computes locations inline.
    reset_all_voxel2loc();
}

void PathMap::sync_tiles_from_soa(){
    if(!soa_valid_) return;
    long n_voxel = get_number_of_voxels();
    calc_all_voxel2loc();
    const bool have_pen = penalty_soa_.size() == (size_t) n_voxel;
    for(long i = 0; i < n_voxel; i++){
        tiles[i].penalty = have_pen ? penalty_soa_[i]
                                    : ((cost[i] < 0.0f) ? TILE_PENALTY_DEFAULT : 0.0f);
        tiles[i].cost = (cost[i] < 0.0f) ? TILE_COST_DEFAULT : cost[i];
        tiles[i].density = density_soa_[i] * path_weight(i);
        tiles[i].previous = nullptr;
    }
    // `visited` mirrors the reached set for the tile-based readers
    visited.assign(n_voxel, 0);
    for(long i = 0; i < n_voxel; i++){
        if(cost[i] >= 0.0f && cost[i] < TILE_COST_DEFAULT) visited[i] = 1;
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

void PathMap::set_particles(const IMP::ParticlesTemp &ps) {
    ps_ = ps;
    refresh_spheres_from_particles();
}

void PathMap::refresh_spheres_from_particles() {
    // The obstacles are re-read from the particles every time they are
    // sampled, not once when they are set. IMP's `SampledDensityMap` held
    // `core::XYZR` *decorators* -- views onto live particles -- so a sample
    // taken after `load_frame` moved the atoms saw where they had moved to.
    // A value copy taken at `set_particles` does not: the legacy trajectory
    // path sets the particles once and samples every frame, and with a stale
    // copy every frame was sampled against frame 0 -- AV mean positions
    // fifteen angstroms from IMP's, on the trajectory tests and nowhere else.
    if (ps_.empty()) return;
    GridSpheres spheres;
    spheres.reserve(ps_.size());
    for (std::size_t i = 0; i < ps_.size(); ++i) {
        IMP::core::XYZR x(ps_[i]);
        spheres.push_back(GridSphere(x.get_coordinates(), x.get_radius()));
    }
    set_spheres(spheres);
}

void PathMap::sample_obstacles(double extra_radius){
    refresh_spheres_from_particles();
    set_origin(pathMapHeader_.get_origin());

    // The radius each atom obstructs with: its own, or the override's when one
    // is set. **Zero is transparent** -- it is not inflated, and a binarized
    // sphere of radius zero contains no point, so the atom blocks nothing
    // while keeping its place in the list.
    const bool overridden = !obstacle_radii_.empty();
    IMP_USAGE_CHECK(!overridden || obstacle_radii_.size() == xyzr_.size(),
                    "PathMap: " << obstacle_radii_.size() << " obstacle radii "
                    "for " << xyzr_.size() << " particles");

    std::vector<double> radii_original;
    radii_original.reserve(xyzr_.size());
    for(size_t i = 0; i < xyzr_.size(); i++){
        const double own = xyzr_[i].get_radius();
        radii_original.emplace_back(own);
        const double r = overridden ? obstacle_radii_[i] : own;
        xyzr_[i].set_radius(r > 0.0 ? r + extra_radius : 0.0);
    }

    // 2. Rasterise the obstacles into the lattice
    DensityGrid::resample();

    // Restore radii
    for (size_t i = 0; i < radii_original.size(); i++) {
        xyzr_[i].set_radius(radii_original[i]);
    }
}

std::vector<IMP::algebra::Vector4D> PathMap::get_xyz_density(){
    long n_voxel = get_number_of_voxels();
    float linker_length = pathMapHeader_.get_max_path_length();
    float grid_spacing = pathMapHeader_.get_simulation_grid_resolution();

    std::vector<IMP::algebra::Vector4D> v;
    if(soa_valid_){
        // same arithmetic as PathMapTile::get_value(PM_TILE_ACCESSIBLE_DENSITY);
        // reached tiles carry a real cost (blocked ones BLOCKED_COST, open
        // unreached ones TILE_COST_DEFAULT) and lie among the ball candidates
        const float sp = header_.get_spacing();
        const float ox = header_.get_xorigin(), oy = header_.get_yorigin(), oz = header_.get_zorigin();
        const bool have_cand = ballc_source_ >= 0 && ballc_nx_ == header_.get_nx() &&
                               ballc_ny_ == header_.get_ny() && ballc_nz_ == header_.get_nz();
        if(have_cand){
            v.reserve(ball_cand_.size() / 4 + 16);
            const float *cp = cost.data();
            for(const BallRun &run : ball_runs_){
                // contiguous run: a reached tile has 0 <= cost < default
                const float *rc = cp + run.idx0;
                const int len = run.len;
                int k = 0;
#ifdef __ARM_NEON
                const float32x4_t zero = vdupq_n_f32(0.0f);
                const float32x4_t dflt = vdupq_n_f32(TILE_COST_DEFAULT);
                for(; k + 4 <= len; k += 4){
                    const float32x4_t cc = vld1q_f32(rc + k);
                    const uint32x4_t reached = vandq_u32(vcgeq_f32(cc, zero), vcltq_f32(cc, dflt));
                    if(vmaxvq_u32(reached) == 0) continue;
                    for(int q = k; q < k + 4; q++){
                        const long i = run.idx0 + q;
                        if(!(cp[i] >= 0.0f && cp[i] < TILE_COST_DEFAULT)) continue;
                        float c = cp[i] * grid_spacing;
                        float density = (c >= 0.0f && c < linker_length) ? 1.0 : 0.0f;
                        density *= density_soa_[i] * path_weight(i);
                        if(density > 0){
                            const float xf = (run.ix0 + q) * sp + ox, yf = run.iy * sp + oy, zf = run.iz * sp + oz;
                            v.emplace_back((double) xf, (double) yf, (double) zf, (double) density);
                        }
                    }
                }
#endif
                for(; k < len; k++){
                    const long i = run.idx0 + k;
                    if(!(cp[i] >= 0.0f && cp[i] < TILE_COST_DEFAULT)) continue;
                    float c = cp[i] * grid_spacing;
                    float density = (c >= 0.0f && c < linker_length) ? 1.0 : 0.0f;
                    density *= density_soa_[i] * path_weight(i);
                    if(density > 0){
                        const float xf = (run.ix0 + k) * sp + ox, yf = run.iy * sp + oy, zf = run.iz * sp + oz;
                        v.emplace_back((double) xf, (double) yf, (double) zf, (double) density);
                    }
                }
            }
            return v;
        }
        calc_all_voxel2loc();
        const float *xl = x_loc_.get();
        const float *yl = y_loc_.get();
        const float *zl = z_loc_.get();
        for(long i = 0; i < n_voxel; i++){
            if(!(cost[i] >= 0.0f && cost[i] < TILE_COST_DEFAULT)) continue;
            float c = cost[i] * grid_spacing;
            float density = (c >= 0.0f && c < linker_length) ? 1.0 : 0.0f;
            density *= density_soa_[i] * path_weight(i);
            if(density > 0){
                v.emplace_back((double) xl[i], (double) yl[i], (double) zl[i], (double) density);
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
            v.emplace_back(r[0], r[1], r[2], (double) density);
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

void PathMap::get_xyz_density_soa(std::vector<float> &vx, std::vector<float> &vy,
                                  std::vector<float> &vz, std::vector<float> &vw){
    vx.clear(); vy.clear(); vz.clear(); vw.clear();
    long n_voxel = get_number_of_voxels();
    float linker_length = pathMapHeader_.get_max_path_length();
    float grid_spacing = pathMapHeader_.get_simulation_grid_resolution();
    const bool have_cand = soa_valid_ && ballc_source_ >= 0 && ballc_nx_ == header_.get_nx() &&
                           ballc_ny_ == header_.get_ny() && ballc_nz_ == header_.get_nz();
    if(have_cand){
        // same tiles, order and arithmetic as get_xyz_density()
        const float sp = header_.get_spacing();
        const float ox = header_.get_xorigin(), oy = header_.get_yorigin(), oz = header_.get_zorigin();
        const float *cp = cost.data();
        for(const BallRun &run : ball_runs_){
            const float *rc = cp + run.idx0;
            const int len = run.len;
            int k = 0;
#ifdef __ARM_NEON
            const float32x4_t zero = vdupq_n_f32(0.0f);
            const float32x4_t dflt = vdupq_n_f32(TILE_COST_DEFAULT);
            for(; k + 4 <= len; k += 4){
                const float32x4_t cc = vld1q_f32(rc + k);
                const uint32x4_t reached = vandq_u32(vcgeq_f32(cc, zero), vcltq_f32(cc, dflt));
                if(vmaxvq_u32(reached) == 0) continue;
                for(int q = k; q < k + 4; q++){
                    const long i = run.idx0 + q;
                    if(!(cp[i] >= 0.0f && cp[i] < TILE_COST_DEFAULT)) continue;
                    float c = cp[i] * grid_spacing;
                    float density = (c >= 0.0f && c < linker_length) ? 1.0 : 0.0f;
                    density *= density_soa_[i] * path_weight(i);
                    if(density > 0){
                        vx.push_back((run.ix0 + q) * sp + ox); vy.push_back(run.iy * sp + oy);
                        vz.push_back(run.iz * sp + oz); vw.push_back(density);
                    }
                }
            }
#endif
            for(; k < len; k++){
                const long i = run.idx0 + k;
                if(!(cp[i] >= 0.0f && cp[i] < TILE_COST_DEFAULT)) continue;
                float c = cp[i] * grid_spacing;
                float density = (c >= 0.0f && c < linker_length) ? 1.0 : 0.0f;
                density *= density_soa_[i] * path_weight(i);
                if(density > 0){
                    vx.push_back((run.ix0 + k) * sp + ox); vy.push_back(run.iy * sp + oy);
                    vz.push_back(run.iz * sp + oz); vw.push_back(density);
                }
            }
        }
        return;
    }
    (void) n_voxel;
    std::vector<IMP::algebra::Vector4D> v = get_xyz_density();
    vx.reserve(v.size()); vy.reserve(v.size()); vz.reserve(v.size()); vw.reserve(v.size());
    for(const auto &p : v){
        vx.push_back((float) p[0]); vy.push_back((float) p[1]);
        vz.push_back((float) p[2]); vw.push_back((float) p[3]);
    }
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

namespace {

//! IMP's header for this grid: the seven fields the map writers read.
IMP::em::DensityHeader em_header_of(const GridHeader &gh) {
    IMP::em::DensityHeader eh;
    eh.update_map_dimensions(gh.get_nx(), gh.get_ny(), gh.get_nz());
    eh.Objectpixelsize_ = gh.get_spacing();
    eh.set_xorigin(gh.get_xorigin());
    eh.set_yorigin(gh.get_yorigin());
    eh.set_zorigin(gh.get_zorigin());
    eh.set_resolution(gh.get_resolution());
    eh.compute_xyz_top(true);
    return eh;
}

}  // namespace

IMP::em::DensityMap* PathMap::create_density_map() const {
    const GridHeader *gh = get_header();
    IMP_NEW(IMP::em::DensityMap, dm, ());
    // set_void_map allocates and zeroes; the spacing and origin go on after,
    // in the order IMP's own readers use, so the tops come out consistent.
    dm->set_void_map(gh->get_nx(), gh->get_ny(), gh->get_nz());
    dm->update_voxel_size(gh->get_spacing());
    dm->set_origin(gh->get_xorigin(), gh->get_yorigin(), gh->get_zorigin());
    dm->get_header_writable()->set_resolution(gh->get_resolution());
    std::copy(data_.begin(), data_.end(), dm->get_data());
    return dm.release();
}

void write_map_feature(
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

    // The lattice is this module's own (`GridHeader`); IMP's map writers want
    // IMP's header, so one is filled in here. This is the last thing in the
    // path-map family that still speaks `IMP.em`, and it is a leaf: it writes
    // a file and returns. It moves to the connection layer with the rest of
    // the IMP-facing surface -- the point of the change underneath it was to
    // get `PathMap` itself off `SampledDensityMap`, not to reimplement four
    // volume formats today.
    const IMP::em::DensityHeader eh = em_header_of(*d->get_header());
    rw->write(name.c_str(), f_data.data(), eh);
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

void PathMap::resize(long nvox){
    // The grid does the allocation *and* drops the location caches, so a
    // resize can never leave `x_loc_` sized for the old shape. Previously the
    // only caller (`set_path_map_header`) reset them by hand afterwards; the
    // invariant belongs where the size changes, not in every caller.
    DensityGrid::resize(nvox);
    reached_valid_ = false;
    ball_source_ = -1;
    ballc_source_ = -1;
    soa_valid_ = false;
    penalty_soa_.clear();
    density_soa_.clear();

    edge_computed.resize(0);
    edge_computed.resize(nvox, false);

    tiles.resize(0);
    for(long i = 0; i < nvox; i++){
        auto tile = PathMapTile(static_cast<int>(i));
        tiles.emplace_back(tile);
    }

}

std::vector<PathMapTile>& PathMap::get_tiles(){
    sync_tiles_from_soa();
    return tiles;
}


IMPBFF_END_NAMESPACE
