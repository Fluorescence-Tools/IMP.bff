/**
 *  \file IMP/bff/PathMapTile.h
 *  \brief Tile used in path search by PathMap
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */
#include <IMP/bff/PathMapTile.h>

IMPBFF_BEGIN_NAMESPACE


void PathMapTile::update_edges(
        PathMap* av,
        std::vector<PathMapTile> &tiles,
        const double nr,
        const float tile_penalty_threshold
){
    const IMP::em::DensityHeader* header = av->get_header();
    const int nn = ceil(nr);
    const double nr2 = nr * nr;

    int nx = header->get_nx();
    int ny = header->get_ny();
    int nz = header->get_nz();

    int x0 = av->get_dim_index_by_voxel(idx, 0);
    int y0 = av->get_dim_index_by_voxel(idx, 1);
    int z0 = av->get_dim_index_by_voxel(idx, 2);

    int za = std::max(0,  z0 - nn);
    int ze = std::min(nz, z0 + nn);
    int ya = std::max(0,  y0 - nn);
    int ye = std::min(ny, y0 + nn);
    int xa = std::max(0,  x0 - nn);
    int xe = std::min(nx, x0 + nn);
    int nx_ny = nx * ny;

    edges.clear();
    for(int z = za; z < ze; z++) {  // z slowest
        int oz = z * nx_ny;
        int dz = (z - z0);
        double dz2 = dz * dz;
        for(int y = ya; y < ye; y++) {
            int oy = y * nx;
            int dy = (y - y0);
            double dz2_dy2 = dz2 + dy * dy;
            for(int x = xa; x < xe; x++) {
                int ox = x;
                int dx = (x - x0);
                double dz2_dy2_dx2 = dz2_dy2 + dx * dx;
                if(dz2_dy2_dx2 <= nr2){
                    int tile_idx = oz + oy + ox;
                    std::cout << idx << ":"<< tile_idx << ":" << idx - tile_idx << std::endl;
                    PathMapTile* tile = &tiles[tile_idx];
                    if(tile->penalty < tile_penalty_threshold){
                        float edge_cost = std::sqrt(dz2_dy2_dx2);
                        edges.emplace_back(
                                PathMapTileEdge(tile_idx, edge_cost)
                        );
                    }
                }
            }
        }
    }
}


void PathMapTile::update_edges_2(
    int nx, int ny, int nz,
    std::vector<PathMapTile> &tiles,
    std::vector<int> neighbor_idxs,
    float tile_penalty_threshold
){

    int x0 = idx % nx;
    int y0 = idx / nx % ny;
    int z0 = idx / (nx * ny);    

    edges.clear();
    for(size_t i = 0; i < neighbor_idxs.size(); i += 5){
        int iz = z0 + neighbor_idxs[i + 0];
        int iy = y0 + neighbor_idxs[i + 1];
        int ix = x0 + neighbor_idxs[i + 2];

        if(iz >= nz || iz < 0) continue;
        if(iy >= nz || iy < 0) continue;
        if(ix >= nz || ix < 0) continue;

        int tile_idx = idx + neighbor_idxs[i + 3];
        PathMapTile* tile = &tiles[tile_idx];
        if(tile->penalty < tile_penalty_threshold){
            // edge_cost is a float stored in an 32bit int
            float edge_cost = *(float*)&neighbor_idxs[i + 4];
            edges.emplace_back(PathMapTileEdge(tile_idx, edge_cost));
        }

    }
}

std::vector<long> PathMapTile::backtrack_to_path(){
    PathMapTile* current = this;
    std::vector<long> path;
    while(current != nullptr){
        path.emplace_back(current->idx);
        current = current->previous;
    }
    std::reverse(path.begin(), path.end());
    return path;
}

float PathMapTile::get_value(
        const int value_type,
        std::pair<float, float> bounds,
        const std::string &name,
        float grid_spacing
){
    float value, c;
    auto clamp = [](float v, std::pair<float, float>b)
        {return std::min(std::max(v, b.first), b.second);};
    switch (value_type) {
        case PM_TILE_PENALTY:
            value = clamp(penalty, bounds);
            break;
        case PM_TILE_DENSITY:
            value = clamp(density, bounds);
            break;
        case PM_TILE_COST_DENSITY:
            c = (cost >= bounds.first && cost < bounds.second) ? cost :0.0f;
            value = c * density;
            value = clamp(value, bounds);
            break;
        case PM_TILE_PATH_LENGTH:
            c = cost * grid_spacing;
            value = (c >= bounds.first && c < bounds.second) ? c : -1.0f;
            break;
        case PM_TILE_PATH_LENGTH_DENSITY:
            c = cost * grid_spacing;
            value = (c >= bounds.first && c < bounds.second) ? c : -1.0f;
            value *= density;
            break;
        case PM_TILE_ACCESSIBLE_DENSITY:
            c = cost * grid_spacing;
            value = (c >= bounds.first && c < bounds.second) ? 1.0 : 0.0f;
            value *= density;
            break;
        case PM_TILE_FEATURE:
            value = clamp(features[name], bounds);
            break;
        case PM_TILE_ACCESSIBLE_FEATURE:
            c = cost * grid_spacing;
            value = (c >= bounds.first && c < bounds.second) ? 1.0 : -1.0f;
            value *= clamp(features[name], bounds);
            break;
        case PM_TILE_COST:
        default:
            value = clamp(cost, bounds);
            break;
    }
    return value;
}

void PathMapTile::set_value(int value_type, float value, const std::string &name){
    IMP_USAGE_CHECK(value_type != PM_TILE_COST_DENSITY, "Cannot set combined features.");
    switch (value_type) {
        case PM_TILE_PENALTY:
            penalty = value;
            break;
        case PM_TILE_COST:
            cost = value;
            break;
        case PM_TILE_FEATURE:
            features.insert({name, value});
            break;
        default:
            cost = value;
            break;
    }
}


IMPBFF_END_NAMESPACE
