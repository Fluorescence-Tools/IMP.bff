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
        const int nn,
        const float tile_penalty_threshold
){
    const IMP::em::DensityHeader* header = av->get_header();

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
        int dz2 = dz * dz;
        for(int y = ya; y < ye; y++) {
            int oy = y * nx;
            int dy = (y - y0);
            int dz2_dy2 = dz2 + dy * dy;
            for(int x = xa; x < xe; x++) {
                int ox = x;
                int dx = (x - x0);
                int dz2_dy2_dx2 = dz2_dy2 + dx * dx;
                PathMapTile* tile = &tiles[oz + oy + ox];
                if(tile->penalty < tile_penalty_threshold){
                    float edge_cost = std::sqrt(dz2_dy2_dx2);
                    edges.emplace_back(
                            PathMapTileEdge(tile, edge_cost)
                    );
                }
            }
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
            c = (cost >= bounds.first && cost < bounds.second) ? cost : 0.0f;
            value = c * density;
            value = clamp(value, bounds);
            break;
        case PM_TILE_PATH_LENGTH:
            c = cost * grid_spacing;
            value = (c >= bounds.first && c < bounds.second) ? c : 0.0f;
            break;
        case PM_TILE_PATH_LENGTH_DENSITY:
            c = cost * grid_spacing;
            value = (c >= bounds.first && c < bounds.second) ? c : 0.0f;
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
            value = (c >= bounds.first && c < bounds.second) ? 1.0 : 0.0f;
            value *= clamp(features[name], bounds);
            break;
        default: // case PM_TILE_COST:
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
