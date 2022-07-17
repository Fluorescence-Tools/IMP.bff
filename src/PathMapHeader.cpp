/**
 *  \file IMP/bff/PathMapHeader.h
 *  \brief Header class for path search class PathMap
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */
#include <IMP/bff/PathMapHeader.h>

IMPBFF_BEGIN_NAMESPACE


PathMapHeader::PathMapHeader(
        double max_path_length,
        double grid_spacing,
        int neighbor_radius,
        double obstacle_threshold
) :
        grid_spacing_(grid_spacing),
        max_path_length_(max_path_length),
        neighbor_radius_(neighbor_radius),
        obstacle_threshold_(obstacle_threshold)
{

    density_header_ = IMP::em::DensityHeader();
    density_header_.Objectpixelsize_ = grid_spacing;
    update_map_dimensions();
}

double PathMapHeader::get_grid_edge_length(){
    double ll = get_max_path_length();
    return 2. * ll;
}

void PathMapHeader::update_map_dimensions(int nx, int ny, int nz){
    if(nx < 0 || ny <0 || nz < 0){
        double grid_edge_length = get_grid_edge_length();
        double spacing = get_simulation_grid_resolution();
        int extent = std::ceil(grid_edge_length / spacing);
        nx = ny = nz = extent;
    }
    density_header_.update_map_dimensions(nx, ny, nz);
    density_header_.compute_xyz_top(true);
}


IMPBFF_END_NAMESPACE
