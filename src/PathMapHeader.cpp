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
        double neighbor_radius,
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

double PathMapHeader::get_simulation_grid_resolution() {
    grid_spacing_ = (double) density_header_.get_spacing();
    return grid_spacing_;
}

void PathMapHeader::set_obstacle_threshold(
        double obstacle_threshold){
    obstacle_threshold_ = obstacle_threshold;
}

int PathMapHeader::get_neighbor_box_size() const{
    return ceil(neighbor_radius_);
}

void PathMapHeader::set_neighbor_radius(double neighbor_radius){
    neighbor_radius_ = neighbor_radius;
}

void PathMapHeader::set_origin(float x, float y, float z){
    density_header_.set_xorigin(x);
    density_header_.set_yorigin(y);
    density_header_.set_zorigin(z);
}

void PathMapHeader::set_path_origin(const IMP::algebra::Vector3D &v){
    path_origin_ = v;
    double lr = get_grid_edge_length() / 2;
    density_header_.set_xorigin(v[0] - lr);
    density_header_.set_yorigin(v[1] - lr);
    density_header_.set_zorigin(v[2] - lr);
}

IMP::algebra::Vector3D PathMapHeader::get_origin() const {
    return algebra::Vector3D(
            (double) density_header_.get_xorigin(),
            (double) density_header_.get_yorigin(),
            (double) density_header_.get_zorigin()
    );
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
