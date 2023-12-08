/**
 * \file IMP/bff/AV.h
 * \brief Simple Accessible Volume decorator.
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */
#include <IMP/bff/AV.h>

IMPBFF_BEGIN_NAMESPACE


void AV::show(std::ostream &out) const {
  out << "(" << IMP::algebra::commas_io(get_parameter()) << ")";
}


std::string AVPairDistanceMeasurement::get_json(){
    // create an empty structure (null)
    nlohmann::json j;
    j["position1_name"] = position_1;
    j["position2_name"] = position_2;
    j["distance"] = distance;
    j["error_neg"] = error_neg;
    j["error_pos"] = error_pos;
    j["Forster_radius"] = forster_radius;
    j["distance_type"] = distance_type;
    return j.dump();
}


double AVPairDistanceMeasurement::score_model(double model){
    auto ev = [](auto f, auto m, auto en, auto ep){
        auto dev = m - f;
        auto w = (dev < 0) ? 1. / en : 1. / ep;
        return .5 * algebra::get_squared(dev * w);
    };
    if(std::isnan(model)){
        return std::numeric_limits<double>::infinity();
    }
    else{
        return 0.5 * ev(model, distance, error_neg, error_pos);
    }
}


double av_distance(
        const IMP::bff::AV& av1,
        const IMP::bff::AV& av2,
        double forster_radius,
        int distance_type,
        int n_samples
){
    // Draw points using Inverse transform sampling
    auto el3getter = [](const IMP::algebra::Vector4D &p) { return p[3]; };
    using points_type = std::vector<IMP::algebra::Vector4D>;
    auto m1 = av1.get_map();
    auto m2 = av2.get_map();
    auto p1 = m1->get_xyz_density();
    auto p2 = m2->get_xyz_density();
    if(!p1.empty() && !p2.empty()){
        InverseSampler<points_type> sampler1(p1, el3getter);
        InverseSampler<points_type> sampler2(p2, el3getter);
        double val = 0.;
        switch(distance_type){
            case DYE_PAIR_EFFICIENCY: {
                for (int s = 0; s < n_samples; s++) {
                    auto tmp = sampler1.get_random() - sampler2.get_random();
                    tmp[3] = 0.0;
                    val += fret_efficiency<double>(tmp.get_magnitude(), forster_radius);
                }
                return val / n_samples;
            }
            case DYE_PAIR_DISTANCE_E: {
                double fret_eff = av_distance(av1, av2, forster_radius, DYE_PAIR_EFFICIENCY, n_samples);
                return distance_fret<double>(fret_eff, forster_radius);
            }
            case DYE_PAIR_DISTANCE_MP: {
                IMP::algebra::Vector3D mp1 = av1.get_mean_position();
                IMP::algebra::Vector3D mp2 = av2.get_mean_position();
                return get_l2_norm((mp1 - mp2));
            }
            case DYE_PAIR_XYZ_DISTANCE: {
                IMP::Particle* p1 = av1.get_particle();
                IMP::Particle* p2 = av2.get_particle();
                IMP::algebra::Vector3D mp1 = IMP::core::XYZ(p1).get_coordinates();
                IMP::algebra::Vector3D mp2 = IMP::core::XYZ(p2).get_coordinates();
                return get_l2_norm((mp1 - mp2));
            }
            case DYE_PAIR_DISTANCE_MEAN:
            default: {
                for (int s = 0; s < n_samples; s++) {
                    auto tmp = sampler1.get_random() - sampler2.get_random();
                    tmp[3] = 0.0;
                    val += tmp.get_magnitude();
                }
                return val / n_samples;
            }
        }
    } else {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

IMP::bff::PathMap* AV::get_map() const{
    // get_map needs to be const
    if(av_map_ == nullptr){
        // cast away const to init map ¯\_(ツ)_/¯
        AV* ptr = (AV*)(this);
        ptr->init_path_map();
        ptr->resample();
    }
    return av_map_;
}

IMP::algebra::Vector3D AV::get_mean_position(bool include_source) const{
    IMP::algebra::Vector3D r = {0.0, 0.0, 0.0};
    double sum = 1.0;
    if(include_source){
        r += get_source_coordinates();
        sum += 1.0;
    }
    auto xyzd = get_map()->get_xyz_density();
    for(auto &a: xyzd){
        if(a[3] <= 0.0f) continue;
        sum += a[3];
        r[0] += a[0] * a[3];
        r[1] += a[1] * a[3];
        r[2] += a[2] * a[3];
    }
    return r /= sum;
}

IMP::Particle* AV::get_source() const{
    return get_model()->get_particle(get_particle_index(0));
}

IMP::algebra::Vector3D AV::get_source_coordinates() const{
    auto xyz = IMP::core::XYZ(get_source());
    return xyz.get_coordinates();
}

void AV::init_path_map(){
    auto path_map_header = create_path_map_header();
    av_map_ = new IMP::bff::PathMap(path_map_header);
    IMP::Particle* parent = get_model()->get_particle(get_particle_index(0));

    auto h = IMP::atom::Hierarchy(get_model(), parent->get_index());
    auto root = IMP::atom::get_root(h);
    av_map_->set_particles(get_leaves(root));
}

void AV::resample(bool shift_xyz){
    auto map = get_map();

    // Update parameters of path map
    if(get_parameters_are_optimized()){
        auto path_map_header = create_path_map_header();
        map->set_path_map_header(path_map_header);
    }

    // Update path_map origin
    IMP::algebra::Vector3D source = get_source_coordinates();
    auto header = get_map()->get_path_map_header_writable();
    header->set_path_origin(source);
    av_map_->set_origin(source);

    // 1.1 Sample obstacles
    map->sample_obstacles(get_linker_width() * 0.5);

    // 2. Block voxels further away from source than linker length
    double critical_radius;
    critical_radius = get_linker_length();
    map->fill_sphere(source, critical_radius, TILE_PENALTY_THRESHOLD, true);

    // 3.1 Unblock voxels in initial sphere
    critical_radius = get_allowed_sphere_radius();
    map->fill_sphere(source, critical_radius, 0, false);

    // 4. Find a path from source to other tiles
    map->update_tiles(); // Update tiles to assure that the nodes are updated
    long source_idx = map->get_voxel_by_location(source);
    map->find_path_dijkstra(source_idx, -1);

    // 5. Remove tiles closer to obstacles than dye radius
    double r = get_radius1();
    map->sample_obstacles(r);
    auto obstacle = map->get_data();
    long nvox = map->get_number_of_voxels();
    for(long i=0; i<nvox; i++){
        if(obstacle[i] > TILE_OBSTACLE_THRESHOLD){
            // Tile density is by default 1 (or another user specified number)
            // Setting the tile density to zeros effectively removes the tile.
            map->tiles[i].density *= 0.0;
        }
    }

    // Shift XYZ to mean AV position
    if(shift_xyz){
        set_coordinates(get_mean_position());
    }
}

void AV::set_av_parameter(const nlohmann::json &j){
    set_linker_length(j.value("linker_length", 20.0));
    algebra::Vector3D r = {j.value("radius1", 3.0),
                           j.value("radius2", 0.0),
                           j.value("radius3", 0.0)};
    set_radius1(r[0]);
    set_radius2(r[0]);
    set_radius3(r[0]);
    set_linker_width(j.value("linker_width", 0.5));
    set_allowed_sphere_radius(j.value("allowed_sphere_radius", 1.5));
    set_contact_volume_thickness(j.value("contact_volume_thickness", 0.0));
    set_contact_volume_trapped_fraction(j.value("contact_volume_trapped_fraction", -1));
    set_simulation_grid_resolution(j.value("simulation_grid_resolution", 1.5));
}

IMP::bff::PathMapHeader AV::create_path_map_header(){
    // create PathMapHeader
    double ll = get_linker_length();
    double dg = get_simulation_grid_resolution();
    IMP::bff::PathMapHeader path_map_header(ll, dg);
    return path_map_header;
}

IMP::ParticleIndex search_labeling_site(
        const IMP::core::Hierarchy& hier,
        std::string json_str,
        const nlohmann::json &json_data
){
    nlohmann::json js = json_data;
    if(json_str != ""){
        std::cout << json_str << std::endl;
        js = nlohmann::json::parse(json_str);
    }

    auto sel = IMP::atom::Selection(hier);
    std::string chain_id = js.value("chain_identifier", "");
    int residue_seq_number = js.value("residue_seq_number", -1);
    std::string atom_name = js.value("atom_name", "");

    std::string s = "search_labeling_site:" + std::string(chain_id) + ":" + std::to_string(residue_seq_number) + ":" + std::string(atom_name);
    IMP::add_to_log(PROGRESS, s);

    if(!chain_id.empty()){
        sel.set_chain_id(chain_id);
    }
    auto p_chain = sel.get_selected_particles(false);

    if(residue_seq_number > 0){
        sel.set_residue_index(residue_seq_number);
    }
    auto p_residue = sel.get_selected_particles(false);

    if(!atom_name.empty()){
        sel.set_atom_type(IMP::atom::AtomType(atom_name));
    }
    auto p_atom = sel.get_selected_particles(false);
    std::clog << p_chain << p_residue << p_atom << std::endl;

    auto p = p_residue[0];
    IMP_USAGE_CHECK(p_residue.size() == 1, "No or ambiguous AV attachment site selected");
    if(p_atom.empty()){
        IMP_WARN("AV attached to atom but residue level resolution\n")
    } else{
        p = p_atom[0];
    }
    return p->get_index();
}

//! Random sampling over AV
    void get_xyz_density();

std::vector<double> av_random_points(const AV& av, int n_samples){
    auto m = av.get_map();
    auto d = m->get_xyz_density();
    std::vector<double> data; 
    if(!d.empty()){
        // Draw points using Inverse transform sampling
        using points_type = std::vector<IMP::algebra::Vector4D>;
        auto el3getter = [](const IMP::algebra::Vector4D &p) { return p[3]; };
        InverseSampler<points_type> sampler(d, el3getter);
        data.reserve(4 * n_samples);
        for (int s = 0; s < n_samples; s++) {
            auto v = sampler.get_random();
            data.emplace_back(v[0]);
            data.emplace_back(v[1]);
            data.emplace_back(v[2]);
            data.emplace_back(v[3]);
        }
    }

    return data;    
}

std::vector<double> av_random_distances(
        const AV& av1,
        const AV& av2,
        int n_samples
){
    auto p1 = av_random_points(av1, n_samples);
    auto p2 = av_random_points(av2, n_samples);
    std::vector<double> data; data.reserve(n_samples);
    for (int s = 0; s < n_samples; s++) {
        auto dx = p1[s * 4 + 0] - p2[s * 4 + 0];
        auto dy = p1[s * 4 + 1] - p2[s * 4 + 1];
        auto dz = p1[s * 4 + 2] - p2[s * 4 + 2];
        auto d2 = dx*dx + dy*dy + dz*dz;
        data.emplace_back(sqrt(d2));
    }
    return data;
}

std::vector<double> av_distance_distribution(
        const AV& av1,
        const AV& av2,
        //double start, double stop, int n_bins,
        std::vector<double> axis,
        int n_samples
){
    auto data = av_random_distances(av1, av2, n_samples);
 
//    // For future versions (requires C++14)
//    using namespace boost::histogram; // strip the boost::histogram prefix
//    auto hist_axis = axis::regular<>(n_bins, start, stop);
//    auto h = make_histogram(hist_axis);
//    std::for_each(data.begin(), data.end(), std::ref(h));
//    auto hist = std::vector<double>(n_bins);
//    for(int i = 0; i < n_bins; i++){
//        hist[i] = h[i];
//    }

    auto hist = std::vector<double>(axis.size(), 0);
    IMP::bff::histogram1D<double>(
            data.data(), data.size(),  // sampled distances
            nullptr, 0,                          // weights
            axis.data(), axis.size(),            // axis
            hist.data(), hist.size(),            // data
            IMP::bff::AXIS_LIN,                  // axis type
            false
    );

    return hist;
}


IMPBFF_END_NAMESPACE
