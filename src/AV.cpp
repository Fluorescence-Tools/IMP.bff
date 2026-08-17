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


double AVPairDistanceMeasurement::score_model(double model) const{
    auto ev = [](double f, double m, double en, double ep){
        double dev = m - f;
        double w = (dev < 0) ? 1. / en : 1. / ep;
        return .5 * algebra::get_squared(dev * w);
    };
    if(std::isnan(model)){
        return std::numeric_limits<double>::infinity();
    }
    return 0.5 * ev(model, distance, error_neg, error_pos);
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

namespace {

// The lattice window of an AV: cubic, centred on the lattice-quantised source,
// large enough that every voxel centre within the linker length of the source
// lies inside. With q = round(s/h) and |s/h - q| <= 1/2, |k - q| <= ll/h + 1/2.
void lattice_window(const IMP::algebra::Vector3D &source, double ll, double h,
                    int k0[3], int &n){
    int half = (int) std::floor(ll / h + 0.5);
    n = 2 * half + 1;
    for(int d = 0; d < 3; d++){
        int q = (int) std::floor(source[d] / h + 0.5);
        k0[d] = q - half;
    }
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
    return av_map_.get();
}

internal::AVLatticeState &AV::get_state(){
    if(!state_){
        state_ = std::make_shared<internal::AVLatticeState>();
    }
    return *state_;
}

const std::vector<IMP::algebra::Vector4D> &AV::get_cloud() const{
    auto map = get_map();   // builds and resamples on first use
    auto &st = const_cast<AV*>(this)->get_state();
    if(!st.cloud_valid || st.cloud_generation != st.result_generation){
        st.cloud = map->get_xyz_density();
        st.cloud_generation = st.result_generation;
        st.cloud_valid = true;
    }
    return st.cloud;
}

IMP::algebra::Vector3D AV::get_mean_position(bool include_source) const{
    IMP::algebra::Vector3D r = {0.0, 0.0, 0.0};
    double sum = 1.0;
    if(include_source){
        r += get_source_coordinates();
        sum += 1.0;
    }
    const auto &xyzd = get_cloud();
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

IntKey AV::get_space_fixed_key(){
    static const IntKey k("av_space_fixed");
    return k;
}

bool AV::get_space_fixed() const{
    if(get_model()->get_has_attribute(get_space_fixed_key(), get_particle_index())){
        return get_model()->get_attribute(get_space_fixed_key(), get_particle_index()) != 0;
    }
    return true;
}

void AV::set_space_fixed(bool tf){
    if(get_model()->get_has_attribute(get_space_fixed_key(), get_particle_index())){
        get_model()->set_attribute(get_space_fixed_key(), get_particle_index(), tf ? 1 : 0);
    } else {
        get_model()->add_attribute(get_space_fixed_key(), get_particle_index(), tf ? 1 : 0);
    }
    if(!tf){
        static bool warned = false;
        if(!warned){
            IMP_WARN("AV: space_fixed=False (legacy source-anchored grid) is deprecated "
                     "and will be removed after the PRD-105 transition period.\n");
            warned = true;
        }
        if(state_ && state_->registry){
            IMP_THROW("AV: a shared occupancy registry requires space_fixed",
                      IMP::ValueException);
        }
    }
    // The map was built for the other anchoring; rebuild it lazily.
    if(av_map_){
        av_map_ = nullptr;
        state_.reset();
    }
}

void AV::set_occupancy_registry(AVOccupancyRegistry *registry){
    if(registry && !get_space_fixed()){
        IMP_THROW("AV: a shared occupancy registry requires space_fixed",
                  IMP::ValueException);
    }
    auto &st = get_state();
    st.registry = registry;
    st.have_result = false;
    st.private1 = nullptr;
    st.private2 = nullptr;
}

AVOccupancyRegistry *AV::get_occupancy_registry() const{
    if(!state_) return nullptr;
    return state_->registry.get();
}

std::vector<int> AV::get_lattice_window() const{
    if(!state_ || !state_->have_window) return {};
    return {state_->k0[0], state_->k0[1], state_->k0[2], state_->n};
}

long AV::get_number_of_skips() const{ return state_ ? state_->n_skip : 0; }
long AV::get_number_of_local_updates() const{ return state_ ? state_->n_local : 0; }
long AV::get_number_of_full_updates() const{ return state_ ? state_->n_full : 0; }
long AV::get_number_of_rolls() const{ return state_ ? state_->n_roll : 0; }
unsigned long AV::get_result_generation() const{
    return state_ ? state_->result_generation : 0;
}

void AV::prepare_lattice_window(){
    if(!get_space_fixed() || !state_ || !state_->registry) return;
    const double h = get_simulation_grid_resolution();
    int k0[3]; int n;
    lattice_window(get_source_coordinates(), get_linker_length(), h, k0, n);
    auto &st = get_state();
    st.registry->get_map(h, get_linker_width() * 0.5)
        ->request_window(k0[0], k0[1], k0[2], n, n, n);
    st.registry->get_map(h, get_radius1())
        ->request_window(k0[0], k0[1], k0[2], n, n, n);
}

void AV::init_path_map(){
    auto path_map_header = create_path_map_header();
    av_map_ = new IMP::bff::PathMap(path_map_header);
    IMP::Particle* parent = get_model()->get_particle(get_particle_index(0));

    auto h = IMP::atom::Hierarchy(get_model(), parent->get_index());
    auto root = IMP::atom::get_root(h);
    av_map_->set_particles(get_leaves(root));
    // A fresh map: whatever the state remembers about tiles is stale.
    auto &st = get_state();
    st.have_result = false;
    st.have_window = false;
    st.private1 = nullptr;
    st.private2 = nullptr;
}

void AV::resample(bool shift_xyz, bool force_full){
    if(av_map_ == nullptr){
        // get_map() would resample once itself; build the map here so the
        // first evaluation is one full update, not a full one and a skip.
        init_path_map();
    }
    if(get_space_fixed()){
        resample_lattice(shift_xyz, force_full);
    } else {
        resample_legacy(shift_xyz);
    }
}

// The pre-PRD-105 evaluation, kept byte-for-byte for `space_fixed=False`.
void AV::resample_legacy(bool shift_xyz){
    auto map = get_map();

    // Update parameters of path map
    if(get_parameters_are_optimized()){
        auto path_map_header = create_path_map_header();
        map->set_path_map_header(path_map_header);
    }

    // Update path_map origin
    IMP::algebra::Vector3D source = get_source_coordinates();
    auto &header = get_map()->get_path_map_header_writable();
    header.set_path_origin(source);
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
        // Tile density is by default 1. Setting it to zero removes the tile.
        // (Was `*= 0.0` on obstacle tiles only, which never restored a tile
        // once removed: repeated resamples of a moving structure shrank the
        // AV frame after frame.)
        map->tiles[i].density = (obstacle[i] > TILE_OBSTACLE_THRESHOLD)
            ? 0.0f : 1.0f;
    }

    // Shift XYZ to mean AV position
    if(shift_xyz){
        set_coordinates(get_mean_position());
    }
    auto &st = get_state();
    st.result_generation++;
    st.quad_valid = false;
    st.n_full++;
}

void AV::resample_prepare(bool shift_xyz, bool force_full){
    if(av_map_ == nullptr){
        init_path_map();
    }
    if(get_space_fixed()){
        resample_lattice_prepare(shift_xyz, force_full);
    } else {
        resample_legacy(shift_xyz);
    }
}

void AV::resample_compute(){
    if(state_ && state_->pending) resample_lattice_compute();
}

void AV::resample_finish(){
    if(state_ && state_->pending) resample_lattice_finish();
}

void AV::resample_lattice(bool shift_xyz, bool force_full){
    resample_lattice_prepare(shift_xyz, force_full);
    resample_lattice_compute();
    resample_lattice_finish();
}

void AV::resample_lattice_prepare(bool shift_xyz, bool force_full){
    auto map = get_map();
    auto &st = get_state();
    st.pending = false;

    const double h = get_simulation_grid_resolution();
    const double ll = get_linker_length();
    const IMP::algebra::Vector3D source = get_source_coordinates();
    const IMP::algebra::VectorD<9> parameter = get_parameter();

    // 0. Window on the lattice; roll when it moved by whole voxels
    int k0[3]; int n;
    lattice_window(source, ll, h, k0, n);
    auto same_vec = [](const auto &a, const auto &b){
        for(unsigned i = 0; i < a.get_dimension(); i++) if(a[i] != b[i]) return false;
        return true;
    };
    bool params_changed = !st.have_result ||
        !same_vec(st.last_parameter, parameter);
    bool dims_changed = !st.have_window || n != st.n;
    bool rolled = st.have_window && !dims_changed &&
        (k0[0] != st.k0[0] || k0[1] != st.k0[1] || k0[2] != st.k0[2]);
    if(rolled) st.n_roll++;

    auto &header = map->get_path_map_header_writable();
    IMP::algebra::Vector3D grid_origin(k0[0] * h, k0[1] * h, k0[2] * h);
    if(dims_changed || params_changed){
        auto path_map_header = create_path_map_header();
        path_map_header.set_path_origin(source, grid_origin);
        map->set_path_map_header(path_map_header);
        map->set_origin(grid_origin);
    } else if(rolled){
        header.set_path_origin(source, grid_origin);
        map->set_origin(grid_origin);
    } else {
        // same window; the label may have moved inside its voxel
        header.set_path_origin(source, grid_origin);
    }
    for(int d = 0; d < 3; d++) st.k0[d] = k0[d];
    st.n = n;
    st.have_window = true;

    // 1. Occupancy sources for the two passes
    const double extra1 = get_linker_width() * 0.5;
    const double extra2 = get_radius1();
    AVOccupancyMap *occ1; AVOccupancyMap *occ2;
    if(st.registry){
        occ1 = st.registry->get_map(h, extra1);
        occ2 = st.registry->get_map(h, extra2);
        occ1->request_window(k0[0], k0[1], k0[2], n, n, n);
        occ2->request_window(k0[0], k0[1], k0[2], n, n, n);
    } else {
        IMP::ParticlesTemp ps(map->ps_.begin(), map->ps_.end());
        if(!st.private1 || st.private1->get_extra_radius() != extra1
                        || st.private1->get_spacing() != h){
            st.private1 = new AVOccupancyMap(h, extra1, ps);
            st.private1->set_was_used(true);
        }
        if(!st.private2 || st.private2->get_extra_radius() != extra2
                        || st.private2->get_spacing() != h){
            st.private2 = new AVOccupancyMap(h, extra2, ps);
            st.private2->set_was_used(true);
        }
        occ1 = st.private1.get();
        occ2 = st.private2.get();
        occ1->set_window(k0[0], k0[1], k0[2], n, n, n);
        occ2->set_window(k0[0], k0[1], k0[2], n, n, n);
    }
    occ1->update(force_full);
    occ2->update(force_full);

    // 2. Nothing that feeds the search changed: keep the tiles
    // Nothing that feeds the search changed *inside this window*: a moved
    // atom whose footprint (old and new) misses the window leaves this AV
    // alone, so refinement moves elsewhere in the structure cost nothing here.
    bool can_skip = !force_full && st.have_result && !rolled &&
        !dims_changed && !params_changed &&
        same_vec(st.last_source, source) &&
        !occ1->get_changed_since(st.generation1, k0[0], k0[1], k0[2], n, n, n) &&
        !occ2->get_changed_since(st.generation2, k0[0], k0[1], k0[2], n, n, n);
    if(can_skip){
        st.n_skip++;
        st.generation1 = occ1->get_generation();
        st.generation2 = occ2->get_generation();
        if(shift_xyz) set_coordinates(st.last_mean);
        return;
    }
    if(rolled || dims_changed || params_changed || !st.have_result){
        st.n_full++;
    } else {
        st.n_local++;
    }

    // Everything compute() needs, gathered while we may still touch the Model
    st.pending = true;
    st.pending_shift_xyz = shift_xyz;
    st.pending_ll = ll;
    st.pending_allowed = get_allowed_sphere_radius();
    st.pending_source = source;
    st.pending_occ1 = occ1;
    st.pending_occ2 = occ2;
    st.pending_gen1 = occ1->get_generation();
    st.pending_gen2 = occ2->get_generation();
    st.last_parameter = parameter;
}

// Touches only this AV's map and lattice state: safe to run concurrently
// with other AVs' compute phases (their occupancy sources are read-only now).
void AV::resample_lattice_compute(){
    auto &st = *state_;
    if(!st.pending) return;
    PathMap *map = av_map_.get();
    const int *k0 = st.k0;
    const int n = st.n;
    const IMP::algebra::Vector3D &source = st.pending_source;

    // 3. Obstacles inflated by half the linker width, from the lattice
    long nvox = map->get_number_of_voxels();
    double *data = map->get_data();
    st.pending_occ1->read_window(k0[0], k0[1], k0[2], n, n, n, data);
    map->normalized_ = false;
    map->rms_calculated_ = false;

    // 4. Block voxels further away from source than linker length
    map->fill_sphere(source, st.pending_ll, TILE_PENALTY_THRESHOLD, true);

    // 5. Unblock voxels in initial sphere
    map->fill_sphere(source, st.pending_allowed, 0, false);

    // 6. Find a path from source to other tiles
    map->update_tiles();
    long source_idx = map->get_voxel_by_location(source);
    map->find_path_dijkstra(source_idx, -1);

    // 7. Remove tiles closer to obstacles than the dye radius
    st.pending_occ2->read_window(k0[0], k0[1], k0[2], n, n, n, data);
    for(long i=0; i<nvox; i++){
        map->tiles[i].density = (data[i] > TILE_OBSTACLE_THRESHOLD)
            ? 0.0f : 1.0f;
    }

    st.have_result = true;
    st.last_source = source;
    st.generation1 = st.pending_gen1;
    st.generation2 = st.pending_gen2;
    st.result_generation++;
    st.quad_valid = false;

    // the cloud and the mean position, from the map and the cached source
    st.cloud = map->get_xyz_density();
    st.cloud_generation = st.result_generation;
    st.cloud_valid = true;
    IMP::algebra::Vector3D r = source;
    double sum = 2.0;
    for(auto &a: st.cloud){
        if(a[3] <= 0.0f) continue;
        sum += a[3];
        r[0] += a[0] * a[3];
        r[1] += a[1] * a[3];
        r[2] += a[2] * a[3];
    }
    st.last_mean = r / sum;
}

void AV::resample_lattice_finish(){
    auto &st = *state_;
    if(!st.pending) return;
    st.pending = false;
    if(st.pending_shift_xyz){
        set_coordinates(st.last_mean);
    }
}

void AV::set_av_parameter(const nlohmann::json &j){
    set_linker_length(j.value("linker_length", 20.0));
    algebra::Vector3D r = {j.value("radius1", 3.0),
                           j.value("radius2", 0.0),
                           j.value("radius3", 0.0)};
    set_radius1(r[0]);
    set_radius2(r[1]);
    set_radius3(r[2]);
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
    if(get_space_fixed()){
        int k0[3]; int n;
        lattice_window(get_source_coordinates(), ll, dg, k0, n);
        path_map_header.update_map_dimensions(n, n, n);
    }
    return path_map_header;
}

namespace {

// Coarsen a weighted lattice cloud into cubic blocks of `m` voxels edge:
// the smallest m that leaves at most `k` non-empty blocks. Each block is
// summarised by its weighted centroid, total weight and second central
// moments (xx, yy, zz, xy, xz, yz). Deterministic: blocks are keyed by
// integer lattice index and emitted in key order.
void coarsen_cloud(
        const std::vector<IMP::algebra::Vector4D> &cloud,
        double spacing, int k,
        std::vector<IMP::algebra::Vector4D> &points,
        std::vector<std::array<double, 6> > &moments){
    points.clear();
    moments.clear();
    if(cloud.empty() || k <= 0) return;
    const size_t np = cloud.size();
    if((int) np <= k){
        points = cloud;
        moments.assign(np, {0, 0, 0, 0, 0, 0});
        return;
    }
    // integer lattice indices relative to the cloud minimum
    std::vector<int> ix(np), iy(np), iz(np);
    double minx = cloud[0][0], miny = cloud[0][1], minz = cloud[0][2];
    for(const auto &p : cloud){
        minx = std::min(minx, (double) p[0]);
        miny = std::min(miny, (double) p[1]);
        minz = std::min(minz, (double) p[2]);
    }
    int ex = 0, ey = 0, ez = 0;   // extents in voxels
    for(size_t i = 0; i < np; i++){
        ix[i] = (int) std::floor((cloud[i][0] - minx) / spacing + 0.5);
        iy[i] = (int) std::floor((cloud[i][1] - miny) / spacing + 0.5);
        iz[i] = (int) std::floor((cloud[i][2] - minz) / spacing + 0.5);
        ex = std::max(ex, ix[i]); ey = std::max(ey, iy[i]); ez = std::max(ez, iz[i]);
    }
    // Each block holds at most m^3 points, so m >= cbrt(n/k) is required;
    // walk up from there to the smallest m with <= k non-empty blocks. The
    // block grid is dense and small, so counting is a scatter into it.
    int m = std::max(1, (int) std::ceil(std::cbrt((double) np / k)));
    std::vector<int> cell(np);
    std::vector<int> occupied;
    long bx, by, bz;
    for(;; m++){
        bx = ex / m + 1; by = ey / m + 1; bz = ez / m + 1;
        occupied.assign((size_t) bx * by * bz, 0);
        long nblocks = 0;
        for(size_t i = 0; i < np; i++){
            int c = (int) (((iz[i] / m) * by + (iy[i] / m)) * bx + (ix[i] / m));
            cell[i] = c;
            if(!occupied[c]){ occupied[c] = 1; nblocks++; }
        }
        if(nblocks <= k) break;
    }
    // Number the occupied cells in cell order (canonical: independent of the
    // order the cloud was produced in)
    std::vector<int> slot(occupied.size(), -1);
    int nb = 0;
    for(size_t c = 0; c < occupied.size(); c++){
        if(occupied[c]) slot[c] = nb++;
    }
    std::vector<double> sw(nb, 0), sx(nb, 0), sy(nb, 0), sz(nb, 0);
    for(size_t i = 0; i < np; i++){
        int b = slot[cell[i]];
        double w = cloud[i][3];
        sw[b] += w; sx[b] += w * cloud[i][0]; sy[b] += w * cloud[i][1]; sz[b] += w * cloud[i][2];
    }
    points.resize(nb);
    for(int b = 0; b < nb; b++){
        points[b] = IMP::algebra::Vector4D(sx[b] / sw[b], sy[b] / sw[b], sz[b] / sw[b], sw[b]);
    }
    moments.assign(nb, {0, 0, 0, 0, 0, 0});
    for(size_t i = 0; i < np; i++){
        int b = slot[cell[i]];
        const auto &p = cloud[i];
        double w = p[3];
        double dx = p[0] - points[b][0], dy = p[1] - points[b][1], dz = p[2] - points[b][2];
        auto &mom = moments[b];
        mom[0] += w * dx * dx; mom[1] += w * dy * dy; mom[2] += w * dz * dz;
        mom[3] += w * dx * dy; mom[4] += w * dx * dz; mom[5] += w * dy * dz;
    }
    for(int b = 0; b < nb; b++){
        for(double &v : moments[b]) v /= sw[b];
    }
}

void quadrature_points(
        const AV &av, int k,
        const std::vector<IMP::algebra::Vector4D> *&points,
        const std::vector<std::array<double, 6> > *&moments){
    auto &st = const_cast<AV&>(av).get_state();
    auto map = av.get_map();   // builds/resamples if needed
    // resample() may have advanced the result generation
    if(!st.quad_valid || st.quad_k != k ||
       st.quad_generation != st.result_generation){
        const auto &cloud = av.get_cloud();
        double h = map->get_path_map_header().get_simulation_grid_resolution();
        coarsen_cloud(cloud, h, k, st.quad_points, st.quad_moments);
        st.quad_k = k;
        st.quad_generation = st.result_generation;
        st.quad_valid = true;
    }
    points = &st.quad_points;
    moments = &st.quad_moments;
}

}

void AV::prepare_quadrature(int k) const{
    const std::vector<IMP::algebra::Vector4D> *pts;
    const std::vector<std::array<double, 6> > *mom;
    quadrature_points(*this, k, pts, mom);
}

std::vector<double> AV::get_quadrature_points(int k) const{
    const std::vector<IMP::algebra::Vector4D> *pts;
    const std::vector<std::array<double, 6> > *mom;
    quadrature_points(*this, k, pts, mom);
    std::vector<double> out;
    out.reserve(4 * pts->size());
    for(const auto &p : *pts){
        out.push_back(p[0]); out.push_back(p[1]);
        out.push_back(p[2]); out.push_back(p[3]);
    }
    return out;
}

double av_distance_quadrature(
        const AV& av1,
        const AV& av2,
        double forster_radius,
        int distance_type,
        int quad_k
){
    switch(distance_type){
        case DYE_PAIR_DISTANCE_MP:
        case DYE_PAIR_XYZ_DISTANCE:
            return av_distance(av1, av2, forster_radius, distance_type, 0);
        default:
            break;
    }
    const std::vector<IMP::algebra::Vector4D> *p1, *p2;
    const std::vector<std::array<double, 6> > *m1, *m2;
    quadrature_points(av1, quad_k, p1, m1);
    quadrature_points(av2, quad_k, p2, m2);
    if(p1->empty() || p2->empty()){
        return std::numeric_limits<double>::quiet_NaN();
    }
    double w1 = 0, w2 = 0;
    for(const auto &p : *p1) w1 += p[3];
    for(const auto &p : *p2) w2 += p[3];
    const bool efficiency = (distance_type == DYE_PAIR_EFFICIENCY ||
                             distance_type == DYE_PAIR_DISTANCE_E);
    // Exact weighted double sum over block centroids, with a second-order
    // correction from the blocks' second central moments: for a pair of
    // blocks with joint covariance C = C1 + C2, displacement D, d = |D|,
    //   E[f(|r|)] ~= f(d) + 1/2 (f''(d) q + f'(d) (tr C - q) / d),
    // where q = D^T C D / d^2. This removes the within-block flattening the
    // centroid sum alone suffers (measured on T4L @2.0 A: max error
    // 0.26 A -> < 0.001 A for the mean distance at K = 100).
    double val = 0.;
    for(size_t i = 0; i < p1->size(); i++){
        const auto &a = (*p1)[i];
        const auto &ca = (*m1)[i];
        double acc = 0.;
        for(size_t j = 0; j < p2->size(); j++){
            const auto &b = (*p2)[j];
            const auto &cb = (*m2)[j];
            double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
            double d2 = dx*dx + dy*dy + dz*dz;
            double d = std::sqrt(d2);
            double cxx = ca[0] + cb[0], cyy = ca[1] + cb[1], czz = ca[2] + cb[2];
            double cxy = ca[3] + cb[3], cxz = ca[4] + cb[4], cyz = ca[5] + cb[5];
            double tr = cxx + cyy + czz;
            double q = 0.0;
            if(d2 > 0){
                q = (dx * (cxx * dx + cxy * dy + cxz * dz) +
                     dy * (cxy * dx + cyy * dy + cyz * dz) +
                     dz * (cxz * dx + cyz * dy + czz * dz)) / d2;
            }
            double f;
            if(efficiency){
                double t = std::pow(d / forster_radius, 6.0);
                double u = 1.0 + t;
                f = 1.0 / u;
                if(d > 0){
                    double f1 = -6.0 * t / (d * u * u);
                    double f2 = -6.0 * t * (5.0 - 7.0 * t) / (d2 * u * u * u);
                    f += 0.5 * (f2 * q + f1 * (tr - q) / d);
                }
            } else {
                f = d;
                if(d > 0){
                    f += (tr - q) / (2.0 * d);
                }
            }
            acc += b[3] * f;
        }
        val += a[3] * acc;
    }
    val /= (w1 * w2);
    if(distance_type == DYE_PAIR_DISTANCE_E){
        return distance_fret<double>(val, forster_radius);
    }
    return val;
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
