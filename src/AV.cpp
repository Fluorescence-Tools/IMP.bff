/**
 * \file AV.cpp
 * \brief Simple Accessible Volume decorator.
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */
#include <IMP/bff/AV.h>
#include <IMP/bff/SelectionExpression.h>
#include <IMP/bff/AVDistance.h>
#include <IMP/bff/StructureIO.h>
#include <IMP/bff/internal/OutputView.h>

#include <chrono>

IMPBFF_BEGIN_NAMESPACE

//! The volume's full point cloud, flat: x, y, z, weight. Defined below.
static std::vector<double> av_cloud(const AV& av);




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
    // One implementation of the asymmetric chi2, `chi2_score`. What stood
    // here selected the *opposite* error bar -- it measured the deviation as
    // data minus model and then read the branch as if it were model minus
    // data, so a model that was too large was judged against error_neg -- and
    // with asymmetric errors this restraint therefore disagreed with
    // `chi2_score` and with the FPS tables by up to two orders of magnitude
    // (experiment 50 A, errors (1, 10): 25.0 against 1.0). Nothing asserted
    // either number; see okf/validation/two_chi2_conventions.md.
    //
    // The factor is 0.5, which is what a Gaussian restraint is worth:
    // -log L = (dev/sigma)^2 / 2, and it is what `IMP::core::Harmonic` scores
    // for k = 1/sigma^2. Halve it here and the restraint is worth half an
    // IMP harmonic, so it mixes with other terms at the wrong weight.
    if(std::isnan(model)){
        return std::numeric_limits<double>::infinity();
    }
    return 0.5 * chi2_score(model, distance, error_neg, error_pos);
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
            case PROBE_PAIR_EFFICIENCY: {
                for (int s = 0; s < n_samples; s++) {
                    auto tmp = sampler1.get_random() - sampler2.get_random();
                    tmp[3] = 0.0;
                    val += fret_efficiency<double>(tmp.get_magnitude(), forster_radius);
                }
                return val / n_samples;
            }
            case PROBE_PAIR_DISTANCE_E: {
                double fret_eff = av_distance(av1, av2, forster_radius, PROBE_PAIR_EFFICIENCY, n_samples);
                return distance_fret<double>(fret_eff, forster_radius);
            }
            case PROBE_PAIR_DISTANCE_MP: {
                IMP::algebra::Vector3D mp1 = av1.get_mean_position();
                IMP::algebra::Vector3D mp2 = av2.get_mean_position();
                return get_l2_norm((mp1 - mp2));
            }
            case PROBE_PAIR_DISTANCE_MIN: {
                // Exact, and not from `n_samples`: a minimum estimated from a
                // sample is biased high, so sampling it would answer a
                // different question than the one the type names.
                return minimum_distance(av_cloud(av1), av_cloud(av2));
            }
            case PROBE_PAIR_XYZ_DISTANCE: {
                IMP::Particle* p1 = av1.get_particle();
                IMP::Particle* p2 = av2.get_particle();
                IMP::algebra::Vector3D mp1 = IMP::core::XYZ(p1).get_coordinates();
                IMP::algebra::Vector3D mp2 = IMP::core::XYZ(p2).get_coordinates();
                return get_l2_norm((mp1 - mp2));
            }
            case PROBE_PAIR_DISTANCE_MEAN:
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

void AV::refresh_cloud_soa() const{
    auto map = get_map();   // builds and resamples on first use
    auto &st = const_cast<AV*>(this)->get_state();
    if(!st.cloud_valid || st.cloud_generation != st.result_generation){
        map->get_xyz_density_soa(st.cloud_x, st.cloud_y, st.cloud_z, st.cloud_w);
        st.cloud_generation = st.result_generation;
        st.cloud_valid = true;
    }
}

const std::vector<IMP::algebra::Vector4D> &AV::get_cloud() const{
    refresh_cloud_soa();
    auto &st = const_cast<AV*>(this)->get_state();
    if(!st.cloud_aos_valid || st.cloud_aos_generation != st.cloud_generation){
        st.cloud.clear();
        st.cloud.reserve(st.cloud_x.size());
        for(size_t i = 0; i < st.cloud_x.size(); i++){
            st.cloud.emplace_back((double) st.cloud_x[i], (double) st.cloud_y[i],
                                  (double) st.cloud_z[i], (double) st.cloud_w[i]);
        }
        st.cloud_aos_generation = st.cloud_generation;
        st.cloud_aos_valid = true;
    }
    return st.cloud;
}

IMP::algebra::Vector3D AV::get_mean_position(bool include_source) const{
    IMP::algebra::Vector3D r = {0.0, 0.0, 0.0};
    /* Zero, not one. `sum` is the total weight the numerator was built from,
       and starting it at 1.0 added a unit of weight belonging to no point:
       every mean position was pulled toward the origin by (1+W)/(2+W), and an
       **empty** volume returned exactly half the source coordinate -- a
       plausible-looking number for a volume that does not exist. */
    double sum = 0.0;
    if(include_source){
        r += get_source_coordinates();
        sum += 1.0;
    }
    refresh_cloud_soa();
    const auto &st = *state_;
    // the same double arithmetic as over Vector4D entries (float values
    // widened on read)
    for(size_t i = 0; i < st.cloud_x.size(); i++){
        const double w = st.cloud_w[i];
        if(w <= 0.0f) continue;
        sum += w;
        r[0] += (double) st.cloud_x[i] * w;
        r[1] += (double) st.cloud_y[i] * w;
        r[2] += (double) st.cloud_z[i] * w;
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

IntKey AV::get_search_grid_factor_key(){
    static const IntKey k("av_search_grid_factor");
    return k;
}

int AV::get_search_grid_factor() const{
    if(get_model()->get_has_attribute(get_search_grid_factor_key(), get_particle_index())){
        return std::max(1, get_model()->get_attribute(get_search_grid_factor_key(), get_particle_index()));
    }
    return 1;
}

void AV::set_search_grid_factor(int f){
    IMP_USAGE_CHECK(f >= 1, "AV: search_grid_factor must be >= 1");
    if(get_model()->get_has_attribute(get_search_grid_factor_key(), get_particle_index())){
        get_model()->set_attribute(get_search_grid_factor_key(), get_particle_index(), f);
    } else {
        get_model()->add_attribute(get_search_grid_factor_key(), get_particle_index(), f);
    }
    if(state_){
        state_->have_result = false;   // the factor is part of the result
        state_->coarse_map = nullptr;
    }
}

IntKey AV::get_search_stencil_key(){
    static const IntKey k("av_search_stencil");
    return k;
}

int AV::get_search_stencil() const{
    if(get_model()->get_has_attribute(get_search_stencil_key(), get_particle_index())){
        return get_model()->get_attribute(get_search_stencil_key(), get_particle_index());
    }
    // **74 by default: the reference metric.** LabelLib's `essentialNeighbours()`
    // (`FlexLabel/src/FlexLabel.cxx:149`) keeps the shells with squared offset
    // length {1,2,3,5,6}; bff reaches the same path lengths with a sqrt(6)
    // neighbour radius. Measured with no obstacles at all, where the accessible
    // volume must be a sphere: 74 gives 30682 voxels for a 20 A linker against
    // LabelLib's 30688 and an analytic 33510 -- 0.02 %. The 26 stencil
    // ({1,2,3}) gives 26146, and the 15 % it loses is not at the rim but in
    // every shell, because its isopath surface is cubic rather than spherical.
    //
    // **26 is the speed option**, and it is a real one: on T4L site 22 at 1.0 A
    // a re-resample is 1.8 ms against 3.4 ms, so ~1.9x, for ~20 % less volume.
    // Choose it where throughput matters more than the metric.
    //
    // The cost of the default: 74's longest jump is sqrt(6) ~ 2.45 voxels, so it
    // **can tunnel through walls thinner than that**, which 26 (sqrt(3) ~ 1.73)
    // cannot. LabelLib has the same property -- this is what the reference
    // metric is -- and it is why a grid coarse enough to make an obstacle layer
    // thinner than ~2.5 voxels is not safe with either. `test_empty_av_gives_nan`
    // records a concrete instance at 1.5 A spacing with `linker_width` 0.5.
    return 74;
}

void AV::set_search_stencil(int stencil){
    IMP_USAGE_CHECK(stencil == 26 || stencil == 30 || stencil == 74,
                    "AV: search stencil must be 74 (default, the LabelLib "
                    "reference metric), 26 (the speed option, ~1.9x faster for "
                    "~20 % less volume) or 30 (historical)");
    if(get_model()->get_has_attribute(get_search_stencil_key(), get_particle_index())){
        get_model()->set_attribute(get_search_stencil_key(), get_particle_index(), stencil);
    } else {
        get_model()->add_attribute(get_search_stencil_key(), get_particle_index(), stencil);
    }
    if(av_map_){
        av_map_ = nullptr;      // header radius and offsets change: rebuild
        if(state_){ state_->have_result = false; state_->coarse_map = nullptr; }
    }
}

IntKey AV::get_compensate_stencil_key(){
    static const IntKey k("av_compensate_stencil");
    return k;
}

bool AV::get_compensate_stencil() const{
    if(get_model()->get_has_attribute(get_compensate_stencil_key(), get_particle_index())){
        return get_model()->get_attribute(get_compensate_stencil_key(), get_particle_index()) != 0;
    }
    // **Opt-in, not default.** Asking for stencil 26 or 30 should give that
    // stencil's own answer -- a caller pinning the historical 30 metric gets it
    // unchanged. Compensation is for the case where 26 is chosen *for speed* and
    // the reference volume is still wanted; it is then an explicit request.
    return false;
}

void AV::set_compensate_stencil(bool tf){
    int v = tf ? 1 : 0;
    if(get_model()->get_has_attribute(get_compensate_stencil_key(), get_particle_index())){
        get_model()->set_attribute(get_compensate_stencil_key(), get_particle_index(), v);
    } else {
        get_model()->add_attribute(get_compensate_stencil_key(), get_particle_index(), v);
    }
    if(av_map_){
        av_map_ = nullptr;              // the path threshold changes: rebuild
        if(state_){ state_->have_result = false; state_->coarse_map = nullptr; }
    }
}

IntKey AV::get_radii_source_key(){
    static const IntKey k("av_radii_source");
    return k;
}

std::string AV::get_radii_source() const{
    // Absent attribute = the default, which is the radii the particles carry
    // (IMP's own) -- so that the volume and `clash_container` size an atom the
    // same way. Stored as an int for the same reason `av_search_mode` is: a
    // Model IntKey is cheap and the string is the API.
    if(get_model()->get_has_attribute(get_radii_source_key(), get_particle_index())){
        return av_radii_source_to_string((AVRadiiSource)
            get_model()->get_attribute(get_radii_source_key(), get_particle_index()));
    }
    return av_radii_source_to_string(AV_RADII_IMP);
}

void AV::set_radii_source(std::string source){
    const int v = (int) av_radii_source_from_string(source);  // validates
    if(get_model()->get_has_attribute(get_radii_source_key(), get_particle_index())){
        if(get_model()->get_attribute(get_radii_source_key(), get_particle_index()) == v){
            return;
        }
        get_model()->set_attribute(get_radii_source_key(), get_particle_index(), v);
    } else {
        if(v == (int) AV_RADII_IMP) return;  // already the default
        get_model()->add_attribute(get_radii_source_key(), get_particle_index(), v);
    }
    // The obstacle set is built from it, exactly as it is from the strip mask,
    // so a new source invalidates the map.
    av_map_ = nullptr;
}

IntKey AV::get_search_mode_key(){
    static const IntKey k("av_search_mode");
    return k;
}

std::string AV::get_search_mode() const{
    if(get_model()->get_has_attribute(get_search_mode_key(), get_particle_index())){
        return get_model()->get_attribute(get_search_mode_key(), get_particle_index()) == 1
               ? "euclidean" : "dijkstra";
    }
    return "dijkstra";
}

void AV::set_search_mode(std::string mode){
    int v;
    if(mode == "dijkstra") v = 0;
    else if(mode == "euclidean") v = 1;
    else IMP_THROW("AV: search mode must be \"dijkstra\" or \"euclidean\"", IMP::ValueException);
    if(get_model()->get_has_attribute(get_search_mode_key(), get_particle_index())){
        get_model()->set_attribute(get_search_mode_key(), get_particle_index(), v);
    } else {
        get_model()->add_attribute(get_search_mode_key(), get_particle_index(), v);
    }
    if(av_map_) av_map_->set_euclidean_search(v == 1);
    if(state_){
        state_->have_result = false;
        if(state_->coarse_map) state_->coarse_map->set_euclidean_search(v == 1);
    }
}

void AV::set_occupancy_registry(AVOccupancyRegistry *registry){
    if(registry && !get_space_fixed()){
        IMP_THROW("AV: a shared occupancy registry requires space_fixed",
                  IMP::ValueException);
    }
    auto &st = get_state();
    // A shared occupancy grid *is* one obstacle set, rasterised once for every
    // volume in it. A `strip_mask` gives this volume an obstacle set of its
    // own -- its labelling site's side chain is absent for it and present for
    // its neighbours -- so it cannot be in that grid, and takes a private map
    // instead. Sharing stays for the volumes that strip nothing.
    if(registry && !get_strip_mask().empty()) registry = nullptr;
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
    lattice_window(get_source_coordinates(), get_effective_linker_length(), h, k0, n);
    auto &st = get_state();
    st.registry->get_map(h, get_linker_width() * 0.5)
        ->request_window(k0[0], k0[1], k0[2], n, n, n);
    st.registry->get_map(h, get_radius1())
        ->request_window(k0[0], k0[1], k0[2], n, n, n);
}

namespace {

//! The radius each atom obstructs this volume with; the mask's atoms get zero.
/*! A stripped atom keeps its place in the list and loses its **size**. That is
    all "removing an obstacle" means to a rasteriser that tests
    `distance < radius`, and it is cheaper and safer than handing each volume a
    list of its own: every volume indexes the same particles whatever it
    strips, so the incremental window machinery needs no per-volume mapping,
    and nothing is mutated on a model that another thread is reading.

    The source atom keeps its radius whatever the mask says -- a mask must not
    be able to delete the anchor. It is dropped from the obstacle set all the
    same, as FPS drops it (`av_routines.cpp:50`), but by
    drop_source_obstruction() rather than here: that is one mechanism on every
    path, where doing it here for some paths and there for others let a shared
    raster and a private one disagree about the same volume. An empty mask
    returns an empty override, which means "use the particles' own radii". A
    mask that cannot be read raises rather than being ignored: computing
    against obstacles the position said to remove is the failure this
    guards.

    The **radii set** enters here too, and for the same reason: the vector this
    returns is what the raster inflates by, so "which van der Waals radii" and
    "which atoms are transparent" are one answer, computed once. Under
    `AV_RADII_IMP` (the default) each atom keeps whatever the particle carries
    -- IMP's united-atom set after `read_pdb`, and the same radius
    `clash_container` measures overlap with; under `AV_RADII_OLGA` it takes
    Olga's name-keyed radius (#IMP::bff::olga_vdw_radius). An empty override --
    "use the particles' own radii" -- is therefore returned only for
    `AV_RADII_IMP` with no mask, and that is the point of the default: no
    override, so no second opinion anywhere about how big an atom is. */
std::vector<double> obstacle_radii(IMP::atom::Hierarchy root,
                                   IMP::Particle* source,
                                   const IMP::ParticlesTemp& all,
                                   const std::string& mask,
                                   AVRadiiSource radii_source) {
    if (mask.empty() && radii_source == AV_RADII_IMP) {
        return std::vector<double>();
    }

    std::set<IMP::ParticleIndex> transparent;
    if (!mask.empty()) {
        const IMP::atom::Selection selection =
                selection_from_expression(root, mask);
        const IMP::ParticleIndexes stripped =
                selection.get_selected_particle_indexes();
        transparent.insert(stripped.begin(), stripped.end());
        transparent.erase(source->get_index());
    }

    std::vector<double> radii(all.size(), 0.0);
    for (std::size_t i = 0; i < all.size(); ++i) {
        if (transparent.count(all[i]->get_index())) continue;  // no size
        /* Olga's table is keyed by *atom name*, so it has an opinion only
           about particles that are atoms. A particle that is not one -- a
           coarse-grained bead, a synthetic obstacle built by a test, anything
           a caller put in the hierarchy -- is outside the table's domain, and
           giving it the 1.50 A unknown-name fallback would silently shrink an
           obstacle set the caller sized deliberately. It keeps its own radius.
           An atom whose *name* the table misses does take the fallback: that
           case is inside Olga's domain and is exactly what Olga does. */
        if (radii_source == AV_RADII_OLGA &&
            IMP::atom::Atom::get_is_setup(all[i])) {
            radii[i] = olga_vdw_particle_radius(all[i]);
        } else {
            radii[i] = IMP::core::XYZR::get_is_setup(all[i])
                               ? IMP::core::XYZR(all[i]).get_radius()
                               : 0.0;
        }
    }
    return radii;
}

}  // namespace

namespace {

//! Subtract the attachment atom's own obstruction from an occupancy window.
/*!
    FPS drops the attachment atom before it rasterises anything
    (`av_routines.cpp:50`: `if (i == atom_i) continue;`), and it has to be
    dropped: an atom cannot block the linker that is tied to it. `IMP.bff` kept
    it, so the source sat inside its own inflated sphere and the volume came
    back empty unless the clearance was raised past it -- a knob standing in
    for a missing rule. Measured on FPS's own reference cloud, dropping it
    takes `p_1bp` from 0.88 to **0.99** of FPS's volume (PRD-121).

    It is done here rather than through the radii override because the
    occupancy raster is **shared** between every volume of the same (spacing,
    extra-radius) class (PRD-105), and each volume drops a *different* atom.
    Giving each its own raster would cost that sharing. The raster stores a
    per-voxel atom **count**, so one atom's contribution subtracts exactly: a
    voxel only this atom covered falls to zero and opens, a voxel any other
    atom covers stays blocked.

    \param[in,out] counts the window, as read from the occupancy
    \param[in] k0,n the window's lattice origin and edge
    \param[in] spacing the lattice spacing
    \param[in] c the attachment atom's position
    \param[in] radius its own radius plus this occupancy's extra radius --
               the same `r + extra_radius` the raster used, so the voxel set
               subtracted is exactly the one that was added
*/
template <typename T>
void drop_source_obstruction(T *counts, const IMP::algebra::Vector3D &origin,
                             int n, double spacing,
                             const IMP::algebra::Vector3D &c,
                             double source_radius, double extra) {
    /* Nothing to subtract when the attachment atom contributed nothing. A
       radius of zero means it is not an obstacle -- the array door builds its
       obstacles from a caller's list and the source is not in it -- and
       subtracting a sphere that was never added would open voxels no rule
       opened. */
    if (source_radius <= 0.0) return;
    const double radius = source_radius + extra;
    if (radius <= 0.0) return;
    const double r2 = radius * radius;
    const long nxy = (long) n * n;
    for (int iz = 0; iz < n; iz++) {
        const double dz = origin[2] + iz * spacing - c[2];
        const double dz2 = dz * dz;
        if (dz2 >= r2) continue;
        for (int iy = 0; iy < n; iy++) {
            const double dy = origin[1] + iy * spacing - c[1];
            const double dyz2 = dz2 + dy * dy;
            if (dyz2 >= r2) continue;
            const long row = (long) iz * nxy + (long) iy * n;
            for (int ix = 0; ix < n; ix++) {
                const double dx = origin[0] + ix * spacing - c[0];
                // strict `<`, as AVOccupancyMap::add_sphere tests
                if (dyz2 + dx * dx < r2) {
                    T &v = counts[row + ix];
                    if (v > 0) v -= 1;
                }
            }
        }
    }
}

}  // namespace

void AV::init_path_map(){
    auto path_map_header = create_path_map_header();
    av_map_ = new IMP::bff::PathMap(path_map_header);
    if(get_space_fixed() && (get_search_stencil() == 26
                            || get_search_stencil() == 74)){
        av_map_->set_symmetric_stencil(true);
    }
    if(get_space_fixed()) av_map_->set_euclidean_search(get_search_mode() == "euclidean");
    IMP::Particle* parent = get_model()->get_particle(get_particle_index(0));

    auto h = IMP::atom::Hierarchy(get_model(), parent->get_index());
    auto root = IMP::atom::get_root(h);
    const IMP::ParticlesTemp leaves = get_leaves(root);
    av_map_->set_particles(leaves);
    /* The legacy anchoring keeps the legacy radii, and refuses Olga's the way
       it refuses the accessible contact volume (see resample_legacy step 5b).
       `space_fixed=False` exists for exactly one thing -- reproducing the
       pre-PRD-105 numbers byte for byte, which
       `references/prd105_legacy_pins.json` is there to prove -- and those
       numbers were computed against the radii the particles carry. Said out
       loud once per handle rather than left to be discovered. Since the
       default is those radii again, this now fires only for a volume that
       explicitly asked for Olga's. */
    AVRadiiSource radii_source =
            av_radii_source_from_string(get_radii_source());
    if(!get_space_fixed() && radii_source != AV_RADII_IMP){
        if(!get_state().warned_radii_source){
            get_state().warned_radii_source = true;
            IMP_WARN("AV " << get_particle()->get_name() << ": radii_source=\""
                     << av_radii_source_to_string(radii_source)
                     << "\" is ignored under space_fixed=False. The legacy"
                        " anchoring is frozen on the pre-PRD-105 numbers, which"
                        " were computed with the radii the particles carry; use"
                        " the default anchoring to get Olga's table."
                     << std::endl);
        }
        radii_source = AV_RADII_IMP;
    }
    const std::vector<double> orad =
            obstacle_radii(root, parent, leaves, get_strip_mask(),
                           radii_source);
    av_map_->set_obstacle_radii(orad);
    if(get_chain_weighting()){
        av_map_->set_linker_weighting(linker_weighting(get_linker_length()));
    }
    // A fresh map: whatever the state remembers about tiles is stale.
    auto &st = get_state();
    st.have_result = false;
    st.have_window = false;
    st.private1 = nullptr;
    st.private2 = nullptr;
    /* The attachment atom's radius *as this map's raster uses it*, found once
       here rather than scanned for on every frame. Only an override answers
       it; with none, the raster reads the particle and so does prepare(). */
    st.source_obstacle_radius = -1.0;
    if(!orad.empty()){
        for(std::size_t i = 0; i < leaves.size(); ++i){
            if(leaves[i] == parent){ st.source_obstacle_radius = orad[i]; break; }
        }
    }
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
    /* An empty volume is a result, not an error -- a site can genuinely be
       buried -- but it must not pass for a computed one. It used to: nothing
       was logged, and the mean position came back as a plausible coordinate
       near the attachment atom. The clearance is named because it is the cause
       in every case seen so far, and because it is the one parameter a caller
       would not think to look at.
       Only the single-call path warns; the network's split
       prepare/compute/finish path evaluates thousands of volumes per frame and
       would drown a run in messages. */
    if(get_map() != nullptr && get_map()->get_xyz_density().empty()){
        IMP_WARN("AV " << get_particle()->get_name() << ": no accessible voxel."
                 << " linker_length=" << get_linker_length()
                 << " linker_width=" << get_linker_width()
                 << " radius1=" << get_radius1()
                 << " grid=" << get_simulation_grid_resolution()
                 << " clearance=" << get_effective_allowed_sphere_radius()
                 << ". The search inflates obstacles by half the linker width,"
                    " so a clearance below that walls the source in."
                 << std::endl);
    }
}

// The pre-PRD-105 evaluation, kept byte-for-byte for `space_fixed=False`.
/* The pre-PRD-105 evaluation, and deliberately *not* given the attachment-atom
   drop the lattice path has. `space_fixed=False` exists to reproduce the old
   numbers -- `references/prd105_legacy_pins.json` is there to prove it still
   does -- so changing its physics would remove the only thing it is for. The
   consequence is real and worth knowing: a volume computed with
   `space_fixed=False` grows out of a source that obstructs itself, and is
   smaller than the same volume on the default path. */
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
    critical_radius = get_effective_linker_length();
    map->fill_sphere(source, critical_radius, TILE_PENALTY_THRESHOLD, true);

    // 3.1 Unblock voxels in initial sphere
    critical_radius = get_effective_allowed_sphere_radius();
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

    /* 5b. The accessible contact volume is **not** applied here, and this is
       the one place in the module where an ACV request is honoured by being
       refused. `space_fixed=False` exists for exactly one thing -- reproducing
       the pre-PRD-105 numbers byte for byte, which
       `references/prd105_legacy_pins.json` is there to prove -- and those
       numbers were computed with the weighting inert. Adding it here would
       change them and remove the only reason the path exists -- and it would
       need a second copy of the rule, against the tiles rather than the SoA
       raster, which is one more place for the two to drift apart. Said out
       loud once per handle rather than left to be discovered. */
    if(get_contact_volume_thickness() > 0.0 &&
       get_contact_volume_trapped_fraction() >= 0.0 &&
       !get_state().warned_contact_volume){
        get_state().warned_contact_volume = true;
        IMP_WARN("AV " << get_particle()->get_name()
                 << ": contact_volume_thickness="
                 << get_contact_volume_thickness()
                 << " is ignored under space_fixed=False. The legacy anchoring"
                    " is frozen on the pre-PRD-105 numbers, which were computed"
                    " without it; use the default anchoring to get the"
                    " accessible contact volume." << std::endl);
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
    if(state_ && state_->pending) resample_lattice_compute(false);
}

void AV::resample_compute_search(){
    if(state_ && state_->pending) resample_lattice_compute(true);
}

void AV::resample_compute_carve(){
    if(state_ && state_->pending) resample_lattice_compute_carve();
}

void AV::resample_finish(){
    if(state_ && state_->pending) resample_lattice_finish();
}

void AV::resample_lattice(bool shift_xyz, bool force_full){
    resample_lattice_prepare(shift_xyz, force_full);
    resample_lattice_compute(false);
    resample_lattice_finish();
}

void AV::resample_lattice_prepare(bool shift_xyz, bool force_full){
    auto map = get_map();
    auto &st = get_state();
    st.pending = false;

    const double h = get_simulation_grid_resolution();
    const double ll = get_effective_linker_length();
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
    // The voxel locations (set_origin) are recomputed in the compute phase,
    // per map and therefore on the threads.
    bool set_origin_needed = false;
    if(dims_changed || params_changed){
        auto path_map_header = create_path_map_header();
        path_map_header.set_path_origin(source, grid_origin);
        map->set_path_map_header(path_map_header);
        set_origin_needed = true;
    } else if(rolled){
        header.set_path_origin(source, grid_origin);
        set_origin_needed = true;
    } else {
        // same window; the label may have moved inside its voxel
        header.set_path_origin(source, grid_origin);
    }
    for(int d = 0; d < 3; d++) st.k0[d] = k0[d];
    st.n = n;
    st.have_window = true;

    // 1. Occupancy sources for the two passes
    //
    // The second pass carves by the dye radius. AV3 carries three of them, so
    // it needs one occupancy source per radius: the density is the fraction of
    // probes that fit, which is what LabelLib/Olga's dyeDensityAV3 computes.
    // AV1 has a single active radius and takes exactly the path it always did.
    const double extra1 = get_linker_width() * 0.5;
    const std::vector<double> dye_radii = get_active_radii();
    const double extra2 = dye_radii[0];
    st.pending_extra1 = extra1;
    st.pending_extra2 = extra2;
    st.pending_extra_dye = dye_radii;
    AVOccupancyMap *occ1; AVOccupancyMap *occ2;
    st.pending_occ_dye.clear();
    if(st.registry){
        /* The radii set is a property of the shared raster, not of this
           volume: hand it over before the first map is created, so the maps
           are built with it rather than re-rastered afterwards. O(1) after the
           first volume -- the registry compares the *source*, not the vector.
           A strip mask never reaches here (set_occupancy_registry drops a
           masked volume to a private raster), so this vector is the pure
           radii-set answer. */
        st.registry->adopt_obstacle_radii(
                av_radii_source_from_string(get_radii_source()),
                map->get_obstacle_radii());
        // (a registry requires space_fixed, so no legacy override here)
        occ1 = st.registry->get_map(h, extra1);
        occ2 = st.registry->get_map(h, extra2);
        occ1->request_window(k0[0], k0[1], k0[2], n, n, n);
        occ2->request_window(k0[0], k0[1], k0[2], n, n, n);
        st.pending_occ_dye.push_back(occ2);
        for(size_t i = 1; i < dye_radii.size(); i++){
            AVOccupancyMap *o = st.registry->get_map(h, dye_radii[i]);
            o->request_window(k0[0], k0[1], k0[2], n, n, n);
            st.pending_occ_dye.push_back(o);
        }
    } else {
        IMP::ParticlesTemp ps(map->ps_.begin(), map->ps_.end());
        // The volume's own radii: the mask's atoms carry zero, so they take
        // part in the raster's indexing and block nothing. A private map is
        // the only place this can be honoured -- a shared raster is one
        // obstacle set for every volume in it.
        const std::vector<double> &radii = map->get_obstacle_radii();
        if(!st.private1 || st.private1->get_extra_radius() != extra1
                        || st.private1->get_spacing() != h){
            st.private1 = new AVOccupancyMap(h, extra1, ps);
            st.private1->set_obstacle_radii(radii);
            st.private1->set_was_used(true);
        }
        if(!st.private2 || st.private2->get_extra_radius() != extra2
                        || st.private2->get_spacing() != h){
            st.private2 = new AVOccupancyMap(h, extra2, ps);
            st.private2->set_obstacle_radii(radii);
            st.private2->set_was_used(true);
        }
        occ1 = st.private1.get();
        occ2 = st.private2.get();
        occ1->set_window(k0[0], k0[1], k0[2], n, n, n);
        occ2->set_window(k0[0], k0[1], k0[2], n, n, n);
        st.pending_occ_dye.push_back(occ2);
        // AV3 without a shared registry: one private source per extra radius
        st.private_dye.resize(dye_radii.size() > 1 ? dye_radii.size() - 1 : 0);
        for(size_t i = 1; i < dye_radii.size(); i++){
            IMP::Pointer<AVOccupancyMap> &pm = st.private_dye[i - 1];
            if(!pm || pm->get_extra_radius() != dye_radii[i]
                   || pm->get_spacing() != h){
                pm = new AVOccupancyMap(h, dye_radii[i], ps);
                pm->set_obstacle_radii(radii);
                pm->set_was_used(true);
            }
            pm->set_window(k0[0], k0[1], k0[2], n, n, n);
            st.pending_occ_dye.push_back(pm.get());
        }
    }
    if(!(st.registry && st.registry_driven_externally)){
        occ1->update(force_full);
        for(size_t i = 0; i < st.pending_occ_dye.size(); i++){
            st.pending_occ_dye[i]->update(force_full);
        }
    }

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
        st.generation1 = occ1->get_generation_after_pending();
        st.generation2 = occ2->get_generation_after_pending();
        if(shift_xyz) set_coordinates(st.last_mean);
        return;
    }
    if(rolled || dims_changed || params_changed || !st.have_result){
        st.n_full++;
    } else {
        st.n_local++;
    }

    // Coarse search grid: lattice points with index = multiple of f, covering
    // the fine window. Built here (Model access: none, but the header/resize
    // are cheap and serial keeps it simple); the coarse occupancy is read in
    // compute().
    const int factor = get_search_grid_factor();
    st.pending_factor = factor;
    if(factor > 1){
        auto floor_div = [](int a, int b){ return (a >= 0) ? a / b : -((-a + b - 1) / b); };
        int c0[3], cn[3];
        for(int d = 0; d < 3; d++){
            c0[d] = floor_div(k0[d], factor);
            int c1 = floor_div(k0[d] + n - 1, factor);
            cn[d] = c1 - c0[d] + 1;
        }
        bool rebuild = !st.coarse_map ||
            cn[0] != st.coarse_n[0] || cn[1] != st.coarse_n[1] || cn[2] != st.coarse_n[2];
        bool moved = !rebuild &&
            (c0[0] != st.coarse_k0[0] || c0[1] != st.coarse_k0[1] || c0[2] != st.coarse_k0[2]);
        if(rebuild || params_changed){
            IMP::bff::PathMapHeader ch(ll, h * factor);
            ch.update_map_dimensions(cn[0], cn[1], cn[2]);
            if(get_search_stencil() == 26) ch.set_neighbor_radius(std::sqrt(3.0) + 1e-6);
            ch.set_path_origin(source, IMP::algebra::Vector3D(
                c0[0] * factor * h, c0[1] * factor * h, c0[2] * factor * h));
            if(!st.coarse_map){
                st.coarse_map = new IMP::bff::PathMap(ch, "AVCoarsePathMap%1%");
                st.coarse_map->set_was_used(true);
                if(get_search_stencil() == 26) st.coarse_map->set_symmetric_stencil(true);
            } else {
                st.coarse_map->set_path_map_header(ch);
            }
            moved = true;
        }
        for(int d = 0; d < 3; d++){ st.coarse_k0[d] = c0[d]; st.coarse_n[d] = cn[d]; }
        if(moved){
            st.coarse_map->set_origin(IMP::algebra::Vector3D(
                c0[0] * factor * h, c0[1] * factor * h, c0[2] * factor * h));
        }
    }

    // Everything compute() needs, gathered while we may still touch the Model
    st.pending = true;
    st.pending_stage = 0;
    st.pending_shift_xyz = shift_xyz;
    st.pending_ll = ll;
    st.pending_allowed = get_effective_allowed_sphere_radius();
    /* Unconditional, shared raster or private: the same subtraction on both,
       so the two cannot compute different volumes for the same position.
       It must be the radius the *raster* used for this atom, not the one the
       Model carries -- under `radii_source = "olga"` those differ, and
       subtracting a sphere of the wrong size would leave a ring of blocked
       voxels around the anchor (too small) or open voxels no rule opened
       (too large). obstacle_radii() is the one place that decides; this reads
       the answer back out of it. */
    st.pending_source_radius =
            st.source_obstacle_radius >= 0.0
                    ? st.source_obstacle_radius
                    : (IMP::core::XYZR::get_is_setup(get_model(),
                                                     get_particle_index(0))
                               ? IMP::core::XYZR(get_model(),
                                                 get_particle_index(0))
                                         .get_radius()
                               : 0.0);
    st.pending_source = source;
    st.pending_set_origin = set_origin_needed;
    st.pending_grid_origin = grid_origin;
    st.pending_occ1 = occ1;
    st.pending_occ2 = occ2;
    if(st.pending_occ_dye.empty()) st.pending_occ_dye.push_back(occ2);
    st.pending_gen1 = occ1->get_generation_after_pending();
    st.pending_gen2 = occ2->get_generation_after_pending();
    st.last_parameter = parameter;
}

// Touches only this AV's map and lattice state: safe to run concurrently
// with other AVs' compute phases (their occupancy sources are read-only now).
void AV::resample_lattice_compute(bool split_stages){
    auto &st = *state_;
    if(!st.pending || st.pending_stage != 0) return;
    st.compute_t0 = std::chrono::steady_clock::now();
    PathMap *map = av_map_.get();
    const int *k0 = st.k0;
    const int n = st.n;
    const IMP::algebra::Vector3D &source = st.pending_source;

    if(st.pending_set_origin){
        map->set_origin_fast(st.pending_grid_origin);
        st.pending_set_origin = false;
    }

    // 3. Obstacles inflated by half the linker width, from the lattice
    //    (integer counts straight into the penalty pass)
    long nvox = map->get_number_of_voxels();
    double *data = map->get_data();
    if(st.window_counts.size() != (size_t) nvox) st.window_counts.resize(nvox);
    int32_t *counts = st.window_counts.data();
    st.pending_occ1->read_window_counts(k0[0], k0[1], k0[2], n, n, n, counts);
    // The attachment atom does not obstruct the linker tied to it.
    {
        const double h_ = map->get_path_map_header()
                                  .get_simulation_grid_resolution();
        drop_source_obstruction(
                counts,
                IMP::algebra::Vector3D(k0[0] * h_, k0[1] * h_, k0[2] * h_), n,
                h_, source, st.pending_source_radius, st.pending_extra1);
    }

    // 4./5./6. Block voxels further away from the source than the linker
    //    length, open the allowed sphere, and search -- in one pass over
    //    the window (search_lattice with spheres). Only tiles with
    //    cost * spacing < linker length carry density, so the exact search
    //    stops there; the source tile keeps the default cost as the
    //    historical search left it (its voxel is not part of the cloud).
    long source_idx = map->get_voxel_by_location(source);
    const float h = (float) map->get_path_map_header().get_simulation_grid_resolution();
    const float ll = (float) map->get_path_map_header().get_max_path_length();
    // smallest float c with c * h >= ll (the density test is c * h < ll)
    float bound = ll / h;
    while(bound * h < ll) bound = std::nextafter(bound, std::numeric_limits<float>::infinity());
    if(st.pending_factor <= 1){
        map->search_lattice(source_idx, bound, source, st.pending_ll, st.pending_allowed, counts);
    } else {
        // the coarse branch reads the fine occupancy from the map data
        for(long i = 0; i < nvox; i++) data[i] = (double) counts[i];
        // Coarse search on the points with fine index = multiple of f, then
        // every fine tile takes min over its (up to 8) surrounding coarse
        // points of (coarse cost * f + straight hop), in fine voxel units.
        const int f = st.pending_factor;
        PathMap *cm = st.coarse_map.get();
        const int *c0 = st.coarse_k0;
        const int *cn = st.coarse_n;
        const long cnvox = (long) cn[0] * cn[1] * cn[2];
        double *cdata = cm->get_data();
        st.pending_occ1->read_window_strided(c0[0] * f, c0[1] * f, c0[2] * f,
                                             cn[0], cn[1], cn[2], f, cdata);
        cm->fill_sphere(source, st.pending_ll, TILE_PENALTY_THRESHOLD, true);
        cm->fill_sphere(source, st.pending_allowed, 0, false);
        long csrc = cm->get_voxel_by_location(source);
        cdata[csrc] = 0.0;   // the source's coarse tile is always open
        cm->search_lattice(csrc, bound / f);
        cm->cost[csrc] = 0.0f;   // its own cost is real here (extension below)
        // extend to the fine grid
        const long nxy = (long) n * n;
        const long cnx = cn[0], cnxy = (long) cn[0] * cn[1];
        map->cost.assign(nvox, TILE_COST_DEFAULT);
        map->visited.assign(nvox, false);
        map->penalty_soa_.assign(nvox, 0.0f);
        for(int iz = 0; iz < n; iz++){
            const int gz = k0[2] + iz;                 // fine global index
            const int czf = (gz >= 0) ? gz / f : -((-gz + f - 1) / f);   // floor(gz/f)
            const int cz_lo = czf - c0[2];
            for(int iy = 0; iy < n; iy++){
                const int gy = k0[1] + iy;
                const int cyf = (gy >= 0) ? gy / f : -((-gy + f - 1) / f);
                const int cy_lo = cyf - c0[1];
                for(int ix = 0; ix < n; ix++){
                    const int gx = k0[0] + ix;
                    const int cxf = (gx >= 0) ? gx / f : -((-gx + f - 1) / f);
                    const int cx_lo = cxf - c0[0];
                    float best = TILE_COST_DEFAULT;
                    for(int dz = 0; dz < 2; dz++){
                        int cz = cz_lo + dz; if(cz < 0 || cz >= cn[2]) continue;
                        float fz = (float) (gz - (c0[2] + cz) * f);
                        for(int dy = 0; dy < 2; dy++){
                            int cy = cy_lo + dy; if(cy < 0 || cy >= cn[1]) continue;
                            float fy = (float) (gy - (c0[1] + cy) * f);
                            for(int dx = 0; dx < 2; dx++){
                                int cx = cx_lo + dx; if(cx < 0 || cx >= cn[0]) continue;
                                long ci = cz * cnxy + cy * cnx + cx;
                                float cc = cm->cost[ci];
                                if(!(cc < TILE_COST_DEFAULT)) continue;
                                float fx = (float) (gx - (c0[0] + cx) * f);
                                float hop = std::sqrt(fx*fx + fy*fy + fz*fz);
                                float c = cc * (float) f + hop;
                                if(c < best) best = c;
                            }
                        }
                    }
                    if(best < TILE_COST_DEFAULT){
                        long i = iz * nxy + iy * n + ix;
                        map->cost[i] = best;
                        map->visited[i] = true;
                    }
                }
            }
        }
        // the fine source tile keeps the default cost, as in the fine path
        map->cost[source_idx] = TILE_COST_DEFAULT;
        map->reached_valid_ = true;
        map->soa_valid_ = true;
        (void) cnvox;
    }

    st.pending_stage = 1;
    if(!split_stages) resample_lattice_compute_carve();
}

// Second half of the compute phase: the dye-radius carve, the cloud and the
// mean. Separate so a scheduler can interleave it with other AVs' searches.
void AV::resample_lattice_compute_carve(){
    auto &st = *state_;
    if(!st.pending || st.pending_stage != 1) return;
    PathMap *map = av_map_.get();
    const int *k0 = st.k0;
    const int n = st.n;
    const IMP::algebra::Vector3D &source = st.pending_source;
    double *data = map->get_data();

    // 7. Remove tiles closer to obstacles than the dye radius
    long nvox = map->get_number_of_voxels();
    if(st.window_counts.size() != (size_t) nvox) st.window_counts.resize(nvox);
    int32_t *counts = st.window_counts.data();
    st.pending_occ2->read_window_counts(k0[0], k0[1], k0[2], n, n, n, counts);
    const double spacing =
            map->get_path_map_header().get_simulation_grid_resolution();
    // ...nor does it obstruct the dye, for the same reason.
    const IMP::algebra::Vector3D win_origin(k0[0] * spacing, k0[1] * spacing,
                                            k0[2] * spacing);
    drop_source_obstruction(counts, win_origin, n, spacing, source,
                            st.pending_source_radius, st.pending_extra2);
    const size_t n_dye = st.pending_occ_dye.empty() ? 1 : st.pending_occ_dye.size();
    if(n_dye <= 1){
        map->carve_lattice(counts);
    } else {
        // AV3: the density is the fraction of dye radii that fit, which is
        // LabelLib/Olga's dyeDensityAV3. `counts` is already the first radius.
        if(st.window_counts_dye.size() != (n_dye - 1) * (size_t) nvox){
            st.window_counts_dye.resize((n_dye - 1) * (size_t) nvox);
        }
        std::vector<const int32_t *> src(n_dye);
        src[0] = counts;
        for(size_t i = 1; i < n_dye; i++){
            int32_t *c = st.window_counts_dye.data() + (i - 1) * (size_t) nvox;
            st.pending_occ_dye[i]->read_window_counts(k0[0], k0[1], k0[2],
                                                     n, n, n, c);
            drop_source_obstruction(c, win_origin, n, spacing, source,
                                    st.pending_source_radius,
                                    st.pending_extra_dye[i]);
            src[i] = c;
        }
        map->carve_lattice_fractional(src.data(), (int) n_dye);
    }
    /* The accessible *contact* volume (PRD-121 G9). Applied here, between the
       carve and the cloud, because the trapped share is a share of the points
       that reach the cloud: the search has run, the densities are final, and
       nothing has read them yet. `counts` is the dye-radius occupancy the
       carve just used -- Olga's `occupancyVdWDye` -- so the surface the layer
       is measured from is the one the dye actually cannot enter.
       For AV3 the classification uses the first radius's raster, the one
       carve_lattice_fractional() already treats as canonical; Olga instead
       classifies each radius against its own and concatenates three clouds,
       which its per-radius normalisation makes depend on the *free* volumes
       alone. One surface and one share is the statement the parameter makes. */
    map->apply_contact_weighting(counts, get_contact_volume_thickness(),
                                 get_contact_volume_trapped_fraction());
    (void) data;

    st.have_result = true;
    st.last_source = source;
    st.generation1 = st.pending_gen1;
    st.generation2 = st.pending_gen2;
    st.result_generation++;
    st.quad_valid = false;

    // the cloud and the mean position, from the map and the cached source
    map->get_xyz_density_soa(st.cloud_x, st.cloud_y, st.cloud_z, st.cloud_w);
    st.cloud_generation = st.result_generation;
    st.cloud_valid = true;
    IMP::algebra::Vector3D r = source;
    /* One, not two: the source is in the numerator once, so it must be in the
       denominator once. The same phantom unit of weight as in
       get_mean_position(), in a second copy that sets the AV's coordinates. */
    double sum = 1.0;
    for(size_t i = 0; i < st.cloud_x.size(); i++){
        const double w = st.cloud_w[i];
        if(w <= 0.0f) continue;
        sum += w;
        r[0] += (double) st.cloud_x[i] * w;
        r[1] += (double) st.cloud_y[i] * w;
        r[2] += (double) st.cloud_z[i] * w;
    }
    st.last_mean = r / sum;
    st.last_compute_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - st.compute_t0).count();
    st.pending_stage = 2;
}

void AV::resample_lattice_finish(){
    auto &st = *state_;
    if(!st.pending) return;
    st.pending = false;
    if(st.pending_shift_xyz){
        set_coordinates(st.last_mean);
    }
}

void AV::set_av_parameter(const std::string &json_text){
    const nlohmann::json j = nlohmann::json::parse(json_text);
    // Only AV1/AV3 (and XYZ, which the reader turns into a fixed point) are
    // AV models. Rotamer-ensemble positions ("R1", PRD-108) are Python-only
    // today: an fps.json handed to ProbeNetworkRestraint must be filtered with
    // IMP.bff.fps_positions_for_docking() first -- warn instead of silently
    // building an AV1 with default parameters for them.
    const std::string stype = j.value("simulation_type", std::string("AV1"));
    if(stype != "AV1" && stype != "AV3" && stype != "XYZ"){
        IMP_WARN("AV: fps.json position has simulation_type '" << stype
                 << "', which is not an AV model (AV1/AV3/XYZ); it is scored "
                 << "as AV1 with its AV parameters. Rotamer-ensemble positions "
                 << "(R1) are Python-only: filter them with "
                 << "IMP.bff.fps_positions_for_docking() before docking.\n");
    }
    set_linker_length(j.value("linker_length", 20.0));
    algebra::Vector3D r = {j.value("radius1", 3.0),
                           j.value("radius2", 0.0),
                           j.value("radius3", 0.0)};
    set_radius1(r[0]);
    set_radius2(r[1]);
    set_radius3(r[2]);
    set_linker_width(j.value("linker_width", 0.5));
    // Absent means *derive* (see get_effective_allowed_sphere_radius); a
    // flat 1.5 here is what made this door return empty volumes at FPS's
    // standard linker width of 4.5 A.
    set_allowed_sphere_radius(j.value("allowed_sphere_radius", -1.0));
    set_contact_volume_thickness(j.value("contact_volume_thickness", 0.0));
    /* `-1.0`, not `-1`. nlohmann deduces the value type from the default, so
       an `int` default read `"contact_volume_trapped_fraction": 0.328` back as
       **0** -- every fps.json fraction truncated to zero or one. Invisible
       while the ACV was inert (PRD-121 G9); the first thing the wiring found. */
    set_contact_volume_trapped_fraction(
            j.value("contact_volume_trapped_fraction", -1.0));
    set_simulation_grid_resolution(j.value("simulation_grid_resolution", 1.5));
    // The position's parameters are calibrated against a structure whose
    // labelling site has been stripped, so the mask is part of the position,
    // not a convenience of whoever wrote the file.
    set_strip_mask(j.value("strip_mask", std::string()));
    // Off unless the file asks for it: turning it on changes every number the
    // volume reports, so it is opt-in and not a default.
    set_chain_weighting(j.value("chain_weighting", false));
    /* Which van der Waals radii the obstacles are inflated by. Absent means
       the module default, `"imp"`: the radii the particles carry, which is
       what `clash_container` reads, so a volume and the clash term beside it
       cannot disagree about how big an atom is. A file reproducing Olga-era
       numbers -- the fitted `contact_volume_trapped_fraction` above was fitted
       against Olga's obstacle set, as were the published <R_DA> this module is
       measured against -- says `"radii_source": "olga"`. It is a property of
       the position and not a run flag for the same reason the strip mask is:
       the rest of the position's parameters were calibrated against one radii
       set, and a file that names it says which. */
    set_radii_source(j.value("radii_source",
                             av_radii_source_to_string(AV_RADII_IMP)));
}

IMP::bff::PathMapHeader AV::create_path_map_header(){
    // create PathMapHeader
    double ll = get_effective_linker_length();
    double dg = get_simulation_grid_resolution();
    IMP::bff::PathMapHeader path_map_header(ll, dg);
    if(get_space_fixed()){
        int k0[3]; int n;
        lattice_window(get_source_coordinates(), ll, dg, k0, n);
        path_map_header.update_map_dimensions(n, n, n);
        if(get_search_stencil() == 26){
            // face, edge, corner neighbours: no length-2 axis jumps, so a
            // path cannot tunnel through a one-voxel wall
            path_map_header.set_neighbor_radius(std::sqrt(3.0) + 1e-6);
        } else if(get_search_stencil() == 74){
            // The reference metric. LabelLib's `essentialNeighbours()`
            // (`FlexLabel/src/FlexLabel.cxx:149`) keeps the shells with squared
            // offset length {1, 2, 3, 5, 6} -- 74 edges -- and notes that a
            // shorter list gives "isopath surfaces that are cubic instead of
            // spherical". The 26 stencil is {1, 2, 3}, and it costs real volume:
            // with no obstacles at all, where the answer must be a sphere, 26
            // returns 26146 voxels against LabelLib's 30688 for a 20 A linker.
            //
            // A radius of sqrt(6) admits {1, 2, 3, 4, 5, 6} = 80 edges. The six
            // extra are the squared-length-4 axis jumps (2,0,0), which are
            // exactly redundant -- their cost 2.0 equals two unit steps -- so
            // Dijkstra returns the same path lengths as LabelLib's 74. They are
            // left in rather than special-cased: the stencil stays a radius
            // rule, and the cost is six offsets.
            path_map_header.set_neighbor_radius(std::sqrt(6.0) + 1e-6);
        }
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
        const std::vector<float> &cx, const std::vector<float> &cy,
        const std::vector<float> &cz, const std::vector<float> &cw,
        double spacing, int k,
        std::vector<IMP::algebra::Vector4D> &points,
        std::vector<std::array<double, 6> > &moments){
    points.clear();
    moments.clear();
    if(cx.empty() || k <= 0) return;
    const size_t np = cx.size();
    if((int) np <= k){
        points.reserve(np);
        for(size_t i = 0; i < np; i++){
            points.emplace_back((double) cx[i], (double) cy[i], (double) cz[i], (double) cw[i]);
        }
        moments.assign(np, {0, 0, 0, 0, 0, 0});
        return;
    }
    // integer lattice indices relative to the cloud minimum
    std::vector<int> ix(np), iy(np), iz(np);
    double minx = cx[0], miny = cy[0], minz = cz[0];
    for(size_t i = 0; i < np; i++){
        minx = std::min(minx, (double) cx[i]);
        miny = std::min(miny, (double) cy[i]);
        minz = std::min(minz, (double) cz[i]);
    }
    int ex = 0, ey = 0, ez = 0;   // extents in voxels
    for(size_t i = 0; i < np; i++){
        ix[i] = (int) std::floor(((double) cx[i] - minx) / spacing + 0.5);
        iy[i] = (int) std::floor(((double) cy[i] - miny) / spacing + 0.5);
        iz[i] = (int) std::floor(((double) cz[i] - minz) / spacing + 0.5);
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
        double w = cw[i];
        sw[b] += w; sx[b] += w * (double) cx[i]; sy[b] += w * (double) cy[i]; sz[b] += w * (double) cz[i];
    }
    points.resize(nb);
    for(int b = 0; b < nb; b++){
        points[b] = IMP::algebra::Vector4D(sx[b] / sw[b], sy[b] / sw[b], sz[b] / sw[b], sw[b]);
    }
    moments.assign(nb, {0, 0, 0, 0, 0, 0});
    for(size_t i = 0; i < np; i++){
        int b = slot[cell[i]];
        double w = cw[i];
        double dx = (double) cx[i] - points[b][0], dy = (double) cy[i] - points[b][1], dz = (double) cz[i] - points[b][2];
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
        av.refresh_cloud_soa();
        double h = map->get_path_map_header().get_simulation_grid_resolution();
        coarsen_cloud(st.cloud_x, st.cloud_y, st.cloud_z, st.cloud_w, h, k, st.quad_points, st.quad_moments);
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
        case PROBE_PAIR_DISTANCE_MP:
        case PROBE_PAIR_XYZ_DISTANCE:
        case PROBE_PAIR_DISTANCE_MIN:
            // All three are exact already; there is nothing for a quadrature
            // to approximate, and a coarsened cloud would give a *wrong*
            // minimum rather than a cheaper one.
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
    const bool efficiency = (distance_type == PROBE_PAIR_EFFICIENCY ||
                             distance_type == PROBE_PAIR_DISTANCE_E);
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
                double x = d / forster_radius;
                double x2 = x * x;
                double t = x2 * x2 * x2;
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
    if(distance_type == PROBE_PAIR_DISTANCE_E){
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
    IMP_LOG_VERBOSE("labelling site " << chain_id << ":" << residue_seq_number
                    << ":" << atom_name << " -> " << p_chain.size()
                    << " chain, " << p_residue.size() << " residue, "
                    << p_atom.size() << " atom matches" << std::endl);

    /* Emptiness first. This read `p_residue[0]` and *then* asked whether
       `p_residue` had exactly one element -- and IMP_USAGE_CHECK compiles out
       of a release build, so a site the structure does not contain indexed
       element zero of an empty vector. Screening a library segfaulted on the
       first structure that was missing a chain, which is the ordinary case a
       library is screened for. Throwing is what the callers already expect:
       IMP::bff::screen_structures catches IMP::Exception and records the
       structure with a NaN rather than ending the run. */
    if(p_residue.empty()){
        IMP_THROW("no labelling site " << (chain_id.empty() ? "*" : chain_id)
                  << ":" << residue_seq_number
                  << (atom_name.empty() ? "" : ":" + atom_name)
                  << " in this structure", ValueException);
    }
    if(p_residue.size() > 1){
        IMP_THROW("ambiguous labelling site " << (chain_id.empty() ? "*" : chain_id)
                  << ":" << residue_seq_number << ": " << p_residue.size()
                  << " residues match", ValueException);
    }
    auto p = p_residue[0];
    if(p_atom.empty()){
        IMP_WARN("labelling site " << chain_id << ":" << residue_seq_number
                 << " has no atom " << atom_name
                 << "; attaching at residue resolution" << std::endl);
    } else{
        p = p_atom[0];
    }
    return p->get_index();
}

//! Random sampling over AV
    void get_xyz_density();

namespace {
// The samplers, as plain vectors. The exported entry points below publish a
// numpy view over a copy; these two are what the other kernels in this file
// call, and they must not allocate a view only to free it.
std::vector<double> sample_points(const AV& av, int n_samples){
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

std::vector<double> sample_distances(
        const AV& av1,
        const AV& av2,
        int n_samples
){
    auto p1 = sample_points(av1, n_samples);
    auto p2 = sample_points(av2, n_samples);
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
}  // namespace

static std::vector<double> av_cloud(const AV& av) {
    std::vector<double> flat;
    auto m = av.get_map();
    if (!m) return flat;
    const auto d = m->get_xyz_density();
    flat.reserve(d.size() * 4);
    for (const auto& v : d) {
        flat.push_back(v[0]); flat.push_back(v[1]);
        flat.push_back(v[2]); flat.push_back(v[3]);
    }
    return flat;
}

void write_av(const AV& av, const std::string& path) {
    const std::size_t dot = path.rfind('.');
    std::string ext = dot == std::string::npos ? "" : path.substr(dot + 1);
    for (auto& c : ext) c = (char) std::tolower((unsigned char) c);

    if (ext == "xyz") {
        write_points_xyz(path, av_cloud(av), "He",
                         "accessible volume, column 5 is the weight");
        return;
    }
    if (ext == "pqr") {
        write_points_pqr(path, av_cloud(av), 1.0);
        return;
    }
    auto m = av.get_map();
    if (!m) IMP_THROW("the volume has no map to write", ValueException);
    if (ext == "dx") {
        const auto* h = m->get_header();
        const int nx = h->get_nx(), ny = h->get_ny(), nz = h->get_nz();
        // The *accessible* density, which is what an accessible volume is:
        // PathMap derives from em::SampledDensityMap and its own `get_value`
        // is the obstacle raster it was built from, not the volume. This is
        // the same feature `write_map_feature` and AVBuilder read, so it also
        // carries the linker weighting when one is set.
        const std::vector<float> tile = m->get_tile_values(
                PM_TILE_ACCESSIBLE_DENSITY);
        std::vector<double> density((std::size_t) nx * ny * nz, 0.0);
        // IMP's voxels run x fastest; OpenDX wants z fastest, so this is the
        // transposition, done once here rather than by every caller.
        for (int ix = 0; ix < nx; ++ix) {
            for (int iy = 0; iy < ny; ++iy) {
                for (int iz = 0; iz < nz; ++iz) {
                    const long v = m->xyz_ind2voxel(ix, iy, iz);
                    if (v < 0 || (std::size_t) v >= tile.size()) continue;
                    density[((std::size_t) ix * ny + iy) * nz + iz] = tile[v];
                }
            }
        }
        std::vector<double> origin(3);
        origin[0] = h->get_xorigin();
        origin[1] = h->get_yorigin();
        origin[2] = h->get_zorigin();
        write_opendx(path, density, nx, ny, nz, origin, h->get_spacing());
        return;
    }
    if (ext == "mrc" || ext == "map" || ext == "ccp4") {
        write_map_feature(m, path, PM_TILE_DENSITY);
        return;
    }
    IMP_THROW("cannot write " << path << ": the extension must be one of "
                              << "xyz, pqr, dx, mrc, map, ccp4",
              ValueException);
}

double rmp_from_measurement(const AV& a, const AV& b,
                            const AVPairDistanceMeasurement& measurement,
                            double accuracy, int n_samples, int seed) {
    return rmp_from_model_distance(av_cloud(a), av_cloud(b),
                                   measurement.distance,
                                   measurement.distance_type,
                                   measurement.forster_radius, accuracy,
                                   n_samples, seed);
}

std::vector<double> rmp_flat_bottom_bounds(
        const AV& a, const AV& b, const AVPairDistanceMeasurement& measurement,
        double accuracy, int n_samples, int seed) {
    const std::vector<double> ca = av_cloud(a), cb = av_cloud(b);
    const int type = measurement.distance_type;
    const double r0 = measurement.forster_radius;
    // An efficiency falls with separation, so its low error bar is the high
    // R_mp bound. Converting each end separately and sorting afterwards keeps
    // that straight without a special case per convention.
    const double lo_obs = measurement.distance -
                          std::max(0.0, measurement.error_neg);
    const double hi_obs = measurement.distance +
                          std::max(0.0, measurement.error_pos);
    auto convert = [&](double target, double fallback) {
        try {
            return rmp_from_model_distance(ca, cb, target, type, r0, accuracy,
                                           n_samples, seed);
        } catch (const IMP::ValueException&) {
            // An error bar running past what the volumes can reach is a
            // one-sided restraint, not a failure.
            return fallback;
        }
    };
    const double centre = rmp_from_model_distance(ca, cb, measurement.distance,
                                                  type, r0, accuracy,
                                                  n_samples, seed);
    double lo = convert(lo_obs, 0.0);
    double hi = convert(hi_obs, std::numeric_limits<double>::infinity());
    if (lo > hi) std::swap(lo, hi);
    std::vector<double> out(3);
    out[0] = lo; out[1] = centre; out[2] = hi;
    return out;
}

double av_overlap(const AV& av, const IMP::atom::Hierarchy& hierarchy,
                  const std::string& selection, double radius) {
    IMP::ParticleIndexes chosen;
    if (selection.empty()) {
        for (auto h : IMP::atom::get_leaves(hierarchy)) {
            chosen.push_back(h.get_particle_index());
        }
    } else {
        chosen = selection_from_expression(hierarchy, selection)
                         .get_selected_particle_indexes();
    }
    IMP::Model* m = hierarchy.get_model();
    std::vector<double> refs;
    refs.reserve(chosen.size() * 3);
    for (auto pi : chosen) {
        if (!IMP::core::XYZ::get_is_setup(m, pi)) continue;
        const auto c = IMP::core::XYZ(m, pi).get_coordinates();
        refs.push_back(c[0]); refs.push_back(c[1]); refs.push_back(c[2]);
    }
    return cloud_overlap(av_cloud(av), refs, radius);
}

void av_random_points(const AV& av, double** out_view, int* n_out_view,
                      int n_samples){
    internal::copy_to_view(sample_points(av, n_samples), out_view, n_out_view);
}

void av_random_distances(const AV& av1, const AV& av2,
                         double** out_view, int* n_out_view, int n_samples){
    internal::copy_to_view(sample_distances(av1, av2, n_samples),
                           out_view, n_out_view);
}

std::vector<double> av_distance_distribution(
        const AV& av1,
        const AV& av2,
        //double start, double stop, int n_bins,
        std::vector<double> axis,
        int n_samples
){
    auto data = sample_distances(av1, av2, n_samples);
 
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
