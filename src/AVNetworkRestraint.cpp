/**
 *  \file IMP/bff/AVNetworkRestraint.h
 *  \brief Simple restraint for networks of accessible volumes.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */
 #include <IMP/bff/AVNetworkRestraint.h>

IMPBFF_BEGIN_NAMESPACE

AVNetworkRestraint::AVNetworkRestraint(
        const IMP::core::Hierarchy &hier,
        std::string fps_json_fn,
        std::string name,
        std::string score_set,
        int n_samples,
        bool space_fixed,
        bool shared_map
) : IMP::Restraint(hier.get_model(), name), n_samples(n_samples),
    space_fixed_(space_fixed), shared_map_(shared_map){
    if(shared_map && !space_fixed){
        IMP_THROW("AVNetworkRestraint: shared_map=True requires space_fixed=True "
                  "(sharing needs commensurate lattice windows)",
                  IMP::ValueException);
    }
    if(!space_fixed){
        IMP_WARN("AVNetworkRestraint: space_fixed=False (legacy source-anchored "
                 "grids) is deprecated and will be removed after the PRD-105 "
                 "transition period.\n");
    }
    auto fps_reader = IMP::bff::FPSReaderWriter(fps_json_fn, score_set);

    distances_ = fps_reader.get_distances();
    auto used_positions = fps_reader.get_used_positions();

    avs_ = create_av_decorated_particles(used_positions, hier);
    for(auto &av: avs_){
        av_pi_.emplace_back(av.second->get_particle_index());
    }

    for(IMP::core::Hierarchy &h : IMP::core::get_leaves(hier)){
        model_ps_.emplace_back(h.get_particle_index());
    }
    configure_avs();
}

void AVNetworkRestraint::configure_avs(){
    if(space_fixed_ && shared_map_ && !registry_ && !avs_.empty()){
        // The obstacle set every AV rasterises: all leaves of the root of
        // the first source (identical for every AV of one hierarchy).
        IMP::bff::AV *first = avs_.begin()->second.get();
        IMP::Particle* parent = first->get_source();
        auto h = IMP::atom::Hierarchy(get_model(), parent->get_index());
        auto root = IMP::atom::get_root(h);
        registry_ = new AVOccupancyRegistry(IMP::atom::get_leaves(root));
        registry_->set_was_used(true);
    }
    for(auto &av: avs_){
        av.second->set_space_fixed(space_fixed_);
        av.second->set_occupancy_registry(
            (space_fixed_ && shared_map_) ? registry_.get() : nullptr);
    }
}

const IMP::bff::AVs AVNetworkRestraint::get_used_avs(){
    IMP::bff::AVs out;
    for(IMP::ParticleIndex pi: av_pi_){
        out.emplace_back(get_model(), pi);
    }
    return out;
}

IMP::ModelObjectsTemp AVNetworkRestraint::do_get_inputs() const {
    IMP::ModelObjectsTemp ret;
    for (size_t i = 0; i < model_ps_.size(); i++) {
        ret.push_back(get_model()->get_particle(model_ps_[i]));
    }
    return ret;
}

std::map<std::string, std::unique_ptr<IMP::bff::AV> > AVNetworkRestraint::create_av_decorated_particles(
        nlohmann::json used_positions,
        const IMP::core::Hierarchy &hier
){
    IMP::Model* model = get_model();
    std::map<std::string, std::unique_ptr<IMP::bff::AV> > avs{};

    for(nlohmann::json::iterator it = used_positions.begin();
            it != used_positions.end(); ++it){
        nlohmann::json position = it.value();
        std::string position_name = it.key();

        // Create new Particle for AV
        IMP_NEW(IMP::Particle, av_particle, (model));
        av_particle->set_name(position_name);

        // Search for labeling site
        IMP::ParticleIndex parent_particle_idx =
                IMP::bff::search_labeling_site(hier, "", position);

        // Store AV particle index
        IMP::ParticleIndex av_index = av_particle->get_index();

        // Decorate AV particle
        IMP::bff::AV::do_setup_particle(model, av_index, parent_particle_idx);
        auto av = new IMP::bff::AV(av_particle);  // ownership passes to avs below
        av->set_av_parameter(position);

        avs[position_name].reset(av);
    }
    return avs;
}

IMP::bff::AV AVNetworkRestraint::get_used_av(std::string name) const{
    IMP::bff::AV* av = get_av(name);
    IMP_USAGE_CHECK(av != nullptr, "AVNetworkRestraint: no AV named " << name);
    return *av;
}

IMP::bff::AV* AVNetworkRestraint::get_av(std::string name) const{
    for (const auto& n : avs_)
        if(n.first == name){
            return n.second.get();
        }
    IMP_WARN("AV not found in AVNetworkRestraint");
    return nullptr;
}

double AVNetworkRestraint::unprotected_evaluate(
        IMP::DerivativeAccumulator *accum) const {
    double score = 0.0;
    n_evaluations_++;
    for(auto &av: avs_){
        av.second->prepare_lattice_window();
    }
    for(auto &av: avs_){
        av.second->resample();
    }
    for(const auto & it : distances_){
        auto distance = it.second;
        double model = get_model_distance(
                distance.position_1,
                distance.position_2,
                distance.forster_radius,
                distance.distance_type
        );
        score += distance.score_model(model);
    }
    return score;
}

double AVNetworkRestraint::get_model_distance(
        std::string position1_name,
        std::string position2_name,
        double forster_radius,
        int distance_type
) const {
    auto av1 = get_av(position1_name);
    auto av2 = get_av(position2_name);
    return av_distance(*av1, *av2, forster_radius,distance_type, n_samples);
}

std::string AVNetworkRestraint::get_diagnostics_json() const{
    nlohmann::json j;
    j["space_fixed"] = space_fixed_;
    j["shared_map"] = shared_map_;
    j["n_samples"] = n_samples;
    j["evaluations"] = n_evaluations_;
    nlohmann::json avs = nlohmann::json::object();
    long skip = 0, local = 0, full = 0, rolls = 0;
    for(const auto &kv : avs_){
        const IMP::bff::AV &av = *kv.second;
        nlohmann::json a;
        a["skip"] = av.get_number_of_skips();
        a["local"] = av.get_number_of_local_updates();
        a["full"] = av.get_number_of_full_updates();
        a["rolls"] = av.get_number_of_rolls();
        a["window"] = av.get_lattice_window();
        skip += av.get_number_of_skips();
        local += av.get_number_of_local_updates();
        full += av.get_number_of_full_updates();
        rolls += av.get_number_of_rolls();
        avs[kv.first] = a;
    }
    j["avs"] = avs;
    j["av_totals"] = {{"skip", skip}, {"local", local}, {"full", full}, {"rolls", rolls}};
    nlohmann::json maps = nlohmann::json::array();
    auto add_map = [&](const AVOccupancyMap *m){
        nlohmann::json o;
        o["spacing"] = m->get_spacing();
        o["extra_radius"] = m->get_extra_radius();
        o["skip"] = m->get_number_of_skips();
        o["local"] = m->get_number_of_local_updates();
        o["full"] = m->get_number_of_full_updates();
        o["grow"] = m->get_number_of_grows();
        o["rolls"] = m->get_number_of_rolls();
        o["moved_last"] = m->get_number_of_moved_last();
        o["moved_total"] = m->get_number_of_moved_total();
        o["extent"] = m->get_extent();
        o["voxels"] = m->get_number_of_voxels();
        maps.push_back(o);
    };
    if(registry_){
        for(const auto &m : registry_->get_maps()) add_map(m.get());
    }
    j["shared_maps"] = maps;
    j["shared_map_classes"] = maps.size();
    return j.dump();
}

IMP_OBJECT_SERIALIZE_IMPL(IMP::bff::AVNetworkRestraint);

IMPBFF_END_NAMESPACE
