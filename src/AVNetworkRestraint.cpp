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
        int n_samples
) : IMP::Restraint(hier.get_model(), name), n_samples(n_samples){
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

std::map<std::string, IMP::bff::AV*> AVNetworkRestraint::create_av_decorated_particles(
        nlohmann::json used_positions,
        const IMP::core::Hierarchy &hier
){
    IMP::Model* model = get_model();
    std::map<std::string, IMP::bff::AV*> avs{};

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
        auto av = new IMP::bff::AV(av_particle);
        av->set_av_parameter(position);

        avs[position_name] = av;
    }
    return avs;
}

IMP::bff::AV* AVNetworkRestraint::get_av(std::string name) const{
    for (const auto& n : avs_)
        if(n.first == name){
            return n.second;
        }
    IMP_WARN("AV not found in AVNetworkRestraint");
    return nullptr;
}

double AVNetworkRestraint::unprotected_evaluate(
        IMP::DerivativeAccumulator *accum) const {
    double score = 0.0;
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

IMPBFF_END_NAMESPACE
