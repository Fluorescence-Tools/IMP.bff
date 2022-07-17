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
        std::string score_set
) : IMP::Restraint(hier.get_model(), name), hier_(hier) {
    auto fps_reader = IMP::bff::FPSReaderWriter(fps_json_fn, score_set);

    distances_ = fps_reader.get_distances();
    auto used_positions = fps_reader.get_used_positions();
    create_av_decorated_particles(used_positions);

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

void AVNetworkRestraint::create_av_decorated_particles(
        nlohmann::json used_positions
){
    std::cout << "create_av_decorated_particles" << std::endl;
    IMP::Model* model = get_model();

    for(nlohmann::json::iterator it = used_positions.begin();
            it != used_positions.end(); ++it){
        nlohmann::json position = it.value();
        std::string position_name = it.key();

        std::cout << position << std::endl;

        // Create new Particle for AV
        IMP_NEW(IMP::Particle, av_particle, (model));
        av_particle->set_name(position_name);
        IMP::ParticleIndex av_index = av_particle->get_index();

        // Search for labeling site
        IMP::ParticleIndex parent_particle_idx =
                IMP::bff::search_labeling_site(hier_, "", position);

        // Store AV particle index
        av_pi_.emplace_back(av_index);

        // Decorate AV particle
        IMP::bff::AV::do_setup_particle(model, av_index, parent_particle_idx);
        auto av = new IMP::bff::AV(av_particle);
        av->set_av_parameter(position);

        avs_.emplace_back(av);
    }
}


IMP::bff::AV* AVNetworkRestraint::get_av(std::string name) const{
    for(auto av: avs_){
        if(av->get_particle()->get_name() == name){
            return av;
        }
    }
    return nullptr;
}


double AVNetworkRestraint::unprotected_evaluate(
        IMP::DerivativeAccumulator *accum) const{
    double score = 0.0;
    for(auto &av: avs_){
        av->resample();
    }
    for(auto it = distances_.begin(); it != distances_.end(); it++){
        auto distance = it->second;
        auto av1 = get_av(distance.position_1);
        auto av2 = get_av(distance.position_2);
        double feature = av_distance(
                *av1, *av2, distance.forster_radius,
                distance.distance_type);
        score += distance.score_model(feature);
    }
    return score;
}

IMPBFF_END_NAMESPACE
