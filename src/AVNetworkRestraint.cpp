/**
 *  \file IMP/bff/AVNetworkRestraint.h
 *  \brief Simple restraint for networks of accessible volumes.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */
 #include <IMP/bff/AVNetworkRestraint.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <thread>
#include <vector>
#include <memory>
#include <functional>

IMPBFF_BEGIN_NAMESPACE

AVNetworkRestraint::AVNetworkRestraint(
        const IMP::core::Hierarchy &hier,
        std::string fps_json_fn,
        std::string name,
        std::string score_set,
        int n_samples,
        bool space_fixed,
        bool shared_map,
        std::string distance,
        int quad_k,
        int search_grid_factor,
        int search_stencil,
        std::string search_mode
) : IMP::Restraint(hier.get_model(), name), n_samples(n_samples),
    space_fixed_(space_fixed), shared_map_(shared_map),
    distance_(distance), quad_k_(quad_k), search_grid_factor_(search_grid_factor),
    search_stencil_(search_stencil), search_mode_(search_mode){
    if(search_mode != "dijkstra" && search_mode != "euclidean"){
        IMP_THROW("AVNetworkRestraint: search_mode must be \"dijkstra\" or \"euclidean\"",
                  IMP::ValueException);
    }
    if(search_mode == "euclidean" && !space_fixed){
        IMP_THROW("AVNetworkRestraint: search_mode=\"euclidean\" requires space_fixed=True",
                  IMP::ValueException);
    }
    if(search_stencil != 26 && search_stencil != 30){
        IMP_THROW("AVNetworkRestraint: search_stencil must be 26 or 30", IMP::ValueException);
    }
    if(search_grid_factor < 1){
        IMP_THROW("AVNetworkRestraint: search_grid_factor must be >= 1", IMP::ValueException);
    }
    if(search_grid_factor > 1 && !space_fixed){
        IMP_THROW("AVNetworkRestraint: search_grid_factor > 1 requires space_fixed=True",
                  IMP::ValueException);
    }
    if(shared_map && !space_fixed){
        IMP_THROW("AVNetworkRestraint: shared_map=True requires space_fixed=True "
                  "(sharing needs commensurate lattice windows)",
                  IMP::ValueException);
    }
    if(distance != "quad" && distance != "mc"){
        IMP_THROW("AVNetworkRestraint: distance must be \"quad\" or \"mc\", got \""
                  << distance << "\"", IMP::ValueException);
    }
    if(quad_k < 1){
        IMP_THROW("AVNetworkRestraint: quad_k must be >= 1", IMP::ValueException);
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
        av.second->set_search_grid_factor(search_grid_factor_);
        av.second->set_search_stencil(search_stencil_);
        av.second->set_search_mode(search_mode_);
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

namespace internal {
struct AVEvalJob {
    struct RasterTask { AVOccupancyMap *m; int z_lo, z_hi; bool local; };
    std::vector<RasterTask> rtasks;
    AVOccupancyMaps maps;
    std::vector<IMP::bff::AV*> all;
    std::vector<const AVPairDistanceMeasurement*> pairs;
    std::vector<double> model;
    std::vector<std::pair<size_t, size_t> > pair_slots;
    bool pipelined = false;
    bool quad = true;
    int qk = 50;
    bool has_model_reading_pairs = false;
    std::unique_ptr<std::thread> runner;   // evaluate_async: the pool's caller thread
    std::chrono::steady_clock::time_point t0;
    double t_registry = 0, t_prepare = 0, t_compute = 0;
    ~AVEvalJob() {
        // an abandoned evaluate_async(): let it finish rather than terminate
        if (runner && runner->joinable()) runner->join();
    }
};
}

std::shared_ptr<internal::AVEvalJob> AVNetworkRestraint::begin_evaluation() const {
    auto job = std::make_shared<internal::AVEvalJob>();
    n_evaluations_++;
    typedef std::chrono::steady_clock clk;
    auto lap = [](clk::time_point &t){ auto n = clk::now(); double d = std::chrono::duration<double>(n - t).count(); t = n; return d; };
    clk::time_point t = clk::now();
    for(auto &av: avs_){
        av.second->prepare_lattice_window();
    }
    // The shared class rasters are independent of each other: they are
    // classified here (Model reads) and rasterised in the pool run.
    const int nthreads = get_number_of_threads();
    job->quad = (distance_ == "quad");
    job->qk = quad_k_;
    job->pipelined = (nthreads > 1);
    if(registry_){
        registry_->refresh_snapshot();
        job->maps = registry_->get_maps();
        if(job->pipelined){
            for(auto &mp : job->maps){
                int action = mp->begin_update();
                if(action == 1){
                    job->rtasks.push_back({mp.get(), 0, 0, true});
                } else if(action == 2){
                    std::vector<int> ext = mp->get_extent();
                    int z0 = ext[2], nz = ext[5];
                    int nslab = std::max(1, std::min(nthreads, nz));
                    for(int sIdx = 0; sIdx < nslab; sIdx++){
                        int a = z0 + (int)((long) nz * sIdx / nslab);
                        int b = z0 + (int)((long) nz * (sIdx + 1) / nslab) - 1;
                        if(b >= a) job->rtasks.push_back({mp.get(), a, b, false});
                    }
                }
            }
        } else {
            for(auto &m : job->maps) m->update();
        }
    }
    job->t_registry = lap(t);
    // prepare (serial: touches the Model; with a pipelined registry it only
    // reads the pending classification)
    job->all.reserve(avs_.size());
    for(auto &av: avs_){
        av.second->set_registry_driven_externally(job->pipelined && registry_);
        av.second->resample_prepare();
        job->all.push_back(av.second.get());
    }
    job->t_prepare = lap(t);
    std::stable_sort(job->all.begin(), job->all.end(), [](const IMP::bff::AV *a, const IMP::bff::AV *b){
        bool pa = a->get_has_pending_compute(), pb = b->get_has_pending_compute();
        if(pa != pb) return pa;
        return a->get_last_compute_seconds() > b->get_last_compute_seconds();
    });
    job->pairs.reserve(distances_.size());
    for(const auto &it : distances_){
        job->pairs.push_back(&it.second);
        if(it.second.distance_type == DYE_PAIR_DISTANCE_MP ||
           it.second.distance_type == DYE_PAIR_XYZ_DISTANCE){
            job->has_model_reading_pairs = true;
        }
    }
    job->model.assign(job->pairs.size(), 0.0);
    std::map<const IMP::bff::AV*, size_t> av_slot;
    for(size_t i = 0; i < job->all.size(); i++) av_slot[job->all[i]] = i;
    job->pair_slots.resize(job->pairs.size());
    for(size_t j = 0; j < job->pairs.size(); j++){
        job->pair_slots[j] = std::make_pair(av_slot[get_av(job->pairs[j]->position_1)],
                                            av_slot[get_av(job->pairs[j]->position_2)]);
    }
    job->t0 = t;
    return job;
}

void AVNetworkRestraint::run_evaluation(internal::AVEvalJob &job) const {
    // Rasters, searches (after all rasters), carves + quadrature (after the
    // own search), pair sums (after both AVs) -- one pool run with explicit
    // dependencies; nothing here reads the Model.
    const bool quad = job.quad;
    const int qk = job.qk;
    auto eval_pair = [&](size_t j){
        job.model[j] = get_model_distance(job.pairs[j]->position_1, job.pairs[j]->position_2,
                                          job.pairs[j]->forster_radius, job.pairs[j]->distance_type);
    };
    if(job.pipelined){
        const size_t nr = job.rtasks.size(), na = job.all.size(), np = job.pairs.size();
        std::atomic<size_t> rasters_left(nr);
        std::vector<std::atomic<int> > stage(na);   // 0 none, 1 searched, 2 ready
        for(auto &f : stage) f.store(0);
        auto spin_until = [](const std::function<bool()> &ok){
            while(!ok()) std::this_thread::yield();
        };
        auto work = [&](size_t t){
            if(t < nr){
                const internal::AVEvalJob::RasterTask &rt = job.rtasks[t];
                if(rt.local) rt.m->apply_local(); else rt.m->raster_slab(rt.z_lo, rt.z_hi);
                rasters_left.fetch_sub(1, std::memory_order_acq_rel);
            } else if(t < nr + na){
                size_t i = t - nr;
                spin_until([&]{ return rasters_left.load(std::memory_order_acquire) == 0; });
                job.all[i]->resample_compute_search();
                stage[i].store(1, std::memory_order_release);
            } else if(t < nr + 2 * na){
                size_t i = t - nr - na;
                spin_until([&]{ return stage[i].load(std::memory_order_acquire) >= 1; });
                job.all[i]->resample_compute_carve();
                if(quad) job.all[i]->prepare_quadrature(qk);
                stage[i].store(2, std::memory_order_release);
            } else {
                size_t j = t - nr - 2 * na;
                size_t a = job.pair_slots[j].first, b = job.pair_slots[j].second;
                spin_until([&]{ return stage[a].load(std::memory_order_acquire) == 2 &&
                                       stage[b].load(std::memory_order_acquire) == 2; });
                if(quad) eval_pair(j);
            }
        };
        // Order matters: rasters, searches, carves, pairs -- no worker waits
        // on a task that an earlier worker has not started.
        get_pool().run(nr + 2 * na + (quad ? np : 0), work);
        if(!quad){
            for(size_t j = 0; j < np; j++) eval_pair(j);
        }
    } else {
        for(size_t i = 0; i < job.all.size(); i++){
            job.all[i]->resample_compute();
            if(quad) job.all[i]->prepare_quadrature(qk);
        }
        for(size_t j = 0; j < job.pairs.size(); j++) eval_pair(j);
    }
}

double AVNetworkRestraint::finish_evaluation(internal::AVEvalJob &job) const {
    typedef std::chrono::steady_clock clk;
    for(auto &mp : job.maps) mp->end_update();
    for(auto *av : job.all) av->resample_finish();
    double score = 0.0;
    for(size_t i = 0; i < job.pairs.size(); i++){
        score += job.pairs[i]->score_model(job.model[i]);
    }
    t_registry_ += job.t_registry;
    t_prepare_ += job.t_prepare;
    t_compute_ += std::chrono::duration<double>(clk::now() - job.t0).count();
    return score;
}

double AVNetworkRestraint::unprotected_evaluate(
        IMP::DerivativeAccumulator *accum) const {
    if(job_){
        // an evaluate_async() nobody waited for: finish it first
        wait_score();
    }
    std::shared_ptr<internal::AVEvalJob> job = begin_evaluation();
    run_evaluation(*job);
    return finish_evaluation(*job);
}

void AVNetworkRestraint::evaluate_async() const {
    if(job_) wait_score();
    std::shared_ptr<internal::AVEvalJob> job = begin_evaluation();
    // Pairs that read the Model at pair time (MP / XYZ types) cannot overlap a
    // Model change: run those synchronously here. The pool's caller must be
    // a thread of its own for the overlap; the workers are the pool's.
    if(!job->pipelined || job->has_model_reading_pairs){
        run_evaluation(*job);
        job_ = job;
        return;
    }
    internal::AVEvalJob *raw = job.get();
    job->runner.reset(new std::thread([this, raw](){ run_evaluation(*raw); }));
    job_ = job;
}

double AVNetworkRestraint::wait_score() const {
    IMP_USAGE_CHECK(job_, "AVNetworkRestraint::wait_score without evaluate_async");
    std::shared_ptr<internal::AVEvalJob> job = job_;
    if(job->runner){
        job->runner->join();
        job->runner.reset();
    }
    double score = finish_evaluation(*job);
    job_.reset();
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
    if(distance_ == "quad"){
        return av_distance_quadrature(*av1, *av2, forster_radius, distance_type, quad_k_);
    }
    return av_distance(*av1, *av2, forster_radius,distance_type, n_samples);
}

internal::ThreadPool &AVNetworkRestraint::get_pool() const{
    int n = get_number_of_threads();
    if(!pool_ || pool_->size() != n){
        pool_ = std::make_shared<internal::ThreadPool>(n);
    }
    return *pool_;
}

int AVNetworkRestraint::get_number_of_threads() const{
    if(n_threads_ > 0) return n_threads_;
    unsigned hw = std::thread::hardware_concurrency();
    return hw > 0 ? (int) hw : 1;
}

std::string AVNetworkRestraint::get_diagnostics_json() const{
    nlohmann::json j;
    j["space_fixed"] = space_fixed_;
    j["shared_map"] = shared_map_;
    j["distance"] = distance_;
    j["quad_k"] = quad_k_;
    j["search_grid_factor"] = search_grid_factor_;
    j["search_stencil"] = search_stencil_;
    j["search_mode"] = search_mode_;
    j["n_samples"] = n_samples;
    j["evaluations"] = n_evaluations_;
    j["threads"] = get_number_of_threads();
    j["timing_ms_total"] = {{"registry", t_registry_ * 1e3}, {"prepare", t_prepare_ * 1e3},
                            {"compute", t_compute_ * 1e3}, {"pairs", t_pairs_ * 1e3}};
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

double AVNetworkRestraint::get_quad_error_estimate(int reference_k) const{
    double worst = 0.0;
    for(const auto &it : distances_){
        const auto &d = it.second;
        auto av1 = get_av(d.position_1);
        auto av2 = get_av(d.position_2);
        double a = av_distance_quadrature(*av1, *av2, d.forster_radius,
                                          d.distance_type, quad_k_);
        double b = av_distance_quadrature(*av1, *av2, d.forster_radius,
                                          d.distance_type, reference_k);
        if(std::isnan(a) || std::isnan(b)) continue;
        worst = std::max(worst, std::fabs(a - b));
    }
    return worst;
}

IMP_OBJECT_SERIALIZE_IMPL(IMP::bff::AVNetworkRestraint);

IMPBFF_END_NAMESPACE
