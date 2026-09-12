/** \file IMP/bff/ModelSearch.cpp */
#include <IMP/bff/ModelSearch.h>
#include <IMP/bff/FitMinimizer.h>
#include <IMP/bff/GraphNode.h>
#include <IMP/bff/GraphPort.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <utility>

IMPBFF_BEGIN_NAMESPACE

namespace {

void require_key(const std::string& key, const char* what) {
  if (key.empty()) {
    throw ModelSearchConfigurationError(std::string(what) +
                                           " key must not be empty");
  }
}

void require_state(const ModelSearchState& state, const char* what) {
  require_key(state.get_key(), what);
  require_key(state.get_structure_key(), "state structure");
  if (!std::isfinite(state.get_reward())) {
    throw ModelSearchConfigurationError(std::string(what) +
                                           " reward must be finite");
  }
}

struct TreeNode {
  explicit TreeNode(double p = 1.0, TreeNode* up = nullptr)
      : prior(p), parent(up) {}

  ModelSearchState state;
  ModelSearchAction parent_action;
  double prior = 1.0;
  std::vector<std::unique_ptr<TreeNode> > children;
  int visits = 0;
  double value_sum = 0.0;
  bool expanded = false;
  bool evaluated = false;
  bool terminal = false;
  TreeNode* parent = nullptr;
  std::string structure_key;
};

double reward_value(const ModelSearchState& state, double root_reward,
                    double scale) {
  return std::tanh((state.get_reward() - root_reward) / scale);
}

double puct(const TreeNode& child, int parent_visits, double c_puct) {
  const double q = child.value_sum / (child.visits + 1.0);
  const double u = c_puct * child.prior *
                   std::sqrt(static_cast<double>(parent_visits)) /
                   (1.0 + child.visits);
  return q + u;
}

}  // namespace

ModelSearchState::ModelSearchState()
    : reward_(0.0), acceptable_(false) {}

ModelSearchState::ModelSearchState(const std::string& key, double reward,
                                         bool acceptable)
    : key_(key), structure_key_(key), reward_(reward),
      acceptable_(acceptable) {}

ModelSearchState::ModelSearchState(
    const std::string& key, const std::string& structure_key, double reward,
    bool acceptable)
    : key_(key), structure_key_(structure_key), reward_(reward),
      acceptable_(acceptable) {}

const std::string& ModelSearchState::get_key() const { return key_; }
const std::string& ModelSearchState::get_structure_key() const {
  return structure_key_;
}
double ModelSearchState::get_reward() const { return reward_; }
bool ModelSearchState::get_acceptable() const { return acceptable_; }

ModelSearchAction::ModelSearchAction()
    : prior_(1.0), terminal_(false) {}

ModelSearchAction::ModelSearchAction(
    const std::string& key, const std::string& predicted_state_key,
    double prior, bool terminal)
    : key_(key), predicted_state_key_(predicted_state_key), prior_(prior),
      terminal_(terminal) {}

const std::string& ModelSearchAction::get_key() const { return key_; }
const std::string& ModelSearchAction::get_predicted_state_key() const {
  return predicted_state_key_;
}
double ModelSearchAction::get_prior() const { return prior_; }
bool ModelSearchAction::get_terminal() const { return terminal_; }

ModelSearchProblem::~ModelSearchProblem() {}
void ModelSearchProblem::request_cancel() {}
void ModelSearchProblem::clear_cancel() {}
void ModelSearchProblem::activate_state(const ModelSearchState&) {}

struct TabularModelSearchProblem::Impl {
  std::map<std::string, ModelSearchState> states;
  std::map<std::string, ModelSearchActions> actions;
  std::string initial_key;
  int evaluations = 0;
};

TabularModelSearchProblem::TabularModelSearchProblem()
    : impl_(new Impl) {}

void TabularModelSearchProblem::add_state(const std::string& key,
                                           double reward,
                                           bool acceptable) {
  require_key(key, "state");
  if (!std::isfinite(reward)) {
    throw ModelSearchConfigurationError("state reward must be finite");
  }
  impl_->states[key] = ModelSearchState(key, reward, acceptable);
}

void TabularModelSearchProblem::set_initial_state(const std::string& key) {
  require_key(key, "initial state");
  impl_->initial_key = key;
}

void TabularModelSearchProblem::add_action(
    const std::string& parent_key, const std::string& action_key,
    const std::string& result_key, double prior, bool terminal) {
  require_key(parent_key, "parent state");
  require_key(action_key, "action");
  require_key(result_key, "result state");
  if (!std::isfinite(prior) || prior < 0.0) {
    throw ModelSearchConfigurationError(
        "action prior must be finite and non-negative");
  }
  ModelSearchActions& actions = impl_->actions[parent_key];
  for (std::size_t i = 0; i < actions.size(); ++i) {
    if (actions[i].get_key() == action_key) {
      throw ModelSearchConfigurationError(
          "action keys must be unique within a state");
    }
  }
  actions.push_back(
      ModelSearchAction(action_key, result_key, prior, terminal));
}

ModelSearchState TabularModelSearchProblem::get_initial_state() {
  const std::map<std::string, ModelSearchState>::const_iterator found =
      impl_->states.find(impl_->initial_key);
  if (found == impl_->states.end()) {
    throw ModelSearchConfigurationError(
        "initial state has not been added to the graph problem");
  }
  return found->second;
}

ModelSearchActions TabularModelSearchProblem::get_actions(
    const ModelSearchState& state) {
  const std::map<std::string, ModelSearchActions>::const_iterator
      found = impl_->actions.find(state.get_key());
  if (found == impl_->actions.end()) return ModelSearchActions();
  return found->second;
}

ModelSearchState TabularModelSearchProblem::evaluate(
    const ModelSearchState&, const ModelSearchAction& action) {
  const std::map<std::string, ModelSearchState>::const_iterator found =
      impl_->states.find(action.get_predicted_state_key());
  if (found == impl_->states.end()) {
    throw ModelSearchConfigurationError(
        "action result state has not been added to the graph problem");
  }
  ++impl_->evaluations;
  return found->second;
}

int TabularModelSearchProblem::get_number_of_evaluations() const {
  return impl_->evaluations;
}

void TabularModelSearchProblem::reset_number_of_evaluations() {
  impl_->evaluations = 0;
}

namespace {

class FitSearchCancelObserver : public FitMinimizerObserver {
 public:
  explicit FitSearchCancelObserver(std::atomic<bool>* cancelled)
      : FitMinimizerObserver("FitSearchCancelObserver%1%"),
        cancelled_(cancelled) {}
  bool report(int, int, double) override { return !cancelled_->load(); }
  IMP_OBJECT_METHODS(FitSearchCancelObserver);

 private:
  std::atomic<bool>* cancelled_;
};

struct FitSearchSnapshot {
  std::vector<double> values;
  std::vector<int> fixed;
};

struct FitSearchGroup {
  std::vector<std::shared_ptr<GraphPort> > ports;
  std::vector<double> seeds;
};

struct FitSearchStructure {
  std::set<std::string> free_groups;
};

}  // namespace

struct FittingModelSearchProblem::Impl {
  std::shared_ptr<GraphNode> objective;
  std::string residual_key = "residuals";
  std::string score_output;
  std::string acceptable_output;
  double complexity_penalty = 0.0;
  std::map<std::string, FitSearchGroup> groups;
  std::vector<std::string> group_order;
  std::vector<std::shared_ptr<GraphPort> > ports;
  std::map<std::string, FitSearchStructure> structures;
  std::map<std::string, ModelSearchActions> actions;
  std::string initial_structure;
  std::map<std::string, FitSearchSnapshot> snapshots;
  std::map<std::string, std::string> snapshot_structures;
  unsigned long long next_snapshot = 1;
  std::atomic<bool> cancelled;
  int last_status = 0;
  std::string last_failure;

  Impl() : cancelled(false) {}

  FitSearchSnapshot capture() const {
    FitSearchSnapshot snapshot;
    snapshot.values.reserve(ports.size());
    snapshot.fixed.reserve(ports.size());
    for (std::size_t i = 0; i < ports.size(); ++i) {
      snapshot.values.push_back(ports[i]->get_value());
      snapshot.fixed.push_back(ports[i]->get_fixed() ? 1 : 0);
    }
    return snapshot;
  }

  void restore(const FitSearchSnapshot& snapshot) {
    if (snapshot.values.size() != ports.size() ||
        snapshot.fixed.size() != ports.size()) {
      throw ModelSearchConfigurationError("invalid cached fit snapshot");
    }
    for (std::size_t i = 0; i < ports.size(); ++i) {
      ports[i]->set_fixed(false);
      ports[i]->set_value(snapshot.values[i]);
      ports[i]->set_fixed(snapshot.fixed[i] != 0);
    }
    if (objective) objective->update();
  }

  int free_port_count(const std::string& structure_key) const {
    const std::map<std::string, FitSearchStructure>::const_iterator structure =
        structures.find(structure_key);
    if (structure == structures.end()) {
      throw ModelSearchConfigurationError("unknown fit structure '" +
                                             structure_key + "'");
    }
    int count = 0;
    for (std::set<std::string>::const_iterator key =
             structure->second.free_groups.begin();
         key != structure->second.free_groups.end(); ++key) {
      count += static_cast<int>(groups.find(*key)->second.ports.size());
    }
    return count;
  }

  std::vector<std::shared_ptr<GraphPort> > apply_structure(
      const std::string& parent_structure,
      const std::string& target_structure) {
    const FitSearchStructure& parent = structures.find(parent_structure)->second;
    const FitSearchStructure& target = structures.find(target_structure)->second;
    std::vector<std::shared_ptr<GraphPort> > free_ports;
    for (std::size_t gi = 0; gi < group_order.size(); ++gi) {
      const std::string& key = group_order[gi];
      FitSearchGroup& group = groups.find(key)->second;
      const bool was_free = parent.free_groups.count(key) != 0;
      const bool is_free = target.free_groups.count(key) != 0;
      if (is_free && !was_free && !group.seeds.empty()) {
        for (std::size_t pi = 0; pi < group.ports.size(); ++pi) {
          group.ports[pi]->set_fixed(false);
          group.ports[pi]->set_value(group.seeds[pi]);
        }
      }
      for (std::size_t pi = 0; pi < group.ports.size(); ++pi) {
        group.ports[pi]->set_fixed(!is_free);
        if (is_free) free_ports.push_back(group.ports[pi]);
      }
    }
    return free_ports;
  }

  std::pair<double, bool> score_current(int n_free) {
    objective->update();
    double reward = 0.0;
    if (!score_output.empty()) {
      const std::shared_ptr<GraphPort> score =
          objective->get_output_port(score_output);
      if (!score) {
        throw ModelSearchConfigurationError(
            "objective has no score output '" + score_output + "'");
      }
      reward = score->get_value();
    } else {
      const std::shared_ptr<GraphPort> residual =
          objective->get_output_port(residual_key);
      if (!residual) {
        throw ModelSearchConfigurationError(
            "objective has no residual output '" + residual_key + "'");
      }
      const std::vector<double>& values = residual->get_values_ref();
      double chi2 = 0.0;
      for (std::size_t i = 0; i < values.size(); ++i) {
        chi2 += values[i] * values[i];
      }
      reward = -0.5 * chi2 - complexity_penalty * n_free;
    }
    if (!std::isfinite(reward)) {
      throw ModelSearchConfigurationError("fit reward is not finite");
    }
    bool acceptable = false;
    if (!acceptable_output.empty()) {
      const std::shared_ptr<GraphPort> port =
          objective->get_output_port(acceptable_output);
      if (!port) {
        throw ModelSearchConfigurationError(
            "objective has no acceptable output '" + acceptable_output + "'");
      }
      acceptable = port->get_value_bool();
    }
    return std::make_pair(reward, acceptable);
  }
};

FittingModelSearchProblem::FittingModelSearchProblem() : impl_(new Impl) {}

FittingModelSearchProblem::FittingModelSearchProblem(
    std::shared_ptr<GraphNode> objective, const std::string& residual_key)
    : impl_(new Impl) {
  set_objective(std::move(objective), residual_key);
}

FittingModelSearchProblem::~FittingModelSearchProblem() {}

void FittingModelSearchProblem::set_objective(
    std::shared_ptr<GraphNode> objective, const std::string& residual_key) {
  if (!objective) {
    throw ModelSearchConfigurationError("fit problem objective is null");
  }
  require_key(residual_key, "residual output");
  impl_->objective = std::move(objective);
  impl_->residual_key = residual_key;
  impl_->snapshots.clear();
  impl_->snapshot_structures.clear();
}

std::shared_ptr<GraphNode> FittingModelSearchProblem::get_objective() const {
  return impl_->objective;
}

void FittingModelSearchProblem::add_parameter_group(
    const std::string& key,
    const std::vector<std::shared_ptr<GraphPort> >& ports,
    const std::vector<double>& enable_values) {
  require_key(key, "parameter group");
  if (impl_->groups.count(key)) {
    throw ModelSearchConfigurationError("duplicate parameter group '" + key +
                                           "'");
  }
  if (ports.empty()) {
    throw ModelSearchConfigurationError("parameter group is empty");
  }
  if (!enable_values.empty() && enable_values.size() != ports.size()) {
    throw ModelSearchConfigurationError(
        "enable values must match the parameter group size");
  }
  std::set<GraphPort*> existing;
  for (std::size_t i = 0; i < impl_->ports.size(); ++i) {
    existing.insert(impl_->ports[i].get());
  }
  for (std::size_t i = 0; i < ports.size(); ++i) {
    if (!ports[i]) {
      throw ModelSearchConfigurationError("parameter group has a null port");
    }
    if (ports[i]->get_is_vector()) {
      throw ModelSearchConfigurationError(
          "fit model search supports scalar parameter ports only");
    }
    if (ports[i]->get_link()) {
      throw ModelSearchConfigurationError(
          "fit model search requires canonical owner ports, not linked followers");
    }
    if (!existing.insert(ports[i].get()).second) {
      throw ModelSearchConfigurationError(
          "a parameter port may belong to only one group");
    }
  }
  FitSearchGroup group;
  group.ports = ports;
  group.seeds = enable_values;
  impl_->groups[key] = group;
  impl_->group_order.push_back(key);
  impl_->ports.insert(impl_->ports.end(), ports.begin(), ports.end());
}

std::vector<std::string>
FittingModelSearchProblem::get_parameter_group_keys() const {
  return impl_->group_order;
}

void FittingModelSearchProblem::add_structure(
    const std::string& key, const std::vector<std::string>& free_groups) {
  require_key(key, "fit structure");
  FitSearchStructure structure;
  for (std::size_t i = 0; i < free_groups.size(); ++i) {
    if (!impl_->groups.count(free_groups[i])) {
      throw ModelSearchConfigurationError("unknown parameter group '" +
                                             free_groups[i] + "'");
    }
    structure.free_groups.insert(free_groups[i]);
  }
  impl_->structures[key] = structure;
}

void FittingModelSearchProblem::set_initial_structure(const std::string& key) {
  require_key(key, "initial structure");
  impl_->initial_structure = key;
}

void FittingModelSearchProblem::add_action(
    const std::string& parent_structure, const std::string& action_key,
    const std::string& result_structure, double prior, bool terminal) {
  require_key(parent_structure, "parent structure");
  require_key(action_key, "action");
  require_key(result_structure, "result structure");
  if (!std::isfinite(prior) || prior < 0.0) {
    throw ModelSearchConfigurationError(
        "action prior must be finite and non-negative");
  }
  ModelSearchActions& actions =
      impl_->actions[parent_structure];
  for (std::size_t i = 0; i < actions.size(); ++i) {
    if (actions[i].get_key() == action_key) {
      throw ModelSearchConfigurationError(
          "action keys must be unique within a structure");
    }
  }
  actions.push_back(ModelSearchAction(action_key, result_structure, prior,
                                         terminal));
}

void FittingModelSearchProblem::set_score_output(const std::string& key) {
  require_key(key, "score output");
  impl_->score_output = key;
}
void FittingModelSearchProblem::clear_score_output() {
  impl_->score_output.clear();
}
const std::string& FittingModelSearchProblem::get_score_output() const {
  return impl_->score_output;
}
void FittingModelSearchProblem::set_complexity_penalty(double value) {
  if (!std::isfinite(value) || value < 0.0) {
    throw ModelSearchConfigurationError(
        "complexity penalty must be finite and non-negative");
  }
  impl_->complexity_penalty = value;
}
double FittingModelSearchProblem::get_complexity_penalty() const {
  return impl_->complexity_penalty;
}
void FittingModelSearchProblem::set_acceptable_output(const std::string& key) {
  require_key(key, "acceptable output");
  impl_->acceptable_output = key;
}
void FittingModelSearchProblem::clear_acceptable_output() {
  impl_->acceptable_output.clear();
}

ModelSearchState FittingModelSearchProblem::get_initial_state() {
  if (!impl_->objective) {
    throw ModelSearchConfigurationError("fit problem has no objective");
  }
  const std::map<std::string, FitSearchStructure>::const_iterator structure =
      impl_->structures.find(impl_->initial_structure);
  if (structure == impl_->structures.end()) {
    throw ModelSearchConfigurationError(
        "initial fit structure has not been declared");
  }
  for (std::size_t gi = 0; gi < impl_->group_order.size(); ++gi) {
    const std::string& key = impl_->group_order[gi];
    const bool should_be_free = structure->second.free_groups.count(key) != 0;
    const FitSearchGroup& group = impl_->groups.find(key)->second;
    for (std::size_t pi = 0; pi < group.ports.size(); ++pi) {
      if (group.ports[pi]->get_fixed() == should_be_free) {
        throw ModelSearchConfigurationError(
            "initial structure does not match the live port fixed mask");
      }
    }
  }
  const FitSearchSnapshot root = impl_->capture();
  try {
    const std::pair<double, bool> score = impl_->score_current(
        impl_->free_port_count(impl_->initial_structure));
    impl_->snapshots[impl_->initial_structure] = root;
    impl_->snapshot_structures[impl_->initial_structure] =
        impl_->initial_structure;
    return ModelSearchState(impl_->initial_structure,
                               impl_->initial_structure, score.first,
                               score.second);
  } catch (...) {
    impl_->restore(root);
    throw;
  }
}

ModelSearchActions FittingModelSearchProblem::get_actions(
    const ModelSearchState& state) {
  const std::map<std::string, ModelSearchActions>::const_iterator
      found = impl_->actions.find(state.get_structure_key());
  if (found == impl_->actions.end()) return ModelSearchActions();
  return found->second;
}

ModelSearchState FittingModelSearchProblem::evaluate(
    const ModelSearchState& parent, const ModelSearchAction& action) {
  impl_->last_status = 0;
  impl_->last_failure.clear();
  const std::map<std::string, FitSearchSnapshot>::const_iterator parent_snapshot =
      impl_->snapshots.find(parent.get_key());
  if (parent_snapshot == impl_->snapshots.end()) {
    throw ModelSearchConfigurationError("parent fit snapshot is not cached");
  }
  const std::string target = action.get_predicted_state_key();
  if (!impl_->structures.count(target)) {
    throw ModelSearchConfigurationError("unknown result structure '" + target +
                                           "'");
  }
  try {
    impl_->restore(parent_snapshot->second);
    std::vector<std::shared_ptr<GraphPort> > free_ports =
        impl_->apply_structure(parent.get_structure_key(), target);
    if (impl_->cancelled.load()) {
      impl_->restore(parent_snapshot->second);
      impl_->last_status = -1;
      impl_->last_failure = "cancelled before minimization";
      return parent;
    }
    if (!free_ports.empty()) {
      FitMinimizer minimizer;
      minimizer.set_parameter_ports(free_ports);
      minimizer.set_objective(impl_->objective, impl_->residual_key);
      IMP::Pointer<FitSearchCancelObserver> observer(
          new FitSearchCancelObserver(&impl_->cancelled));
      minimizer.set_observer(observer.get());
      impl_->last_status = minimizer.run();
      if (minimizer.get_cancelled() || impl_->cancelled.load()) {
        impl_->restore(parent_snapshot->second);
        impl_->last_status = -1;
        impl_->last_failure = "minimization cancelled";
        return parent;
      }
      if (impl_->last_status < 1 || impl_->last_status > 4) {
        std::ostringstream message;
        message << "minimizer did not converge (status " << impl_->last_status
                << ")";
        impl_->last_failure = message.str();
        impl_->restore(parent_snapshot->second);
        return parent;
      }
    } else {
      impl_->objective->update();
      impl_->last_status = 1;
    }
    const std::pair<double, bool> score =
        impl_->score_current(static_cast<int>(free_ports.size()));
    std::ostringstream state_key;
    state_key << target << "@" << impl_->next_snapshot++;
    const std::string key = state_key.str();
    impl_->snapshots[key] = impl_->capture();
    impl_->snapshot_structures[key] = target;
    return ModelSearchState(key, target, score.first, score.second);
  } catch (const std::exception& error) {
    impl_->last_failure = error.what();
    impl_->last_status = 0;
    impl_->restore(parent_snapshot->second);
    return parent;
  }
}

void FittingModelSearchProblem::request_cancel() {
  impl_->cancelled.store(true);
}
void FittingModelSearchProblem::clear_cancel() {
  impl_->cancelled.store(false);
}
void FittingModelSearchProblem::activate_state(
    const ModelSearchState& state) {
  restore_state(state.get_key());
}
bool FittingModelSearchProblem::has_cached_state(
    const std::string& state_key) const {
  return impl_->snapshots.count(state_key) != 0;
}
std::vector<double> FittingModelSearchProblem::get_cached_values(
    const std::string& state_key) const {
  const std::map<std::string, FitSearchSnapshot>::const_iterator found =
      impl_->snapshots.find(state_key);
  if (found == impl_->snapshots.end()) return std::vector<double>();
  return found->second.values;
}
std::vector<int> FittingModelSearchProblem::get_cached_fixed(
    const std::string& state_key) const {
  const std::map<std::string, FitSearchSnapshot>::const_iterator found =
      impl_->snapshots.find(state_key);
  if (found == impl_->snapshots.end()) return std::vector<int>();
  return found->second.fixed;
}
void FittingModelSearchProblem::restore_state(const std::string& state_key) {
  const std::map<std::string, FitSearchSnapshot>::const_iterator found =
      impl_->snapshots.find(state_key);
  if (found == impl_->snapshots.end()) {
    throw ModelSearchConfigurationError("unknown cached fit state '" +
                                           state_key + "'");
  }
  impl_->restore(found->second);
}
int FittingModelSearchProblem::get_last_fit_status() const {
  return impl_->last_status;
}
const std::string& FittingModelSearchProblem::get_last_failure() const {
  return impl_->last_failure;
}

ModelSearchConfig::ModelSearchConfig()
    : n_simulations_(200), c_puct_(1.5), reward_scale_(4.0),
      dirichlet_alpha_(0.15), dirichlet_fraction_(0.25), seed_(42u) {}

void ModelSearchConfig::set_number_of_simulations(int value) {
  n_simulations_ = value;
}
int ModelSearchConfig::get_number_of_simulations() const {
  return n_simulations_;
}
void ModelSearchConfig::set_c_puct(double value) { c_puct_ = value; }
double ModelSearchConfig::get_c_puct() const { return c_puct_; }
void ModelSearchConfig::set_reward_scale(double value) {
  reward_scale_ = value;
}
double ModelSearchConfig::get_reward_scale() const { return reward_scale_; }
void ModelSearchConfig::set_dirichlet_alpha(double value) {
  dirichlet_alpha_ = value;
}
double ModelSearchConfig::get_dirichlet_alpha() const {
  return dirichlet_alpha_;
}
void ModelSearchConfig::set_dirichlet_fraction(double value) {
  dirichlet_fraction_ = value;
}
double ModelSearchConfig::get_dirichlet_fraction() const {
  return dirichlet_fraction_;
}
void ModelSearchConfig::set_seed(unsigned int value) { seed_ = value; }
unsigned int ModelSearchConfig::get_seed() const { return seed_; }

ModelSearchResult::ModelSearchResult()
    : n_states_evaluated_(0), n_simulations_(0), improvement_(0.0),
      acceptable_(false), cancelled_(false) {}

const ModelSearchState& ModelSearchResult::get_root_state() const {
  return root_state_;
}
const ModelSearchState& ModelSearchResult::get_best_state() const {
  return best_state_;
}
const std::vector<std::string>& ModelSearchResult::get_best_path() const {
  return best_path_;
}
int ModelSearchResult::get_number_of_states_evaluated() const {
  return n_states_evaluated_;
}
int ModelSearchResult::get_number_of_simulations() const {
  return n_simulations_;
}
double ModelSearchResult::get_improvement() const { return improvement_; }
bool ModelSearchResult::get_acceptable() const { return acceptable_; }
bool ModelSearchResult::get_cancelled() const { return cancelled_; }

struct ModelSearch::Impl {
  std::shared_ptr<ModelSearchProblem> problem;
  ModelSearchConfig config;
  std::atomic<bool> cancel_requested;
  Impl() : cancel_requested(false) {}
};

ModelSearch::ModelSearch() : impl_(new Impl) {}
ModelSearch::ModelSearch(std::shared_ptr<ModelSearchProblem> problem)
    : impl_(new Impl) {
  impl_->problem = std::move(problem);
}
ModelSearch::~ModelSearch() {}

void ModelSearch::set_problem(
    std::shared_ptr<ModelSearchProblem> problem) {
  impl_->problem = std::move(problem);
}
std::shared_ptr<ModelSearchProblem> ModelSearch::get_problem() const {
  return impl_->problem;
}
void ModelSearch::set_config(const ModelSearchConfig& config) {
  impl_->config = config;
}
ModelSearchConfig ModelSearch::get_config() const { return impl_->config; }
void ModelSearch::request_cancel() {
  impl_->cancel_requested.store(true);
  if (impl_->problem) impl_->problem->request_cancel();
}
void ModelSearch::clear_cancel() {
  impl_->cancel_requested.store(false);
  if (impl_->problem) impl_->problem->clear_cancel();
}
bool ModelSearch::get_cancel_requested() const {
  return impl_->cancel_requested.load();
}

ModelSearchResult ModelSearch::run() {
  if (!impl_->problem) {
    throw ModelSearchConfigurationError("model search has no problem");
  }
  const ModelSearchConfig cfg = impl_->config;
  if (cfg.get_number_of_simulations() < 0 ||
      !std::isfinite(cfg.get_c_puct()) || cfg.get_c_puct() < 0.0 ||
      !std::isfinite(cfg.get_reward_scale()) ||
      cfg.get_reward_scale() <= 0.0 ||
      !std::isfinite(cfg.get_dirichlet_fraction()) ||
      cfg.get_dirichlet_fraction() < 0.0 ||
      cfg.get_dirichlet_fraction() > 1.0 ||
      (cfg.get_dirichlet_fraction() > 0.0 &&
       (!std::isfinite(cfg.get_dirichlet_alpha()) ||
        cfg.get_dirichlet_alpha() <= 0.0))) {
    throw ModelSearchConfigurationError("invalid model-search controls");
  }

  const ModelSearchState initial = impl_->problem->get_initial_state();
  require_state(initial, "initial state");
  const double root_reward = initial.get_reward();
  std::mt19937_64 rng(cfg.get_seed());
  int evaluator_calls = 0;

  std::unique_ptr<TreeNode> root(new TreeNode);
  root->state = initial;
  root->structure_key = initial.get_structure_key();
  root->evaluated = true;
  root->visits = 1;

  const auto expand = [&](TreeNode* node) {
    ModelSearchActions actions =
        impl_->problem->get_actions(node->state);
    std::set<std::string> action_keys;
    double prior_sum = 0.0;
    for (std::size_t i = 0; i < actions.size(); ++i) {
      require_key(actions[i].get_key(), "action");
      require_key(actions[i].get_predicted_state_key(), "predicted state");
      if (!action_keys.insert(actions[i].get_key()).second) {
        throw ModelSearchConfigurationError(
            "action keys must be unique within a state");
      }
      if (!std::isfinite(actions[i].get_prior()) ||
          actions[i].get_prior() < 0.0) {
        throw ModelSearchConfigurationError(
            "action prior must be finite and non-negative");
      }
      prior_sum += actions[i].get_prior();
    }
    const bool uniform = !actions.empty() && prior_sum == 0.0;

    std::vector<double> priors(actions.size(), 0.0);
    for (std::size_t i = 0; i < actions.size(); ++i) {
      priors[i] = uniform ? 1.0 / actions.size()
                          : actions[i].get_prior() / prior_sum;
    }
    if (node == root.get() && actions.size() > 1 &&
        cfg.get_dirichlet_fraction() > 0.0) {
      std::gamma_distribution<double> gamma(cfg.get_dirichlet_alpha(), 1.0);
      std::vector<double> noise(actions.size());
      double noise_sum = 0.0;
      for (std::size_t i = 0; i < noise.size(); ++i) {
        noise[i] = gamma(rng);
        noise_sum += noise[i];
      }
      if (noise_sum > 0.0) {
        for (std::size_t i = 0; i < priors.size(); ++i) {
          priors[i] = (1.0 - cfg.get_dirichlet_fraction()) * priors[i] +
                      cfg.get_dirichlet_fraction() * noise[i] / noise_sum;
        }
      }
    }

    std::set<std::string> ancestor_keys;
    for (TreeNode* ancestor = node->parent; ancestor;
         ancestor = ancestor->parent) {
      ancestor_keys.insert(ancestor->structure_key);
    }
    for (std::size_t i = 0; i < actions.size(); ++i) {
      if (ancestor_keys.count(actions[i].get_predicted_state_key())) continue;
      std::unique_ptr<TreeNode> child(new TreeNode(priors[i], node));
      child->parent_action = actions[i];
      child->terminal = actions[i].get_terminal();
      child->structure_key = actions[i].get_predicted_state_key();
      node->children.push_back(std::move(child));
    }
    node->expanded = true;
  };

  expand(root.get());
  int done = 0;
  while (done < cfg.get_number_of_simulations() &&
         !impl_->cancel_requested.load()) {
    std::vector<TreeNode*> path(1, root.get());
    TreeNode* node = root.get();
    while (node->expanded && !node->children.empty()) {
      TreeNode* best = node->children[0].get();
      double best_score = puct(*best, node->visits, cfg.get_c_puct());
      for (std::size_t i = 1; i < node->children.size(); ++i) {
        TreeNode* candidate = node->children[i].get();
        const double score = puct(*candidate, node->visits,
                                  cfg.get_c_puct());
        if (score > best_score) {
          best = candidate;
          best_score = score;
        }
      }
      node = best;
      path.push_back(node);
      if (node->terminal) break;
    }

    double value = 0.0;
    if (node->terminal) {
      node->state = node->parent->state;
      node->evaluated = true;
      value = reward_value(node->state, root_reward, cfg.get_reward_scale());
    } else if (!node->evaluated) {
      node->state = impl_->problem->evaluate(node->parent->state,
                                             node->parent_action);
      ++evaluator_calls;
      require_state(node->state, "evaluated state");
      node->structure_key = node->state.get_structure_key();
      node->evaluated = true;
      for (TreeNode* ancestor = node->parent; ancestor;
           ancestor = ancestor->parent) {
        if (ancestor->structure_key == node->structure_key) {
          node->terminal = true;
          node->expanded = true;
          break;
        }
      }
      value = reward_value(node->state, root_reward, cfg.get_reward_scale());
    } else if (node->expanded && node->children.empty()) {
      value = reward_value(node->state, root_reward, cfg.get_reward_scale());
    } else {
      expand(node);
      value = reward_value(node->state, root_reward, cfg.get_reward_scale());
    }

    for (std::size_t i = 0; i < path.size(); ++i) {
      ++path[i]->visits;
      path[i]->value_sum += value;
    }
    ++done;
  }

  TreeNode* best = root.get();
  std::vector<TreeNode*> stack(1, root.get());
  while (!stack.empty()) {
    TreeNode* node = stack.back();
    stack.pop_back();
    if (node->evaluated &&
        node->state.get_reward() > best->state.get_reward()) {
      best = node;
    }
    for (std::size_t i = 0; i < node->children.size(); ++i) {
      stack.push_back(node->children[i].get());
    }
  }

  ModelSearchResult result;
  result.root_state_ = initial;
  result.best_state_ = best->state;
  result.n_states_evaluated_ = 1 + evaluator_calls;
  result.n_simulations_ = done;
  result.improvement_ = best->state.get_reward() - root_reward;
  result.acceptable_ = best->state.get_acceptable();
  result.cancelled_ = done < cfg.get_number_of_simulations() &&
                      impl_->cancel_requested.load();

  TreeNode* path_node = root.get();
  while (path_node->expanded && !path_node->children.empty()) {
    TreeNode* most_visited = path_node->children[0].get();
    for (std::size_t i = 1; i < path_node->children.size(); ++i) {
      if (path_node->children[i]->visits > most_visited->visits) {
        most_visited = path_node->children[i].get();
      }
    }
    if (most_visited->terminal || most_visited->visits == 0) break;
    result.best_path_.push_back(most_visited->parent_action.get_key());
    if (!most_visited->evaluated) break;
    path_node = most_visited;
  }
  impl_->problem->activate_state(result.cancelled_ ? result.root_state_
                                                  : result.best_state_);
  return result;
}

IMPBFF_END_NAMESPACE
