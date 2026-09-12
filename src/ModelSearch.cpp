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
    std::vector<std::shared_ptr<GraphPort> > free_ports;
    for (std::size_t gi = 0; gi < impl_->group_order.size(); ++gi) {
      const FitSearchGroup& group = impl_->groups.find(impl_->group_order[gi])->second;
      for (std::size_t pi = 0; pi < group.ports.size(); ++pi) {
        if (!group.ports[pi]->get_fixed()) free_ports.push_back(group.ports[pi]);
      }
    }
    if (!free_ports.empty() && !impl_->cancelled.load()) {
      FitMinimizer minimizer;
      minimizer.set_parameter_ports(free_ports);
      minimizer.set_objective(impl_->objective, impl_->residual_key);
      impl_->last_status = minimizer.run();
      if (impl_->last_status < 1 || impl_->last_status > 4) {
        throw ModelSearchConfigurationError(
            "initial fit structure did not converge");
      }
    } else {
      impl_->objective->update();
      impl_->last_status = impl_->cancelled.load() ? -1 : 1;
    }
    const std::pair<double, bool> score = impl_->score_current(
        static_cast<int>(free_ports.size()));
    impl_->snapshots[impl_->initial_structure] = impl_->capture();
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

namespace {

struct MultiStructureRecord {
  std::shared_ptr<GraphNode> objective;
  std::vector<std::shared_ptr<GraphNode> > graph_nodes;
  std::string residual_key;
  std::string score_output;
  std::string acceptable_output;
  std::vector<double> initial_values;
  std::vector<int> fixed;
  double effective_sample_size = 0.0;
  double complexity = 0.0;
  bool use_bic = false;
};

}  // namespace

struct MultiStructureModelSearchProblem::Impl {
  std::map<std::string, std::shared_ptr<GraphPort> > parameters;
  std::vector<std::string> parameter_order;
  std::map<std::string, MultiStructureRecord> structures;
  std::vector<std::string> structure_order;
  std::map<std::string, ModelSearchActions> actions;
  std::string initial_structure;
  std::string active_structure;
  std::map<std::string, FitSearchSnapshot> snapshots;
  std::map<std::string, std::string> snapshot_structures;
  unsigned long long next_snapshot = 1;
  std::atomic<bool> cancelled;
  int last_status = 0;
  std::string last_failure;

  Impl() : cancelled(false) {}

  void validate_registry() const {
    std::set<GraphPort*> owners;
    if (parameter_order.size() != parameters.size()) {
      throw ModelSearchConfigurationError(
          "canonical parameter registry order is inconsistent");
    }
    for (std::size_t i = 0; i < parameter_order.size(); ++i) {
      const std::map<std::string, std::shared_ptr<GraphPort> >::const_iterator
          found = parameters.find(parameter_order[i]);
      if (found == parameters.end() || !found->second) {
        throw ModelSearchConfigurationError(
            "canonical parameter registry contains a missing owner");
      }
      if (found->second->get_is_vector() || found->second->get_link()) {
        throw ModelSearchConfigurationError(
            "canonical parameters must remain scalar unlinked owner ports");
      }
      if (!owners.insert(found->second.get()).second) {
        throw ModelSearchConfigurationError(
            "canonical parameter registry contains a duplicate owner");
      }
    }
  }

  const MultiStructureRecord& structure(const std::string& key) const {
    const std::map<std::string, MultiStructureRecord>::const_iterator found =
        structures.find(key);
    if (found == structures.end()) {
      throw ModelSearchConfigurationError("unknown model structure '" + key +
                                           "'");
    }
    return found->second;
  }

  MultiStructureRecord& structure(const std::string& key) {
    const std::map<std::string, MultiStructureRecord>::iterator found =
        structures.find(key);
    if (found == structures.end()) {
      throw ModelSearchConfigurationError("unknown model structure '" + key +
                                           "'");
    }
    return found->second;
  }

  FitSearchSnapshot capture() const {
    FitSearchSnapshot snapshot;
    snapshot.values.reserve(parameter_order.size());
    snapshot.fixed.reserve(parameter_order.size());
    for (std::size_t i = 0; i < parameter_order.size(); ++i) {
      const std::shared_ptr<GraphPort>& port =
          parameters.find(parameter_order[i])->second;
      snapshot.values.push_back(port->get_value());
      snapshot.fixed.push_back(port->get_fixed() ? 1 : 0);
    }
    return snapshot;
  }

  void restore_values(const FitSearchSnapshot& snapshot) {
    if (snapshot.values.size() != parameter_order.size() ||
        snapshot.fixed.size() != parameter_order.size()) {
      throw ModelSearchConfigurationError("invalid cached model snapshot");
    }
    for (std::size_t i = 0; i < parameter_order.size(); ++i) {
      const std::shared_ptr<GraphPort>& port =
          parameters.find(parameter_order[i])->second;
      port->set_fixed(false);
      port->set_value(snapshot.values[i]);
      port->set_fixed(snapshot.fixed[i] != 0);
    }
  }

  void select_and_update(const std::string& structure_key) {
    MultiStructureRecord& selected = structure(structure_key);
    active_structure = structure_key;
    selected.objective->update();
  }

  void restore(const FitSearchSnapshot& snapshot,
               const std::string& structure_key) {
    restore_values(snapshot);
    select_and_update(structure_key);
  }

  void apply_initial(const MultiStructureRecord& target) {
    for (std::size_t i = 0; i < parameter_order.size(); ++i) {
      const std::shared_ptr<GraphPort>& port =
          parameters.find(parameter_order[i])->second;
      port->set_fixed(false);
      port->set_value(target.initial_values[i]);
      port->set_fixed(target.fixed[i] != 0);
    }
  }

  std::vector<std::shared_ptr<GraphPort> > apply_transition(
      const MultiStructureRecord& parent,
      const MultiStructureRecord& target) {
    std::vector<std::shared_ptr<GraphPort> > free_ports;
    for (std::size_t i = 0; i < parameter_order.size(); ++i) {
      const std::shared_ptr<GraphPort>& port =
          parameters.find(parameter_order[i])->second;
      const bool was_free = parent.fixed[i] == 0;
      const bool is_free = target.fixed[i] == 0;
      port->set_fixed(false);
      if (!is_free || !was_free) {
        port->set_value(target.initial_values[i]);
      }
      port->set_fixed(!is_free);
      if (is_free) free_ports.push_back(port);
    }
    return free_ports;
  }

  std::pair<double, bool> score_current(
      const MultiStructureRecord& selected) {
    selected.objective->update();
    double reward = 0.0;
    if (!selected.score_output.empty()) {
      const std::shared_ptr<GraphPort> score =
          selected.objective->get_output_port(selected.score_output);
      if (!score) {
        throw ModelSearchConfigurationError(
            "objective has no score output '" + selected.score_output + "'");
      }
      reward = score->get_value();
    } else {
      if (selected.residual_key.empty()) {
        throw ModelSearchConfigurationError(
            "structure has neither a residual nor score output");
      }
      const std::shared_ptr<GraphPort> residual =
          selected.objective->get_output_port(selected.residual_key);
      if (!residual) {
        throw ModelSearchConfigurationError(
            "objective has no residual output '" + selected.residual_key +
            "'");
      }
      const std::vector<double>& values = residual->get_values_ref();
      double chi2 = 0.0;
      for (std::size_t i = 0; i < values.size(); ++i) {
        chi2 += values[i] * values[i];
      }
      reward = -0.5 * chi2;
      if (selected.use_bic) {
        reward -= 0.5 * selected.complexity *
                  std::log(selected.effective_sample_size);
      }
    }
    if (!std::isfinite(reward)) {
      throw ModelSearchConfigurationError("model-search reward is not finite");
    }
    bool acceptable = false;
    if (!selected.acceptable_output.empty()) {
      const std::shared_ptr<GraphPort> port =
          selected.objective->get_output_port(selected.acceptable_output);
      if (!port) {
        throw ModelSearchConfigurationError(
            "objective has no acceptable output '" +
            selected.acceptable_output + "'");
      }
      acceptable = port->get_value_bool();
    }
    return std::make_pair(reward, acceptable);
  }

  void rollback(const FitSearchSnapshot& snapshot,
                const std::string& structure_key) {
    restore_values(snapshot);
    active_structure = structure_key;
    try {
      structure(structure_key).objective->update();
    } catch (...) {
      // Preserve the original transition failure. Values and masks are the
      // transaction boundary; a subsequent explicit restore will re-evaluate.
    }
  }
};

MultiStructureModelSearchProblem::MultiStructureModelSearchProblem()
    : impl_(new Impl) {}

MultiStructureModelSearchProblem::~MultiStructureModelSearchProblem() {}

void MultiStructureModelSearchProblem::add_parameter(
    const std::string& canonical_id, std::shared_ptr<GraphPort> owner) {
  require_key(canonical_id, "canonical parameter");
  if (!impl_->structures.empty()) {
    throw ModelSearchConfigurationError(
        "canonical parameters must be registered before structures");
  }
  if (impl_->parameters.count(canonical_id)) {
    throw ModelSearchConfigurationError("duplicate canonical parameter id '" +
                                         canonical_id + "'");
  }
  if (!owner) {
    throw ModelSearchConfigurationError("canonical parameter is null");
  }
  if (owner->get_is_vector()) {
    throw ModelSearchConfigurationError(
        "model search supports scalar canonical parameters only");
  }
  if (owner->get_link()) {
    throw ModelSearchConfigurationError(
        "canonical parameters must be owner ports, not linked followers");
  }
  for (std::map<std::string, std::shared_ptr<GraphPort> >::const_iterator it =
           impl_->parameters.begin();
       it != impl_->parameters.end(); ++it) {
    if (it->second.get() == owner.get()) {
      throw ModelSearchConfigurationError(
          "a canonical parameter port cannot have more than one id");
    }
  }
  impl_->parameters[canonical_id] = owner;
  impl_->parameter_order.push_back(canonical_id);
}

std::vector<std::string>
MultiStructureModelSearchProblem::get_parameter_ids() const {
  return impl_->parameter_order;
}

std::shared_ptr<GraphPort> MultiStructureModelSearchProblem::get_parameter(
    const std::string& canonical_id) const {
  const std::map<std::string, std::shared_ptr<GraphPort> >::const_iterator found =
      impl_->parameters.find(canonical_id);
  if (found == impl_->parameters.end()) return std::shared_ptr<GraphPort>();
  return found->second;
}

void MultiStructureModelSearchProblem::add_structure(
    const std::string& key, std::shared_ptr<GraphNode> objective,
    const std::vector<std::string>& parameter_ids,
    const std::vector<std::shared_ptr<GraphPort> >& parameter_ports,
    const std::vector<double>& initial_values,
    const std::vector<int>& fixed_mask, const std::string& residual_key) {
  require_key(key, "model structure");
  if (impl_->structures.count(key)) {
    throw ModelSearchConfigurationError("duplicate model structure '" + key +
                                         "'");
  }
  if (!objective) {
    throw ModelSearchConfigurationError("model structure objective is null");
  }
  impl_->validate_registry();
  const std::size_t count = impl_->parameter_order.size();
  if (parameter_ids.size() != count || parameter_ports.size() != count ||
      initial_values.size() != count || fixed_mask.size() != count) {
    throw ModelSearchConfigurationError(
        "structure parameter state must cover the complete canonical registry");
  }
  std::set<std::string> seen_ids;
  std::set<GraphPort*> seen_ports;
  MultiStructureRecord record;
  record.objective = std::move(objective);
  record.residual_key = residual_key;
  record.initial_values.resize(count);
  record.fixed.resize(count);
  for (std::size_t i = 0; i < count; ++i) {
    require_key(parameter_ids[i], "structure parameter");
    if (!seen_ids.insert(parameter_ids[i]).second) {
      throw ModelSearchConfigurationError(
          "duplicate canonical parameter id in structure '" + key + "'");
    }
    if (!parameter_ports[i] ||
        !seen_ports.insert(parameter_ports[i].get()).second) {
      throw ModelSearchConfigurationError(
          "duplicate or null canonical parameter port in structure '" + key +
          "'");
    }
    const std::map<std::string, std::shared_ptr<GraphPort> >::const_iterator
        registered = impl_->parameters.find(parameter_ids[i]);
    if (registered == impl_->parameters.end()) {
      throw ModelSearchConfigurationError("unknown canonical parameter id '" +
                                           parameter_ids[i] + "'");
    }
    if (registered->second.get() != parameter_ports[i].get()) {
      throw ModelSearchConfigurationError(
          "structure parameter is not the registered canonical owner for '" +
          parameter_ids[i] + "'");
    }
    const std::vector<std::string>::const_iterator ordered = std::find(
        impl_->parameter_order.cbegin(), impl_->parameter_order.cend(),
        parameter_ids[i]);
    const std::size_t index = static_cast<std::size_t>(
        std::distance(impl_->parameter_order.cbegin(), ordered));
    if (!std::isfinite(initial_values[i])) {
      throw ModelSearchConfigurationError(
          "structure initial parameter values must be finite");
    }
    if (fixed_mask[i] != 0 && fixed_mask[i] != 1) {
      throw ModelSearchConfigurationError(
          "structure fixed mask values must be zero or one");
    }
    record.initial_values[index] = initial_values[i];
    record.fixed[index] = fixed_mask[i];
  }
  if (seen_ids.size() != impl_->parameters.size()) {
    throw ModelSearchConfigurationError(
        "structure is missing canonical parameter ids");
  }
  if (!residual_key.empty() &&
      !record.objective->get_output_port(residual_key)) {
    throw ModelSearchConfigurationError("objective has no residual output '" +
                                         residual_key + "'");
  }
  impl_->structures[key] = record;
  impl_->structure_order.push_back(key);
}

std::vector<std::string>
MultiStructureModelSearchProblem::get_structure_keys() const {
  return impl_->structure_order;
}

std::shared_ptr<GraphNode>
MultiStructureModelSearchProblem::get_structure_objective(
    const std::string& key) const {
  return impl_->structure(key).objective;
}

void MultiStructureModelSearchProblem::add_structure_node(
    const std::string& structure_key, std::shared_ptr<GraphNode> node) {
  if (!node) {
    throw ModelSearchConfigurationError("structure graph node is null");
  }
  MultiStructureRecord& selected = impl_->structure(structure_key);
  if (selected.objective.get() == node.get()) {
    throw ModelSearchConfigurationError(
        "the structure objective is retained automatically");
  }
  for (std::size_t i = 0; i < selected.graph_nodes.size(); ++i) {
    if (selected.graph_nodes[i].get() == node.get()) {
      throw ModelSearchConfigurationError(
          "duplicate node in structure graph");
    }
  }
  selected.graph_nodes.push_back(std::move(node));
}

void MultiStructureModelSearchProblem::set_initial_structure(
    const std::string& key) {
  require_key(key, "initial structure");
  if (!impl_->structures.count(key)) {
    throw ModelSearchConfigurationError("unknown initial structure '" + key +
                                         "'");
  }
  impl_->initial_structure = key;
}

void MultiStructureModelSearchProblem::add_action(
    const std::string& parent_structure, const std::string& action_key,
    const std::string& result_structure, double prior, bool terminal) {
  require_key(action_key, "action");
  impl_->structure(parent_structure);
  impl_->structure(result_structure);
  if (!std::isfinite(prior) || prior < 0.0) {
    throw ModelSearchConfigurationError(
        "action prior must be finite and non-negative");
  }
  ModelSearchActions& actions = impl_->actions[parent_structure];
  for (std::size_t i = 0; i < actions.size(); ++i) {
    if (actions[i].get_key() == action_key) {
      throw ModelSearchConfigurationError(
          "action keys must be unique within a structure");
    }
  }
  actions.push_back(ModelSearchAction(action_key, result_structure, prior,
                                      terminal));
}

void MultiStructureModelSearchProblem::set_structure_score_output(
    const std::string& structure_key, const std::string& output_key) {
  require_key(output_key, "score output");
  MultiStructureRecord& selected = impl_->structure(structure_key);
  if (!selected.objective->get_output_port(output_key)) {
    throw ModelSearchConfigurationError("objective has no score output '" +
                                         output_key + "'");
  }
  selected.score_output = output_key;
}

void MultiStructureModelSearchProblem::clear_structure_score_output(
    const std::string& structure_key) {
  impl_->structure(structure_key).score_output.clear();
}

void MultiStructureModelSearchProblem::set_structure_acceptable_output(
    const std::string& structure_key, const std::string& output_key) {
  require_key(output_key, "acceptable output");
  MultiStructureRecord& selected = impl_->structure(structure_key);
  if (!selected.objective->get_output_port(output_key)) {
    throw ModelSearchConfigurationError(
        "objective has no acceptable output '" + output_key + "'");
  }
  selected.acceptable_output = output_key;
}

void MultiStructureModelSearchProblem::clear_structure_acceptable_output(
    const std::string& structure_key) {
  impl_->structure(structure_key).acceptable_output.clear();
}

void MultiStructureModelSearchProblem::set_structure_bic_metadata(
    const std::string& structure_key, double effective_sample_size,
    double complexity) {
  if (!std::isfinite(effective_sample_size) || effective_sample_size <= 0.0) {
    throw ModelSearchConfigurationError(
        "effective sample size must be finite and positive");
  }
  if (!std::isfinite(complexity) || complexity < 0.0) {
    throw ModelSearchConfigurationError(
        "model complexity must be finite and non-negative");
  }
  MultiStructureRecord& selected = impl_->structure(structure_key);
  selected.effective_sample_size = effective_sample_size;
  selected.complexity = complexity;
  selected.use_bic = true;
}

void MultiStructureModelSearchProblem::clear_structure_bic_metadata(
    const std::string& structure_key) {
  MultiStructureRecord& selected = impl_->structure(structure_key);
  selected.effective_sample_size = 0.0;
  selected.complexity = 0.0;
  selected.use_bic = false;
}

ModelSearchState MultiStructureModelSearchProblem::get_initial_state() {
  impl_->validate_registry();
  const MultiStructureRecord& selected =
      impl_->structure(impl_->initial_structure);
  const FitSearchSnapshot previous = impl_->capture();
  const std::string previous_structure = impl_->active_structure;
  try {
    impl_->apply_initial(selected);
    impl_->active_structure = impl_->initial_structure;
    std::vector<std::shared_ptr<GraphPort> > free_ports;
    for (std::size_t i = 0; i < impl_->parameter_order.size(); ++i) {
      const std::shared_ptr<GraphPort>& port =
          impl_->parameters.find(impl_->parameter_order[i])->second;
      if (!port->get_fixed()) free_ports.push_back(port);
    }
    if (!free_ports.empty() && !impl_->cancelled.load()) {
      if (selected.residual_key.empty()) {
        throw ModelSearchConfigurationError(
            "initial structure with free parameters requires residuals");
      }
      FitMinimizer minimizer;
      minimizer.set_parameter_ports(free_ports);
      minimizer.set_objective(selected.objective, selected.residual_key);
      impl_->last_status = minimizer.run();
      if (impl_->last_status < 1 || impl_->last_status > 4) {
        throw ModelSearchConfigurationError(
            "initial model-search structure did not converge");
      }
    } else {
      selected.objective->update();
      impl_->last_status = impl_->cancelled.load() ? -1 : 1;
    }
    const std::pair<double, bool> score = impl_->score_current(selected);
    impl_->snapshots[impl_->initial_structure] = impl_->capture();
    impl_->snapshot_structures[impl_->initial_structure] =
        impl_->initial_structure;
    return ModelSearchState(impl_->initial_structure,
                            impl_->initial_structure, score.first,
                            score.second);
  } catch (...) {
    impl_->restore_values(previous);
    impl_->active_structure = previous_structure;
    throw;
  }
}

ModelSearchActions MultiStructureModelSearchProblem::get_actions(
    const ModelSearchState& state) {
  const std::map<std::string, ModelSearchActions>::const_iterator found =
      impl_->actions.find(state.get_structure_key());
  if (found == impl_->actions.end()) return ModelSearchActions();
  return found->second;
}

ModelSearchState MultiStructureModelSearchProblem::evaluate(
    const ModelSearchState& parent, const ModelSearchAction& action) {
  impl_->validate_registry();
  impl_->last_status = 0;
  impl_->last_failure.clear();
  const std::map<std::string, FitSearchSnapshot>::const_iterator parent_snapshot =
      impl_->snapshots.find(parent.get_key());
  if (parent_snapshot == impl_->snapshots.end()) {
    throw ModelSearchConfigurationError(
        "parent model-search snapshot is not cached");
  }
  const std::map<std::string, std::string>::const_iterator parent_structure =
      impl_->snapshot_structures.find(parent.get_key());
  if (parent_structure == impl_->snapshot_structures.end() ||
      parent_structure->second != parent.get_structure_key()) {
    throw ModelSearchConfigurationError(
        "parent state structure does not match its cached snapshot");
  }
  const std::string target_key = action.get_predicted_state_key();
  const MultiStructureRecord& source =
      impl_->structure(parent.get_structure_key());
  const MultiStructureRecord& target = impl_->structure(target_key);
  try {
    impl_->restore_values(parent_snapshot->second);
    std::vector<std::shared_ptr<GraphPort> > free_ports =
        impl_->apply_transition(source, target);
    impl_->active_structure = target_key;
    if (impl_->cancelled.load()) {
      impl_->last_status = -1;
      impl_->last_failure = "cancelled before minimization";
      impl_->rollback(parent_snapshot->second, parent.get_structure_key());
      return parent;
    }
    if (!free_ports.empty()) {
      if (target.residual_key.empty()) {
        throw ModelSearchConfigurationError(
            "a structure with free parameters requires a residual output");
      }
      FitMinimizer minimizer;
      minimizer.set_parameter_ports(free_ports);
      minimizer.set_objective(target.objective, target.residual_key);
      IMP::Pointer<FitSearchCancelObserver> observer(
          new FitSearchCancelObserver(&impl_->cancelled));
      minimizer.set_observer(observer.get());
      impl_->last_status = minimizer.run();
      if (minimizer.get_cancelled() || impl_->cancelled.load()) {
        impl_->last_status = -1;
        impl_->last_failure = "minimization cancelled";
        impl_->rollback(parent_snapshot->second, parent.get_structure_key());
        return parent;
      }
      if (impl_->last_status < 1 || impl_->last_status > 4) {
        std::ostringstream message;
        message << "minimizer did not converge (status " << impl_->last_status
                << ")";
        impl_->last_failure = message.str();
        impl_->rollback(parent_snapshot->second, parent.get_structure_key());
        return parent;
      }
    } else {
      target.objective->update();
      impl_->last_status = 1;
    }
    const std::pair<double, bool> score = impl_->score_current(target);
    std::ostringstream state_key;
    state_key << target_key << "@" << impl_->next_snapshot++;
    const std::string key = state_key.str();
    impl_->snapshots[key] = impl_->capture();
    impl_->snapshot_structures[key] = target_key;
    return ModelSearchState(key, target_key, score.first, score.second);
  } catch (const std::exception& error) {
    impl_->last_failure = error.what();
    impl_->last_status = 0;
    impl_->rollback(parent_snapshot->second, parent.get_structure_key());
    return parent;
  }
}

void MultiStructureModelSearchProblem::request_cancel() {
  impl_->cancelled.store(true);
}

void MultiStructureModelSearchProblem::clear_cancel() {
  impl_->cancelled.store(false);
}

void MultiStructureModelSearchProblem::activate_state(
    const ModelSearchState& state) {
  const std::map<std::string, std::string>::const_iterator found =
      impl_->snapshot_structures.find(state.get_key());
  if (found == impl_->snapshot_structures.end() ||
      found->second != state.get_structure_key()) {
    throw ModelSearchConfigurationError(
        "model-search state does not match a cached structure");
  }
  restore_state(state.get_key());
}

bool MultiStructureModelSearchProblem::has_cached_state(
    const std::string& state_key) const {
  return impl_->snapshots.count(state_key) != 0;
}

std::vector<double> MultiStructureModelSearchProblem::get_cached_values(
    const std::string& state_key) const {
  const std::map<std::string, FitSearchSnapshot>::const_iterator found =
      impl_->snapshots.find(state_key);
  if (found == impl_->snapshots.end()) return std::vector<double>();
  return found->second.values;
}

std::vector<int> MultiStructureModelSearchProblem::get_cached_fixed(
    const std::string& state_key) const {
  const std::map<std::string, FitSearchSnapshot>::const_iterator found =
      impl_->snapshots.find(state_key);
  if (found == impl_->snapshots.end()) return std::vector<int>();
  return found->second.fixed;
}

void MultiStructureModelSearchProblem::restore_state(
    const std::string& state_key) {
  const std::map<std::string, FitSearchSnapshot>::const_iterator found =
      impl_->snapshots.find(state_key);
  const std::map<std::string, std::string>::const_iterator structure =
      impl_->snapshot_structures.find(state_key);
  if (found == impl_->snapshots.end() ||
      structure == impl_->snapshot_structures.end()) {
    throw ModelSearchConfigurationError("unknown cached model state '" +
                                         state_key + "'");
  }
  impl_->restore(found->second, structure->second);
}

const std::string& MultiStructureModelSearchProblem::get_active_structure()
    const {
  return impl_->active_structure;
}

std::shared_ptr<GraphNode>
MultiStructureModelSearchProblem::get_active_objective() const {
  if (impl_->active_structure.empty()) return std::shared_ptr<GraphNode>();
  return impl_->structure(impl_->active_structure).objective;
}

int MultiStructureModelSearchProblem::get_last_fit_status() const {
  return impl_->last_status;
}

const std::string& MultiStructureModelSearchProblem::get_last_failure() const {
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
