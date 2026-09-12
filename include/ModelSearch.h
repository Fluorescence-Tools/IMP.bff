/**
 *  \file IMP/bff/ModelSearch.h
 *  \brief Model-independent Monte Carlo tree search for fitted structures.
 *
 *  The tree knows only opaque state/action keys and scalar rewards.  The
 *  problem owns every model-specific object and performs each transition in
 *  C++; this keeps fitting data and the inner optimiser on the same side of
 *  the Python boundary.  A future adapter may therefore keep one state as a
 *  single fit and another as a heterogeneous global fit without teaching the
 *  search either representation.
 */
#ifndef IMPBFF_MODELSEARCH_H
#define IMPBFF_MODELSEARCH_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/IMPCompatibility.h>

#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

class IMPBFFEXPORT ModelSearchConfigurationError
    : public std::domain_error {
 public:
  explicit ModelSearchConfigurationError(const std::string& what)
      : std::domain_error(what) {}
};

//! One fully evaluated model state; its key is opaque to the tree.
class IMPBFFEXPORT ModelSearchState {
 public:
  ModelSearchState();
  ModelSearchState(const std::string& key, double reward,
                      bool acceptable = false);
  ModelSearchState(const std::string& key,
                      const std::string& structure_key, double reward,
                      bool acceptable = false);
  const std::string& get_key() const;
  //! Canonical structural identity; several fitted snapshots may share it.
  const std::string& get_structure_key() const;
  double get_reward() const;
  bool get_acceptable() const;
  IMP_SHOWABLE_INLINE(ModelSearchState,
                      out << "ModelSearchState(" << key_ << ", "
                          << reward_ << ")");

 private:
  std::string key_;
  std::string structure_key_;
  double reward_;
  bool acceptable_;
};
IMP_VALUES(ModelSearchState, ModelSearchStates);

//! A structural move offered by a problem at an evaluated state.
class IMPBFFEXPORT ModelSearchAction {
 public:
  ModelSearchAction();
  ModelSearchAction(const std::string& key,
                       const std::string& predicted_state_key,
                       double prior = 1.0, bool terminal = false);
  const std::string& get_key() const;
  const std::string& get_predicted_state_key() const;
  double get_prior() const;
  bool get_terminal() const;
  IMP_SHOWABLE_INLINE(ModelSearchAction,
                      out << "ModelSearchAction(" << key_ << " -> "
                          << predicted_state_key_ << ")");

 private:
  std::string key_;
  std::string predicted_state_key_;
  double prior_;
  bool terminal_;
};
IMP_VALUES(ModelSearchAction, ModelSearchActions);

//! Model-specific topology and evaluation, implemented entirely in C++.
/*!
  `evaluate()` may return a key different from the action's predicted key.
  This is how an optimiser reports that a nominally larger model collapsed
  onto a canonical ancestor.  The tree records its reward, then makes that
  node a dead end so the collapsed structure is not searched repeatedly.
*/
class IMPBFFEXPORT ModelSearchProblem {
 public:
  virtual ~ModelSearchProblem();
  virtual ModelSearchState get_initial_state() = 0;
  virtual ModelSearchActions get_actions(
      const ModelSearchState& state) = 0;
  virtual ModelSearchState evaluate(
      const ModelSearchState& parent,
      const ModelSearchAction& action) = 0;
  //! Propagate a cooperative cancellation request into an active evaluator.
  virtual void request_cancel();
  virtual void clear_cancel();
  //! Make a cached state current after search (the root when cancelled).
  virtual void activate_state(const ModelSearchState& state);
};

//! A callback-free finite problem, useful for persisted/pre-scored graphs.
/*!
  This is also the executable contract fixture for the SWIG surface.  Real
  fitting adapters subclass #ModelSearchProblem in C++ and keep their
  heterogeneous model states behind the same opaque keys.
*/
class IMPBFFEXPORT TabularModelSearchProblem : public ModelSearchProblem {
 public:
  TabularModelSearchProblem();
  void add_state(const std::string& key, double reward,
                 bool acceptable = false);
  void set_initial_state(const std::string& key);
  void add_action(const std::string& parent_key, const std::string& action_key,
                  const std::string& result_key, double prior = 1.0,
                  bool terminal = false);
  ModelSearchState get_initial_state() override;
  ModelSearchActions get_actions(
      const ModelSearchState& state) override;
  ModelSearchState evaluate(const ModelSearchState& parent,
                               const ModelSearchAction& action) override;
  int get_number_of_evaluations() const;
  void reset_number_of_evaluations();

 private:
  struct Impl;
  std::shared_ptr<Impl> impl_;
};

class GraphNode;
class GraphPort;

//! A live, callback-free adapter from declarative structures to FitMinimizer.
/*!
  Structures name the parameter groups that are free.  Transitions restore
  their parent's cached snapshot, apply the target fix/free mask and any seed
  for a newly enabled group, then run #FitMinimizer on the resulting free
  owner ports.  Every successful state caches its own values and fixed mask.

  The graph is shared and is therefore evaluated sequentially.  There is no
  graph clone and no hidden parallelism: each transition is a transaction over
  one graph, and the next transition begins by restoring its own parent.  On
  cancellation, non-convergence or an exception the parent snapshot is
  restored and the transition collapses to that parent.

  Ports must be scalar, unlinked canonical owners.  A caller with linked ports
  supplies each link target once; accepting followers as independent state
  would make snapshot restoration order-dependent.
*/
class IMPBFFEXPORT FittingModelSearchProblem : public ModelSearchProblem {
 public:
  FittingModelSearchProblem();
  FittingModelSearchProblem(std::shared_ptr<GraphNode> objective,
                           const std::string& residual_key = "residuals");
  ~FittingModelSearchProblem();

  void set_objective(std::shared_ptr<GraphNode> objective,
                     const std::string& residual_key = "residuals");
  std::shared_ptr<GraphNode> get_objective() const;

  //! Declare one independently fixable group and optional enable-time seeds.
  void add_parameter_group(
      const std::string& key,
      const std::vector<std::shared_ptr<GraphPort> >& ports,
      const std::vector<double>& enable_values = std::vector<double>());
  std::vector<std::string> get_parameter_group_keys() const;

  //! Declare a structure by the groups that are free in it.
  void add_structure(const std::string& key,
                     const std::vector<std::string>& free_groups);
  void set_initial_structure(const std::string& key);
  void add_action(const std::string& parent_structure,
                  const std::string& action_key,
                  const std::string& result_structure, double prior = 1.0,
                  bool terminal = false);

  //! Use this scalar objective output as reward (higher is better).
  void set_score_output(const std::string& key);
  void clear_score_output();
  const std::string& get_score_output() const;
  //! Otherwise reward is -chi2/2 - this penalty times the free-port count.
  void set_complexity_penalty(double value);
  double get_complexity_penalty() const;
  //! Optional scalar/bool output defining result.acceptable.
  void set_acceptable_output(const std::string& key);
  void clear_acceptable_output();

  ModelSearchState get_initial_state() override;
  ModelSearchActions get_actions(
      const ModelSearchState& state) override;
  ModelSearchState evaluate(const ModelSearchState& parent,
                               const ModelSearchAction& action) override;
  void request_cancel() override;
  void clear_cancel() override;
  void activate_state(const ModelSearchState& state) override;

  //! Snapshot inspection and explicit winner application.
  bool has_cached_state(const std::string& state_key) const;
  std::vector<double> get_cached_values(const std::string& state_key) const;
  std::vector<int> get_cached_fixed(const std::string& state_key) const;
  void restore_state(const std::string& state_key);
  int get_last_fit_status() const;
  const std::string& get_last_failure() const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
  FittingModelSearchProblem(const FittingModelSearchProblem&) = delete;
  FittingModelSearchProblem& operator=(const FittingModelSearchProblem&) = delete;
};

//! Search controls, independent of any fitting model or data family.
class IMPBFFEXPORT ModelSearchConfig {
 public:
  ModelSearchConfig();
  void set_number_of_simulations(int value);
  int get_number_of_simulations() const;
  void set_c_puct(double value);
  double get_c_puct() const;
  void set_reward_scale(double value);
  double get_reward_scale() const;
  void set_dirichlet_alpha(double value);
  double get_dirichlet_alpha() const;
  void set_dirichlet_fraction(double value);
  double get_dirichlet_fraction() const;
  void set_seed(unsigned int value);
  unsigned int get_seed() const;
  IMP_SHOWABLE_INLINE(ModelSearchConfig,
                      out << "ModelSearchConfig(" << n_simulations_
                          << " simulations)");

 private:
  int n_simulations_;
  double c_puct_;
  double reward_scale_;
  double dirichlet_alpha_;
  double dirichlet_fraction_;
  unsigned int seed_;
};
IMP_VALUES(ModelSearchConfig, ModelSearchConfigs);

//! The best evaluated state and an auditable summary of one run.
class IMPBFFEXPORT ModelSearchResult {
 public:
  ModelSearchResult();
  const ModelSearchState& get_root_state() const;
  const ModelSearchState& get_best_state() const;
  const std::vector<std::string>& get_best_path() const;
  int get_number_of_states_evaluated() const;
  int get_number_of_simulations() const;
  double get_improvement() const;
  bool get_acceptable() const;
  bool get_cancelled() const;
  IMP_SHOWABLE_INLINE(ModelSearchResult,
                      out << "ModelSearchResult(" << n_simulations_
                          << " simulations, best=" << best_state_.get_key()
                          << ")");

 private:
  friend class ModelSearch;
  ModelSearchState root_state_;
  ModelSearchState best_state_;
  std::vector<std::string> best_path_;
  int n_states_evaluated_;
  int n_simulations_;
  double improvement_;
  bool acceptable_;
  bool cancelled_;
};
IMP_VALUES(ModelSearchResult, ModelSearchResults);

//! PUCT traversal, lazy expansion, cancellation, and bookkeeping in C++.
class IMPBFFEXPORT ModelSearch {
 public:
  ModelSearch();
  explicit ModelSearch(std::shared_ptr<ModelSearchProblem> problem);
  ~ModelSearch();
  void set_problem(std::shared_ptr<ModelSearchProblem> problem);
  std::shared_ptr<ModelSearchProblem> get_problem() const;
  void set_config(const ModelSearchConfig& config);
  ModelSearchConfig get_config() const;

  //! Thread-safe cooperative cancellation, checked between simulations.
  void request_cancel();
  void clear_cancel();
  bool get_cancel_requested() const;

  //! Run from the problem's evaluated initial state.
  ModelSearchResult run();

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
  ModelSearch(const ModelSearch&) = delete;
  ModelSearch& operator=(const ModelSearch&) = delete;
};

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_MODELSEARCH_H
