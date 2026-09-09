/**
 *  \file IMP/bff/Sampler.h
 *  \brief ChiSurf's MCMC samplers over a bff Port/Node model, in C++.
 *
 *  Phase 6 of removing chinet from chisurf: the samplers of
 *  chisurf/core/fitting/sample.py and ensemble.py ported 1:1 onto the bff
 *  Port/Node runtime, so that a Markov-chain step never crosses the
 *  SWIG boundary. Phases 1-5 moved the parameter runtime (Port, Node,
 *  Session, FactorGraph) into bff and chisurf onto it; what stayed slow
 *  was the sampler loop itself -- measured at 1.67 us per port set+get
 *  across the boundary against 0.06 us for a Python attribute, paid once
 *  per proposal per parameter. This class writes the walker into the
 *  input ports, evaluates the node graph and accepts, all in C++.
 *
 *  Four algorithms behind one interface. The interface is the point:
 *  a caller names a backend and everything else -- the objective, the
 *  bounds, the walkers, the chain, the diagnostics -- is the same, so
 *  swapping a sampler is a string and not a rewrite.
 *
 *  - "stretch": the affine-invariant ensemble stretch move of
 *    Goodman & Weare as chisurf's EnsembleSampler implements it. The
 *    walkers are split into two randomly assigned halves; a walker of
 *    the active half is proposed along the line to a walker of the
 *    complementary half, stretched by z drawn from g(z) ~ 1/sqrt(z) on
 *    [1/a, a] (a the stretch scale, 2.0 by default), and accepted with
 *    probability min(1, z^(n-1) p(y)/p(x)). One log-posterior evaluation
 *    per walker per step.
 *
 *  - "de": Differential-Evolution MCMC (ter Braak) with chisurf's two
 *    refinements: every tenth generation proposes with gamma = 1 (a
 *    direct mode-to-mode jump), and a fraction of proposals are snooker
 *    updates along the line to another chain. The population is seeded
 *    and jittered exactly as sample_differential_evolution seeds it.
 *
 *  - "slice": the ensemble slice sampler (Karamanis & Beutler; what
 *    `zeus` implements). Like "stretch" the walkers are split into two
 *    halves and a walker moves along the line to one of the other half,
 *    but the move is a *slice* rather than a Metropolis proposal: a
 *    height is drawn under the density, the interval is stepped out
 *    until both ends are below it, and points are drawn and the interval
 *    shrunk until one lands above. Every step is accepted, which is what
 *    it buys -- no tuning of an acceptance rate, and no rejected
 *    evaluations. The direction scale mu is tuned during warm-up from
 *    the ratio of expansions to contractions and frozen afterwards.
 *
 *    **The stepping-out must be allowed to overshoot.** Capping it
 *    tightly looks harmless -- every draw is still inside the slice --
 *    and it silently truncates the tails: measured on a Gaussian target,
 *    a cap that binds gives a posterior some 30% too narrow, with no
 *    diagnostic saying so. #get_slice_truncations counts the times the
 *    cap bound, and it should be zero on a converged run.
 *
 *  - "metropolis": the blocked random-walk Metropolis of
 *    walk_mcmc_blocked. Blocks come from an attached FactorGraph's
 *    sampling blocks (two parameters share a block when the same set of
 *    datasets depends on both) or from an explicit partition, and the
 *    fallback is a single block covering every parameter -- an ordinary
 *    full-vector random walk, chisurf's own degenerate case. Proposals
 *    are full-covariance per block, seeded from a diagonal scaled by the
 *    parameter values (chisurf's fallback -- its curvature seed reads
 *    Fit.covariance_matrix, which has no bff counterpart), adapted
 *    during warm-up by Nesterov dual averaging with Stan's windowed
 *    schedule, then frozen for the recorded chain.
 *
 *  Semantics ported exactly:
 *
 *  - the log-posterior is chisurf's lnprob_parts: a box test against the
 *    parameter bounds that rejects (-inf) before anything is evaluated,
 *    then the sum of the ports' prior log-densities (the JSON
 *    specification mirrored onto each port, chisurf's _smooth_prior
 *    contract), then -chi2/2 of the node-graph output with the chi2max
 *    cutoff. Proposals out of bounds are REJECTED, never clipped --
 *    chisurf's lnprior short-circuits to -inf.
 *  - acceptance in log space with a strict >, as every chisurf sampler
 *    does; non-finite log-posteriors reject.
 *  - every tunable keeps chisurf's default (stretch scale 2.0, DE jitter
 *    1e-4 and snooker fraction 0.1, warm-up lengths, per-block target
 *    acceptance rates, the 2.38/sqrt(d) optimal scaling).
 *
 *  The random source is std::mt19937_64, reseedable; runs under a fixed
 *  seed are reproducible. Stream parity with chisurf's numpy generators
 *    is neither claimed nor required -- the contract is statistical.
 *
 *  What is deliberately NOT here:
 *
 *  - walk_mcmc's Robbins-Monro variant; the blocked sampler with one
 *    block is chisurf's own fallback for an unstructured fit and covers
 *    that role.
 *  - Callable/Product priors: chisurf keeps runtime-only priors on the
 *    Python side and leaves the port's prior empty; the C++ sampler sees
 *    exactly the serialisable specifications, and such priors must be
 *    folded into the objective node.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_SAMPLER_H
#define IMPBFF_SAMPLER_H

#include <IMP/bff/bff_config.h>

#include <functional>
#include <limits>
#include <stdexcept>
#include <memory>
#include <random>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

class Port;
class Node;
class FactorGraph;

//! Raised for a misconfigured sampler (no objective, fixed ports, a
//! degenerate ensemble, an unknown algorithm...).
/*!
    std::domain_error so the wrapper maps it to a ValueError -- the
    contract chisurf's samplers raise under. LinkCycleError (Port.h)
    derives domain_error for the same reason: IMP's wrapper handler maps
    it to IMP.ValueException, a ValueError, which is what chinet's and
    chisurf's callers catch.
*/
class IMPBFFEXPORT SamplerConfigurationError : public std::domain_error {
 public:
  explicit SamplerConfigurationError(const std::string& what_arg)
      : std::domain_error(what_arg) {}
};

//! ChiSurf's samplers (stretch, differential evolution, blocked
//! Metropolis) over free-parameter ports and a node-graph objective.
class IMPBFFEXPORT Sampler {
 public:
  //! Configure the algorithm and the seed; everything else is set later.
  /*!
      \param[in] algorithm "stretch", "de" (differential evolution) or
                 "metropolis" (blocked random walk). The names chisurf's
                 resolve_sampler accepts are aliases of these.
      \param[in] seed the RNG seed; a fixed seed gives a reproducible run
  */
  explicit Sampler(const std::string& algorithm = "stretch",
                   unsigned int seed = 42u);
  ~Sampler();

  // ------------------------------------------------------------ parameters

  //! The free parameters, as the ports the walker is written into.
  /*!
      Values start from the ports' current values, bounds from the ports
      that enforce them (a port without enforcement is unbounded) and
      priors from the ports' prior specifications. A fixed port cannot be
      sampled and is refused here.
  */
  void set_parameter_ports(
      const std::vector<std::shared_ptr<Port> >& parameters);
  //! The parameter ports (empty in plain-vector mode).
  std::vector<std::shared_ptr<Port> > get_parameter_ports() const;
  //! The parameter names (the ports' names, or x0, x1, ... without ports).
  std::vector<std::string> get_parameter_names() const;
  //! Override the starting values (the ports' current values by default).
  void set_initial_values(const std::vector<double>& values);
  //! The starting values the next run initialises from.
  std::vector<double> get_initial_values() const;
  //! Explicit box bounds (NaN for infinite); overrides the ports' bounds.
  void set_bounds(const std::vector<double>& lower,
                  const std::vector<double>& upper);

  // ------------------------------------------------------------- objective

  //! The primary objective: a node graph evaluated entirely in C++.
  /*!
      \param[in] node the node whose update() computes the objective
      \param[in] output_port the output port read after the update; its
                 value is chi^2 unless set_output_is_log_likelihood()
                 says otherwise
  */
  void set_objective(std::shared_ptr<Node> node,
                     const std::string& output_port = "chi2");
  //! The objective node, or a null pointer.
  std::shared_ptr<Node> get_objective() const;
  //! Read the output port as a log-likelihood instead of a chi^2.
  /*!
      The chi^2 reading is chisurf's (lnlike = -chi2/2, -inf above
      chi2max); the log-likelihood reading is for a node that computes
      one directly. Reported chi^2 is then -2 lnlike, as chisurf's
      ensemble_result reports it when no blobs are available.
  */
  void set_output_is_log_likelihood(bool v);
  //! Whether the output port is read as a log-likelihood.
  bool get_output_is_log_likelihood() const;

  //! The secondary objective: a plain log-posterior function.
  /*!
      The function returns the log-posterior (up to a constant) of a
      parameter vector; priors and bounds are then the caller's business,
      exactly as they are for a log_prob_fn handed to chisurf's ensemble
      samplers. Replaces the node objective when set.
  */
  void set_objective_function(
      std::function<double(const std::vector<double>&)> objective);
  //! Whether an objective (node or function) has been set.
  bool has_objective() const;

  // -------------------------------------------------------------- blocking

  //! Attach the fit's factor graph; its sampling blocks become the blocks.
  /*!
      Keys the graph cannot place in the parameter vector are dropped, as
      chisurf's _default_blocks drops them; a cover that misses parameters
      falls back to one block over everything, with the same warning
      chisurf logs. Only the "metropolis" algorithm uses blocks.

      The graph is borrowed (FactorGraph is a plain class, not
      shared_ptr-owned -- the caller keeps it alive for the sampler's
      lifetime, as chisurf's fit keeps its graph).
  */
  void set_factor_graph(FactorGraph* graph);
  //! The attached factor graph, or a null pointer.
  FactorGraph* get_factor_graph() const;
  //! An explicit block partition: flat indices and one size per block.
  void set_blocks(const std::vector<int>& flat_indices,
                  const std::vector<int>& block_sizes);
  //! The block partition the metropolis sweep proposes over.
  std::vector<std::vector<int> > get_blocks() const;
  //! One size per block (the SWIG-friendly read of get_blocks()).
  std::vector<int> get_block_sizes() const;

  // ------------------------------------------------------------ tunables

  //! Select the algorithm (aliases as in the constructor).
  void set_algorithm(const std::string& algorithm);
  //! The canonical algorithm name: "stretch", "de" or "metropolis".
  const std::string& get_algorithm() const;
  //! Reseed the RNG (also what a fixed-seed constructor started from).
  void set_seed(unsigned int seed);
  //! The seed.
  unsigned int get_seed() const;

  //! The number of walkers (stretch); 0 keeps chisurf's default.
  /*!
      chisurf defaults to max(2*ndim + 2, 10) and refuses fewer than
      2*ndim walkers (and fewer than 4 outright) unless
      live_dangerously is set.
  */
  void set_number_of_walkers(int n);
  //! The number of walkers the next stretch run uses.
  int get_number_of_walkers() const;
  //! The stretch scale a of the stretch move (2.0, chisurf's default).
  void set_stretch_scale(double a);
  //! The stretch scale.
  double get_stretch_scale() const;

  //! The slice sampler's direction scale (tuned during warm-up).
  double get_slice_mu() const;
  //! Set it explicitly, and stop tuning it.
  void set_slice_mu(double mu);
  //! Cap on stepping-out expansions per side; 0 restores the default.
  /*! Raise it, never lower it, unless a run is known to be pathological:
      a cap that binds truncates the tails silently. */
  void set_slice_max_steps(int n);
  int get_slice_max_steps() const;
  //! How often the stepping-out cap bound. Non-zero means the chain is
  //! truncated and its width is not to be believed.
  long get_slice_truncations() const;
  //! Skip the nwalkers >= 2*ndim requirement (chisurf's flag).
  void set_live_dangerously(bool v);
  //! Whether the walker-count requirement is skipped.
  bool get_live_dangerously() const;

  //! The DE population size; 0 keeps chisurf's max(8, 2*ndim).
  void set_number_of_chains(int n);
  //! The DE population size the next run uses.
  int get_number_of_chains() const;
  //! The DE jitter (relative width of the additive noise, 1e-4).
  void set_jitter(double jitter);
  //! The DE jitter.
  double get_jitter() const;
  //! The fraction of DE proposals that are snooker updates (0.1).
  void set_snooker(double fraction);
  //! The snooker fraction.
  double get_snooker() const;

  //! The relative proposal width of the diagonal covariance seed (0.1).
  void set_step_size(double step_size);
  //! The relative proposal width.
  double get_step_size() const;

  //! Seed the blocked proposal from a full covariance matrix.
  /*!
      chisurf's _seed_block_covariances: the curvature of the objective at
      the optimum is very nearly the ideal preconditioner, and the fit
      computes it for the error bars anyway -- the caller passes it once,
      at begin. Each block takes its submatrix; a block whose submatrix is
      not finite and positive definite falls back to the diagonal seed.
      A curvature-seeded block keeps its *shape* through warm-up (only the
      scale adapts): an exact curvature cannot be improved on by a short
      chain, a diagonal guess can. Only "metropolis" uses it.

      \param[in] cov (ndim, ndim) rows over the full parameter vector; an
                 empty vector clears the seed.
  */
  void set_proposal_covariance(
      const std::vector<std::vector<double> >& cov);
  //! The full-covariance seed, or an empty vector when unset.
  std::vector<std::vector<double> > get_proposal_covariance() const;

  //! Warm-up steps, discarded; -1 keeps chisurf's per-algorithm default.
  /*!
      None for stretch (chisurf's ensemble samplers have no warm-up), DE's
      min(500, max(50, n_steps/4)), the blocked sampler's
      clip(n_steps/20, 100, 500).
  */
  void set_n_adapt(int n);
  //! The warm-up length (-1 for the default).
  int get_n_adapt() const;

  //! The sampling temperature; above one flattens the posterior.
  void set_temp(double temp);
  //! The sampling temperature.
  double get_temp() const;
  //! The hard chi^2 cutoff; above it the log-likelihood is -inf.
  void set_chi2max(double chi2max);
  //! The chi^2 cutoff.
  double get_chi2max() const;

  //! The relative spread of the initial walker cloud (1e-3, chisurf's).
  /*!
      The scale of _ensemble_walker_start: the bounded range times 1e-4
      where the parameter has finite bounds, |value| * std otherwise, and
      the absolute std wherever that would leave a direction with no
      spread at all (a parameter the ensemble could then never move).
  */
  void set_walker_start_std(double std);
  //! The relative spread of the initial walker cloud.
  double get_walker_start_std() const;

  //! Explicit initial walker positions (rows; one per walker/chain).
  /*!
      Used as given, as a state handed to run_mcmc is; the degeneracy
      check still applies. Overrides the automatic spread.
  */
  void set_walker_start(const std::vector<std::vector<double> >& start);
  //! The explicit start, or an empty vector when the spread is automatic.
  std::vector<std::vector<double> > get_walker_start() const;

  // ----------------------------------------------------------------- run

  //! Run the chain: n_steps steps, every thin-th state recorded.
  /*!
      n_steps counts the steps each walker takes, so n_steps/thin states
      per walker are recorded -- chisurf's counting for every sampler
      here. A second call continues from where the first left off; a
      fresh chain wants reset(). DE restores the starting values into the
      parameter ports afterwards, as chisurf's DE sampler does.

      \param[in] n_steps steps per walker
      \param[in] thin record only every thin-th step (1 records all)
  */
  void run(int n_steps, int thin = 1);
  //! Take one step and record it (run(1, 1)).
  void step();
  //! Discard the chain, the bookkeeping and the ensemble; reseed.
  void reset();

  // -------------------------------------------------------------- results

  //! The recorded chain, flattened: one row per (step, walker).
  /*!
      (n_recorded * n_walkers, ndim) in chisurf's flat ordering -- the
      reshape of its (n_steps, nwalkers, ndim) chain under flat=True.
  */
  std::vector<std::vector<double> > get_chain() const;
  //! One walker's recorded chain: (n_recorded, ndim).
  std::vector<std::vector<double> > get_chain_of_walker(int walker) const;
  //! The current walker positions: one row per walker.
  std::vector<std::vector<double> > get_walkers() const;
  //! The recorded log-posteriors, flat, one per chain row.
  std::vector<double> get_log_prob() const;
  //! The recorded log-priors, flat (0 in function-objective mode).
  std::vector<double> get_lnprior() const;
  //! The recorded chi^2, flat (-2 lnpost in function-objective mode).
  std::vector<double> get_chi2() const;
  //! The overall acceptance rate of the recorded chain.
  double get_acceptance_rate() const;
  //! The per-walker acceptance fractions (a single rate otherwise).
  std::vector<double> get_acceptance_fractions() const;
  //! The per-block acceptance rates of the blocked sampler.
  std::vector<double> get_block_acceptance_rates() const;
  //! The number of log-posterior evaluations performed.
  unsigned int get_number_of_evaluations() const;
  //! The number of recorded states.
  unsigned int get_iteration() const;
  //! The number of sampled parameters.
  unsigned int get_number_of_parameters() const;
  //! A one-line summary: algorithm, dimension, walkers, acceptance.
  std::string describe() const;

  //! A per-step progress hook: observer(steps_done, steps_total).
  /*!
      C++-only (not SWIG-wrapped); called once per step of run().
  */
  void set_observer(std::function<void(int, int)> observer);

 private:
  //! chisurf's lnprob_parts: bounds box, priors, then the objective.
  struct Parts {
    //! log-posterior (lnlike + lnprior)
    double lnpost;
    //! log-prior alone (0 when only a function objective is set)
    double lnprior;
    //! chi^2 of the data misfit (inf when the box or a prior rejected)
    double chi2;
  };

  //! A prior specification parsed off a port's prior JSON.
  struct PriorSpec {
    std::string kind;
    //! The numeric fields, by name (mu, sigma, lb, ub, ...).
    std::vector<std::pair<std::string, double> > numbers;
    bool empty() const { return kind.empty(); }
    double get(const std::string& name, double fallback) const;
  };

  //! Nesterov dual averaging of a log step size (_DualAveraging).
  struct DualAveraging {
    double target = 0.234;
    double gamma = 0.05;
    double t0 = 10.0;
    double kappa = 0.75;
    double mu = 0.0;
    double log_eps = 0.0;
    double log_eps_bar = 0.0;
    double h_bar = 0.0;
    int counter = 0;
    void restart(double log_eps);
    double update(double alpha);
    double averaged() const { return log_eps_bar; }
  };

  //! Everything the blocked sampler keeps per block.
  struct BlockState {
    std::vector<int> indices;             // positions in the parameter vector
    std::vector<double> factor;           // flat lower-triangular Cholesky
    double log_scale = 0.0;
    DualAveraging adapter;
    long accepted = 0;                     // recorded-phase counters
    long proposed = 0;
    //! Seeded from the caller's curvature: warm-up adapts the scale only,
    //! never replaces the shape (chisurf's from_curvature flag).
    bool from_curvature = false;
  };

  //! The explicit partition set through set_blocks (empty when derived).
  std::vector<std::vector<int> > explicit_blocks_;
  //! The DE jitter widths, seeded from the starting values.
  std::vector<double> de_noise_;

  // ------------------------------------------------------------- internals
  //! Validate the configuration; throw SamplerConfigurationError.
  void validate() const;
  //! Derive bounds, priors and names from the ports (idempotent).
  void configure_from_ports();
  //! (Re)build the ensemble from the start, per algorithm.
  void initialize_ensemble();
  //! chisurf's _ensemble_walker_start: a spread, bounded walker cloud.
  std::vector<std::vector<double> > spread_walkers(int n) const;
  //! chisurf's walkers_independent: refuse a degenerate ensemble.
  bool walkers_independent(
      const std::vector<std::vector<double> >& coords) const;
  //! Evaluate the log-posterior parts of one parameter vector.
  Parts evaluate(const std::vector<double>& x);
  //! The log-prior alone: box, then the parsed priors (chisurf lnprior).
  double log_prior(const std::vector<double>& x) const;
  //! One prior's lnpdf (chisurf priors.py, kind for kind).
  static double prior_lnpdf(const PriorSpec& spec, double x);
  //! Parse {"kind": ..., name: number, ...} off a port's prior string.
  static PriorSpec parse_prior(const std::string& json);

  //! One ensemble step of the stretch move (EnsembleSampler._step).
  void stretch_step();
  //! One ensemble slice sweep: both halves, every walker, always accepted.
  void slice_step();
  //! Slice-sample one walker along \p direction from \p x; returns the new point.
  std::vector<double> slice_along(const std::vector<double>& x,
                                  const std::vector<double>& direction,
                                  double log_p_x, int* expansions,
                                  int* contractions, bool* truncated);
  //! One DE generation (sample_differential_evolution._generation).
  void de_generation(long generation_index);
  //! One blocked sweep (walk_mcmc_blocked._sweep).
  void blocked_sweep(bool adapt);
  //! Rebuild the block partition (factor graph, explicit, or single).
  void rebuild_blocks();
  //! Record the current ensemble into the chain.
  void record_state();
  //! Seed the block covariances and adapters (the warm-up's start).
  void seed_blocks();

  // ------------------------------------------------------------- the model
  std::vector<std::shared_ptr<Port> > parameters_;
  std::shared_ptr<Node> objective_node_;
  std::shared_ptr<Port> output_port_;
  bool output_is_log_likelihood_ = false;
  std::function<double(const std::vector<double>&)> objective_function_;
  FactorGraph* factor_graph_ = nullptr;  //!< borrowed

  unsigned int ndim_ = 0;
  std::vector<double> initial_values_;
  std::vector<double> lower_, upper_;      //!< NaN/inf mean infinite
  std::vector<PriorSpec> priors_;
  std::vector<std::string> names_;
  bool bounds_explicit_ = false;

  // -------------------------------------------------------- the tunables
  std::string algorithm_ = "stretch";
  unsigned int seed_ = 42u;
  int n_walkers_setting_ = 0;              //!< 0: chisurf default
  double stretch_scale_ = 2.0;
  //! The slice sampler's direction scale, tuned during warm-up.
  double slice_mu_ = 1.0;
  //! How many stepping-out expansions the cap cut short. Should stay 0.
  long slice_truncations_ = 0;
  long slice_expansions_ = 0;
  long slice_contractions_ = 0;
  //! Cap on stepping-out expansions per side. Generous on purpose.
  int slice_max_steps_ = 10000;
  bool slice_tuning_ = true;
  bool live_dangerously_ = false;
  int n_chains_setting_ = 0;               //!< 0: chisurf default
  double jitter_ = 1e-4;
  double snooker_ = 0.1;
  double step_size_ = 0.1;
  //! The full-covariance proposal seed (empty: diagonal fallback).
  std::vector<std::vector<double> > proposal_covariance_;
  int n_adapt_setting_ = -1;               //!< -1: chisurf default
  double temp_ = 1.0;
  double chi2max_ = std::numeric_limits<double>::infinity();
  double walker_start_std_ = 1e-3;
  std::vector<std::vector<double> > walker_start_override_;

  // ------------------------------------------------------- the run state
  mutable std::mt19937_64 rng_;
  std::vector<std::vector<double> > walkers_;   //!< rows: walkers/chains
  std::vector<Parts> walker_parts_;            //!< their log-posteriors
  std::vector<long> accepted_;                  //!< per walker, recorded
  std::vector<long> substep_accepted_;          //!< per walker, this substep
  std::vector<double> acceptance_fractions_;    //!< cached at record time
  std::vector<BlockState> blocks_;
  long de_accepted_ = 0, de_proposed_ = 0;      //!< recorded-phase DE
  std::vector<std::vector<double> > chain_;     //!< flat, per record
  std::vector<double> log_prob_, ln_prior_, chi2_;
  unsigned int iteration_ = 0;
  unsigned int n_evaluations_ = 0;
  bool initialized_ = false;
  int thin_ = 1;
  std::function<void(int, int)> observer_;
};

IMPBFF_END_NAMESPACE

#endif // IMPBFF_SAMPLER_H
