/**
 *  \file IMP/bff/FitMinimizer.h
 *  \brief ChiSurf's bounded least-squares optimiser over a bff GraphPort/GraphNode model, in C++.
 *
 *  The deterministic counterpart of `MCMCSampler.h`, and the last piece of a fit
 *  that still returned to Python on every iteration. `MCMCSampler` moved the
 *  Markov walk into C++; this moves `fit.run()` -- the Levenberg-Marquardt
 *  minimisation chisurf's `core/math/optimization/leastsqbound.py` performs
 *  through `scipy.optimize._minpack._lmdif`.
 *
 *  **Why, given that MINPACK is already compiled.** Measured on this tree
 *  (okf/log.md, 2026-09-01 (6)): of a 512-point three-parameter fit, MINPACK's
 *  own arithmetic is about 4% and the residual *callback* is 57%. Replacing
 *  one compiled optimiser with another buys nothing on its own -- and that
 *  measurement is the reason this class exists in the shape it does. The cost
 *  is the boundary, not the algorithm: a Python-driven loop crosses once to
 *  write each parameter, once to evaluate the model, once to form residuals,
 *  every iteration. What removes that is not a faster optimiser but an
 *  optimiser on the *same side of the boundary as the objective*, so that the
 *  whole loop -- parameters, model, residuals, chi-square -- runs in C++ and
 *  the caller crosses once per `run()` rather than once per part per
 *  iteration.
 *
 *  So the point of the port is composition, not speed in isolation:
 *  `GraphExpression -> FitChiSquared -> FitMinimizer` is a fit that never re-enters the
 *  interpreter. A model bff cannot represent still gains, but less: a Python
 *  `GraphNode` director costs one crossing per residual evaluation instead of the
 *  four or five the numpy path pays.
 *
 *  **The algorithm is MINPACK's, ported rather than re-derived.** `lmdif`
 *  with its `enorm`, `fdjac2`, `qrfac`, `qrsolv`, `lmpar` and `covar`
 *  (More, Garbow & Hillstrom, ANL-80-74; public domain), so that a fit that
 *  converged under scipy converges here to the same answer from the same
 *  start. Every default is chisurf's: `ftol = xtol = 1.49012e-8`,
 *  `gtol = 0`, `factor = 100`, `maxfev = 200 * (n + 1)` for the
 *  forward-difference Jacobian.
 *
 *  **Bounds are chisurf's transform, not a projection.** `leastsqbound` maps
 *  each bounded parameter to an unconstrained internal coordinate --
 *  `lower + (upper - lower)/2 * (sin(xi) + 1)` for a two-sided bound,
 *  `lower - 1 + sqrt(xi^2 + 1)` for a lower bound alone, the mirror of it for
 *  an upper -- optimises there, and maps back. A bound that is `None`,
 *  infinite or NaN means "no constraint" (chisurf's `_is_unbounded`, which
 *  fixed a real defect: a parameter declared `(-inf, inf)` used to take the
 *  two-sided branch and poison the internal vector with NaN).
 *
 *  **The covariance comes out in external coordinates.** `covar()` is run on
 *  the final `R` after its columns are divided by the transform's gradient,
 *  which is exactly what `leastsqbound`'s `full_output` branch does with
 *  `_internal2external_grad`. Reporting the internal covariance would be a
 *  silently wrong error bar, and the whole reason chisurf recomputes a
 *  Jacobian after every fit; here the Jacobian the optimiser already built is
 *  reused, correctly.
 *
 *  \see MCMCSampler, FitChiSquared, GraphExpression, GraphNode
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_FITMINIMIZER_H
#define IMPBFF_FITMINIMIZER_H

#include <IMP/bff/bff_config.h>
// IMP_OBJECT_METHODS, used by FitMinimizerObserver below. IMPCompatibility.h is what
// resolves the IMP macros in *either* configuration -- IMP's own headers in
// the module build, the standalone definitions otherwise -- so a header that
// uses one has to include it rather than rely on a neighbour having done so.
#include <IMP/bff/IMPCompatibility.h>

#include <IMP/Object.h>
#include <IMP/Pointer.h>

#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

class GraphPort;
class GraphNode;

//! Raised for a misconfigured minimiser (no objective, fixed ports, a
//! bounds list that does not match the parameters...).
/*!
    `std::domain_error` so the wrapper maps it to a `ValueError`, which is
    the contract chisurf's optimiser raises under -- the same choice
    `MCMCSamplerConfigurationError` and `GraphLinkCycleError` make.
*/
class IMPBFFEXPORT FitMinimizerConfigurationError : public std::domain_error {
 public:
  explicit FitMinimizerConfigurationError(const std::string& what_arg)
      : std::domain_error(what_arg) {}
};

//! What a caller is told while a minimisation runs, and how it stops one.
/*!
    A SWIG director, so a Python caller -- a GUI progress dialog, most of the
    time -- subclasses it and overrides `report()`. This is the same shape
    `RRTCollision` uses for the growth loop: the C++ loop asks, the Python
    subclass answers.

    `report()` is called once per residual evaluation. **Returning false
    cancels the run**, which is how chisurf's `OptimizationCancelled` travels
    without an exception having to cross the boundary in the middle of a
    Fortran-shaped loop. The optimiser then stops at the *last accepted*
    parameter vector -- not the trial point it was evaluating -- so a
    cancelled fit leaves the model somewhere it has actually been.
*/
class IMPBFFEXPORT FitMinimizerObserver : public IMP::Object {
 public:
  explicit FitMinimizerObserver(std::string name = "FitMinimizerObserver%1%")
      : IMP::Object(name) {}

  //! One residual evaluation happened; return false to cancel.
  /*!
      \param[in] n_evaluations residual evaluations completed so far
      \param[in] total the evaluation budget to report against, from
                 `minimizer_reported_total()` -- an estimate that always
                 moves and never arrives, not a limit
      \param[in] chi2 the sum of squared residuals at this evaluation
  */
  virtual bool report(int n_evaluations, int total, double chi2) {
    return true;
  }

  IMP_OBJECT_METHODS(FitMinimizerObserver);
};

//! Levenberg-Marquardt iterations a well-posed fit typically needs.
/*! chisurf's `_EXPECTED_ITERATIONS`. Each iteration costs `n + 1` residual
    evaluations -- `n` for the forward-difference Jacobian and one for the
    trial step. */
IMPBFFEXPORT extern const int MINIMIZER_EXPECTED_ITERATIONS;

//! An evaluation budget a progress bar can usefully report against.
/*! chisurf's `_expected_evaluations`. MINPACK's `200 * (n + 1)` is a
    give-up limit, not an expectation, and using it as a denominator makes a
    fast fit look broken. */
IMPBFFEXPORT int minimizer_expected_evaluations(int n, int maxfev);

//! The total to report `nfev` against, so a progress bar always moves.
/*! chisurf's `_reported_total`. The reported fraction is asymptotic,
    `0.99 * (1 - 2^(-nfev / (2 * expected)))`, because a linear bar against a
    guessed total saturates on exactly the fits long enough for anyone to be
    watching -- a real MFD fit takes 450 evaluations against an estimate of
    42, and sat at 98-99% for 91% of its runtime. */
IMPBFFEXPORT int minimizer_reported_total(int nfev, int expected);

//! ChiSurf's bounded Levenberg-Marquardt over free-parameter ports and a
//! node-graph objective.
class IMPBFFEXPORT FitMinimizer {
 public:
  //! Configure the algorithm; everything else is set later.
  /*!
      \param[in] algorithm "leastsq" -- MINPACK's `lmdif` with chisurf's
                 bounds transform. The only algorithm today; the name is
                 taken now so that adding one later is not an API change.
  */
  explicit FitMinimizer(const std::string& algorithm = "leastsq");
  ~FitMinimizer();

  // ------------------------------------------------------------ parameters

  //! The free parameters, as the ports the trial vector is written into.
  /*!
      Values start from the ports' current values and bounds from the ports
      that enforce them (`GraphPort::get_is_bounded`); a port without enforcement
      is unbounded in that direction. A fixed port cannot be optimised and
      is refused here.
  */
  void set_parameter_ports(
      const std::vector<std::shared_ptr<GraphPort> >& parameters);
  std::vector<std::shared_ptr<GraphPort> > get_parameter_ports() const;
  std::vector<std::string> get_parameter_names() const;

  //! The starting point. Defaults to the ports' current values.
  void set_initial_values(const std::vector<double>& values);
  std::vector<double> get_initial_values() const;

  //! Bounds, overriding whatever the ports carry.
  /*! Non-finite entries mean "no constraint in that direction", which is
      chisurf's `_is_unbounded` -- `None`, `+/-inf` and `nan` alike. */
  void set_bounds(const std::vector<double>& lower,
                  const std::vector<double>& upper);
  std::vector<double> get_lower_bounds() const;
  std::vector<double> get_upper_bounds() const;

  // ------------------------------------------------------------- objective

  //! The node graph producing the residual vector.
  /*!
      \param[in] node the node to evaluate; `FitChiSquared` is the intended
                 one, but any node works -- including a Python `GraphNode`
                 director wrapping a model bff cannot represent
      \param[in] residual_key the node's output port carrying the residual
                 vector; empty keeps the current key (default "residuals")

      The optimiser writes the trial vector into the parameter ports, calls
      `GraphNode::update()`, and reads the residuals off that port. Nothing else
      crosses.
  */
  void set_objective(std::shared_ptr<GraphNode> node,
                     const std::string& residual_key = "");
  std::shared_ptr<GraphNode> get_objective() const;
  void set_residual_port_key(const std::string& key);
  const std::string& get_residual_port_key() const;

  //! A C++ residual function instead of a graph. Not SWIG-wrapped.
  /*! The module carries no `std_function.i`, and a Python callable is a
      `GraphNode` director anyway -- the same division `MCMCSampler` makes. */
  void set_residual_function(
      std::function<std::vector<double>(const std::vector<double>&)> f);
  bool has_objective() const;

  //! Where progress is reported and cancellation is asked for.
  /*! Null (the default) means neither. */
  void set_observer(FitMinimizerObserver* observer);
  FitMinimizerObserver* get_observer() const;
  //! Stop reporting; equivalent to `set_observer(nullptr)` from Python,
  //! where a null pointer is awkward to spell.
  void clear_observer();

  // --------------------------------------------------------------- options

  void set_algorithm(const std::string& algorithm);
  const std::string& get_algorithm() const;

  //! Relative error desired in the sum of squares. chisurf: 1.49012e-8.
  void set_ftol(double v);
  double get_ftol() const;
  //! Relative error desired in the approximate solution. chisurf: 1.49012e-8.
  void set_xtol(double v);
  double get_xtol() const;
  //! Orthogonality desired between the residuals and the Jacobian columns.
  void set_gtol(double v);
  double get_gtol() const;
  //! Hard limit on residual evaluations; 0 means `200 * (n + 1)`.
  void set_maxfev(int v);
  int get_maxfev() const;
  //! Step length for the forward-difference Jacobian; 0 means machine eps.
  void set_epsfcn(double v);
  double get_epsfcn() const;
  //! Initial step bound, `factor * ||diag * x||`. In (0.1, 100]; chisurf: 100.
  void set_factor(double v);
  double get_factor() const;
  //! Per-variable scale factors; empty means MINPACK's internal scaling.
  void set_diag(const std::vector<double>& diag);
  std::vector<double> get_diag() const;

  // ------------------------------------------------------------------- run

  //! Minimise, and return MINPACK's `info`.
  /*!
      1, 2, 3 and 4 are convergence (ftol, xtol, both, gtol); 5 is
      `maxfev` reached; 6, 7 and 8 say a tolerance is too small to make
      further progress; 0 is improper input. -1 means the observer
      cancelled.

      The parameter ports are left holding the solution, and the objective
      graph is left evaluated there -- so a caller reads the fitted model
      curve off the graph without evaluating anything again.
  */
  int run();

  //! Reset the results, keeping the configuration.
  void reset();

  //! Score many candidate parameter vectors in one crossing of the boundary.
  /*!
      \param[in] in_candidates the candidates, row-major `n_rows x n_cols`
      \param[in] n_rows how many candidates
      \param[in] n_cols must equal the number of free parameters
      \param[out] out_view `n_rows` objective values, allocated here
      \param[out] n_out_view how many, i.e. \p n_rows

      The objective is the sum of the squared residuals -- the same number
      `get_chi2()` reports -- and a NaN becomes `+inf`, so a sampler rejects
      the candidate rather than propagating the NaN into a posterior. That is
      `FitChiSquared`'s convention and `FitJointChiSquared`'s.

      **Why this exists.** A chi-square surface, a support-plane interval, a
      population sampler and a random restart all ask the same question of a
      few hundred to a few million points, and asked one at a time each point
      pays a boundary crossing for arithmetic that is often shorter than the
      crossing. This walks the candidates in C++ and hands back one array.
      The graph is evaluated exactly as `run()` evaluates it, so the numbers
      are the optimiser's own, not a second implementation of them.

      **Not parallel, and it cannot be.** The candidates are independent but
      the graph is not: they are scored by writing the shared parameter ports,
      so two threads would race over one set of ports. Parallelism here needs
      a graph per thread, which is a different feature (a `GraphNode` deep copy)
      and not this one.

      The parameter ports are restored to what they held on entry, and the
      graph is left evaluated at those values -- scoring a surface must not
      quietly move the fit somebody is holding.
   */
  void compute_objective_batch(double* in_candidates, int n_rows, int n_cols,
                               double** out_view, int* n_out_view);

  // --------------------------------------------------------------- results

  //! The solution, in external (bounded) coordinates.
  std::vector<double> get_x() const;
  //! Residuals at the solution.
  std::vector<double> get_residuals() const;
  //! Sum of squared residuals at the solution.
  double get_chi2() const;
  //! `chi2 / (n_residuals - n_free - 1)`, chisurf's reduced chi-square.
  double get_chi2r() const;
  //! Residual evaluations used.
  int get_number_of_evaluations() const;
  //! MINPACK's `info` from the last run; -1 if the observer cancelled.
  int get_status() const;
  //! Whether the last run stopped because the observer said to.
  bool get_cancelled() const;
  //! MINPACK's message for the last `info`.
  std::string get_message() const;

  //! Parameter covariance at the solution, `n x n` row-major.
  /*!
      `(J^T J)^-1` in **external** coordinates, built by MINPACK's `covar`
      from the `R` the optimiser already computed -- no second Jacobian.

      This is the covariance of a *weighted* least-squares problem, so it is
      the parameter covariance only when the weights are real standard
      deviations. It must be scaled by the residual variance when they are
      not -- the defect okf/log.md 2026-09-01 (7) records, where unweighted
      data reported error bars 22x too large.

      Empty when the run did not converge, or when `R` is singular.
  */
  std::vector<double> get_covariance() const;

  //! Parameter standard errors: the square roots of the covariance diagonal.
  /*! Empty when there is no covariance. */
  std::vector<double> get_errors() const;

  //! The covariance differenced afresh at the solution, over the graph.
  /*!
      `get_covariance()` is free but not always usable: `lmdif` differences
      at a step **relative to the parameter** (`h_j = sqrt(epsfcn)|x_j|`), so
      a parameter that converged near zero is differenced at a step near
      zero, its Jacobian column is round-off, and `covar` reports an enormous
      variance for it. That is what a relative step means and it is not a
      defect -- scipy's `leastsqbound` reports the same thing for the same
      fit -- but it makes the matrix unusable for every real TCSPC model,
      whose free vector runs from a scatter fraction of 1e-5 to a lifetime of
      3.15.

      chisurf's answer was to rebuild the Jacobian in numpy, which costs
      `p + 1` trips through `Model.update_model()` and was **32% of a TCSPC
      fit** -- the largest thing still pulling a fit that optimises entirely
      in C++ back across the boundary. This does the same arithmetic here,
      where the graph and the data already are.

      **The step rule is chisurf's `approx_grad`, deliberately, and not
      `lmdif`'s.** `h_j = epsilon * max(|x_j|, floor)`, rounded to an exactly
      representable difference (`(x + h) - x`) so the divisor is the step the
      model actually saw. The floor is what makes it resolve: a parameter at
      1e-5 is differenced at `epsilon`, not at `epsilon * 1e-5`. The
      optimiser's `epsfcn` is *not* reused and must not be -- it is tuned for
      convergence (chisurf's 1e-6 is worth 23 more converged fits in 88) and
      this is tuned for resolution. Two different jobs, two different steps.

      Columns whose partial derivative is identically zero are **dropped**
      rather than inverted, which `get_covariance_parameters()` then reports.
      A zero column is legitimate: chisurf's `E_FRET` is a free parameter
      whose value comes from a callable, so nothing the optimiser writes
      moves it and the graph gives it a deliberately dangling port. Inverting
      a singular `alpha` instead would hand every FRET fit garbage error
      bars.

      \param[in] epsilon relative step; 0 means `sqrt(machine eps)`, which is
                 chisurf's `FINITE_DIFFERENCE_STEP`
      \param[in] floor absolute floor on the step scale; 0 means 1.0, which
                 is `approx_grad`'s `max(|x|, 1)`
      \return the covariance of the kept parameters, `k x k` row-major, in
              **external** coordinates. Empty before `run()`, or when no
              parameter moves the objective.

      Costs `p + 1` evaluations of the objective: one per parameter, plus one
      that puts the ports back where the solution left them.
  */
  std::vector<double> compute_covariance(double epsilon = 0.0,
                                         double floor = 0.0);

  //! The same covariance, at a parameter vector of the caller's choosing.
  /*!
      compute_covariance() reads the solution and reuses the residuals the
      optimiser left behind; this one takes the point and evaluates `f0`
      itself, so it works without a fit having been run. That is what
      chisurf's `covariance_matrix(fit, model=...)` is -- the curvature at
      wherever the model currently sits, asked for by a posterior view, an
      error-propagation onto a derived quantity, or a sampler looking for a
      preconditioner.

      Costs `p + 2` evaluations rather than `p + 1`: `f0` cannot be assumed.
  */
  std::vector<double> compute_covariance_at(const std::vector<double>& x,
                                            double epsilon = 0.0,
                                            double floor = 0.0);

  //! The forward-difference Jacobian at `x`, `n x m` row-major.
  /*!
      chisurf's `approx_grad`, in C++ and over the graph: one row per
      parameter, one column per residual, at the same step rule
      compute_covariance() uses. Rows are **not** dropped here -- a parameter
      that does not move the objective gets a row of zeros, which is what
      `approx_grad` returns and what its callers filter on.

      Costs `p + 2` evaluations.
  */
  std::vector<double> compute_jacobian(const std::vector<double>& x,
                                       double epsilon = 0.0,
                                       double floor = 0.0);

  //! Which parameters the last compute_covariance() kept, in order.
  /*! Indices into the parameter list. Shorter than it when a parameter does
      not move the objective; see compute_covariance(). */
  std::vector<int> get_covariance_parameters() const;

  //! A one-line summary: algorithm, parameters, status, chi-square.
  std::string describe() const;

  //! Evaluate the objective at an external (bounded) parameter vector.
  /*!
      Writes the parameter ports, updates the graph and returns the residual
      port -- one crossing, whatever the graph behind it costs. Public
      because the covariance is *not* the optimiser's business: a routine
      that differences the objective at its own step rule is written on top
      of this without touching `lmdif` or `epsfcn`, which is what
      compute_covariance() is. It is also the honest way for a test to ask
      "what does the graph say here?" without going through a fit.

      Unlike the loop inside run(), this does not count the evaluation and
      does not report it to the observer -- a progress bar measures the
      minimisation, not everything that ever touched the objective.
  */
  std::vector<double> evaluate_external(const std::vector<double>& xe);

 private:
  //! `epsilon = 0 -> sqrt(machine eps)`, `floor = 0 -> 1.0`, in one place.
  static void covariance_defaults(double* epsilon, double* floor);
  //! The forward-difference loop, given `f0`; restores the ports at the end.
  std::vector<double> jacobian_impl(const std::vector<double>& x,
                                    const std::vector<double>& f0,
                                    double epsilon, double floor);
  //! `pinvh(J J')` over the rows that are not identically zero.
  /*! Sets covariance_parameters_ to the rows it kept. */
  std::vector<double> covariance_from_jacobian(const std::vector<double>& jac,
                                               std::size_t m);
  //! Evaluate at an internal vector, counting it and reporting it.
  bool evaluate_internal(const std::vector<double>& xi,
                         std::vector<double>* fvec);
  void validate() const;
  void configure_from_ports();

  //! chisurf's `_internal2external_lambda`, one coordinate.
  double to_external(double xi, unsigned int i) const;
  //! chisurf's `_external2internal_lambda`, one coordinate.
  double to_internal(double xe, unsigned int i) const;
  //! chisurf's `_internal2external_grad`, one coordinate.
  double external_gradient(double xi, unsigned int i) const;
  //! The forward-difference step in *internal* coordinates for parameter
  //! `i`, given MINPACK's relative-step fraction `eps` (PRD-120).
  double fdjac2_step(double xi, unsigned int i, double eps) const;

  int lmdif(std::vector<double>* x, std::vector<double>* fvec);

  std::string algorithm_ = "leastsq";

  std::vector<std::shared_ptr<GraphPort> > parameters_;
  std::shared_ptr<GraphNode> objective_node_;
  std::string residual_key_ = "residuals";
  std::function<std::vector<double>(const std::vector<double>&)>
      residual_function_;
  //! Owned, because it outlives the call that set it: a Python observer
  //! held only by a raw pointer would be collected between iterations.
  IMP::Pointer<FitMinimizerObserver> observer_;

  std::vector<double> initial_values_;
  std::vector<double> lower_;
  std::vector<double> upper_;
  std::vector<double> diag_;
  unsigned int ndim_ = 0;

  double ftol_ = 1.49012e-8;
  double xtol_ = 1.49012e-8;
  double gtol_ = 0.0;
  double epsfcn_ = 0.0;
  double factor_ = 100.0;
  int maxfev_ = 0;

  std::vector<double> x_;
  std::vector<double> fvec_;
  //! The final `R`, `n x n` row-major, already in external coordinates.
  std::vector<double> r_;
  //! `covar`'s column permutation, 0-based (MINPACK's `ipvt` less one).
  std::vector<int> ipvt_;
  std::vector<double> covariance_;
  //! Which parameters the last compute_covariance() kept.
  std::vector<int> covariance_parameters_;
  double chi2_ = 0.0;
  int n_evaluations_ = 0;
  int expected_evaluations_ = 1;
  int status_ = 0;
  bool cancelled_ = false;
};

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_FITMINIMIZER_H
