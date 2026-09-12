/**
 *  \file IMP/bff/Sampler.cpp
 *  \brief ChiSurf's MCMC samplers over a bff GraphPort/GraphNode model, in C++.
 *
 *  The implementation notes below name the chisurf function each piece is
 *  ported from; the algorithms are those functions, not re-derivations of
 *  them. Numeric helpers (Cholesky, a Jacobi eigensolver for the walker
 *  degeneracy test, a JSON object reader for the port prior specs) live in
 *  the sampler_detail namespace because this file is compiled into IMP's
 *  unity build, where two anonymous namespaces are the same namespace and
 *  same-named helpers would collide.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#include <IMP/bff/Sampler.h>

#include <IMP/bff/FactorGraph.h>
#include <IMP/bff/GraphNode.h>
#include <IMP/bff/GraphPort.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <sstream>
#include <utility>

IMPBFF_BEGIN_NAMESPACE

namespace sampler_detail {

//! A dense matrix, flat row-major (kept out of the header on purpose).
typedef std::vector<std::vector<double> > Matrix;

//! Lower-triangular Cholesky factor of A, or an empty vector on failure.
Matrix cholesky(const Matrix& a) {
  const std::size_t n = a.size();
  Matrix l(n, std::vector<double>(n, 0.0));
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      double sum = a[i][j];
      for (std::size_t k = 0; k < j; ++k) sum -= l[i][k] * l[j][k];
      if (i == j) {
        if (!(sum > 0.0)) return Matrix();
        l[i][j] = std::sqrt(sum);
      } else {
        l[i][j] = sum / l[j][j];
      }
    }
  }
  return l;
}

//! chisurf's _cholesky_or_diagonal: ridge, then drop the correlations.
Matrix cholesky_or_diagonal(const Matrix& cov) {
  const std::size_t k = cov.size();
  double scale = 0.0;
  for (std::size_t i = 0; i < k; ++i) scale += cov[i][i];
  scale /= (k > 0 ? static_cast<double>(k) : 1.0);
  const double ridges[4] = {0.0, 1e-10, 1e-6, 1e-3};
  for (int r = 0; r < 4; ++r) {
    Matrix ridged = cov;
    for (std::size_t i = 0; i < k; ++i) ridged[i][i] += ridges[r] * scale;
    Matrix factor = cholesky(ridged);
    if (!factor.empty()) return factor;
  }
  Matrix diag(k, std::vector<double>(k, 0.0));
  for (std::size_t i = 0; i < k; ++i)
    diag[i][i] = std::sqrt(std::max(cov[i][i], 1e-30));
  return diag;
}

//! chisurf's _regularised_covariance: a shrunk sample covariance.
/*!  Returns an empty matrix when the window is unusable (too short, not
     finite, or a zero variance), as chisurf returns None. */
Matrix regularised_covariance(const std::vector<const double*>& draws,
                              std::size_t n_draws, const std::vector<int>& idx) {
  const std::size_t k = idx.size();
  if (n_draws < 3 || k == 0) return Matrix();
  Matrix cov(k, std::vector<double>(k, 0.0));
  std::vector<double> mean(k, 0.0);
  for (std::size_t d = 0; d < n_draws; ++d)
    for (std::size_t j = 0; j < k; ++j)
      mean[j] += draws[d][static_cast<std::size_t>(idx[j])];
  for (std::size_t j = 0; j < k; ++j) mean[j] /= static_cast<double>(n_draws);
  for (std::size_t d = 0; d < n_draws; ++d)
    for (std::size_t i = 0; i < k; ++i)
      for (std::size_t j = 0; j < k; ++j)
        cov[i][j] += (draws[d][static_cast<std::size_t>(idx[i])] - mean[i]) *
                     (draws[d][static_cast<std::size_t>(idx[j])] - mean[j]);
  const double denom = static_cast<double>(n_draws) - 1.0;
  for (std::size_t i = 0; i < k; ++i)
    for (std::size_t j = 0; j < k; ++j) cov[i][j] /= denom;
  for (std::size_t i = 0; i < k; ++i)
    for (std::size_t j = 0; j < k; ++j)
      if (!std::isfinite(cov[i][j])) return Matrix();
  for (std::size_t i = 0; i < k; ++i)
    if (!(cov[i][i] > 0.0)) return Matrix();
  const double weight = static_cast<double>(n_draws) /
                        (static_cast<double>(n_draws) + 5.0);
  for (std::size_t i = 0; i < k; ++i)
    for (std::size_t j = 0; j < k; ++j)
      cov[i][j] = weight * cov[i][j] + (1.0 - weight) * (i == j ? cov[i][i] : 0.0);
  return cov;
}

//! A lower-triangular matrix flattened row-major (BlockState::factor's form).
std::vector<double> flatten_factor(const Matrix& m) {
  std::vector<double> flat;
  flat.reserve(m.size() * m.size());
  for (std::size_t i = 0; i < m.size(); ++i)
    for (std::size_t j = 0; j < m[i].size(); ++j) flat.push_back(m[i][j]);
  return flat;
}

//! k distinct entries of a population, uniformly (a partial Fisher-Yates,
//! standing in for numpy's rng.choice(..., replace=False)).
std::vector<int> sample_without_replacement(const std::vector<int>& population,
                                            std::size_t k, std::mt19937_64& rng) {
  std::vector<int> pool = population;
  if (k > pool.size()) k = pool.size();
  for (std::size_t i = 0; i < k; ++i) {
    std::uniform_int_distribution<std::size_t> pick(i, pool.size() - 1);
    const std::size_t j = pick(rng);
    std::swap(pool[i], pool[j]);
  }
  pool.resize(k);
  return pool;
}

//! chisurf's _adaptation_windows: Stan's init/doubling/term schedule.
void adaptation_windows(int n_adapt, int& init_buffer, int& term_buffer,
                        std::vector<int>& window_ends) {
  init_buffer = 75;
  term_buffer = 50;
  int base_window = 25;
  window_ends.clear();
  n_adapt = int(n_adapt);
  if (n_adapt < 20) {
    init_buffer = n_adapt;
    term_buffer = 0;
    return;
  }
  if (init_buffer + base_window + term_buffer > n_adapt) {
    init_buffer = int(std::lround(0.15 * n_adapt));
    term_buffer = int(std::lround(0.10 * n_adapt));
    base_window = n_adapt - init_buffer - term_buffer;
    if (base_window < 2) {
      init_buffer = n_adapt;
      term_buffer = 0;
      return;
    }
  }
  int start = init_buffer, window = base_window;
  const int last = n_adapt - term_buffer;
  while (start + window <= last) {
    int end = start + window;
    if (end + 2 * window > last) end = last;
    window_ends.push_back(end);
    start = end;
    window *= 2;
  }
}

//! Eigenvalues of a symmetric matrix by cyclic Jacobi rotations.
void symmetric_eigenvalues(Matrix a, std::vector<double>& values) {
  const std::size_t n = a.size();
  values.assign(n, 0.0);
  if (n == 0) return;
  for (int sweep = 0; sweep < 100; ++sweep) {
    double off = 0.0;
    for (std::size_t p = 0; p < n; ++p)
      for (std::size_t q = p + 1; q < n; ++q) off += a[p][q] * a[p][q];
    if (off <= 1e-300) break;
    for (std::size_t p = 0; p < n; ++p) {
      for (std::size_t q = p + 1; q < n; ++q) {
        if (std::fabs(a[p][q]) <= 1e-300) continue;
        const double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
        const double t =
            (theta >= 0.0 ? 1.0 : -1.0) /
            (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
        const double c = 1.0 / std::sqrt(t * t + 1.0);
        const double s = t * c;
        for (std::size_t k = 0; k < n; ++k) {
          const double akp = a[k][p], akq = a[k][q];
          a[k][p] = c * akp - s * akq;
          a[k][q] = s * akp + c * akq;
        }
        for (std::size_t k = 0; k < n; ++k) {
          const double apk = a[p][k], aqk = a[q][k];
          a[p][k] = c * apk - s * aqk;
          a[q][k] = s * apk + c * aqk;
        }
      }
    }
  }
  for (std::size_t i = 0; i < n; ++i) values[i] = a[i][i];
}

//! Skip whitespace; \return the character at the position (0 at the end).
char json_peek(const std::string& s, std::size_t& i) {
  while (i < s.size() && (s[i] == ' ' || s[i] == '\t' || s[i] == '\n' ||
                          s[i] == '\r'))
    ++i;
  return i < s.size() ? s[i] : '\0';
}

//! Read a JSON string (the opening quote is at i); \return the contents.
std::string json_string(const std::string& s, std::size_t& i) {
  std::string out;
  ++i;  // the opening quote
  while (i < s.size() && s[i] != '"') {
    if (s[i] == '\\' && i + 1 < s.size()) {
      ++i;
      if (s[i] == 'n') {
        out += '\n';
      } else if (s[i] == 't') {
        out += '\t';
      } else {
        out += s[i];
      }
    } else {
      out += s[i];
    }
    ++i;
  }
  ++i;  // the closing quote
  return out;
}

}  // namespace sampler_detail

// ---------------------------------------------------------------- lifecycle

Sampler::Sampler(const std::string& algorithm, unsigned int seed)
    : algorithm_("stretch"), seed_(seed), rng_(seed) {
  set_algorithm(algorithm);
}

Sampler::~Sampler() {}

void Sampler::set_algorithm(const std::string& algorithm) {
  std::string canonical;
  if (algorithm == "stretch" || algorithm == "ensemble" ||
      algorithm == "sample_ensemble" || algorithm == "affine" ||
      algorithm == "emcee") {
    canonical = "stretch";
  } else if (algorithm == "de" || algorithm == "differential_evolution" ||
             algorithm == "sample_differential_evolution") {
    canonical = "de";
  } else if (algorithm == "slice" || algorithm == "zeus" ||
             algorithm == "ensemble_slice" || algorithm == "sample_slice") {
    canonical = "slice";
  } else if (algorithm == "metropolis" || algorithm == "walk" ||
             algorithm == "walk_mcmc" || algorithm == "blocked" ||
             algorithm == "walk_mcmc_blocked") {
    canonical = "metropolis";
  } else {
    throw SamplerConfigurationError(
        "unknown sampler algorithm '" + algorithm +
        "'; expected 'stretch', 'slice', 'de' or 'metropolis'");
  }
  if (canonical != algorithm_) {
    algorithm_ = canonical;
    initialized_ = false;  // the ensemble belongs to the old algorithm
  }
}

const std::string& Sampler::get_algorithm() const { return algorithm_; }

void Sampler::set_seed(unsigned int seed) {
  seed_ = seed;
  rng_.seed(seed);
}

unsigned int Sampler::get_seed() const { return seed_; }

// --------------------------------------------------------------- parameters

void Sampler::set_parameter_ports(
    const std::vector<std::shared_ptr<GraphPort> >& parameters) {
  for (std::size_t i = 0; i < parameters.size(); ++i) {
    if (!parameters[i])
      throw SamplerConfigurationError("parameter " + std::to_string(i) +
                               " is a null port");
    if (parameters[i]->get_fixed())
      throw SamplerConfigurationError(
          "parameter '" + parameters[i]->get_name() +
          "' is fixed; a fixed port cannot be sampled (chisurf samples "
          "the free parameters)");
  }
  parameters_ = parameters;
  ndim_ = static_cast<unsigned int>(parameters.size());
  initialized_ = false;
}

std::vector<std::shared_ptr<GraphPort> > Sampler::get_parameter_ports() const {
  return parameters_;
}

std::vector<std::string> Sampler::get_parameter_names() const {
  if (!names_.empty()) return names_;
  std::vector<std::string> out;
  for (std::size_t i = 0; i < parameters_.size(); ++i)
    out.push_back(parameters_[i]->get_name().empty()
                      ? "x" + std::to_string(i)
                      : parameters_[i]->get_name());
  return out;
}

void Sampler::set_initial_values(const std::vector<double>& values) {
  initial_values_ = values;
  if (ndim_ == 0) ndim_ = static_cast<unsigned int>(values.size());
  initialized_ = false;
}

std::vector<double> Sampler::get_initial_values() const { return initial_values_; }

void Sampler::set_bounds(const std::vector<double>& lower,
                         const std::vector<double>& upper) {
  if (lower.size() != upper.size())
    throw SamplerConfigurationError(
        "set_bounds: lower and upper must have the same length");
  if (ndim_ != 0 && lower.size() != ndim_)
    throw SamplerConfigurationError("set_bounds: expected " +
                             std::to_string(ndim_) + " bounds, got " +
                             std::to_string(lower.size()));
  lower_ = lower;
  upper_ = upper;
  bounds_explicit_ = true;
  initialized_ = false;
}

// ---------------------------------------------------------------- objective

void Sampler::set_objective(std::shared_ptr<GraphNode> node,
                            const std::string& output_port) {
  if (!node)
    throw SamplerConfigurationError("set_objective: the node is a null pointer");
  if (!node->get_output_port(output_port))
    throw SamplerConfigurationError(
        "set_objective: node '" + node->get_name() +
        "' has no output port '" + output_port + "'");
  objective_node_ = node;
  output_port_ = node->get_output_port(output_port);
  objective_function_ = nullptr;
  initialized_ = false;
}

std::shared_ptr<GraphNode> Sampler::get_objective() const { return objective_node_; }

void Sampler::set_output_is_log_likelihood(bool v) {
  output_is_log_likelihood_ = v;
}

bool Sampler::get_output_is_log_likelihood() const {
  return output_is_log_likelihood_;
}

void Sampler::set_objective_function(
    std::function<double(const std::vector<double>&)> objective) {
  if (!objective)
    throw SamplerConfigurationError("set_objective_function: a null function");
  objective_function_ = objective;
  objective_node_.reset();
  output_port_.reset();
  initialized_ = false;
}

bool Sampler::has_objective() const {
  return objective_node_ != nullptr || objective_function_ != nullptr;
}

// ----------------------------------------------------------------- blocking

void Sampler::set_factor_graph(FactorGraph* graph) {
  factor_graph_ = graph;
  initialized_ = false;
}

FactorGraph* Sampler::get_factor_graph() const {
  return factor_graph_;
}

void Sampler::set_blocks(const std::vector<int>& flat_indices,
                         const std::vector<int>& block_sizes) {
  explicit_blocks_.clear();
  std::size_t offset = 0;
  for (std::size_t b = 0; b < block_sizes.size(); ++b) {
    const int size = block_sizes[b];
    if (size <= 0)
      throw SamplerConfigurationError("set_blocks: block sizes must be positive");
    if (offset + static_cast<std::size_t>(size) > flat_indices.size())
      throw SamplerConfigurationError("set_blocks: more block entries than indices");
    std::vector<int> block(
        flat_indices.begin() + offset,
        flat_indices.begin() + offset + static_cast<std::size_t>(size));
    for (std::size_t i = 0; i < block.size(); ++i) {
      if (ndim_ != 0 &&
          (block[i] < 0 ||
           block[i] >= static_cast<int>(ndim_)))
        throw SamplerConfigurationError("set_blocks: index out of range");
    }
    explicit_blocks_.push_back(block);
    offset += static_cast<std::size_t>(size);
  }
  initialized_ = false;
}

std::vector<std::vector<int> > Sampler::get_blocks() const {
  std::vector<std::vector<int> > out;
  for (std::size_t b = 0; b < blocks_.size(); ++b)
    out.push_back(blocks_[b].indices);
  return out;
}

std::vector<int> Sampler::get_block_sizes() const {
  std::vector<int> out;
  for (std::size_t b = 0; b < blocks_.size(); ++b)
    out.push_back(static_cast<int>(blocks_[b].indices.size()));
  return out;
}

// ---------------------------------------------------------------- tunables

void Sampler::set_number_of_walkers(int n) {
  n_walkers_setting_ = n;
  initialized_ = false;
}

int Sampler::get_number_of_walkers() const {
  if (algorithm_ != "stretch" && algorithm_ != "slice") return 1;
  if (n_walkers_setting_ > 0) return n_walkers_setting_;
  const int d = static_cast<int>(ndim_);
  return std::max(2 * d + 2, 10);
}

void Sampler::set_stretch_scale(double a) {
  if (!(a > 1.0))
    throw SamplerConfigurationError("the stretch scale must exceed one");
  stretch_scale_ = a;
}

double Sampler::get_stretch_scale() const { return stretch_scale_; }

double Sampler::get_slice_mu() const { return slice_mu_; }

void Sampler::set_slice_mu(double mu) {
  if (!(mu > 0.0))
    throw SamplerConfigurationError(
        "the slice direction scale must be positive");
  slice_mu_ = mu;
  // An explicit scale is a decision; tuning would overwrite it.
  slice_tuning_ = false;
}

void Sampler::set_slice_max_steps(int n) {
  slice_max_steps_ = n > 0 ? n : 10000;
}

int Sampler::get_slice_max_steps() const { return slice_max_steps_; }

long Sampler::get_slice_truncations() const { return slice_truncations_; }

void Sampler::set_live_dangerously(bool v) { live_dangerously_ = v; }

bool Sampler::get_live_dangerously() const { return live_dangerously_; }

void Sampler::set_number_of_chains(int n) {
  n_chains_setting_ = n;
  initialized_ = false;
}

int Sampler::get_number_of_chains() const {
  if (algorithm_ != "de") return 1;
  int n = n_chains_setting_ > 0
              ? n_chains_setting_
              : std::max(8, 2 * static_cast<int>(ndim_));
  return std::max(4, n);
}

void Sampler::set_jitter(double jitter) {
  if (!(jitter >= 0.0))
    throw SamplerConfigurationError("the jitter must not be negative");
  jitter_ = jitter;
}

double Sampler::get_jitter() const { return jitter_; }

void Sampler::set_snooker(double fraction) {
  if (!(fraction >= 0.0 && fraction <= 1.0))
    throw SamplerConfigurationError("the snooker fraction must lie in [0, 1]");
  snooker_ = fraction;
}

double Sampler::get_snooker() const { return snooker_; }

void Sampler::set_step_size(double step_size) {
  if (!(step_size > 0.0))
    throw SamplerConfigurationError("the step size must be positive");
  step_size_ = step_size;
}

double Sampler::get_step_size() const { return step_size_; }

void Sampler::set_proposal_covariance(
    const std::vector<std::vector<double> >& cov) {
  if (!cov.empty()) {
    const std::size_t n = cov.size();
    for (const std::vector<double>& row : cov) {
      if (row.size() != n) {
        throw SamplerConfigurationError(
            "set_proposal_covariance: the covariance must be square");
      }
    }
  }
  proposal_covariance_ = cov;
}

std::vector<std::vector<double> > Sampler::get_proposal_covariance() const {
  return proposal_covariance_;
}

void Sampler::set_n_adapt(int n) { n_adapt_setting_ = n; }

int Sampler::get_n_adapt() const { return n_adapt_setting_; }

void Sampler::set_temp(double temp) {
  if (!(temp > 0.0))
    throw SamplerConfigurationError("the temperature must be positive");
  temp_ = temp;
}

double Sampler::get_temp() const { return temp_; }

void Sampler::set_chi2max(double chi2max) { chi2max_ = chi2max; }

double Sampler::get_chi2max() const { return chi2max_; }

void Sampler::set_walker_start_std(double std) {
  if (!(std > 0.0))
    throw SamplerConfigurationError("the walker start spread must be positive");
  walker_start_std_ = std;
}

double Sampler::get_walker_start_std() const { return walker_start_std_; }

void Sampler::set_walker_start(
    const std::vector<std::vector<double> >& start) {
  walker_start_override_ = start;
  initialized_ = false;
}

std::vector<std::vector<double> > Sampler::get_walker_start() const {
  return walker_start_override_;
}

void Sampler::set_observer(std::function<void(int, int)> observer) {
  observer_ = observer;
}

// --------------------------------------------------------------- internals

double Sampler::PriorSpec::get(const std::string& name,
                               double fallback) const {
  for (std::size_t i = 0; i < numbers.size(); ++i)
    if (numbers[i].first == name) return numbers[i].second;
  return fallback;
}

//! chisurf priors.py, kind for kind; unknown kinds contribute nothing.
double Sampler::prior_lnpdf(const PriorSpec& spec, double x) {
  const double inf = std::numeric_limits<double>::infinity();
  const double k2pi = 6.283185307179586476925286766559;
  const double half_ln_2pi = 0.918938533204672741780329736406;
  if (spec.kind == "uniform") {
    const double lb = spec.get("lb", -inf), ub = spec.get("ub", inf);
    if (x < lb || x > ub) return -inf;
    const double width = ub - lb;
    if (std::isfinite(width) && width > 0.0) return -std::log(width);
    return 0.0;
  }
  if (spec.kind == "normal") {
    const double mu = spec.get("mu", 0.0), sigma = spec.get("sigma", 1.0);
    const double z = (x - mu) / sigma;
    return -0.5 * z * z - std::log(sigma) - half_ln_2pi;
  }
  if (spec.kind == "truncated_normal") {
    const double lb = spec.get("lb", -inf), ub = spec.get("ub", inf);
    if (x < lb || x > ub) return -inf;
    const double mu = spec.get("mu", 0.0), sigma = spec.get("sigma", 1.0);
    const double z = (x - mu) / sigma;
    return -0.5 * z * z - std::log(sigma) - half_ln_2pi;
  }
  if (spec.kind == "half_normal") {
    const double sigma = spec.get("sigma", 1.0), loc = spec.get("loc", 0.0);
    if (x < loc) return -inf;
    const double z = (x - loc) / sigma;
    return -0.5 * z * z - std::log(sigma) + 0.5 * std::log(2.0 / M_PI);
  }
  if (spec.kind == "lognormal") {
    const double mu = spec.get("mu", 0.0), sigma = spec.get("sigma", 1.0);
    if (x <= 0.0) return -inf;
    const double lx = std::log(x);
    const double z = (lx - mu) / sigma;
    return -0.5 * z * z - lx - std::log(sigma) - half_ln_2pi;
  }
  if (spec.kind == "exponential") {
    const double scale = spec.get("scale", 1.0), loc = spec.get("loc", 0.0);
    if (x < loc) return -inf;
    return -(x - loc) / scale - std::log(scale);
  }
  if (spec.kind == "gamma") {
    const double alpha = spec.get("alpha", 1.0), beta = spec.get("beta", 1.0),
                 loc = spec.get("loc", 0.0);
    const double t = x - loc;
    if (t <= 0.0) return -inf;
    return (alpha - 1.0) * std::log(t) - beta * t + alpha * std::log(beta) -
           std::lgamma(alpha);
  }
  if (spec.kind == "beta") {
    const double alpha = spec.get("alpha", 1.0), beta = spec.get("beta", 1.0);
    if (x <= 0.0 || x >= 1.0) return -inf;
    const double log_b =
        std::lgamma(alpha) + std::lgamma(beta) - std::lgamma(alpha + beta);
    return (alpha - 1.0) * std::log(x) + (beta - 1.0) * std::log1p(-x) - log_b;
  }
  return 0.0;  // an unrecognised kind is no prior (prior_from_state -> None)
}

Sampler::PriorSpec Sampler::parse_prior(const std::string& json) {
  PriorSpec spec;
  if (json.empty()) return spec;
  std::size_t i = 0;
  if (sampler_detail::json_peek(json, i) != '{') return spec;
  ++i;
  while (true) {
    char c = sampler_detail::json_peek(json, i);
    if (c == '\0') return PriorSpec();  // malformed: no prior
    if (c == '}') return spec;
    if (c == ',') {
      ++i;
      continue;
    }
    if (c != '"') return PriorSpec();
    const std::string key = sampler_detail::json_string(json, i);
    if (sampler_detail::json_peek(json, i) != ':') return PriorSpec();
    ++i;
    c = sampler_detail::json_peek(json, i);
    if (c == '"') {
      const std::string value = sampler_detail::json_string(json, i);
      if (key == "kind") spec.kind = value;
    } else if (c == 'n' || c == 't' || c == 'f') {
      // null / true / false: skipped, as nothing numeric is keyed there
      while (i < json.size() && json[i] != ',' && json[i] != '}') ++i;
    } else {
      const char* begin = json.c_str() + i;
      char* end = nullptr;
      const double value = std::strtod(begin, &end);
      if (end == begin) return PriorSpec();
      i += static_cast<std::size_t>(end - begin);
      spec.numbers.push_back(std::make_pair(key, value));
    }
  }
}

void Sampler::configure_from_ports() {
  if (parameters_.empty()) {
    // Plain-vector mode: no ports, no priors, bounds only when explicit.
    if (!bounds_explicit_) {
      lower_.assign(ndim_, -std::numeric_limits<double>::infinity());
      upper_.assign(ndim_, std::numeric_limits<double>::infinity());
    }
    priors_.assign(ndim_, PriorSpec());
    names_.clear();
    for (unsigned int i = 0; i < ndim_; ++i)
      names_.push_back("x" + std::to_string(i));
    if (initial_values_.empty()) initial_values_.assign(ndim_, 0.0);
    return;
  }
  if (initial_values_.empty()) {
    initial_values_.resize(ndim_);
    for (unsigned int i = 0; i < ndim_; ++i)
      initial_values_[i] = parameters_[i]->get_value();
  }
  if (!bounds_explicit_) {
    // A port that enforces bounds gives its box; anything else is
    // unbounded (NaN stored on a port stands for no bound, as (nan, nan)
    // is what the Python surface reports for enforcement off).
    lower_.resize(ndim_);
    upper_.resize(ndim_);
    for (unsigned int i = 0; i < ndim_; ++i) {
      const std::shared_ptr<GraphPort>& p = parameters_[i];
      double lb = -std::numeric_limits<double>::infinity();
      double ub = std::numeric_limits<double>::infinity();
      if (p->get_is_bounded()) {
        if (std::isnan(p->get_lower_bound()))
          lb = -std::numeric_limits<double>::infinity();
        else
          lb = p->get_lower_bound();
        if (std::isnan(p->get_upper_bound()))
          ub = std::numeric_limits<double>::infinity();
        else
          ub = p->get_upper_bound();
      }
      lower_[i] = lb;
      upper_[i] = ub;
    }
  }
  priors_.resize(ndim_);
  for (unsigned int i = 0; i < ndim_; ++i)
    priors_[i] = parse_prior(parameters_[i]->get_prior());
  names_.clear();
  for (unsigned int i = 0; i < ndim_; ++i) {
    const std::string& n = parameters_[i]->get_name();
    names_.push_back(n.empty() ? "x" + std::to_string(i) : n);
  }
}

//! chisurf fit.lnprior with explicit bounds: the box first, then the priors.
double Sampler::log_prior(const std::vector<double>& x) const {
  const double inf = std::numeric_limits<double>::infinity();
  for (unsigned int i = 0; i < ndim_; ++i) {
    if (x[i] < lower_[i] || x[i] > upper_[i]) return -inf;
  }
  double lp = 0.0;
  for (unsigned int i = 0; i < ndim_; ++i) {
    if (priors_[i].empty()) continue;
    lp += prior_lnpdf(priors_[i], x[i]);
    if (!std::isfinite(lp)) return -inf;
  }
  return lp;
}

//! chisurf fit.lnprob_parts: the prior short-circuits the evaluation.
Sampler::Parts Sampler::evaluate(const std::vector<double>& x) {
  ++n_evaluations_;
  const double inf = std::numeric_limits<double>::infinity();
  Parts parts;
  parts.lnprior = log_prior(x);
  if (!std::isfinite(parts.lnprior)) {
    parts.lnpost = -inf;
    parts.chi2 = inf;
    return parts;
  }
  double lnlike;
  if (objective_function_) {
    lnlike = objective_function_(x);
    if (std::isnan(lnlike))
      throw SamplerConfigurationError("the log-probability returned NaN");
    parts.chi2 = -2.0 * lnlike;
  } else {
    for (unsigned int i = 0; i < ndim_; ++i)
      parameters_[i]->set_value(x[i]);
    objective_node_->update();
    const double value = output_port_->get_value();
    if (std::isnan(value))
      throw SamplerConfigurationError("the log-probability returned NaN");
    if (output_is_log_likelihood_) {
      lnlike = value;
      parts.chi2 = -2.0 * value;
    } else {
      parts.chi2 = value;
      lnlike = (value < chi2max_) ? -0.5 * value : -inf;
    }
  }
  parts.lnpost = lnlike + parts.lnprior;
  if (std::isnan(parts.lnpost)) {
    // -inf likelihood plus a -inf prior is a rejection; NaN is chisurf's
    // hard error. A finite prior cannot make -inf finite again, so this
    // branch is only ever the function objective's doing.
    parts.lnpost = -inf;
  }
  return parts;
}

std::vector<std::vector<double> > Sampler::spread_walkers(int n) const {
  // chisurf's _ensemble_walker_start: the bounded range where there is
  // one, the value (floored by the absolute std) where there is not, so
  // no direction is ever left without spread.
  std::vector<double> spread(ndim_);
  for (unsigned int i = 0; i < ndim_; ++i) {
    const bool lo = std::isfinite(lower_[i]), hi = std::isfinite(upper_[i]);
    if (lo && hi) {
      spread[i] = (upper_[i] - lower_[i]) * 1e-4;
    } else {
      spread[i] = std::fabs(initial_values_[i]) > 1e-15
                      ? std::fabs(initial_values_[i]) * walker_start_std_
                      : walker_start_std_;
    }
    if (!(spread[i] > 0.0)) spread[i] = walker_start_std_;
  }
  std::vector<std::vector<double> > start(
      static_cast<std::size_t>(n), std::vector<double>(ndim_));
  std::normal_distribution<double> normal(0.0, 1.0);
  for (int w = 0; w < n; ++w)
    for (unsigned int i = 0; i < ndim_; ++i) {
      double v = initial_values_[i] + spread[i] * normal(rng_);
      start[static_cast<std::size_t>(w)][i] =
          std::min(std::max(v, lower_[i]), upper_[i]);
    }
  return start;
}

//! chisurf's walkers_independent, through the Gram matrix eigenvalues.
bool Sampler::walkers_independent(
    const std::vector<std::vector<double> >& coords) const {
  const std::size_t n = coords.size();
  if (n < 2 || ndim_ == 0) return false;
  std::vector<double> mean(ndim_, 0.0);
  for (std::size_t w = 0; w < n; ++w)
    for (unsigned int i = 0; i < ndim_; ++i) {
      if (!std::isfinite(coords[w][i])) return false;
      mean[i] += coords[w][i];
    }
  for (unsigned int i = 0; i < ndim_; ++i) mean[i] /= static_cast<double>(n);
  sampler_detail::Matrix c(
      n, std::vector<double>(ndim_, 0.0));
  std::vector<double> col_max(ndim_, 0.0);
  for (std::size_t w = 0; w < n; ++w)
    for (unsigned int i = 0; i < ndim_; ++i) {
      c[w][i] = coords[w][i] - mean[i];
      col_max[i] = std::max(col_max[i], std::fabs(c[w][i]));
    }
  for (unsigned int i = 0; i < ndim_; ++i)
    if (!(col_max[i] > 0.0)) return false;
  std::vector<double> col_norm(ndim_, 0.0);
  for (std::size_t w = 0; w < n; ++w)
    for (unsigned int i = 0; i < ndim_; ++i) {
      c[w][i] /= col_max[i];
      col_norm[i] += c[w][i] * c[w][i];
    }
  for (unsigned int i = 0; i < ndim_; ++i) {
    if (!(col_norm[i] > 0.0)) return false;
    col_norm[i] = std::sqrt(col_norm[i]);
  }
  for (std::size_t w = 0; w < n; ++w)
    for (unsigned int i = 0; i < ndim_; ++i) c[w][i] /= col_norm[i];
  // cond(c) = sqrt(largest/smallest eigenvalue of c^T c)
  sampler_detail::Matrix gram(ndim_, std::vector<double>(ndim_, 0.0));
  for (std::size_t w = 0; w < n; ++w)
    for (unsigned int i = 0; i < ndim_; ++i)
      for (unsigned int j = 0; j < ndim_; ++j)
        gram[i][j] += c[w][i] * c[w][j];
  std::vector<double> values;
  sampler_detail::symmetric_eigenvalues(gram, values);
  double lmin = values[0], lmax = values[0];
  for (std::size_t i = 1; i < values.size(); ++i) {
    lmin = std::min(lmin, values[i]);
    lmax = std::max(lmax, values[i]);
  }
  if (!(lmin > 0.0)) return false;
  return std::sqrt(lmax / lmin) <= 1e8;
}

void Sampler::validate() const {
  if (ndim_ == 0)
    throw SamplerConfigurationError(
        "no parameters: call set_parameter_ports() or set_initial_values()");
  if (!has_objective())
    throw SamplerConfigurationError(
        "no objective: call set_objective() or set_objective_function()");
  if (!initial_values_.empty() && initial_values_.size() != ndim_)
    throw SamplerConfigurationError("the initial values do not match the parameters");
  if (bounds_explicit_ &&
      (lower_.size() != ndim_ || upper_.size() != ndim_))
    throw SamplerConfigurationError("the bounds do not match the parameters");
}

void Sampler::rebuild_blocks() {
  blocks_.clear();
  std::vector<std::vector<int> > partition;
  if (!explicit_blocks_.empty()) {
    partition = explicit_blocks_;
  } else if (factor_graph_) {
    // chisurf's _default_blocks: the graph's partition, a parameter the
    // graph cannot place dropped, and one block over everything when the
    // cover misses anything.
    const std::vector<std::vector<std::string> > graph_blocks =
        factor_graph_->get_sampling_blocks();
    std::vector<int> covered;
    for (std::size_t b = 0; b < graph_blocks.size(); ++b) {
      std::vector<int> block;
      for (std::size_t k = 0; k < graph_blocks[b].size(); ++k) {
        const int idx = factor_graph_->index_of(graph_blocks[b][k]);
        if (idx >= 0 && idx < static_cast<int>(ndim_)) block.push_back(idx);
      }
      if (block.empty()) continue;
      std::sort(block.begin(), block.end());
      partition.push_back(block);
      covered.insert(covered.end(), block.begin(), block.end());
    }
    std::sort(covered.begin(), covered.end());
    bool full = covered.size() == ndim_;
    for (unsigned int i = 0; full && i < ndim_; ++i)
      full = covered[i] == static_cast<int>(i);
    if (!full) {
      partition.clear();
      std::vector<int> all(ndim_);
      for (unsigned int i = 0; i < ndim_; ++i)
        all[i] = static_cast<int>(i);
      partition.push_back(all);
    }
  } else {
    std::vector<int> all(ndim_);
    for (unsigned int i = 0; i < ndim_; ++i) all[i] = static_cast<int>(i);
    partition.push_back(all);
  }
  for (std::size_t b = 0; b < partition.size(); ++b) {
    if (partition[b].empty()) continue;
    BlockState state;
    state.indices = partition[b];
    blocks_.push_back(state);
  }
}

//! chisurf's _seed_block_covariances: the caller's curvature where it is
//! finite and positive definite for a block, the diagonal otherwise.
void Sampler::seed_blocks() {
  const bool have_full =
      proposal_covariance_.size() == static_cast<std::size_t>(ndim_);
  for (std::size_t b = 0; b < blocks_.size(); ++b) {
    BlockState& block = blocks_[b];
    const std::size_t k = block.indices.size();
    block.from_curvature = false;
    sampler_detail::Matrix cov;
    if (have_full) {
      // The block's submatrix of the caller's full covariance.
      sampler_detail::Matrix candidate(k, std::vector<double>(k, 0.0));
      bool usable = true;
      for (std::size_t i = 0; i < k && usable; ++i) {
        for (std::size_t j = 0; j < k; ++j) {
          const double v = proposal_covariance_
              [static_cast<std::size_t>(block.indices[i])]
              [static_cast<std::size_t>(block.indices[j])];
          if (!std::isfinite(v)) { usable = false; break; }
          candidate[i][j] = v;
        }
        if (usable && !(candidate[i][i] > 0.0)) usable = false;
      }
      if (usable && !sampler_detail::cholesky(candidate).empty()) {
        cov = candidate;
        block.from_curvature = true;
      }
    }
    if (cov.empty()) {
      cov.assign(k, std::vector<double>(k, 0.0));
      for (std::size_t i = 0; i < k; ++i) {
        double scale =
            std::fabs(initial_values_[static_cast<std::size_t>(
                          block.indices[i])]) *
            step_size_;
        if (scale < 1e-15) scale = step_size_;
        cov[i][i] = scale * scale;
      }
    }
    block.factor = sampler_detail::flatten_factor(sampler_detail::cholesky_or_diagonal(cov));
    // chisurf's per-block target: 0.44 for a singleton, 0.44/sqrt(k)
    // (floored at 0.234) otherwise, and the theoretically optimal start
    // 2.38/sqrt(k) for the log scale -- OPTIMAL_RWM_SCALING.
    const double target = k == 1
                              ? 0.44
                              : std::max(0.234,
                                         0.44 / std::sqrt(static_cast<double>(k)));
    const double log_scale =
        std::log(2.38 / std::sqrt(static_cast<double>(std::max<std::size_t>(1, k))));
    block.adapter.target = target;
    block.adapter.restart(log_scale);
    block.log_scale = log_scale;
    block.accepted = 0;
    block.proposed = 0;
  }
}

void Sampler::DualAveraging::restart(double log_eps) {
  mu = log_eps;
  log_eps = log_eps;
  log_eps_bar = log_eps;
  h_bar = 0.0;
  counter = 0;
}

double Sampler::DualAveraging::update(double alpha) {
  ++counter;
  const double eta = 1.0 / (static_cast<double>(counter) + t0);
  h_bar = (1.0 - eta) * h_bar + eta * (target - alpha);
  log_eps = mu - std::sqrt(static_cast<double>(counter)) / gamma * h_bar;
  const double weight = std::pow(static_cast<double>(counter), -kappa);
  log_eps_bar = weight * log_eps + (1.0 - weight) * log_eps_bar;
  return log_eps;
}

void Sampler::initialize_ensemble() {
  configure_from_ports();
  walkers_.clear();
  walker_parts_.clear();
  chain_.clear();
  log_prob_.clear();
  ln_prior_.clear();
  chi2_.clear();
  iteration_ = 0;
  acceptance_fractions_.clear();
  de_accepted_ = 0;
  de_proposed_ = 0;

  if (algorithm_ == "stretch" || algorithm_ == "slice") {
    const int n = get_number_of_walkers();
    if (n < 4)
      throw SamplerConfigurationError(
          "an ensemble of " + std::to_string(n) +
          " walkers is too small to be split into two halves that propose "
          "from each other; use at least 4");
    if (n < 2 * static_cast<int>(ndim_) && !live_dangerously_)
      throw SamplerConfigurationError(
          "an ensemble of " + std::to_string(n) + " walkers cannot span " +
          std::to_string(ndim_) + " dimensions; use at least " +
          std::to_string(2 * ndim_) + " walkers");
    walkers_ = walker_start_override_.empty() ? spread_walkers(n)
                                              : walker_start_override_;
    if (walkers_.size() < 4)
      throw SamplerConfigurationError(
          "the explicit walker start must hold at least 4 walkers");
    if (!walkers_independent(walkers_))
      throw SamplerConfigurationError(
          "Initial state has a large condition number. The walkers span "
          "less than the full parameter space, so the chain cannot explore "
          "it -- spread them out.");
    accepted_.assign(walkers_.size(), 0);
  } else if (algorithm_ == "de") {
    // chisurf's seeding: a population spread around the start, big enough
    // that the first difference vectors mean something, clipped to the
    // bounds, with member 0 left exactly at the start.
    const int n = get_number_of_chains();
    walkers_.assign(static_cast<std::size_t>(n),
                    std::vector<double>(ndim_, 0.0));
    const double floor_scale = std::max(jitter_, 1e-3);
    std::normal_distribution<double> normal(0.0, 1.0);
    for (int c = 0; c < n; ++c)
      for (unsigned int i = 0; i < ndim_; ++i) {
        double scale =
            std::fabs(initial_values_[i]) * std::max(jitter_, 1e-3) * 10.0;
        if (scale < 1e-12) scale = floor_scale;
        double v = initial_values_[i] + normal(rng_) * scale;
        walkers_[static_cast<std::size_t>(c)][i] =
            std::min(std::max(v, lower_[i]), upper_[i]);
      }
    for (unsigned int i = 0; i < ndim_; ++i)
      walkers_[0][i] = initial_values_[i];
    de_noise_.resize(ndim_);
    for (unsigned int i = 0; i < ndim_; ++i) {
      de_noise_[i] = std::fabs(initial_values_[i]) * jitter_;
      if (de_noise_[i] < 1e-15) de_noise_[i] = jitter_;
    }
    de_accepted_ = 0;
    de_proposed_ = 0;
  } else {  // metropolis
    rebuild_blocks();
    seed_blocks();
    walkers_.assign(1, initial_values_);
  }
  accepted_.assign(walkers_.size(), 0);
  walker_parts_.resize(walkers_.size());
  for (std::size_t w = 0; w < walkers_.size(); ++w)
    walker_parts_[w] = evaluate(walkers_[w]);
  initialized_ = true;
}

void Sampler::record_state() {
  for (std::size_t w = 0; w < walkers_.size(); ++w) {
    chain_.push_back(walkers_[w]);
    log_prob_.push_back(walker_parts_[w].lnpost);
    ln_prior_.push_back(walker_parts_[w].lnprior);
    chi2_.push_back(walker_parts_[w].chi2);
  }
  ++iteration_;
  acceptance_fractions_.assign(walkers_.size(), 0.0);
  for (std::size_t w = 0; w < walkers_.size(); ++w)
    acceptance_fractions_[w] =
        static_cast<double>(accepted_[w]) / static_cast<double>(iteration_);
}

// ------------------------------------------------------------------- moves

std::vector<double> Sampler::slice_along(const std::vector<double>& x,
                                         const std::vector<double>& direction,
                                         double log_p_x, int* expansions,
                                         int* contractions, bool* truncated) {
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  // The slice height, in log space: y = log p(x) + log u, u ~ U(0, 1).
  const double threshold = log_p_x + std::log(uniform(rng_));

  // The initial interval straddles x at a random offset, so the move is
  // reversible: [L, R) of unit length with x somewhere inside it.
  double left = -uniform(rng_);
  double right = left + 1.0;

  std::vector<double> probe(x.size());
  const std::size_t n = x.size();
  // A point on the line, as a lambda would be if this file were C++14.
  struct At {
    static void fill(std::vector<double>* out, const std::vector<double>& x0,
                     const std::vector<double>& d, double t, std::size_t n) {
      for (std::size_t i = 0; i < n; ++i) (*out)[i] = x0[i] + t * d[i];
    }
  };

  // Stepping out. This is the half that must be allowed to overshoot: a cap
  // that binds leaves the interval short of the slice and every draw comes
  // from a truncated line, which narrows the posterior without ever
  // producing a rejected point to notice.
  int steps = 0;
  while (steps < slice_max_steps_) {
    At::fill(&probe, x, direction, left, n);
    if (evaluate(probe).lnpost <= threshold) break;
    left -= 1.0;
    ++steps;
    ++(*expansions);
  }
  if (steps >= slice_max_steps_) *truncated = true;
  steps = 0;
  while (steps < slice_max_steps_) {
    At::fill(&probe, x, direction, right, n);
    if (evaluate(probe).lnpost <= threshold) break;
    right += 1.0;
    ++steps;
    ++(*expansions);
  }
  if (steps >= slice_max_steps_) *truncated = true;

  // Shrinkage. A draw outside the slice replaces the end it came from, so
  // the interval closes on the slice and the walk is exact.
  for (int attempt = 0; attempt < slice_max_steps_; ++attempt) {
    const double t = left + uniform(rng_) * (right - left);
    At::fill(&probe, x, direction, t, n);
    if (evaluate(probe).lnpost > threshold) return probe;
    if (t < 0.0) {
      left = t;
    } else {
      right = t;
    }
    ++(*contractions);
  }
  // The interval collapsed onto x without a point above the threshold, which
  // happens when the density is flat to numerical precision. Staying is the
  // correct answer: the walker is already in the slice.
  return x;
}

void Sampler::slice_step() {
  // Karamanis & Beutler's ensemble slice sampler, which is what zeus runs:
  // the halves and the differential direction of "stretch" above, and a
  // slice along that direction instead of a Metropolis proposal.
  const int n = static_cast<int>(walkers_.size());
  const int half = n / 2;
  std::vector<int> order(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) order[static_cast<std::size_t>(i)] = i;
  std::shuffle(order.begin(), order.end(), rng_);
  // Every move is accepted; the tally is filled so that a caller reading
  // acceptance across backends is not handed zeros for this one.
  substep_accepted_.assign(walkers_.size(), 1);

  long expansions = 0, contractions = 0;
  std::vector<double> direction(ndim_);
  for (int pass = 0; pass < 2; ++pass) {
    for (int s = 0; s < (pass == 0 ? half : n - half); ++s) {
      const int a = pass == 0 ? order[static_cast<std::size_t>(s)]
                              : order[static_cast<std::size_t>(half + s)];
      const int comp_offset = pass == 0 ? half : 0;
      const int comp_size = pass == 0 ? n - half : half;
      // The direction is the difference of **two other** walkers, not the
      // line from this one to a partner. That is not a detail: a direction
      // built from the current point makes the move's geometry depend on
      // where the walker is, the slice is then not a slice of a fixed line,
      // and the chain is not the target. Measured on the correlated Gaussian
      // the (partner - self) version samples 20% narrow with no truncation
      // and no other symptom.
      if (comp_size < 2) continue;
      std::uniform_int_distribution<int> comp_pick(0, comp_size - 1);
      const int i1 = comp_pick(rng_);
      int i2 = comp_pick(rng_);
      while (i2 == i1) i2 = comp_pick(rng_);
      const int p1 = order[static_cast<std::size_t>(comp_offset + i1)];
      const int p2 = order[static_cast<std::size_t>(comp_offset + i2)];
      for (unsigned int d = 0; d < ndim_; ++d) {
        direction[d] = slice_mu_ *
                       (walkers_[static_cast<std::size_t>(p1)][d] -
                        walkers_[static_cast<std::size_t>(p2)][d]);
      }
      int e = 0, c = 0;
      bool truncated = false;
      const std::vector<double> moved = slice_along(
          walkers_[static_cast<std::size_t>(a)], direction,
          walker_parts_[static_cast<std::size_t>(a)].lnpost, &e, &c,
          &truncated);
      walkers_[static_cast<std::size_t>(a)] = moved;
      walker_parts_[static_cast<std::size_t>(a)] = evaluate(moved);
      expansions += e;
      contractions += c;
      if (truncated) ++slice_truncations_;
    }
  }
  slice_expansions_ += expansions;
  slice_contractions_ += contractions;

  // Tuning: the scale is right when the interval neither has to be grown nor
  // shrunk much, so mu is moved toward the ratio the two counts imply. Frozen
  // once warm-up is over, because a scale that keeps adapting on the recorded
  // chain makes the chain non-Markovian.
  if (slice_tuning_ && (expansions + contractions) > 0) {
    const double ratio = 2.0 * static_cast<double>(expansions) /
                         static_cast<double>(expansions + contractions);
    slice_mu_ *= std::pow(ratio > 0.0 ? ratio : 0.5, 0.25);
    if (!(slice_mu_ > 1e-8)) slice_mu_ = 1e-8;
    if (slice_mu_ > 1e8) slice_mu_ = 1e8;
  }
}

void Sampler::stretch_step() {
  // EnsembleSampler._step: two randomly assigned halves, a stretch per
  // active walker along the line to a complementary walker, accepted in
  // log space against (ndim - 1) ln z.
  const int n = static_cast<int>(walkers_.size());
  const int half = n / 2;
  std::vector<int> order(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) order[static_cast<std::size_t>(i)] = i;
  std::shuffle(order.begin(), order.end(), rng_);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  std::vector<double> q(ndim_);
  // chisurf's EnsembleSampler.sample() carries the per-walker acceptance
  // flags of the LAST substep of a stored step into its accepted tally;
  // substep_accepted_ is those flags, folded in by run().
  substep_accepted_.assign(walkers_.size(), 0);
  for (int pass = 0; pass < 2; ++pass) {
    // pass 0: order[0..half) active, the rest the complement; pass 1: swap
    for (int s = 0; s < (pass == 0 ? half : n - half); ++s) {
      const int a = pass == 0 ? order[static_cast<std::size_t>(s)]
                              : order[static_cast<std::size_t>(half + s)];
      const int comp_offset = pass == 0 ? half : 0;
      const int comp_size = pass == 0 ? n - half : half;
      std::uniform_int_distribution<int> comp_pick(0, comp_size - 1);
      const int p =
          order[static_cast<std::size_t>(comp_offset + comp_pick(rng_))];
      // z ~ g(z) propto 1/sqrt(z) on [1/a, a], by inverting its CDF
      const double u = uniform(rng_);
      const double z = ((stretch_scale_ - 1.0) * u + 1.0) *
                       ((stretch_scale_ - 1.0) * u + 1.0) / stretch_scale_;
      const double factor = (static_cast<double>(ndim_) - 1.0) * std::log(z);
      for (unsigned int d = 0; d < ndim_; ++d) {
        const double partner = walkers_[static_cast<std::size_t>(p)][d];
        const double s_v = walkers_[static_cast<std::size_t>(a)][d];
        q[d] = partner - (partner - s_v) * z;
      }
      const Parts proposal = evaluate(q);
      const double delta =
          factor + proposal.lnpost - walker_parts_[static_cast<std::size_t>(a)].lnpost;
      if (delta > std::log(uniform(rng_))) {
        walkers_[static_cast<std::size_t>(a)] = q;
        walker_parts_[static_cast<std::size_t>(a)] = proposal;
        ++substep_accepted_[static_cast<std::size_t>(a)];
      }
    }
  }
}

void Sampler::de_generation(long generation_index) {
  // sample_differential_evolution._generation, including the every-tenth
  // gamma = 1 jump and the snooker updates with their line Jacobian.
  const double gamma0 =
      2.38 / std::sqrt(2.0 * static_cast<double>(std::max(1, static_cast<int>(ndim_))));
  const double gamma = (generation_index % 10 == 9) ? 1.0 : gamma0;
  const int n = static_cast<int>(walkers_.size());
  std::vector<int> order(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) order[static_cast<std::size_t>(i)] = i;
  std::shuffle(order.begin(), order.end(), rng_);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  std::normal_distribution<double> normal(0.0, 1.0);
  std::vector<double> proposal(ndim_, 0.0), direction(ndim_, 0.0);
  for (int oi = 0; oi < n; ++oi) {
    const int c = order[static_cast<std::size_t>(oi)];
    // Three (or two) distinct others, as rng.choice(..., replace=False)
    std::vector<int> others;
    others.reserve(static_cast<std::size_t>(n - 1));
    for (int o = 0; o < n; ++o)
      if (o != c) others.push_back(o);
    double log_jacobian = 0.0;
    std::vector<int> picked;
    if (uniform(rng_) < snooker_ && n >= 4) {
      picked = sampler_detail::sample_without_replacement(others, 3, rng_);
      const int z = picked[0], j = picked[1], k = picked[2];
      for (unsigned int d = 0; d < ndim_; ++d)
        direction[d] = walkers_[static_cast<std::size_t>(c)][d] -
                       walkers_[static_cast<std::size_t>(z)][d];
      double norm = 0.0;
      for (unsigned int d = 0; d < ndim_; ++d) norm += direction[d] * direction[d];
      if (!(norm > 0.0)) continue;
      double proj = 0.0;
      for (unsigned int d = 0; d < ndim_; ++d)
        proj += (walkers_[static_cast<std::size_t>(j)][d] -
                 walkers_[static_cast<std::size_t>(k)][d]) *
                direction[d];
      proj /= norm;
      const double w = 1.2 + uniform(rng_) * (2.2 - 1.2);
      for (unsigned int d = 0; d < ndim_; ++d)
        proposal[d] = walkers_[static_cast<std::size_t>(c)][d] +
                      w * proj * direction[d];
      double new_norm = 0.0;
      for (unsigned int d = 0; d < ndim_; ++d) {
        const double r = proposal[d] - walkers_[static_cast<std::size_t>(z)][d];
        new_norm += r * r;
      }
      if (!(new_norm > 0.0)) continue;
      log_jacobian = 0.5 * (static_cast<double>(ndim_) - 1.0) *
                     (std::log(new_norm) - std::log(norm));
    } else {
      picked = sampler_detail::sample_without_replacement(others, 2, rng_);
      const int j = picked[0], k = picked[1];
      for (unsigned int d = 0; d < ndim_; ++d)
        proposal[d] = walkers_[static_cast<std::size_t>(c)][d] +
                      gamma * (walkers_[static_cast<std::size_t>(j)][d] -
                               walkers_[static_cast<std::size_t>(k)][d]) +
                      normal(rng_) * de_noise_[d];
    }
    const Parts candidate = evaluate(proposal);
    ++de_proposed_;
    if (!std::isfinite(candidate.lnpost)) continue;
    const double delta = (candidate.lnpost -
                          walker_parts_[static_cast<std::size_t>(c)].lnpost) /
                             temp_ +
                         log_jacobian;
    if (delta > std::log(uniform(rng_))) {
      walkers_[static_cast<std::size_t>(c)] = proposal;
      walker_parts_[static_cast<std::size_t>(c)] = candidate;
      ++de_accepted_;
    }
  }
}

void Sampler::blocked_sweep(bool adapt) {
  // walk_mcmc_blocked._sweep: a correlated proposal per block, accepted in
  // log space, with dual averaging of the per-block log scale in warm-up.
  std::vector<double> current = walkers_[0];
  Parts parts = walker_parts_[0];
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  std::normal_distribution<double> normal(0.0, 1.0);
  for (std::size_t b = 0; b < blocks_.size(); ++b) {
    BlockState& block = blocks_[b];
    const std::size_t k = block.indices.size();
    std::vector<double> trial = current;
    std::vector<double> draw(k, 0.0);
    // One standard-normal vector per proposal, shared across the rows of
    // L: draw = L z. Drawing a fresh normal per matrix *entry* -- the bug
    // this replaces -- yields a diagonal proposal with the right marginal
    // widths and zero correlation, which on a strongly correlated
    // posterior mixes orders of magnitude worse and was invisible for as
    // long as the seed itself was diagonal (found 2026-09-02 by the
    // curvature seed, the first non-diagonal factor to arrive).
    std::vector<double> z(k);
    for (std::size_t j = 0; j < k; ++j) z[j] = normal(rng_);
    for (std::size_t i = 0; i < k; ++i) {
      double value = 0.0;
      for (std::size_t j = 0; j <= i; ++j)
        value += block.factor[i * k + j] * z[j];
      draw[i] = value;
    }
    const double scale = std::exp(block.log_scale);
    for (std::size_t i = 0; i < k; ++i)
      trial[static_cast<std::size_t>(block.indices[i])] =
          current[static_cast<std::size_t>(block.indices[i])] + scale * draw[i];
    const Parts trial_parts = evaluate(trial);
    ++block.proposed;
    double alpha = 0.0;
    if (std::isfinite(trial_parts.lnpost)) {
      const double delta =
          (trial_parts.lnpost - parts.lnpost) / temp_;
      alpha = delta >= 0.0 ? 1.0 : std::exp(delta);
      if (delta > std::log(uniform(rng_))) {
        current = trial;
        parts = trial_parts;
        ++block.accepted;
      }
    }
    if (adapt) block.log_scale = block.adapter.update(alpha);
  }
  walkers_[0] = current;
  walker_parts_[0] = parts;
}

// -------------------------------------------------------------------- run

void Sampler::run(int n_steps, int thin) {
  thin_ = std::max(1, thin);
  validate();
  if (!initialized_) {
    initialize_ensemble();
    // The warm-ups chisurf runs before recording: none for the ensemble
    // samplers, a spread-out population for DE, the windowed covariance
    // and scale adaptation for the blocked walk.
    const int n_samples = std::max(1, n_steps / thin_);
    int n_adapt = n_adapt_setting_;
    if (algorithm_ == "de") {
      if (n_adapt < 0)
        n_adapt = std::min(500, std::max(50, (n_samples * thin_) / 4));
      for (int g = 0; g < n_adapt; ++g) de_generation(g);
      de_accepted_ = 0;
      de_proposed_ = 0;
    } else if (algorithm_ == "slice") {
      // Warm-up here tunes one number, the direction scale, and it is tuned
      // by *running* -- the expansion/contraction counts are the signal.
      // Frozen afterwards: a scale that keeps adapting on the recorded chain
      // makes the chain non-Markovian.
      if (n_adapt < 0)
        n_adapt = std::min(200, std::max(20, (n_samples * thin_) / 20));
      const bool was_tuning = slice_tuning_;
      for (int i = 0; i < n_adapt; ++i) slice_step();
      slice_tuning_ = false;
      // The truncation counter answers a question about the recorded chain,
      // so warm-up's expansions -- taken with an untuned scale, where a long
      // stepping-out is expected -- do not count against it.
      slice_truncations_ = 0;
      slice_expansions_ = 0;
      slice_contractions_ = 0;
      (void)was_tuning;
    } else if (algorithm_ == "metropolis") {
      if (n_adapt < 0)
        n_adapt = std::min(500, std::max(100, (n_samples * thin_) / 20));
      if (n_adapt > 0) {
        int init_buffer = 0, term_buffer = 0;
        std::vector<int> window_ends;
        sampler_detail::adaptation_windows(n_adapt, init_buffer, term_buffer,
                                           window_ends);
        // chisurf's growing windows: every estimate uses all draws since
        // the end of the initial buffer, never a sliding one.
        std::vector<std::vector<double> > warmup(
            static_cast<std::size_t>(n_adapt),
            std::vector<double>(ndim_, 0.0));
        const std::size_t window_start = static_cast<std::size_t>(init_buffer);
        for (int i = 0; i < n_adapt; ++i) {
          blocked_sweep(true);
          warmup[static_cast<std::size_t>(i)] = walkers_[0];
          const bool is_end = std::find(window_ends.begin(), window_ends.end(),
                                        i + 1) != window_ends.end();
          if (is_end) {
            std::vector<const double*> visited;
            for (std::size_t d = window_start;
                 d <= static_cast<std::size_t>(i); ++d)
              visited.push_back(&warmup[d][0]);
            for (std::size_t b = 0; b < blocks_.size(); ++b) {
              BlockState& block = blocks_[b];
              // An exact curvature seed cannot be improved on by a short
              // warm-up chain; only its scale adapts (chisurf's
              // from_curvature contract).
              if (block.from_curvature) continue;
              const std::size_t k = block.indices.size();
              sampler_detail::Matrix empirical = sampler_detail::
                  regularised_covariance(
                      visited, visited.size(), block.indices);
              if (empirical.empty()) continue;
              sampler_detail::Matrix current_cov(
                  k, std::vector<double>(k, 0.0));
              double size = 0.0, empirical_size = 0.0;
              for (std::size_t r = 0; r < k; ++r)
                for (std::size_t cc = 0; cc < k; ++cc) {
                  double dot = 0.0;
                  for (std::size_t j = 0; j < k; ++j)
                    dot += block.factor[r * k + j] * block.factor[cc * k + j];
                  current_cov[r][cc] =
                      std::exp(2.0 * block.log_scale) * dot;
                }
              for (std::size_t r = 0; r < k; ++r) {
                size += current_cov[r][r];
                empirical_size += empirical[r][r];
              }
              if (!(empirical_size > 0.0 && std::isfinite(size) && size > 0.0))
                continue;
              for (std::size_t r = 0; r < k; ++r)
                for (std::size_t cc = 0; cc < k; ++cc)
                  empirical[r][cc] *= size / empirical_size;
              block.factor = sampler_detail::flatten_factor(
                  sampler_detail::cholesky_or_diagonal(empirical));
              block.log_scale = 0.0;
              block.adapter.restart(0.0);
            }
          }
        }
        for (std::size_t b = 0; b < blocks_.size(); ++b)
          blocks_[b].log_scale = blocks_[b].adapter.averaged();
        for (std::size_t b = 0; b < blocks_.size(); ++b) {
          blocks_[b].accepted = 0;
          blocks_[b].proposed = 0;
        }
      }
    }
  }

  const int n_stored = std::max(1, n_steps / thin_);
  const int total = n_stored * thin_;
  int done = 0;
  for (int s = 0; s < n_stored; ++s) {
    for (int t = 0; t < thin_; ++t) {
      if (algorithm_ == "stretch") {
        stretch_step();
        if (t == thin_ - 1)
          for (std::size_t w = 0; w < walkers_.size(); ++w)
            accepted_[w] += substep_accepted_[w];
      } else if (algorithm_ == "slice") {
        slice_step();
        // Every slice move is accepted by construction, so the tally that
        // means something for "stretch" would read 100% here and say
        // nothing. It is still filled, so a caller reading acceptance
        // across backends is not handed zeros.
        if (t == thin_ - 1)
          for (std::size_t w = 0; w < walkers_.size(); ++w)
            accepted_[w] += substep_accepted_[w];
      } else if (algorithm_ == "de") {
        de_generation(t + s * thin_ + 1);  // chisurf numbers recordings from 1
      } else {
        blocked_sweep(false);
      }
      ++done;
      if (observer_) observer_(done, total);
    }
    record_state();
  }

  if (algorithm_ == "de" && !parameters_.empty()) {
    // chisurf's DE sampler leaves the model where it started it.
    for (unsigned int i = 0; i < ndim_; ++i)
      parameters_[i]->set_value(initial_values_[i]);
    objective_node_->update();
  }
}

void Sampler::step() { run(1, 1); }

void Sampler::reset() {
  walkers_.clear();
  walker_parts_.clear();
  chain_.clear();
  log_prob_.clear();
  ln_prior_.clear();
  chi2_.clear();
  acceptance_fractions_.clear();
  iteration_ = 0;
  de_accepted_ = 0;
  de_proposed_ = 0;
  initialized_ = false;
  rng_.seed(seed_);
}

// ----------------------------------------------------------------- results

std::vector<std::vector<double> > Sampler::get_chain() const {
  return chain_;
}

std::vector<std::vector<double> > Sampler::get_chain_of_walker(
    int walker) const {
  std::vector<std::vector<double> > out;
  const std::size_t n = walkers_.size();
  if (n == 0) return out;
  if (walker < 0 || static_cast<std::size_t>(walker) >= n)
    throw SamplerConfigurationError("no walker " + std::to_string(walker));
  for (unsigned int s = 0; s < iteration_; ++s)
    out.push_back(chain_[static_cast<std::size_t>(s) * n +
                         static_cast<std::size_t>(walker)]);
  return out;
}

std::vector<std::vector<double> > Sampler::get_walkers() const {
  return walkers_;
}

std::vector<double> Sampler::get_log_prob() const { return log_prob_; }

std::vector<double> Sampler::get_lnprior() const { return ln_prior_; }

std::vector<double> Sampler::get_chi2() const { return chi2_; }

double Sampler::get_acceptance_rate() const {
  if (iteration_ == 0) return std::nan("");
  if (algorithm_ == "stretch") {
    double sum = 0.0;
    for (std::size_t w = 0; w < acceptance_fractions_.size(); ++w)
      sum += acceptance_fractions_[w];
    return acceptance_fractions_.empty()
               ? std::nan("")
               : sum / static_cast<double>(acceptance_fractions_.size());
  }
  if (algorithm_ == "de") {
    return static_cast<double>(de_accepted_) /
           static_cast<double>(std::max<long>(1, de_proposed_));
  }
  long accepted = 0, proposed = 0;
  for (std::size_t b = 0; b < blocks_.size(); ++b) {
    accepted += blocks_[b].accepted;
    proposed += blocks_[b].proposed;
  }
  return static_cast<double>(accepted) /
         static_cast<double>(std::max<long>(1, proposed));
}

std::vector<double> Sampler::get_acceptance_fractions() const {
  if (algorithm_ == "stretch") return acceptance_fractions_;
  return std::vector<double>(1, get_acceptance_rate());
}

std::vector<double> Sampler::get_block_acceptance_rates() const {
  std::vector<double> out;
  for (std::size_t b = 0; b < blocks_.size(); ++b)
    out.push_back(blocks_[b].proposed > 0
                      ? static_cast<double>(blocks_[b].accepted) /
                            static_cast<double>(blocks_[b].proposed)
                      : std::nan(""));
  return out;
}

unsigned int Sampler::get_number_of_evaluations() const {
  return n_evaluations_;
}

unsigned int Sampler::get_iteration() const { return iteration_; }

unsigned int Sampler::get_number_of_parameters() const { return ndim_; }

std::string Sampler::describe() const {
  std::ostringstream out;
  out << "Sampler(algorithm='" << algorithm_ << "', n_parameters=" << ndim_
      << ", n_walkers=" << walkers_.size() << ", iteration=" << iteration_
      << ", acceptance_rate=" << get_acceptance_rate() << ")";
  return out.str();
}

IMPBFF_END_NAMESPACE
