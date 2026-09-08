/**
 * \file ChiSquared.cpp
 * \brief The data misfit of a model curve, as a node in the model graph.
 *
 * Copyright 2007-2023 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/ChiSquared.h>

#include <algorithm>
#include <new>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

IMPBFF_BEGIN_NAMESPACE

namespace {

//! Signed Poisson deviance residual for one point.
/** ChiSurf's ``deviance_residuals``: the model is floored at the smallest
    positive double so its logarithm stays finite, the ``y log(y/mu)`` term
    is zero where there are no counts, and the deviance is clipped at zero
    before the square root so rounding cannot produce a NaN. */
double deviance_residual(double y, double mu) {
  const double tiny = std::numeric_limits<double>::min();
  if (!(mu > tiny)) mu = tiny;
  const double ylog = (y > 0.0) ? y * (std::log(y) - std::log(mu)) : 0.0;
  double dev = 2.0 * (mu - y + ylog);
  if (!(dev > 0.0)) dev = 0.0;
  const double sign = (mu > y) ? 1.0 : ((mu < y) ? -1.0 : 0.0);
  return sign * std::sqrt(dev);
}

}  // namespace

ChiSquared::ChiSquared(const std::string& name) : Node(name) {}

void ChiSquared::set_data(const std::vector<double>& y,
                          const std::vector<double>& ey) {
  if (!ey.empty() && ey.size() != y.size()) {
    throw std::domain_error(
        "ChiSquared::set_data: the errors must be as long as the data");
  }
  data_y_ = y;
  data_ey_ = ey;
  set_valid(false);
}

void ChiSquared::set_fit_range(int xmin, int xmax) {
  if (xmax >= 0 && xmax < xmin) {
    throw std::domain_error(
        "ChiSquared::set_fit_range: the range end precedes its start");
  }
  xmin_ = xmin;
  xmax_ = xmax;
  set_valid(false);
}

void ChiSquared::set_mask(const std::vector<double>& mask) {
  mask_ = mask;
  set_valid(false);
}

void ChiSquared::set_noise_model(NoiseModel model) {
  noise_model_ = model;
  set_valid(false);
}

void ChiSquared::set_noise_model_name(const std::string& name) {
  if (name == "poisson" || name == "mle" || name == "ml") {
    set_noise_model(NOISE_POISSON);
  } else if (name == "default" || name == "neyman" || name == "wls" ||
             name == "chi2") {
    set_noise_model(NOISE_NEYMAN);
  } else {
    throw std::domain_error(
        "ChiSquared: unknown noise model '" + name +
        "'; expected 'default' or 'poisson'");
  }
}

std::string ChiSquared::get_noise_model_name() const {
  return noise_model_ == NOISE_POISSON ? "poisson" : "default";
}

void ChiSquared::resolve_window(int* begin, int* end) const {
  const int n = static_cast<int>(data_y_.size());
  int lo = std::max(0, xmin_);
  int hi = (xmax_ < 0) ? n : std::min(n, xmax_);
  if (hi < lo) hi = lo;
  *begin = lo;
  *end = hi;
}

void ChiSquared::set_data_arrays(double* in_data_y, int n_data_y,
                                 double* in_data_ey, int n_data_ey) {
  data_y_.assign(in_data_y, in_data_y + std::max(0, n_data_y));
  data_ey_.assign(in_data_ey, in_data_ey + std::max(0, n_data_ey));
  // `set_valid(false)`, as the std::vector overload does. Without it the
  // numpy door -- the one a per-iteration caller is told to use -- left the
  // node claiming to be valid, so the next `update()` skipped `evaluate()`
  // and served residuals against the *previous* data. The two doors must
  // agree about more than the numbers.
  set_valid(false);
}

void ChiSquared::set_mask_array(double* in_mask_a, int n_mask_a) {
  mask_.assign(in_mask_a, in_mask_a + std::max(0, n_mask_a));
  set_valid(false);
}

void ChiSquared::compute_weighted_residuals_array(double* in_model_y,
                                                  int n_model_y,
                                                  double** out_wres,
                                                  int* n_out_wres) const {
  // Wrapping the caller's buffer rather than copying it: the vector overload
  // exists for callers who already hold one, and going through it here would
  // add exactly the copy this entry point is for.
  int begin = 0, end = 0;
  resolve_window(&begin, &end);
  const int window = end - begin;
  const int available = std::max(0, n_model_y) - begin;
  const int n = std::max(0, std::min(window, available));

  double* out = static_cast<double*>(
      std::malloc(std::max<std::size_t>(static_cast<std::size_t>(n), 1) *
                  sizeof(double)));
  if (out == nullptr) throw std::bad_alloc();

  for (int i = 0; i < n; ++i) {
    const std::size_t k = static_cast<std::size_t>(begin + i);
    const double y = data_y_[k];
    const double mu = in_model_y[k];
    if (noise_model_ == NOISE_POISSON) {
      out[i] = deviance_residual(y, mu);
    } else {
      const double e = data_ey_.empty() ? 1.0 : data_ey_[k];
      out[i] = (y - mu) / e;
    }
  }

  if (!mask_.empty() && n > 0) {
    const int m = static_cast<int>(mask_.size());
    const int mlo = std::max(0, xmin_);
    if (mlo + n <= m) {
      for (int i = 0; i < n; ++i) {
        out[i] *= mask_[static_cast<std::size_t>(mlo + i)];
      }
    }
  }

  *out_wres = out;
  *n_out_wres = n;
}

void ChiSquared::set_dataset(const Dataset& dataset) {
  dataset_ = dataset;
  has_dataset_ = true;
  set_valid(false);
}

std::vector<double> ChiSquared::compute_weighted_residuals(
    const std::vector<double>& model_y) const {
  if (has_dataset_) {
    if (xmin_ != 0 || xmax_ >= 0) {
      IMP_THROW("this ChiSquared scores a Dataset, which carries its own mask, "
                "and an index window was also set. Two masking mechanisms that "
                "disagree exclude the wrong points quietly, so pick one: mask "
                "the dataset, or drop the window",
                IMP::ValueException);
    }
    ResidualKind kind = residual_kind_;
    if (!has_residual_kind_) {
      // The family's own: the deviance is the likelihood's residual for
      // counts, Pearson is the natural one for a measured variance.
      kind = (dataset_.get_noise_family() == NOISE_FAMILY_POISSON)
                 ? RESIDUAL_DEVIANCE
                 : RESIDUAL_PEARSON;
    }
    double* view = nullptr;
    int n_view = 0;
    dataset_.residuals(model_y, kind, &view, &n_view);
    std::vector<double> out(view, view + n_view);
    std::free(view);
    return out;
  }
  int begin = 0, end = 0;
  resolve_window(&begin, &end);

  // ChiSurf slices data and model by the same window and then truncates to
  // whichever came out shorter, so a model that stops early is not an error.
  const int window = end - begin;
  const int model_available = static_cast<int>(model_y.size()) - begin;
  const int n = std::max(0, std::min(window, model_available));

  std::vector<double> wres(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    const double y = data_y_[static_cast<std::size_t>(begin + i)];
    const double mu = model_y[static_cast<std::size_t>(begin + i)];
    if (noise_model_ == NOISE_POISSON) {
      wres[static_cast<std::size_t>(i)] = deviance_residual(y, mu);
    } else {
      const double e = data_ey_.empty()
                           ? 1.0
                           : data_ey_[static_cast<std::size_t>(begin + i)];
      wres[static_cast<std::size_t>(i)] = (y - mu) / e;
    }
  }

  // The mask is indexed like the data and applied over the same window.
  // ChiSurf skips it outright when the window no longer matches the
  // residuals, rather than silently masking the wrong points.
  if (!mask_.empty() && n > 0) {
    const int m = static_cast<int>(mask_.size());
    const int mlo = std::max(0, xmin_);
    const int mhi = std::min(m, std::max(mlo, xmax_ < 0 ? m : xmax_));
    if (mhi - mlo == n) {
      for (int i = 0; i < n; ++i) {
        wres[static_cast<std::size_t>(i)] *=
            mask_[static_cast<std::size_t>(mlo + i)];
      }
    }
  }
  return wres;
}

double ChiSquared::compute_chi2(const std::vector<double>& model_y) const {
  const std::vector<double> wres = compute_weighted_residuals(model_y);
  double chi2 = 0.0;
  for (double r : wres) chi2 += r * r;
  // ChiSurf turns a NaN misfit into an infinite one so a sampler rejects it
  // instead of propagating the NaN into the posterior.
  if (std::isnan(chi2)) return std::numeric_limits<double>::infinity();
  return chi2;
}

double ChiSquared::get_chi2r(int n_free) const {
  const double dof =
      static_cast<double>(wres_.size()) - static_cast<double>(n_free) - 1.0;
  return chi2_ / dof;
}

void ChiSquared::update() {
  const std::shared_ptr<Port> model_port = get_input_port(model_key_);
  if (model_port) model_port->set_sanitize(false);
  const std::shared_ptr<Port> res = get_output_port(residuals_key_);
  if (res) res->set_sanitize(false);
  Node::update();
}

void ChiSquared::evaluate() {
  const std::shared_ptr<Port> model_port = get_input_port(model_key_);
  if (!model_port) {
    throw std::domain_error(
        "ChiSquared '" + get_name() + "': no input port '" + model_key_ +
        "' carrying the model curve");
  }
  wres_ = compute_weighted_residuals(model_port->get_values_ref());
  chi2_ = 0.0;
  for (double r : wres_) chi2_ += r * r;
  if (std::isnan(chi2_)) chi2_ = std::numeric_limits<double>::infinity();

  const std::shared_ptr<Port> out = get_output_port(get_name());
  if (!out) {
    throw std::domain_error(
        "ChiSquared '" + get_name() +
        "' writes chi-square to the output port keyed by its own name, "
        "which this node does not have");
  }
  out->set_value(chi2_);
  // The residuals themselves, when the graph asked for them. Absent by
  // default: a `Sampler` wants the scalar and would otherwise pay for a
  // copy of the whole residual vector on every move.
  const std::shared_ptr<Port> res = get_output_port(residuals_key_);
  if (res) res->set_value_vector(wres_);
  set_valid(true);
}

std::string ChiSquared::describe() const {
  std::ostringstream out;
  out << "points         : " << data_y_.size() << "\n"
      << "noise model    : " << get_noise_model_name() << "\n"
      << "fit range      : [" << xmin_ << ", "
      << (xmax_ < 0 ? static_cast<int>(data_y_.size()) : xmax_) << ")\n"
      << "residuals      : " << wres_.size() << "\n"
      << "chi2           : " << chi2_ << "\n";
  return out.str();
}

void weighted_residuals(
    double* in_data_y, int n_data_y,
    double* in_data_ey, int n_data_ey,
    double* in_model_y, int n_model_y,
    int xmin, int xmax, const std::string& noise_model,
    double** out_wres, int* n_out_wres) {
  const int n_data = std::max(0, n_data_y);
  const int lo = std::max(0, xmin);
  int hi = (xmax < 0) ? n_data : std::min(n_data, xmax);
  if (hi < lo) hi = lo;
  const int available = std::max(0, n_model_y) - lo;
  const int n = std::max(0, std::min(hi - lo, available));

  double* out = static_cast<double*>(
      std::malloc(std::max<std::size_t>(static_cast<std::size_t>(n), 1) *
                  sizeof(double)));
  if (out == nullptr) throw std::bad_alloc();

  const bool poisson = (noise_model == "poisson");
  const bool have_ey = n_data_ey > 0;
  for (int i = 0; i < n; ++i) {
    const int k = lo + i;
    const double y = in_data_y[k];
    const double mu = in_model_y[k];
    if (poisson) {
      out[i] = deviance_residual(y, mu);
    } else {
      out[i] = (y - mu) / (have_ey ? in_data_ey[k] : 1.0);
    }
  }
  *out_wres = out;
  *n_out_wres = n;
}

IMPBFF_END_NAMESPACE
