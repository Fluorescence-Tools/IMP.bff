/**
 * \file Dataset.cpp
 * \brief Measured values of any rank, and the noise family that goes with them.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/Dataset.h>

#include <IMP/bff/internal/OutputView.h>

#include <cmath>
#include <limits>
#include <sstream>

IMPBFF_BEGIN_NAMESPACE

namespace {

//! Signed square root of a point's Poisson deviance.
/*! `2*(mu - y + y*log(y/mu))` is non-negative and vanishes at `y == mu`; the
    sign is the direction of the miss, so a run of like signs still means what
    it means. `y == 0` drops the logarithm, which is its limit. */
double deviance_term(double y, double mu) {
  // Spelt exactly as ChiSquared::deviance_residual spells it -- the model
  // floored at the smallest positive double so the logarithm stays finite,
  // and `y*(log y - log mu)` rather than `y*log(y/mu)`. The two differ only
  // in rounding, but after a square root and a near-cancellation that showed
  // up at 1e-8, and "the objective is unchanged" should mean it.
  const double tiny = std::numeric_limits<double>::min();
  if (!(mu > tiny)) mu = tiny;
  const double ylog = (y > 0.0) ? y * (std::log(y) - std::log(mu)) : 0.0;
  double d = 2.0 * (mu - y + ylog);
  if (!(d > 0.0)) d = 0.0;
  return d;
}

}  // namespace

Dataset::Dataset() {}

void Dataset::set_values(const std::vector<double>& values,
                         const std::vector<int>& shape) {
  if (shape.empty()) {
    IMP_THROW("a dataset needs a shape; a curve's shape is its length",
              IMP::ValueException);
  }
  std::size_t n = 1;
  for (int e : shape) {
    if (e <= 0) IMP_THROW("every extent must be positive", IMP::ValueException);
    n *= static_cast<std::size_t>(e);
  }
  if (n != values.size()) {
    std::ostringstream s;
    s << "the shape describes " << n << " values and " << values.size()
      << " were given";
    IMP_THROW(s.str(), IMP::ValueException);
  }
  values_ = values;
  shape_ = shape;
  mask_.clear();
}

void Dataset::set_values(const std::vector<double>& values) {
  set_values(values, std::vector<int>(1, static_cast<int>(values.size())));
}

void Dataset::set_mask(const std::vector<double>& mask) {
  if (!mask.empty() && mask.size() != values_.size()) {
    IMP_THROW("a mask is empty or as long as the data", IMP::ValueException);
  }
  mask_ = mask;
}

void Dataset::set_stored_variance(const std::vector<double>& variance) {
  if (!variance.empty() && variance.size() != values_.size()) {
    IMP_THROW("a stored variance is as long as the data", IMP::ValueException);
  }
  stored_variance_ = variance;
}

void Dataset::set_constant_variance(double variance) {
  if (!(variance > 0.0)) {
    IMP_THROW("a variance must be positive", IMP::ValueException);
  }
  constant_variance_ = variance;
}

void Dataset::check_model(const std::vector<double>& model) const {
  if (model.size() != values_.size()) {
    std::ostringstream s;
    s << "the model has " << model.size() << " values and the data has "
      << values_.size();
    IMP_THROW(s.str(), IMP::ValueException);
  }
  if (family_ == NOISE_FAMILY_STORED && stored_variance_.empty()) {
    // Deliberately not a fall back to ones. Unweighted least squares on counts
    // spanning four decades fits the peak and ignores the tail, silently.
    IMP_THROW("this dataset says its variance was measured and stored, and "
              "none was stored; set_stored_variance() or choose another "
              "noise family",
              IMP::ValueException);
  }
}

void Dataset::variance(const std::vector<double>& model, double** out_view,
                       int* n_out_view) const {
  check_model(model);
  std::vector<double> v(values_.size(), 0.0);
  for (std::size_t i = 0; i < values_.size(); ++i) {
    switch (family_) {
      case NOISE_FAMILY_POISSON:
        // The model's prediction, not the datum: this is the whole point.
        v[i] = model[i];
        break;
      case NOISE_FAMILY_STORED:
        v[i] = stored_variance_[i];
        break;
      case NOISE_FAMILY_CONSTANT:
        v[i] = constant_variance_;
        break;
    }
  }
  internal::copy_to_view(v, out_view, n_out_view);
}

void Dataset::residuals(const std::vector<double>& model, ResidualKind kind,
                        double** out_view, int* n_out_view) const {
  check_model(model);
  std::vector<double> r(values_.size(), 0.0);
  for (std::size_t i = 0; i < values_.size(); ++i) {
    if (!mask_.empty() && mask_[i] == 0.0) continue;
    const double y = values_[i], mu = model[i];
    switch (kind) {
      case RESIDUAL_PEARSON: {
        double var = constant_variance_;
        if (family_ == NOISE_FAMILY_POISSON) var = mu;
        else if (family_ == NOISE_FAMILY_STORED) var = stored_variance_[i];
        r[i] = (var > 0.0) ? (y - mu) / std::sqrt(var) : 0.0;
        break;
      }
      case RESIDUAL_DEVIANCE:
        r[i] = (y >= mu ? 1.0 : -1.0) * std::sqrt(deviance_term(y, mu));
        break;
      case RESIDUAL_NEYMAN:
        r[i] = (y > 0.0) ? (y - mu) / std::sqrt(y) : (y - mu);
        break;
    }
  }
  internal::copy_to_view(r, out_view, n_out_view);
}

double Dataset::objective(const std::vector<double>& model) const {
  check_model(model);
  double total = 0.0;
  for (std::size_t i = 0; i < values_.size(); ++i) {
    if (!mask_.empty() && mask_[i] == 0.0) continue;
    const double y = values_[i], mu = model[i];
    if (family_ == NOISE_FAMILY_POISSON) {
      total += deviance_term(y, mu);
    } else {
      const double var = (family_ == NOISE_FAMILY_STORED) ? stored_variance_[i]
                                                          : constant_variance_;
      if (var > 0.0) {
        const double d = y - mu;
        total += d * d / var;
      }
    }
  }
  return total;
}

int Dataset::get_number_of_active_points() const {
  if (mask_.empty()) return static_cast<int>(values_.size());
  int n = 0;
  for (double m : mask_) {
    if (m != 0.0) ++n;
  }
  return n;
}

std::string Dataset::describe() const {
  std::ostringstream s;
  s << "Dataset(shape=[";
  for (std::size_t i = 0; i < shape_.size(); ++i) {
    if (i) s << ", ";
    s << shape_[i];
  }
  s << "], noise=";
  switch (family_) {
    case NOISE_FAMILY_POISSON: s << "poisson (variance = model)"; break;
    case NOISE_FAMILY_STORED: s << "stored"; break;
    case NOISE_FAMILY_CONSTANT: s << "constant"; break;
  }
  s << ", active=" << get_number_of_active_points() << "/" << get_size() << ")";
  return s.str();
}

IMPBFF_END_NAMESPACE
