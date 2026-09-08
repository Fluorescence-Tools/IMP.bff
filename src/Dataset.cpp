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
  coordinates_.clear();   // the shape may have moved under them
  coordinate_names_.clear();
}

void Dataset::set_values(const std::vector<double>& values) {
  set_values(values, std::vector<int>(1, static_cast<int>(values.size())));
}

void Dataset::set_values_array(double* in_values, int n_values) {
  // The array path, because a curve's values change on every write and the
  // list path walks them one Python object at a time: measured at 2 ms for
  // 117k against 85 us. A sync that costs more than the arithmetic will be
  // skipped, and then the two halves drift.
  set_values(std::vector<double>(in_values, in_values + std::max(0, n_values)));
}

void Dataset::set_stored_variance_array(double* in_variance, int n_variance) {
  set_stored_variance(
      std::vector<double>(in_variance, in_variance + std::max(0, n_variance)));
}

void Dataset::set_mask_array(double* in_mask, int n_mask) {
  set_mask(std::vector<double>(in_mask, in_mask + std::max(0, n_mask)));
}

void Dataset::set_coordinate_array(int index, const std::string& name,
                                   double* in_coordinate, int n_coordinate) {
  set_coordinate(index, name,
                 std::vector<double>(in_coordinate,
                                     in_coordinate + std::max(0, n_coordinate)));
}

void Dataset::set_coordinate(int index, const std::string& name,
                             const std::vector<double>& values) {
  if (index < 0) IMP_THROW("a coordinate index is not negative", IMP::ValueException);
  if (values.size() != values_.size()) {
    IMP_THROW("a coordinate has one value per point: " << values_.size()
              << " expected, " << values.size() << " given. Every y has an x",
              IMP::ValueException);
  }
  if (index >= static_cast<int>(coordinates_.size())) {
    coordinates_.resize(static_cast<std::size_t>(index) + 1);
    coordinate_names_.resize(static_cast<std::size_t>(index) + 1);
  }
  coordinates_[static_cast<std::size_t>(index)] = values;
  coordinate_names_[static_cast<std::size_t>(index)] = name;
}

void Dataset::set_grid_axis(int dimension, const std::string& name,
                            const std::vector<double>& values) {
  if (dimension < 0 || dimension >= get_rank()) {
    IMP_THROW("this dataset has " << get_rank() << " dimensions and no "
              << dimension, IMP::ValueException);
  }
  if (static_cast<int>(values.size()) != shape_[dimension]) {
    IMP_THROW("dimension " << dimension << " has extent " << shape_[dimension]
              << " and " << values.size() << " coordinates were given",
              IMP::ValueException);
  }
  // Expanded rather than stored compactly: the general form is one coordinate
  // per value, and a second storage shape would be a second thing to get
  // wrong for the sake of a few kilobytes.
  std::size_t inner = 1;
  for (std::size_t d = static_cast<std::size_t>(dimension) + 1;
       d < shape_.size(); ++d) {
    inner *= static_cast<std::size_t>(shape_[d]);
  }
  std::vector<double> full(values_.size(), 0.0);
  for (std::size_t i = 0; i < values_.size(); ++i) {
    const std::size_t along =
        (i / inner) % static_cast<std::size_t>(shape_[dimension]);
    full[i] = values[along];
  }
  set_coordinate(dimension, name, full);
}

const std::vector<double>& Dataset::get_coordinate(int index) const {
  static const std::vector<double> none;
  if (index < 0 || index >= static_cast<int>(coordinates_.size())) return none;
  return coordinates_[static_cast<std::size_t>(index)];
}

std::string Dataset::get_coordinate_name(int index) const {
  if (index < 0 || index >= static_cast<int>(coordinate_names_.size())) {
    return std::string();
  }
  return coordinate_names_[static_cast<std::size_t>(index)];
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
  if (family_ == NOISE_FAMILY_PROPAGATED && propagated_variance_.empty()) {
    IMP_THROW("this dataset says its variance was propagated and none was; "
              "build it with set_linear_combination()",
              IMP::ValueException);
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

std::vector<double> Dataset::own_variance() const {
  std::vector<double> v(values_.size(), 0.0);
  for (std::size_t i = 0; i < values_.size(); ++i) {
    switch (family_) {
      case NOISE_FAMILY_POISSON: v[i] = values_[i]; break;
      case NOISE_FAMILY_STORED: v[i] = stored_variance_.empty() ? 0.0 : stored_variance_[i]; break;
      case NOISE_FAMILY_CONSTANT: v[i] = constant_variance_; break;
      case NOISE_FAMILY_PROPAGATED: v[i] = propagated_variance_.empty() ? 0.0 : propagated_variance_[i]; break;
    }
  }
  return v;
}

void Dataset::set_as_source(const std::string& name) {
  source_names_.assign(1, name);
  source_variance_.assign(1, own_variance());
  derivative_.assign(1, std::vector<double>(values_.size(), 1.0));
}

void Dataset::finish_propagation(const std::vector<double>& values,
                                 const std::string& provenance) {
  const std::size_t n = values.size();
  std::vector<double> var(n, 0.0);
  for (std::size_t k = 0; k < source_names_.size(); ++k) {
    for (std::size_t i = 0; i < n; ++i) {
      const double d = derivative_[k][i];
      var[i] += d * d * source_variance_[k][i];
    }
  }
  const std::vector<std::string> names = source_names_;
  const std::vector<std::vector<double> > svar = source_variance_, der = derivative_;
  set_values(values);                       // clears the mask, keeps nothing else
  source_names_ = names;
  source_variance_ = svar;
  derivative_ = der;
  family_ = NOISE_FAMILY_PROPAGATED;
  propagated_variance_ = var;
  stored_variance_.clear();
  provenance_ = provenance;
}

//! op: 0 add, 1 subtract, 2 multiply, 3 divide.
Dataset Dataset::binary(const Dataset& a, const Dataset& b, int op,
                        const char* symbol) {
  const std::size_t n = a.values_.size();
  if (b.values_.size() != n) {
    IMP_THROW("both operands must have the same size", IMP::ValueException);
  }
  if (a.source_names_.empty() || b.source_names_.empty()) {
    IMP_THROW("both operands must carry uncertainty; call set_as_source() on a "
              "measured channel, or build it from ones that do. An untracked "
              "operand would contribute no variance and say nothing",
              IMP::ValueException);
  }
  Dataset out;
  // The union of the sources, so that an operand appearing on both sides --
  // which is what makes anisotropy correlated -- is one source with the two
  // derivatives added, not two independent ones.
  std::vector<std::string> names = a.source_names_;
  std::vector<std::vector<double> > svar = a.source_variance_;
  std::vector<std::vector<double> > der(names.size(),
                                        std::vector<double>(n, 0.0));
  std::vector<double> values(n, 0.0), da(n, 0.0), db(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const double x = a.values_[i], y = b.values_[i];
    switch (op) {
      case 0: values[i] = x + y; da[i] = 1.0;      db[i] = 1.0;      break;
      case 1: values[i] = x - y; da[i] = 1.0;      db[i] = -1.0;     break;
      case 2: values[i] = x * y; da[i] = y;        db[i] = x;        break;
      default:
        if (y == 0.0) {
          IMP_THROW("a division by zero at point " << i, IMP::ValueException);
        }
        values[i] = x / y; da[i] = 1.0 / y; db[i] = -x / (y * y);
        break;
    }
  }
  for (std::size_t k = 0; k < names.size(); ++k) {
    for (std::size_t i = 0; i < n; ++i) der[k][i] = da[i] * a.derivative_[k][i];
  }
  for (std::size_t k = 0; k < b.source_names_.size(); ++k) {
    std::size_t at = names.size();
    for (std::size_t j = 0; j < names.size(); ++j) {
      if (names[j] == b.source_names_[k]) { at = j; break; }
    }
    if (at == names.size()) {
      names.push_back(b.source_names_[k]);
      svar.push_back(b.source_variance_[k]);
      der.push_back(std::vector<double>(n, 0.0));
    }
    for (std::size_t i = 0; i < n; ++i) {
      der[at][i] += db[i] * b.derivative_[k][i];
    }
  }
  out.source_names_ = names;
  out.source_variance_ = svar;
  out.derivative_ = der;
  out.finish_propagation(values, "(a " + std::string(symbol) + " b)");
  return out;
}

Dataset Dataset::add(const Dataset& a, const Dataset& b) { return binary(a, b, 0, "+"); }
Dataset Dataset::subtract(const Dataset& a, const Dataset& b) { return binary(a, b, 1, "-"); }
Dataset Dataset::multiply(const Dataset& a, const Dataset& b) { return binary(a, b, 2, "*"); }
Dataset Dataset::divide(const Dataset& a, const Dataset& b) { return binary(a, b, 3, "/"); }

Dataset Dataset::affine(const Dataset& a, double scale, double offset) {
  if (a.source_names_.empty()) {
    IMP_THROW("the operand must carry uncertainty", IMP::ValueException);
  }
  Dataset out;
  const std::size_t n = a.values_.size();
  std::vector<double> values(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) values[i] = scale * a.values_[i] + offset;
  out.source_names_ = a.source_names_;
  out.source_variance_ = a.source_variance_;
  out.derivative_ = a.derivative_;
  for (auto& d : out.derivative_) {
    for (double& e : d) e *= scale;         // squares in the variance
  }
  std::ostringstream p;
  p << scale << "*a + " << offset;
  out.finish_propagation(values, p.str());
  return out;
}

Dataset Dataset::transform(const Dataset& a, const std::vector<double>& values,
                           const std::vector<double>& derivative) {
  const std::size_t n = a.values_.size();
  if (values.size() != n || derivative.size() != n) {
    IMP_THROW("f(a) and f'(a) are as long as a", IMP::ValueException);
  }
  if (a.source_names_.empty()) {
    IMP_THROW("the operand must carry uncertainty", IMP::ValueException);
  }
  Dataset out;
  out.source_names_ = a.source_names_;
  out.source_variance_ = a.source_variance_;
  out.derivative_ = a.derivative_;
  for (auto& d : out.derivative_) {
    for (std::size_t i = 0; i < n; ++i) d[i] *= derivative[i];
  }
  out.finish_propagation(values, "f(a)");
  return out;
}

void Dataset::set_linear_combination(const std::vector<Dataset>& sources,
                                     const std::vector<double>& coefficients) {
  if (sources.empty()) {
    IMP_THROW("a combination needs at least one source", IMP::ValueException);
  }
  if (sources.size() != coefficients.size()) {
    IMP_THROW("a combination needs one coefficient per source",
              IMP::ValueException);
  }
  const std::size_t n = sources[0].get_values().size();
  for (const auto& s : sources) {
    if (s.get_values().size() != n) {
      IMP_THROW("every source of a combination has the same size",
                IMP::ValueException);
    }
  }
  std::vector<double> v(n, 0.0), var(n, 0.0);
  std::ostringstream prov;
  for (std::size_t k = 0; k < sources.size(); ++k) {
    const double c = coefficients[k];
    const std::vector<double>& y = sources[k].get_values();
    for (std::size_t i = 0; i < n; ++i) {
      v[i] += c * y[i];
      // c^2, not c. This is the whole point: the coefficient squares in the
      // variance and a Poisson reading of the result does not know that.
      double vi = 0.0;
      switch (sources[k].get_noise_family()) {
        case NOISE_FAMILY_POISSON:
          // The measured value, because the combination cannot be decomposed
          // into per-source predictions. See the header.
          vi = y[i];
          break;
        case NOISE_FAMILY_STORED:
        case NOISE_FAMILY_PROPAGATED:
          vi = sources[k].stored_variance_.empty()
                   ? sources[k].propagated_variance_[i]
                   : sources[k].stored_variance_[i];
          break;
        case NOISE_FAMILY_CONSTANT:
          vi = sources[k].constant_variance_;
          break;
      }
      var[i] += c * c * vi;
    }
    if (k) prov << " + ";
    prov << coefficients[k] << "*s" << k;
  }
  set_values(v);
  family_ = NOISE_FAMILY_PROPAGATED;
  propagated_variance_ = var;
  stored_variance_.clear();
  provenance_ = prov.str();
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
      case NOISE_FAMILY_PROPAGATED:
        v[i] = propagated_variance_[i];
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
        else if (family_ == NOISE_FAMILY_PROPAGATED) var = propagated_variance_[i];
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
      const double var = (family_ == NOISE_FAMILY_STORED)
                             ? stored_variance_[i]
                             : (family_ == NOISE_FAMILY_PROPAGATED
                                    ? propagated_variance_[i]
                                    : constant_variance_);
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
    case NOISE_FAMILY_PROPAGATED:
      s << "propagated (" << provenance_ << ")";
      break;
  }
  s << ", active=" << get_number_of_active_points() << "/" << get_size() << ")";
  return s.str();
}

IMPBFF_END_NAMESPACE
