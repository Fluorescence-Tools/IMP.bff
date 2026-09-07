/**
 *  \file IMP/bff/CrosstalkMatrix.cpp
 *  \brief The excitation and emission crosstalk matrices: implementation.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/CrosstalkMatrix.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>

#include <Eigen/Dense>

IMPBFF_BEGIN_NAMESPACE

// --------------------------------------------------------- CrosstalkMatrix

CrosstalkMatrix::CrosstalkMatrix() {}

CrosstalkMatrix::CrosstalkMatrix(const std::vector<std::string>& rows,
                                 const std::vector<std::string>& columns,
                                 const std::vector<double>& values)
    : rows_(rows), columns_(columns), values_(values) {
  const std::size_t expected =
      static_cast<std::size_t>(rows.size()) * columns.size();
  if (values.size() != expected) {
    throw std::invalid_argument(
        "CrosstalkMatrix: got " + std::to_string(values.size()) +
        " values for a " + std::to_string(rows.size()) + " x " +
        std::to_string(columns.size()) + " matrix (expected " +
        std::to_string(expected) + ")");
  }
}

int CrosstalkMatrix::index_of(const std::vector<std::string>& labels,
                              const std::string& label) const {
  for (std::size_t i = 0; i < labels.size(); ++i) {
    if (labels[i] == label) return static_cast<int>(i);
  }
  return -1;
}

double CrosstalkMatrix::get(int row, int column) const {
  if (row < 0 || row >= get_n_rows() || column < 0 ||
      column >= get_n_columns()) {
    throw std::out_of_range("CrosstalkMatrix::get: index out of range");
  }
  return values_[static_cast<std::size_t>(row) * columns_.size() + column];
}

double CrosstalkMatrix::get(const std::string& row,
                            const std::string& column) const {
  const int i = index_of(rows_, row);
  const int j = index_of(columns_, column);
  if (i < 0 || j < 0) return 0.0;
  return get(i, j);
}

CrosstalkMatrix CrosstalkMatrix::select(
        const std::vector<std::string>& rows,
        const std::vector<std::string>& columns) const {
  std::vector<double> out(rows.size() * columns.size(), 0.0);
  for (std::size_t i = 0; i < rows.size(); ++i) {
    const int si = index_of(rows_, rows[i]);
    if (si < 0) continue;
    for (std::size_t j = 0; j < columns.size(); ++j) {
      const int sj = index_of(columns_, columns[j]);
      if (sj < 0) continue;
      out[i * columns.size() + j] = get(si, sj);
    }
  }
  return CrosstalkMatrix(rows, columns, out);
}

void CrosstalkMatrix::scale_row(const std::string& label, double factor) {
  const int i = index_of(rows_, label);
  if (i < 0) {
    throw std::out_of_range("CrosstalkMatrix::scale_row: no row " + label);
  }
  for (int j = 0; j < get_n_columns(); ++j) {
    values_[static_cast<std::size_t>(i) * columns_.size() + j] *= factor;
  }
}

void CrosstalkMatrix::scale_column(const std::string& label, double factor) {
  const int j = index_of(columns_, label);
  if (j < 0) {
    throw std::out_of_range("CrosstalkMatrix::scale_column: no column " +
                            label);
  }
  for (int i = 0; i < get_n_rows(); ++i) {
    values_[static_cast<std::size_t>(i) * columns_.size() + j] *= factor;
  }
}

// ------------------------------------------------------------- the algebra

namespace {

//! A row-major ``(rows, columns)`` view of a flat buffer.
Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                               Eigen::RowMajor>>
flat_matrix(const std::vector<double>& values, int rows, int columns,
            const char* what) {
  const std::size_t expected =
      static_cast<std::size_t>(rows) * columns;
  if (values.size() != expected) {
    throw std::invalid_argument(std::string(what) + ": expected " +
                                std::to_string(expected) + " values, got " +
                                std::to_string(values.size()));
  }
  return Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                        Eigen::RowMajor>>(
      values.data(), rows, columns);
}

//! One column of non-negative least squares (Lawson-Hanson active set).
/**
    The classical algorithm: grow the passive set by the largest positive
    gradient, solve the passive subproblem unconstrained, and step back along
    the line to the first coefficient that would turn negative until the
    subproblem is feasible. `maxiter = 3 n` matches the reference
    implementation's default.
 */
Eigen::VectorXd nnls_column(const Eigen::MatrixXd& a,
                            const Eigen::VectorXd& y) {
  const int n = static_cast<int>(a.cols());
  const double eps = std::numeric_limits<double>::epsilon();
  const double a_max = a.size() > 0 ? a.cwiseAbs().maxCoeff() : 0.0;
  const double tol = 10.0 * eps * a_max * static_cast<double>(
      std::max(a.rows(), Eigen::Index(1)));
  const int max_iter = 3 * n;

  Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
  std::vector<int> passive;
  passive.reserve(n);

  Eigen::VectorXd w = a.transpose() * y;
  for (int iter = 0; iter <= max_iter; ++iter) {
    // the entering coefficient: the largest gradient off the passive set
    int enter = -1;
    double best = tol;
    for (int i = 0; i < n; ++i) {
      if (std::find(passive.begin(), passive.end(), i) != passive.end()) {
        continue;
      }
      if (w[i] > best) {
        best = w[i];
        enter = i;
      }
    }
    if (enter < 0) break;  // KKT satisfied: no improving inactive direction
    passive.push_back(enter);
    std::sort(passive.begin(), passive.end());

    // make the passive subproblem feasible
    for (;;) {
      if (passive.empty()) {
        x.setZero();
        break;
      }
      const int p = static_cast<int>(passive.size());
      Eigen::MatrixXd ap(a.rows(), p);
      for (int c = 0; c < p; ++c) ap.col(c) = a.col(passive[c]);
      Eigen::VectorXd z = ap.colPivHouseholderQr().solve(y);

      int blocking = -1;
      double alpha = std::numeric_limits<double>::infinity();
      for (int c = 0; c < p; ++c) {
        if (z[c] <= 0.0) {
          const double denom = x[passive[c]] - z[c];
          const double t =
              denom > 0.0 ? x[passive[c]] / denom
                          : 0.0;  // already at/below zero: leave immediately
          if (t < alpha) {
            alpha = t;
            blocking = c;
          }
        }
      }
      if (alpha >= 1.0) {
        // feasible: accept, move everything passive to x
        for (int c = 0; c < p; ++c) x[passive[c]] = z[c];
        break;
      }
      for (int c = 0; c < p; ++c) {
        x[passive[c]] += alpha * (z[c] - x[passive[c]]);
      }
      // drop every coefficient that hit zero, then re-solve
      std::vector<int> still;
      for (int c = 0; c < p; ++c) {
        if (x[passive[c]] <= tol) {
          x[passive[c]] = 0.0;
        } else {
          still.push_back(passive[c]);
        }
      }
      passive = still;
      if (blocking >= 0 && passive.empty()) break;
    }

    w = a.transpose() * (y - a * x);
  }
  // numerical slack: nothing may leave negative
  return x.cwiseMax(0.0);
}

//! The solve behind `crosstalk_invert_mixing` for a whole block of items.
/**
    `a` is the transposed mixing matrix (detectors x sources); `y` holds the
    measured signals row-major as (detectors x items). One solve per item --
    non-negative least squares does not factor, so per-item solves are also
    what the regularised and unconstrained paths share.
 */
Eigen::MatrixXd invert_block(const Eigen::MatrixXd& a,
                             const Eigen::MatrixXd& y, bool nonneg,
                             double ridge) {
  const int n_src = static_cast<int>(a.cols());
  const int n_items = static_cast<int>(y.cols());
  Eigen::MatrixXd out(n_src, n_items);

  if (nonneg) {
    // the ridge as an augmented system: rows sqrt(lambda) I with zero
    // targets, the same construction the Python side used
    Eigen::MatrixXd a_solve = a;
    if (ridge > 0.0) {
      a_solve.conservativeResize(a.rows() + n_src, n_src);
      a_solve.bottomRows(n_src).setZero();
      a_solve.bottomRows(n_src).setIdentity();
      a_solve.bottomRows(n_src) *= std::sqrt(ridge);
    }
    // zero once: the augmentation rows stay zero, the measurement head is
    // overwritten per item -- zeroing inside the loop would wipe it when
    // there is no augmentation at all and the tail *is* the head
    Eigen::VectorXd y_aug = Eigen::VectorXd::Zero(a_solve.rows());
    for (int item = 0; item < n_items; ++item) {
      y_aug.head(y.rows()) = y.col(item);
      out.col(item) = nnls_column(a_solve, y_aug);
    }
    return out;
  }
  if (ridge > 0.0) {
    // Tikhonov closed form: x = (A^T A + lambda I)^-1 A^T y
    Eigen::MatrixXd gram = a.transpose() * a;
    gram.diagonal().array() += ridge;
    out = gram.partialPivLu().solve(a.transpose() * y);
    return out;
  }
  // minimum-norm least squares: the pseudo-inverse path a numpy caller
  // means, via a complete orthogonal decomposition
  Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(a);
  out = cod.pseudoInverse() * y;
  return out;
}

}  // namespace

void crosstalk_apply_mixing(const std::vector<double>& matrix, int n_sources,
                            int n_detectors,
                            const std::vector<double>& sources,
                            double** out_view, int* n_out_view) {
  const Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                       Eigen::RowMajor>>
      m = flat_matrix(matrix, n_sources, n_detectors, "crosstalk_apply_mixing");
  if (n_sources == 0 ||
      sources.size() % static_cast<std::size_t>(n_sources) != 0) {
    throw std::invalid_argument(
        "crosstalk_apply_mixing: sources size is not a multiple of "
        "n_sources");
  }
  const std::size_t n_items = sources.size() / n_sources;
  const Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                       Eigen::RowMajor>>
      s(sources.data(), n_sources, static_cast<int>(n_items));
  Eigen::MatrixXd out = m.transpose() * s;

  double* buf = static_cast<double*>(
      std::malloc(sizeof(double) * out.size()));
  if (buf == nullptr) throw std::bad_alloc();
  Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                           Eigen::RowMajor>>(buf, out.rows(), out.cols()) =
      out;
  *out_view = buf;
  *n_out_view = static_cast<int>(out.size());
}

void crosstalk_invert_mixing(const std::vector<double>& matrix, int n_sources,
                             int n_detectors,
                             const std::vector<double>& measured, bool nonneg,
                             double ridge, double** out_view,
                             int* n_out_view) {
  const Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                       Eigen::RowMajor>>
      m = flat_matrix(matrix, n_sources, n_detectors,
                      "crosstalk_invert_mixing");
  if (n_detectors == 0 ||
      measured.size() % static_cast<std::size_t>(n_detectors) != 0) {
    throw std::invalid_argument(
        "crosstalk_invert_mixing: measured size is not a multiple of "
        "n_detectors");
  }
  const std::size_t n_items = measured.size() / n_detectors;
  const Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                       Eigen::RowMajor>>
      y(measured.data(), n_detectors, static_cast<int>(n_items));

  Eigen::MatrixXd out = invert_block(m.transpose(), y, nonneg, ridge);

  double* buf = static_cast<double*>(
      std::malloc(sizeof(double) * out.size()));
  if (buf == nullptr) throw std::bad_alloc();
  Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                           Eigen::RowMajor>>(buf, out.rows(), out.cols()) =
      out;
  *out_view = buf;
  *n_out_view = static_cast<int>(out.size());
}

void crosstalk_shuffle_unmix(const std::vector<double>& matrix, int n_sources,
                             int n_detectors,
                             const std::vector<double>& counts,
                             const std::vector<double>& abundances,
                             unsigned long long seed, double** out_view,
                             int* n_out_view) {
  const Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                       Eigen::RowMajor>>
      m = flat_matrix(matrix, n_sources, n_detectors,
                      "crosstalk_shuffle_unmix");
  if (n_detectors == 0 ||
      counts.size() % static_cast<std::size_t>(n_detectors) != 0) {
    throw std::invalid_argument(
        "crosstalk_shuffle_unmix: counts size is not a multiple of "
        "n_detectors");
  }
  const std::size_t n_items = counts.size() / n_detectors;
  if (!abundances.empty() &&
      abundances.size() != n_sources * n_items) {
    throw std::invalid_argument(
        "crosstalk_shuffle_unmix: abundances must be empty or "
        "n_sources * n_items");
  }

  // spectral shape B[k, d] = P(detector d | photon from source k):
  // each row of M normalised to a distribution over detectors
  Eigen::VectorXd row_mass = m.rowwise().sum();
  Eigen::MatrixXd b(n_sources, n_detectors);
  for (int k = 0; k < n_sources; ++k) {
    for (int d = 0; d < n_detectors; ++d) {
      b(k, d) = row_mass[k] > 0.0 ? m(k, d) / row_mass[k] : 0.0;
    }
  }

  const Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                       Eigen::RowMajor>>
      y(counts.data(), n_detectors, static_cast<int>(n_items));

  // the abundance prior, estimated where it was not supplied
  Eigen::MatrixXd a(n_sources, static_cast<int>(n_items));
  if (abundances.empty()) {
    a = invert_block(m.transpose(), y, /*nonneg=*/true, /*ridge=*/0.0);
  } else {
    a = Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                       Eigen::RowMajor>>(
        abundances.data(), n_sources, static_cast<int>(n_items));
  }

  std::mt19937_64 generator(seed);
  Eigen::MatrixXd out = Eigen::MatrixXd::Zero(n_sources,
                                              static_cast<int>(n_items));
  for (int det = 0; det < n_detectors; ++det) {
    // posterior over sources for photons in this detector: p[k] ∝ a_k B[k,d]
    // (an array broadcast, not `a * b.col(det).asDiagonal()`: Eigen 5.0.1
    // mis-evaluates that product when the diagonal wraps a Block -- it
    // silently degrades to a * b(0, det))
    Eigen::MatrixXd w =
        a.array().colwise() * b.col(det).array();
    Eigen::VectorXd spectral = b.col(det);
    const double sp_tot = spectral.sum();
    for (int item = 0; item < static_cast<int>(n_items); ++item) {
      const double tot = w.col(item).sum();
      for (int k = 0; k < n_sources; ++k) {
        w(k, item) = tot > 0.0
                         ? w(k, item) / tot
                         : (sp_tot > 0.0 ? spectral[k] / sp_tot : 0.0);
      }
    }
    // split the integer detector counts across sources by a binomial chain;
    // the last source takes the rest, so the total is preserved exactly
    for (int item = 0; item < static_cast<int>(n_items); ++item) {
      long long remaining =
          static_cast<long long>(std::llround(y(det, item)));
      double cumulative = 0.0;
      for (int k = 0; k < n_sources; ++k) {
        if (k == n_sources - 1) {
          out(k, item) += static_cast<double>(remaining);
          break;
        }
        const double denom = 1.0 - cumulative;
        double cond = denom > 1e-12 ? w(k, item) / denom : 0.0;
        cond = std::min(1.0, std::max(0.0, cond));
        std::binomial_distribution<long long> draw(remaining, cond);
        const long long assigned = draw(generator);
        out(k, item) += static_cast<double>(assigned);
        remaining -= assigned;
        cumulative += w(k, item);
      }
    }
  }

  double* buf = static_cast<double*>(
      std::malloc(sizeof(double) * out.size()));
  if (buf == nullptr) throw std::bad_alloc();
  Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                           Eigen::RowMajor>>(buf, out.rows(), out.cols()) =
      out;
  *out_view = buf;
  *n_out_view = static_cast<int>(out.size());
}

IMPBFF_END_NAMESPACE
