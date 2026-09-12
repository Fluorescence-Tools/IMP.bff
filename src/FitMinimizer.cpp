/**
 *  \file FitMinimizer.cpp
 *  \brief ChiSurf's bounded least-squares optimiser, in C++.
 *
 *  The Levenberg-Marquardt core is MINPACK's `lmdif` and its helpers
 *  (More, Garbow & Hillstrom, ANL-80-74; public domain), transcribed rather
 *  than re-derived so that a fit which converged under
 *  `scipy.optimize.leastsq` converges here to the same answer from the same
 *  start. The bounds transform, the tolerances and the progress arithmetic
 *  around it are chisurf's `leastsqbound`.
 *
 *  Matrices follow MINPACK's column-major convention -- element (i, j) of an
 *  `m x n` matrix with leading dimension `lda` is `a[j * lda + i]` -- because
 *  the transcription is line by line and re-indexing it would be the one
 *  place a transcription can silently go wrong.
 *
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/FitMinimizer.h>

#include <IMP/bff/GraphNode.h>
#include <IMP/bff/GraphPort.h>
#include <IMP/bff/internal/OutputView.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <utility>

IMPBFF_BEGIN_NAMESPACE

const int MINIMIZER_EXPECTED_ITERATIONS = 6;

namespace {

//! chisurf's `_BAR_HALF_LIFE`: how many multiples of the estimate the bar
//! takes to cover its first half.
const double BAR_HALF_LIFE = 2.0;

//! chisurf's `_MAX_RUNNING_RATIO`: a full bar is reserved for completion.
const double MAX_RUNNING_RATIO = 0.99;

const double MACHINE_EPS = std::numeric_limits<double>::epsilon();

//! MINPACK's `dpmpar(2)`, the smallest positive magnitude.
const double DWARF = std::numeric_limits<double>::min();

//! chisurf's `_is_unbounded`: `None`, infinite and NaN all mean "no bound".
/*! Only `None` used to be recognised, so a parameter declared `(-inf, inf)`
    fell through to the two-sided branch and was mapped with
    `arcsin(inf/inf - 1)` -> NaN, silently poisoning the internal vector. */
bool is_unbounded(double v) { return !std::isfinite(v); }

//! Eigen-decomposition of a symmetric matrix by cyclic Jacobi rotations.
/*!
    `a` is `n x n` row-major and is consumed; the eigenvalues come back in
    `values` and the eigenvectors as the *columns* of `vectors`, which is
    LAPACK's convention and therefore `scipy.linalg.eigh`'s.

    Jacobi rather than a tridiagonal reduction because the matrices here are
    the size of a fit's free vector -- three to a dozen -- and Jacobi is
    thirty lines with no workspace, no pivoting and no LAPACK dependency. It
    is also the more accurate of the two for small matrices: the rotations
    are orthogonal to working precision, so a nearly singular `J'J` keeps the
    relative accuracy of its small eigenvalues, which is exactly the corner
    this is used in. `Sampler.cpp` carries the eigenvalue-only sibling of
    this routine; the covariance needs the vectors too, so it cannot share.
*/
void symmetric_eigen(std::vector<double>& a, int n,
                     std::vector<double>* values,
                     std::vector<double>* vectors) {
  values->assign(n, 0.0);
  vectors->assign(static_cast<std::size_t>(n) * n, 0.0);
  for (int i = 0; i < n; ++i) (*vectors)[i * n + i] = 1.0;
  if (n <= 0) return;
  for (int sweep = 0; sweep < 100; ++sweep) {
    double off = 0.0;
    for (int p = 0; p < n; ++p)
      for (int q = p + 1; q < n; ++q) off += a[p * n + q] * a[p * n + q];
    if (!(off > 0.0)) break;
    for (int p = 0; p < n; ++p) {
      for (int q = p + 1; q < n; ++q) {
        const double apq = a[p * n + q];
        if (apq == 0.0) continue;
        const double theta = (a[q * n + q] - a[p * n + p]) / (2.0 * apq);
        const double t = (theta >= 0.0 ? 1.0 : -1.0) /
                         (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
        const double c = 1.0 / std::sqrt(t * t + 1.0);
        const double s = t * c;
        // A' = G' A G, applied as the two one-sided sweeps rather than a
        // matrix product: only rows p, q and columns p, q change.
        for (int k = 0; k < n; ++k) {
          const double akp = a[k * n + p], akq = a[k * n + q];
          a[k * n + p] = c * akp - s * akq;
          a[k * n + q] = s * akp + c * akq;
        }
        for (int k = 0; k < n; ++k) {
          const double apk = a[p * n + k], aqk = a[q * n + k];
          a[p * n + k] = c * apk - s * aqk;
          a[q * n + k] = s * apk + c * aqk;
        }
        for (int k = 0; k < n; ++k) {
          const double vkp = (*vectors)[k * n + p];
          const double vkq = (*vectors)[k * n + q];
          (*vectors)[k * n + p] = c * vkp - s * vkq;
          (*vectors)[k * n + q] = s * vkp + c * vkq;
        }
      }
    }
  }
  for (int i = 0; i < n; ++i) (*values)[i] = a[i * n + i];
}

//! `scipy.linalg.pinvh` for a symmetric matrix, `n x n` row-major.
/*!
    The Moore-Penrose pseudo-inverse through the eigen-decomposition, with
    scipy's own default cutoff: an eigenvalue is inverted when its magnitude
    exceeds `max|w| * n * machine eps`, and dropped otherwise. Written to
    match rather than to improve on it, because chisurf's error bars are
    pinned to `pinvh`'s answer and the two must not disagree on a
    rank-deficient fit -- which is the case the cutoff exists for.

    A true inverse would do for a full-rank `J'J`, and this reduces to one
    there. It is the degenerate fit -- two parameters the data cannot tell
    apart -- where the two differ, and that fit is common enough that
    chisurf's numpy path has always used `pinvh`.
*/
std::vector<double> pseudo_inverse_symmetric(std::vector<double> a, int n) {
  std::vector<double> w, v;
  symmetric_eigen(a, n, &w, &v);
  double wmax = 0.0;
  for (int i = 0; i < n; ++i) wmax = std::max(wmax, std::fabs(w[i]));
  const double cutoff = wmax * static_cast<double>(n) * MACHINE_EPS;
  std::vector<double> out(static_cast<std::size_t>(n) * n, 0.0);
  for (int k = 0; k < n; ++k) {
    if (!(std::fabs(w[k]) > cutoff)) continue;
    const double inv = 1.0 / w[k];
    for (int i = 0; i < n; ++i) {
      const double vik = v[i * n + k];
      if (vik == 0.0) continue;
      for (int j = 0; j < n; ++j) out[i * n + j] += vik * inv * v[j * n + k];
    }
  }
  return out;
}

//! MINPACK's `enorm`: the Euclidean norm, scaled to avoid over/underflow.
double enorm(int n, const double* x) {
  const double rdwarf = 3.834e-20, rgiant = 1.304e19;
  double s1 = 0.0, s2 = 0.0, s3 = 0.0;
  double x1max = 0.0, x3max = 0.0;
  const double agiant = rgiant / static_cast<double>(n > 0 ? n : 1);
  for (int i = 0; i < n; ++i) {
    const double xabs = std::fabs(x[i]);
    if (xabs > rdwarf && xabs < agiant) {
      s2 += xabs * xabs;  // intermediate components
    } else if (xabs > rdwarf) {
      if (xabs > x1max) {  // large components
        const double ratio = x1max / xabs;
        s1 = 1.0 + s1 * ratio * ratio;
        x1max = xabs;
      } else {
        const double ratio = xabs / x1max;
        s1 += ratio * ratio;
      }
    } else {
      if (xabs > x3max) {  // small components
        const double ratio = x3max / xabs;
        s3 = 1.0 + s3 * ratio * ratio;
        x3max = xabs;
      } else if (xabs != 0.0) {
        const double ratio = xabs / x3max;
        s3 += ratio * ratio;
      }
    }
  }
  if (s1 != 0.0) return x1max * std::sqrt(s1 + (s2 / x1max) / x1max);
  if (s2 != 0.0) {
    if (s2 >= x3max)
      return std::sqrt(s2 * (1.0 + (x3max / s2) * (x3max * s3)));
    return std::sqrt(x3max * ((s2 / x3max) + (x3max * s3)));
  }
  return x3max * std::sqrt(s3);
}

//! MINPACK's `qrfac`: Householder QR with column pivoting.
/*! \param[in,out] a `m x n`, overwritten with the factorisation
    \param[out] ipvt the column permutation, 0-based here (MINPACK's is 1-based)
    \param[out] rdiag the diagonal of R
    \param[out] acnorm the original column norms
    \param[in] wa scratch, length n */
void qrfac(int m, int n, double* a, int lda, int* ipvt, double* rdiag,
           double* acnorm, double* wa) {
  const double p05 = 0.05;
  for (int j = 0; j < n; ++j) {
    acnorm[j] = enorm(m, &a[j * lda]);
    rdiag[j] = acnorm[j];
    wa[j] = rdiag[j];
    ipvt[j] = j;
  }
  const int minmn = std::min(m, n);
  for (int j = 0; j < minmn; ++j) {
    // Bring the column of largest norm into the pivot position.
    int kmax = j;
    for (int k = j; k < n; ++k)
      if (rdiag[k] > rdiag[kmax]) kmax = k;
    if (kmax != j) {
      for (int i = 0; i < m; ++i)
        std::swap(a[j * lda + i], a[kmax * lda + i]);
      rdiag[kmax] = rdiag[j];
      wa[kmax] = wa[j];
      std::swap(ipvt[j], ipvt[kmax]);
    }
    // The Householder transformation reducing column j to a multiple of
    // the j-th unit vector.
    double ajnorm = enorm(m - j, &a[j * lda + j]);
    if (ajnorm != 0.0) {
      if (a[j * lda + j] < 0.0) ajnorm = -ajnorm;
      for (int i = j; i < m; ++i) a[j * lda + i] /= ajnorm;
      a[j * lda + j] += 1.0;
      // Apply it to the remaining columns and update the norms.
      for (int k = j + 1; k < n; ++k) {
        double sum = 0.0;
        for (int i = j; i < m; ++i) sum += a[j * lda + i] * a[k * lda + i];
        const double temp = sum / a[j * lda + j];
        for (int i = j; i < m; ++i) a[k * lda + i] -= temp * a[j * lda + i];
        if (rdiag[k] != 0.0) {
          double t = a[k * lda + j] / rdiag[k];
          rdiag[k] *= std::sqrt(std::max(0.0, 1.0 - t * t));
          t = rdiag[k] / wa[k];
          if (p05 * t * t <= MACHINE_EPS) {
            rdiag[k] = enorm(m - j - 1, &a[k * lda + j + 1]);
            wa[k] = rdiag[k];
          }
        }
      }
    }
    rdiag[j] = -ajnorm;
  }
}

//! MINPACK's `qrsolv`: solve `(R^T R + D^T D) x = R^T Q^T b` by Givens.
void qrsolv(int n, double* r, int ldr, const int* ipvt, const double* diag,
            const double* qtb, double* x, double* sdiag, double* wa) {
  const double p5 = 0.5, p25 = 0.25;
  // Copy R and Q^T b, saving R's diagonal in x.
  for (int j = 0; j < n; ++j) {
    for (int i = j; i < n; ++i) r[j * ldr + i] = r[i * ldr + j];
    x[j] = r[j * ldr + j];
    wa[j] = qtb[j];
  }
  for (int j = 0; j < n; ++j) {
    // Prepare the row of D to be eliminated, locating the diagonal element
    // using the permutation from the QR factorisation.
    const int l = ipvt[j];
    if (diag[l] != 0.0) {
      for (int k = j; k < n; ++k) sdiag[k] = 0.0;
      sdiag[j] = diag[l];
      // The transformations modify only a single element of (Q^T b, 0)
      // beyond the first n, which is initially zero.
      double qtbpj = 0.0;
      for (int k = j; k < n; ++k) {
        if (sdiag[k] == 0.0) continue;
        double sn, cs;
        if (std::fabs(r[k * ldr + k]) < std::fabs(sdiag[k])) {
          const double cotan = r[k * ldr + k] / sdiag[k];
          sn = p5 / std::sqrt(p25 + p25 * cotan * cotan);
          cs = sn * cotan;
        } else {
          const double tn = sdiag[k] / r[k * ldr + k];
          cs = p5 / std::sqrt(p25 + p25 * tn * tn);
          sn = cs * tn;
        }
        r[k * ldr + k] = cs * r[k * ldr + k] + sn * sdiag[k];
        const double temp = cs * wa[k] + sn * qtbpj;
        qtbpj = -sn * wa[k] + cs * qtbpj;
        wa[k] = temp;
        for (int i = k + 1; i < n; ++i) {
          const double t = cs * r[k * ldr + i] + sn * sdiag[i];
          sdiag[i] = -sn * r[k * ldr + i] + cs * sdiag[i];
          r[k * ldr + i] = t;
        }
      }
    }
    // Store S's diagonal element and restore R's.
    sdiag[j] = r[j * ldr + j];
    r[j * ldr + j] = x[j];
  }
  // Solve the triangular system; a singular one gets a least-squares answer.
  int nsing = n;
  for (int j = 0; j < n; ++j) {
    if (sdiag[j] == 0.0 && nsing == n) nsing = j;
    if (nsing < n) wa[j] = 0.0;
  }
  for (int k = nsing - 1; k >= 0; --k) {
    double sum = 0.0;
    for (int i = k + 1; i < nsing; ++i) sum += r[k * ldr + i] * wa[i];
    wa[k] = (wa[k] - sum) / sdiag[k];
  }
  for (int j = 0; j < n; ++j) x[ipvt[j]] = wa[j];
}

//! MINPACK's `lmpar`: the Levenberg-Marquardt parameter and its step.
void lmpar(int n, double* r, int ldr, const int* ipvt, const double* diag,
           const double* qtb, double delta, double* par, double* x,
           double* sdiag, double* wa1, double* wa2) {
  const double p1 = 0.1, p001 = 0.001;
  // The Gauss-Newton direction; a rank-deficient Jacobian gets a
  // least-squares solution.
  int nsing = n;
  for (int j = 0; j < n; ++j) {
    wa1[j] = qtb[j];
    if (r[j * ldr + j] == 0.0 && nsing == n) nsing = j;
    if (nsing < n) wa1[j] = 0.0;
  }
  for (int k = nsing - 1; k >= 0; --k) {
    wa1[k] /= r[k * ldr + k];
    const double temp = wa1[k];
    for (int i = 0; i < k; ++i) wa1[i] -= r[k * ldr + i] * temp;
  }
  for (int j = 0; j < n; ++j) x[ipvt[j]] = wa1[j];
  // Evaluate the function at the origin and test the Gauss-Newton direction.
  int iter = 0;
  for (int j = 0; j < n; ++j) wa2[j] = diag[j] * x[j];
  double dxnorm = enorm(n, wa2);
  double fp = dxnorm - delta;
  if (fp <= p1 * delta) {
    *par = 0.0;
    return;
  }
  // A full-rank Jacobian gives a lower bound for the zero of the function.
  double parl = 0.0;
  if (nsing >= n) {
    for (int j = 0; j < n; ++j) {
      const int l = ipvt[j];
      wa1[j] = diag[l] * (wa2[l] / dxnorm);
    }
    for (int j = 0; j < n; ++j) {
      double sum = 0.0;
      for (int i = 0; i < j; ++i) sum += r[j * ldr + i] * wa1[i];
      wa1[j] = (wa1[j] - sum) / r[j * ldr + j];
    }
    const double temp = enorm(n, wa1);
    parl = ((fp / delta) / temp) / temp;
  }
  // An upper bound for the zero of the function.
  for (int j = 0; j < n; ++j) {
    double sum = 0.0;
    for (int i = 0; i <= j; ++i) sum += r[j * ldr + i] * qtb[i];
    wa1[j] = sum / diag[ipvt[j]];
  }
  const double gnorm = enorm(n, wa1);
  double paru = gnorm / delta;
  if (paru == 0.0) paru = DWARF / std::min(delta, p1);
  // Move an input `par` outside (parl, paru) to the nearer endpoint.
  *par = std::max(*par, parl);
  *par = std::min(*par, paru);
  if (*par == 0.0) *par = gnorm / dxnorm;

  for (;;) {
    ++iter;
    if (*par == 0.0) *par = std::max(DWARF, p001 * paru);
    double temp = std::sqrt(*par);
    for (int j = 0; j < n; ++j) wa1[j] = temp * diag[j];
    qrsolv(n, r, ldr, ipvt, wa1, qtb, x, sdiag, wa2);
    for (int j = 0; j < n; ++j) wa2[j] = diag[j] * x[j];
    dxnorm = enorm(n, wa2);
    temp = fp;
    fp = dxnorm - delta;
    if (std::fabs(fp) <= p1 * delta ||
        (parl == 0.0 && fp <= temp && temp < 0.0) || iter == 10) {
      break;
    }
    // The Newton correction.
    for (int j = 0; j < n; ++j) {
      const int l = ipvt[j];
      wa1[j] = diag[l] * (wa2[l] / dxnorm);
    }
    for (int j = 0; j < n; ++j) {
      wa1[j] /= sdiag[j];
      const double t = wa1[j];
      for (int i = j + 1; i < n; ++i) wa1[i] -= r[j * ldr + i] * t;
    }
    temp = enorm(n, wa1);
    const double parc = ((fp / delta) / temp) / temp;
    if (fp > 0.0) parl = std::max(parl, *par);
    if (fp < 0.0) paru = std::min(paru, *par);
    *par = std::max(parl, *par + parc);
  }
  if (iter == 0) *par = 0.0;
}

//! MINPACK's `covar`: `(J^T J)^-1` from the R of the final QR.
/*! \param[in,out] r on entry the `n x n` upper triangle of the factorisation,
                     on exit the symmetric covariance, permuted back
    \param[in] ipvt the 0-based column permutation
    \param[in] tol the rank tolerance
    \param[in] wa scratch, length n
    \return the rank found; 0 means the whole matrix is singular */
int covar(int n, double* r, int ldr, const int* ipvt, double tol, double* wa) {
  // Invert R in place, in its own upper triangle.
  const double tolr = tol * std::fabs(r[0]);
  int l = 0;
  bool full = true;
  for (int k = 0; k < n && full; ++k) {
    if (std::fabs(r[k * ldr + k]) <= tolr) {
      full = false;
      break;
    }
    r[k * ldr + k] = 1.0 / r[k * ldr + k];
    for (int j = 0; j < k; ++j) {
      const double temp = r[k * ldr + k] * r[k * ldr + j];
      r[k * ldr + j] = 0.0;
      for (int i = 0; i <= j; ++i) r[k * ldr + i] -= temp * r[j * ldr + i];
    }
    l = k + 1;
  }
  // Form the upper triangle of the inverse of R^T R.
  for (int k = 0; k < l; ++k) {
    for (int j = 0; j < k; ++j) {
      const double temp = r[k * ldr + j];
      for (int i = 0; i <= j; ++i) r[j * ldr + i] += temp * r[k * ldr + i];
    }
    const double temp = r[k * ldr + k];
    for (int i = 0; i <= k; ++i) r[k * ldr + i] *= temp;
  }
  // Permute back: the lower triangle of the covariance, plus its diagonal
  // in `wa`.
  for (int j = 0; j < n; ++j) {
    const int jj = ipvt[j];
    const bool sing = j >= l;
    for (int i = 0; i <= j; ++i) {
      if (sing) r[j * ldr + i] = 0.0;
      const int ii = ipvt[i];
      if (ii > jj) r[jj * ldr + ii] = r[j * ldr + i];
      if (ii < jj) r[ii * ldr + jj] = r[j * ldr + i];
    }
    wa[jj] = r[j * ldr + j];
  }
  // Symmetrise.
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i <= j; ++i) r[j * ldr + i] = r[i * ldr + j];
    r[j * ldr + j] = wa[j];
  }
  return l;
}

}  // namespace

// ------------------------------------------------------------ progress bar

int minimizer_expected_evaluations(int n, int maxfev) {
  const int limit = (maxfev > 0) ? maxfev : 200 * (n + 1);
  const int estimate = MINIMIZER_EXPECTED_ITERATIONS * (n + 1);
  return std::max(1, std::min(limit, estimate));
}

int minimizer_reported_total(int nfev, int expected) {
  const double scale =
      std::max(1.0, BAR_HALF_LIFE * static_cast<double>(std::max(1, expected)));
  // The ceiling is folded into the curve rather than clamped on afterwards:
  // a constant fraction over integer totals oscillates, and the reported
  // fraction stepped backwards once every hundred evaluations.
  double fraction =
      MAX_RUNNING_RATIO * (1.0 - std::pow(2.0, -static_cast<double>(nfev) / scale));
  fraction = std::max(fraction, 1e-6);
  const double total = std::ceil(static_cast<double>(nfev) / fraction);
  return std::max(nfev + 1, static_cast<int>(total));
}

// -------------------------------------------------------------- the class

FitMinimizer::FitMinimizer(const std::string& algorithm) { set_algorithm(algorithm); }

FitMinimizer::~FitMinimizer() {}

void FitMinimizer::set_algorithm(const std::string& algorithm) {
  if (algorithm != "leastsq") {
    throw FitMinimizerConfigurationError(
        "FitMinimizer: unknown algorithm '" + algorithm + "'; the only one is "
        "'leastsq' (MINPACK's lmdif with chisurf's bounds transform)");
  }
  algorithm_ = algorithm;
}

const std::string& FitMinimizer::get_algorithm() const { return algorithm_; }

void FitMinimizer::set_parameter_ports(
    const std::vector<std::shared_ptr<GraphPort> >& parameters) {
  for (std::size_t i = 0; i < parameters.size(); ++i) {
    if (!parameters[i])
      throw FitMinimizerConfigurationError(
          "set_parameter_ports: a null port in the list");
    if (parameters[i]->get_fixed())
      throw FitMinimizerConfigurationError(
          "set_parameter_ports: port '" + parameters[i]->get_name() +
          "' is fixed and cannot be optimised");
  }
  parameters_ = parameters;
  configure_from_ports();
}

std::vector<std::shared_ptr<GraphPort> > FitMinimizer::get_parameter_ports() const {
  return parameters_;
}

std::vector<std::string> FitMinimizer::get_parameter_names() const {
  std::vector<std::string> names;
  names.reserve(parameters_.size());
  for (std::size_t i = 0; i < parameters_.size(); ++i)
    names.push_back(parameters_[i]->get_name());
  return names;
}

void FitMinimizer::configure_from_ports() {
  const double inf = std::numeric_limits<double>::infinity();
  ndim_ = static_cast<unsigned int>(parameters_.size());
  initial_values_.assign(ndim_, 0.0);
  lower_.assign(ndim_, -inf);
  upper_.assign(ndim_, inf);
  for (unsigned int i = 0; i < ndim_; ++i) {
    initial_values_[i] = parameters_[i]->get_value();
    if (parameters_[i]->get_is_bounded()) {
      lower_[i] = parameters_[i]->get_lower_bound();
      upper_[i] = parameters_[i]->get_upper_bound();
    }
  }
}

void FitMinimizer::set_initial_values(const std::vector<double>& values) {
  if (ndim_ != 0 && values.size() != ndim_)
    throw FitMinimizerConfigurationError(
        "set_initial_values: as many values as parameters are needed");
  initial_values_ = values;
  if (ndim_ == 0) {
    const double inf = std::numeric_limits<double>::infinity();
    ndim_ = static_cast<unsigned int>(values.size());
    lower_.assign(ndim_, -inf);
    upper_.assign(ndim_, inf);
  }
}

std::vector<double> FitMinimizer::get_initial_values() const {
  return initial_values_;
}

void FitMinimizer::set_bounds(const std::vector<double>& lower,
                           const std::vector<double>& upper) {
  if (lower.size() != upper.size())
    throw FitMinimizerConfigurationError(
        "set_bounds: the lower and upper lists differ in length");
  if (ndim_ != 0 && lower.size() != ndim_)
    throw FitMinimizerConfigurationError(
        "set_bounds: as many bounds as parameters are needed");
  for (std::size_t i = 0; i < lower.size(); ++i) {
    if (!is_unbounded(lower[i]) && !is_unbounded(upper[i]) &&
        !(upper[i] > lower[i])) {
      std::ostringstream m;
      m << "set_bounds: parameter " << i << " has an upper bound (" << upper[i]
        << ") that does not exceed its lower bound (" << lower[i] << ")";
      throw FitMinimizerConfigurationError(m.str());
    }
  }
  lower_ = lower;
  upper_ = upper;
  ndim_ = static_cast<unsigned int>(lower.size());
  if (initial_values_.size() != ndim_) initial_values_.assign(ndim_, 0.0);
}

std::vector<double> FitMinimizer::get_lower_bounds() const { return lower_; }
std::vector<double> FitMinimizer::get_upper_bounds() const { return upper_; }

void FitMinimizer::set_objective(std::shared_ptr<GraphNode> node,
                              const std::string& residual_key) {
  if (!node)
    throw FitMinimizerConfigurationError("set_objective: the node is a null pointer");
  if (!residual_key.empty()) residual_key_ = residual_key;
  if (!node->get_output_port(residual_key_)) {
    throw FitMinimizerConfigurationError(
        "set_objective: node '" + node->get_name() + "' has no output port '" +
        residual_key_ + "' to carry the residual vector");
  }
  objective_node_ = node;
  residual_function_ = nullptr;
}

std::shared_ptr<GraphNode> FitMinimizer::get_objective() const {
  return objective_node_;
}

void FitMinimizer::set_residual_port_key(const std::string& key) {
  if (key.empty())
    throw FitMinimizerConfigurationError(
        "set_residual_port_key: the key is empty");
  residual_key_ = key;
}

const std::string& FitMinimizer::get_residual_port_key() const {
  return residual_key_;
}

void FitMinimizer::set_residual_function(
    std::function<std::vector<double>(const std::vector<double>&)> f) {
  if (!f)
    throw FitMinimizerConfigurationError("set_residual_function: a null function");
  residual_function_ = f;
  objective_node_.reset();
}

bool FitMinimizer::has_objective() const {
  return objective_node_ != nullptr || residual_function_ != nullptr;
}

void FitMinimizer::set_observer(FitMinimizerObserver* observer) {
  observer_ = observer;
}

FitMinimizerObserver* FitMinimizer::get_observer() const { return observer_.get(); }

void FitMinimizer::clear_observer() { observer_ = nullptr; }

void FitMinimizer::set_ftol(double v) { ftol_ = v; }
double FitMinimizer::get_ftol() const { return ftol_; }
void FitMinimizer::set_xtol(double v) { xtol_ = v; }
double FitMinimizer::get_xtol() const { return xtol_; }
void FitMinimizer::set_gtol(double v) { gtol_ = v; }
double FitMinimizer::get_gtol() const { return gtol_; }
void FitMinimizer::set_maxfev(int v) { maxfev_ = v; }
int FitMinimizer::get_maxfev() const { return maxfev_; }
void FitMinimizer::set_epsfcn(double v) { epsfcn_ = v; }
double FitMinimizer::get_epsfcn() const { return epsfcn_; }

void FitMinimizer::set_factor(double v) {
  if (!(v > 0.0))
    throw FitMinimizerConfigurationError("set_factor: the factor must be positive");
  factor_ = v;
}
double FitMinimizer::get_factor() const { return factor_; }

void FitMinimizer::set_diag(const std::vector<double>& diag) {
  for (std::size_t i = 0; i < diag.size(); ++i)
    if (!(diag[i] > 0.0))
      throw FitMinimizerConfigurationError(
          "set_diag: every scale factor must be positive");
  diag_ = diag;
}
std::vector<double> FitMinimizer::get_diag() const { return diag_; }

// ------------------------------------------------------- bounds transform

double FitMinimizer::to_external(double xi, unsigned int i) const {
  const bool lo_free = is_unbounded(lower_[i]);
  const bool up_free = is_unbounded(upper_[i]);
  if (lo_free && up_free) return xi;
  if (up_free) return lower_[i] - 1.0 + std::sqrt(xi * xi + 1.0);
  if (lo_free) return upper_[i] + 1.0 - std::sqrt(xi * xi + 1.0);
  return lower_[i] + ((upper_[i] - lower_[i]) / 2.0) * (std::sin(xi) + 1.0);
}

double FitMinimizer::to_internal(double xe, unsigned int i) const {
  const bool lo_free = is_unbounded(lower_[i]);
  const bool up_free = is_unbounded(upper_[i]);
  if (lo_free && up_free) return xe;
  if (up_free) {
    // A start below its own lower bound would take the square root of a
    // negative number; the port would have clipped it, so clip it here too
    // rather than hand the optimiser a NaN it cannot recover from.
    const double v = std::max(xe, lower_[i]) - lower_[i] + 1.0;
    return std::sqrt(v * v - 1.0);
  }
  if (lo_free) {
    const double v = upper_[i] - std::min(xe, upper_[i]) + 1.0;
    return std::sqrt(v * v - 1.0);
  }
  double t = 2.0 * (xe - lower_[i]) / (upper_[i] - lower_[i]) - 1.0;
  t = std::min(1.0, std::max(-1.0, t));
  return std::asin(t);
}

double FitMinimizer::external_gradient(double xi, unsigned int i) const {
  const bool lo_free = is_unbounded(lower_[i]);
  const bool up_free = is_unbounded(upper_[i]);
  // chisurf's `_internal2external_grad` tests `is None` here while its
  // transform lambdas test `_is_unbounded`, so a parameter declared
  // `(-inf, inf)` takes the two-sided branch there and its gradient comes
  // out `(inf - -inf) * cos(v) / 2` = inf -- an infinite column scaling on
  // a covariance the caller then reads. The two tests agree here.
  if (lo_free && up_free) return 1.0;
  if (up_free) return xi / std::sqrt(xi * xi + 1.0);
  if (lo_free) return -xi / std::sqrt(xi * xi + 1.0);
  return (upper_[i] - lower_[i]) * std::cos(xi) / 2.0;
}

//! PRD-120: the two-sided transform's derivative is `(upper-lower)/2 *
//! cos(xi)`, which for a decade-spanning interval (`ub` defensively left at
//! `1e9`) is on the order of the *half-width* almost everywhere -- huge --
//! except right at the two edges, where it collapses to zero. `fdjac2`'s
//! usual step, `eps * |xi|` in the transformed coordinate, ignores that: for
//! a parameter sitting well inside a wide box (deep in the interior, but
//! close to one edge of `[-pi/2, pi/2]` simply because the box dwarfs the
//! parameter's own scale -- the ICS `A0` case, bound `(0, 1e9)`, value
//! `0.4`) its probe lands kilometres away in physical units, the
//! finite-difference column is dominated by the transform's own curvature
//! rather than the objective's, and the linear model LM builds from it
//! stops predicting anything -- the fit stalls at a wrong, plausible-looking
//! point (measured: chi2r ~608 against ~1.0 unbounded on the ICS 2D-Gaussian
//! recovery, see okf/prds/prd-120.md).
//!
//! The fix asks the question `fdjac2` actually needs -- "how far do I step
//! the *physical* parameter" -- first, in external units, exactly as an
//! unbounded parameter would (`eps * |xe|`), then inverts that through the
//! transform to find the internal step that produces it.
//!
//! That externally-relative step is not uniformly better, though: right at
//! a genuine edge (a parameter whose *true* optimum sits near its bound --
//! a scatter fraction pinned near 0) the *original* `eps * |xi|` step is
//! already the safer one, because the transform's own derivative is
//! shrinking to zero there and the plain internal step shrinks with it; an
//! externally-relative probe instead falls back to a fixed absolute `eps`
//! once `|xe|` is tiny, overshooting into the transform's curvature from the
//! other side. Measured while closing this out: on a TCSPC lifetime fit
//! whose scatter fraction converges near its `(0, 100)` lower bound, always
//! preferring the external step moved chi2r from 3.10 (scipy `trf`
//! agrees to 3.09) to 3.25 -- a regression this PRD's own "answers must not
//! move" gate is designed to catch.
//!
//! So take whichever of the two candidates is the *smaller* step: the
//! conservative choice either way, since a smaller forward-difference probe
//! is never a worse local linear approximation than a larger one for a
//! curve this smooth. Deep in an oversized box (ICS) the internal step is
//! the one that overshoots and the external one wins; right at an edge
//! (TCSPC) the external step is the one that overshoots and the internal
//! one -- already correctly shrunk by the transform -- wins.
//!
//! Deliberately scoped to the two-sided branch only: one-sided
//! (`lo_free` xor `up_free`) and unbounded parameters keep the original
//! `eps * |xi|` step verbatim, since the decade-spanning failure is a
//! two-sided (`sin`) phenomenon -- the one-sided (`sqrt`) transform's
//! derivative does not collapse the same way -- and every already-pinned
//! fit with a one-sided or absent bound must not move.
double FitMinimizer::fdjac2_step(double xi, unsigned int i, double eps) const {
  const bool lo_free = is_unbounded(lower_[i]);
  const bool up_free = is_unbounded(upper_[i]);
  double h = eps * std::fabs(xi);
  if (h == 0.0) h = eps;
  if (lo_free || up_free) return h;  // unchanged: one-sided or unbounded.

  const double xe = to_external(xi, i);
  double h_ext = eps * std::fabs(xe);
  if (h_ext == 0.0) h_ext = eps;
  const double xi_perturbed = to_internal(xe + h_ext, i);
  const double h_candidate = xi_perturbed - xi;
  if (h_candidate == 0.0) return h;
  return (std::fabs(h_candidate) < std::fabs(h)) ? h_candidate : h;
}

// ------------------------------------------------------------- evaluation

std::vector<double> FitMinimizer::evaluate_external(
    const std::vector<double>& xe) {
  if (residual_function_) return residual_function_(xe);
  for (unsigned int i = 0; i < ndim_; ++i) parameters_[i]->set_value(xe[i]);
  objective_node_->update();
  const std::shared_ptr<GraphPort> out =
      objective_node_->get_output_port(residual_key_);
  if (!out) {
    throw FitMinimizerConfigurationError(
        "FitMinimizer: the objective node lost its output port '" + residual_key_ +
        "' while running");
  }
  return out->get_values_ref();
}

bool FitMinimizer::evaluate_internal(const std::vector<double>& xi,
                                  std::vector<double>* fvec) {
  std::vector<double> xe(ndim_);
  for (unsigned int i = 0; i < ndim_; ++i) xe[i] = to_external(xi[i], i);
  std::vector<double> res = evaluate_external(xe);
  if (!fvec->empty() && res.size() != fvec->size()) {
    std::ostringstream m;
    m << "FitMinimizer: the objective returned " << res.size()
      << " residuals after returning " << fvec->size()
      << "; the residual length must not depend on the parameters";
    throw FitMinimizerConfigurationError(m.str());
  }
  *fvec = res;
  ++n_evaluations_;
  if (observer_) {
    double chi2 = 0.0;
    for (std::size_t i = 0; i < fvec->size(); ++i)
      chi2 += (*fvec)[i] * (*fvec)[i];
    const int total = minimizer_reported_total(n_evaluations_,
                                               expected_evaluations_);
    if (!observer_->report(n_evaluations_, total, chi2)) {
      cancelled_ = true;
      return false;
    }
  }
  return true;
}

// -------------------------------------------------------------------- run

void FitMinimizer::validate() const {
  if (ndim_ == 0)
    throw FitMinimizerConfigurationError(
        "no parameters: call set_parameter_ports() or set_initial_values()");
  if (!has_objective())
    throw FitMinimizerConfigurationError(
        "no objective: call set_objective() or set_residual_function()");
  if (initial_values_.size() != ndim_)
    throw FitMinimizerConfigurationError(
        "the starting values and the parameters differ in number");
  if (lower_.size() != ndim_ || upper_.size() != ndim_)
    throw FitMinimizerConfigurationError(
        "the bounds and the parameters differ in number");
  if (!diag_.empty() && diag_.size() != ndim_)
    throw FitMinimizerConfigurationError(
        "the scale factors and the parameters differ in number");
}

void FitMinimizer::compute_objective_batch(double* in_candidates, int n_rows,
                                       int n_cols, double** out_view,
                                       int* n_out_view) {
  if (!has_objective()) {
    throw FitMinimizerConfigurationError(
        "no objective: call set_objective() or set_residual_function()");
  }
  if (ndim_ == 0) {
    throw FitMinimizerConfigurationError(
        "no parameters: call set_parameter_ports() or set_initial_values()");
  }
  if (n_cols != static_cast<int>(ndim_)) {
    std::ostringstream m;
    m << "compute_objective_batch: the candidates have " << n_cols
      << " columns and there are " << ndim_
      << " free parameters; one column per parameter, one row per candidate";
    throw FitMinimizerConfigurationError(m.str());
  }
  const int rows = std::max(0, n_rows);
  if (rows == 0 || in_candidates == nullptr) {
    internal::new_double_view(0, out_view, n_out_view);
    return;
  }

  // What the ports hold now, so a surface scan does not move somebody's fit.
  // Not needed for the residual-function path, which never touches them.
  std::vector<double> saved;
  if (!residual_function_) {
    saved.resize(ndim_);
    for (unsigned int i = 0; i < ndim_; ++i) saved[i] = parameters_[i]->get_value();
  }

  double* out = internal::new_double_view(static_cast<std::size_t>(rows),
                                          out_view, n_out_view);
  if (out == nullptr) return;

  std::vector<double> xe(ndim_);
  for (int r = 0; r < rows; ++r) {
    const double* row = in_candidates + static_cast<std::size_t>(r) * ndim_;
    for (unsigned int i = 0; i < ndim_; ++i) xe[i] = row[i];
    // evaluate_external, not evaluate_internal: the bounds transform belongs
    // to the optimiser's internal coordinates, and a caller scanning a
    // surface gave the values it means, in the parameters' own units.
    const std::vector<double>& res = evaluate_external(xe);
    double chi2 = 0.0;
    for (double v : res) chi2 += v * v;
    out[r] = std::isnan(chi2) ? std::numeric_limits<double>::infinity() : chi2;
  }

  // Put the ports back and leave the graph consistent with them, so the model
  // curve a caller reads off the graph is the one the ports say it is.
  if (!residual_function_) {
    for (unsigned int i = 0; i < ndim_; ++i) parameters_[i]->set_value(saved[i]);
    objective_node_->update();
  }
}

void FitMinimizer::reset() {
  x_.clear();
  fvec_.clear();
  r_.clear();
  ipvt_.clear();
  covariance_.clear();
  covariance_parameters_.clear();
  chi2_ = 0.0;
  n_evaluations_ = 0;
  status_ = 0;
  cancelled_ = false;
}

int FitMinimizer::run() {
  validate();
  reset();
  const int n = static_cast<int>(ndim_);
  expected_evaluations_ = minimizer_expected_evaluations(n, maxfev_);

  std::vector<double> xi(ndim_);
  for (unsigned int i = 0; i < ndim_; ++i)
    xi[i] = to_internal(initial_values_[i], i);

  std::vector<double> fvec;
  status_ = lmdif(&xi, &fvec);

  x_.assign(ndim_, 0.0);
  for (unsigned int i = 0; i < ndim_; ++i) x_[i] = to_external(xi[i], i);
  fvec_ = fvec;
  chi2_ = 0.0;
  for (std::size_t i = 0; i < fvec_.size(); ++i) chi2_ += fvec_[i] * fvec_[i];
  if (std::isnan(chi2_)) chi2_ = std::numeric_limits<double>::infinity();

  // Leave the ports -- and the graph hanging off them -- at the solution, so
  // that a caller reads the fitted curve without evaluating anything again.
  if (!residual_function_) {
    for (unsigned int i = 0; i < ndim_; ++i) parameters_[i]->set_value(x_[i]);
    objective_node_->update();
  }

  // The covariance, from the R the optimiser already has, with the columns
  // taken back to external coordinates. `leastsqbound` only builds this for a
  // converged fit, and so does this: an R from a run that stopped on maxfev
  // describes a point the optimiser was passing through.
  covariance_.clear();
  if (!r_.empty() && (status_ >= 1 && status_ <= 4)) {
    std::vector<double> r = r_;
    for (int j = 0; j < n; ++j) {
      const double g = external_gradient(xi[ipvt_[j]], ipvt_[j]);
      if (!(std::fabs(g) > 0.0)) {
        // A two-sided parameter sitting exactly on a bound has cos(xi) = 0,
        // so its internal derivative vanishes and the external covariance is
        // unbounded. Say so with a non-finite column rather than dividing by
        // zero silently.
        for (int i = 0; i <= j && i < n; ++i)
          r[j * n + i] = std::numeric_limits<double>::infinity();
      } else {
        for (int i = 0; i <= j && i < n; ++i) r[j * n + i] /= g;
      }
    }
    std::vector<double> wa(ndim_, 0.0);
    const int rank = covar(n, &r[0], n, &ipvt_[0], MACHINE_EPS, &wa[0]);
    if (rank > 0) {
      // Row-major on the way out; the matrix is symmetric, so the transpose
      // of MINPACK's column-major answer is the answer.
      covariance_.assign(static_cast<std::size_t>(n) * n, 0.0);
      for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) covariance_[i * n + j] = r[j * n + i];
    }
  }

  if (observer_ && !cancelled_) {
    // A converged fit fills its own bar; `leastsqbound` reports (nfev, nfev).
    observer_->report(n_evaluations_, n_evaluations_, chi2_);
  }
  return status_;
}

//! MINPACK's `lmdif`, over the internal (unconstrained) coordinates.
int FitMinimizer::lmdif(std::vector<double>* xv, std::vector<double>* fvecv) {
  const double p1 = 0.1, p5 = 0.5, p25 = 0.25, p75 = 0.75, p0001 = 1.0e-4;
  const int n = static_cast<int>(ndim_);
  double* x = &(*xv)[0];

  int maxfev = maxfev_;
  if (maxfev <= 0) maxfev = 200 * (n + 1);

  // Evaluate at the starting point; its length is what `m` is.
  std::vector<double> fvec;
  if (!evaluate_internal(*xv, &fvec)) {
    *fvecv = fvec;
    return -1;
  }
  const int m = static_cast<int>(fvec.size());
  if (m < n) {
    std::ostringstream msg;
    msg << "FitMinimizer: " << n << " parameters against " << m
        << " residuals; a least-squares fit needs at least as many residuals "
           "as parameters";
    throw FitMinimizerConfigurationError(msg.str());
  }
  const int ldfjac = m;

  const bool internal_scaling = diag_.empty();
  std::vector<double> diag = internal_scaling ? std::vector<double>(n, 1.0)
                                              : diag_;
  std::vector<double> fjac(static_cast<std::size_t>(ldfjac) * n, 0.0);
  std::vector<double> qtf(n, 0.0);
  std::vector<double> wa1(n), wa2(n), wa3(n), wa4(m);
  std::vector<double> trial(n);
  ipvt_.assign(n, 0);

  double fnorm = enorm(m, &fvec[0]);
  double par = 0.0, delta = 0.0, xnorm = 0.0;
  int iter = 1;
  int info = 0;

  for (;;) {  // outer loop
    // The Jacobian, by forward differences (MINPACK's `fdjac2`).
    {
      const double eps = std::sqrt(std::max(epsfcn_, MACHINE_EPS));
      for (int j = 0; j < n; ++j) {
        const double temp = x[j];
        // PRD-120: chosen so the probe is a sane relative step in the
        // *physical* parameter, not the transformed coordinate -- see
        // fdjac2_step's docstring.
        const double h = fdjac2_step(temp, static_cast<unsigned int>(j), eps);
        x[j] = temp + h;
        std::vector<double> f = fvec;
        if (!evaluate_internal(*xv, &f)) {
          x[j] = temp;
          *fvecv = fvec;
          return -1;
        }
        x[j] = temp;
        for (int i = 0; i < m; ++i)
          fjac[static_cast<std::size_t>(j) * ldfjac + i] =
              (f[i] - fvec[i]) / h;
      }
    }

    qrfac(m, n, &fjac[0], ldfjac, &ipvt_[0], &wa1[0], &wa2[0], &wa3[0]);

    if (iter == 1) {
      if (internal_scaling) {
        for (int j = 0; j < n; ++j)
          diag[j] = (wa2[j] == 0.0) ? 1.0 : wa2[j];
      }
      for (int j = 0; j < n; ++j) wa3[j] = diag[j] * x[j];
      xnorm = enorm(n, &wa3[0]);
      delta = factor_ * xnorm;
      if (delta == 0.0) delta = factor_;
    }

    // Q^T fvec, its first n components into qtf.
    for (int i = 0; i < m; ++i) wa4[i] = fvec[i];
    for (int j = 0; j < n; ++j) {
      const std::size_t col = static_cast<std::size_t>(j) * ldfjac;
      if (fjac[col + j] != 0.0) {
        double sum = 0.0;
        for (int i = j; i < m; ++i) sum += fjac[col + i] * wa4[i];
        const double temp = -sum / fjac[col + j];
        for (int i = j; i < m; ++i) wa4[i] += fjac[col + i] * temp;
      }
      fjac[col + j] = wa1[j];
      qtf[j] = wa4[j];
    }

    // The norm of the scaled gradient.
    double gnorm = 0.0;
    if (fnorm != 0.0) {
      for (int j = 0; j < n; ++j) {
        const int l = ipvt_[j];
        if (wa2[l] == 0.0) continue;
        double sum = 0.0;
        for (int i = 0; i <= j; ++i)
          sum += fjac[static_cast<std::size_t>(j) * ldfjac + i] *
                 (qtf[i] / fnorm);
        gnorm = std::max(gnorm, std::fabs(sum / wa2[l]));
      }
    }
    if (gnorm <= gtol_) info = 4;
    if (info != 0) break;

    if (internal_scaling)
      for (int j = 0; j < n; ++j) diag[j] = std::max(diag[j], wa2[j]);

    double actred = 0.0, prered = 0.0, ratio = 0.0;
    for (;;) {  // inner loop
      lmpar(n, &fjac[0], ldfjac, &ipvt_[0], &diag[0], &qtf[0], delta, &par,
            &wa1[0], &wa2[0], &wa3[0], &wa4[0]);
      // The direction p, the trial point x + p, and the norm of p.
      for (int j = 0; j < n; ++j) {
        wa1[j] = -wa1[j];
        wa2[j] = x[j] + wa1[j];
        wa3[j] = diag[j] * wa1[j];
      }
      const double pnorm = enorm(n, &wa3[0]);
      if (iter == 1) delta = std::min(delta, pnorm);

      for (int j = 0; j < n; ++j) trial[j] = wa2[j];
      std::vector<double> ftrial = fvec;
      if (!evaluate_internal(trial, &ftrial)) {
        *fvecv = fvec;
        return -1;
      }
      const double fnorm1 = enorm(m, &ftrial[0]);

      actred = -1.0;
      if (p1 * fnorm1 < fnorm) {
        const double t = fnorm1 / fnorm;
        actred = 1.0 - t * t;
      }
      // The scaled predicted reduction and directional derivative.
      for (int j = 0; j < n; ++j) wa3[j] = 0.0;
      for (int j = 0; j < n; ++j) {
        const double temp = wa1[ipvt_[j]];
        for (int i = 0; i <= j; ++i)
          wa3[i] += fjac[static_cast<std::size_t>(j) * ldfjac + i] * temp;
      }
      const double temp1 = enorm(n, &wa3[0]) / fnorm;
      const double temp2 = (std::sqrt(par) * pnorm) / fnorm;
      prered = temp1 * temp1 + temp2 * temp2 / p5;
      const double dirder = -(temp1 * temp1 + temp2 * temp2);
      ratio = (prered != 0.0) ? actred / prered : 0.0;

      // Update the step bound.
      if (ratio <= p25) {
        double temp;
        if (actred >= 0.0) {
          temp = p5;
        } else {
          temp = p5 * dirder / (dirder + p5 * actred);
        }
        if (p1 * fnorm1 >= fnorm || temp < p1) temp = p1;
        delta = temp * std::min(delta, pnorm / p1);
        par /= temp;
      } else if (par == 0.0 || ratio >= p75) {
        delta = pnorm / p5;
        par = p5 * par;
      }

      if (ratio >= p0001) {  // successful iteration
        for (int j = 0; j < n; ++j) {
          x[j] = wa2[j];
          wa2[j] = diag[j] * x[j];
        }
        fvec = ftrial;
        xnorm = enorm(n, &wa2[0]);
        fnorm = fnorm1;
        ++iter;
      }

      // Convergence.
      if (std::fabs(actred) <= ftol_ && prered <= ftol_ && p5 * ratio <= 1.0)
        info = 1;
      if (delta <= xtol_ * xnorm) info = 2;
      if (std::fabs(actred) <= ftol_ && prered <= ftol_ &&
          p5 * ratio <= 1.0 && info == 2)
        info = 3;
      if (info != 0) break;
      // Termination, and tolerances too stringent to make progress.
      if (n_evaluations_ >= maxfev) info = 5;
      if (std::fabs(actred) <= MACHINE_EPS && prered <= MACHINE_EPS &&
          p5 * ratio <= 1.0)
        info = 6;
      if (delta <= MACHINE_EPS * xnorm) info = 7;
      if (gnorm <= MACHINE_EPS) info = 8;
      if (info != 0) break;
      if (ratio >= p0001) break;  // repeat the outer loop
    }
    if (info != 0) break;
  }

  // R, in external-free form: the upper n x n triangle of fjac, which is
  // what `covar` consumes.
  r_.assign(static_cast<std::size_t>(n) * n, 0.0);
  for (int j = 0; j < n; ++j)
    for (int i = 0; i <= j && i < n; ++i)
      r_[static_cast<std::size_t>(j) * n + i] =
          fjac[static_cast<std::size_t>(j) * ldfjac + i];

  *fvecv = fvec;
  return info;
}

// ---------------------------------------------------------------- results

std::vector<double> FitMinimizer::get_x() const { return x_; }
std::vector<double> FitMinimizer::get_residuals() const { return fvec_; }
double FitMinimizer::get_chi2() const { return chi2_; }

double FitMinimizer::get_chi2r() const {
  const double dof = static_cast<double>(fvec_.size()) -
                     static_cast<double>(ndim_) - 1.0;
  return chi2_ / dof;
}

int FitMinimizer::get_number_of_evaluations() const { return n_evaluations_; }
int FitMinimizer::get_status() const { return status_; }
bool FitMinimizer::get_cancelled() const { return cancelled_; }

std::string FitMinimizer::get_message() const {
  std::ostringstream m;
  switch (status_) {
    case -1:
      return "The optimisation was cancelled by its observer.";
    case 0:
      return "Improper input parameters.";
    case 1:
      m << "Both actual and predicted relative reductions in the sum of "
           "squares are at most " << ftol_;
      return m.str();
    case 2:
      m << "The relative error between two consecutive iterates is at most "
        << xtol_;
      return m.str();
    case 3:
      m << "Both actual and predicted relative reductions in the sum of "
           "squares are at most " << ftol_
        << " and the relative error between two consecutive iterates is at "
           "most " << xtol_;
      return m.str();
    case 4:
      m << "The cosine of the angle between the residuals and any column of "
           "the Jacobian is at most " << gtol_ << " in absolute value";
      return m.str();
    case 5:
      m << "Number of calls to the objective has reached maxfev = "
        << (maxfev_ > 0 ? maxfev_ : 200 * (static_cast<int>(ndim_) + 1)) << ".";
      return m.str();
    case 6:
      m << "ftol=" << ftol_
        << " is too small, no further reduction in the sum of squares is "
           "possible.";
      return m.str();
    case 7:
      m << "xtol=" << xtol_
        << " is too small, no further improvement in the approximate solution "
           "is possible.";
      return m.str();
    case 8:
      m << "gtol=" << gtol_
        << " is too small, the residuals are orthogonal to the columns of the "
           "Jacobian to machine precision.";
      return m.str();
    default:
      return "Unknown error.";
  }
}

std::vector<double> FitMinimizer::get_covariance() const { return covariance_; }

std::vector<double> FitMinimizer::get_errors() const {
  std::vector<double> errors;
  if (covariance_.empty()) return errors;
  const std::size_t n = ndim_;
  errors.assign(n, std::numeric_limits<double>::quiet_NaN());
  for (std::size_t i = 0; i < n; ++i) {
    const double v = covariance_[i * n + i];
    errors[i] = (v >= 0.0) ? std::sqrt(v)
                           : std::numeric_limits<double>::quiet_NaN();
  }
  return errors;
}

std::vector<double> FitMinimizer::jacobian_impl(const std::vector<double>& x,
                                             const std::vector<double>& f0,
                                             double epsilon, double floor) {
  const std::size_t m = f0.size();
  const int n = static_cast<int>(ndim_);
  std::vector<double> jac(static_cast<std::size_t>(n) * m, 0.0);
  std::vector<double> xp = x;
  for (int k = 0; k < n; ++k) {
    double step = epsilon * std::max(std::fabs(x[k]), floor);
    // Round the step to an exactly representable difference, so the divisor
    // below is the step the objective actually saw. `approx_grad` does the
    // same, and without it the gradient is wrong in the last digits for
    // every parameter whose binary expansion does not end where the step
    // begins -- which is most of them.
    step = (x[k] + step) - x[k];
    if (step == 0.0) continue;
    xp[k] = x[k] + step;
    const std::vector<double> f = evaluate_external(xp);
    xp[k] = x[k];
    if (f.size() != m) {
      std::ostringstream msg;
      msg << "FitMinimizer: the objective returned " << f.size()
          << " residuals after returning " << m
          << "; the residual length must not depend on the parameters";
      throw FitMinimizerConfigurationError(msg.str());
    }
    for (std::size_t i = 0; i < m; ++i)
      jac[k * m + i] = (f[i] - f0[i]) / step;
  }
  // Put the ports back where they were, and the graph with them: the last
  // thing evaluated above was a perturbation, and a caller reading the model
  // curve off an output port would otherwise get it. `approx_grad` restores
  // for the same reason, and learned it the hard way -- a cache keyed on
  // "did a value change" then serves residuals belonging to a different
  // parameter vector.
  evaluate_external(x);
  return jac;
}

void FitMinimizer::covariance_defaults(double* epsilon, double* floor) {
  // chisurf's `FINITE_DIFFERENCE_STEP` and `approx_grad`'s `max(|x|, 1.0)`.
  // Defaulted here rather than at the call site so that a caller who asks
  // for "the covariance" gets the one the error bars are pinned to.
  if (!(*epsilon > 0.0)) *epsilon = std::sqrt(MACHINE_EPS);
  if (!(*floor > 0.0)) *floor = 1.0;
}

std::vector<double> FitMinimizer::compute_jacobian(const std::vector<double>& x,
                                                double epsilon, double floor) {
  if (ndim_ == 0 || !has_objective() || x.size() != ndim_)
    return std::vector<double>();
  covariance_defaults(&epsilon, &floor);
  const std::vector<double> f0 = evaluate_external(x);
  if (f0.empty()) return std::vector<double>();
  return jacobian_impl(x, f0, epsilon, floor);
}

std::vector<double> FitMinimizer::covariance_from_jacobian(
    const std::vector<double>& jac, std::size_t m) {
  const int n = static_cast<int>(ndim_);
  covariance_parameters_.clear();
  std::vector<int> kept;
  for (int k = 0; k < n; ++k) {
    double ss = 0.0;
    for (std::size_t i = 0; i < m; ++i) {
      const double v = jac[k * m + i];
      ss += v * v;
    }
    // Exactly `covariance_matrix`'s `important_parameters` test. A row of
    // zeros is a parameter that does not move the objective -- chisurf's
    // `E_FRET` is one by construction -- and it is dropped rather than
    // inverted; see the header.
    if (ss > 0.0) kept.push_back(k);
  }
  const int p = static_cast<int>(kept.size());
  if (p == 0) return std::vector<double>();
  covariance_parameters_ = kept;

  std::vector<double> alpha(static_cast<std::size_t>(p) * p, 0.0);
  for (int a = 0; a < p; ++a) {
    for (int b = a; b < p; ++b) {
      double s = 0.0;
      const double* ja = &jac[kept[a] * m];
      const double* jb = &jac[kept[b] * m];
      for (std::size_t i = 0; i < m; ++i) s += ja[i] * jb[i];
      alpha[a * p + b] = s;
      alpha[b * p + a] = s;
    }
  }
  return pseudo_inverse_symmetric(alpha, p);
}

std::vector<double> FitMinimizer::compute_covariance_at(
    const std::vector<double>& x, double epsilon, double floor) {
  covariance_parameters_.clear();
  if (ndim_ == 0 || !has_objective() || x.size() != ndim_)
    return std::vector<double>();
  covariance_defaults(&epsilon, &floor);
  const std::vector<double> f0 = evaluate_external(x);
  if (f0.empty()) return std::vector<double>();
  return covariance_from_jacobian(jacobian_impl(x, f0, epsilon, floor),
                                  f0.size());
}

std::vector<double> FitMinimizer::compute_covariance(double epsilon,
                                                  double floor) {
  covariance_parameters_.clear();
  if (ndim_ == 0 || !has_objective() || x_.size() != ndim_)
    return std::vector<double>();
  covariance_defaults(&epsilon, &floor);

  // `run()` leaves the ports at the solution and the graph evaluated there,
  // so `fvec_` *is* f0 -- the same numbers, not an approximation of them.
  // That is one of the p + 2 evaluations saved, and it is only sound because
  // nothing between `run()` and here may have written a port. A caller who
  // cannot promise that wants compute_covariance_at() instead.
  std::vector<double> f0 = fvec_;
  if (f0.empty()) f0 = evaluate_external(x_);
  if (f0.empty()) return std::vector<double>();
  return covariance_from_jacobian(jacobian_impl(x_, f0, epsilon, floor),
                                  f0.size());
}

std::vector<int> FitMinimizer::get_covariance_parameters() const {
  return covariance_parameters_;
}

std::string FitMinimizer::describe() const {
  std::ostringstream out;
  out << "algorithm      : " << algorithm_ << "\n"
      << "parameters     : " << ndim_ << "\n"
      << "residuals      : " << fvec_.size() << "\n"
      << "evaluations    : " << n_evaluations_ << "\n"
      << "status         : " << status_ << " (" << get_message() << ")\n"
      << "chi2           : " << chi2_ << "\n";
  return out.str();
}

IMPBFF_END_NAMESPACE
