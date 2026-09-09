/**
 * \file IMP/bff/LineSearch.h
 * \brief Given a curvature and a gradient, find a step that actually improves
 *        the objective.
 *
 * Copyright 2007-2023 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_LINESEARCH_H
#define IMPBFF_LINESEARCH_H

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! \name Damped Newton steps
//! @{

/**
 * \brief In-place Cholesky solve of `A x = b`, `A` symmetric, row-major.
 *
 * Returns false if `A` is not positive definite, which is information rather
 * than an error: a Newton or scoring matrix that fails here is telling the
 * caller the quadratic model is not a bowl at this point, and the damping
 * loop below responds by adding to the diagonal until it is.
 */
/**
 * \brief A Cholesky factor kept, so `n` right-hand sides cost one factorisation.
 *
 * `cholesky_solve` factorises on every call, which is right for a stepper that
 * changes its matrix every time and wrong for anything that solves repeatedly
 * against the same curvature -- inverting an `n x n` information matrix column
 * by column, for instance, where it turns `O(n^3)` into `O(n^4)`.
 */
class CholeskyFactor {
 public:
  //! Factor `A` (`n x n`, row-major, symmetric). False if it is not positive definite.
  bool factor(const double* A, std::size_t n_) {
    n = n_;
    L.assign(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j <= i; ++j) {
        double s = A[i * n + j];
        for (std::size_t k = 0; k < j; ++k) s -= L[i * n + k] * L[j * n + k];
        if (i == j) {
          if (!(s > 0.0) || !std::isfinite(s)) return false;
          L[i * n + j] = std::sqrt(s);
        } else {
          L[i * n + j] = s / L[j * n + j];
        }
      }
    }
    ok = true;
    return true;
  }
  //! Solve `A x = b` with the kept factor.
  bool solve(const double* b, double* x) const {
    if (!ok) return false;
    std::vector<double> y(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
      double s = b[i];
      for (std::size_t k = 0; k < i; ++k) s -= L[i * n + k] * y[k];
      y[i] = s / L[i * n + i];
    }
    for (std::size_t ii = n; ii-- > 0;) {
      double s = y[ii];
      for (std::size_t k = ii + 1; k < n; ++k) s -= L[k * n + ii] * x[k];
      x[ii] = s / L[ii * n + ii];
      if (ii == 0) break;
    }
    for (std::size_t i = 0; i < n; ++i) if (!std::isfinite(x[i])) return false;
    return true;
  }
  //! `sum log L_ii`, so `log det A = 2 *` this.
  double log_det_half() const {
    double s = 0.0;
    for (std::size_t i = 0; i < n; ++i) s += std::log(L[i * n + i]);
    return s;
  }

 private:
  std::size_t n = 0;
  bool ok = false;
  std::vector<double> L;
};

inline bool cholesky_solve(const double* A, const double* b, std::size_t n, double* x) {
  std::vector<double> L(n * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      double s = A[i * n + j];
      for (std::size_t k = 0; k < j; ++k) s -= L[i * n + k] * L[j * n + k];
      if (i == j) {
        if (!(s > 0.0) || !std::isfinite(s)) return false;
        L[i * n + j] = std::sqrt(s);
      } else {
        L[i * n + j] = s / L[j * n + j];
      }
    }
  }
  std::vector<double> y(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    double s = b[i];
    for (std::size_t k = 0; k < i; ++k) s -= L[i * n + k] * y[k];
    y[i] = s / L[i * n + i];
  }
  for (std::size_t ii = n; ii-- > 0;) {
    double s = y[ii];
    for (std::size_t k = ii + 1; k < n; ++k) s -= L[k * n + ii] * x[k];
    x[ii] = s / L[ii * n + ii];
    if (ii == 0) break;
  }
  for (std::size_t i = 0; i < n; ++i) if (!std::isfinite(x[i])) return false;
  return true;
}

//! What one step attempt did.
struct StepResult {
  bool accepted = false;        //!< the objective decreased
  double f_new = 0.0;           //!< its value at the accepted point
  double decrement = 0.0;       //!< `grad . step`, the scale-free convergence measure
  double mu = 0.0;              //!< the damping the accepted step needed
  int n_eval = 0;               //!< objective evaluations spent
  double alpha = 0.0;           //!< the fraction of the damped step actually taken
  bool gradient_fallback = false;  //!< the curvature direction failed and a gradient step was used
};

/**
 * \brief A damped Newton (or Fisher scoring) step with a gradient fallback.
 *
 * **The step.** Solve `(A + mu diag|A|) s = g` and move to `theta + s`, where
 * `A` is minus the curvature of the objective and `g` its gradient. With
 * `mu = 0` this is Newton's method; the damping interpolates toward a scaled
 * gradient step, which is Levenberg-Marquardt's idea and the reason it
 * survives a curvature that is wrong, indefinite or merely far from the mode.
 *
 * **The schedule.** Multiply `mu` by `mu_up` whenever the solve fails or the
 * objective does not decrease; divide by `mu_down` on success, with a floor.
 * Growing faster than shrinking is deliberate: a rejected step means the
 * quadratic model is untrustworthy right now and the cost of finding that out
 * again is another objective evaluation, while a successful step costs
 * nothing to be slightly over-damped.
 *
 * **Backtrack along the direction before re-damping it.** When the full step
 * overshoots, `theta + alpha s` for `alpha < 1` is tried first, accepted on
 * Armijo's sufficient-decrease condition (Nocedal & Wright, *Numerical
 * Optimization* 2nd ed., §3.1, Alg. 3.1), and only if the whole ray fails does
 * `mu` grow. The two shrink the step in different ways and it matters which:
 * raising `mu` shrinks EVERY direction, including the ones the quadratic model
 * describes perfectly, while backtracking keeps the Newton direction and moves
 * less far along it. On a Poisson model whose mean is exponential in some
 * coordinates and linear in others -- the usual case -- the first behaviour
 * costs an order of magnitude in iterations. Measured on the s89 posterior
 * (114 coordinates, 1694 dof): the pure damping schedule reached its stopping
 * rule in 50 steps and 122 objective evaluations, this one in 15 and 27.
 *
 * **The fallback matters more than it looks.** If no damping makes the
 * curvature direction descend, the code tries a plain backtracking gradient
 * step. If THAT fails too, the gradient and the objective disagree -- a
 * derivative bug, or an objective that is not differentiable where it is being
 * evaluated -- and the caller should report a failure rather than a converged
 * fit. Silently returning the current point is how a wrong derivative gets
 * reported as an answer.
 *
 * **The decrement, not the gradient norm, is the convergence measure.**
 * `grad . step` is the predicted improvement in the objective's own units, so
 * it is comparable across parameters with wildly different scales, which a
 * gradient norm is not. It is also what distinguishes a fit that has stalled
 * from one that is simply approaching slowly.
 *
 * `Objective` is any callable `double(const double* theta)`. A non-finite
 * value is treated as a rejection, so a caller may return infinity for a
 * point outside the model's domain rather than guarding every step.
 */
template <typename Objective>
class DampedNewton {
 public:
  double mu = 1e-3;             //!< current damping, carried between steps
  double mu_up = 4.0;           //!< growth on rejection
  double mu_down = 10.0;        //!< shrink on acceptance
  double mu_min = 1e-12;        //!< floor
  int max_damping = 25;         //!< attempts before the gradient fallback
  int max_backtrack = 40;       //!< halvings of the gradient step
  int max_line = 12;            //!< halvings of the damped step, before re-damping
  double c1 = 1e-4;             //!< Armijo's sufficient-decrease constant

  /**
   * \brief One step. On acceptance `theta` is advanced in place.
   *
   * \param A minus the curvature, `n x n` row-major and symmetric
   * \param grad the gradient of the objective being MAXIMISED
   * \param f_old its current value, as a quantity being MINIMISED (`-log p`)
   * \param theta the parameters, advanced in place on acceptance
   */
  StepResult step(const double* A, const double* grad, std::size_t n,
                  double f_old, double* theta, Objective f) {
    StepResult r;
    std::vector<double> damped(n * n), s(n, 0.0), cand(n);
    std::vector<double> diag(n);
    for (std::size_t i = 0; i < n; ++i)
      diag[i] = std::max(std::fabs(A[i * n + i]), 1e-10);

    for (int attempt = 0; attempt < max_damping; ++attempt) {
      for (std::size_t i = 0; i < n * n; ++i) damped[i] = A[i];
      for (std::size_t i = 0; i < n; ++i) damped[i * n + i] += mu * diag[i];
      if (!cholesky_solve(damped.data(), grad, n, s.data())) { mu *= mu_up; continue; }
      const double gs = dot(grad, s.data(), n);      // predicted improvement
      if (!(gs > 0.0)) { mu *= mu_up; continue; }    // not an ascent direction
      double alpha = 1.0;
      for (int k = 0; k < max_line; ++k) {
        for (std::size_t i = 0; i < n; ++i) cand[i] = theta[i] + alpha * s[i];
        const double f_new = f(cand.data());
        ++r.n_eval;
        if (std::isfinite(f_new) && f_new <= f_old - c1 * alpha * gs) {
          for (std::size_t i = 0; i < n; ++i) theta[i] = cand[i];
          r.accepted = true; r.f_new = f_new;
          r.decrement = alpha * gs;
          r.alpha = alpha;
          //  a full step means the quadratic model was good; a backtracked one
          //  means it was not, and the damping should not be relaxed for it
          if (alpha == 1.0) mu = std::max(mu / mu_down, mu_min);
          r.mu = mu;
          return r;
        }
        //  Halving, and not the textbook quadratic interpolation of Nocedal &
        //  Wright Sec. 3.5, which was tried and is WORSE here: the interpolated
        //  minimiser hit its own lower safeguard on every backtrack, whatever
        //  the upper one was set to, because the objective along the ray has a
        //  barrier near the full step that a quadratic cannot represent. The
        //  grid it is used on went from 51 steps to 97 (measured 2026-09-09).
        alpha *= 0.5;
      }
      mu *= mu_up;
    }

    //  the curvature direction failed: a plain gradient step must still descend
    double gn = 0.0;
    for (std::size_t i = 0; i < n; ++i) gn += grad[i] * grad[i];
    gn = std::sqrt(gn);
    double t = 1.0 / std::max(gn, 1e-12);
    for (int k = 0; k < max_backtrack; ++k) {
      for (std::size_t i = 0; i < n; ++i) cand[i] = theta[i] + t * grad[i];
      const double f_new = f(cand.data());
      ++r.n_eval;
      if (std::isfinite(f_new) && f_new < f_old) {
        for (std::size_t i = 0; i < n; ++i) { s[i] = t * grad[i]; theta[i] = cand[i]; }
        r.accepted = true; r.f_new = f_new; r.gradient_fallback = true;
        r.decrement = dot(grad, s.data(), n);
        r.mu = mu;
        return r;
      }
      t *= 0.5;
    }
    r.mu = mu;
    return r;                   //!< not accepted: the caller must not call this converged
  }

 private:
  static double dot(const double* a, const double* b, std::size_t n) {
    double s = 0.0;
    for (std::size_t i = 0; i < n; ++i) s += a[i] * b[i];
    return s;
  }
};

//! @}

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_LINESEARCH_H
