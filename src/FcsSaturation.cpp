/**
 *  \file FcsSaturation.cpp
 *  \brief Saturated FCS forward model: steady-state photophysics on a grid.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/FcsSaturation.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <complex>
#include <cstdlib>
#include <map>
#include <mutex>
#include <stdexcept>

IMPBFF_BEGIN_NAMESPACE

namespace {

const double kSatPi = 3.14159265358979323846;
const double kSatPlanck = 6.62607015e-34;
const double kSatCLight = 2.99792458e8;
const double kSatAvogadro = 6.02214076e23;

//! The (K_dark, K_exc) generators with derived diagonals: zero the diagonal,
//! set it to -sum(column), so K = K_d + k_exc * K_e is a proper generator.
void generator_matrices(const std::vector<double>& dark,
                        const std::vector<double>& exc, int n,
                        Eigen::MatrixXd* k_d, Eigen::MatrixXd* k_e) {
  if (n < 1 || dark.size() != static_cast<std::size_t>(n) * n ||
      exc.size() != static_cast<std::size_t>(n) * n) {
    throw std::invalid_argument(
        "fcs saturation: the scheme matrices must be flat n_states^2");
  }
  *k_d = Eigen::Map<const Eigen::Matrix<
      double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
      dark.data(), n, n);
  *k_e = Eigen::Map<const Eigen::Matrix<
      double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
      exc.data(), n, n);
  k_d->diagonal().setZero();
  k_e->diagonal().setZero();
  k_d->diagonal() = -k_d->colwise().sum();
  k_e->diagonal() = -k_e->colwise().sum();
}

//! The normalised stationary distribution of K, or false when degenerate.
bool stationary(const Eigen::MatrixXd& K, Eigen::VectorXd* p_eq) {
  const int n = static_cast<int>(K.rows());
  Eigen::MatrixXd K_eq = K;
  K_eq.row(n - 1).setOnes();
  Eigen::VectorXd b = Eigen::VectorXd::Zero(n);
  b(n - 1) = 1.0;
  Eigen::PartialPivLU<Eigen::MatrixXd> lu(K_eq);
  Eigen::VectorXd p = lu.solve(b);
  if (!p.allFinite()) {
    p = K_eq.colPivHouseholderQr().solve(b);
    if (!p.allFinite()) return false;
  }
  p = p.cwiseMax(0.0);
  const double total = p.sum();
  if (!(total > 0.0)) return false;
  *p_eq = p / total;
  return true;
}

//! Uniform-grid trapezoid over the leading (r) axis with the 2*pi*r element,
//! then over z: the cylindrical volume integral both integrals below take.
double cyl_integral(const Eigen::MatrixXd& f, const std::vector<double>& r,
                    double dr, double dz) {
  const int nr = static_cast<int>(f.rows());
  const int nz = static_cast<int>(f.cols());
  std::vector<double> along_z(nz, 0.0);
  for (int j = 0; j < nz; ++j) {
    double s = 0.0;
    for (int i = 0; i < nr; ++i) {
      const double w = (i == 0 || i == nr - 1) ? 0.5 : 1.0;
      s += w * f(i, j) * 2.0 * kSatPi * r[static_cast<std::size_t>(i)];
    }
    along_z[static_cast<std::size_t>(j)] = s * dr;
  }
  double out = 0.0;
  for (int j = 0; j < nz; ++j) {
    const double w = (j == 0 || j == nz - 1) ? 0.5 : 1.0;
    out += w * along_z[static_cast<std::size_t>(j)];
  }
  return out * dz;
}

//! The cached 0th-order Hankel quadrature matrix 2*pi*r*J0(kr*r)*dr —
//! thousands of Bessel evaluations that depend on nothing but the grids,
//! exactly the Python lru_cache's contract (and its performance).
const Eigen::MatrixXd& hankel_matrix(int n_r, double r_max, double kr_max,
                                     int n_kr) {
  struct Key {
    int n_r; double r_max; double kr_max; int n_kr;
    bool operator<(const Key& o) const {
      if (n_r != o.n_r) return n_r < o.n_r;
      if (r_max != o.r_max) return r_max < o.r_max;
      if (kr_max != o.kr_max) return kr_max < o.kr_max;
      return n_kr < o.n_kr;
    }
  };
  static std::map<Key, Eigen::MatrixXd> cache;
  static std::mutex mutex;
  std::lock_guard<std::mutex> lock(mutex);
  const Key key{n_r, r_max, kr_max, n_kr};
  auto it = cache.find(key);
  if (it != cache.end()) return it->second;
  if (cache.size() > 16) cache.clear();

  Eigen::MatrixXd m(n_kr, n_r);
  const double dr = n_r > 1 ? r_max / (n_r - 1) : 1.0;
  for (int a = 0; a < n_kr; ++a) {
    const double kr = kr_max * a / (n_kr > 1 ? n_kr - 1 : 1);
    for (int i = 0; i < n_r; ++i) {
      const double r = r_max * i / (n_r > 1 ? n_r - 1 : 1);
      const double w = (i == 0 || i == n_r - 1) ? 0.5 : 1.0;
      m(a, i) = 2.0 * kSatPi * r * w * dr * bessel_j0(kr * r);
    }
  }
  return cache.emplace(key, std::move(m)).first->second;
}

double* publish(const std::vector<double>& values, double** out_view,
                int* n_out_view) {
  double* out = static_cast<double*>(
      std::malloc(sizeof(double) * values.size()));
  if (out == nullptr) throw std::bad_alloc();
  for (std::size_t i = 0; i < values.size(); ++i) out[i] = values[i];
  *out_view = out;
  *n_out_view = static_cast<int>(values.size());
  return out;
}

}  // namespace

double bessel_j0(double x) {
  // Periodic trapezoid on (1/pi) * int_0^pi cos(x sin t) dt: geometric
  // convergence for a smooth periodic integrand; 64 nodes give ~1e-15 up
  // to |x| ~ 40, which is beyond the kr_max * r_max = 30 this model uses.
  const int n = 64;
  double s = 0.5 * (std::cos(0.0) + std::cos(x * std::sin(kSatPi)));
  for (int i = 1; i < n; ++i) {
    s += std::cos(x * std::sin(kSatPi * i / n));
  }
  return s / n;
}

double fcs_excitation_rate_peak(double power_w, double extinction, double w0,
                                double wavelength_m) {
  const double photon_energy = kSatPlanck * kSatCLight / wavelength_m;
  const double flux = power_w / photon_energy;
  const double sigma_cm2 = 1000.0 * std::log(10.0) * extinction / kSatAvogadro;
  const double sigma_m2 = sigma_cm2 * 1e-4;
  return sigma_m2 * 2.0 * flux / (kSatPi * w0 * w0);
}

void fcs_gaussian_g_diff(const std::vector<double>& tau, double w0, double z0,
                         double diffusion, double** out_view,
                         int* n_out_view) {
  std::vector<double> g(tau.size());
  for (std::size_t i = 0; i < tau.size(); ++i) {
    const double lateral = 1.0 / (1.0 + 4.0 * diffusion * tau[i] / (w0 * w0));
    const double axial =
        1.0 / std::sqrt(1.0 + 4.0 * diffusion * tau[i] / (z0 * z0));
    g[i] = lateral * axial;
  }
  publish(g, out_view, n_out_view);
}

void fcs_bunching_factor(const std::vector<double>& tau, double k_exc_0,
                         const std::vector<double>& dark_matrix,
                         const std::vector<double>& exc_matrix, int n_states,
                         const std::vector<double>& brightness,
                         const std::vector<double>& brightness_b,
                         double** out_view, int* n_out_view) {
  std::vector<double> ones(tau.size(), 1.0);
  Eigen::MatrixXd k_d, k_e;
  generator_matrices(dark_matrix, exc_matrix, n_states, &k_d, &k_e);
  const Eigen::MatrixXd K =
      k_d + std::max(0.0, k_exc_0) * k_e;

  Eigen::VectorXd p_eq;
  if (!stationary(K, &p_eq)) {
    publish(ones, out_view, n_out_view);
    return;
  }
  Eigen::Map<const Eigen::VectorXd> q_a(brightness.data(), n_states);
  const bool cross = !brightness_b.empty();
  Eigen::VectorXd q_b =
      cross ? Eigen::Map<const Eigen::VectorXd>(brightness_b.data(), n_states)
                  .eval()
            : q_a.eval();
  if (static_cast<int>(brightness.size()) != n_states ||
      (cross && static_cast<int>(brightness_b.size()) != n_states)) {
    throw std::invalid_argument(
        "fcs_bunching_factor: brightness length must equal n_states");
  }
  const double avg_a = q_a.dot(p_eq);
  const double avg_b = q_b.dot(p_eq);
  if (!(avg_a > 0.0) || !(avg_b > 0.0)) {
    publish(ones, out_view, n_out_view);
    return;
  }
  const double norm = avg_a * avg_b;

  // The relaxation modes: eigenvalues clipped to Re <= 0 (numerical noise
  // on the stationary one must not blow exp up), amplitudes
  // c_m = (q_a^T V) o (V^-1 (q_b o p_eq)) -- complex for a cyclic scheme,
  // the sum real.
  Eigen::EigenSolver<Eigen::MatrixXd> solver(K);
  if (solver.info() != Eigen::Success) {
    publish(ones, out_view, n_out_view);
    return;
  }
  Eigen::VectorXcd evals = solver.eigenvalues();
  for (int m = 0; m < n_states; ++m) {
    evals(m) = std::complex<double>(std::min(evals(m).real(), 0.0),
                                    evals(m).imag());
  }
  const Eigen::MatrixXcd evecs = solver.eigenvectors();
  const Eigen::VectorXcd rhs =
      (q_b.array() * p_eq.array()).matrix().cast<std::complex<double>>();
  const Eigen::VectorXcd left = evecs.transpose() * q_a.cast<std::complex<double>>();
  const Eigen::VectorXcd right = evecs.partialPivLu().solve(rhs);
  const Eigen::VectorXcd c_m = left.array() * right.array();

  std::vector<double> x(tau.size());
  for (std::size_t t = 0; t < tau.size(); ++t) {
    std::complex<double> s(0.0, 0.0);
    for (int m = 0; m < n_states; ++m) {
      s += c_m(m) * std::exp(evals(m) * tau[t]);
    }
    x[t] = s.real() / norm;
  }
  publish(x, out_view, n_out_view);
}

void fcs_saturated_curve_shape(
    const std::vector<double>& tau, double power_w, double extinction,
    const std::vector<double>& dark_matrix,
    const std::vector<double>& exc_matrix, int n_states,
    const std::vector<double>& brightness, double w0, double z0,
    double diffusion, bool include_bunching, int n_r, int n_z,
    double wavelength_m, const std::vector<double>& brightness_b,
    double** out_view, int* n_out_view) {
  const std::size_t nt = tau.size();
  const double k_exc_0 =
      power_w > 0.0
          ? fcs_excitation_rate_peak(power_w, extinction, w0, wavelength_m)
          : 0.0;

  std::vector<double> g(nt, 0.0);
  if (k_exc_0 <= 0.0) {
    // The exact zero-power limit: the unsaturated Gaussian at amplitude 1.
    for (std::size_t i = 0; i < nt; ++i) {
      const double lateral = 1.0 / (1.0 + 4.0 * diffusion * tau[i] / (w0 * w0));
      const double axial =
          1.0 / std::sqrt(1.0 + 4.0 * diffusion * tau[i] / (z0 * z0));
      g[i] = lateral * axial;
    }
  } else {
    const double extent = 5.0;  // GRID_EXTENT_WAISTS
    const double r_max = extent * w0;
    const double dr = n_r > 1 ? r_max / (n_r - 1) : 1.0;
    const double dz = n_z > 1 ? 2.0 * extent * z0 / (n_z - 1) : 1.0;
    std::vector<double> r(static_cast<std::size_t>(n_r));
    for (int i = 0; i < n_r; ++i) r[static_cast<std::size_t>(i)] = dr * i;

    // The excitation grid and the per-point steady state. The generator at
    // each point is K_d + k_exc(r, z) K_e with the last row replaced by the
    // normalisation constraint; N is tiny, PartialPivLU per point is fine.
    Eigen::MatrixXd k_d, k_e;
    generator_matrices(dark_matrix, exc_matrix, n_states, &k_d, &k_e);
    Eigen::Map<const Eigen::VectorXd> q_a(brightness.data(), n_states);
    const bool cross = !brightness_b.empty();
    Eigen::VectorXd q_b =
        cross
            ? Eigen::Map<const Eigen::VectorXd>(brightness_b.data(), n_states)
                  .eval()
            : q_a.eval();

    Eigen::MatrixXd profile(n_r, n_z);
    Eigen::MatrixXd profile_b(n_r, n_z);
    Eigen::MatrixXd K(n_states, n_states);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(n_states);
    b(n_states - 1) = 1.0;
    for (int i = 0; i < n_r; ++i) {
      const double ri = r[static_cast<std::size_t>(i)];
      for (int j = 0; j < n_z; ++j) {
        const double zj = -extent * z0 + dz * j;
        const double psf = std::exp(-2.0 * ri * ri / (w0 * w0)) *
                           std::exp(-2.0 * zj * zj / (z0 * z0));
        const double k_exc = k_exc_0 * psf;
        K = k_d + k_exc * k_e;
        K.row(n_states - 1).setOnes();
        Eigen::VectorXd p = K.partialPivLu().solve(b);
        if (!p.allFinite()) p = K.colPivHouseholderQr().solve(b);
        profile(i, j) = q_a.dot(p);
        if (cross) profile_b(i, j) = q_b.dot(p);
      }
    }
    const Eigen::MatrixXd& prof_b = cross ? profile_b : profile;

    const double int_a = cyl_integral(profile, r, dr, dz);
    const double int_b =
        cross ? cyl_integral(prof_b, r, dr, dz) : int_a;
    const Eigen::MatrixXd prod = profile.cwiseProduct(prof_b);
    const double int_ab = cyl_integral(prod, r, dr, dz);

    if (int_a > 0.0 && int_b > 0.0 && int_ab != 0.0) {
      // Reciprocal space: Hankel along r (cached quadrature matrix), real
      // DFT along z. nz is tens, so the direct O(nz^2) transform is
      // microseconds and buys freedom from any FFT dependency.
      const double kr_max = 30.0 / r_max;
      const int n_kr = std::min(n_r, 64);
      const Eigen::MatrixXd& j0m = hankel_matrix(n_r, r_max, kr_max, n_kr);

      const int n_kz = n_z / 2 + 1;
      // f_zk[j][k] = sum_z profile(:, z) e^{-i 2pi k z / nz} * dz, done as
      // two real transforms; then the Hankel product.
      Eigen::MatrixXd fz_re(n_r, n_kz), fz_im(n_r, n_kz);
      Eigen::MatrixXd gz_re, gz_im;
      auto rdft = [&](const Eigen::MatrixXd& src, Eigen::MatrixXd* re,
                      Eigen::MatrixXd* im) {
        re->resize(n_r, n_kz);
        im->resize(n_r, n_kz);
        for (int k = 0; k < n_kz; ++k) {
          for (int i = 0; i < n_r; ++i) {
            double sr = 0.0, si = 0.0;
            for (int j = 0; j < n_z; ++j) {
              const double angle = -2.0 * kSatPi * k * j / n_z;
              sr += src(i, j) * std::cos(angle);
              si += src(i, j) * std::sin(angle);
            }
            (*re)(i, k) = sr * dz;
            (*im)(i, k) = si * dz;
          }
        }
      };
      rdft(profile, &fz_re, &fz_im);
      Eigen::MatrixXd f_re = j0m * fz_re;
      Eigen::MatrixXd f_im = j0m * fz_im;
      Eigen::MatrixXd csd;
      if (cross) {
        rdft(prof_b, &gz_re, &gz_im);
        const Eigen::MatrixXd g_re = j0m * gz_re;
        const Eigen::MatrixXd g_im = j0m * gz_im;
        csd = f_re.cwiseProduct(g_re) + f_im.cwiseProduct(g_im);
      } else {
        csd = f_re.cwiseProduct(f_re) + f_im.cwiseProduct(f_im);
      }

      // kr weight and the rfft mirror restoration (double every paired
      // bin; DC and, for even nz, Nyquist have no partner).
      std::vector<double> kz(static_cast<std::size_t>(n_kz));
      for (int k = 0; k < n_kz; ++k) {
        kz[static_cast<std::size_t>(k)] = 2.0 * kSatPi * k / (n_z * dz);
      }
      double total = 0.0;
      for (int a = 0; a < n_kr; ++a) {
        const double kr = kr_max * a / (n_kr > 1 ? n_kr - 1 : 1);
        for (int k = 0; k < n_kz; ++k) {
          double mirror = 2.0;
          if (k == 0 || (n_z % 2 == 0 && k == n_kz - 1)) mirror = 1.0;
          csd(a, k) *= kr * mirror;
          total += csd(a, k);
        }
      }

      if (total != 0.0) {
        const double v0 = std::pow(kSatPi, 1.5) * w0 * w0 * z0;
        const double amplitude = v0 * int_ab / (int_a * int_b);
        // The separable propagator, factored: exponentials are
        // (n_kr + n_kz) * n_tau, the rest two matrix products.
        for (std::size_t t = 0; t < nt; ++t) {
          double s = 0.0;
          for (int a = 0; a < n_kr; ++a) {
            const double kr = kr_max * a / (n_kr > 1 ? n_kr - 1 : 1);
            const double dec_r = std::exp(-diffusion * kr * kr * tau[t]);
            if (dec_r == 0.0) continue;
            double inner = 0.0;
            for (int k = 0; k < n_kz; ++k) {
              inner += csd(a, k) *
                       std::exp(-diffusion *
                                kz[static_cast<std::size_t>(k)] *
                                kz[static_cast<std::size_t>(k)] * tau[t]);
            }
            s += dec_r * inner;
          }
          g[t] = s * amplitude / total;
        }
      }
    }
  }

  if (include_bunching) {
    double* x_view = nullptr;
    int n_x = 0;
    fcs_bunching_factor(tau, k_exc_0, dark_matrix, exc_matrix, n_states,
                        brightness, brightness_b, &x_view, &n_x);
    for (std::size_t i = 0; i < nt; ++i) g[i] *= x_view[i];
    std::free(x_view);
  }
  publish(g, out_view, n_out_view);
}

IMPBFF_END_NAMESPACE
