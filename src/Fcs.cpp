/**
 * \file Fcs.cpp
 * \brief FCS forward models beyond the closed forms: the MDF and saturation shapes.
 *
 * Sections in the order of IMP/bff/Fcs.h; each is marked with the file it
 * came from.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

// -------- from FcsMdf.cpp --------
/**
 *  (formerly FcsMdf.cpp, now a section of this file)
 *  \brief Enderlein molecule-detection-function FCS forward model.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/Fcs.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

IMPBFF_BEGIN_NAMESPACE

namespace {

const double kPi = 3.14159265358979323846;

//! w(z): axial dependence of the lateral 1/e^2 beam radius.
inline double beam_w(double z, double w0, double lambda_ex, double n) {
  const double q = lambda_ex * z / (kPi * w0 * w0 * n);
  return w0 * std::sqrt(1.0 + q * q);
}

//! kappa(z): axial collection efficiency of the imaged pinhole.
inline double kappa(double z, double r0, double lambda_em, double n,
                    double a) {
  const double q = lambda_em * z / (kPi * r0 * r0 * n);
  const double R2 = r0 * r0 * (1.0 + q * q);
  return 1.0 - std::exp(-2.0 * a * a / R2);
}

//! Uniform-grid trapezoid: h * (sum - (f0 + fn)/2).
inline double trapz_uniform(const std::vector<double>& f, double h) {
  double s = 0.0;
  for (double v : f) s += v;
  s -= 0.5 * (f.front() + f.back());
  return s * h;
}

//! The symmetric axial grid spanning `span` Rayleigh ranges.
void axial_grid(double w0, double r0, double lambda_ex, double lambda_em,
                double n, int n_grid, double span, std::vector<double>& z,
                double* h) {
  const double z_r = kPi * std::max(w0 * w0, r0 * r0) * n /
                     std::min(lambda_ex, lambda_em);
  const double z_max = span * z_r;
  z.resize(static_cast<std::size_t>(n_grid));
  const double step =
      n_grid > 1 ? 2.0 * z_max / static_cast<double>(n_grid - 1) : 0.0;
  for (int i = 0; i < n_grid; ++i) {
    z[static_cast<std::size_t>(i)] = -z_max + step * i;
  }
  *h = step;
}

}  // namespace

void hermgauss(int n, std::vector<double>& nodes,
               std::vector<double>& weights) {
  if (n < 1) throw std::invalid_argument("hermgauss: n must be positive");
  nodes.assign(static_cast<std::size_t>(n), 0.0);
  weights.assign(static_cast<std::size_t>(n), 0.0);
  // Roots come in +/- pairs; find the positive half by Newton over the
  // *orthonormal* recurrence h_{k+1} = x sqrt(2/(k+1)) h_k
  // - sqrt(k/(k+1)) h_{k-1}, which stays O(1) where the plain physicists'
  // recurrence overflows past n ~ 20. Weight: w_i = 2 / h'_n(x_i)^2.
  const int m = (n + 1) / 2;
  const double pim4 = 0.7511255444649425;  // pi^(-1/4), h_0's value
  double x = 0.0;
  for (int i = 0; i < m; ++i) {
    // Standard initial guesses (largest root first, then stepping down).
    if (i == 0) {
      x = std::sqrt(2.0 * n + 1.0) -
          1.85575 * std::pow(2.0 * n + 1.0, -1.0 / 6.0);
    } else if (i == 1) {
      x -= 1.14 * std::pow(static_cast<double>(n), 0.426) / x;
    } else if (i == 2) {
      x = 1.86 * x - 0.86 * nodes[static_cast<std::size_t>(n - 1)];
    } else if (i == 3) {
      x = 1.91 * x - 0.91 * nodes[static_cast<std::size_t>(n - 2)];
    } else {
      x = 2.0 * x - nodes[static_cast<std::size_t>(n - i + 1)];
    }
    double pp = 0.0;
    for (int iteration = 0; iteration < 100; ++iteration) {
      double p1 = pim4;
      double p2 = 0.0;
      for (int j = 0; j < n; ++j) {
        const double p3 = p2;
        p2 = p1;
        p1 = x * std::sqrt(2.0 / (j + 1.0)) * p2 -
             std::sqrt(static_cast<double>(j) / (j + 1.0)) * p3;
      }
      pp = std::sqrt(2.0 * n) * p2;  // h'_n(x) via the recurrence identity
      const double dx = p1 / pp;
      x -= dx;
      if (std::fabs(dx) <= 1e-15 * std::max(1.0, std::fabs(x))) break;
    }
    // Ascending order, mirrored: the i-th largest positive root sits at the
    // top, its negative twin at the bottom (numpy's ordering).
    nodes[static_cast<std::size_t>(n - 1 - i)] = x;
    nodes[static_cast<std::size_t>(i)] = -x;
    const double w = 2.0 / (pp * pp);
    weights[static_cast<std::size_t>(n - 1 - i)] = w;
    weights[static_cast<std::size_t>(i)] = w;
  }
  if (n % 2 == 1) {
    // The middle root of an odd order is exactly zero; the Newton loop
    // above found it from a nonzero guess, but pin it regardless.
    nodes[static_cast<std::size_t>(n / 2)] = 0.0;
  }
}

double fcs_mdf_effective_volume(double w0, double r0,
                                double excitation_wavelength,
                                double emission_wavelength,
                                double refractive_index,
                                double pinhole_radius, int n_grid,
                                double span) {
  std::vector<double> z;
  double h = 0.0;
  axial_grid(w0, r0, excitation_wavelength, emission_wavelength,
             refractive_index, n_grid, span, z, &h);
  std::vector<double> k(z.size()), k2w2(z.size());
  for (std::size_t i = 0; i < z.size(); ++i) {
    const double ki =
        kappa(z[i], r0, emission_wavelength, refractive_index, pinhole_radius);
    const double wi = beam_w(z[i], w0, excitation_wavelength, refractive_index);
    k[i] = ki;
    k2w2[i] = ki * ki / (wi * wi);
  }
  const double num = trapz_uniform(k, h);
  const double den = trapz_uniform(k2w2, h);
  return kPi * num * num / den;
}

void fcs_mdf_axial_profiles(double w0, double r0, double excitation_wavelength,
                            double emission_wavelength, double refractive_index,
                            double pinhole_radius, int n_grid, double span,
                            std::vector<double>& z, double& h,
                            std::vector<double>& kappa_z,
                            std::vector<double>& w2_z) {
  axial_grid(w0, r0, excitation_wavelength, emission_wavelength,
             refractive_index, n_grid, span, z, &h);
  kappa_z.resize(z.size());
  w2_z.resize(z.size());
  for (std::size_t i = 0; i < z.size(); ++i) {
    kappa_z[i] = kappa(z[i], r0, emission_wavelength, refractive_index,
                       pinhole_radius);
    const double wi = beam_w(z[i], w0, excitation_wavelength, refractive_index);
    w2_z[i] = wi * wi;
  }
}

double fcs_mdf_g_raw(double t, double sep2, double diffusion,
                     const std::vector<double>& z, double h,
                     const std::vector<double>& k,
                     const std::vector<double>& w2,
                     const std::vector<double>& xi,
                     const std::vector<double>& hq, double w0, double r0,
                     double lambda_ex, double lambda_em, double refr,
                     double a) {
  t = std::max(t, 1e-18);
  const double s = 4.0 * diffusion * t;
  const double sqrt_s = std::sqrt(s);
  std::vector<double> integrand(z.size());
  for (std::size_t i = 0; i < z.size(); ++i) {
    double inner = 0.0;
    for (std::size_t j = 0; j < xi.size(); ++j) {
      const double zp = z[i] + sqrt_s * xi[j];
      const double kzp = kappa(zp, r0, lambda_em, refr, a);
      const double wzp = beam_w(zp, w0, lambda_ex, refr);
      const double w_sum = w2[i] + wzp * wzp;
      double g_lat = 1.0 / (4.0 + w_sum / (2.0 * diffusion * t));
      if (sep2 > 0.0) {
        // Two-focus cross-correlation: the lateral overlap of foci a
        // distance d apart decays with the lag, which is what yields an
        // absolute D from the known d.
        g_lat *= std::exp(-sep2 / (4.0 * diffusion * t + 0.5 * w_sum));
      }
      inner += hq[j] * kzp * g_lat;
    }
    integrand[i] = k[i] * inner;
  }
  return trapz_uniform(integrand, h) / s;
}

void fcs_mdf_g_diff(const std::vector<double>& tau, double w0, double r0,
                    double diffusion, double excitation_wavelength,
                    double emission_wavelength, double refractive_index,
                    double pinhole_radius, int n_grid, double span,
                    bool normalize, int n_herm, double separation,
                    double** out_view, int* n_out_view) {
  std::vector<double> z, k, w2;
  double h = 0.0;
  fcs_mdf_axial_profiles(w0, r0, excitation_wavelength, emission_wavelength,
                         refractive_index, pinhole_radius, n_grid, span, z, h,
                         k, w2);
  std::vector<double> xi, hq;
  hermgauss(n_herm, xi, hq);

  const double d2 = separation * separation;
  // The small-lag plateau, always without the separation term: an
  // autocorrelation's g(0), and the cross-correlation's normaliser.
  const double num0 =
      fcs_mdf_g_raw(1e-15, 0.0, diffusion, z, h, k, w2, xi, hq, w0, r0,
                    excitation_wavelength, emission_wavelength,
                    refractive_index, pinhole_radius);
  const double veff =
      normalize ? 1.0
                : fcs_mdf_effective_volume(w0, r0, excitation_wavelength,
                                           emission_wavelength,
                                           refractive_index, pinhole_radius);

  double* out =
      static_cast<double*>(std::malloc(sizeof(double) * tau.size()));
  if (out == nullptr) throw std::bad_alloc();
  for (std::size_t i = 0; i < tau.size(); ++i) {
    const double t = tau[i];
    const double raw =
        t <= 0.0
            ? num0
            : fcs_mdf_g_raw(t, d2, diffusion, z, h, k, w2, xi, hq, w0, r0,
                            excitation_wavelength, emission_wavelength,
                            refractive_index, pinhole_radius);
    out[i] = raw / num0 / veff;
  }
  *out_view = out;
  *n_out_view = static_cast<int>(tau.size());
}

// --------------------------------------------------------- FcsMdfCurve

// Empty for the reason `FretSpectrum`'s is: `add_port` reaches for
// `shared_from_this()`, so a node cannot own ports until something owns it.
FcsMdfCurve::FcsMdfCurve(const std::string& name) : Node(name) {}

namespace {

//! A scalar input port, remembered by raw pointer (the node owns it).
/*! Named for this file: the library builds as one translation unit, so an
    anonymous namespace does not keep two helpers of the same name apart. */
void add_mdf_port(Node* node, const std::string& key, double value,
                  Port** slot) {
  std::shared_ptr<Port> port(new Port(value));
  node->add_input_port(key, port);
  *slot = port.get();
}

}  // namespace

void FcsMdfCurve::build_ports() {
  if (w0_port_ != nullptr) return;  // built once; ports cannot be removed
  add_mdf_port(this, "w0", 0.25, &w0_port_);
  add_mdf_port(this, "wem", 0.25, &wem_port_);
  add_mdf_port(this, "D", 300.0, &d_port_);
  add_mdf_port(this, "diam", 0.0, &diam_port_);
  set_valid(false);
}

void FcsMdfCurve::set_axis(const std::vector<double>& tau) {
  tau_ = tau;
  curve_built_ = false;
  set_valid(false);
}

void FcsMdfCurve::set_axis_array(double* in_axis, int n_axis) {
  tau_.assign(in_axis, in_axis + n_axis);
  curve_built_ = false;
  set_valid(false);
}

void FcsMdfCurve::set_optics(double excitation_wavelength,
                             double emission_wavelength,
                             double refractive_index, double pinhole_radius) {
  excitation_wavelength_ = excitation_wavelength;
  emission_wavelength_ = emission_wavelength;
  refractive_index_ = refractive_index;
  pinhole_radius_ = pinhole_radius;
  curve_built_ = false;
  set_valid(false);
}

void FcsMdfCurve::set_quadrature(int n_grid, double span, int n_herm) {
  if (n_grid < 2 || n_herm < 1) {
    throw std::domain_error(
        "FcsMdfCurve::set_quadrature: the axial grid needs at least two "
        "points and the Gauss-Hermite rule at least one node");
  }
  n_grid_ = n_grid;
  span_ = span;
  n_herm_ = n_herm;
  curve_built_ = false;
  set_valid(false);
}

void FcsMdfCurve::set_length_scale(double micrometres_per_unit) {
  length_scale_ = micrometres_per_unit;
  curve_built_ = false;
  set_valid(false);
}

void FcsMdfCurve::set_normalize(bool v) {
  normalize_ = v;
  curve_built_ = false;
  set_valid(false);
}

void FcsMdfCurve::refresh_quadrature() {
  if (quadrature_n_ == n_herm_) return;
  hermgauss(n_herm_, herm_nodes_, herm_weights_);
  quadrature_n_ = n_herm_;
  ++quadrature_builds_;
}

void FcsMdfCurve::refresh_profiles(double w0, double r0) {
  if (profiles_built_ && profile_w0_ == w0 && profile_r0_ == r0 &&
      profile_lambda_ex_ == excitation_wavelength_ &&
      profile_lambda_em_ == emission_wavelength_ &&
      profile_n_ == refractive_index_ && profile_a_ == pinhole_radius_ &&
      profile_n_grid_ == n_grid_ && profile_span_ == span_) {
    return;
  }
  fcs_mdf_axial_profiles(w0, r0, excitation_wavelength_, emission_wavelength_,
                         refractive_index_, pinhole_radius_, n_grid_, span_,
                         z_, h_, kappa_z_, w2_z_);
  profiles_built_ = true;
  profile_w0_ = w0;
  profile_r0_ = r0;
  profile_lambda_ex_ = excitation_wavelength_;
  profile_lambda_em_ = emission_wavelength_;
  profile_n_ = refractive_index_;
  profile_a_ = pinhole_radius_;
  profile_n_grid_ = n_grid_;
  profile_span_ = span_;
  ++profile_builds_;
}

void FcsMdfCurve::evaluate() {
  if (w0_port_ == nullptr) {
    throw std::domain_error(
        "FcsMdfCurve::evaluate: build_ports() has not been called, so the "
        "node has no waists to integrate over");
  }
  const std::shared_ptr<Port> out = get_output_port(get_name());
  if (!out) {
    throw std::domain_error(
        "'" + get_name() +
        "' writes its shape to the output port keyed by its own name, which "
        "this node does not have");
  }
  const double w0 = w0_port_->get_value() * length_scale_;
  const double r0 = wem_port_->get_value() * length_scale_;
  const double diffusion = d_port_->get_value();
  const double separation = diam_port_->get_value() * length_scale_;

  // Bit-identical ports mean a bit-identical shape, so re-publishing the one
  // already computed *is* the model. This is what makes the graph route
  // cheaper than the director rather than merely differently arranged: a
  // Levenberg-Marquardt Jacobian differences one parameter at a time, and
  // the columns belonging to N, to the baseline or to a bunching term move
  // nothing this node reads.
  if (!(curve_built_ && curve_w0_ == w0 && curve_r0_ == r0 &&
        curve_d_ == diffusion && curve_sep_ == separation)) {
    refresh_quadrature();
    refresh_profiles(w0, r0);
    const double d2 = separation * separation;
    // The small-lag plateau, always without the separation term: an
    // autocorrelation's g(0), and the cross-correlation's normaliser.
    const double num0 = fcs_mdf_g_raw(
        1e-15, 0.0, diffusion, z_, h_, kappa_z_, w2_z_, herm_nodes_,
        herm_weights_, w0, r0, excitation_wavelength_, emission_wavelength_,
        refractive_index_, pinhole_radius_);
    const double veff =
        normalize_ ? 1.0
                   : fcs_mdf_effective_volume(
                         w0, r0, excitation_wavelength_, emission_wavelength_,
                         refractive_index_, pinhole_radius_);
    curve_.resize(tau_.size());
    for (std::size_t i = 0; i < tau_.size(); ++i) {
      const double t = tau_[i];
      const double raw =
          t <= 0.0 ? num0
                   : fcs_mdf_g_raw(t, d2, diffusion, z_, h_, kappa_z_, w2_z_,
                                   herm_nodes_, herm_weights_, w0, r0,
                                   excitation_wavelength_,
                                   emission_wavelength_, refractive_index_,
                                   pinhole_radius_);
      curve_[i] = raw / num0 / veff;
    }
    curve_built_ = true;
    curve_w0_ = w0;
    curve_r0_ = r0;
    curve_d_ = diffusion;
    curve_sep_ = separation;
    ++kernel_evaluations_;
  }

  // Fit transport, so it does not sanitise: a non-finite shape (a waist the
  // optimiser walked to zero) has to reach `ChiSquared` and make the misfit
  // infinite rather than be floored, which in a fit reads as a *good* fit.
  out->set_sanitize(false);
  out->set_value_vector(curve_);
  set_valid(true);
}

std::string FcsMdfCurve::describe() const {
  std::ostringstream out;
  out << "FcsMdfCurve(name='" << get_name() << "', lags=" << tau_.size()
      << ", n_grid=" << n_grid_ << ", span=" << span_
      << ", n_herm=" << n_herm_ << ", normalize=" << (normalize_ ? 1 : 0)
      << ", quadrature_builds=" << quadrature_builds_
      << ", profile_builds=" << profile_builds_
      << ", kernel_evaluations=" << kernel_evaluations_ << ")";
  return out.str();
}

IMPBFF_END_NAMESPACE

// -------- from FcsSaturation.cpp --------
/**
 *  (formerly FcsSaturation.cpp, now a section of this file)
 *  \brief Saturated FCS forward model: steady-state photophysics on a grid.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <complex>
#include <map>
#include <mutex>

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

// -------- from FcsSaturationCurve.cpp --------
/**
 *  (formerly FcsSaturationCurve.cpp, now a section of this file)
 *  \brief Saturated FCS forward model as a graph node.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */



IMPBFF_BEGIN_NAMESPACE

FcsSaturationCurve::FcsSaturationCurve(const std::string& name)
    : Node(name) {}

namespace {

//! A scalar input port, remembered by raw pointer (the node owns it).
void add_sat_port(Node* node, const std::string& key, double value,
                  Port** slot) {
  std::shared_ptr<Port> port(new Port(value));
  node->add_input_port(key, port);
  *slot = port.get();
}

}  // namespace

void FcsSaturationCurve::build_ports() {
  if (power_port_ != nullptr) return;  // built once; ports cannot be removed
  // Defaults are the model's own (mW, nm, um^2/s) -- they are overwritten
  // at build time, but a sane value keeps a bare node from crashing.
  add_sat_port(this, "power", 0.0, &power_port_);
  add_sat_port(this, "extinction", 80000.0, &extinction_port_);
  add_sat_port(this, "w0", 200.0, &w0_port_);
  add_sat_port(this, "z0", 1000.0, &z0_port_);
  add_sat_port(this, "D", 400.0, &d_port_);
  add_sat_port(this, "N", 1.0, &n_port_);
  add_sat_port(this, "b", 1.0, &b_port_);
  add_sat_port(this, "bg", 0.0, &bg_port_);
  set_valid(false);
}

void FcsSaturationCurve::set_axis(const std::vector<double>& tau) {
  tau_ = tau;
  curve_built_ = false;
  set_valid(false);
}

void FcsSaturationCurve::set_axis_array(double* in_axis, int n_axis) {
  tau_.assign(in_axis, in_axis + n_axis);
  curve_built_ = false;
  set_valid(false);
}

void FcsSaturationCurve::set_scheme(const std::vector<double>& dark_matrix,
                                     const std::vector<double>& exc_matrix,
                                     int n_states,
                                     const std::vector<double>& brightness) {
  dark_matrix_ = dark_matrix;
  exc_matrix_ = exc_matrix;
  n_states_ = n_states;
  brightness_ = brightness;
  curve_built_ = false;
  set_valid(false);
}

void FcsSaturationCurve::set_quadrature(int n_r, int n_z) {
  n_r_ = n_r;
  n_z_ = n_z;
  curve_built_ = false;
  set_valid(false);
}

void FcsSaturationCurve::set_wavelength(double wavelength_m) {
  wavelength_m_ = wavelength_m;
  curve_built_ = false;
  set_valid(false);
}

void FcsSaturationCurve::set_include_bunching(bool v) {
  include_bunching_ = v;
  curve_built_ = false;
  set_valid(false);
}

void FcsSaturationCurve::evaluate() {
  if (power_port_ == nullptr) {
    throw std::domain_error(
        "FcsSaturationCurve::evaluate: build_ports() has not been called, "
        "so the node has no ports to read");
  }
  const std::shared_ptr<Port> out = get_output_port(get_name());
  if (!out) {
    throw std::domain_error(
        "'" + get_name() +
        "' writes its shape to the output port keyed by its own name, which "
        "this node does not have");
  }
  if (n_states_ < 2) {
    throw std::domain_error(
        "FcsSaturationCurve::evaluate: the photokinetic scheme has not been "
        "set (call set_scheme first)");
  }

  // Read ports in the model's units and scale to SI, exactly as the Python
  // wrapper `saturated_curve_shape` does.
  const double power_W = power_port_->get_value() * 1e-3;    // mW -> W
  const double extinction = extinction_port_->get_value();  // M^-1 cm^-1
  const double w0_m = w0_port_->get_value() * 1e-9;          // nm -> m
  const double z0_m = z0_port_->get_value() * 1e-9;          // nm -> m
  const double D_m2s = d_port_->get_value() * 1e-12;         // um^2/s -> m^2/s

  // Bit-identical ports mean a bit-identical shape.  A Levenberg-Marquardt
  // Jacobian differences one parameter at a time, and the columns belonging
  // to N, b or bg move nothing this node reads.
  if (!(curve_built_ && curve_power_ == power_W &&
        curve_extinction_ == extinction && curve_w0_ == w0_m &&
        curve_z0_ == z0_m && curve_d_ == D_m2s)) {
    // fcs_saturated_curve_shape publishes the result through an owned
    // malloc'd buffer (see publish() in FcsSaturation.cpp); we adopt it
    // and free our copy.
    double* g_ptr = nullptr;
    int n_g = 0;
    fcs_saturated_curve_shape(
        tau_, power_W, extinction, dark_matrix_, exc_matrix_, n_states_,
        brightness_, w0_m, z0_m, D_m2s, include_bunching_, n_r_, n_z_,
        wavelength_m_, /*brightness_b=*/std::vector<double>(),
        &g_ptr, &n_g);
    curve_.assign(g_ptr, g_ptr + static_cast<std::size_t>(n_g));
    std::free(g_ptr);

    curve_built_ = true;
    curve_power_ = power_W;
    curve_extinction_ = extinction;
    curve_w0_ = w0_m;
    curve_z0_ = z0_m;
    curve_d_ = D_m2s;
  }

  // Fit transport: a non-finite shape must reach ChiSquared and make the
  // misfit infinite, not be silently floored.
  out->set_sanitize(false);
  out->set_value_vector(curve_);
  set_valid(true);
}

std::string FcsSaturationCurve::describe() const {
  std::ostringstream out;
  out << "FcsSaturationCurve(name='" << get_name()
      << "', lags=" << tau_.size()
      << ", n_states=" << n_states_
      << ", n_r=" << n_r_ << ", n_z=" << n_z_
      << ", wavelength=" << wavelength_m_
      << ", include_bunching=" << (include_bunching_ ? 1 : 0) << ")";
  return out.str();
}

IMPBFF_END_NAMESPACE
