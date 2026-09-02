/**
 *  \file FcsMdf.cpp
 *  \brief Enderlein molecule-detection-function FCS forward model.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/FcsMdf.h>

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
