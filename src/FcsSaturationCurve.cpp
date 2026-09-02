/**
 *  \file FcsSaturationCurve.cpp
 *  \brief Saturated FCS forward model as a graph node.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/FcsSaturationCurve.h>

#include <IMP/bff/FcsSaturation.h>

#include <cmath>
#include <sstream>
#include <stdexcept>

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
