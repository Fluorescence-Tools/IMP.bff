/**\file FRETSpectrumNode.cpp
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/FRETSpectrumNode.h>
#include "internal/SpectrumNodeHelpers.h"
#include <IMP/bff/PhotophysicsLifetimeSpectrum.h>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

IMPBFF_BEGIN_NAMESPACE

// ------------------------------------------------------------ FRETSpectrumNode

FRETSpectrumNode::FRETSpectrumNode(const std::string& name) : GraphNode(name) {}

void FRETSpectrumNode::build_ports() {
  if (donor_port_ != nullptr) return;
  std::shared_ptr<GraphPort> donor(new GraphPort(std::vector<double>{1.0, 4.0}));
  donor->set_sanitize(false);
  add_input_port(donor_port_key(), donor);
  donor_port_ = donor.get();

  std::shared_ptr<GraphPort> distances(new GraphPort(std::vector<double>{1.0, 52.0}));
  distances->set_sanitize(false);
  add_input_port(distance_port_key(), distances);
  distance_port_ = distances.get();

  spectrum_node_detail::add_scalar_port(this, "x_donly", 0.0, &x_donly_port_);
  spectrum_node_detail::add_scalar_port(this, "forster_radius", 52.0, &forster_radius_port_);
  spectrum_node_detail::add_scalar_port(this, "tau0", 4.0, &tau0_port_);
  spectrum_node_detail::add_scalar_port(this, "kappa2", 2.0 / 3.0, &kappa2_port_);
  set_valid(false);
}

void FRETSpectrumNode::evaluate() {
  if (donor_port_ == nullptr) {
    throw std::domain_error("FRETSpectrumNode '" + get_name() +
                            "' has no ports; call build_ports() first");
  }
  // By key, not through the cached pointer, for the reason
  // `PhotophysicsAnisotropySpectrumNode` gives: these two are the ports a caller replaces.
  const std::shared_ptr<GraphPort> donor_in = get_input_port(donor_port_key());
  const std::shared_ptr<GraphPort> distance_in = get_input_port(distance_port_key());
  if (!donor_in || !distance_in) {
    throw std::domain_error("FRETSpectrumNode '" + get_name() +
                            "' is missing one of its two spectrum ports");
  }
  const std::vector<double>& donor = donor_in->get_values_ref();
  const std::vector<double>& distances = distance_in->get_values_ref();
  spectrum_node_detail::check_interleaved("FRETSpectrumNode '" + get_name() + "' (donor)", donor);
  spectrum_node_detail::check_interleaved("FRETSpectrumNode '" + get_name() + "' (distances)",
                    distances);

  const double x_donly = std::fabs(x_donly_port_->get_value());
  const double forster_radius = forster_radius_port_->get_value();
  const double tau0 = tau0_port_->get_value();
  const double kappa2 = kappa2_port_->get_value();
  if (!(tau0 > 0.0)) {
    throw std::domain_error("FRETSpectrumNode '" + get_name() +
                            "': tau0 is not positive");
  }

  // The donor as *rates*, which is the space the two spectra combine in.
  const std::size_t n_donor = donor.size() / 2;
  std::vector<double> donor_rates(2 * n_donor);
  for (std::size_t i = 0; i < n_donor; ++i) {
    donor_rates[2 * i] = donor[2 * i];
    donor_rates[2 * i + 1] = 1.0 / donor[2 * i + 1];
  }

  // Each distance as a transfer rate: k = 3/2 kappa^2 / tau0 * (R0/r)^6.
  const std::size_t n_distances = distances.size() / 2;
  std::vector<double> fret_rates(2 * n_distances);
  const double prefactor = 1.5 * kappa2 / tau0;
  for (std::size_t i = 0; i < n_distances; ++i) {
    const double ratio = forster_radius / distances[2 * i + 1];
    const double ratio3 = ratio * ratio * ratio;
    fret_rates[2 * i] = distances[2 * i];
    fret_rates[2 * i + 1] = prefactor * ratio3 * ratio3;
  }

  // Rates add: a quenched species decays at its own rate plus the transfer
  // rate, so the quenched spectrum is the Cartesian product with the
  // amplitudes multiplied. FRET-major, which is the order ChiSurf's
  // `ere2(fret, donor)` produces.
  std::vector<double> combined;
  combined.reserve(2 * (n_distances * n_donor + n_donor));
  for (std::size_t i = 0; i < n_distances; ++i) {
    for (std::size_t j = 0; j < n_donor; ++j) {
      combined.push_back(fret_rates[2 * i] * donor_rates[2 * j] *
                         (1.0 - x_donly));
      combined.push_back(fret_rates[2 * i + 1] + donor_rates[2 * j + 1]);
    }
  }
  // The donor-only fraction, appended unconditionally -- including at
  // `x_donly == 0`, where it carries zero amplitude. A node that sometimes
  // returns a shorter spectrum is a second code path for no gain, and the
  // length is observable.
  for (std::size_t j = 0; j < n_donor; ++j) {
    combined.push_back(donor_rates[2 * j] * x_donly);
    combined.push_back(donor_rates[2 * j + 1]);
  }

  // Back to lifetimes, which is what the instrument node reconvolves.
  spectrum_.resize(combined.size());
  for (std::size_t i = 0; i < combined.size() / 2; ++i) {
    spectrum_[2 * i] = combined[2 * i];
    spectrum_[2 * i + 1] = 1.0 / combined[2 * i + 1];
  }
  spectrum_node_detail::publish(this, spectrum_);
  set_valid(true);
}

std::string FRETSpectrumNode::describe() const {
  std::ostringstream out;
  out << "species out    : " << spectrum_.size() / 2 << "\n"
      << "R0 / tau0      : " << forster_radius_port_->get_value() << " / "
      << tau0_port_->get_value() << "\n";
  return out.str();
}

IMPBFF_END_NAMESPACE
