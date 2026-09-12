/**\file PhotophysicsLifetimeSpectrumNode.cpp
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/PhotophysicsLifetimeSpectrumNode.h>
#include "internal/SpectrumNodeHelpers.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

IMPBFF_BEGIN_NAMESPACE

// ------------------------------------------------------ PhotophysicsLifetimeSpectrumNode

// Empty for the reason `TCSPCDecay`'s is: `add_port` reaches for
// `shared_from_this()`, so a node cannot own ports until something owns it.
PhotophysicsLifetimeSpectrumNode::PhotophysicsLifetimeSpectrumNode(const std::string& name)
    : GraphNode(name) {}

void PhotophysicsLifetimeSpectrumNode::set_number_of_lifetimes(int n) {
  if (n <= 0) {
    throw std::domain_error(
        "PhotophysicsLifetimeSpectrumNode::set_number_of_lifetimes: a spectrum has at "
        "least one species");
  }
  lifetime_ports_.clear();
  lifetime_ports_.reserve(static_cast<std::size_t>(2 * n));
  for (int i = 0; i < n; ++i) {
    std::ostringstream a, t;
    a << "a" << i;
    t << "t" << i;
    GraphPort* amplitude = nullptr;
    GraphPort* lifetime = nullptr;
    // Ports cannot be removed from a node, so a second call with a smaller
    // count leaves the surplus ports visible but unread -- same caveat, and
    // same remedy (build the node once), as `TCSPCDecay`.
    if (get_input_port(a.str())) {
      amplitude = get_input_port(a.str()).get();
      lifetime = get_input_port(t.str()).get();
    } else {
      spectrum_node_detail::add_scalar_port(this, a.str(), 1.0, &amplitude);
      spectrum_node_detail::add_scalar_port(this, t.str(), 1.0, &lifetime);
    }
    lifetime_ports_.push_back(amplitude);
    lifetime_ports_.push_back(lifetime);
  }
  n_lifetimes_ = n;
  spectrum_.assign(static_cast<std::size_t>(2 * n), 0.0);
  set_valid(false);
}

void PhotophysicsLifetimeSpectrumNode::evaluate() {
  if (n_lifetimes_ <= 0) {
    throw std::domain_error("PhotophysicsLifetimeSpectrumNode '" + get_name() +
                            "' has no species; call "
                            "set_number_of_lifetimes() first");
  }
  spectrum_.resize(static_cast<std::size_t>(2 * n_lifetimes_));
  double sum = 0.0;
  for (int i = 0; i < n_lifetimes_; ++i) {
    const std::size_t k = static_cast<std::size_t>(2 * i);
    double amplitude = lifetime_ports_[k]->get_value();
    if (absolute_amplitudes_) amplitude = std::fabs(amplitude);
    spectrum_[k] = amplitude;
    sum += amplitude;
    // Absolute, for the reason `TCSPCDecay` takes it: a lifetime walked
    // through zero is a growing exponential, not a decay.
    spectrum_[k + 1] = std::fabs(lifetime_ports_[k + 1]->get_value());
  }
  if (normalize_amplitudes_) {
    const double scale = std::fabs(sum);
    for (int i = 0; i < n_lifetimes_; ++i) {
      spectrum_[static_cast<std::size_t>(2 * i)] /= scale;
    }
  }
  spectrum_node_detail::publish(this, spectrum_);
  set_valid(true);
}

IMPBFF_END_NAMESPACE
