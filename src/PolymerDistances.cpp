/**\file PolymerDistances.cpp
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/PolymerDistances.h>
#include "internal/SpectrumNodeHelpers.h"
#include <IMP/bff/PolymerChain.h>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

IMPBFF_BEGIN_NAMESPACE

// ------------------------------------------------------- PolymerDistances

PolymerDistances::PolymerDistances(const std::string& name) : GraphNode(name) {}

void PolymerDistances::set_mode(const std::string& mode) {
  if (!mode_.empty()) {
    throw std::domain_error("PolymerDistances '" + get_name() +
                            "': the mode is set once");
  }
  parameter_ports_.clear();
  GraphPort* p = nullptr;
  if (mode == "worm_like_chain" || mode == "worm_like_chain_linker") {
    spectrum_node_detail::add_scalar_port(this, "chain_length", 100.0, &p);
    parameter_ports_.push_back(p);
    spectrum_node_detail::add_scalar_port(this, "persistence_length", 30.0, &p);
    parameter_ports_.push_back(p);
    // Both modes carry the linker width. Without the linker it is inert --
    // evaluate() never reads it -- which mirrors chisurf exactly: the
    // model's `w` is a fitting parameter whether or not the linker is on,
    // and with it off the numpy path fits an inert parameter too. A port
    // the optimiser can claim is what keeps a free-but-inert `w` from
    // refusing the whole graph.
    spectrum_node_detail::add_scalar_port(this, "sigma_linker", 6.0, &p);
    parameter_ports_.push_back(p);
  } else if (mode == "saw_nu") {
    spectrum_node_detail::add_scalar_port(this, "r_rms", 50.0, &p);
    parameter_ports_.push_back(p);
    spectrum_node_detail::add_scalar_port(this, "nu", 0.588, &p);
    parameter_ports_.push_back(p);
  } else if (mode == "ising_chain") {
    const char* keys[] = {"n_residues", "b_structured", "b_unstructured",
                          "coupling", "field"};
    const double defaults[] = {50.0, 4.0, 8.0, 1.5, 0.0};
    for (int i = 0; i < 5; ++i) {
      spectrum_node_detail::add_scalar_port(this, keys[i], defaults[i], &p);
      parameter_ports_.push_back(p);
    }
  } else {
    throw std::domain_error("PolymerDistances '" + get_name() +
                            "': unknown mode '" + mode + "'");
  }
  mode_ = mode;
  set_valid(false);
}

void PolymerDistances::set_axis(const std::vector<double>& axis) {
  if (axis.size() < 2) {
    throw std::domain_error("PolymerDistances '" + get_name() +
                            "': a distance axis needs at least two points");
  }
  axis_ = axis;
  set_valid(false);
}

void PolymerDistances::set_axis_array(double* in_axis, int n_axis) {
  set_axis(std::vector<double>(in_axis, in_axis + n_axis));
}

void PolymerDistances::evaluate() {
  if (mode_.empty()) {
    throw std::domain_error("PolymerDistances '" + get_name() +
                            "' has no mode; call set_mode() first");
  }
  if (axis_.empty()) {
    throw std::domain_error("PolymerDistances '" + get_name() +
                            "' has no distance axis");
  }
  // The kernels publish malloc'd views (their numpy contract); copy and
  // free. They are the same functions chisurf's rdf.py forwarders call, so
  // the graph path and the Python path evaluate one implementation.
  double* view = nullptr;
  int n_view = 0;
  if (mode_ == "worm_like_chain" || mode_ == "worm_like_chain_linker") {
    const double chain_length = parameter_ports_[0]->get_value();
    const double persistence_length = parameter_ports_[1]->get_value();
    if (!(chain_length > 0.0)) {
      throw std::domain_error("PolymerDistances '" + get_name() +
                              "': chain_length is not positive");
    }
    // The kernel takes the dimensionless kappa; the ports carry what the
    // model fits. Same derivation as chisurf's model property.
    const double kappa = persistence_length / chain_length;
    if (mode_ == "worm_like_chain") {
      // `distance = false`, always: chisurf's forwarder has never applied
      // the r^2 factor its signature advertises -- the flag was accepted
      // and dropped on the floor, every fit in the stack was made against
      // that behaviour, and preserving the answer is the port's contract
      // (rdf.py says the same over the same call; owner decision pending
      // in PRD-105 on the flag itself).
      worm_like_chain(axis_, kappa, chain_length, true, false, &view, &n_view);
    } else {
      worm_like_chain_linker(axis_, kappa, chain_length,
                             parameter_ports_[2]->get_value(), true,
                             &view, &n_view);
    }
  } else if (mode_ == "saw_nu") {
    saw_nu(axis_, parameter_ports_[0]->get_value(),
           parameter_ports_[1]->get_value(), 1.1615, &view, &n_view);
  } else {  // ising_chain
    const int n_residues = static_cast<int>(
        std::lround(parameter_ports_[0]->get_value()));
    ising_chain(axis_, n_residues, parameter_ports_[1]->get_value(),
                parameter_ports_[2]->get_value(),
                parameter_ports_[3]->get_value(),
                parameter_ports_[4]->get_value(), n_k_, &view, &n_view);
  }
  if (view == nullptr || n_view != static_cast<int>(axis_.size())) {
    std::free(view);
    throw std::domain_error("PolymerDistances '" + get_name() +
                            "': the kernel returned no distribution");
  }
  spectrum_.resize(2 * axis_.size());
  for (std::size_t j = 0; j < axis_.size(); ++j) {
    spectrum_[2 * j] = view[j];
    spectrum_[2 * j + 1] = axis_[j];
  }
  std::free(view);
  spectrum_node_detail::publish(this, spectrum_);
  set_valid(true);
}

IMPBFF_END_NAMESPACE
