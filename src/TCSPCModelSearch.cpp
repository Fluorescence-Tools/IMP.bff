/**
 * \file TCSPCModelSearch.cpp
 * \brief Native TCSPC lifetime model families for structural search.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/TCSPCModelSearch.h>

#include <IMP/bff/FitChiSquared.h>
#include <IMP/bff/GraphPort.h>
#include <IMP/bff/ModelSearch.h>
#include <IMP/bff/TCSPCDecay.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>

IMPBFF_BEGIN_NAMESPACE

namespace {

const char* const kScatter = "instrument.scatter";
const char* const kBackground = "instrument.background";
const char* const kN0 = "instrument.n0";
const char* const kTimeshift = "instrument.timeshift";

std::string component_id(const char* field, int index) {
  std::ostringstream out;
  out << "lifetime." << field << "." << index;
  return out.str();
}

std::string structure_key(int components) {
  std::ostringstream out;
  out << "lifetime.components." << components;
  return out.str();
}

std::shared_ptr<GraphPort> make_owner(double value, double lower,
                                      double upper) {
  if (!std::isfinite(value) || !std::isfinite(lower) ||
      !std::isfinite(upper) || lower > upper) {
    throw std::domain_error(
        "TCSPC lifetime search parameter bounds and values must be finite");
  }
  std::shared_ptr<GraphPort> port(new GraphPort(value));
  port->set_bounds(lower, upper);
  port->set_is_bounded(true);
  return port;
}

int effective_sample_size(const FitDataset& dataset) {
  const std::vector<double>& mask = dataset.get_mask();
  if (mask.empty()) return dataset.get_size();
  return static_cast<int>(std::count_if(
      mask.begin(), mask.end(), [](double value) { return value != 0.0; }));
}

}  // namespace

TCSPCLifetimeSearchSpace::TCSPCLifetimeSearchSpace() {}
TCSPCLifetimeSearchSpace::~TCSPCLifetimeSearchSpace() {}

std::shared_ptr<MultiStructureModelSearchProblem>
TCSPCLifetimeSearchSpace::get_problem() const {
  return problem_;
}

std::vector<std::string> TCSPCLifetimeSearchSpace::get_parameter_ids() const {
  return parameter_ids_;
}

std::shared_ptr<GraphPort> TCSPCLifetimeSearchSpace::get_parameter(
    const std::string& id) const {
  const std::map<std::string, std::shared_ptr<GraphPort> >::const_iterator it =
      parameters_.find(id);
  return it == parameters_.end() ? std::shared_ptr<GraphPort>() : it->second;
}

std::vector<std::string> TCSPCLifetimeSearchSpace::get_structure_keys() const {
  return problem_->get_structure_keys();
}

std::shared_ptr<TCSPCDecay> TCSPCLifetimeSearchSpace::get_decay(
    const std::string& structure) const {
  const std::map<std::string, std::shared_ptr<TCSPCDecay> >::const_iterator it =
      decays_.find(structure);
  return it == decays_.end() ? std::shared_ptr<TCSPCDecay>() : it->second;
}

std::shared_ptr<FitChiSquared> TCSPCLifetimeSearchSpace::get_objective(
    const std::string& structure) const {
  const std::map<std::string, std::shared_ptr<FitChiSquared> >::const_iterator
      it = objectives_.find(structure);
  return it == objectives_.end() ? std::shared_ptr<FitChiSquared>()
                                 : it->second;
}

TCSPCLifetimeSearchFactory::TCSPCLifetimeSearchFactory()
    : has_dataset_(false),
      channel_width_(0.0),
      excitation_period_(0.0),
      minimum_components_(1),
      maximum_components_(4),
      convolution_stop_(-1),
      stop_(-1),
      scale_start_(0),
      scale_stop_(-1),
      normalize_amplitudes_(true),
      absolute_amplitudes_(true) {}

void TCSPCLifetimeSearchFactory::set_dataset(const FitDataset& dataset) {
  dataset_ = dataset;
  has_dataset_ = true;
}

void TCSPCLifetimeSearchFactory::set_response(
    const std::vector<double>& response) {
  response_ = response;
}

void TCSPCLifetimeSearchFactory::set_timing(double channel_width,
                                            double excitation_period) {
  channel_width_ = channel_width;
  excitation_period_ = excitation_period;
}

void TCSPCLifetimeSearchFactory::set_component_range(int minimum,
                                                     int maximum) {
  minimum_components_ = minimum;
  maximum_components_ = maximum;
}

void TCSPCLifetimeSearchFactory::set_convolution_range(int convolution_stop,
                                                       int stop) {
  convolution_stop_ = convolution_stop;
  stop_ = stop;
}

void TCSPCLifetimeSearchFactory::set_scale_range(int start, int stop) {
  scale_start_ = start;
  scale_stop_ = stop;
}

void TCSPCLifetimeSearchFactory::set_normalize_amplitudes(bool value) {
  normalize_amplitudes_ = value;
}

void TCSPCLifetimeSearchFactory::set_absolute_amplitudes(bool value) {
  absolute_amplitudes_ = value;
}

void TCSPCLifetimeSearchFactory::set_parameter(const std::string& id,
                                               double initial, bool free,
                                               double lower, double upper) {
  if (id.empty()) {
    throw std::domain_error("TCSPC lifetime search parameter id is empty");
  }
  if (!std::isfinite(initial) || !std::isfinite(lower) ||
      !std::isfinite(upper) || lower > upper || initial < lower ||
      initial > upper) {
    throw std::domain_error(
        "TCSPC lifetime search parameter has invalid value or bounds");
  }
  ParameterSetting setting = {initial, free, lower, upper};
  settings_[id] = setting;
}

std::shared_ptr<TCSPCLifetimeSearchSpace>
TCSPCLifetimeSearchFactory::build() const {
  if (!has_dataset_ || dataset_.get_rank() != 1 || dataset_.get_size() <= 0) {
    throw std::domain_error(
        "TCSPC lifetime search requires a non-empty rank-1 dataset");
  }
  if (response_.size() != static_cast<std::size_t>(dataset_.get_size())) {
    throw std::domain_error(
        "TCSPC lifetime search response must match the dataset length");
  }
  if (!(channel_width_ > 0.0) || !(excitation_period_ > 0.0) ||
      !std::isfinite(channel_width_) || !std::isfinite(excitation_period_)) {
    throw std::domain_error(
        "TCSPC lifetime search timing must be finite and positive");
  }
  if (minimum_components_ < 1 || maximum_components_ < minimum_components_) {
    throw std::domain_error("invalid TCSPC lifetime component range");
  }
  for (std::size_t i = 0; i < dataset_.get_values().size(); ++i) {
    if (!std::isfinite(dataset_.get_values()[i])) {
      throw std::domain_error(
          "TCSPC lifetime search dataset values must be finite");
    }
  }
  const int samples = effective_sample_size(dataset_);
  if (samples <= 0) {
    throw std::domain_error("TCSPC lifetime search dataset mask is empty");
  }

  std::shared_ptr<TCSPCLifetimeSearchSpace> space(
      new TCSPCLifetimeSearchSpace);
  space->problem_.reset(new MultiStructureModelSearchProblem);

  std::map<std::string, ParameterSetting> defaults;
  for (int i = 0; i < maximum_components_; ++i) {
    const double amplitude = 1.0 / static_cast<double>(i + 1);
    const double tau = std::min(excitation_period_,
                                std::max(channel_width_, 1.0 + i));
    defaults[component_id("amplitude", i)] =
        ParameterSetting{amplitude, true, 0.0, 1.0};
    defaults[component_id("tau", i)] = ParameterSetting{
        tau, true, std::max(channel_width_ * 0.01, 1e-12),
        excitation_period_};
  }
  defaults[kScatter] = ParameterSetting{0.0, false, 0.0, 1.0};
  const double data_max = std::max(
      1.0, *std::max_element(dataset_.get_values().begin(),
                             dataset_.get_values().end()));
  const double data_sum = std::max(
      data_max, std::accumulate(dataset_.get_values().begin(),
                                dataset_.get_values().end(), 0.0));
  defaults[kBackground] =
      ParameterSetting{0.0, true, 0.0, data_max};
  defaults[kN0] = ParameterSetting{
      data_sum, true, 0.0, 100.0 * data_sum};
  defaults[kTimeshift] = ParameterSetting{
      0.0, false, -0.5 * response_.size(), 0.5 * response_.size()};

  for (std::map<std::string, ParameterSetting>::const_iterator it =
           settings_.begin();
       it != settings_.end(); ++it) {
    if (!defaults.count(it->first)) {
      throw std::domain_error("unknown TCSPC lifetime parameter id '" +
                              it->first + "'");
    }
    defaults[it->first] = it->second;
  }

  for (int i = 0; i < maximum_components_; ++i) {
    space->parameter_ids_.push_back(component_id("amplitude", i));
    space->parameter_ids_.push_back(component_id("tau", i));
  }
  space->parameter_ids_.push_back(kScatter);
  space->parameter_ids_.push_back(kBackground);
  space->parameter_ids_.push_back(kN0);
  space->parameter_ids_.push_back(kTimeshift);
  std::vector<std::shared_ptr<GraphPort> > owner_ports;
  owner_ports.reserve(space->parameter_ids_.size());
  for (std::size_t i = 0; i < space->parameter_ids_.size(); ++i) {
    const std::string& id = space->parameter_ids_[i];
    const ParameterSetting& setting = defaults[id];
    space->parameters_[id] =
        make_owner(setting.initial, setting.lower, setting.upper);
    owner_ports.push_back(space->parameters_[id]);
    space->problem_->add_parameter(id, space->parameters_[id]);
  }

  for (int n = minimum_components_; n <= maximum_components_; ++n) {
    const std::string key = structure_key(n);
    std::shared_ptr<TCSPCDecay> decay(new TCSPCDecay(key + ".decay"));
    decay->set_number_of_lifetimes(n);
    decay->add_output_port(decay->get_name(),
                           std::shared_ptr<GraphPort>(new GraphPort(
                               std::vector<double>(1, 0.0), false, true)));
    decay->set_response(response_);
    decay->set_timing(channel_width_, excitation_period_);
    decay->set_convolution_range(
        convolution_stop_ < 0 ? dataset_.get_size() : convolution_stop_,
        stop_ < 0 ? dataset_.get_size() : stop_);
    decay->set_scale_range(scale_start_,
                           scale_stop_ < 0 ? dataset_.get_size() : scale_stop_);
    decay->set_normalize_amplitudes(normalize_amplitudes_);
    decay->set_absolute_amplitudes(absolute_amplitudes_);

    for (int i = 0; i < n; ++i) {
      decay->get_input_port("a" + std::to_string(i))->set_link(
          space->parameters_[component_id("amplitude", i)]);
      decay->get_input_port("t" + std::to_string(i))->set_link(
          space->parameters_[component_id("tau", i)]);
    }
    decay->get_input_port("scatter")->set_link(space->parameters_[kScatter]);
    decay->get_input_port("background")->set_link(
        space->parameters_[kBackground]);
    decay->get_input_port("n0")->set_link(space->parameters_[kN0]);
    decay->get_input_port("timeshift")->set_link(
        space->parameters_[kTimeshift]);

    std::shared_ptr<FitChiSquared> objective(
        new FitChiSquared(key + ".objective"));
    objective->set_dataset(dataset_);
    std::shared_ptr<GraphPort> model(new GraphPort(std::vector<double>(1, 0.0)));
    model->set_link(decay->get_output_port(decay->get_name()));
    objective->add_input_port("model", model);
    objective->add_output_port(
        objective->get_name(),
        std::shared_ptr<GraphPort>(new GraphPort(0.0, false, true)));
    objective->add_output_port(
        "residuals", std::shared_ptr<GraphPort>(new GraphPort(
                         std::vector<double>(1, 0.0), false, true)));

    std::vector<double> initial;
    std::vector<int> fixed;
    initial.reserve(space->parameter_ids_.size());
    fixed.reserve(space->parameter_ids_.size());
    int complexity = 0;
    for (std::size_t i = 0; i < space->parameter_ids_.size(); ++i) {
      const std::string& id = space->parameter_ids_[i];
      ParameterSetting setting = defaults[id];
      bool active = setting.free;
      if (id.compare(0, 19, "lifetime.amplitude.") == 0) {
        const int index = std::stoi(id.substr(19));
        active = setting.free && index < n &&
                 !(normalize_amplitudes_ && index == 0);
        if (!settings_.count(id)) setting.initial = 1.0 / n;
      } else if (id.compare(0, 13, "lifetime.tau.") == 0) {
        const int index = std::stoi(id.substr(13));
        active = setting.free && index < n;
      }
      initial.push_back(setting.initial);
      fixed.push_back(active ? 0 : 1);
      if (active) ++complexity;
    }

    space->problem_->add_structure(key, objective, space->parameter_ids_,
                                   owner_ports, initial, fixed, "residuals");
    space->problem_->set_structure_bic_metadata(key, samples, complexity);
    space->decays_[key] = decay;
    space->objectives_[key] = objective;
  }

  for (int n = minimum_components_; n <= maximum_components_; ++n) {
    const std::string key = structure_key(n);
    space->problem_->add_action(key, "stop", key, 0.25, true);
    if (n < maximum_components_) {
      space->problem_->add_action(key, "add-component", structure_key(n + 1),
                                  0.5);
    }
    if (n > minimum_components_) {
      space->problem_->add_action(key, "remove-component",
                                  structure_key(n - 1), 0.25);
    }
  }
  space->problem_->set_initial_structure(structure_key(minimum_components_));
  return space;
}

IMPBFF_END_NAMESPACE
