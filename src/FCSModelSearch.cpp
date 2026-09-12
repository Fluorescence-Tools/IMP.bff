/**
 * \file FCSModelSearch.cpp
 * \brief Native analytical-FCS model families for model search.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/FCSModelSearch.h>

#include <IMP/bff/GraphExpression.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

IMPBFF_BEGIN_NAMESPACE

namespace {

const char* kParameterIds[] = {
    "fcs.N", "fcs.baseline", "fcs.structure_parameter",
    "fcs.diffusion_time.1", "fcs.diffusion_fraction.1",
    "fcs.diffusion_time.2", "fcs.relaxation_amplitude.1",
    "fcs.relaxation_time.1"};

const std::size_t kParameterCount =
    sizeof(kParameterIds) / sizeof(kParameterIds[0]);

std::shared_ptr<GraphPort> make_parameter(double value, double lower,
                                          double upper,
                                          const std::string& name) {
  return std::make_shared<GraphPort>(value, false, false, false, true, lower,
                                     upper, GRAPH_PORT_FLOAT, name);
}

void require_finite(const std::vector<double>& values,
                    const std::string& what) {
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (!std::isfinite(values[i])) {
      throw std::invalid_argument("FCS model search: " + what +
                                  " must contain only finite values");
    }
  }
}

double median(std::vector<double> values) {
  if (values.empty()) return 0.0;
  const std::size_t middle = values.size() / 2;
  std::nth_element(values.begin(), values.begin() + middle, values.end());
  const double upper = values[middle];
  if (values.size() % 2 != 0) return upper;
  std::nth_element(values.begin(), values.begin() + middle - 1, values.end());
  return 0.5 * (values[middle - 1] + upper);
}

std::string structure_key(int dimensions, int components, int relaxations) {
  std::ostringstream out;
  out << "fcs." << dimensions << "d." << components << "diff."
      << relaxations << "relax";
  return out.str();
}

std::string diffusion_term(const std::string& time, int dimensions) {
  std::string term = "(1+x/abs(" + time + "))**(-1)";
  if (dimensions == 3) {
    term += "*(1+x/(s*s*abs(" + time + ")))**(-0.5)";
  }
  return term;
}

std::string analytical_expression(int dimensions, int components,
                                  int relaxations) {
  std::string diffusion = diffusion_term("td1", dimensions);
  if (components == 2) {
    diffusion = "a1*(" + diffusion + ")+(1-a1)*(" +
                diffusion_term("td2", dimensions) + ")";
  }
  std::string expression = "b+1/abs(N)*(" + diffusion + ")";
  if (relaxations == 1) {
    expression += "*(1-ba1+ba1*exp(-x/abs(bt1)))";
  }
  return expression;
}

std::string canonical_id_for_variable(const std::string& variable) {
  if (variable == "N") return "fcs.N";
  if (variable == "b") return "fcs.baseline";
  if (variable == "s") return "fcs.structure_parameter";
  if (variable == "td1") return "fcs.diffusion_time.1";
  if (variable == "a1") return "fcs.diffusion_fraction.1";
  if (variable == "td2") return "fcs.diffusion_time.2";
  if (variable == "ba1") return "fcs.relaxation_amplitude.1";
  if (variable == "bt1") return "fcs.relaxation_time.1";
  throw std::logic_error("FCS model search: unregistered expression variable '" +
                         variable + "'");
}

}  // namespace

FCSModelSearchConfig::FCSModelSearchConfig() {}

void FCSModelSearchConfig::set_max_diffusion_components(int value) {
  if (value < 1 || value > 2) {
    throw std::invalid_argument(
        "FCS model search supports one or two analytical diffusion "
        "components");
  }
  max_diffusion_components_ = value;
}

void FCSModelSearchConfig::set_max_relaxation_terms(int value) {
  if (value < 0 || value > 1) {
    throw std::invalid_argument(
        "FCS model search supports zero or one analytical relaxation term");
  }
  max_relaxation_terms_ = value;
}

void FCSModelSearchConfig::set_noise_model(FitNoiseModel value) {
  if (value != FIT_NOISE_NEYMAN && value != FIT_NOISE_POISSON) {
    throw std::invalid_argument("FCS model search: unknown noise model");
  }
  noise_model_ = value;
}

std::shared_ptr<MultiStructureModelSearchProblem>
FCSModelSearchFactory::create_analytical(
    const std::vector<double>& axis, const std::vector<double>& data,
    const std::vector<double>& errors, const FCSModelSearchConfig& config,
    const std::vector<double>& mask, int fit_begin, int fit_end) {
  if (axis.empty() || axis.size() != data.size()) {
    throw std::invalid_argument(
        "FCS model search: axis and data must be non-empty and equally long");
  }
  if (!errors.empty() && errors.size() != data.size()) {
    throw std::invalid_argument(
        "FCS model search: errors must be empty or as long as the data");
  }
  if (!mask.empty() && mask.size() != data.size()) {
    throw std::invalid_argument(
        "FCS model search: mask must be empty or as long as the data");
  }
  if (!config.get_include_2d() && !config.get_include_3d()) {
    throw std::invalid_argument(
        "FCS model search: at least one diffusion mode must be enabled");
  }
  if (fit_begin < 0 || fit_begin >= static_cast<int>(data.size()) ||
      (fit_end >= 0 &&
       (fit_end <= fit_begin || fit_end > static_cast<int>(data.size())))) {
    throw std::invalid_argument("FCS model search: invalid fit range");
  }
  require_finite(axis, "axis");
  require_finite(data, "data");
  require_finite(errors, "errors");
  require_finite(mask, "mask");
  if (config.get_noise_model() == FIT_NOISE_NEYMAN) {
    for (std::size_t i = 0; i < errors.size(); ++i) {
      if (!(errors[i] > 0.0)) {
        throw std::invalid_argument(
            "FCS model search: weighted errors must be positive");
      }
    }
  }

  const int end = fit_end < 0 ? static_cast<int>(data.size()) : fit_end;
  std::vector<double> fit_axis(axis.begin() + fit_begin, axis.begin() + end);
  std::vector<double> positive_axis;
  for (std::size_t i = 0; i < fit_axis.size(); ++i) {
    if (fit_axis[i] > 0.0) positive_axis.push_back(fit_axis[i]);
  }
  if (positive_axis.empty()) {
    throw std::invalid_argument(
        "FCS model search: the fit range needs a positive lag");
  }
  // The median of the last quarter still contains substantial correlation
  // for slow curves and biases both baseline and N.  The final lag is the
  // least contaminated baseline estimate available without a second model;
  // the half-amplitude crossing then gives the analytical one-component td.
  const double baseline = data[static_cast<std::size_t>(end - 1)];
  const double amplitude = data[static_cast<std::size_t>(fit_begin)] - baseline;
  const double particles = amplitude > 1e-12 ? 1.0 / amplitude : 1.0;
  double td1 = std::max(median(positive_axis), 1e-12);
  if (amplitude > 1e-12) {
    const double half = baseline + 0.5 * amplitude;
    double best = std::numeric_limits<double>::infinity();
    for (int i = fit_begin; i < end; ++i) {
      const double distance = std::fabs(data[static_cast<std::size_t>(i)] - half);
      if (distance < best && axis[static_cast<std::size_t>(i)] > 0.0) {
        best = distance;
        td1 = axis[static_cast<std::size_t>(i)];
      }
    }
  }
  const double td2 = std::max(4.0 * td1, td1 + 1e-12);

  std::shared_ptr<MultiStructureModelSearchProblem> problem =
      std::make_shared<MultiStructureModelSearchProblem>();
  std::map<std::string, std::shared_ptr<GraphPort> > parameters;
  parameters["fcs.N"] = make_parameter(particles, 1e-9, 1e12, "N");
  parameters["fcs.baseline"] =
      make_parameter(baseline, -1e6, 1e6, "b");
  parameters["fcs.structure_parameter"] =
      make_parameter(3.5, 1e-3, 1e4, "s");
  parameters["fcs.diffusion_time.1"] =
      make_parameter(td1, 1e-12, 1e12, "td1");
  parameters["fcs.diffusion_fraction.1"] =
      make_parameter(0.5, 0.0, 1.0, "a1");
  parameters["fcs.diffusion_time.2"] =
      make_parameter(td2, 1e-12, 1e12, "td2");
  parameters["fcs.relaxation_amplitude.1"] =
      make_parameter(0.1, 0.0, 1.0, "ba1");
  parameters["fcs.relaxation_time.1"] =
      make_parameter(std::max(td1 / 20.0, 1e-12), 1e-12, 1e12, "bt1");

  std::vector<std::string> parameter_ids;
  std::vector<std::shared_ptr<GraphPort> > parameter_ports;
  for (std::size_t i = 0; i < kParameterCount; ++i) {
    const std::string id(kParameterIds[i]);
    parameter_ids.push_back(id);
    parameter_ports.push_back(parameters[id]);
    problem->add_parameter(id, parameters[id]);
  }
  const std::vector<double> initial_values = {
      particles, baseline, 3.5, td1, 0.5, td2, 0.1,
      std::max(td1 / 20.0, 1e-12)};

  std::vector<int> dimensions;
  if (config.get_include_2d()) dimensions.push_back(2);
  if (config.get_include_3d()) dimensions.push_back(3);

  for (std::size_t d = 0; d < dimensions.size(); ++d) {
    for (int components = 1;
         components <= config.get_max_diffusion_components(); ++components) {
      for (int relaxations = 0;
           relaxations <= config.get_max_relaxation_terms(); ++relaxations) {
        const std::string key =
            structure_key(dimensions[d], components, relaxations);
        std::shared_ptr<GraphExpression> model =
            std::make_shared<GraphExpression>(key + ".model");
        model->set_expression(
            analytical_expression(dimensions[d], components, relaxations));
        const std::vector<std::string> variables = model->get_variable_names();
        for (std::size_t i = 0; i < variables.size(); ++i) {
          if (variables[i] == "x") {
            model->add_input_port(
                "x", std::make_shared<GraphPort>(axis, true, false, false,
                                                  false, 0.0, 0.0,
                                                  GRAPH_PORT_FLOAT_VECTOR,
                                                  "x"));
          } else {
            const std::string id = canonical_id_for_variable(variables[i]);
            std::shared_ptr<GraphPort> follower =
                std::make_shared<GraphPort>(parameters[id]->get_value());
            follower->set_link(parameters[id]);
            model->add_input_port(variables[i], follower);
          }
        }
        std::shared_ptr<GraphPort> curve = std::make_shared<GraphPort>(
            std::vector<double>(axis.size(), 0.0), false, true, false, false,
            0.0, 0.0, GRAPH_PORT_FLOAT_VECTOR, key + ".curve");
        model->add_output_port(model->get_name(), curve);

        std::shared_ptr<FitChiSquared> objective =
            std::make_shared<FitChiSquared>(key + ".chi2");
        objective->set_data(data, errors);
        objective->set_noise_model(config.get_noise_model());
        objective->set_fit_range(fit_begin, fit_end);
        objective->set_mask(mask);
        std::shared_ptr<GraphPort> model_input =
            std::make_shared<GraphPort>(std::vector<double>(axis.size(), 0.0));
        model_input->set_link(curve);
        objective->add_input_port("model", model_input);
        objective->add_output_port(
            objective->get_name(),
            std::make_shared<GraphPort>(0.0, false, true));
        objective->add_output_port(
            "residuals", std::make_shared<GraphPort>(
                             std::vector<double>(axis.size(), 0.0), false,
                             true));

        std::vector<int> fixed(kParameterCount, 1);
        fixed[0] = 0;  // N
        fixed[1] = 0;  // baseline
        fixed[3] = 0;  // first diffusion time
        if (dimensions[d] == 3) fixed[2] = 0;
        if (components == 2) {
          fixed[4] = 0;
          fixed[5] = 0;
        }
        if (relaxations == 1) {
          fixed[6] = 0;
          fixed[7] = 0;
        }
        problem->add_structure(
            key, objective, parameter_ids, parameter_ports, initial_values,
            fixed, "residuals");
        problem->add_structure_node(key, model);
        const double complexity = static_cast<double>(
            3 + (dimensions[d] == 3 ? 1 : 0) +
            (components == 2 ? 2 : 0) + (relaxations == 1 ? 2 : 0));
        problem->set_structure_bic_metadata(
            key, static_cast<double>(end - fit_begin), complexity);
      }
    }
  }

  const int initial_dimensions = config.get_include_2d() ? 2 : 3;
  problem->set_initial_structure(
      structure_key(initial_dimensions, 1, 0));

  for (std::size_t d = 0; d < dimensions.size(); ++d) {
    for (int components = 1;
         components <= config.get_max_diffusion_components(); ++components) {
      for (int relaxations = 0;
           relaxations <= config.get_max_relaxation_terms(); ++relaxations) {
        const std::string key =
            structure_key(dimensions[d], components, relaxations);
        problem->add_action(key, "fit-and-stop", key, 0.2, true);
        if (dimensions[d] == 2 && config.get_include_3d()) {
          problem->add_action(
              key, "use-3d-diffusion",
              structure_key(3, components, relaxations), 1.0);
        }
        if (components < config.get_max_diffusion_components()) {
          problem->add_action(
              key, "add-diffusion-component",
              structure_key(dimensions[d], components + 1, relaxations), 1.0);
        }
        if (relaxations < config.get_max_relaxation_terms()) {
          problem->add_action(
              key, "add-relaxation",
              structure_key(dimensions[d], components, relaxations + 1), 1.0);
        }
      }
    }
  }
  return problem;
}

IMPBFF_END_NAMESPACE
