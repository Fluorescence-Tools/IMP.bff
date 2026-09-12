#ifndef IMPBFF_FCS_MODEL_SEARCH_H
#define IMPBFF_FCS_MODEL_SEARCH_H

/**
 * \file IMP/bff/FCSModelSearch.h
 * \brief Native analytical-FCS model families for model search.
 *
 * The fitting state belongs to bff: one canonical parameter registry is
 * shared by every complete objective topology.  Selecting another FCS model
 * therefore selects another graph over the same owner ports; it never copies
 * values between graphs and never asks a Python callback to evaluate a model.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/bff_config.h>
#include <IMP/bff/FitChiSquared.h>
#include <IMP/bff/GraphNode.h>
#include <IMP/bff/GraphPort.h>
#include <IMP/bff/ModelSearch.h>

#include <memory>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Controls which analytical FCS structures a factory builds.
class IMPBFFEXPORT FCSModelSearchConfig {
 public:
  FCSModelSearchConfig();

  void set_include_2d(bool value) { include_2d_ = value; }
  bool get_include_2d() const { return include_2d_; }
  void set_include_3d(bool value) { include_3d_ = value; }
  bool get_include_3d() const { return include_3d_; }

  //! One or two analytical diffusion components.
  void set_max_diffusion_components(int value);
  int get_max_diffusion_components() const {
    return max_diffusion_components_;
  }

  //! Zero or one exponential relaxation term in the initial implementation.
  void set_max_relaxation_terms(int value);
  int get_max_relaxation_terms() const { return max_relaxation_terms_; }

  void set_noise_model(FitNoiseModel value);
  FitNoiseModel get_noise_model() const { return noise_model_; }
  IMP_SHOWABLE_INLINE(
      FCSModelSearchConfig,
      out << "FCSModelSearchConfig(" << max_diffusion_components_
          << " diffusion, " << max_relaxation_terms_ << " relaxation)");

 private:
  bool include_2d_ = true;
  bool include_3d_ = true;
  int max_diffusion_components_ = 2;
  int max_relaxation_terms_ = 1;
  FitNoiseModel noise_model_ = FIT_NOISE_NEYMAN;
};

IMP_VALUES(FCSModelSearchConfig, FCSModelSearchConfigs);

//! Builds complete callback-free FCS search graphs from measured data.
/**
 * The analytical family currently spans 2-D and 3-D Gaussian diffusion, one
 * or two diffusion components, and zero or one exponential relaxation term.
 * The returned problem owns its canonical parameters, every forward-model
 * node, and every objective. MDF and photokinetic curves can therefore join
 * later as additional complete topologies without a second search engine or
 * parameter state.
 */
class IMPBFFEXPORT FCSModelSearchFactory {
 public:
  //! Build the analytical FCS family for one correlation curve.
  /**
   * `axis`, `data`, and `errors` use the caller's units.  The canonical
   * diffusion times use the same units as `axis`. Empty `errors` means unit
   * weights. `mask` and `[fit_begin, fit_end)` have the FitChiSquared
   * semantics; a negative end means the end of the data.
   */
  static std::shared_ptr<MultiStructureModelSearchProblem> create_analytical(
      const std::vector<double>& axis, const std::vector<double>& data,
      const std::vector<double>& errors,
      const FCSModelSearchConfig& config = FCSModelSearchConfig(),
      const std::vector<double>& mask = std::vector<double>(),
      int fit_begin = 0, int fit_end = -1);
};

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_FCS_MODEL_SEARCH_H
