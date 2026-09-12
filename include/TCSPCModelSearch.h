/**
 * \file IMP/bff/TCSPCModelSearch.h
 * \brief Native TCSPC lifetime model families for structural search.
 *
 * A search changes the number of lifetime components, and therefore changes
 * the objective graph.  Parameter identity must not change with it.  This
 * module builds one canonical registry of owner ports and links every
 * candidate graph to those owners.  There is no graph-to-graph parameter
 * copying and no application callback in either fitting or scoring.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_TCSPCMODELSEARCH_H
#define IMPBFF_TCSPCMODELSEARCH_H

#include <IMP/bff/bff_config.h>

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <IMP/bff/FitDataset.h>

IMPBFF_BEGIN_NAMESPACE

class FitChiSquared;
class GraphPort;
class MultiStructureModelSearchProblem;
class TCSPCDecay;

//! A complete, callback-free TCSPC lifetime search space.
/**
 * The space owns the search problem, canonical parameter ports, and every
 * candidate graph.  Keeping those together is intentional: graph ports hold
 * their nodes weakly to avoid cycles, while a model family must keep all its
 * alternative topologies alive for the duration of a search.
 */
class IMPBFFEXPORT TCSPCLifetimeSearchSpace {
 public:
  ~TCSPCLifetimeSearchSpace();

  std::shared_ptr<MultiStructureModelSearchProblem> get_problem() const;
  std::vector<std::string> get_parameter_ids() const;
  std::shared_ptr<GraphPort> get_parameter(const std::string& id) const;
  std::vector<std::string> get_structure_keys() const;
  std::shared_ptr<TCSPCDecay> get_decay(const std::string& structure) const;
  std::shared_ptr<FitChiSquared> get_objective(
      const std::string& structure) const;

 private:
  friend class TCSPCLifetimeSearchFactory;
  TCSPCLifetimeSearchSpace();

  std::shared_ptr<MultiStructureModelSearchProblem> problem_;
  std::vector<std::string> parameter_ids_;
  std::map<std::string, std::shared_ptr<GraphPort> > parameters_;
  std::map<std::string, std::shared_ptr<TCSPCDecay> > decays_;
  std::map<std::string, std::shared_ptr<FitChiSquared> > objectives_;
};

//! Builds native multi-exponential TCSPC candidate graphs.
/**
 * The first model family is deliberately narrow: a plain lifetime spectrum
 * observed through a measured IRF.  The instrument parameters are already
 * canonical registry entries, so adding alternative IRF or nuisance
 * topologies later does not require changing parameter identity or the
 * search machinery.
 */
class IMPBFFEXPORT TCSPCLifetimeSearchFactory {
 public:
  TCSPCLifetimeSearchFactory();

  void set_dataset(const FitDataset& dataset);
  void set_response(const std::vector<double>& response);
  void set_timing(double channel_width, double excitation_period);
  void set_component_range(int minimum, int maximum);
  void set_convolution_range(int convolution_stop, int stop);
  void set_scale_range(int start, int stop);

  void set_normalize_amplitudes(bool value);
  void set_absolute_amplitudes(bool value);

  //! Configure one canonical parameter before build().
  /**
   * Known ids are `lifetime.amplitude.N`, `lifetime.tau.N`,
   * `instrument.scatter`, `instrument.background`, `instrument.n0`, and
   * `instrument.timeshift`.  Component ids are zero based and must be below
   * the configured maximum component count.
   */
  void set_parameter(const std::string& id, double initial, bool free,
                     double lower, double upper);

  //! Build all component-count topologies and their native BIC-scored problem.
  std::shared_ptr<TCSPCLifetimeSearchSpace> build() const;

 private:
  struct ParameterSetting {
    double initial;
    bool free;
    double lower;
    double upper;
  };

  FitDataset dataset_;
  bool has_dataset_;
  std::vector<double> response_;
  double channel_width_;
  double excitation_period_;
  int minimum_components_;
  int maximum_components_;
  int convolution_stop_;
  int stop_;
  int scale_start_;
  int scale_stop_;
  bool normalize_amplitudes_;
  bool absolute_amplitudes_;
  std::map<std::string, ParameterSetting> settings_;
};

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_TCSPCMODELSEARCH_H
