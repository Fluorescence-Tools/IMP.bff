/**\file IMP/bff/GaussianDistances.h
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#ifndef IMPBFF_GAUSSIANDISTANCES_H
#define IMPBFF_GAUSSIANDISTANCES_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/GraphNode.h>
#include <IMP/bff/GraphPort.h>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! A distance distribution built from a sum of generalised normals.
/**
 * The producer for the FRET model everyone actually fits: a distance
 * distribution as a weighted sum of (possibly skewed) normal distributions,
 * evaluated on a fixed distance axis. Downstream it becomes a rate
 * distribution and then a lifetime spectrum, so the whole chain from a
 * fitted mean distance to a model curve is
 *
 *     GaussianDistances -> FRETSpectrumNode -> TCSPCDecay -> FitChiSquared
 *
 * and nothing in it returns to the caller.
 *
 * \par The axis is data, the components are ports
 * The distance axis is a setting -- its range and resolution do not move
 * during a fit -- so the caller hands it over once. The mean, width, skew
 * and weight of each component are what an optimiser writes, so they are
 * ports.
 *
 * \par Two normalisations, in this order, and they are not interchangeable
 * Each component is normalised to unit sum on the axis *before* it is
 * weighted, and the weighted sum is normalised again. Normalising only at
 * the end would let a wide component contribute less than its weight says
 * merely because more of it falls off the end of the axis.
 *
 * Ports, created by set_number_of_components():
 *
 * | port | what it is |
 * |---|---|
 * | `mean0`, `sigma0`, `shape0`, `amplitude0`, ... | one distance component |
 *
 * The output is interleaved `(p0, r0, p1, r1, ...)` -- weight and distance,
 * the layout `FRETSpectrumNode` reads.
 */
class IMPBFFEXPORT GaussianDistances : public GraphNode {
 public:
  explicit GaussianDistances(const std::string& name = "distances");

  //! Build `4 * n` ports named `mean0`, `sigma0`, `shape0`, `amplitude0`,
  //! `mean1`, ... -- the index is a suffix with no separator.
  void set_number_of_components(int n);
  int get_number_of_components() const { return n_components_; }

  //! The distance axis the distribution is evaluated on, in Angstrom.
  void set_axis(const std::vector<double>& axis);
  const std::vector<double>& get_axis() const { return axis_; }
  //! Numpy in, for the axis a caller holds as an array.
  void set_axis_array(double* in_axis, int n_axis);

  //! The interleaved `(p, r)` distribution from the last evaluation.
  const std::vector<double>& get_distribution() const { return spectrum_; }

  //! Use the distance-between-two-Gaussians kernel instead of the skewed one.
  /*! These are different probability densities, not two parameterisations of
      one: the two-cloud form carries an extra `r/mean` factor and an
      antisymmetric second term. The builder used to refuse a model with this
      set, which left `GaussianModel` without a graph whenever it was on. */
  void set_distance_between_gaussians(bool value) {
    distance_between_gaussians_ = value;
    set_valid(false);
  }
  bool get_distance_between_gaussians() const {
    return distance_between_gaussians_;
  }

  void evaluate() override;

 private:
  std::vector<double> axis_;
  std::vector<double> density_;
  std::vector<double> spectrum_;
  std::vector<GraphPort*> component_ports_;
  int n_components_ = 0;
  bool distance_between_gaussians_ = false;
};

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_GAUSSIANDISTANCES_H
