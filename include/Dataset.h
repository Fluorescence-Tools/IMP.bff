/**
 *  \file IMP/bff/Dataset.h
 *  \brief Measured values of any rank, and the noise family that goes with them.
 *
 *  A decay is 1-D, an image is 2-D, a density is 3-D, and the only thing that
 *  differs is the shape. So there is no curve class here: a curve is the
 *  rank-1 case.
 *
 *  \par What a dataset owns, and what it does not
 *  It owns the **noise family** -- "these are counts", "these have measured
 *  standard deviations" -- and not the per-point weights. For counting data
 *  the weight is not a property of the data at all: weighting by the data,
 *  \f$1/y\f$, is Neyman and is *biased*, pulling a fit toward the low bins in
 *  a way that does not vanish as counts grow; weighting by the model,
 *  \f$1/\lambda\f$, is Pearson, and it is what the expected Fisher
 *  information uses. The model changes every iteration, so the weights do.
 *
 *  Hence #variance takes the model. A dataset that owned per-point sigmas
 *  would be right for a Gaussian measurement and quietly wrong for counts,
 *  and the failure is the worst kind: the fit still converges, to the wrong
 *  place, and a reduced chi-square near one hides it.
 *
 *  \par Everything is named for its weighting
 *  #ResidualKind says Pearson, deviance or Neyman, never "weighted". One
 *  repository was found computing "weighted residuals" two different ways,
 *  only one of which was what its own rule meant; the phrase without the
 *  qualifier is what allowed it.
 *
 *  \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_DATASET_H
#define IMPBFF_DATASET_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/Base.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! How the variance of a measured value is known.
enum NoiseFamily {
  //! Counts. \f$\mathrm{var} = \lambda\f$, the model's own prediction.
  NOISE_FAMILY_POISSON = 0,
  //! A measurement with its own per-point variance, supplied with the data.
  NOISE_FAMILY_STORED = 1,
  //! Gaussian with one variance for every point.
  NOISE_FAMILY_CONSTANT = 2
};

//! Which residual, named for how it is weighted rather than that it is.
enum ResidualKind {
  //! \f$(y-\mu)/\sqrt{\mathrm{var}(\mu)}\f$ -- weighted by the *model*.
  RESIDUAL_PEARSON = 0,
  //! Signed square root of the Poisson deviance; the likelihood's own residual.
  RESIDUAL_DEVIANCE = 1,
  //! \f$(y-\mu)/\sqrt{y}\f$ -- weighted by the *data*, and biased.
  /*! Offered because it is what a great deal of existing analysis used, and
      named so that nobody reaches for it without meaning to. */
  RESIDUAL_NEYMAN = 2
};

//! Measured values of any rank, with the noise family they were measured under.
class IMPBFFEXPORT Dataset {
 public:
  Dataset();

  //! The values, row-major, and their shape.
  /*! \throws IMP::ValueException unless the shape's product is the number of
      values, and every extent is positive. */
  void set_values(const std::vector<double>& values,
                  const std::vector<int>& shape);
  //! A rank-1 dataset: the shape is its length.
  void set_values(const std::vector<double>& values);

  //! The values, row-major.
  const std::vector<double>& get_values() const { return values_; }
  //! Extents, slowest axis first.
  std::vector<int> get_shape() const { return shape_; }
  //! How many axes: 1 for a curve, 2 for an image, 3 for a density.
  int get_rank() const { return static_cast<int>(shape_.size()); }
  //! How many values.
  int get_size() const { return static_cast<int>(values_.size()); }

  //! Which points count. Empty means all of them.
  /*! \throws IMP::ValueException if it is neither empty nor the data's size. */
  void set_mask(const std::vector<double>& mask);
  const std::vector<double>& get_mask() const { return mask_; }

  //! How the variance is known.
  void set_noise_family(NoiseFamily family) { family_ = family; }
  NoiseFamily get_noise_family() const { return family_; }

  //! Per-point variance, for #NOISE_FAMILY_STORED.
  void set_stored_variance(const std::vector<double>& variance);
  //! The single variance, for #NOISE_FAMILY_CONSTANT.
  void set_constant_variance(double variance);

  //! The variance at every point, given what the model predicts there.
  /*!
      Poisson takes it from \p model, which is why this is a function and not
      an array: the weights move with the fit.

      \throws IMP::ValueException if the family is #NOISE_FAMILY_STORED and no
              variance was stored. **Not** a silent fall back to ones: on
              counts spanning four decades, unweighted least squares fits the
              peak and ignores the tail, and the tail is usually where the
              information is.
  */
  void variance(const std::vector<double>& model,
                double** out_view, int* n_out_view) const;

  //! Residuals of \p model against these values, weighted as \p kind says.
  /*!
      A Poisson dataset answers for every kind, because one analysis wants
      more than one: a deviance per degree of freedom to report, Pearson signs
      for a runs test. Choosing once, at construction, is a choice a real
      analysis does not make.
  */
  void residuals(const std::vector<double>& model, ResidualKind kind,
                 double** out_view, int* n_out_view) const;

  //! The family's own objective: Poisson deviance, or chi-square.
  /*! Minimised in both cases, and for Poisson it is twice the difference from
      the saturated log-likelihood -- so minimising it maximises the
      likelihood, and the dropped \f$-\log y!\f$ cancels in any ratio. */
  double objective(const std::vector<double>& model) const;

  //! Points that are not masked out.
  int get_number_of_active_points() const;

  //! A line naming the shape, the family and the mask.
  std::string describe() const;

 private:
  std::vector<double> values_, mask_, stored_variance_;
  std::vector<int> shape_;
  NoiseFamily family_ = NOISE_FAMILY_POISSON;
  double constant_variance_ = 1.0;

  void check_model(const std::vector<double>& model) const;
};

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_DATASET_H
