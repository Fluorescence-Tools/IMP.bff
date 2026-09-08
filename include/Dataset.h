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
  NOISE_FAMILY_CONSTANT = 2,
  //! A variance carried here by propagation from other datasets.
  /*! Set by #set_linear_combination, and stored rather than derived from the
      model for the reason given there: the components' own predictions are
      not recoverable from the combination. */
  NOISE_FAMILY_PROPAGATED = 3
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

  //! The coordinate along one axis: the `x` of a curve, a pixel edge, a lag.
  /*!
      Optional, and per dimension, because a curve without its x is not a
      curve and a 2-D dataset may need one axis and not the other.

      \param[in] dimension which axis, 0 for the slowest
      \param[in] values as long as that dimension's extent
      \throws IMP::ValueException if the dimension does not exist or the
              length is wrong.
  */
  void set_axis(int dimension, const std::vector<double>& values);
  //! The coordinate along one axis, empty if none was given.
  const std::vector<double>& get_axis(int dimension) const;
  //! Whether an axis was given for this dimension.
  bool get_has_axis(int dimension) const;

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

  //! Declare this dataset an independent source of uncertainty.
  /*!
      A measured channel. Its variance is its family's -- counts for Poisson,
      the stored array otherwise -- and every dataset derived from it tracks
      \f$\partial\,\mathrm{value}/\partial\,\mathrm{this}\f$, so that a
      later combination knows what it was built from.

      \param[in] name what it is called in #describe and in error messages
  */
  void set_as_source(const std::string& name);

  //! `a + b`, `a - b`, `a * b`, `a / b`, elementwise, with the variance carried through.
  /*!
      First-order propagation over the **independent sources**, not over the
      operands:
      \f[ \mathrm{Var}(z) = \sum_s \left(\frac{\partial z}{\partial s}\right)^2
          \mathrm{Var}(s) \f]
      with the derivatives composed by the chain rule as the expression is
      built. Propagating over the operands instead is the standard mistake and
      it is wrong whenever they share a source.

      **Anisotropy is exactly that case.**
      \f$r = (VV - G\,VH)/(VV + 2G\,VH)\f$ has the same two channels above
      and below the line, so numerator and denominator are strongly
      correlated; treating them as independent gets the variance wrong in a
      way no amount of care about each half will fix. Tracing back to
      \f$VV\f$ and \f$VH\f$ gets it right without the caller thinking about
      it.

      \throws IMP::ValueException if the operands differ in size, or if an
              operand carries no uncertainty at all -- an untracked operand
              would silently contribute none.
  */
  static Dataset add(const Dataset& a, const Dataset& b);
  static Dataset subtract(const Dataset& a, const Dataset& b);
  static Dataset multiply(const Dataset& a, const Dataset& b);
  static Dataset divide(const Dataset& a, const Dataset& b);
  //! `scale * a + offset`; the offset is exact and the scale is a constant.
  static Dataset affine(const Dataset& a, double scale, double offset = 0.0);
  //! `f(a)` for a function whose derivative the caller supplies per point.
  /*! The general case, so that a transform this class does not know about --
      a logarithm, a power, an instrument linearisation -- still carries its
      uncertainty. \p values and \p derivative are `f(a_i)` and
      `f'(a_i)`. */
  static Dataset transform(const Dataset& a, const std::vector<double>& values,
                           const std::vector<double>& derivative);

  //! Names of the independent sources this dataset was built from.
  std::vector<std::string> get_source_names() const { return source_names_; }

  //! Build this dataset as `sum_i coefficient_i * source_i`, propagating the variance.
  /*!
      The case that motivates it is the constructed magic-angle decay, which
      is not measured at 54.7 degrees but assembled from a parallel and a
      perpendicular channel:

          I = I_VV + 2 G I_VH

      **The result is not Poisson, and treating it as though it were is
      wrong.** Each channel is counts, so \f$\mathrm{Var} = \sum_i c_i^2
      \mathrm{Var}_i\f$ and the combination's variance is
      \f$I_{VV} + 4G^2 I_{VH}\f$ -- while a Poisson reading of the result
      would say \f$I_{VV} + 2G I_{VH}\f$. The perpendicular channel's
      contribution is understated by a factor of \f$2G\f$, and the misfit is
      pushed wherever the perpendicular counts are largest. The coefficient
      squares; the Poisson assumption does not know that.

      A Poisson source contributes its **measured** values as its variance,
      because the combination does not let a model be decomposed back into
      per-channel predictions. That is data weighting, with the bias described
      on #NOISE_FAMILY_POISSON, and it is the price of constructing the
      combination at all. **The unbiased alternative is not to construct it**:
      keep the channels as two datasets, give each its Poisson family, and
      let one objective hold both. That the abstraction supports the honest
      route is the argument for preferring it.

      \throws IMP::ValueException unless there is at least one source, the
              counts match, and every source has the same size.
  */
  void set_linear_combination(const std::vector<Dataset>& sources,
                              const std::vector<double>& coefficients);

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
  std::vector<double> values_, mask_, stored_variance_, propagated_variance_;
  std::string provenance_;
  //! The independent sources, their per-point variance, and d(value)/d(source).
  std::vector<std::vector<double> > axes_;
  std::vector<std::string> source_names_;
  std::vector<std::vector<double> > source_variance_, derivative_;

  //! Own variance, whatever the family says, for use in propagation.
  std::vector<double> own_variance() const;
  //! Recompute values_/propagated_variance_ from the sources.
  void finish_propagation(const std::vector<double>& values,
                          const std::string& provenance);
  static Dataset binary(const Dataset& a, const Dataset& b, int op,
                        const char* symbol);
  std::vector<int> shape_;
  NoiseFamily family_ = NOISE_FAMILY_POISSON;
  double constant_variance_ = 1.0;

  void check_model(const std::vector<double>& model) const;
};

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_DATASET_H
