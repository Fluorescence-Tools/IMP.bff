/**
 * \file IMP/bff/FitChiSquared.h
 * \brief The data misfit of a model curve, as a node in the model graph.
 *
 * Copyright 2007-2023 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_FITCHISQUARED_H
#define IMPBFF_FITCHISQUARED_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/FitDataset.h>

#include <string>
#include <vector>

#include <IMP/bff/GraphNode.h>
#include <IMP/bff/GraphPort.h>

IMPBFF_BEGIN_NAMESPACE

//! How the residual between data and model is weighted.
enum FitNoiseModel {
  //! Weighted least squares, ``(data - model) / error``. ChiSurf's
  //! ``"default"`` noise model, i.e. Neyman chi-square.
  FIT_NOISE_NEYMAN = 0,
  //! Signed Poisson deviance residuals, the maximum-likelihood estimator
  //! for low counts. ChiSurf's ``"poisson"`` noise model.
  FIT_NOISE_POISSON = 1
};

//! The chi-square of a model curve against measured data.
/**
 * A fit's objective, evaluated entirely in C++ so that a sampler never has
 * to cross into Python to learn how badly a parameter vector fits.
 *
 * The node takes the model curve on an input port and writes chi-square to
 * an output port, which is exactly the shape `MCMCSampler::set_objective()`
 * consumes: a whole fit becomes one node graph, and a move costs no
 * interpreter at all.
 *
 * The arithmetic is ChiSurf's, kept identical on purpose:
 * `calculate_weighted_residuals` over the index window ``[xmin, xmax)``,
 * truncated to the shorter of data and model, optionally multiplied by a
 * fit mask, then squared and summed with NaN counting as an infinite
 * misfit. Priors are deliberately *not* part of it -- a reported chi-square
 * is the data misfit alone, and `MCMCSampler` adds the prior separately.
 *
 * Argument errors are thrown as ``std::domain_error``, which is what bff
 * surfaces to Python as ``ValueError``.
 *
 * \see MCMCSampler, GraphNode
 */
class IMPBFFEXPORT FitChiSquared : public GraphNode {
 public:
  explicit FitChiSquared(const std::string& name = "chi2");

  //! Set the measured curve and its per-point errors.
  /** \param y measured values
      \param ey per-point errors; only used by the Neyman noise model */
  //! Score against a FitDataset, which carries its own noise family.
  /*!
      The alternative to #set_noise_model, and the one that composes: the
      dataset says whether it is counts or a measurement with stored
      variances, and the residuals follow from that rather than from a
      setting on the objective. A #FitJointChiSquared whose members carry
      datasets can hold a Poisson decay and a Gaussian correlation curve at
      once, which it cannot when each member is *told* its noise model by
      whoever built it.

      The dataset's own mask applies. The index window does not, and setting
      both is refused rather than resolved: two masking mechanisms that
      disagree is how the wrong points get excluded quietly.

      \note Residuals from this path follow the dataset's sign convention --
      positive where the data exceeds the model, for every kind. The
      \p noise_model path's `"poisson"` residuals have the opposite sign to
      its own `"default"` ones; that is a defect of long standing which
      chisurf plots today, so it is left alone rather than flipped underneath
      it. See PRD-140.
  */
  void set_dataset(const FitDataset& dataset);
  //! Whether a dataset was given.
  bool get_has_dataset() const { return has_dataset_; }
  //! The dataset, if one was given.
  const FitDataset& get_dataset() const { return dataset_; }
  //! Which residual the dataset path reports; the family's own by default.
  void set_residual_kind(FitResidualKind kind) { residual_kind_ = kind; has_residual_kind_ = true; }

  void set_data(const std::vector<double>& y, const std::vector<double>& ey);

  const std::vector<double>& get_data() const { return data_y_; }
  const std::vector<double>& get_errors() const { return data_ey_; }

  //! Restrict the fit to the index window ``[xmin, xmax)``.
  void set_fit_range(int xmin, int xmax);
  int get_fit_range_min() const { return xmin_; }
  int get_fit_range_max() const { return xmax_; }

  //! Multiplicative per-point weights, in full data indexing.
  /** Applied only when the window length matches the residual length,
      which is ChiSurf's guard against a mask that no longer fits. An
      empty mask disables masking. */
  void set_mask(const std::vector<double>& mask);
  const std::vector<double>& get_mask() const { return mask_; }

  void set_noise_model(FitNoiseModel model);
  FitNoiseModel get_noise_model() const { return noise_model_; }

  //! Set the noise model by ChiSurf's name: "default" or "poisson".
  void set_noise_model_name(const std::string& name);
  std::string get_noise_model_name() const;

  //! The key of the input port carrying the model curve.
  void set_model_port_key(const std::string& key) { model_key_ = key; }
  const std::string& get_model_port_key() const { return model_key_; }

  //! The key of the output port the residual vector is written to.
  /** Written by `evaluate()` **only when the node has such a port**, so
      nothing changes for a graph that only wants chi-square. It exists
      because `FitMinimizer` needs the residuals rather than their sum, and
      needs them from *any* objective node -- a `FitChiSquared`, or a Python
      `GraphNode` director wrapping a model this library cannot represent. A
      port is the one thing both can present. */
  void set_residuals_port_key(const std::string& key) { residuals_key_ = key; }
  const std::string& get_residuals_port_key() const { return residuals_key_; }

  //! Weighted residuals from the last evaluation.
  const std::vector<double>& get_weighted_residuals() const { return wres_; }

  //! Chi-square from the last evaluation.
  double get_chi2() const { return chi2_; }

  //! Reduced chi-square, ``chi2 / (n_points - n_free - 1)``.
  double get_chi2r(int n_free) const;

  //! Number of residuals the last evaluation produced.
  unsigned int get_number_of_residuals() const {
    return static_cast<unsigned int>(wres_.size());
  }

  //! Compute the residuals of an explicit model curve without a graph.
  /** The same arithmetic `evaluate()` runs, for callers that hold a model
      curve directly rather than on a port. */
  std::vector<double> compute_weighted_residuals(
      const std::vector<double>& model_y) const;

  //! Set the measured curve and its errors from numpy, copying nothing extra.
  /** \param in_data_y the measured curve
      \param n_data_y its length
      \param in_data_ey the per-point errors, same length
      \param n_data_ey its length

      The `std::vector` overload above means a Python caller builds two lists
      of a few thousand floats before any arithmetic happens, which for a
      per-iteration path costs more than the residual it is setting up for.
      Here the buffers are read where they already lie. */
  void set_data_arrays(double* in_data_y, int n_data_y,
                       double* in_data_ey, int n_data_ey);

  //! Set the multiplicative per-point mask from numpy. Empty clears it.
  void set_mask_array(double* in_mask_a, int n_mask_a);

  //! Weighted residuals of a model curve handed in as numpy.
  /** \param in_model_y the model curve, in full data indexing
      \param n_model_y its length
      \param out_wres the residuals over the fit window
      \param n_out_wres how many

      The point of the whole class for ChiSurf: this is model-agnostic. It
      needs the model's *curve*, not the model, so every model -- whether its
      curve was computed by this library, by numpy, or by anything else --
      gets the residual formed in one pass with a single allocation, instead
      of the slice, subtract, divide, copy and mask that each allocate again.

      \see compute_weighted_residuals, which is the same arithmetic for a
      caller that already holds `std::vector`. */
  void compute_weighted_residuals_array(double* in_model_y, int n_model_y,
                                        double** out_wres,
                                        int* n_out_wres) const;

  //! Compute chi-square of an explicit model curve without a graph.
  double compute_chi2(const std::vector<double>& model_y) const;

  //! Read the model curve from the input port, write chi-square out.
  void evaluate() override;

  //! Clear sanitising on this node's numeric transport, then pull as usual.
  /** The model input is written by `GraphNode::update()` *before* `evaluate()`
      runs, so the flag cannot be set there: a NaN curve would already have
      been floored by the time this node saw it. */
  void update() override;

  std::string describe() const;

 private:
  FitDataset dataset_;
  bool has_dataset_ = false;
  FitResidualKind residual_kind_ = RESIDUAL_PEARSON;
  bool has_residual_kind_ = false;
  std::vector<double> data_y_;
  std::vector<double> data_ey_;
  std::vector<double> mask_;
  std::vector<double> wres_;
  std::string model_key_ = "model";
  std::string residuals_key_ = "residuals";
  double chi2_ = 0.0;
  int xmin_ = 0;
  int xmax_ = -1;  //!< -1 means "to the end of the data"
  FitNoiseModel noise_model_ = FIT_NOISE_NEYMAN;

  //! The half-open window actually used, clamped to the data.
  void resolve_window(int* begin, int* end) const;
};

//! Weighted residuals of a model curve against data, in one pass.
/** \param in_data_y the measured curve
    \param n_data_y its length
    \param in_data_ey per-point errors; may be empty, then treated as 1
    \param n_data_ey its length
    \param in_model_y the model curve, in the same indexing as the data
    \param n_model_y its length
    \param xmin first index of the fit window
    \param xmax one past the last; negative means "to the end of the data"
    \param noise_model "default" (weighted least squares) or "poisson"
    \param out_wres the residuals
    \param n_out_wres how many

    **Stateless, and deliberately so.** The obvious design is a `FitChiSquared`
    holding the data and reused across iterations, but that copies the data
    into the node, and a caller editing its curve in place would then be
    served residuals against a stale copy. Reading every buffer where it lies
    removes that hazard, and it turned out to be the faster arrangement too:
    the cache lookup guarding such a node cost several times the arithmetic
    it protected.

    **Model-agnostic.** It needs the model's *curve*, not the model, so it
    serves every model there is -- whatever computed the curve, the residual
    afterwards is the same operation.

    Semantics follow ChiSurf's `calculate_weighted_residuals` exactly,
    including that a model which stops early is truncated against rather than
    treated as an error.
 */
IMPBFFEXPORT void fit_weighted_residuals(
    double* in_data_y, int n_data_y,
    double* in_data_ey, int n_data_ey,
    double* in_model_y, int n_model_y,
    int xmin, int xmax, const std::string& noise_model,
    double** out_wres, int* n_out_wres);

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_FITCHISQUARED_H
