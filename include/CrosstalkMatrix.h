/**
 *  \file IMP/bff/CrosstalkMatrix.h
 *  \brief The excitation and emission crosstalk matrices and their algebra.
 *
 *  The two matrices a light path is described by — the **excitation** matrix
 *  ``excitation[l, k]``, the relative rate at which laser ``l`` excites
 *  chromophore ``k``, and the **emission** matrix ``emission[k, m]``, the
 *  detected brightness of chromophore ``k`` in detection channel ``m``
 *  (emission spectrum x filter transmission x detector efficiency x quantum
 *  yield) — together carry everything the scalar Hellenkamp factors
 *  (``alpha`` leakage, ``beta`` excitation flux, ``gamma`` detection/QY,
 *  ``delta`` direct excitation) summarise in two colours, and everything they
 *  cannot express: any number of chromophores and lasers, and arbitrary
 *  inter-channel bleed including acceptor-to-acceptor leakage.
 *
 *  Moved here from ChiSurf's ``core/fluorescence/crosstalk.py`` (owner,
 *  2026-09-04: the excitation and emission crosstalk matrix definition lives
 *  in bff, not in chisurf) — this is the same placement the FCS forward model
 *  took (FcsMdf.h). The matrix *definition* is the part that two
 *  implementations disagree on silently: the label convention, what a row and
 *  a column mean, how a payload is ordered into values, and what the forward
 *  and inverse mixing do. ChiSurf keeps the light-path calculator that
 *  *builds* a payload and the views that *display* corrected signals; the
 *  algebra in between is this.
 *
 *  Convention. A crosstalk matrix maps **sources** (rows: lasers for the
 *  excitation matrix, chromophores for the emission matrix) to **detectors**
 *  (columns). ``M[i, j]`` is the contribution of source ``i`` to detector
 *  ``j``; the forward model is ``measured_j = sum_i source_i * M[i, j]``
 *  (``measured = M^T @ sources``) and the inverse recovers the sources from
 *  the measured detector signals. Rows and columns are *labelled* — the
 *  labels are how a payload built against one instrument description is
 *  ordered against another consumer's expectation, and a requested label the
 *  matrix does not carry contributes zeros rather than an error, because a
 *  missing element of a light path is a dark element, not a broken one.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#ifndef IMPBFF_CROSSTALKMATRIX_H
#define IMPBFF_CROSSTALKMATRIX_H

#include <IMP/bff/bff_config.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! A labelled crosstalk matrix: sources as rows, detectors as columns.
/**
    The one type both crosstalk matrices share — an excitation matrix and an
    emission matrix are the same thing at different ends of the light path,
    and giving them one definition is what stops a payload built by one tool
    and consumed by another from agreeing on nothing but the shape.

    The values are row-major, ``values[i * n_columns() + j]`` the contribution
    of row ``i`` to column ``j``. Labels name the rows and columns; \select
    re-orders and sub-samples against them, filling a requested-but-missing
    label with zeros.
 */
class IMPBFFEXPORT CrosstalkMatrix {
 public:
  //! The empty matrix.
  CrosstalkMatrix();

  //! ``n_rows x n_columns`` values, row-major; throws on a size mismatch.
  CrosstalkMatrix(const std::vector<std::string>& rows,
                  const std::vector<std::string>& columns,
                  const std::vector<double>& values);

  int get_n_rows() const { return static_cast<int>(rows_.size()); }
  int get_n_columns() const { return static_cast<int>(columns_.size()); }

  const std::vector<std::string>& get_rows() const { return rows_; }
  const std::vector<std::string>& get_columns() const { return columns_; }
  const std::vector<double>& get_values() const { return values_; }

  //! The contribution of row \p row to column \p column, by position.
  double get(int row, int column) const;

  //! The same, by label; a label the matrix does not carry contributes 0.
  double get(const std::string& row, const std::string& column) const;

  //! The sub-matrix in the requested label order.
  /**
      Rows and columns are re-ordered to match, and a requested label this
      matrix does not carry contributes an all-zero row or column — the same
      rule `get` follows, because sub-setting against a fuller instrument
      description is how one payload serves several consumers.
   */
  CrosstalkMatrix select(const std::vector<std::string>& rows,
                         const std::vector<std::string>& columns) const;

  //! Multiply one row (a source) by \p factor.
  /**
      How a per-source quantum yield folds into the emission matrix: the
      detected brightness of a chromophore scales with what it emits.
   */
  void scale_row(const std::string& label, double factor);

  //! Multiply one column (a detector) by \p factor.
  /**
      How a per-detector efficiency folds into the emission matrix: every
      chromophore's brightness in that channel scales with how well the
      channel detects.
   */
  void scale_column(const std::string& label, double factor);

 private:
  //! Position of \p label, or -1.
  int index_of(const std::vector<std::string>& labels,
               const std::string& label) const;

  std::vector<std::string> rows_;
  std::vector<std::string> columns_;
  std::vector<double> values_;
};

//! Forward mixing: predict detector signals from source signals.
/**
    ``out[d * n_items + ...] = sum_i sources[i] * M[i, d]`` — the measured =
    M^T @ sources of the convention above, evaluated over \p n_items
    independent items (bursts, pixels) stored row-major as
    ``sources[i * n_items + item]``. With one item this is a plain
    matrix-vector product.

    \param[in] matrix the ``(n_sources, n_detectors)`` crosstalk matrix,
               row-major
    \param[in] n_sources,n_detectors its shape; ``n_sources * n_detectors``
               must equal \p matrix 's size
    \param[in] sources ``n_sources * n_items`` source signals; ``n_items`` is
               inferred and must be an exact division
    \param[out] out_view,n_out_view ``n_detectors * n_items`` detector
               signals, allocated here (a managed view — the wrapper owns it)
 */
IMPBFFEXPORT void crosstalk_apply_mixing(
        const std::vector<double>& matrix, int n_sources, int n_detectors,
        const std::vector<double>& sources, double** out_view,
        int* n_out_view);

//! Inverse mixing: recover source signals from measured detector signals.
/**
    The inverse of `crosstalk_apply_mixing`: given the detector signals,
    estimate the sources that produced them, over \p n_items independent
    items stored row-major as ``measured[d * n_items + item]``.

    Three solves share this one implementation, because they are three
    regularisations of one problem:

    - the plain least-squares solution (minimum-norm, via a complete
      orthogonal decomposition — the pseudo-inverse a numpy caller would
      reach for), which is fast but returns negative sources and amplifies
      noise when the matrix is ill-conditioned (strong spectral overlap);
    - ``ridge > 0`` adds Tikhonov regularisation ``x = (A^T A + lambda I)^-1
      A^T y``, damping that amplification at the cost of a small bias;
    - ``nonneg`` solves non-negative least squares per item (Lawson-Hanson),
      the physically-constrained unmixing — sources cannot be negative —
      that stays stable where the plain solve blows up. ``ridge`` applies to
      it too, as an augmented system.

    \param[in] matrix the ``(n_sources, n_detectors)`` crosstalk matrix,
               row-major
    \param[in] n_sources,n_detectors its shape
    \param[in] measured ``n_detectors * n_items`` measured signals
    \param[in] nonneg constrain the solution non-negative (NNLS)
    \param[in] ridge Tikhonov strength; 0 is unregularised
    \param[out] out_view,n_out_view ``n_sources * n_items`` recovered source
               signals, allocated here (a managed view)
 */
IMPBFFEXPORT void crosstalk_invert_mixing(
        const std::vector<double>& matrix, int n_sources, int n_detectors,
        const std::vector<double>& measured, bool nonneg, double ridge,
        double** out_view, int* n_out_view);

//! Integer, statistics-preserving unmixing by photon reassignment.
/**
    The least-squares inverses return *fractional* — and unconstrained,
    negative — source estimates, which destroys the integer, Poisson nature
    of photon-counting data. This instead **reassigns each detected photon to
    a source**: every photon counted in detector ``d`` is attributed to
    exactly one source by a multinomial draw, so the output is a non-negative
    integer per-source stream whose total equals the input's, elementwise.

    The probability that a photon in detector ``d`` came from source ``k`` is
    P(k | d) proportional to ``a_k * B[k, d]`` with ``B`` the row-normalised
    matrix (each source's spectral shape as a distribution over detectors)
    and ``a_k`` the source abundance — estimated per item by non-negative
    least squares when not supplied. A multinomial thinning of a Poisson
    count is itself Poisson, so burst-variance and maximum-likelihood
    analyses downstream see genuine photon-counting data, and the expected
    assignment equals the soft (Richardson-Lucy / EM) unmixing.

    \param[in] matrix the ``(n_sources, n_detectors)`` crosstalk matrix,
               row-major
    \param[in] n_sources,n_detectors its shape
    \param[in] counts ``n_detectors * n_items`` **integer** photon counts per
               detector; non-integral values are truncated, as a counter
               would truncate them
    \param[in] abundances ``n_sources * n_items`` abundance prior, or empty
               to estimate them by non-negative least squares
    \param[in] seed seed of the random generator; identical seeds give
               identical shuffles
    \param[out] out_view,n_out_view ``n_sources * n_items`` non-negative
               integer counts, allocated here (a managed view)
 */
IMPBFFEXPORT void crosstalk_shuffle_unmix(
        const std::vector<double>& matrix, int n_sources, int n_detectors,
        const std::vector<double>& counts, const std::vector<double>& abundances,
        unsigned long long seed, double** out_view, int* n_out_view);

IMPBFF_END_NAMESPACE

#endif /* IMPBFF_CROSSTALKMATRIX_H */
