/**
 * \file IMP/bff/GraphExpression.h
 * \brief A vectorised arithmetic expression, compiled once and evaluated as a node.
 *
 * Copyright 2007-2023 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_GRAPHEXPRESSION_H
#define IMPBFF_GRAPHEXPRESSION_H

#include <IMP/bff/bff_config.h>

#include <memory>
#include <string>
#include <vector>

#include <IMP/bff/GraphNode.h>
#include <IMP/bff/GraphPort.h>

IMPBFF_BEGIN_NAMESPACE

//! An arithmetic model equation, compiled once and evaluated in C++.
/**
 * ChiSurf's parse models are equation strings such as
 * ``b+1/abs(N)*(1+x/td)**(-1)/sqrt(1+1/s**2*x/td)``, evaluated per step with
 * Python's ``eval`` over numpy arrays. That costs an interpreter round trip
 * on every sampler move, which for a cheap model is most of the move.
 *
 * This compiles the same string once into reverse Polish form and evaluates
 * it over vectors, so the equation becomes an ordinary node in the model
 * graph: its free names are input ports, its curve is the output port, and
 * a `FitChiSquared` node downstream turns it into a fit that never leaves C++.
 *
 * Semantics follow numpy, because that is what the equations were written
 * against: scalars broadcast against vectors, ``**`` binds tighter than
 * unary minus and is right-associative, and division by zero yields
 * infinity rather than throwing.
 *
 * The evaluator itself is **not here**. It is ``pto::ExpressionEngine``
 * from ptolib (https://github.com/tpeulen/ptolib), vendored verbatim as
 * ``include/internal/ptolib.h``, and this class is a `GraphNode` wrapper over
 * it. That is not an accident of history: tttrlib's `DataStore` gates burst
 * columns with the same evaluator, neither library may depend on the other,
 * so the one implementation lives in the header both carry. imp.bff had
 * its own copy until 2026-08-31 and a stale byte-copy of the engine until
 * 2026-09-08; two evaluators of the same language is one too many.
 *
 * An equation the engine cannot compile is refused as ``std::domain_error``
 * (which bff surfaces to Python as ``ValueError``) so a caller can fall
 * back rather than get a wrong curve. There is no second evaluator behind
 * this one to fall through to any more -- which is the point. The old
 * fallback did not merely duplicate the engine, it answered
 * *wrongly*: ExprTk evaluated a multi-argument function at element 0 and
 * broadcast the result, so ``hypot(x,y)`` came back constant, silently.
 *
 * \see FitChiSquared, GraphNode, MCMCSampler
 */

class IMPBFFEXPORT GraphExpression : public GraphNode {
 public:
  explicit GraphExpression(const std::string& name = "expression");

  //! Compile an equation string; throws if it cannot be represented.
  void set_expression(const std::string& expression);
  const std::string& get_expression() const { return expression_; }

  //! The equation in the engine's spelling, after the Python rewrites.
  /** ``**`` becomes ``^``, and ``&``, ``|``, ``~`` become ``and``, ``or``,
      ``not``. Produced by ``ExpressionEngine::normalise``, so a caller that
      holds a second evaluator hands it the same string and the two cannot
      drift. */
  const std::string& get_translated_expression() const { return translated_; }

  //! The free variable names the equation uses, in first-appearance order.
  /** These are exactly the input ports the node expects. */
  std::vector<std::string> get_variable_names() const { return variables_; }

  //! Evaluate with variables supplied directly, bypassing the ports.
  /** Vectors and scalars may be mixed; a scalar broadcasts, exactly as in
      numpy. One call evaluates the whole batch, so a caller with a column
      of a million rows pays one call rather than a million. */
  std::vector<double> compute(
      const std::vector<std::string>& names,
      const std::vector<std::vector<double> >& values) const;

  //! Evaluate many parameter sets at once, one result row per set.
  /** \param names the variable names each row supplies, in order
      \param rows one vector of values per row, each as long as \p names
      \return one result per row

      For callers that sweep a parameter grid or score a table of
      candidates: the equation is compiled once and every row reuses it. */
  std::vector<double> compute_batch(
      const std::vector<std::string>& names,
      const std::vector<std::vector<double> >& rows) const;

  //! Evaluate over columns handed in as one contiguous array.
  /** \param names one name per column, in row order of \p in_columns
      \param in_columns ``n_vars x n_rows``, row *i* being column
             \c names[i]
      \param n_vars number of columns supplied
      \param n_rows length of each column
      \param out_values the result, one entry per row
      \param n_out_values length of the result

      The entry point for whole-table work. ``compute()`` takes
      ``std::vector``, which from Python means building a list of a million
      floats per column before any arithmetic happens -- for a table query
      that conversion costs an order of magnitude more than the evaluation.
      Here the buffer is read where it already lies. */
  void compute_columns(const std::vector<std::string>& names,
                       double* in_columns, int n_vars, int n_rows,
                       double** out_values, int* n_out_values) const;

  //! Evaluate a boolean query over columns, answering with a mask.
  /** \param names one name per column, in row order of \p in_columns
      \param in_columns ``n_vars x n_rows``
      \param n_vars number of columns supplied
      \param n_rows length of each column
      \param out_mask one byte per row, 1 where the query holds
      \param n_out_mask length of the mask

      The answer to a gate is a boolean, and returning it as an array of
      doubles costs eight bytes a row to say one bit. This is the form a
      selection actually wants. */
  void compute_mask(const std::vector<std::string>& names,
                    double* in_columns, int n_vars, int n_rows,
                    unsigned char** out_mask, int* n_out_mask) const;

  //! Evaluate a model curve: scalar parameters and one axis, nothing copied.
  /** \param parameter_names one name per scalar, in order of \p in_parameters
      \param in_parameters one value per scalar parameter
      \param n_parameters how many scalars
      \param axis_name the name the equation gives the curve axis
      \param in_axis the axis itself
      \param n_axis its length, and the length of the result
      \param out_values the curve
      \param n_out_values its length

      This is the shape a fit actually has, and neither of the other entry
      points fits it. ``compute()`` takes ``std::vector``, so from Python
      every call converts the axis to a list -- measured at 130 us against
      numpy's 30 for a 4096-point curve, which is four times *slower* than
      the interpreter it exists to replace. ``compute_columns()`` reads numpy
      directly but wants every operand the same length, so a scalar parameter
      has to be materialised as a full column, and the memcpy per operand
      costs more than the arithmetic: also slower than numpy from about two
      thousand points up.

      Here a scalar stays a scalar. The engine already broadcasts one -- that
      is what its `is_vector` flag is for -- so this simply stops lying to it
      about the shape of the data. */
  void compute_curve(const std::vector<std::string>& parameter_names,
                     double* in_parameters, int n_parameters,
                     const std::string& axis_name,
                     double* in_axis, int n_axis,
                     double** out_values, int* n_out_values) const;

  //! Fix which parameter slot each variable reads, once, for a whole fit.
  /** \param parameter_names the order \ref compute_curve_bound will supply
      \param axis_name the name the equation gives the curve axis

      A fit re-evaluates the same equation with the same variables thousands
      of times, changing only their values. `compute_curve()` takes the names
      on every call, which means SWIG builds a `std::vector<std::string>` --
      allocating and copying every name -- per iteration, and the engine then
      looks each one up again. Measured at about 0.1 us per name, so a
      four-parameter model spends ~13% of a 512-point evaluation re-passing
      strings that cannot have changed.

      This resolves the binding once. Call it after `set_expression()` and
      whenever the parameter order changes; then evaluate with
      \ref compute_curve_bound, across which no string passes at all. */
  void bind_parameters(const std::vector<std::string>& parameter_names,
                       const std::string& axis_name);

  //! Evaluate a curve against the binding \ref bind_parameters established.
  /** \param in_parameters one value per name given to \ref bind_parameters,
             in that order
      \param n_parameters how many, checked against the binding
      \param in_axis the curve axis
      \param n_axis its length, and the length of the result
      \param out_values the curve
      \param n_out_values its length

      Throws `std::domain_error` if no binding has been established, or if
      \p n_parameters disagrees with it -- a silently mismatched binding
      would evaluate the right equation against the wrong values. */
  void compute_curve_bound(double* in_parameters, int n_parameters,
                           double* in_axis, int n_axis,
                           double** out_values, int* n_out_values) const;

  //! Whether \ref bind_parameters has been called for the current equation.
  bool has_parameter_binding() const;

  //! Whether a compiled plan is currently cached.
  /** Evaluation compiles on the first call and on a change of batch shape;
      during a fit neither happens, so the equation is parsed once and the
      sampler only ever re-binds values. */
  bool has_compiled_plan() const;

  //! Number of times the equation has been handed to the parser.
  /** A fit that re-parsed per step would show this climbing; it is pinned
      by a test precisely so that cannot regress unnoticed. */
  unsigned int get_number_of_compilations() const;

  //! Read the variables from the input ports, write the curve out.
  void evaluate() override;

  //! Whether an expression string can be compiled at all.
  static bool is_supported(const std::string& expression);

  std::string describe() const;

 private:
  //! Holds the ptolib evaluator, so its header stays in the implementation.
  /*! It must stay out of this header for two reasons: SWIG parses this file
      and would try to parse that one too, and every translation unit that
      includes an IMP header should not acquire ptolib's. */
  struct Impl;

  std::string expression_;   //!< as the caller wrote it, in Python syntax
  std::string translated_;   //!< the same equation in the engine's spelling
  std::vector<std::string> variables_;
  std::shared_ptr<Impl> impl_;

  void compile(const std::string& expression);
};

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_GRAPHEXPRESSION_H
