/**
 * \file Expression.cpp
 * \brief A model equation as a node, evaluated by ptolib's engine.
 *
 * This file used to be ~1900 lines: a tokeniser, a shunting-yard parser, an
 * RPN compiler with constant folding and common-subexpression elimination, a
 * block-vectorised evaluator with NEON kernels, and a vendored 1.6 MB copy of
 * ExprTk behind it as a fallback. All of that now lives in
 * `pto::ExpressionEngine` (ptolib), and what remains here is the part that is
 * genuinely imp.bff's: presenting it as a `Node` with ports, and the
 * numpy-facing entry points SWIG wraps.
 *
 * Why the engine moved rather than being shared some other way: `DataStore`
 * gates burst columns with the same language, tttrlib cannot depend on
 * imp.bff, and `AGENTS.md` places photons and curves in tttrlib. One
 * implementation had to sit on that side of the layering. See T-20260831-12.
 *
 * Copyright 2007-2023 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/Expression.h>

// The engine is ptolib's (https://github.com/tpeulen/ptolib), carried here as
// the verbatim copy `include/internal/ptolib.h` and compiled once in
// `src/Pto.cpp`. One implementation serves this node, tttrlib's DataStore
// gating and every other consumer; `test/test_vendored_headers.py` fails on
// drift from the sibling checkout. ptolib owns it; changes go there.
#include <IMP/bff/internal/ptolib.h>

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <new>
#include <stdexcept>
#include <string>

IMPBFF_BEGIN_NAMESPACE

using pto::ExprColumn;
using pto::ExprScalarType;

//! The compiled engine, kept out of the header.
struct Expression::Impl {
  pto::ExpressionEngine engine;
  //! Counts trips to the parser, so a fit that re-parsed per step is visible.
  unsigned int compilations = 0;

  //! Where each of the engine's variables reads its value from.
  /*! One entry per engine variable, in the engine's own order: \c kAxisSlot
      means "the curve axis", anything else is an index into the parameter
      array. Resolved once by bind_parameters() so a fit's inner loop passes
      no names and does no lookups. Empty until then. */
  static constexpr int kAxisSlot = -1;  // constexpr: ODR-used (push_back), so it needs a definition -- inline in C++17
  std::vector<int> binding;
  std::size_t bound_parameters = 0;

  //! Scratch for `evaluate()`, so a fit's inner loop allocates nothing.
  /*! `evaluate()` runs once per residual evaluation for the whole length of
      the curve. Rebuilding these per call was two heap allocations and a
      full-length copy on a path that runs tens of times per fit iteration. */
  std::vector<ExprColumn> columns;
  std::vector<double> curve;
};

namespace {

//! One column of doubles, or one value broadcast along the batch.
ExprColumn column_of(const double* data, bool is_vector) {
  ExprColumn c;
  c.data = data;
  c.type = ExprScalarType::Float64;
  c.is_vector = is_vector;
  return c;
}

//! Map the caller's names onto the engine's variables, by NAME.
/*! Not by position. `variables()` is in first-appearance order, which for
    `(g-b)/(r-b)` is `g, b, r` -- not the order a caller holds its columns
    in. Binding positionally silently evaluates the right expression over the
    wrong columns, which is a wrong answer rather than an error. */
template <typename T>
std::vector<T> bind_by_name(const std::vector<std::string>& engine_vars,
                            const std::vector<std::string>& names,
                            const std::vector<T>& supplied,
                            const char* where) {
  std::vector<T> bound;
  bound.reserve(engine_vars.size());
  for (const std::string& want : engine_vars) {
    std::vector<std::string>::const_iterator it =
        std::find(names.begin(), names.end(), want);
    if (it == names.end()) {
      throw std::domain_error(std::string(where) +
                              ": no column for variable '" + want + "'");
    }
    bound.push_back(supplied[static_cast<std::size_t>(it - names.begin())]);
  }
  return bound;
}

}  // namespace

Expression::Expression(const std::string& name)
    : Node(name), impl_(std::make_shared<Impl>()) {}

void Expression::compile(const std::string& expression) {
  // Compiled into a probe first, and only committed once it succeeds. A
  // failed `set_expression` must leave the node exactly as it was -- an
  // equation typed wrongly in a GUI should not destroy the fit that was
  // already running. Compiling in place would clear the program and leave
  // every later `compute()` throwing "no compiled program".
  pto::ExpressionEngine probe;
  if (!probe.compile(expression)) {
    // Refused, not approximated. There is no second evaluator behind this
    // one, deliberately: the one that used to be there answered
    // multi-argument functions wrongly rather than declining them.
    throw std::domain_error("Expression: cannot compile '" + expression +
                            "'");
  }
  Impl& impl = *impl_;
  impl.engine = probe;
  // A new equation has new variables, so any binding for the old one is
  // meaningless. Dropped rather than repaired: silently reusing it would
  // evaluate the new equation against the old equation's slots.
  impl.binding.clear();
  impl.bound_parameters = 0;
  // One, not one more. This counts trips to the parser *for the equation
  // currently held*, so a fit that re-parsed on every sampler move would
  // show it climbing; a caller that simply sets a new equation would not.
  impl.compilations = 1;
  expression_ = expression;
  translated_ = pto::ExpressionEngine::normalise(expression);
  variables_ = impl.engine.variables();
}

void Expression::set_expression(const std::string& expression) {
  compile(expression);
  set_valid(false);
}

bool Expression::is_supported(const std::string& expression) {
  // Never throws: callers use this to decide whether to offer an equation at
  // all, and an exception escaping a predicate is its own bug.
  try {
    pto::ExpressionEngine probe;
    return probe.compile(expression);
  } catch (...) {
    return false;
  }
}

bool Expression::has_compiled_plan() const { return impl_->engine.ready(); }

unsigned int Expression::get_number_of_compilations() const {
  return impl_->compilations;
}

std::vector<double> Expression::compute(
    const std::vector<std::string>& names,
    const std::vector<std::vector<double> >& values) const {
  if (names.size() != values.size()) {
    throw std::domain_error("Expression::compute: one value list per name");
  }
  std::vector<const std::vector<double>*> supplied;
  supplied.reserve(values.size());
  for (const std::vector<double>& v : values) supplied.push_back(&v);

  const std::vector<const std::vector<double>*> bound =
      bind_by_name(variables_, names, supplied, "Expression");

  // numpy broadcasting: the result is as long as the longest operand and
  // every scalar repeats along it. An empty column means an empty answer, as
  // in numpy -- checked before the length rule, or a zero-length operand
  // looks like a mismatch.
  for (const std::vector<double>* v : bound) {
    if (v->empty()) return std::vector<double>();
  }
  std::size_t n = 1;
  for (const std::vector<double>* v : bound) n = std::max(n, v->size());
  for (const std::vector<double>* v : bound) {
    if (v->size() != 1 && v->size() != n) {
      throw std::domain_error("Expression: operands of incompatible lengths");
    }
  }

  std::vector<ExprColumn> columns;
  columns.reserve(bound.size());
  for (const std::vector<double>* v : bound) {
    columns.push_back(column_of(v->data(), v->size() > 1));
  }
  std::vector<double> result(n);
  impl_->engine.compute_values(columns, n, result.data());
  return result;
}

std::vector<double> Expression::compute_batch(
    const std::vector<std::string>& names,
    const std::vector<std::vector<double> >& rows) const {
  std::vector<double> out;
  out.reserve(rows.size());
  std::vector<std::vector<double> > one(names.size());
  for (const std::vector<double>& row : rows) {
    if (row.size() != names.size()) {
      throw std::domain_error(
          "Expression::compute_batch: every row needs one value per name");
    }
    for (std::size_t i = 0; i < names.size(); ++i) one[i].assign(1, row[i]);
    const std::vector<double> r = compute(names, one);
    out.push_back(r.empty() ? 0.0 : r[0]);
  }
  return out;
}

//! Bind the rows of an `n_vars x n_rows` block to the engine's variables.
static std::vector<const double*> bind_block(
    const std::vector<std::string>& variables,
    const std::vector<std::string>& names, double* in_columns,
    std::size_t n_rows, const char* where) {
  std::vector<const double*> supplied;
  supplied.reserve(names.size());
  for (std::size_t i = 0; i < names.size(); ++i) {
    supplied.push_back(in_columns + i * n_rows);
  }
  return bind_by_name(variables, names, supplied, where);
}

void Expression::compute_columns(const std::vector<std::string>& names,
                                 double* in_columns, int n_vars, int n_rows,
                                 double** out_values,
                                 int* n_out_values) const {
  if (static_cast<int>(names.size()) != n_vars) {
    throw std::domain_error(
        "Expression::compute_columns: one name per supplied column");
  }
  if (n_rows < 0) {
    throw std::domain_error("Expression::compute_columns: negative length");
  }
  const std::size_t n = static_cast<std::size_t>(n_rows);
  const std::vector<const double*> bound =
      bind_block(variables_, names, in_columns, n, "Expression");

  double* out = static_cast<double*>(
      std::malloc(std::max<std::size_t>(n, 1) * sizeof(double)));
  if (out == nullptr) throw std::bad_alloc();

  std::vector<ExprColumn> columns;
  columns.reserve(bound.size());
  for (const double* p : bound) columns.push_back(column_of(p, n > 1));
  try {
    impl_->engine.compute_values(columns, n, out);
  } catch (...) {
    std::free(out);
    throw;
  }
  *out_values = out;
  *n_out_values = static_cast<int>(n);
}

void Expression::bind_parameters(
    const std::vector<std::string>& parameter_names,
    const std::string& axis_name) {
  Impl& impl = *impl_;
  std::vector<int> resolved;
  resolved.reserve(variables_.size());
  for (const std::string& want : variables_) {
    if (want == axis_name) {
      resolved.push_back(Impl::kAxisSlot);
      continue;
    }
    std::vector<std::string>::const_iterator it =
        std::find(parameter_names.begin(), parameter_names.end(), want);
    if (it == parameter_names.end()) {
      throw std::domain_error(
          "Expression::bind_parameters: no parameter or axis named '" + want +
          "'");
    }
    resolved.push_back(static_cast<int>(it - parameter_names.begin()));
  }
  impl.binding.swap(resolved);
  impl.bound_parameters = parameter_names.size();
}

bool Expression::has_parameter_binding() const {
  return !impl_->binding.empty() || variables_.empty();
}

void Expression::compute_curve_bound(double* in_parameters, int n_parameters,
                                     double* in_axis, int n_axis,
                                     double** out_values,
                                     int* n_out_values) const {
  const Impl& impl = *impl_;
  if (impl.binding.size() != variables_.size()) {
    throw std::domain_error(
        "Expression::compute_curve_bound: bind_parameters() has not been "
        "called for this equation");
  }
  if (static_cast<std::size_t>(n_parameters) != impl.bound_parameters) {
    // Refused rather than truncated: a binding that no longer matches the
    // caller's array would read the right equation from the wrong slots.
    throw std::domain_error(
        "Expression::compute_curve_bound: bound to " +
        std::to_string(impl.bound_parameters) + " parameters but given " +
        std::to_string(n_parameters));
  }
  if (n_axis < 0) {
    throw std::domain_error(
        "Expression::compute_curve_bound: negative axis length");
  }
  const std::size_t n = static_cast<std::size_t>(n_axis);

  std::vector<ExprColumn> columns;
  columns.reserve(impl.binding.size());
  for (int slot : impl.binding) {
    if (slot == Impl::kAxisSlot) {
      columns.push_back(column_of(in_axis, true));
    } else {
      columns.push_back(column_of(in_parameters + slot, false));
    }
  }

  double* out = static_cast<double*>(
      std::malloc(std::max<std::size_t>(n, 1) * sizeof(double)));
  if (out == nullptr) throw std::bad_alloc();
  try {
    impl.engine.compute_values(columns, n, out);
  } catch (...) {
    std::free(out);
    throw;
  }
  *out_values = out;
  *n_out_values = static_cast<int>(n);
}

void Expression::compute_curve(
    const std::vector<std::string>& parameter_names, double* in_parameters,
    int n_parameters, const std::string& axis_name, double* in_axis,
    int n_axis, double** out_values, int* n_out_values) const {
  if (static_cast<int>(parameter_names.size()) != n_parameters) {
    throw std::domain_error(
        "Expression::compute_curve: one name per parameter value");
  }
  if (n_axis < 0) {
    throw std::domain_error("Expression::compute_curve: negative axis length");
  }
  const std::size_t n = static_cast<std::size_t>(n_axis);

  // One ExprColumn per variable: the axis is a vector, every parameter is a
  // single value the engine broadcasts. Nothing is widened and nothing is
  // copied -- the engine reads the caller's numpy buffers where they lie.
  std::vector<ExprColumn> columns;
  columns.reserve(variables_.size());
  for (const std::string& want : variables_) {
    if (want == axis_name) {
      columns.push_back(column_of(in_axis, true));
      continue;
    }
    std::vector<std::string>::const_iterator it =
        std::find(parameter_names.begin(), parameter_names.end(), want);
    if (it == parameter_names.end()) {
      throw std::domain_error(
          "Expression::compute_curve: no parameter or axis named '" + want +
          "'");
    }
    columns.push_back(column_of(
        in_parameters + (it - parameter_names.begin()), false));
  }

  double* out = static_cast<double*>(
      std::malloc(std::max<std::size_t>(n, 1) * sizeof(double)));
  if (out == nullptr) throw std::bad_alloc();
  try {
    impl_->engine.compute_values(columns, n, out);
  } catch (...) {
    std::free(out);
    throw;
  }
  *out_values = out;
  *n_out_values = static_cast<int>(n);
}

void Expression::compute_mask(const std::vector<std::string>& names,
                              double* in_columns, int n_vars, int n_rows,
                              unsigned char** out_mask,
                              int* n_out_mask) const {
  if (static_cast<int>(names.size()) != n_vars || n_rows < 0) {
    throw std::domain_error("Expression::compute_mask: one name per column");
  }
  const std::size_t n = static_cast<std::size_t>(n_rows);
  const std::vector<const double*> bound =
      bind_block(variables_, names, in_columns, n, "Expression");

  unsigned char* out = static_cast<unsigned char*>(
      std::malloc(std::max<std::size_t>(n, 1)));
  if (out == nullptr) throw std::bad_alloc();

  std::vector<ExprColumn> columns;
  columns.reserve(bound.size());
  for (const double* p : bound) columns.push_back(column_of(p, n > 1));

  // The engine answers a gate in bits, which is the right unit for a
  // selection over millions of rows. This entry point predates that and
  // publishes one byte per row through numpy, so the words are unpacked
  // here rather than the engine being asked for a weaker form.
  try {
    std::vector<std::uint64_t> words((n + 63) / 64, 0);
    impl_->engine.compute_mask(columns, n, words.empty() ? nullptr
                                                         : words.data());
    for (std::size_t i = 0; i < n; ++i) {
      out[i] = (words[i >> 6] >> (i & 63)) & 1ULL ? 1 : 0;
    }
  } catch (...) {
    std::free(out);
    throw;
  }
  *out_mask = out;
  *n_out_mask = static_cast<int>(n);
}

void Expression::evaluate() {
  // The inner loop of a fit that lives entirely in C++, so it reads the
  // ports **where they lie**. The obvious spelling --
  // `values.push_back(p->get_value_vector())` and then `compute(names,
  // values)` -- copied every operand per call, and one of the operands is
  // the curve axis: several kilobytes of memcpy per residual evaluation,
  // for a buffer nothing here writes to. It also re-resolved every variable
  // name against the engine on each call, which `bind_parameters()` exists
  // to make unnecessary.
  Impl& impl = *impl_;
  impl.columns.clear();
  impl.columns.reserve(variables_.size());
  std::size_t n = 1;
  for (const std::string& name : variables_) {
    const std::shared_ptr<Port> p = get_input_port(name);
    if (!p) {
      throw std::domain_error("Expression '" + get_name() +
                              "': no input port for variable '" + name + "'");
    }
    const std::vector<double>& v = p->get_values_ref();
    // A scalar broadcasts; the batch length is the longest operand, which is
    // numpy's rule and the one the equations were written against.
    if (v.size() > n) n = v.size();
    impl.columns.push_back(column_of(v.data(), v.size() > 1));
  }
  impl.curve.resize(n);
  impl.engine.compute_values(impl.columns, n, impl.curve.data());

  // The curve goes to the output port keyed by the node's *own* name, which
  // is the graph's convention here rather than a fixed "value".
  const std::shared_ptr<Port> out = get_output_port(get_name());
  if (!out) {
    throw std::domain_error(
        "Expression '" + get_name() +
        "' writes its curve to the output port keyed by its own name, "
        "which this node does not have");
  }
  // A NaN curve must reach the misfit as a NaN: floored to `tiny` it would
  // read as a model that is zero everywhere, which is a plausible fit.
  out->set_sanitize(false);
  out->set_value_vector(impl.curve);
  set_valid(true);
}

std::string Expression::describe() const {
  std::string s = "Expression '" + get_name() + "'\n";
  s += "  equation:  " + expression_ + "\n";
  s += "  evaluated: pto::ExpressionEngine\n";
  s += "  variables: ";
  for (std::size_t i = 0; i < variables_.size(); ++i) {
    s += (i ? ", " : "") + variables_[i];
  }
  return s + "\n";
}

IMPBFF_END_NAMESPACE
