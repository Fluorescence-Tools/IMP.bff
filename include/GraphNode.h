/**
 *  \file IMP/bff/GraphNode.h
 *  \brief A lazily-evaluated function over named ports: chinet's GraphNode.
 *
 *  A GraphNode owns named ports (inputs and outputs) and computes its outputs
 *  from its inputs on demand. Evaluation is lazy and invalidation-driven:
 *  writing an input port invalidates the node, a reactive input evaluates
 *  it immediately, and evaluate() of an upstream node pushes its result
 *  through the links into the followers, invalidating (or, when they are
 *  reactive, re-evaluating) the nodes downstream. A node is valid after it
 *  has been evaluated, invalid from construction, and is_valid() is
 *  additionally false while any node an input is linked to is invalid.
 *
 *  Three callback mechanisms, mirroring chinet:
 *
 *  - string operators (set_callback("multiply_double", "C")), the four
 *    arithmetic operators of chinet's C-callback protocol;
 *  - a std::function taking the input and output port maps
 *    (set_callback_function), the C++ spelling of chinet's callback_class;
 *  - a virtual evaluate(), which a C++ subclass can override outright --
 *    and, since GraphNode is a SWIG director, a Python subclass too: chinet's
 *    set_python_callback_function (signature introspection, ports from
 *    parameters, outputs from a dict return) lives in the Python layers
 *    on top (chisurf/core/nodes.py), which override evaluate() to call
 *    the Python callable over the port maps.
 *
 *  Ported from chinet's chinet/node.py (phase 1 of removing chinet from
 *  chisurf). Standalone, like FactorGraph: no IMP particles, restraints
 *  or decorators. Ports are shared_ptr-owned by the node; a GraphNode must be
 *  shared_ptr-owned too (its ports hold it weakly), which every
 *  Python-wrapped node is.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_GRAPHNODE_H
#define IMPBFF_GRAPHNODE_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/GraphPort.h>

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! A function over named ports, evaluated lazily when its inputs settle.
class IMPBFFEXPORT GraphNode : public GraphObject,
                          public std::enable_shared_from_this<GraphNode> {
 public:
  //! Ports by name; the map every port lookup answers from.
  typedef std::map<std::string, std::shared_ptr<GraphPort> > GraphPortMap;
  //! A C++ callback over the port maps (chinet's callback_class.run).
  /*!
      \param[in] inputs the node's input ports
      \param[out] outputs the node's output ports, to be written
  */
  typedef std::function<void(const GraphPortMap& inputs, GraphPortMap& outputs)>
      Callback;

  //! An empty node, invalid until evaluated.
  explicit GraphNode(const std::string& name = "");
  virtual ~GraphNode();

  //! A node with ports (chinet's GraphNode(ports=...)): port names follow keys.
  /*!
      A factory rather than a constructor because ports attach to their
      node through a weak_ptr, and a node can only hand those out once a
      shared_ptr owns it -- which is not yet true inside a constructor.
  */
  static std::shared_ptr<GraphNode> make_graph_node(const std::string& name,
                                         const GraphPortMap& ports);

  //! Adopt a port map (chinet's set_ports): each port keeps its own
  //! is_output flag and takes its dict key as its name.
  void set_ports(const GraphPortMap& ports);
  //! Add a port under a key, as input or output; replaces an existing key.
  void add_port(const std::string& key, std::shared_ptr<GraphPort> port,
                bool is_output);
  //! Add an input port.
  void add_input_port(const std::string& key, std::shared_ptr<GraphPort> port);
  //! Add an output port.
  void add_output_port(const std::string& key, std::shared_ptr<GraphPort> port);

  //! A port by name, or a null pointer.
  std::shared_ptr<GraphPort> get_port(const std::string& name) const;
  //! An input port by name, or a null pointer.
  std::shared_ptr<GraphPort> get_input_port(const std::string& name) const;
  //! An output port by name, or a null pointer.
  std::shared_ptr<GraphPort> get_output_port(const std::string& name) const;
  //! All ports.
  GraphPortMap get_ports() const;
  //! The input ports.
  GraphPortMap get_input_ports() const;
  //! The output ports.
  GraphPortMap get_output_ports() const;

  //! Set a string callback (chinet's set_callback).
  /*!
      \param[in] callback operator name: addition_double, addition_int,
                  multiply_double or multiply_int -- the result is written
                  to the output port keyed by the node's name
      \param[in] callback_type "C" (built-in operator) or
                  "CLASS"/"CPP"/"C++" (a callback object, set separately);
                  anything else disables the node
  */
  void set_callback(const std::string& callback,
                    const std::string& callback_type);
  //! Set a C++ callback (chinet's callback_class). Not SWIG-wrapped.
  void set_callback_function(Callback callback);
  //! The operator name last set.
  const std::string& get_callback() const;
  //! The callback kind last set: 0 operator, 1 object, -1 none.
  int get_callback_type() const;
  //! The callback-type tag verbatim ("C", "CLASS", ...), chinet's
  //! callback_type document field; set_callback stores it unmodified.
  const std::string& get_callback_type_string() const;
  //! The raw validity flag (chinet's node_valid_), which is what a saved
  //! document carries -- is_valid() additionally consults the inputs.
  bool get_node_valid() const;
  //! GraphPort keys in insertion order (chinet's ports dict order, which the
  //! map of get_ports() does not preserve).
  const std::vector<std::string>& get_port_order() const;

  //! How many times this node has evaluated, ever.
  /*!
      Monotonic, never reset. The invalidation contract -- that `update()`
      re-evaluates only what a change reached -- is only assertable by
      counting, and a consumer should not have to monkey-patch its own
      callbacks to do it. Counts evaluations driven by `update()` and by a
      reactive port; a direct call to `evaluate()` is the caller's own and is
      not counted, because the caller already knows it made it.
  */
  unsigned long long get_evaluation_count() const { return eval_count_; }

#ifndef SWIG
  //! About to evaluate. Called by GraphNode and GraphPort; not for callers.
  void note_evaluation() { ++eval_count_; }
#endif

  //! Skip `evaluate()` when every input holds what it held last time.
  /*!
      Off by default, and deliberately opt-in: it is only sound for a node
      that is a *function* of its inputs. An `GraphExpression` is one; a node that
      draws a random number, reads a clock, or accumulates across calls is
      not, and memoising it would freeze its output. The class cannot tell
      which it has -- a Python director is opaque -- so the caller says.

      **What it is for.** A finite-difference Jacobian evaluates the graph
      once per free parameter, and a joint fit's members are mostly untouched
      by any one of them: perturbing dataset 1's amplitude leaves dataset 2's
      model exactly where it was, yet `update()` recomputes it, because
      invalidation travels by *reachability* and not by whether anything
      moved. With this on, the untouched member compares its inputs, finds
      them identical and keeps its outputs. The saving grows with the number
      of members and the number of free parameters, which is exactly the
      direction a joint fit gets expensive in.

      **What it costs when it does not hit.** One comparison of the input
      values, which stops at the first difference, plus a copy of them after
      an evaluation that did happen. For scalar parameters that is nothing.
      For a node whose input is a long vector the comparison is O(n) against
      an evaluation that is at least O(n) and usually far more per element --
      a compare against a call to `exp` is not a close race -- but a node
      whose evaluation is cheaper than reading its own inputs should leave
      this off.

      Measured on `GraphExpression("a*exp(-x/t)") -> ChiSquared`, 3 000 updates,
      2026-09-08:

      | points  | value rewritten | every step new |
      |---------|-----------------|----------------|
      | 256     | 2.2x faster     | 4.9% slower    |
      | 4 096   | 5.8x faster     | 4.4% slower    |
      | 65 536  | 10.0x faster    | 5.6% slower    |

      The hit grows with the array because the comparison is the cheap half
      of what it replaces; the miss is a flat few per cent, which is what
      makes the default `off` a choice rather than a hedge.

      Turning it off discards the record, so a caller may switch it per phase
      (on for a Jacobian, off for a stochastic pass) and the next phase starts
      from a real evaluation. Keeping the record across an off phase would be
      wrong, not merely wasteful: the node goes on evaluating while off
      without updating it, so the record would describe inputs its outputs no
      longer came from.
   */
  void set_memoize(bool v);
  bool get_memoize() const { return memoize_; }

  //! How many times `update()` skipped `evaluate()` because nothing changed.
  /*! Monotonic, never reset, and zero while #set_memoize is off. The saving
      is only assertable by counting, the same argument as
      #get_evaluation_count. */
  unsigned long long get_memo_hit_count() const { return memo_hits_; }

  //! False while any node a linked input follows is invalid.
  bool inputs_valid() const;
  //! Valid from evaluation, false while inputs are unsettled or changed.
  bool is_valid() const;
  //! Force the validity flag (invalidates nothing; chinet semantics).
  /*! Clears any #set_memoize record: a caller invalidating a node by hand
      means something the inputs do not show, so the node must really run. */
  void set_valid(bool v);

#ifndef SWIG
  //! A port of this node was written. Called by GraphPort; not for callers.
  /*!
      Counts writes so that update() can tell a node that *computed* from one
      that has nothing to do. It cannot ask whether evaluate() was overridden
      -- a director subclass in Python is indistinguishable from C++ -- but it
      can ask whether anything came out, which is the question it actually
      means.
  */
  void note_port_write() { ++write_epoch_; }

  //! A port of this node was written: invalidate, but keep the memo record.
  /*!
      The one invalidation that #set_memoize survives, and the reason it is
      separate from `set_valid(false)`. A port write is exactly the change
      the memo is there to *filter* -- it says "something upstream moved",
      and comparing the inputs is how the node finds out whether it moved to
      somewhere new. Every other invalidation means something the inputs
      cannot show (a dataset swapped, an internal coefficient changed, a
      caller saying "do it again"), and those clear the record so the node
      really does evaluate. Called by GraphPort; not for callers.
   */
  void invalidate_from_port() {
    node_valid_ = false;
    ++write_epoch_;
  }
#endif

  //! Compute the outputs from the inputs; marks the node valid.
  /*!
      No-op when the node has neither a callback object nor an operator
      (chinet: callback_class None and callback_type < 0). After the
      callback runs, every other node sharing an output port is
      invalidated.
  */
  virtual void evaluate();
  //! Pull linked inputs to their sources' values, evaluating upstream first.
  /*!
      Recursively updates any upstream node that is invalid, copies the
      sources' values into this node's input ports, and evaluates itself
      if still invalid. The link graph is acyclic, so this terminates.
      Virtual so a director subclass can intercept the lazy pull, as
      evaluate() above is.
  */
  virtual void update();

  //! A one-line summary: name, ports, validity.
  std::string describe() const;

  //! Invalidate the cached execution plan.
  //!
  //! Every structural mutator calls this, and so does a GraphPort whose link
  //! changed: linking is what decides which inputs the plan has to pull.
  void touch_structure();

 private:
  //! Rebuild the in_/out_ lookups from ports_ in insertion order.
  void fill_input_output_port_lookups();

  //! The operator resolved from its name once, instead of per evaluation.
  enum Operator { OP_NONE = 0, OP_ADD, OP_MUL };

  //! Writes to this node's ports, so update() can see that evaluate() worked.
  unsigned long long write_epoch_ = 0;
  //! Evaluations, for #get_evaluation_count.
  unsigned long long eval_count_ = 0;

  //! Everything ``evaluate()`` and ``update()`` would otherwise re-derive
  //! from the string-keyed port maps on every call.
  //!
  //! Sampling calls ``update()`` millions of times over a graph whose shape
  //! never changes, so the map lookups, string comparisons and shared_ptr
  //! refcount traffic dominated the arithmetic. The plan is rebuilt lazily
  //! whenever ``structure_epoch_`` moves, which every structural mutator
  //! bumps. Raw pointers are safe here: the node owns its ports through
  //! ``ports_``, so a planned port outlives the plan.
  struct Plan {
    long epoch = -1;
    //! ``set_name`` is not virtual on GraphObject, so a rename cannot bump
    //! the epoch. The operator writes to the output port keyed by the
    //! node's name, so the plan re-resolves when the name moves.
    std::string name;
    std::vector<GraphPort*> linked_inputs;
    Operator op = OP_NONE;
    bool as_int = false;
    GraphPort* operand_a = nullptr;
    GraphPort* operand_b = nullptr;
    GraphPort* op_output = nullptr;
    std::vector<GraphPort*> outputs;
  };

  //! Bring plan_ up to date with the current structure if it is stale.
  void refresh_plan() const;

  mutable Plan plan_;
  long structure_epoch_ = 0;

  GraphPortMap ports_;
  GraphPortMap in_;
  GraphPortMap out_;
  //! Insertion order of ports_ (std::map sorts; chinet's dicts do not, and
  //! the operator callbacks take the first two inputs by insertion order).
  std::vector<std::string> port_order_;
  std::string callback_;
  std::string callback_type_string_;
  int callback_type_ = -1;
  Callback callback_class_;
  bool node_valid_ = false;
  //! #set_memoize, and the inputs as of the last evaluation under it.
  bool memoize_ = false;
  bool memo_valid_ = false;
  std::vector<double> memo_inputs_;
  unsigned long long memo_hits_ = 0;
#ifndef SWIG
  //! Do the inputs hold exactly what they held at the last evaluation?
  bool memo_matches() const;
  //! Record the inputs as of an evaluation that just happened.
  void memo_record();
#endif
};

IMPBFF_END_NAMESPACE

#endif // IMPBFF_GRAPHNODE_H
