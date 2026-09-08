/**
 *  \file IMP/bff/EvaluationGraph.h
 *  \brief A graph of nodes that is assembled, then run on demand.
 *
 *  `Node` already evaluates lazily: assembling computes nothing, `update()`
 *  pulls what a node needs, and a node that is still valid is not run again.
 *  What it does not have is a handle for the *graph* -- somewhere to say what
 *  the interesting results are called and to ask for them by name.
 *
 *  That is what this is. Outputs are registered under **labels**, each naming
 *  a (node, port) pair, and `run()` evaluates exactly what the labels ask for.
 *  A label rather than a node name on purpose: a settings file or a
 *  colleague's notebook names `p_R`, and a node should be renameable or
 *  re-plumbable without breaking either.
 *
 *  A label may name **any** port, not only a terminal one. Naming an
 *  intermediate is the normal mode of use rather than a corner case -- an
 *  instrument basis that several histograms share and that a weight scan
 *  should not touch; the expected counts, which one wants to compare against
 *  the data for two milliseconds' work rather than compute a posterior for
 *  twenty seconds; a mode that is a result for one purpose and an input for
 *  another. Because `update()` walks upstream and never down, naming an
 *  intermediate stops the run there by construction.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_EVALUATIONGRAPH_H
#define IMPBFF_EVALUATIONGRAPH_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/Node.h>
#include <IMP/bff/Port.h>

#include <map>
#include <memory>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! What a run did. Not what it produced -- the values stay on the ports.
struct IMPBFFEXPORT RunReport {
  //! Nodes upstream of what was asked for, this one included.
  int nodes_visited = 0;
  //! Of those, the ones that actually evaluated.
  int nodes_evaluated = 0;
  //! Wall-clock seconds the run took.
  double seconds = 0.0;
  std::string describe() const;
};

//! A graph of nodes with named outputs, evaluated on demand.
class IMPBFFEXPORT EvaluationGraph {
 public:
  EvaluationGraph();

  //! Name a port so a run can ask for it.
  /*!
      \param[in] label what callers and settings files know it by
      \param[in] node the node the port belongs to
      \param[in] port_name the port's name on that node
      \param[in] provenance an opaque string identifying the code behind the
                 node, stored and round-tripped and **never interpreted**.
                 `IMP.bff.node_provenance()` builds one from a Python node's
                 module, class and source hash; anything else the caller finds
                 more meaningful works equally well, because this class only
                 hands it back.
      \throws IMP::ValueException if the node has no such port, or the label
              is already taken -- silently rebinding a label is how a saved
              settings file comes to mean something else.
  */
  void add_output(const std::string& label, std::shared_ptr<Node> node,
                  const std::string& port_name,
                  const std::string& provenance = "");

  //! Forget a label. The node and its port are untouched.
  void remove_output(const std::string& label);

  //! Registered labels, in the order they were added.
  std::vector<std::string> get_output_labels() const;

  //! The port a label names, or null if there is no such label.
  std::shared_ptr<Port> get_output_port(const std::string& label) const;

  //! The node a label names, or null.
  std::shared_ptr<Node> get_output_node(const std::string& label) const;

  //! The provenance recorded with a label; empty if none or no such label.
  std::string get_output_provenance(const std::string& label) const;

  //! Evaluate what every registered label needs.
  RunReport run();

  //! Evaluate what these labels need, and no more.
  /*!
      Two labels sharing an upstream subgraph evaluate it once: the first
      leaves it valid and the second finds it so.

      \param[in] outputs labels to produce
      \throws IMP::ValueException if a label is not registered -- a
              misspelling that quietly produced nothing would be worse.
  */
  RunReport run(const std::vector<std::string>& outputs);

  //! The nodes a label depends on, itself included, upstream first.
  std::vector<std::string> get_dependencies(const std::string& label) const;

  //! Labels and the ports they name, as JSON.
  /*!
      **The nodes are not written**, and not for want of trying: a node's work
      is its callback, and a callback is code. Saving its *source* was
      measured and does not work. `inspect.getsource` fails outright for a
      class defined in a notebook cell or by `exec` ("is a built-in class"),
      which is where this kind of model is usually written; and where it does
      succeed, re-executing the text in a fresh namespace fails on the first
      name it does not carry -- the source is not the closure, and a callback
      closes over its imports, its module constants and whatever else it
      calls. Capturing all of that is what `cloudpickle` and `dill` do, at the
      price of breaking across library versions.

      There is a second reason not to, which would stand even if it worked:
      a document that carries code is a document that executes when loaded,
      and these are shared between colleagues and copied off servers.

      So what round-trips is the naming -- which label points at which port of
      which node, by name -- plus whatever #add_output was given as
      provenance. Rebuilding the nodes is the caller's, and #from_json
      re-attaches the labels to them.
  */
  std::string to_json() const;

  //! Re-attach labels to nodes the caller has rebuilt.
  /*!
      \param[in] json a #to_json document
      \param[in] nodes the rebuilt nodes, by the name they had when written
      \throws IMP::ValueException if the document names a node that is not
              among \p nodes, or a port that node does not have.
  */
  void from_json(const std::string& json,
                 const std::map<std::string, std::shared_ptr<Node> >& nodes);

  //! How many labels are registered.
  unsigned int get_number_of_outputs() const;

 private:
  struct Output {
    std::string label;
    std::shared_ptr<Node> node;
    std::string port_name;
    std::string provenance;
  };
  std::vector<Output> outputs_;
  std::map<std::string, int> index_of_;

  const Output* find(const std::string& label) const;
  RunReport run_nodes(const std::vector<std::shared_ptr<Node> >& roots);
};

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_EVALUATIONGRAPH_H
