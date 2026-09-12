/**
 *  \file IMP/bff/InferenceFactorGraph.h
 *  \brief The factor structure of a model's posterior, as an explicit graph.
 *
 *  A posterior factorises exactly:
 *
 *      p(theta | D)  proportional to  prod_k L_k(theta_Sk) * prod_i pi_i(theta_i)
 *
 *  where Sk is the set of free variables dataset k's model actually reads and
 *  pi_i is variable i's prior. This class materialises that factorisation:
 *  variables are free parameters, factors are per-dataset likelihoods and
 *  per-variable priors.
 *
 *  What it answers:
 *
 *  - Relevance: which factors (and which local fits) must be re-evaluated
 *    when a given set of variables changes. This is the query that turns a
 *    global objective from O(N datasets) into O(what moved).
 *  - Structure: moralising the graph, eliminating variables in a greedy
 *    min-fill order and building a junction tree exposes the cliques, the
 *    separators and the treewidth. Cliques are the blocks a sampler or a scan
 *    should move jointly; the separator is the set of variables the datasets
 *    actually share.
 *  - Identifiability: treewidth and the connected components are results,
 *    not implementation details -- "these datasets are conditionally
 *    independent given {R0}" is the statement a global fit exists to make.
 *
 *  The design follows the architecture of mature probabilistic
 *  graphical-model toolkits (a model object separate from any inference
 *  engine, moralisation plus triangulation to expose blocks, relevance
 *  pruning per query). Their discrete sum-product kernels do not transfer to
 *  a continuous posterior; the structural machinery does, and needs nothing
 *  beyond the C++ standard library.
 *
 *  Ported from ChiSurf's chisurf.core.fitting.factorgraph (PRD-68). This
 *  layer is deliberately standalone: no IMP particles, restraints or
 *  decorators participate in the graph. Structural state is expected to
 *  enter as plain values at the model boundary.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_INFERENCEFACTORGRAPH_H
#define IMPBFF_INFERENCEFACTORGRAPH_H

#include <IMP/bff/bff_config.h>

#include <ostream>
#include <set>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! What a factor is: a dataset likelihood or a per-variable prior.
enum InferenceFactorKind {
  INFERENCE_FACTOR_PRIOR = 0,
  INFERENCE_FACTOR_LIKELIHOOD = 1,
  //! A factor over a hyperparameter and the variables it governs.
  /*!
      Not a prior on one constant: a roughness penalty couples every
      coefficient of the curve it smooths *and* the weight that scales it, and
      integrating that weight out rather than optimising it is what makes such
      an analysis work. Declared as a INFERENCE_FACTOR_PRIOR it is indistinguishable from a
      Gaussian on a scalar, and the elimination order that follows is not the
      one the analysis uses. Scope semantics are the others': it couples
      everything it names.
  */
  INFERENCE_FACTOR_HYPER = 2
};

//! One edge of a junction (clique) tree, joining two maximal cliques.
/*!
    The edge carries the separator: the variables the two cliques share.
    Conditioning on a separator makes the two sides independent.
*/
struct InferenceJunctionTreeEdge {
  //! Index of the first clique (into the clique list the tree is over).
  int first;
  //! Index of the second clique.
  int second;
  //! Sorted variable keys the two cliques share.
  std::vector<std::string> separator;
};

//! Stream a junction-tree edge (clique pair and its separator).
inline std::ostream& operator<<(std::ostream& out,
                                const InferenceJunctionTreeEdge& e) {
  out << "(" << e.first << ", " << e.second << ")[";
  for (std::size_t i = 0; i < e.separator.size(); ++i) {
    if (i) out << ", ";
    out << e.separator[i];
  }
  out << "]";
  return out;
}

//! The variables and factors of a model's posterior, with graph queries.
/*!
    Build by calling add_variable() for every free parameter and add_factor()
    for every likelihood and prior. Variables carry the position of the
    parameter in the model's flat free-parameter vector so the graph and a
    parameter array stay aligned; a fit_index groups likelihood factors by
    the local dataset they belong to (-1 for globals, priors, or an
    ungrouped model).

    Derived structure (the moral graph, elimination orders and the cliques)
    is cached for the lifetime of an unmutated graph; adding a node or a
    factor drops the caches. All queries are deterministic: ties break on
    the variable's flat-vector index.
*/
class IMPBFFEXPORT InferenceFactorGraph {
 public:
  InferenceFactorGraph();

  //! Add a free parameter (a variable of the posterior).
  /*!
      \param[in] key stable identifier of the variable
      \param[in] name human-readable parameter name
      \param[in] index position in the model's flat free-parameter vector
      \param[in] fit_index local fit the parameter belongs to, or -1 for a
                           global parameter
      \param[in] size how many free numbers the variable holds -- 1 for a
                 scalar, 24 for a curve on a 24-dimensional basis. It is the
                 *free* dimension that matters, not the constrained one: a
                 25-coefficient spline constrained to sum to zero has 24, and
                 24 is what enters #block_cost, elimination and treewidth.
      \param[in] role what kind of quantity it is -- the caller's vocabulary,
                 stored and returned unchanged. A consumer that keeps a
                 parallel dictionary of roles is keeping it because this
                 argument did not exist.
  */
  void add_variable(const std::string& key, const std::string& name,
                    int index, int fit_index = -1, int size = 1,
                    const std::string& role = "");

  //! Add a factor (a dataset likelihood or a per-variable prior).
  /*!
      Every variable named in scope must already exist. A factor couples all
      the variables it reads: it contributes a clique over its scope to the
      moral graph.

      \param[in] key stable identifier of the factor
      \param[in] kind INFERENCE_FACTOR_PRIOR or INFERENCE_FACTOR_LIKELIHOOD
      \param[in] scope keys of the variables this factor depends on
      \param[in] fit_index local fit of a likelihood factor, or -1
      \param[in] size number of residuals the factor contributes
  */
  void add_factor(const std::string& key, InferenceFactorKind kind,
                  const std::vector<std::string>& scope,
                  int fit_index = -1, int size = 0);

  //! Drop every cached derived structure (moral graph, orders, cliques).
  void invalidate();

  //! Number of variables (free parameters).
  unsigned int get_number_of_variables() const;
  //! Number of factors, priors and likelihoods together.
  unsigned int get_number_of_factors() const;
  //! Number of likelihood factors, ordered by local-fit index.
  unsigned int get_number_of_likelihood_factors() const;

  //! Variable keys ordered by their flat-parameter-vector index.
  std::vector<std::string> get_variable_keys() const;
  //! Factor keys in the order they were added.
  std::vector<std::string> get_factor_keys() const;
  //! The work a block's movement dirties: the sizes of the factors it touches.
  /*!
      #block_cost answers "how many local fits does this force", which is what
      it has always answered and what a scheduler over datasets wants. It is
      deliberately silent about a variable that forces no fit at all -- a
      hyperparameter reaching the data only through its own penalty costs
      zero there, correctly and unhelpfully.

      This is the other question, kept separate rather than folded in: the
      residuals that have to be recomputed, summed over every factor the block
      touches, priors and hyper factors included. For a hyperparameter over a
      24-coefficient curve it is 24 + 1 rather than 0 -- it forces no fit, and
      it dirties the roughness factor.
  */
  int factor_cost(const std::vector<std::string>& block) const;

  //! Free numbers a variable holds; 1 unless it was given a size, 0 if absent.
  int get_variable_size(const std::string& key) const;
  //! A variable's role as the caller declared it; empty if none or absent.
  std::string get_variable_role(const std::string& key) const;
  //! A factor's kind; INFERENCE_FACTOR_PRIOR for an unknown key, so check with factor_keys().
  InferenceFactorKind get_factor_kind(const std::string& factor_key) const;
  //! Factors of one kind.
  unsigned int get_number_of_factors_of_kind(InferenceFactorKind kind) const;

  //! The graph as a JSON document: every variable and factor, nothing derived.
  /*!
      Round trips exactly through #from_json. Structure only -- the moral
      graph, the orders and the cliques are recomputed on demand and are not
      written, because they are answers rather than state.
  */
  std::string to_json() const;
  //! Replace this graph with one read from #to_json's output.
  void from_json(const std::string& json);
  //! #to_json to a file.
  void save(const std::string& path) const;
  //! #from_json from a file.
  void load(const std::string& path);

  //! Flat-parameter-vector position of a variable key, or -1 if absent.
  int index_of(const std::string& key) const;
  //! Variable key at a position of the flat parameter vector, "" if absent.
  std::string key_at(int index) const;
  //! Keys of the factors a variable appears in.
  std::vector<std::string> factors_of(const std::string& key) const;
  //! Scope (variable keys) of a factor, empty if the factor is unknown.
  std::vector<std::string> variables_of(const std::string& factor_key) const;

  //! The independent sub-problems: components of the moral graph.
  /*!
      Variables in different components share no factor, so their posteriors
      are independent and can be optimised, scanned or sampled separately.
      Returned largest first.
  */
  std::vector<std::vector<std::string> > connected_components() const;

  //! A greedy variable-elimination order.
  /*!
      \param[in] heuristic "min_fill" (default) or "min_degree". min_fill
             repeatedly eliminates the variable whose elimination adds the
             fewest new edges; min_degree the one with the fewest neighbours.
             Ties break on the flat-vector index, so the order is
             deterministic.

      \return variable keys in elimination order
  */
  std::vector<std::string> get_elimination_order(
      const std::string& heuristic = "min_fill") const;

  //! The maximal cliques induced by the default elimination order.
  /*!
      Each eliminated variable together with its then-remaining neighbours
      forms a clique of the triangulated graph; non-maximal cliques are
      dropped. Cliques are the blocks a sampler or a scan should move
      jointly. Largest first.
  */
  std::vector<std::vector<std::string> > get_cliques() const;

  //! max clique size - 1: the model's structural difficulty.
  /*!
      The cost of exact marginalisation is exponential in this number, and
      it is the dimension a blocked sampler would have to move jointly. A
      star-shaped global fit (many datasets, a few shared globals) has a
      small treewidth however many datasets it holds; a single dataset whose
      likelihood couples all its parameters has n_free - 1.
  */
  int get_treewidth() const;

  //! The junction (clique) tree edges over get_cliques().
  /*!
      Each edge carries the separator of the two cliques it joins. The tree
      is the maximum-weight spanning tree of the clique graph weighted by
      separator size, the standard construction guaranteeing the
      running-intersection property. A disconnected fit gives a forest.
  */
  std::vector<InferenceJunctionTreeEdge> get_junction_tree_edges() const;

  //! The distinct separators of the junction tree, longest first.
  /*!
      The separator of a global fit is the set of parameters its datasets
      genuinely share: condition on it and the datasets become independent.
  */
  std::vector<std::vector<std::string> > get_separators() const;

  //! A partition of the variables for block-wise sampling.
  /*!
      Two variables land in the same block exactly when the same set of
      datasets depends on both. Unlike get_cliques() this is a true
      partition, which is what a block sampler needs. Cheapest blocks
      (fewest datasets) come first so a sweep front-loads the inexpensive
      moves; variables no likelihood touches come last.
  */
  std::vector<std::vector<std::string> > get_sampling_blocks() const;

  //! How many local fits a move of block must recompute.
  int block_cost(const std::vector<std::string>& block) const;

  //! The factors that must be re-evaluated for a set of changed variables.
  std::vector<std::string> affected_factors(
      const std::vector<std::string>& changed) const;

  //! The local fits that must be recomputed for a set of changed variables.
  /*!
      This is the query that makes a global objective proportional to what
      actually moved rather than to the number of datasets. Returns sorted
      fit indices.
  */
  std::vector<int> affected_fits(
      const std::vector<std::string>& changed) const;

  //! The variables no likelihood factor depends on.
  /*!
      Such a variable is free but, as far as the graph can tell, reaches no
      data. A caller using affected_fits() to skip work must treat a change
      to one of these conservatively and recompute everything.
  */
  std::vector<std::string> get_unexplained_variables() const;

  //! A short human-readable structure report.
  /*!
      Renders the counts, the treewidth, the separators and the independent
      components -- the identifiability statement a global fit exists to
      make.
  */
  std::string describe() const;

 private:
  struct Variable {
    std::string key, name;
    int index;
    int fit_index;
    int size = 1;
    std::string role;
  };
  struct Factor {
    std::string key;
    InferenceFactorKind kind;
    std::vector<int> scope;  // variable positions
    int fit_index;
    int size;
  };

  //! Adjacency over variable positions; builds and caches the moral graph.
  const std::vector<std::set<int> >& moral_adjacency() const;
  bool is_complete() const;
  std::vector<std::vector<std::string> > compute_cliques(
      const std::vector<int>& order) const;

  std::vector<Variable> variables_;
  std::vector<Factor> factors_;
  std::map<std::string, int> variable_index_of_;  // key -> position
  std::map<std::string, int> factor_index_of_;    // key -> position
  std::vector<std::vector<int> > incidence_;      // var pos -> factor keys

  mutable bool has_moral_ = false;
  mutable std::vector<std::set<int> > moral_;
  mutable std::vector<int> order_min_fill_, order_min_degree_;
  mutable bool has_order_min_fill_ = false, has_order_min_degree_ = false;
  mutable std::vector<std::vector<std::string> > cliques_;
  mutable bool has_cliques_ = false;
};

IMPBFF_END_NAMESPACE

#endif // IMPBFF_INFERENCEFACTORGRAPH_H
