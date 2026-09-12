/**
 *  \file InferenceFactorGraph.cpp
 *  \brief The factor structure of a model's posterior (see InferenceFactorGraph.h).
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/bff/InferenceFactorGraph.h>
#include <IMP/bff/IMPCompatibility.h>

#include <IMP/bff/internal/json.h>

#include <algorithm>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

IMPBFF_BEGIN_NAMESPACE

namespace {

//! Union-find over clique indices, for Kruskal's spanning tree.
struct DisjointSets {
  explicit DisjointSets(int n) : parent_(n) {
    for (int i = 0; i < n; ++i) parent_[i] = i;
  }
  int find(int i) {
    while (parent_[i] != i) {
      parent_[i] = parent_[parent_[i]];
      i = parent_[i];
    }
    return i;
  }
  void unite(int i, int j) { parent_[find(i)] = find(j); }
  std::vector<int> parent_;
};

}  // namespace

InferenceFactorGraph::InferenceFactorGraph() = default;

void InferenceFactorGraph::add_variable(const std::string& key, const std::string& name,
                               int index, int fit_index, int size,
                               const std::string& role) {
  if (variable_index_of_.count(key)) {
    throw std::invalid_argument("duplicate variable key: " + key);
  }
  if (size < 1) {
    IMP_THROW("variable '" << key << "' must hold at least one free number",
              IMP::ValueException);
  }
  Variable v;
  v.key = key;
  v.name = name;
  v.index = index;
  v.fit_index = fit_index;
  v.size = size;
  v.role = role;
  variables_.push_back(v);
  variable_index_of_[key] = static_cast<int>(variables_.size()) - 1;
  incidence_.push_back({});
  invalidate();
}

void InferenceFactorGraph::add_factor(const std::string& key, InferenceFactorKind kind,
                             const std::vector<std::string>& scope,
                             int fit_index, int size) {
  if (factor_index_of_.count(key)) {
    throw std::invalid_argument("duplicate factor key: " + key);
  }
  Factor f;
  f.key = key;
  f.kind = kind;
  f.fit_index = fit_index;
  f.size = size;
  for (const auto& k : scope) {
    auto it = variable_index_of_.find(k);
    if (it == variable_index_of_.end()) {
      throw std::invalid_argument("unknown variable in scope of " + key +
                                  ": " + k);
    }
    f.scope.push_back(it->second);
  }
  factors_.push_back(f);
  factor_index_of_[key] = static_cast<int>(factors_.size()) - 1;
  const int fpos = static_cast<int>(factors_.size()) - 1;
  for (int v : f.scope) incidence_[v].push_back(fpos);
  invalidate();
}

void InferenceFactorGraph::invalidate() {
  has_moral_ = false;
  moral_.clear();
  has_order_min_fill_ = false;
  has_order_min_degree_ = false;
  order_min_fill_.clear();
  order_min_degree_.clear();
  has_cliques_ = false;
  cliques_.clear();
}

unsigned int InferenceFactorGraph::get_number_of_variables() const {
  return static_cast<unsigned int>(variables_.size());
}

unsigned int InferenceFactorGraph::get_number_of_factors() const {
  return static_cast<unsigned int>(factors_.size());
}

unsigned int InferenceFactorGraph::get_number_of_likelihood_factors() const {
  unsigned int n = 0;
  for (const auto& f : factors_) {
    if (f.kind == INFERENCE_FACTOR_LIKELIHOOD) ++n;
  }
  return n;
}

std::vector<std::string> InferenceFactorGraph::get_variable_keys() const {
  std::vector<const Variable*> vs;
  vs.reserve(variables_.size());
  for (const auto& v : variables_) vs.push_back(&v);
  std::stable_sort(vs.begin(), vs.end(),
                   [](const Variable* a, const Variable* b) {
                     return a->index < b->index;
                   });
  std::vector<std::string> out;
  out.reserve(vs.size());
  for (const auto* v : vs) out.push_back(v->key);
  return out;
}

std::vector<std::string> InferenceFactorGraph::get_factor_keys() const {
  std::vector<std::string> out;
  out.reserve(factors_.size());
  for (const auto& f : factors_) out.push_back(f.key);
  return out;
}

int InferenceFactorGraph::index_of(const std::string& key) const {
  auto it = variable_index_of_.find(key);
  if (it == variable_index_of_.end()) return -1;
  return variables_[it->second].index;
}

std::string InferenceFactorGraph::key_at(int index) const {
  for (const auto& v : variables_) {
    if (v.index == index) return v.key;
  }
  return "";
}

std::vector<std::string> InferenceFactorGraph::factors_of(
    const std::string& key) const {
  auto it = variable_index_of_.find(key);
  if (it == variable_index_of_.end()) return {};
  std::vector<std::string> out;
  for (int f : incidence_[it->second]) out.push_back(factors_[f].key);
  return out;
}

std::vector<std::string> InferenceFactorGraph::variables_of(
    const std::string& factor_key) const {
  auto it = factor_index_of_.find(factor_key);
  if (it == factor_index_of_.end()) return {};
  std::vector<std::string> out;
  for (int v : factors_[it->second].scope) out.push_back(variables_[v].key);
  return out;
}

const std::vector<std::set<int> >& InferenceFactorGraph::moral_adjacency() const {
  if (has_moral_) return moral_;
  moral_.assign(variables_.size(), {});
  for (const auto& f : factors_) {
    std::vector<int> scope = f.scope;
    std::sort(scope.begin(), scope.end());
    for (std::size_t a = 0; a < scope.size(); ++a) {
      for (std::size_t b = a + 1; b < scope.size(); ++b) {
        moral_[scope[a]].insert(scope[b]);
        moral_[scope[b]].insert(scope[a]);
      }
    }
  }
  has_moral_ = true;
  return moral_;
}

bool InferenceFactorGraph::is_complete() const {
  const auto& adj = moral_adjacency();
  const std::size_t n = adj.size();
  // The adjacency is symmetric, so this sum is *twice* the edge count. It
  // used to be compared against n(n-1)/2 -- half of what it should be --
  // which called a graph complete whenever |E| = n(n-1)/4. A four-variable
  // star has three edges and hits that exactly, and a star is the shape of
  // every global fit: several datasets around one shared parameter. It then
  // took the complete-graph shortcut below, reporting one clique over
  // everything and a treewidth of n-1 for a graph whose treewidth is 1.
  std::size_t directed = 0;
  for (const auto& s : adj) directed += s.size();
  return directed == n * (n - 1);
}

std::vector<std::vector<std::string> > InferenceFactorGraph::connected_components()
    const {
  const auto& adj = moral_adjacency();
  const int n = static_cast<int>(adj.size());
  std::vector<int> comp(n, -1);
  std::vector<std::vector<int> > components;
  for (int start = 0; start < n; ++start) {
    if (comp[start] != -1) continue;
    std::vector<int> stack{start};
    comp[start] = static_cast<int>(components.size());
    std::vector<int> members{start};
    while (!stack.empty()) {
      int v = stack.back();
      stack.pop_back();
      for (int u : adj[v]) {
        if (comp[u] == -1) {
          comp[u] = comp[start];
          members.push_back(u);
          stack.push_back(u);
        }
      }
    }
    components.push_back(members);
  }
  std::stable_sort(components.begin(), components.end(),
                   [](const std::vector<int>& a, const std::vector<int>& b) {
                     return a.size() > b.size();
                   });
  std::vector<std::vector<std::string> > out;
  out.reserve(components.size());
  for (const auto& members : components) {
    std::vector<int> sorted = members;
    std::sort(sorted.begin(), sorted.end(), [&](int a, int b) {
      return variables_[a].index < variables_[b].index;
    });
    std::vector<std::string> keys;
    keys.reserve(sorted.size());
    for (int v : sorted) keys.push_back(variables_[v].key);
    out.push_back(keys);
  }
  return out;
}

std::vector<std::string> InferenceFactorGraph::get_elimination_order(
    const std::string& heuristic) const {
  bool use_min_fill;
  if (heuristic == "min_fill") {
    use_min_fill = true;
  } else if (heuristic == "min_degree") {
    use_min_fill = false;
  } else {
    throw std::invalid_argument(
        "unknown elimination heuristic '" + heuristic +
        "'; expected 'min_fill' or 'min_degree'");
  }
  if (use_min_fill && has_order_min_fill_) {
    std::vector<std::string> out;
    out.reserve(order_min_fill_.size());
    for (int v : order_min_fill_) out.push_back(variables_[v].key);
    return out;
  }
  if (!use_min_fill && has_order_min_degree_) {
    std::vector<std::string> out;
    out.reserve(order_min_degree_.size());
    for (int v : order_min_degree_) out.push_back(variables_[v].key);
    return out;
  }

  const int n = static_cast<int>(variables_.size());
  std::vector<int> order;
  if (n == 0) {
    // empty: nothing to eliminate
  } else if (is_complete()) {
    // Every node has the same cost at every step; the tie-break on the
    // flat-vector index already decides, so skip the O(n^3) greedy loop.
    order.reserve(n);
    for (int v = 0; v < n; ++v) order.push_back(v);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
      return variables_[a].index < variables_[b].index;
    });
  } else {
    std::vector<std::set<int> > adj = moral_adjacency();
    std::vector<bool> alive(n, true);
    order.reserve(n);
    for (int step = 0; step < n; ++step) {
      int best = -1;
      long best_cost = -1;
      int best_tie = -1;
      for (int v = 0; v < n; ++v) {
        if (!alive[v]) continue;
        long cost;
        if (use_min_fill) {
          cost = 0;
          std::vector<int> nb(adj[v].begin(), adj[v].end());
          for (std::size_t a = 0; a < nb.size(); ++a) {
            for (std::size_t b = a + 1; b < nb.size(); ++b) {
              if (!adj[nb[a]].count(nb[b])) ++cost;
            }
          }
        } else {
          cost = static_cast<long>(adj[v].size());
        }
        int tie = variables_[v].index;
        if (best == -1 || cost < best_cost ||
            (cost == best_cost && tie < best_tie)) {
          best = v;
          best_cost = cost;
          best_tie = tie;
        }
      }
      std::vector<int> nb(adj[best].begin(), adj[best].end());
      // Eliminating a variable marries its neighbours: whatever it linked
      // remains coupled once it is summed/maximised out.
      for (std::size_t a = 0; a < nb.size(); ++a) {
        for (std::size_t b = a + 1; b < nb.size(); ++b) {
          adj[nb[a]].insert(nb[b]);
          adj[nb[b]].insert(nb[a]);
        }
      }
      for (int u : nb) adj[u].erase(best);
      adj[best].clear();
      alive[best] = false;
      order.push_back(best);
    }
  }
  if (use_min_fill) {
    order_min_fill_ = order;
    has_order_min_fill_ = true;
  } else {
    order_min_degree_ = order;
    has_order_min_degree_ = true;
  }
  std::vector<std::string> out;
  out.reserve(order.size());
  for (int v : order) out.push_back(variables_[v].key);
  return out;
}

std::vector<std::vector<std::string> > InferenceFactorGraph::compute_cliques(
    const std::vector<int>& order) const {
  // Simulate the elimination: each eliminated variable plus its
  // then-remaining neighbours is a clique of the triangulated graph.
  std::vector<std::set<int> > adj = moral_adjacency();
  std::vector<bool> alive(variables_.size(), true);
  std::vector<std::vector<int> > raw;
  for (int v : order) {
    if (!alive[v]) continue;
    std::vector<int> nb(adj[v].begin(), adj[v].end());
    std::vector<int> clique = nb;
    clique.push_back(v);
    raw.push_back(clique);
    for (std::size_t a = 0; a < nb.size(); ++a) {
      for (std::size_t b = a + 1; b < nb.size(); ++b) {
        adj[nb[a]].insert(nb[b]);
        adj[nb[b]].insert(nb[a]);
      }
    }
    for (int u : nb) adj[u].erase(v);
    adj[v].clear();
    alive[v] = false;
  }
  // Drop non-maximal cliques (those contained in another); largest first.
  std::sort(raw.begin(), raw.end(),
            [](const std::vector<int>& a, const std::vector<int>& b) {
              return a.size() > b.size();
            });
  std::vector<std::set<int> > maximal;
  for (const auto& clique : raw) {
    std::set<int> cs(clique.begin(), clique.end());
    bool subsumed = false;
    for (const auto& kept : maximal) {
      if (std::includes(kept.begin(), kept.end(), cs.begin(), cs.end())) {
        subsumed = true;
        break;
      }
    }
    if (!subsumed) maximal.push_back(cs);
  }
  std::vector<std::vector<std::string> > out;
  out.reserve(maximal.size());
  for (const auto& cs : maximal) {
    std::vector<int> members(cs.begin(), cs.end());
    std::sort(members.begin(), members.end(), [&](int a, int b) {
      return variables_[a].index < variables_[b].index;
    });
    std::vector<std::string> keys;
    keys.reserve(members.size());
    for (int v : members) keys.push_back(variables_[v].key);
    out.push_back(keys);
  }
  return out;
}

std::vector<std::vector<std::string> > InferenceFactorGraph::get_cliques() const {
  if (has_cliques_) return cliques_;
  if (variables_.empty()) {
    cliques_.clear();
    has_cliques_ = true;
    return cliques_;
  }
  if (is_complete()) {
    // A complete graph has a single maximal clique whatever the order.
    std::vector<std::string> keys = get_variable_keys();
    cliques_.assign(1, keys);
    has_cliques_ = true;
    return cliques_;
  }
  // Elimination order over variable positions.
  std::vector<std::string> key_order = get_elimination_order();
  std::vector<int> order;
  order.reserve(key_order.size());
  for (const auto& k : key_order) order.push_back(variable_index_of_.at(k));
  cliques_ = compute_cliques(order);
  has_cliques_ = true;
  return cliques_;
}

int InferenceFactorGraph::get_treewidth() const {
  const auto& cs = get_cliques();
  if (cs.empty()) return 0;
  std::size_t widest = 0;
  for (const auto& c : cs) widest = std::max(widest, c.size());
  return static_cast<int>(widest) - 1;
}

std::vector<InferenceJunctionTreeEdge> InferenceFactorGraph::get_junction_tree_edges() const {
  const auto& cliques = get_cliques();
  const int n = static_cast<int>(cliques.size());
  struct Candidate {
    int i, j, weight;
  };
  std::vector<Candidate> candidates;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      std::set<std::string> a(cliques[i].begin(), cliques[i].end());
      int shared = 0;
      for (const auto& k : cliques[j]) {
        if (a.count(k)) ++shared;
      }
      if (shared > 0) candidates.push_back({i, j, shared});
    }
  }
  // Maximum-weight spanning tree (Kruskal), weights = separator size;
  // a disconnected fit keeps a forest.
  std::stable_sort(candidates.begin(), candidates.end(),
                   [](const Candidate& a, const Candidate& b) {
                     return a.weight > b.weight;
                   });
  DisjointSets ds(n);
  std::vector<InferenceJunctionTreeEdge> edges;
  for (const auto& c : candidates) {
    if (ds.find(c.i) != ds.find(c.j)) {
      ds.unite(c.i, c.j);
      std::set<std::string> a(cliques[c.i].begin(), cliques[c.i].end());
      std::vector<std::string> sep;
      for (const auto& k : cliques[c.j]) {
        if (a.count(k)) sep.push_back(k);
      }
      std::sort(sep.begin(), sep.end(), [&](const std::string& x,
                                            const std::string& y) {
        return index_of(x) < index_of(y);
      });
      InferenceJunctionTreeEdge e;
      e.first = c.i;
      e.second = c.j;
      e.separator = sep;
      edges.push_back(e);
    }
  }
  std::sort(edges.begin(), edges.end(),
            [](const InferenceJunctionTreeEdge& a, const InferenceJunctionTreeEdge& b) {
              if (a.first != b.first) return a.first < b.first;
              return a.second < b.second;
            });
  return edges;
}

std::vector<std::vector<std::string> > InferenceFactorGraph::get_separators() const {
  std::vector<InferenceJunctionTreeEdge> edges = get_junction_tree_edges();
  std::set<std::vector<std::string> > seen;
  for (const auto& e : edges) {
    if (!e.separator.empty()) seen.insert(e.separator);
  }
  std::vector<std::vector<std::string> > out(seen.begin(), seen.end());
  std::sort(out.begin(), out.end(),
            [](const std::vector<std::string>& a,
               const std::vector<std::string>& b) {
              return a.size() > b.size();
            });
  return out;
}

std::vector<std::vector<std::string> > InferenceFactorGraph::get_sampling_blocks()
    const {
  // Group variables by their likelihood-factor neighbourhood: two variables
  // land in the same block exactly when the same set of datasets depends
  // on both.
  const unsigned int n_likelihood = get_number_of_likelihood_factors();
  std::map<std::vector<int>, std::vector<int> > groups;
  for (int v = 0; v < static_cast<int>(variables_.size()); ++v) {
    std::set<int> nb;
    for (int f : incidence_[v]) {
      if (factors_[f].kind == INFERENCE_FACTOR_LIKELIHOOD && factors_[f].fit_index >= 0) {
        nb.insert(factors_[f].fit_index);
      }
    }
    std::vector<int> key(nb.begin(), nb.end());
    groups[key].push_back(v);
  }
  // Cheap blocks first; ties on the smallest flat-vector index so the
  // partition is deterministic. An empty neighbourhood costs nothing to
  // evaluate but explains nothing either, so it sorts last.
  std::vector<std::pair<std::vector<int>, std::vector<int> > > ordered(
      groups.begin(), groups.end());
  std::stable_sort(ordered.begin(), ordered.end(),
                   [&](const auto& a, const auto& b) {
                     const unsigned int cost_a =
                         a.first.empty() ? n_likelihood + 1
                                         : static_cast<unsigned int>(
                                               a.first.size());
                     const unsigned int cost_b =
                         b.first.empty() ? n_likelihood + 1
                                         : static_cast<unsigned int>(
                                               b.first.size());
                     if (cost_a != cost_b) return cost_a < cost_b;
                     int min_a = variables_[a.second.front()].index;
                     int min_b = variables_[b.second.front()].index;
                     return min_a < min_b;
                   });
  std::vector<std::vector<std::string> > out;
  out.reserve(ordered.size());
  for (const auto& kv : ordered) {
    std::vector<int> members = kv.second;
    std::sort(members.begin(), members.end(), [&](int a, int b) {
      return variables_[a].index < variables_[b].index;
    });
    std::vector<std::string> keys;
    keys.reserve(members.size());
    for (int v : members) keys.push_back(variables_[v].key);
    out.push_back(keys);
  }
  return out;
}

int InferenceFactorGraph::get_variable_size(const std::string& key) const {
  auto it = variable_index_of_.find(key);
  return it == variable_index_of_.end() ? 0 : variables_[it->second].size;
}

std::string InferenceFactorGraph::get_variable_role(const std::string& key) const {
  auto it = variable_index_of_.find(key);
  return it == variable_index_of_.end() ? std::string() : variables_[it->second].role;
}

InferenceFactorKind InferenceFactorGraph::get_factor_kind(const std::string& factor_key) const {
  auto it = factor_index_of_.find(factor_key);
  return it == factor_index_of_.end() ? INFERENCE_FACTOR_PRIOR : factors_[it->second].kind;
}

unsigned int InferenceFactorGraph::get_number_of_factors_of_kind(InferenceFactorKind kind) const {
  unsigned int n = 0;
  for (const auto& f : factors_) {
    if (f.kind == kind) ++n;
  }
  return n;
}

int InferenceFactorGraph::factor_cost(const std::vector<std::string>& block) const {
  long long total = 0;
  for (const auto& key : affected_factors(block)) {
    auto it = factor_index_of_.find(key);
    if (it != factor_index_of_.end()) total += factors_[it->second].size;
  }
  const long long cap = std::numeric_limits<int>::max();
  return static_cast<int>(total > cap ? cap : total);
}

int InferenceFactorGraph::block_cost(const std::vector<std::string>& block) const {
  // Two factors, and both are the point of the number. How many local fits
  // have to be redone when this block moves -- which is what it has always
  // been -- times how big the block is, which is the product of its
  // variables' free dimensions.
  //
  // The product, not the sum: this is the size of the block's own state, the
  // same quantity a discrete graphical model calls a clique table. A block of
  // scalars gives exactly the old answer, so callers written before variables
  // had sizes are unaffected.
  const long long fits = static_cast<long long>(affected_fits(block).size());
  long long dim = 1;
  const long long cap = std::numeric_limits<int>::max();
  for (const auto& key : block) {
    auto it = variable_index_of_.find(key);
    if (it == variable_index_of_.end()) continue;
    const long long s = variables_[it->second].size;
    // Saturate rather than wrap: a block this large is beyond sampling
    // anyway, and a negative cost would be read as cheap.
    if (s != 0 && dim > cap / s) return static_cast<int>(cap);
    dim *= s;
  }
  const long long cost = fits * dim;
  return static_cast<int>(cost > cap ? cap : cost);
}

std::string InferenceFactorGraph::to_json() const {
  nlohmann::json j;
  j["format"] = "imp.bff.factorgraph";
  j["version"] = 1;
  nlohmann::json vars = nlohmann::json::array();
  for (const auto& v : variables_) {
    nlohmann::json o;
    o["key"] = v.key;
    o["name"] = v.name;
    o["index"] = v.index;
    o["fit_index"] = v.fit_index;
    o["size"] = v.size;
    o["role"] = v.role;
    vars.push_back(o);
  }
  j["variables"] = vars;
  nlohmann::json facs = nlohmann::json::array();
  for (const auto& f : factors_) {
    nlohmann::json o;
    o["key"] = f.key;
    o["kind"] = static_cast<int>(f.kind);
    nlohmann::json scope = nlohmann::json::array();
    // By key, not by position: a document that referred to positions would
    // be readable only against the graph that wrote it.
    for (int p : f.scope) scope.push_back(variables_[p].key);
    o["scope"] = scope;
    o["fit_index"] = f.fit_index;
    o["size"] = f.size;
    facs.push_back(o);
  }
  j["factors"] = facs;
  return j.dump(2);
}

void InferenceFactorGraph::from_json(const std::string& json) {
  nlohmann::json j = nlohmann::json::parse(json, nullptr, false);
  if (j.is_discarded()) {
    IMP_THROW("InferenceFactorGraph::from_json: not JSON", IMP::ValueException);
  }
  if (j.contains("format") &&
      j.at("format").get<std::string>() != "imp.bff.factorgraph") {
    IMP_THROW("InferenceFactorGraph::from_json: unexpected format '"
                      << j.at("format").get<std::string>() << "'",
              IMP::ValueException);
  }
  if (!j.contains("variables") || !j.contains("factors")) {
    IMP_THROW("InferenceFactorGraph::from_json: the document needs 'variables' and "
              "'factors'", IMP::ValueException);
  }
  InferenceFactorGraph fresh;
  for (const auto& v : j.at("variables")) {
    fresh.add_variable(v.at("key").get<std::string>(),
                       v.value("name", std::string()),
                       v.value("index", -1), v.value("fit_index", -1),
                       v.value("size", 1), v.value("role", std::string()));
  }
  for (const auto& f : j.at("factors")) {
    std::vector<std::string> scope;
    for (const auto& s : f.at("scope")) scope.push_back(s.get<std::string>());
    fresh.add_factor(f.at("key").get<std::string>(),
                     static_cast<InferenceFactorKind>(f.value("kind", 0)), scope,
                     f.value("fit_index", -1), f.value("size", 0));
  }
  *this = fresh;
}

void InferenceFactorGraph::save(const std::string& path) const {
  std::ofstream out(path.c_str());
  if (!out) IMP_THROW("cannot write " << path, IMP::ValueException);
  out << to_json();
  if (!out) IMP_THROW("failed while writing " << path, IMP::ValueException);
}

void InferenceFactorGraph::load(const std::string& path) {
  std::ifstream in(path.c_str());
  if (!in) IMP_THROW("cannot read " << path, IMP::ValueException);
  std::ostringstream all;
  all << in.rdbuf();
  from_json(all.str());
}

std::vector<std::string> InferenceFactorGraph::affected_factors(
    const std::vector<std::string>& changed) const {
  std::set<int> touched;
  for (const auto& k : changed) {
    auto it = variable_index_of_.find(k);
    if (it == variable_index_of_.end()) continue;
    for (int f : incidence_[it->second]) touched.insert(f);
  }
  std::vector<std::string> out;
  out.reserve(touched.size());
  for (int f : touched) out.push_back(factors_[f].key);
  return out;
}

std::vector<int> InferenceFactorGraph::affected_fits(
    const std::vector<std::string>& changed) const {
  std::set<int> out;
  for (const auto& k : changed) {
    auto it = variable_index_of_.find(k);
    if (it == variable_index_of_.end()) continue;
    for (int f : incidence_[it->second]) {
      if (factors_[f].kind == INFERENCE_FACTOR_LIKELIHOOD && factors_[f].fit_index >= 0) {
        out.insert(factors_[f].fit_index);
      }
    }
  }
  return std::vector<int>(out.begin(), out.end());
}

std::vector<std::string> InferenceFactorGraph::get_unexplained_variables() const {
  std::vector<bool> covered(variables_.size(), false);
  for (const auto& f : factors_) {
    if (f.kind != INFERENCE_FACTOR_LIKELIHOOD) continue;
    for (int v : f.scope) covered[v] = true;
  }
  std::vector<const Variable*> unexplained;
  for (std::size_t v = 0; v < variables_.size(); ++v) {
    if (!covered[v]) unexplained.push_back(&variables_[v]);
  }
  std::stable_sort(unexplained.begin(), unexplained.end(),
                   [](const Variable* a, const Variable* b) {
                     return a->index < b->index;
                   });
  std::vector<std::string> out;
  out.reserve(unexplained.size());
  for (const auto* v : unexplained) out.push_back(v->key);
  return out;
}

std::string InferenceFactorGraph::describe() const {
  std::ostringstream lines;
  long long free_numbers = 0;
  for (const auto& v : variables_) free_numbers += v.size;
  lines << "variables      : " << variables_.size() << " ("
        << free_numbers << " free numbers)\n";
  lines << "likelihoods    : " << get_number_of_likelihood_factors() << "\n";
  lines << "priors         : " << get_number_of_factors_of_kind(INFERENCE_FACTOR_PRIOR) << "\n";
  if (get_number_of_factors_of_kind(INFERENCE_FACTOR_HYPER)) {
    lines << "hyper factors  : " << get_number_of_factors_of_kind(INFERENCE_FACTOR_HYPER) << "\n";
  }
  lines << "treewidth      : " << get_treewidth() << "\n";
  lines << "components     : " << connected_components().size() << "\n";
  std::vector<std::vector<std::string> > seps;
  for (const auto& s : get_separators()) {
    if (!s.empty()) seps.push_back(s);
  }
  if (!seps.empty()) {
    lines << "separators     : ";
    for (std::size_t i = 0; i < seps.size(); ++i) {
      if (i) lines << ", ";
      lines << "{";
      for (std::size_t j = 0; j < seps[i].size(); ++j) {
        if (j) lines << ", ";
        const auto& key = seps[i][j];
        auto it = variable_index_of_.find(key);
        const std::string& name =
            (it != variable_index_of_.end() && !variables_[it->second].name.empty())
                ? variables_[it->second].name
                : key;
        lines << name;
      }
      lines << "}";
    }
    lines << "\n";
  } else {
    lines << "separators     : (none - no shared parameters)\n";
  }
  return lines.str();
}

IMPBFF_END_NAMESPACE
