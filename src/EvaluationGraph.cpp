/**
 * \file EvaluationGraph.cpp
 * \brief A graph of nodes that is assembled, then run on demand.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/EvaluationGraph.h>

#include <IMP/bff/internal/json.h>

#include <chrono>
#include <set>
#include <sstream>

IMPBFF_BEGIN_NAMESPACE

namespace {

//! Every node `n` depends on, `n` included, sources first.
/*!
    Walks input links only. That is what makes a label on an intermediate stop
    a run there: nothing downstream is reachable this way, so nothing
    downstream is asked to evaluate.
*/
void collect_upstream(const std::shared_ptr<Node>& n,
                      std::set<const Node*>& seen,
                      std::vector<std::shared_ptr<Node> >& out) {
    if (!n || seen.count(n.get())) return;
    seen.insert(n.get());
    for (const auto& kv : n->get_input_ports()) {
        const std::shared_ptr<Port>& p = kv.second;
        if (!p) continue;
        const std::shared_ptr<Port>& src = p->get_link_ref();
        if (!src) continue;
        const std::shared_ptr<Node> up = src->get_node();
        // A port linked to another port of the *same* node is a shared
        // parameter, not a dependency; recursing would not terminate.
        if (up && up.get() != n.get()) collect_upstream(up, seen, out);
    }
    out.push_back(n);
}

}  // namespace

std::string RunReport::describe() const {
    std::ostringstream s;
    s << nodes_evaluated << " of " << nodes_visited << " nodes evaluated in "
      << seconds << " s";
    return s.str();
}

EvaluationGraph::EvaluationGraph() {}

const EvaluationGraph::Output* EvaluationGraph::find(
        const std::string& label) const {
    auto it = index_of_.find(label);
    return it == index_of_.end() ? nullptr : &outputs_[it->second];
}

void EvaluationGraph::add_output(const std::string& label,
                                 std::shared_ptr<Node> node,
                                 const std::string& port_name) {
    if (label.empty()) {
        IMP_THROW("an output needs a label", IMP::ValueException);
    }
    if (index_of_.count(label)) {
        IMP_THROW("the label '" << label << "' is already taken; remove it "
                  "first if you mean to move it", IMP::ValueException);
    }
    if (!node) {
        IMP_THROW("output '" << label << "' has no node", IMP::ValueException);
    }
    if (!node->get_port(port_name)) {
        IMP_THROW("node '" << node->get_name() << "' has no port '"
                  << port_name << "' for output '" << label << "'",
                  IMP::ValueException);
    }
    Output o;
    o.label = label;
    o.node = node;
    o.port_name = port_name;
    outputs_.push_back(o);
    index_of_[label] = static_cast<int>(outputs_.size()) - 1;
}

void EvaluationGraph::remove_output(const std::string& label) {
    auto it = index_of_.find(label);
    if (it == index_of_.end()) return;
    outputs_.erase(outputs_.begin() + it->second);
    index_of_.clear();
    for (std::size_t i = 0; i < outputs_.size(); ++i) {
        index_of_[outputs_[i].label] = static_cast<int>(i);
    }
}

std::vector<std::string> EvaluationGraph::get_output_labels() const {
    std::vector<std::string> out;
    out.reserve(outputs_.size());
    for (const auto& o : outputs_) out.push_back(o.label);
    return out;
}

std::shared_ptr<Port> EvaluationGraph::get_output_port(
        const std::string& label) const {
    const Output* o = find(label);
    return o ? o->node->get_port(o->port_name) : std::shared_ptr<Port>();
}

std::shared_ptr<Node> EvaluationGraph::get_output_node(
        const std::string& label) const {
    const Output* o = find(label);
    return o ? o->node : std::shared_ptr<Node>();
}

unsigned int EvaluationGraph::get_number_of_outputs() const {
    return static_cast<unsigned int>(outputs_.size());
}

std::vector<std::string> EvaluationGraph::get_dependencies(
        const std::string& label) const {
    const Output* o = find(label);
    if (!o) {
        IMP_THROW("no output labelled '" << label << "'", IMP::ValueException);
    }
    std::set<const Node*> seen;
    std::vector<std::shared_ptr<Node> > order;
    collect_upstream(o->node, seen, order);
    std::vector<std::string> names;
    names.reserve(order.size());
    for (const auto& n : order) names.push_back(n->get_name());
    return names;
}

RunReport EvaluationGraph::run_nodes(
        const std::vector<std::shared_ptr<Node> >& roots) {
    // Everything the roots depend on, gathered first so that the report can
    // say how much was in play and how much of it actually ran. Gathering is
    // cheap next to evaluating and it is the only way to answer honestly.
    std::set<const Node*> seen;
    std::vector<std::shared_ptr<Node> > order;
    for (const auto& r : roots) collect_upstream(r, seen, order);

    std::vector<unsigned long long> before;
    before.reserve(order.size());
    for (const auto& n : order) before.push_back(n->get_evaluation_count());

    const auto t0 = std::chrono::steady_clock::now();
    // update() on each root, not on every node in the order: update() pulls
    // its own upstream, and a node left valid by an earlier root is found
    // valid by a later one. That is what makes two labels sharing a subgraph
    // evaluate it once.
    for (const auto& r : roots) r->update();
    const auto t1 = std::chrono::steady_clock::now();

    RunReport rep;
    rep.nodes_visited = static_cast<int>(order.size());
    for (std::size_t i = 0; i < order.size(); ++i) {
        if (order[i]->get_evaluation_count() != before[i]) ++rep.nodes_evaluated;
    }
    rep.seconds = std::chrono::duration<double>(t1 - t0).count();
    return rep;
}

RunReport EvaluationGraph::run() {
    std::vector<std::shared_ptr<Node> > roots;
    roots.reserve(outputs_.size());
    for (const auto& o : outputs_) roots.push_back(o.node);
    return run_nodes(roots);
}

RunReport EvaluationGraph::run(const std::vector<std::string>& outputs) {
    std::vector<std::shared_ptr<Node> > roots;
    roots.reserve(outputs.size());
    for (const auto& label : outputs) {
        const Output* o = find(label);
        if (!o) {
            // Not a silent no-op: a misspelt label that produced nothing
            // would look exactly like a graph that had nothing to do.
            IMP_THROW("no output labelled '" << label << "'",
                      IMP::ValueException);
        }
        roots.push_back(o->node);
    }
    return run_nodes(roots);
}

std::string EvaluationGraph::to_json() const {
    nlohmann::json j;
    j["format"] = "imp.bff.evaluationgraph";
    j["version"] = 1;
    nlohmann::json outs = nlohmann::json::array();
    for (const auto& o : outputs_) {
        nlohmann::json e;
        e["label"] = o.label;
        e["node"] = o.node->get_name();
        e["port"] = o.port_name;
        outs.push_back(e);
    }
    j["outputs"] = outs;
    return j.dump(2);
}

void EvaluationGraph::from_json(
        const std::string& json,
        const std::map<std::string, std::shared_ptr<Node> >& nodes) {
    nlohmann::json j = nlohmann::json::parse(json, nullptr, false);
    if (j.is_discarded()) {
        IMP_THROW("EvaluationGraph::from_json: not JSON", IMP::ValueException);
    }
    if (j.contains("format") &&
        j.at("format").get<std::string>() != "imp.bff.evaluationgraph") {
        IMP_THROW("EvaluationGraph::from_json: unexpected format '"
                  << j.at("format").get<std::string>() << "'",
                  IMP::ValueException);
    }
    if (!j.contains("outputs")) {
        IMP_THROW("EvaluationGraph::from_json: the document needs 'outputs'",
                  IMP::ValueException);
    }
    EvaluationGraph fresh;
    for (const auto& e : j.at("outputs")) {
        const std::string node_name = e.at("node").get<std::string>();
        auto it = nodes.find(node_name);
        if (it == nodes.end()) {
            IMP_THROW("the document names node '" << node_name
                      << "', which is not among the nodes given; a graph's "
                         "callbacks are the caller's to rebuild",
                      IMP::ValueException);
        }
        fresh.add_output(e.at("label").get<std::string>(), it->second,
                         e.at("port").get<std::string>());
    }
    *this = fresh;
}

IMPBFF_END_NAMESPACE
