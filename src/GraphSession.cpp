/**
 *  \file IMP/bff/GraphSession.cpp
 *  \brief Save and load node graphs in chinet's session format (see
 *         GraphSession.h).
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/bff/GraphSession.h>

#include <IMP/bff/internal/json.h>
#include <IMP/bff/internal/ordered_map.h>

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

namespace {

//! Documents whose object key order must survive a parse: the session
//! format's "ports"/"nodes" maps restore insertion order on load (see
//! internal/ordered_map.h for why this is not nlohmann's default json).
using OrderedJson = nlohmann::basic_json<nlohmann::ordered_map>;

//! Quoted-and-escaped JSON text for one string.
std::string quoted(const std::string& s) { return nlohmann::json(s).dump(); }

//! A JSON value in python json.dumps' default spacing: ", " between
//! elements, ": " after keys. chinet saves with json.dumps, so a re-saved
//! session lines up against a chinet-written one at the byte level for
//! every layout this file controls (an object nested in a prior iterates
//! alphabetically here and in insertion order there; readers do not care).
std::string emit_value(const nlohmann::json& v) {
  switch (v.type()) {
    case nlohmann::json::value_t::array: {
      std::string out = "[";
      for (nlohmann::json::const_iterator it = v.begin(); it != v.end();
           ++it) {
        if (it != v.begin()) out += ", ";
        out += emit_value(*it);
      }
      return out + "]";
    }
    case nlohmann::json::value_t::object: {
      std::string out = "{";
      for (nlohmann::json::const_iterator it = v.begin(); it != v.end();
           ++it) {
        if (it != v.begin()) out += ", ";
        out += quoted(it.key()) + ": " + emit_value(it.value());
      }
      return out + "}";
    }
    default:
      return v.dump();  // null, bool, number, string: no nesting
  }
}

//! One document line from pre-emitted fields, in the order given: chinet's
//! field order, so the files diff against chinet's own writes.
std::string emit_line(
    const std::vector<std::pair<std::string, std::string> >& fields) {
  std::string out = "{";
  for (std::size_t i = 0; i < fields.size(); ++i) {
    if (i > 0) out += ", ";
    out += quoted(fields[i].first) + ": " + fields[i].second;
  }
  return out + "}";
}

//! A {key: uid} map field ("nodes", "ports") in the caller's order.
std::string emit_id_map(
    const std::vector<std::pair<std::string, std::string> >& entries) {
  std::vector<std::pair<std::string, std::string> > fields;
  fields.reserve(entries.size());
  for (std::size_t i = 0; i < entries.size(); ++i) {
    fields.push_back(std::make_pair(entries[i].first,
                                    quoted(entries[i].second)));
  }
  return emit_line(fields);
}

//! json's spelling of a python bool.
std::string emit_bool(bool v) { return v ? "true" : "false"; }

//! chinet's int/float distinction on the wire: a value-type code of 0 or 2
//! is an int port, and numpy's int64 came through json as integers.
//! Anything non-integral in an int port still writes as a float rather
//! than being silently rounded.
std::string emit_number(double e, bool as_int) {
  const double kInt64Max = 9223372036854775807.0;
  if (as_int && e == std::floor(e) && e >= -kInt64Max && e <= kInt64Max) {
    return nlohmann::json(static_cast<std::int64_t>(e)).dump();
  }
  return nlohmann::json(e).dump();
}

//! A port's value field: a JSON scalar for a scalar port, a JSON array
//! for a vector port (empty data writes as []). Link-aware, as chinet's
//! .value is: a linked port saves the value it reads.
std::string emit_port_value(const GraphPort& p) {
  // A bool port writes JSON `true`/`false`, not 1/0: the document should say
  // what the port means, and a reader that restores the value before the
  // value_type (which the session loader does) would otherwise see a plain
  // integer. Codes 0-3 keep chinet's rendering exactly.
  const int element = p.get_element_type();
  if (element == GRAPH_PORT_BOOL) {
    if (p.get_is_vector()) {
      const std::vector<double> values = p.get_value_vector();
      std::string out = "[";
      for (std::size_t i = 0; i < values.size(); ++i) {
        if (i > 0) out += ", ";
        out += (values[i] != 0.0) ? "true" : "false";
      }
      return out + "]";
    }
    if (p.current_size() == 0) return "[]";
    return p.get_value_bool() ? "true" : "false";
  }
  const bool as_int = element == GRAPH_PORT_INT;
  if (p.get_is_vector()) {
    std::vector<std::string> elements;
    if (as_int) {
      // From the exact store, not from get_value_vector(): the doubles are a
      // mirror and a value past 2^53 is already rounded in them.
      const std::vector<long long> exact = p.get_value_vector_int();
      for (std::size_t i = 0; i < exact.size(); ++i) {
        elements.push_back(nlohmann::json(exact[i]).dump());
      }
    } else {
      const std::vector<double> values = p.get_value_vector();
      for (std::size_t i = 0; i < values.size(); ++i) {
        elements.push_back(emit_number(values[i], as_int));
      }
    }
    std::string out = "[";
    for (std::size_t i = 0; i < elements.size(); ++i) {
      if (i > 0) out += ", ";
      out += elements[i];
    }
    return out + "]";
  }
  if (p.current_size() == 0) return "[]";  // scalar flag, no data written yet
  if (as_int) return nlohmann::json(p.get_value_int()).dump();
  return emit_number(p.get_value(), as_int);
}

//! A port's prior, as a JSON object. chinet's slot is dict-or-null; a
//! prior string here that does not parse to an object is refused rather
//! than quietly saved as null (that would be data loss at a distance).
nlohmann::json port_prior_json(const GraphPort& p) {
  if (p.get_prior().empty()) return nlohmann::json();
  nlohmann::json prior =
      nlohmann::json::parse(p.get_prior(), nullptr, false);
  if (prior.is_discarded() || !prior.is_object()) {
    throw std::invalid_argument(
        "IMP::bff::GraphSession::save: the prior of port '" + p.get_name() +
        "' is not JSON object text");
  }
  return prior;
}

//! One port line, chinet's field order.
std::string port_document(const GraphPort& p) {
  std::vector<std::pair<std::string, std::string> > fields;
  fields.push_back(std::make_pair("_id", quoted(p.get_uid())));
  fields.push_back(std::make_pair("name", quoted(p.get_name())));
  fields.push_back(std::make_pair("type", quoted("port")));
  fields.push_back(std::make_pair("precursor", quoted(p.get_precursor())));
  fields.push_back(
      std::make_pair("death", nlohmann::json(p.get_death()).dump()));
  fields.push_back(std::make_pair("fixed", emit_bool(p.get_fixed())));
  fields.push_back(std::make_pair("is_output", emit_bool(p.get_is_output())));
  fields.push_back(
      std::make_pair("is_reactive", emit_bool(p.get_is_reactive())));
  fields.push_back(
      std::make_pair("is_bounded", emit_bool(p.get_is_bounded())));
  fields.push_back(std::make_pair("value", emit_port_value(p)));
  fields.push_back(std::make_pair(
      "bounds",
      "[" + nlohmann::json(p.get_lower_bound()).dump() + ", " +
          nlohmann::json(p.get_upper_bound()).dump() + "]"));
  fields.push_back(std::make_pair(
      "link", p.get_link() ? quoted(p.get_link()->get_uid())
                           : std::string("null")));
  fields.push_back(
      std::make_pair("value_type", nlohmann::json(p.get_value_type()).dump()));
  fields.push_back(std::make_pair("prior", emit_value(port_prior_json(p))));
  return emit_line(fields);
}

//! One node line, chinet's field order. "valid" is the raw flag and
//! "ports" follows the node's own insertion order.
std::string node_document(const GraphNode& n) {
  std::vector<std::pair<std::string, std::string> > port_entries;
  const std::vector<std::string> order = n.get_port_order();
  for (std::size_t i = 0; i < order.size(); ++i) {
    std::shared_ptr<GraphPort> p = n.get_port(order[i]);
    if (p) port_entries.push_back(std::make_pair(order[i], p->get_uid()));
  }
  std::vector<std::pair<std::string, std::string> > fields;
  fields.push_back(std::make_pair("_id", quoted(n.get_uid())));
  fields.push_back(std::make_pair("name", quoted(n.get_name())));
  fields.push_back(std::make_pair("type", quoted("node")));
  fields.push_back(std::make_pair("precursor", quoted(n.get_precursor())));
  fields.push_back(
      std::make_pair("death", nlohmann::json(n.get_death()).dump()));
  fields.push_back(std::make_pair("callback", quoted(n.get_callback())));
  fields.push_back(
      std::make_pair("callback_type", quoted(n.get_callback_type_string())));
  fields.push_back(std::make_pair("valid", emit_bool(n.get_node_valid())));
  fields.push_back(std::make_pair("ports", emit_id_map(port_entries)));
  return emit_line(fields);
}

//! The session line: chinet's field order, nodes in the order added.
std::string session_document(
    const GraphSession& s,
    const std::vector<std::pair<std::string, std::string> >& node_entries) {
  std::vector<std::pair<std::string, std::string> > fields;
  fields.push_back(std::make_pair("_id", quoted(s.get_uid())));
  fields.push_back(std::make_pair("name", quoted(s.get_name())));
  fields.push_back(std::make_pair("type", quoted("session")));
  fields.push_back(std::make_pair("precursor", quoted(s.get_precursor())));
  fields.push_back(
      std::make_pair("death", nlohmann::json(s.get_death()).dump()));
  fields.push_back(std::make_pair("nodes", emit_id_map(node_entries)));
  return emit_line(fields);
}

//! A field of a loaded document, or null if absent.
//!
//! The load path reads documents as OrderedJson, whose objects
//! keep the order the file wrote -- exactly what Python's json.loads gives
//! chinet (a dict), and what the load order depends on: a node's "ports"
//! map restores the node's port insertion order, and that order is
//! observable (operator callbacks compute over the first two inputs in
//! insertion order). nlohmann's default json sorts object keys, so it must
//! not be used where document order is read.
const OrderedJson* doc_find(const OrderedJson& d,
                                       const char* key) {
  OrderedJson::const_iterator it = d.find(key);
  return it == d.end() ? nullptr : &it.value();
}

std::string doc_string(const OrderedJson& d, const char* key,
                       const std::string& fallback) {
  const OrderedJson* v = doc_find(d, key);
  return (v != nullptr && v->is_string()) ? v->get<std::string>() : fallback;
}

bool doc_bool(const OrderedJson& d, const char* key,
              bool fallback) {
  const OrderedJson* v = doc_find(d, key);
  return (v != nullptr && v->is_boolean()) ? v->get<bool>() : fallback;
}

int doc_int(const OrderedJson& d, const char* key, int fallback) {
  const OrderedJson* v = doc_find(d, key);
  return (v != nullptr && v->is_number()) ? v->get<int>() : fallback;
}

//! Restore a port's state from its document, chinet's set_document.
/*!
    Bounds go on before the value (enforcement is still off, so setting
    them does not clip) and the type code after it (the write promotes an
    int port for float input; the document's code is authoritative, as in
    chinet). Enabling is_bounded afterwards clips the loaded value into
    the bounds -- chinet's load did not clip, but chinet-written values
    are always within bounds, so this only bites hand-edited files.
*/
void apply_port_document(const std::shared_ptr<GraphPort>& p,
                         const OrderedJson& d) {
  const OrderedJson* v = doc_find(d, "_id");
  if (v != nullptr && v->is_string()) p->set_uid(v->get<std::string>());
  v = doc_find(d, "name");
  if (v != nullptr && v->is_string()) p->set_name(v->get<std::string>());
  v = doc_find(d, "precursor");
  if (v != nullptr && v->is_string()) {
    p->set_precursor(v->get<std::string>());
  }
  p->set_death(doc_int(d, "death", 0));

  v = doc_find(d, "bounds");
  if (v != nullptr && v->is_array() && v->size() >= 2 && (*v)[0].is_number() &&
      (*v)[1].is_number()) {
    p->set_bounds((*v)[0].get<double>(), (*v)[1].get<double>());
  }
  // The type is declared BEFORE the value is written, because a write no
  // longer retypes the port -- it coerces to whatever the port already is.
  // Restoring the value first would coerce it to the default type and lose
  // it: a saved 2.5 would come back 2, and the type restore that followed
  // could not put the half back.
  v = doc_find(d, "value_type");
  if (v != nullptr && v->is_number()) p->set_value_type(v->get<int>());
  v = doc_find(d, "value");
  if (v != nullptr) {
    // JSON `true`/`false` is not a number, so a bool port's value would be
    // dropped here and then read back as `false` when value_type restored it.
    // The value is restored before the type (see GraphSession.h), so this branch
    // has to accept booleans itself rather than leave them to the conversion.
    if (v->is_array()) {
      // An all-integer array is restored through the exact entry point; a
      // wide value would already be rounded by the time it reached a double.
      bool all_integral = true;
      for (OrderedJson::const_iterator it = v->begin(); it != v->end(); ++it) {
        if (!it->is_number_integer() && !it->is_boolean()) all_integral = false;
      }
      if (all_integral) {
        std::vector<long long> exact;
        for (OrderedJson::const_iterator it = v->begin(); it != v->end();
             ++it) {
          exact.push_back(it->is_boolean() ? (it->get<bool>() ? 1 : 0)
                                           : it->get<long long>());
        }
        p->set_value_vector_int(exact);
      } else {
        std::vector<double> values;
        for (OrderedJson::const_iterator it = v->begin(); it != v->end();
             ++it) {
          values.push_back(it->is_number() ? it->get<double>() : 0.0);
        }
        p->set_value_vector(values);
      }
    } else if (v->is_boolean()) {
      p->set_value_bool(v->get<bool>());
    } else if (v->is_number_integer()) {
      p->set_value(v->get<long long>());
    } else if (v->is_number()) {
      p->set_value(v->get<double>());
    }
  }

  // Applied a second time, deliberately. The first application declared the
  // element type so the value would coerce correctly; this one restores the
  // document's exact code, including chinet's quirk that a float *vector*
  // built by the constructor reports the scalar code. Re-applying a code with
  // the same element type converts nothing, so it costs no precision.
  v = doc_find(d, "value_type");
  if (v != nullptr && v->is_number()) p->set_value_type(v->get<int>());

  // fixed last of the flags that gate writes: value is already restored.
  p->set_fixed(doc_bool(d, "fixed", false));
  p->set_port_type(doc_bool(d, "is_output", false));
  p->set_is_reactive(doc_bool(d, "is_reactive", false));
  p->set_is_bounded(doc_bool(d, "is_bounded", false));
  v = doc_find(d, "prior");
  p->set_prior((v != nullptr && v->is_object()) ? v->dump() : "");
}

//! Restore a node's own state from its document (its ports are restored
//! by the session pass that resolves the ports map).
void apply_node_document(const std::shared_ptr<GraphNode>& n,
                         const OrderedJson& d) {
  const OrderedJson* v = doc_find(d, "_id");
  if (v != nullptr && v->is_string()) n->set_uid(v->get<std::string>());
  v = doc_find(d, "name");
  if (v != nullptr && v->is_string()) n->set_name(v->get<std::string>());
  v = doc_find(d, "precursor");
  if (v != nullptr && v->is_string()) {
    n->set_precursor(v->get<std::string>());
  }
  n->set_death(doc_int(d, "death", 0));
  n->set_callback(doc_string(d, "callback", ""),
                  doc_string(d, "callback_type", ""));
  n->set_valid(doc_bool(d, "valid", false));
}

//! Rebuild a session from its document and the object documents, in
//! chinet's load order: instantiate every node and port, then -- in
//! document order -- restore links and node port maps, then fill the
//! session's node map, then adopt unclaimed ports as free ports.
std::shared_ptr<GraphSession> rebuild_session(
    const OrderedJson& session_doc,
    const std::vector<OrderedJson>& objects) {
  std::shared_ptr<GraphSession> session = std::make_shared<GraphSession>();
  const OrderedJson* v = doc_find(session_doc, "_id");
  if (v != nullptr && v->is_string()) session->set_uid(v->get<std::string>());
  v = doc_find(session_doc, "name");
  if (v != nullptr && v->is_string()) session->set_name(v->get<std::string>());
  v = doc_find(session_doc, "precursor");
  if (v != nullptr && v->is_string()) {
    session->set_precursor(v->get<std::string>());
  }
  session->set_death(doc_int(session_doc, "death", 0));

  std::map<std::string, std::shared_ptr<GraphPort> > ports_by_uid;
  std::map<std::string, std::shared_ptr<GraphNode> > nodes_by_uid;
  for (std::size_t i = 0; i < objects.size(); ++i) {
    const OrderedJson& d = objects[i];
    const std::string type = doc_string(d, "type", "");
    const std::string oid = doc_string(d, "_id", "");
    if (oid.empty()) continue;
    if (type == "port") {
      std::shared_ptr<GraphPort> p = std::make_shared<GraphPort>();
      apply_port_document(p, d);
      ports_by_uid[oid] = p;
    } else if (type == "node") {
      std::shared_ptr<GraphNode> n = std::make_shared<GraphNode>(doc_string(d, "name", ""));
      apply_node_document(n, d);
      nodes_by_uid[oid] = n;
    }
  }

  // Document order, as chinet's pass 3: a link set after its port's node
  // was attached invalidates that node again, and this order reproduces
  // exactly that.
  std::set<std::string> claimed;
  for (std::size_t i = 0; i < objects.size(); ++i) {
    const OrderedJson& d = objects[i];
    const std::string type = doc_string(d, "type", "");
    if (type == "port") {
      const OrderedJson* link = doc_find(d, "link");
      if (link == nullptr || !link->is_string()) continue;
      std::map<std::string, std::shared_ptr<GraphPort> >::iterator self_it =
          ports_by_uid.find(doc_string(d, "_id", ""));
      std::map<std::string, std::shared_ptr<GraphPort> >::iterator target_it =
          ports_by_uid.find(link->get<std::string>());
      if (self_it != ports_by_uid.end() && target_it != ports_by_uid.end()) {
        self_it->second->set_link(target_it->second);
      }
    } else if (type == "node") {
      const OrderedJson* ports_map = doc_find(d, "ports");
      if (ports_map == nullptr || !ports_map->is_object()) continue;
      std::map<std::string, std::shared_ptr<GraphNode> >::iterator node_it =
          nodes_by_uid.find(doc_string(d, "_id", ""));
      if (node_it == nodes_by_uid.end()) continue;
      for (OrderedJson::const_iterator it = ports_map->begin();
           it != ports_map->end(); ++it) {
        if (!it.value().is_string()) continue;
        const std::string port_oid = it.value().get<std::string>();
        std::map<std::string, std::shared_ptr<GraphPort> >::iterator port_it =
            ports_by_uid.find(port_oid);
        if (port_it == ports_by_uid.end()) continue;
        node_it->second->add_port(it.key(), port_it->second,
                                  port_it->second->get_is_output());
        claimed.insert(port_oid);
      }
    }
  }

  const OrderedJson* nodes_map = doc_find(session_doc, "nodes");
  if (nodes_map != nullptr && nodes_map->is_object()) {
    for (OrderedJson::const_iterator it = nodes_map->begin();
         it != nodes_map->end(); ++it) {
      if (!it.value().is_string()) continue;
      std::map<std::string, std::shared_ptr<GraphNode> >::iterator node_it =
          nodes_by_uid.find(it.value().get<std::string>());
      if (node_it != nodes_by_uid.end()) {
        session->add_node(it.key(), node_it->second);
      }
    }
  }

  // Unclaimed ports become the session's free ports, in document order:
  // what chinet's DB held for them, on this side of the port.
  for (std::size_t i = 0; i < objects.size(); ++i) {
    const OrderedJson& d = objects[i];
    if (doc_string(d, "type", "") != "port") continue;
    const std::string oid = doc_string(d, "_id", "");
    if (oid.empty() || claimed.count(oid) > 0) continue;
    std::map<std::string, std::shared_ptr<GraphPort> >::iterator port_it =
        ports_by_uid.find(oid);
    if (port_it != ports_by_uid.end()) session->add_port(port_it->second);
  }

  return session;
}

//! The file's non-empty lines, carriage returns stripped.
std::vector<std::string> nonempty_lines(const std::string& content) {
  std::vector<std::string> lines;
  std::size_t start = 0;
  while (start <= content.size()) {
    std::size_t end = content.find('\n', start);
    if (end == std::string::npos) end = content.size();
    std::string line = content.substr(start, end - start);
    std::size_t first = line.find_first_not_of(" \t\r");
    if (first != std::string::npos) {
      std::size_t last = line.find_last_not_of(" \t\r");
      lines.push_back(line.substr(first, last - first + 1));
    }
    if (end == content.size()) break;
    start = end + 1;
  }
  return lines;
}

}  // namespace

GraphSession::GraphSession() : GraphObject("session") {}

GraphSession::~GraphSession() = default;

void GraphSession::add_node(const std::string& key, std::shared_ptr<GraphNode> node) {
  if (!node) {
    throw std::invalid_argument("IMP::bff::GraphSession::add_node: null node");
  }
  if (nodes_.find(key) == nodes_.end()) node_order_.push_back(key);
  nodes_[key] = node;
}

GraphSession::GraphNodeMap GraphSession::get_nodes() const { return nodes_; }

std::shared_ptr<GraphNode> GraphSession::get_node(const std::string& key) const {
  GraphNodeMap::const_iterator it = nodes_.find(key);
  return it == nodes_.end() ? std::shared_ptr<GraphNode>() : it->second;
}

void GraphSession::add_port(std::shared_ptr<GraphPort> port) {
  if (!port) {
    throw std::invalid_argument("IMP::bff::GraphSession::add_port: null port");
  }
  for (std::size_t i = 0; i < ports_.size(); ++i) {
    if (ports_[i] == port) return;
  }
  ports_.push_back(port);
}

GraphSession::GraphPortList GraphSession::get_ports() const { return ports_; }

std::shared_ptr<GraphPort> GraphSession::get_port(const std::string& name) const {
  for (std::size_t i = 0; i < ports_.size(); ++i) {
    if (ports_[i]->get_name() == name) return ports_[i];
  }
  return std::shared_ptr<GraphPort>();
}

void GraphSession::clear() {
  nodes_.clear();
  node_order_.clear();
  ports_.clear();
}

unsigned int GraphSession::get_number_of_nodes() const {
  return static_cast<unsigned int>(nodes_.size());
}

unsigned int GraphSession::get_number_of_ports() const {
  std::set<const GraphPort*> seen;
  for (GraphNodeMap::const_iterator it = nodes_.begin(); it != nodes_.end(); ++it) {
    const GraphNode::GraphPortMap ports = it->second->get_ports();
    for (GraphNode::GraphPortMap::const_iterator pit = ports.begin();
         pit != ports.end(); ++pit) {
      seen.insert(pit->second.get());
    }
  }
  for (std::size_t i = 0; i < ports_.size(); ++i) {
    seen.insert(ports_[i].get());
  }
  return static_cast<unsigned int>(seen.size());
}

void GraphSession::save(const std::string& path) const {
  std::ofstream out(path.c_str(), std::ios::binary | std::ios::trunc);
  if (!out) {
    throw std::runtime_error("IMP::bff::GraphSession::save: cannot open '" + path +
                             "' for writing");
  }

  std::vector<std::pair<std::string, std::string> > node_entries;
  for (std::size_t i = 0; i < node_order_.size(); ++i) {
    GraphNodeMap::const_iterator it = nodes_.find(node_order_[i]);
    if (it != nodes_.end()) {
      node_entries.push_back(std::make_pair(it->first, it->second->get_uid()));
    }
  }
  out << session_document(*this, node_entries) << "\n";

  // One document per uid, as chinet's registry had: a port shared by two
  // nodes (or held both ways) is written once, under the first path to it.
  std::set<std::string> written;
  for (std::size_t i = 0; i < node_order_.size(); ++i) {
    GraphNodeMap::const_iterator it = nodes_.find(node_order_[i]);
    if (it == nodes_.end()) continue;
    if (!written.insert(it->second->get_uid()).second) continue;
    out << node_document(*it->second) << "\n";
    const std::vector<std::string> order = it->second->get_port_order();
    for (std::size_t k = 0; k < order.size(); ++k) {
      std::shared_ptr<GraphPort> p = it->second->get_port(order[k]);
      if (!p || !written.insert(p->get_uid()).second) continue;
      out << port_document(*p) << "\n";
    }
  }
  for (std::size_t i = 0; i < ports_.size(); ++i) {
    if (!written.insert(ports_[i]->get_uid()).second) continue;
    out << port_document(*ports_[i]) << "\n";
  }

  out.flush();
  if (!out) {
    throw std::runtime_error("IMP::bff::GraphSession::save: failed writing '" +
                             path + "'");
  }
}

std::shared_ptr<GraphSession> GraphSession::load(const std::string& path) {
  std::ifstream in(path.c_str(), std::ios::binary);
  if (!in) {
    throw std::runtime_error("IMP::bff::GraphSession::load: cannot open '" + path +
                             "'");
  }
  const std::string content((std::istreambuf_iterator<char>(in)),
                            std::istreambuf_iterator<char>());
  const std::vector<std::string> lines = nonempty_lines(content);
  if (lines.empty()) {
    return std::shared_ptr<GraphSession>();  // chinet: an empty file loads as None
  }

  // chinet's format test: a JSONL session line is an object without a
  // "session" key; a single line carrying one is a legacy monolithic
  // document, and so is anything the first line cannot parse (the legacy
  // writer could pretty-print across lines).
  const OrderedJson first =
      OrderedJson::parse(lines[0], nullptr, false);
  bool monolithic;
  if (first.is_discarded() || !first.is_object()) {
    monolithic = true;
  } else {
    monolithic = lines.size() == 1 && first.find("session") != first.end();
  }
  if (monolithic) {
    const OrderedJson whole =
        OrderedJson::parse(content, nullptr, false);
    if (whole.is_discarded() || !whole.is_object()) {
      throw std::runtime_error(
          "IMP::bff::GraphSession::load: '" + path +
          "' is neither a JSONL session nor a legacy session document");
    }
    const OrderedJson* session_doc = doc_find(whole, "session");
    const OrderedJson* objects = doc_find(whole, "objects");
    if (session_doc == nullptr || !session_doc->is_object() ||
        objects == nullptr || !objects->is_array()) {
      throw std::runtime_error(
          "IMP::bff::GraphSession::load: '" + path +
          "' is a legacy document without a session object and an objects "
          "array");
    }
    std::vector<OrderedJson> documents;
    for (OrderedJson::const_iterator it = objects->begin();
         it != objects->end(); ++it) {
      documents.push_back(*it);
    }
    return rebuild_session(*session_doc, documents);
  }

  std::vector<OrderedJson> documents;
  for (std::size_t i = 0; i < lines.size(); ++i) {
    OrderedJson d =
        OrderedJson::parse(lines[i], nullptr, false);
    if (d.is_discarded() || !d.is_object()) {
      throw std::runtime_error("IMP::bff::GraphSession::load: '" + path +
                               "' has a malformed line " +
                               std::to_string(i + 1));
    }
    documents.push_back(d);
  }
  for (std::size_t i = 0; i < documents.size(); ++i) {
    if (doc_string(documents[i], "type", "") == "session") {
      return rebuild_session(documents[i], documents);
    }
  }
  return std::shared_ptr<GraphSession>();  // no session line, chinet's None
}

IMPBFF_END_NAMESPACE
