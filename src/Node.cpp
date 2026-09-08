/**
 *  \file Node.cpp
 *  \brief The lazily-evaluated function over ports (see Node.h).
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/bff/Node.h>

#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

namespace {

//! upper() of the callback-type tag, chinet's set_callback dispatch.
std::string to_upper(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), ::toupper);
  return s;
}

//! Elementwise a+b with numpy's scalar broadcasting; mismatched vector
//! sizes are an error there and here.
std::vector<double> broadcast_add(Port* a, Port* b) {
  const bool av = a->get_is_vector();
  const bool bv = b->get_is_vector();
  if (!av && !bv) return std::vector<double>(1, a->get_value() + b->get_value());
  if (!av) {
    const double s = a->get_value();
    std::vector<double> r = b->get_value_vector();
    for (double& e : r) e = s + e;
    return r;
  }
  if (!bv) {
    const double s = b->get_value();
    std::vector<double> r = a->get_value_vector();
    for (double& e : r) e = e + s;
    return r;
  }
  std::vector<double> ra = a->get_value_vector();
  std::vector<double> rb = b->get_value_vector();
  if (ra.size() != rb.size()) {
    throw std::invalid_argument(
        "operator callback over ports of different lengths");
  }
  for (std::size_t i = 0; i < ra.size(); ++i) ra[i] += rb[i];
  return ra;
}

//! Elementwise a*b with the same broadcasting rules.
std::vector<double> broadcast_mul(Port* a, Port* b) {
  const bool av = a->get_is_vector();
  const bool bv = b->get_is_vector();
  if (!av && !bv) return std::vector<double>(1, a->get_value() * b->get_value());
  if (!av) {
    const double s = a->get_value();
    std::vector<double> r = b->get_value_vector();
    for (double& e : r) e = s * e;
    return r;
  }
  if (!bv) {
    const double s = b->get_value();
    std::vector<double> r = a->get_value_vector();
    for (double& e : r) e = e * s;
    return r;
  }
  std::vector<double> ra = a->get_value_vector();
  std::vector<double> rb = b->get_value_vector();
  if (ra.size() != rb.size()) {
    throw std::invalid_argument(
        "operator callback over ports of different lengths");
  }
  for (std::size_t i = 0; i < ra.size(); ++i) ra[i] *= rb[i];
  return ra;
}

//! Write a computed result to an output port the way chinet's
//! out_[name].value = result does: a scalar result (both inputs scalar)
//! through the scalar write, anything else through the vector write. The
//! int operators compute int(v1 op v2), which numpy refuses over vectors
//! longer than one; so does this.
void write_result(Port* out, const std::vector<double>& result, bool as_int,
                  bool result_is_scalar) {
  if (result_is_scalar) {
    if (as_int) {
      // long long, not int: the port's integer store is 64-bit, and an
      // expression yielding 3e9 used to wrap through a 32-bit cast.
      out->set_value(static_cast<long long>(result[0]));
    } else {
      out->set_value(result[0]);
    }
    return;
  }
  if (as_int) {
    throw std::invalid_argument(
        "int operator callbacks need scalar ports");
  }
  out->set_value_vector(result);
}

}  // namespace

Node::Node(const std::string& name) : BaseObject(name) {}

Node::~Node() = default;

std::shared_ptr<Node> Node::make_node(const std::string& name,
                                      const PortMap& ports) {
  std::shared_ptr<Node> n = std::make_shared<Node>(name);
  n->set_ports(ports);
  return n;
}

void Node::set_ports(const PortMap& ports) {
  touch_structure();
  for (const auto& kv : ports) {
    kv.second->set_name(kv.first);
    add_port(kv.first, kv.second, kv.second->get_is_output());
  }
  fill_input_output_port_lookups();
}

void Node::add_port(const std::string& key, std::shared_ptr<Port> port,
                    bool is_output) {
  touch_structure();
  port->set_port_type(is_output);
  try {
    port->set_node(shared_from_this());
  } catch (const std::bad_weak_ptr&) {
    throw std::runtime_error(
        "a Node that owns ports must itself be owned by a std::shared_ptr "
        "(every Python-wrapped node is)");
  }
  if (ports_.find(key) == ports_.end()) port_order_.push_back(key);
  ports_[key] = port;
  fill_input_output_port_lookups();
}

void Node::add_input_port(const std::string& key,
                          std::shared_ptr<Port> port) {
  add_port(key, port, false);
}

void Node::add_output_port(const std::string& key,
                           std::shared_ptr<Port> port) {
  add_port(key, port, true);
}

std::shared_ptr<Port> Node::get_port(const std::string& name) const {
  PortMap::const_iterator it = ports_.find(name);
  return it == ports_.end() ? std::shared_ptr<Port>() : it->second;
}

std::shared_ptr<Port> Node::get_input_port(const std::string& name) const {
  PortMap::const_iterator it = in_.find(name);
  return it == in_.end() ? std::shared_ptr<Port>() : it->second;
}

std::shared_ptr<Port> Node::get_output_port(const std::string& name) const {
  PortMap::const_iterator it = out_.find(name);
  return it == out_.end() ? std::shared_ptr<Port>() : it->second;
}

Node::PortMap Node::get_ports() const { return ports_; }

Node::PortMap Node::get_input_ports() const { return in_; }

Node::PortMap Node::get_output_ports() const { return out_; }

void Node::set_callback(const std::string& callback,
                        const std::string& callback_type) {
  touch_structure();
  callback_ = callback;
  callback_type_string_ = callback_type;
  const std::string t = to_upper(callback_type);
  if (t == "C") {
    callback_type_ = 0;
  } else if (t == "CLASS" || t == "CPP" || t == "C++") {
    callback_type_ = 1;
  } else {
    callback_type_ = -1;
  }
}

void Node::set_callback_function(Callback callback) {
  touch_structure();
  callback_class_ = callback;
  callback_type_ = 1;
}

const std::string& Node::get_callback() const { return callback_; }

int Node::get_callback_type() const { return callback_type_; }

const std::string& Node::get_callback_type_string() const {
  return callback_type_string_;
}

bool Node::get_node_valid() const { return node_valid_; }

const std::vector<std::string>& Node::get_port_order() const {
  return port_order_;
}

bool Node::inputs_valid() const {
  for (const auto& kv : in_) {
    if (kv.second->is_linked()) {
      std::shared_ptr<Port> output_port = kv.second->get_link();
      std::shared_ptr<Node> output_node = output_port->get_node();
      if (output_node && output_node.get() == this) continue;
      if (output_node && !output_node->is_valid()) return false;
    }
  }
  return true;
}

bool Node::is_valid() const {
  if (in_.empty()) return true;
  if (!inputs_valid()) return false;
  return node_valid_;
}

void Node::set_valid(bool v) { node_valid_ = v; }

void Node::evaluate() {
  // chinet: a node with neither a callback object nor an operator does
  // nothing and stays invalid.
  if (!callback_class_ && callback_type_ < 0) return;

  refresh_plan();

  if (callback_type_ == 0 && plan_.op != OP_NONE && plan_.operand_a &&
      plan_.operand_b) {
    Port* v1 = plan_.operand_a;
    Port* v2 = plan_.operand_b;
    if (plan_.op_output == nullptr) {
      throw std::invalid_argument(
          "operator callback '" + callback_ +
          "' writes to the output port keyed by the node's name ('" +
          get_name() + "'), which this node does not have");
    }
    const bool scalar = !v1->get_is_vector() && !v2->get_is_vector();
    if (scalar) {
      // The hot path: no vector, no allocation, no map lookup. Sampling
      // drives millions of these and the arithmetic is two loads and an op.
      const double a = v1->get_value(), b = v2->get_value();
      const double r = (plan_.op == OP_ADD) ? a + b : a * b;
      if (plan_.as_int) {
        plan_.op_output->set_value(static_cast<long long>(r));
      } else {
        plan_.op_output->set_value(r);
      }
    } else {
      write_result(plan_.op_output,
                   plan_.op == OP_ADD ? broadcast_add(v1, v2)
                                      : broadcast_mul(v1, v2),
                   plan_.as_int, false);
    }
  }

  if (callback_class_) callback_class_(in_, out_);

  // A node sharing one of our output ports now reads stale results.
  for (Port* out : plan_.outputs) {
    if (const std::shared_ptr<Node> n = out->get_node()) {
      if (n.get() != this) n->set_valid(false);
    }
  }
  node_valid_ = true;
}

//! Rebuild the cached plan when the node's structure has moved on.
void Node::refresh_plan() const {
  if (plan_.epoch == structure_epoch_ && plan_.name == get_name()) return;
  Plan plan;
  plan.epoch = structure_epoch_;
  plan.name = get_name();

  for (const auto& kv : in_) {
    if (kv.second->is_linked()) plan.linked_inputs.push_back(kv.second.get());
  }
  for (const auto& kv : out_) plan.outputs.push_back(kv.second.get());

  if (callback_ == "addition_double") {
    plan.op = OP_ADD;
  } else if (callback_ == "addition_int") {
    plan.op = OP_ADD;
    plan.as_int = true;
  } else if (callback_ == "multiply_double") {
    plan.op = OP_MUL;
  } else if (callback_ == "multiply_int") {
    plan.op = OP_MUL;
    plan.as_int = true;
  }

  if (plan.op != OP_NONE) {
    // The operator protocol: the first two inputs in insertion order, the
    // result written to the output port keyed by the node's own name.
    Port* operands[2] = {nullptr, nullptr};
    int found = 0;
    for (const std::string& key : port_order_) {
      PortMap::const_iterator it = in_.find(key);
      if (it == in_.end()) continue;
      operands[found++] = it->second.get();
      if (found == 2) break;
    }
    if (found == 2) {
      plan.operand_a = operands[0];
      plan.operand_b = operands[1];
      PortMap::const_iterator o = out_.find(get_name());
      if (o != out_.end()) plan.op_output = o->second.get();
    } else {
      // Fewer than two inputs: chinet computes nothing at all.
      plan.op = OP_NONE;
    }
  }

  plan_ = plan;
}

void Node::touch_structure() { ++structure_epoch_; }

void Node::update() {
  refresh_plan();
  for (Port* p : plan_.linked_inputs) {
    const std::shared_ptr<Port>& source_port = p->get_link_ref();
    if (!source_port) continue;  // unlinked since the plan was built
    if (const std::shared_ptr<Node> source_node = source_port->get_node()) {
      // Not `this`: an input linked to another port of the *same* node --
      // two variables of one expression that are one number, which is how a
      // parameter shared inside a single model is spelt -- needs no upstream
      // evaluation, it needs the copy below. Recursing would re-enter this
      // node forever and take the stack with it. `inputs_valid()` has always
      // skipped the self-link; this is the same rule, in the place that
      // walks.
      if (source_node.get() != this && !source_node->is_valid()) {
        source_node->update();
      }
    }
    // p.value = source.value: the write follows the source's shape.
    if (!source_port->get_is_vector() && source_port->current_size() == 1) {
      p->set_value(source_port->get_value());
    } else {
      // `get_values_ref`, not `get_value_vector`: this runs once per linked
      // input per update, and a joint objective has one full residual vector
      // per member coming through it on every iteration of a fit.
      p->set_value_vector(source_port->get_values_ref());
    }
  }
  // node_valid_, not is_valid(): the loop above has just brought every
  // linked input's source node up to date, so inputs_valid() is known true
  // and re-deriving it would walk the whole upstream graph a second time
  // for every node in it.
  if (!node_valid_) evaluate();
}

std::string Node::describe() const {
  std::ostringstream out;
  out << "Node(name='" << get_name() << "', uid='" << get_uid() << "'";
  out << ", ports=[";
  bool first = true;
  for (const std::string& key : port_order_) {
    if (!first) out << ", ";
    first = false;
    out << (in_.count(key) ? "-" : "+") << key;
  }
  out << "]";
  out << ", valid=" << (node_valid_ ? "true" : "false");
  out << ", callback='" << callback_ << "'";
  out << ")";
  return out.str();
}

void Node::fill_input_output_port_lookups() {
  touch_structure();
  in_.clear();
  out_.clear();
  for (const std::string& key : port_order_) {
    PortMap::const_iterator it = ports_.find(key);
    if (it == ports_.end()) continue;
    if (it->second->get_is_output()) {
      out_[key] = it->second;
    } else {
      in_[key] = it->second;
    }
  }
}

IMPBFF_END_NAMESPACE
