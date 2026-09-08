/**
 *  \file Port.cpp
 *  \brief The reactive value cell (see Port.h).
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/bff/Port.h>

#include <IMP/bff/Node.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>

IMPBFF_BEGIN_NAMESPACE

namespace {

//! A fresh uuid-shaped id, standing in for chinet's uuid4 oids.
std::string new_uid() {
  static std::random_device rd;
  std::uint64_t a = (static_cast<std::uint64_t>(rd()) << 32) ^ rd();
  std::uint64_t b = (static_cast<std::uint64_t>(rd()) << 32) ^ rd();
  char buf[40];
  std::snprintf(buf, sizeof(buf),
                "%08x-%04x-%04x-%04x-%04x%08x",
                static_cast<unsigned>((a >> 32) & 0xffffffffu),
                static_cast<unsigned>((a >> 16) & 0xffffu),
                static_cast<unsigned>(a & 0xffffu),
                static_cast<unsigned>((b >> 48) & 0xffffu),
                static_cast<unsigned>((b >> 16) & 0xffffu),
                static_cast<unsigned>(b & 0xffffu));
  return std::string(buf);
}

//! chinet sanitises float writes: NaN becomes the smallest normal, +inf the
//! largest double, -inf the smallest. Optimisers write both.
double sanitize(double v) {
  if (std::isnan(v)) return 2.2250738585072014e-308;  // np.finfo(float).tiny
  if (v > 1.7976931348623157e308) return 1.7976931348623157e308;
  if (v < -1.7976931348623157e308) return -1.7976931348623157e308;
  return v;
}

//! A dependency-graph vertex, exactly chinet's _dependency_vertex: a port
//! attached to a node shares that node as its vertex; an unattached port
//! (e.g. the standalone ports backing chisurf fitting parameters) is its
//! own vertex. Held by shared_ptr so the traversal stays alive while the
//! ids it compares by are in use.
struct Vertex {
  void* id;
  std::shared_ptr<Node> node;  // set iff the vertex is a node
  std::shared_ptr<Port> port;  // set iff the vertex is a standalone port
};

Vertex vertex_of(const std::shared_ptr<Port>& p) {
  std::shared_ptr<Node> n = p->get_node();
  if (n) return Vertex{n.get(), n, std::shared_ptr<Port>()};
  return Vertex{p.get(), std::shared_ptr<Node>(), p};
}

//! The vertices v directly depends on through existing links (chinet's
//! _vertex_dependencies): a node depends on the sources of its linked
//! input ports; a standalone port on whatever it is linked to.
void vertex_dependencies(const Vertex& v, std::vector<Vertex>& out) {
  if (v.node) {
    for (const auto& kv : v.node->get_input_ports()) {
      if (kv.second->is_linked()) out.push_back(vertex_of(kv.second->get_link()));
    }
  } else if (v.port && v.port->is_linked()) {
    out.push_back(vertex_of(v.port->get_link()));
  }
}

}  // namespace

BaseObject::BaseObject(const std::string& name)
    : name_(name), uid_(new_uid()), precursor_(uid_), death_(0) {}

BaseObject::~BaseObject() = default;

const std::string& BaseObject::get_name() const { return name_; }

void BaseObject::set_name(const std::string& name) { name_ = name; }

const std::string& BaseObject::get_uid() const { return uid_; }

void BaseObject::set_uid(const std::string& uid) { uid_ = uid; }

const std::string& BaseObject::get_precursor() const { return precursor_; }

void BaseObject::set_precursor(const std::string& precursor) {
  precursor_ = precursor;
}

int BaseObject::get_death() const { return death_; }

void BaseObject::set_death(int death) { death_ = death; }

int port_value_type_element(int t) {
  switch (t) {
    case PORT_FLOAT:
    case PORT_FLOAT_VECTOR: return PORT_FLOAT;
    case PORT_BOOL:
    case PORT_BOOL_VECTOR: return PORT_BOOL;
    default: return PORT_INT;
  }
}

bool port_value_type_is_vector(int t) {
  return t == PORT_INT_VECTOR || t == PORT_FLOAT_VECTOR ||
         t == PORT_BOOL_VECTOR;
}

int port_value_type_of(int element, bool is_vector) {
  switch (port_value_type_element(element)) {
    case PORT_FLOAT: return is_vector ? PORT_FLOAT_VECTOR : PORT_FLOAT;
    case PORT_BOOL: return is_vector ? PORT_BOOL_VECTOR : PORT_BOOL;
    default: return is_vector ? PORT_INT_VECTOR : PORT_INT;
  }
}

std::string port_value_type_name(int t) {
  switch (t) {
    case PORT_INT: return "int";
    case PORT_FLOAT: return "float";
    case PORT_INT_VECTOR: return "int[]";
    case PORT_FLOAT_VECTOR: return "float[]";
    case PORT_BOOL: return "bool";
    case PORT_BOOL_VECTOR: return "bool[]";
    default: return "unknown(" + std::to_string(t) + ")";
  }
}

Port::Port()
    : BaseObject(), buffer_(1, 0.0), value_type_(PORT_FLOAT),
      is_vector_(false) {
  // Float, not int. chinet defaulted this to the integer code because the
  // type was inferred from the first write and so the default never survived
  // contact with data. Now that the declared type is what a write coerces to,
  // an int default would silently truncate the first float anybody stored --
  // and a float port is what almost every caller wants (every fitting
  // parameter, every spectrum, every residual).
}

Port::Port(double value, bool fixed, bool is_output, bool is_reactive,
           bool is_bounded, double lb, double ub, int value_type,
           const std::string& name)
    : BaseObject(name),
      buffer_(1, value),
      value_type_(value_type),
      is_vector_(false),
      fixed_(fixed),
      is_output_(is_output),
      is_reactive_(is_reactive),
      is_bounded_(is_bounded),
      lb_(lb),
      ub_(ub) {
  // chinet's __init__ dtype rules: a float argument is float input, so an
  // int-typed scalar port is promoted (only code 0 -> 1; chinet never
  // touches the vector codes in the constructor). An int argument (the int
  // constructor below) is int input and keeps the port integral.
  if (value_type_ == 0) value_type_ = 1;
  if (port_value_type_element(value_type_) != PORT_FLOAT) {
    // The argument was a double but the port is declared integral: convert the
    // slots once, here, rather than leaving them as raw doubles a later read
    // would misinterpret.
    store_doubles(std::vector<double>(1, buffer_[0]));
  }
  finalize_storage();
}

Port::Port(long long value, bool fixed, bool is_output, bool is_reactive,
           bool is_bounded, double lb, double ub, int value_type,
           const std::string& name)
    : BaseObject(name),
      buffer_(1, int_as_slot(value)),
      value_type_(value_type),
      is_vector_(false),
      fixed_(fixed),
      is_output_(is_output),
      is_reactive_(is_reactive),
      is_bounded_(is_bounded),
      lb_(lb),
      ub_(ub) {
  // int input on an int port keeps the type code (chinet stores int64);
  // a float-declared port must hold the value as a double slot instead.
  if (port_value_type_element(value_type_) == PORT_FLOAT) {
    buffer_.assign(1, static_cast<double>(value));
  }
  finalize_storage();
}

Port::Port(bool value, bool fixed, bool is_output, bool is_reactive,
           bool is_bounded, double lb, double ub, int value_type,
           const std::string& name)
    : BaseObject(name),
      buffer_(1, int_as_slot(value ? 1 : 0)),
      value_type_(value_type),
      is_vector_(false),
      fixed_(fixed),
      is_output_(is_output),
      is_reactive_(is_reactive),
      is_bounded_(is_bounded),
      lb_(lb),
      ub_(ub) {
  // A bool argument declares a flag port; unlike the int constructor there is
  // no promotion rule to apply, because bool does not promote.
  if (port_value_type_element(value_type_) != PORT_BOOL) {
    value_type_ = port_value_type_of(PORT_BOOL, false);
  }
  finalize_storage();
}

Port::Port(const std::vector<double>& values, bool fixed, bool is_output,
           bool is_reactive, bool is_bounded, double lb, double ub,
           int value_type, const std::string& name)
    : BaseObject(name),
      buffer_(values),
      value_type_(value_type),
      is_vector_(true),
      fixed_(fixed),
      is_output_(is_output),
      is_reactive_(is_reactive),
      is_bounded_(is_bounded),
      lb_(lb),
      ub_(ub) {
  // chinet promotes an int-typed port fed float data; note it leaves a
  // vector port at code 1 here (the setter path is what writes 3), and the
  // port keeps that quirk.
  bool is_f = (value_type_ == 1 || value_type_ == 3) ||
              !(buffer_.empty() && (value_type_ == 0 || value_type_ == 2));
  if (value_type_ == 0 && is_f) value_type_ = 1;
  if (!is_f) {
    for (double& e : buffer_) e = std::trunc(e);
  }
  if (is_bounded_) clip_to_bounds();
}

Port::~Port() = default;

double Port::get_value() const {
  const Port* src = link_ ? link_.get() : this;
  if (src->buffer_.empty()) return 0.0;  // chinet's .item() raises on empty;
                                         // callers of an empty port have
                                         // nothing to read
  return src->element_as_double(0);
}

std::vector<double> Port::get_value_vector() const {
  return get_values_ref();
}

const std::vector<double>& Port::get_values_ref() const {
  const Port* src = link_ ? link_.get() : this;
  // A float port hands out the store itself: no copy, no branch downstream.
  // This is the read the evaluation hot path makes every iteration.
  if (port_value_type_element(src->value_type_) == PORT_FLOAT) {
    return src->buffer_;
  }
  // Anything else materialises a view. It is a cache, not a second store --
  // it is never read back, and every write drops it.
  if (!src->double_cache_valid_) {
    src->double_cache_.resize(src->buffer_.size());
    for (std::size_t i = 0; i < src->buffer_.size(); ++i) {
      src->double_cache_[i] = src->element_as_double(i);
    }
    src->double_cache_valid_ = true;
  }
  return src->double_cache_;
}

void Port::get_value_view(double** out_values, int* n_out_values) const {
  const std::vector<double>& d = get_values_ref();
  *n_out_values = static_cast<int>(d.size());
  // numpy.i's ARGOUTVIEWM typemap hands the buffer to a numpy array whose
  // capsule destructor calls free() -- so this is malloc, not new[], and
  // the plain-double buffer keeps the C allocator contract.
  *out_values = static_cast<double*>(std::malloc(d.size() * sizeof(double)));
  if (!*out_values && !d.empty()) throw std::bad_alloc();
  std::copy(d.begin(), d.end(), *out_values);
}

void Port::set_value(double v) { write_value(v, /*input_is_float=*/true); }

void Port::set_value(long long v) {
  if (fixed_) return;
  const int element = port_value_type_element(value_type_);
  if (element != PORT_INT) {
    // A float port takes it as a float, a bool port as a flag. Accepted,
    // coerced, never retyped.
    write_value(static_cast<double>(v), false);
    return;
  }
  value_type_ = PORT_INT;
  is_vector_ = false;
  store_ints(std::vector<long long>(1, v));
  finalize_storage();
  update_attached_node();
  propagate_to_followers();
}

long long Port::get_value_int() const {
  const Port* src = link_ ? link_.get() : this;
  if (src->buffer_.empty()) return 0;
  if (port_value_type_element(src->value_type_) == PORT_FLOAT) {
    return static_cast<long long>(std::trunc(src->buffer_[0]));
  }
  return slot_as_int(src->buffer_[0]);
}

std::vector<long long> Port::get_value_vector_int() const {
  const Port* src = link_ ? link_.get() : this;
  const bool as_float =
      port_value_type_element(src->value_type_) == PORT_FLOAT;
  std::vector<long long> out;
  out.reserve(src->buffer_.size());
  for (std::size_t i = 0; i < src->buffer_.size(); ++i) {
    out.push_back(as_float
                      ? static_cast<long long>(std::trunc(src->buffer_[i]))
                      : slot_as_int(src->buffer_[i]));
  }
  return out;
}

void Port::set_value_bool(bool v) {
  // Writes a truth value into whatever the port already is; it does not make
  // the port a bool port. Declaring a flag is the constructor's job, or
  // set_value_type()'s -- see PortValueType.
  write_value(v ? 1.0 : 0.0, false);
}

bool Port::get_value_bool() const {
  const Port* src = link_ ? link_.get() : this;
  return !src->buffer_.empty() && src->element_as_double(0) != 0.0;
}

std::vector<int> Port::get_value_vector_bool() const {
  const Port* src = link_ ? link_.get() : this;
  std::vector<int> out;
  out.reserve(src->buffer_.size());
  for (std::size_t i = 0; i < src->buffer_.size(); ++i) {
    out.push_back(src->element_as_double(i) != 0.0 ? 1 : 0);
  }
  return out;
}

void Port::set_value_vector_int(const std::vector<long long>& v) {
  if (fixed_) return;
  const int element = port_value_type_element(value_type_);
  if (element != PORT_INT) {
    write_vector(std::vector<double>(v.begin(), v.end()), false);
    return;
  }
  value_type_ = PORT_INT_VECTOR;
  is_vector_ = true;
  store_ints(v);
  finalize_storage();
  update_attached_node();
  propagate_to_followers();
}

void Port::set_value_vector(const std::vector<double>& v) {
  write_vector(v, /*input_is_float=*/true);
}

void Port::set_values_array(double* in_values, int n_values) {
  const int n = std::max(0, n_values);
  write_vector(std::vector<double>(in_values, in_values + n),
               /*input_is_float=*/true);
}


long long Port::slot_as_int(double slot) {
  long long v;
  std::memcpy(&v, &slot, sizeof(v));
  return v;
}

double Port::int_as_slot(long long v) {
  double slot;
  std::memcpy(&slot, &v, sizeof(slot));
  return slot;
}

double Port::element_as_double(std::size_t i) const {
  if (port_value_type_element(value_type_) == PORT_FLOAT) return buffer_[i];
  return static_cast<double>(slot_as_int(buffer_[i]));
}

void Port::store_doubles(const std::vector<double>& v) {
  const int element = port_value_type_element(value_type_);
  if (element == PORT_FLOAT) {
    buffer_ = v;
  } else {
    buffer_.resize(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) {
      const double e = v[i];
      const long long as_int =
          (element == PORT_BOOL) ? ((e != 0.0 && !std::isnan(e)) ? 1 : 0)
                                 : static_cast<long long>(std::trunc(e));
      buffer_[i] = int_as_slot(as_int);
    }
  }
  double_cache_valid_ = false;
}

void Port::store_ints(const std::vector<long long>& v) {
  const int element = port_value_type_element(value_type_);
  buffer_.resize(v.size());
  for (std::size_t i = 0; i < v.size(); ++i) {
    if (element == PORT_FLOAT) {
      buffer_[i] = static_cast<double>(v[i]);
    } else {
      buffer_[i] = int_as_slot((element == PORT_BOOL) ? (v[i] != 0 ? 1 : 0)
                                                      : v[i]);
    }
  }
  double_cache_valid_ = false;
}

void Port::finalize_storage() {
  if (is_bounded_) clip_to_bounds();
  double_cache_valid_ = false;
}

bool Port::is_bool_typed() const {
  return port_value_type_element(value_type_) == PORT_BOOL;
}

void Port::write_value(double v, bool input_is_float) {
  if (fixed_) return;
  // The element type is DECLARED at construction and never changes on a
  // write. A value of another kind is accepted and coerced -- writing 1.5 to
  // an integer port stores 1, it does not turn the port into a float port.
  //
  // This is a deliberate divergence from chinet, which inferred the dtype
  // from every write (numpy's promotion rules) and so let a port silently
  // become something else halfway through a session. A port that can change
  // type is not typed; whatever holds a reference to it had already decided
  // what it was.
  const int element = port_value_type_element(value_type_);
  if (input_is_float && sanitize_) v = sanitize(v);
  value_type_ = port_value_type_of(element, false);
  is_vector_ = false;
  store_doubles(std::vector<double>(1, v));   // coerces to the element type
  finalize_storage();
  update_attached_node();
  propagate_to_followers();
}

void Port::write_vector(const std::vector<double>& v, bool input_is_float) {
  if (fixed_) return;
  // As write_value(): the element type is fixed, only the shape follows the
  // data. Vector-ness is not part of what was declared -- a port that is
  // handed an array becomes an array port of the same element type.
  const int element = port_value_type_element(value_type_);
  std::vector<double> incoming = v;
  if (input_is_float && sanitize_ && element == PORT_FLOAT) {
    for (double& e : incoming) e = sanitize(e);
  }
  value_type_ = port_value_type_of(element, true);
  is_vector_ = true;
  store_doubles(incoming);
  finalize_storage();
  update_attached_node();
  propagate_to_followers();
}

int Port::get_value_type() const { return value_type_; }

void Port::set_value_type(int t) {
  // Convert through the *values*, never through the raw slots: a slot means
  // a double under one element type and an int64 bit pattern under another,
  // so reinterpreting one as the other is nonsense rather than a conversion.
  // Reading both forms out first is also what keeps an integer exact across
  // an int -> int retype, which is the path a document restore takes.
  const bool was_float = port_value_type_element(value_type_) == PORT_FLOAT;
  std::vector<double> as_doubles(buffer_.size());
  std::vector<long long> as_ints(buffer_.size());
  for (std::size_t i = 0; i < buffer_.size(); ++i) {
    as_doubles[i] = element_as_double(i);
    as_ints[i] = was_float
                     ? static_cast<long long>(std::trunc(buffer_[i]))
                     : slot_as_int(buffer_[i]);
  }
  value_type_ = t;
  if (port_value_type_element(value_type_) == PORT_FLOAT) {
    store_doubles(as_doubles);
  } else {
    store_ints(as_ints);   // exact when the values were already integers
  }
}

bool Port::get_is_vector() const { return is_vector_; }

unsigned int Port::current_size() const {
  return static_cast<unsigned int>(buffer_.size());
}

bool Port::is_valid() const { return true; }

bool Port::get_fixed() const { return fixed_; }

void Port::set_fixed(bool v) { fixed_ = v; }

bool Port::get_is_output() const { return is_output_; }

void Port::set_port_type(bool v) { is_output_ = v; }

bool Port::get_is_reactive() const { return is_reactive_; }

void Port::set_is_reactive(bool v) { is_reactive_ = v; }

bool Port::get_is_bounded() const { return is_bounded_; }

void Port::set_is_bounded(bool v) {
  is_bounded_ = v;
  if (is_bounded_) clip_to_bounds();
}

double Port::get_lower_bound() const { return lb_; }

double Port::get_upper_bound() const { return ub_; }

void Port::set_bounds(double lb, double ub) {
  lb_ = lb;
  ub_ = ub;
  if (is_bounded_) clip_to_bounds();
}

void Port::clip_to_bounds() {
  const int element = port_value_type_element(value_type_);
  for (std::size_t i = 0; i < buffer_.size(); ++i) {
    const double clipped = std::min(std::max(element_as_double(i), lb_), ub_);
    if (element == PORT_FLOAT) {
      buffer_[i] = clipped;
    } else if (element == PORT_BOOL) {
      buffer_[i] = int_as_slot((clipped != 0.0) ? 1 : 0);
    } else {
      buffer_[i] = int_as_slot(static_cast<long long>(std::trunc(clipped)));
    }
  }
  double_cache_valid_ = false;
}

const std::string& Port::get_prior() const { return prior_; }

void Port::set_prior(const std::string& json) { prior_ = json; }

std::shared_ptr<Port> Port::get_link() const { return link_; }

void Port::set_link(std::shared_ptr<Port> v) {
  if (v.get() == link_.get()) return;
  if (v && would_create_cycle(v)) {
    throw LinkCycleError(
        "Linking these ports would create a cycle; "
        "the port graph must remain acyclic.");
  }
  unlink();
  if (!v) return;
  std::weak_ptr<Port> self;
  try {
    self = shared_from_this();
  } catch (const std::bad_weak_ptr&) {
    throw std::runtime_error(
        "a Port that follows another Port must be owned by a "
        "std::shared_ptr (every Python-wrapped port is)");
  }
  link_ = v;
  v->linked_to_.push_back(self);
  touch_attached_node_structure();
  update_attached_node();
}

bool Port::unlink() {
  if (link_) {
    // Drop this port (and any followers that died unlinked) from the
    // source's follower list, then break the link. chinet's list remove
    // only matches the port itself; expired followers are skipped by every
    // read there, and skipping them here too keeps the list compact.
    auto& followers = link_->linked_to_;
    followers.erase(
        std::remove_if(followers.begin(), followers.end(),
                       [this](const std::weak_ptr<Port>& w) {
                         std::shared_ptr<Port> p = w.lock();
                         return !p || p.get() == this;
                       }),
        followers.end());
    link_.reset();
    touch_attached_node_structure();
    update_attached_node();
  }
  return true;
}

bool Port::is_linked() const { return link_ != nullptr; }

bool Port::would_create_cycle(const std::shared_ptr<Port>& v) {
  if (!v) return false;

  std::shared_ptr<Port> self;
  try {
    self = shared_from_this();
  } catch (const std::bad_weak_ptr&) {
    throw std::runtime_error(
        "a Port whose links are being checked must be owned by a "
        "std::shared_ptr (every Python-wrapped port is)");
  }
  Vertex sv = vertex_of(self);
  Vertex vv = vertex_of(v);
  // Ports of the same node may be linked (a node may read its own output);
  // a standalone port linked to itself is a degenerate self-cycle.
  if (sv.id == vv.id) return this == v.get();

  // The subgraph reachable from vv following dependency edges, plus sv.
  std::vector<Vertex> nodes;
  nodes.push_back(sv);
  std::vector<void*> seen;
  seen.push_back(sv.id);
  std::vector<Vertex> stack;
  stack.push_back(vv);
  while (!stack.empty()) {
    Vertex cur = stack.back();
    stack.pop_back();
    if (std::find(seen.begin(), seen.end(), cur.id) != seen.end()) continue;
    seen.push_back(cur.id);
    nodes.push_back(cur);
    std::vector<Vertex> deps;
    vertex_dependencies(cur, deps);
    for (const Vertex& dep : deps) {
      if (std::find(seen.begin(), seen.end(), dep.id) == seen.end()) {
        stack.push_back(dep);
      }
    }
  }

  // Adjacency and in-degrees over the subgraph: every existing dependency
  // edge except the one sv itself currently contributes (the new link
  // supersedes it), plus the proposed edge sv -> vv.
  std::map<void*, std::vector<void*> > adjacency;
  std::map<void*, int> in_degree;
  for (const Vertex& n : nodes) {
    adjacency[n.id] = std::vector<void*>();
    in_degree[n.id] = 0;
  }
  auto add_edge = [&](void* a, void* b) {
    if (adjacency.find(b) == adjacency.end()) return;
    adjacency[a].push_back(b);
    ++in_degree[b];
  };
  for (const Vertex& n : nodes) {
    if (n.node) {
      for (const auto& kv : n.node->get_input_ports()) {
        const Port* p = kv.second.get();
        if (p == this || !p->is_linked()) continue;
        add_edge(n.id, vertex_of(p->get_link()).id);
      }
    } else {
      if (n.port.get() == this || !n.port->is_linked()) continue;
      add_edge(n.id, vertex_of(n.port->get_link()).id);
    }
  }
  add_edge(sv.id, vv.id);

  // Kahn's algorithm: repeatedly peel off zero-in-degree vertices.
  std::vector<void*> queue;
  for (const auto& kv : in_degree) {
    if (kv.second == 0) queue.push_back(kv.first);
  }
  std::size_t removed = 0;
  while (!queue.empty()) {
    void* k = queue.back();
    queue.pop_back();
    ++removed;
    for (void* m : adjacency[k]) {
      if (--in_degree[m] == 0) queue.push_back(m);
    }
  }
  // Leftover vertices keep an incoming edge: the subgraph is cyclic.
  return removed != adjacency.size();
}

std::shared_ptr<Node> Port::get_node() const { return node_.lock(); }

void Port::set_node(std::shared_ptr<Node> n) { node_ = n; }

//! Tell the owning node its cached execution plan is stale.
//!
//! Only linking does this. Value writes go through update_attached_node(),
//! which runs on every propagation and must stay cheap -- invalidating the
//! plan there would rebuild it millions of times during sampling.
void Port::touch_attached_node_structure() {
  if (std::shared_ptr<Node> n = node_.lock()) n->touch_structure();
}

void Port::update_attached_node() {
  if (std::shared_ptr<Node> n = node_.lock()) {
    n->set_valid(false);
    if (is_reactive_ && !is_output_) n->evaluate();
  }
}

void Port::propagate_to_followers() {
  // Overwhelmingly the common case during sampling, and it used to pay for
  // a heap-allocated snapshot before discovering it had nobody to give it to.
  if (linked_to_.empty()) return;
  // chinet pushes its own (post-clip) array into every follower's setter --
  // which flips a scalar follower to vector, quirk and all. Reproduced.
  //
  // Typed, so that a wide integer reaches the follower as the integer it is:
  // pushing the double view would round it on the way, which is the same
  // mistake the old mirror made in three other places.
  if (port_value_type_element(value_type_) == PORT_FLOAT) {
    for (const auto& w : linked_to_) {
      if (std::shared_ptr<Port> f = w.lock()) f->set_value_vector(buffer_);
    }
    return;
  }
  std::vector<long long> exact(buffer_.size());
  for (std::size_t i = 0; i < buffer_.size(); ++i) {
    exact[i] = slot_as_int(buffer_[i]);
  }
  for (const auto& w : linked_to_) {
    if (std::shared_ptr<Port> f = w.lock()) f->set_value_vector_int(exact);
  }
}

std::string Port::describe() const {
  std::ostringstream out;
  out << "Port(name='" << get_name() << "', uid='" << get_uid()
      << "', value_type=" << value_type_ << ", value=[";
  const std::vector<double>& d = get_values_ref();
  for (std::size_t i = 0; i < d.size(); ++i) {
    if (i) out << ", ";
    out << d[i];
  }
  out << "]";
  if (fixed_) out << ", fixed";
  if (is_output_) out << ", output";
  if (is_reactive_) out << ", reactive";
  if (is_bounded_) out << ", bounded=[" << lb_ << ", " << ub_ << "]";
  if (link_) out << ", linked";
  out << ")";
  return out.str();
}

IMPBFF_END_NAMESPACE
