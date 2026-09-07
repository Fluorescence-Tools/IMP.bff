/**
 *  \file IMP/bff/Port.h
 *  \brief A reactive value cell: chinet's Port, standalone.
 *
 *  A Port is the unit of state in the chinet computation model chisurf is
 *  built on: a scalar or an array, a fixed flag, optional hard bounds, an
 *  optional prior specification, and a link to another port it follows.
 *  Writing a value marks the port's node invalid and, when the port is
 *  reactive, re-evaluates it; the write also propagates to every port that
 *  follows this one through a link. The link graph must stay acyclic --
 *  linking is the only operation that adds dependency edges, so the DAG
 *  invariant is enforced there, with Kahn's algorithm (LinkCycleError).
 *
 *  Ported from chinet's chinet/port.py (phase 1 of removing chinet from
 *  chisurf: bff absorbs the parameter/node runtime). This layer is
 *  deliberately standalone, like FactorGraph: no IMP particles, restraints
 *  or decorators take part. Persistence is phase 2: Session.h reads and
 *  writes the chinet document from this state (the document dict itself,
 *  chinet's schema.py and db.py/MMFDB backend are not ported; BaseObject
 *  carries the precursor/death document fields so identity follows the
 *  object, as chinet's document did).
 *
 *  Where C++ and numpy differ the divergence is noted in the member
 *  documentation; the two that matter are
 *
 *  - value data is stored as double. chinet stores float64 or int64 numpy
 *    arrays and infers the dtype from every write; here the type code
 *    (get_value_type(), a PortValueType) carries the element type and
 *    integral ports round integral values. An int vector keeps its type
 *    through set_value_vector_int(). The int and bool element types are
 *    stored exactly, in an int64 vector, and `data_` is a double mirror kept
 *    for the zero-copy readers -- so an integer port round-trips to
 *    INT64_MAX, not only to 2^53.
 *  - bool is bff's, not chinet's: PORT_BOOL and PORT_BOOL_VECTOR are codes
 *    4 and 5, so a chinet document (which only ever holds 0-3) still reads
 *    and writes unchanged. Unlike int/float, bool does not promote -- see
 *    PortValueType for why a declared type must not be inferred away.
 *  - bounds report (nan, nan) when enforcement is off, standing in for
 *    chinet's (None, None).
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_PORT_H
#define IMPBFF_PORT_H

#include <IMP/bff/bff_config.h>

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

class Node;

//! Raised when linking two ports would create a cycle in the link graph.
/*!
    The port graph must remain a directed acyclic graph: a node may only
    depend on nodes that do not transitively depend on it. Subclasses
    std::domain_error -- not std::runtime_error -- because IMP's wrapper
    exception handler maps domain_error to IMP.ValueException (a Python
    ValueError), preserving chinet's LinkCycleError(ValueError) catch
    contract; runtime_error would cross as a plain RuntimeError.
*/
class IMPBFFEXPORT LinkCycleError : public std::domain_error {
 public:
  explicit LinkCycleError(const std::string& what_arg)
      : std::domain_error(what_arg) {}
};

//! The identity half of chinet's BaseObject: a name and a uid.
/*!
    chinet's BaseObject also carries a document dict and database
    registration; those are not ported -- Session writes documents from
    the objects' own state, and there is no registry. What persistence
    needs beyond the runtime's name and uid is the precursor and death
    document fields, which live here so they follow the object across
    sessions exactly as chinet's document did.
*/
class IMPBFFEXPORT BaseObject {
 public:
  explicit BaseObject(const std::string& name = "");
  virtual ~BaseObject();

  //! Human-readable name (ports take the key they are added under).
  const std::string& get_name() const;
  //! Set the name.
  void set_name(const std::string& name);
  //! Unique id, generated at construction, empty-string patchable.
  const std::string& get_uid() const;
  //! Overwrite the uid (as chinet's oid setter does).
  void set_uid(const std::string& uid);
  //! The precursor id: chinet's oid_precursor, the document this object
  //! was derived from. Defaults to the construction uid and, as in
  //! chinet, does not follow a later set_uid().
  const std::string& get_precursor() const;
  //! Overwrite the precursor id (what loading a document restores).
  void set_precursor(const std::string& precursor);
  //! Time of death, chinet's death document field; 0 while alive.
  int get_death() const;
  //! Overwrite the time of death (documents round-trip it verbatim).
  void set_death(int death);

 private:
  std::string name_;
  std::string uid_;
  //! chinet's oid_precursor / time_of_death: document identity.
  std::string precursor_;
  int death_ = 0;
};

//! The type of the values a Port holds: element type, plus vector-ness.
/*!
    The numbers are **not** free to change: Session.h writes the code into the
    chinet document verbatim and `test/session/chinet_fixture.jsonl` pins it, so
    0-3 mean what chinet meant by them. Bool is a bff extension and takes fresh
    codes rather than renumbering.

    Element type and vector-ness are tangled in one code because chinet tangled
    them; use port_value_type_is_vector() and port_value_type_element() rather
    than comparing codes by hand, and port_value_type_of() to build one.

    **The element type is declared once and never changes.** It is set by the
    constructor, or by set_value_type(), and nothing else moves it -- a write
    of another kind is *accepted* and coerced, not refused and not promoted.
    Writing 1.5 to an integer port stores 1; writing 3 to a float port stores
    3.0; writing 3.7 to a bool port stores `true`.

    This is a deliberate divergence from chinet, which inferred the dtype from
    every write (numpy's promotion rules), so an integer port became a float
    port the first time anything stored 1.5 in it. A port whose type can
    change is not typed: whatever holds a reference to it has already decided
    what it is, and one stray write should not be able to answer differently
    on their behalf. It also cost exactness -- a single float write turned an
    exact int64 port into a rounding one for the rest of the session.

    Vector-ness is **not** part of the declaration: it is shape, and it
    follows the data. Writing an array to a scalar port gives an array port of
    the same element type.

*/
enum PortValueType {
  PORT_INT = 0,           //!< scalar integer (chinet int64)
  PORT_FLOAT = 1,         //!< scalar double (chinet float64)
  PORT_INT_VECTOR = 2,    //!< integer array
  PORT_FLOAT_VECTOR = 3,  //!< double array
  PORT_BOOL = 4,          //!< scalar boolean; stored as 0.0 / 1.0
  PORT_BOOL_VECTOR = 5    //!< boolean array; stored as 0.0 / 1.0
};

//! The element type a code carries, as the scalar code for it.
/*! PORT_INT, PORT_FLOAT or PORT_BOOL, whether or not \p t is a vector code. */
IMPBFFEXPORT int port_value_type_element(int t);
//! Whether a code is one of the vector codes.
IMPBFFEXPORT bool port_value_type_is_vector(int t);
//! The code for an element type at a given vector-ness.
/*! \param[in] element PORT_INT, PORT_FLOAT or PORT_BOOL
    \param[in] is_vector whether the array form is wanted */
IMPBFFEXPORT int port_value_type_of(int element, bool is_vector);
//! A readable name for a code ("int", "float[]", ...), for messages and repr.
IMPBFFEXPORT std::string port_value_type_name(int t);

//! A value cell with links, bounds and invalidation-driven evaluation.
/*!
    Ports are shared_ptr-owned: a Node holds its ports, a follower holds the
    port it follows, and SWIG hands Python a shared_ptr from every
    constructor. A port that takes part in a link must be shared_ptr-owned
    (std::enable_shared_from_this); a stack port can hold and read a value
    but cannot itself follow another port.

    Two ports of the same node may be linked (a node may read its own
    output); a standalone port linked to itself is a degenerate cycle and is
    rejected.
*/
class IMPBFFEXPORT Port : public BaseObject,
                          public std::enable_shared_from_this<Port> {
 public:
  //! An empty scalar integer port with value 0 (chinet's Port()).
  Port();
  //! A scalar port; value_type 0 by default, promoted to 1 for a double.
  /*!
      \param[in] value initial scalar value
      \param[in] fixed freeze the port against value writes
      \param[in] is_output mark as a node output port
      \param[in] is_reactive re-evaluate the attached node on writes
      \param[in] is_bounded enforce lb/ub on every write
      \param[in] lb lower bound
      \param[in] ub upper bound
      \param[in] value_type 0 int, 1 float, 2 int vector, 3 float vector
      \param[in] name port name
  */
  Port(double value, bool fixed = false, bool is_output = false,
       bool is_reactive = false, bool is_bounded = false, double lb = 0.0,
       double ub = 0.0, int value_type = 0, const std::string& name = "");
  //! The int twin of the scalar constructor: keeps the integer type code.
  /*! `long long`, not `int`: a Python integer wider than 32 bits does not
      convert to `int`, so SWIG fell through to the `double` overload and the
      port came back **typed float, with the value rounded** --
      `Port(value=2**53+1)` reported `float` and lost a digit. */
  Port(long long value, bool fixed = false, bool is_output = false,
       bool is_reactive = false, bool is_bounded = false, double lb = 0.0,
       double ub = 0.0, int value_type = 0, const std::string& name = "");
  //! The bool twin of the scalar constructor: a declared flag port.
  /*! The value is stored as 0.0 or 1.0 and the port is PORT_BOOL. Note that
      in Python `True` is an `int`, so the binding dispatches on
      `isinstance(v, bool)` before it looks for an integer. */
  Port(bool value, bool fixed = false, bool is_output = false,
       bool is_reactive = false, bool is_bounded = false, double lb = 0.0,
       double ub = 0.0, int value_type = PORT_BOOL,
       const std::string& name = "");

  //! A vector port (float-typed; set_value_vector_int() keeps int).
  /*! There is deliberately **no** `std::vector<int>` constructor overload.
      A Python list of floats converts to `std::vector<int>` as happily as to
      `std::vector<double>` -- lossily, by truncation -- so an overload here
      makes SWIG's dispatch decide the element type, and it decided wrong:
      `Port(value=[1.5, 2.5])` came back as `[1, 2]`. The `value_type`
      argument already says what the caller means, so the overload bought
      nothing and cost that. */
  Port(const std::vector<double>& values, bool fixed = false,
       bool is_output = false, bool is_reactive = false,
       bool is_bounded = false, double lb = 0.0, double ub = 0.0,
       int value_type = 0, const std::string& name = "");
  virtual ~Port();

  //! Scalar value; follows the link if linked (chinet's .value).
  /*!
      Returns element 0 of the underlying data. Whether a port reads as a
      scalar or a vector is the port's own vectorness flag
      (get_is_vector()), exactly as in chinet: the flag follows writes, not
      the link.
  */
  double get_value() const;
  //! Write a scalar float value (no-op when fixed); marks the node invalid.
  void set_value(double v);
  //! Write a scalar int value: keeps an integral port integral, and exact.
  /*! Widened from `int` for the reason on the integer constructor. An
      integral port stores this in int64, so it round-trips beyond 2^53. */
  void set_value(long long v);
  //! The scalar value as an exact integer (the int64 store, when there is one).
  long long get_value_int() const;
  //! Write a scalar bool value: makes the port PORT_BOOL and stores 0.0/1.0.
  void set_value_bool(bool v);
  //! Read the scalar value as a truth value (`!= 0`), whatever the type.
  bool get_value_bool() const;
  //! Full value (link-aware), C++ side of the value property.
  std::vector<double> get_value_vector() const;

  //! The values where they lie, with no copy. **Not SWIG-wrapped.**
  /*!
      `get_value_vector()` returns by value, which is right for a Python
      caller and wrong for the inside of a fit: an `Expression` reading its
      axis, a `ChiSquared` reading the model curve and a `Minimizer` reading
      the residuals each copied a full-length vector on **every** residual
      evaluation. At 512 points that was several kilobytes of memcpy per
      iteration, for buffers none of them writes to.

      The reference follows the link, exactly as `get_value_vector()` does,
      and is valid until the port (or its link's) storage is next written.
  */
  const std::vector<double>& get_values_ref() const;

  //! Whether writes replace non-finite values (chinet's behaviour, default).
  /*!
      chinet floors a NaN to `np.finfo(float).tiny` and clamps infinities to
      the largest double, so a stored value is always JSON-serialisable.
      That is right for a document and **wrong for the numeric transport of
      a fit**: a NaN there means "this trial point is impossible", and
      floored to `tiny` it reads as *zero* -- so a model that blows up looks
      like a model that fits perfectly, and the optimiser is attracted to it
      instead of rejecting it. Measured: a NaN model curve came out of
      `ChiSquared` as chi2 = 14 rather than infinity, where the numpy path
      returns NaN and MINPACK rejects the step.

      So the fit's own ports -- an `Expression`'s curve, a `ChiSquared`'s
      model input and residual output, a `JointChiSquared`'s blocks -- turn
      it off. Everything else keeps chinet's behaviour, which is what a
      saved document still needs.
  */
  void set_sanitize(bool v) { sanitize_ = v; }
  bool get_sanitize() const { return sanitize_; }
  //! Write a vector value (no-op when fixed); marks the node invalid.
  void set_value_vector(const std::vector<double>& v);
  //! Write an integer vector: the port becomes PORT_INT_VECTOR, not float.
  /*! Named rather than overloaded, for the dispatch reason on the vector
      constructor above. This is what closes the int-vector divergence from
      chinet the file comment records, and the values are stored exactly. */
  void set_value_vector_int(const std::vector<long long>& v);
  //! The values as exact integers.
  std::vector<long long> get_value_vector_int() const;
  //! The values as truth values, one per element.
  std::vector<int> get_value_vector_bool() const;

  //! Set a vector value from a numpy array, with no Python list in between.
  /*! \param in_values the values
      \param n_values how many

      `set_value_vector` takes a `std::vector<double>`, which from Python
      means building a list of every element before any of it is stored. On
      a per-iteration path that dominates: a Python objective handing back a
      512-point residual spent more time converting it than the optimiser
      spent on the whole step, and the director objective it serves measured
      *slower* than the scipy loop it replaced because of it. */
  void set_values_array(double* in_values, int n_values);
  //! The value as a fresh numpy array (the SWIG-facing read; managed).
  /*!
      \param[out] out_values a malloc'd buffer of the value; numpy.i's
                  ARGOUTVIEWM typemap passes it to a numpy array that
                  free()s it, so it is C-allocated on purpose
      \param[out] n_out_values its length
  */
  void get_value_view(double** out_values, int* n_out_values) const;

  //! Type code; a PortValueType.
  int get_value_type() const;
  //! Set the type code, converting the stored data to match.
  /*! int codes truncate (numpy's `astype(int64)`), bool codes take `!= 0`,
      float codes leave the data alone. This is the only way out of a bool
      port -- see PortValueType on why writing does not do it. */
  void set_value_type(int t);
  //! The element type: PORT_INT, PORT_FLOAT or PORT_BOOL.
  int get_element_type() const { return port_value_type_element(get_value_type()); }
  //! Whether the port reads as a vector (follows writes, as in chinet).
  bool get_is_vector() const;
  //! Number of stored elements.
  unsigned int current_size() const;
  //! Ports are always self-consistent: true (chinet compatibility).
  bool is_valid() const;

  //! Fixed flag: a fixed port ignores value writes.
  bool get_fixed() const;
  //! Set the fixed flag.
  void set_fixed(bool v);
  //! Output-port flag.
  bool get_is_output() const;
  //! Set the output-port flag (chinet's set_port_type).
  void set_port_type(bool v);
  //! Reactive flag: a reactive input re-evaluates its node when written.
  bool get_is_reactive() const;
  //! Set the reactive flag.
  void set_is_reactive(bool v);

  //! Whether lb/ub are enforced on writes.
  bool get_is_bounded() const;
  //! Enable/disable enforcement; enabling clips the current data.
  void set_is_bounded(bool v);
  //! Stored lower bound (reported whether or not enforcement is on).
  double get_lower_bound() const;
  //! Stored upper bound.
  double get_upper_bound() const;
  //! Store both bounds; clips the current data when enforcement is on.
  void set_bounds(double lb, double ub);

  //! Prior specification as JSON ("" for none; chinet stores a dict).
  /*!
      chisurf keeps its Prior objects on its side and mirrors a
      JSON-serialisable state dict onto the port; a raw string is the same
      contract without a JSON dependency in the runtime. Uniform (box)
      priors are the bounds and stay there.
  */
  const std::string& get_prior() const;
  //! Set the prior specification, or clear it with "".
  void set_prior(const std::string& json);

  //! The port this one follows, or a null pointer.
  std::shared_ptr<Port> get_link() const;

  //! The link without a refcount bump, for the evaluation hot path.
  const std::shared_ptr<Port>& get_link_ref() const { return link_; }
  //! Follow v, or unlink when v is null. Throws LinkCycleError on a cycle.
  void set_link(std::shared_ptr<Port> v);
  //! Break the link; marks the attached node invalid. Always true.
  bool unlink();
  //! Whether this port follows another.
  bool is_linked() const;
  //! Whether set_link(v) would break the DAG invariant (side-effect free).
  bool would_create_cycle(const std::shared_ptr<Port>& v);

  //! The node this port belongs to, or a null pointer.
  std::shared_ptr<Node> get_node() const;
  //! Attach the port to a node (nodes call this; shared_ptr-owned).
  void set_node(std::shared_ptr<Node> n);

  //! Mark the attached node invalid; evaluate it if this port is reactive.
  void update_attached_node();

  //! Invalidate the owning node's cached plan (link changes only).
  void touch_attached_node_structure();

  //! A one-line summary: name, type code, value, flags.
  std::string describe() const;

 private:
  //! Whether this port's declared element type is bool (writes coerce).
  bool is_bool_typed() const;
  //! One slot reinterpreted as the integer it holds, and the reverse.
  static long long slot_as_int(double slot);
  static double int_as_slot(long long v);
  //! The i-th element as a double, whatever the element type is.
  double element_as_double(std::size_t i) const;
  //! Replace the buffer, reading the input as doubles / as exact integers.
  void store_doubles(const std::vector<double>& v);
  void store_ints(const std::vector<long long>& v);
  //! Clip if bounded; always drops the derived double view.
  void finalize_storage();
  //! chinet's general write path, shared by every value entry point.
  void write_value(double v, bool input_is_float);
  //! The vector twin of write_value().
  void write_vector(const std::vector<double>& v, bool input_is_float);
  //! Clip data_ into [lb, ub] (np.clip on write, when bounded).
  void clip_to_bounds();
  //! Push the current data into every port following this one.
  void propagate_to_followers();

  //! **The** store: one 8-byte slot per element, read according to the
  //! element type in value_type_. A float slot holds the double itself; an
  //! int or bool slot holds the int64 bit pattern, copied in and out with
  //! memcpy -- so there is no strict-aliasing game and no alignment
  //! question, a double slot being already 8-byte aligned and exactly as
  //! wide as an int64.
  /*!
      One buffer and a tag, rather than a typed vector per type. The previous
      design kept a double array *and* an int64 array and had to answer "which
      one is right" on every path; it answered wrong three times (the type
      conversion, the Python constructor and the session writer each rounded a
      wide integer back through the double). A single store cannot disagree
      with itself, and a fourth element type is now a tag rather than a fourth
      array to keep in step.
  */
  std::vector<double> buffer_;
  //! A double view of a non-float port, built on demand and dropped on every
  //! write. Strictly derived -- nothing reads it back into buffer_ -- so it
  //! cannot become a second authority. Float ports never use it: for them
  //! get_values_ref() hands out buffer_ itself, which is what keeps the
  //! evaluation hot path a zero-copy read.
  mutable std::vector<double> double_cache_;
  mutable bool double_cache_valid_ = false;
  bool sanitize_ = true;
  int value_type_ = 0;
  bool is_vector_ = false;
  bool fixed_ = false;
  bool is_output_ = false;
  bool is_reactive_ = false;
  bool is_bounded_ = false;
  double lb_ = 0.0;
  double ub_ = 0.0;
  std::string prior_;
  std::shared_ptr<Port> link_;
  //! Ports following this one; weak, so a follower that dies just drops out
  //! (chinet's strong list keeps the follower alive instead -- a reference
  //! cycle only Python's GC could break).
  std::vector<std::weak_ptr<Port> > linked_to_;
  std::weak_ptr<Node> node_;
};

IMPBFF_END_NAMESPACE

#endif // IMPBFF_PORT_H
