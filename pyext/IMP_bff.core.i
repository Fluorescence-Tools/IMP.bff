/*
 * The core's Python surface, in the order SWIG needs it: every topic file and
 * loose header that names nothing of IMP's particle world. The IMP-module
 * build includes this and then layer.i (the connection layer); the standalone
 * build includes this alone, after its own preamble (standalone/pyext).
 * PRD-137 step 6c.
 */

/* Vendored from numpy upstream; not ours, not edited. */
%include "numpy.i"

%pythoncode %{
import json
import numpy as np
%}

/* Does this build actually thread? The kernels ask; the compiler may not listen. */
%include "IMP/bff/Parallel.h"

/* STL templates and the numpy typemaps everything below depends on. */
%include "IMP_bff.types.i"

/* Where the shipped cgprobe data lives. */
%include "IMP/bff/DataPaths.h"

/*
 * The factor structure of a model's posterior: variables (free parameters),
 * factors (per-dataset likelihoods and per-variable priors), and the graph
 * queries a structured optimiser needs -- elimination order, cliques,
 * junction tree, treewidth, separators, sampling blocks and relevance.
 * Ported from ChiSurf's PRD-68 factor graph. Deliberately standalone: no
 * IMP particles, restraints or decorators take part in the graph.
 */
%template(VectorString) std::vector<std::string >;
%template(VectorInt) std::vector<int >;
/* Port's exact integer store: a Python int is arbitrary precision and a
   C++ `int` is 32 bits, so an integer port has to speak int64 or large
   values fall through to the double overload and come back rounded. */
%template(VectorInt64) std::vector<long long >;
%include "IMP/bff/FactorGraph.h"
%template(VectorJunctionTreeEdge) std::vector<IMP::bff::JunctionTreeEdge>;

/*
 * The reactive Port/Node runtime ported from chinet (phase 1 of removing
 * chinet from chisurf): value cells with links, bounds and invalidation-
 * driven evaluation, and nodes computing outputs from inputs over those
 * ports. Deliberately standalone, like the factor graph above -- no IMP
 * particles, restraints or decorators take part. Persistence is phase 2:
 * Session (below) reads and writes chinet's session format; chinet's
 * schema.py and db.py/MMFDB backend stay out of bff entirely.
 *
 * Ports and nodes are shared_ptr-owned: a node holds its ports, a follower
 * holds the port it follows, and every Python constructor hands back a
 * shared_ptr. get_value_vector() stays C++-only -- this module cannot name
 * std::vector<double> (see IMP_bff.types.i) -- so the Python surface reads
 * values through get_value_view(), a managed numpy array, and through the
 * `value` property below. Python callables as node callbacks are chinet's
 * callback_class: the SWIG director below lets a Python subclass override
 * evaluate()/update(), and the set_python_callback_function ergonomics
 * (signature introspection, auto ports) are a chisurf-side helper
 * (chisurf/core/nodes.py) -- bff stays C++-minimal.
 */
/* LinkCycleError crosses to Python as IMP.ValueException (a ValueError):
   IMP's handle_imp_exception maps std::domain_error there, which is why
   the C++ class derives domain_error and not runtime_error -- chinet's
   error subclasses ValueError and its callers catch it as one. */

%shared_ptr(IMP::bff::BaseObject);
%shared_ptr(IMP::bff::Port);
%shared_ptr(IMP::bff::Node);
/* A Python subclass overrides evaluate()/update() to run a Python callable
   over the port maps -- chinet's callback_class. The chinet ergonomics of
   set_python_callback_function (inspect the signature, auto-create the
   ports) live in chisurf's core/nodes.py, not here. Directors and
   %shared_ptr cooperate through SWIG's shared_ptr director support.

   The director lifetime (T-20260901-13): a C++ Minimizer or Sampler holds
   a Python Node subclass across the whole run through a shared_ptr, but
   the director keeps only a *weak* pointer back to the Python proxy -- a
   ResidualNode bound to `_` was collected when the next tuple unpacking
   rebound `_`, and Node_update then dispatched into a dead object.
   MinimizerObserver's answer (IMP_SWIG_DIRECTOR) does NOT work here:
   `_director_objects.register` silently refuses anything without IMP's
   `get_ref_count`, and Node is a plain shared_ptr class, not an
   IMP::Object -- the macro would read as protection and protect nothing.
   The real fix mirrors the C++ ownership on the Python side: the
   %pythonappend hooks below set_objective (Minimizer, Sampler) stash the
   node proxy on the wrapper that holds the shared_ptr, so the proxy lives
   exactly as long as the C++ reference does. */
%feature("director") IMP::bff::Node;
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** out_values, int* n_out_values)};
%ignore IMP::bff::Port::get_value_vector;
%ignore IMP::bff::Node::set_callback_function;
/* Zero-copy in for a port's vector value: a Python objective writes its
   residual here once per iteration, and a std::vector<double> conversion
   would build a list of every point first. */
%apply(double* IN_ARRAY1, int DIM1) {(double* in_values, int n_values)};
%ignore IMP::bff::Port::get_values_ref;
%include "IMP/bff/Port.h"
%include "IMP/bff/Node.h"
%template(MapStringPort) std::map<std::string, std::shared_ptr<IMP::bff::Port> >;

/* The chinet surface on top of the accessors: scalar-or-array values that
   follow the port's vectorness flag, dict priors, tuple bounds with NaN
   standing in for chinet's None. */

/* The properties name their accessors through self, not as bare functions:
   SWIG emits an %extend %pythoncode block into the proxy class *before* the
   method definitions, so a bare name is a NameError at class-creation time
   (and no line inside %pythoncode may start with '#' -- SWIG's preprocessor
   reads it as a directive). */
%extend IMP::bff::Port {
    %pythoncode {
        def _set_value(self, v):
            """Write, dispatching on the value's type and the port's.

            `isinstance(v, bool)` is tested FIRST because in Python `bool` is
            a subclass of `int`: `hasattr(True, '__index__')` is true, so the
            integer branch would swallow every flag. A bool port also keeps
            its type when it is handed a number -- see PortValueType -- so
            the int/float branches below are safe on one.
            """
            if isinstance(v, (bool, np.bool_)):
                self.set_value_bool(bool(v))
            elif hasattr(v, '__len__') and not isinstance(v, (str, bytes)):
                arr = np.asarray(v).ravel()
                if arr.dtype.kind in 'iub' and self.get_element_type() != PORT_FLOAT:
                    self.set_value_vector_int([int(x) for x in arr])
                else:
                    self.set_value_vector(arr.astype(np.float64))
            elif hasattr(v, '__index__'):
                self.set_value(int(v))
            else:
                self.set_value(float(v))

        def _get_value(self):
            """Read back in the port's own type, not always as float.

            A typed port that answers every read as a double is only half a
            type: `Port(value=True).value` used to be `1.0`. Scalars come back
            as `bool`/`int`/`float` and arrays with the matching dtype.
            """
            t = self.get_value_type()
            element = self.get_element_type()
            if not self.get_is_vector() and len(self.get_value_view()) == 1:
                if element == PORT_BOOL:
                    return self.get_value_bool()
                if element == PORT_INT:
                    return self.get_value_int()
                return self.get_value()
            values = self.get_value_view()
            if element == PORT_BOOL:
                return values.astype(bool)
            if element == PORT_INT:
                return np.asarray(self.get_value_vector_int(), dtype=np.int64)
            return values
        value = property(_get_value, _set_value)
        """The bare `PORT_*` names above resolve at CALL time, in the module
        globals SWIG emits the enum into -- which is why they are allowed here
        while the note at the top of this block forbids bare names in the
        class-body expressions (those run at class-creation time, before the
        module is finished)."""
        value_type = property(lambda self: self.get_value_type(),
                              lambda self, v: self.set_value_type(int(v)))
        element_type = property(lambda self: self.get_element_type())
        type_name = property(lambda self: port_value_type_name(
            self.get_value_type()))
        fixed = property(lambda self: self.get_fixed(),
                         lambda self, v: self.set_fixed(v))
        is_output = property(lambda self: self.get_is_output(),
                             lambda self, v: self.set_port_type(v))
        is_reactive = property(lambda self: self.get_is_reactive(),
                               lambda self, v: self.set_is_reactive(v))
        reactive = property(lambda self: self.get_is_reactive(),
                            lambda self, v: self.set_is_reactive(v))
        bounded = property(lambda self: self.get_is_bounded(),
                           lambda self, v: self.set_is_bounded(v))
        is_vector = property(lambda self: self.get_is_vector())
        link = property(lambda self: self.get_link(),
                        lambda self, v: self.set_link(v))
        name = property(lambda self: self.get_name(),
                        lambda self, v: self.set_name(v))
        uid = property(lambda self: self.get_uid(),
                       lambda self, v: self.set_uid(v))
        precursor = property(lambda self: self.get_precursor(),
                             lambda self, v: self.set_precursor(v))
        bounds = property(
            lambda self: ((self.get_lower_bound(), self.get_upper_bound())
                          if self.get_is_bounded()
                          else (float('nan'), float('nan'))),
            lambda self, v: self.set_bounds(float(v[0]), float(v[1])))
        prior = property(
            lambda self: (json.loads(self.get_prior())
                          if self.get_prior() else None),
            lambda self, v: self.set_prior(
                json.dumps(v) if v is not None else ""))

        def _document_value(self):
            v = self.value
            if isinstance(v, np.ndarray):
                v = v.tolist()
            element = self.get_element_type()
            if element == PORT_BOOL:
                v = [bool(x) for x in v] if isinstance(v, list) else bool(v)
            elif element == PORT_INT:
                v = [int(x) for x in v] if isinstance(v, list) else int(v)
            return v

        def get_json(self, indent=0):
            """The port's document as JSON, chinet's Port.get_json.

            The fields chinet's _update_doc_from_data wrote (value,
            value_type, bounds, the flags, the prior, the link target), so
            what read a chinet document reads this one; chisurf's
            Parameter pickle path is the consumer that matters. Stored
            bounds are written whether or not enforcement is on, as chinet
            wrote self._bounds, and a link is its target's uid.
            """
            link = self.get_link()
            doc = {
                "type": "port",
                "name": self.get_name(),
                "fixed": self.get_fixed(),
                "is_output": self.get_is_output(),
                "is_reactive": self.get_is_reactive(),
                "is_bounded": self.get_is_bounded(),
                "value": self._document_value(),
                "bounds": [self.get_lower_bound(), self.get_upper_bound()],
                "link": link.get_uid() if link is not None else None,
                "value_type": self.get_value_type(),
                "prior": (json.loads(self.get_prior())
                          if self.get_prior() else None),
            }
            return json.dumps(doc, indent=indent if indent > 0 else None)

        def read_json(self, s):
            """Restore the port from a get_json document, chinet's
            Port.read_json (its set_document).

            Field order mirrors the session loader: value_type FIRST (so the
            value coerces to the right element type), then the value, then
            value_type again (restoring the document's exact code, including
            chinet's quirk that a constructor-built float vector reports the
            scalar code -- a same-element retype converts nothing), then
            bounds, the flags, enforcement (enabling clips
            the value, as on load), the prior. The type leads because a write
            no longer retypes the port -- it coerces -- so restoring the value
            into a default-typed port would round it away before the type
            arrived. The link field is
            written for readers but not restored -- chinet's set_document
            left it alone too; a link needs its target, not a uid string.
            The value is written even when the port is fixed, as chinet's
            set_document wrote _data directly: `fixed` is lifted for the
            write and restored from the document below.
            """
            doc = json.loads(s)
            if "value_type" in doc:
                self.set_value_type(int(doc["value_type"]))
            if "value" in doc:
                was_fixed = self.get_fixed()
                self.set_fixed(False)
                v = doc["value"]
                if isinstance(v, (list, tuple)):
                    self.set_value_vector([float(x) for x in v])
                elif isinstance(v, bool):
                    self.set_value(int(v))
                elif isinstance(v, float):
                    self.set_value(v)
                else:
                    self.set_value(int(v))
                self.set_fixed(was_fixed)
            if "value_type" in doc:
                self.set_value_type(int(doc["value_type"]))
            if "bounds" in doc and isinstance(doc["bounds"], (list, tuple)) \
                    and len(doc["bounds"]) == 2:
                self.set_bounds(float(doc["bounds"][0]),
                                float(doc["bounds"][1]))
            if "fixed" in doc:
                self.set_fixed(bool(doc["fixed"]))
            if "is_output" in doc:
                self.set_port_type(bool(doc["is_output"]))
            if "is_reactive" in doc:
                self.set_is_reactive(bool(doc["is_reactive"]))
            if "is_bounded" in doc:
                self.set_is_bounded(bool(doc["is_bounded"]))
            if "prior" in doc:
                self.prior = doc["prior"] if isinstance(doc["prior"], dict) \
                    else None
            if "name" in doc and doc["name"]:
                self.set_name(str(doc["name"]))
            return True
    }
}

/*
 * chinet's Python constructor surface, which chisurf's Parameter builds its
 * ports through:
 *
 *     Port(value=..., fixed=..., is_output=..., is_reactive=...,
 *          is_bounded=..., lb=..., ub=..., value_type=..., name=...,
 *          prior={...})
 *
 * SWIG cannot generate keyword arguments for an overloaded constructor (it
 * needs every overload keyword-enabled or it falls back to positional), so
 * the keyword form is one honest Python shim over the positional overloads:
 * the value is normalised -- python scalar, numpy scalar, list, tuple or
 * array -- then dispatched to the matching C++ constructor, and the prior
 * dict is serialised onto the port after construction. Positional calls go
 * through untouched.
 *
 * Type rules, as chinet's: a python/numpy float scalar is a float port, an
 * int scalar an int port, and an array-like is vector-valued. chinet kept
 * numpy's int64 for an int array (with its scalar type code); here an
 * integer-kind array is passed value_type 2 -- an int vector by bff's own
 * code conventions -- unless the caller says otherwise.
 */
%pythoncode %{
_Port_init_positional = Port.__init__

def _Port_kwargs_init(self, *args, **kwargs):
    if not kwargs:
        _Port_init_positional(self, *args)
        return
    names = ("value", "fixed", "is_output", "is_reactive", "is_bounded",
             "lb", "ub", "value_type", "name", "prior", "oid")
    unknown = set(kwargs) - set(names)
    if unknown:
        raise TypeError(
            "Port() got an unexpected keyword argument '"
            + "', '".join(sorted(unknown))
            + "'. Accepted: " + ", ".join(names))
    if len(args) > len(names):
        raise TypeError("Port() takes at most " + str(len(names))
                        + " positional arguments")
    for position, argument in enumerate(args):
        key = names[position]
        if key in kwargs:
            raise TypeError(
                "Port() got multiple values for argument '" + key + "'")
        kwargs[key] = argument
    value = kwargs.get("value", 0)
    fixed = bool(kwargs.get("fixed", False))
    is_output = bool(kwargs.get("is_output", False))
    is_reactive = bool(kwargs.get("is_reactive", False))
    is_bounded = bool(kwargs.get("is_bounded", False))
    lb = float(kwargs.get("lb", 0.0))
    ub = float(kwargs.get("ub", 0.0))
    name = str(kwargs.get("name", ""))
    value_type = kwargs.get("value_type", None)
    if value_type is not None:
        value_type = int(value_type)
    want_bool = isinstance(value, (bool, np.bool_))
    if want_bool:
        value = int(value)
        if value_type is None:
            value_type = 4
    if hasattr(value, '__len__') and not isinstance(value, (str, bytes)):
        arr = np.asarray(value).ravel()
        if value_type is None:
            if arr.dtype.kind == 'b':
                value_type = 5
            else:
                value_type = 2 if arr.dtype.kind in 'iu' else 0
        if arr.dtype.kind in 'iub':
            _Port_init_positional(self, [], False, is_output, is_reactive,
                                  is_bounded, lb, ub, value_type, name)
            self.set_value_vector_int([int(x) for x in arr])
            if fixed:
                self.set_fixed(True)
        else:
            _Port_init_positional(self, arr.astype(np.float64).tolist(), fixed,
                                  is_output, is_reactive, is_bounded, lb, ub,
                                  value_type, name)
    else:
        if isinstance(value, (int, np.integer)):
            if value_type is None:
                value_type = 0
            _Port_init_positional(self, int(value), fixed, is_output,
                                  is_reactive, is_bounded, lb, ub,
                                  value_type, name)
        else:
            if value_type is None:
                value_type = 0
            _Port_init_positional(self, float(value), fixed, is_output,
                                  is_reactive, is_bounded, lb, ub,
                                  value_type, name)
    prior = kwargs.get("prior", None)
    if prior is not None:
        if not isinstance(prior, dict):
            raise TypeError("Port(prior=...) takes a dict or None")
        self.set_prior(json.dumps(prior))
    if kwargs.get("oid", None) is not None:
        self.set_uid(str(kwargs["oid"]))

Port.__init__ = _Port_kwargs_init
%}

%extend IMP::bff::Node {
    %pythoncode {
        name = property(lambda self: self.get_name(),
                        lambda self, v: self.set_name(v))
        uid = property(lambda self: self.get_uid(),
                       lambda self, v: self.set_uid(v))
        precursor = property(lambda self: self.get_precursor(),
                             lambda self, v: self.set_precursor(v))
        ports = property(lambda self: self.get_ports())
        inputs = property(lambda self: self.get_input_ports())
        outputs = property(lambda self: self.get_output_ports())
    }
}

/*
 * Session persistence (phase 2): the chinet JSONL session format, field
 * for field, so existing chisurf .csp projects load with no importer and
 * a re-saved file lines up against a chinet-written one. The session is
 * the registry -- there is no hidden global database behind the wrapper,
 * so a port or node reaches a file only through add_port()/add_node(),
 * never by construction alone. chinet's schema.py conversion and db.py /
 * MMFDB backend are not ported; the legacy monolithic
 * {"session": ..., "objects": [...]} format is read, not written.
 */
%shared_ptr(IMP::bff::Session);
%include "IMP/bff/Session.h"
%template(MapStringNode) std::map<std::string, std::shared_ptr<IMP::bff::Node> >;
%template(VectorPort) std::vector<std::shared_ptr<IMP::bff::Port> >;

/* The chinet surface on top of the accessors: dict-of-nodes reads, and
   the document-identity fields as attributes. Properties name their
   accessors through self, as above -- a bare name is a NameError at
   class-creation time. */
%extend IMP::bff::Session {
    %pythoncode {
        nodes = property(lambda self: self.get_nodes())
        name = property(lambda self: self.get_name(),
                        lambda self, v: self.set_name(v))
        uid = property(lambda self: self.get_uid(),
                       lambda self, v: self.set_uid(v))
        precursor = property(lambda self: self.get_precursor(),
                             lambda self, v: self.set_precursor(v))
        death = property(lambda self: self.get_death(),
                         lambda self, v: self.set_death(v))
    }
}

/* chisurf's entry point once chinet is gone: one default session per
   process, made on first use. The name resolves immediately (nothing
   lazy about the attribute itself); only the session inside is made on
   demand, and it stays out of C++ -- a hidden global is not a bff
   object's business. */
%pythoncode %{
_session = None

def get_session():
    """The process's default bff session, made on first call."""
    global _session
    if _session is None:
        _session = Session()
    return _session
%}

/*
 * ChiSurf's samplers (phase 6 of removing chinet from chisurf): the
 * affine-invariant stretch move, differential evolution with snooker
 * updates, and the blocked Metropolis walk of chisurf.core.fitting.sample,
 * ported 1:1 and run entirely in C++ -- the walker is written into the
 * parameter ports, the objective node graph is updated and read, and the
 * accept/reject is decided without a single boundary crossing. What that
 * removes is measured: 1.67 us per port set+get across the wrapper against
 * 0.06 us for a Python attribute, paid once per proposal per parameter.
 *
 * The objective is a node (any node -- including a director subclass whose
 * evaluate() runs Python, which is the existing callback machinery and the
 * way a Python log-posterior plugs in) whose named output port carries chi^2
 * by default; bounds come from the parameter ports that enforce them and
 * priors from the ports' JSON prior specifications, exactly the contract
 * chisurf's lnprior/lnprob_parts implement. The C++ std::function objective
 * and the C++ progress observer are not wrapped: this SWIG carries no
 * std_function.i, and a Python callable is a Node director anyway.
 *
 * get_blocks() is ignored because std::vector<std::vector<int>> is not a
 * template this module may name (see IMP_bff.types.i); get_block_sizes()
 * is the wrapped read of the same partition.
 */
%shared_ptr(IMP::bff::Sampler);
/* Same Node-director lifetime fix as Minimizer::set_objective
   (T-20260901-13): the sampler holds the objective node by shared_ptr for
   its lifetime, so the Python proxy must live as long. */
%pythonappend IMP::bff::Sampler::set_objective %{
        self.__dict__['_objective_node_keepalive'] = (
            args[0] if args else kwargs.get('node'))
%}
%ignore IMP::bff::Sampler::set_objective_function;
%ignore IMP::bff::Sampler::set_observer;
%ignore IMP::bff::Sampler::get_blocks;
%include "IMP/bff/Sampler.h"

/*
 * The data misfit of a model curve, as a node in the model graph. Lets a
 * whole fit -- parameters, model, chi-square -- live in one C++ graph, so
 * Sampler can drive it without crossing into Python per move. Ported from
 * ChiSurf's calculate_weighted_residuals / get_chi2.
 */
/* Numpy in and out for the residual path: it runs once per fit iteration for
   every model, so nothing here may go through a Python list. */
%apply(double* IN_ARRAY1, int DIM1) {(double* in_data_y, int n_data_y)};
%apply(double* IN_ARRAY1, int DIM1) {(double* in_data_ey, int n_data_ey)};
%apply(double* IN_ARRAY1, int DIM1) {(double* in_mask_a, int n_mask_a)};
%apply(double* IN_ARRAY1, int DIM1) {(double* in_model_y, int n_model_y)};
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** out_wres, int* n_out_wres)};
%shared_ptr(IMP::bff::ChiSquared);
/* A dataset of any rank, with its noise family. Its variance and residuals
   are managed views for the same reason the solver's are. */
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** out_view, int* n_out_view)};
%include "IMP/bff/Dataset.h"

%include "IMP/bff/ChiSquared.h"

/*
 * The grouping: one misfit over several datasets. Sharing a parameter
 * between them is already `Port::set_link`; this is the other half, the
 * joint objective, so a `Minimizer` pointed at it moves the shared
 * parameters using every dataset's curvature at once.
 *
 * Deliberately not called a fit group -- ChiSurf's FitGroup is a container
 * with a selection, a history and a run policy, and this is the arithmetic
 * underneath it.
 */
%shared_ptr(IMP::bff::JointChiSquared);
%include "IMP/bff/JointChiSquared.h"

/*
 * A TCSPC decay as a node: the multi-exponential model curve a
 * time-correlated instrument produces, so a lifetime fit joins a parse fit
 * on the graph instead of returning to Python once per iteration.
 *
 * The kernels are tttrlib's, taken header-only from the vendored copy of
 * DecayConvolution.h; this class is the graph around them, exactly as
 * Expression is the graph around tttrlib's expression engine.
 */
%apply(double* IN_ARRAY1, int DIM1) {(double* in_response, int n_response)};
%apply(double* IN_ARRAY1, int DIM1) {(double* in_data_y, int n_data_y)};
%apply(double* IN_ARRAY1, int DIM1) {(double* in_data_ey, int n_data_ey)};
%apply(double* IN_ARRAY1, int DIM1) {(double* in_table, int n_table)};
%shared_ptr(IMP::bff::TcspcDecay);
%include "IMP/bff/TcspcDecay.h"

/*
 * The other half of a decay model: the nodes that *produce* the interleaved
 * spectrum TcspcDecay reconvolves. Kept separate from the instrument on
 * purpose -- every TCSPC model shares the instrument and differs only in the
 * photophysics upstream of it, so folding the physics into TcspcDecay would
 * make that class the union of every model anyone fits, and folding it into
 * the caller would put the deriving back in Python once per iteration.
 */
%apply(double* IN_ARRAY1, int DIM1) {(double* in_axis, int n_axis)};
%shared_ptr(IMP::bff::LifetimeSpectrumNode);
%shared_ptr(IMP::bff::AnisotropySpectrum);
%shared_ptr(IMP::bff::PolymerDistances);
%shared_ptr(IMP::bff::GaussianDistances);
%shared_ptr(IMP::bff::FretSpectrum);
%include "IMP/bff/SpectrumNode.h"

/*
 * A model equation compiled once and evaluated in C++, so ChiSurf's parse
 * models stop paying for a Python eval() on every sampler move.
 */
%apply(double* IN_ARRAY2, int DIM1, int DIM2) {(double* in_columns, int n_vars, int n_rows)};
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** out_values, int* n_out_values)};
%apply(unsigned char** ARGOUTVIEWM_ARRAY1, int* DIM1) {(unsigned char** out_mask, int* n_out_mask)};

/* compute_curve: scalar parameters and one axis, each read where it lies.
   Two separate IN_ARRAY1 pairs, because a fit's parameters and its axis are
   different lengths -- which is the whole point of the entry point. */
%apply(double* IN_ARRAY1, int DIM1) {(double* in_parameters, int n_parameters)};
%apply(double* IN_ARRAY1, int DIM1) {(double* in_axis, int n_axis)};

%shared_ptr(IMP::bff::Expression);
%include "IMP/bff/Expression.h"

/*
 * The excitation and emission crosstalk matrices a FRET correction is stated
 * in, and their algebra. Upstream of the decay nodes: a crosstalk matrix is
 * instrument calibration -- what the light path does to a photon between the
 * laser and the detector -- which the light-path calculator builds a payload
 * of and every corrected model consumes. ChiSurf keeps the calculator and the
 * views; the definition and the algebra are this (see IMP_bff.crosstalk.i).
 */
%include "IMP_bff.crosstalk.i"

/*
 * ChiSurf's bounded Levenberg-Marquardt, over the same Port/Node graph
 * `Sampler` walks. `fit.run()` becomes one crossing instead of one per
 * parameter per iteration; see Minimizer.h for why that, and not the
 * optimiser's own arithmetic, is what the port is for.
 *
 * The observer is an IMP::Object director rather than a std::function, for
 * the reason `RRTCollision` is: the C++ loop has to ask a Python object
 * something, and this module carries no std_function.i. It is owned through
 * an IMP::Pointer because it outlives the call that set it -- a progress
 * dialog held only by a raw pointer would be collected mid-fit.
 *
 * set_residual_function is not wrapped: a Python residual is a Node
 * director, which is the same division Sampler makes.
 */
IMP_SWIG_OBJECT(IMP::bff, MinimizerObserver, MinimizerObservers);
/* IMP_SWIG_DIRECTOR, not a bare %feature("director"): the C++ side holds
   this object across the whole run, and SWIG's director keeps only a *weak*
   pointer back to the Python proxy. A progress dialog passed inline --
   `m.set_observer(MyObserver())` -- therefore lost its proxy to the garbage
   collector and the next report() segfaulted. IMP's macro registers every
   director instance in `_director_objects`, which is the module-wide answer
   to exactly this. (`RRTCollision` has a bare feature because it is a call
   argument and cannot outlive its caller's reference; this one can.) */
IMP_SWIG_DIRECTOR(IMP::bff, MinimizerObserver);
%shared_ptr(IMP::bff::Minimizer);
%ignore IMP::bff::Minimizer::set_residual_function;
/* The Node-director lifetime fix (T-20260901-13; rationale at the Node
   director block above): the C++ side keeps a shared_ptr to the objective
   node for the minimiser's lifetime, so the Python proxy must live as
   long too -- a director whose proxy is collected dispatches into a dead
   object on the next update. Stashing it on the wrapper is the Python
   mirror of the C++ reference. */
%pythonappend IMP::bff::Minimizer::set_objective %{
        self.__dict__['_objective_node_keepalive'] = (
            args[0] if args else kwargs.get('node'))
%}
%include "IMP/bff/Minimizer.h"

%extend IMP::bff::Sampler {
    %pythoncode {
        algorithm = property(lambda self: self.get_algorithm(),
                             lambda self, v: self.set_algorithm(v))
        seed = property(lambda self: self.get_seed(),
                        lambda self, v: self.set_seed(v))
        temp = property(lambda self: self.get_temp(),
                        lambda self, v: self.set_temp(v))
        chi2max = property(lambda self: self.get_chi2max(),
                           lambda self, v: self.set_chi2max(v))
        stretch_scale = property(lambda self: self.get_stretch_scale(),
                                 lambda self, v: self.set_stretch_scale(v))
        jitter = property(lambda self: self.get_jitter(),
                          lambda self, v: self.set_jitter(v))
        snooker = property(lambda self: self.get_snooker(),
                           lambda self, v: self.set_snooker(v))
        step_size = property(lambda self: self.get_step_size(),
                             lambda self, v: self.set_step_size(v))
        n_adapt = property(lambda self: self.get_n_adapt(),
                           lambda self, v: self.set_n_adapt(v))
        live_dangerously = property(
            lambda self: self.get_live_dangerously(),
            lambda self, v: self.set_live_dangerously(v))
        number_of_walkers = property(
            lambda self: self.get_number_of_walkers(),
            lambda self, v: self.set_number_of_walkers(v))
        number_of_chains = property(
            lambda self: self.get_number_of_chains(),
            lambda self, v: self.set_number_of_chains(v))
        walker_start_std = property(
            lambda self: self.get_walker_start_std(),
            lambda self, v: self.set_walker_start_std(v))
        walker_start = property(lambda self: self.get_walker_start(),
                                lambda self, v: self.set_walker_start(v))
        initial_values = property(lambda self: self.get_initial_values(),
                                  lambda self, v: self.set_initial_values(v))
        parameter_ports = property(
            lambda self: self.get_parameter_ports(),
            lambda self, v: self.set_parameter_ports(v))
        parameter_names = property(lambda self: self.get_parameter_names())
        objective = property(lambda self: self.get_objective(),
                             lambda self, v: self.set_objective(v))
        output_is_log_likelihood = property(
            lambda self: self.get_output_is_log_likelihood(),
            lambda self, v: self.set_output_is_log_likelihood(v))
        factor_graph = property(lambda self: self.get_factor_graph(),
                                lambda self, v: self.set_factor_graph(v))
        block_sizes = property(lambda self: self.get_block_sizes())

        chain = property(lambda self: np.asarray(self.get_chain()))
        walkers = property(lambda self: np.asarray(self.get_walkers()))
        log_prob = property(lambda self: np.asarray(self.get_log_prob()))
        lnprior = property(lambda self: np.asarray(self.get_lnprior()))
        chi2 = property(lambda self: np.asarray(self.get_chi2()))
        acceptance_rate = property(lambda self: self.get_acceptance_rate())
        acceptance_fractions = property(
            lambda self: np.asarray(self.get_acceptance_fractions()))
        block_acceptance_rates = property(
            lambda self: np.asarray(self.get_block_acceptance_rates()))
        iteration = property(lambda self: self.get_iteration())
        n_evaluations = property(
            lambda self: self.get_number_of_evaluations())
    }
}

/*
 * The bounds are a pair everywhere in chisurf, so they read as one here too
 * -- `bounds = [(lo, hi), ...]`, which is exactly what `leastsqbound` takes
 * and what `model.parameter_bounds` produces. A `None` in either slot means
 * "no bound in that direction", the spelling chisurf uses; it becomes the
 * infinity the C++ tests for.
 *
 * `get_covariance()` is flat because SWIG may not name a vector-of-vector
 * here (IMP_bff.types.i); the square shape is restored in `covariance`.
 *
 * No `#` comments below: inside `%pythoncode` SWIG reads a line starting
 * with `#` as one of its own preprocessor directives and stops.
 */
%extend IMP::bff::Minimizer {
    %pythoncode {
        algorithm = property(lambda self: self.get_algorithm(),
                             lambda self, v: self.set_algorithm(v))
        parameter_ports = property(
            lambda self: self.get_parameter_ports(),
            lambda self, v: self.set_parameter_ports(v))
        parameter_names = property(lambda self: self.get_parameter_names())
        initial_values = property(lambda self: self.get_initial_values(),
                                  lambda self, v: self.set_initial_values(v))
        objective = property(lambda self: self.get_objective(),
                             lambda self, v: self.set_objective(v))
        residual_port_key = property(
            lambda self: self.get_residual_port_key(),
            lambda self, v: self.set_residual_port_key(v))
        observer = property(lambda self: self.get_observer(),
                            lambda self, v: self.set_observer(v))

        ftol = property(lambda self: self.get_ftol(),
                        lambda self, v: self.set_ftol(v))
        xtol = property(lambda self: self.get_xtol(),
                        lambda self, v: self.set_xtol(v))
        gtol = property(lambda self: self.get_gtol(),
                        lambda self, v: self.set_gtol(v))
        maxfev = property(lambda self: self.get_maxfev(),
                          lambda self, v: self.set_maxfev(v))
        epsfcn = property(lambda self: self.get_epsfcn(),
                          lambda self, v: self.set_epsfcn(v))
        factor = property(lambda self: self.get_factor(),
                          lambda self, v: self.set_factor(v))
        diag = property(lambda self: self.get_diag(),
                        lambda self, v: self.set_diag(v))

        def _get_bounds(self):
            import math
            lo, hi = self.get_lower_bounds(), self.get_upper_bounds()
            return [(None if math.isinf(a) else a,
                     None if math.isinf(b) else b) for a, b in zip(lo, hi)]

        def _set_bounds(self, v):
            lo, hi = [], []
            for a, b in v:
                lo.append(-np.inf if a is None else float(a))
                hi.append(np.inf if b is None else float(b))
            self.set_bounds(lo, hi)

        bounds = property(_get_bounds, _set_bounds)

        x = property(lambda self: np.asarray(self.get_x()))
        residuals = property(lambda self: np.asarray(self.get_residuals()))
        chi2 = property(lambda self: self.get_chi2())
        chi2r = property(lambda self: self.get_chi2r())
        status = property(lambda self: self.get_status())
        cancelled = property(lambda self: self.get_cancelled())
        message = property(lambda self: self.get_message())
        n_evaluations = property(
            lambda self: self.get_number_of_evaluations())
        errors = property(lambda self: np.asarray(self.get_errors()))

        def _covariance(self):
            c = np.asarray(self.get_covariance())
            if c.size == 0:
                return c.reshape(0, 0)
            n = len(self.get_x())
            return c.reshape(n, n)

        covariance = property(_covariance)

        def finite_difference_covariance(self, epsilon=0.0, floor=0.0):
            """The covariance differenced afresh at the solution, in C++.

            Returns ``(cov, used)`` -- a ``k x k`` matrix and the indices of
            the ``k`` parameters it is over. ``k`` is short of the parameter
            count when a parameter does not move the objective, which is
            chisurf's ``important_parameters`` and is why the two come back
            together rather than the matrix alone.

            ``epsilon = 0`` means ``sqrt(machine eps)`` and ``floor = 0``
            means 1.0, i.e. chisurf's ``approx_grad`` step. See
            Minimizer.h for why this is not the optimiser's ``epsfcn``.
            """
            c = np.asarray(self.compute_covariance(epsilon, floor))
            used = [int(k) for k in self.get_covariance_parameters()]
            k = len(used)
            return c.reshape(k, k), used

        def finite_difference_covariance_at(self, x, epsilon=0.0, floor=0.0):
            """The same, at a point of the caller's choosing.

            Takes ``p + 2`` evaluations rather than ``p + 1`` -- the
            residuals at ``x`` cannot be assumed the way they can straight
            after a run. This is the entry point for a curvature wanted
            *without* a fit: a posterior view, an error propagated onto a
            derived quantity, a sampler's preconditioner.
            """
            c = np.asarray(self.compute_covariance_at(
                [float(v) for v in x], epsilon, floor))
            used = [int(k) for k in self.get_covariance_parameters()]
            k = len(used)
            return c.reshape(k, k), used

        def jacobian(self, x, epsilon=0.0, floor=0.0):
            """The forward-difference Jacobian at ``x``, shape (p, m).

            One row per parameter, one column per residual -- the shape
            chisurf's ``approx_grad`` returns, including the rows of zeros it
            keeps for parameters that do not move the objective.
            """
            j = np.asarray(self.compute_jacobian(
                [float(v) for v in x], epsilon, floor))
            n = len(x)
            if j.size == 0:
                return j.reshape(n, 0)
            return j.reshape(n, j.size // n)
    }
}

/*
 * Headers are %include'd directly here, which is IMP's own convention -- see
 * modules/core/pyext/swig.i-in, which does the same for 101 of them. A separate
 * .i file per header earns its place only when it carries real content:
 * IMP_SWIG_OBJECT/DECORATOR declarations for reference-counted types, %ignore
 * for a method that will not marshal, a %feature("shadow"). Two do
 * (IMP_bff.av.i, IMP_bff.pathmap.i, IMP_bff.observables.i); the fourteen
 * one-line wrappers that used
 * to sit beside them did not, and are gone.
 */

/* Reference-counted types: decorators, objects, and the restraint. */
/* The lattice PathMap is built on. It has to be wrapped *before* pathmap.i:
   SWIG only gives a class the methods of a base it has seen, and PathMap's
   whole numeric surface -- get_number_of_voxels, get_value, the voxel/location
   conversions -- is inherited. Before this was here PathMap still wrapped, and
   silently came out with none of it. */
/* GridHeader is wrapped as an ordinary class, deliberately NOT with
   IMP_SWIG_VALUE. That macro installs a typemap forbidding a value type from
   being returned by pointer, and `get_header()` returns one -- as
   `IMP::em::DensityMap::get_header` did before it, wrapped the same way. Using
   the macro and then %ignore-ing the pointer accessors removes them from
   Python altogether (an %extend of the same name is ignored too), and
   `test_av_lattice.py` reads the grid extent through `get_header()`. */
/* The raw voxel array stays in C++: `double*` trips the same value-return
   typemap, and Python already has bounds-keeping doors onto the same data
   (`get_xyz_density`, `get_tile_values`). */
%ignore IMP::bff::DensityGrid::get_data;

IMP_SWIG_VALUE(IMP::bff, GridSphere, GridSpheres);
IMP_SWIG_OBJECT(IMP::bff, DensityGrid, DensityGrids);
// The free write_mrc takes a float buffer; Python writes a grid with
// DensityGrid.write_mrc(path) and a path map's feature with write_map_feature.
%ignore IMP::bff::write_mrc;
%include "IMP/bff/DensityGrid.h"

%include "IMP_bff.pathmap.i"
/* Olga's name-keyed van der Waals radii and the AV's choice of radii set.
   Wrapped before av.i for the same reason: AV.h names AVRadiiSource. The
   `std::map` accessor and the per-particle vector stay C++-only (`#ifndef
   SWIG` in the header) -- `MapStringDouble` is not declared until
   IMP_bff.probe.i, below, and this module cannot name
   `std::vector<double>` at all (IMP_bff.types.i). Python reads the table
   through olga_vdw_radius(), get_olga_vdw_atom_names() and
   olga_vdw_radii_csv(). */
%include "IMP/bff/VdwRadii.h"

/* A label's states, whatever represents them: the cloud, the kernels it is
   measured with, the distances between two. Every representation and every
   consumer of a `States` comes after this -- including the decorator's
   Model doors (av.i), which return an AccessibleVolume. */
%include "IMP_bff.states.i"

/* The accessible volume: one representation of `States`, and the label
   distributions that produce one on demand. */
%include "IMP_bff.avmodel.i"

/* The inflated-sphere raster, core; the layer's particle view derives from it. */
%include "IMP_bff.occupancy.i"

/* fps.json: the definition, the reader, the writer -- and the distance
   measurement record, which the network restraint below takes by the map. */
%include "IMP_bff.fps.i"

/* Reading a frame out of a hierarchy without a SWIG call per atom, and the
   `ProteinFrame` value that reading is for. */
IMP_SWIG_VALUE(IMP::bff, ProteinFrame, ProteinFrames);
// Hidden from SWIG's member wrapping and republished as an `(n, 3)` table:
// a flat `coords` and a shaped `coords` cannot both be called `coords`.
%ignore IMP::bff::ProteinFrame::coords;
%include "IMP/bff/HierarchyFrame.h"
%attribute_np2(IMP::bff::ProteinFrame, std::vector<double>, coords,
               get_coords, 3, set_coords);
%attribute(IMP::bff::ProteinFrame, int, n_atoms, get_n_atoms);
%template(ProteinFrameList) std::vector<IMP::bff::ProteinFrame>;

/* Collisional quenching, short of the model: which moieties quench a dye and
   how hard, the fields a quenched dye lives in (stickiness, PET rate, FRET
   rate), the photon race along a walk, and the solvent-accessible surface that
   says how buried a quencher is. */
%include "IMP_bff.quenching.i"

/* A probe as a species -- dye, fluorescent protein or spin label: what it is,
   its spectra where it has them, the Förster radius they give, its .pto
   library, and the coarse-grained system it is simulated as. */
%include "IMP_bff.probe.i"

/* Distances and orientations over point clouds. */

/* Dynamics: the dye's position, then its excited state. */
%include "IMP/bff/BrownianWalk.h"
/* `diffusion_propagate` reports two arrays -- the population at each reported
   step and the final density -- and both are managed numpy views, so Python
   gets a 2-tuple. The population was a `std::vector<double>&` out-parameter,
   which made a caller construct a wrapped vector of a class this module does
   not own (see the note in `IMP_bff.types.i`). */
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** out_fluorescence, int* n_out_fluorescence)};
%include "IMP/bff/DiffusionSolver.h"
/* A dense network, evaluated in batches. The outputs are one managed view --
   `n_rows * n_outputs` of them -- for the same reason the diffusion solver's
   are: a walked SWIG proxy costs ~340 ns an element. */
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** out_view, int* n_out_view)};
%include "IMP/bff/NeuralNet.h"
/* The evaluation graph: labels naming ports, run on demand. Node is a
   director, so a graph of Python-subclassed nodes works here too. The
   provenance helpers that go with it are Python, so they live in their own
   interface file. */
%include "IMP_bff.evaluationgraph.i"
/* The photon trace returns delay and emitted-flag views, and the decay curve
   accumulates in place into the caller's histogram. The %apply has to sit
   here, before the header that declares them. */
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** out_delays, int* n_delays)};
%apply(unsigned char** ARGOUTVIEWM_ARRAY1, int* DIM1) {(unsigned char** out_emitted, int* n_emitted)};
%apply(double* INPLACE_ARRAY1, int DIM1) {(double* decay, int n_decay)};
%include "IMP/bff/PhotonSimulation.h"

/* Scoring orchestration: CHARMM36, LJ, Boltzmann, AABB, rotamer score. */
%include "IMP_bff.scoring.i"

/* Greedy Olga: which pair to measure next. */
%include "IMP_bff.greedyolga.i"

/* What comes out: (amplitude, rate) pairs. */
/* MOL2 in: the atoms and bonds a force-field component is built from. */
%include "IMP/bff/Mol2IO.h"
%template(Mol2Atoms) std::vector<IMP::bff::Mol2Atom>;
%template(MapIntString) std::map<int, std::string>;

/* Local sequence alignment: which template a fusion's domain matches. */
%include "IMP/bff/SequenceAlignment.h"
%template(AlignedBlocks) std::vector<IMP::bff::AlignedBlock>;

/* Weighted averaging over a rotamer ensemble. */

/* Protein side-chain packing: the 1:1 FASPR port (Dunbrack 2010).
   Attribution and the non-redistributed library note in the header. */

/* Leader clustering: how a conformer stack becomes a rotamer library. */
%include "IMP/bff/Clustering.h"

/* Bond connectivity, and the angles, torsions and rings that follow from it. */
%include "IMP/bff/MolecularGraph.h"
%template(VectorVectorInt) std::vector<std::vector<int> >;
// `LabelledGraph` answers in labels, so its returns need the string forms.
// `VectorVectorString` is instantiated in types.i for the multi-dye kernel.
%template(VectorPairStringString) std::vector<std::pair<std::string, std::string> >;
%template(VectorPairIntInt) std::vector<std::pair<int, int> >;
%template(PairIntInt) std::pair<int, int>;

/* Internal coordinates, general over dyes and side chains: measure a chi
   vector, rebuild coordinates from one (FASPR's construction, IUPAC
   convention). The .drot dye-rotamer store and the side-chain packer share
   this kernel. */
%include "IMP/bff/ZMatrix.h"
%template(VectorPairIntInt) std::vector<std::pair<int, int> >;
%template(VectorVector3D) std::vector<IMP::algebra::Vector3D>;

/* The force-field system in, through the ihm C reader IMP vendors. */
%include "IMP/bff/ForceFieldCIF.h"

/* The force-field system out, template CIF, and rotamer library IO. */
// The rotamer library value -- an ensemble of conformers with a weight
// each -- before the three readers that return one: the numpy/text pair
// (componenttemplate.i), the PDB+trajectory loader (probesampling.i)
// and `read_drot`.
IMP_SWIG_VALUE(IMP::bff, RotamerLibrary, RotamerLibraries);
// Same as `ProteinFrame::coords`: the shaped attributes below replace the raw
// members, so a library answers `(n_rotamers, n_atoms, 3)` and not a flat
// vector a caller has to fold by hand.
%ignore IMP::bff::RotamerLibrary::coords;
%ignore IMP::bff::RotamerLibrary::weights;
/* RotamerLibrary.h now carries the `.drot` reader and writer, and
   `write_drot_with_provenance` takes #MfdbTags, so the profile header has to be
   wrapped before this one rather than in labelizer.i, where it used to sit.
   The value semantics must come with it. Without these, MfdbTags arrives as an opaque
   SwigPyObject with no destructor and the wrapper rejects it by type. */
IMP_SWIG_VALUE(IMP::bff, MfdbTag, MfdbTags);
IMP_SWIG_VALUE(IMP::bff, MfdbColumn, MfdbColumns);
IMP_SWIG_VALUE(IMP::bff, MfdbAttribution, MfdbAttributions);
/* Pto.h is the container too, since PRD-138 folded PtoProfile.h into it. The
   container classes are kept out of the binding for now: nothing in Python
   reads or writes a .pto except through the typed doors (the rotamer library,
   the probe container, the Labelizer store), and whether to publish a raw
   reader/writer is the interface pass's decision, not a consolidation's. */
%ignore IMP::bff::PtoObject;
%ignore IMP::bff::PtoWriter;
%ignore IMP::bff::PtoReader;
%include "IMP/bff/Pto.h"
%include "IMP/bff/RotamerLibrary.h"
%attribute_np3(IMP::bff::RotamerLibrary, std::vector<double>, coords,
               get_coords, n_rotamers, 3, set_coords);
%attribute_np(IMP::bff::RotamerLibrary, std::vector<double>, weights,
              get_weights, set_weights);

%include "IMP_bff.componenttemplate.i"

/* The system layer: Label, Quencher, attachment, backbone frame, strip. */
/* The selection language, before anything that takes a mask. */
%include "IMP_bff.selection.i"

%include "IMP_bff.label.i"

/* The cgprobe topology builder: force-field system assembly, graph helpers. */
%include "IMP_bff.topology.i"

/* The cgprobe samplers: linker Metropolis, rotamer libraries, RRT, kinetics. */
/* where a kernel runs, and how an accelerator is found */
%include "IMP_bff.compute.i"

/* the shape every simulation shares (PRD-139) */
%include "IMP_bff.simulation.i"
%include "IMP_bff.sampling.i"

/* The cgprobe MD runner: IMP MD/MC driver. */
%include "IMP_bff.sim.i"

/* Building an accessible volume: the two doors, the PDB read and the strip. */
%include "IMP_bff.avbuilder.i"

/* The native Labelizer: per-residue label-site scores and the FRET pair score.
   After `avbuilder.i`, whose `PDBAtomRecord` reader it groups into residues,
   and whose accessible volumes the pair layer measures distances over. */
%include "IMP_bff.labelizer.i"

/* Structures and trajectories in and out: PDB, MOL2, mmCIF, BinaryCIF, DCD.
   After `avbuilder.i`, which declares the `PDBAtomRecord` the MOL2 writer takes
   -- SWIG resolves a type at the point of use. */
%include "IMP_bff.structureio.i"

/* The processes that deactivate or depolarise a dye. */
%include "IMP_bff.photophysics.i"

/* The channels that deactivate an excited dye. They consume `States`, so they
   come after the object that is one. */
%include "IMP_bff.interactionterms.i"

/* The output contract, a C++ value with a Python surface. Its `from_states`
   reduction sums the terms, so it comes after them. */
%include "IMP_bff.observables.i"

/* FRET over every pair of two labelled ensembles, and the rate along one
   dye's trajectory. The pair values consume `States`. */
%include "IMP_bff.fret.i"

/* The closed-form distance distributions and the FCS shapes. The distances
   between two labels that used to sit here are in states.i, with the
   States they consume (PRD-138); this header is wrapped where it always was. */
%include "IMP_bff.distributions.i"

/* Chi-squared scoring of labelling data, with and without volumes. */
%include "IMP_bff.proberestraints.i"

/* What repeated docking says about a model's precision (the restraint that
   does the docking is the connection layer's, in layer.i). */
%include "IMP_bff.modelprecision.i"

/* The particle picture of a tethered dye: the walk, the photons, the
   equilibrium, as kernels and as an object. */
%include "IMP_bff.probesampling.i"

/* ...and the field picture. */
%include "IMP_bff.griddiffusion.i"

/* One labelled site's donor decay, both ways. Last of the physics: it consumes
   an accessible volume, the PET tables, the walk and the grid solver. */
%include "IMP_bff.quenchingmodel.i"

/* The third rotamer family in the same container: the Dunbrack-2010
   backbone-dependent side-chain table, which is a distribution over (phi,
   psi) rather than an ensemble and so has kinds of its own. */
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** out_view, int* n_out_view)};

/* The rotamer layer, all of it: RotamerEnsemble subclasses States (a SWIG
   value type) and returns the FRET pair values, so the whole header waits
   until after both. */
%include "IMP_bff.rotamer.i"

/* The parameter tables the coarse-grained potentials read; the potentials
   themselves are IMP restraints, in layer.i. */
%include "IMP_bff.potentialtables.i"

/*
 * The flat user-facing surface: every public name is an attribute of `IMP.bff`
 * itself, a SWIG wrapper over a C++ declaration.
 *
 * Nothing is lazy and there is no `__getattr__` hook: a missing attribute is
 * a missing attribute, and `dir()` is the module's own. `rmf` and `rotamer`
 * are modules of this one (`dependencies.py`), so the readers and writers
 * over them are ordinary C++. A PMI restraint has no C++ spelling, and what
 * one would do apart from PMI's bookkeeping is
 * `probe_network_restraint_set`.
 */

/* where the large data lives, for both builds (PRD-137 6d) */
%include "IMP_bff.data.i"

/* LabelLib's API, name for name (PRD-137 6d) */
%include "IMP_bff.labellib.i"
