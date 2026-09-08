/*
 * What IMP's kernel interface provides to a module's SWIG and the standalone
 * build must provide itself (PRD-137 step 6c): the exception translation,
 * the IMP_SWIG_* declarations core.i is written with (values as Python lists,
 * reference-counted objects, directors, cereal pickling, show() as str/repr),
 * and the algebra vectors as plain sequences. Reduced to what the core uses.
 */

%{
#include <IMP/bff/Base.h>
#include <IMP/algebra/VectorD.h>
#include <cereal/archives/binary.hpp>
#include <sstream>
#include <stdexcept>
%}

/* ---- exceptions: IMP's hierarchy onto Python's ---- */
/* ---- exceptions, as IMP's kernel does them ----
   IMP creates Python exception classes at module init that mirror the C++
   ones (IMP.Exception and its family; ValueException also derives from
   ValueError, IOException from IOError, IndexException from IndexError,
   TypeException from TypeError) and translates every C++ exception in one
   handler, including the std:: ones the core throws directly
   (std::domain_error is a ValueException, std::out_of_range an
   IndexException). Same classes, same mapping, under IMP.bff. */
%{
static PyObject *imp_bff_exception = 0, *imp_bff_internal_exception = 0,
    *imp_bff_model_exception = 0, *imp_bff_usage_exception = 0,
    *imp_bff_index_exception = 0, *imp_bff_io_exception = 0,
    *imp_bff_value_exception = 0, *imp_bff_event_exception = 0,
    *imp_bff_type_exception = 0;

static PyObject *imp_bff_new_exception(PyObject *m, const char *qualname,
                                       const char *name, PyObject *base,
                                       PyObject *pybase) {
    PyObject *cls;
    if (pybase) {
        PyObject *bases = PyTuple_Pack(2, base, pybase);
        cls = PyErr_NewException(qualname, bases, NULL);
        Py_DECREF(bases);
    } else {
        cls = PyErr_NewException(qualname, base, NULL);
    }
    Py_INCREF(cls);
    PyModule_AddObject(m, name, cls);
    return cls;
}

static void imp_bff_handle_exception(void) {
    try {
        throw;
    /* std:: exceptions, mapped as IMP maps them */
    } catch (const std::out_of_range &e) {
        PyErr_SetString(imp_bff_index_exception, e.what());
    } catch (const std::domain_error &e) {
        PyErr_SetString(imp_bff_value_exception, e.what());
    } catch (const std::ios::failure &e) {
        PyErr_SetString(imp_bff_io_exception, e.what());
    } catch (const std::length_error &e) {
        PyErr_SetString(imp_bff_internal_exception, e.what());
    /* the module's own (Base.h's standalone branch) */
    } catch (const IMP::ValueException &e) {
        PyErr_SetString(imp_bff_value_exception, e.what());
    } catch (const IMP::ModelException &e) {
        PyErr_SetString(imp_bff_model_exception, e.what());
    } catch (const IMP::UsageException &e) {
        PyErr_SetString(imp_bff_usage_exception, e.what());
    } catch (const IMP::IOException &e) {
        PyErr_SetString(imp_bff_io_exception, e.what());
    } catch (const IMP::Exception &e) {
        PyErr_SetString(imp_bff_exception, e.what());
    } catch (const std::bad_alloc &e) {
        PyErr_SetString(PyExc_MemoryError, e.what());
    } catch (const std::exception &e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
    } catch (...) {
        PyErr_SetString(PyExc_RuntimeError, "Unknown error caught by Python wrapper");
    }
}
%}

%init %{
    imp_bff_exception = imp_bff_new_exception(m, "_IMP_bff.Exception", "Exception", NULL, NULL);
    imp_bff_internal_exception = imp_bff_new_exception(m, "_IMP_bff.InternalException", "InternalException", imp_bff_exception, NULL);
    imp_bff_model_exception = imp_bff_new_exception(m, "_IMP_bff.ModelException", "ModelException", imp_bff_exception, NULL);
    imp_bff_usage_exception = imp_bff_new_exception(m, "_IMP_bff.UsageException", "UsageException", imp_bff_exception, NULL);
    imp_bff_event_exception = imp_bff_new_exception(m, "_IMP_bff.EventException", "EventException", imp_bff_exception, NULL);
    imp_bff_index_exception = imp_bff_new_exception(m, "_IMP_bff.IndexException", "IndexException", imp_bff_exception, PyExc_IndexError);
    imp_bff_io_exception = imp_bff_new_exception(m, "_IMP_bff.IOException", "IOException", imp_bff_exception, PyExc_IOError);
    imp_bff_value_exception = imp_bff_new_exception(m, "_IMP_bff.ValueException", "ValueException", imp_bff_exception, PyExc_ValueError);
    imp_bff_type_exception = imp_bff_new_exception(m, "_IMP_bff.TypeException", "TypeException", imp_bff_exception, PyExc_TypeError);
%}

%pythoncode %{
Exception = _IMP_bff.Exception
InternalException = _IMP_bff.InternalException
ModelException = _IMP_bff.ModelException
UsageException = _IMP_bff.UsageException
EventException = _IMP_bff.EventException
IndexException = _IMP_bff.IndexException
IOException = _IMP_bff.IOException
ValueException = _IMP_bff.ValueException
TypeException = _IMP_bff.TypeException
%}

%exception {
    try {
        $action
    } catch (...) {
        /* a director method that raised has already set the Python error */
        if (!PyErr_Occurred()) imp_bff_handle_exception();
        SWIG_fail;
    }
}
%feature("director:except") {
    if ($error != NULL) {
        throw Swig::DirectorMethodException();
    }
}

/* ---- show() as __str__ / __repr__ ---- */
%define IMPBFF_SHOWABLE(Namespace, Name)
%extend Namespace::Name {
    std::string __str__() const { std::ostringstream o; self->show(o); return o.str(); }
    std::string __repr__() const { std::ostringstream o; self->show(o); return o.str(); }
}
%enddef

/* ---- values: a std::vector of them crosses as a Python list, as in IMP ---- */
/* IMP spells a value's plural `IMP::Vector<T>`, which derives from
   std::vector<T>; the shim build spells it std::vector<T> outright. SWIG has
   to know the type it will actually write into the wrapper, so both
   spellings get the same typemaps. Declared here for both builds -- where
   nothing uses it, nothing matches it. */
namespace IMP {
template <class T>
class Vector : public std::vector<T> {
 public:
    Vector();
    Vector(const std::vector<T>& o);
};
}

%define IMP_SWIG_VALUE_VECTOR_TYPEMAPS(Namespace, Name, VectorType)
%typemap(out) VectorType {
    // bound to a reference first: SWIG hands a value return through
    // SwigValueWrapper, which is not the vector but converts to it
    const VectorType& imp_bff_out = $1;
    $result = PyList_New(imp_bff_out.size());
    for (size_t i = 0; i < imp_bff_out.size(); ++i) {
        PyList_SET_ITEM($result, i, SWIG_NewPointerObj(new Namespace::Name(imp_bff_out[i]), $descriptor(Namespace::Name*), SWIG_POINTER_OWN));
    }
}
%typemap(out) const VectorType& {
    $result = PyList_New($1->size());
    for (size_t i = 0; i < $1->size(); ++i) {
        PyList_SET_ITEM($result, i, SWIG_NewPointerObj(new Namespace::Name((*$1)[i]), $descriptor(Namespace::Name*), SWIG_POINTER_OWN));
    }
}
%typemap(in) const VectorType& (VectorType imp_bff_tmp) {
    VectorType* imp_bff_direct = 0;
    if (SWIG_IsOK(SWIG_ConvertPtr($input, (void**)&imp_bff_direct, $descriptor(VectorType*), 0)) && imp_bff_direct) {
        $1 = imp_bff_direct;
    } else {
        if (!PySequence_Check($input) || PyUnicode_Check($input)) { SWIG_exception_fail(SWIG_TypeError, "a sequence of " #Name " is needed"); }
        Py_ssize_t n = PySequence_Size($input);
        for (Py_ssize_t i = 0; i < n; ++i) {
            PyObject* item = PySequence_GetItem($input, i);
            Namespace::Name* p = 0;
            int res = SWIG_ConvertPtr(item, (void**)&p, $descriptor(Namespace::Name*), 0);
            if (!SWIG_IsOK(res) || !p) { Py_DECREF(item); SWIG_exception_fail(SWIG_TypeError, "item " #Name " expected"); }
            imp_bff_tmp.push_back(*p);
            Py_DECREF(item);
        }
        $1 = &imp_bff_tmp;
    }
}
%typemap(in) VectorType (VectorType imp_bff_tmp) {
    VectorType* imp_bff_direct = 0;  // NOLINT
    if (SWIG_IsOK(SWIG_ConvertPtr($input, (void**)&imp_bff_direct, $descriptor(VectorType*), 0)) && imp_bff_direct) {
        imp_bff_tmp = *imp_bff_direct;
    } else {
        if (!PySequence_Check($input) || PyUnicode_Check($input)) { SWIG_exception_fail(SWIG_TypeError, "a sequence of " #Name " is needed"); }
        Py_ssize_t n = PySequence_Size($input);
        for (Py_ssize_t i = 0; i < n; ++i) {
            PyObject* item = PySequence_GetItem($input, i);
            Namespace::Name* p = 0;
            int res = SWIG_ConvertPtr(item, (void**)&p, $descriptor(Namespace::Name*), 0);
            if (!SWIG_IsOK(res) || !p) { Py_DECREF(item); SWIG_exception_fail(SWIG_TypeError, "item " #Name " expected"); }
            imp_bff_tmp.push_back(*p);
            Py_DECREF(item);
        }
    }
    $1 = imp_bff_tmp;
}
%typemap(typecheck, precedence=SWIG_TYPECHECK_VECTOR) const VectorType&, VectorType {
    $1 = 0;
    void* imp_bff_direct = 0;
    if (SWIG_IsOK(SWIG_ConvertPtr($input, &imp_bff_direct, $descriptor(VectorType*), 0))) {
        $1 = 1;
    } else if (PySequence_Check($input) && !PyUnicode_Check($input)) {
        if (PySequence_Size($input) == 0) { $1 = 1; }
        else {
            PyObject* item = PySequence_GetItem($input, 0);
            void* p = 0;
            $1 = SWIG_IsOK(SWIG_ConvertPtr(item, &p, $descriptor(Namespace::Name*), 0)) ? 1 : 0;
            Py_DECREF(item);
        }
    }
}
%enddef

%define IMP_SWIG_VALUE(Namespace, Name, PluralName)
IMP_SWIG_VALUE_VECTOR_TYPEMAPS(Namespace, Name, IMP::Vector<Namespace::Name >)
%typemap(out) Namespace::Name const& {
    $result = SWIG_NewPointerObj(new Namespace::Name(*$1), $descriptor(Namespace::Name*), SWIG_POINTER_OWN);
}
%typemap(out) std::vector<Namespace::Name> {
    $result = PyList_New($1.size());
    for (size_t i = 0; i < $1.size(); ++i) {
        PyList_SET_ITEM($result, i, SWIG_NewPointerObj(new Namespace::Name($1[i]), $descriptor(Namespace::Name*), SWIG_POINTER_OWN));
    }
}
%typemap(out) const std::vector<Namespace::Name>& {
    $result = PyList_New($1->size());
    for (size_t i = 0; i < $1->size(); ++i) {
        PyList_SET_ITEM($result, i, SWIG_NewPointerObj(new Namespace::Name((*$1)[i]), $descriptor(Namespace::Name*), SWIG_POINTER_OWN));
    }
}
/* A wrapped std::vector<Name> (the PluralName##List template) is taken as it
   is; any other sequence is copied item by item. The copy happens before the
   item's reference is dropped: a wrapped vector's __getitem__ hands out a
   fresh owning proxy, and dropping it first would leave `p` dangling. */
%typemap(in) const std::vector<Namespace::Name>& (std::vector<Namespace::Name> imp_bff_tmp) {
    std::vector<Namespace::Name>* imp_bff_direct = 0;
    if (SWIG_IsOK(SWIG_ConvertPtr($input, (void**)&imp_bff_direct, $descriptor(std::vector<Namespace::Name>*), 0)) && imp_bff_direct) {
        $1 = imp_bff_direct;
    } else {
        if (!PySequence_Check($input) || PyUnicode_Check($input)) { SWIG_exception_fail(SWIG_TypeError, "a sequence of " #Name " is needed"); }
        Py_ssize_t n = PySequence_Size($input);
        imp_bff_tmp.reserve(n);
        for (Py_ssize_t i = 0; i < n; ++i) {
            PyObject* item = PySequence_GetItem($input, i);
            Namespace::Name* p = 0;
            int res = SWIG_ConvertPtr(item, (void**)&p, $descriptor(Namespace::Name*), 0);
            if (!SWIG_IsOK(res) || !p) { Py_DECREF(item); SWIG_exception_fail(SWIG_TypeError, "item " #Name " expected"); }
            imp_bff_tmp.push_back(*p);
            Py_DECREF(item);
        }
        $1 = &imp_bff_tmp;
    }
}
%typemap(in) std::vector<Namespace::Name> (std::vector<Namespace::Name> imp_bff_tmp) {
    std::vector<Namespace::Name>* imp_bff_direct = 0;
    if (SWIG_IsOK(SWIG_ConvertPtr($input, (void**)&imp_bff_direct, $descriptor(std::vector<Namespace::Name>*), 0)) && imp_bff_direct) {
        imp_bff_tmp = *imp_bff_direct;
    } else {
        if (!PySequence_Check($input) || PyUnicode_Check($input)) { SWIG_exception_fail(SWIG_TypeError, "a sequence of " #Name " is needed"); }
        Py_ssize_t n = PySequence_Size($input);
        for (Py_ssize_t i = 0; i < n; ++i) {
            PyObject* item = PySequence_GetItem($input, i);
            Namespace::Name* p = 0;
            int res = SWIG_ConvertPtr(item, (void**)&p, $descriptor(Namespace::Name*), 0);
            if (!SWIG_IsOK(res) || !p) { Py_DECREF(item); SWIG_exception_fail(SWIG_TypeError, "item " #Name " expected"); }
            imp_bff_tmp.push_back(*p);
            Py_DECREF(item);
        }
    }
    $1 = imp_bff_tmp;
}
%typemap(typecheck, precedence=SWIG_TYPECHECK_VECTOR) const std::vector<Namespace::Name>&, std::vector<Namespace::Name> {
    $1 = 0;
    void* imp_bff_direct = 0;
    if (SWIG_IsOK(SWIG_ConvertPtr($input, &imp_bff_direct, $descriptor(std::vector<Namespace::Name>*), 0))) {
        $1 = 1;
    } else if (PySequence_Check($input) && !PyUnicode_Check($input)) {
        if (PySequence_Size($input) == 0) { $1 = 1; }
        else {
            PyObject* item = PySequence_GetItem($input, 0);
            void* p = 0;
            $1 = SWIG_IsOK(SWIG_ConvertPtr(item, &p, $descriptor(Namespace::Name*), 0)) ? 1 : 0;
            Py_DECREF(item);
        }
    }
}
%pythoncode %{
PluralName = list
%}
IMPBFF_SHOWABLE(Namespace, Name)
%enddef

/* IMP_SWIG_VALUE_INSTANCE and friends are the layer's; the core does not use them. */

/* ---- objects: intrusive reference counting, as IMP's Object ---- */
%define IMP_SWIG_OBJECT(Namespace, Name, PluralName)
%feature("ref") Namespace::Name "$this->ref();"
%feature("unref") Namespace::Name "$this->unref();"
%typemap(out) std::vector<IMP::Pointer<Namespace::Name> >, const std::vector<IMP::Pointer<Namespace::Name> >& {
    const std::vector<IMP::Pointer<Namespace::Name> >& v = *&$1;
    $result = PyList_New(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        Namespace::Name* p = v[i].get();
        PyList_SET_ITEM($result, i, SWIG_NewPointerObj(p, $descriptor(Namespace::Name*), SWIG_POINTER_OWN));
    }
}
%typemap(out) std::vector<Namespace::Name*> {
    $result = PyList_New($1.size());
    for (size_t i = 0; i < $1.size(); ++i) {
        PyList_SET_ITEM($result, i, SWIG_NewPointerObj($1[i], $descriptor(Namespace::Name*), SWIG_POINTER_OWN));
    }
}
%typemap(in) const std::vector<Namespace::Name*>& (std::vector<Namespace::Name*> imp_bff_tmp) {
    if (!PySequence_Check($input)) { SWIG_exception_fail(SWIG_TypeError, "a sequence of " #Name " is needed"); }
    Py_ssize_t n = PySequence_Size($input);
    for (Py_ssize_t i = 0; i < n; ++i) {
        PyObject* item = PySequence_GetItem($input, i);
        Namespace::Name* p = 0;
        int res = SWIG_ConvertPtr(item, (void**)&p, $descriptor(Namespace::Name*), 0);
        if (!SWIG_IsOK(res)) { Py_DECREF(item); SWIG_exception_fail(SWIG_TypeError, "item " #Name " expected"); }
        imp_bff_tmp.push_back(p);
        Py_DECREF(item);
    }
    $1 = &imp_bff_tmp;
}
%pythoncode %{
PluralName = list
%}
IMPBFF_SHOWABLE(Namespace, Name)
%enddef

/* IMP's answer to the director lifetime problem (IMP_kernel.directors.i):
   the C++ side keeps only a weak pointer to the Python instance, so once
   the caller drops its proxy a callback would land on a dead object. Every
   Python subclass instance is therefore held in a registry, released again
   when nothing else -- Python or C++ (the ref count) -- refers to it. */
%pythoncode %{
import sys as _imp_bff_sys

class _DirectorObjects:
    """@internal Keeps Python subclasses of director classes alive while
       C++ still refers to them."""
    def __init__(self):
        self._objects = []

    def register(self, obj):
        # Unconditionally: where IMP is linked, IMP::Object is not wrapped --
        # SWIG is told the macros, not the class -- so `get_ref_count` is not
        # on the proxy, and a registry that asked for it first held nothing at
        # all. The callback then landed on a collected object and the process
        # died. cleanup() below copes with either shape.
        self.cleanup()
        self._objects.append(obj)

    def cleanup(self):
        # three references are ours: the list, `x`, and getrefcount's argument
        def still_needed(x):
            if _imp_bff_sys.getrefcount(x) > 3:
                return True
            counted = getattr(x, "get_ref_count", None)
            # without the C++ count, Python's is all there is to go on
            return counted() > 1 if counted is not None else False

        self._objects = [x for x in self._objects if still_needed(x)]

    def get_object_count(self):
        return len(self._objects)

_director_objects = _DirectorObjects()
%}

%define IMP_SWIG_DIRECTOR(Namespace, Name)
%feature("director") Namespace::Name;
%pythonappend Namespace::Name::Name %{
        if self.__class__ != Name:
            _director_objects.register(self)
%}
%enddef


/* ---- pickling through cereal, as IMP's values do ---- */
%define IMP_SWIG_VALUE_SERIALIZE_IMPL(Namespace, Name)
%extend Namespace::Name {
    PyObject* _get_as_binary() const {
        std::ostringstream oss;
        { cereal::BinaryOutputArchive ba(oss); ba(*self); }
        const std::string s = oss.str();
        return PyBytes_FromStringAndSize(s.data(), s.size());
    }
    void _set_from_binary(PyObject* p) {
        char* buf; Py_ssize_t len;
        if (PyBytes_AsStringAndSize(p, &buf, &len) < 0) { throw IMP::ValueException("bytes expected"); }
        std::istringstream iss(std::string(buf, len));
        cereal::BinaryInputArchive ba(iss);
        ba(*self);
    }
    %pythoncode %{
    def __getstate__(self):
        p = self._get_as_binary()
        if len(self.__dict__) > 1:
            d = self.__dict__.copy()
            del d['this']
            p = (d, p)
        return p
    def __setstate__(self, p):
        if not hasattr(self, 'this'):
            self.__init__()
        if isinstance(p, tuple):
            d, p = p
            self.__dict__.update(d)
        return self._set_from_binary(p)
    %}
}
%enddef

/* ---- IMP::algebra vectors as plain sequences ---- */
%define IMPBFF_VECTOR_TYPEMAPS(N, Type)
%typemap(in) Type (Type imp_bff_v), const Type& (Type imp_bff_v) {
    if (!PySequence_Check($input) || PySequence_Size($input) != N) { SWIG_exception_fail(SWIG_TypeError, "a sequence of " #N " numbers is needed"); }
    for (int i = 0; i < N; ++i) {
        PyObject* item = PySequence_GetItem($input, i);
        double d = PyFloat_AsDouble(item);
        Py_DECREF(item);
        if (d == -1.0 && PyErr_Occurred()) { SWIG_exception_fail(SWIG_TypeError, "a number is needed"); }
        imp_bff_v[i] = d;
    }
    $1 = &imp_bff_v;
}
%typemap(in) Type (Type imp_bff_v) {
    if (!PySequence_Check($input) || PySequence_Size($input) != N) { SWIG_exception_fail(SWIG_TypeError, "a sequence of " #N " numbers is needed"); }
    for (int i = 0; i < N; ++i) {
        PyObject* item = PySequence_GetItem($input, i);
        double d = PyFloat_AsDouble(item);
        Py_DECREF(item);
        if (d == -1.0 && PyErr_Occurred()) { SWIG_exception_fail(SWIG_TypeError, "a number is needed"); }
        imp_bff_v[i] = d;
    }
    $1 = imp_bff_v;
}
%typemap(typecheck, precedence=SWIG_TYPECHECK_DOUBLE_ARRAY) Type, const Type& {
    $1 = (PySequence_Check($input) && !PyUnicode_Check($input) && PySequence_Size($input) == N) ? 1 : 0;
}
%typemap(out) Type {
    const Type& v = $1;
    $result = PyTuple_New(N);
    for (int i = 0; i < N; ++i) PyTuple_SET_ITEM($result, i, PyFloat_FromDouble(v[i]));
}
%typemap(out) const Type&, Type& {
    const Type& v = *$1;
    $result = PyTuple_New(N);
    for (int i = 0; i < N; ++i) PyTuple_SET_ITEM($result, i, PyFloat_FromDouble(v[i]));
}
%typemap(out) std::vector<Type> {
    const std::vector<Type>& v = $1;
    $result = PyList_New(v.size());
    for (size_t k = 0; k < v.size(); ++k) {
        PyObject* t = PyTuple_New(N);
        for (int i = 0; i < N; ++i) PyTuple_SET_ITEM(t, i, PyFloat_FromDouble(v[k][i]));
        PyList_SET_ITEM($result, k, t);
    }
}
%typemap(out) const std::vector<Type>&, std::vector<Type>& {
    const std::vector<Type>& v = *$1;
    $result = PyList_New(v.size());
    for (size_t k = 0; k < v.size(); ++k) {
        PyObject* t = PyTuple_New(N);
        for (int i = 0; i < N; ++i) PyTuple_SET_ITEM(t, i, PyFloat_FromDouble(v[k][i]));
        PyList_SET_ITEM($result, k, t);
    }
}
%typemap(in) const std::vector<Type>& (std::vector<Type> imp_bff_vs) {
    if (!PySequence_Check($input)) { SWIG_exception_fail(SWIG_TypeError, "a sequence of " #N "-vectors is needed"); }
    Py_ssize_t n = PySequence_Size($input);
    for (Py_ssize_t k = 0; k < n; ++k) {
        PyObject* row = PySequence_GetItem($input, k);
        Type v;
        bool ok = PySequence_Check(row) && PySequence_Size(row) == N;
        for (int i = 0; ok && i < N; ++i) {
            PyObject* item = PySequence_GetItem(row, i);
            v[i] = PyFloat_AsDouble(item);
            Py_DECREF(item);
            if (v[i] == -1.0 && PyErr_Occurred()) ok = false;
        }
        Py_DECREF(row);
        if (!ok) { SWIG_exception_fail(SWIG_TypeError, "row " #N " numbers expected"); }
        imp_bff_vs.push_back(v);
    }
    $1 = &imp_bff_vs;
}
%typemap(typecheck, precedence=SWIG_TYPECHECK_VECTOR) const std::vector<Type>& {
    $1 = (PySequence_Check($input) && !PyUnicode_Check($input)) ? 1 : 0;
}
%enddef
// SWIG must know these names exist so that the headers' unqualified
// `algebra::Vector3D` (written inside IMP::bff) resolves to the typemapped
// type instead of being emitted verbatim into the wrapper.
/* IMP declares VectorD in namespace IMP and pulls it into IMP::algebra, and
   headers use both spellings; SWIG matches typemaps on the spelling it sees,
   so both are declared and both get the typemaps below. */
namespace IMP {
template <int D> class VectorD {
 public:
    VectorD();
    VectorD(const VectorD<D>& o);
    double operator[](unsigned int i) const;
};
}

namespace IMP { namespace algebra {
// Enough of the shape for SWIG to treat it as an ordinary value: without a
// default constructor in view it wraps every return in SwigValueWrapper and
// writes code the real class does not support.
using IMP::VectorD;
typedef VectorD<3> Vector3D;
typedef VectorD<4> Vector4D;
typedef IMP::Vector<VectorD<3> > Vector3Ds;
typedef IMP::Vector<VectorD<4> > Vector4Ds;
} }
IMPBFF_VECTOR_TYPEMAPS(3, IMP::algebra::Vector3D)
IMPBFF_VECTOR_TYPEMAPS(4, IMP::algebra::Vector4D)
IMPBFF_VECTOR_TYPEMAPS(3, IMP::VectorD<3>)
IMPBFF_VECTOR_TYPEMAPS(4, IMP::VectorD<4>)
