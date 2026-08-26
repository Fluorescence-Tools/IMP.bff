/*
 * The channels that deactivate an excited dye.
 *
 * These were an `abc.ABC` and three `@dataclass`es in `photophysics.py`, and
 * they are `IMP::Object`s here rather than values -- `total_rate` sums a list
 * of *different* terms, which needs one polymorphic handle rather than a
 * value per type.
 *
 * Python's duck typing let each term take the participants it needed
 * (`rate_constants(states)`, `(states, atoms)`, `(donor, acceptor)`); C++ takes
 * the same two and reports through `arity` which it reads. `PETTerm`'s
 * quenching atoms move into the term, because resolving which of a structure's
 * atoms quench and how hard is a function of names and parameters and does not
 * change from one set of dye states to the next.
 *
 * The two constructors had a `%feature("shadow")` each, for three reasons that
 * are all typemaps or features now: keyword arguments with C++ defaults
 * (`%feature("compactdefaultargs")` -- SWIG's default-argument *overloads* are
 * what disables `kwargs`, and compact ones do not), the conversions of a dict,
 * two name lists and an `(N, 3)` array (std_map, std_vector and numpy.i do
 * those), and `kappa2=None` for "no orientation factor" (a typemap, below,
 * since NaN is how C++ spells it).
 */

IMP_SWIG_OBJECT(IMP::bff, InteractionTerm, InteractionTerms);
IMP_SWIG_OBJECT(IMP::bff, RadiativeTerm, RadiativeTerms);
IMP_SWIG_OBJECT(IMP::bff, PETTerm, PETTerms);
IMP_SWIG_OBJECT(IMP::bff, FRETTerm, FRETTerms);

%attribute_py(IMP::bff::InteractionTerm, int, arity, get_arity);
%attribute_py(IMP::bff::InteractionTerm, bool, needs_orientations,
              get_needs_orientations);
%attribute_py(IMP::bff::RadiativeTerm, double, lifetime, get_lifetime);
%attribute_py(IMP::bff::PETTerm, double, dye_radius, get_dye_radius);
%attribute_py(IMP::bff::FRETTerm, double, forster_radius, get_forster_radius);
%attribute_py(IMP::bff::FRETTerm, bool, used_isotropic_kappa2,
              get_used_isotropic_kappa2);
%attribute_py(IMP::bff::FRETTerm, double, kappa2, get_kappa2);

// The quencher atoms, straight from numpy.
%apply(double* IN_ARRAY2, int DIM1, int DIM2) {
    (double* coords, int n_atoms, int n_dim)
};

// `kappa2=None` is "no orientation factor given", and C++ spells that NaN.
// The typemap is on the argument's name, so it applies to the one parameter
// that means it and to no other double.
%typemap(in) double kappa2 {
    if ($input == Py_None) {
        $1 = std::numeric_limits<double>::quiet_NaN();
    } else {
        double v = PyFloat_AsDouble($input);
        if (v == -1.0 && PyErr_Occurred()) SWIG_fail;
        $1 = v;
    }
}
%typemap(typecheck, precedence=SWIG_TYPECHECK_DOUBLE) double kappa2 {
    $1 = ($input == Py_None || PyNumber_Check($input)) ? 1 : 0;
}

// Keyword arguments *and* C++ default arguments: SWIG expands defaults into
// overloads by default, and it refuses `kwargs` for an overloaded function.
// A compact default keeps one wrapper, so both work.
%feature("compactdefaultargs") IMP::bff::PETTerm::PETTerm;
%feature("compactdefaultargs") IMP::bff::FRETTerm::FRETTerm;
%feature("kwargs") IMP::bff::PETTerm::PETTerm;
%feature("kwargs") IMP::bff::FRETTerm::FRETTerm;

%include "IMP/bff/InteractionTerms.h"

// `get_kQ`/`get_rC` publish managed numpy views (ARGOUTVIEWM_ARRAY1); a
// property over them hands the array straight back, no Python in between.
%attribute_py(IMP::bff::PETTerm, PyObject*, kQ, get_kQ);
%attribute_py(IMP::bff::PETTerm, PyObject*, rC, get_rC);
