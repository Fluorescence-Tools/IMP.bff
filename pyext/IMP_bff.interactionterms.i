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

// Both take every argument by keyword with defaults, which SWIG cannot express
// for a function that has them; the shadows also spell `None` for "no kappa2"
// and "no attenuation", which is NaN on the C++ side.
%feature("shadow") IMP::bff::PETTerm::PETTerm %{
def __init__(self, parameters, res_names, atom_names, coords, dye_radius=3.5):
    _IMP_bff.PETTerm_swiginit(self, _IMP_bff.new_PETTerm(
        dict(parameters), [str(r) for r in res_names],
        [str(a) for a in atom_names],
        np.ascontiguousarray(np.asarray(coords, dtype=np.float64)),
        float(dye_radius)))
%}

%feature("shadow") IMP::bff::FRETTerm::FRETTerm %{
def __init__(self, donor, acceptor, refractive_index=1.4, kappa2=None,
             r_min=7.0):
    _IMP_bff.FRETTerm_swiginit(self, _IMP_bff.new_FRETTerm(
        donor, acceptor, float(refractive_index),
        float("nan") if kappa2 is None else float(kappa2), float(r_min)))
%}

%include "IMP/bff/InteractionTerms.h"

%extend IMP::bff::PETTerm {
    %pythoncode %{
        @property
        def kQ(self):
            """Per-atom quenching rate, 1/ns. Zero where the atom is inert."""
            return _IMP_bff.PETTerm_get_kQ(self)

        @property
        def rC(self):
            """Per-atom attenuation length, A."""
            return _IMP_bff.PETTerm_get_rC(self)
    %}
}

%pythoncode %{
def _term_rates(term, first, second=None):
    """A term's rate constants as a numpy array."""
    return np.asarray(term.rate_constants(first, second if second is not None
                                          else States()), dtype=np.float64)
%}
