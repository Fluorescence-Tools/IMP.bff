/*
 * The output contract as a C++ value, with the Python surface it had as a
 * dataclass.
 *
 * `LifetimeSpectrum` was 112 lines of Python holding two numpy arrays. Moving
 * it to C++ is not about the arithmetic -- the two averages are a dozen
 * multiply-adds. It is about what the object *is*: in `atom` and `core` the
 * Python is a shim because the domain objects are C++ and SWIG exposes them,
 * and no amount of kernel porting produces that shape while the objects
 * themselves are Python holding numpy arrays.
 *
 * What has to survive is the surface, because it is what the tests and every
 * caller use: attribute access rather than getters, `len()`, `repr()`, and a
 * `decay()` that returns the caller's shape.
 */

IMP_SWIG_VALUE(IMP::bff, LifetimeSpectrum, LifetimeSpectrums);

// The class has two constructors -- the real one and a default one, which the
// value-vector template needs -- and SWIG will not generate keyword arguments
// for an overloaded function. Callers say `LifetimeSpectrum(a, k, exact=False)`,
// so the shadow restores them. Same reason and same shape as the one on
// AVNetworkRestraint.
%feature("shadow") IMP::bff::LifetimeSpectrum::LifetimeSpectrum %{
def __init__(self, *args, **kwargs):
    if not args and not kwargs:
        _IMP_bff.LifetimeSpectrum_swiginit(self, _IMP_bff.new_LifetimeSpectrum())
        return
    names = ("amplitudes", "rate_constants", "exact")
    defaults = {"exact": True}
    if len(args) > len(names):
        raise TypeError("LifetimeSpectrum() takes at most %d positional "
                        "arguments (%d given)" % (len(names), len(args)))
    values = dict(zip(names, args))
    for k, v in kwargs.items():
        if k not in names:
            raise TypeError("LifetimeSpectrum() got an unexpected keyword "
                            "argument %r" % k)
        if k in values:
            raise TypeError("LifetimeSpectrum() got multiple values for "
                            "argument %r" % k)
        values[k] = v
    for k in names:
        if k not in values:
            if k not in defaults:
                raise TypeError("LifetimeSpectrum() missing required "
                                "argument %r" % k)
            values[k] = defaults[k]
    _IMP_bff.LifetimeSpectrum_swiginit(
        self, _IMP_bff.new_LifetimeSpectrum(*[values[k] for k in names]))
%}

// Read-only properties. `%attribute_np` wraps the getter in `np.array`, which
// would copy; these getters already publish a managed numpy view over the
// object's own buffer, so `%attribute_py` passes it through untouched.
%attribute_py(IMP::bff::LifetimeSpectrum, PyObject*, amplitudes, get_amplitudes);
%attribute_py(IMP::bff::LifetimeSpectrum, PyObject*, rate_constants, get_rate_constants);
%attribute_py(IMP::bff::LifetimeSpectrum, PyObject*, lifetimes, get_lifetimes);
%attribute_py(IMP::bff::LifetimeSpectrum, bool, exact, get_exact);
%attribute_py(IMP::bff::LifetimeSpectrum, int, n_species, get_n_species);
%attribute_py(IMP::bff::LifetimeSpectrum, double, total_amplitude, get_total_amplitude);
%attribute_py(IMP::bff::LifetimeSpectrum, double, species_averaged_lifetime,
              get_species_averaged_lifetime);
%attribute_py(IMP::bff::LifetimeSpectrum, double, intensity_averaged_lifetime,
              get_intensity_averaged_lifetime);

// `convolve_distance_with_k2_ratio` publishes its centres and weights as two
// managed views; the names are unique to it, so this binds only that call.
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** out_centres, int* n_centres)};
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** out_hist, int* n_hist)};

%include "IMP/bff/LifetimeSpectrum.h"

// `decay` publishes a flat managed view (ARGOUTVIEWM_ARRAY1) whose 1-D shape is
// the kernel's own; the caller's time axis shape is not restated here.
%extend IMP::bff::LifetimeSpectrum {
    int __len__() { return $self->get_n_species(); }
}
