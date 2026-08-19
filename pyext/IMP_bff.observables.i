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

%include "IMP/bff/LifetimeSpectrum.h"

%extend IMP::bff::LifetimeSpectrum {
    %pythoncode %{
        def __len__(self):
            return self.get_n_species()

        def decay(self, time):
            """``F(t) = sum_i a_i exp(-k_i t)`` on the caller's axis. **Unconvolved.**

            The reshape is the whole reason this is not the raw C++ method: a
            caller handing in a ``(3, 4)`` time array expects a ``(3, 4)``
            decay, and a numpy view over a flat buffer is flat. The kernel
            stays 1-D; the shape is the caller's, so restoring it belongs on
            this side of the boundary.
            """
            t = np.asarray(time, dtype=np.float64)
            out = _IMP_bff.LifetimeSpectrum_decay(
                self, np.ascontiguousarray(t.ravel()))
            return out.reshape(t.shape)
    %}
}
