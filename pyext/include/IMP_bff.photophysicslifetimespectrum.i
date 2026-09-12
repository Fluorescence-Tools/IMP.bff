/*
 * The output contract as a C++ value.
 *
 * `PhotophysicsLifetimeSpectrum` is a domain object, not a pair of numpy arrays behind an
 * interface: it is C++ for the same reason `atom` and `core` values are, and
 * SWIG gives it attribute access, `len()`, `repr()` and a `decay()` that
 * returns the caller's shape.
 */

IMP_SWIG_VALUE(IMP::bff, PhotophysicsLifetimeSpectrum, PhotophysicsLifetimeSpectrums);

// Keyword arguments come from SWIG, not from a hand-written shadow: the
// constructor is a single declaration with default arguments now (an
// overloaded one cannot carry them), so %feature("kwargs") generates
// PhotophysicsLifetimeSpectrum(amplitudes, rate_constants, exact=True) and the 28 lines
// of Python that re-derived names, defaults and error messages are gone.
%feature("kwargs") IMP::bff::PhotophysicsLifetimeSpectrum::PhotophysicsLifetimeSpectrum;

// Read-only properties. `%attribute_np` wraps the getter in `np.array`, which
// would copy; these getters already publish a managed numpy view over the
// object's own buffer, so `%attribute_py` passes it through untouched.
%attribute_py(IMP::bff::PhotophysicsLifetimeSpectrum, PyObject*, amplitudes, get_amplitudes);
%attribute_py(IMP::bff::PhotophysicsLifetimeSpectrum, PyObject*, rate_constants, get_rate_constants);
%attribute_py(IMP::bff::PhotophysicsLifetimeSpectrum, PyObject*, lifetimes, get_lifetimes);
%attribute_py(IMP::bff::PhotophysicsLifetimeSpectrum, bool, exact, get_exact);
%attribute_py(IMP::bff::PhotophysicsLifetimeSpectrum, int, n_species, get_n_species);
%attribute_py(IMP::bff::PhotophysicsLifetimeSpectrum, double, total_amplitude, get_total_amplitude);
%attribute_py(IMP::bff::PhotophysicsLifetimeSpectrum, double, species_averaged_lifetime,
              get_species_averaged_lifetime);
%attribute_py(IMP::bff::PhotophysicsLifetimeSpectrum, double, intensity_averaged_lifetime,
              get_intensity_averaged_lifetime);

// `convolve_distance_with_k2_ratio` publishes its centres and weights as two
// managed views; the names are unique to it, so this binds only that call.
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** out_centres, int* n_centres)};
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** out_hist, int* n_hist)};

%include "IMP/bff/PhotophysicsLifetimeSpectrum.h"

// `decay` publishes a flat managed view (ARGOUTVIEWM_ARRAY1) whose 1-D shape is
// the kernel's own; the caller's time axis shape is not restated here.
%extend IMP::bff::PhotophysicsLifetimeSpectrum {
    int __len__() { return $self->get_n_species(); }
}
