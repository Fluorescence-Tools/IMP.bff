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

/*
 * Building a spectrum -- the reductions that turn a forward model's internals
 * into the contract the header states.
 *
 * `IMP.bff.observables` was a module for these three functions and a re-export
 * of the two names that had already become C++. What is left is argument
 * marshalling and a shape check, so it lives here rather than in a file of its
 * own.
 */
%pythoncode %{
def lifetime_spectrum_from_rates(rates, weights=None, exact=True):
    """One species per rate constant.

    The static limit taken literally: if a population of weight ``w_i`` decays
    at ``k_i`` and keeps that rate, then ``F(t) = sum w_i exp(-k_i t)`` is the
    decay, with no fitting and no approximation. Whether that premise holds is
    the caller's to know, and *exact* records the answer.

    :param rates: total deactivation rate per species, 1/ns.
    :param weights: population per species; uniform if omitted.
    :param exact: see :class:`LifetimeSpectrum`.
    """
    k = np.asarray(rates, dtype=np.float64).ravel()
    w = (np.ones_like(k) if weights is None
         else np.asarray(weights, dtype=np.float64).ravel())
    if w.shape != k.shape:
        raise ValueError(f"one weight per rate: {w.shape} against {k.shape}")
    return LifetimeSpectrum(w, k, exact=exact)


def rate_constants(terms, *participants, **kwargs):
    """Total deactivation rate per state, in 1/ns.

    A thin, named pass-through to
    :func:`IMP.bff.photophysics.total_rate`, kept here because *the rate
    constants are themselves an observable* -- one of the four this package
    emits. A caller that wants FRET rates rather than a decay should be able to
    ask for them without going through a spectrum.

    The channels add because they are parallel.
    """
    from IMP.bff.photophysics import total_rate
    return total_rate(terms, *participants, **kwargs)


def lifetime_spectrum_from_states(terms, *participants, weights=None,
                                  exact=True, **kwargs):
    """The spectrum of a state ensemble under a set of interaction terms.

    One species per state -- per accessible-volume point, per rotamer, per
    conformer -- carrying the summed rate of every channel acting on it. This is
    the reduction that connects :mod:`IMP.bff.representation` and
    :mod:`IMP.bff.photophysics` to an experiment-neutral answer, and it is
    representation-agnostic for the same reason the terms are: it consumes
    states.

    .. warning::
       This is the **static** limit. It is exact when each state holds its rate
       for the whole excited-state lifetime, and wrong when the dye reorganises
       fast enough to average over rates -- then the decay is not a sum of
       exponentials at all and the population has to be propagated instead
       (:class:`~IMP.bff.GridDiffusionSolver`, or the Brownian walk). Pass ``exact=False`` when using it outside that limit, so the
       spectrum says what it is.

    :param terms: the interaction terms to sum.
    :param participants: states, in the terms' arity order -- the donor's
        states first, then an acceptor's or the quencher atoms.
    :param weights: population per state. Taken from the first participant's
        ``weights`` when omitted, which is what an accessible volume's
        occupancy already is.
    """
    k = np.asarray(rate_constants(terms, *participants, **kwargs),
                   dtype=np.float64).ravel()
    if weights is None and participants:
        w = getattr(participants[0], "weights", None)
        weights = None if w is None else np.asarray(w, dtype=np.float64).ravel()
    return lifetime_spectrum_from_rates(k, weights, exact=exact)
%}
