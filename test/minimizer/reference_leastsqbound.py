"""ChiSurf's old bounded least-squares optimiser, kept **only** as a reference.

Nothing in any product imports this. It is the implementation
`IMP::bff::FitMinimizer` was transcribed from -- scipy's `leastsq` under
ChiSurf's bounds transform -- and it is here so that the parity tests in
`test_minimizer.py` can still assert the port is 1:1: same answer, same
number of function evaluations, same evaluation *points*, same covariance,
from the same start.

**Why a frozen copy rather than the live module.** It used to be imported
from the sibling `chisurf` checkout, which was honest while ChiSurf still
ran on it. On 2026-09-01 ChiSurf's copy was deleted -- there is to be one
implementation of this algorithm in the stack and it is bff's -- and the
parity tests would then have skipped themselves silently. A skipped parity
test on a 1:1 port is worth nothing; it is the only thing that says the
answer did not move, and it is how the permuted-covariance defect
(`okf/log.md` 2026-09-01) was found. So the reference outlived the
implementation, which is the right way round.

Do not fix bugs in this file. If it disagrees with `FitMinimizer`, one of the
two is wrong and the question is *which* -- that is the whole point of
keeping it. Do not import it from `src/` or `pyext/` either; it is test
scaffolding, and `AGENTS.md`'s language rule does not make an exception for
a reference implementation that something ships against.

Original header: "Constrained multivariate least-squares optimization".
"""

import math
import warnings

import numpy as np

from scipy.optimize import _minpack, leastsq
import scipy
try:
    from scipy.optimize.minpack import _check_func
except ImportError:
    from scipy.optimize._minpack_py import _check_func

class OptimizationCancelled(Exception):
    """Signal that a least-squares optimization was cancelled by the caller.

    This exception is intended to be raised from a user-provided
    ``progress_callback`` (for example, when a GUI progress dialog's Cancel
    button is pressed). ``leastsqbound`` treats this exception specially and
    will propagate it to the caller instead of swallowing it together with
    other callback errors.
    """

    pass

def _internal2external_grad(
        xi,
        bounds
):
    """
Calculate the internal (unconstrained) to external (constained)
parameter gradiants.
"""
    grad = np.empty_like(xi)
    for i, (v, bound) in enumerate(zip(xi, bounds)):
        lower, upper = bound
        if lower is None and upper is None:  # No constraints
            grad[i] = 1.0
        elif upper is None:  # only lower bound
            grad[i] = v / np.sqrt(v * v + 1.)
        elif lower is None:  # only upper bound
            grad[i] = -v / np.sqrt(v * v + 1.)
        else:  # lower and upper bounds
            grad[i] = (upper - lower) * np.cos(v) / 2.
    return grad

def _internal2external_func(bounds):
    """
Make a function which converts between internal (unconstrained) and
external (constrained) parameters.
"""
    ls = [_internal2external_lambda(b) for b in bounds]

    def convert_i2e(xi):
        """Convert internal (unconstrained) parameters to external (constrained).

        Parameters
        ----------
        xi : np.ndarray
            Internal (unconstrained) parameter vector.

        Returns
        -------
        np.ndarray
            External (constrained) parameter vector.
        """
        xe = np.empty_like(xi)
        xe[:] = [l(p) for l, p in zip(ls, xi)]
        return xe

    return convert_i2e

def _is_unbounded(v) -> bool:
    """Whether a bound value means "no constraint".

    Callers express an absent bound either as ``None`` or as a non-finite value
    (``+/-inf``, ``nan``). Only ``None`` used to be recognised, so a parameter
    declared ``(-inf, inf)`` fell through to the *bounded* branch and was mapped
    with ``arcsin((x - lower) / (upper - lower) - 1)`` -> ``arcsin(inf/inf - 1)``
    -> **nan**, silently poisoning the optimiser's internal parameter vector.
    """
    return v is None or not np.isfinite(v)


def _internal2external_lambda(bound):
    """
Make a lambda function which converts a single internal (uncontrained)
parameter to a external (constrained) parameter.
"""
    lower, upper = bound
    lo_free, up_free = _is_unbounded(lower), _is_unbounded(upper)
    if lo_free and up_free:  # no constraints
        return lambda x: x
    elif up_free:  # only lower bound
        return lambda x: lower - 1. + np.sqrt(x * x + 1.)
    elif lo_free:  # only upper bound
        return lambda x: upper + 1. - np.sqrt(x * x + 1.)
    else:
        return lambda x: lower + ((upper - lower) / 2.) * (np.sin(x) + 1.)

def _external2internal_func(bounds):
    """
Make a function which converts between external (constrained) and
internal (unconstrained) parameters.
"""
    ls = [_external2internal_lambda(b) for b in bounds]

    def convert_e2i(xe):
        """Convert external (constrained) parameters to internal (unconstrained).

        Parameters
        ----------
        xe : np.ndarray
            External (constrained) parameter vector.

        Returns
        -------
        np.ndarray
            Internal (unconstrained) parameter vector.
        """
        xi = np.empty_like(xe)
        xi[:] = [l(p) for l, p in zip(ls, xe)]
        return xi

    return convert_e2i

def _external2internal_lambda(bound):
    """
Make a lambda function which converts an single external (constrained)
parameter to a internal (unconstrained) parameter.
"""
    lower, upper = bound
    lo_free, up_free = _is_unbounded(lower), _is_unbounded(upper)
    if lo_free and up_free:  # no constraints
        return lambda x: x
    elif up_free:  # only lower bound
        return lambda x: np.sqrt((x - lower + 1.) ** 2 - 1)
    elif lo_free:  # only upper bound
        return lambda x: np.sqrt((upper - x + 1.) ** 2 - 1)
    else:
        return lambda x: np.arcsin((2. * (x - lower) / (upper - lower)) - 1.)


#: Levenberg-Marquardt iterations a well-posed fit typically needs. Each costs
#: ``n + 1`` residual evaluations -- ``n`` for the forward-difference Jacobian and
#: one for the trial step -- so the expected budget is ``_EXPECTED_ITERATIONS *
#: (n + 1)``.
#:
#: An estimate only -- fits routinely run four to ten times longer, and
#: :func:`_reported_total` is built so that overrunning it degrades smoothly
#: instead of pinning the bar.
_EXPECTED_ITERATIONS = 6

#: How many multiples of the estimate the bar takes to cover its first half.
#:
#: This is the whole pacing decision. A bar that reads ``done = nfev / estimate``
#: saturates the instant the estimate is passed: on a real MFD fit -- 450
#: evaluations against an estimate of 42 -- it reached 98% after ten seconds and
#: then sat at 98-99% for **91% of the fit**. Reporting 450 times is no use when
#: 410 of those reports are the same number.
#:
#: So the fraction is asymptotic rather than linear (see :func:`_reported_total`):
#: it always moves, always increases, and never arrives. Two means the estimate
#: buys the first ~29% of the bar, which keeps that same fit visibly moving to
#: the end instead of dying at ten seconds.
_BAR_HALF_LIFE = 2.0

#: The most a *running* fit may report. A full bar is reserved for the
#: completion report, so "converged" stays distinguishable from "outran the
#: estimate and is still going" -- which is the difference between a bar the
#: user can trust and one they learn to ignore.
_MAX_RUNNING_RATIO = 0.99


def _expected_evaluations(n: int, maxfev: int) -> int:
    """Return an evaluation budget a progress bar can usefully report against.

    **MINPACK's ``200 * (n + 1)`` is a give-up limit, not an expectation.** Using
    it as the denominator made a four-parameter fit -- which converges in about
    fifty evaluations -- creep to five percent and then jump straight to done,
    which reads as a broken progress bar rather than a fast fit.

    The estimate here is what the algorithm actually costs: one Jacobian
    (``n`` evaluations) plus one trial step per iteration, over the number of
    iterations a well-posed problem needs. It is deliberately an *estimate*, and
    :func:`_reported_total` is what makes overrunning it harmless.

    Parameters
    ----------
    n : int
        Number of free parameters.
    maxfev : int
        Hard evaluation limit; the estimate is never larger than this.

    Returns
    -------
    int
    """
    limit = int(maxfev) if maxfev and maxfev > 0 else 200 * (int(n) + 1)
    return int(max(1, min(limit, _EXPECTED_ITERATIONS * (int(n) + 1))))


def _reported_total(nfev: int, expected: int) -> int:
    """Return the total to report ``nfev`` against, so the bar always moves.

    Nobody can know how many evaluations a fit needs -- that is the thing being
    discovered. A bar that divides by a *guess* is therefore wrong twice over: it
    finishes early when the guess is high, and it saturates when the guess is
    low. Saturating is much the worse of the two, because it happens on exactly
    the fits long enough for anyone to be watching. A real MFD fit takes 450
    evaluations against an estimate of 42, and a linear bar sat at 98-99% for
    **91% of its runtime** -- reporting 450 times, 410 of them the same number.

    So the reported fraction is asymptotic instead::

        done / total = 1 - 2 ** (-nfev / (_BAR_HALF_LIFE * expected))

    which is strictly increasing in ``nfev``, never reaches one, and degrades
    gracefully in both directions: a fit that beats the estimate ends low and
    snaps to full on completion (over in a blink, so nobody sees it), and a fit
    that takes ten times the estimate is still visibly moving at the end.

    The returned number is a genuine current estimate of the total -- ``nfev``
    divided by the fraction -- so ``(done, total)`` keeps meaning what it says
    and every consumer of the callback keeps working unchanged.

    Parameters
    ----------
    nfev : int
        Evaluations so far.
    expected : int
        The a-priori estimate, from :func:`_expected_evaluations`.

    Returns
    -------
    int
        Total to report against; always greater than ``nfev``.
    """
    scale = max(1.0, _BAR_HALF_LIFE * float(max(1, expected)))
    # The ceiling is folded into the curve rather than clamped on afterwards. A
    # clamp creates a *constant* region, and a constant fraction over integer
    # totals oscillates: ``ceil(594/0.99)`` is exactly 600 while ``ceil(595/0.99)``
    # is 602, so the reported fraction stepped backwards once every hundred
    # evaluations. Multiplying instead keeps the curve strictly increasing.
    fraction = _MAX_RUNNING_RATIO * (1.0 - 2.0 ** (-float(nfev) / scale))
    fraction = max(fraction, 1e-6)
    total = math.ceil(float(nfev) / fraction)
    return int(max(int(nfev) + 1, total))


def _report_progress(callback, nfev: int, total: int, chi2, chi2r) -> None:
    """Report one evaluation to *callback*, whatever signature it has.

    The extra ``chi2`` / ``chi2r`` keywords are a convenience for callers that
    want to show the objective in the label -- but the *documented* signature is
    ``callback(evaluated, total)``, and a callback written to it raises
    ``TypeError`` on the keywords. That landed in a blanket ``except Exception``
    and was swallowed, so a correctly written callback silently never fired and
    the progress bar sat at zero for the whole fit. Offer the extras, fall back
    to the documented call.

    Parameters
    ----------
    callback : callable
        The progress sink.
    nfev : int
        Evaluations completed.
    total : int
        Estimated evaluation budget.
    chi2, chi2r : float or None
        Objective values, when they could be computed.

    Raises
    ------
    OptimizationCancelled
        If the callback raises it -- that is how a caller aborts a fit, so it
        must not be treated as a callback failure.
    """
    try:
        try:
            callback(nfev, total, chi2=chi2, chi2r=chi2r)
        except TypeError:
            # Either the callback takes only (done, total), or it raised a
            # TypeError of its own; the retry distinguishes them by outcome.
            callback(nfev, total)
    except OptimizationCancelled:
        raise
    except Exception:
        # A broken progress bar must not take the optimisation down with it.
        pass



def _pivot_indices(ipvt) -> np.ndarray:
    """Return MINPACK's column permutation as 0-based indices.

    **SciPy changed the convention and this code did not.** MINPACK's ``ipvt``
    is 1-based in Fortran, and every caller therefore wrote ``ipvt - 1``;
    modern SciPy returns it already 0-based (1.18's ``_minpack_py.py`` builds
    its own ``cov_x`` with a bare ``perm = retval[1]['ipvt']``, and the array
    contains a literal ``0``). Subtracting one from that produces ``-1``,
    which numpy reads as *the last element* -- so the permutation stops being
    a permutation and the covariance comes back with its rows and columns
    cyclically shifted. Measured against ``scipy.optimize.leastsq``'s own
    ``cov_x`` on a four-parameter fit: every element wrong, by up to a factor
    of nine.

    This was latent rather than harmful only because
    :meth:`~chisurf.core.fitting.fit.Fit.update_error_estimates` takes its own
    finite-difference Jacobian and never reads this matrix. Anything that did
    read it -- and reusing it is the obvious way to remove a third of a fit's
    runtime -- would have shipped silently permuted error bars.

    Detected rather than assumed, so both conventions keep working: a 1-based
    permutation of ``n`` columns contains ``n`` and not ``0``.

    Parameters
    ----------
    ipvt : array_like
        MINPACK's permutation, in either convention.

    Returns
    -------
    numpy.ndarray
        0-based column indices.
    """
    ipvt = np.asarray(ipvt, dtype=int)
    return ipvt - 1 if ipvt.min() > 0 else ipvt


def leastsqbound(
        func, x0,
        args=(),
        bounds=None,
        Dfun=None,
        full_output=0,
        col_deriv=0,
        ftol=1.49012e-8,
        xtol=1.49012e-8,
        gtol=0.0,
        maxfev=0,
        epsfcn=0.0,
        factor=100,
        diag=None,
        progress_callback=None,
        progress_total=None
):
    """
Bounded minimization of the sum of squares of a set of equations.

::

x = arg min(sum(func(y)**2,axis=0))
y

Parameters
----------
func : callable
should take at least one (possibly length N vector) argument and
returns M floating point numbers.
x0 : ndarray
The starting estimate for the minimization.
args : tuple
Any extra arguments to func are placed in this tuple.
bounds : list
``(min, max)`` pairs for each element in ``x``, defining
the bounds on that parameter. Use None for one of ``min`` or
``max`` when there is no bound in that direction.
Dfun : callable
A function or method to compute the Jacobian of func with derivatives
across the rows. If this is None, the Jacobian will be estimated.
full_output : bool
non-zero to return all optional outputs.
col_deriv : bool
non-zero to specify that the Jacobian function computes derivatives
down the columns (faster, because there is no transpose operation).
ftol : float
Relative error desired in the sum of squares.
xtol : float
Relative error desired in the approximate solution.
gtol : float
Orthogonality desired between the function vector and the columns of
the Jacobian.
maxfev : int
The maximum number of calls to the function. If zero, then 100*(N+1) is
the maximum where N is the number of elements in x0.
epsfcn : float
A suitable step length for the forward-difference approximation of the
Jacobian (for Dfun=None). If epsfcn is less than the machine precision,
it is assumed that the relative errors in the functions are of the
order of the machine precision.
factor : float
A parameter determining the initial step bound
(``factor * || diag * x||``). Should be in interval ``(0.1, 100)``.
diag : sequence
N positive entries that serve as a scale factors for the variables.

Returns
-------
x : ndarray
The solution (or the result of the last iteration for an unsuccessful
call).
cov_x : ndarray
Uses the fjac and ipvt optional outputs to construct an
estimate of the jacobian around the solution. ``None`` if a
singular matrix encountered (indicates very flat curvature in
some direction). This matrix must be multiplied by the
residual standard deviation to get the covariance of the
parameter estimates -- see curve_fit.
infodict : dict
a dictionary of optional outputs with the key s::

- 'nfev' : the number of function calls
- 'fvec' : the function evaluated at the output
- 'fjac' : A permutation of the R matrix of a QR
factorization of the final approximate
Jacobian matrix, stored column wise.
Together with ipvt, the covariance of the
estimate can be approximated.
- 'ipvt' : an integer array of length N which defines
a permutation matrix, p, such that
fjac*p = q*r, where r is upper triangular
with diagonal elements of nonincreasing
magnitude. Column j of p is column ipvt(j)
of the identity matrix.
- 'qtf' : the vector (np.transpose(q) * fvec).

mesg : str
A string message giving information about the cause of failure.
ier : int
An integer flag. If it is equal to 1, 2, 3 or 4, the solution was
found. Otherwise, the solution was not found. In either case, the
optional output variable 'mesg' gives more information.

Notes
-----
"leastsq" is a wrapper around MINPACK's lmdif and lmder algorithms.

cov_x is a Jacobian approximation to the Hessian of the least squares
objective function.
This approximation assumes that the objective function is based on the
difference between some observed target data (ydata) and a (non-linear)
function of the parameters `f(xdata, params)` ::

func(params) = ydata - f(xdata, params)

so that the objective function is ::

min sum((ydata - f(xdata, params))**2, axis=0)
params

Contraints on the parameters are enforced using an internal parameter list
with appropiate transformations such that these internal parameters can be
optimized without constraints. The transfomation between a given internal
parameter, p_i, and a external parameter, p_e, are as follows:

With ``min`` and ``max`` bounds defined ::

p_i = np.arcsin((2 * (p_e - min) / (max - min)) - 1.)
p_e = min + ((max - min) / 2.) * (np.sin(p_i) + 1.)

With only ``max`` defined ::

p_i = np.sqrt((max - p_e + 1.)**2 - 1.)
p_e = max + 1. - np.sqrt(p_i**2 + 1.)

With only ``min`` defined ::

p_i = np.sqrt((p_e - min + 1.)**2 - 1.)
p_e = min - 1. + np.sqrt(p_i**2 + 1.)

These transfomations are used in the MINUIT package, and described in
detail in the section 1.3.1 of the MINUIT User's Guide.

To Do
-----
Currently the ``factor`` and ``diag`` parameters scale the
internal parameter list, but should scale the external parameter list.

The `qtf` vector in the infodic dictionary reflects internal parameter
list, it should be correct to reflect the external parameter list.

References
----------
* F. James and M. Winkler. MINUIT User's Guide, July 16, 2004.

"""
    # Optional progress tracking state (best-effort only). When no callback
    # is provided, these variables remain unused and the behavior matches
    # the original implementation.
    nfev = 0
    # The a-priori estimate; ``_reported_total`` turns it into a fraction that
    # keeps moving whether the fit beats it or runs many times past it.
    expected = progress_total if progress_total is not None else None

    def _compute_objective(residuals, f_args):
        """Return (chi2, chi2r) for the current residual vector.

        chi2 is the sum of squared residuals. For chi2r we try to infer a
        model instance from ``f_args`` that exposes ``n_points`` and
        ``n_free`` and, if successful, apply the usual
        ``chi2 / (n_points - n_free - 1)`` formula.
        """

        try:
            r = np.array(residuals, ndmin=1).ravel()
        except Exception:
            return None, None

        try:
            chi2_val = float(np.dot(r, r))
        except Exception:
            chi2_val = None

        chi2r_val = None
        if chi2_val is not None:
            model = None
            for arg in f_args:
                if hasattr(arg, "n_points") and hasattr(arg, "n_free"):
                    model = arg
                    break
            if model is not None:
                try:
                    dof = float(model.n_points - getattr(model, "n_free", 0) - 1)
                    if dof > 0:
                        chi2r_val = chi2_val / dof
                except Exception:
                    pass

        return chi2_val, chi2r_val

    # use leastsq if no bounds are present
    if bounds is None:
        if progress_callback is not None:
            # Estimate a total evaluation budget similar to MINPACK defaults
            # so that the callback can report a normalized progress fraction.
            x0_arr = np.array(x0, ndmin=1)
            n = len(x0_arr)
            if expected is None:
                expected = _expected_evaluations(n, maxfev)

            def _wrapped_func(x, *f_args):
                """Wrapper around the objective function with progress reporting.

                Calls the original ``func``, increments the call counter, and
                invokes ``progress_callback`` after each evaluation. Used when
                no bounds are given.

                Parameters
                ----------
                x : np.ndarray
                    Parameter vector.
                f_args : tuple
                    Additional arguments passed to ``func``.

                Returns
                -------
                np.ndarray
                    Residual vector from ``func``.
                """
                nonlocal nfev
                res = func(x, *f_args)
                nfev += 1
                total = _reported_total(nfev, expected)
                chi2_val, chi2r_val = _compute_objective(res, f_args)
                _report_progress(
                    progress_callback, nfev, total, chi2_val, chi2r_val
                )
                return res

            result = leastsq(
                _wrapped_func,
                x0,
                args,
                Dfun,
                full_output,
                col_deriv,
                ftol,
                xtol,
                gtol,
                maxfev,
                epsfcn,
                factor,
                diag,
            )
            # The estimate is deliberately generous, so a fit that converges
            # early would otherwise leave the bar stranded partway. Reporting
            # the budget as met at the end is not a fudge: the fit is finished.
            _report_progress(progress_callback, nfev, nfev, None, None)
            return result
        return leastsq(
            func,
            x0,
            args,
            Dfun,
            full_output,
            col_deriv,
            ftol,
            xtol,
            gtol,
            maxfev,
            epsfcn,
            factor,
            diag,
        )

    # create function which convert between internal and external parameters
    i2e = _internal2external_func(bounds)
    e2i = _external2internal_func(bounds)

    x0 = np.array(x0, ndmin=1)
    i0 = e2i(x0)
    n = len(x0)
    if len(bounds) != n:
        raise ValueError('length of x0 != length of bounds')
    if not isinstance(args, tuple):
        args = (args,)
    m = _check_func('leastsq', 'func', func, x0, args, n)[0]
    if n > m[0]:
        raise TypeError('Improper input: N=%s must not exceed M=%s' % (n, m))

    # define a wrapped func which accepts internal parameters, converts them
    # to external parameters and calls func. Progress reporting is layered on
    # top of this base wrapper when a callback is provided.
    def _base_wfunc(x, *f_args):
        """Base wrapper: convert internal params to external and call func.

        Parameters
        ----------
        x : np.ndarray
            Internal (unconstrained) parameter vector.
        f_args : tuple
            Additional arguments passed to ``func``.

        Returns
        -------
        np.ndarray
            Residual vector from ``func`` evaluated at external params.
        """
        return func(i2e(x), *f_args)

    if Dfun is None:
        if (maxfev == 0):
            maxfev = 200 * (n + 1)
        if progress_callback is not None and expected is None:
            expected = _expected_evaluations(n, maxfev)

        if progress_callback is not None:
            def wfunc(x, *f_args):
                """Wrapper with progress reporting when bounds are given.

                Converts internal params to external, calls ``_base_wfunc``,
                and reports progress via ``progress_callback``.

                Parameters
                ----------
                x : np.ndarray
                    Internal (unconstrained) parameter vector.
                f_args : tuple
                    Additional arguments passed to ``func``.

                Returns
                -------
                np.ndarray
                    Residual vector from ``func``.
                """
                nonlocal nfev
                res = _base_wfunc(x, *f_args)
                nfev += 1
                total = _reported_total(nfev, expected)
                chi2_val, chi2r_val = _compute_objective(res, f_args)
                _report_progress(
                    progress_callback, nfev, total, chi2_val, chi2r_val
                )
                return res
        else:
            wfunc = _base_wfunc

        retval = _minpack._lmdif(
            wfunc,
            i0,
            args,
            full_output,
            ftol,
            xtol,
            gtol,
            maxfev,
            epsfcn,
            factor,
            diag,
        )
        if progress_callback is not None:
            # See the unbounded branch: a converged fit fills its own bar.
            _report_progress(progress_callback, nfev, nfev, None, None)
    else:
        if col_deriv:
            _check_func('leastsq', 'Dfun', Dfun, x0, args, n, (n, m))
        else:
            _check_func('leastsq', 'Dfun', Dfun, x0, args, n, (m, n))
        if (maxfev == 0):
            maxfev = 100 * (n + 1)

        def wDfun(x, *args):  # wrapped Dfun
            """Wrapper around the Jacobian function.

            Converts internal parameters to external before calling ``Dfun``.

            Parameters
            ----------
            x : np.ndarray
                Internal (unconstrained) parameter vector.
            *args : tuple
                Additional arguments passed to ``Dfun``.

            Returns
            -------
            np.ndarray
                Jacobian matrix evaluated at external parameters.
            """
            return Dfun(i2e(x), *args)

        retval = _minpack._lmder(
            func,
            wDfun,
            i0,
            args,
            full_output,
            col_deriv,
            ftol,
            xtol,
            gtol,
            maxfev,
            factor,
            diag
        )

    errors = {0: ["Improper input parameters.", TypeError],
              1: ["Both actual and predicted relative reductions "
                  "in the sum of squares\n are at most %f" % ftol, None],
              2: ["The relative error between two consecutive "
                  "iterates is at most %f" % xtol, None],
              3: ["Both actual and predicted relative reductions in "
                  "the sum of squares\n are at most %f and the "
                  "relative error between two consecutive "
                  "iterates is at \n most %f" % (ftol, xtol), None],
              4: ["The cosine of the angle between func(x) and any "
                  "column of the\n Jacobian is at most %f in "
                  "absolute value" % gtol, None],
              5: ["Number of calls to function has reached "
                  "maxfev = %d." % maxfev, ValueError],
              6: ["ftol=%f is too small, no further reduction "
                  "in the sum of squares\n is possible.""" % ftol, ValueError],
              7: ["xtol=%f is too small, no further improvement in "
                  "the approximate\n solution is possible." % xtol, ValueError],
              8: ["gtol=%f is too small, func(x) is orthogonal to the "
                  "columns of\n the Jacobian to machine "
                  "precision." % gtol, ValueError],
              'unknown': ["Unknown error.", TypeError]}

    info = retval[-1]  # The FORTRAN return value

    if (info not in [1, 2, 3, 4] and not full_output):
        if info in [5, 6, 7, 8]:
            warnings.warn(errors[info][0], RuntimeWarning)
        else:
            try:
                raise errors[info][1](errors[info][0])
            except KeyError:
                raise errors['unknown'][1](errors['unknown'][0])

    mesg = errors[info][0]
    x = i2e(retval[0])  # internal params to external params

    if full_output:
        # convert fjac from internal params to external
        ipvt = _pivot_indices(retval[1]['ipvt'])
        grad = _internal2external_grad(retval[0], bounds)
        retval[1]['fjac'] = (retval[1]['fjac'].T / np.take(grad, ipvt)).T
        cov_x = None
        if info in [1, 2, 3, 4]:
            # from numpy.dual import pinv
            from numpy.linalg import LinAlgError, pinv

            perm = np.take(np.eye(n), ipvt, 0)
            r = np.triu(np.transpose(retval[1]['fjac'])[:n, :])
            R = np.dot(r, perm)
            try:
                cov_x = pinv(np.dot(np.transpose(R), R))
            except LinAlgError:
                pass
        return (x, cov_x) + retval[1:-1] + (mesg, info)
    else:
        return (x, info)
