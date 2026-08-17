"""``njit``/``jit`` that fall back to plain Python when numba is absent.

The numerical kernels in ``distance_metrics``, ``distributions``, ``polymer``
and ``av._kernels`` are written as numba-jitted loops. numba is not part of what
IMP brings, and IMP.bff ships through conda-forge as part of IMP where the
runtime dependency list is a public contract -- so it is used when present and
done without when not, rather than being required.

The cost of doing without is real and worth stating: these are explicit Python
loops, so an un-jitted kernel is orders of magnitude slower, not a few percent.
The fallback exists so that ``import IMP.bff`` works and the answers are the
same, not so that the slow path is a reasonable place to live. Anything doing
serious AV or polymer work should install numba.

Usage mirrors numba's, so the call sites read unchanged::

    from IMP.bff._jit import njit

    @njit
    def kernel(...): ...
"""

from __future__ import annotations

__all__ = ["njit", "jit", "prange", "get_num_threads", "HAS_NUMBA"]

try:
    import numba as _nb

    HAS_NUMBA = True
    njit = _nb.njit
    jit = _nb.jit
    prange = _nb.prange
    get_num_threads = _nb.get_num_threads

except ImportError:  # pragma: no cover - exercised only where numba is absent
    HAS_NUMBA = False

    #: ``prange`` degrades to ``range``: a ``parallel=True`` kernel that is not
    #: being compiled is just a loop, and the results are identical because
    #: every such kernel here is written race-free (each iteration owns its
    #: output slice).
    prange = range

    def get_num_threads():
        """One "thread" without numba -- the RNG-seeding loops then run once."""
        return 1

    def _passthrough(*args, **kwargs):
        """Accept every spelling numba's decorators accept, and do nothing.

        numba is used here both bare (``@njit``) and called
        (``@jit(nopython=True, nogil=True)``), so the shim has to cope with
        being handed the function directly or being handed options first.
        """
        if len(args) == 1 and not kwargs and callable(args[0]):
            return args[0]

        def decorator(func):
            return func

        return decorator

    njit = _passthrough
    jit = _passthrough
