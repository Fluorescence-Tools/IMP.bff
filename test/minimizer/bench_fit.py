"""A/B a whole ChiSurf fit across the three paths it can take.

Run it directly; it is a benchmark, not a test, so nothing here asserts.

    $E/bin/python test/minimizer/bench_fit.py

The three paths:

``scipy``
    ``leastsqbound`` -- MINPACK's ``lmdif`` in Fortran, calling back into
    Python for every residual. **This is history, not a path anything
    takes.** ChiSurf's copy was deleted on 2026-09-01 -- one implementation
    of the algorithm in the stack, and it is bff's -- so this row is measured
    against the frozen `reference_leastsqbound.py` beside this file, which
    exists to keep the parity tests honest and to keep this baseline
    meaningful.
``director``
    ``IMP.bff.Minimizer`` driving the *same* Python residual through a
    ``Node`` director. **This is the fallback that ships**, for every model
    the graph cannot represent. It was rejected as the fallback when it
    measured 1.54 ms against scipy's 1.33; timing `minimize` alone on the
    current build it is 1.11x on the parse fit and 1.06x on the decay --
    still slower, by a few per cent, and paid deliberately so that one
    implementation of the algorithm exists rather than three.

    **The whole-run rows below flatter it less than that.** A scipy run
    leaves no covariance behind, so its error estimate goes through
    `covariance_matrix` -- which is itself in C++ over the graph now --
    while a director run differences the *director*, in Python. Read the
    ratios as a whole-fit comparison, not as an optimiser one.
``graph``
    ``Expression -> ChiSquared -> Minimizer``. The parameters are ports the
    optimiser writes in C++, the curve is computed in C++, the data live in
    the node; nothing crosses the SWIG boundary per iteration.

A third table does the same for a **TCSPC lifetime fit** -- the model
anyone actually fits. Its curve is `IMP.bff.TcspcDecay`, which reconvolves
the lifetime spectrum with the measured response using tttrlib's own
kernels; there is no ``director`` row for it either, because before the node
existed a lifetime model simply fell back to scipy.

A second table does the same for a ``FitGroup`` -- which is what the GUI
actually builds, ``global_optimize_local_first`` shipping ``false`` so a
group is exactly one optimisation over ``GlobalFitModel``. Its graph is one
``Expression -> ChiSquared`` per member under a ``JointChiSquared``, with the
shared parameter as a ``Port`` link. There is no ``director`` row: before
this existed a group simply fell back to scipy, so scipy *is* the before.

A final table splits one ``fit.run()`` into optimise / error estimate /
rest and counts the Python ``update_model()`` calls inside it. That count is
the honest metric -- it does not move with the machine -- and it is what says
whether the fit stays in C++ *after* the optimiser stops. It used to be 7 for
a TCSPC fit, 5 of them rebuilding a Jacobian for the error bars; it is 2.

**Interleaved and min-of-many, deliberately.** A plain before/after on this
machine is useless -- load drifts by more than the effect (okf/log.md
2026-09-01 (5) learned that the expensive way). The three paths alternate
within one process and the best time of many is reported, which is the
arrangement that survived being wrong twice.

Expect roughly 1.26 / 1.26 / 0.50 ms at 512 points and three free
parameters, and for the group roughly 4.0 / 4.5 / 2.2 ms over four members.
The ratios are stable; the absolute numbers are not, and drop by a third if a
test suite is running alongside.
"""
import sys
import time
from pathlib import Path

import numpy as np

# ChiSurf is a sibling checkout (AGENTS.md), not an installed package.
for candidate in (Path(__file__).resolve().parents[3] / "chisurf",):
    if candidate.is_dir() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import chisurf as cs                                          # noqa: E402
import chisurf.core.data                                      # noqa: E402
import chisurf.core.fitting.fit as F                          # noqa: E402
import chisurf.core.fitting.minimizer as M                    # noqa: E402
import chisurf.core.math.optimization as O                    # noqa: E402
from chisurf.core.models.parse import ParseModel              # noqa: E402
from chisurf.core.models.tcspc.lifetime import LifetimeModel  # noqa: E402

N = 512


def make_fit(n=N, seed=3):
    x = np.linspace(0.1, 20.0, n)
    y = 2.5 * np.exp(-x / 3.1) + 0.4
    y = y + np.random.default_rng(seed).normal(0, 0.02, n)
    data = cs.core.data.DataCurve(x=x, y=y, ey=np.full(n, 0.02))
    fit = F.Fit(model_class=ParseModel, data=data)
    fit.model.func = "a*exp(-x/t)+b"
    for name, value in (("a", 1.0), ("t", 1.0), ("b", 0.0)):
        fit.model.parameters_all_dict[name].value = value
    fit.xmin, fit.xmax = 0, n - 1
    return fit


class forced:
    """Force ``fit.run()`` down one of the three paths, for the length of a block.

    ``scipy`` is the only one that needs explaining. ChiSurf's own bounded
    Levenberg-Marquardt was deleted on 2026-09-01 -- one implementation of
    this algorithm in the stack, and it is bff's -- so the row is measured
    against the *frozen* copy in `reference_leastsqbound.py`, which sits
    beside this file and is imported by nothing else. It is kept as a
    baseline and as the parity reference, not as a path anything takes.
    """

    def __init__(self, mode):
        self.mode = mode

    def __enter__(self):
        self._minimize, self._graph = M.minimize, M.graph_objective
        if self.mode == "scipy":
            from reference_leastsqbound import leastsqbound

            def fallback(func, x0, args=(), bounds=None,
                         progress_callback=None, n_free=0, fit=None,
                         model=None, **kw):
                return leastsqbound(func, x0, args=args, bounds=bounds,
                                    progress_callback=progress_callback, **kw)
            M.minimize = fallback
        elif self.mode == "director":
            # A refused graph is what sends `minimize` to its fallback.
            M.graph_objective = lambda *a, **k: None
        return self

    def __exit__(self, *exc):
        M.minimize, M.graph_objective = self._minimize, self._graph
        return False


def _best_of(make, mode, reps, run=lambda fit: fit.run()):
    """Best of *reps* whole ``run()`` calls down one path."""
    with forced(mode):
        best, result = float("inf"), None
        for _ in range(reps):
            fit = make()
            start = time.perf_counter()
            run(fit)
            best = min(best, time.perf_counter() - start)
            model = getattr(fit, "_model", None) or fit.model
            result = (fit.chi2r, [p.value for p in model.parameters])
        return best, result


def timed(mode, reps):
    """Best of *reps* whole ``fit.run()`` calls down one path."""
    return _best_of(make_fit, mode, reps)


MEMBERS = 4


def make_group(n_members=MEMBERS, n=N, seed=5):
    """A group of decays sharing one lifetime -- the shape the GUI builds.

    Each member keeps its own amplitude and its own data; the lifetime is one
    number for the whole group, linked the way ChiSurf links it (the
    follower's ``Parameter.link``, which is a ``Port`` link underneath). So
    the free vector is one amplitude per member plus the shared lifetime.
    """
    rng = np.random.default_rng(seed)
    x = np.linspace(0.1, 20.0, n)
    curves = []
    for k in range(n_members):
        y = (2.0 + 0.5 * k) * np.exp(-x / 3.1) + rng.normal(0, 0.02, n)
        curves.append(cs.core.data.DataCurve(x=x, y=y, ey=np.full(n, 0.02)))
    group = F.FitGroup(data=cs.core.data.DataGroup(curves),
                       model_class=ParseModel)
    for member in group:
        member.model.func = "a*exp(-x/t)"
        member.xmin, member.xmax = 0, n - 1
        member.model.find_parameters()
    master = group[0].model.parameters_all_dict["t"]
    for member in list(group)[1:]:
        member.model.parameters_all_dict["t"].link = master
    for member in group:
        member.model.find_parameters()
    group._model.find_parameters()
    for member in group:
        member.model.parameters_all_dict["a"].value = 1.0
    master.value = 1.0
    return group


def timed_group(mode, reps):
    """Best of *reps* whole ``FitGroup.run()`` calls down one path."""
    return _best_of(make_group, mode, reps,
                    run=lambda fit: fit.run(local_first=False))


def main_group(rounds=3, reps=3):
    best = {}
    for _ in range(rounds):
        for mode in ("scipy", "director", "graph"):
            elapsed, result = timed_group(mode, reps)
            if mode not in best or elapsed < best[mode][0]:
                best[mode] = (elapsed, result)
    print("")
    print("FitGroup, %d members, %d free parameters:" % (MEMBERS, MEMBERS + 1))
    for mode in ("scipy", "director", "graph"):
        elapsed, (chi2r, values) = best[mode]
        print("%-9s %.3f ms   chi2r %.5f   x %s"
              % (mode, elapsed * 1e3, chi2r, np.round(values, 6)))
    print("graph vs director %.2fx   (historical: vs scipy %.2fx)"
          % (best["director"][0] / best["graph"][0],
             best["scipy"][0] / best["graph"][0]))
    # A group is one optimisation over one objective, so the graph must land
    # where the numpy path lands -- including on the shared parameter, which
    # is the only thing making it a group rather than N fits.
    for mode in ("director", "graph"):
        np.testing.assert_allclose(best[mode][1][1], best["scipy"][1][1],
                                   rtol=1e-5)
    print("all three paths agree to 1e-5")


def main(rounds=4, reps=4):
    best = {}
    for _ in range(rounds):
        for mode in ("scipy", "director", "graph"):
            elapsed, result = timed(mode, reps)
            if mode not in best or elapsed < best[mode][0]:
                best[mode] = (elapsed, result)
    for mode in ("scipy", "director", "graph"):
        elapsed, (chi2r, values) = best[mode]
        print("%-9s %.3f ms   chi2r %.5f   x %s"
              % (mode, elapsed * 1e3, chi2r, np.round(values, 6)))
    print("graph vs director %.2fx   (historical: vs scipy %.2fx)"
          % (best["director"][0] / best["graph"][0],
             best["scipy"][0] / best["graph"][0]))
    # The answer must not depend on the path; a speed-up that moved the
    # minimum would not be one.
    reference = best["scipy"][1][1]
    for mode in ("director", "graph"):
        np.testing.assert_allclose(best[mode][1][1], reference, rtol=1e-6)
    print("all three paths agree to 1e-6")


DECAY_N = 512
DECAY_DT = 0.032
DECAY_REP = 80.0


def make_decay_fit(n=DECAY_N, dt=DECAY_DT, rep_rate=DECAY_REP, seed=7):
    """A single-exponential TCSPC fit of Poisson counts, ChiSurf's defaults.

    Free: the scatter fraction, the constant background, the lifetime and the
    timeshift -- which is the set a default `LifetimeModel` hands the
    optimiser. The amplitude is redundant (`normalize_amplitudes`) and `n0`
    is autoscaled, i.e. computed from the data rather than fitted.
    """
    import tttrlib
    x = np.arange(n) * dt
    irf = 1000.0 * np.exp(-0.5 * ((x - 1.0) / 0.08) ** 2)
    irf_y = irf / irf.sum()
    clean = np.zeros(n)
    tttrlib.fconv_per_cs(clean, irf_y, np.array([1.0, 3.1]),
                         1000.0 / rep_rate, n - 1, n - 1, dt)
    clean = clean / clean.max() * 4000.0 + 2.0
    y = np.random.default_rng(seed).poisson(clean).astype(float)
    data = cs.core.data.DataCurve(x=x, y=y, ey=np.sqrt(np.maximum(y, 1.0)))
    fit = F.Fit(model_class=LifetimeModel, data=data)
    fit.xmin, fit.xmax = 0, n
    model = fit.model
    model.convolve._irf = cs.core.curve.Curve(x=x, y=irf.copy())
    model.convolve.dt = dt
    model.convolve.rep_rate = rep_rate
    model.convolve.stop = n * dt
    model.find_parameters()
    model.parameters_all_dict["tL1"].value = 4.0
    return fit


def timed_decay(mode, reps):
    """Best of *reps* whole ``fit.run()`` calls down one path."""
    return _best_of(make_decay_fit, mode, reps)


def main_decay(rounds=3, reps=3):
    best = {}
    for _ in range(rounds):
        for mode in ("scipy", "director", "graph"):
            elapsed, result = timed_decay(mode, reps)
            if mode not in best or elapsed < best[mode][0]:
                best[mode] = (elapsed, result)
    print("")
    print("TCSPC lifetime fit, %d channels, 4 free parameters:" % DECAY_N)
    for mode in ("scipy", "director", "graph"):
        elapsed, (chi2r, values) = best[mode]
        print("%-9s %.3f ms   chi2r %.5f   x %s"
              % (mode, elapsed * 1e3, chi2r, np.round(values, 6)))
    print("graph vs director %.2fx   (historical: vs scipy %.2fx)"
          % (best["director"][0] / best["graph"][0],
             best["scipy"][0] / best["graph"][0]))
    # Compared against the uncertainty rather than by a blanket relative
    # tolerance: the background here is determined to about 0.8, so the two
    # optimisers stopping 1e-4 apart on it is agreement.
    b = np.asarray(best["scipy"][1][1])
    scale = np.maximum(np.abs(b), 1e-3 * np.max(np.abs(b)))
    for mode in ("director", "graph"):
        a = np.asarray(best[mode][1][1])
        assert np.all(np.abs(a - b) < 1e-2 * scale), (mode, a, b)
    print("all three paths agree to well inside the uncertainty")


def make_polarised_fit(polarization="vv", rho=2.5, seed=13, **kwargs):
    """The same instrument, one node further upstream.

    A VV decay is `LifetimeSpectrumNode -> AnisotropySpectrum ->
    TcspcDecay -> ChiSquared`: the polarisation is a *spectrum* transform
    (the product of two sums of exponentials is a sum of exponentials), so
    it costs one more node and no more crossings. Before it existed, a
    polarised model was refused outright and fell back to numpy, so scipy
    *is* the before.

    The data are generated from the polarised model rather than reused from
    the magic-angle fixture: a rotational correlation time fitted to a decay
    that carries no depolarisation is a flat direction, and a benchmark on a
    flat direction measures how long two optimisers take to stop, not how
    long they take to agree.
    """
    fit = make_decay_fit(**kwargs)
    model = fit.model
    model.anisotropy.polarization_type = polarization
    model.anisotropy.add_rotation(b=0.2, rho=rho)
    # Redundant against `r0` -- the amplitudes are normalised to it -- so
    # fitting it is fitting a flat direction.
    for parameter in model.anisotropy._bs:
        parameter.fixed = True
    model.anisotropy._g.value = 1.3
    model.find_parameters()
    model.update_model()
    clean = np.maximum(np.asarray(model.y, dtype=float), 1e-9)
    counts = np.random.default_rng(seed).poisson(clean).astype(float)
    fit.data.y = counts
    fit.data.ey = np.sqrt(np.maximum(counts, 1.0))
    model.parameters_all_dict["tL1"].value = 4.0
    return fit


def timed_polarised(mode, reps):
    """Best of *reps* whole ``fit.run()`` calls down one path."""
    return _best_of(make_polarised_fit, mode, reps)


def main_polarised(rounds=3, reps=3):
    best = {}
    for _ in range(rounds):
        for mode in ("scipy", "director", "graph"):
            elapsed, result = timed_polarised(mode, reps)
            if mode not in best or elapsed < best[mode][0]:
                best[mode] = (elapsed, result)
    print("")
    print("VV TCSPC fit, %d channels, one rotation, 5 free parameters:"
          % DECAY_N)
    for mode in ("scipy", "director", "graph"):
        elapsed, (chi2r, values) = best[mode]
        print("%-9s %.3f ms   chi2r %.5f   x %s"
              % (mode, elapsed * 1e3, chi2r, np.round(values, 6)))
    print("graph vs director %.2fx   (historical: vs scipy %.2fx)"
          % (best["director"][0] / best["graph"][0],
             best["scipy"][0] / best["graph"][0]))
    b = np.asarray(best["scipy"][1][1])
    scale = np.maximum(np.abs(b), 1e-3 * np.max(np.abs(b)))
    for mode in ("director", "graph"):
        a = np.asarray(best[mode][1][1])
        assert np.all(np.abs(a - b) < 1e-2 * scale), (mode, a, b)
    print("all three paths agree to well inside the uncertainty")


def make_fret_fit(distance=45.0, sigma=6.0, x_donly=0.2, seed=17, **kwargs):
    """A Gaussian-distance FRET decay -- the chain, end to end.

    `GaussianDistances -> FretSpectrum -> TcspcDecay -> ChiSquared`: three
    producers and an instrument, and the caller crosses once. The numpy path
    rebuilds the whole of that per iteration -- a 96-point distance
    distribution, 96 transfer rates, their product with the donor spectrum
    and the donor-only mixture -- which is why this is the widest gap of the
    four tables.
    """
    from chisurf.core.models.tcspc.fret import GaussianModel
    fit = make_decay_fit(**kwargs)
    model_data = fit.data
    fit = F.Fit(model_class=GaussianModel, data=model_data)
    fit.xmin, fit.xmax = 0, DECAY_N
    model = fit.model
    x = np.arange(DECAY_N) * DECAY_DT
    irf = 1000.0 * np.exp(-0.5 * ((x - 1.0) / 0.08) ** 2)
    model.convolve._irf = cs.core.curve.Curve(x=x, y=irf)
    model.convolve.dt = DECAY_DT
    model.convolve.rep_rate = DECAY_REP
    model.convolve.stop = DECAY_N * DECAY_DT
    model.fret_parameters.xDOnly = x_donly
    model.lifetimes._lifetimes[0].value = 4.0
    model.lifetimes._lifetimes[0].fixed = True
    model.gaussians._gaussianMeans[0].value = distance
    model.gaussians._gaussianSigma[0].value = sigma
    model.gaussians._gaussianSigma[0].fixed = True
    # Redundant with one component: the weights are normalised to sum to one.
    for parameter in model.gaussians._gaussianAmplitudes:
        parameter.fixed = True
    model.find_parameters()
    model.update_model()
    clean = np.maximum(np.asarray(model.y, dtype=float), 1e-9)
    counts = np.random.default_rng(seed).poisson(clean).astype(float)
    fit.data.y = counts
    fit.data.ey = np.sqrt(np.maximum(counts, 1.0))
    model.gaussians._gaussianMeans[0].value = 55.0
    return fit


def timed_fret(mode, reps):
    """Best of *reps* whole ``fit.run()`` calls down one path."""
    return _best_of(make_fret_fit, mode, reps)


def main_fret(rounds=3, reps=3):
    best = {}
    for _ in range(rounds):
        for mode in ("scipy", "director", "graph"):
            elapsed, result = timed_fret(mode, reps)
            if mode not in best or elapsed < best[mode][0]:
                best[mode] = (elapsed, result)
    print("")
    print("FRET fit (Gaussian distances), %d channels, %d-point "
          "distribution:" % (DECAY_N, 96))
    for mode in ("scipy", "director", "graph"):
        elapsed, (chi2r, values) = best[mode]
        print("%-9s %.3f ms   chi2r %.5f   x %s"
              % (mode, elapsed * 1e3, chi2r, np.round(values, 6)))
    print("graph vs director %.2fx   (historical: vs scipy %.2fx)"
          % (best["director"][0] / best["graph"][0],
             best["scipy"][0] / best["graph"][0]))
    b = np.asarray(best["scipy"][1][1])
    scale = np.maximum(np.abs(b), 1e-3 * np.max(np.abs(b)))
    for mode in ("director", "graph"):
        a = np.asarray(best[mode][1][1])
        assert np.all(np.abs(a - b) < 5e-2 * scale), (mode, a, b)
    print("all three paths agree to well inside the uncertainty")


# ------------------------------------------- where a fit's time actually goes

FIXTURES = (
    ("parse", make_fit),
    ("tcspc", make_decay_fit),
    ("tcspc VV", make_polarised_fit),
    ("FRET", make_fret_fit),
)


def split(make, reps=6):
    """One ``fit.run()``, timed in three parts, best-of-*reps*.

    The parts are the optimisation (:func:`minimize`), the error estimate
    (:meth:`Fit.update_error_estimates`) and everything else -- the parameter
    discovery, the model update, the result snapshot. Best-of rather than
    mean because load on this laptop drifts by more than the effect.
    """
    best = None
    for _ in range(reps):
        fit = make()
        timing = {"minimize": 0.0, "errors": 0.0}
        original_minimize = M.minimize
        original_errors = F.Fit.update_error_estimates

        def timed_minimize(*a, **kw):
            start = time.perf_counter()
            try:
                return original_minimize(*a, **kw)
            finally:
                timing["minimize"] += time.perf_counter() - start

        def timed_errors(self):
            start = time.perf_counter()
            try:
                return original_errors(self)
            finally:
                timing["errors"] += time.perf_counter() - start

        M.minimize = timed_minimize
        F.Fit.update_error_estimates = timed_errors
        try:
            start = time.perf_counter()
            fit.run()
            total = time.perf_counter() - start
        finally:
            M.minimize = original_minimize
            F.Fit.update_error_estimates = original_errors
        if best is None or total < best[0]:
            best = (total, timing["minimize"], timing["errors"])
    return best


def evaluations(make):
    """Python ``update_model()`` calls in one ``run()``, and where they fall.

    **The honest metric.** The clock moves with the machine; this does not.
    A model evaluation inside the error estimate is a crossing the graph was
    built to remove, and counting them says whether it was removed -- which
    is a different question from whether the fit got faster.

    Also records whether a C++ covariance was *offered* to
    `update_error_estimates`, separately from whether it was used: the guard
    on MINPACK's own matrix is applied when the stash is written, not when it
    is read, so the two were once different answers.
    """
    fit = make()
    model_class = type(fit.model)
    original_update = model_class.update_model
    original_errors = F.Fit.update_error_estimates
    state = {"total": 0, "errors": 0, "inside": False, "offered": None}

    def counting(self, *a, **kw):
        state["total"] += 1
        if state["inside"]:
            state["errors"] += 1
        return original_update(self, *a, **kw)

    def watched(self):
        state["offered"] = "_cpp_covariance" in self.__dict__
        state["inside"] = True
        try:
            return original_errors(self)
        finally:
            state["inside"] = False

    model_class.update_model = counting
    F.Fit.update_error_estimates = watched
    try:
        fit.run()
    finally:
        model_class.update_model = original_update
        F.Fit.update_error_estimates = original_errors
    return state


def main_split():
    """The table that says whether the fit still leaves C++ once it stops.

    A fit that optimises entirely in C++ and then calls the Python model
    ``p + 1`` more times to rebuild a Jacobian for its error bars has handed
    a third of the saving back. Before 2026-09-01 that is what every TCSPC
    fit did: 32% of the run, five model evaluations, because MINPACK's own
    covariance was *correctly* refused -- it differences at a step relative
    to the parameter, and a scatter fraction of 1.5e-5 beside a lifetime of
    3.15 gets a column of round-off. The remedy is not a better matrix in
    numpy but the same differences taken over the graph, at a step with an
    absolute floor.

    So the number to watch is the last column: a `run()` should evaluate the
    Python model **twice**, whatever the model is, and none of it for the
    covariance.
    """
    print("")
    print("Where a fit's time goes, and what it costs in model evaluations:")
    print("%-10s %8s %9s %9s %7s  %s"
          % ("fit", "total", "optimise", "errors", "rest", "update_model()"))
    for name, make in FIXTURES:
        total, minimise, errors = split(make)
        counted = evaluations(make)
        rest = total - minimise - errors
        print("%-10s %6.2f ms %8.0f%% %8.0f%% %6.0f%%  %d total, %d in errors"
              % (name, total * 1e3, 100 * minimise / total,
                 100 * errors / total, 100 * rest / total,
                 counted["total"], counted["errors"]))
    print("(a C++ covariance was offered in every row; the question is "
          "whether the error estimate still evaluated the Python model)")


if __name__ == "__main__":
    main()
    main_group()
    main_decay()
    main_polarised()
    main_fret()
    main_split()
