"""How much compute each ChiSurf model spends, and whether it spends it in bff.

Run it directly; it is a benchmark, not a test, so nothing here asserts.

    $E/bin/python test/minimizer/bench_models.py

`census_models.py` answers *whether* a model becomes a C++ graph. That is
half of what is needed to decide what to port next, and on its own it is the
misleading half: it ranks a model nobody fits the same as the one every
session runs. This script measures the other half -- what one evaluation of
the model curve costs -- so the order of work is set by measurement rather
than by which class looked interesting.

**Cost alone is the trap this script exists to avoid, and it caught one.**
The most expensive model here by a factor of three, `MaxEntLifetimeModel` at
~190 ms a curve, has *nothing to move*: 97% of that is already inside
`tttrlib.solve_tcspc_mem_lifetime`, in C++, where the layering rule
(`AGENTS.md`: photons/curves -> tttrlib) says it belongs. Ranking by time
alone would have sent a session to port it. So each row is also profiled and
split three ways:

``native``
    inside a `tttrlib` or `IMP` extension call. **Already where it should
    be** -- porting this buys nothing.
``numpy``
    inside a numpy built-in. Vectorised, so it is not interpreter overhead --
    but the arrays live in Python and are marshalled to reach it, which is
    the half of the standing rule (*the data stay where the computation is*)
    that a "it is already vectorised" reading misses.
``py``
    interpreter-level, in ChiSurf's own `.py` files.

`numpy + py` is the work; `native` is done. The ranking column is
`us/LM-iter` restricted to that share.

**Two things the split gets wrong, both in the same direction.** It reads a
function by the file it is defined in, so a numpy routine whose *wrapper* is
Python -- `np.linalg.inv` is the one that matters here, and `np.linalg.solve`
and friends behave the same -- is charged to `py` even though the work is in
LAPACK. `DeerTikhonovModel` is the row this distorts: it reads 99% `py` and a
real share of that is BLAS. Treat a high `py%` on a linear-algebra model as a
question, not a verdict, and profile it before porting. The `native` column
has no such problem in the other direction: `tttrlib`/`IMP` calls are
extension objects and are always attributed correctly, so a row that reads
`native` really is done.

`DeerTikhonovModel` is also **the noisiest row in the table** -- it has been
seen between 46 and 232 ms a curve on an otherwise idle machine, without
growing within a process. Its rank is meaningful; its absolute number is
not.

**The unit that matters is not one `update_model()`.** Levenberg-Marquardt
differences the objective once per free parameter per iteration, so a model
with `n` free parameters costs `(n + 1)` evaluations per iteration, and the
column that decides anything is that product. A cheap curve with fifteen
free parameters can outweigh an expensive one with two.

The fixture is `census_models.py`'s, deliberately the same one: a model that
cannot be built there is not measurable here either, and the two tables are
meant to be read side by side.

**min-of-many, interleaved with nothing.** Unlike `bench_fit.py` there is no
A/B here -- each row is measured on its own -- so the defence against a
drifting machine is only the minimum. Read the ranks, not the absolute
microseconds; the ranks are what survive a loaded laptop.
"""
import cProfile
import pstats
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import census_models as C                                # noqa: E402
import chisurf.core.fitting.minimizer as M               # noqa: E402

REPEATS = 7
BATCH = 5


def time_update(model):
    """Best per-call time of `model.update_model()`, in microseconds."""
    model.update_model()                       # warm: first call builds caches
    best = float("inf")
    for _ in range(REPEATS):
        start = time.perf_counter()
        for _ in range(BATCH):
            model.update_model()
        best = min(best, (time.perf_counter() - start) / BATCH)
    return best * 1e6


def where_it_runs(model):
    """Fraction of one `update_model()` spent native / in numpy / in Python.

    `cProfile`'s own overhead inflates the interpreter share, so these are
    proportions to rank by, not a budget to quote. What they are reliable
    for is the distinction that matters: a call into a `tttrlib` or `IMP`
    extension is attributed to that extension, so a model that is already a
    C++ kernel behind a thin wrapper reads as `native` no matter how long
    the kernel takes.
    """
    profiler = cProfile.Profile()
    profiler.enable()
    for _ in range(3):
        model.update_model()
    profiler.disable()
    buckets = {"native": 0.0, "numpy": 0.0, "py": 0.0}
    for (filename, _, funcname), (_, _, tottime, _, _) in \
            pstats.Stats(profiler).stats.items():
        if filename.endswith(".py"):
            bucket = "py"
        elif "tttrlib" in funcname or "IMP" in funcname or "_imp" in funcname:
            bucket = "native"
        elif "numpy" in funcname or "scipy" in funcname:
            bucket = "numpy"
        else:
            bucket = "py"
        buckets[bucket] += tottime
    total = sum(buckets.values()) or 1.0
    return {k: v / total for k, v in buckets.items()}


def in_bff(fit, model):
    """`census_models.verdict` without re-timing it: yes / no / error."""
    try:
        built = M.graph_objective(fit, model)
    except Exception as exc:
        return "raised %s" % type(exc).__name__
    return "yes" if built is not None else "no"


def row(name, model_class):
    try:
        fit = C.make_fit(model_class, "vm")
    except Exception as exc:
        return name, None, None, None, "%s" % type(exc).__name__
    if fit is None:
        return name, None, None, None, "-"
    model = fit.model
    free = len(getattr(model, "parameters", []))
    try:
        micros = time_update(model)
        share = where_it_runs(model)
    except Exception as exc:
        return name, None, free, None, "update raised %s" % type(exc).__name__
    return name, micros, free, (in_bff(fit, model), share), None


def main():
    classes = C.model_classes()
    rows = []
    for name in sorted(classes):
        short = name.replace("chisurf.core.models.", "")
        rows.append(row(short, classes[name]))

    measured = [r for r in rows if r[1] is not None]
    skipped = [r for r in rows if r[1] is None]

    def movable(r):
        """us per LM iteration that is *not* already in a C++ extension.

        One curve per free parameter per iteration, plus the undifferenced
        one -- times the share that is not `native`. This is the ranking,
        and it is the only number here that answers "port what next".
        """
        share = r[3][1]
        return r[1] * (r[2] + 1) * (share["numpy"] + share["py"])

    measured.sort(key=lambda r: -movable(r))

    width = max(len(r[0]) for r in measured)
    print("%-*s  %9s %4s %11s %11s  %5s %5s %5s  %s"
          % (width, "model", "us/curve", "free", "us/LM-iter", "movable",
             "nat%", "np%", "py%", "graph"))
    print("-" * (width + 66))
    movable_total = 0.0
    for r in measured:
        name, micros, free, (bff, share), _ = r
        movable_total += movable(r)
        print("%-*s  %9.1f %4d %11.1f %11.1f  %5.0f %5.0f %5.0f  %s"
              % (width, name, micros, free, micros * (free + 1), movable(r),
                 100 * share["native"], 100 * share["numpy"],
                 100 * share["py"], bff))

    print("")
    grand = sum(r[1] * (r[2] + 1) for r in measured)
    print("%d models measured. %.1f ms per LM iteration in total; %.1f ms "
          "(%.0f%%) of that is already inside a C++ extension, and %.1f ms "
          "(%.0f%%) is the compute that could move."
          % (len(measured), grand / 1e3, (grand - movable_total) / 1e3,
             100.0 * (grand - movable_total) / grand, movable_total / 1e3,
             100.0 * movable_total / grand))
    if skipped:
        print("")
        print("not measurable with this fixture (see census_models.py):")
        for name, _, _, _, why in skipped:
            print("  %-*s  %s" % (width, name, why))


if __name__ == "__main__":
    main()
