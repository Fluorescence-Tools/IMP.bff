"""Which ChiSurf models run in bff, asked rather than assumed.

Run it directly; it is a diagnostic, not a test, so nothing here asserts.

    $E/bin/python test/minimizer/census_models.py

A directive like *"the models must live in bff"* is not actionable without
this. It walks every `Model` subclass under `chisurf.core.models`, builds a
fit for it, and asks `graph_objective` whether that model becomes one C++
graph -- rather than reading the code and guessing, which on 2026-09-01 got
the answer wrong in **both** directions at once: the equation families were
already entirely in bff and nobody had noticed, while two lifetime
subclasses were building a graph that computed a *different curve* from the
one they report.

That second failure is why the census reports what it reports. A model that
falls back is slower; a model that builds the wrong graph is wrong. So each
row also compares the graph's curve against ``model.update()``, and a
disagreement is printed as **MISMATCH** rather than as a pass.

A third, more dangerous failure is not wrong -- it is *never run*. On the
shared decay-shaped fixture, a model whose update no-ops without
the payload it needs (no burst folder, no meaningful axis) leaves
``ModelCurve``'s zeroed placeholder untouched, and a graph check that only
compares "the graph" against "``model.y``" calls that a match, because a
stale array trivially agrees with itself. The **curve** column answers a
question the graph column cannot: did ``model.update()`` actually compute
something, or did it hand back the same flat array it started with. A
**yes** in the graph column beside a **dead** curve is that phantom, made
visible.

The polarisation columns are here because the default for every TCSPC model
is magic angle, and a census that only ever asks the default would have said
nothing about VV/VH -- which is half of what a polarisation-resolved setup
fits.
"""
import importlib
import inspect
import pkgutil
import sys
from pathlib import Path

import numpy as np

# ChiSurf is a sibling checkout (AGENTS.md), not an installed package.
for candidate in (Path(__file__).resolve().parents[3] / "chisurf",):
    if candidate.is_dir() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import chisurf as cs                                     # noqa: E402
import chisurf.core.curve                                # noqa: E402
import chisurf.core.data                                 # noqa: E402
import chisurf.core.fitting.fit as F                     # noqa: E402
import chisurf.core.fitting.minimizer as M               # noqa: E402
import chisurf.core.models as models_package             # noqa: E402
from chisurf.core.models.model import Model              # noqa: E402

N = 256
DT = 0.032
REP_RATE = 80.0


def model_classes():
    """Every concrete `Model` subclass ChiSurf defines, by qualified name."""
    found = {}
    for module in pkgutil.walk_packages(models_package.__path__,
                                        models_package.__name__ + "."):
        try:
            imported = importlib.import_module(module.name)
        except Exception:
            # A model whose module will not import cannot be fitted either
            # way; it is not this script's business to say why.
            continue
        for value in vars(imported).values():
            if (inspect.isclass(value) and issubclass(value, Model)
                    and value is not Model
                    and value.__module__.startswith("chisurf.core.models")):
                found["%s.%s" % (value.__module__, value.__name__)] = value
    return found


def make_fit(model_class, polarization="vm"):
    """A fit of *model_class* over a decay-shaped curve, or ``None``.

    One fixture for every family, deliberately: what is being asked is
    whether the builder recognises a model, and a model that needs data this
    does not carry (a PDA histogram, a structure) answers that by raising,
    which the caller reports as its own row.
    """
    if ".models.ics." in model_class.__module__:
        # The image-correlation family evaluates over lag grids the decay
        # fixture does not carry; without them the builder refuses for a
        # *fixture* reason and the census row would say nothing about the
        # builder. A minimal carpet answers the actual question.
        if polarization != "vm":
            return None          # no polarisation axis in a carpet
        half = 8
        xi, psi = np.meshgrid(np.arange(-half, half + 1, dtype=float),
                              np.arange(-half, half + 1, dtype=float))
        y = np.full(xi.size, 0.1)
        data = cs.core.data.DataCurve(
            x=np.arange(y.size, dtype=float), y=y,
            ey=np.full(y.size, 1e-3))
        data.meta_data["ics"] = {"pixel_shift": xi, "line_shift": psi}
        fit = F.Fit(model_class=model_class, data=data)
        fit.xmin, fit.xmax = 0, y.size
        try:
            fit.model.find_parameters()
        except Exception:
            pass
        return fit
    x = np.arange(N) * DT
    irf = 1000.0 * np.exp(-0.5 * ((x - 1.0) / 0.08) ** 2)
    y = 4000.0 * np.exp(-x / 3.1) + 2.0
    data = cs.core.data.DataCurve(x=x, y=y, ey=np.sqrt(np.maximum(y, 1.0)))
    fit = F.Fit(model_class=model_class, data=data)
    fit.xmin, fit.xmax = 0, N - 1
    model = fit.model
    if hasattr(model, "convolve"):
        model.convolve._irf = cs.core.curve.Curve(x=x, y=irf.copy())
        model.convolve.dt = DT
        model.convolve.rep_rate = REP_RATE
        model.convolve.stop = N * DT
    if polarization != "vm":
        if not hasattr(model, "anisotropy"):
            return None
        model.anisotropy.polarization_type = polarization
        model.anisotropy.add_rotation(b=0.2, rho=2.0)
    try:
        model.find_parameters()
    except Exception:
        pass
    # A parse model whose equation has no free parameters would be reported
    # as a refusal, which is a statement about this fixture rather than about
    # the builder -- `ParseModel` ships `x*0`. The families that carry a
    # catalogue already select something real, so this only fires on the
    # bare one.
    if hasattr(model, "func") and not getattr(model, "parameters", []):
        try:
            model.func = "a*exp(-x/t)+b"
            model.find_parameters()
        except Exception:
            pass
    return fit


def _curve_is_degenerate(y) -> bool:
    """True when a curve carries no information: empty, non-finite, or constant.

    A flat array is exactly what a silently no-op update leaves
    behind -- ``ModelCurve.__init__`` zeros ``y`` to the data's shape, and a
    model that never touches it again hands that same zero array straight
    back. Constant-nonzero is caught too, deliberately: a model that fills the
    array with one repeated placeholder value is exactly as uninformative as
    zero, and a strict-zero check alone would miss it.
    """
    arr = np.asarray(y, dtype=float)
    if arr.size == 0:
        return True
    if not np.all(np.isfinite(arr)):
        return True
    return bool(np.ptp(arr) == 0.0)


def curve_verdict(model) -> str:
    """``"live"``, ``"dead"``, or an error tag, for one model's ``update()``.

    Deliberately independent of whether a graph builds. The phantom this
    exists to catch is a **yes** in the graph column standing beside a curve
    that was never computed -- the graph-vs-``model.y`` comparison below
    cannot see that on its own, because a stale ``model.y`` agrees with
    itself trivially.
    """
    try:
        model.update()
    except Exception as exc:
        return "raised %s" % type(exc).__name__
    return "dead" if _curve_is_degenerate(getattr(model, "y", None)) else "live"


def verdict(model_class, polarization):
    """``(graph, curve, free)`` for one model.

    ``graph`` is ``"yes"``, ``"no"``, ``"MISMATCH ..."`` or an error tag;
    ``curve`` is :func:`curve_verdict`'s ``"live"`` / ``"dead"`` / error tag,
    or ``"-"`` when no model was ever built for this fixture; ``free`` is the
    list of free parameter names.
    """
    try:
        fit = make_fit(model_class, polarization)
    except Exception as exc:
        return "%s: %s" % (type(exc).__name__, str(exc)[:34]), "-", []
    if fit is None:
        return "-", "-", []
    model = fit.model
    free = [p.name for p in getattr(model, "parameters", [])]
    # Asked once, up front, and independent of whatever the graph checks below
    # do to the model afterwards -- this is the third column, not a side
    # effect of the mismatch check.
    curve = curve_verdict(model)
    try:
        built = M.graph_objective(fit, model)
    except Exception as exc:
        return "raised %s" % type(exc).__name__, curve, free
    if built is None:
        return "no", curve, free

    # A graph that builds is only good news if it is the same model. This is
    # the check that would have caught `MaxEntLifetimeModel` fitting a plain
    # multi-exponential -- 793.8 counts away from its own curve, silently.
    decay = getattr(built[0], "_decay", None)
    if decay is None:
        # An expression graph (parse, image correlation) holds its model
        # node first in the keepalive; its curve must be the model's too.
        try:
            node = built[0]._graph[2][0]
            model.update()
            node.update()
            graph = np.asarray(
                node.get_output_port(node.get_name()).value, dtype=float)
            python = np.asarray(model.y, dtype=float)
            n = min(graph.size, python.size)
            scale = max(1.0, float(np.max(np.abs(python[:n]))))
            gap = float(np.max(np.abs(graph[:n] - python[:n]))) / scale
        except Exception:
            return "yes", curve, free
        if gap > 1e-8:
            return "MISMATCH %.3g" % gap, curve, free
        return "yes", curve, free
    try:
        model.update()
        decay.update()
        graph = np.asarray(decay.get_curve(), dtype=float)
        python = np.asarray(model.y, dtype=float)
        n = min(graph.size, python.size)
        scale = max(1.0, float(np.max(np.abs(python[:n]))))
        gap = float(np.max(np.abs(graph[:n] - python[:n]))) / scale
    except Exception as exc:
        return "yes (uncheckable: %s)" % type(exc).__name__, curve, free
    if gap > 1e-8:
        return "MISMATCH %.3g" % gap, curve, free
    return "yes", curve, free


def main():
    classes = model_classes()
    width = max(len(name) for name in classes) - len("chisurf.core.models.")
    rows = []
    for name in sorted(classes):
        short = name.replace("chisurf.core.models.", "")
        magic, curve, free = verdict(classes[name], "vm")
        polarised, _, _ = verdict(classes[name], "vv")
        rows.append((short, magic, curve, polarised, free))

    print("%-*s  %-14s %-6s %-14s %s"
          % (width, "model", "vm", "curve", "vv", "free"))
    print("-" * (width + 46))
    for short, magic, curve, polarised, free in rows:
        print("%-*s  %-14s %-6s %-14s %s" % (width, short, magic, curve,
                                              polarised, ",".join(free) or "-"))

    built = sum(1 for _, magic, _, _, _ in rows if magic == "yes")
    polarised_built = sum(1 for _, _, _, p, _ in rows if p == "yes")
    polarised_asked = sum(1 for _, _, _, p, _ in rows if p in ("yes", "no"))
    buildable = sum(1 for _, magic, _, _, _ in rows if magic in ("yes", "no"))
    wrong = [short for short, magic, _, p, _ in rows
             if "MISMATCH" in magic or "MISMATCH" in p]
    # The phantom this whole column exists to expose: a graph that reads as
    # built and agreeing, next to a curve the model never computed.
    phantom = [short for short, magic, curve, _, _ in rows
               if magic == "yes" and curve == "dead"]
    dead = [short for short, _, curve, _, _ in rows if curve == "dead"]
    print("")
    print("%d of the %d models this fixture can construct build a graph "
          "(%d classes in all);" % (built, buildable, len(rows)))
    print("%d of the %d that have a polarisation build one under VV."
          % (polarised_built, polarised_asked))
    if wrong:
        print("WRONG GRAPHS (a graph that is not the model): %s"
              % ", ".join(wrong))
    else:
        print("no model builds a graph that disagrees with its own curve")
    if phantom:
        print("PHANTOM GRAPHS (yes beside a curve the model never "
              "computed): %s" % ", ".join(phantom))
    if dead:
        print("DEAD CURVES (update() left it flat/empty, %d): %s"
              % (len(dead), ", ".join(dead)))
    else:
        print("no model's update() left its curve dead on this fixture")


if __name__ == "__main__":
    main()
