"""A whole fit as one C++ graph, timed against ChiSurf's Python evaluation.

The model is ChiSurf's own collinear test case, ``c + a*x + b*x**2`` over 96
points. On the bff side the polynomial and its chi-square are nodes, so a
sampler move touches no interpreter at all; on the ChiSurf side the same
move goes through ParseModel's ``eval`` and the parameter wrappers.
"""

import time

import numpy as np

from IMP import bff

SIGMA = 0.1
N_WALKERS = 16
N_STEPS = 4000

_KEEP_ALIVE = []


def _op(name, a, b, op):
    """One operator node, with its result on a port named after the node."""
    node = bff.GraphNode(name)
    node.add_input_port("a", a)
    node.add_input_port("b", b)
    node.add_output_port(name, bff.GraphPort(0.0, False, True))
    node.set_callback(op, "C")
    _KEEP_ALIVE.append(node)
    return node


def _follow(port):
    p = bff.GraphPort(0.0)
    p.set_link(port)
    return p


def build_graph(x, y, ey):
    """c + a*x + b*x**2, then its chi-square, as bff nodes."""
    a, b, c = bff.GraphPort(2.0, name="a"), bff.GraphPort(0.5, name="b"), bff.GraphPort(1.0, name="c")
    xp = bff.GraphPort(list(x), name="x")

    ax = _op("ax", xp, a, "multiply_double")
    x2 = _op("x2", xp, bff.GraphPort(list(x)), "multiply_double")
    bx2 = _op("bx2", _follow(x2.get_output_port("x2")), b, "multiply_double")
    s1 = _op("s1", _follow(ax.get_output_port("ax")),
             _follow(bx2.get_output_port("bx2")), "addition_double")
    model = _op("model", _follow(s1.get_output_port("s1")), c, "addition_double")

    chi2 = bff.FitChiSquared("chi2")
    chi2.add_input_port("model", _follow(model.get_output_port("model")))
    chi2.add_output_port("chi2", bff.GraphPort(0.0, name="chi2"))
    chi2.set_data(list(y), list(ey))
    _KEEP_ALIVE.append(chi2)
    return [a, b, c], chi2


def main():
    rng = np.random.default_rng(0)
    x = np.linspace(1.0, 2.0, 96)
    y = 1.0 + 2.0 * x + 0.5 * x ** 2 + rng.normal(0.0, SIGMA, x.size)
    ey = np.ones_like(y) * SIGMA

    params, objective = build_graph(x, y, ey)
    sampler = bff.Sampler("stretch", 42)
    sampler.set_parameter_ports(params)
    sampler.set_objective(objective, "chi2")
    sampler.set_number_of_walkers(N_WALKERS)
    sampler.run(200)
    t0 = time.perf_counter()
    sampler.run(N_STEPS)
    dt = time.perf_counter() - t0
    evals = N_STEPS * N_WALKERS
    print(f"bff: whole fit in C++            : {evals / dt:10.0f} evals/s "
          f"({dt / evals * 1e6:6.2f} us/eval)")

    chain = np.asarray(sampler.chain).reshape(-1, len(params))
    mean = chain.mean(axis=0)
    print(f"     posterior mean (a, b, c)    : "
          f"{mean[0]:.3f} {mean[1]:.3f} {mean[2]:.3f}   (truth 2.0 0.5 1.0)")

    try:
        import sys
        sys.path.insert(0, "/Users/tpeulen/dev/chisurf/test/fitting")
        from test_blocked_sampler import _collinear_fit
        from chisurf.core.fitting.fit import get_chi2
        fit = _collinear_fit()
        model_obj = fit.model
        n = 2000
        t0 = time.perf_counter()
        for i in range(n):
            model_obj.parameters[0].value = 1.0 + 1e-6 * i
            model_obj.update_model()
            get_chi2([], model=model_obj, reduced=False)
        dt = time.perf_counter() - t0
        print(f"chisurf: same fit in Python      : {n / dt:10.0f} evals/s "
              f"({dt / n * 1e6:6.2f} us/eval)")
    except Exception as exc:  # pragma: no cover - chisurf is optional here
        print(f"chisurf comparison unavailable: {exc}")


if __name__ == "__main__":
    main()
