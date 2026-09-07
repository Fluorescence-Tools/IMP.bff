"""The speed evidence for phase 6: bff.Sampler against the same samplers
in Python.

Not collected by pytest -- the name is bench_, not test_. Run it directly:

    PYTHONPATH=<build>/lib python imp.bff/test/sampler/bench_sampler.py

Three samplers on the same toy problem (the 2D correlated Gaussian of
test_sampler.py), N steps of the stretch move per walker, 16 walkers:

(i)   bff.Sampler with the node-graph objective -- the C++ path: the
      walker is written into the parameter ports, the node graph is
      updated and read, and the accept/reject is decided, all in C++;
(ii)  the same stretch move in Python (chisurf's EnsembleSampler._step
      loop shape) driving the *same* node graph through Python -- every
      port write, node update and output read crosses the wrapper;
(iii) the same Python sampler with a pure-Python objective (a numpy
      chi^2, no ports at all) -- what the boundary costs on top of a
      Python objective.

The number that justified this port: 1.67 us per port set+get across the
wrapper against 0.06 us for a Python attribute. On a two-parameter model
that is the better part of the per-move tax; it grows with the number of
free parameters.
"""
import time

import numpy as np

from IMP.bff import Node, Port, Sampler

MU = np.array([1.0, -0.5])
SIGMA = np.array([[1.0, 0.8], [0.8, 2.0]])
P = np.linalg.inv(SIGMA)

N_STEPS = 10000
N_WALKERS = 16


def const_port(value):
    return Port(float(value))


def linked_port(source):
    p = Port(0.0)
    p.set_link(source)
    return p


#: Every node kept alive: a downstream node holds its upstream ports
#: strongly but their nodes only weakly (chinet ownership -- the session
#: or the caller owns nodes), and a garbage-collected intermediate node
#: makes update() silently read stale values.
_KEEP_ALIVE = []


def operator_node(name, a, b, op):
    node = Node(name)
    node.add_input_port("a", a)
    node.add_input_port("b", b)
    node.add_output_port(name, Port(0.0, False, True))
    node.set_callback(op, "C")
    _KEEP_ALIVE.append(node)
    return node


def gaussian_graph():
    x1, x2 = Port(0.0, name="x1"), Port(0.0, name="x2")
    d1 = operator_node("d1", x1, const_port(-MU[0]), "addition_double")
    d2 = operator_node("d2", x2, const_port(-MU[1]), "addition_double")
    s11 = operator_node("s11", linked_port(d1.get_output_port("d1")),
                        linked_port(d1.get_output_port("d1")), "multiply_double")
    s12 = operator_node("s12", linked_port(d1.get_output_port("d1")),
                        linked_port(d2.get_output_port("d2")), "multiply_double")
    s22 = operator_node("s22", linked_port(d2.get_output_port("d2")),
                        linked_port(d2.get_output_port("d2")), "multiply_double")
    t11 = operator_node("t11", linked_port(s11.get_output_port("s11")),
                        const_port(P[0, 0]), "multiply_double")
    t12 = operator_node("t12", linked_port(s12.get_output_port("s12")),
                        const_port(2.0 * P[0, 1]), "multiply_double")
    t22 = operator_node("t22", linked_port(s22.get_output_port("s22")),
                        const_port(P[1, 1]), "multiply_double")
    u = operator_node("u", linked_port(t11.get_output_port("t11")),
                      linked_port(t12.get_output_port("t12")), "addition_double")
    chi2 = operator_node("chi2", linked_port(u.get_output_port("u")),
                         linked_port(t22.get_output_port("t22")), "addition_double")
    return [x1, x2], chi2


# ---------------------------------------------------------------- (i) C++

def bench_cpp_sampler():
    params, objective = gaussian_graph()
    params[0].set_value(MU[0])
    params[1].set_value(MU[1])
    sampler = Sampler("stretch", 42)
    sampler.set_parameter_ports(params)
    sampler.set_objective(objective, "chi2")
    sampler.set_number_of_walkers(N_WALKERS)
    sampler.run(200)  # warm-up, out of the timing
    t0 = time.perf_counter()
    sampler.run(N_STEPS)
    dt = time.perf_counter() - t0
    return N_STEPS * N_WALKERS / dt, N_STEPS / dt


# ------------------------------------------- (ii)/(iii) the Python samplers

def python_stretch(walkers, log_prob_fn, n_steps, a=2.0, rng=None):
    """chisurf's EnsembleSampler._step loop shape, in Python."""
    if rng is None:
        rng = np.random.default_rng(42)
    nwalkers, ndim = walkers.shape
    lp = np.array([log_prob_fn(w) for w in walkers])
    for _ in range(n_steps):
        order = rng.permutation(nwalkers)
        half = nwalkers // 2
        for active, complement in ((order[:half], order[half:]),
                                   (order[half:], order[:half])):
            s = walkers[active]
            z = ((a - 1.0) * rng.random(len(active)) + 1.0) ** 2 / a
            partners = walkers[complement][
                rng.integers(len(complement), size=len(active))]
            q = partners - (partners - s) * z[:, None]
            factors = (ndim - 1.0) * np.log(z)
            new_lp = np.array([log_prob_fn(x) for x in q])
            take = (factors + new_lp - lp[active]) > np.log(
                rng.random(len(active)))
            idx = active[take]
            walkers[idx] = q[take]
            lp[idx] = new_lp[take]
    return walkers


def bench_python_graph_sampler():
    params, objective = gaussian_graph()
    params[0].set_value(MU[0])
    params[1].set_value(MU[1])
    p0, p1 = params
    output = objective.get_output_port("chi2")

    def log_prob(x):
        # the same graph, driven from Python: three boundary crossings per
        # evaluation (two writes, one read) plus the node update
        p0.set_value(x[0])
        p1.set_value(x[1])
        objective.update()
        chi2 = output.get_value()
        return -0.5 * chi2

    walkers = MU + 1e-2 * np.random.default_rng(1).standard_normal(
        (N_WALKERS, 2))
    python_stretch(walkers.copy(), log_prob, 20)  # warm-up
    t0 = time.perf_counter()
    python_stretch(walkers, log_prob, N_STEPS)
    dt = time.perf_counter() - t0
    return N_STEPS * N_WALKERS / dt, N_STEPS / dt


def bench_python_objective_sampler():
    def log_prob(x):
        d = x - MU
        return -0.5 * float(d @ P @ d)

    walkers = MU + 1e-2 * np.random.default_rng(1).standard_normal(
        (N_WALKERS, 2))
    python_stretch(walkers.copy(), log_prob, 20)  # warm-up
    t0 = time.perf_counter()
    python_stretch(walkers, log_prob, N_STEPS)
    dt = time.perf_counter() - t0
    return N_STEPS * N_WALKERS / dt, N_STEPS / dt


def main():
    print(f"toy: 2D correlated Gaussian, stretch move, {N_STEPS} steps, "
          f"{N_WALKERS} walkers")
    cpp_evals, cpp_steps = bench_cpp_sampler()
    pyg_evals, pyg_steps = bench_python_graph_sampler()
    pyo_evals, pyo_steps = bench_python_objective_sampler()
    print()
    print(f"(i)   bff.Sampler, node-graph objective   : "
          f"{cpp_steps:10.0f} steps/s  ({cpp_evals:10.0f} walker-evals/s)")
    print(f"(ii)  python sampler, same graph via bff  : "
          f"{pyg_steps:10.0f} steps/s  ({pyg_evals:10.0f} walker-evals/s)")
    print(f"(iii) python sampler, python objective    : "
          f"{pyo_steps:10.0f} steps/s  ({pyo_evals:10.0f} walker-evals/s)")
    print()
    print(f"(i) / (ii): {cpp_evals / pyg_evals:6.1f}x  "
          f"(the SWIG boundary per move, removed)")
    print(f"(i) / (iii): {cpp_evals / pyo_evals:6.1f}x "
          f"(boundary removed, objective also in C++)")


if __name__ == "__main__":
    main()
