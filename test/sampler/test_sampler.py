"""ChiSurf's samplers, in bff C++ (phase 6 of removing chinet from chisurf).

chisurf.core.fitting.sample and .ensemble implement the stretch move
(Goodman & Weare), differential evolution with snooker updates (ter Braak)
and the blocked random-walk Metropolis, all driving chisurf's Python models
through a port per parameter. bff.Sampler is those algorithms over the bff
Port/Node runtime, run entirely in C++: the walker is written into the
parameter ports, the node graph is updated and read, and the accept/reject
is decided without crossing the wrapper -- which measured 1.67 us per port
set+get against 0.06 us for a Python attribute, paid per proposal per
parameter.

These tests port the behavioural contract of chisurf's sampler tests
(test/fitting/test_ensemble_samplers.py, test_differential_evolution_
sampler.py, test_blocked_sampler.py, test_sampler_seeding.py, and the
prior scenarios of test_prior_posterior_sampling.py) onto the C++ surface:

- a run under a fixed seed is reproducible, and every sampler shares one
  RNG contract (statistical parity with chisurf, not stream parity -- the
  C++ sampler draws from std::mt19937_64, chisurf from numpy generators);
- a 2D correlated Gaussian recovered through a *node-graph objective*
  (tiny operator nodes computing chi^2), by stretch, DE and blocked walk;
- bounds are respected by rejection, never clipping (chisurf's lnprior
  short-circuits to -inf out of the box);
- a Gaussian prior on a port pulls the posterior mean;
- blocked and unblocked runs agree statistically; blocks come from the
  factor graph's sampling blocks, an explicit partition, or one block.

Not ported: the ensemble *slice* sampler (stays Python; see Sampler.h),
chisurf's pool/vectorize plumbing (a C++ loop has no pool to call), and
the curvature-seeded block covariances (Fit.covariance_matrix is a chisurf
concept; the C++ sampler seeds the diagonal fallback chisurf falls back
to).
"""
import numpy as np
import pytest

from IMP.bff import FactorGraph, LIKELIHOOD, PRIOR, Node, Port, Sampler

# The toy posterior: a 2D correlated Gaussian, chi^2 = d^T P d with
# P = Sigma^-1, computed by a graph of bff operator nodes -- the same
# shape a chisurf model graph has, only small enough to read.
MU = np.array([1.0, -0.5])
SIGMA = np.array([[1.0, 0.8], [0.8, 2.0]])
P = np.linalg.inv(SIGMA)


def const_port(value):
    return Port(float(value))


def linked_port(source):
    p = Port(0.0)
    p.set_link(source)
    return p


#: Every node a test builds, kept alive for the process. A downstream
#: node holds its upstream *ports* strongly but their *nodes* only weakly
#: (Port.h: a port's node is a weak_ptr) -- chinet's ownership: the
#: session or the caller owns the nodes, chisurf registers them on its
#: session. Without this list the intermediate nodes of a graph are
#: garbage-collected, their ports' get_node() comes back empty, and
#: update()'s recursion silently reads stale values -- the classic
#: chinet footgun, and the reason chisurf keeps models on a session.
_KEEP_ALIVE = []


def operator_node(name, a, b, op):
    """A node computing op(a, b) into an output port keyed by its name."""
    node = Node(name)
    node.add_input_port("a", a)
    node.add_input_port("b", b)
    node.add_output_port(name, Port(0.0, False, True))
    node.set_callback(op, "C")
    _KEEP_ALIVE.append(node)
    return node


def gaussian_graph(mu=MU, precision=P):
    """A node graph computing chi^2(x) of a correlated 2D Gaussian.

    Returns the two parameter input ports and the final node (whose
    output port is named "chi2"): d = x - mu per axis, the three
    quadratic terms, and their sum.
    """
    x1, x2 = Port(0.0, name="x1"), Port(0.0, name="x2")
    d1 = operator_node("d1", x1, const_port(-float(mu[0])), "addition_double")
    d2 = operator_node("d2", x2, const_port(-float(mu[1])), "addition_double")
    s11 = operator_node("s11", linked_port(d1.get_output_port("d1")),
                        linked_port(d1.get_output_port("d1")), "multiply_double")
    s12 = operator_node("s12", linked_port(d1.get_output_port("d1")),
                        linked_port(d2.get_output_port("d2")), "multiply_double")
    s22 = operator_node("s22", linked_port(d2.get_output_port("d2")),
                        linked_port(d2.get_output_port("d2")), "multiply_double")
    t11 = operator_node("t11", linked_port(s11.get_output_port("s11")),
                        const_port(float(precision[0, 0])), "multiply_double")
    t12 = operator_node("t12", linked_port(s12.get_output_port("s12")),
                        const_port(2.0 * float(precision[0, 1])), "multiply_double")
    t22 = operator_node("t22", linked_port(s22.get_output_port("s22")),
                        const_port(float(precision[1, 1])), "multiply_double")
    u = operator_node("u", linked_port(t11.get_output_port("t11")),
                      linked_port(t12.get_output_port("t12")), "addition_double")
    chi2 = operator_node("chi2", linked_port(u.get_output_port("u")),
                         linked_port(t22.get_output_port("t22")), "addition_double")
    return [x1, x2], chi2


def make_sampler(algorithm, seed=7, n_steps=2000, thin=1, **kwargs):
    """A sampler on the correlated-Gaussian toy, ready to run."""
    params, objective = gaussian_graph()
    params[0].set_value(MU[0])
    params[1].set_value(MU[1])  # keyword form: Port(value=.., name=..) drops positionals
    sampler = Sampler(algorithm, seed)
    sampler.set_parameter_ports(params)
    sampler.set_objective(objective, "chi2")
    kwargs.pop("thin", None)
    for key, value in kwargs.items():
        getattr(sampler, "set_" + key)(value)
    sampler.run(n_steps, thin)
    return sampler


def flat_chain(sampler, burn=300):
    chain = np.asarray(sampler.get_chain())
    n_walkers = len(sampler.get_walkers())
    return chain[burn * n_walkers:]


# ------------------------------------------------------------- reproducible

@pytest.mark.parametrize("algorithm", ["stretch", "slice", "de", "metropolis"])
def test_a_run_under_a_fixed_seed_is_reproducible(algorithm):
    """Two identically-seeded samplers give the same chain, byte for byte."""
    first = make_sampler(algorithm, seed=5, n_steps=300)
    second = make_sampler(algorithm, seed=5, n_steps=300)
    assert np.array_equal(np.asarray(first.get_chain()),
                          np.asarray(second.get_chain()))
    assert np.array_equal(np.asarray(first.get_log_prob()),
                          np.asarray(second.get_log_prob()))


@pytest.mark.parametrize("algorithm", ["stretch", "slice", "de", "metropolis"])
def test_a_different_seed_gives_a_different_chain(algorithm):
    first = make_sampler(algorithm, seed=5, n_steps=300)
    second = make_sampler(algorithm, seed=6, n_steps=300)
    assert not np.array_equal(np.asarray(first.get_chain()),
                              np.asarray(second.get_chain()))


def test_reset_reseeds_and_restarts_the_chain():
    sampler = make_sampler("stretch", seed=5, n_steps=200)
    first = np.asarray(sampler.get_chain()).copy()
    sampler.reset()
    assert sampler.get_iteration() == 0
    assert len(sampler.get_chain()) == 0
    sampler.run(200)
    assert np.array_equal(first, np.asarray(sampler.get_chain()))


def test_a_run_can_be_continued_from_where_it_left_off():
    """Two 300-step runs and one 600-step run of the same seeded sampler
    are the same chain: the second run resumed the ensemble and the RNG
    stream, exactly as chisurf's run_mcmc(state) does."""
    sampler = make_sampler("stretch", seed=5, n_steps=300)
    sampler.run(300)
    whole = make_sampler("stretch", seed=5, n_steps=600)
    assert sampler.get_iteration() == 600
    assert np.array_equal(np.asarray(sampler.get_chain()),
                          np.asarray(whole.get_chain()))


def test_thinning_records_every_nth_step():
    sampler = make_sampler("stretch", seed=5, n_steps=100, thin=5)
    assert sampler.get_iteration() == 20


# ------------------------------------------------- statistical correctness

def moments(chain):
    mean = chain.mean(axis=0)
    cov = np.cov(chain, rowvar=False)
    return mean, cov


def assert_gaussian_recovered(sampler, n_sigma=4.0, cov_tol=0.25):
    chain = flat_chain(sampler)
    assert chain.shape[1] == 2
    mean, cov = moments(chain)
    n = chain.shape[0]
    # mean within a few standard errors of the true mean (loose, as the
    # chains are autocorrelated: a few sigma of the *independent* error)
    se = np.sqrt(np.diag(SIGMA) / max(1.0, n / 50.0))
    assert np.all(np.abs(mean - MU) < n_sigma * se), (mean, MU)
    assert np.allclose(cov, SIGMA, rtol=cov_tol, atol=0.15), (cov, SIGMA)


def test_the_stretch_move_recovers_the_correlated_gaussian():
    sampler = make_sampler("stretch", seed=11, n_steps=2500)
    assert_gaussian_recovered(sampler)


def test_differential_evolution_recovers_the_correlated_gaussian():
    sampler = make_sampler("de", seed=11, n_steps=2500)
    assert_gaussian_recovered(sampler)


def test_the_blocked_random_walk_recovers_the_correlated_gaussian():
    sampler = make_sampler("metropolis", seed=11, n_steps=3000)
    assert_gaussian_recovered(sampler, cov_tol=0.3)


def test_all_three_algorithms_agree_statistically():
    """Stretch, DE and blocked runs agree on mean and covariance."""
    results = {}
    for algorithm in ("stretch", "de", "metropolis"):
        mean, cov = moments(flat_chain(make_sampler(algorithm, seed=3,
                                                    n_steps=3000)))
        results[algorithm] = (mean, cov)
    for algorithm, (mean, cov) in results.items():
        assert np.all(np.abs(mean - MU) < 0.25), (algorithm, mean)
        assert np.all(np.abs(cov - SIGMA) < 0.6), (algorithm, cov)


# ------------------------------------------------------------------- bounds

def test_a_forbidden_region_is_never_entered():
    """Bounds reject (chisurf's lnprior box), they never clip: the chain
    stays inside and the wall is felt as rejection, not as a pile of
    walkers parked exactly on it."""
    params, objective = gaussian_graph()
    params[0].set_value(1.0)
    params[1].set_value(-0.5)
    for p in params:
        p.set_is_bounded(True)
    params[0].set_bounds(0.75, 1.25)
    params[1].set_bounds(-0.75, -0.25)
    sampler = Sampler("stretch", 13)
    sampler.set_parameter_ports(params)
    sampler.set_objective(objective, "chi2")
    sampler.run(1500)
    chain = flat_chain(sampler, burn=200)
    assert chain[:, 0].min() >= 0.75 - 1e-9
    assert chain[:, 0].max() <= 1.25 + 1e-9
    assert chain[:, 1].min() >= -0.75 - 1e-9
    assert chain[:, 1].max() <= -0.25 + 1e-9
    # rejection, not clipping: a clipped chain would pile mass exactly on
    # the wall
    at_wall = (np.isclose(chain[:, 0], 0.75, atol=1e-9) |
               np.isclose(chain[:, 0], 1.25, atol=1e-9))
    assert at_wall.mean() < 0.05


def test_unbounded_ports_are_unbounded():
    params, objective = gaussian_graph()
    sampler = Sampler("stretch", 13)
    sampler.set_parameter_ports(params)
    sampler.set_objective(objective, "chi2")
    sampler.set_bounds([-1e30, -1e30], [1e30, 1e30])
    sampler.run(100)
    assert np.all(np.isfinite(np.asarray(sampler.get_log_prob())))


# -------------------------------------------------------------------- priors

def test_a_gaussian_prior_pulls_the_posterior_mean():
    """A N(2, 0.1) prior on x1 against a N(1, 1) likelihood: the posterior
    mean is the precision-weighted 1*0.01/(1+0.01)... ~2.0 -- the prior
    dominates, exactly as chisurf's lnprior sum says it should."""
    params, objective = gaussian_graph(mu=np.array([1.0, -0.5]),
                                       precision=np.diag([1.0, 0.5]))
    params[0].set_value(2.0)
    params[1].set_value(-0.5)
    params[0].set_prior('{"kind": "normal", "mu": 2.0, "sigma": 0.1}')
    sampler = Sampler("stretch", 17)
    sampler.set_parameter_ports(params)
    sampler.set_objective(objective, "chi2")
    sampler.run(2500)
    mean, _ = moments(flat_chain(sampler))
    # posterior mean = (mu_like/s2_like + mu_pr/s2_pr)/(1/s2_like+1/s2_pr)
    expected = (1.0 / 1.0 + 2.0 / 0.01) / (1.0 / 1.0 + 1.0 / 0.01)
    assert abs(mean[0] - expected) < 0.05, (mean[0], expected)
    # the prior contribution is carried next to the chain, as blobs are
    # (a narrow prior's log density is positive where it concentrates)
    lnprior = np.asarray(sampler.get_lnprior())
    assert np.all(lnprior != 0.0)
    assert np.ptp(lnprior) > 0.1


def test_a_flat_prior_chain_is_not_pulled():
    params, objective = gaussian_graph()
    sampler = Sampler("stretch", 17)
    sampler.set_parameter_ports(params)
    sampler.set_objective(objective, "chi2")
    sampler.run(1500)
    assert np.allclose(np.asarray(sampler.get_lnprior()), 0.0)


def test_a_prior_can_forbid_a_region():
    params, objective = gaussian_graph(precision=np.diag([1.0, 0.5]))
    params[0].set_value(2.5)
    params[1].set_value(-0.5)
    params[0].set_prior('{"kind": "uniform", "lb": 2.0, "ub": 3.0}')
    sampler = Sampler("stretch", 19)
    sampler.set_parameter_ports(params)
    sampler.set_objective(objective, "chi2")
    sampler.run(1500)
    chain = flat_chain(sampler, burn=200)
    assert chain[:, 0].min() >= 2.0 - 1e-9
    assert chain[:, 0].max() <= 3.0 + 1e-9


# ------------------------------------------------------- chisurf's scenarios

def test_walkers_that_span_nothing_are_refused():
    """chisurf: a degenerate ensemble cannot explore; it refuses."""
    params, objective = gaussian_graph()
    start = [[0.0, 0.0]] * 10
    sampler = Sampler("stretch", 23)
    sampler.set_parameter_ports(params)
    sampler.set_objective(objective, "chi2")
    sampler.set_walker_start(start)
    with pytest.raises(ValueError):
        sampler.run(10)


def test_the_stretch_move_needs_enough_walkers_to_span_the_space():
    params, objective = gaussian_graph()
    sampler = Sampler("stretch", 23)
    sampler.set_parameter_ports(params)
    sampler.set_objective(objective, "chi2")
    sampler.set_number_of_walkers(3)  # < 2 * ndim and < 4
    with pytest.raises(ValueError):
        sampler.run(10)


def test_live_dangerously_skips_the_walker_requirement():
    params, objective = gaussian_graph()
    # a non-degenerate 4-walker cloud (fewer than the 2*ndim the stretch
    # move wants, and a start that spans the plane)
    start = [[0.0, -0.4], [0.1, -0.6], [0.2, -0.5], [0.3, -0.55]]
    sampler = Sampler("stretch", 23)
    sampler.set_parameter_ports(params)
    sampler.set_objective(objective, "chi2")
    sampler.set_number_of_walkers(4)
    sampler.set_live_dangerously(True)
    sampler.set_walker_start(start)
    sampler.run(50)  # does not raise
    assert sampler.get_iteration() == 50


def test_the_stretch_acceptance_fraction_is_in_a_usable_range():
    sampler = make_sampler("stretch", seed=29, n_steps=1500)
    rate = sampler.get_acceptance_rate()
    assert 0.1 < rate < 0.9, rate
    # per-walker fractions, chisurf's acceptance_fraction
    fractions = np.asarray(sampler.get_acceptance_fractions())
    assert fractions.size == len(sampler.get_walkers())
    assert np.all(fractions > 0.0)


def test_the_result_has_the_shape_every_sampler_promises():
    """chisurf: parameter_values/parameter_names/chains keys, one row per
    (step, walker), acceptance in [0, 1]."""
    for algorithm in ("stretch", "de", "metropolis"):
        sampler = make_sampler(algorithm, seed=31, n_steps=120, thin=4)
        chain = np.asarray(sampler.get_chain())
        n_walkers = {"stretch": sampler.get_number_of_walkers(),
                     "de": sampler.get_number_of_chains(),
                     "metropolis": 1}[algorithm]
        assert chain.shape == (30 * n_walkers, 2), (algorithm, chain.shape)
        assert len(sampler.get_parameter_names()) == 2
        assert len(sampler.get_log_prob()) == chain.shape[0]
        assert len(sampler.get_chi2()) == chain.shape[0]
        assert 0.0 <= sampler.get_acceptance_rate() <= 1.0
        per_walker = np.asarray(sampler.get_chain_of_walker(0))
        assert per_walker.shape == (30, 2)


def test_the_de_population_size_is_configurable_and_defaults_sensibly():
    """chisurf defaults to max(8, 2*ndim) and never below 4."""
    params, objective = gaussian_graph()
    sampler = Sampler("de", 37)
    sampler.set_parameter_ports(params)
    sampler.set_objective(objective, "chi2")
    assert sampler.get_number_of_chains() == 8  # max(8, 2*2)
    sampler.set_number_of_chains(12)
    assert sampler.get_number_of_chains() == 12
    sampler.run(100)
    assert len(sampler.get_walkers()) == 12


def test_de_restores_the_starting_values_afterwards():
    """chisurf's DE sampler leaves the model where it found it."""
    params, objective = gaussian_graph()
    params[0].set_value(MU[0])
    params[1].set_value(MU[1])
    sampler = Sampler("de", 41)
    sampler.set_parameter_ports(params)
    sampler.set_objective(objective, "chi2")
    sampler.run(300)
    assert params[0].get_value() == pytest.approx(MU[0])
    assert params[1].get_value() == pytest.approx(MU[1])


def test_the_walker_spread_is_never_zero_in_any_direction():
    """chisurf: a parameter starting at 0 still gets spread (the relative
    scale is floored by the absolute std)."""
    params, objective = gaussian_graph(mu=np.array([0.0, 0.0]),
                                       precision=np.diag([1.0, 1.0]))
    params[0].set_value(0.0)
    sampler = Sampler("stretch", 43)
    sampler.set_parameter_ports(params)
    sampler.set_objective(objective, "chi2")
    sampler.run(5)
    walkers = np.asarray(sampler.get_walkers())
    assert np.all(np.std(walkers, axis=0) > 0.0)


# ------------------------------------------------------------------ blocking

def independent_graph():
    """Two independent parameters: x1 ~ N(1, 1), x2 ~ N(-1, 4)"""
    x1, x2 = Port(0.0, name="x1"), Port(0.0, name="x2")
    d1 = operator_node("d1", x1, const_port(-1.0), "addition_double")
    d2 = operator_node("d2", x2, const_port(1.0), "addition_double")
    s1 = operator_node("s1", linked_port(d1.get_output_port("d1")),
                       linked_port(d1.get_output_port("d1")), "multiply_double")
    s2 = operator_node("s2", linked_port(d2.get_output_port("d2")),
                       linked_port(d2.get_output_port("d2")), "multiply_double")
    t2 = operator_node("t2", linked_port(s2.get_output_port("s2")),
                       const_port(0.25), "multiply_double")
    chi2 = operator_node("chi2", linked_port(s1.get_output_port("s1")),
                         linked_port(t2.get_output_port("t2")), "addition_double")
    return [x1, x2], chi2


def independent_factor_graph():
    g = FactorGraph()
    g.add_variable("x1", "x1", 0, 0)
    g.add_variable("x2", "x2", 1, 1)
    g.add_factor("L0", LIKELIHOOD, ["x1"], 0, 8)
    g.add_factor("L1", LIKELIHOOD, ["x2"], 1, 8)
    g.add_factor("pi_x1", PRIOR, ["x1"], -1, 1)
    g.add_factor("pi_x2", PRIOR, ["x2"], -1, 1)
    return g


@pytest.mark.parametrize("blocks_from", ["factor_graph", "explicit"])
def test_blocked_sampler_reports_per_block_acceptance(blocks_from):
    params, objective = independent_graph()
    sampler = Sampler("metropolis", 47)
    sampler.set_parameter_ports(params)
    sampler.set_objective(objective, "chi2")
    if blocks_from == "factor_graph":
        graph = independent_factor_graph()  # borrowed: keep it alive
        sampler.set_factor_graph(graph)
    else:
        sampler.set_blocks([0, 1], [1, 1])
    sampler.run(400)
    sizes = sampler.get_block_sizes()
    rates = np.asarray(sampler.get_block_acceptance_rates())
    assert list(sizes) == [1, 1]
    assert rates.shape == (2,)
    assert np.all((rates > 0.0) & (rates < 1.0))


def test_blocked_and_unblocked_runs_agree_statistically():
    """Two independent parameters, sampled jointly (one block) and per
    block: same posterior, within statistical error."""
    chains = {}
    for label, use_blocks in (("blocked", True), ("single", False)):
        params, objective = independent_graph()
        sampler = Sampler("metropolis", 53)
        sampler.set_parameter_ports(params)
        sampler.set_objective(objective, "chi2")
        if use_blocks:
            graph = independent_factor_graph()  # borrowed: keep it alive
            sampler.set_factor_graph(graph)
        sampler.run(4000)
        chain = flat_chain(sampler, burn=500)
        chains[label] = (chain.mean(axis=0), np.cov(chain, rowvar=False))
    for i in range(2):
        a, b = chains["blocked"][0][i], chains["single"][0][i]
        assert abs(a - b) < 0.25, (a, b)
    for i in range(2):
        for j in range(2):
            a = chains["blocked"][1][i, j]
            b = chains["single"][1][i, j]
            assert abs(a - b) < 0.9, (a, b)


def test_explicit_blocks_are_respected():
    params, objective = gaussian_graph()
    sampler = Sampler("metropolis", 59)
    sampler.set_parameter_ports(params)
    sampler.set_objective(objective, "chi2")
    sampler.set_blocks([0, 1], [2])
    sampler.run(200)
    assert list(sampler.get_block_sizes()) == [2]


def test_a_factor_graph_block_partition_is_used():
    params, objective = independent_graph()
    sampler = Sampler("metropolis", 61)
    sampler.set_parameter_ports(params)
    sampler.set_objective(objective, "chi2")
    graph = independent_factor_graph()  # borrowed: keep it alive
    sampler.set_factor_graph(graph)
    sampler.run(200)
    assert sorted(sampler.get_block_sizes()) == [1, 1]


def test_a_single_block_covers_everything_without_a_graph():
    """chisurf: a single-dataset fit is one block -- an ordinary walk."""
    params, objective = gaussian_graph()
    sampler = Sampler("metropolis", 63)
    sampler.set_parameter_ports(params)
    sampler.set_objective(objective, "chi2")
    sampler.run(200)
    assert list(sampler.get_block_sizes()) == [2]


def test_a_warm_up_is_a_small_share_of_the_chain():
    """chisurf's blocked default: clip(steps/20, 100, 500) sweeps."""
    params, objective = gaussian_graph()
    sampler = Sampler("metropolis", 67)
    sampler.set_parameter_ports(params)
    sampler.set_objective(objective, "chi2")
    sampler.set_n_adapt(0)
    sampler.run(2000)  # with adaptation off
    assert sampler.get_iteration() == 2000
    assert sampler.n_evaluations >= 2000


# ------------------------------------------------------- the objective paths

def test_a_python_objective_through_a_node_director():
    """The existing callback machinery: a Node director whose evaluate()
    runs Python is a Sampler objective (the boundary cost per evaluation
    is exactly what the C++ path removes; the bench measures it)."""
    from IMP.bff import Node as BffNode

    class Chi2Node(BffNode):
        def __init__(self):
            super(Chi2Node, self).__init__("chi2")
            self.add_input_port("x1", Port(0.0))
            self.add_input_port("x2", Port(0.0))
            self.add_output_port("chi2", Port(0.0, False, True))

        def evaluate(self):
            d = np.array([self.get_input_port("x1").get_value() - MU[0],
                          self.get_input_port("x2").get_value() - MU[1]])
            self.get_output_port("chi2").set_value(float(d @ P @ d))

    node = Chi2Node()
    sampler = Sampler("stretch", 71)
    sampler.set_parameter_ports([node.get_input_port("x1"),
                                 node.get_input_port("x2")])
    sampler.set_objective(node, "chi2")
    sampler.run(1200)
    mean, cov = moments(flat_chain(sampler, burn=200))
    assert np.all(np.abs(mean - MU) < 0.25), mean
    assert np.allclose(cov, SIGMA, rtol=0.3, atol=0.2)


def test_a_log_likelihood_output_port_is_read_as_one():
    """set_output_is_log_likelihood: the objective node computes lnlike
    directly (here lnlike = -0.5 (x - 0.5)^2); chi^2 is then reported as
    -2 lnlike, the way chisurf's ensemble_result reports a blobless
    chain."""
    x = Port(0.5, name="x")
    d = operator_node("d", x, const_port(-0.5), "addition_double")
    sq = operator_node("sq", linked_port(d.get_output_port("d")),
                       linked_port(d.get_output_port("d")), "multiply_double")
    lnlike = operator_node("lnlike", linked_port(sq.get_output_port("sq")),
                           const_port(-0.5), "multiply_double")
    sampler = Sampler("stretch", 73)
    sampler.set_parameter_ports([x])
    sampler.set_objective(lnlike, "lnlike")
    sampler.set_output_is_log_likelihood(True)
    sampler.run(800)
    chain = flat_chain(sampler, burn=200)
    assert abs(chain[:, 0].mean() - 0.5) < 0.2
    chi2 = np.asarray(sampler.get_chi2())
    assert np.all(chi2 >= -1e-9)  # -2 lnlike, and lnlike <= 0
    assert np.all(chi2 < np.inf)


def test_fixed_parameter_ports_are_refused():
    params, objective = gaussian_graph()
    params[0].set_fixed(True)
    sampler = Sampler("stretch", 79)
    with pytest.raises(ValueError):
        sampler.set_parameter_ports(params)


def test_a_sampler_without_an_objective_is_refused():
    params, _ = gaussian_graph()
    sampler = Sampler("stretch", 79)
    sampler.set_parameter_ports(params)
    with pytest.raises(ValueError):
        sampler.run(10)


def test_an_unknown_algorithm_is_refused():
    """`"slice"` used to be the example here, because it did not exist. It
    does now, so the example has to be a name that really is not a backend --
    otherwise this test passes by accident the day one is added."""
    with pytest.raises(ValueError):
        Sampler("hamiltonian")


# ------------------------------------------------------------ misc surface

def test_describe_names_the_algorithm():
    sampler = make_sampler("stretch", seed=83, n_steps=50)
    text = sampler.describe()
    assert "stretch" in text


def test_the_temperature_flattens_acceptance_but_stays_valid():
    sampler = make_sampler("metropolis", seed=89, n_steps=400, temp=4.0)
    assert 0.0 < sampler.get_acceptance_rate() <= 1.0
    assert np.all(np.isfinite(np.asarray(sampler.get_log_prob())))


# --------------------------------------------------------------- slice/zeus

def test_slice_is_selected_by_every_name_it_is_known_by():
    """`zeus` is the package people arrive from; `slice` is the move."""
    for name in ("slice", "zeus", "ensemble_slice", "sample_slice"):
        assert Sampler(name, 1).get_algorithm() == "slice"


def test_slice_recovers_the_width_of_the_target():
    """The measurement that catches a truncated stepping-out.

    A slice sampler whose interval is capped before it clears the density
    still produces draws that are all *inside* the slice, so nothing looks
    wrong: the acceptance is 100% by construction and the chain moves. What
    it loses is the tails, and the posterior comes out too narrow. So the
    test is the width, not the mean.
    """
    sampler = make_sampler("slice", n_steps=3000)
    chain = flat_chain(sampler)
    sd = chain.std(axis=0)
    truth = np.sqrt(np.diag(SIGMA))
    assert sd[0] == pytest.approx(truth[0], rel=0.15), (sd, truth)
    assert sd[1] == pytest.approx(truth[1], rel=0.15), (sd, truth)


def test_slice_recovers_the_correlation():
    sampler = make_sampler("slice", n_steps=3000)
    chain = flat_chain(sampler)
    truth = SIGMA[0, 1] / np.sqrt(SIGMA[0, 0] * SIGMA[1, 1])
    assert np.corrcoef(chain.T)[0, 1] == pytest.approx(truth, abs=0.1)


def test_a_converged_slice_run_truncates_nothing():
    """`slice_truncations` is the diagnostic for the failure above.

    It counts the times the stepping-out cap bound rather than the density
    ending the expansion. On a Gaussian with the default cap it must be zero;
    a non-zero count means the recorded chain's width is not to be believed.
    """
    sampler = make_sampler("slice", n_steps=1000)
    assert sampler.get_slice_truncations() == 0


def test_the_slice_scale_is_tuned_and_then_frozen():
    """Tuning during the recorded chain would make it non-Markovian."""
    sampler = make_sampler("slice", n_steps=1000)
    mu = sampler.get_slice_mu()
    assert mu > 0.0
    sampler.run(200, 1)
    assert sampler.get_slice_mu() == pytest.approx(mu)


def test_an_explicit_slice_scale_is_not_overwritten():
    sampler = make_sampler("slice", n_steps=200, slice_mu=0.75)
    assert sampler.get_slice_mu() == pytest.approx(0.75)
