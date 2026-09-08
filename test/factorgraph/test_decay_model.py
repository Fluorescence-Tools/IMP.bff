"""PRD-139: the factor graph carrying a Bayesian decay analysis.

The fixture is not a caricature written here. `graph_spec_ensemble.json` was
emitted by the consumer (`../ucfret`, commit 52e616e) **by measurement** --
perturb each variable, recompute each histogram's expected counts, record what
actually changed -- so the scopes are what the model has rather than what
anyone believed it had. Its own first draft said every likelihood reads every
global variable; measuring said otherwise, and the tests below assert the
measured shape.

What that shape is, and why each part needs something the class did not have:

* **Sizes.** 46 variables holding 113 free numbers. The distance distribution
  is 25 spline coefficients constrained to sum to zero, so its *free*
  dimension is 24 -- and 24, not 25, is what enters a block's cost. Without a
  size every block looked the same width.
* **Roles.** distribution, physics, calibration, instrument, hyper. The
  consumer kept a parallel dictionary because the class had nowhere to put
  them.
* **A HYPER factor.** `log10_lam` appears in no likelihood's scope at all: it
  reaches the data only through the roughness factor that couples it to all 24
  distribution coefficients. That is what makes it a hyperparameter rather
  than a parameter, and as a PRIOR it would be indistinguishable from a
  Gaussian on one constant.
"""

import json
from pathlib import Path

import numpy as np
import pytest

import IMP.bff

SPEC = Path(__file__).resolve().parent / "graph_spec_ensemble.json"


@pytest.fixture(scope="module")
def spec():
    return json.loads(SPEC.read_text(encoding="utf-8"))


def _build(spec):
    g = IMP.bff.FactorGraph()
    offset = 0
    for v in spec["variables"]:
        g.add_variable(v["key"], v.get("name", v["key"]), offset, -1,
                       v.get("size", 1), v.get("role", ""))
        offset += v.get("size", 1)
    kinds = {"PRIOR": IMP.bff.PRIOR, "LIKELIHOOD": IMP.bff.LIKELIHOOD,
             "HYPER": IMP.bff.HYPER}
    for f in spec["factors"]:
        g.add_factor(f["key"], kinds[f["kind"]], f["scope"],
                     f.get("fit_index", -1), f.get("size", 0))
    return g


@pytest.fixture(scope="module")
def graph(spec):
    return _build(spec)


def test_the_model_goes_in_whole(graph, spec):
    assert graph.get_number_of_variables() == spec["n_variables"] == 46
    assert graph.get_number_of_factors() == len(spec["factors"]) == 54
    assert graph.get_number_of_likelihood_factors() == spec["n_histograms"] == 8


def test_the_free_numbers_add_up(graph, spec):
    total = sum(graph.get_variable_size(v["key"]) for v in spec["variables"])
    assert total == spec["unconstrained_dim"] == 113


def test_the_distribution_is_24_free_numbers_not_25(graph):
    """25 spline coefficients, constrained to sum to zero. The free dimension
    is what the graph must carry: it is what block cost, elimination and
    treewidth are about. The sum-to-zero basis is the consumer's transform and
    bff knows nothing of it."""
    assert graph.get_variable_size("c") == 24
    assert graph.get_variable_role("c") == "distribution"


def test_every_variable_carries_its_role(graph, spec):
    for v in spec["variables"]:
        assert graph.get_variable_role(v["key"]) == v["role"]
    roles = {graph.get_variable_role(v["key"]) for v in spec["variables"]}
    assert roles == {"distribution", "physics", "calibration", "instrument", "hyper"}


def test_a_block_costs_its_width(graph):
    """The check that could not fail before, because it could not differ.

    Cost is how many local fits the block's movement forces, times how wide
    the block is. The distribution touches four histograms and is 24 numbers
    across, so its block costs 24 times what a scalar in the same factors
    would."""
    fits = len(graph.affected_fits(["c"]))
    assert fits == 4
    assert graph.block_cost(["c"]) == 24 * fits


def test_the_cost_law_on_a_block_whose_width_is_the_only_difference(graph):
    """The same claim without the real graph's incidental structure: two
    graphs alike but for one variable's size."""
    def one(size):
        g = IMP.bff.FactorGraph()
        g.add_variable("wide", "wide", 0, -1, size)
        g.add_variable("s", "s", size, -1, 1)
        g.add_factor("L", IMP.bff.LIKELIHOOD, ["wide", "s"], 0, 100)
        return g.block_cost(["wide", "s"])
    assert one(24) == 24 * one(1)


def test_a_hyperparameter_forces_no_local_fit(graph):
    """`log10_lam` reaches no histogram, so moving it re-runs no local fit and
    its block cost is zero. That is what the number has always meant -- fits
    forced, not work done -- and it is worth pinning because a sampler reading
    it as "free" would be drawing the wrong conclusion from the right
    number."""
    assert len(graph.affected_fits(["log10_lam"])) == 0
    assert graph.block_cost(["log10_lam"]) == 0


def test_the_hyperparameter_is_not_a_prior(graph):
    assert graph.get_factor_kind("rough:c") == IMP.bff.HYPER
    assert graph.get_number_of_factors_of_kind(IMP.bff.HYPER) == 1
    assert "hyper factors  : 1" in graph.describe()


def test_the_hyperparameter_reaches_the_data_only_through_that_factor(graph, spec):
    """`log10_lam` is in no likelihood's scope. It couples to the data solely
    through the roughness factor -- which is the whole reason it is integrated
    out rather than optimised."""
    for f in spec["factors"]:
        if f["kind"] == "LIKELIHOOD":
            assert "log10_lam" not in f["scope"]
    touching = set(graph.factors_of("log10_lam"))
    assert touching == {"rough:c", "pi:log10_lam"}


def test_the_distribution_is_read_by_one_sample_only(graph, spec):
    """The measured structure, and the reason this is neither a star nor a
    complete graph: the donor-only and acceptor-only histograms say nothing
    about the distance distribution directly."""
    reading = [f["key"] for f in spec["factors"]
               if f["kind"] == "LIKELIHOOD" and "c" in f["scope"]]
    assert len(reading) == 4
    assert len(reading) < spec["n_histograms"]


def test_the_scopes_really_do_differ(spec):
    """No two of the eight histograms read the same set. A graph in which they
    did would decompose quite differently."""
    scopes = [frozenset(f["scope"]) for f in spec["factors"]
              if f["kind"] == "LIKELIHOOD"]
    assert len(set(scopes)) == len(scopes) == 8


def test_the_decomposition_is_computed_over_the_real_scopes(graph):
    comps = graph.connected_components()
    assert len(comps) >= 1
    assert sum(len(c) for c in comps) == graph.get_number_of_variables()
    assert graph.get_treewidth() > 0


def test_a_round_trip_reproduces_the_graph(graph, tmp_path):
    text = graph.to_json()
    back = IMP.bff.FactorGraph()
    back.from_json(text)
    assert back.to_json() == text
    assert back.describe() == graph.describe()
    for key in graph.get_variable_keys():
        assert back.get_variable_size(key) == graph.get_variable_size(key)
        assert back.get_variable_role(key) == graph.get_variable_role(key)
        assert back.index_of(key) == graph.index_of(key)
    for key in graph.get_factor_keys():
        assert back.get_factor_kind(key) == graph.get_factor_kind(key)
        assert back.variables_of(key) == graph.variables_of(key)

    path = tmp_path / "graph.json"
    graph.save(str(path))
    loaded = IMP.bff.FactorGraph()
    loaded.load(str(path))
    assert loaded.to_json() == text


def test_loading_something_else_is_refused():
    g = IMP.bff.FactorGraph()
    with pytest.raises(ValueError):
        g.from_json("{not json")
    with pytest.raises(ValueError):
        g.from_json(json.dumps({"format": "something.else", "variables": [],
                                "factors": []}))


def test_the_defaults_are_the_old_behaviour():
    """Calls written before variables had sizes must mean what they meant."""
    g = IMP.bff.FactorGraph()
    g.add_variable("a", "a", 0)
    g.add_variable("b", "b", 1)
    g.add_factor("L", IMP.bff.LIKELIHOOD, ["a", "b"], 0, 10)
    assert g.get_variable_size("a") == 1
    assert g.get_variable_role("a") == ""
    assert g.block_cost(["a", "b"]) == len(g.affected_fits(["a", "b"]))


def test_work_dirtied_is_a_separate_number_from_fits_forced(graph):
    """The consumer's addition rather than a redefinition: block_cost keeps
    meaning fits forced, and factor_cost answers the other question -- the
    residuals a move dirties, summed over every factor it touches, priors and
    hyper factors included.

    For the hyperparameter that is 24 + 1 rather than 0: it forces no fit, and
    it dirties the roughness factor over the 24 distribution coefficients and
    its own prior."""
    assert graph.block_cost(["log10_lam"]) == 0
    assert graph.factor_cost(["log10_lam"]) == 25

    # and for a variable that does reach the data, it is the histograms it
    # touches plus its own prior
    assert graph.factor_cost(["c"]) > graph.block_cost(["c"])
