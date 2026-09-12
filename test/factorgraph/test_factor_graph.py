"""The factor structure of a model's posterior, in bff (ported from PRD-68).

ChiSurf's factor graph made a global fit's factorisation explicit: which
datasets a parameter touches (relevance), how the fit decomposes into blocks
and separators, and what its treewidth is. These tests pin the same behaviour
for the bff C++ implementation -- the layer that lets a complex MDL
(spectroscopy bound to structures) carry its own structure, standalone, with
no IMP particles taking part in the graph.
"""
import pytest

from IMP.bff import InferenceFactorGraph, INFERENCE_FACTOR_LIKELIHOOD, INFERENCE_FACTOR_PRIOR


def single_dataset_graph(n_locals=3):
    """One dataset whose likelihood reads every parameter: a complete graph."""
    g = InferenceFactorGraph()
    keys = []
    for i in range(n_locals):
        key = f"x{i}"
        keys.append(key)
        g.add_variable(key, key, i, 0)
    g.add_factor("L0", INFERENCE_FACTOR_LIKELIHOOD, keys, 0, 32)
    for key in keys:
        g.add_factor(f"pi_{key}", INFERENCE_FACTOR_PRIOR, [key], -1, 1)
    return g, keys


def star_graph(n_datasets=4, n_globals=2, n_locals=2):
    """Datasets with local parameters sharing a few globals: a star."""
    g = InferenceFactorGraph()
    index = 0
    local_keys = {}
    for k in range(n_datasets):
        for j in range(n_locals):
            key = f"fit{k}:x{j}"
            g.add_variable(key, f"x{j}", index, k)
            local_keys.setdefault(k, []).append(key)
            index += 1
    global_keys = []
    for j in range(n_globals):
        key = f"g:g{j}"
        g.add_variable(key, f"g{j}", index, -1)
        global_keys.append(key)
        index += 1
    for k in range(n_datasets):
        g.add_factor(f"L{k}", INFERENCE_FACTOR_LIKELIHOOD, local_keys[k] + global_keys, k, 32)
    for keys in local_keys.values():
        for key in keys:
            g.add_factor(f"pi_{key}", INFERENCE_FACTOR_PRIOR, [key], -1, 1)
    return g, local_keys, global_keys


def unlinked_group(n_datasets=3, n_locals=2):
    """Datasets that share nothing: n independent sub-problems."""
    g = InferenceFactorGraph()
    index = 0
    local_keys = {}
    for k in range(n_datasets):
        for j in range(n_locals):
            key = f"fit{k}:x{j}"
            g.add_variable(key, f"x{j}", index, k)
            local_keys.setdefault(k, []).append(key)
            index += 1
        g.add_factor(f"L{k}", INFERENCE_FACTOR_LIKELIHOOD, local_keys[k], k, 32)
    return g, local_keys


def test_single_dataset_is_one_clique():
    g, keys = single_dataset_graph(3)
    assert g.get_number_of_likelihood_factors() == 1
    cliques = g.get_cliques()
    assert len(cliques) == 1
    assert set(cliques[0]) == set(keys)
    # every parameter coupled by one likelihood: n_free - 1
    assert g.get_treewidth() == len(keys) - 1


def test_complete_graph_shortcut_agrees_with_the_result():
    # A complete graph has a single maximal clique whatever the order.
    g, keys = single_dataset_graph(5)
    assert g.get_cliques() == (tuple(sorted(keys, key=g.index_of)),)
    assert g.get_elimination_order("min_degree") == g.get_elimination_order(
        "min_fill")


def test_empty_graph_has_no_cliques():
    g = InferenceFactorGraph()
    assert g.get_cliques() == ()
    assert g.get_treewidth() == 0
    assert g.get_number_of_variables() == 0


def test_unknown_heuristic_raises():
    g, _ = single_dataset_graph()
    with pytest.raises(Exception):
        g.get_elimination_order("max_entropy")


def test_duplicate_keys_and_unknown_scope_raise():
    g = InferenceFactorGraph()
    g.add_variable("x", "x", 0, 0)
    with pytest.raises(Exception):
        g.add_variable("x", "x", 1, 0)
    with pytest.raises(Exception):
        g.add_factor("L", INFERENCE_FACTOR_LIKELIHOOD, ["nope"], 0, 1)


def test_linking_a_parameter_creates_a_star_with_one_separator():
    g, local_keys, global_keys = star_graph(n_datasets=4)
    # one clique per dataset: locals + the globals they share
    cliques = g.get_cliques()
    assert len(cliques) == 4
    for k in range(4):
        assert set(cliques[k]) == set(local_keys[k] + global_keys)
    # the separator is exactly the shared globals
    separators = g.get_separators()
    assert len(separators) == 1
    assert set(separators[0]) == set(global_keys)
    # treewidth stays small however many datasets
    assert g.get_treewidth() == 2 + len(global_keys) - 1


def test_unlinked_group_decomposes_into_independent_datasets():
    g, local_keys = unlinked_group(n_datasets=3)
    components = g.connected_components()
    assert len(components) == 3
    for k in range(3):
        assert set(components[k]) == set(local_keys[k])
    assert g.get_separators() == ()


def test_relevance_maps_parameters_to_the_datasets_they_touch():
    g, local_keys, global_keys = star_graph(n_datasets=3)
    # a local parameter reaches exactly its own dataset
    assert list(g.affected_fits(local_keys[1])) == [1]
    # a global reaches every dataset
    assert list(g.affected_fits(global_keys)) == [0, 1, 2]
    # a prior touching one dataset's variable does not add a fit
    assert list(g.affected_fits(local_keys[0][:1])) == [0]
    # block cost is the number of fits a move recomputes
    assert g.block_cost(local_keys[0]) == 1
    assert g.block_cost(global_keys) == 3


def test_relevance_of_an_unlinked_group_is_one_dataset_each():
    g, local_keys = unlinked_group(n_datasets=3)
    for k in range(3):
        assert list(g.affected_fits(local_keys[k])) == [k]


def test_affected_factors_names_the_factors_whose_scope_intersected():
    g, _, global_keys = star_graph(n_datasets=2)
    factors = set(g.affected_factors(global_keys))
    assert "L0" in factors and "L1" in factors
    assert g.variables_of("L0") == g.variables_of("L0")


def test_sampling_blocks_partition_the_variables_by_neighbourhood():
    g, local_keys, global_keys = star_graph(n_datasets=3)
    blocks = g.get_sampling_blocks()
    # disjoint cover
    seen = [v for b in blocks for v in b]
    assert sorted(seen) == sorted(local_keys[0] + local_keys[1] +
                                  local_keys[2] + global_keys)
    # locals per dataset (cheap, one fit each), globals last (all fits)
    assert set(blocks[0]) == set(local_keys[0])
    assert set(blocks[-1]) == set(global_keys)
    # a plain fit has one neighbourhood and so one block
    g1, _ = single_dataset_graph(3)
    assert len(g1.get_sampling_blocks()) == 1


def test_priors_become_factors_but_do_not_couple_variables():
    g, keys = single_dataset_graph(2)
    # single-variable scopes add no moral edges: the graph stays complete
    assert g.get_cliques() == (tuple(sorted(keys, key=g.index_of)),)
    assert g.get_number_of_factors() == 3  # one likelihood, two priors


def test_elimination_orders_are_deterministic_and_repeatable():
    g, _, _ = star_graph(n_datasets=4)
    order = g.get_elimination_order("min_fill")
    assert order == g.get_elimination_order("min_fill")
    assert sorted(order) == sorted(g.get_variable_keys())
    # min_degree is also valid and deterministic (may differ from min_fill)
    order_d = g.get_elimination_order("min_degree")
    assert order_d == g.get_elimination_order("min_degree")
    assert sorted(order_d) == sorted(g.get_variable_keys())


def test_junction_tree_edges_carry_their_separator():
    g, local_keys, global_keys = star_graph(n_datasets=3)
    edges = g.get_junction_tree_edges()
    cliques = g.get_cliques()
    # a tree over 3 cliques: 2 edges, each separator the shared globals
    assert len(edges) == 2
    for e in edges:
        shared = set(cliques[e.first]) & set(cliques[e.second])
        assert set(e.separator) == shared
        assert shared == set(global_keys)


def test_unexplained_variables_are_reported():
    g = InferenceFactorGraph()
    g.add_variable("seen", "seen", 0, 0)
    g.add_variable("blind", "blind", 1, 0)
    g.add_factor("L0", INFERENCE_FACTOR_LIKELIHOOD, ["seen"], 0, 8)
    assert g.get_unexplained_variables() == ("blind",)


def test_describe_reports_the_identifiability_statement():
    g, local_keys, global_keys = star_graph(n_datasets=4, n_globals=2)
    text = g.describe()
    assert "variables      : 10" in text
    assert "likelihoods    : 4" in text
    assert "treewidth      : 3" in text
    assert "components     : 1" in text
    assert "separators     : {g0, g1}" in text


def test_mutation_invalidates_cached_structure():
    g, _, global_keys = star_graph(n_datasets=2)
    treewidth_before = g.get_treewidth()
    # add a variable no factor reads: an unexplained variable joins the graph
    g.add_variable("ghost", "ghost", 99, -1)
    assert g.get_unexplained_variables() == ("ghost",)
    # the ghost is its own component now
    assert len(g.connected_components()) == 2
    assert g.get_treewidth() == treewidth_before


def test_a_four_variable_star_is_not_complete():
    """The completeness shortcut double-counted, and a star hit it exactly.

    `is_complete()` summed the symmetric adjacency -- which is *twice* the
    edge count -- and compared it against ``n(n-1)/2``. That is true whenever
    ``|E| == n(n-1)/4``, and three datasets around one shared parameter is
    four variables with three edges: ``2*3 == 4*3/2``. The graph was then
    taken for a clique, so ``get_cliques()`` returned one block over
    everything and ``get_treewidth()`` said 3 where the answer is 1 -- for
    the shape a global fit *has*.

    The existing star fixture above has ten variables and sixteen edges, so
    it never hit the coincidence; that is why this needs its own test rather
    than a bigger one.
    """
    g = InferenceFactorGraph()
    g.add_variable("shared", "a", 0, -1)
    for k in range(3):
        g.add_variable(f"local{k}", "c", k + 1, k)
        g.add_factor(f"L{k}", INFERENCE_FACTOR_LIKELIHOOD, ["shared", f"local{k}"], k, 32)

    cliques = [tuple(c) for c in g.get_cliques()]
    assert len(cliques) == 3, cliques
    assert all(len(c) == 2 for c in cliques)
    assert g.get_treewidth() == 1
    assert len(g.get_separators()) == 1
    assert list(g.get_separators()[0]) == ["shared"]


def test_a_genuinely_complete_four_variable_graph_still_is():
    """The guard must not overshoot: six edges over four variables is one
    clique, and that is what the shortcut exists for."""
    g, keys = single_dataset_graph(4)
    assert g.get_treewidth() == 3
    assert [tuple(c) for c in g.get_cliques()] == [tuple(keys)]
