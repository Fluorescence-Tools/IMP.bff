"""Focused contract tests for the model-independent native MCTS core."""

import numpy as np

import IMP.bff as bff


def _problem():
    problem = bff.TabularModelSearchProblem()
    for key, reward, acceptable in (
        ("root", 0.0, False),
        ("small", 2.0, False),
        ("large", 1.0, False),
        ("winner", 5.0, True),
    ):
        problem.add_state(key, reward, acceptable)
    problem.set_initial_state("root")
    problem.add_action("root", "grow-small", "small", 0.5)
    problem.add_action("root", "grow-large", "large", 0.5)
    problem.add_action("small", "refine", "winner", 0.9)
    problem.add_action("small", "undo", "root", 0.1)
    problem.add_action("large", "stop", "large", 1.0, True)
    problem.add_action("winner", "undo", "root", 0.5)
    problem.add_action("winner", "stop", "winner", 0.5, True)
    return problem


def _search(problem, simulations=80):
    config = bff.ModelSearchConfig()
    config.set_number_of_simulations(simulations)
    config.set_reward_scale(1.0)
    config.set_dirichlet_fraction(0.0)
    config.set_seed(7)
    search = bff.ModelSearch(problem)
    search.set_config(config)
    return search


def test_native_search_is_lazy_reproducible_and_selects_best_reward():
    first_problem = _problem()
    first = _search(first_problem).run()
    second = _search(_problem()).run()

    assert first.get_best_state().get_key() == "winner"
    assert first.get_improvement() == 5.0
    assert first.get_acceptable()
    assert first.get_number_of_simulations() == 80
    assert list(first.get_best_path()) == list(second.get_best_path())
    # Only three non-root states exist. Traversals revisit cached nodes rather
    # than re-running the evaluator on every simulation.
    assert first_problem.get_number_of_evaluations() <= 3
    assert first.get_number_of_states_evaluated() == (
        1 + first_problem.get_number_of_evaluations()
    )


def test_cancel_before_run_returns_the_accepted_root_without_evaluation():
    problem = _problem()
    search = _search(problem)
    search.request_cancel()
    result = search.run()

    assert result.get_cancelled()
    assert result.get_number_of_simulations() == 0
    assert result.get_best_state().get_key() == "root"
    assert problem.get_number_of_evaluations() == 0


def _linear_fit_problem():
    x = np.linspace(-2.0, 2.0, 21)
    y = 2.0 * x + 3.0

    model = bff.GraphExpression("linear_model")
    model.set_expression("a*x+b")
    a = bff.GraphPort(2.0)
    b = bff.GraphPort(0.0, True)
    axis = bff.GraphPort(list(x))
    model.add_input_port("a", a)
    model.add_input_port("b", b)
    model.add_input_port("x", axis)
    curve = bff.GraphPort([0.0], False, True)
    model.add_output_port("linear_model", curve)

    objective = bff.FitChiSquared("linear_fit")
    objective.set_data_arrays(np.ascontiguousarray(y), np.ones(y.size))
    model_input = bff.GraphPort([0.0])
    model_input.link = curve
    objective.add_input_port("model", model_input)
    objective.add_output_port("linear_fit", bff.GraphPort(0.0, False, True))
    objective.add_output_port("residuals", bff.GraphPort([0.0], False, True))
    objective._graph = (model, curve, model_input, axis)

    problem = bff.FittingModelSearchProblem(objective, "residuals")
    problem.add_parameter_group("slope", [a])
    problem.add_parameter_group("intercept", [b], [1.0])
    problem.add_structure("root", ["slope"])
    problem.add_structure("with-intercept", ["slope", "intercept"])
    problem.set_initial_structure("root")
    problem.add_action("root", "enable-intercept", "with-intercept", 0.9)
    problem.add_action("root", "stop", "root", 0.1, True)
    problem.add_action("with-intercept", "disable-intercept", "root", 0.1)
    problem.add_action("with-intercept", "stop", "with-intercept", 0.9, True)
    problem.set_complexity_penalty(1.0)
    problem._keepalive = (objective, model, a, b)
    return problem, a, b


def test_live_fit_adapter_caches_masks_and_activates_the_winner():
    problem, a, b = _linear_fit_problem()
    search = _search(problem, simulations=40)
    result = search.run()

    assert result.get_best_state().get_structure_key() == "with-intercept"
    assert abs(a.value - 2.0) < 1e-8
    assert abs(b.value - 3.0) < 1e-8
    assert not a.fixed and not b.fixed
    assert list(problem.get_cached_fixed("root")) == [0, 1]
    assert list(problem.get_cached_fixed(result.get_best_state().get_key())) == [0, 0]
    problem.restore_state("root")
    assert a.value == 2.0 and b.value == 0.0
    assert not a.fixed and b.fixed


def test_live_fit_adapter_cancellation_restores_the_root_snapshot():
    problem, a, b = _linear_fit_problem()
    search = _search(problem, simulations=40)
    search.request_cancel()
    result = search.run()

    assert result.get_cancelled()
    assert result.get_best_state().get_key() == "root"
    assert a.value == 2.0 and b.value == 0.0
    assert not a.fixed and b.fixed


def test_live_fit_adapter_rolls_back_a_failed_candidate_to_its_parent():
    problem, a, b = _linear_fit_problem()
    extras = [bff.GraphPort(0.0, True) for _ in range(22)]
    problem.add_parameter_group("too-many", extras, [0.0] * len(extras))
    problem.add_structure("overparameterized", ["slope", "too-many"])
    problem.add_action("root", "overparameterize", "overparameterized", 1.0)
    root = problem.get_initial_state()
    action = next(
        action for action in problem.get_actions(root)
        if action.get_key() == "overparameterize"
    )

    failed = problem.evaluate(root, action)

    assert failed.get_key() == root.get_key()
    assert problem.get_last_fit_status() == 0
    assert problem.get_last_failure()
    assert a.value == 2.0 and b.value == 0.0
    assert not a.fixed and b.fixed
    assert all(port.fixed for port in extras)


def _joint_fit_problem():
    x = np.linspace(0.1, 8.0, 48)
    ys = (2.0 * np.exp(-x / 3.0), 5.0 * np.exp(-x / 3.0))
    members = []
    parameters = []
    for index, y in enumerate(ys):
        name = "joint_%d" % index
        model = bff.GraphExpression(name + "_model")
        model.set_expression("a*exp(-x/t)")
        amplitude = bff.GraphPort(float((2.0, 5.0)[index]))
        lifetime = bff.GraphPort(1.0, index == 0)
        axis = bff.GraphPort(list(x))
        model.add_input_port("a", amplitude)
        model.add_input_port("t", lifetime)
        model.add_input_port("x", axis)
        curve = bff.GraphPort([0.0], False, True)
        model.add_output_port(name + "_model", curve)
        objective = bff.FitChiSquared(name)
        objective.set_data_arrays(np.ascontiguousarray(y), np.ones(y.size))
        model_input = bff.GraphPort([0.0])
        model_input.link = curve
        objective.add_input_port("model", model_input)
        objective.add_output_port(name, bff.GraphPort(0.0, False, True))
        objective.add_output_port("residuals", bff.GraphPort([0.0], False, True))
        objective._graph = (model, curve, model_input, axis)
        members.append(objective)
        parameters.append((amplitude, lifetime))
    parameters[1][1].link = parameters[0][1]

    joint = bff.FitJointChiSquared("joint")
    joint.add_output_port("joint", bff.GraphPort(0.0, False, True))
    joint.add_output_port("residuals", bff.GraphPort([0.0], False, True))
    for member in members:
        joint.add_member(member, "residuals")

    problem = bff.FittingModelSearchProblem(joint, "residuals")
    problem.add_parameter_group("amplitudes", [parameters[0][0], parameters[1][0]])
    problem.add_parameter_group("lifetime", [parameters[0][1]], [2.0])
    problem.add_structure("fixed-lifetime", ["amplitudes"])
    problem.add_structure("free-lifetime", ["amplitudes", "lifetime"])
    problem.set_initial_structure("fixed-lifetime")
    problem.add_action("fixed-lifetime", "free-lifetime", "free-lifetime", 1.0)
    problem.add_action("free-lifetime", "stop", "free-lifetime", 1.0, True)
    problem.set_complexity_penalty(0.01)
    problem._keepalive = (joint, members, parameters)
    return problem, parameters


def test_live_adapter_treats_a_joint_objective_like_any_other_graph():
    problem, parameters = _joint_fit_problem()
    result = _search(problem, simulations=20).run()

    assert result.get_best_state().get_structure_key() == "free-lifetime"
    assert abs(parameters[0][1].value - 3.0) < 1e-5
    assert parameters[1][1].value == parameters[0][1].value


def _expression_objective(name, expression, x, y, owners):
    model = bff.GraphExpression(name + "_model")
    model.set_expression(expression)
    followers = []
    for key, owner in owners.items():
        follower = bff.GraphPort(owner.value)
        follower.link = owner
        model.add_input_port(key, follower)
        followers.append(follower)
    axis = bff.GraphPort(list(x))
    model.add_input_port("x", axis)
    curve = bff.GraphPort([0.0], False, True)
    model.add_output_port(name + "_model", curve)

    objective = bff.FitChiSquared(name)
    objective.set_data_arrays(np.ascontiguousarray(y), np.ones(len(y)))
    model_input = bff.GraphPort([0.0])
    model_input.link = curve
    objective.add_input_port("model", model_input)
    objective.add_output_port(name, bff.GraphPort(0.0, False, True))
    objective.add_output_port("residuals", bff.GraphPort([0.0], False, True))
    objective._graph = (model, followers, axis, curve, model_input)
    return objective


def _multi_structure_linear_problem():
    x = np.linspace(-2.0, 2.0, 25)
    y = 2.0 * x + 3.0
    slope = bff.GraphPort(2.0)
    intercept = bff.GraphPort(0.0, True)
    root_objective = _expression_objective(
        "through_origin", "a*x", x, y, {"a": slope}
    )
    expanded_objective = _expression_objective(
        "with_intercept", "a*x+b", x, y, {"a": slope, "b": intercept}
    )

    problem = bff.MultiStructureModelSearchProblem()
    problem.add_parameter("slope", slope)
    problem.add_parameter("intercept", intercept)
    ids = ["slope", "intercept"]
    ports = [slope, intercept]
    problem.add_structure(
        "root", root_objective, ids, ports, [2.0, 0.0], [0, 1], "residuals"
    )
    # The deliberately bad slope seed proves that a parameter which remains
    # free is warm-started in the shared registry rather than reset per graph.
    problem.add_structure(
        "expanded",
        expanded_objective,
        ids,
        ports,
        [-100.0, 1.0],
        [0, 0],
        "residuals",
    )
    problem.add_structure_node("root", root_objective._graph[0])
    problem.add_structure_node("expanded", expanded_objective._graph[0])
    del root_objective._graph
    del expanded_objective._graph
    problem.set_structure_bic_metadata("root", len(x), 1.0)
    problem.set_structure_bic_metadata("expanded", len(x), 2.0)
    problem.set_initial_structure("root")
    problem.add_action("root", "add-intercept", "expanded", 1.0)
    problem.add_action("expanded", "stop", "expanded", 1.0, True)
    problem._keepalive = (root_objective, expanded_objective, slope, intercept)
    return problem, slope, intercept, root_objective, expanded_objective


def test_multi_structure_problem_rejects_noncanonical_or_incomplete_state():
    owner = bff.GraphPort(1.0)
    other = bff.GraphPort(1.0)
    objective = _expression_objective(
        "identity_contract", "a*x", [1.0, 2.0], [1.0, 2.0], {"a": owner}
    )
    problem = bff.MultiStructureModelSearchProblem()
    problem.add_parameter("a", owner)

    with np.testing.assert_raises(ValueError):
        problem.add_structure(
            "missing", objective, [], [], [], [], "residuals"
        )
    with np.testing.assert_raises(ValueError):
        problem.add_structure(
            "copied-owner", objective, ["a"], [other], [1.0], [0], "residuals"
        )
    with np.testing.assert_raises(ValueError):
        problem.add_parameter("second-name", owner)


def test_multi_structure_search_switches_graphs_and_activates_canonical_winner():
    problem, slope, intercept, root_objective, expanded_objective = (
        _multi_structure_linear_problem()
    )
    root = problem.get_initial_state()
    assert list(problem.get_cached_values(root.get_key())) == [2.0, 0.0]
    assert problem.get_active_objective().get_uid() == root_objective.get_uid()

    result = _search(problem, simulations=20).run()

    assert result.get_best_state().get_structure_key() == "expanded"
    assert problem.get_active_structure() == "expanded"
    assert problem.get_active_objective().get_uid() == expanded_objective.get_uid()
    assert abs(slope.value - 2.0) < 1e-8
    assert abs(intercept.value - 3.0) < 1e-8
    assert not slope.fixed and not intercept.fixed

    problem.restore_state("root")
    assert problem.get_active_structure() == "root"
    assert slope.value == 2.0 and intercept.value == 0.0
    assert not slope.fixed and intercept.fixed


def test_multi_structure_cancel_rolls_back_registry_and_active_graph():
    problem, slope, intercept, root_objective, _ = (
        _multi_structure_linear_problem()
    )
    root = problem.get_initial_state()
    action = problem.get_actions(root)[0]
    problem.request_cancel()

    collapsed = problem.evaluate(root, action)

    assert collapsed.get_key() == root.get_key()
    assert problem.get_last_fit_status() == -1
    assert problem.get_last_failure()
    assert problem.get_active_structure() == "root"
    assert problem.get_active_objective().get_uid() == root_objective.get_uid()
    assert slope.value == 2.0 and intercept.value == 0.0
    assert not slope.fixed and intercept.fixed


def test_multi_structure_failed_fit_rolls_back_registry_and_active_graph():
    x = np.array([1.0])
    slope = bff.GraphPort(1.0)
    intercept = bff.GraphPort(0.0, True)
    root_objective = _expression_objective(
        "failure_root", "a*x", x, x, {"a": slope}
    )
    invalid_objective = _expression_objective(
        "failure_target", "a*x+b", x, x, {"a": slope, "b": intercept}
    )
    problem = bff.MultiStructureModelSearchProblem()
    problem.add_parameter("slope", slope)
    problem.add_parameter("intercept", intercept)
    ids = ["slope", "intercept"]
    ports = [slope, intercept]
    problem.add_structure(
        "root", root_objective, ids, ports, [1.0, 0.0], [0, 1]
    )
    # Two free parameters cannot be fitted to one residual.
    problem.add_structure(
        "underdetermined", invalid_objective, ids, ports, [50.0, 20.0], [0, 0]
    )
    problem.set_initial_structure("root")
    problem.add_action("root", "fail", "underdetermined", 1.0)
    problem._keepalive = (root_objective, invalid_objective, slope, intercept)
    root = problem.get_initial_state()

    collapsed = problem.evaluate(root, problem.get_actions(root)[0])

    assert collapsed.get_key() == root.get_key()
    assert problem.get_last_fit_status() == 0
    assert problem.get_last_failure()
    assert problem.get_active_structure() == "root"
    assert slope.value == 1.0 and intercept.value == 0.0
    assert not slope.fixed and intercept.fixed


def test_multi_structure_problem_accepts_a_joint_target_objective():
    x = np.linspace(0.1, 8.0, 48)
    lifetime = bff.GraphPort(1.0, True)
    amplitudes = [bff.GraphPort(2.0), bff.GraphPort(5.0)]
    members = []
    for index, amplitude in enumerate(amplitudes):
        y = amplitude.value * np.exp(-x / 3.0)
        members.append(
            _expression_objective(
                "multi_joint_%d" % index,
                "a*exp(-x/t)",
                x,
                y,
                {"a": amplitude, "t": lifetime},
            )
        )
    joint = bff.FitJointChiSquared("multi_joint")
    joint.add_output_port("multi_joint", bff.GraphPort(0.0, False, True))
    joint.add_output_port("residuals", bff.GraphPort([0.0], False, True))
    for member in members:
        joint.add_member(member, "residuals")

    problem = bff.MultiStructureModelSearchProblem()
    problem.add_parameter("amplitude-0", amplitudes[0])
    problem.add_parameter("amplitude-1", amplitudes[1])
    problem.add_parameter("lifetime", lifetime)
    ids = ["amplitude-0", "amplitude-1", "lifetime"]
    ports = [amplitudes[0], amplitudes[1], lifetime]
    problem.add_structure(
        "fixed-lifetime", members[0], ids, ports, [2.0, 5.0, 1.0], [0, 1, 1]
    )
    problem.add_structure(
        "joint-free-lifetime", joint, ids, ports, [2.0, 5.0, 2.0], [0, 0, 0]
    )
    problem.add_structure_node("fixed-lifetime", members[0]._graph[0])
    for member in members:
        problem.add_structure_node("joint-free-lifetime", member._graph[0])
        del member._graph
    problem.set_initial_structure("fixed-lifetime")
    problem.add_action(
        "fixed-lifetime", "fit-joint", "joint-free-lifetime", 1.0
    )
    problem.add_action(
        "joint-free-lifetime", "stop", "joint-free-lifetime", 1.0, True
    )
    problem._keepalive = (joint, members, amplitudes, lifetime)

    result = _search(problem, simulations=20).run()

    assert result.get_best_state().get_structure_key() == "joint-free-lifetime"
    assert abs(lifetime.value - 3.0) < 1e-5
    assert problem.get_active_objective().get_uid() == joint.get_uid()
