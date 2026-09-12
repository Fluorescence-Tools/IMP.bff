"""Native FCS topology and canonical-registry contracts for model search."""

import numpy as np
import pytest

import IMP.bff as bff


def _curve_2d(axis, n=2.5, baseline=1.0, td=0.7):
    return baseline + (1.0 / n) / (1.0 + axis / td)


def test_factory_builds_complete_analytical_family_over_one_registry():
    axis = np.geomspace(1.0e-3, 20.0, 80)
    data = _curve_2d(axis)

    problem = bff.FCSModelSearchFactory.create_analytical(
        list(axis), list(data), [0.01] * len(axis)
    )

    assert list(problem.get_parameter_ids()) == [
        "fcs.N",
        "fcs.baseline",
        "fcs.structure_parameter",
        "fcs.diffusion_time.1",
        "fcs.diffusion_fraction.1",
        "fcs.diffusion_time.2",
        "fcs.relaxation_amplitude.1",
        "fcs.relaxation_time.1",
    ]
    assert set(problem.get_structure_keys()) == {
        f"fcs.{dimension}d.{components}diff.{relaxations}relax"
        for dimension in (2, 3)
        for components in (1, 2)
        for relaxations in (0, 1)
    }

    # Each topology reads the exact same owner port through a follower. A
    # topology switch therefore cannot transfer or name-match parameters.
    owner = problem.get_parameter("fcs.N")
    for key in problem.get_structure_keys():
        objective = problem.get_structure_objective(key)
        model = objective.get_input_port("model").link.get_node()
        assert model is not None  # the problem owns the complete graph
        follower = model.get_input_port("N")
        assert follower.link.uid == owner.uid


def test_single_2d_topology_fits_without_a_python_model_callback():
    axis = np.geomspace(1.0e-3, 20.0, 100)
    data = _curve_2d(axis, n=3.2, baseline=0.97, td=0.42)
    config = bff.FCSModelSearchConfig()
    config.set_include_3d(False)
    config.set_max_diffusion_components(1)
    config.set_max_relaxation_terms(0)
    problem = bff.FCSModelSearchFactory.create_analytical(
        list(axis), list(data), [0.002] * len(axis), config
    )
    search_config = bff.ModelSearchConfig()
    search_config.set_number_of_simulations(12)
    search_config.set_dirichlet_fraction(0.0)
    search_config.set_seed(11)
    search = bff.ModelSearch(problem)
    search.set_config(search_config)

    result = search.run()

    assert result.get_best_state().get_structure_key() == "fcs.2d.1diff.0relax"
    assert abs(problem.get_parameter("fcs.N").value - 3.2) < 5.0e-3
    assert abs(problem.get_parameter("fcs.baseline").value - 0.97) < 3.0e-4
    assert abs(problem.get_parameter("fcs.diffusion_time.1").value - 0.42) < 1.0e-3
    objective = problem.get_active_objective()
    objective.update()
    residuals = np.asarray(objective.get_output_port("residuals").value)
    assert np.max(np.abs(residuals * 0.002)) < 1.5e-4


def test_factory_rejects_incomplete_data_and_an_empty_mode_set():
    with pytest.raises((ValueError, RuntimeError)):
        bff.FCSModelSearchFactory.create_analytical([1.0], [], [])

    config = bff.FCSModelSearchConfig()
    config.set_include_2d(False)
    config.set_include_3d(False)
    with pytest.raises((ValueError, RuntimeError)):
        bff.FCSModelSearchFactory.create_analytical(
            [1.0, 2.0], [1.2, 1.1], [0.1, 0.1], config
        )
