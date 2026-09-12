"""Native TCSPC model-family construction and structural search."""

import numpy as np

import IMP.bff as bff


def _synthetic_two_component_data():
    n = 192
    dt = 0.05
    period = 12.5
    x = np.arange(n) * dt
    irf = np.exp(-0.5 * ((x - 0.8) / 0.08) ** 2)
    source = bff.TCSPCDecay("source")
    source.set_number_of_lifetimes(2)
    source.add_output_port("source", bff.GraphPort([0.0], False, True))
    source.set_response_array(np.ascontiguousarray(irf))
    source.set_timing(dt, period)
    source.set_convolution_range(n, n)
    source.get_input_port("a0").value = 0.65
    source.get_input_port("t0").value = 0.7
    source.get_input_port("a1").value = 0.35
    source.get_input_port("t1").value = 3.8
    source.get_input_port("n0").value = 20000.0
    source.set_normalize_amplitudes(True)
    source.update()
    clean = np.asarray(source.get_output_port("source").value)
    data = bff.FitDataset()
    data.set_values_array(np.ascontiguousarray(clean))
    data.set_noise_family(bff.FIT_NOISE_FAMILY_POISSON)
    return data, irf, dt, period


def _space(maximum=3):
    data, irf, dt, period = _synthetic_two_component_data()
    factory = bff.TCSPCLifetimeSearchFactory()
    factory.set_dataset(data)
    factory.set_response(list(irf))
    factory.set_timing(dt, period)
    factory.set_component_range(1, maximum)
    factory.set_parameter("instrument.n0", 15000.0, True, 0.0, 1e6)
    factory.set_parameter("instrument.background", 0.0, False, 0.0, 1e5)
    return factory.build()


def test_every_topology_reads_one_canonical_parameter_registry():
    space = _space()
    assert list(space.get_structure_keys()) == [
        "lifetime.components.1",
        "lifetime.components.2",
        "lifetime.components.3",
    ]
    tau0 = space.get_parameter("lifetime.tau.0")
    first = space.get_decay("lifetime.components.1")
    third = space.get_decay("lifetime.components.3")
    assert first.get_input_port("t0").link.uid == tau0.uid
    assert third.get_input_port("t0").link.uid == tau0.uid
    assert first.get_output_port(first.name).get_node().uid == first.uid


def test_search_uses_native_objectives_and_activates_a_complete_topology():
    space = _space(maximum=2)
    search = bff.ModelSearch(space.get_problem())
    config = bff.ModelSearchConfig()
    config.set_number_of_simulations(40)
    config.set_dirichlet_fraction(0.0)
    config.set_seed(3)
    search.set_config(config)
    result = search.run()

    assert result.get_best_state().get_structure_key() == "lifetime.components.2"
    active = space.get_problem().get_active_objective()
    assert active is not None
    active.update()
    assert np.isfinite(np.asarray(active.get_output_port("residuals").value)).all()


def test_factory_rejects_unknown_parameters_instead_of_falling_back():
    data, irf, dt, period = _synthetic_two_component_data()
    factory = bff.TCSPCLifetimeSearchFactory()
    factory.set_dataset(data)
    factory.set_response(list(irf))
    factory.set_timing(dt, period)
    factory.set_parameter("chisurf.callback", 1.0, True, 0.0, 2.0)
    with np.testing.assert_raises(ValueError):
        factory.build()
