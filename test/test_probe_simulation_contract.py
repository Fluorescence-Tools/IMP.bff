"""The simulation/trajectory binding contract in either build lane."""

import json

import numpy as np
import IMP.bff as bff


def test_grid_simulation_returns_the_shared_trajectory_in_femtoseconds():
    simulation = bff.ProbeDiffusionSimulation([1] * 343, 1.0, [0.0] * 3)
    simulation.set_parameters(json.dumps({
        "diffusion_coefficient": 0.0, "t_step": 0.1,
        "n_trajectories": 1, "seed": 7,
    }))
    trajectory = simulation.run(8, 2)
    assert isinstance(simulation, bff.ProbeSimulation)
    assert isinstance(trajectory, bff.ProbeSimulationTrajectory)
    assert trajectory.n_frames == 4
    assert trajectory.n_atoms == 1
    coordinates = np.asarray(trajectory.coordinates)
    assert coordinates.shape == (4, 1, 3)
    # Zero diffusion must leave the chosen starting voxel unchanged.
    np.testing.assert_array_equal(coordinates, np.repeat(coordinates[:1], 4, axis=0))
    np.testing.assert_allclose(trajectory.times_fs, [0, 200000, 400000, 600000])
    assert np.asarray(trajectory.potential_energy).size == 0
    assert np.asarray(trajectory.kinetic_energy).size == 0
