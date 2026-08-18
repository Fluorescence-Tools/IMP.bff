"""Photon and decay-curve Monte-Carlo against a per-frame quenching rate.

Two entry points, both driven by a quenching rate sampled along a trajectory:
:func:`simulate_photon_trace` samples individual photons (giving a photon
histogram with the counting statistics of a real measurement), and
:func:`simulate_quenched_decay` integrates the excited-state population directly
(giving a smooth curve for the same model, much faster, no shot noise).

Moved here from QuEst (``quest/core/photon.py``) by PRD-109.
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

import IMP.bff

__all__ = ["simulate_photon_trace", "simulate_quenched_decay"]

def simulate_photon_trace(
    n_ph: int,
    k_quench,
    t_step: float = 0.01,
    tau0: float = 0.25,
    random_seed: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Monte-Carlo photon trace against a per-frame quenching rate.

    :param n_ph: number of excitation events to simulate.
    :param k_quench: ``(n_frames,)`` quenching rate in 1/ns along the trajectory.
    :param t_step: trajectory time step in ns.
    :param tau0: intrinsic (unquenched) lifetime in ns.
    :param random_seed: makes the trace reproducible, at the cost of running
        single-threaded (see :func:`_trace_seeded`).
    :returns: ``(delay_times, emitted)`` -- the sampled delay time of each event
        and a ``uint8`` flag saying whether a photon was emitted or the dye was
        quenched first. Delay time is 0 for quenched events.
    """
    emitted = IMP.bff.VectorInt()
    dts = IMP.bff.photon_trace(
        int(n_ph), np.ascontiguousarray(np.asarray(k_quench, dtype=np.float64).ravel()),
        float(t_step), float(tau0),
        -1 if random_seed is None else int(random_seed), emitted)
    return (np.asarray(dts, dtype=np.float64),
            np.asarray(emitted, dtype=np.uint8))


def simulate_quenched_decay(
    n_curves: int,
    decay: np.ndarray,
    dt_tac: float,
    k_quench,
    t_step: float,
    tau0: float,
    random_seed: Optional[int] = None,
) -> None:
    """Fill *decay* in place with a trajectory-driven fluorescence decay curve.

    For each of *n_curves* virtual trajectories, a random starting frame is
    picked along *k_quench* and the excited-state population is propagated in
    steps of *t_step* under a total rate ``1/tau0 + k_quench[frame]``, with the
    emitted intensity accumulated into TAC bins of width *dt_tac*.

    The result is the smooth decay of the same model
    :func:`simulate_photon_trace` samples photon-by-photon -- no shot noise, and
    far cheaper when only the curve is wanted.

    :param n_curves: number of virtual decay curves to average over.
    :param decay: 1-D array filled in place.
    :param dt_tac: TAC bin width in ns.
    :param k_quench: ``(n_frames,)`` quenching rate in 1/ns per trajectory frame.
    :param t_step: trajectory time step in ns.
    :param tau0: intrinsic excited-state lifetime in ns.
    :param random_seed: makes the curve reproducible. The numba this replaced
        had no seed at all -- the starting frames were drawn from the global
        state -- so the function could not be pinned in a test.
    """
    decay += np.asarray(IMP.bff.quenched_decay(
        int(n_curves), int(np.asarray(decay).shape[0]), float(dt_tac),
        np.ascontiguousarray(np.asarray(k_quench, dtype=np.float64).ravel()),
        float(t_step), float(tau0),
        -1 if random_seed is None else int(random_seed)), dtype=np.float64)
