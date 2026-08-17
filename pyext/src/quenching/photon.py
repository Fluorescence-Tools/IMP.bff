"""Photon and decay-curve Monte-Carlo against a per-frame quenching rate.

Two entry points, both driven by a quenching rate sampled along a trajectory:
:func:`simulate_photon_trace` samples individual photons (giving a photon
histogram with the counting statistics of a real measurement), and
:func:`simulate_quenched_decay` integrates the excited-state population directly
(giving a smooth curve for the same model, much faster, no shot noise).

Moved here from QuEst (``quest/core/photon.py``) by PRD-109.
"""

from __future__ import annotations

import os
from typing import Optional, Tuple

import numpy as np

from .._jit import njit, prange, get_num_threads

__all__ = ["simulate_photon_trace", "simulate_quenched_decay"]

@njit(cache=True)
def _ranf() -> float:
    return np.random.random()


@njit(cache=True)
def _log_exp_wait(tau0: float) -> float:
    """An exponentially distributed waiting time with mean *tau0*.

    Inverse-transform sampling: ``-tau0 * ln(u)`` for ``u`` uniform on (0, 1].
    ``np.random.random()`` is uniform on **[0, 1)**, so the draw is ``1 - u``,
    which lands in (0, 1] -- never 0, so the log is finite, and never above 1,
    so the result is never negative.

    QuEst computed ``log(1 / (u + EPS)) * tau0`` with ``EPS = 2.4414062e-4``,
    inherited from a C ancestor where the epsilon guarded ``log(1/0)``. But it
    also pushed the argument **below 1** whenever ``u > 1 - EPS``, making the
    logarithm negative: a photon emitted before it was excited. That happened to
    2.4e-4 of all photons -- 8 in a 40 000-photon trace, measured -- and each
    one was then silently dropped by any histogram starting at 0, so the decay
    curve quietly held fewer photons than the trace said were emitted. Fixed on
    the move (PRD-109); the epsilon is gone because the reformulation does not
    need it.
    """
    return -np.log(1.0 - _ranf()) * tau0


@njit(cache=True, nogil=True)
def _photon_walk(kq, t_step, tau0, shift_nbr, n_frames):
    """Emit one photon: sample its lifetime, then race it against quenching."""
    dt = _log_exp_wait(tau0)
    n_step = int(dt / t_step)
    for j in range(shift_nbr, shift_nbr + n_step):
        if _ranf() < kq[j % n_frames] * t_step:
            return 0.0, np.uint8(0)
    return dt, np.uint8(1)


@njit(cache=True, parallel=True)
def _seed_worker_threads(base_seed: int) -> None:
    """Give every numba worker thread a fresh RNG state.

    Numba keeps one RNG state per thread and initialises it deterministically,
    so an *unseeded* run repeats itself for the whole life of the process --
    four Monte-Carlo runs in a row return byte-identical results, which reads as
    precision that is not there. Seeding one iteration per thread reaches each
    thread's own state.
    """
    for index in prange(get_num_threads()):
        np.random.seed(base_seed + index)


@njit(cache=True, parallel=True)
def _trace(n_ph: int, kq, t_step: float, tau0: float):
    dts = np.zeros(n_ph, dtype=np.float64)
    phs = np.zeros(n_ph, dtype=np.uint8)
    n_frames = kq.shape[0]
    if n_frames == 0:
        return dts, phs
    shift_nbrs = np.random.randint(0, n_frames if n_frames > 1 else 1, n_ph)

    for i in prange(n_ph):
        dts[i], phs[i] = _photon_walk(kq, t_step, tau0, shift_nbrs[i], n_frames)
    return dts, phs


# Seeded counterpart of `_trace`. It has to stay sequential: under
# `parallel=True` numba keeps one RNG state per thread and hands iterations to
# whichever thread is free, so the draws a photon sees would depend on the
# scheduler and the seed would not pin the result down.
@njit(cache=True, nogil=True)
def _trace_seeded(n_ph: int, kq, t_step: float, tau0: float, random_seed: int):
    np.random.seed(random_seed)
    dts = np.zeros(n_ph, dtype=np.float64)
    phs = np.zeros(n_ph, dtype=np.uint8)
    n_frames = kq.shape[0]
    if n_frames == 0:
        return dts, phs
    shift_nbrs = np.random.randint(0, n_frames if n_frames > 1 else 1, n_ph)

    for i in range(n_ph):
        dts[i], phs[i] = _photon_walk(kq, t_step, tau0, shift_nbrs[i], n_frames)
    return dts, phs


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
    kq = np.asarray(k_quench, dtype=np.float32)
    if random_seed is None:
        # Draw the thread seeds from the OS so consecutive unseeded runs are
        # genuinely independent samples, not the same one four times.
        _seed_worker_threads(int.from_bytes(os.urandom(4), "little"))
        return _trace(int(n_ph), kq, float(t_step), float(tau0))
    return _trace_seeded(int(n_ph), kq, float(t_step), float(tau0), int(random_seed))


# Curves are processed in `n_blocks` contiguous blocks, each accumulating into
# its own row of `partial`, which is then reduced in block order.
#
# **This is not an optimisation, it is the correctness fix.** The pre-move
# kernel accumulated straight into `decay[bin_idx]` from inside the `prange`.
# numba can privatise an array reduction only when it can prove the index is
# loop-invariant, and `bin_idx` is data-dependent -- so that was a genuine data
# race with lost updates, not merely a non-deterministic summation order.
# Measured on 8 threads (constant rate, 200 curves): sums scattered over
# 196.73 - 200.61 against a stable single-thread 200.612864255, i.e. up to **2%
# of the emitted intensity silently dropped**, varying run to run. Blocking also
# makes the reduction order fixed, so the result is now reproducible for a given
# `n_blocks`.
@njit(cache=True, parallel=True)
def _decay(n_curves: int, decay, dt_tac: float, k_quench, t_step: float,
           tau0: float, n_blocks: int):
    n_frames = k_quench.shape[0]
    n_bins = decay.shape[0]
    if n_frames == 0 or n_bins == 0 or n_curves <= 0:
        return

    intrinsic_rate = 1.0 / tau0 if tau0 > 0.0 else 0.0
    max_time = dt_tac * n_bins

    shifts = np.random.randint(0, n_frames // 2 if n_frames > 1 else 1, n_curves)
    partial = np.zeros((n_blocks, n_bins), dtype=np.float64)

    for block in prange(n_blocks):
        lo = block * n_curves // n_blocks
        hi = (block + 1) * n_curves // n_blocks
        for i in range(lo, hi):
            shift = shifts[i]
            intensity = 1.0
            t = 0.0
            frame = 0

            while t < max_time and intensity > 1e-6 and frame < n_frames:
                idx = (shift + frame) % n_frames
                total_rate = intrinsic_rate + float(k_quench[idx])
                dt = t_step
                if t + dt > max_time:
                    dt = max_time - t
                if dt <= 0.0:
                    break

                emitted = intensity * total_rate * dt
                bin_idx = int(t / dt_tac)
                if 0 <= bin_idx < n_bins:
                    partial[block, bin_idx] += emitted

                intensity *= np.exp(-total_rate * dt)
                t += dt
                frame += 1

    for block in range(n_blocks):
        for bin_idx in range(n_bins):
            decay[bin_idx] += partial[block, bin_idx]


def simulate_quenched_decay(
    n_curves: int,
    decay: np.ndarray,
    dt_tac: float,
    k_quench,
    t_step: float,
    tau0: float,
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
    """
    # `np.asarray(decay, dtype=np.float64)` *copies* whenever `decay` is not
    # already a float64 C-contiguous array, and the kernel would then fill the
    # copy -- an in-place function silently doing nothing. Fill a buffer we own
    # and write it back, so the contract holds for any array the caller passes.
    buffer = np.ascontiguousarray(decay, dtype=np.float64)
    n_curves = int(n_curves)
    # One block per thread, but never more blocks than curves -- an empty block
    # is harmless, yet `n_blocks > n_curves` would make most of them empty and
    # the reduction pointlessly wide.
    n_blocks = max(1, min(int(get_num_threads()), n_curves))
    _decay(
        n_curves,
        buffer,
        float(dt_tac),
        np.asarray(k_quench, dtype=np.float32),
        float(t_step),
        float(tau0),
        n_blocks,
    )
    if buffer is not decay:
        decay[...] = buffer
