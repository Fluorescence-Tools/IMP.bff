"""Utilities for kinetic trajectory reconstruction from rotamer libraries."""

from __future__ import annotations

import random
import numpy as np


def rotamer_transition_matrix(transition_counts: list[list[int]]) -> np.ndarray:
    """Convert raw transition counts to probabilities.
    
    Returns matrix P where P[i,j] is the probability of jumping from state i to j.
    Rows sum to 1.
    """
    counts = np.array(transition_counts, dtype=float)
    row_sums = counts.sum(axis=1)
    
    # Handle states with no outgoing transitions (shouldn't happen in long walks)
    # by making them sink states (jump to self)
    for i in range(len(row_sums)):
        if row_sums[i] == 0:
            counts[i, i] = 1.0
            row_sums[i] = 1.0
            
    return counts / row_sums[:, np.newaxis]


def rotamer_correlation_times(transition_counts: list[list[int]], timestep: float) -> np.ndarray:
    """Calculate relaxation times from the transition matrix.
    
    Args:
        transition_counts: Matrix of transition counts.
        timestep: Time between frames in the sampling walk.
        
    Returns:
        Array of relaxation times.
    """
    p = rotamer_transition_matrix(transition_counts)
    vals, _ = np.linalg.eig(p.T)
    # Sort eigenvalues by magnitude
    vals = np.sort(np.abs(vals))[::-1]
    
    # Relaxation times: t_i = -timestep / ln(lambda_i)
    # lambda_1 should be 1.0 (stationary state)
    times = []
    for v in vals[1:]:
        if v > 1e-10 and v < 0.99999999:
            times.append(-timestep / np.log(v))
        else:
            times.append(0.0)
    return np.array(times)


def rotamer_rotational_correlation_time(transition_counts: list[list[int]], timestep: float) -> float:
    """Estimate the slowest rotational correlation time.
    
    Args:
        transition_counts: Matrix of transition counts.
        timestep: Time between frames in the sampling walk.
        
    Returns:
        The slowest correlation time.
    """
    times = rotamer_correlation_times(transition_counts, timestep)
    if len(times) > 0:
        return float(np.max(times))
    return 0.0


def reconstruct_rotamer_trajectory(
    lib: dict, 
    n_frames: int, 
    start_index: int | None = None,
    seed: int | None = None
) -> list[int]:
    """Generate a sequence of rotamer indices using transition probabilities.
    
    Args:
        lib: Rotamer library dict with 'transitions'
        n_frames: Length of desired trajectory
        start_index: Initial rotamer index (0-based)
        
    Returns:
        List of rotamer indices (0-based)
    """
    if "transitions" not in lib or lib["transitions"] is None:
        raise ValueError("Library does not contain transition data")
        
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
        
    p_matrix = rotamer_transition_matrix(lib["transitions"])
    n_states = p_matrix.shape[0]
    
    weights = np.array(lib["weight"], dtype=float)
    weights /= weights.sum() # Ensure exact sum to 1.0
    
    if start_index is None:
        # Sample starting state from Boltzmann weights
        start_index = np.random.choice(n_states, p=weights)
        
    traj = [start_index]
    current = start_index
    
    for _ in range(n_frames - 1):
        # Sample next state based on current row of P
        probs = p_matrix[current]
        probs /= probs.sum() # Ensure exact sum to 1.0
        nxt = np.random.choice(n_states, p=probs)
        traj.append(nxt)
        current = nxt
        
    return traj
