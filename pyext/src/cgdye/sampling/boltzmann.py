"""Boltzmann scoring and statistical weight calculation."""

from __future__ import annotations

import numpy as np


def boltzmann_weights(energies: np.ndarray, temperature: float = 298.15) -> np.ndarray:
    """Compute normalized Boltzmann weights from energies.
    
    Args:
        energies: Array of energies in kcal/mol (or consistent units with kT)
        temperature: Temperature in Kelvin
        
    Returns:
        Normalized weights summing to 1.0.
    """
    # kB in kcal/(mol*K)
    KB = 0.0019872041
    kt = KB * temperature
    
    # Use log-sum-exp trick for numerical stability
    e_min = np.min(energies)
    shifted_energies = energies - e_min
    unnormalized_weights = np.exp(-shifted_energies / kt)
    
    return unnormalized_weights / np.sum(unnormalized_weights)


def rotamer_cluster_weights(assignments: np.ndarray, frame_weights: np.ndarray, n_clusters: int) -> np.ndarray:
    """Aggregate frame weights into cluster weights.
    
    Args:
        assignments: Cluster index for each frame (n_frames,)
        frame_weights: Normalized weight for each frame (n_frames,)
        n_clusters: Total number of clusters
        
    Returns:
        Array of shape (n_clusters,) with normalized cluster weights.
    """
    weights = np.zeros(n_clusters)
    for i in range(len(assignments)):
        weights[assignments[i]] += frame_weights[i]
    
    return weights / np.sum(weights)
