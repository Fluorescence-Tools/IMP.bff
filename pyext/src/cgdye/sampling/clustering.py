"""Conformational clustering for rotamer library generation."""

from __future__ import annotations

import numpy as np


def rmsd_no_align(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """Compute RMSD between two sets of coordinates without alignment.
    
    Assumes atoms are in the same order.
    """
    diff = coords1 - coords2
    return np.sqrt(np.mean(np.sum(diff**2, axis=-1)))


def cluster_frames_leader(coords: np.ndarray, threshold: float) -> list[int]:
    """Greedy leader clustering algorithm.
    
    Returns a list of cluster center indices.
    
    Args:
        coords: Array of shape (n_frames, n_atoms, 3)
        threshold: RMSD threshold for a new cluster (in Angstroms)
    """
    if coords.ndim != 3:
        raise ValueError("coords must be (n_frames, n_atoms, 3)")
    
    n_frames = coords.shape[0]
    if n_frames == 0:
        return []
    
    centers = [0]
    for i in range(1, n_frames):
        # Compute RMSD to all existing centers
        is_new = True
        for c_idx in centers:
            rmsd = rmsd_no_align(coords[i], coords[c_idx])
            if rmsd < threshold:
                is_new = False
                break
        if is_new:
            centers.append(i)
            
    return centers


def assign_frames_to_clusters(coords: np.ndarray, centers: list[int]) -> np.ndarray:
    """Assign each frame to the nearest cluster center.
    
    Returns array of shape (n_frames,) containing center indices.
    """
    n_frames = coords.shape[0]
    assignments = np.zeros(n_frames, dtype=int)
    center_coords = coords[centers]
    
    for i in range(n_frames):
        # Vectorized RMSD to all centers
        diffs = center_coords - coords[i]
        rmsds = np.sqrt(np.mean(np.sum(diffs**2, axis=-1), axis=-1))
        assignments[i] = np.argmin(rmsds)
        
    return assignments
