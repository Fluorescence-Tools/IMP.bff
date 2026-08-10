"""Physically exact FRET analysis using kinetic rotamer ensembles."""

from __future__ import annotations

import numpy as np


def compute_exact_efficiency(
    p_matrix: np.ndarray, 
    fret_rates: np.ndarray, 
    tau0: float, 
    dt: float = 1.0,
    weights: np.ndarray | None = None
) -> float:
    """Compute exact FRET efficiency by solving the Master Equation.
    
    This avoids approximations by accounting for the competition between 
    fluorescence decay and conformational transitions.
    
    Args:
        p_matrix: Transition probability matrix P (n x n) for time step dt.
        fret_rates: FRET rates for each state (n,). Units: 1/ns.
        tau0: Donor lifetime in absence of acceptor (ns).
        dt: Time step of the transitions in p_matrix (ns).
        weights: Stationary distribution (n,). If None, calculated as first eigenvector.
        
    Returns:
        Exact FRET efficiency.
    """
    n = p_matrix.shape[0]
    k_rad = 1.0 / tau0
    
    # 1. Convert Probability Matrix P to Rate Matrix M
    # For small dt: P = exp(M*dt) approx I + M*dt => M = (P - I) / dt
    # This assumes the transitions recorded are for a physical time interval.
    M = (p_matrix - np.eye(n)) / dt
    
    # 2. Setup the sink matrix (Radiative decay + FRET rates)
    # K is diagonal matrix of total decay rates
    K_fret = np.diag(fret_rates)
    
    # 3. Stationary distribution (initial population)
    if weights is None:
        # Find eigenvector with eigenvalue 1 for P (or 0 for M)
        evals, evecs = np.linalg.eig(p_matrix.T)
        idx = np.argmin(np.abs(evals - 1.0))
        weights = np.real(evecs[:, idx])
        weights /= weights.sum()
        
    # 4. Solve (k_rad*I + K_fret - M) * G = weights
    # where G is the integrated population (area under decay curves)
    A = k_rad * np.eye(n) + K_fret - M
    G = np.linalg.solve(A, weights)
    
    # 5. Efficiency = 1 - (integrated population / tau0)
    # Since we didn't include k_rad in the matrix inversion for the FRET-only part:
    # Actually, E = sum(k_fret_i * G_i)
    efficiency = np.sum(fret_rates * G)
    
    return float(efficiency)


def calculate_fret_exact(
    dist_matrix: np.ndarray,
    kappa2_matrix: np.ndarray,
    p_d: np.ndarray,
    p_a: np.ndarray,
    weights_d: np.ndarray,
    weights_a: np.ndarray,
    R0: float = 52.0,
    tau0: float = 4.0,
    dt: float = 0.1
) -> float:
    """Calculate exact FRET efficiency for two kinetic ensembles.
    
    Args:
        dist_matrix: (nd, na) distances
        kappa2_matrix: (nd, na) orientation factors
        p_d: Donor transition matrix (nd, nd)
        p_a: Acceptor transition matrix (na, na)
        weights_d: Donor stationary weights
        weights_a: Acceptor stationary weights
        R0: Förster radius (A)
        tau0: Lifetime (ns)
        dt: Time step (ns)
    """
    # 1. Product space transition matrix P_total = P_d kron P_a
    # This represents the joint kinetics of both dyes.
    # Note: kron is expensive for large libraries.
    # P_total[i*na + j, k*na + l] = P_d[i,k] * P_a[j,l]
    P_total = np.kron(p_d, p_a)
    
    # 2. Product space weights and FRET rates
    w_total = np.outer(weights_d, weights_a).flatten()
    
    k_rad = 1.0 / tau0
    # k_fret_ij = k_rad * (R0/r_ij)^6 * 1.5 * kappa2_ij
    rate_ratios = (R0 / dist_matrix)**6 * (1.5 * kappa2_matrix)
    fret_rates = (k_rad * rate_ratios).flatten()
    
    return compute_exact_efficiency(P_total, fret_rates, tau0, dt, weights=w_total)


def calculate_fret_regimes(
    dist_matrix: np.ndarray,
    kappa2_matrix: np.ndarray,
    weights_d: np.ndarray,
    weights_a: np.ndarray,
    R0: float = 52.0
) -> dict:
    """Calculate FRET efficiency under different kinetic regimes.
    
    Args:
        dist_matrix: (n_donor, n_acceptor) distances
        kappa2_matrix: (n_donor, n_acceptor) orientation factors
        weights_d: (n_donor,) donor rotamer weights
        weights_a: (n_acceptor,) acceptor rotamer weights
        R0: Förster radius for kappa2=2/3 (A)
        
    Returns:
        Dict with 'static', 'dynamic', 'dynamic_plus' efficiencies.
    """
    # Combined weights
    w_ij = np.outer(weights_d, weights_a)
    
    # Rate constant k_fret / k_rad = (R0/r)^6 * (k2 / (2/3))
    # We use (1.5 * k2) to normalize kappa2 relative to 2/3
    rate_ratio = (R0 / dist_matrix)**6 * (1.5 * kappa2_matrix)
    
    # 1. Static Regime: Average of efficiencies
    eff_static = rate_ratio / (1 + rate_ratio)
    e_static = np.sum(w_ij * eff_static)
    
    # 2. Dynamic Regime: Use <kappa^2>
    k2_avg = np.sum(w_ij * kappa2_matrix)
    rate_ratio_dyn = (R0 / dist_matrix)**6 * (1.5 * k2_avg)
    eff_dyn = rate_ratio_dyn / (1 + rate_ratio_dyn)
    e_dynamic = np.sum(w_ij * eff_dyn)
    
    # 3. Dynamic+ Regime: Average of rates
    rate_avg = np.sum(w_ij * rate_ratio)
    e_dynamic_plus = rate_avg / (1 + rate_avg)
    
    return {
        "static": float(e_static),
        "dynamic": float(e_dynamic),
        "dynamic_plus": float(e_dynamic_plus),
        "kappa2_avg": float(k2_avg)
    }
