"""Mean-field rotamer probability update for dye rotamer libraries.

Algorithmic source: IMP RotamerCalculator.cpp, RotamerCalculator::transform(),
equations (39)–(42) from Bauer et al. / Dunbrack mean-field algorithm.

The IMP C++ implementation operates on amino-acid side chains using a hardcoded
topology table (chi-axis atoms: N-CA-CB-CG etc.).  For the dye linker — a
non-canonical flexible chain — we cannot reuse that table, so we reimplement
the *algorithm* (equations 39–42) in Python/NumPy against our own coordinate
arrays and LJ scoring backend.

What this module does
---------------------
Given:
  - n_clusters rotamer cluster coordinate sets for the dye  (n_clusters, n_atoms, 3)
  - initial Boltzmann weights  q[i]  for each cluster  (n_clusters,)
  - static protein coordinates                          (n_prot_atoms, 3)
  - LJ pair parameters between dye atoms and protein atoms

It iterates (Eq. 41/42):

  E_P[i] = E_bb[i] + Σ_{j≠i} q[j] * E_SC[ij]          (Eq. 41)
  q[i]   ∝ P_boltz[i] * exp(-K * E_P[i])               (Eq. 42)

until convergence or max_iter, returning updated cluster weights.

For a *single* dye (no dye–dye cross-terms) this reduces to a pure
protein–dye interaction reweighting:

  E_P[i] = E_bb[i]  (Eq. 39 only)
  q[i]   ∝ P_boltz[i] * exp(-K * E_bb[i])

For *multi-dye* systems the E_SC cross term between the two dyes' rotamers
is computed as the pairwise LJ interaction between rotamer i of dye A and
rotamer j of dye B, then averaged with the current q of the partner dye.
"""

from __future__ import annotations

import numpy as np

from IMP.bff.cgdye.topology.dye import lj_cross, lj_energy


# ---------------------------------------------------------------------------
# LJ helpers (mirror IMP's PairScore, but in NumPy)
# ---------------------------------------------------------------------------

def _lj_energy_pairs(
    coords_a: np.ndarray,   # (n_a, 3)  dye atoms in one rotamer
    coords_b: np.ndarray,   # (n_b, 3)  reference atoms (protein or other dye)
    rmin: np.ndarray,       # (n_a * n_b,) or None for cut-off only
    eps: np.ndarray,        # (n_a * n_b,)
    r_cutoff: float = 12.0,
) -> float:
    """Compute total LJ interaction energy between two atom sets.

    Vectorized over all n_a × n_b pairs.
    """
    # (n_a, n_b, 3)
    diff = coords_a[:, None, :] - coords_b[None, :, :]
    r = np.linalg.norm(diff, axis=-1)           # (n_a, n_b)
    r_flat = r.ravel()

    # rmin / eps broadcast: if None assume uniform zeros → no interaction
    if rmin is None or eps is None:
        return 0.0

    # Repulsive-only LJ within the cutoff, through the shared kernel
    lj = lj_energy(r_flat, rmin, eps, repulsive_only=True, cutoff=r_cutoff)
    return float(lj.sum())


def _build_cross_lj_params(
    dye_elements: list[str],
    ref_elements: list[str],
) -> tuple[np.ndarray, np.ndarray]:
    """Build flattened LJ rmin / eps arrays for all dye × ref atom pairs.

    Returns (rmin, eps) arrays of shape (n_dye * n_ref,).
    """
    n_d = len(dye_elements)
    n_r = len(ref_elements)
    rmin_arr = np.zeros(n_d * n_r, dtype=np.float64)
    eps_arr  = np.zeros(n_d * n_r, dtype=np.float64)

    k = 0
    for ea in dye_elements:
        for eb in ref_elements:
            r, e = lj_cross(ea, eb)
            rmin_arr[k] = r
            eps_arr[k]  = e
            k += 1
    return rmin_arr, eps_arr


# ---------------------------------------------------------------------------
# Trick 1 (AABB) used here too — import from scoring
# ---------------------------------------------------------------------------

def _aabb_overlap(
    dye_box: np.ndarray,   # (6,)  dye rotamer box
    ref_box: np.ndarray,   # (6,)  protein box
) -> bool:
    """Return True if the two AABBs overlap.  Port of IMP::intersect()."""
    return (
        dye_box[3] >= ref_box[0] and dye_box[0] <= ref_box[3] and
        dye_box[4] >= ref_box[1] and dye_box[1] <= ref_box[4] and
        dye_box[5] >= ref_box[2] and dye_box[2] <= ref_box[5]
    )


# ---------------------------------------------------------------------------
# Main public function — single-dye mean-field reweighting (Eq. 39 + 42)
# ---------------------------------------------------------------------------

def rotamer_mean_field_weights(
    rotamer_coords: np.ndarray,         # (n_clusters, n_dye_atoms, 3)
    initial_weights: np.ndarray,        # (n_clusters,)  Boltzmann prior
    protein_coords: np.ndarray,         # (n_prot_atoms, 3)
    dye_elements: list[str],            # element per dye atom, len = n_dye_atoms
    protein_elements: list[str],        # element per protein atom, len = n_prot_atoms
    K: float = 1.0,
    n_iter: int = 10,
    aabb_pad: float = 3.5,
    r_cutoff: float = 12.0,
) -> np.ndarray:
    """Iterative mean-field update of dye rotamer cluster weights.

    Implements equations (39) and (42) from IMP RotamerCalculator.cpp for a
    *single* dye (no dye–dye cross-terms).  The protein acts as the fixed
    backbone environment (E_bb), and q[i] is updated each iteration.

    For a single dye E_P[i] = E_bb[i] (the Σ q_lm E_SC term vanishes because
    there is only one dye).  The loop still converges in one iteration but is
    kept general for the multi-dye extension.

    Args:
        rotamer_coords:   Cluster centre coordinates, shape (n_clusters, n_dye_atoms, 3).
        initial_weights:  Initial Boltzmann weights, shape (n_clusters,), sum ≈ 1.
        protein_coords:   Static protein atom positions, shape (n_prot_atoms, 3).
        dye_elements:     Element symbol for each dye atom.
        protein_elements: Element symbol for each protein atom.
        K:                Boltzmann scaling factor (IMP's ``K`` in Eq. 42).
                          K=1 ≈ kT-scale energies; increase to sharpen distribution.
        n_iter:           Number of mean-field iterations (default 10, usually
                          converges in 3–5).
        aabb_pad:         AABB padding in Å (same 3.5 Å default as IMP).
        r_cutoff:         Atom-pair distance cut-off in Å.

    Returns:
        weights: shape (n_clusters,), updated and renormalized.
    """
    n_clusters = rotamer_coords.shape[0]
    q = initial_weights.copy().astype(np.float64)

    # Pre-build LJ cross params for dye × protein (Eq. 39 — E_bb)
    rmin_dp, eps_dp = _build_cross_lj_params(dye_elements, protein_elements)

    # AABB for protein (static reference) — Trick 1 pre-filter
    prot_box = np.concatenate([
        protein_coords.min(axis=0) - aabb_pad,
        protein_coords.max(axis=0) + aabb_pad,
    ])

    # Pre-compute bounding boxes per rotamer cluster — Trick 1
    rot_boxes = np.concatenate([
        rotamer_coords.min(axis=1) - aabb_pad,    # (n_clusters, 3)
        rotamer_coords.max(axis=1) + aabb_pad,    # (n_clusters, 3)
    ], axis=1)                                     # (n_clusters, 6)

    # Pre-compute E_bb[i] = protein–dye interaction energy for each cluster
    # (Eq. 39).  This is constant across iterations for a static protein.
    E_bb = np.zeros(n_clusters, dtype=np.float64)
    for i in range(n_clusters):
        # Trick 1: skip if bounding boxes don't overlap
        if not _aabb_overlap(rot_boxes[i], prot_box):
            continue
        E_bb[i] = _lj_energy_pairs(
            rotamer_coords[i],     # dye rotamer i
            protein_coords,
            rmin_dp, eps_dp,
            r_cutoff=r_cutoff,
        )

    # Iterative mean-field update — Eq. 41 / 42
    # For single dye: E_P[i] = E_bb[i]  (no SC-SC cross term)
    for _iteration in range(n_iter):
        E_P = E_bb.copy()  # Eq. 41 simplified for single dye

        # Eq. 42: q[i] ∝ P_boltz[i] * exp(-K * E_P[i])
        log_q = np.log(np.maximum(q, 1e-300)) - K * E_P
        log_q -= log_q.max()               # numerical stability (log-sum-exp)
        q_new = np.exp(log_q)
        denom = q_new.sum()
        if denom > 0:
            q = q_new / denom
        else:
            # Degenerate: fall back to uniform
            q = np.ones(n_clusters, dtype=np.float64) / n_clusters

    return q


# ---------------------------------------------------------------------------
# Multi-dye extension — Eq. 40/41/42 with dye–dye SC cross term
# ---------------------------------------------------------------------------

def rotamer_mean_field_weights_multi_dye(
    rotamer_coords_list: list[np.ndarray],    # list of (n_i, n_atoms_i, 3)
    initial_weights_list: list[np.ndarray],   # list of (n_i,)
    protein_coords: np.ndarray,               # (n_prot_atoms, 3)
    dye_elements_list: list[list[str]],
    protein_elements: list[str],
    K: float = 1.0,
    n_iter: int = 10,
    aabb_pad: float = 3.5,
    r_cutoff: float = 12.0,
) -> list[np.ndarray]:
    """Multi-dye mean-field weight update including dye–dye cross-interactions.

    Implements equations (39)–(42) from IMP RotamerCalculator.cpp for
    multiple dyes whose rotamer clouds may interact with each other.

    Args:
        rotamer_coords_list:  One array per dye, shape (n_i, n_atoms_i, 3).
        initial_weights_list: Initial weights per dye, each shape (n_i,).
        protein_coords:       Static protein atoms, shape (n_prot_atoms, 3).
        dye_elements_list:    Element lists per dye.
        protein_elements:     Element list for protein.
        K:                    Boltzmann scaling (Eq. 42).
        n_iter:               Mean-field iterations.
        aabb_pad:             AABB padding in Å.
        r_cutoff:             Atom-pair distance cut-off in Å.

    Returns:
        Updated weight arrays, one per dye.
    """
    n_dyes = len(rotamer_coords_list)
    # q[d] holds current weights for dye d
    q_list = [w.copy().astype(np.float64) for w in initial_weights_list]

    # AABB for protein
    prot_box = np.concatenate([
        protein_coords.min(axis=0) - aabb_pad,
        protein_coords.max(axis=0) + aabb_pad,
    ])

    # Pre-compute E_bb[d][i] — protein interaction per dye d, cluster i  (Eq. 39)
    E_bb_all = []
    rmin_dp_all = []
    eps_dp_all  = []
    rot_boxes_all = []

    for d in range(n_dyes):
        coords_d = rotamer_coords_list[d]
        n_clusters_d = coords_d.shape[0]
        rmin_dp, eps_dp = _build_cross_lj_params(
            dye_elements_list[d], protein_elements
        )
        rmin_dp_all.append(rmin_dp)
        eps_dp_all.append(eps_dp)

        rot_boxes_d = np.concatenate([
            coords_d.min(axis=1) - aabb_pad,
            coords_d.max(axis=1) + aabb_pad,
        ], axis=1)
        rot_boxes_all.append(rot_boxes_d)

        E_bb_d = np.zeros(n_clusters_d, dtype=np.float64)
        for i in range(n_clusters_d):
            if not _aabb_overlap(rot_boxes_d[i], prot_box):
                continue
            E_bb_d[i] = _lj_energy_pairs(
                coords_d[i], protein_coords,
                rmin_dp, eps_dp,
                r_cutoff=r_cutoff,
            )
        E_bb_all.append(E_bb_d)

    # Pre-compute cross LJ params for each dye pair
    cross_rmin = {}
    cross_eps  = {}
    for d1 in range(n_dyes):
        for d2 in range(d1 + 1, n_dyes):
            rp, ep = _build_cross_lj_params(
                dye_elements_list[d1], dye_elements_list[d2]
            )
            cross_rmin[(d1, d2)] = rp
            cross_eps[(d1, d2)]  = ep

    # Mean-field iteration — Eq. 41 / 42
    for _iteration in range(n_iter):
        E_P_all = [E_bb_all[d].copy() for d in range(n_dyes)]

        # Add dye–dye SC cross terms (Eq. 40 / 41)
        for d1 in range(n_dyes):
            n_i = rotamer_coords_list[d1].shape[0]
            for d2 in range(n_dyes):
                if d1 == d2:
                    continue
                d_lo, d_hi = (d1, d2) if d1 < d2 else (d2, d1)
                rmin_cc = cross_rmin[(d_lo, d_hi)]
                eps_cc  = cross_eps[(d_lo, d_hi)]

                n_j = rotamer_coords_list[d2].shape[0]
                for i in range(n_i):
                    for j in range(n_j):
                        # AABB pre-filter — Trick 1
                        if not _aabb_overlap(
                            rot_boxes_all[d1][i],
                            rot_boxes_all[d2][j],
                        ):
                            continue
                        E_sc_ij = _lj_energy_pairs(
                            rotamer_coords_list[d1][i],
                            rotamer_coords_list[d2][j],
                            rmin_cc, eps_cc,
                            r_cutoff=r_cutoff,
                        )
                        # Eq. 41: weighted by partner probability
                        E_P_all[d1][i] += q_list[d2][j] * E_sc_ij

        # Eq. 42: renormalize each dye independently
        for d in range(n_dyes):
            log_q = np.log(np.maximum(q_list[d], 1e-300)) - K * E_P_all[d]
            log_q -= log_q.max()
            q_new = np.exp(log_q)
            denom = q_new.sum()
            q_list[d] = q_new / denom if denom > 0 else np.ones_like(q_new) / len(q_new)

    return q_list
