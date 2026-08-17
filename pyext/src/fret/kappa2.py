"""Orientation factor κ² of donor/acceptor transition dipoles.

Lives in ``fret`` (label-pair physics, no coordinate model attached); cgdye's
rotamer scoring and the ensemble pair kernels import it. Resolves the PRD-47
κ² placement on the bff side: κ² lives in imp.bff.
"""

from __future__ import annotations

import numpy as np


def kappa2_from_dipoles(mu_donor: np.ndarray, mu_acceptor: np.ndarray, r_vectors: np.ndarray) -> np.ndarray:
    """Compute orientation factors from dipole and distance vectors.

    Parameters
    ----------
    mu_donor : numpy.ndarray
        Donor transition-dipole vectors with shape ``(n_donor, 3)``.
    mu_acceptor : numpy.ndarray
        Acceptor transition-dipole vectors with shape ``(n_acceptor, 3)``.
    r_vectors : numpy.ndarray
        Donor-to-acceptor vectors with shape ``(n_donor, n_acceptor, 3)``.

    Returns
    -------
    numpy.ndarray
        ``kappa^2`` matrix.
    """
    r_vectors = np.asarray(r_vectors, dtype=np.float64)
    r_norm = np.linalg.norm(r_vectors, axis=2, keepdims=True)
    r_unit = np.divide(r_vectors, r_norm, out=np.zeros_like(r_vectors), where=r_norm > 0.0)
    cos_da = np.einsum("ik,jk->ij", mu_donor, mu_acceptor)
    cos_dr = np.einsum("ik,ijk->ij", mu_donor, r_unit)
    cos_ar = np.einsum("jk,ijk->ij", mu_acceptor, r_unit)
    return np.power(cos_da - 3.0 * cos_dr * cos_ar, 2)


def kappa2_isotropic() -> float:
    """⟨κ²⟩ for isotropically and rapidly reorienting dipoles (2/3)."""
    return 2.0 / 3.0


__all__ = ["kappa2_from_dipoles", "kappa2_isotropic"]
