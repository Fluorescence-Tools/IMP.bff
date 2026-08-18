"""The Lennard-Jones term, and the CHARMM36 parameters it reads.

Sterics: how hard two atoms push each other apart. Nothing here is about
fluorescence, which is why it is a tool rather than a domain -- and why it sits
here rather than inside either of the two things that need it. A rotamer's score
against a protein is a clash energy
(:mod:`IMP.bff.representation.rotamer.scoring`); so is an explicit dye's
internal energy (:mod:`IMP.bff.cgdye.topology.dye`). Those are a core
representation and a legacy package respectively, and the first must not depend
on the second to compute a steric term.

Lifted out of ``cgdye/topology/dye.py`` in the 2026-08-18 cleanup, when
``rotamer/`` moved into ``representation/`` and left that edge pointing the
wrong way.
"""

from __future__ import annotations

import numpy as np

__all__ = ["CHARMM36_LJ", "lj_params", "lj_cross", "lj_parameter_arrays",
           "lj_energy"]

CHARMM36_LJ = {
    "C": {"rmin_half": 2.02446316, "epsilon": -0.06394724},
    "N": {"rmin_half": 1.89285714, "epsilon": -0.15428571},
    "O": {"rmin_half": 1.693, "epsilon": -0.12642017},
    "S": {"rmin_half": 2.1, "epsilon": -0.47},
    "H": {"rmin_half": 0.98357778, "epsilon": -0.03466645},
}


def lj_params(element):
    return CHARMM36_LJ.get(element, CHARMM36_LJ["C"])


def lj_cross(elem_i, elem_j):
    """Lorentz-Berthelot cross parameters ``(rmin, eps)`` for two elements.

    ``rmin = rmin_half_i + rmin_half_j`` (Angstrom), ``eps = sqrt(eps_i eps_j)``
    (kcal/mol, positive well depth).
    """
    pi = lj_params(elem_i)
    pj = lj_params(elem_j)
    rmin = pi["rmin_half"] + pj["rmin_half"]
    eps = (pi["epsilon"] * pj["epsilon"]) ** 0.5
    return rmin, eps


def lj_parameter_arrays(elements):
    """``(rmin_half, epsilon)`` arrays for a sequence of element symbols.

    Unknown elements fall back to carbon, as ``lj_params`` does.
    """
    params = [lj_params(e) for e in elements]
    rmin_half = np.array([p["rmin_half"] for p in params], dtype=np.float64)
    epsilon = np.array([p["epsilon"] for p in params], dtype=np.float64)
    return rmin_half, epsilon


def lj_energy(r, rmin, eps, *, repulsive_only=False, cutoff=None, r_floor=0.01):
    """12-6 Lennard-Jones energy ``eps * ((rmin/r)^12 - 2 (rmin/r)^6)``.

    Vectorised over numpy arrays (``r``, ``rmin``, ``eps`` broadcast against
    each other); scalars work too. ``r`` is clamped to ``r_floor`` to keep the
    energy finite. ``repulsive_only=True`` zeroes the attractive tail
    (``r >= rmin``), ``cutoff`` zeroes pairs beyond that distance. With
    positive ``eps`` (the ``lj_cross`` convention) the well depth at ``rmin``
    is ``-eps``.
    """
    r_arr = np.asarray(r, dtype=np.float64)
    r_safe = np.maximum(r_arr, r_floor)
    ratio6 = np.power(np.asarray(rmin, dtype=np.float64) / r_safe, 6)
    energy = np.asarray(eps, dtype=np.float64) * (ratio6 * ratio6 - 2.0 * ratio6)
    if repulsive_only:
        energy = np.where(r_arr < rmin, energy, 0.0)
    if cutoff is not None:
        energy = np.where(r_arr < cutoff, energy, 0.0)
    if np.ndim(energy) == 0:
        return float(energy)
    return energy
