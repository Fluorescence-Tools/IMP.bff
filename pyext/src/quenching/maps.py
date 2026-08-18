"""Rate and mobility maps on an accessible-volume grid.

The *field* formulation of the quenched-dye model: instead of walking a dye and
reading rates along its trajectory (:mod:`IMP.bff.dynamics.brownian`), assign
every voxel a diffusion coefficient and a decay rate, and propagate the excited
state on the grid (:mod:`IMP.bff.dynamics.smoluchowski`). Deterministic, no shot
noise, and the natural home for a *distance-dependent* quenching law.

Moved here from ChiSurf (``chisurf/core/structure/av/functions.py``) by PRD-109.

Two quenching laws, deliberately both kept:

:func:`quenching_rate_map`
    ``k(r) = 1/tau0 + sum_a kQ_a * exp(-(|r - r_a| - r_dye) / rC_a)`` -- a
    Dexter-like exponential in the dye-to-atom distance, per **atom**.
:func:`IMP.bff.quenching_rate_grid`
    A hard contact sphere per quencher: full rate inside ``quench_radius``,
    zero outside, rates adding where spheres overlap. This is QuEst's law.

They are different models, not two spellings of one. The exponential has no
sharp contact radius and reaches further at low amplitude; the step model is
what QuEst's published parameters were calibrated against.
"""

from __future__ import annotations

from collections import OrderedDict

import numpy as np

import IMP.bff


__all__ = [
    "grid_axis",
    "atomic_quenching_parameters",
    "diffusion_coefficient_map",
    "radial_diffusion_map",
    "slow_diffusion_near_atoms",
    "quenching_rate_map",
    "fret_rate_map",
]


def grid_axis(ng: int, dg: float) -> np.ndarray:
    """Voxel-centre offsets from the grid anchor, in Angstrom.

    ``(i - (ng - 1) // 2) * dg`` -- the same integer centre offset as
    :func:`IMP.bff.grid_center_index`, so a map built here indexes the way the
    Brownian walk samples it.
    """
    return (np.arange(int(ng), dtype=np.float64) - (int(ng) - 1) // 2) * float(dg)


def atomic_quenching_parameters(atoms, quencher):
    """Per-atom quenching rate ``kQ`` and characteristic distance ``rC``.

    :param atoms: structured array with ``res_name`` and ``atom_name`` fields.
    :param quencher: ``{res_name: {atom_name: (kQ, rC)}}``. Atoms not named in
        it do not quench.
    :returns: ``(kQ, rC)``, each ``(n_atoms,)`` float64, zero where absent.
    """
    n = len(atoms)
    kQ = np.zeros(n, dtype=np.float64)
    rC = np.zeros(n, dtype=np.float64)
    res_names = atoms["res_name"]
    atom_names = atoms["atom_name"]
    for i in range(n):
        by_atom = quencher.get(str(res_names[i]))
        if not by_atom:
            continue
        entry = by_atom.get(str(atom_names[i]))
        if entry is None:
            continue
        kQ[i], rC[i] = float(entry[0]), float(entry[1])
    return kQ, rC


# Every kernel below iterates the FULL grid, `range(ng)`.
#
# ChiSurf's originals ran `for ix in range(-npm, npm)` with `npm = (ng - 1) // 2`
# -- one short on every axis, so the outermost slab was never written. In
# `assign_diffusion_to_grid_1` the output was `np.empty_like(...)`, so that slab
# was **uninitialised memory** read back as a diffusion coefficient, and that
# variant is the one `DynamicAV` used. In the rate maps the output was zeroed,
# so the outer face silently reported no quenching and no FRET. For an even
# `ng` -- which is what the IMP AV path produces -- it missed two slabs, not
# one. Fixed on the move (PRD-109).


def _slow_near_atoms(d_map, density, axis, r0, atoms_xyz, min_distance_sq, factor):
    """Slow the mobility wherever the dye contacts an atom. **C++.**

    The factor is applied *once per contacting atom*, so it compounds: a voxel
    touching 45 atoms at 0.985 ends at 0.5, and at 0.9 at 9e-3. Only values
    within a whisker of 1.0 are meaningful -- a property of the model, not a
    taste.
    """
    ng = int(density.shape[0])
    return np.asarray(IMP.bff.slow_near_atoms(
        np.ascontiguousarray(d_map, dtype=np.float64).ravel(),
        np.ascontiguousarray(density, dtype=np.float64).ravel(),
        np.ascontiguousarray(axis, dtype=np.float64).ravel(),
        np.ascontiguousarray(r0, dtype=np.float64).ravel(),
        np.ascontiguousarray(atoms_xyz, dtype=np.float64).ravel(),
        float(min_distance_sq), float(factor)),
        dtype=np.float64).reshape(ng, ng, ng)


def slow_diffusion_near_atoms(
    d_map, density, r0, dg, atoms_xyz, min_distance, slow_factor
) -> np.ndarray:
    """Scale a diffusion map down once per atom in contact with each voxel.

    The more atoms crowd a voxel the slower the dye moves there; the factor is
    applied once per contacting atom, so it compounds. Voxels outside the
    accessible volume are set to zero -- the dye cannot be there, and a non-zero
    coefficient would let the solver leak population into the protein.

    J. Chem. Phys. 126, 044707 (2007), eq. (3).

    :param d_map: base diffusion coefficient per voxel (A^2/ns).
    :param density: binary occupancy of the accessible volume.
    :param r0: grid anchor (the attachment point).
    :param dg: voxel edge in Angstrom.
    :param atoms_xyz: ``(n, 3)`` atom coordinates.
    :param min_distance: contact distance, roughly dye radius + 2 * vdW.
    :param slow_factor: factor applied per contacting atom, in [0, 1].
    """
    density = np.ascontiguousarray(density, dtype=np.float64)
    ng = density.shape[0]
    return _slow_near_atoms(
        np.ascontiguousarray(d_map, dtype=np.float64),
        density,
        grid_axis(ng, dg),
        np.asarray(r0, dtype=np.float64),
        np.ascontiguousarray(atoms_xyz, dtype=np.float64),
        float(min_distance) ** 2,
        float(slow_factor),
    )


def radial_diffusion_map(density, dg, f) -> np.ndarray:
    """A diffusion map from a radial profile ``f(distance from the anchor)``.

    Used to give the dye a mobility that depends on how far the linker is
    stretched, before :func:`slow_diffusion_near_atoms` adds the local crowding.
    Vectorised, so *f* must accept an array.
    """
    density = np.asarray(density)
    ng = density.shape[0]
    axis = grid_axis(ng, dg)
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    return np.asarray(f(np.sqrt(x * x + y * y + z * z)), dtype=np.float64)


def diffusion_coefficient_map(
    density, r0, dg, atoms_xyz, *, free_diffusion, min_distance, slow_factor,
    radial_profile=None,
) -> np.ndarray:
    """The mobility field: a base coefficient, slowed by nearby atoms.

    :param free_diffusion: the unhindered dye diffusion coefficient (A^2/ns),
        used where *radial_profile* is not given.
    :param radial_profile: optional callable of the distance from the anchor,
        replacing the constant base coefficient.
    """
    density = np.asarray(density)
    if radial_profile is None:
        base = np.full(density.shape, float(free_diffusion), dtype=np.float64)
    else:
        base = radial_diffusion_map(density, dg, radial_profile)
    return slow_diffusion_near_atoms(
        base, density, r0, dg, atoms_xyz, min_distance, slow_factor
    )


def _quenching_map(density, axis, r0, atoms_xyz, kQ, rC, dye_radius, inv_tau0):
    """PET rate field. **C++.**

    Voxels outside the accessible volume stay at **zero**, not at ``1/tau0``:
    the dye cannot be there, so it has no decay rate there.
    """
    ng = int(density.shape[0])
    return np.asarray(IMP.bff.quenching_map(
        np.ascontiguousarray(density, dtype=np.float64).ravel(),
        np.ascontiguousarray(axis, dtype=np.float64).ravel(),
        np.ascontiguousarray(r0, dtype=np.float64).ravel(),
        np.ascontiguousarray(atoms_xyz, dtype=np.float64).ravel(),
        np.ascontiguousarray(kQ, dtype=np.float64).ravel(),
        np.ascontiguousarray(rC, dtype=np.float64).ravel(),
        float(dye_radius), float(inv_tau0)),
        dtype=np.float64).reshape(ng, ng, ng)


def quenching_rate_map(
    density, r0, dg, atoms_xyz, kQ, rC, *, tau0, dye_radius
) -> np.ndarray:
    """Total decay rate per voxel: intrinsic plus exponential PET.

    ``k(r) = 1/tau0 + sum_a kQ_a * exp(-(|r - r_a| - r_dye) / rC_a)``

    The distance is measured from the dye **surface** (the atom distance less
    the dye radius), and the sum runs over every atom with a non-zero ``kQ`` and
    ``rC``. Unlike the contact-sphere law in
    :func:`IMP.bff.quenching_rate_grid` there is no cut-off: a distant atom
    contributes exponentially little rather than nothing.

    Voxels outside the accessible volume stay at zero -- **not** at ``1/tau0``.
    The solver reads this as the rate field on a domain masked by the same
    occupancy, so an unreachable voxel carries no rate at all.
    """
    density = np.ascontiguousarray(density, dtype=np.float64)
    ng = density.shape[0]
    return _quenching_map(
        density,
        grid_axis(ng, dg),
        np.asarray(r0, dtype=np.float64),
        np.ascontiguousarray(atoms_xyz, dtype=np.float64),
        np.ascontiguousarray(kQ, dtype=np.float64),
        np.ascontiguousarray(rC, dtype=np.float64),
        float(dye_radius),
        1.0 / float(tau0) if tau0 > 0.0 else 0.0,
    )


def _fret_map(density_d, density_a, axis_d, axis_a, r0_d, r0_a, r0_6, kf, step):
    """FRET rate field, donor volume against acceptor volume. **C++.**

    **The harmonic mean**: the acceptor-weighted mean transfer *time* is
    accumulated and inverted, which is the static limit -- the acceptor does not
    move within the donor's excited-state lifetime. ``fret_rate_trace``
    accumulates the arithmetic mean of *rates* along a trajectory, which is the
    fast-exchange limit. Different physics, not variants.
    """
    ng = int(density_d.shape[0])
    return np.asarray(IMP.bff.fret_map(
        np.ascontiguousarray(density_d, dtype=np.float64).ravel(),
        np.ascontiguousarray(density_a, dtype=np.float64).ravel(),
        np.ascontiguousarray(axis_d, dtype=np.float64).ravel(),
        np.ascontiguousarray(axis_a, dtype=np.float64).ravel(),
        np.ascontiguousarray(r0_d, dtype=np.float64).ravel(),
        np.ascontiguousarray(r0_a, dtype=np.float64).ravel(),
        float(r0_6), float(kf), int(step)),
        dtype=np.float64).reshape(ng, ng, ng)


def fret_rate_map(
    density_donor,
    density_acceptor,
    r0_donor,
    r0_acceptor,
    dg_donor,
    dg_acceptor,
    forster_radius: float,
    kf: float,
    acceptor_step: int = 2,
) -> np.ndarray:
    """An effective FRET rate for every donor voxel, from the acceptor cloud.

    A fixed donor position does not have *a* FRET rate: it has a distribution of
    them, one per accessible acceptor position, and the FRET-induced donor decay
    is ``sum_i x_i exp(-k_i t)``. This approximates that sum by a single
    exponential whose rate is the reciprocal of the **mean transfer time**::

        k_eff = 1 / <1 / k_i>_A

    which is the *harmonic* mean of the rates, not the arithmetic one. The two
    differ, and the choice matters: the harmonic mean is dominated by the slow
    (distant) acceptor positions, the arithmetic mean by the fast (close) ones.
    :func:`IMP.bff.fret_rate_trace` takes the arithmetic mean instead, because
    it models a fast-exchanging acceptor where the *rates* average.

    :param kf: the donor's radiative rate ``1/tau0``.
    :param acceptor_step: stride over the acceptor grid; the cost is the product
        of the two grids' sizes.
    """
    density_donor = np.ascontiguousarray(density_donor, dtype=np.float64)
    density_acceptor = np.ascontiguousarray(density_acceptor, dtype=np.float64)
    if float(kf) <= 0.0:
        raise ValueError("kf (the donor's radiative rate) must be positive.")
    return _fret_map(
        density_donor,
        density_acceptor,
        grid_axis(density_donor.shape[0], dg_donor),
        grid_axis(density_acceptor.shape[0], dg_acceptor),
        np.asarray(r0_donor, dtype=np.float64),
        np.asarray(r0_acceptor, dtype=np.float64),
        float(forster_radius) ** 6,
        float(kf),
        max(1, int(acceptor_step)),
    )
