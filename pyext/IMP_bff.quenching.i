/*
 * The fields a quenched dye lives in: stickiness, PET rate, FRET rate.
 *
 * Four headers under one shim because the Python that used to wrap them was one
 * module: every function here is a *named composition* of a kernel two headers
 * down -- `slow_factor_grid` is `stamp_spheres` with the multiplying combine,
 * `quenching_rate_map` is `quenching_map` with the axis built for it -- and the
 * names are the ones the physics uses. The renames below exist only so a flat
 * `ng^3` buffer can come back as an `(ng, ng, ng)` cube; nothing here computes.
 *
 * `radial_diffusion_map` is the one function that stays Python, because it takes
 * a *callable* of the distance from the anchor. A Python function is a Python
 * object; there is no C++ spelling of it that does not amount to calling back
 * into the interpreter per voxel.
 */

%rename(_sphere_points) IMP::bff::sphere_points;
%rename(_solvent_accessible_surface) IMP::bff::solvent_accessible_surface;
%rename(_quenching_rate_per_frame) IMP::bff::quenching_rate_per_frame;
%rename(_slow_factor_grid) IMP::bff::slow_factor_grid;
%rename(_quenching_rate_grid) IMP::bff::quenching_rate_grid;
%rename(_av_contact_mask) IMP::bff::av_contact_mask;
%rename(_grid_axis) IMP::bff::grid_axis;
%rename(_diffusion_coefficient_map) IMP::bff::diffusion_coefficient_map;
%rename(_quenching_rate_map) IMP::bff::quenching_rate_map;
%rename(_fret_rate_map) IMP::bff::fret_rate_map;
%rename(_fret_rate_trace) IMP::bff::fret_rate_trace;
%rename(_fret_rate_pair_trace) IMP::bff::fret_rate_pair_trace;

%include "IMP/bff/SolventAccessibleSurface.h"
%include "IMP/bff/QuenchingGrid.h"
%include "IMP/bff/QuenchingMap.h"
%include "IMP/bff/FRETRateTrace.h"

%pythoncode %{
def _flat(a, dtype=np.float64):
    return np.ascontiguousarray(np.asarray(a, dtype=dtype)).ravel()


def _kappa2(value):
    """The orientation factor to pass down; ``None`` is the isotropic 2/3.

    Resolved here rather than by a C++ sentinel: `kappa2` is physically in
    [0, 4] and NaN is an *error* there -- one arriving from a failed orientation
    calculation must raise rather than be answered with the isotropic average --
    so there is no value left over to mean "unspecified".
    """
    return float(kappa2_isotropic()) if value is None else float(value)


def _cube(flat):
    """A flat grid as ``(ng, ng, ng)``, taking ``ng`` from its own length."""
    n = int(np.asarray(flat).size)
    ng = int(round(n ** (1.0 / 3.0)))
    for cand in (ng, ng - 1, ng + 1):
        if cand > 0 and cand ** 3 == n:
            return np.asarray(flat, dtype=np.float64).reshape(cand, cand, cand)
    raise ValueError("a grid of %d values is not a cube" % n)


def sphere_points(n):
    """*n* roughly equidistant points on the unit sphere (golden spiral)."""
    return np.asarray(_IMP_bff._sphere_points(int(n)),
                      dtype=np.float64).reshape(-1, 3)


def solvent_accessible_surface(xyz, vdw, probe_atom_indices, points=None,
                               probe=1.0, radius=2.5):
    """Accessible surface area (A^2) of each atom in *probe_atom_indices*.

    :param xyz: ``(n, 3)`` atom coordinates.
    :param vdw: ``(n,)`` van der Waals radii.
    :param probe_atom_indices: indices of the atoms to measure.
    :param points: unit-sphere sample points; :func:`sphere_points` by default.
    :param probe: probe radius.
    :param radius: radius at which the sphere points are placed.
    """
    if points is None or np.asarray(points).ndim != 2:
        points = np.zeros(0)
    return np.asarray(_IMP_bff._solvent_accessible_surface(
        _flat(xyz), _flat(vdw),
        [int(i) for i in np.asarray(probe_atom_indices).ravel()],
        _flat(points), float(probe), float(radius)), dtype=np.float64)


def quenching_rate_per_frame(collided, k_quench):
    """Sum the rates of the quenching atoms the dye touched, frame by frame.

    The **per-atom** route to a quenching trace, as against sampling a stamped
    rate grid along the trajectory. Use it when the contact flags are what you
    have -- from a distance calculation against explicit atoms rather than from
    a voxel map.

    :param collided: ``(n_frames, n_atoms)`` flags, non-zero where the dye was
        within the critical distance of that atom in that frame.
    :param k_quench: ``(n_atoms,)`` quenching rate per atom, in 1/ns.
    :returns: ``(n_frames,)`` total quenching rate.
    """
    collided = np.ascontiguousarray(collided)
    k_quench = _flat(k_quench)
    if collided.ndim != 2:
        raise ValueError("`collided` must be (n_frames, n_atoms).")
    if collided.shape[1] != k_quench.shape[0]:
        raise ValueError(
            "`collided` has %d atoms but `k_quench` has %d."
            % (collided.shape[1], k_quench.shape[0]))
    return np.asarray(_IMP_bff._quenching_rate_per_frame(
        [int(v) for v in collided.astype(np.int32).ravel()],
        int(collided.shape[0]), k_quench), dtype=np.float64)


def slow_factor_grid(density, ng, dg, slow_radius, rs, r0, slow_fact):
    """Per-voxel diffusion scaling from overlapping sticky spheres.

    Stickiness **multiplies** where spheres overlap. Voxels outside the AV keep
    1.0, so the factor is only meaningful where the walk can go.
    """
    return _cube(_IMP_bff._slow_factor_grid(
        _flat(density), int(ng), float(dg), _flat(slow_radius), _flat(rs),
        _flat(r0), _flat(slow_fact)))


def quenching_rate_grid(density, ng, dg, radius, rs, r0, values):
    """Per-voxel quenching rate (1/ns) from overlapping quencher spheres.

    Rates **add up** where contact spheres overlap, which is the physical
    composition rule for independent PET channels. Same geometry and the same
    indexing convention as :func:`slow_factor_grid`; only the accumulator
    differs.
    """
    return _cube(_IMP_bff._quenching_rate_grid(
        _flat(density), int(ng), float(dg), _flat(radius), _flat(rs),
        _flat(r0), _flat(values)))


def av_contact_mask(density, ng, dg, slow_radius, rs, r0):
    """The contact ("slow") part of an AV density grid, as a binary mask."""
    ng = int(ng)
    mask = _IMP_bff._av_contact_mask(
        _flat(density), ng, float(dg), _flat(slow_radius), _flat(rs), _flat(r0))
    return np.ascontiguousarray(mask, dtype=np.uint8).reshape(ng, ng, ng)


def grid_axis(ng, dg):
    """Voxel-centre offsets from the grid anchor, in Angstrom."""
    return np.asarray(_IMP_bff._grid_axis(int(ng), float(dg)),
                      dtype=np.float64)


def radial_diffusion_map(density, dg, f):
    """A diffusion map from a radial profile ``f(distance from the anchor)``.

    Used to give the dye a mobility that depends on how far the linker is
    stretched, before :func:`slow_diffusion_near_atoms` adds the local crowding.
    Vectorised, so *f* must accept an array.

    Python, and it has to be: *f* is a Python callable, and there is no C++
    spelling of one that does not call back into the interpreter per voxel.
    """
    ng = np.asarray(density).shape[0]
    axis = grid_axis(ng, dg)
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    return np.asarray(f(np.sqrt(x * x + y * y + z * z)), dtype=np.float64)


def slow_diffusion_near_atoms(d_map, density, r0, dg, atoms_xyz, min_distance,
                              slow_factor):
    """Scale a diffusion map down once per atom in contact with each voxel.

    The more atoms crowd a voxel the slower the dye moves there; the factor is
    applied once per contacting atom, so it compounds. Voxels outside the
    accessible volume are set to zero -- the dye cannot be there, and a non-zero
    coefficient would let the solver leak population into the protein.

    J. Chem. Phys. 126, 044707 (2007), eq. (3).
    """
    ng = np.asarray(density).shape[0]
    return _cube(_IMP_bff.slow_near_atoms(
        _flat(d_map), _flat(density), grid_axis(ng, dg), _flat(r0),
        _flat(atoms_xyz), float(min_distance) ** 2, float(slow_factor)))


def diffusion_coefficient_map(density, r0, dg, atoms_xyz, free_diffusion,
                              min_distance, slow_factor, radial_profile=None):
    """The mobility field: a base coefficient, slowed by nearby atoms.

    :param free_diffusion: the unhindered dye diffusion coefficient (A^2/ns),
        used where *radial_profile* is not given.
    :param radial_profile: optional callable of the distance from the anchor,
        replacing the constant base coefficient.
    """
    base = (np.zeros(0) if radial_profile is None
            else radial_diffusion_map(density, dg, radial_profile))
    return _cube(_IMP_bff._diffusion_coefficient_map(
        _flat(density), _flat(r0), float(dg), _flat(atoms_xyz),
        float(free_diffusion), float(min_distance), float(slow_factor),
        _flat(base)))


def quenching_rate_map(density, r0, dg, atoms_xyz, kQ, rC, tau0, dye_radius):
    """Total decay rate per voxel: intrinsic plus exponential PET.

    ``k(r) = 1/tau0 + sum_a kQ_a * exp(-(|r - r_a| - r_dye) / rC_a)``

    The distance is measured from the dye **surface**. Unlike the contact-sphere
    law in :func:`quenching_rate_grid` there is no cut-off: a distant atom
    contributes exponentially little rather than nothing. Voxels outside the
    accessible volume stay at zero -- **not** at ``1/tau0``.
    """
    return _cube(_IMP_bff._quenching_rate_map(
        _flat(density), _flat(r0), float(dg), _flat(atoms_xyz), _flat(kQ),
        _flat(rC), float(tau0), float(dye_radius)))


def fret_rate_map(density_donor, density_acceptor, r0_donor, r0_acceptor,
                  dg_donor, dg_acceptor, forster_radius, kf, acceptor_step=2):
    """An effective FRET rate for every donor voxel, from the acceptor cloud.

    A fixed donor position does not have *a* FRET rate: it has a distribution of
    them, one per accessible acceptor position. This approximates that sum by a
    single exponential whose rate is the reciprocal of the **mean transfer
    time** -- the *harmonic* mean of the rates, not the arithmetic one, so it is
    dominated by the slow (distant) acceptor positions.
    :func:`fret_rate_trace` takes the arithmetic mean instead, because it models
    a fast-exchanging acceptor where the *rates* average.

    :param kf: the donor's radiative rate ``1/tau0``.
    :param acceptor_step: stride over the acceptor grid; the cost is the product
        of the two grids' sizes.
    """
    return _cube(_IMP_bff._fret_rate_map(
        _flat(density_donor), _flat(density_acceptor), _flat(r0_donor),
        _flat(r0_acceptor), float(dg_donor), float(dg_acceptor),
        float(forster_radius), float(kf), int(acceptor_step)))


def fret_rate_trace(trajectory, acceptor_points, R0, tau0, r_min=7.0,
                    max_acceptor_points=512, kappa2=None):
    """Per-frame FRET rate (1/ns) along a donor trajectory.

    For every donor position the transfer rate is averaged over the acceptor's
    accessible volume::

        k_FRET(t) = (1/tau0) * < (R0 / |r_D(t) - r_A|)^6 >_A * (kappa2 / (2/3))

    This resolves the **donor** in time -- the point of simulating its diffusion
    at all -- while treating the acceptor as sampling its own volume quickly
    compared with the donor's excited-state lifetime. That is the *fast-acceptor
    limit*, and it is not valid if the acceptor is immobilised; use
    :func:`fret_rate_pair_trace` then.

    :param kappa2: orientation factor; the isotropic 2/3 by default.
    """
    return _IMP_bff._fret_rate_trace(
        _flat(trajectory), _flat(np.asarray(acceptor_points)[..., :3]),
        float(R0), float(tau0), float(r_min), int(max_acceptor_points),
        _kappa2(kappa2))


def fret_rate_pair_trace(donor_trajectory, acceptor_trajectory, R0, tau0,
                         r_min=7.0, kappa2=None):
    """Per-frame FRET rate (1/ns) from *two* trajectories.

    Both dyes are resolved in time::

        k_FRET(t) = (1/tau0) * (R0 / |r_D(t) - r_A(t)|)^6 * (kappa2 / (2/3))

    This is the general case; :func:`fret_rate_trace` averages the acceptor over
    its accessible volume instead. The two do not agree, and the reason is
    **occupancy**, not the averaging order: a cloud weights every accessible
    voxel the same, while a trajectory is weighted by where the dye actually
    spends its time.

    Unequal frame counts raise rather than truncating -- see the C++ note.
    """
    return _IMP_bff._fret_rate_pair_trace(
        _flat(donor_trajectory), _flat(acceptor_trajectory), float(R0),
        float(tau0), float(r_min),
        _kappa2(kappa2))
%}
