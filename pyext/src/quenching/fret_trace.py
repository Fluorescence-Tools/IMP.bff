"""Time-resolved FRET rate along a dye trajectory.

:mod:`IMP.bff.representation.distance` owns the *static* pair metrics -- ensemble averages
over two accessible-volume clouds. These are their per-frame counterparts, and
they live here because they need a **trajectory**, which only the diffusion
model in this package produces. They are also the hot inner loop of a run
(``n_frames * n_acceptor``, routinely 10^9 distance evaluations), so unlike the
Monte-Carlo samplers in ``fret/distance.py`` they are jitted rather than
vectorised.

Moved here from QuEst (``quest/core/av.py``) by PRD-109. QuEst's own OKF had
already classified these as belonging upstream.
"""

from __future__ import annotations

import numpy as np

import IMP.bff

from IMP.bff.fret.kappa2 import kappa2_isotropic

__all__ = [
    "MAX_KAPPA2",
    "fret_rate_trace",
    "fret_rate_pair_trace",
]

#: kappa2 ranges over [0, 4]: 0 for perpendicular dipoles, 4 for collinear ones.
MAX_KAPPA2 = 4.0


def _rate_trace(trajectory, acceptor_points, R0, tau0, r_min2, kappa2_scale):
    """FRET rate per frame against a static acceptor cloud. **C++.**

    The **arithmetic mean of rates** over the cloud -- the fast-exchange limit,
    where the acceptor re-randomises within the donor's excited-state lifetime.
    ``fret_rate_map`` accumulates the mean transfer *time* and inverts it, which
    is the static limit. Different physics, not two spellings.
    """
    return np.asarray(IMP.bff.fret_rate_trace_kernel(
        np.ascontiguousarray(trajectory, dtype=np.float64).ravel(),
        np.ascontiguousarray(acceptor_points, dtype=np.float64).ravel(),
        float(R0), float(tau0), float(r_min2), float(kappa2_scale)),
        dtype=np.float64)


def _rate_pair_trace(donor, acceptor, R0, tau0, r_min2, kappa2_scale):
    """FRET rate per frame from two trajectories, paired frame by frame. **C++.**"""
    return np.asarray(IMP.bff.fret_rate_pair_trace_kernel(
        np.ascontiguousarray(donor, dtype=np.float64).ravel(),
        np.ascontiguousarray(acceptor, dtype=np.float64).ravel(),
        float(R0), float(tau0), float(r_min2), float(kappa2_scale)),
        dtype=np.float64)


def _kappa2_scale(kappa2) -> float:
    """Turn an orientation factor into a multiplier on ``R0**6``.

    A published Forster radius is quoted *at* the isotropic 2/3, so the transfer
    rate carries the ratio ``kappa2 / (2/3)`` -- equivalently the ``1.5 *
    kappa2`` of the literature -- and passing 2/3 leaves the rate as it was.

    **This model does not compute kappa2, and cannot**: its dye is a
    structureless sphere diffusing in a volume, with no transition dipole to
    orient. What this does is make the *assumption* explicit and adjustable
    instead of silently baked into R0. Modelling it needs a rotational degree of
    freedom per dye and a correlation time to go with it -- a change to the
    physics, not a call to a library. :func:`IMP.bff.kappa2_from_dipoles` does
    compute it, but it takes the two dipole vectors, which is exactly what a
    point in a volume has not got. Use :class:`IMP.bff.RotamerEnsemble` when
    dipoles matter.
    """
    value = float(kappa2_isotropic() if kappa2 is None else kappa2)
    if not np.isfinite(value) or value < 0.0 or value > MAX_KAPPA2:
        raise ValueError(
            f"kappa2 must be in [0, {MAX_KAPPA2:g}] (0 = perpendicular dipoles, "
            f"4 = collinear); got {value!r}."
        )
    return value / kappa2_isotropic()


def fret_rate_trace(
    trajectory,
    acceptor_points,
    R0: float,
    tau0: float,
    r_min: float = 7.0,
    max_acceptor_points: int = 512,
    kappa2=None,
) -> np.ndarray:
    """Per-frame FRET rate (1/ns) along a donor trajectory.

    For every donor position the transfer rate is averaged over the acceptor's
    accessible volume::

        k_FRET(t) = (1/tau0) * < (R0 / |r_D(t) - r_A|)^6 >_A * (kappa2 / (2/3))

    This resolves the **donor** in time -- the point of simulating its diffusion
    at all -- while treating the acceptor as sampling its own volume quickly
    compared with the donor's excited-state lifetime. That is the *fast-acceptor
    limit*, and it is not valid if the acceptor is immobilised; use
    :func:`fret_rate_pair_trace` then.

    :param trajectory: ``(n_frames, 3)`` donor positions in Angstrom.
    :param acceptor_points: ``(n, 3)`` acceptor AV cloud in Angstrom.
    :param R0: Forster radius in Angstrom.
    :param tau0: unquenched donor lifetime in ns.
    :param r_min: closest approach of the two dye centres, in Angstrom. Both
        clouds are centre positions and can overlap in space, where
        ``(R0/r)^6`` would diverge; two dye spheres cannot interpenetrate, so
        the distance is floored at the sum of their radii.
    :param max_acceptor_points: cap on the acceptor points averaged over, taken
        by an even stride. The average converges long before the full cloud, and
        the cost is ``n_frames * n_acceptor``.
    :param kappa2: orientation factor; the isotropic 2/3 by default.
    """
    trajectory = np.ascontiguousarray(trajectory, dtype=np.float64).reshape((-1, 3))
    points = np.ascontiguousarray(acceptor_points, dtype=np.float64).reshape((-1, 3))
    if points.shape[0] == 0:
        raise ValueError("Acceptor accessible volume is empty; cannot compute FRET rates.")
    if tau0 <= 0.0:
        raise ValueError("tau0 must be positive to derive a FRET rate.")
    if points.shape[0] > max_acceptor_points:
        stride = int(np.ceil(points.shape[0] / float(max_acceptor_points)))
        points = np.ascontiguousarray(points[::stride])
    r_min = max(float(r_min), 1e-6)
    return _rate_trace(
        trajectory, points, float(R0), float(tau0), r_min * r_min, _kappa2_scale(kappa2)
    )


def fret_rate_pair_trace(
    donor_trajectory,
    acceptor_trajectory,
    R0: float,
    tau0: float,
    r_min: float = 7.0,
    kappa2=None,
) -> np.ndarray:
    """Per-frame FRET rate (1/ns) from *two* trajectories.

    Both dyes are resolved in time::

        k_FRET(t) = (1/tau0) * (R0 / |r_D(t) - r_A(t)|)^6 * (kappa2 / (2/3))

    This is the general case. :func:`fret_rate_trace` averages the acceptor over
    its accessible volume instead, which is the fast-acceptor limit.

    The two do not agree, and the reason is **occupancy**, not the averaging
    order. A cloud is weighted uniformly -- every accessible voxel counts the
    same -- while a trajectory is weighted by where the dye actually spends its
    time, and a sticky dye dwells near the protein surface rather than out in
    the free volume. Which way that moves the efficiency depends on the
    geometry; it is not a fixed bias. Measured on 148l E15->E90 with a sticky
    acceptor (D = 2.5 A^2/ns, slow radius 11 A) it raised E from 0.73 to 0.78,
    because the sticky region there lies toward the donor.

    The two walks are independent, so pairing frame *i* with frame *i* samples
    the joint distribution correctly **provided both were simulated on the same
    time base** -- same ``t_step`` and the same number of frames. A mismatch
    raises rather than truncating: a trajectory is a *concatenation of one walk
    per excitation*, so cutting it at an arbitrary frame splits a walk and pairs
    the tail of one excitation against the head of another.

    **Unequal lengths are the normal case, not an error case.** The acceptor sits
    at a different site with different residues around it and is therefore
    quenched differently, so the two dyes do not need the same number of
    excitations to collect the same number of photons. Matching them is the
    *caller's* job and the lever is the **number of walks**, not ``t_max`` or
    ``t_step`` -- those are usually already shared and changing them will not
    help.

    :param donor_trajectory: ``(n_frames, 3)`` donor centres in Angstrom.
    :param acceptor_trajectory: ``(n_frames, 3)`` acceptor centres, frame-aligned.
    :param R0: Forster radius in Angstrom.
    :param tau0: unquenched donor lifetime in ns.
    :param r_min: closest approach of the two dye centres, in Angstrom.
    :param kappa2: orientation factor; the isotropic 2/3 by default.
    """
    donor = np.ascontiguousarray(donor_trajectory, dtype=np.float64).reshape((-1, 3))
    acceptor = np.ascontiguousarray(acceptor_trajectory, dtype=np.float64).reshape((-1, 3))
    if donor.shape[0] == 0 or acceptor.shape[0] == 0:
        raise ValueError("Both trajectories must have frames to pair.")
    if tau0 <= 0.0:
        raise ValueError("tau0 must be positive to derive a FRET rate.")
    if donor.shape[0] != acceptor.shape[0]:
        raise ValueError(
            "Donor and acceptor trajectories must have the same number of "
            f"frames to pair ({donor.shape[0]} vs {acceptor.shape[0]}). This is "
            "expected whenever the two dyes are quenched differently -- they sit "
            "at different sites -- so they need different numbers of excitations "
            "for the same photon count. Simulate the same number of *walks* for "
            "both dyes; t_max and t_step are usually already shared and are not "
            "the lever. Truncating here is refused because a trajectory "
            "concatenates one walk per excitation, so an arbitrary cut splits a "
            "walk and pairs one excitation's tail against another's head."
        )
    r_min = max(float(r_min), 1e-6)
    return _rate_pair_trace(
        donor, acceptor, float(R0), float(tau0), r_min * r_min, _kappa2_scale(kappa2)
    )
