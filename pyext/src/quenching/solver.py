"""Explicit propagation of the dye's excited state on an AV grid.

The field counterpart of the Brownian walk: instead of sampling trajectories and
histogramming photons, integrate

    dp/dt = div(D grad p) - k p

on the accessible-volume grid, with ``D`` the mobility field
(:func:`IMP.bff.diffusion_coefficient_map`) and ``k`` the decay-rate field
(:func:`IMP.bff.quenching_rate_map`, plus a FRET map if there is an acceptor).
Deterministic and free of shot noise, at the cost of a stability limit on the
time step.

Two things it gives that the walk does not: the **equilibrium** occupancy of the
volume (propagate with ``k = 0`` until the distribution stops moving -- the
sticky regions fill up, which a uniform AV density does not represent), and the
donor decay as a smooth curve in one pass.

Moved here from ChiSurf (``chisurf/core/structure/av/functions.py``,
``DiffusionIterator``) by PRD-109. ChiSurf carried an OpenCL backend beside the
Python one; it is not brought across -- ``IMP.bff`` ships through conda-forge as
part of IMP and cannot take a ``pyopencl`` dependency, and the kernel here is
jitted.
"""

from __future__ import annotations

from typing import NamedTuple, Optional

import numpy as np

from .._jit import njit, prange

__all__ = [
    "GridDiffusionResult",
    "GridDiffusionSolver",
    "diffusion_stability_limit",
    "equilibrium_occupancy",
]


class GridDiffusionResult(NamedTuple):
    """Time axis, surviving excited-state fraction, and the final density."""

    time: np.ndarray
    fluorescence: np.ndarray
    density: np.ndarray


def diffusion_stability_limit(d_max: float, dg: float, k_max: float = 0.0) -> float:
    """Largest stable time step for the explicit scheme.

    The updated voxel keeps the coefficient ``1 - 6 D dt/dg^2 - k dt``, so both
    terms constrain the step:

        dt <= 1 / (6 D / dg^2 + k_max)

    With ``k_max = 0`` this is the familiar ``dt <= dg^2 / (6 D)``.

    **The rate term used to be left out, and it dominates.** On T4L site 19 at
    2.5 A the diffusion term contributes 0.16 to that coefficient and the
    quenching term 2.01 -- twelve times more -- so a step this function called
    safe diverged. Worse, the divergence is not always visible: site 124 broke
    the criterion by the same margin and still returned a smooth, finite,
    entirely plausible decay. Since :class:`GridDiffusionSolver` now integrates
    the rate term exponentially it is unconditionally stable in ``k`` and only
    the diffusion term binds, but the honest criterion is kept here for callers
    sizing a step for the plain scheme, and because a caller who asks with a
    ``k_max`` deserves the answer that accounts for it.
    """
    d_max = float(d_max)
    k_max = float(k_max)
    denominator = 6.0 * d_max / float(dg) ** 2 + k_max
    if denominator <= 0.0:
        return float("inf")
    return 1.0 / denominator


def equilibrium_occupancy(
    diffusion_map, bounds, flux_form: str = "smoluchowski"
) -> np.ndarray:
    """The stationary occupancy of the volume, in closed form.

    Which closed form depends on how the flux is discretised, and the two
    disagree completely:

    ``"smoluchowski"`` (default, and the physics)
        ``p_eq ∝ 1`` on the accessible domain -- **independent of D**.
    ``"ito"`` (the inherited ChiSurf behaviour)
        ``p_eq ∝ 1/D``, so the dye piles up wherever it moves slowly.

    **The default is Smoluchowski because equilibrium is thermodynamics and
    mobility is kinetics.** A dye slowed by friction near the protein surface,
    with no attractive interaction, must still be found uniformly across its
    accessible volume at equilibrium -- it merely takes longer to get around.
    Letting a friction field set the distribution asserts a potential that was
    never specified.

    This is also the form the field's own canonical treatment uses: the
    **Haas-Steinberg** equation for diffusion-modulated FRET,

        ∂N(r,t)/∂t = -[1/tau_D + k_T(r)] N + D d/dr [ p(r) d/dr ( N / p(r) ) ]

    is written with exactly this structure so that its stationary state is the
    *given* distance distribution ``p(r)`` for any ``D``. ``p(r)`` comes from the
    chain statistics -- here, from the accessible volume -- and ``D`` is a
    separate kinetic parameter fitted against it. QuEst's notebook
    ``04_diffusion_modulated_fret.ipynb`` integrates that equation, but with a
    constant ``D`` and a uniform ``p(r)``, where both discretisations coincide
    and neither is tested.

    A genuinely *sticky* dye -- one with an attractive interaction, not merely a
    slower one -- is represented by a non-uniform ``p_eq``, which is a separate
    physical input this signature is shaped to accept later. It is not the same
    thing as a mobility field, and conflating them is what the inherited
    ``"ito"`` behaviour did.

    The ``"ito"`` branch below documents what it is:

    :func:`_step` propagates the flux as ``d[i]*p[i] - d[j]*p[j]``, which
    discretises

        ∂p/∂t = ∇²(D p)                    (the Itô / divergence form)

    **not** ``∇·(D ∇p)``. The stationary state of ``∇²(Dp)`` is ``D p = const``,
    i.e. ``p_eq ∝ 1/D`` -- so the dye accumulates where it moves slowly. It is
    not a small effect: on T4L site 132 the ratio of peak to mean occupancy is
    ~75. Verified against the iterative solver to a maximum relative deviation
    of **1.1e-13** (``test_quenching_field.py``).

    Either way, use this rather than :meth:`GridDiffusionSolver.equilibrium`,
    which spends tens of thousands of iterations converging to it -- and in the
    ``"ito"`` form on a real site, with ``D`` varying by orders of magnitude
    through the compounding slow factor, may not converge in any reasonable
    number at all.
    """
    diffusion_map = np.ascontiguousarray(diffusion_map, dtype=np.float64)
    mask = np.asarray(bounds) > 0
    occupancy = np.zeros_like(diffusion_map)
    if flux_form == "smoluchowski":
        usable = mask
    elif flux_form == "ito":
        usable = mask & (diffusion_map > 0.0)
    else:
        raise ValueError(
            f"flux_form must be 'smoluchowski' or 'ito', not {flux_form!r}")
    if not usable.any():
        return occupancy
    if flux_form == "smoluchowski":
        occupancy[usable] = 1.0
    else:
        occupancy[usable] = 1.0 / diffusion_map[usable]
    total = occupancy.sum()
    return occupancy / total if total else occupancy


@njit(cache=True, parallel=True)
def _step_smoluchowski(nxt, cur, d, k, bounds):
    """One explicit Euler step of the **Smoluchowski** form.

    The flux between two voxels is ``D_ij (p_i - p_j)`` with ``D_ij`` the
    interface mobility (arithmetic mean), rather than ``D_i p_i - D_j p_j``.
    The difference is the whole physics: this flux vanishes when ``p`` is
    uniform *whatever* ``D`` does in space, so a spatially varying mobility
    changes how fast the dye redistributes and **not where it ends up**. See
    :func:`equilibrium_occupancy` for why that is the right requirement.
    """
    ng = cur.shape[0]
    for ix in prange(ng):
        if ix == 0 or ix == ng - 1:
            for iy in range(ng):
                for iz in range(ng):
                    nxt[ix, iy, iz] = 0.0
            continue
        for iy in range(ng):
            nxt[ix, iy, 0] = 0.0
            nxt[ix, iy, ng - 1] = 0.0
        for iz in range(ng):
            nxt[ix, 0, iz] = 0.0
            nxt[ix, ng - 1, iz] = 0.0
    for ix in prange(1, ng - 1):
        for iy in range(1, ng - 1):
            for iz in range(1, ng - 1):
                if bounds[ix, iy, iz] == 0:
                    nxt[ix, iy, iz] = 0.0
                    continue
                p0 = cur[ix, iy, iz]
                d0 = d[ix, iy, iz]
                flux = (
                    0.5 * (d0 + d[ix - 1, iy, iz]) * (p0 - cur[ix - 1, iy, iz]) * bounds[ix - 1, iy, iz]
                    + 0.5 * (d0 + d[ix + 1, iy, iz]) * (p0 - cur[ix + 1, iy, iz]) * bounds[ix + 1, iy, iz]
                    + 0.5 * (d0 + d[ix, iy - 1, iz]) * (p0 - cur[ix, iy - 1, iz]) * bounds[ix, iy - 1, iz]
                    + 0.5 * (d0 + d[ix, iy + 1, iz]) * (p0 - cur[ix, iy + 1, iz]) * bounds[ix, iy + 1, iz]
                    + 0.5 * (d0 + d[ix, iy, iz - 1]) * (p0 - cur[ix, iy, iz - 1]) * bounds[ix, iy, iz - 1]
                    + 0.5 * (d0 + d[ix, iy, iz + 1]) * (p0 - cur[ix, iy, iz + 1]) * bounds[ix, iy, iz + 1]
                )
                nxt[ix, iy, iz] = (p0 - flux) * k[ix, iy, iz]
    return nxt


@njit(cache=True, parallel=True)
def _step(nxt, cur, d, k, bounds):
    """One explicit Euler step of ``dp/dt = div(D grad p) - k p``.

    *d* carries ``D * dt / dg^2`` and *k* carries ``exp(-k * dt)``, folded in by
    the caller so the inner loop is pure arithmetic. The rate term is applied as
    a factor rather than subtracted: ``exp(-k dt)`` is the *exact* solution of
    ``dp/dt = -k p`` over the step, so the decay contributes no stability
    constraint at all, where ``1 - k dt`` goes negative and diverges once
    ``k dt > 1`` -- which on a strongly quenched site it does. Flux to a voxel outside the
    volume is dropped by the ``bounds`` factors, which is a no-flux (reflecting)
    wall -- the dye cannot enter the protein.
    """
    ng = cur.shape[0]
    # The 7-point stencil cannot be evaluated on the outer shell, so that shell
    # is written to zero rather than left alone. Leaving it alone is what
    # ChiSurf did, and with a ping-pong pair of buffers those voxels then hold
    # the state from *two* steps ago -- stale data that never decays and never
    # diffuses, silently added into every population sum. It only stayed
    # invisible because an accessible volume never reaches the grid edge, which
    # `GridDiffusionSolver` now checks instead of assuming.
    for ix in prange(ng):
        if ix == 0 or ix == ng - 1:
            for iy in range(ng):
                for iz in range(ng):
                    nxt[ix, iy, iz] = 0.0
            continue
        for iy in range(ng):
            nxt[ix, iy, 0] = 0.0
            nxt[ix, iy, ng - 1] = 0.0
        for iz in range(ng):
            nxt[ix, 0, iz] = 0.0
            nxt[ix, ng - 1, iz] = 0.0
    for ix in prange(1, ng - 1):
        for iy in range(1, ng - 1):
            for iz in range(1, ng - 1):
                if bounds[ix, iy, iz] == 0:
                    nxt[ix, iy, iz] = 0.0
                    continue
                dp = d[ix, iy, iz] * cur[ix, iy, iz]
                flux = (
                    (dp - d[ix - 1, iy, iz] * cur[ix - 1, iy, iz]) * bounds[ix - 1, iy, iz]
                    + (dp - d[ix + 1, iy, iz] * cur[ix + 1, iy, iz]) * bounds[ix + 1, iy, iz]
                    + (dp - d[ix, iy - 1, iz] * cur[ix, iy - 1, iz]) * bounds[ix, iy - 1, iz]
                    + (dp - d[ix, iy + 1, iz] * cur[ix, iy + 1, iz]) * bounds[ix, iy + 1, iz]
                    + (dp - d[ix, iy, iz - 1] * cur[ix, iy, iz - 1]) * bounds[ix, iy, iz - 1]
                    + (dp - d[ix, iy, iz + 1] * cur[ix, iy, iz + 1]) * bounds[ix, iy, iz + 1]
                )
                nxt[ix, iy, iz] = (cur[ix, iy, iz] - flux) * k[ix, iy, iz]
    return nxt


class GridDiffusionSolver:
    """Propagate an excited-state density on a masked 3-D grid.

    :param diffusion_map: ``D`` per voxel, in A^2/ns.
    :param bounds: non-zero where the dye may be; the domain mask.
    :param density: initial population; normalised on the first step.
    :param rate_map: ``k`` per voxel, in 1/ns. Zero (the default) propagates
        pure diffusion, which is how the equilibrium occupancy is obtained.
    :param t_step: integration step in ns. Must satisfy
        :func:`diffusion_stability_limit`; a larger one raises rather than
        returning a diverged field.
    :param dg: voxel edge in Angstrom.
    """

    def __init__(
        self,
        diffusion_map,
        bounds,
        density,
        rate_map=None,
        t_step: float = 1.0,
        dg: float = 1.0,
        check_stability: bool = True,
        flux_form: str = "smoluchowski",
    ):
        self.diffusion_map = np.ascontiguousarray(diffusion_map, dtype=np.float64)
        self.bounds = np.ascontiguousarray(
            np.asarray(bounds) > 0, dtype=np.float64
        )
        self.density = np.ascontiguousarray(density, dtype=np.float64)
        if rate_map is None:
            rate_map = np.zeros_like(self.diffusion_map)
        self.rate_map = np.ascontiguousarray(rate_map, dtype=np.float64)
        self.dg = float(dg)
        self.n_iterations = 0
        self.check_stability = bool(check_stability)
        self.t_step = float(t_step)
        if flux_form not in ("smoluchowski", "ito"):
            raise ValueError(
                f"flux_form must be 'smoluchowski' or 'ito', not {flux_form!r}")
        #: How the flux between voxels is discretised, which decides where the
        #: dye sits at equilibrium. ``"smoluchowski"`` keeps the equilibrium
        #: independent of ``D`` (the Haas-Steinberg structure); ``"ito"`` is the
        #: inherited ChiSurf behaviour, whose equilibrium is ``p ∝ 1/D``. They
        #: coincide when ``D`` is uniform.
        self.flux_form = flux_form
        self._kernel = _step if flux_form == "ito" else _step_smoluchowski

    def _validate(self):
        # The stencil cannot reach the outer shell, so a domain touching it
        # would lose population there with no warning.
        edge = self.bounds.copy()
        edge[1:-1, 1:-1, 1:-1] = 0.0
        if edge.any():
            raise ValueError(
                "The domain reaches the outer shell of the grid, where the "
                "7-point stencil cannot be evaluated and population is "
                "discarded. Enlarge the grid so the accessible volume is "
                "surrounded by at least one empty voxel."
            )
        if not self.check_stability:
            return
        limit = diffusion_stability_limit(float(self.diffusion_map.max()), self.dg)
        if self.t_step > limit:
            raise ValueError(
                f"t_step {self.t_step:g} ns exceeds the explicit-scheme "
                f"stability limit dg^2 / (6 D_max) = {limit:g} ns for "
                f"dg = {self.dg:g} A and D_max = {self.diffusion_map.max():g} "
                "A^2/ns. The scheme would diverge rather than lose accuracy."
            )

    def run(self, n_steps: int, n_out: int = 10) -> GridDiffusionResult:
        """Integrate *n_steps* steps, reporting every *n_out*.

        :returns: the time axis, the surviving excited-state fraction
            ``sum(p) / sum(p_0)`` at each reported step -- the donor decay -- and
            the final density.
        """
        self._validate()
        n_steps = int(n_steps)
        n_out = max(1, int(n_out))

        cur = self.density * self.bounds
        total = cur.sum()
        if total > 0.0:
            cur = cur / total
        nxt = np.zeros_like(cur)

        # Fold dt into the coefficients so the inner loop is pure arithmetic.
        d = self.diffusion_map * (self.t_step / self.dg ** 2)
        # Exact over the step, and unconditionally stable: see `_step`.
        k = np.exp(-self.rate_map * self.t_step)

        n_reports = n_steps // n_out + 1
        time = np.arange(n_reports, dtype=np.float64) * self.t_step * n_out
        fluorescence = np.zeros(n_reports, dtype=np.float64)

        i_out = 0
        for step in range(n_steps):
            if step % n_out == 0 and i_out < n_reports:
                fluorescence[i_out] = cur.sum()
                i_out += 1
            nxt = self._kernel(nxt, cur, d, k, self.bounds)
            # Swap **every** step. ChiSurf swapped only on odd steps
            # (`if time_i % 2 > 0`) while always computing `n <- f(p)`, so every
            # even step recomputed the previous one from a stale buffer and half
            # the evolution was thrown away -- the field advanced at roughly
            # half the requested rate. Fixed on the move (PRD-109); pinned by
            # the free-diffusion variance test, which is what caught it.
            cur, nxt = nxt, cur
            self.n_iterations += 1

        while i_out < n_reports:
            fluorescence[i_out] = cur.sum()
            i_out += 1

        self.density = cur
        return GridDiffusionResult(time, fluorescence, cur)

    def equilibrium(
        self, n_steps: int = 20000, tolerance: float = 1e-8, n_check: int = 100
    ) -> np.ndarray:
        """Propagate with no decay until the occupancy stops moving.

        **Prefer :func:`equilibrium_occupancy`**, which is the same answer in
        closed form (``p ∝ 1/D``, agreeing to 1.1e-13). This iterates toward it,
        which on a real site — where ``D`` spans orders of magnitude through the
        compounding slow factor — can fail to converge in any practical number
        of steps: on T4L site 132 the distribution was still drifting after
        40 000 iterations. Kept as the independent check that the closed form is
        right, and for the case where someone changes the flux discretisation
        and the closed form no longer holds.

        :returns: the normalised equilibrium density.
        """
        self._validate()
        cur = self.density * self.bounds
        total = cur.sum()
        if total > 0.0:
            cur = cur / total
        nxt = np.zeros_like(cur)
        zero_rate = np.ones_like(cur)   # exp(-0 * dt): no decay
        d = self.diffusion_map * (self.t_step / self.dg ** 2)

        previous = cur.copy()
        for step in range(int(n_steps)):
            nxt = self._kernel(nxt, cur, d, zero_rate, self.bounds)
            cur, nxt = nxt, cur
            self.n_iterations += 1
            if (step + 1) % int(n_check) == 0:
                drift = float(np.abs(cur - previous).sum())
                if drift < tolerance:
                    break
                previous = cur.copy()

        total = cur.sum()
        if total > 0.0:
            cur = cur / total
        self.density = cur
        return cur
