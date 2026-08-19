"""PRD-115 stage 0: the adjoint of the lattice field solver.

``GridDiffusionSolver.gradient`` returns ``dL/dD``, ``dL/dk`` and ``dL/dp0``
for every voxel from one backward sweep. A wrong adjoint does not crash: it
returns a plausible field that a fit then follows to the wrong minimum. So the
check is the dot-product identity against the *forward* solver, which shares
no code with the reverse sweep --

    <dL/dtheta, v>  (adjoint)   ==   d/dh L(theta + h v)  (central difference of run)

-- for random directions ``v`` in each of ``D``, ``k`` and ``p0``, both flux
forms (the Smoluchowski stencil is self-adjoint up to the domain mask, the Ito
one is not), and step counts that do and do not divide the checkpoint length.
The transpose, the ``bounds`` masking, the report seeding and the checkpoint
restart all have to be right for these to agree to roundoff.

Then the same on the five parameters PRD-111 fits, chained through the field formulas,
against the finite-difference Jacobian the identifiability benchmark uses.
"""

import numpy as np

import IMP
import IMP.test

from IMP.bff.sampling.smoluchowski import (
    GridDiffusionGradient,
    GridDiffusionSolver,
    diffusion_stability_limit,
)


def sphere_with_bite(ng=15, radius=5.5):
    """A domain inside the shell: a sphere with a corner missing, so the mask
    has both a smooth and a ragged boundary."""
    c = (ng - 1) / 2.0
    i, j, k = np.mgrid[0:ng, 0:ng, 0:ng]
    bounds = ((i - c) ** 2 + (j - c) ** 2 + (k - c) ** 2 < radius ** 2).astype(float)
    bounds[:6, :6, :] = 0.0
    bounds[0, :, :] = bounds[-1, :, :] = 0.0
    bounds[:, 0, :] = bounds[:, -1, :] = 0.0
    bounds[:, :, 0] = bounds[:, :, -1] = 0.0
    return bounds


def make_solver(rng, ng=15, flux_form="smoluchowski", dg=1.5):
    bounds = sphere_with_bite(ng)
    D = 5.0 + 20.0 * rng.uniform(size=(ng,) * 3)          # A^2/ns
    k = 0.25 + 1.5 * rng.uniform(size=(ng,) * 3)          # 1/ns
    p0 = bounds * (0.5 + rng.uniform(size=(ng,) * 3))
    dt = 0.5 * diffusion_stability_limit(D.max(), dg)
    return GridDiffusionSolver(D, bounds, p0, rate_map=k, t_step=dt, dg=dg,
                               flux_form=flux_form), p0


def loss(solver, p0, w, n_steps, n_out):
    """L = <w, F> for a fresh run from p0 (run() mutates the density)."""
    solver.density = p0.copy()
    return float(np.dot(w, solver.run(n_steps, n_out).fluorescence))


class DotProductIdentityTests(IMP.test.TestCase):

    def _check(self, flux_form, n_steps, n_out, seed):
        rng = np.random.default_rng(seed)
        solver, p0 = make_solver(rng, flux_form=flux_form)
        n_reports = n_steps // n_out + 1
        w = rng.uniform(-1, 1, size=n_reports)

        g = solver.gradient(w, n_steps, n_out, density=p0)
        self.assertIsInstance(g, GridDiffusionGradient)
        for name in ("d_diffusion", "d_rate", "d_density"):
            self.assertEqual(getattr(g, name).shape, p0.shape)

        h = 1e-4
        for name, field, attr, scale in (
            ("D", solver.diffusion_map, "diffusion_map", 1.0),
            ("k", solver.rate_map, "rate_map", 0.1),
            ("p0", p0, None, 1.0),
        ):
            v = rng.uniform(-1, 1, size=field.shape) * scale
            if attr is None:
                v = v * solver.bounds
            base = field.copy()
            if attr is not None:
                setattr(solver, attr, base + h * v)
                lp = loss(solver, p0, w, n_steps, n_out)
                setattr(solver, attr, base - h * v)
                lm = loss(solver, p0, w, n_steps, n_out)
                setattr(solver, attr, base)
            else:
                lp = loss(solver, p0 + h * v, w, n_steps, n_out)
                lm = loss(solver, p0 - h * v, w, n_steps, n_out)
            fd = (lp - lm) / (2 * h)
            adj = float(np.sum(getattr(g, {"D": "d_diffusion", "k": "d_rate", "p0": "d_density"}[name]) * v))
            self.assertAlmostEqual(
                adj, fd, delta=1e-7 * (1.0 + abs(fd)),
                msg=f"{flux_form} {name}: adjoint {adj:.12g} vs FD {fd:.12g} "
                    f"(n_steps={n_steps}, n_out={n_out})")

    def test_smoluchowski_steps_divide_the_checkpoint_length(self):
        self._check("smoluchowski", 100, 10, seed=1)

    def test_smoluchowski_steps_do_not_divide(self):
        self._check("smoluchowski", 47, 5, seed=2)

    def test_ito(self):
        self._check("ito", 64, 8, seed=3)

    def test_ito_odd_lengths(self):
        self._check("ito", 33, 4, seed=4)

    def test_wrong_dL_dF_length_raises(self):
        rng = np.random.default_rng(5)
        solver, p0 = make_solver(rng)
        with self.assertRaises(ValueError):
            solver.gradient(np.ones(3), 100, 10, density=p0)


class ParameterJacobianTests(IMP.test.TestCase):
    """The five PRD-111 parameters through the map builders: the adjoint
    chained through numpy must reproduce the finite-difference Jacobian the
    identifiability benchmark computes."""

    def test_matches_central_differences_on_a_toy_site(self):
        rng = np.random.default_rng(11)
        ng, dg = 17, 1.5
        bounds = sphere_with_bite(ng, radius=6.5)
        c = (ng - 1) // 2
        ax = (np.arange(ng) - c) * dg
        X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
        # a handful of "atoms" around the domain, some of them quenchers
        atoms = np.array([[7.5, 0.0, 0.0], [-2.0, 6.0, 1.0], [0.0, -7.0, -2.0], [3.0, 3.0, 7.0]])
        kq_atom = np.array([1.5, 0.0, 2.5, 0.0])
        rc_atom = np.array([1.2, 1.2, 1.2, 1.2])
        p0 = bounds * (1.0 + 0.2 * rng.uniform(size=(ng,) * 3))
        # theta = (free_diffusion, slow_factor, contact_distance, kQ_scale, rC_scale),
        # the PRD-111 vector; the fields are the benchmark's formulas in numpy
        # (D = D_free * slow_factor^(#atoms within contact_distance),
        #  k = 1/tau0 + sum_a kQ_a exp(-(|r - r_a| - r_dye)/rC_a)) so this test
        # does not depend on the C++ map kernels -- the solver adjoint is the
        # part under test.
        theta0 = np.array([12.0, 0.95, 4.0, 1.0, 1.0])
        tau0, r_dye = 4.0, 1.0
        n_steps, n_out = 200, 10
        dt = 0.5 * diffusion_stability_limit(40.0, dg)   # theta-independent, as PRD-111 requires
        dist = np.stack([np.sqrt((X - a[0]) ** 2 + (Y - a[1]) ** 2 + (Z - a[2]) ** 2) for a in atoms])

        def fields(theta):
            fd_, sf, cd, kqs, rcs = theta
            # smooth "count" so contact_distance is differentiable in the test
            count = np.sum(1.0 / (1.0 + np.exp((dist - cd) / 0.3)), axis=0)
            D = fd_ * sf ** count
            k = 1.0 / tau0 + np.sum(kq_atom[:, None, None, None] * kqs
                                    * np.exp(-(dist - r_dye) / (rc_atom[:, None, None, None] * rcs)), axis=0)
            return D, k

        def decay(theta):
            D, k = fields(theta)
            s = GridDiffusionSolver(D, bounds, p0.copy(), rate_map=k, t_step=dt, dg=dg)
            return s.run(n_steps, n_out).fluorescence

        F0 = decay(theta0)
        w = rng.uniform(-1, 1, size=F0.shape)

        # adjoint on the fields, then chain rule through the map builders by
        # finite differences of the *maps* (cheap: no solve) -- the solver
        # gradient is the expensive, tested part.
        D0, k0 = fields(theta0)
        s = GridDiffusionSolver(D0, bounds, p0.copy(), rate_map=k0, t_step=dt, dg=dg)
        g = s.gradient(w, n_steps, n_out, density=p0)
        grad_adj = np.zeros(5)
        for i in range(5):
            h = 1e-6 * max(1.0, abs(theta0[i]))
            tp, tm = theta0.copy(), theta0.copy()
            tp[i] += h
            tm[i] -= h
            Dp, kp = fields(tp)
            Dm, km = fields(tm)
            grad_adj[i] = (np.sum(g.d_diffusion * (Dp - Dm)) + np.sum(g.d_rate * (kp - km))) / (2 * h)

        # the benchmark's way: central differences of the whole solve
        grad_fd = np.zeros(5)
        for i in range(5):
            h = 1e-4 * max(1.0, abs(theta0[i]))
            tp, tm = theta0.copy(), theta0.copy()
            tp[i] += h
            tm[i] -= h
            grad_fd[i] = (np.dot(w, decay(tp)) - np.dot(w, decay(tm))) / (2 * h)

        scale = np.abs(grad_fd).max()
        np.testing.assert_allclose(grad_adj, grad_fd, rtol=1e-4, atol=1e-6 * scale)


if __name__ == "__main__":
    IMP.test.main()
