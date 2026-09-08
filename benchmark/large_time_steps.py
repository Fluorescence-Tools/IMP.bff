"""Can this propagation take larger time steps, and what does each method cost?

`diffusion_propagate` runs four thousand explicit Euler sweeps because the
stability bound `dt <= dg^2/(6D)` says it must. That bound is what makes the
solver cost `dg^-5`, and it is the reason `quench_pinn` exists as a surrogate.
This script asks what can be done about the bound itself.

Two properties of the operator decide the answer, and both are worth stating
because a general PDE has neither:

* **It does not change over a run.** `d`, `decay` and `bounds` are fixed, so
  the whole propagation is one matrix exponential, `p(t) = exp(-L t) p0`,
  rather than a sequence of different problems.
* **It is symmetrisable.** The Smoluchowski flux form is symmetric as it
  stands. The Ito form is not, but `D^(1/2) A D^(-1/2)` is symmetric to
  machine precision -- the change of variables `q = sqrt(d) p` costs a
  diagonal multiply. So the symmetric Krylov machinery serves both.

Everything is measured against `scipy.sparse.linalg.expm_multiply`, i.e.
against the exact solution of the semi-discrete equation, not against the
shipped scheme. That matters: the shipped scheme has an error of its own
(first-order Lie splitting plus forward Euler) and it is the budget every
other method has to fit inside, not a target to reproduce.

    python benchmark/large_time_steps.py        # needs scipy, not IMP.bff
"""
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

NG, NSTEP, NOUT = 41, 4000, 100


# --------------------------------------------------------------------------
# The fixture, as a sparse operator
# --------------------------------------------------------------------------
def fields(ng, seed=3):
    """The same fields as `gpu_diffusion_wgsl.py`, so the numbers compare."""
    rng = np.random.default_rng(seed)
    i = np.arange(ng) - ng // 2
    X, Y, Z = np.meshgrid(i, i, i, indexing="ij")
    inside = (X**2 + Y**2 + Z**2) <= (ng // 2 - 2) ** 2
    cur = (inside / inside.sum()).astype(np.float64)
    d = (0.10 + 0.05 * rng.random((ng, ng, ng))) * inside
    decay = np.exp(-(0.02 + 0.20 * (X**2 + Y**2 + Z**2) / (ng**2)) * 0.005) * inside + (~inside)
    return cur, d, decay, inside.astype(np.float64)


def operator(ng, flux="smoluchowski", start="smooth", quench="mild", seed=3):
    """A, K and p0 on the active voxels.

    The shipped step is `p <- exp(-K) (I - A) p`; the exact one is
    `p <- exp(-(A+K)) p`. `d` already carries `D dt/dg^2` and `decay` is
    `exp(-k dt)`, so time is measured in shipped steps throughout.
    """
    cur, D, decay, B = fields(ng, seed)
    interior = np.zeros((ng, ng, ng), bool)
    interior[1:-1, 1:-1, 1:-1] = True
    idx = np.flatnonzero((B != 0) & interior)
    pos = -np.ones(ng ** 3, np.int64)
    pos[idx] = np.arange(len(idx))
    ix, iy, iz = np.unravel_index(idx, (ng, ng, ng))
    d0 = D[ix, iy, iz]
    rows, cols, vals = [], [], []
    diag = np.zeros(len(idx))
    for ox, oy, oz in [(-1,0,0),(1,0,0),(0,-1,0),(0,1,0),(0,0,-1),(0,0,1)]:
        jx, jy, jz = ix + ox, iy + oy, iz + oz
        b, dm = B[jx, jy, jz], D[jx, jy, jz]
        if flux == "smoluchowski":
            a = 0.5 * (d0 + dm) * b
            diag += a; off = -a
        else:                                    # Ito: (d_c p_c - d_m p_m) b_m
            diag += d0 * b; off = -dm * b
        j = pos[jx * ng * ng + jy * ng + jz]
        m = j >= 0
        rows.append(np.arange(len(idx))[m]); cols.append(j[m]); vals.append(off[m])
    rows.append(np.arange(len(idx))); cols.append(np.arange(len(idx))); vals.append(diag)
    A = sp.csr_matrix((np.concatenate(vals),
                       (np.concatenate(rows), np.concatenate(cols))),
                      shape=(len(idx),) * 2)
    if flux != "smoluchowski":
        s = np.sqrt(d0)                          # q = sqrt(d) p macht A symmetrisch
        A = sp.diags(s) @ A @ sp.diags(1.0 / s)
    k = -np.log(decay.ravel()[idx])
    if quench == "strong":                       # ein Quencher statt eines glatten Feldes
        c = (ng // 2 + ng // 5, ng // 2, ng // 2)
        r2 = (ix - c[0])**2 + (iy - c[1])**2 + (iz - c[2])**2
        k = k + 300.0 * k.mean() * np.exp(-r2 / 4.0)
    if start == "smooth":
        p0 = cur.ravel()[idx].copy()
    else:                                        # Punktquelle: ein Voxel
        p0 = np.zeros(len(idx))
        p0[np.argmin((ix - ng // 2 - ng // 4)**2 + (iy - ng // 2)**2 + (iz - ng // 2)**2)] = 1.0
    return A.tocsr(), sp.diags(k), k, p0


# --------------------------------------------------------------------------
# The methods
# --------------------------------------------------------------------------
def shipped(A, k, p0, times):
    """What runs today: forward Euler on the diffusion, exact on the decay."""
    ex = np.exp(-k)
    p = p0.copy(); out = [p.sum()]; nv = 0
    for s in range(int(times[-1])):
        p = ex * (p - A @ p); nv += 1
        if (s + 1) % NOUT == 0: out.append(p.sum())
    return np.array(out), nv


def rkl2(L, p0, s, tau, T):
    """Runge-Kutta-Legendre, second order: an explicit method whose stability
    polynomial is a shifted Legendre polynomial, so its real-axis stability
    interval grows as s^2 rather than staying at 2. `s` stencil applications
    advance the solution by up to `s(s+1)/4` explicit steps. Nothing else
    changes -- no solver, no new data structure, the same kernel in the loop.
    Meyer, Balsara & Aslam, J. Comput. Phys. 257 (2014) 594."""
    b = np.empty(s + 1); b[0] = b[1] = 1.0 / 3.0
    for j in range(2, s + 1):
        b[j] = (j * j + j - 2.0) / (2.0 * j * (j + 1.0))
    if s >= 2: b[2] = 1.0 / 3.0
    a = 1.0 - b
    w1 = 4.0 / (s * s + s - 2.0)
    F = lambda y: -(L @ y)
    p = p0.copy(); t = 0.0; nv = 0; out = [(0.0, p.sum())]
    while t < T - 1e-9:
        h = min(tau, T - t)
        Y0 = p; M0 = F(Y0); nv += 1
        Ym1 = Y0 + b[1] * w1 * h * M0
        Ym2 = Y0
        for j in range(2, s + 1):
            mu = (2.0 * j - 1.0) / j * b[j] / b[j - 1]
            nu = -(j - 1.0) / j * b[j] / b[j - 2]
            mut = mu * w1
            Ym2, Ym1 = Ym1, (mu * Ym1 + nu * Ym2 + (1.0 - mu - nu) * Y0
                             + mut * h * F(Ym1) - a[j - 1] * mut * h * M0)
            nv += 1
        p = Ym1; t += h; out.append((t, p.sum()))
    return np.array(out), nv


def lanczos(L, p0, m, times, dtype=np.float64):
    """The whole propagation as one matrix exponential.

    Three vectors of memory and no stored basis: the observable is a plain
    sum, so the only thing the recursion has to keep from each Krylov vector
    is its column sum `u_j = 1^T v_j`, one scalar. The trace at *every* output
    time then comes out of the same m x m tridiagonal:

        F(t) = beta * u^T exp(-T_m t) e_1 = beta * sum_i c_i exp(-lam_i t)

    -- which is a sum of m exponentials, i.e. already the form a fluorescence
    decay is fitted with. Reorthogonalisation is deliberately absent: it costs
    m^2 n flops, more than the steps it replaces, and it is measurably
    unnecessary here (see the note).

    The final density needs the basis after all: either store it, or run the
    recursion a second time with the coefficients known. Not done here."""
    v = p0.astype(dtype); beta = float(np.linalg.norm(v))
    vj = v / dtype(beta); vp = np.zeros_like(vj)
    al = np.zeros(m); be = np.zeros(max(m - 1, 1)); u = np.zeros(m)
    mm = m
    for j in range(m):
        u[j] = vj.sum()
        w = L @ vj
        al[j] = float(vj @ w)
        w = w - dtype(al[j]) * vj - (dtype(be[j - 1]) * vp if j else 0)
        if j < m - 1:
            be[j] = float(np.linalg.norm(w))
            if be[j] < 1e-12: mm = j + 1; break
            vp, vj = vj, w / dtype(be[j])
    al, be, u = al[:mm], be[:mm - 1], u[:mm]
    lam, Q = np.linalg.eigh(np.diag(al) + np.diag(be, 1) + np.diag(be, -1))
    y = (np.exp(-np.outer(times, lam)) * Q[0]) @ Q.T
    return beta * (y @ u), mm


# --------------------------------------------------------------------------
def main():
    times = np.arange(0, NSTEP + 1, NOUT, dtype=float)
    A, K, k, p0 = operator(NG)
    L = (A + K).tocsr()
    rho = spl.eigsh(L, k=1, which="LM", return_eigenvectors=False)[0]
    print("ng=%d, %d active voxels, %d steps" % (NG, len(p0), NSTEP))
    print("rho(A+K) = %.3f; forward Euler is stable to 2, so the shipped step "
          "sits at %.0f%% of the limit" % (rho, 100 * rho / 2))
    print("stiffness ||L||*T = %.0f\n" % (rho * NSTEP))

    t0 = time.perf_counter()
    ref = spl.expm_multiply(-L, p0, start=0.0, stop=float(NSTEP),
                            num=len(times), endpoint=True).sum(axis=1)
    print("reference: expm_multiply, %.1f s" % (time.perf_counter() - t0))
    err = lambda tr: np.abs(np.asarray(tr) - ref).max() / np.abs(ref).max()

    tr, nv = shipped(A, k, p0, times)
    budget = err(tr)
    print("shipped scheme:            %5d matvecs   error %.2e"
          "   <- the budget every other method has to fit inside\n" % (nv, budget))

    print("--- RKL2 super-time-stepping (explicit, same kernel) ---")
    print("%-4s %10s %9s %9s %11s" % ("s", "steps/cyc", "matvecs", "vs 4000", "error"))
    for s in (4, 6, 8, 12, 16, 24, 32):
        tau = 0.95 * (2.0 / rho) * (s * s + s - 2.0) / 4.0
        out, nv = rkl2(L, p0, s, tau, float(NSTEP))
        e = err(np.interp(times, out[:, 0], out[:, 1]))
        print("%-4d %10.1f %9d %8.1fx %11.2e%s"
              % (s, tau, nv, NSTEP / nv, e, "  <- over budget" if e > budget else ""))

    print("\n--- Krylov (Lanczos) on exp(-Lt): one basis, all output times ---")
    print("%-4s %9s %9s %13s %13s" % ("m", "matvecs", "vs 4000", "error f64", "error f32"))
    for m in (32, 64, 96, 128, 192, 256):
        e64 = err(lanczos(L, p0, m, times)[0])
        e32 = err(lanczos(L.astype(np.float32), p0, m, times, np.float32)[0])
        print("%-4d %9d %8.1fx %13.2e %13.2e" % (m, m, NSTEP / m, e64, e32))

    print("\n--- how hard the start is decides m ---")
    print("%-34s %9s %13s" % ("case", "m", "error f32"))
    for flux in ("smoluchowski", "ito"):
        for start in ("smooth", "point"):
            for quench in ("mild", "strong"):
                A2, K2, k2, q0 = operator(NG, flux, start, quench)
                L2 = (A2 + K2).tocsr()
                r2 = spl.expm_multiply(-L2, q0, start=0.0, stop=float(NSTEP),
                                       num=len(times), endpoint=True).sum(axis=1)
                for m in (64, 128, 192, 256):
                    tr, _ = lanczos(L2.astype(np.float32), q0, m, times, np.float32)
                    e = np.abs(tr - r2).max() / np.abs(r2).max()
                    if e < budget:
                        print("%-34s %9d %13.2e"
                              % ("%s / %s / %s" % (flux, start, quench), m, e))
                        break
                else:
                    print("%-34s %9s %13s" % ("%s / %s / %s" % (flux, start, quench),
                                              ">256", "-"))


if __name__ == "__main__":
    main()
