"""The same propagation in PyTorch, against the C++ core and the WGSL prototype.

torch is the obvious thing to reach for: it is already installed wherever the
`quench_pinn` surrogate is worked on, it has MPS on this machine and CUDA on
the compute box, and it needs no plugin, no shader and no C ABI. So before
writing a wgpu backend it is worth knowing what it would be competing with.

Both routes are implemented here on the same fields as
`benchmark/gpu_diffusion_wgsl.py`, so every number in this file compares
directly with every number in that one:

* the **stepping** route -- 4 000 explicit sweeps, written as seven shifted
  multiply-adds with the stencil weights precomputed, which is the same
  arithmetic the tuned WGSL kernel does;
* the **Krylov** route -- the whole propagation as one matrix exponential,
  the method `src/KrylovDiffusion.cpp` implements.

Both eager and under `torch.compile`, because the difference is not small: the
step is seven elementwise operations on a grid, and eager torch launches each
of them separately with a temporary in between, four thousand times over.
Compiling fuses them into one kernel, which is what the C++ and the WGSL
kernel do by construction. Measured at 5.3-13.2x, and it is the difference
between torch being uncompetitive here and torch being competitive.

**Run it in its own process.** torch bundles its own libomp, and a bff built
against another one aborts on import with OMP error #15. `KMP_DUPLICATE_LIB_OK`
silences that, but it is documented as unsafe and it is not something to
measure through -- so this script imports no IMP.bff, and the C++ figures it
is compared against are taken from a separate run.

    python benchmark/torch_diffusion.py
"""
import time

import numpy as np
import torch

NSTEP, NOUT = 4000, 100


def fields(ng, seed=3):
    rng = np.random.default_rng(seed)
    i = np.arange(ng) - ng // 2
    X, Y, Z = np.meshgrid(i, i, i, indexing="ij")
    inside = (X**2 + Y**2 + Z**2) <= (ng // 2 - 2) ** 2
    cur = (inside / inside.sum()).astype(np.float64)
    d = (0.10 + 0.05 * rng.random((ng, ng, ng))) * inside
    decay = np.exp(-(0.02 + 0.20 * (X**2 + Y**2 + Z**2) / ng**2) * 0.005) * inside + (~inside)
    return cur, d, decay, inside.astype(np.float64)


# The six neighbour slices, in the order the weights are built.
SL = [(slice(0, -2), slice(1, -1), slice(1, -1)), (slice(2, None), slice(1, -1), slice(1, -1)),
      (slice(1, -1), slice(0, -2), slice(1, -1)), (slice(1, -1), slice(2, None), slice(1, -1)),
      (slice(1, -1), slice(1, -1), slice(0, -2)), (slice(1, -1), slice(1, -1), slice(2, None))]
IN = (slice(1, -1), slice(1, -1), slice(1, -1))


def coefficients(d, decay, bnds, generator):
    """Seven weights per interior voxel: the self term and the six neighbours.

    `generator=False` gives the step operator `decay*(I-A)`, `generator=True`
    the generator `L = A + K` whose exponential the Krylov route takes."""
    d0, k0, b0 = d[IN], decay[IN], bnds[IN]
    a = [0.5 * (d0 + d[s]) * bnds[s] for s in SL]
    asum = sum(a)
    if generator:
        with np.errstate(divide="ignore"):
            k = np.where(b0 != 0, -np.log(np.where(k0 > 0, k0, 1.0)), 0.0)
        self_ = (asum + k) * b0
        off = [-x * b0 for x in a]
    else:
        self_ = k0 * (1.0 - asum) * b0
        off = [k0 * x * b0 for x in a]
    return self_, off


def to(x, dev, dt):
    """`torch.tensor`, not `torch.as_tensor`: the latter *shares* memory with the
    numpy array, and the stepping loop below writes into its buffers. Sharing
    meant the first run overwrote the caller's initial density with the
    propagated one, so every later run on the same fields silently started from
    a decayed field -- and reported a trace 2x too small rather than an error."""
    return torch.tensor(np.ascontiguousarray(x), dtype=dt, device=dev)


def make_apply(self_, off, compile_=False):
    def apply_(p, out):
        acc = self_ * p[IN]
        for q in range(6):
            acc = acc + off[q] * p[SL[q]]
        out[IN] = acc
        return out
    return torch.compile(apply_, dynamic=False) if compile_ else apply_


def stepping(ng, cur, d, decay, bnds, dev, dt, compile_=False):
    self_, off = coefficients(d, decay, bnds, generator=False)
    self_ = to(self_, dev, dt); off = [to(x, dev, dt) for x in off]
    a = to(cur, dev, dt); b = torch.zeros_like(a)
    step = make_apply(self_, off, compile_)
    if compile_:                      # compile outside the clock
        step(a.clone(), b.clone()); sync(dev)
    sync(dev); t0 = time.perf_counter()
    trace = []
    for s in range(NSTEP):
        if s % NOUT == 0:
            trace.append(a.sum())
        a, b = step(a, b), a
    trace.append(a.sum())
    tr = torch.stack(trace).cpu().numpy()
    sync(dev)
    return tr, time.perf_counter() - t0


def krylov(ng, cur, d, decay, bnds, dev, dt, m, compile_=False):
    self_, off = coefficients(d, decay, bnds, generator=True)
    self_ = to(self_, dev, dt); off = [to(x, dev, dt) for x in off]
    apply_ = make_apply(self_, off, compile_)
    if compile_:
        z = to(cur, dev, dt); apply_(z, torch.zeros_like(z)); sync(dev)
    sync(dev); t0 = time.perf_counter()
    v0 = to(cur, dev, dt)
    beta = v0.norm()
    vj = v0 / beta
    vp = torch.zeros_like(vj); w = torch.zeros_like(vj)
    al, be, u = [], [], []
    # Everything stays on the device: pulling alpha and beta back each step
    # would be two synchronisations per matvec, which on MPS costs more than
    # the matvec.
    for j in range(m):
        apply_(vj, w)
        al.append((vj * w).sum()); u.append(vj.sum())
        if j + 1 >= m:
            break
        w = w - al[j] * vj - (be[j - 1] * vp if j else 0.0)
        be.append(w.norm())
        vp, vj = vj, w / be[j]
        w = torch.zeros_like(vj)
    # .cpu() before .double(): MPS has no float64, so the cast has to happen
    # on the host side of the transfer.
    A = torch.stack(al).cpu().double().numpy()
    B = torch.stack(be).cpu().double().numpy() if be else np.zeros(0)
    U = torch.stack(u).cpu().double().numpy()
    beta = float(beta)
    lam, Q = np.linalg.eigh(np.diag(A) + np.diag(B, 1) + np.diag(B, -1))
    lam = np.maximum(lam, 0.0)
    times = np.arange(0, NSTEP + 1, NOUT, dtype=float)
    amp = beta * (U @ Q) * Q[0]
    tr = (amp[None, :] * np.exp(-np.outer(times, lam))).sum(axis=1)
    sync(dev)
    return tr, time.perf_counter() - t0


def sync(dev):
    if dev == "mps":
        torch.mps.synchronize()


def main(reps=3):
    print("torch %s | mps %s | cuda %s | %d cpu threads"
          % (torch.__version__, torch.backends.mps.is_available(),
             torch.cuda.is_available(), torch.get_num_threads()))
    devices = [("cpu", torch.float64), ("cpu", torch.float32)]
    if torch.backends.mps.is_available():
        devices.append(("mps", torch.float32))   # MPS has no f64
    if torch.cuda.is_available():
        devices += [("cuda", torch.float32), ("cuda", torch.float64)]
    print("\n%-6s %-18s %11s %12s %12s %11s"
          % ("grid", "device", "stepping", "krylov m=64", "krylov m=128", "deviation"))
    for ng in (41, 61, 81):
        cur, d, decay, bnds = fields(ng)
        ref = None
        for dev, dt in devices:
            for compile_ in (False, True):
                def best(fn):
                    out = None; ts = []
                    for _ in range(reps):
                        out, t = fn(); ts.append(t)
                    return out, min(ts) * 1e3
                try:
                    tr_s, t_s = best(lambda: stepping(ng, cur, d, decay, bnds,
                                                      dev, dt, compile_))
                    if ref is None:
                        ref = tr_s
                    _, t_64 = best(lambda: krylov(ng, cur, d, decay, bnds,
                                                  dev, dt, 64, compile_))
                    tr_k, t_128 = best(lambda: krylov(ng, cur, d, decay, bnds,
                                                      dev, dt, 128, compile_))
                except Exception as e:
                    print("ng=%-3d %-18s failed: %s"
                          % (ng, "%s/%s%s" % (dev, str(dt).replace("torch.", ""),
                                              " +compile" if compile_ else ""),
                             str(e)[:60]))
                    continue
                dev_rel = np.abs(tr_k - ref).max() / np.abs(ref).max()
                print("ng=%-3d %-18s %10.1f %11.1f %11.1f %11.1e"
                      % (ng, "%s/%s%s" % (dev, str(dt).replace("torch.", ""),
                                          " +compile" if compile_ else ""),
                         t_s, t_64, t_128, dev_rel))


if __name__ == "__main__":
    main()
