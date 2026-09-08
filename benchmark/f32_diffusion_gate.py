"""Das Tor: hält der Stencil in f32?

WGSL kennt kein f64. Der Löser rechnet durchgehend in double. Bevor
irgendeine GPU-Klempnerei geschrieben wird, muss dieselbe Iteration in
float32 gegen die f64-Referenz gehalten werden -- auf den Gittern, die die
Tests benutzen, und über so viele Schritte wie ein echter Lauf.
"""
import numpy as np, IMP.bff as b

def sweep(cur, d, decay, bounds, smol, dtype):
    """Ein Schritt, genau wie internal::lattice_sweep, in der gegebenen Genauigkeit."""
    c, dd, de, bo = (a.astype(dtype) for a in (cur, d, decay, bounds))
    nxt = np.zeros_like(c)
    s = (slice(1, -1),) * 3
    p0, d0 = c[s], dd[s]
    flux = np.zeros_like(p0)
    for ax in range(3):
        for off in (-1, 1):
            sl = [slice(1, -1)] * 3
            sl[ax] = slice(1 + off, (-1 + off) or None)
            m = tuple(sl)
            if smol:
                flux = flux + dtype(0.5) * (d0 + dd[m]) * (p0 - c[m]) * bo[m]
            else:
                flux = flux + (d0 * p0 - dd[m] * c[m]) * bo[m]
    inner = (p0 - flux) * de[s]
    nxt[s] = np.where(bo[s] != 0, inner, 0)
    return nxt

def run(ng, n_steps, n_out, dtype, seed=3):
    rng = np.random.default_rng(seed)
    i = np.arange(ng) - ng // 2
    X, Y, Z = np.meshgrid(i, i, i, indexing="ij")
    inside = (X**2 + Y**2 + Z**2) <= (ng // 2 - 2) ** 2
    cur = (inside / inside.sum()).astype(np.float64)
    # ortsabhängige Mobilität und Löschung, wie ein echtes Feld
    d = (0.10 + 0.05 * rng.random((ng, ng, ng))) * inside
    decay = np.exp(-(0.02 + 0.20 * (X**2 + Y**2 + Z**2) / (ng**2)) * 0.005) * inside + (~inside)
    bounds = inside.astype(np.float64)
    trace = []
    c = cur.astype(dtype)
    for step in range(n_steps):
        if step % n_out == 0:
            trace.append(float(c.sum()))
        c = sweep(c, d, decay, bounds, True, dtype)
    trace.append(float(c.sum()))
    return np.array(trace), c

for ng, n_steps in ((21, 2000), (41, 4000)):
    t64, c64 = run(ng, n_steps, 100, np.float64)
    t32, c32 = run(ng, n_steps, 100, np.float32)
    rel = np.abs(t32 - t64) / np.maximum(np.abs(t64), 1e-30)
    dens = np.abs(c32.astype(np.float64) - c64) / max(np.abs(c64).max(), 1e-30)
    print("ng=%2d, %5d Schritte:  Fluoreszenz max. rel. Abweichung %.3e | "
          "letzte Dichte max. rel. %.3e | Restpopulation %.4f"
          % (ng, n_steps, rel.max(), dens.max(), t64[-1]))
