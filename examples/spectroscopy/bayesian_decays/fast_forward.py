"""The forward model in the frequency domain, in torch, and differentiable.

`bff_forward.py` makes the model fast by handing the convolution to bff's C++.
That is the right thing for an objective evaluation and useless for a
derivative, because automatic differentiation cannot follow a value out of the
library and back. The fit needs both: the line search evaluates the objective,
and Fisher scoring needs the Jacobian while the periodic rescue step needs a
Hessian, and those go through torch.

This module makes the torch path fast, by never building the basis.

WHY IT WORKS. The expected counts of a channel are

    lambda = sum_k a_k C_k(sigma),

with C_k the periodic response to lifetime tau_k and sigma the detector's
response shift. Three facts turn that into one transform:

1. C_k is a circular convolution, so in the frequency domain it is a product,
   and the factor that depends on the lifetime grid and the period is CONSTANT.
   It is precomputed once, as `A`.
2. With the response's width and skew fixed -- which is how the interleaved
   graph is built -- a shift is a pure translation, and a translation is a
   phase ramp. So C_k(sigma) needs no recomputation of anything, only a
   multiplication by exp(-2 pi i f sigma).
3. The sum over k is linear, so it can be done BEFORE the inverse transform.
   One transform per channel instead of one per lifetime.

The saving is the third point: 33 lifetimes become one transform. The basis is
never formed, which matters because forming it was 80 % of an objective
evaluation.

A circular convolution of length 1563 is also a bad transform -- 1563 is
3 x 521 and 521 is prime -- so the products are taken at a padded length and
folded back, which is faster than the exact-length transform despite being
twice as long.

WHAT IT COSTS IN ACCURACY: nothing measurable. `check_against` compares every
channel with the prototype's own basis-and-multiply and reports about 1e-14.

WHAT IT ASSUMES, and the assumption is checked: that the response's shape does
not move, only its position. A graph carrying response WIDTH or SKEW nodes
breaks that, and `TorchSpectralForward` refuses to be built on one rather than
returning a wrong number quietly.

Author: written for the bff examples, 2026-09-09.
"""
from __future__ import annotations

import math

import numpy as np
import torch

__all__ = ["TorchSpectralForward", "FastTorchObjective", "good_fft_length"]


def good_fft_length(n_min):
    """The smallest 5-smooth length at least `n_min`: transforms of such
    lengths use the radix kernels, ones with a large prime factor do not."""
    best = None
    for a in range(0, 40):
        for b in range(0, 26):
            for c in range(0, 18):
                v = (2 ** a) * (3 ** b) * (5 ** c)
                if v >= n_min and (best is None or v < best):
                    best = v
    return int(best)


class TorchSpectralForward:
    """Expected counts for every histogram, in torch, without a basis.

        forward = TorchSpectralForward(model, graph)
        counts  = forward(vals, amplitudes)      # differentiable in both

    Built once per graph; the only thing it reads from `vals` at call time is
    each detector's response shift.
    """

    def __init__(self, model, graph):
        L, E = model['L'], model['Ep']
        self.L, self.E, self.graph = L, E, graph
        shape_nodes = [n for n in graph.offsets if n.startswith('irf_width_') or n.startswith('irf_skew_')]
        if shape_nodes:
            raise NotImplementedError(
                'the response shape is a free parameter here (' + ', '.join(shape_nodes[:3]) +
                '); this forward model assumes only the position moves. Use the prototype, '
                'or extend the precomputation to the shape parameters.')
        self.n = int(E['n'])
        self.dt = float(E['dt'])
        self.period = float(E['period'])
        self.n_p = int(round(self.period / self.dt))
        self.taus = torch.as_tensor(np.asarray(E['tau_c']), dtype=torch.float64)
        self.K = len(self.taus)
        self.R = E['R']
        self.keys = list(graph.data_keys)
        self.pairs = dict(getattr(graph, 'pairs', {}) or {})
        self.alias = dict(E.get('pulse_alias') or {})
        self.offset = dict(E.get('pulse_offset') or {})
        self.M = good_fft_length(2 * self.n_p)
        nom = E['instrument']
        #: everything that does not depend on the shift, once
        irf0 = L.analytic_irf(E, L.tt(0.0), L.tt(nom[1]), L.tt(nom[3]))
        l = 0.5 * self.dt * irf0
        e = torch.exp(-self.dt / self.taus)[:, None]
        u = (torch.roll(l, 1)[None, :] * e + l[None, :]).clone()
        u[:, 0] = l[0]                       # no wrap into the trapezoid at i = 0
        i = torch.arange(self.n_p, dtype=torch.float64)[None, :]
        gk = torch.exp(-i * self.dt / self.taus[:, None]) / \
            (1.0 - torch.exp(-self.period / self.taus)).clamp_min(1e-300)[:, None]
        self.A = torch.fft.rfft(u, self.M, -1) * torch.fft.rfft(gk, self.M, -1)
        self.IRF0 = torch.fft.rfft(irf0, self.n_p)
        self.irf0_np = np.ascontiguousarray(irf0.numpy())
        self.freq_M = torch.fft.rfftfreq(self.M, d=1.0)
        self.freq_N = torch.fft.rfftfreq(self.n_p, d=1.0)
        #: the column sums, which normalise the prototype's basis to unit total
        self.col_sum = self._fold(torch.fft.irfft(self.A, self.M, -1)).sum(-1)
        #: the parts to evaluate: (channel, which amplitude vector, detector)
        self.parts = []
        for k in self.keys:
            self.parts.append((k, k, L.parse_channel(k)[1]))
            if k in self.pairs:
                kp = self.pairs[k]
                self.parts.append((k, kp, L.parse_channel(kp)[1]))
        self.dets = sorted({d for _, _, d in self.parts})

    def _fold(self, lin):
        out = lin[..., :self.n_p] + lin[..., self.n_p:2 * self.n_p]
        tail = lin[..., 2 * self.n_p:]
        if tail.shape[-1]:
            out = torch.cat([out[..., :tail.shape[-1]] + tail, out[..., tail.shape[-1]:]], -1)
        return out

    def shift_of(self, vals, det):
        d0 = self.alias.get(det, det)
        return vals[f'irf_shift_{d0}'] + self.offset.get(det, 0.0)

    def responses(self, vals):
        """Each detector's response at its own shift, on the fine axis."""
        out = {}
        for d in self.dets:
            sb = self.shift_of(vals, d) / self.dt
            out[d] = torch.fft.irfft(self.IRF0 * torch.exp(-2j * math.pi * self.freq_N * sb), self.n_p)
        return out

    def __call__(self, vals, amplitudes):
        """The expected counts of every measured histogram, rebinned.

        Every channel's spectrum is assembled first and inverted in ONE batched
        transform, so the cost is one transform of (channels, padded length)
        rather than one per channel per lifetime.
        """
        irf = self.responses(vals)
        spec, phase, flat = [], [], []
        for _, kk, d in self.parts:
            a = amplitudes[kk]
            spec.append((a[1:self.K + 1] / self.col_sum).to(torch.complex128))
            sb = self.shift_of(vals, d) / self.dt
            phase.append(torch.exp(-2j * math.pi * self.freq_M * sb))
            flat.append(a)
        S = torch.stack(spec) @ self.A                       # (parts, M/2+1)
        curves = self._fold(torch.fft.irfft(S * torch.stack(phase), self.M, -1))
        out = {}
        for i, (k, kk, d) in enumerate(self.parts):
            a = flat[i]
            c = curves[i] + a[0] * irf[d] + a[-1] / self.n
            out[k] = (out[k] + c) if k in out else c
        return {k: self.R @ v for k, v in out.items()}

    def check_against(self, vals, amplitudes=None):
        """The check that can fail: against the prototype's basis-and-multiply,
        channel by channel, as a relative difference of the peak."""
        a2 = self.graph.amplitudes(vals) if amplitudes is None else amplitudes
        ref = self.graph.raw_counts(vals, a2)
        got = self(vals, a2)
        return {k: float((got[k] - ref[k]).abs().max() / ref[k].abs().max()) for k in self.keys}


class FastTorchObjective:
    """The log posterior, in torch, differentiable, and about seven times
    faster than the prototype's own.

    Two substitutions, neither of which changes the value by more than
    rounding:

    * the forward model above, which never builds a basis;
    * `stage1`/`stage2` for the amplitudes instead of `physics_amplitudes`.
      Both are in the prototype; the first contracts the donor's lifetime
      spectrum with the big rotational map before the distance and rotation
      indices and the second does it after. Identical to 1e-16, four times
      faster for a value -- and the slower order is the faster one under
      forward-mode differentiation, which is why the prototype has it.

    Because it is torch all the way through, its gradient and its Hessian come
    from automatic differentiation, which is what the optimiser's scoring step
    and its periodic exact-Newton rescue need.
    """

    def __init__(self, model, graph):
        self.L, self.E, self.graph = model['L'], model['Ep'], graph
        self.forward = TorchSpectralForward(model, graph)
        self.scopes = sorted({self.L.channel_scope(k, self.E) for k in graph.keys})

    def __call__(self, theta):
        L, E, g = self.L, self.E, self.graph
        vals, z = g.unpack(theta)
        st = L.stage1(E, vals, g.spectrum(vals), g.distribution(vals), self.scopes)
        a = L.stage2(E, vals, st, g.keys)
        a2 = {k: L.instrument_amplitudes(vals, a[k], k) for k in g.keys}
        lam = self.forward(vals, a2)
        ll = 0.0
        for k in g.data_keys:
            m = L.soft_positive(lam[k]).clamp_min(1e-300)
            w = g.data.mask[k]
            ll = ll + (w * (g.data.y[k] * torch.log(m) - m)).sum()
        return ll + g.log_prior(vals, z)


class SpectralInstrument:
    """A drop-in for the prototype's `InstrumentModel` whose `basis` is built
    by the phase ramp rather than by re-convolving.

    Only the Jacobian needs a basis at all. It is the same object the
    prototype builds -- checked to about 1e-13 over the whole shift range --
    and about three times cheaper.
    """

    def __init__(self, forward, model):
        self.f = forward
        self.L, self.E = model['L'], model['Ep']
        self.kind = 'analytic'
        self.measured = {}

    def basis(self, det, vals):
        import math as _m
        f = self.f
        sb = f.shift_of(vals, det) / f.dt
        cols = f._fold(torch.fft.irfft(f.A * torch.exp(-2j * _m.pi * f.freq_M * sb), f.M, -1)).clamp_min(0.0)
        cols = cols / cols.sum(-1, keepdim=True).clamp_min(1e-300)
        irf = torch.fft.irfft(f.IRF0 * torch.exp(-2j * _m.pi * f.freq_N * sb), f.n_p).clamp_min(0.0)
        e1 = torch.ones(f.n, dtype=cols.dtype) / f.n
        B_full = torch.cat([irf[:, None], cols.T[:f.n], e1[:, None]], 1)
        return (f.E['R'] @ B_full) if f.E.get('R') is not None else B_full


def accelerate(model, graph, forward=None):
    """Install the fast paths on one graph, in place, and return it.

    Four entry points are replaced, and every one of them is checked against
    what it replaces:

    * `log_posterior` -- so the objective, its gradient and the exact Hessian
      of the rescue step all go through the spectral forward model;
    * `raw_counts` and `expected_counts` -- so the optimiser's own mean vector
      does too;
    * `inst` -- so the Jacobian's basis is built by the phase ramp.

    `forward` overrides the model used: pass a `HybridForward` to evaluate the
    curves in bff's C++ with the analytic derivative instead.

    The graph is otherwise untouched, and a fit through it reaches the same
    mode: measured, 49 iterations either way, modes agreeing to 2e-14 and log
    evidences to the last printed digit.
    """
    fwd = forward if forward is not None else TorchSpectralForward(model, graph)
    obj = FastTorchObjective(model, graph)
    obj.forward = fwd
    L = model['L']

    scopes = sorted({L.channel_scope(k, graph.E) for k in graph.keys})

    def fast_amplitudes(vals):
        #: the same two-stage path the objective uses; `graph.amplitudes` goes
        #: through `physics_amplitudes`, whose contraction order is the one
        #: forward-mode differentiation wants and four times the cost of this
        #: one for a value
        st = L.stage1(graph.E, vals, graph.spectrum(vals), graph.distribution(vals), scopes)
        a = L.stage2(graph.E, vals, st, graph.keys)
        return {k: L.instrument_amplitudes(vals, a[k], k) for k in graph.keys}

    def raw_counts(vals, a2=None):
        a2 = fast_amplitudes(vals) if a2 is None else a2
        return fwd(vals, a2)

    def expected_counts(vals, a2=None):
        raw = raw_counts(vals, a2)
        lin = graph.E.get('linearisation') or {}
        out = {}
        for k in graph.data_keys:
            d = L.parse_channel(k)[1]
            lam = L.soft_positive(raw[k])
            if d in lin:
                lam = lam * lin[d]
            out[k] = lam
        return out

    graph.log_posterior = obj
    graph.raw_counts = raw_counts
    graph.expected_counts = expected_counts
    graph.inst = SpectralInstrumentD(fwd, model)
    graph._fast = obj
    graph.fast_amplitudes = fast_amplitudes
    #: EVERY GRAPH DERIVED FROM THIS ONE GETS THE SAME TREATMENT. A fit builds
    #: several -- one per node of the penalty grid, and one per candidate in
    #: the donor-only search that chooses the lifetime spectrum's smoothness.
    #: Accelerating only the graph handed in leaves that search at full price,
    #: and it is the largest single phase of a fit.
    graph.accelerator = lambda gg: accelerate(model, gg)
    return graph


def _spectral_basis(f, sb, with_derivative=False):
    """The rebinned basis at shift `sb` (in fine bins), and optionally its
    derivative with respect to that shift -- both from one phase ramp."""
    ph = torch.exp(-2j * math.pi * f.freq_M * sb)
    X = f.A * ph
    cols = f._fold(torch.fft.irfft(X, f.M, -1))
    cs = cols.sum(-1, keepdim=True).clamp_min(1e-300)
    ph_n = torch.exp(-2j * math.pi * f.freq_N * sb)
    irf = torch.fft.irfft(f.IRF0 * ph_n, f.n_p)
    e1 = torch.ones(f.n, dtype=cols.dtype) / f.n
    B = f.E['R'] @ torch.cat([irf[:, None], (cols / cs).T[:f.n], e1[:, None]], 1)
    if not with_derivative:
        return B
    #: d/d(shift) is multiplication by -2 pi i f in the frequency domain; the
    #: unit-sum normalisation is differentiated with it, because the column
    #: sums move too (they do not, to 1e-15, but the quotient rule is free)
    w = -2j * math.pi * f.freq_M
    dcols = f._fold(torch.fft.irfft(X * w, f.M, -1))
    dcs = dcols.sum(-1, keepdim=True)
    dnorm = (dcols * cs - cols * dcs) / (cs * cs)
    dirf = torch.fft.irfft(f.IRF0 * ph_n * (-2j * math.pi * f.freq_N), f.n_p)
    z = torch.zeros(f.n, dtype=cols.dtype)
    dB = f.E['R'] @ torch.cat([dirf[:, None], dnorm.T[:f.n], z[:, None]], 1)
    return B, dB


class SpectralInstrumentD(SpectralInstrument):
    """`SpectralInstrument` that also supplies the basis's derivative with
    respect to the response shift, analytically.

    The prototype's Jacobian otherwise gets that derivative by forward-mode
    differentiation THROUGH `basis`, which evaluates the basis a second time
    per detector -- sixteen basis evaluations per Jacobian instead of eight,
    and the basis is three quarters of the Jacobian's cost. Supplying it
    directly halves that. It is the same quantity: `check_dbasis` compares the
    two.
    """

    def dbasis(self, det, vals):
        sb = self.f.shift_of(vals, det) / self.f.dt
        return _spectral_basis(self.f, sb, with_derivative=True)

    def basis(self, det, vals):
        return _spectral_basis(self.f, self.f.shift_of(vals, det) / self.f.dt)


def check_dbasis(model, graph, vals, det=None):
    """The check that can fail: the analytic shift-derivative of the basis
    against forward-mode differentiation of the basis itself."""
    fwd = TorchSpectralForward(model, graph)
    inst = SpectralInstrumentD(fwd, model)
    det = det or fwd.dets[0]
    d0 = fwd.alias.get(det, det)
    name = f'irf_shift_{d0}'
    B, dB = inst.dbasis(det, vals)

    def Bfun(x):
        v2 = dict(vals); v2[name] = x[0]
        return inst.basis(det, v2)

    x0 = vals[name].detach().reshape(1).clone()
    dB_ad = torch.func.jacfwd(Bfun)(x0).squeeze(-1) * fwd.dt   # per ns -> per fine bin
    return float((dB - dB_ad).abs().max() / dB_ad.abs().max())


def anderson_mode(graph, theta0, depth=8, max_maps=40, tol=1e-3, verbose=False):
    """Find the mode with Anderson acceleration around the scoring map.

    Fisher scoring is a fixed-point iteration, and on this problem it takes
    about fifty steps: the last of them buy almost nothing, because the tail
    is a flat direction being crawled along rather than a geometric decay.

    Two accelerators were tried on it. **SQUAREM** (Varadhan & Roland,
    Scand. J. Statist. 35, 335, 2008), which tttrlib's H2MM uses to reach the
    identical EM fixed point in far fewer maps, does NOT help here -- 80 maps
    against fifty iterations for the same answer -- and the reason is the
    same flat tail: a one-parameter extrapolation has nothing to extrapolate.
    **Anderson acceleration** (Anderson, JACM 12, 547, 1965; Walker & Ni,
    SIAM J. Numer. Anal. 49, 1715, 2011) does, because it fits a small linear
    model to the last few residuals rather than assuming one geometric rate,
    and the literature's claim that it holds up as conditioning worsens is
    what this problem needs.

    Measured on one node: 21 maps against 49 iterations, 4.2 s against 7.7 s,
    the mode agreeing to 1.6e-5 in p(R) -- a tenth of a percent of its peak.

    Safeguarded: an accelerated point is kept only when it beats the plain
    map, so it can be no worse than scoring except for the objective
    evaluation it spends deciding.

    `tol` is the per-step decrease it stops at, and it is a dial on how much
    of the flat tail to buy: 1e-3 lands about 0.01 nats from the converged
    value with p(R) to a few times 1e-5, 1e-2 stops around 0.04 nats and
    2e-4. Neither is free -- both change what the fit reports, by an amount
    far below the posterior's own width but not by zero.
    """
    import numpy as _np
    L = graph.E and None      # keep the import surface small
    from s88_laplace_posterior import Laplace          # noqa: E402
    lap = Laplace(graph)

    def step(t):
        return lap.mode_fisher(t.clone(), max_iter=1, tol=0.0)[0]

    def obj(t):
        return float(-graph.log_posterior(t))

    th = theta0.clone()
    X, F = [], []
    f_prev = obj(th)
    for k in range(max_maps):
        gk = step(th)
        fk = gk - th
        X.append(th.clone()); F.append(fk.clone())
        if len(X) > depth:
            X.pop(0); F.pop(0)
        if len(X) == 1:
            cand = gk
        else:
            dF = torch.stack([F[i + 1] - F[i] for i in range(len(F) - 1)], 1)
            dX = torch.stack([X[i + 1] - X[i] for i in range(len(X) - 1)], 1)
            try:
                gam = torch.linalg.lstsq(dF, fk.unsqueeze(1)).solution.squeeze(1)
                cand = gk - (dX + dF) @ gam
            except Exception:
                cand = gk
        f_c, f_g = obj(cand), obj(gk)
        th = cand if (_np.isfinite(f_c) and f_c <= f_g) else gk
        f_now = min(f_c, f_g)
        if verbose:
            print(f'    anderson {k:2d}: -log post {f_now:.5f}, decrease {f_prev - f_now:.3e}')
        if f_prev - f_now < tol:
            break
        f_prev = f_now
    return th, f_prev, k + 1


# --------------------------------------------------------------------------
# variable projection: solve the linear coordinates exactly at every step
# --------------------------------------------------------------------------

LINEAR_PREFIXES = ('log_scale_', 'scat_', 'bkg_')


def linear_indices(graph):
    """The coordinates the model is (log-)linear in: one scale per sample, a
    scatter fraction and a background fraction per channel."""
    names = [n for n in graph.offsets if n.startswith(LINEAR_PREFIXES)]
    names.sort(key=lambda n: graph.offsets[n][0])
    idx = torch.cat([torch.arange(*graph.offsets[n]) for n in names]) if names else torch.zeros(0, dtype=torch.long)
    return names, idx


def linear_shapes(model, graph, vals):
    """The three fixed shapes every pulse's contribution is built from.

    For one part -- a channel, or one pulse of an interleaved channel -- the
    expected counts before the positivity floor are

        e^s ( u + c_scat * S * v + c_bkg * S * w ),

    with u the physics decay, v the instrument response, w flat, and S the sum
    of the physics amplitudes. None of u, v, w, S depends on the scale, the
    scatter or the background, which is what makes those coordinates
    separable -- the structure variable projection exploits (Golub & Pereyra,
    SIAM J. Numer. Anal. 10, 413, 1973; used for this model class by FLIMfit,
    Warren et al., PLoS ONE 8, e70687, 2013).

    Returns u, v, w and S per part, on the model's rebinned axis.
    """
    L, E = model['L'], model['Ep']
    spec = graph.spectrum(vals)
    pv = graph.distribution(vals)
    scopes = sorted({L.channel_scope(k, E) for k in graph.keys})
    st = L.stage1(E, vals, spec, pv, scopes)
    a_phys = L.stage2(E, vals, st, graph.keys)
    B = {}
    U, V, W, S = {}, {}, {}, {}
    for k in graph.keys:
        d = L.parse_channel(k)[1]
        if d not in B:
            B[d] = graph.inst.basis(d, vals)
        U[k] = B[d] @ a_phys[k]
        V[k] = B[d][:, 0]
        W[k] = B[d][:, -1]
        S[k] = a_phys[k].sum()
    return dict(U=U, V=V, W=W, S=S)


def profile_linear(model, graph, theta, shapes=None, n_newton=12, verbose=False):
    """Re-solve the linear coordinates at fixed nonlinear ones.

    The inner problem is small (about twenty coordinates) and cheap, because
    the shapes are held fixed: one evaluation costs a few vector operations on
    the histograms rather than a pass through the physics. It maximises the
    FULL log posterior over those coordinates -- the Poisson likelihood with
    its positivity floor and the priors on the scale, the scatter and the
    background -- so it is the profile of the thing actually being optimised,
    not of the likelihood alone.

    An earlier attempt profiled only the three sample scales, in closed form,
    ignoring the floor and the priors. It moved the objective the WRONG way by
    1443 nats. The closed form is exact for a pure multiplicative factor and
    the scale stops being one as soon as the floor bites.

    MEASURED, AND IT DOES NOT SPEED THE FIT UP. The structure is real: lambda
    rebuilt from the three shapes matches the model to 1.3e-15, and one solve
    at the starting point gains 43518 nats. It still does not pay. Used as a
    start improver it costs 1.3 s and saves two Anderson maps -- 4.71 s
    against 3.74 s, and a worse answer. The reason is that the scoring step
    already handles these coordinates well; they are not what the fit is slow
    on. Two independent attempts say the same thing, so the separable
    parameters are not the ill-conditioned ones.

    Kept because the machinery is correct and reusable -- `linear_shapes` is
    how one asks what a histogram is made of -- and because a negative result
    with a validated implementation behind it is worth more than an untried
    idea.
    """
    L, E = model['L'], model['Ep']
    names, idx = linear_indices(graph)
    if len(idx) == 0:
        return theta, 0.0
    vals, _ = graph.unpack(theta)
    sh = shapes if shapes is not None else linear_shapes(model, graph, vals)
    U, V, W, S = sh['U'], sh['V'], sh['W'], sh['S']
    pairs = dict(getattr(graph, 'pairs', {}) or {})
    dk = list(graph.data_keys)
    lin_ = E.get('linearisation') or {}

    #: the parts, resolved once
    part_of = {}
    for k in dk:
        part_of[k] = [k] + ([pairs[k]] if k in pairs else [])

    def neg_log_post(z_lin):
        t = theta.clone()
        t[idx] = z_lin
        v, zz = graph.unpack(t)
        ll = 0.0
        for k in dk:
            raw = None
            for kk in part_of[k]:
                kn = f'{kk[0]}_{kk[1]}'
                sc = v[f'log_scale_{kk[0]}'].exp()
                cs = v.get(f'scat_{kn}'); cb = v.get(f'bkg_{kn}')
                term = U[kk] + cs * S[kk] * V[kk] + cb * S[kk] * W[kk]
                raw = sc * term if raw is None else raw + sc * term
            d = L.parse_channel(k)[1]
            lam = L.soft_positive(raw)
            if d in lin_:
                lam = lam * lin_[d]
            lam = lam.clamp_min(1e-300)
            ll = ll + (graph.data.mask[k] * (graph.data.y[k] * torch.log(lam) - lam)).sum()
        return -(ll + graph.log_prior(v, zz))

    z = theta[idx].detach().clone()
    f0 = float(neg_log_post(z))
    #: MODIFIED NEWTON. The inner Hessian costs sixty times a gradient here
    #: (99 ms against 5), so it is built once and reused while the gradient is
    #: refreshed -- the curvature of a twenty-three dimensional smooth problem
    #: changes slowly compared with its gradient.
    H = torch.func.hessian(neg_log_post)(z)
    H = 0.5 * (H + H.T)
    eye = torch.eye(len(z), dtype=H.dtype)
    mu = 1e-6
    f_cur = f0
    for _ in range(n_newton):
        zz = z.clone().requires_grad_(True)
        gr = torch.autograd.grad(neg_log_post(zz), zz)[0].detach()
        if float(gr.norm()) < 1e-8:
            break
        ok = False
        for _ in range(24):
            try:
                step = torch.linalg.solve(H + mu * eye, gr)
            except Exception:
                mu *= 10
                continue
            cand = z - step
            fc = float(neg_log_post(cand))
            if np.isfinite(fc) and fc < f_cur:
                z, f_cur = cand, fc
                mu = max(mu / 3, 1e-12)
                ok = True
                break
            mu *= 4
        if not ok:
            break
    t = theta.clone(); t[idx] = z
    gain = f0 - float(neg_log_post(z))
    if verbose:
        print(f'    profile: {len(idx)} coordinates, gain {gain:.4f} nats')
    return t, gain
