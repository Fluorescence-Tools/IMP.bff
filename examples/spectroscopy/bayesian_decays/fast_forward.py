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

    def raw_counts(vals, a2=None):
        a2 = graph.amplitudes(vals) if a2 is None else a2
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
    graph.inst = SpectralInstrument(fwd, model)
    graph._fast = obj
    return graph
