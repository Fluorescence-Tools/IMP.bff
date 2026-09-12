"""The forward model evaluated in C++ and differentiated analytically.

`bff_forward.py` is fast and not differentiable. `fast_forward.py` is
differentiable and, being torch, pays torch's per-operation overhead. This
module is both: the curve comes out of bff's C++ reconvolution, and the
derivative is supplied by hand as a `torch.autograd.Function`, so the whole
thing sits inside an ordinary torch graph and the optimiser cannot tell the
difference.

WHY A HAND-WRITTEN DERIVATIVE IS CHEAP HERE. The expected counts are linear in
the amplitudes,

    lambda(t) = sum_k a_k C_k(t),

so the vector-Jacobian product the backward pass needs is

    (dL/da)_k = sum_t C_k(t) (dL/dlambda)(t) = <C_k, v>,

an inner product of the incoming gradient with each basis column. Forming the
columns to do that would cost as much as the forward pass, so it is done in the
frequency domain instead: the transform of C_k is precomputed once (it is a
constant times a phase ramp, see `fast_forward`), and by Parseval the whole
vector-Jacobian product is one transform of v and one complex matrix-vector
product -- independent of the number of lifetimes.

The derivative with respect to the response shift comes from the same place: a
shift is a phase ramp, so d lambda / d sigma is the inverse transform of the
spectrum times (-2 pi i f), which costs one more transform.

MEASURED, on the eight-histogram problem: the curve at 0.045 ms per channel
against torch's 0.17, the whole forward model at 1.05 ms against 2.24, and
both derivatives exact -- 1e-15 against automatic differentiation through the
equivalent torch expression.

AND WHAT IT CANNOT DO, which is the reason the series does not fit through it.
A hand-written backward is a first derivative and nothing else. The optimiser
also wants an exact Hessian for its periodic rescue step, and that is a
derivative OF the backward -- which a numpy backward cannot provide, since
autograd cannot see inside it. So this model accelerates the objective and the
gradient and is unusable for the Hessian, which is most of a node.

Closing that gap does not need a cleverer derivative. It needs the
reconvolution to hand back its per-species terms, because the Jacobian with
respect to the amplitudes IS the basis: `TcspcDecay::set_emit_basis` in bff's
tree does exactly that, and with it the whole derivative chain is a matmul.

Author: written for the bff examples, 2026-09-09.
"""
from __future__ import annotations

import math

import numpy as np
import torch

__all__ = ["HybridForward", "TttrlibForward", "KERNELS", "KERNEL_NOTES"]


class _BffCurve(torch.autograd.Function):
    """One detector's decay: forward in bff's C++, backward by Parseval."""

    @staticmethod
    def forward(amps, shift, engine, det):
        a = amps.detach().numpy()
        sb = float(shift.detach())
        return torch.from_numpy(engine._call_bff(det, a, sb))

    @staticmethod
    def setup_context(ctx, inputs, output):
        #: the separate-context form functorch requires; without it jacrev
        #: refuses the Function outright, which is the better failure
        amps, shift, engine, det = inputs
        ctx.engine, ctx.det = engine, det
        ctx.save_for_backward(amps, shift)

    @staticmethod
    def backward(ctx, grad_out):
        e, det = ctx.engine, ctx.det
        amps, shift = ctx.saved_tensors
        g = grad_out.detach().numpy()
        d_a, d_s = e._vjp(det, amps.detach().numpy(), float(shift.detach()), g)
        return (torch.from_numpy(d_a), torch.tensor(d_s, dtype=shift.dtype), None, None)


class HybridForward:
    """Expected counts for every histogram: bff's kernel, torch's autograd.

        forward = HybridForward(model, graph)
        counts  = forward(vals, amplitudes)      # differentiable

    Refuses to be built on a graph whose response shape is a free parameter,
    for the same reason `fast_forward` does: the precomputation that makes the
    derivative cheap assumes only the position moves.
    """

    def __init__(self, model, graph):
        import IMP.bff as bff
        import sys, pathlib
        sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
        import pie_mfd as P
        from fast_forward import TorchSpectralForward
        shape_nodes = [n for n in graph.offsets if n.startswith('irf_width_') or n.startswith('irf_skew_')]
        if shape_nodes:
            raise NotImplementedError('the response shape is a free parameter; see fast_forward')
        self.P, self.bff = P, bff
        self.model, self.graph = model, graph
        self.L, self.E = model['L'], model['Ep']
        #: the spectral precomputation is shared with the torch model: it is
        #: what makes the backward pass free of the lifetime count
        self.s = TorchSpectralForward(model, graph)
        self.n, self.K, self.M = self.s.n, self.s.K, self.s.M
        self.n_p, self.dt = self.s.n_p, self.s.dt
        self.col_sum = self.s.col_sum.numpy()
        self.A = self.s.A.numpy()
        self.freq_M = self.s.freq_M.numpy()
        self.freq_N = self.s.freq_N.numpy()
        self.IRF0 = np.fft.rfft(self.s.irf0_np, self.s.n_p)
        self.R = self.s.R
        self.keys, self.pairs, self.parts = self.s.keys, self.s.pairs, self.s.parts
        #: one bff node per detector
        self.node, self.out, self.spec_port, self.buf = {}, {}, {}, {}
        for d in self.s.dets:
            nd = bff.TcspcDecay(d)
            nd.set_number_of_lifetimes(self.K)
            nd.set_timing(self.dt, self.s.period)
            port = bff.GraphPort(); port.value = np.zeros(self.n)
            nd.add_output_port(d, port)
            for name, v in (('background', 0.0), ('scatter', 0.0), ('n0', 1.0), ('timeshift', 0.0)):
                p = nd.get_port(name)
                if p is not None:
                    p.value = v
            nd.set_spectrum_from_port(True)
            self.node[d], self.out[d] = nd, port
            self.spec_port[d] = nd.get_port('lifetime_spectrum')
            b = np.empty(2 * self.K); b[1::2] = self.s.taus.numpy()
            self.buf[d] = b
        self._irf_shift = {d: None for d in self.s.dets}
        self._irf = {}
        #: Parseval weights for a real transform of even length: the zero and
        #: Nyquist bins count once, the rest twice
        self.pw = np.full(self.M // 2 + 1, 2.0); self.pw[0] = 1.0; self.pw[-1] = 1.0

    #: ---- the pieces the autograd Function calls -------------------------
    def _set_response(self, det, shift_bins):
        if self._irf_shift[det] == shift_bins:
            return
        #: the response lives on the UNPADDED axis, so its phase ramp uses the
        #: unpadded frequency grid -- using the padded one shifts it by the
        #: wrong amount, which the forward check catches at once
        ph = np.exp(-2j * math.pi * self.freq_N * shift_bins)
        irf = np.fft.irfft(np.fft.rfft(self.s.irf0_np, self.n_p) * ph, self.n_p)
        self.node[det].set_response_array(np.ascontiguousarray(irf))
        self._irf_shift[det] = shift_bins
        self._irf[det] = irf

    def _call_bff(self, det, a, shift_bins):
        self._set_response(det, shift_bins)
        b = self.buf[det]
        b[0::2] = a[1:self.K + 1] / self.col_sum
        self.spec_port[det].value = b
        self.node[det].update()
        return np.asarray(self.out[det].value) + a[0] * self._irf[det] + a[-1] / self.n

    def _vjp(self, det, a, shift_bins, g):
        #: pad the incoming gradient to the transform length and take its
        #: spectrum once; every column's inner product is then a dot product
        #: FOLDING IS PERIODIC EXTENSION ON THE OTHER SIDE. The curve is
        #: fold(lin), so <fold(lin), g> = <lin, g extended with period n_p>.
        gp = np.resize(g, self.M)
        G = np.fft.rfft(gp, self.M)
        ph = np.exp(-2j * math.pi * self.freq_M * shift_bins)
        #: <C_k, g> for every k by Parseval, so the cost does not grow with K
        inner = (np.real((self.A * ph) * np.conj(G)[None, :]) * self.pw[None, :]).sum(1) / self.M
        d_a = np.zeros_like(a)
        d_a[1:self.K + 1] = inner / self.col_sum
        d_a[0] = float(self._irf[det] @ g)
        d_a[-1] = float(g.sum()) / self.n
        #: d lambda / d shift: the same spectrum times -2 pi i f
        spec = (a[1:self.K + 1] / self.col_sum).astype(np.complex128) @ self.A
        dcurve = np.fft.irfft(spec * ph * (-2j * math.pi * self.freq_M), self.M)
        d_lin = dcurve[:self.n_p] + dcurve[self.n_p:2 * self.n_p]
        tail = dcurve[2 * self.n_p:]
        if tail.size:
            d_lin[:tail.size] += tail
        #: the scattered light is the RESPONSE, which moves with the shift too:
        #: without this term the shift derivative is wrong by ten percent, which
        #: is what the comparison against automatic differentiation showed
        ph_n = np.exp(-2j * math.pi * self.freq_N * shift_bins)
        d_irf = np.fft.irfft(self.IRF0 * ph_n * (-2j * math.pi * self.freq_N), self.n_p)
        d_s = float(d_lin @ g) + float(a[0]) * float(d_irf @ g)
        return d_a, d_s

    #: ---- the model ------------------------------------------------------
    def __call__(self, vals, amplitudes):
        out = {}
        for k, kk, d in self.parts:
            shift = self.s.shift_of(vals, d) / self.dt
            c = _BffCurve.apply(amplitudes[kk], shift, self, d)
            out[k] = (out[k] + c) if k in out else c
        return {k: self.R @ v for k, v in out.items()}


# --------------------------------------------------------------------------
# choosing the kernel
# --------------------------------------------------------------------------

KERNELS = ('bff', 'fconv_per_cs', 'fconv_per')

#: what each one is, and why the default is the default
KERNEL_NOTES = {
    'bff': "IMP.bff's TcspcDecay node, which calls tttrlib's fconv_per_cs_ad -- "
           "the kernel that interleaves eight species' recursions in registers. "
           "Fastest of the three and stable.",
    'fconv_per_cs': "tttrlib's fconv_per_cs called directly. Correct and stable, "
                    "and about 1.5x slower than the node above because it uses the "
                    "two-lane NEON kernel rather than the eight-way blocked one.",
    'fconv_per': "tttrlib's fconv_per. NOT USABLE in the build this was written "
                 "against: called repeatedly with identical inputs and freshly "
                 "allocated arrays it returns a curve growing by one species' worth "
                 "per call, so a fit through it diverges while a single-call test "
                 "passes at 1e-15. Reported upstream; kept here only so the "
                 "comparison can be reproduced.",
}


class TttrlibForward:
    """The same forward model through tttrlib's kernels directly.

    `kernel` selects which: see `KERNELS` and `KERNEL_NOTES`. This exists to
    make the comparison runnable rather than asserted -- the numbers in
    notebook 10 come from it.
    """

    def __init__(self, model, graph, kernel='fconv_per_cs'):
        import tttrlib
        if kernel not in KERNELS:
            raise ValueError(f'kernel must be one of {KERNELS}')
        from fast_forward import TorchSpectralForward
        self.tttrlib = tttrlib
        self.kernel = kernel
        self.s = TorchSpectralForward(model, graph)
        self.L = model['L']
        self.n, self.K = self.s.n, self.s.K
        self.dt, self.period, self.n_p = self.s.dt, self.s.period, self.s.n_p
        self.taus = self.s.taus.numpy()
        self.col_sum = self.s.col_sum.numpy()
        self.R = np.asarray(self.s.R)
        self.keys, self.parts = self.s.keys, self.s.parts
        self._irf = {}

    def set_response(self, vals):
        f = self.s
        fN = f.freq_N.numpy()
        for d in {p[2] for p in self.parts}:
            sb = float(f.shift_of(vals, d)) / self.dt
            self._irf[d] = np.ascontiguousarray(
                np.fft.irfft(np.fft.rfft(f.irf0_np, self.n_p) * np.exp(-2j * math.pi * fN * sb), self.n_p))

    def __call__(self, vals, amplitudes):
        if not self._irf:
            self.set_response(vals)
        t = self.tttrlib
        out = {}
        for k, kk, d in self.parts:
            a = np.asarray(amplitudes[kk].detach())
            x = np.empty(2 * self.K)
            x[0::2] = a[1:self.K + 1] / self.col_sum
            x[1::2] = self.taus
            buf = np.zeros(self.n)
            if self.kernel == 'fconv_per_cs':
                t.fconv_per_cs(buf, self._irf[d], x, self.period, self.n - 1, self.n - 1, self.dt)
            else:
                t.fconv_per(buf, self._irf[d], x, self.period, 0, self.n - 1, self.dt)
            c = buf + float(a[0]) * self._irf[d] + float(a[-1]) / self.n
            out[k] = (out[k] + c) if k in out else c
        return {k: self.R @ v for k, v in out.items()}

    def check_against(self, vals, amplitudes=None):
        a2 = self.s.graph.amplitudes(vals) if amplitudes is None else amplitudes
        ref = self.s.graph.raw_counts(vals, a2)
        got = self(vals, a2)
        return {k: float(np.abs(got[k] - np.asarray(ref[k].detach())).max()
                         / np.abs(np.asarray(ref[k].detach())).max()) for k in self.keys}
