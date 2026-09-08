"""The decay model as a bff node graph.

The expensive part of fitting a distance distribution is not the linear
algebra and not the statistics: it is evaluating the forward model over and
over. Measured on the eight-histogram problem of this series, one evaluation
of the log posterior costs about 55 ms, and 41 ms of that is rebuilding the
periodic reconvolution of 33 lifetimes for each of eight detectors.

bff does that reconvolution in C++ and does it per curve rather than per
basis. `TcspcDecay` takes a response, a lifetime spectrum and a period, and
returns the convolved, periodic decay -- which is exactly the quantity the
likelihood needs. Asking it for the finished curve instead of building a basis
and multiplying is the whole speedup, and it is a factor of fifteen.

    forward = BffForward(model, graph, vals)
    counts  = forward(amplitudes)          # {channel: expected counts}

The class is deliberately small: one `TcspcDecay` node per detector, its ports
wired once, and an evaluation that writes amplitudes into ports and reads the
curve out. Nothing is cached that the caller could invalidate without saying
so; `set_response` is the one call that must be made when the response
parameters move.

WHAT THIS DOES NOT DO. It computes the model, not its derivative. The fit's
Jacobian with respect to the amplitudes is the basis itself, and getting that
out of one call needs `TcspcDecay::set_emit_basis`, which exists in bff's tree
but not in the build this example runs against. Until then the derivative
comes from the prototype's own automatic differentiation, and this class
accelerates the objective evaluations -- which is where most of the time is.

Author: written for the bff examples, 2026-09-08.
"""
from __future__ import annotations

import numpy as np

__all__ = ["BffForward", "structural_graph"]


class BffForward:
    """The expected counts of every histogram, computed by bff.

    Parameters
    ----------
    model : the environment dict from `bd.build_model`
    graph : the prototype's factor graph, which supplies the channel list, the
            interleaved pairing and the rebinning matrix
    vals  : a value dict, used once to place each detector's response

    The amplitude convention is the prototype's: for channel k, `a[0]` is
    scattered light, `a[1:K+1]` are the amplitudes of the lifetime grid and
    `a[-1]` is background, all against columns normalised to unit sum. bff's
    amplitudes are absolute, so the lifetime amplitudes are divided by the
    column sums -- which are the same for any response, because a circular
    convolution with a unit-sum kernel preserves the integral over the period.
    That identity is checked in `check_against`, not assumed.
    """

    def __init__(self, model, graph, vals):
        import IMP.bff as bff
        import sys
        sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
        import pie_mfd as P
        self.bff = bff
        self.model, self.graph, self.P = model, graph, P
        L, Ep = model['L'], model['Ep']
        self.L, self.Ep = L, Ep
        self.n = int(Ep['n'])
        self.dt = float(Ep['dt'])
        self.period = float(Ep['period'])
        self.taus = np.asarray(Ep['tau_c'], float)
        self.K = len(self.taus)
        self.R = np.asarray(Ep['R'])
        self.keys = list(graph.data_keys)
        self.pairs = dict(getattr(graph, 'pairs', {}) or {})
        self.dets = sorted({L.parse_channel(k)[1] for k in graph.keys})
        #: the sum of each periodic column, once: it does not depend on the response
        nominal = L.analytic_irf(Ep, L.tt(0.0), L.tt(Ep['instrument'][1]), L.tt(Ep['instrument'][3]))
        self.col_sum = L.periodic_columns(Ep['dec'], nominal, L.tt(self.taus), self.period).numpy().sum(1)
        self.nodes, self.out, self.irf, self.spec_port = {}, {}, {}, {}
        #: the interleaved (amplitude, lifetime) buffer the vector port takes
        self._buf = np.empty(2 * self.K)
        self._buf[1::2] = self.taus
        for d in self.dets:
            node = bff.TcspcDecay(d)
            node.set_number_of_lifetimes(self.K)
            node.set_timing(self.dt, self.period)
            port = bff.Port()
            port.value = np.zeros(self.n)
            node.add_output_port(d, port)
            for name, v in (('background', 0.0), ('scatter', 0.0), ('n0', 1.0), ('timeshift', 0.0)):
                p = node.get_port(name)
                if p is not None:
                    p.value = v
            #: THE SPECTRUM GOES IN AS ONE VECTOR, not as 2K scalar ports.
            #: Measured: writing 33 scalar amplitude ports costs 0.05 ms, which
            #: is MORE than the convolution it feeds (0.046 ms). Through the
            #: vector port the same curve costs 0.046 ms in total, and the
            #: result is bit-identical.
            node.set_spectrum_from_port(True)
            self.spec_port[d] = node.get_port('lifetime_spectrum')
            self.nodes[d], self.out[d] = node, port
        self.set_response(vals)

    def set_response(self, vals):
        """Place each detector's instrument response. Call this whenever the
        response parameters move; nothing else here depends on them."""
        for d in self.dets:
            r = np.ascontiguousarray(np.asarray(self.P.fine_response(self.model, d, vals), dtype=float))
            self.irf[d] = r
            self.nodes[d].set_response_array(r)

    def curve(self, det, a):
        """One detector's decay on the fine axis, for one amplitude vector."""
        node, port = self.nodes[det], self.out[det]
        self._buf[0::2] = a[1:self.K + 1] / self.col_sum
        self.spec_port[det].value = self._buf
        node.update()
        return np.asarray(port.value) + float(a[0]) * self.irf[det] + float(a[-1]) / self.n

    def __call__(self, amplitudes):
        """The expected counts of every measured histogram, on the model's
        rebinned axis. A histogram of the labelled sample holds both pulses, so
        its partner scope's curve is added to it."""
        out = {}
        for k in self.keys:
            det = self.L.parse_channel(k)[1]
            fine = self.curve(det, np.asarray(amplitudes[k].detach()))
            if k in self.pairs:
                kp = self.pairs[k]
                fine = fine + self.curve(self.L.parse_channel(kp)[1], np.asarray(amplitudes[kp].detach()))
            out[k] = self.R @ fine
        return out

    def check_against(self, vals, amplitudes=None):
        """The check that can fail: bff's forward model against the
        prototype's, channel by channel, as a relative difference of the peak.
        Two independent implementations of the same periodic reconvolution --
        one in C++, one in torch — have to agree, or one of them is wrong."""
        a2 = self.graph.amplitudes(vals) if amplitudes is None else amplitudes
        ref = self.graph.raw_counts(vals, a2)
        got = self(a2)
        return {k: float(np.abs(got[k] - np.asarray(ref[k].detach())).max()
                         / max(float(np.abs(np.asarray(ref[k].detach())).max()), 1e-300))
                for k in self.keys}


def structural_graph(graph, model):
    """The model's structure as a `bff.FactorGraph`: which variables exist,
    which factors touch which of them, and what that implies for how it can be
    decomposed.

    This is the graph in the graphical-model sense, and it is separate from the
    dataflow that computes the curves. It knows nothing about numbers; it knows
    that the donor-only fraction and the spline coefficients appear together in
    every likelihood factor, and therefore cannot be updated independently.
    """
    import IMP.bff as bff
    L = model['L']
    fg = bff.FactorGraph()
    idx = 0
    for name, (a, b) in sorted(graph.offsets.items(), key=lambda kv: kv[1][0]):
        fg.add_variable(name, name, idx)
        idx += 1
    names = set(graph.offsets)
    #: the priors: one factor per variable it constrains
    #: factor kinds are integers on this build: 0 a prior, 1 a likelihood
    #: (`get_number_of_likelihood_factors` counts the second)
    for name in sorted(names):
        fg.add_factor(f'prior_{name}', 0, [name])
    #: the likelihoods: each histogram sees the physics, the calibration, its
    #: own detector's response and its own scale, scatter and background
    shared = [n for n in ('c', 'x_d0', 'spec_eps', 'w_a', 'w_rho', 'w_rho_a', 'r0_d', 'r0_a',
                          'g', 'g_r', 'l1', 'l2', 'C_GD', 'C_GA', 'C_RD', 'C_RA',
                          'G_GREEN', 'G_RED', 'QY_D', 'QY_A', 'EX_AG', 'EX_DR', 'log10_lam')
              if n in names]
    for k in graph.data_keys:
        samp, det, _ = L.parse_channel(k)
        d0 = (model['Ep'].get('pulse_alias') or {}).get(det, det)
        scope = list(shared)
        for n in (f'irf_shift_{d0}', f'irf_width_{d0}', f'irf_skew_{d0}',
                  f'log_scale_{samp}', f'scat_{samp}_{k[1]}', f'bkg_{samp}_{k[1]}'):
            if n in names:
                scope.append(n)
        if k in (getattr(graph, 'pairs', {}) or {}):
            kp = graph.pairs[k]
            for n in (f'log_scale_{kp[0]}', f'scat_{kp[0]}_{kp[1]}'):
                if n in names:
                    scope.append(n)
        fg.add_factor(f'counts_{k[0]}_{k[1]}', 1, sorted(set(scope)))
    return fg


# --------------------------------------------------------------------------
# drawing the factor graph
# --------------------------------------------------------------------------

#: what each parameter is, for the drawing: the group it belongs to and a
#: label a reader who has never seen a TCSPC histogram can follow
GROUPS = {
    'distance': ('c', 'log10_lam', 'x_d0'),
    'photophysics': ('spec_eps', 'w_a', 'w_rho', 'w_rho_a', 'r0_d', 'r0_a'),
    'detection': ('g', 'g_r', 'l1', 'l2', 'G_GREEN', 'G_RED'),
    'crosstalk': ('C_GD', 'C_GA', 'C_RD', 'C_RA', 'EX_AG', 'EX_DR', 'QY_D', 'QY_A'),
    'response': None,          # irf_*
    'per sample': None,        # log_scale_*
    'per channel': None,       # scat_*, bkg_*
}

PRETTY = {
    'c': 'p(R/R0)\nspline coefficients', 'log10_lam': 'roughness\nweight',
    'x_d0': 'donor-only\nfraction', 'spec_eps': 'donor lifetime\nspectrum',
    'w_a': 'acceptor\nlifetimes', 'w_rho': 'donor rotation', 'w_rho_a': 'acceptor rotation',
    'r0_d': 'donor r0', 'r0_a': 'acceptor r0', 'g': 'g factor', 'g_r': 'g factor (red)',
    'l1': 'l1', 'l2': 'l2', 'G_GREEN': 'green\nefficiency', 'G_RED': 'red\nefficiency',
    'C_GD': 'donor into\ngreen', 'C_GA': 'acceptor into\ngreen', 'C_RD': 'donor into\nred',
    'C_RA': 'acceptor into\nred', 'EX_AG': 'green excites\nacceptor',
    'EX_DR': 'red excites\ndonor', 'QY_D': 'donor\nquantum yield', 'QY_A': 'acceptor\nquantum yield',
}


def group_of(name):
    for g, members in GROUPS.items():
        if members and name in members:
            return g
    if name.startswith('irf_'):
        return 'response'
    if name.startswith('log_scale_'):
        return 'per sample'
    if name.startswith('scat_') or name.startswith('bkg_'):
        return 'per channel'
    return 'other'


def factor_graph_nx(graph, model):
    """The model as a bipartite `networkx` graph: one node per parameter, one
    per factor, an edge wherever a factor's scope contains a parameter."""
    import networkx as nx
    L = model['L']
    G = nx.Graph()
    for name in graph.offsets:
        G.add_node(name, kind='variable', group=group_of(name),
                   size=int(graph.offsets[name][1] - graph.offsets[name][0]))
    names = set(graph.offsets)
    shared = [n for n in ('log10_lam', 'c', 'x_d0', 'spec_eps', 'w_a', 'w_rho', 'w_rho_a',
                          'r0_d', 'r0_a', 'g', 'g_r', 'l1', 'l2', 'C_GD', 'C_GA', 'C_RD', 'C_RA',
                          'G_GREEN', 'G_RED', 'QY_D', 'QY_A', 'EX_AG', 'EX_DR') if n in names]
    for k in graph.data_keys:
        samp, det, _ = L.parse_channel(k)
        d0 = (model['Ep'].get('pulse_alias') or {}).get(det, det)
        f = f'{samp} {k[1]}'
        G.add_node(f, kind='factor', group='likelihood')
        scope = list(shared)
        for n in (f'irf_shift_{d0}', f'irf_width_{d0}', f'irf_skew_{d0}',
                  f'log_scale_{samp}', f'scat_{samp}_{k[1]}', f'bkg_{samp}_{k[1]}'):
            if n in names:
                scope.append(n)
        if k in (getattr(graph, 'pairs', {}) or {}):
            kp = graph.pairs[k]
            for n in (f'log_scale_{kp[0]}', f'scat_{kp[0]}_{kp[1]}'):
                if n in names:
                    scope.append(n)
        for n in sorted(set(scope)):
            G.add_edge(f, n)
    return G


def plot_factor_graph(G, title='', figsize=(15, 9)):
    """Draw it so the structure is readable: the eight likelihood factors down
    the middle, the parameters that every one of them touches on the left, and
    the ones local to a detector or a channel on the right.

    The point of the picture is the left-hand column. Everything there appears
    in every likelihood, which is what a treewidth of 26 means in practice and
    why the degeneracies of notebook 5 are possible.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    colour = {'distance': '#c0392b', 'photophysics': '#e67e22', 'detection': '#16a085',
              'crosstalk': '#2980b9', 'response': '#8e44ad', 'per sample': '#7f8c8d',
              'per channel': '#95a5a6', 'other': '#bdc3c7', 'likelihood': '#2c3e50'}
    facs = [n for n, d in G.nodes(data=True) if d['kind'] == 'factor']
    vars_ = [n for n, d in G.nodes(data=True) if d['kind'] == 'variable']
    shared = [v for v in vars_ if G.degree(v) == len(facs)]
    local = [v for v in vars_ if G.degree(v) < len(facs)]
    order = ['distance', 'photophysics', 'detection', 'crosstalk']
    shared.sort(key=lambda v: (order.index(G.nodes[v]['group']) if G.nodes[v]['group'] in order else 9, v))
    local.sort(key=lambda v: (G.nodes[v]['group'], v))
    pos = {}
    for i, v in enumerate(shared):
        pos[v] = (0.0, 1.0 - 2.0 * i / max(len(shared) - 1, 1))
    for i, f in enumerate(facs):
        pos[f] = (1.0, 0.85 - 1.7 * i / max(len(facs) - 1, 1))
    for i, v in enumerate(local):
        pos[v] = (2.0, 1.0 - 2.0 * i / max(len(local) - 1, 1))
    fig, ax = plt.subplots(figsize=figsize)
    for u, v in G.edges():
        f, w = (u, v) if G.nodes[u]['kind'] == 'factor' else (v, u)
        c = colour[G.nodes[w]['group']]
        ax.plot(*zip(pos[f], pos[w]), color=c, lw=0.4, alpha=0.35, zorder=1)
    for v in shared + local:
        d = G.nodes[v]
        x, y = pos[v]
        ax.scatter([x], [y], s=90 + 26 * min(d['size'], 12), color=colour[d['group']],
                   edgecolor='white', linewidth=0.7, zorder=3)
        ax.text(x - 0.045 if x == 0 else x + 0.045, y, PRETTY.get(v, v.replace('_', ' ')),
                ha='right' if x == 0 else 'left', va='center', fontsize=6.6,
                color='0.15', zorder=4)
    for f in facs:
        x, y = pos[f]
        ax.scatter([x], [y], s=340, marker='s', color=colour['likelihood'], zorder=3)
        ax.text(x, y, f.replace(' ', '\n'), ha='center', va='center', fontsize=6,
                color='white', zorder=4)
    ax.set_xlim(-0.62, 2.62); ax.set_ylim(-1.15, 1.15); ax.axis('off')
    ax.text(0.0, 1.09, 'in every likelihood', ha='center', fontsize=9, weight='bold')
    ax.text(1.0, 1.09, 'the eight histograms', ha='center', fontsize=9, weight='bold')
    ax.text(2.0, 1.09, 'local to a detector or channel', ha='center', fontsize=9, weight='bold')
    seen, handles = [], []
    for g in ['distance', 'photophysics', 'detection', 'crosstalk', 'response', 'per sample', 'per channel']:
        if any(G.nodes[v]['group'] == g for v in vars_):
            handles.append(plt.Line2D([], [], marker='o', ls='', color=colour[g], label=g))
    ax.legend(handles=handles, loc='lower center', ncol=len(handles), frameon=False, fontsize=7.5,
              bbox_to_anchor=(0.5, -0.06))
    ax.set_title(title or 'the factor graph: every parameter, and which histograms see it', fontsize=11)
    fig.tight_layout()
    return fig


# --------------------------------------------------------------------------
# the objective, assembled from the fast parts
# --------------------------------------------------------------------------

class FastObjective:
    """The log posterior, evaluated the quick way.

    Three changes against the prototype's own `log_posterior`, none of which
    changes the value by more than rounding:

    1. **bff computes the curves.** The prototype builds a basis -- one column
       per lifetime per detector -- and multiplies. It needs the basis for the
       Jacobian, but an objective evaluation does not, and building it is most
       of the cost.
    2. **The spectrum goes into bff as one vector.** Writing 33 scalar
       amplitude ports costs more than the convolution they feed.
    3. **The amplitudes come from `stage1`/`stage2`.** The prototype already
       contains that path -- it is what the analytic Jacobian uses -- and it
       contracts the lifetime index of the big rotational map before the
       distance and rotation indices instead of after. Same numbers, a quarter
       of the time. The slower order in `physics_amplitudes` is not an
       oversight: it is faster under forward-mode differentiation with a
       hundred tangents, which is a different job.

    Use it for objective evaluations. The gradient and the Jacobian still come
    from the prototype's automatic differentiation.
    """

    def __init__(self, model, graph, theta):
        import numpy as _np
        self.model, self.graph = model, graph
        self.L, self.E = model['L'], model['Ep']
        self.forward = BffForward(model, graph, graph.unpack(theta)[0])
        self.scopes = sorted({self.L.channel_scope(k, self.E) for k in graph.keys})
        self.y = {k: _np.asarray(graph.data.y[k]) for k in graph.data_keys}
        self.mask = {k: _np.asarray(graph.data.mask[k]) for k in graph.data_keys}
        self.soft = float(self.L.SOFT)

    def _soft_positive(self, x):
        #: SOFT * softplus(x / SOFT), the prototype's positivity floor, in numpy
        z = x / self.soft
        return self.soft * np.where(z > 30.0, z, np.log1p(np.exp(np.minimum(z, 30.0))))

    def __call__(self, theta):
        L, E, g = self.L, self.E, self.graph
        vals, z = g.unpack(theta)
        st = L.stage1(E, vals, g.spectrum(vals), g.distribution(vals), self.scopes)
        a = L.stage2(E, vals, st, g.keys)
        a2 = {k: L.instrument_amplitudes(vals, a[k], k) for k in g.keys}
        self.forward.set_response(vals)
        lam = self.forward(a2)
        ll = 0.0
        for k in g.data_keys:
            m = np.maximum(self._soft_positive(lam[k]), 1e-300)
            ll += float(np.sum(self.mask[k] * (self.y[k] * np.log(m) - m)))
        return ll + float(g.log_prior(vals, z))
