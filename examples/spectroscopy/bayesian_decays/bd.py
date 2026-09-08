"""Shared helpers for the `bayesian_decays` notebook series.

The series develops, one notebook at a time, the model that ended up in
`../smfret_pie_mfd`: what the measurement is, why an unregularised inversion
of it does not work, what prior fixes that and what the prior's weight does,
how the posterior is computed without a sampler, whether it is calibrated,
what the second laser buys, and what changes when the photons arrive one
molecule at a time.

`THEORY.md` beside this file carries the probability calculus; every notebook
points into it by section.

The physics and the posterior are the ucfret prototype
`s88_laplace_posterior` (investigation `pinn_pR_anisotropy`), imported
through `pie_mfd` in the sibling example folder so that there is one
implementation and not two.  Point `UCFRET_S88` at the prototype directory if
it is not at `~/dev/ucfret/investigation/pinn_pR_anisotropy`.

Author: written for the bff examples, 2026-09-08.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
_SIB = HERE.parent / 'smfret_pie_mfd'
if str(_SIB) not in sys.path:
    sys.path.insert(0, str(_SIB))

import pie_mfd as P                                     # noqa: E402
from pie_mfd import (default_settings, build_model, truth_distribution, physics_tables,   # noqa: F401,E402
                     simulate_stream, to_tttr, burst_search, burst_decays, assemble_data,
                     burst_observables, burst_truth, background_fraction, poisson_reference,
                     summary_table, node_table, plot_decays, plot_es, plot_pR, plot_groups,
                     classify_bursts, group_indices, group_summary, label_groups,
                     hdbscan_labels, check_hdbscan_selection, constant_values, fit)

__all__ = ['P', 'default_settings', 'build_model', 'truth_distribution', 'physics_tables',
           'simulate_stream', 'to_tttr', 'burst_search', 'burst_decays', 'assemble_data',
           'burst_observables', 'burst_truth', 'background_fraction', 'poisson_reference',
           'summary_table', 'node_table', 'plot_decays', 'plot_es', 'plot_pR', 'plot_groups',
           'classify_bursts', 'group_indices', 'group_summary', 'label_groups',
           'hdbscan_labels', 'check_hdbscan_selection', 'constant_values', 'fit',
           'gaussian_truth', 'expected_counts', 'simulate_ensemble', 'fit_at_lambda',
           'fit_ensemble', 'window', 'figure', 'mix_nodes', 'joint_mode', 'keep_inline', 'assert_inline']


# --------------------------------------------------------------------------
# truths and ensemble data (no bursts: the idealised measurement)
# --------------------------------------------------------------------------

def gaussian_truth(model, centres, weights=None, widths=0.06):
    """A distance distribution as a mixture of Gaussians in R/R0, on the
    model's own grid, so that the truth and what is compared against it are
    the same array."""
    rel = model['rel']
    centres = np.atleast_1d(np.asarray(centres, float))
    w = np.ones_like(centres) if weights is None else np.asarray(weights, float)
    sd = np.broadcast_to(np.asarray(widths, float), centres.shape)
    p = np.zeros_like(rel)
    for c, wi, s in zip(centres, w, sd):
        p += wi * np.exp(-0.5 * ((rel - c) / max(s, 1e-6)) ** 2)
    return p / p.sum()


def expected_counts(model, st, p_true, x_d0, targets, bkg=0.01, scat=0.02):
    """The expected counts of the twelve PIE channels for one truth.

    `targets`: the counts wanted in each sample's reference channel, which is
    how a measurement's length enters. Returns (expected counts per channel,
    the value dict, the graph they came from).
    """
    import torch
    L = model['L']; Ep = model['Ep']; keys = model['keys']
    vals = dict(P.constant_values(model, st))
    vals['x_d0'] = L.tt(float(x_d0))
    raw = L.physics_amplitudes(Ep, vals, Ep['cD'], L.tt(p_true), keys)
    for samp, ref in (('D0', ('D0', 'gv_vv')), ('DA', ('DA', 'gv_vv')), ('A0', ('A0', 'r2v_vv'))):
        vals[f'log_scale_{samp}'] = L.tt(math.log(targets[samp] / float(raw[ref].sum())))
    for k in keys:
        kk = f'{k[0]}_{k[1]}'; det = L.parse_channel(k)[1]
        vals[f'scat_{kk}'] = L.tt(scat if (det.startswith('g') and '2' not in det) else 1e-6)
        vals[f'bkg_{kk}'] = L.tt(1e-9 if k[0] == 'A0' else max(float(bkg), 1e-9))
    g0 = L.FactorGraph(Ep, keys, L.default_variables(Ep, keys, n_coef=model['spl'].shape[1], irf_shape=False),
                       L.PoissonCountsFactor({k: torch.ones(model['n_bin']) for k in keys}),
                       L.InstrumentModel(Ep, 'analytic'), model['spl'],
                       L.PSplineFactor(model['spl'].shape[1], spl=model['spl']))
    g0.rel = model['rel']
    a2 = {k: L.instrument_amplitudes(vals, raw[k], k) for k in keys}
    g0.assert_physical(vals, a2, where='truth')
    return g0.expected_counts(vals, a2), vals, g0


def simulate_ensemble(model, st, p_true, x_d0=0.2, photons=3e5, seed=0, bkg=0.01, scat=0.02):
    """One Poisson realisation of the eight histograms: no molecules, no
    bursts, no selection -- the data are exactly a realisation of the model
    being fitted. This is the idealised measurement the series works with
    until notebook 07."""
    import torch
    L = model['L']
    targets = {'D0': float(photons), 'DA': float(photons), 'A0': float(photons) * 0.6}
    lam, vals, g0 = expected_counts(model, st, p_true, x_d0, targets, bkg, scat)
    gen = torch.Generator().manual_seed(int(seed))
    y, lam_h = L.simulate_pie(lam, gen)
    return y, lam_h, vals


# --------------------------------------------------------------------------
# fitting
# --------------------------------------------------------------------------

def _route_b_graph(model, y):
    """The graph the whole series fits: twelve physics channels, eight
    histograms, each labelled histogram carrying both pulses."""
    import torch
    L = model['L']; Ep = dict(model['Ep']); keys = model['keys']; spl = model['spl']
    kind = {'gv': 'vv', 'gh': 'vh', 'rv': 'vv', 'rh': 'vh'}
    pairs = {('DA', f'{d}_{kind[d]}'): ('A0', f'{L.PIE_PARTNER[d]}_{kind[d]}') for d in L.PIE_DETS}
    Ep['pie_pairs'] = pairs
    win = {k: (torch.ones(model['n_bin']) if k[0] == 'D0' else
               (Ep['win_green'] if k[0] == 'DA' else Ep['win_red'])) for k in keys}
    data_keys = [k for k in keys if k[0] != 'A0']
    med = {'D0': max(float(y[('D0', 'gv_vv')].sum()), 1.0),
           'DA': max(float(y[('DA', 'gv_vv')].sum()), 1.0),
           'A0': max(float((y[('DA', 'rv_vv')] * Ep['win_red']).sum()), 1.0)}
    n_coef = spl.shape[1]
    g = L.FactorGraph(Ep, keys,
                      L.default_variables(Ep, keys, n_coef=n_coef, scale_medians=med, irf_shape=False),
                      L.PoissonCountsFactor({k: y[k] for k in data_keys},
                                            mask={k: torch.ones(model['n_bin']) for k in data_keys}),
                      L.InstrumentModel(Ep, 'analytic'), spl, L.PSplineFactor(n_coef, spl=spl))
    g.rel = model['rel']
    g = g.with_fixed(**{f'bkg_{kp[0]}_{kp[1]}': L.tt([math.log(1e-9)]) for kp in pairs.values()})
    g.start_y = {k: y[k] * win[k] for k in keys}
    return g


def fit_at_lambda(model, y, log10_lam, verbose=False, hessian='fisher'):
    """One Laplace approximation at one value of the roughness weight.

    Returns the node dict the prototype produces: the mode, the covariance,
    the log evidence, the expected counts, and Rule 0's numbers.
    """
    L = model['L']
    g = _route_b_graph(model, y)
    tr = g.index['log10_lam'].transform
    gi = g.with_fixed(log10_lam=tr.to_unconstrained(L.tt([float(log10_lam)])))
    th0, _ = L.start_from_data(gi, y, verbose=verbose)
    r = L.laplace_at(gi, th0, optimiser='fisher', verbose=verbose, hessian=hessian)
    vals, _ = gi.unpack(r['theta'])
    lam = gi.expected_counts(vals)
    rows, dev, dof = L.rule0(gi, y, lam)
    r.update(graph=gi, log10_lam=float(log10_lam), dev=dev, dof=dof, rows=rows,
             lam={k: lam[k].detach() for k in gi.data_keys}, p=gi.distribution(vals).detach().numpy())
    return r


def fit_ensemble(model, y, lam_nodes=(1.0, 0.0, -1.0), seed=0, verbose=False, refine=0.0):
    """The whole procedure: a Laplace at each node of the roughness grid, the
    nodes mixed by their evidences."""
    import torch
    L = model['L']
    keep = L.LAM_REFINE_STEP
    L.LAM_REFINE_STEP = float(refine)
    try:
        g = _route_b_graph(model, y)
        gen = torch.Generator().manual_seed(int(seed))
        post = L.fit_sample(g, y, gen, model['rel'], verbose=verbose, start='mem', optimiser='fisher',
                            lam_nodes=tuple(float(x) for x in lam_nodes), hessian='fisher')
    finally:
        L.LAM_REFINE_STEP = keep
    return g, post


# --------------------------------------------------------------------------
# small utilities the notebooks share
# --------------------------------------------------------------------------

def window(dE=0.01):
    """The R/R0 interval outside which the transfer efficiency is within dE of
    zero or one, and the distance therefore unbounded (THEORY.md 7.3)."""
    return ((dE / (1 - dE)) ** (1 / 6), ((1 - dE) / dE) ** (1 / 6))


def figure(name, fig=None, caption=None):
    """Save a figure beside the notebooks with its caption, the way the
    ucfret investigation does it: one file per result, named after the step
    that made it, and never overwritten silently."""
    import matplotlib.pyplot as plt
    from datetime import datetime
    d = HERE / 'figures'
    d.mkdir(exist_ok=True)
    stamp = datetime.now().strftime('%Y%m%dT%H%M%S')
    p = d / f'{stamp}_{name}.png'
    (fig or plt.gcf()).savefig(p, dpi=150, bbox_inches='tight')
    if caption:
        p.with_suffix('.txt').write_text(caption)
    return p


def mix_nodes(nodes, prior_slope=0.0):
    """Assemble the evidence-weighted mixture from a list of Laplace nodes.

    The posterior of everything except the roughness weight is
    sum_l w_l N(theta_l, Sigma_l) with w_l proportional to the evidence at
    node l (THEORY.md 5.3). This packs a list of `fit_at_lambda` results into
    the dictionary the prototype's plotting and summary functions expect.
    """
    ev = np.array([n['evidence'] for n in nodes], float)
    lg = np.array([n['log10_lam'] for n in nodes], float)
    t = ev + prior_slope * lg
    w = np.exp(t - np.nanmax(t[np.isfinite(t)]))
    w = np.where(np.isfinite(t), w, 0.0); w = w / w.sum()
    best = int(np.nanargmax(np.where(np.isfinite(t), t, -np.inf)))
    b = nodes[best]
    d = {float(n['log10_lam']): n for n in nodes}
    post = dict(nodes=d, lam_nodes=[float(x) for x in lg], weights=w, best_lam=float(lg[best]),
                ev=ev, theta=b['theta'], Sigma=b['Sigma'], H=b['H'], graph=b['graph'],
                lam=b['lam'], rows=b['rows'], dev=b['dev'], dof=b['dof'],
                evidence_weights=w.copy(), evidence_best=float(lg[best]), lam_prior_slope=prior_slope)
    return post


def joint_mode(model, y, verbose=False):
    """The mode of the JOINT posterior over everything including the roughness
    weight -- what one gets by treating a variance component as a coordinate.
    THEORY.md 5.2 says why this is wrong; this returns it so it can be shown."""
    L = model['L']
    g = _route_b_graph(model, y)
    th0, _ = L.start_from_data(g, y, verbose=verbose)
    r = L.laplace_at(g, th0, optimiser='fisher', verbose=verbose, hessian='fisher')
    vals, _ = g.unpack(r['theta'])
    lam = g.expected_counts(vals)
    rows, dev, dof = L.rule0(g, y, lam)
    r.update(graph=g, dev=dev, dof=dof, rows=rows, p=g.distribution(vals).detach().numpy(),
             log10_lam=float(vals['log10_lam']), lam={k: lam[k].detach() for k in g.data_keys})
    return r


def keep_inline(verbose=False):
    """Stop the figures from silently disappearing.

    Some of the prototype's modules call `matplotlib.use('Agg')` when they are
    imported -- they were written for headless scripts. In a notebook that is
    a trap: after the switch, a cell that draws a figure produces NOTHING, and
    says nothing about it either. Three figures went missing that way before
    this was found.

    Called once in the first cell, this restores the inline backend and makes
    later switches to a file backend no-ops for the rest of the session.
    `assert_inline` is the check to run after anything heavy has been
    imported.
    """
    try:
        from IPython import get_ipython
        ip = get_ipython()
        if ip is None:
            return False
        import matplotlib
        ip.run_line_magic('matplotlib', 'inline')
        if not getattr(matplotlib, '_bd_use_patched', False):
            real_use = matplotlib.use

            def _use(backend, *a, **k):
                if str(backend).lower() in ('agg', 'pdf', 'ps', 'svg', 'cairo', 'template'):
                    return
                return real_use(backend, *a, **k)
            matplotlib.use = _use
            matplotlib._bd_use_patched = True
        if verbose:
            print('backend:', matplotlib.get_backend())
        return True
    except Exception:
        return False


def assert_inline():
    """A check that can fail: is this session still going to draw anything?"""
    import matplotlib
    b = matplotlib.get_backend()
    assert 'inline' in b.lower(), f'the figure backend is {b}; figures would be drawn to nothing'
    return b


def fisher_sd(model, graph, theta, mask=None):
    """The posterior standard deviations a DESIGN implies, before any data are
    fitted: the expected information plus the prior curvature, inverted.

        I(theta) = J^T diag(w / lambda) J,   J = d lambda / d theta,

    with w a per-bin weight of 1 or 0 that says which bins the design
    measures. Changing `mask` and nothing else asks what a different
    experiment would have been able to say. Returns {name: standard deviation}
    in the model's own coordinates plus the whole matrix.
    """
    import torch
    L = model['L']
    dk = list(graph.data_keys)

    def mean_vec(t):
        v, _ = graph.unpack(t)
        lam = graph.expected_counts(v)
        return torch.cat([lam[k] for k in dk])

    J = torch.func.jacfwd(mean_vec)(theta)
    m = mean_vec(theta).clamp_min(1e-12)
    if mask is None:
        w = torch.cat([graph.data.mask[k] for k in dk])
    else:
        w = torch.cat([mask[k] for k in dk])
    F = J.T @ (J * (w / m)[:, None])
    Hp = L.chunked_hessian(lambda t: -graph.log_prior(*graph.unpack(t)), theta)
    H = F + Hp
    H = 0.5 * (H + H.T)
    Sig = torch.linalg.inv(H + 1e-10 * torch.eye(H.shape[0], dtype=H.dtype))
    sd = torch.sqrt(torch.diagonal(Sig).clamp_min(0))
    out = {}
    for n, (a, b) in graph.offsets.items():
        out[n] = float(sd[a]) if b - a == 1 else float(sd[a:b].mean())
    return out, H, Sig
