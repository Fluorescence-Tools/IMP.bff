"""Photon sorting: what a single-molecule measurement carries that a pooled
decay does not, and where sorting stops working.

Notebooks 12 to 15 fit pooled histograms. That is right for a cuvette and it
throws away what a single-molecule experiment is for: the photons arrive in
bursts, each burst is one molecule, and a burst can be measured before its
photons are added to anyone else's.

Two rulers, not one
-------------------
The fraction of a burst's photons that came out red says how much energy
transferred. The mean arrival time of its donor photons says how much the donor
was quenched. For one molecule at one distance these are the same measurement
and the population sits on a line. **When they disagree the molecule was not at
one distance while it was being watched** -- and the size of the disagreement is
how the disagreement is read.

Both reference curves here are computed from the model's own per-molecule
tables, not from a single-exponential formula. The donor in this model has a
lifetime spectrum, so `E = 1 - tau/tau_0` is an approximation; taking the locus
straight out of the tables avoids the question of how good an approximation, and
`line_error()` measures it for the notebook that wants to quote it.

Exchange
--------
`exchange_stream` gives each molecule a two-state Markov jump process and emits
each photon from whichever state the molecule was in at that instant. At an
exchange time much longer than a burst it must reproduce the existing static
simulation, and `check_slow_limit()` is that gate.

Author: written for the bff examples, 2026-09-09.
See `okf/prd-single-molecule-sorting.md` in the ucfret repository for why.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bd import P                                            # noqa: E402

__all__ = ['two_state_truth', 'per_molecule', 'line_table', 'static_line', 'dynamic_line',
           'line_curvature', 'above_static', 'exchange_stream', 'check_slow_limit', 'observables',
           'line_slope', 'scatter_about_line', 'donor_only_from_s', 'plot_e_tau', 'plot_e_s', 'segment', 'group_photons', 'burst_state_fraction']

#: the four physics channels of the labelled sample under the green pulse, in
#: the order `physics_tables` builds them: green parallel, green perpendicular,
#: red parallel, red perpendicular
GREEN = (0, 1)
RED = (2, 3)


# --------------------------------------------------------------------------
# the truth: two states on the model's own grid
# --------------------------------------------------------------------------

def two_state_truth(model, r1: float, r2: float, f1: float = 0.5):
    """Two states at `r1` and `r2` in units of R0, with equilibrium fractions
    `f1` and `1 - f1`, placed on the nearest grid points.

    Both states are single grid points, not Gaussians: a state that a molecule
    occupies has one distance, and the width a fit recovers for it is then
    entirely the method's, which is what the sweeps here are measuring.
    """
    rel = model['rel']
    j1 = int(np.argmin(np.abs(rel - r1))); j2 = int(np.argmin(np.abs(rel - r2)))
    if j1 == j2:
        raise ValueError(f'r1={r1} and r2={r2} land on the same grid point '
                         f'({rel[j1]:.3f}); the grid cannot tell them apart')
    p = np.zeros_like(rel)
    p[j1] = f1; p[j2] = 1.0 - f1
    return p, (j1, j2)


# --------------------------------------------------------------------------
# the two axes, from the model's own tables
# --------------------------------------------------------------------------

def per_molecule(model, tab, j: int):
    """`(E*, tau)` for a molecule sitting at grid point `j` for the whole burst.

    `E*` is the proximity ratio -- of the photons the green pulse produced, the
    fraction that came out red. It is not corrected for crosstalk, direct
    excitation or detection efficiency, which is deliberate: it is what the
    detectors report, and every correction is a place to be wrong.

    `tau` is the intensity-weighted mean arrival time of the donor photons in
    nanoseconds, including the instrument's offset, because that is also what
    the detectors report.
    """
    w = tab['w_da'][j]
    n_fine = model['n_fine']; dt = model['dt']
    t = (np.arange(n_fine) + 0.5) * dt
    num = den = 0.0
    for d in GREEN:
        pdf = np.diff(np.concatenate(([0.0], tab['cdf_da'][j][d])))
        num += w[d] * float((t * pdf).sum()); den += w[d]
    tau = num / max(den, 1e-300)
    return float(w[list(RED)].sum() / max(w.sum(), 1e-300)), float(tau)


def line_table(model, st, js=None, stride=4, extra=(), verbose=False):
    """The per-molecule tables over MANY grid points, which the reference lines
    need and a simulation does not.

    `physics_tables` only fills the grid points a truth actually occupies --
    two, for a two-state sample -- because that is all the photons come from.
    The static line is a locus over the whole grid, so it gets its own table.
    Cached on the model dict, since it costs a few seconds and never changes.
    `extra` adds grid points that must be present whatever the stride -- the two
    states of a dynamic line, for instance, which the stride will usually miss.
    """
    key = ('_line_table', stride if js is None else tuple(js), tuple(sorted(extra)))
    if key in model:
        return model[key]
    rel = model['rel']
    if js is None:
        js = list(range(0, len(rel), int(stride)))
    js = sorted(set(list(js) + [int(e) for e in extra]))
    p = np.zeros_like(rel)
    p[list(js)] = 1.0
    p = p / p.sum()
    tab = P.physics_tables(model, st, p, verbose=verbose)
    model[key] = tab
    return tab


def static_line(model, tab, used=None):
    """The locus every molecule that does NOT move must lie on: `(tau, E*)`
    over the distance grid.

    This is not a formula -- it is the same table the simulation draws photons
    from, read out one grid point at a time, so it carries the multi-exponential
    donor, the instrument response and the crosstalk exactly as the data do.
    Pass the table from `line_table`, not the one the simulation uses: that one
    holds only the grid points the truth occupies.
    """
    js = list(tab['used']) if used is None else list(used)
    if len(js) < 3:
        raise ValueError(f'the static line needs the tables over many grid points and this '
                         f'one has {len(js)}; build it with line_table(model, st)')
    out = np.array([per_molecule(model, tab, int(j))[::-1] for j in js])   # (tau, E)
    o = np.argsort(out[:, 0])
    return out[o], np.asarray(js)[o]


def dynamic_line(model, tab, j1: int, j2: int, n: int = 101):
    """The locus of a molecule that interconverts between two fixed states
    FASTER than its burst, as the fraction of time in state 1 runs from 0 to 1.

    It lies above the static line, and the reason is worth stating: the
    efficiency averages the two states weighted by how many photons each
    produced, while the arrival time averages them weighted by photons TIMES
    lifetime. A long-lived state therefore pulls the clock further than it pulls
    the intensity, and the two rulers disagree by an amount that depends only on
    how far apart the states are.
    """
    n_fine = model['n_fine']; dt = model['dt']
    t = (np.arange(n_fine) + 0.5) * dt
    w1, w2 = tab['w_da'][j1], tab['w_da'][j2]
    m1 = np.array([float((t * np.diff(np.concatenate(([0.0], tab['cdf_da'][j1][d])))).sum())
                   for d in GREEN])
    m2 = np.array([float((t * np.diff(np.concatenate(([0.0], tab['cdf_da'][j2][d])))).sum())
                   for d in GREEN])
    out = []
    for f in np.linspace(0.0, 1.0, n):
        w = f * w1 + (1.0 - f) * w2
        num = float((f * w1[list(GREEN)] * m1 + (1.0 - f) * w2[list(GREEN)] * m2).sum())
        den = float((f * w1[list(GREEN)] + (1.0 - f) * w2[list(GREEN)]).sum())
        out.append((num / max(den, 1e-300), float(w[list(RED)].sum() / max(w.sum(), 1e-300))))
    return np.array(out)


def line_curvature(model, tab):
    """How far the static locus is from being a straight line.

    The textbook static FRET line is straight, because it assumes one donor
    lifetime, no instrument offset in the arrival time and a proximity ratio
    already corrected to a transfer efficiency. None of the three holds here:
    the donor has a lifetime spectrum, the mean arrival time carries the
    response's own position, and `E*` is what the detectors report. So the
    locus this model produces is a CURVE, and quoting how far it is from the
    straight line through its ends is the honest way to say by how much.

    Returns `(max deviation in E*, the two end points)`.
    """
    line, _ = static_line(model, tab)
    tau, E = line[:, 0], line[:, 1]
    chord = np.interp(tau, [tau[0], tau[-1]], [E[0], E[-1]])
    return float(np.max(np.abs(E - chord))), ((tau[0], E[0]), (tau[-1], E[-1]))


def above_static(model, tab, j1, j2, n=51):
    """How far the dynamic line rises above the static one, in `E*`.

    This is the signal the whole diagram is read for: a molecule that
    interconverted during its burst sits above the line a static one must lie
    on, and the gap is what says it moved. Zero would mean the diagram carries
    no information about exchange between these two states.
    """
    line, _ = static_line(model, tab)
    d = dynamic_line(model, tab, j1, j2, n)
    stat = np.interp(d[:, 0], line[:, 0], line[:, 1])
    return d[:, 0], d[:, 1] - stat


# --------------------------------------------------------------------------
# exchange
# --------------------------------------------------------------------------

def exchange_stream(model, st, tab, states, t_exchange_ms, rng=None,
                    n_bursts=None, sample='DA'):
    """Bursts of molecules that interconvert between two states.

    `states` is `((j1, f1), (j2, f2))` -- grid points and equilibrium
    fractions. `t_exchange_ms` is the relaxation time of the two-state process,
    `1 / (k12 + k21)`; the individual rates follow from detailed balance, so a
    state that is rarely occupied is also short-lived. `numpy.inf` means no
    exchange at all, which must reproduce the static simulation.

    A molecule's trajectory is a continuous-time Markov chain over the burst.
    Each dwell interval emits its own photons, at a rate proportional to that
    state's brightness, with micro times from that state's decay. Nothing about
    the photon physics changes -- only which state's table a photon comes from.
    """
    rng = rng or np.random.default_rng(int(st['seed']))
    n_bursts = int(st['n_bursts']) if n_bursts is None else int(n_bursts)
    (j1, f1), (j2, f2) = states
    tot_f = f1 + f2
    f1, f2 = f1 / tot_f, f2 / tot_f
    period_s = model['period'] * 1e-9
    dur = max(int(st['burst_duration_ms'] * 1e-3 / period_s), 1)
    ref = tab['w_d0'].sum()
    #: detailed balance: k12 f1 = k21 f2 and k12 + k21 = 1 / t_exchange
    if np.isfinite(t_exchange_ms) and t_exchange_ms > 0:
        k = 1.0 / (t_exchange_ms * 1e-3 / period_s)       # per pulse
        k12, k21 = k * f2, k * f1
    else:
        k12 = k21 = 0.0

    macro, micro, chan, flag, per_burst = [], [], [], [], []
    frac1 = []
    t = 0
    x_d0 = st['donor_only_fraction']
    for _ in range(n_bursts):
        t += int(rng.exponential(st['burst_spacing_ms'] * 1e-3 / period_s))
        donor_only = rng.random() < x_d0
        if donor_only:
            segs = [(-1, 0, dur)]
        else:
            #: the trajectory
            s = 0 if rng.random() < f1 else 1
            segs, u = [], 0
            while u < dur:
                rate = (k12 if s == 0 else k21)
                dwell = dur - u if rate <= 0 else rng.exponential(1.0 / rate)
                v = min(u + dwell, float(dur))
                segs.append((j1 if s == 0 else j2, u, v))
                u = v; s = 1 - s
        n_tot = 0; in1 = 0.0
        for j, u, v in segs:
            dtn = (v - u) / dur
            if j >= 0 and j == j1:
                in1 += dtn
            w = tab['w_d0'] if j < 0 else tab['w_da'][j]
            cdf = tab['cdf_d0'] if j < 0 else tab['cdf_da'][j]
            parts = [(w, cdf, rng.poisson(st['photons_per_burst'] * w.sum() / ref * dtn), 0)]
            if j >= 0:
                parts.append((tab['w_a0'], tab['cdf_a0'],
                              rng.poisson(st['photons_per_burst'] * st['red_green_ratio'] * dtn), 4))
            for w_, cdf_, n_, base in parts:
                if n_ <= 0:
                    continue
                c = rng.choice(len(w_), size=int(n_), p=w_ / w_.sum())
                for i in range(len(w_)):
                    m = int((c == i).sum())
                    if not m:
                        continue
                    micro.append(P._draw_micro(cdf_[i], m, rng))
                    chan.append(np.full(m, (base + i) % 4, dtype=np.int8))
                    macro.append(t + rng.integers(int(u), max(int(v), int(u) + 1), m, dtype=np.int64))
                    flag.append(np.zeros(m, dtype=bool)); n_tot += m
        n_s = rng.poisson(st['scatter_fraction'] * max(n_tot, 1))
        if n_s > 0:
            d = int(rng.integers(0, 2))
            det = tab['dets'][tab['keys']['DA'][d]]
            micro.append(P._draw_micro(tab['scatter'][det], n_s, rng))
            chan.append(np.full(n_s, d, dtype=np.int8))
            macro.append(t + rng.integers(0, dur, n_s, dtype=np.int64))
            flag.append(np.zeros(n_s, dtype=bool))
        per_burst.append((t, n_tot, -1 if donor_only else j1))
        frac1.append(np.nan if donor_only else in1)
        t += dur
    t_end = max(t, 1)
    n_bg = int(st['background_cps'] * 4 * t_end * period_s)
    if n_bg > 0:
        macro.append(rng.integers(0, t_end, n_bg, dtype=np.int64))
        micro.append(rng.integers(0, model['n_fine'], n_bg).astype(np.uint16))
        chan.append(rng.integers(0, 4, n_bg).astype(np.int8))
        flag.append(np.ones(n_bg, dtype=bool))
    macro = np.concatenate(macro); micro = np.concatenate(micro)
    chan = np.concatenate(chan); flag = np.concatenate(flag)
    o = np.argsort(macro, kind='stable')
    return dict(macro=macro[o].astype(np.uint64), micro=micro[o], chan=chan[o],
                background=flag[o], bursts=per_burst, duration_s=t_end * period_s,
                n_background=n_bg, fraction_state1=np.array(frac1))


def burst_state_fraction(stream):
    """the fraction of each burst spent in state 1 (NaN for donor-only)"""
    return stream['fraction_state1']


def check_slow_limit(model, st, tab, states, n_bursts=1500, seed=7, verbose=True):
    """M3's gate: with no exchange, the new simulation must agree with the one
    already in the example.

    The comparison is a two-sample Kolmogorov-Smirnov test on the per-burst
    proximity ratio and on the mean donor arrival time. A p-value is reported
    rather than a pass mark on a tolerance nobody chose.
    """
    from scipy.stats import ks_2samp
    (j1, f1), (j2, f2) = states
    p = np.zeros_like(model['rel']); p[j1] = f1; p[j2] = f2; p = p / p.sum()
    a = P._sample_bursts(model, st, tab, np.random.default_rng(seed), n_bursts, p, 'DA')
    b = exchange_stream(model, st, tab, states, np.inf,
                        np.random.default_rng(seed + 1), n_bursts=n_bursts)
    out = {}
    for name, s in (('existing simulation', a), ('exchange code, no exchange', b)):
        tt = P.to_tttr(model, s)
        _, keep = P.burst_search(model, tt, st)
        out[name] = P.burst_observables(tt, keep, model)
    A, B = out['existing simulation'], out['exchange code, no exchange']
    res = {}
    for i, nm in ((0, 'proximity ratio'), (3, 'mean donor arrival time')):
        k = ks_2samp(A[:, i], B[:, i])
        res[nm] = dict(statistic=float(k.statistic), p_value=float(k.pvalue))
        if verbose:
            print(f'  {nm:<26} KS {k.statistic:.4f}  p = {k.pvalue:.3f}  '
                  f'({len(A)} vs {len(B)} bursts)')
    if verbose:
        print('  a small p would say the two simulations differ; they should not')
    return res


# --------------------------------------------------------------------------
# reading the bursts
# --------------------------------------------------------------------------

#: the columns of the feature table this module builds. The first four are the
#: ones `pie_mfd.burst_observables` produces, in the same order, so anything
#: written against that table still works; the rest are what the scatter check
#: needs and that table does not carry.
COLUMNS = ('E*', 'S', 'size', 'tau', 'n_donor', 'n_sensitised', 'n_direct', 'tau_sd')


def observables(model, st, stream):
    """`(tttr, selection, bursts, features)` -- the photon object, the indices
    of the burst photons, the burst boundaries, and the per-burst table.

    The columns are `COLUMNS`. The first four reproduce
    `pie_mfd.burst_observables` exactly; the last four are the photon counts the
    two ratios were formed from and the spread of the donor arrival times within
    the burst, which is what turns "the bursts scatter about the line" into a
    statement with a prediction attached.
    """
    tt = P.to_tttr(model, stream)
    sel, keep = P.burst_search(model, tt, st)
    mi = np.asarray(tt.micro_times); ch = np.asarray(tt.routing_channels)
    split = int(round(model['t_pie'] / model['dt']))
    green_pulse = mi < split
    dt = model['dt']
    rows = []
    for a, b in keep:
        sl = slice(a, b + 1)
        g, r, gp = ch[sl] < 2, ch[sl] >= 2, green_pulse[sl]
        n_gg = int((g & gp).sum()); n_gr = int((r & gp).sum()); n_rr = int((r & ~gp).sum())
        tot = n_gg + n_gr + n_rr
        if tot < 1:
            continue
        dm = mi[sl][g & gp].astype(float) * dt
        tau = float(dm.mean()) if len(dm) else np.nan
        tsd = float(dm.std()) if len(dm) > 1 else np.nan
        rows.append((n_gr / max(n_gg + n_gr, 1), (n_gg + n_gr) / tot, tot, tau,
                     n_gg, n_gr, n_rr, tsd))
    return tt, sel, keep, np.array(rows)


def segment(feat, columns=(0, 3), min_cluster_size=None, method='leaf'):
    """HDBSCAN on the burst observables.

    The default reads the two INDEPENDENT rulers -- the proximity ratio and the
    mean donor arrival time -- and nothing else, so a group is a group because
    the intensity and the clock agree about it, not because the bursts happened
    to be bright.
    """
    X = np.asarray(feat)[:, list(columns)]
    X = (X - X.mean(0)) / np.maximum(X.std(0), 1e-12)
    n = max(int(0.02 * len(X)), 20) if min_cluster_size is None else int(min_cluster_size)
    return P.hdbscan_labels(X, min_cluster_size=n, method=method)


def group_photons(bursts, labels, group):
    """the photon indices of one group's bursts"""
    return P.group_indices(bursts, labels, group)


def line_slope(line, tau, window=0.4):
    """`dE*/dtau` of the static line at each of `tau`, by a local straight-line
    fit over a window of that width in nanoseconds.

    A finite difference between neighbouring points on the line is useless
    here: the line is sampled at distance-grid points, and near its ends the
    efficiency changes while the arrival time barely does, so the difference
    quotient blows up. Fitting locally instead gives the slope that actually
    projects an error in the clock onto the efficiency axis.
    """
    t, e = np.asarray(line[:, 0]), np.asarray(line[:, 1])
    out = np.zeros_like(np.asarray(tau, float))
    for i, x in enumerate(np.atleast_1d(tau)):
        m = np.abs(t - x) <= window
        if m.sum() < 3:
            k = np.argsort(np.abs(t - x))[:3]
            m = np.zeros_like(t, bool); m[k] = True
        A = np.vstack([t[m] - x, np.ones(m.sum())]).T
        out[i] = float(np.linalg.lstsq(A, e[m], rcond=None)[0][0])
    return out


def scatter_about_line(feat, line):
    """How far the bursts sit from the static line, against what counting
    statistics alone predict.

    A static sample must sit ON the line, and the spread around it is then
    arithmetic rather than physics. **Two terms, and the second is the larger
    one**, which is the trap this function exists to avoid: the proximity ratio
    of a burst of `n` green-pulse photons is uncertain by `sqrt(E(1-E)/n)`; but
    the burst's mean arrival time is uncertain too, by the spread of its donor
    arrival times over the root of their number, and the line is steep, so that
    uncertainty appears as a vertical residual as well. Leaving the second term
    out makes a perfectly ordinary static sample look four times noisier than it
    should be -- which is what the first version of this check reported.
    """
    f = np.asarray(feat)
    tau, E = f[:, 3], f[:, 0]
    n_gr = f[:, 4] + f[:, 5]                 # the photons E* is formed from
    n_d, tsd = f[:, 4], f[:, 7]
    ok = np.isfinite(tau) & np.isfinite(E) & (n_gr > 0) & np.isfinite(tsd) & (n_d > 1)
    tau, E, n_gr, n_d, tsd = tau[ok], E[ok], n_gr[ok], n_d[ok], tsd[ok]
    on = np.interp(tau, line[:, 0], line[:, 1])
    resid = E - on
    sl = line_slope(line, tau)
    var_E = np.clip(on * (1.0 - on), 0, None) / np.maximum(n_gr, 1)
    var_tau = tsd ** 2 / np.maximum(n_d, 1)
    expect = np.sqrt(var_E + sl ** 2 * var_tau)
    pred = float(np.sqrt((expect ** 2).mean()))
    return dict(residual_sd=float(resid.std()), predicted_sd=pred,
                from_counting_E=float(np.sqrt(var_E.mean())),
                from_the_clock=float(np.sqrt((sl ** 2 * var_tau).mean())),
                ratio=float(resid.std() / max(pred, 1e-30)),
                median_burst_size=float(np.median(f[ok, 2])), n_bursts=int(len(resid)),
                residual=resid)


def donor_only_from_s(feat, threshold=0.85):
    """The fraction of bursts whose stoichiometry says the acceptor is not
    there. This is the quantity the second laser measures directly and the
    donor-excitation-only geometry has to infer."""
    f = np.asarray(feat)
    ok = np.isfinite(f[:, 1])
    return float((f[ok, 1] > threshold).mean()), int(ok.sum())


# --------------------------------------------------------------------------
# the diagrams
# --------------------------------------------------------------------------

def plot_e_tau(feat, model=None, tab=None, states=None, ax=None, bins=70,
               title='', tau_range=None, e_range=(-0.05, 1.05)):
    """`E*` against the mean donor arrival time, as a two-dimensional histogram,
    with the static line every non-moving molecule must lie on and, if two
    states are named, the line a fast-exchanging molecule lies on instead."""
    import matplotlib.pyplot as plt
    if ax is None:
        _, ax = plt.subplots(figsize=(5.2, 4.2))
    f = np.asarray(feat)
    tau, E = f[:, 3], f[:, 0]
    ok = np.isfinite(tau) & np.isfinite(E)
    tr = tau_range or (float(np.nanpercentile(tau[ok], 0.5)) - 0.2,
                       float(np.nanpercentile(tau[ok], 99.5)) + 0.2)
    ax.hist2d(tau[ok], E[ok], bins=bins, range=[tr, e_range], cmap='viridis',
              cmin=1)
    if model is not None and tab is not None:
        line, _ = static_line(model, tab)
        ax.plot(line[:, 0], line[:, 1], 'w-', lw=2.0)
        ax.plot(line[:, 0], line[:, 1], 'k-', lw=1.0, label='static line')
        if states is not None:
            d = dynamic_line(model, tab, states[0][0], states[1][0])
            ax.plot(d[:, 0], d[:, 1], 'w--', lw=2.0)
            ax.plot(d[:, 0], d[:, 1], 'C3--', lw=1.2, label='dynamic line')
    ax.set_xlim(*tr); ax.set_ylim(*e_range)
    ax.set_xlabel('mean donor arrival time / ns')
    ax.set_ylabel('proximity ratio $E^*$')
    ax.legend(fontsize=7, loc='lower left', frameon=True)
    if title:
        ax.set_title(title, fontsize=9)
    return ax


def plot_e_s(feat, ax=None, bins=70, title='', s_range=(-0.05, 1.05)):
    """stoichiometry against `E*`: what the second laser buys, as a picture."""
    import matplotlib.pyplot as plt
    if ax is None:
        _, ax = plt.subplots(figsize=(5.2, 4.2))
    f = np.asarray(feat)
    ax.hist2d(f[:, 0], f[:, 1], bins=bins, range=[(-0.05, 1.05), s_range],
              cmap='viridis', cmin=1)
    ax.set_xlabel('proximity ratio $E^*$'); ax.set_ylabel('stoichiometry $S$')
    if title:
        ax.set_title(title, fontsize=9)
    return ax
