"""A single-molecule FRET experiment with pulsed interleaved excitation and
multi-parameter detection: simulated as a photon stream, searched for bursts
with tttrlib, and analysed for the distance distribution p(R/R0).

WHAT THE EXPERIMENT IS.  Molecules carrying a donor and an acceptor dye
diffuse through a confocal spot.  While one is inside, photons arrive in a
burst.  Two lasers alternate within one 50 ns period: the green pulse excites
the donor (and, a little, the acceptor directly), and 25 ns later the red
pulse excites the acceptor alone -- this is PIE, pulsed interleaved
excitation, and it is what tells a molecule that has lost its acceptor from
one that simply transfers nothing.  Four detectors record the two colours and
the two polarisations.  Every photon carries a MACRO time (which laser pulse
it belongs to, so: when the burst happened) and a MICRO time (its delay after
that pulse, a few hundred picoseconds of resolution).  That pair is what a
TCSPC card writes and what `tttrlib.TTTR` holds.

WHAT IS INFERRED.  The micro times of the burst photons form eight decay
histograms -- four detectors for a donor-only reference sample, four for the
labelled sample, each of the latter holding BOTH pulses in one histogram.
From those, the posterior of the donor-acceptor distance distribution
p(R/R0), of the fraction of molecules that have no acceptor, of the
anisotropies and of the calibration constants, all at once and analytically:
a Laplace approximation at each node of a grid over the roughness penalty
weight, mixed by their evidences.  There is no sampler anywhere.

THE PHYSICS IS THE PROTOTYPE'S.  Both the simulation and the fit call the
same forward model (`s88_laplace_posterior` in the ucfret investigation), so
this example measures the inference, not a disagreement between two
independently written physics implementations.  What the simulation adds on
top, and the fit does not know, is the experiment: molecules drawn one at a
time, Poisson burst sizes, uncorrelated background, scattered light, and a
burst search that decides which photons are used at all.

Point `UCFRET_S88` at the prototype directory if it is not at
`~/dev/ucfret/investigation/pinn_pR_anisotropy`.

Author: written for the bff examples, 2026-09-08.
"""
from __future__ import annotations

import math
import os
import sys
from pathlib import Path

import numpy as np

__all__ = ["default_settings", "load_prototype", "build_model", "truth_distribution",
           "physics_tables", "simulate_stream", "to_tttr", "burst_search",
           "burst_decays", "fit", "rule0_table"]


# --------------------------------------------------------------------------
# 0. settings -- everything a user of the notebook edits
# --------------------------------------------------------------------------

def default_settings() -> dict:
    """One dictionary, so the notebook's first cell is the only thing to edit."""
    return dict(
        # -- the sample: a simple mixture of FRET molecules ------------------
        R0=52.0,                      # Foerster radius, Angstrom
        species=[dict(R=42.0, fraction=0.5, width=3.0, name='high FRET'),
                 dict(R=62.0, fraction=0.5, width=4.0, name='low FRET')],
        donor_only_fraction=0.20,     # molecules whose acceptor is missing or bleached
        rho_donor=1.0,                # rotational correlation times, ns
        rho_acceptor=1.0,
        r0_donor=0.36,                # fundamental anisotropies
        r0_acceptor=0.18,
        # -- the experiment --------------------------------------------------
        photons_per_burst=110,        # mean detected photons of an unquenched donor, green pulse
        n_bursts=6000,                # labelled sample
        n_bursts_donor=2500,          # donor-only reference sample
        red_green_ratio=0.55,         # red-pulse photons per green-pulse photon of a labelled molecule
        burst_duration_ms=1.2,
        burst_spacing_ms=25.0,        # mean time between molecules entering the spot
        background_cps=350.0,         # per detector, uncorrelated
        scatter_fraction=0.02,        # of the burst photons, at the response itself
        # -- the burst search -------------------------------------------------
        search_window_ms=0.5,         # sliding time window
        search_min_photons=12,        # photons in that window for it to be a burst
        burst_min_photons=40,         # photons a burst must have to be kept
        burst_max_duration_ms=6.0,
        # -- the analysis -----------------------------------------------------
        knots=25,                     # spline coefficients of log p(R/R0)
        lam_nodes=(2.0, 1.0, 0.0, -1.0),   # log10 roughness weights the evidence mixes over
        threads=4,
        seed=2026,
    )


# --------------------------------------------------------------------------
# 1. the prototype
# --------------------------------------------------------------------------

def prototype_dir() -> Path:
    cand = []
    if os.environ.get('UCFRET_S88'):
        cand.append(Path(os.environ['UCFRET_S88']))
    cand.append(Path.home() / 'dev' / 'ucfret' / 'investigation' / 'pinn_pR_anisotropy')
    here = Path(__file__).resolve()
    cand += [p / 'ucfret' / 'investigation' / 'pinn_pR_anisotropy' for p in here.parents]
    for c in cand:
        if (c / 's88_laplace_posterior.py').exists():
            return c
    raise RuntimeError('the ucfret prototype s88_laplace_posterior.py was not found; '
                       'set UCFRET_S88 to the directory that holds it')


def load_prototype(threads: int = 4):
    """Import the prototype (and the physics maps it caches) once."""
    os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')
    d = prototype_dir()
    if str(d) not in sys.path:
        sys.path.insert(0, str(d))
    cwd = Path.cwd()
    os.chdir(d)                      # the maps are cached beside the prototype
    try:
        import torch
        torch.set_num_threads(int(threads))
        import s88_laplace_posterior as L
    finally:
        os.chdir(cwd)
    return L


def build_model(st: dict) -> dict:
    """The environment: the distance grid, the spline basis, the interleaved
    time axis, and the twelve physics channels of an MFD PIE measurement."""
    L = load_prototype(st.get('threads', 4))
    d = prototype_dir()
    cwd = Path.cwd(); os.chdir(d)
    try:
        E, rel, edges, spl, dx = L.environment(st['knots'], instrument=0, homogeneous=True)
    finally:
        os.chdir(cwd)
    Ep = L.pie_environment(E)
    keys = L.pie_keys()
    return dict(L=L, E=E, Ep=Ep, rel=np.asarray(rel), spl=spl, edges=edges, keys=keys,
                dt=float(E['dt']), n_fine=int(E['n']), n_bin=int(Ep['n_bin']),
                period=float(E['period']), t_pie=float(Ep['t_pie']))


# --------------------------------------------------------------------------
# 2. the truth: a mixture of FRET molecules on the model's distance grid
# --------------------------------------------------------------------------

def truth_distribution(model: dict, st: dict) -> np.ndarray:
    """p(R/R0) on the grid, from the species table.  Each species is a Gaussian
    in Angstrom of the given width; the grid is in units of R0, and a molecule
    is drawn by picking one of its points, so this array is both the truth
    plotted later and the distribution actually sampled from."""
    rel = model['rel']; R0 = st['R0']
    #: the cells of the grid, so a Gaussian is integrated rather than sampled
    mid = 0.5 * (rel[1:] + rel[:-1])
    lo = np.concatenate([[rel[0] - (mid[0] - rel[0])], mid])
    hi = np.concatenate([mid, [rel[-1] + (rel[-1] - mid[-1])]])
    from scipy.stats import norm
    p = np.zeros_like(rel)
    tot = sum(s['fraction'] for s in st['species'])
    for s in st['species']:
        mu, sd = s['R'] / R0, max(s['width'], 1e-6) / R0
        p += (s['fraction'] / tot) * (norm.cdf(hi, mu, sd) - norm.cdf(lo, mu, sd))
    if p.sum() <= 0:
        raise ValueError('the species lie outside the distance grid')
    return p / p.sum()


def constant_values(model: dict, st: dict):
    """The physics constants of the simulated setup, as the prototype's value
    dict: crosstalks, quantum yields, direct excitation, the g factor, the
    polarisation mixing, the anisotropies and the rotational times, and a
    nominal response for every detector."""
    L = model['L']; S = L.S; Ep = model['Ep']
    tt = L.tt
    w_rho = L.log_bump(L.M.RHO_GRID, st['rho_donor'], 0.25)
    w_rho_a = L.log_bump(L.M.RHO_A_GRID, st['rho_acceptor'], 0.25)
    vals = dict(x_d0=tt(0.0), w_a=tt(S.W_A_TRUE), w_rho=tt(w_rho), w_rho_a=tt(w_rho_a),
                r0_d=tt(st['r0_donor']), r0_a=tt(st['r0_acceptor']),
                g=tt(S.G_TRUE), l1=tt(S.L1), l2=tt(S.L2),
                C_GD=tt(S.C_GD), C_GA=tt(S.C_GA), C_RD=tt(S.C_RD), C_RA=tt(S.C_RA),
                G_GREEN=tt(S.G_GREEN), G_RED=tt(S.G_RED),
                QY_D=tt(S.QY_D), QY_A=tt(S.QY_A), EX_AG=tt(S.EX_AG))
    nom = Ep['instrument']
    for d in sorted({L.parse_channel(k)[1] for k in model['keys']}):
        d0 = (Ep.get('pulse_alias') or {}).get(d, d)
        vals[f'irf_shift_{d0}'] = tt(nom[2]); vals[f'irf_width_{d0}'] = tt(nom[1]); vals[f'irf_skew_{d0}'] = tt(nom[3])
    return vals


def fine_response(model: dict, det: str, vals) -> "torch.Tensor":
    """The response of one detector on the FINE time axis (32 ps bins over the
    50 ns period).  A detector seen from the red pulse ('r2v' is 'rv' 25 ns
    later) is the same response displaced by the pulse delay."""
    L = model['L']; Ep = model['Ep']
    off = (Ep.get('pulse_offset') or {}).get(det, 0.0)
    d0 = (Ep.get('pulse_alias') or {}).get(det, det)
    return L.analytic_irf(Ep, vals[f'irf_shift_{d0}'] + off, vals[f'irf_width_{d0}'], vals[f'irf_skew_{d0}'])


def physics_tables(model: dict, st: dict, p_true: np.ndarray, verbose=True):
    """For every point of the distance grid, and for a donor-only molecule,
    the expected micro-time distribution in each detector and the relative
    brightness of that detector.

    This is the whole forward model evaluated one molecule at a time: the
    amplitudes are linear in the distance distribution, so putting all the
    weight on one grid point gives that single molecule's decay, and the
    mixture of them over p(R/R0) is exactly the ensemble decay the fit sees.
    """
    import torch
    L = model['L']; Ep = model['Ep']; keys = model['keys']
    vals = constant_values(model, st)
    spec = Ep['cD']                                  # the donor's lifetime spectrum
    dets = {k: L.parse_channel(k)[1] for k in keys}
    B = {d: L.basis_from_irf(Ep, fine_response(model, d, vals), rebin=False) for d in set(dets.values())}
    used = np.nonzero(p_true > 1e-9)[0]
    nR = len(model['rel'])
    da_keys = [k for k in keys if k[0] == 'DA']
    a0_keys = [k for k in keys if k[0] == 'A0']
    d0_keys = [k for k in keys if k[0] == 'D0']

    def decays(kk, a):
        w, cdf = [], []
        for k in kk:
            v = (B[dets[k]] @ a[k]).clamp_min(0.0).numpy()
            s = v.sum(); w.append(s)
            cdf.append(np.cumsum(v / max(s, 1e-300)))
        return np.array(w), np.array(cdf)

    #: the donor-only molecule (no acceptor at all): the D0 scope
    v0 = dict(vals); v0['x_d0'] = L.tt(1.0)
    a_d0 = L.physics_amplitudes(Ep, v0, spec, L.tt(p_true), d0_keys)
    w_d0, cdf_d0 = decays(d0_keys, a_d0)
    #: the acceptor under its own laser: the same for every labelled molecule
    a_a0 = L.physics_amplitudes(Ep, vals, spec, L.tt(p_true), a0_keys)
    w_a0, cdf_a0 = decays(a0_keys, a_a0)
    #: one molecule at each populated distance, under the green pulse
    w_da = np.zeros((nR, len(da_keys))); cdf_da = np.zeros((nR, len(da_keys), model['n_fine']))
    onehot = np.zeros(nR)
    for j in used:
        onehot[:] = 0.0; onehot[j] = 1.0
        a = L.physics_amplitudes(Ep, vals, spec, L.tt(onehot), da_keys)
        w_da[j], cdf_da[j] = decays(da_keys, a)
        if verbose and (j == used[0] or j == used[-1]):
            print(f'  grid point {j:3d}: R/R0 {model["rel"][j]:.3f}, '
                  f'green fraction {w_da[j][:2].sum() / w_da[j].sum():.3f}')
    #: the scatter: photons at the response itself, in the green detectors
    scat = {d: np.cumsum(B[d][:, 0].numpy() / max(float(B[d][:, 0].sum()), 1e-300)) for d in set(dets.values())}
    return dict(vals=vals, keys=dict(D0=d0_keys, DA=da_keys, A0=a0_keys), dets=dets,
                w_d0=w_d0, cdf_d0=cdf_d0, w_a0=w_a0, cdf_a0=cdf_a0, w_da=w_da, cdf_da=cdf_da,
                scatter=scat, used=used)


# --------------------------------------------------------------------------
# 3. the photon stream
# --------------------------------------------------------------------------

def _draw_micro(cdf: np.ndarray, n: int, rng) -> np.ndarray:
    """n micro times from a cumulative distribution over the fine time axis."""
    return np.minimum(np.searchsorted(cdf, rng.random(n)), len(cdf) - 1).astype(np.uint16)


def _sample_bursts(model, st, tab, rng, n_bursts, p_true, sample):
    """One measurement: molecules crossing the spot one at a time.

    `sample` is 'DA' (the labelled sample, both pulses) or 'D0' (the
    donor-only reference, green pulse only).  Returns the photon arrays and,
    for the labelled sample, which distance every burst had.
    """
    macro, micro, chan, flag, per_burst = [], [], [], [], []
    period_s = model['period'] * 1e-9
    dur_pulses = int(st['burst_duration_ms'] * 1e-3 / period_s)
    ref = tab['w_d0'].sum()                       # an unquenched donor's detected photons
    t = 0
    x_d0 = st['donor_only_fraction'] if sample == 'DA' else 1.0
    for b in range(n_bursts):
        t += int(rng.exponential(st['burst_spacing_ms'] * 1e-3 / period_s))
        donor_only = rng.random() < x_d0
        j = -1 if donor_only else int(rng.choice(len(p_true), p=p_true))
        w = tab['w_d0'] if donor_only else tab['w_da'][j]
        cdf = tab['cdf_d0'] if donor_only else tab['cdf_da'][j]
        n_g = rng.poisson(st['photons_per_burst'] * w.sum() / ref)
        parts = [(w, cdf, n_g, 0)]
        if sample == 'DA' and not donor_only:
            n_r = rng.poisson(st['photons_per_burst'] * st['red_green_ratio'])
            parts.append((tab['w_a0'], tab['cdf_a0'], n_r, 4))
        n_tot = 0
        for w_, cdf_, n_, base in parts:
            if n_ <= 0:
                continue
            c = rng.choice(len(w_), size=int(n_), p=w_ / w_.sum())
            for i in range(len(w_)):
                m = int((c == i).sum())
                if m:
                    micro.append(_draw_micro(cdf_[i], m, rng))
                    chan.append(np.full(m, (base + i) % 4, dtype=np.int8))
                    macro.append(t + rng.integers(0, max(dur_pulses, 1), m, dtype=np.int64))
                    flag.append(np.zeros(m, dtype=bool)); n_tot += m
        #: scattered excitation light, at the response itself
        n_s = rng.poisson(st['scatter_fraction'] * max(n_tot, 1))
        if n_s > 0:
            d = int(rng.integers(0, 2))                   # a green detector
            det = tab['dets'][tab['keys']['DA'][d]]
            micro.append(_draw_micro(tab['scatter'][det], n_s, rng))
            chan.append(np.full(n_s, d, dtype=np.int8))
            macro.append(t + rng.integers(0, max(dur_pulses, 1), n_s, dtype=np.int64))
            flag.append(np.zeros(n_s, dtype=bool))
        per_burst.append((t, n_tot, j))
        t += dur_pulses
    #: uncorrelated background over the whole record, flat in the micro time
    t_end = max(t, 1)
    n_bg = int(st['background_cps'] * 4 * t_end * period_s)
    if n_bg > 0:
        macro.append(rng.integers(0, t_end, n_bg, dtype=np.int64))
        micro.append(rng.integers(0, model['n_fine'], n_bg).astype(np.uint16))
        chan.append(rng.integers(0, 4, n_bg).astype(np.int8))
        flag.append(np.ones(n_bg, dtype=bool))
    macro = np.concatenate(macro); micro = np.concatenate(micro); chan = np.concatenate(chan)
    flag = np.concatenate(flag)
    o = np.argsort(macro, kind='stable')
    return dict(macro=macro[o].astype(np.uint64), micro=micro[o], chan=chan[o],
                background=flag[o], bursts=per_burst, duration_s=t_end * period_s, n_background=n_bg)


def simulate_stream(model: dict, st: dict, p_true: np.ndarray, tab: dict, rng=None) -> dict:
    """Both measurements of the experiment: the labelled sample under both
    pulses, and the donor-only reference sample."""
    rng = rng or np.random.default_rng(st['seed'])
    out = {}
    out['DA'] = _sample_bursts(model, st, tab, rng, st['n_bursts'], p_true, 'DA')
    out['D0'] = _sample_bursts(model, st, tab, rng, st['n_bursts_donor'], p_true, 'D0')
    return out


def to_tttr(model: dict, stream: dict):
    """The photon arrays as a `tttrlib.TTTR` object -- what a measurement file
    would give.  The macro time unit is the laser period; the micro time is in
    fine TCSPC channels over that period."""
    import tttrlib
    t = tttrlib.TTTR()
    n = len(stream['macro'])
    t.append_events(stream['macro'], stream['micro'], stream['chan'],
                    np.zeros(n, dtype=np.int8), shift_macro_time=False)
    h = t.get_header()
    h.set_macro_time_resolution(model['period'] * 1e-9)     # seconds per macro tick
    h.set_micro_time_resolution(model['dt'] * 1e-9)
    h.set_number_of_micro_time_channels(model['n_fine'])
    t.set_header(h)
    return t


# --------------------------------------------------------------------------
# 4. the burst search -- tttrlib
# --------------------------------------------------------------------------

def burst_search(model: dict, tttr, st: dict):
    """Bursts as tttrlib finds them: a sliding time window that must hold at
    least `search_min_photons`, and the resulting stretch of photons kept only
    if it is short enough and bright enough.

    Returns (index array of the burst photons, list of (start, stop) index
    pairs, one per burst).
    """
    import tttrlib
    cal = model['period'] * 1e-9                     # seconds per macro tick
    mt = np.asarray(tttr.macro_times, dtype=np.uint64)
    ranges = tttrlib.ranges_by_time_window(
        mt, st['search_window_ms'] * 1e-3, -1.0, int(st['search_min_photons']), -1, cal, False)
    ranges = np.asarray(ranges).reshape(-1, 2)
    #: consecutive windows that both passed belong to one molecule: merge them,
    #: so that a burst is the whole crossing and not the sliding window
    merged = []
    n_ph = len(mt)
    for a, b in ranges:
        b = min(int(b), n_ph - 1)
        if merged and a <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], int(b))
        else:
            merged.append([int(a), int(b)])
    keep, idx = [], []
    for a, b in merged:
        if b - a + 1 < st['burst_min_photons']:
            continue
        if (mt[b] - mt[a]) * cal * 1e3 > st['burst_max_duration_ms']:
            continue
        keep.append((int(a), int(b))); idx.append(np.arange(a, b + 1))
    sel = np.concatenate(idx) if idx else np.zeros(0, dtype=int)
    return sel, keep


def burst_observables(tttr, bursts, model):
    """Per burst, the numbers a single-molecule experiment is read through.

    Column 0, the **proximity ratio** E*: of the photons the green pulse
    produced, the fraction that came out red.  It rises with transfer, so it
    rises as the dyes come closer.

    Column 1, the **stoichiometry** S: the green pulse's share of all the
    photons.  A molecule with a working acceptor sits near one half; one that
    stays dark under the red pulse goes to one.  This is what PIE buys.

    Column 2, the burst size.

    Column 3, the **mean arrival time of the donor photons** in nanoseconds:
    the average micro time in the green detectors after the green pulse.  It
    is the donor's fluorescence lifetime plus the instrument's offset, and it
    falls as transfer quenches the donor.  It is an INDEPENDENT measure of the
    same distance -- the intensity says one thing and the clock says another,
    and for a static population the two agree.

    The pulse a photon belongs to is read from its micro time, since the red
    pulse sits half a period after the green one.
    """
    mi = np.asarray(tttr.micro_times); ch = np.asarray(tttr.routing_channels)
    split = int(round(model['t_pie'] / model['dt']))
    green_pulse = mi < split
    rows = []
    for a, b in bursts:
        s = slice(a, b + 1)
        g, r, gp = ch[s] < 2, ch[s] >= 2, green_pulse[s]
        n_gg = int((g & gp).sum()); n_gr = int((r & gp).sum()); n_rr = int((r & ~gp).sum())
        tot = n_gg + n_gr + n_rr
        if tot < 1:
            continue
        dm = mi[s][g & gp]
        tau = float(dm.mean()) * model['dt'] if len(dm) else float('nan')
        rows.append((n_gr / max(n_gg + n_gr, 1), (n_gg + n_gr) / tot, tot, tau))
    return np.array(rows)


# --------------------------------------------------------------------------
# 5. the decays of the burst photons
# --------------------------------------------------------------------------

def burst_decays(model: dict, tttr, sel, sample: str):
    """The micro times of the selected photons, histogrammed on the fine axis
    and then rebinned onto the model's axis -- fine right after a pulse, where
    the decay changes fastest, coarse before the next one."""
    import torch
    mi = np.asarray(tttr.micro_times)[sel]
    ch = np.asarray(tttr.routing_channels)[sel]
    L = model['L']; Ep = model['Ep']
    dets = ('gv', 'gh', 'rv', 'rh')
    kind = {'gv': 'vv', 'gh': 'vh', 'rv': 'vv', 'rh': 'vh'}
    out = {}
    for i, d in enumerate(dets):
        fine = np.bincount(mi[ch == i], minlength=model['n_fine']).astype(float)
        out[(sample, f'{d}_{kind[d]}')] = Ep['R'] @ torch.tensor(fine, dtype=torch.float64)
    return out


def assemble_data(model: dict, y_d0: dict, y_da: dict):
    """The eight measured histograms in the twelve-channel layout the model
    uses: the donor-only sample's four, the labelled sample's four -- each of
    which holds BOTH pulses -- and the partner channels that name the red
    pulse's view of the same histogram."""
    L = model['L']
    y = dict(y_d0); y.update(y_da)
    for d in L.PIE_DETS:
        ka = ('DA', next(k[1] for k in y_da if k[1].startswith(d + '_')))
        y[('A0', f'{L.PIE_PARTNER[d]}_{ka[1].split("_")[1]}')] = y[ka]
    return y


# --------------------------------------------------------------------------
# 6. the posterior
# --------------------------------------------------------------------------

def fit(model: dict, y: dict, st: dict, verbose=True):
    """The Laplace posterior of the whole model on the eight histograms.

    Route B of the prototype: a labelled histogram is fitted unmasked with
    both pulses in its mean, which is what the data are; the two pulses only
    split the histogram for the starting values.  The roughness weight of the
    P-spline prior on log p(R/R0) is not maximised but integrated out: a
    Laplace approximation at each node of a grid, the nodes mixed by their
    evidences (INLA; Rue, Martino & Chopin 2009).
    """
    import torch
    L = model['L']; Ep = dict(model['Ep']); keys = model['keys']; spl = model['spl']
    kind = {'gv': 'vv', 'gh': 'vh', 'rv': 'vv', 'rh': 'vh'}
    pairs = {('DA', f'{d}_{kind[d]}'): ('A0', f'{L.PIE_PARTNER[d]}_{kind[d]}') for d in L.PIE_DETS}
    Ep['pie_pairs'] = pairs
    win = {k: (torch.ones(model['n_bin']) if k[0] == 'D0' else
               (Ep['win_green'] if k[0] == 'DA' else Ep['win_red'])) for k in keys}
    data_keys = [k for k in keys if k[0] != 'A0']
    masks = {k: torch.ones(model['n_bin']) for k in data_keys}
    n_coef = spl.shape[1]
    med = {'D0': max(float(y[('D0', 'gv_vv')].sum()), 1.0),
           'DA': max(float(y[('DA', 'gv_vv')].sum()), 1.0),
           'A0': max(float((y[('DA', 'rv_vv')] * Ep['win_red']).sum()), 1.0)}
    V = L.default_variables(Ep, keys, n_coef=n_coef, scale_medians=med, link='softmax', irf_shape=False)
    ps = L.PSplineFactor(n_coef, link='softmax', spl=spl)
    g = L.FactorGraph(Ep, keys, V, L.PoissonCountsFactor({k: y[k] for k in data_keys}, mask=masks),
                      L.InstrumentModel(Ep, 'analytic'), spl, ps)
    g.rel = model['rel']
    g = g.with_fixed(**{f'bkg_{kp[0]}_{kp[1]}': L.tt([math.log(1e-9)]) for kp in pairs.values()})
    g.start_y = {k: y[k] * win[k] for k in keys}
    gen = torch.Generator().manual_seed(int(st['seed']))
    post = L.fit_sample(g, y, gen, model['rel'], verbose=verbose, start='mem', optimiser='fisher',
                        lam_nodes=tuple(float(x) for x in st['lam_nodes']), hessian='fisher')
    return g, post


def rule0_table(model, g, post, y):
    """Rule 0: the Poisson deviance per degree of freedom of every histogram
    and a runs test on its weighted residuals.  A number about a fitted
    quantity is not reportable without these."""
    rows = []
    for k, r in post['rows'].items():
        rows.append(dict(channel=f'{k[0]} {k[1]}', dpd=r['dpd'], runs_p=r['runs_p']))
    rows.append(dict(channel='all', dpd=post['dev'] / post['dof'], runs_p=float('nan')))
    return rows


def burst_truth(model: dict, stream: dict, tttr, bursts):
    """Which molecule every found burst was, by matching its photons back to
    the simulated crossings.

    This is what makes the comparison honest.  A burst search keeps the bright
    crossings, and brightness depends on the distance: a molecule at short
    distance sends most of its photons to the acceptor and, with the
    acceptor's own quantum yield and the red detector's efficiency, is not as
    bright as an unquenched donor.  So the distribution of the molecules
    BEHIND THE SELECTED BURSTS is not the distribution in the cuvette, and it
    is the former that the decays measure.  Returns (p_selected, x_d0_selected,
    counts per grid point).
    """
    mt = np.asarray(tttr.macro_times)
    starts = np.array([b[0] for b in stream['bursts']], dtype=np.int64)
    dur = np.diff(np.append(starts, starts[-1] + 1))
    j_of = np.array([b[2] for b in stream['bursts']], dtype=np.int64)
    nR = len(model['rel'])
    cnt = np.zeros(nR); n_d0 = 0; n_hit = 0
    for a, b in bursts:
        t = int(np.median(mt[a:b + 1]))
        i = int(np.searchsorted(starts, t, side='right')) - 1
        if i < 0 or t - starts[i] > dur[i]:
            continue                      # a stretch of background alone
        n_hit += 1
        if j_of[i] < 0:
            n_d0 += 1
        else:
            cnt[j_of[i]] += 1
    p = cnt / max(cnt.sum(), 1.0)
    return p, n_d0 / max(n_hit, 1), cnt


def poisson_reference(post, n_draw=200, seed=0):
    """What the deviance would be if the fitted means were the truth.

    The Poisson deviance of a CORRECT model is not the number of degrees of
    freedom: for counts of this size it is close to the number of BINS, and
    the reported dof subtracts the full parameter dimension even though a
    penalised spline uses far fewer of them.  So the reference is measured
    rather than assumed -- Poisson draws at the fitted means, their deviance
    against those same means, per degree of freedom.  A fit is good when
    D/dof sits inside this reference's own scatter.
    """
    import torch
    rng = np.random.default_rng(seed)
    lam = {k: np.asarray(v) for k, v in post['lam'].items()}
    sel = {k: post['rows'][k]['sel'] for k in lam}
    d = []
    for _ in range(n_draw):
        tot = 0.0
        for k, l in lam.items():
            m = l[sel[k]]
            yk = rng.poisson(m)
            with np.errstate(divide='ignore', invalid='ignore'):
                t = np.where(yk > 0, yk * np.log(np.maximum(yk, 1e-300) / m), 0.0)
            tot += 2.0 * float((t - (yk - m)).sum())
        d.append(tot)
    d = np.array(d) / post['dof']
    return float(d.mean()), float(d.std())


# --------------------------------------------------------------------------
# 7. figures
# --------------------------------------------------------------------------

def _inline():
    """Restore the notebook's inline figure backend.

    Importing parts of the prototype pulls in a module that calls
    `matplotlib.use('Agg')`, and after that a notebook cell draws nothing at
    all -- silently, which is the worst way for a figure to go missing.  Every
    plotting function here puts the backend back first.
    """
    try:
        from IPython import get_ipython
        ip = get_ipython()
        if ip is None:
            return
        import matplotlib
        if 'inline' not in matplotlib.get_backend():
            ip.run_line_magic('matplotlib', 'inline')
    except Exception:
        pass


def plot_decays(model, y, post, title=''):
    """The eight histograms with the fitted model, and the weighted residuals
    (data minus model over the square root of the model, the Pearson residual
    of Poisson counts) under each one."""
    _inline()
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    tc = model['Ep']['tc']; wid = model['Ep']['wid']
    keys = [k for k in sorted(post['lam'], key=str)]
    fig = plt.figure(figsize=(13, 7.5))
    gs = GridSpec(4, 4, height_ratios=[3, 1.15, 3, 1.15], hspace=0.05, wspace=0.22)
    for i, k in enumerate(keys):
        r, c = divmod(i, 4)
        ax = fig.add_subplot(gs[2 * r, c]); ar = fig.add_subplot(gs[2 * r + 1, c], sharex=ax)
        yk = np.asarray(y[k]); lk = np.asarray(post['lam'][k])
        ax.step(tc, yk / wid, where='mid', color='0.6', lw=0.8, label='counts')
        ax.step(tc, lk / wid, where='mid', color='C3', lw=1.0, label='model')
        ax.set_yscale('log'); ax.set_ylim(max(0.3, np.min(lk / wid) * 0.5), np.max(yk / wid) * 2)
        ax.set_title(f'{k[0]} {k[1]}   D/dof {post["rows"][k]["dpd"]:.2f}, runs p {post["rows"][k]["runs_p"]:.2f}', fontsize=8)
        ax.tick_params(labelbottom=False, labelsize=7)
        w = (yk - lk) / np.sqrt(np.maximum(lk, 1e-12))
        ar.axhline(0, color='0.7', lw=0.6)
        ar.step(tc, w, where='mid', color='C0', lw=0.6)
        ar.set_ylim(-4.5, 4.5); ar.tick_params(labelsize=7)
        if r == 1:
            ar.set_xlabel('micro time / ns', fontsize=8)
        if c == 0:
            ax.set_ylabel('counts per ns', fontsize=8); ar.set_ylabel('w. res.', fontsize=8)
        if i == 0:
            ax.legend(fontsize=6, frameon=False)
    fig.suptitle(title or 'burst decays and weighted residuals', fontsize=10)
    return fig


def plot_es(rows, title=''):
    """The burst-wise proximity ratio against the stoichiometry: the standard
    PIE plot.  Molecules with no acceptor have S near one; the FRET
    populations sit at S near one half."""
    _inline()
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.4), gridspec_kw=dict(width_ratios=[2, 1]))
    ax[0].hexbin(rows[:, 0], rows[:, 1], gridsize=45, extent=(-0.1, 1.1, 0, 1.05), cmap='viridis', mincnt=1)
    ax[0].set_xlabel('proximity ratio E*'); ax[0].set_ylabel('stoichiometry S'); ax[0].set_title(title or 'bursts', fontsize=9)
    ax[1].hist(rows[rows[:, 1] < 0.8, 0], bins=40, range=(-0.1, 1.1), color='C0')
    ax[1].set_xlabel('E* of the acceptor-carrying bursts'); ax[1].set_ylabel('bursts')
    fig.tight_layout()
    return fig


def plot_pR(model, post, p_true=None, p_selected=None, dE=0.01, title=''):
    """The posterior of p(R/R0) with its band, against the truth.

    The band is the delta method at each grid point over the evidence mixture
    of the penalty nodes, clipped at zero.  The grey regions are where the
    transfer efficiency is within dE of zero or one: no measurement of the
    decay can place a distance there, and a fit that draws a line through
    them is drawing its prior.
    """
    _inline()
    import matplotlib.pyplot as plt
    L = model['L']; rel = model['rel']
    m, lo1, hi1, lo2, hi2 = L.delta_bands(None, post, rel, space='linear')
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    ax.fill_between(rel, lo2, hi2, color='C0', alpha=0.18, lw=0, label='posterior 2 sd')
    ax.fill_between(rel, lo1, hi1, color='C0', alpha=0.35, lw=0, label='posterior 1 sd')
    ax.plot(rel, m, color='C0', lw=1.6, label='posterior mean')
    if p_true is not None:
        ax.plot(rel, p_true, color='k', lw=1.2, ls='--', label='truth (in the cuvette)')
    if p_selected is not None:
        ax.plot(rel, p_selected, color='C3', lw=1.2, label='truth (the selected bursts)')
    ax.set_xlabel('R / R0'); ax.set_ylabel('p(R/R0) per grid point')
    ax.set_xlim(rel[0], rel[-1]); ax.set_ylim(0, None)
    L.shade_window(ax, dE)
    ax.legend(fontsize=7, frameon=False)
    ax.set_title(title or 'the distance distribution', fontsize=10)
    fig.tight_layout()
    return fig


def summary_table(model, post, truth: dict):
    """Posterior mean, sd and 95 % interval of the quantities the experiment
    is for, beside their true values."""
    L = model['L']
    mom = L.mixture_summary_moments(post['graph'], post, model['rel'])
    rows = []
    for n in ('x_d0', 'mean', 'sd', 'tail', 'r0_d', 'r0_a', 'g'):
        w = sum(wi for wi, _ in mom.values())
        mu = sum(wi * dm[n][0] for wi, dm in mom.values()) / w
        sd = math.sqrt(sum(wi * (dm[n][1] ** 2 + (dm[n][0] - mu) ** 2) for wi, dm in mom.values()) / w)
        q = [L.mixture_quantile(mom, n, q) for q in (0.025, 0.975)]
        t = truth.get(n, float('nan'))
        rows.append(dict(quantity=n, truth=t, mean=mu, sd=sd, lo=q[0], hi=q[1],
                         inside=bool(q[0] <= t <= q[1]) if t == t else None))
    return rows


# --------------------------------------------------------------------------
# 8. does the interval mean what it says?
# --------------------------------------------------------------------------

def run_once(model: dict, st: dict, tab: dict, p_true: np.ndarray, seed: int, verbose=False):
    """The whole pipeline once, from photons to posterior, at one seed."""
    rng = np.random.default_rng(seed)
    s2 = dict(st); s2['seed'] = seed
    stream = simulate_stream(model, s2, p_true, tab, rng)
    tttr = {s: to_tttr(model, d) for s, d in stream.items()}
    sel, bursts = {}, {}
    for s in ('DA', 'D0'):
        sel[s], bursts[s] = burst_search(model, tttr[s], s2)
    p_sel, x_sel, _ = burst_truth(model, stream['DA'], tttr['DA'], bursts['DA'])
    y = assemble_data(model, burst_decays(model, tttr['D0'], sel['D0'], 'D0'),
                      burst_decays(model, tttr['DA'], sel['DA'], 'DA'))
    graph, post = fit(model, y, s2, verbose=verbose)
    rel = model['rel']; mu = float((p_sel * rel).sum())
    truth = dict(x_d0=x_sel, mean=mu, sd=float(np.sqrt((p_sel * (rel - mu) ** 2).sum())),
                 tail=float(p_sel[rel > 1.3].sum()), r0_d=st['r0_donor'], r0_a=st['r0_acceptor'],
                 g=model['L'].S.G_TRUE)
    return dict(y=y, graph=graph, post=post, truth=truth, p_selected=p_sel, stream=stream,
                tttr=tttr, sel=sel, bursts=bursts)


def repeat_experiment(model, st, tab, p_true, seeds, verbose=False):
    """The same experiment at several seeds, so that the intervals can be
    checked against something that could contradict them.

    A single posterior covering its truth is not evidence that the interval is
    honest, and one that misses is not evidence that it is not: only the
    spread over repeated realisations says anything.  The pull -- (posterior
    mean minus truth) divided by the posterior sd -- should scatter about zero
    with a spread of one.
    """
    names = ('x_d0', 'mean', 'sd', 'tail', 'r0_d', 'r0_a', 'g')
    out = {n: [] for n in names}
    stats = []
    for s in seeds:
        r = run_once(model, st, tab, p_true, s, verbose=verbose)
        ref, ref_sd = poisson_reference(r['post'], 100, s)
        stats.append(dict(seed=s, dpd=r['post']['dev'] / r['post']['dof'], ref=ref, ref_sd=ref_sd))
        for row in summary_table(model, r['post'], r['truth']):
            out[row['quantity']].append((row['truth'], row['mean'], row['sd']))
    return out, stats


def plot_pulls(out, title=''):
    """The pull per quantity over the repeated experiments."""
    _inline()
    import matplotlib.pyplot as plt
    names = list(out)
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    for i, n in enumerate(names):
        z = [(m - t) / s for t, m, s in out[n]]
        ax.plot(np.full(len(z), i) + np.linspace(-0.12, 0.12, len(z)), z, 'o', ms=5, color='C0')
        ax.plot([i - 0.25, i + 0.25], [np.mean(z)] * 2, '-', color='C3', lw=2)
    ax.axhline(0, color='0.6', lw=0.8)
    ax.axhspan(-1, 1, color='0.85', zorder=0, lw=0)
    ax.set_xticks(range(len(names))); ax.set_xticklabels(names)
    ax.set_ylabel('(posterior mean - truth) / posterior sd')
    ax.set_ylim(-4, 4)
    ax.set_title(title or 'pulls over repeated experiments (grey: one posterior sd; red: the mean pull)', fontsize=9)
    fig.tight_layout()
    return fig


def node_table(post):
    """The grid over the roughness weight: what each node fitted, how much
    evidence it carries, and therefore its weight in the posterior mixture."""
    rows = []
    for l, w in zip(post['lam_nodes'], post['weights']):
        nd = post['nodes'][l]
        rows.append(dict(log10_lam=l, evidence=nd['evidence'], weight=float(w),
                         dpd=nd['dev'] / nd['dof'], iterations=nd['n_iter']))
    rows.sort(key=lambda r: r['log10_lam'])
    return rows


def background_fraction(stream, sel):
    """What fraction of the SELECTED photons of each detector is uncorrelated
    background.  The simulation knows which photons it drew as background; a
    real measurement would estimate it, and the fit carries it as a free
    parameter per histogram either way."""
    bg = stream['background'][sel]; ch = stream['chan'][sel]
    return {int(c): (float(bg[ch == c].mean()) if int((ch == c).sum()) else 0.0) for c in range(4)}


# --------------------------------------------------------------------------
# 9. two steps: classify the bursts, then analyse each group
# --------------------------------------------------------------------------

def _hdbscan_select(parents, children, lambdas, sizes, n_points, method='eom'):
    """Which nodes of a condensed tree are the clusters.

    tttrlib builds the mutual-reachability tree and condenses it, and it can
    label the points once told which nodes to keep -- but nothing in it
    computes that selection, which is the whole of HDBSCAN's cluster
    extraction.  This is the excess-of-mass rule of Campello, Moulavi &
    Sander, "Density-Based Clustering Based on Hierarchical Density
    Estimates", PAKDD 2013, section 4, as implemented in `hdbscan` (McInnes,
    Healy & Astels, J. Open Source Software 2, 205, 2017) and in
    `sklearn.cluster.HDBSCAN`:

        stability(C) = sum over the points p leaving C of (lambda_p - lambda_birth(C)),
        lambda = 1 / distance,

    and a cluster is kept only if its own stability exceeds the total
    stability of the selected clusters below it; otherwise those descendants
    are kept instead. `method='leaf'` keeps every leaf of the tree, which
    separates populations that excess of mass merges.

    `check_hdbscan_selection` compares what this produces against
    scikit-learn's labels on data where the two can disagree. That check can
    fail, and it did: the first version let the root compete with its
    children, and since the root holds every point from lambda 0 its
    stability beats any set of descendants -- one cluster containing the whole
    sample. The reference implementation walks every node except the root.
    """
    parents = np.asarray(parents, np.int64); children = np.asarray(children, np.int64)
    lam = np.asarray(lambdas, float); size = np.asarray(sizes, np.int64)
    nodes = np.unique(parents)
    root = int(nodes.min())
    birth = {int(ch): float(lm) for ch, lm in zip(children, lam) if ch >= n_points}
    birth[root] = 0.0
    stab = {int(c): 0.0 for c in nodes}
    for pa, lm, sz in zip(parents, lam, size):
        stab[int(pa)] += float(sz) * (float(lm) - birth[int(pa)])
    kids = {int(c): [] for c in nodes}
    for pa, ch in zip(parents, children):
        if ch >= n_points:
            kids[int(pa)].append(int(ch))
    if method == 'leaf':
        sel = {int(c): (len(kids[int(c)]) == 0) for c in nodes}
    else:
        sel = {int(c): True for c in nodes}
        for c in [int(x) for x in sorted(nodes, reverse=True) if int(x) != root]:
            if not kids[c]:
                continue
            below, stack = 0.0, list(kids[c])       # the selected frontier under c
            while stack:
                k = stack.pop()
                if sel[k]:
                    below += stab[k]
                else:
                    stack.extend(kids[k])
            if stab[c] < below:
                sel[c] = False
            else:
                stack = list(kids[c])
                while stack:
                    k = stack.pop(); sel[k] = False; stack.extend(kids[k])
    sel[root] = False
    is_sel = np.zeros(int(max(int(nodes.max()), int(children.max()))) + 1, dtype=np.uint8)
    for c, v in sel.items():
        is_sel[c] = 1 if v else 0
    return is_sel


def hdbscan_labels(X, min_cluster_size, min_samples=None, method='eom', alpha=1.0):
    """HDBSCAN through tttrlib: the mutual-reachability minimum spanning tree,
    the condensed tree, the cluster selection above, and the labels.

    tttrlib assigns every point to a cluster; scikit-learn calls points
    outside the selected clusters noise. Labels are renumbered 0..k-1.
    """
    import tttrlib
    X = np.ascontiguousarray(np.asarray(X, float))
    ms = int(min_samples or min_cluster_size)
    #: tttrlib gained a complete `hdbscan` -- selection, labels and membership
    #: strengths -- after this example was written (reported from here, 2026-09-08).
    #: Prefer it where it exists; the code below is the fallback, and is what
    #: produced the outputs stored in the notebooks.
    if hasattr(tttrlib, 'hdbscan'):
        r = tttrlib.hdbscan(X, int(min_cluster_size), ms, float(alpha), method, False, 0.0, 0)
        lab_new = np.asarray(r.labels)
        out_new = np.full(len(X), -1, dtype=int)
        for i, c in enumerate(sorted({int(v) for v in lab_new if v >= 0})):
            out_new[lab_new == c] = i
        return out_new
    mst = np.asarray(tttrlib.mutual_reachability_mst(X, ms, float(alpha)))
    #: linkage is order-dependent and tttrlib refuses unsorted edges
    mst = mst[np.argsort(mst[:, 2], kind='stable')]
    src = np.ascontiguousarray(mst[:, 0].astype(np.int64))
    dst = np.ascontiguousarray(mst[:, 1].astype(np.int64))
    w = np.ascontiguousarray(mst[:, 2].astype(float))
    parents, children, lambdas, sizes = tttrlib.hdbscan_condensed_tree(src, dst, w, int(min_cluster_size))
    is_sel = _hdbscan_select(parents, children, lambdas, sizes, len(X), method)
    lab = np.asarray(tttrlib.hdbscan_label_points(np.ascontiguousarray(np.asarray(parents, np.int64)),
                                                  np.ascontiguousarray(np.asarray(children, np.int64)),
                                                  is_sel, len(X)))
    out = np.full(len(X), -1, dtype=int)
    for i, c in enumerate(sorted({int(v) for v in lab if v >= 0})):
        out[lab == c] = i
    return out


def check_hdbscan_selection(seeds=(0, 1, 2, 3, 4), verbose=True):
    """The check that could fail: tttrlib's HDBSCAN with the selection above,
    against scikit-learn's, on data where the two could disagree.

    Well-separated blobs agree trivially and prove nothing, so the cases here
    include blobs that OVERLAP -- where excess of mass has a real choice
    between a parent and its children, and a wrong stability shows up as a
    different number of clusters -- and several minimum cluster sizes, which
    move the condensed tree. Agreement is the adjusted Rand index, 1 only if
    the partitions match up to relabelling. Because tttrlib has no noise
    label, the comparison is also reported over the points scikit-learn does
    cluster.
    """
    from sklearn.cluster import HDBSCAN
    from sklearn.metrics import adjusted_rand_score
    rows = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        for sep, sd in ((3.0, 0.25), (1.5, 0.45), (1.0, 0.5)):
            X = np.vstack([rng.normal([0, 0], sd, (200, 2)), rng.normal([sep, sep], sd, (200, 2)),
                           rng.normal([0, sep], sd, (150, 2))])
            for mcs in (15, 25, 60):
                for method in ('eom', 'leaf'):
                    a = hdbscan_labels(X, min_cluster_size=mcs, method=method)
                    b = HDBSCAN(min_cluster_size=mcs, cluster_selection_method=method, copy=True).fit_predict(X)
                    keep = b >= 0
                    rows.append(dict(case=f'seed {seed}, blobs {sep}/{sd}, min size {mcs}', method=method,
                                     ari=float(adjusted_rand_score(a, b)),
                                     ari_core=float(adjusted_rand_score(a[keep], b[keep])) if keep.sum() > 1 else float('nan'),
                                     n_tttrlib=len(set(a[a >= 0])), n_sklearn=len(set(b[b >= 0])),
                                     noise_sklearn=int((b < 0).sum())))
    if verbose:
        print('  tttrlib labels every point; scikit-learn calls points outside the selected clusters noise,')
        print('  so "core" compares only the points scikit-learn does cluster.')
        for m in ('eom', 'leaf'):
            r = [x for x in rows if x['method'] == m]
            core = np.array([x['ari_core'] for x in r]); ari = np.array([x['ari'] for x in r])
            same = sum(1 for x in r if x['n_tttrlib'] == x['n_sklearn'])
            print(f'  {m}: {len(r)} cases; core Rand min {np.nanmin(core):.4f}, mean {np.nanmean(core):.4f}; '
                  f'all points {ari.mean():.4f}; cluster count agrees in {same}/{len(r)}')
            for x in sorted(r, key=lambda x: x['ari_core'])[:2]:
                if x['ari_core'] < 0.999:
                    print(f"     weakest: {x['case']} -- core Rand {x['ari_core']:.3f}, "
                          f"{x['n_tttrlib']} vs {x['n_sklearn']} clusters, sklearn noise {x['noise_sklearn']}")
    return rows


def classify_bursts(feat, min_cluster_size=None, columns=(0, 1, 3), method='leaf'):
    """Group the bursts by what each one looks like on its own.

    A pooled decay of everything is a mixture, and the distance distribution
    then has to carry every population at once. But each burst already says
    roughly which population it belongs to: its proximity ratio, its
    stoichiometry, and the mean arrival time of its donor photons. Cluster on
    those three, pool each cluster, and each pooled decay is a much narrower
    problem.

    HDBSCAN through tttrlib. It is preferred to k-means because the number of
    populations is what one wants to learn rather than declare.

    `method`: 'leaf' keeps every leaf of the condensed tree, 'eom' the
    excess-of-mass clusters. Excess of mass prefers the parent whenever a
    split is not a deep valley in density, and two FRET populations whose
    per-burst observables overlap by a few standard deviations are exactly
    that case -- it returns them merged. 'leaf' is the default for that
    reason, and the choice is visible in the figure: if the leaves split one
    population in two, the plot shows it.

    The features are standardised first, so that a nanosecond and a unit of
    proximity ratio count the same; without that the clustering is a statement
    about the units.
    """
    X = np.asarray(feat, float)[:, list(columns)]
    ok = np.isfinite(X).all(1)
    Z = (X - np.nanmean(X[ok], 0)) / np.nanstd(X[ok], 0)
    mcs = int(min_cluster_size or max(25, int(ok.sum()) // 25))
    lab = np.full(len(X), -1, dtype=int)
    lab[ok] = hdbscan_labels(Z[ok], min_cluster_size=mcs, method=method)
    return lab


def group_indices(bursts, labels, group):
    """The photon indices of every burst in one group."""
    idx = [np.arange(a, b + 1) for (a, b), l in zip(bursts, labels) if l == group]
    return np.concatenate(idx) if idx else np.zeros(0, dtype=int)


def group_summary(feat, labels):
    """What each group is, in the observables that made it."""
    rows = []
    for g in sorted({int(l) for l in labels}):
        m = labels == g
        rows.append(dict(group=g, n=int(m.sum()), photons=float(feat[m, 2].sum()),
                         E=float(np.nanmean(feat[m, 0])), E_sd=float(np.nanstd(feat[m, 0])),
                         S=float(np.nanmean(feat[m, 1])), tau=float(np.nanmean(feat[m, 3]))))
    return rows


def label_groups(rows, s_donor_only=0.85):
    """Name the groups: a high stoichiometry means no working acceptor; the
    rest are FRET populations, numbered by proximity ratio."""
    names = {}
    fret = [r for r in rows if r['group'] >= 0 and r['S'] < s_donor_only]
    for r in rows:
        if r['group'] < 0:
            names[r['group']] = 'unclustered'
        elif r['S'] >= s_donor_only:
            names[r['group']] = 'donor only'
    for i, r in enumerate(sorted(fret, key=lambda r: r['E'])):
        names[r['group']] = f'FRET {i + 1}'
    return names


def plot_groups(feat, labels, names=None, title=''):
    """The bursts in the plane that made the groups."""
    _inline()
    import matplotlib.pyplot as plt
    names = names or {}
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.8))
    for g in sorted({int(l) for l in labels}):
        m = labels == g
        c = '0.75' if g < 0 else f'C{g % 10}'
        ax[0].plot(feat[m, 0], feat[m, 1], '.', ms=3, color=c, alpha=0.6,
                   label=f"{names.get(g, g)} ({int(m.sum())})")
        ax[1].plot(feat[m, 0], feat[m, 3], '.', ms=3, color=c, alpha=0.6)
    ax[0].set_xlabel('proximity ratio E*'); ax[0].set_ylabel('stoichiometry S')
    ax[0].set_xlim(-0.1, 1.1); ax[0].set_ylim(0, 1.05)
    ax[0].legend(fontsize=6.5, frameon=False, markerscale=2.5, loc='lower left')
    ax[1].set_xlabel('proximity ratio E*'); ax[1].set_ylabel('mean donor arrival time / ns')
    ax[1].set_xlim(-0.1, 1.1)
    fig.suptitle(title or 'step 1: the bursts classified by what each one looks like', fontsize=10)
    fig.tight_layout()
    return fig
