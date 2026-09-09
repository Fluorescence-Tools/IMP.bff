"""What was measured, as data rather than as code.

A time-resolved FRET measurement is described by four lists and nothing else:

* the **samples** -- what was in the cuvette, or which molecules the burst
  belongs to.  A sample carries a donor, an acceptor, or both.
* the **excitations** -- the laser pulses.  Each excites some chromophores and
  arrives at some delay within the period.
* the **detectors** -- each sees one colour and has its own instrument
  response and timing shift.
* the **channels** -- the histograms somebody actually recorded.  One channel
  is a sample, an excitation, a detector and a polarisation.

Everything the model needs follows from those four: which physics each
histogram carries, which unknowns exist, which of them the data can move, and
what the factor graph looks like.  Before this module the channel list was a
literal in the source and the physics each channel carries was a six-entry
dictionary written out by hand, so a colleague with four histograms instead of
twelve, or one laser instead of two, had to edit the model.

    exp = mfd_pie()                     # or build one field by field
    exp.to_json('my_experiment.json')
    exp = Experiment.from_json('my_experiment.json')
    exp.channel_keys()                  # the histograms, in the model's grammar
    exp.scope_table()                   # what physics each one carries

**THE SCOPES ARE DERIVED, NOT LISTED.**  A histogram carries a directly excited
donor if its sample has a donor and its pulse reaches the donor; it carries
FRET if the sample has both chromophores and the pulse reaches the donor; and
so on.  The six cases that interleaved excitation produces are a CONSEQUENCE of
that rule, and `check_against_prototype()` asserts that the rule reproduces the
hand-written table it replaces, key by key, without having been shown it.

Author: written for the bff examples, 2026-09-09.
See `okf/prd-experiment-specification.md` in the ucfret repository for why.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Sequence

__all__ = ['Detector', 'Excitation', 'Sample', 'Channel', 'Scope', 'Experiment',
           'build', 'graph', 'simulate', 'truth_values', 'fit', 'rule0',
           'identifiability', 'unidentified', 'factor_scopes', 'structural_graph',
           'mfd_donor_excitation', 'mfd_pie', 'separate_measurements',
           'magic_angle_minimal', 'GEOMETRIES', 'check_against_prototype']

#: the two chromophores this model knows about
CHROMOPHORES = ('donor', 'acceptor')

#: how a polarisation is written in the model's channel keys.  The spec accepts
#: the readable names and stores those; the codes are what the physics uses.
POLARISATION = {'parallel': 'vv', 'perpendicular': 'vh', 'magic': 'ma',
                'vv': 'vv', 'vh': 'vh', 'ma': 'ma'}
_READABLE = {'vv': 'parallel', 'vh': 'perpendicular', 'ma': 'magic'}


# --------------------------------------------------------------------------
# the four things an experiment is made of
# --------------------------------------------------------------------------

@dataclass
class Detector:
    """One detection channel with its own instrument response.

    `name` keys the response and the timing shift, so two polarisations
    recorded through the SAME response share a detector (the ensemble
    geometry's `g` and `r`) while two that are calibrated apart do not (the
    MFD geometry's `gv`, `gh`, `rv`, `rh`).

    `colour` is which emitter the detector mostly sees, `green` for donor
    emission and `red` for acceptor emission.  It decides which row of the
    crosstalk matrix the channel is on -- not that the channel sees only that
    emitter, which is what crosstalk means.
    """
    name: str
    colour: str

    def __post_init__(self):
        if self.colour not in ('green', 'red'):
            raise ValueError(f"detector {self.name!r}: colour must be 'green' or 'red', "
                             f"not {self.colour!r}")


@dataclass
class Excitation:
    """One laser pulse, and HOW STRONGLY it reaches each chromophore.

    `excites` maps a chromophore to the relative probability that this pulse
    puts it in the excited state.  The strength is either a number -- `1.0` for
    the chromophore the wavelength was chosen for -- or the NAME of a model
    parameter, for a path whose strength is not known and has to be fitted::

        Excitation('green', {'donor': 1.0, 'acceptor': 'EX_AG'})
        Excitation('red',   {'acceptor': 1.0, 'donor': 'EX_DR'})

    Saying it this way is the difference between direct excitation being a term
    in the model and being an assumption.  It is also what separates the two
    pulses: a green and a red pulse reach the same two chromophores, and only
    the strengths tell them apart -- which is why a scope cannot be four
    booleans, as the first version of this file had it and the gate caught.

    `delay_ns` is where the pulse sits in the laser period; the first is at zero.
    """
    name: str
    excites: dict
    delay_ns: float = 0.0

    def __post_init__(self):
        if not isinstance(self.excites, dict):
            #: a bare list is read as "reaches these, at full strength"
            self.excites = {c: 1.0 for c in self.excites}
        bad = [c for c in self.excites if c not in CHROMOPHORES]
        if bad:
            raise ValueError(f"excitation {self.name!r}: unknown chromophore(s) {bad}")
        if not self.excites:
            raise ValueError(f"excitation {self.name!r} excites nothing; leave it out instead")
        for c, v in self.excites.items():
            if not isinstance(v, (int, float, str)):
                raise ValueError(f"excitation {self.name!r}: the strength for {c!r} must be a "
                                 f"number or the name of a parameter, not {v!r}")

    def strength(self, chromophore):
        """the relative excitation probability, `0` if the pulse misses it"""
        return self.excites.get(chromophore, 0.0)

    def reaches(self, chromophore) -> bool:
        return chromophore in self.excites


@dataclass
class Sample:
    """What the molecules carry.

    A donor-only reference carries `('donor',)`, an acceptor-only reference
    `('acceptor',)`, a labelled sample both.  A labelled sample is never
    entirely labelled -- some molecules have lost their acceptor -- and that
    fraction is a free parameter of the model, not of the description.
    """
    name: str
    carries: tuple

    def __post_init__(self):
        self.carries = tuple(self.carries)
        bad = [c for c in self.carries if c not in CHROMOPHORES]
        if bad:
            raise ValueError(f"sample {self.name!r}: unknown chromophore(s) {bad}")
        if not self.carries:
            raise ValueError(f"sample {self.name!r} carries no chromophore; it cannot emit")


@dataclass
class Channel:
    """One recorded histogram."""
    sample: str
    excitation: str
    detector: str
    polarisation: str = 'magic'

    def __post_init__(self):
        p = POLARISATION.get(self.polarisation)
        if p is None:
            raise ValueError(f"channel {self.sample}/{self.excitation}/{self.detector}: "
                             f"polarisation must be one of {sorted(set(POLARISATION))}, "
                             f"not {self.polarisation!r}")
        self.polarisation = _READABLE[p]

    @property
    def code(self) -> str:
        """the model's polarisation code: `vv`, `vh` or `ma`"""
        return POLARISATION[self.polarisation]


@dataclass(frozen=True)
class Scope:
    """The physics one histogram carries, derived from its sample and pulse.

    * `donor` -- the strength with which this pulse excites the donor of these
      molecules: `0.0` if there is no donor or the pulse misses it, `1.0` if
      the pulse was chosen for it, and otherwise the name of the parameter that
      says how weakly it does so.
    * `acceptor` -- the same for direct excitation of the acceptor.
    * `fret` -- both chromophores are present and the donor is excited, so
      energy transfer happens and the histogram carries the DISTANCE.
    * `sensitised` -- the acceptor lights up through that transfer.

    Two histograms carry the same physics exactly when these four agree.  The
    strengths have to be here and not just a pair of booleans: under
    interleaved excitation both pulses reach both chromophores, and it is only
    the strengths that tell the green pulse from the red one.
    """
    donor: object
    acceptor: object
    fret: bool
    sensitised: bool

    def describe(self) -> str:
        bits = []
        if self.donor:
            bits.append('donor' if self.donor == 1.0 else f'donor ({self.donor})')
        if self.fret:
            bits.append('FRET')
        if self.sensitised:
            bits.append('sensitised acceptor')
        if self.acceptor:
            bits.append('acceptor directly' if self.acceptor == 1.0
                        else f'acceptor directly ({self.acceptor})')
        return ' + '.join(bits) if bits else 'nothing emits'


def scope_of(sample: Sample, excitation: Excitation) -> Scope:
    """THE RULE.  Six lines, in place of a table of six cases written out by
    hand -- and it produces the cases of geometries that were never tabulated."""
    has_d = 'donor' in sample.carries
    has_a = 'acceptor' in sample.carries
    to_d = excitation.strength('donor') if has_d else 0.0
    to_a = excitation.strength('acceptor') if has_a else 0.0
    transfer = bool(has_d and has_a and to_d)
    #: the acceptor sensitised by a donor that was itself only weakly excited
    #: is second order and the model drops it; that is a modelling choice and
    #: it belongs here, where it can be read
    return Scope(donor=to_d, acceptor=to_a, fret=transfer,
                 sensitised=bool(transfer and to_d == 1.0))


# --------------------------------------------------------------------------
# the experiment
# --------------------------------------------------------------------------

@dataclass
class Experiment:
    """Samples, pulses, detectors and the histograms that were recorded."""
    name: str
    samples: list
    excitations: list
    detectors: list
    channels: list
    #: True when the pulses share one laser period, so that one histogram holds
    #: both windows and the model must separate them.  Two lasers used in two
    #: separate acquisitions are NOT interleaved, and the description has to say
    #: which it is: it changes what a histogram is.
    interleaved: bool = False
    notes: str = ''

    # -- lookups -----------------------------------------------------------
    def sample(self, name) -> Sample:
        return self._one(self.samples, name, 'sample')

    def excitation(self, name) -> Excitation:
        return self._one(self.excitations, name, 'excitation')

    def detector(self, name) -> Detector:
        return self._one(self.detectors, name, 'detector')

    @staticmethod
    def _one(items, name, what):
        for it in items:
            if it.name == name:
                return it
        raise KeyError(f'no {what} named {name!r}; have '
                       f'{[i.name for i in items]}')

    # -- validation --------------------------------------------------------
    def validate(self) -> 'Experiment':
        """Refuse a description that cannot be true.  Returns self, so it can
        be chained onto a constructor."""
        for what, items in (('sample', self.samples), ('excitation', self.excitations),
                            ('detector', self.detectors)):
            names = [i.name for i in items]
            dup = {n for n in names if names.count(n) > 1}
            if dup:
                raise ValueError(f'duplicate {what} name(s) {sorted(dup)}')
            if not items:
                raise ValueError(f'the experiment has no {what}s')
        if not self.channels:
            raise ValueError('the experiment records no histograms')
        seen = set()
        for c in self.channels:
            self.sample(c.sample); self.excitation(c.excitation); self.detector(c.detector)
            sig = (c.sample, c.excitation, c.detector, c.polarisation)
            if sig in seen:
                raise ValueError(f'the histogram {sig} is listed twice')
            seen.add(sig)
        if self.interleaved:
            if len(self.excitations) != 2:
                raise ValueError(
                    f'{len(self.excitations)} interleaved pulses: the model addresses a pulse by '
                    'giving each detector a second alias, so it handles exactly two. Three would '
                    'need a third alias and a third window, which is not implemented.')
            if self.excitations[0].delay_ns == self.excitations[1].delay_ns:
                raise ValueError('two interleaved pulses at the same delay are one pulse; '
                                 'give the second its position in the period')
        #: a histogram that carries nothing is a description error, not a
        #: measurement: it says the sample cannot emit under that pulse
        for c in self.channels:
            s = scope_of(self.sample(c.sample), self.excitation(c.excitation))
            if not (s.donor or s.acceptor or s.sensitised):
                raise ValueError(
                    f'the histogram {c.sample}/{c.excitation}/{c.detector} carries nothing: '
                    f'sample {c.sample!r} has {self.sample(c.sample).carries} and pulse '
                    f'{c.excitation!r} reaches {self.excitation(c.excitation).excites}')
        return self

    # -- derived ------------------------------------------------------------
    def pulse_index(self, name) -> int:
        return [e.name for e in self.excitations].index(name)

    def detector_alias(self, detector: str, excitation: str) -> str:
        """The label the model uses for (detector, pulse).

        One physical detector sees both pulses, and the model separates them by
        giving the second pulse an alias of the detector -- `gv` under the
        first pulse, `g2v` under the second.  The alias shares the detector's
        response and shift, which is the point: they are the same hardware.
        """
        if not self.interleaved or self.pulse_index(excitation) == 0:
            return detector
        return detector[0] + '2' + detector[1:]

    def channel_keys(self) -> list:
        """The histograms in the model's own grammar: `(sample, 'det_kind')`."""
        return [(c.sample, f'{self.detector_alias(c.detector, c.excitation)}_{c.code}')
                for c in self.channels]

    def scope_table(self) -> dict:
        """`{channel key: Scope}` -- what each histogram carries."""
        return {k: scope_of(self.sample(c.sample), self.excitation(c.excitation))
                for k, c in zip(self.channel_keys(), self.channels)}

    def detectors_used(self) -> list:
        return sorted({c.detector for c in self.channels})

    def samples_used(self) -> list:
        return sorted({c.sample for c in self.channels})

    def summary(self) -> str:
        """A paragraph a reader can check against their own setup."""
        out = [f'{self.name}: {len(self.channels)} histograms']
        out.append('  samples     ' + ', '.join(
            f'{s.name} ({"+".join(s.carries)})' for s in self.samples))
        out.append('  excitations ' + ', '.join(
            f'{e.name} at {e.delay_ns:g} ns, reaches '
            + ', '.join(f'{c} x {v}' for c, v in e.excites.items())
            for e in self.excitations))
        out.append('  detectors   ' + ', '.join(f'{d.name} ({d.colour})' for d in self.detectors))
        out.append('  histograms')
        for c, (k, s) in zip(self.channels, self.scope_table().items()):
            out.append(f'    {k[0]:<3} {k[1]:<8} {c.polarisation:<13} {s.describe()}')
        if self.notes:
            out.append('  ' + self.notes)
        return '\n'.join(out)

    # -- JSON ---------------------------------------------------------------
    def to_dict(self) -> dict:
        return dict(name=self.name, notes=self.notes, interleaved=bool(self.interleaved),
                    samples=[asdict(s) for s in self.samples],
                    excitations=[asdict(e) for e in self.excitations],
                    detectors=[asdict(d) for d in self.detectors],
                    channels=[asdict(c) for c in self.channels])

    @classmethod
    def from_dict(cls, d: dict) -> 'Experiment':
        return cls(name=d['name'], notes=d.get('notes', ''),
                   interleaved=bool(d.get('interleaved', False)),
                   samples=[Sample(**s) for s in d['samples']],
                   excitations=[Excitation(**e) for e in d['excitations']],
                   detectors=[Detector(**x) for x in d['detectors']],
                   channels=[Channel(**c) for c in d['channels']]).validate()

    def to_json(self, path=None, indent=2) -> str:
        text = json.dumps(self.to_dict(), indent=indent)
        if path is not None:
            Path(path).write_text(text)
        return text

    @classmethod
    def from_json(cls, path_or_text) -> 'Experiment':
        text = str(path_or_text)
        if not text.lstrip().startswith('{'):          # a path, not the JSON itself
            text = Path(text).read_text()
        return cls.from_dict(json.loads(text))


# --------------------------------------------------------------------------
# the geometries, as descriptions rather than as code paths
# --------------------------------------------------------------------------

#: the four detectors of a multiparameter (MFD) setup: two colours, each split
#: by a polarising beam splitter, each calibrated on its own
_MFD_DETECTORS = [Detector('gv', 'green'), Detector('gh', 'green'),
                  Detector('rv', 'red'), Detector('rh', 'red')]
_MFD_POL = {'gv': 'parallel', 'gh': 'perpendicular',
            'rv': 'parallel', 'rh': 'perpendicular'}

#: the green pulse reaches the acceptor a little; that path is what the model
#: calls direct excitation, and it is measurable rather than assumed away
_GREEN = Excitation('green', {'donor': 1.0, 'acceptor': 'EX_AG'}, 0.0)
#: and the red pulse reaches the donor a little, which is what makes the
#: acceptor's direct excitation identifiable when both pulses are present
_RED = Excitation('red', {'acceptor': 1.0, 'donor': 'EX_DR'}, 25.0)


def mfd_donor_excitation() -> Experiment:
    """One laser, four detectors, the labelled sample alone.

    The commonest single-molecule FRET measurement, and the hardest case for
    this model: there is no donor-only reference, so the donor's unquenched
    lifetime spectrum has to come out of the same histograms that carry the
    distance.
    """
    return Experiment(
        name='MFD, donor excitation only',
        samples=[Sample('DA', ('donor', 'acceptor'))],
        excitations=[_GREEN],
        detectors=list(_MFD_DETECTORS),
        channels=[Channel('DA', 'green', d.name, _MFD_POL[d.name]) for d in _MFD_DETECTORS],
        notes='no donor-only and no acceptor-only reference: see the identifiability table',
    ).validate()


def mfd_pie(full: bool = True) -> Experiment:
    """Two lasers interleaved in the period, three samples, four detectors.

    `full=True` describes every (sample, pulse) as its own physics channel --
    twenty-four of them in twelve histograms, which is what the measurement
    is.  `full=False` is the earlier route in which only the labelled sample's
    red-pulse window is modelled.
    """
    samples = [Sample('D0', ('donor',)), Sample('A0', ('acceptor',)),
               Sample('DA', ('donor', 'acceptor'))]
    chans = []
    if full:
        for s in samples:
            for d in _MFD_DETECTORS:
                for e in ('green', 'red'):
                    chans.append(Channel(s.name, e, d.name, _MFD_POL[d.name]))
    else:
        for d in _MFD_DETECTORS:
            chans.append(Channel('D0', 'green', d.name, _MFD_POL[d.name]))
        for d in _MFD_DETECTORS:
            chans.append(Channel('DA', 'green', d.name, _MFD_POL[d.name]))
        for d in _MFD_DETECTORS:
            chans.append(Channel('A0', 'red', d.name, _MFD_POL[d.name]))
        samples = [Sample('D0', ('donor',)), Sample('DA', ('donor', 'acceptor')),
                   Sample('A0', ('acceptor',))]
    return Experiment(
        name='MFD-PIE' + ('' if full else ' (labelled sample only under the red pulse)'),
        samples=samples, excitations=[_GREEN, _RED],
        detectors=list(_MFD_DETECTORS), channels=chans, interleaved=True,
        notes='the two windows of one sample share its concentration and counting time, '
              'so their ratio measures the excitation crosstalk',
    ).validate()


def separate_measurements() -> Experiment:
    """Three cuvettes measured on their own: donor-only, labelled,
    acceptor-only.

    Two detectors, each recorded parallel and perpendicular through the same
    response.  The acceptor-only sample is excited by its own laser, in its own
    acquisition -- which is why its pulse is named `red` although no red pulse
    is interleaved with anything.
    """
    dets = [Detector('g', 'green'), Detector('r', 'red')]
    ch = []
    for pol in ('parallel', 'perpendicular'):
        ch.append(Channel('D0', 'green', 'g', pol))
    for pol in ('parallel', 'perpendicular'):
        ch.append(Channel('DA', 'green', 'g', pol))
    for pol in ('parallel', 'perpendicular'):
        ch.append(Channel('DA', 'green', 'r', pol))
    for pol in ('parallel', 'perpendicular'):
        ch.append(Channel('A0', 'acceptor laser', 'r', pol))
    return Experiment(
        name='separate donor, FRET and acceptor measurements',
        samples=[Sample('D0', ('donor',)), Sample('DA', ('donor', 'acceptor')),
                 Sample('A0', ('acceptor',))],
        #: two lasers, but in two SEPARATE acquisitions rather than interleaved
        #: in one period -- the acceptor-only cuvette is measured on its own,
        #: with its own laser at full strength
        excitations=[Excitation('green', {'donor': 1.0, 'acceptor': 'EX_AG'}, 0.0),
                     Excitation('acceptor laser', {'acceptor': 1.0}, 0.0)],
        detectors=dets, channels=ch, interleaved=False,
        notes="the donor-only sample's red channel is left out: it constrains the "
              'crosstalk, which is taken as calibrated here',
    ).validate()


def magic_angle_minimal() -> Experiment:
    """Two histograms and no polarisation: the smallest experiment this model
    accepts.  At the magic angle the anisotropy term vanishes, so the
    rotational parameters are not measured -- they are returned from the
    prior."""
    dets = [Detector('g', 'green'), Detector('r', 'red')]
    return Experiment(
        name='magic angle, labelled sample only',
        samples=[Sample('DA', ('donor', 'acceptor'))],
        excitations=[Excitation('green', {'donor': 1.0, 'acceptor': 'EX_AG'}, 0.0)],
        detectors=dets,
        channels=[Channel('DA', 'green', 'g', 'magic'), Channel('DA', 'green', 'r', 'magic')],
        notes='no polarisation: the anisotropy is not measured here',
    ).validate()


GEOMETRIES = {
    'mfd_donor_excitation': mfd_donor_excitation,
    'mfd_pie': mfd_pie,
    'separate_measurements': separate_measurements,
    'magic_angle_minimal': magic_angle_minimal,
}


# --------------------------------------------------------------------------
# from the description to a model that can be fitted
# --------------------------------------------------------------------------

def build(spec: Experiment, knots: int = 25, threads: int = 4) -> dict:
    """The environment the description implies: the distance grid, the spline
    basis, the time axis, and -- when the pulses are interleaved -- the second
    window and the detector aliases that address it."""
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from bd import P                                       # noqa: E402
    import numpy as np
    L = P.load_prototype(threads=threads)
    import os
    d = P.prototype_dir(); cwd = Path.cwd(); os.chdir(d)
    try:
        E, rel, edges, spl, dx = L.environment(knots, instrument=0, homogeneous=True)
    finally:
        os.chdir(cwd)
    if spec.interleaved:
        Eu = L.pie_environment(E)
        Eu['pie_full'] = True
    else:
        Eu = E
    #: `Ep` is the name the older helpers in this folder use for the
    #: environment the model is actually fitted in; kept as an alias so
    #: `bff_forward.structural_graph` and `bd`'s plotting take this model too
    return dict(L=L, E=E, Eu=Eu, Ep=Eu, rel=np.asarray(rel), spl=spl, edges=edges,
                spec=spec, keys=spec.channel_keys(), n_coef=int(spl.shape[1]),
                n_bin=int(Eu['n_bin']), dt=float(E['dt']), period=float(E['period']))


def reference_channel(spec: Experiment, sample: str):
    """The histogram a sample's acquisition scale is quoted in: its first, in
    the order the description lists them.  Which one it is does not matter --
    the scale is defined by it -- but it has to be SAID, or two runs of the
    same description mean different things."""
    for k, c in zip(spec.channel_keys(), spec.channels):
        if c.sample == sample:
            return k
    raise KeyError(f'sample {sample!r} records no histogram')


def truth_values(spec: Experiment, model: dict, p_true, photons=3e5, scat=0.02, bkg=0.01):
    """Every constant of the forward model at its nominal value, the scales set
    so that each sample's reference channel collects `photons`.

    This is the prototype's `truth_values` with the reference channel taken
    from the description instead of being written in, which is the whole
    difference: the prototype's version names `('A0', 'r_vv')` and `g_vv`, and
    neither exists in three of the four geometries here.
    """
    L = model['L']; Eu = model['Eu']; keys = model['keys']
    S = L.S; tt = L.tt
    vals = dict(x_d0=tt(S.X_D0_TRUE), w_a=tt(S.W_A_TRUE),
                w_rho=tt(L.log_bump(L.M.RHO_GRID, S.RHO_TRUE, 0.25)),
                w_rho_a=tt(L.log_bump(L.M.RHO_A_GRID, S.RHO_TRUE, 0.25)),
                r0_d=tt(S.R0_ANISO), r0_a=tt(0.5 * S.R0_ANISO), g=tt(S.G_TRUE),
                l1=tt(S.L1), l2=tt(S.L2), C_GD=tt(S.C_GD), C_GA=tt(S.C_GA),
                C_RD=tt(S.C_RD), C_RA=tt(S.C_RA), G_GREEN=tt(S.G_GREEN), G_RED=tt(S.G_RED),
                QY_D=tt(S.QY_D), QY_A=tt(S.QY_A), EX_AG=tt(S.EX_AG))
    if spec.interleaved:
        vals['g_r'] = tt(S.G_TRUE); vals['EX_DR'] = tt(0.002)
    for d in sorted({L.parse_channel(k)[1] for k in keys}):
        vals[f'irf_shift_{d}'] = tt(Eu['instrument'][2])
        vals[f'irf_width_{d}'] = tt(Eu['instrument'][1])
        vals[f'irf_skew_{d}'] = tt(Eu['instrument'][3])
    spec_d = Eu['cD']
    raw = L.physics_amplitudes(Eu, vals, spec_d, tt(p_true), keys)
    import math
    for samp in spec.samples_used():
        ref = reference_channel(spec, samp)
        vals[f'log_scale_{samp}'] = tt(math.log(photons / max(float(raw[ref].sum()), 1e-30)))
    for k in keys:
        kk = f'{k[0]}_{k[1]}'
        vals[f'scat_{kk}'] = tt(scat); vals[f'bkg_{kk}'] = tt(bkg)
    return vals, spec_d


def graph(spec: Experiment, model: dict, y=None, photons=3e5, irf_shape=False):
    """The factor graph the description implies.

    Which variables exist is not a choice made here: `default_variables` reads
    the channel list and adds a scatter and a background per histogram, a
    response and a shift per detector actually used, and one acquisition scale
    per sample.  The description decides the channel list; everything else
    follows.
    """
    import torch
    L = model['L']; Eu = dict(model['Eu']); keys = model['keys']; spl = model['spl']
    n_coef = model['n_coef']
    med = {}
    if y is not None:
        for samp in spec.samples_used():
            med[samp] = max(float(y[reference_channel(spec, samp)].sum()), 1.0)
    V = L.default_variables(Eu, keys, n_coef=n_coef, photons=photons,
                            scale_medians=(med or None), irf_shape=irf_shape)
    ps = L.PSplineFactor(n_coef, spl=spl)
    if y is None:
        y = {k: torch.ones(model['n_bin']) for k in keys}
    g = L.FactorGraph(Eu, keys, V, L.PoissonCountsFactor({k: y[k] for k in keys}),
                      L.InstrumentModel(Eu, 'analytic'), spl, ps)
    g.rel = model['rel']
    return g


def simulate(spec: Experiment, model: dict, p_true, photons=3e5, seed=0, scat=0.02, bkg=0.01):
    """One Poisson realisation of exactly the histograms the description lists."""
    import torch
    L = model['L']
    g = graph(spec, model, None, photons=photons)
    vals, spec_d = truth_values(spec, model, p_true, photons, scat, bkg)
    a2 = {k: L.instrument_amplitudes(vals, L.physics_amplitudes(
        model['Eu'], vals, spec_d, L.tt(p_true), model['keys'])[k], k) for k in model['keys']}
    lam = g.expected_counts(vals, a2)
    gen = torch.Generator().manual_seed(int(seed))
    y = {k: torch.poisson(lam[k], generator=gen) for k in model['keys']}
    return y, lam, vals


def fit(spec: Experiment, model: dict, y, lam_nodes=(1.0, 0.0, -1.0), seed=0,
        verbose=False, accelerate=True):
    """The Laplace posterior of the whole model on exactly these histograms.

    The roughness weight of the P-spline prior on log p(R/R0) is not maximised
    but integrated out: a Laplace approximation at each node of a grid, the
    nodes mixed by their evidences (INLA; Rue, Martino & Chopin 2009).
    """
    import torch
    L = model['L']
    g = graph(spec, model, y)
    if accelerate:
        try:
            from fast_forward import accelerate as _acc
            _acc(g)
        except Exception:
            pass
    g.start_y = dict(y)
    gen = torch.Generator().manual_seed(int(seed))
    post = L.fit_sample(g, y, gen, model['rel'], verbose=verbose, start='mem',
                        optimiser='fisher', hessian='fisher',
                        lam_nodes=tuple(float(x) for x in lam_nodes))
    return g, post


def rule0(model, g, post, y):
    """Rule 0, per histogram: the Poisson deviance per degree of freedom and a
    runs test on the weighted residuals.  A number about a fitted quantity is
    not reportable without them."""
    rows = [dict(channel=f'{k[0]} {k[1]}', counts=float(y[k].sum()),
                 dpd=r['dpd'], runs_p=r['runs_p']) for k, r in post['rows'].items()]
    rows.append(dict(channel='all', counts=float(sum(float(y[k].sum()) for k in g.data_keys)),
                     dpd=post['dev'] / post['dof'], runs_p=float('nan')))
    return rows


# --------------------------------------------------------------------------
# what the experiment can and cannot determine
# --------------------------------------------------------------------------

def identifiability(spec: Experiment, model: dict, p_true=None, photons=3e5, step=1e-4):
    """For every unknown: how much of what the posterior will say about it
    comes from the DATA rather than from the prior, and which histograms carry
    that information.

    The measure is Fisher's, in the coordinate the variable is optimised in.
    Moving a coordinate by `dz` changes the expected counts of histogram `k` by
    `dlam`, and a Poisson measurement notices that to the extent of
    `sum_b dlam_b^2 / lam_b`.  Dividing by `dz^2` gives the information the
    histogram carries about the coordinate; adding the prior's own precision
    gives the posterior's, and the ratio is the data's share.

    **What this is not.** It is the DIAGONAL, so it answers "could this
    histogram see this coordinate move on its own", not "can the fit separate
    it from everything else".  A coordinate with a large share here can still
    be undetermined through a correlation, which is what the full covariance is
    for.  A share near zero, though, is conclusive: nothing in the data moves
    when the coordinate does, so whatever the fit reports for it is the prior.
    """
    import numpy as np, torch
    L = model['L']
    if p_true is None:
        p_true = np.exp(-0.5 * ((model['rel'] - 0.9) / 0.08) ** 2)
        p_true = p_true / p_true.sum()
    g = graph(spec, model, None, photons=photons)
    vals, spec_d = truth_values(spec, model, p_true, photons)
    #: the values at the truth, in the coordinates the fit works in
    zd = L.prior_median_z(g)
    for v in g.free:
        if v.name in vals:
            zd[v.name] = v.transform.to_unconstrained(L.tt(vals[v.name]).reshape(-1))
    zd['c'] = g.index['c'].transform.to_unconstrained(
        L.coefficients_of(g, p_true, float(vals['x_d0'])).reshape(-1))
    th = g.pack(zd)
    keys = model['keys']

    def counts(theta):
        v, _ = g.unpack(theta)
        return {k: g.expected_counts(v)[k].detach().numpy() for k in keys}

    base = counts(th)
    prior_prec = {}
    for var in g.free:
        sd = getattr(var.prior, 'sd', None)
        if sd is None:
            lo = getattr(var.prior, 'lo', None); hi = getattr(var.prior, 'hi', None)
            sd = (float(hi) - float(lo)) / math.sqrt(12.0) if lo is not None else 1.0
        prior_prec[var.name] = 1.0 / np.asarray(np.broadcast_to(np.asarray(sd, float),
                                                               (var.size,))) ** 2
    rows = []
    for var in g.free:
        info = np.zeros((var.size, len(keys)))
        a0 = g.offsets[var.name][0]
        for i in range(var.size):
            t2 = th.clone(); h = step * max(1.0, abs(float(th[a0 + i])))
            t2[a0 + i] = th[a0 + i] + h
            pert = counts(t2)
            for j, k in enumerate(keys):
                d = (pert[k] - base[k]) / h
                lam = np.maximum(base[k], 1e-12)
                info[i, j] = float((d * d / lam).sum())
        tot = info.sum(axis=1)
        share = tot / (tot + prior_prec[var.name])
        best = int(np.argmax(info.sum(axis=0))) if len(keys) else -1
        carried = [f'{keys[j][0]} {keys[j][1]}' for j in np.argsort(-info.sum(axis=0))[:3]
                   if info.sum(axis=0)[j] > 0.01 * max(info.sum(), 1e-30)]
        rows.append(dict(variable=var.name, size=var.size, role=var.group,
                         data_share=float(share.max()), data_share_min=float(share.min()),
                         information=float(tot.max()),
                         carried_by=', '.join(carried) if carried else 'nothing'))
    return rows


def factor_scopes(spec: Experiment, model: dict, p_true=None, photons=3e5,
                  threshold=1e-6) -> dict:
    """Which unknowns each histogram's likelihood factor actually touches --
    MEASURED, by moving each coordinate and seeing whether that histogram's
    expected counts move.

    The hand-written version of this list is in `bff_forward.structural_graph`
    and it is a reasonable guess. A guess is the wrong thing to build a graph
    decomposition on, though: the elimination order, the treewidth and the
    separators are all statements about which variables share a factor, so a
    scope that is wrong by one variable is a decomposition that is wrong.
    """
    import numpy as np
    rows = _per_channel_information(spec, model, p_true, photons)
    out = {}
    for k, info in rows.items():
        tot = max(sum(info.values()), 1e-300)
        out[k] = sorted(n for n, v in info.items() if v > threshold * tot)
    return out


def _per_channel_information(spec, model, p_true=None, photons=3e5, step=1e-4):
    """`{channel: {variable: Fisher information this channel carries}}`."""
    import numpy as np
    L = model['L']
    if p_true is None:
        p_true = np.exp(-0.5 * ((model['rel'] - 0.9) / 0.08) ** 2)
        p_true = p_true / p_true.sum()
    g = graph(spec, model, None, photons=photons)
    vals, _ = truth_values(spec, model, p_true, photons)
    zd = L.prior_median_z(g)
    for v in g.free:
        if v.name in vals:
            zd[v.name] = v.transform.to_unconstrained(L.tt(vals[v.name]).reshape(-1))
    zd['c'] = g.index['c'].transform.to_unconstrained(
        L.coefficients_of(g, p_true, float(vals['x_d0'])).reshape(-1))
    th = g.pack(zd)
    keys = model['keys']

    def counts(theta):
        v, _ = g.unpack(theta)
        c = g.expected_counts(v)
        return {k: c[k].detach().numpy() for k in keys}

    base = counts(th)
    out = {k: {} for k in keys}
    for var in g.free:
        a0 = g.offsets[var.name][0]
        acc = {k: 0.0 for k in keys}
        for i in range(var.size):
            h = step * max(1.0, abs(float(th[a0 + i])))
            t2 = th.clone(); t2[a0 + i] = th[a0 + i] + h
            pert = counts(t2)
            for k in keys:
                d = (pert[k] - base[k]) / h
                acc[k] += float((d * d / np.maximum(base[k], 1e-12)).sum())
        for k in keys:
            out[k][var.name] = acc[k]
    return out


def structural_graph(spec: Experiment, model: dict, scopes=None, p_true=None, photons=3e5):
    """The model as an `IMP.bff.FactorGraph`: which unknowns exist, which
    factors touch which of them, and what that implies for how the posterior
    can be decomposed.

    This is the graph in the graphical-model sense and it knows nothing about
    numbers -- it knows that the spline coefficients and the donor-only
    fraction appear together in every likelihood factor and therefore cannot be
    updated independently. `scopes` defaults to the MEASURED ones.
    """
    import IMP.bff as bff
    g = graph(spec, model, None, photons=photons)
    scopes = factor_scopes(spec, model, p_true, photons) if scopes is None else scopes
    fg = bff.FactorGraph()
    for i, name in enumerate(sorted(g.offsets, key=lambda n: g.offsets[n][0])):
        fg.add_variable(name, name, i)
    for name in sorted(g.offsets):
        fg.add_factor(f'prior_{name}', 0, [name])          # 0 = a prior
    for k in g.data_keys:
        sc = [n for n in scopes.get(k, []) if n in g.offsets]
        if sc:
            fg.add_factor(f'counts_{k[0]}_{k[1]}', 1, sorted(sc))   # 1 = a likelihood
    return fg, scopes


def unidentified(rows, threshold=0.05):
    """the unknowns whose posterior is the prior, whatever the fit prints"""
    return [r for r in rows if r['data_share'] < threshold]


# --------------------------------------------------------------------------
# the gate: the rule must reproduce the table it replaces
# --------------------------------------------------------------------------

def check_against_prototype(verbose=True) -> dict:
    """E2 of the PRD.  The derived channel keys must equal the hand-written
    lists, and the derived scopes must partition the channels exactly as the
    hand-written `channel_scope` does.

    The partition is what is compared, not the names: the rule does not know
    the prototype's six labels and is not asked to invent them.  Two channels
    must land in the same derived scope if and only if the prototype gives them
    the same name -- which is a check that can fail, and did while the alias
    rule was wrong.
    """
    import sys, os
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from bd import P                                     # noqa: E402
    L = P.load_prototype(threads=1)

    report = {}

    def partition(keys, scopes, ref):
        got, want = {}, {}
        for k in keys:
            got.setdefault(scopes[k], []).append(k)
            want.setdefault(ref[k], []).append(k)
        return (sorted(sorted(v) for v in got.values()),
                sorted(sorted(v) for v in want.values()))

    # -- the interleaved geometry, all twenty-four physics channels
    exp = mfd_pie(full=True)
    keys = exp.channel_keys()
    ref_keys = L.pie_full_keys()
    assert sorted(keys) == sorted(ref_keys), (
        f'derived keys differ from pie_full_keys():\n  only derived: '
        f'{sorted(set(keys) - set(ref_keys))}\n  only prototype: {sorted(set(ref_keys) - set(keys))}')
    ref = {k: L.channel_scope(k, {'pie_full': True}) for k in keys}
    got, want = partition(keys, exp.scope_table(), ref)
    assert got == want, f'the derived scopes partition differently:\n{got}\n{want}'
    assert len(set(ref.values())) == 6, 'the interleaved geometry should have six scopes'
    report['mfd_pie(full=True)'] = dict(channels=len(keys), scopes=len(got))

    # -- the earlier interleaved route
    exp = mfd_pie(full=False)
    assert sorted(exp.channel_keys()) == sorted(L.pie_keys()), 'pie_keys() not reproduced'
    report['mfd_pie(full=False)'] = dict(channels=len(exp.channels))

    # -- three cuvettes: the ensemble channel list
    exp = separate_measurements()
    assert sorted(exp.channel_keys()) == sorted(L.CHANNELS_FULL), (
        f'CHANNELS_FULL not reproduced: {sorted(exp.channel_keys())}')
    ref = {k: L.channel_scope(k, {}) for k in exp.channel_keys()}
    got, want = partition(exp.channel_keys(), exp.scope_table(), ref)
    assert got == want, 'the ensemble scopes partition differently'
    report['separate_measurements'] = dict(channels=len(exp.channels), scopes=len(got))

    # -- the two geometries that have no hand-written list to compare against
    for nm in ('mfd_donor_excitation', 'magic_angle_minimal'):
        exp = GEOMETRIES[nm]()
        report[nm] = dict(channels=len(exp.channels),
                          scopes=len(set(exp.scope_table().values())))

    # -- the JSON round trip, every geometry
    for nm, make in GEOMETRIES.items():
        a = make()
        b = Experiment.from_json(a.to_json())
        assert a.to_dict() == b.to_dict(), f'{nm} does not survive JSON'
        assert a.channel_keys() == b.channel_keys(), f'{nm} keys change through JSON'
    report['json'] = dict(round_trip='exact for every geometry')

    if verbose:
        for k, v in report.items():
            print(f'  {k:<26} {v}')
        print('PASS: the derived channels and scopes reproduce the hand-written ones')
    return report


def check_refusals(verbose=True) -> list:
    """The four descriptions that cannot be true, and the messages that name
    the offending item.  A validator that never refuses anything is not a
    validator."""
    cases = []

    def bad(what, make):
        try:
            make()
        except (ValueError, KeyError) as e:
            cases.append((what, str(e).split('\n')[0]))
            return
        raise AssertionError(f'{what}: was accepted and should not have been')

    bad('a channel naming a detector that does not exist', lambda: Experiment(
        'x', [Sample('DA', ('donor', 'acceptor'))], [_GREEN], [Detector('g', 'green')],
        [Channel('DA', 'green', 'nope', 'magic')]).validate())
    bad('two detectors with the same name', lambda: Experiment(
        'x', [Sample('DA', ('donor', 'acceptor'))], [_GREEN],
        [Detector('g', 'green'), Detector('g', 'red')],
        [Channel('DA', 'green', 'g', 'magic')]).validate())
    bad('an excitation that excites nothing', lambda: Excitation('dark', {}))
    bad('a histogram that carries nothing', lambda: Experiment(
        'x', [Sample('A0', ('acceptor',))], [Excitation('green', {'donor': 1.0})],
        [Detector('g', 'green')], [Channel('A0', 'green', 'g', 'magic')]).validate())
    bad('three interleaved pulses', lambda: Experiment(
        'x', [Sample('DA', ('donor', 'acceptor'))],
        [_GREEN, _RED, Excitation('ir', {'acceptor': 1.0}, 40.0)],
        [Detector('g', 'green')], [Channel('DA', 'green', 'g', 'magic')],
        interleaved=True).validate())

    if verbose:
        for what, msg in cases:
            print(f'  refused: {what}\n           {msg}')
    return cases


if __name__ == '__main__':
    print('the geometries')
    for nm, make in GEOMETRIES.items():
        print()
        print(make().summary())
    print()
    print('the gate')
    check_against_prototype()
    print()
    print('what the validator refuses')
    check_refusals()
