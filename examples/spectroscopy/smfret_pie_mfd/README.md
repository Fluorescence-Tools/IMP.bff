# smFRET PIE-MFD: from a photon stream to a distance distribution

*2026-09-08*

`smfret_pie_mfd.ipynb` simulates a single-molecule FRET experiment with pulsed
interleaved excitation and multi-parameter detection as a **photon stream**,
finds the bursts with `tttrlib`, builds the fluorescence decays out of the burst
photons, and infers the donor–acceptor distance distribution p(R/R0) from them —
analytically, with no sampler.

Open the notebook and edit its first cell; everything else follows from it.

## Files

| file | what it is |
|---|---|
| `smfret_pie_mfd.ipynb` | the example, in order, with its figures |
| `pie_mfd.py` | the functions the notebook calls: the simulation, the burst search, the decays, the fit, the figures |

## What it needs

* `tttrlib` ≥ 0.27 — the photon stream object and the burst search
* the ucfret prototype `s88_laplace_posterior.py`, which carries the forward
  model and the Laplace posterior. It is looked for at
  `~/dev/ucfret/investigation/pinn_pR_anisotropy`; set `UCFRET_S88` to point
  elsewhere. Its physics maps are cached beside it on first use.
* `numpy`, `scipy`, `torch`, `matplotlib`

## The measurement, in one paragraph

Molecules with a donor and an acceptor dye cross a confocal spot one at a time.
Two lasers alternate inside one 50 ns period: the green pulse excites the donor,
and 25 ns later the red pulse excites the acceptor alone — which is what
distinguishes a molecule that transfers nothing from one whose acceptor is dead.
Four detectors record green and red, each in two polarisations, so the dyes'
rotation is measured alongside the transfer. Every photon carries a macro time
(which pulse, hence when) and a micro time (its delay after that pulse, ~30 ps
resolution). Micro times histogrammed give the decays; macro times give the
bursts.

## What is inferred

Eight histograms — four detectors on a donor-only reference sample, four on the
labelled sample, the latter holding both pulses — give, in one fit: p(R/R0), the
fraction of molecules with no acceptor, both anisotropies and rotational times,
the donor's lifetime spectrum, the g factor, the crosstalks, the quantum yields,
the direct excitation, the response shifts, and a scale, a scatter fraction and a
background per histogram.

The roughness weight of the spline prior on log p(R/R0) is not maximised — that
flattens the distribution — but integrated out: a Laplace approximation at each
node of a grid over the weight, the nodes mixed by their evidences (INLA; Rue,
Martino & Chopin, *JRSS-B* **71**, 319, 2009).

## What the notebook checks, and what it does not

It reports weighted residuals under every histogram with a runs test, and the
Poisson deviance per degree of freedom against a **measured** reference — Poisson
draws at the fitted means, because for counts of this size the deviance of a
correct model is not the number of degrees of freedom. It then repeats the whole
experiment at a few seeds and plots the pulls, which is the only thing that says
whether the intervals mean what they claim.

It is not a calibration study. A single run's coverage is not evidence, and the
notebook says so where it reports one.

## Why simulation and fit share a forward model

Both call the same prototype. The example therefore measures the inference, not a
disagreement between two hand-written physics implementations. What the
simulation adds, and the fit never learns, is the experiment itself: molecules
drawn one at a time, Poisson burst sizes, uncorrelated background, scattered
laser light, and a burst search that decides which photons are used at all — a
selection that, as the notebook shows, moves the distribution the decays measure
away from the one in the cuvette.
