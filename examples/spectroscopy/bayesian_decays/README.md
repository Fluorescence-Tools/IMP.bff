# Bayesian fluorescence decays — how the model was arrived at

*2026-09-08, extended 2026-09-09*

Eleven notebooks that build one model, in the order it was actually built,
five that show how to point it at an experiment of your own, and three on what a
single-molecule measurement carries that a pooled decay does not. The end of
the road is `../smfret_pie_mfd`, which is the whole thing as a single runnable
example; this series exists so that nothing in it looks arbitrary.

Each notebook states the mathematics of the step it adds, runs it, and shows
the figure. Several also record what was tried and **abandoned**, because a
model is only understandable together with the alternatives that were rejected
and the measurement that rejected them.

**`THEORY.md`** beside this file is the probability calculus in full — every
symbol, every prior, every approximation with what it costs. The notebooks
point into it by section.

## The order

| | notebook | what it establishes |
|---|---|---|
| 0 | `00_what_came_before.ipynb` | the amortised neural network this project tried first: what it measured, what it got right, and the calibration gate it never passed — which is why the rest of the series exists |
| 1 | `01_the_measurement.ipynb` | what the data are; the forward model; the ruler $E(r)=1/(1+r^6)$ and the window where it stops reading; that the expected counts are **linear in the distance distribution** |
| 2 | `02_the_inverse_problem.ipynb` | that inverting a decay for $p(R)$ is ill posed even with every nuisance known exactly — 15 of 128 directions above the noise — and what a roughness prior is for |
| 3 | `03_the_penalty_weight.ipynb` | that the prior's weight **cannot** be fitted jointly (the joint mode flattens two populations into one), and the evidence mixture that replaces it |
| 4 | `04_the_posterior.ipynb` | why drawing from the posterior and pushing it through the softmax overshoots, and the delta method that replaces it; the band's space; the gauge |
| 5 | `05_calibration_and_limits.ipynb` | whether the intervals mean what they say; the deviance's true reference; the bias and the five explanations excluded by measurement; the resolution rule |
| 6 | `06_the_experiment.ipynb` | that the bias is a property of the measurement, and what a second laser buys — computed from the information alone, then measured |
| 7 | `07_single_molecules.ipynb` | the photon stream, `tttrlib` bursts, and an A/B showing the single-molecule experiment costs the inference almost nothing |
| 8 | `08_classify_then_pool.ipynb` | classify the bursts by what each one looks like, pool the groups, analyse each — and separate populations the pooled analysis cannot |
| 9 | `09_what_failed.ipynb` | every mistake serious enough to change a result, with the measurement that caught it — wrong answers that fit, silent losses, and checks that could not fail |
| 10 | `10_the_factor_graph_in_bff.ipynb` | the model as a factor graph in bff, drawn with `networkx`, and made fast: the forward model 38x, the whole objective 8x, agreeing to 1e-15 |
| 11 | `11_configuring_the_network.ipynb` | **the experiment as data**: samples, pulses, detectors and the histograms recorded, in Python and in JSON, from which the channels, the physics, the unknowns and the factor graph are derived |
| 12 | `12_mfd_donor_excitation.ipynb` | one laser, four detectors, the labelled sample alone — the commonest smFRET measurement and the hardest case for this model |
| 13 | `13_mfd_pie.ipynb` | two lasers interleaved, three samples, twelve histograms, and what the second pulse makes measurable |
| 14 | `14_separate_measurements.ipynb` | three cuvettes measured on their own — the ensemble geometry the first ten notebooks use, now as a description |
| 15 | `15_magic_angle_minimal.ipynb` | two histograms and no polarisation: the smallest experiment the model accepts, and the seven numbers it must not be asked for |
| 16 | `16_the_two_dimensional_diagrams.ipynb` | **photon sorting**: $E^*$ against the mean donor arrival time with its static and dynamic lines, stoichiometry against $E^*$, and how much of the scatter is arithmetic |
| 17 | `17_exchange_between_states.ipynb` | molecules that interconvert, over five exchange times: where burst sorting recovers both states, and where it returns a confident answer at a distance no molecule ever had |
| 18 | `18_how_close_is_too_close.ipynb` | the two states brought together until the diagram, the clustering and the sub-ensemble fits each stop separating them, against the pooled analysis on the same photons |

Notebooks 0, 1, 2, 4, 6–10 run in a few minutes each; 3 and 5 do several fits
and take longer. Every one was executed end to end before being committed, and
carries its outputs.

**Notebooks 16–18 are the single-molecule analysis proper.** 1–15 fit pooled
histograms, which is what a cuvette gives you; 16–18 sort the photons burst by
burst first, and are mostly about establishing where that stops working. They
need `single_molecule.py`.

**Notebooks 11–15 are the ones to read if you have your own measurement.** 11 is
about describing it; 12–15 each take one description, say what it can and cannot
determine *before* fitting, simulate data from it, fit, and report Rule 0
against a measured deviance reference. The four differ only in their
description.

## What it needs

* `tttrlib` ≥ 0.27 — the photon stream, the burst search, and the HDBSCAN used
  in notebook 8
* the ucfret prototype `s88_laplace_posterior.py`, which carries the physics
  and the posterior. Looked for at
  `~/dev/ucfret/investigation/pinn_pR_anisotropy`; set `UCFRET_S88` to point
  elsewhere.
* `IMP.bff` for notebook 10
* `numpy`, `scipy`, `torch`, `matplotlib`, `networkx`, and `scikit-learn` for
  the cross-checks in notebooks 8 and 9

`single_molecule.py` is the sorting layer: the two per-burst rulers, the static
and dynamic reference lines taken from the model's own tables, a two-state
Markov exchange in the burst simulation with the gate that it reproduces the
static one, and the sub-ensemble fit of a group of bursts.

`experiment.py` is the description layer: the four dataclasses, the rule that
derives what each histogram carries, the JSON round trip, the identifiability
measurement and the bff factor graph. Running it directly
(`python experiment.py`) prints every geometry, checks the derived channels and
scopes against the literals in the prototype, and shows what the validator
refuses.

`bd.py` is the shared code. It imports `pie_mfd` from `../smfret_pie_mfd` so
that the series and the finished example run the same implementation.
`bff_forward.py` is notebook 10's: the decay model as `IMP.bff` nodes, the
structural factor graph, and the drawing.

## The short version of the argument

A fluorescence decay measures a distance through the transfer efficiency
$E=1/(1+(R/R_0)^6)$, which is smooth and saturating. So a decay determines a
*smooth functional* of the distance distribution and is blind to the rest —
measurably so: of 128 grid points, about fifteen directions move the data by
more than the noise. Any method that returns a distribution is therefore
returning its prior wherever the data are silent, and the only honest thing to
do is to say which prior, why, and how much of the answer is it.

That is what this series is: the prior, its weight treated as a
hyperparameter and integrated out rather than chosen, the posterior
approximated without a sampler, and the whole thing put through checks that
could have failed — a fit statistic against a *measured* reference, residual
runs tests, rank calibration, and a resolution test that says how close two
populations may be before the answer stops being about them.

## Attribution

Written for the bff examples, 2026-09-08. The forward model, the
Laplace-on-a-factor-graph posterior and the evidence mixture over the penalty
weight are the ucfret prototype `s88_laplace_posterior` (investigation
`pinn_pR_anisotropy`); the figures in `figures/` that are dated before
2026-09-08 are that investigation's own results, reproduced here with their
captions. The notebooks, the photon-stream simulation and the burst analysis
are the material around it.
