# Bayesian fluorescence decays — how the model was arrived at

*2026-09-08*

Eight notebooks that build one model, in the order it was actually built. The
end of the road is `../smfret_pie_mfd`, which is the whole thing as a single
runnable example; this series exists so that nothing in it looks arbitrary.

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
| 1 | `01_the_measurement.ipynb` | what the data are; the forward model; the ruler $E(r)=1/(1+r^6)$ and the window where it stops reading; that the expected counts are **linear in the distance distribution** |
| 2 | `02_the_inverse_problem.ipynb` | that inverting a decay for $p(R)$ is ill posed even with every nuisance known exactly — 15 of 128 directions above the noise — and what a roughness prior is for |
| 3 | `03_the_penalty_weight.ipynb` | that the prior's weight **cannot** be fitted jointly (the joint mode flattens two populations into one), and the evidence mixture that replaces it |
| 4 | `04_the_posterior.ipynb` | why drawing from the posterior and pushing it through the softmax overshoots, and the delta method that replaces it; the band's space; the gauge |
| 5 | `05_calibration_and_limits.ipynb` | whether the intervals mean what they say; the deviance's true reference; the bias and the five explanations excluded by measurement; the resolution rule |
| 6 | `06_the_experiment.ipynb` | that the bias is a property of the measurement, and what a second laser buys — computed from the information alone, then measured |
| 7 | `07_single_molecules.ipynb` | the photon stream, `tttrlib` bursts, and an A/B showing the single-molecule experiment costs the inference almost nothing |
| 8 | `08_classify_then_pool.ipynb` | classify the bursts by what each one looks like, pool the groups, analyse each — and separate populations the pooled analysis cannot |

Notebooks 1, 2 and 4–8 run in a few minutes each; 3 and 5 do several fits and
take longer. Every one was executed end to end before being committed, and
carries its outputs.

## What it needs

* `tttrlib` ≥ 0.27 — the photon stream, the burst search, and the HDBSCAN used
  in notebook 8
* the ucfret prototype `s88_laplace_posterior.py`, which carries the physics
  and the posterior. Looked for at
  `~/dev/ucfret/investigation/pinn_pR_anisotropy`; set `UCFRET_S88` to point
  elsewhere.
* `numpy`, `scipy`, `torch`, `matplotlib`, and `scikit-learn` for one
  cross-check in notebook 8

`bd.py` is the shared code. It imports `pie_mfd` from `../smfret_pie_mfd` so
that the series and the finished example run the same implementation.

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
