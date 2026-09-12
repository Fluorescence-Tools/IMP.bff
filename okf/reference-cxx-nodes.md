---
type: reference
title: "The C++ nodes bff already has, and which of a decay analysis's chain they cover"
description: Eleven Node subclasses exist in C++ -- lifetime spectra, anisotropy, FRET spectra, distance distributions, the TCSPC instrument, chi-squared with a Poisson option, FCS curves and the expression engine. Written so a consumer building an evaluation graph reaches for them instead of writing a Python callback per node. Maps them onto the Bayesian decay chain of reference-bayesian-decay-model.md section 5, row by row: what is already there, what is missing and generic enough that bff should grow it, and what is genuinely the consumer's.
resource: /Users/tpeulen/dev/imp.bff
tags: [reference, imp.bff, nodes, evaluation-graph, prd-139]
timestamp: '2026-09-08T00:00:00Z'
---
# What is already C++, so it need not be written again in Python

`Node` can be subclassed in Python -- chinet allowed it and the SWIG director
still does, and `EvaluationGraph` works with such nodes. But **the callback is
the wrong default**, and this repository has the numbers rather than an
opinion about it:

- a fit optimised entirely in C++ that then called the Python model five more
  times to rebuild a Jacobian for its error bars spent **32% of a TCSPC fit**
  doing it (`okf/log.md` 2026-09-01);
- a per-element SWIG director measured as a regression against plain scipy,
  1.54 ms against 1.33 (same log).

`SpectrumNode.h` states the rule this file catalogues: putting a model's
physics in the caller "would put the deriving back in Python once per
iteration, which is the arrangement that measured as a regression. A node per
producer is the third option and the only one that composes."

## The nodes that exist

| class | header | produces |
|---|---|---|
| `PhotophysicsLifetimeSpectrumNode` | `SpectrumNode.h` | an interleaved (amplitude, lifetime) spectrum from its parts |
| `PhotophysicsAnisotropySpectrumNode` | `SpectrumNode.h` | the parallel/perpendicular spectra from a spectrum plus rotational parameters |
| `FRETSpectrumNode` | `SpectrumNode.h` | a spectrum from a donor spectrum and a distance distribution |
| `PolymerDistances` | `SpectrumNode.h` | a distance distribution from a polymer model |
| `GaussianDistances` | `SpectrumNode.h` | a distance distribution from Gaussian components |
| `TCSPCDecay` | `TCSPCDecay.h` | the model curve through the instrument |
| `ChiSquared` | `ChiSquared.h` | the objective for one curve |
| `JointChiSquared` | `JointChiSquared.h` | the objective over several |
| `FCSMdfCurve`, `FCSSaturationCurve` | `FCS.h` | correlation curves |
| `Expression` | `Expression.h` | an arbitrary equation over ports |

`TCSPCDecay` is worth reading closely, because it is the largest single piece
of a decay analysis that is already done. Its ports are the amplitudes and
lifetimes (or a whole spectrum through `lifetime_spectrum`), `scatter`,
`background`, `n0` -- which it *writes* when autoscaling -- and `timeshift`.
So reconvolution, the scatter fraction, the constant offset, the scale to the
data and the shift against the response are all inside it. `ChiSquared` takes
`"default"` or `"poisson"`, the latter being signed Poisson deviance
residuals, the maximum-likelihood estimator for low counts.

## Mapped onto the Bayesian decay chain

Rows are `reference-bayesian-decay-model.md` §5, with its measured costs.

| chain node | cost | status |
|---|---|---|
| `transforms` | negligible | **missing, generic.** Unconstrained → constrained (logit, log, softmax, sum-to-zero). No spectroscopy in it; bff's, if it is wanted. |
| `spectrum` | negligible | **`PhotophysicsLifetimeSpectrumNode`.** |
| `distribution` | negligible | **partial.** `PolymerDistances` and `GaussianDistances` produce distance distributions from *their* parameterisations; a spline-basis one does not exist. Generic enough to belong here. |
| `physics amplitudes` | ~6 ms | **`FRETSpectrumNode` + `PhotophysicsAnisotropySpectrumNode`** cover the FRET and anisotropy parts. Crosstalk and the transfer maps are the consumer's. |
| `instrument basis` | ~3 ms each | **`TCSPCDecay`.** |
| `expected counts` | ~2 ms | **`TCSPCDecay`** (scatter, background, scale are its ports). |
| `log likelihood` | negligible | **`ChiSquared` / `JointChiSquared`**, with `"poisson"`. |
| `log prior` | negligible | **missing, generic.** Gaussian, uniform, half-normal on named variables. |
| `jacobian` | 15 ms analytic, 100 ms by AD | **missing, and the one that matters most.** This is exactly the node that measured as 32% of a TCSPC fit when it crossed into Python. |
| `fisher + prior curvature` | ~10 ms | **missing, generic.** `Jᵀ W J` plus the prior's Hessian: linear algebra. |
| `mode` | 5–20 s | **`Minimizer`** exists; whether it fits this loop is unchecked. |
| `laplace` | ~30 ms | **missing, generic.** Covariance, log determinant, evidence from a curvature matrix. |
| `mixture`, `densities` | ~50 ms | **missing, mostly generic.** Combining per-node Laplaces and summarising into a curve with a band. |

Read the table as a recommendation, not an inventory: **the front two thirds
of that chain already exists in C++.** A consumer writing Python callbacks for
the spectrum, the instrument and the objective would be reimplementing tested
code and paying a director crossing per iteration for the privilege.

What is missing divides cleanly. `transforms`, `log prior`, `fisher`,
`laplace`, `mixture` and a spline-basis distribution are **generic** -- no
lifetime, no anisotropy, no crosstalk in any of them -- and belong here or in
tttrlib. The transfer maps, the crosstalk and the calibration are the
consumer's and should enter as data or as its own nodes.

## What a graph document carries, and what it cannot

Related, and asked in the same breath, so it is answered here too.

`EvaluationGraph::to_json` writes the **naming** -- which label points at
which port of which node -- plus an opaque provenance string per output. It
does not write the nodes, and saving their Python source instead was measured
and does not work: `inspect.getsource` raises "is a built-in class" for a
class defined in a notebook cell or by `exec`, and where it does succeed,
executing the text in a fresh namespace fails on the first name it does not
carry, because the source is not the closure. A document that carried code
would also execute when loaded, which is not a property to give a file that
is shared between colleagues.

**This is a second argument for C++ nodes.** A graph of `TCSPCDecay`,
`ChiSquared` and their kin is describable by name and by parameter, so it
*can* be saved and reloaded whole. A graph of Python callbacks can only be
saved as naming plus a fingerprint, with the caller rebuilding the code. The
more of a model lives in C++ nodes, the more of it a document can carry.

## Where this came from

Written 2026-09-08 while building `EvaluationGraph` for PRD-139, after the
owner observed that "in old chinet it was possible to set a python callback,
but it is nicer if the cpp already offers most that is needed, either in
tttrlib or in bff". The consumer's side is `../ucfret`; the target model is
`reference-bayesian-decay-model.md`.
