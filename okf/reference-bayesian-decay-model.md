---
okf_version: "0.2"
status: reference
---

# Reference: what a Bayesian decay analysis looks like, for the graph infrastructure

**Audience**: whoever designs the `FactorGraph` and `Node`/`Port` work of
[PRD-139](prds/prd-139.md). This file is the target, not a proposal: it
describes a model that exists, runs, and is measured, so that the
infrastructure can be designed against something concrete rather than an
imagined use.

**The model lives in `../ucfret`** and stays there: the physics, the priors
and the inference are that project's business. What is asked of bff is the
graph infrastructure underneath -- a structural factor graph that can carry
this model faithfully, and an evaluation graph the user assembles and then
runs. Nothing about spectroscopy needs to move into bff.

Sources: `../ucfret/investigation/pinn_pR_anisotropy/s88_laplace_posterior.py`
(the prototype, ~5900 lines, works and is measured);
`../ucfret/okf/writeup-laplace-posterior.md` (the method and its outcomes);
`../ucfret/investigation/pinn_pR_anisotropy/S88_DEVELOPMENT.md` (how it is put
together); `../ucfret/okf/prd-laplace-reference-implementation.md` (the
consumer plan, step R3 is the one that calls into bff).

## 1. The analysis in one paragraph

A sample carries a donor dye and, on some molecules, an acceptor. A pulsed
laser excites the donor; every emitted photon is timed against the pulse and
sorted by detector (colour and polarisation). The donor's decay is shortened
where an acceptor is near, by an amount that depends on the distance to the
sixth power, so the shape of the decay carries the distribution of
donor-acceptor distances. **The object of the analysis is the posterior of
that distribution, p(R_DA | data), with every calibration constant treated as
an unknown carrying a prior rather than a fixed number.** It is computed
analytically: a Gaussian (Laplace) approximation at the mode for each value of
a smoothing weight, and the posterior as the evidence-weighted mixture over
those weights. There is no sampler, by requirement.

## 2. The data

Between 8 and 12 **histograms** (the prototype's two configurations). Each is
one sample measured on one detector: 130 to 226 bins of photon counts, widening
after each excitation pulse, 1e5 to 1e7 photons in total. Under pulsed
interleaved excitation a single histogram contains TWO excitation pulses and
their wrap-around tails, so one histogram is fed by two model channels whose
means are added.

Sizes to design against: 8 histograms x 130 bins = 1040 counts (the ensemble
configuration), 12 x 226 = 2712 (the interleaved one). Data are `double`
arrays; the likelihood is Poisson.

## 3. The parameters -- what a variable is here

The ensemble configuration has **46 variables holding 113 unconstrained
numbers**. This is the payload for "a variable needs a size and a role".

| role | count | numbers | the variables |
|---|---|---|---|
| `distribution` | 1 | **24** | `c`: 25 spline coefficients constrained to sum to zero, so 24 free numbers. The size bff needs is the FREE one, 24; the basis (Helmert contrasts) is the consumer's transform and bff need not know it |
| `physics` | 7 | 51 | `spec_eps` (33, the donor lifetime spectrum), `w_rho` (8) and `w_rho_a` (4) and `w_a` (3) (rotational and acceptor-lifetime weights, on simplices), `x_d0`, `r0_d`, `r0_a` |
| `calibration` | 12 | 12 | `g`, `l1`, `l2`, four crosstalk constants, two detection efficiencies, two quantum yields, one excitation crosstalk |
| `instrument` | 25 | 25 | per detector: response shift, width, skew (6); per sample: a log scale (3); **per histogram: scatter and background (16)** |
| `hyper` | 1 | 1 | `log10_lam`, the weight of the roughness prior |

Two facts the current interface cannot express, and both change an answer:

- **`c` is 24 numbers and `bkg_D0_g_vv` is one.** A sampling block containing
  the distribution is not the same size as one containing a background, and
  `block_cost` cannot see the difference today.
- **The roles above are not decoration.** They decide what may be held fixed
  in a diagnostic, what is reported, what a block should contain, and what a
  user interface shows. The consumer keeps a parallel dictionary for them
  today.

The interleaved configuration adds two calibration constants and one more
response set, and grows the per-histogram nuisances to 24: 121 numbers over
12 histograms.

**Machine-readable**: `../ucfret/investigation/pinn_pR_anisotropy/s88_scratch/graph_spec_ensemble.json`
and `graph_spec_pie_full.json` -- every variable with its key, size, role,
transform and prior family, and every factor with its kind and measured scope.
Regenerate with `emit_graph_spec.py`. The interleaved model is 86 variables
holding 153 numbers over 12 histograms.

## 4. The factors -- what reads what

**Likelihood, one per histogram** (8 or 12 of them). Its scope is NOT the
same for every histogram, and the difference is the physics rather than an
implementation detail. **The scopes below were measured, not declared**: every
variable was perturbed and every histogram's expected counts recomputed
(`../ucfret/investigation/pinn_pR_anisotropy/s88_scratch/emit_graph_spec.py`,
whose JSON output is the machine-readable version of this section). An earlier
draft of this file said "every likelihood reads everything global"; that is
wrong, and the measurement is why this paragraph exists.

| histogram | variables in its scope | numbers | reads the distance distribution? |
|---|---|---|---|
| donor-only, green parallel | 13 | 52 | **no** |
| donor-only, green perpendicular | 14 | 53 | no |
| labelled, green parallel | 21 | 88 | **yes** |
| labelled, red perpendicular | 21 | 88 | yes |
| acceptor-only, red parallel | 12 | 21 | no |

In the ensemble configuration all 8 histograms have DIFFERENT scopes (8
distinct patterns), and in the interleaved one all 12 do. **No variable is in
every scope.** The structure is:

- the **distance distribution `c` (24 numbers) is read by the four labelled-sample
  histograms only** -- the donor-only and acceptor-only samples have nothing to
  say about it directly, which is exactly why they are measured: they pin the
  nuisances that would otherwise be confounded with it;
- the donor's spectrum, its rotational weights and `r0_d` are read by the
  donor-only and labelled histograms, not by the acceptor-only ones;
- the acceptor's variables (`w_a`, `w_rho_a`, `r0_a`, `QY_A`, `EX_AG`, the red
  crosstalks) are read by the labelled and acceptor-only histograms;
- each histogram reads its own detector's three response numbers, its sample's
  scale, and its own scatter and background -- and, under interleaved
  excitation, the scatter and background of BOTH its pulses (its partner
  channel's nodes), which is why an interleaved histogram has four local
  nuisances rather than two;
- `g` appears only in perpendicular channels, `l1` only in parallel ones,
  `l2` only in perpendicular ones.

So this is not a complete graph and not a simple star: it is three overlapping
groups sharing subsets of a global block, and that is precisely the
factorisation worth having bff analyse.

**What is NOT a factor**: the physics between the parameters and the expected
counts (transfer maps, spectra, anisotropy, crosstalk, convolution with the
instrument response). That is a deterministic computation, and it belongs to
the evaluation graph of §5, not to the factor graph of §4.

## 5. The evaluation graph -- what actually computes

The chain from inputs to posterior densities, which is what the user should be
able to assemble from nodes and then run. Each row is a natural node; the
sizes are for the ensemble configuration.

| node | inputs | output | cost |
|---|---|---|---|
| `transforms` | the unconstrained vector (113) | the constrained values by name | negligible |
| `spectrum` | `spec_eps` | the donor lifetime spectrum (33) | negligible |
| `distribution` | `c` | p(R_DA) on a 128-point grid | negligible |
| `physics amplitudes` | the above, the physics and calibration variables, the transfer maps (fixed, ~24 MB) | abstract amplitudes per channel (8 x 35) | ~6 ms |
| `instrument basis` | the response variables | one (130 x 35) matrix per detector (2 to 4 of them) | ~3 ms each |
| `expected counts` | amplitudes, basis, scales, scatter, background | the mean of every histogram (8 x 130) | ~2 ms |
| `log likelihood` | expected counts, the data | one number | negligible |
| `log prior` | the transformed values, the hyper factor | one number | negligible |
| `jacobian` | as `expected counts` | d(mean)/d(theta), (1040 x 113) | 15 ms analytic, 100 ms by AD |
| `fisher + prior curvature` | the Jacobian, the means, the prior's Hessian | (113 x 113) | ~10 ms |
| `mode` | all of the above, iterated | theta*, 10 to 30 iterations | 5 to 20 s |
| `laplace` | the curvature at the mode | covariance, log determinant, evidence | ~30 ms |
| `mixture` | the per-weight Laplaces | weights, the posterior mixture | negligible |
| `densities` | the mixture | **p(R_DA) with its band, the summaries with their moments, the evidence per node** | ~50 ms |

Two structural points for the design. The **basis is per detector, not per
channel**: several channels share one, and rebuilding it per channel was worth
a factor of four. The **mode search is the only expensive node**, and it is a
loop over the cheap ones, so a graph that recomputes only what a changed input
touches is worth having: changing the data changes everything, changing one
prior width changes only the prior and the mode, changing the smoothing weight
changes only its own factor and the mode.

## 6. What the infrastructure has to do

1. **Assemble without evaluating.** Building the graph must compute nothing.
   Ports exist, links exist, values do not. (`Port::set_is_reactive(false)` is
   presumably how; the requirement is the owner's, 2026-09-08: "should have a
   'run' option, so that it is not eval automatically".)
2. **`run()`** walks the graph once and produces the posterior densities.
3. **Recompute only what changed.** After a `run()`, writing one input port
   and running again must re-evaluate only the nodes downstream of it. The
   consumer will assert this with a counter on the node callbacks.
4. **Save and load.** A graph, with its structure, its port values and the
   paths of the data it reads, written to JSON and read back to give the same
   posterior densities. The data are referenced, not embedded.
5. **Carry the structure faithfully** (§3, §4): variable sizes and roles, a
   hyperparameter factor kind, so that bff's components, elimination order,
   treewidth, separators and sampling blocks are about this model and not
   about a caricature of it.

## 7. Numbers, so the design is grounded

- 113 parameters, 46 variables, 8 histograms, 1040 counts (ensemble);
  121 parameters, 12 histograms, 2712 counts (interleaved).
- One fit: 5 penalty nodes x (10 to 30 Fisher iterations) = 25 s ensemble,
  20 minutes interleaved (the latter only because its Jacobian falls back to
  automatic differentiation; the fix is a consumer-side step).
- Memory: under 1 GB per fit after the Hessian was chunked; the transfer maps
  are 24 MB and are shared and read-only.
- A saved fit (the consumer's own format, for comparison) is 617 kB and loads
  in a millisecond: per weight node the mode, the covariance, the evidence,
  the expected counts and the fit statistics, plus the data and the truth.

## 8. What would make this a bad fit for bff, so it can be said early

If carrying this model would require bff to know about lifetimes, anisotropy
or crosstalk, the boundary is in the wrong place: those belong to the consumer
and enter the graph as node callbacks it supplies. What bff owns is the graph
-- structure, laziness, invalidation, serialisation, and the analysis of the
factorisation. If the callbacks must be C++ to satisfy the compute/display
line, that is a separate conversation with the consumer, whose modelling layer
is Python and whose arrays are torch tensors.

*Written 2026-09-08 from the ucfret side, at the owner's request ("the bff
impl agent needs a reference file to have context to the target"), by Claude
Code in the `bayesian.decays` session.*
