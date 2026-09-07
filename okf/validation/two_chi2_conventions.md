---
type: validation
title: "Three implementations of one asymmetric chi2 — the restraint scorer had the error bars swapped"
description: Three implementations of the asymmetric chi2 disagreed - the C++ restraint scorer picked the opposite error bar from the C++ kernel and the Python diagnostics, up to 100x apart, with no test asserting any of them. Resolved 2026-08-24 on the model-minus-data residual (a too-large model is a positive deviation): one implementation, four pinned restraint values moved, and the branch is now pinned in both directions.
resource: /Users/tpeulen/dev/imp.bff
tags: [validation, imp.bff, chi2, restraint, docking, fps, asymmetric-errors, prd-117]
timestamp: '2026-08-24T00:00:00Z'
---

# Three implementations of one asymmetric chi2

An experimental FRET distance carries **asymmetric** errors: `error_neg` and
`error_pos` (FLR-CIF `distance_error_minus` / `distance_error_plus`). Which of
the two applies depends on which side of the experimental value the model
lands. This repository decided that in three places -- `chi2_score` in
`AVDistance.cpp`, `AVPairDistanceMeasurement::score_model` in `AV.cpp`, and
`PairDistance.chi2` in `docking.i` -- and the middle one, which is what every
optimisation minimises, disagreed with the other two.

**C++** — `AVPairDistanceMeasurement::score_model` (`src/AV.cpp:36`), which is
what `AVNetworkRestraint`, `AVMeanDistanceRestraint` and the docking engine
(`docking.i`, through `meas.score_model`) actually minimise:

```cpp
auto ev = [](double f, double m, double en, double ep){
    double dev = m - f;
    double w = (dev < 0) ? 1. / en : 1. / ep;
    return .5 * algebra::get_squared(dev * w);
};
return 0.5 * ev(model, distance, error_neg, error_pos);
```

Called as `ev(model, distance, ...)` the parameters bind `f = model`,
`m = distance`, so `dev = experiment - model` -- the deviation is measured
one way and the branch below it reads as if it were measured the other. A
model that is **too large** gives `dev < 0` and is scored against
**`error_neg`**.

**Python** — `PairDistance.chi2` (`pyext/IMP_bff.docking.i`), what the
FPS-style diagnostic tables report to the user:

```python
d = self.distance_model - self.distance_exp
err = self.error_pos if d > 0 else self.error_neg
return (d / err) ** 2
```

A model **above** the experiment is scored against **`error_pos`**.

## Measured

Experiment 50 Å, `error_neg` 1 Å, `error_pos` 10 Å:

| model | C++ `score_model` | Python `chi2` | ratio |
|---|---|---|---|
| 60 Å (10 Å above) | 25.0 | 1.0 | 25x |
| 40 Å (10 Å below) | 0.25 | 100.0 | 1/400 |

The restraint scorer punishes the direction the kernel and the diagnostic
forgive, and forgives the direction they punish. There is also a second, smaller
difference: `ev` already multiplies by `0.5` and the caller multiplies by
`0.5` again, so the C++ value is `0.25 * (dev/err)^2` — neither a chi2
(`1.0 *`) nor a log-likelihood (`0.5 *`).

## What is not covered

`test/test_AVNetworkRestraint.py:105` calls `score_model` and **discards the
result**; line 74-76 asserts only that `unprotected_evaluate` returns the same
number twice. No test in the tree asserts what either implementation computes,
which is why a 100x disagreement sat in it.

## The convention, and what it changed

**The residual is model minus data** (maintainer's ruling, 2026-08-24): a model
distance that is too large is a positive deviation, and a positive deviation is
judged against the error on that side, `error_pos`. A model that is too small
is negative and judged against `error_neg`.

Under that rule the three implementations resolve like this:

- `chi2_score` **was already right** -- residual `model - data`, branch on its
  sign -- and is untouched but for its comment. It is the one implementation
  now.
- `AVPairDistanceMeasurement::score_model` **changed**: it computed
  `data - model` and then took the branch as if that were `model - data`, so
  the error bars came out swapped. It delegates to `chi2_score` now (times the
  0.25 it has always carried). This is the defect above, fixed.
- `PairDistance.chi2` / `.efficiency_*` in `docking.i` **stopped being a
  second copy**: they call `chi2_score` and `fret_efficiency`. Its `residual`
  was already `model - data` and stays that way.

**The 0.25 is gone too** (same day, once it was clear the package is
pre-release): `ev` multiplied by 0.5 and its caller multiplied by 0.5 again, so
the restraint was worth a quarter of `chi2_score` where a Gaussian is worth a
half -- \(-\log L = z^2/2\), which is exactly what `IMP::core::Harmonic`
scores for \(k = 1/\sigma^2\). At 0.25 this restraint entered a scoring
function at half the weight of every IMP harmonic beside it, which is not
something a scale factor may do. `score_model` is `0.5 * chi2_score` now, and
the pins below doubled with it.

### Four pinned values moved

The AV network restraint's own regression pins -- not parity against another
program -- carried the swapped branch:

| test | before | after the branch fix | after the 0.5 |
|---|---|---|---|
| `test_AVNetworkRestraint.py` (quadrature) | 13.505098587465483 | 11.349675582043005 | 22.69935116408601 |
| `test_av_lattice.py` legacy MC | 11.918 | 10.0 | 20.0 (places=0) |
| `test_av_lattice.py` stencil 30 | 13.079781979252157 | 11.066107280330783 | 22.132214560661566 |
| `test_av_lattice.py` stencil 26 | 13.508167976104748 | 11.352757271766384 | 22.705514543532768 |

The model *distances* in those tests are unchanged and still assert -- what
moved is only the score they are weighed with, which is what a branch fix
should look like.

### And it is pinned now

`test_distance_conventions.py::TestChi2Score` asserts which error bar divides
in each direction, that `score_model` is `0.25 * chi2_score` across a sweep of
models, and that the asymmetry runs the way the errors say. The reason this
divergence survived is that the one test calling `score_model` threw the
result away.
