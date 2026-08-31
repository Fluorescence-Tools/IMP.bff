# Greedy pair selection against Olga

**Date: 2026-08-31.** The `select_informative_pairs` family had been pinned to
the Python it replaced and, for the chi-squared tail, to
`scipy.special.gammaincc`. Neither is Olga. This is the comparison against
**Olga's own selector**, compiled from `../chisurf/junk/olga/src/best_dist.h`
and run on the same two matrices.

## Result

| case | selection | expected-RMSD decay |
|---|---|---|
| **T4L docking ensemble**, 50 frames x 33 pairs, err 0.06, 10 selected | **identical, all ten, in order** | 1.6e-4 A max over steps 2-10; **3.2e-3 A at step 1** |
| **synthetic**, 60 frames x 12 well-separated candidates, err 0.05, 8 selected | **identical, all eight, in order** | < 5e-5 A over steps 2-8; 2.3e-3 A at step 1 |
| **synthetic**, 60 frames x 25 candidates cut from one cloth | same eight pairs, **two adjacent ranks swapped** | see "the tie" below |

`test/restraints/test_greedy_olga_ab.py` holds all three, and skips itself
wherever the Olga sources are not checked out.

## Why they do not agree to machine precision

Olga carries `Eigen::MatrixXf` — **float32** — and does not evaluate the
chi-squared right tail. It fits it with a **64-piece quadratic spline**
(`chiSqRTSpline`), rebuilt at every greedy step for that step's `ndof`. This
module carries float64 and evaluates the tail.

The spline's accuracy, measured over the domain it is fitted on:

| ndof | worst absolute error in probability | at chisq |
|---|---|---|
| **1** | **1.125e-1** | 0.02 |
| 2 | 1.489e-2 | 0.03 |
| 3 | 1.370e-3 | 0.03 |
| 5 | 8.280e-5 | 0.47 |
| 9 | 3.229e-5 | 1.79 |

At one degree of freedom the fit is off by **0.11 in probability** near the
origin — the tail has an infinite derivative there and a quadratic cannot
follow it. Olga's own source knows: `chiSqRTSpline` ends with

```cpp
// assert(fabs(diff)<tolerance || ndof==1);
```

an accuracy assert that is commented out and carries an `|| ndof==1` escape.

This is not a corner case in the greedy. `bestPair` uses
`nDof = max(int(selPairs.size()) - 1, 1)`, so **ndof is 1 for the first two
selections** — the two that matter most, and the only steps where the two
implementations' decays differ by more than two ten-thousandths of an
angstrom.
Step 1 of the T4L run: Olga 3.087162 A, this module 3.083970 A.

## The tie

The third case was built to fail and does: 25 candidates all drawn from the
same family over the same three-cluster ensemble, so several are
interchangeable. At the second selection the best two candidates are
**2.6e-4 A apart out of 1.70 A** — a relative gap of 1.5e-4, well inside what a
float32 spline carrying 1e-1 error at ndof=1 can reorder. This module takes
pair 11, Olga takes pair 2, and one step later they swap back: the same eight
pairs come out, in a different order.

The test asserts the only thing that can be asserted rigorously — that a
disagreement happens **only on a tie**. Up to the first divergence both
selectors have added the same pairs, so the accumulated chi-squared is
identical and the two picks are directly comparable; after it the histories
differ and a score under one says nothing about the other. At that first
divergence the gap must be under a thousandth of an angstrom, and this module's
pick must be the minimum of its own scoring.

Where the candidates have genuinely different separating power — which is the
case a real experiment plan is, and what both the T4L ensemble and the
well-separated synthetic are — the selections are identical.

## Reproducing

Olga is a Qt/pteros application; the selector is not. Two edits make
`best_dist.h` compile standalone, neither semantic, both applied to a copy:

* drop the `pteros/pteros.h`, `theobald_rmsd.h` and `center.h` includes and
  everything after `sys2xyz` that needs them. What is kept is the whole
  selection: `chiSqRTSpline`, `rmsdColMean`, `rmsdMeanMeanAdd`, `rmsdMeanMean`,
  `chiSquared`, `bestPair`, `greedySelection`, `precisionDecay`.
* `spline.hpp:165` needs the `template` disambiguator clang requires on a
  dependent `.cast<int>()`. GCC accepted it as written.

The test does both in a temporary directory and never writes to the checkout.

## What this does not cover

The efficiencies and the RMSD matrix themselves. Both implementations were
handed the *same* matrices, computed here by
`ProbeNetworkRestraint::get_pair_efficiencies` and `pairwise_rmsd`. Whether
this module's accessible volumes agree with Olga's is a separate question, and
`okf/validation/av_vs_rotamer.md` is the nearest thing to an answer.
