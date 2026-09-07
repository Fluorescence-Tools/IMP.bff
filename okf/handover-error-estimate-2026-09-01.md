# Handover: the optimiser, and everything it needs, is bff's

> ## DONE, 2026-09-01 20:10 — read `okf/log.md` (17) for what actually happened
>
> All three stages landed. A `fit.run()` now makes **2** Python
> `update_model()` calls whatever the model is (was 7 / 8 / 9), the error
> estimate is 0–2% of a fit rather than 32%, and the two other
> Levenberg-Marquardt implementations are deleted. imp.bff **1830 passed**
> (was 1816), none failed; chisurf `test/fitting` **1003 passed** with the
> same two pre-existing failures.
>
> **This page was the plan and it is kept as the plan** — the reasoning about
> *why* the refusal was correct, what the step rule had to be, and the ten
> traps are all still the best account of the problem. Two things in it turned
> out to be wrong, and both are corrected below where they appear:
>
> 1. **"The director is now a wash."** It is not: timing `minimize` alone
>    with the covariance switched off, the director is **1.11x** scipy on the
>    parse fit and **1.06x** on a TCSPC decay. The measurement that said
>    otherwise timed whole `run()` calls, by which point one side's error
>    estimate had already moved into C++ and the other's had not. See
>    *The objection...* below.
> 2. **`evaluate_external` was private, not public.** It is public now, and
>    documented as the entry point a covariance is written on top of.
>
> One defect was found on the way and is *not* fixed: a Python `Node` held by
> a `Minimizer` can be garbage-collected and segfault it —
> **`T-20260901-13`**. It predates this ticket but now sits on a path users
> take, because the fallback is a `ResidualNode`.

Written 2026-09-01. Everything below is measured on this machine (Apple
silicon, `arm64` conda env, Python 3.12) unless it says otherwise.

**Ticket: `T-20260901-10`** (board: `okf/agent-board.md`, Resolved).
Adjacent, deliberately *not* this ticket: `T-20260901-11`.
Background: `okf/log.md` entries **2026-09-01 (9)**, **(12)** and **(15)**;
the rule this serves is
[`../chisurf/okf/architecture/compute-display-line.md`](../../chisurf/okf/architecture/compute-display-line.md).

> **Scope widened by the owner, 2026-09-01 17:37: *"lmdif etc must be also in
> bff."*** This page was first written as "compute the covariance in C++" and
> proposed keeping ChiSurf's numpy optimiser as a documented fallback. That
> is superseded. There is to be **one** implementation of the
> Levenberg-Marquardt machinery and it is bff's; ChiSurf calls it.

---

## The one-paragraph version

`lmdif` and its helpers are *already* in bff — `src/Minimizer.cpp`, 1041
lines, MINPACK transcribed line by line. What is not in bff is everything
around it: ChiSurf still carries a **second** bounded Levenberg-Marquardt
(794 lines of `leastsqbound.py`, still the fallback for every model that
cannot build a graph), a **third** copy inside a plugin, a second
finite-difference Jacobian with a *different step rule*, and a second
covariance with six consumers. The immediate symptom is measurable — a TCSPC
fit optimises entirely in C++ and then calls the Python model **five more
times** to rebuild a Jacobian for the error bars, 32% of the fit — but the
reason to do it is that *two implementations of one algorithm do not average
out, they disagree, and the disagreement is silent.* This stack has already
paid for that once, with a treewidth of 3 where the answer is 1.

---

## One algorithm, three implementations

| where | lines | what it is | still used? |
|---|---|---|---|
| `imp.bff/src/Minimizer.cpp` | 1041 | MINPACK `lmdif` + `enorm`, `qrfac`, `qrsolv`, `lmpar`, `covar`, `fdjac2`, plus ChiSurf's bounds transform | yes — every model that builds a graph |
| `chisurf/core/math/optimization/leastsqbound.py` | 794 | scipy `leastsq` + the *same* bounds transform | yes — the fallback for every model that does not |
| `chisurf/plugins/fluorescence_decay/lltf/core/optimization/leastsqbound.py` | 365 | a plugin's own copy, reached from `lltf/core/fitter.py:968` | yes, independently |

And the pieces the optimiser needs, doubled the same way:

| where | what it is | note |
|---|---|---|
| `Minimizer.cpp` `fdjac2` (~line 786) | forward-difference Jacobian | `h = eps*|x_j|`, floored only at exactly 0 |
| `chisurf/core/fitting/fit.py:2748` `approx_grad` | forward-difference Jacobian | `step = eps*max(|x_k|, 1.0)` — **a different step rule** |
| `Minimizer.cpp` `covar` (~line 308) | `(J'J)^-1` from the final R | published as `get_covariance()` |
| `chisurf/core/fitting/fit.py:2875` `covariance_matrix` | `(J'J)^-1` from `approx_grad` | **six consumers besides the error estimate** — see *Scope surprise* |

The two Jacobians differ in both the multiplier and the floor:

| | multiplier | floor |
|---|---|---|
| `lmdif` / `fdjac2` | `sqrt(epsfcn)` = **1e-3** at `epsfcn = 1e-6` | only at exactly 0 |
| `approx_grad` | `sqrt(machine eps)` = **1.49e-8** | at **1.0** |

That is not a detail. It is the whole reason the C++ covariance gets refused,
and the numpy step rule is the one the error bars are currently pinned to.

---

## What the numbers are

`fit.run()` timed, and Python `update_model()` calls counted inside it
(512 channels; see *How to reproduce the measurement*):

| fit | total | optimise | error estimate | rest | `update_model()` per `run()` |
|---|---|---|---|---|---|
| parse | 0.50 ms | 60% | 4% | 36% | **2** |
| TCSPC lifetime | 3.26 ms | 41% | **32%** | 27% | **7** — 5 of them for the errors |
| TCSPC VV | 4.51 ms | 44% | **34%** | 22% | **8** — 6 for the errors |
| FRET (Gaussian) | 24.3 ms | 80% | 12% | 8% | **9** — 7 for the errors |

The parse row is what compliance looks like: two model evaluations in a whole
`run()`, none of them for the covariance. **The target is to make the TCSPC
row look like the parse row.**

---

## The objection that used to block this no longer holds — measured

`okf/log.md` 2026-09-01 (9) recorded that routing a Python residual through
a C++ `Node` director was a **regression** (1.33 ms against scipy's 1.32),
and that is why `minimize` falls back to `leastsqbound` rather than to
`ResidualNode`. Both the log and the module docstring in
`chisurf/core/fitting/minimizer.py` still say so.

**Re-measured on the current build, three independent runs of
`test/minimizer/bench_fit.py` (interleaved, min-of-many):**

```
run 1:   scipy 1.307 ms    director 1.301 ms
run 2:   scipy 1.293 ms    director 1.274 ms
run 3:   scipy 1.299 ms    director 1.289 ms
```

The director is now a **wash or marginally faster** — the gap is well inside
the noise of a measurement that resolves 2.7x cleanly.

> ### ✗ Wrong, and corrected while doing the work
>
> **That comparison is confounded.** It times whole `fit.run()` calls, and a
> `run()` includes the error estimate — which by then had already moved into
> C++ for the scipy row (it leaves no covariance behind, so
> `covariance_matrix` runs, and stage 2 put *that* on the graph) while the
> director row differences the director, in Python. The two rows were not
> doing the same work.
>
> **Timing `minimize` alone, with the covariance switched off so both do the
> same work:**
>
> ```
>                 scipy      director    ratio
> parse           0.863 ms   0.961 ms    1.11x
> tcspc          11.74  ms  12.47  ms    1.06x
> ```
>
> So the director is **6–11% slower**, not equal. It was 1.29x / 1.22x until
> `ResidualNode` stopped calling `get_input_port` once per parameter per
> residual evaluation — for an answer settled in its constructor.
>
> **This does not change the conclusion, only the argument for it.** The
> reason to have one implementation was never speed: three copies of one
> algorithm and two of its Jacobian, with *different step rules*, disagree
> quietly, and this stack has already paid for that once. A few per cent on
> the path taken by models that build no graph is a cheap price, and
> `T-20260901-08`/`-09` are how to stop paying it.

Correct `okf/log.md` (9) and the `minimizer.py` docstring while you are there
— they assert something false either way, and the next person to read them
will conclude this ticket is a bad idea. *(Done: both now carry the measured
figures.)*

---

## Why the covariance falls back, which is not a bug anywhere

Read this before touching anything: the obvious fix ("just use the C++
covariance") is wrong and the second-most-obvious ("lower `epsfcn`") is
worse.

`Minimizer` ends holding a QR factorisation of its Jacobian and publishes
`(J'J)^-1` from it (`Minimizer::get_covariance`, external coordinates, via
MINPACK's `covar`). `Fit._optimiser_covariance` prefers it, and
`chisurf/core/fitting/minimizer.py:_forward_differences_resolved` decides
whether it may be trusted. It refuses whenever the free vector spans orders
of magnitude, because:

* **`lmdif` differences at a step relative to the parameter.**
  `src/Minimizer.cpp`, in the `fdjac2` block (~line 786):

  ```cpp
  const double eps = std::sqrt(std::max(epsfcn_, MACHINE_EPS));
  double h = eps * std::fabs(temp);
  if (h == 0.0) h = eps;          // a floor only at *exactly* zero
  ```

* So a parameter that converged near zero is differenced at a step near
  zero, its Jacobian column is round-off, and `covar` reports an enormous
  variance for it. **Nothing about that is wrong in the optimiser** — it is
  what a relative step means, and *scipy's own* `leastsqbound` covariance for
  the same fit is no better (0.094 on the same parameter and exactly zero for
  two others). The two implementations fail together, which is how it is
  known to be the estimator's property rather than the port's.

* **Every real TCSPC model hits it.** The default free vector is a scatter
  fraction of `1.5e-5` beside a lifetime of `3.15` — five orders apart.

* **Lowering `epsfcn` is not the way out.** `fit.py:DEFAULT_EPSFCN = 1e-6`
  exists because `epsfcn = 0` (machine epsilon, ~1.5e-8 relative) is *below
  the noise floor* of a recursive convolution over ~1024 channels: over 88
  randomised fits across four conditions, fits reaching `chi2r < 1.1` went
  **60/88 at `epsfcn = 0` to 83/88 at `1e-6`**. That constant is load-bearing
  for convergence and must not be moved to fix the covariance.

The resolution is that **the optimiser's step rule and the covariance's step
rule are different quantities** and should stop sharing a number. One is
tuned for convergence, the other for resolution. `approx_grad`'s rule —
`eps*max(|x|, 1.0)` with `eps = sqrt(machine eps)` — is the specification for
the second, because it is what the current error bars are pinned to.

---

## The task, in three stages

Each stage is independently landable and independently testable. Do them in
this order; stage 1 is the measured 32% and stage 3 is the largest diff.

### Stage 1 — the covariance, in C++, over the graph

`Minimizer` already has the entry point — though it was **private**, not
public as written here; it was made public as part of this work, and
documented as what a covariance routine is built on top of:

```cpp
std::vector<double> evaluate_external(const std::vector<double>& xe);
```

It evaluates the objective at an external parameter vector and returns the
residuals, so this needs **no change to `lmdif` or `epsfcn` at all**.

Add something like
`std::vector<double> Minimizer::compute_covariance(double epsilon, double floor)`
which:

- reads the solution (`get_x()`) and either evaluates once for `f0` or reuses
  `get_residuals()` — check it is current for `x` before trusting it;
- steps each parameter by `epsilon * max(|x_j|, floor)`, rounded to an
  exactly representable difference the way `approx_grad` does
  (`step = (x + step) - x`), and forms the forward difference through
  `evaluate_external`;
- drops columns whose partial derivative is identically zero, exactly as
  `covariance_matrix` does with `important_parameters` — see trap 3;
- forms `alpha = J'J` and returns its pseudo-inverse. The Python side uses
  `scipy.linalg.pinvh`; Eigen's
  `completeOrthogonalDecomposition().pseudoInverse()` is the counterpart, and
  the choice is observable on rank-deficient columns, so pin it against the
  numpy answer rather than assuming;
- publishes in **external** coordinates — `Minimizer.h` documents that
  `get_covariance()` runs `covar` on the external Jacobian
  (`_internal2external_grad`); do the same or every bounded fit disagrees.

Then `Fit.update_error_estimates` asks the minimizer instead of calling
`covariance_matrix(fit, f0=...)`. Keep `_forward_differences_resolved` as the
guard on the **QR** covariance — that matrix is free when it is trustworthy —
and make this the fallback beneath it rather than numpy.

### Stage 2 — the other six consumers of `covariance_matrix`

This is the scope surprise and it is why stage 1 does not finish the job.
Beyond the error estimate, `covariance_matrix` is called by:

| caller | what for |
|---|---|
| `chisurf/core/fitting/graphview.py:328` | the posterior view |
| `chisurf/core/fitting/derived.py:280` | error propagation onto derived quantities |
| `chisurf/core/fitting/engine.py:551, 660` | sampling preconditioner |
| `chisurf/core/fitting/sample.py:313` | MCMC preconditioner |
| `chisurf/core/fitting/fit.py:610` | the `gradient` property |
| `chisurf/core/fitting/fit.py:623` | the `covariance_matrix` property |

Each of these builds a `p+1`-evaluation numpy Jacobian over the Python model.
Several run *per sampler move*. Route them through the same C++ routine —
but note the signature difference: they pass `model=` explicitly, because for
a `FitGroup` `fit.model` is the *selected member's* model rather than the
global one. Do not lose that.

### Stage 3 — delete the second and third optimisers

- `minimize`'s fallback becomes `ResidualNode` + `bff.Minimizer` instead of
  `leastsqbound`. The director measurement above says this is free.
- `chisurf/core/math/optimization/leastsqbound.py` (794 lines) goes, along
  with its re-export in `core/math/optimization/__init__.py`. Watch two
  live users: `OptimizationCancelled` is imported by `fit.py:32`, and
  `sample.py:1297` calls `leastsqbound` directly.
- `chisurf/plugins/fluorescence_decay/lltf/core/optimization/leastsqbound.py`
  (365 lines) is a **third** copy with its own caller
  (`lltf/core/fitter.py:968`). It is a plugin, so it may reasonably land
  after the core — but it is the same algorithm and the same argument, so it
  should not be left unmentioned.
- `approx_grad` and `covariance_matrix` go once stage 2 has moved their
  callers. **Do not delete them before that** — they are the reference the
  C++ routine is pinned against.

### Done when

- A TCSPC `run()` makes **at most 2** Python `update_model()` calls, as a
  parse fit already does. Count them; do not infer it from the clock.
- No `leastsqbound` remains under `chisurf/core/`, and `grep -rn "scipy.optimize"`
  over the fitting package comes back empty (or each hit has a written
  reason).
- `test_the_error_estimates_do_not_move` passes in all three graph-fit files
  **after being made to actually test something** — trap 1.
- `chisurf test/fitting` still **1000 passed** plus the two known
  pre-existing failures (`test_fit_state`, `test_pcf_experiment`; both fail
  on `HEAD` too — verify with `git show HEAD:<file>`, **not** `git stash`).
- `imp.bff` still **1816 passed**, plus tests for the new method.
- `test/minimizer/bench_fit.py` grows the split measurement, and the "the
  director is a regression" claim is corrected in `okf/log.md` (9) and in the
  `minimizer.py` module docstring.

---

## Traps

**1. The test that looks like it protects you does not, for TCSPC.**
`test_the_error_estimates_do_not_move` compares the graph path against the
scipy path — but for a TCSPC fit **both currently fall back to the same numpy
`covariance_matrix`**, so it is comparing numpy against numpy and passes
trivially. It will not catch a wrong C++ Jacobian. Before changing anything,
make it compare the *new* C++ matrix against `covariance_matrix`'s answer
explicitly, and watch it fail. (It genuinely discriminates for the parse fit,
where the QR covariance is accepted — which is why nobody noticed.)

**2. `test_a_covariance_from_an_unresolved_step_is_not_reused` needs
rewriting, not deleting.** The property it defends is still true: *a
round-off column must not become an error bar.* What changes is the remedy —
today it refuses the matrix, tomorrow it differences at a step that resolves.
Keep a test that fails if a near-zero parameter comes back with a
plausible-looking error bar it has not earned.

**3. A column of exact zeros is legitimate and must survive.** `E_FRET` is a
free parameter neither optimiser can move — its value comes from a callable
and its setter is ignored — so its Jacobian column is exactly zero **by
design**, and the graph gives it a deliberately dangling port.
`covariance_matrix` handles this by dropping such parameters from
`important_parameters`, and `_propagate_redundant_error_estimates` then fills
them in. Invert a singular `alpha` instead and every FRET fit gets garbage
error bars.

**4. `_optimiser_covariance` *pops* the stash.** `self.__dict__.pop(...)`, so
it is consumed on first read and a second call silently falls back. If you
add a second consumer, that is where the bug will be.

**5. The guard is applied at stash time, not at read time.**
`minimizer.py:~1258`, inside `minimize`, before `fit._cpp_covariance` is set.
So "the C++ covariance was offered" and "the C++ covariance was used" are
different questions; the instrumentation below answers the second.

**6. `covariance_matrix` takes an explicit `model=`, and for a group it
matters.** `fit.model` on a `FitGroup` is the *selected member's* model, not
the global one. Four of the six stage-2 callers pass `model=` for exactly
this reason.

**7. Two `ninja` runs in `/Users/tpeulen/dev/imp/cmake-build-arm64` corrupt
it, and the symptom lies.** You get `ld: library 'IMP.algebra-lib' not found`
and a half-linked `_IMP_bff.so`. Recovery: `$E/bin/cmake .` then a single
`ninja -j8`. A `ninja` that timed out into the background is **still
running**; wait for it. A new *header* also needs `$E/bin/cmake .` first, or
SWIG reports `Unable to find 'IMP/bff/<X>.h'`; a new *method* on an existing
class does not.

**8. Do not `git stash` in these checkouts.** Both carry pre-existing stashes
and long-lived working-tree modifications that are not yours — including
whole directories that are untracked (`chisurf/core/fitting/` is `??`). To
A/B against `HEAD`, copy the file aside and `git show HEAD:<path> > <path>`,
then copy back.

**9. `epsfcn` is not yours to tune.** See above: `1e-6` is worth 23 more
converged fits out of 88. If the covariance needs a different step, it needs
its *own* step.

**10. SWIG.** No `#` comments inside `%pythoncode` — SWIG reads them as its
own preprocessor directives. A director held *across* calls needs
`IMP_SWIG_DIRECTOR`, not a bare `%feature("director")`; the latter keeps only
a weak pointer to the Python proxy and segfaults when it is collected. This
matters for stage 3, where `ResidualNode` becomes the fallback and will be
held for the length of a fit.

---

## How to reproduce the measurement

Both were run from `/Users/tpeulen/dev/chisurf` with `PYTHONPATH=.`,
importing `bench_fit` from `/Users/tpeulen/dev/imp.bff/test/minimizer` for
its four fixtures (`make_fit`, `make_decay_fit`, `make_polarised_fit`,
`make_fret_fit`).

**The time split** — wrap `M.minimize` and `Fit.update_error_estimates` in
timers, run `fit.run()`, report best-of-6 for the total. Do not report a
single run; load on this laptop drifts by more than the effect
(`okf/log.md` 2026-09-01 (5) learned that the expensive way).

**The evaluation count, which is the honest metric** — monkeypatch
`type(fit.model).update_model` with a counting wrapper, and record
`"_cpp_covariance" in self.__dict__` at the top of `update_error_estimates`
to separate *offered* from *used*:

```
parse      C++ covariance offered=True | update_model() during errors: 0 | total 2
tcspc      C++ covariance offered=True | update_model() during errors: 5 | total 7
tcspc VV   C++ covariance offered=True | update_model() during errors: 6 | total 8
FRET       C++ covariance offered=True | update_model() during errors: 7 | total 9
```

`offered=True` everywhere — the stash is written and then refused. That is
the whole of stage 1 in one line.

Both scripts were scratch and are not checked in; they are short enough to
rewrite from this description, and if you want them permanent,
`test/minimizer/bench_fit.py` is where they belong.

---

## Also open, and deliberately not this ticket

**`T-20260901-11` — the model curve is recomputed for display.** After a fit
the answer is written back through the ordinary setters and
`model.update_model()` re-runs ChiSurf's own pipeline, for a curve the graph
already holds on an output port. The obstacle is real and documented: the
graph is **private to the fit** because `Parameter.value` keeps a
`_frozen_value` cache inside `factorgraph.frozen_structure()`, so a port
written from C++ is not seen by the next Python read. But a *curve* is not a
`Parameter`, so publishing it may be much cheaper than publishing parameters
would be. After this ticket, not with it — this one is bounded, that one is
exploratory.

**The data are duplicated.** `DataCurve` holds `x`/`y`/`ey` as numpy and
`ChiSquared.set_data_arrays` copies them into the node. One copy per fit is
cheap; "the data stay local" is still not literally true.

**The models that never build a graph at all.** `T-20260901-08` (the FRET
distance distributions still without a node) and `T-20260901-09` (DEER, ICS,
PCH/FIDA, PDA2C/3C, MFD). The scoreboard is
`test/minimizer/census_models.py`, currently reading *"7 of the 42
constructible model classes build a graph; 3 of the 14 with a polarisation
build one under VV; no model builds a graph that disagrees with its own
curve."* Keep that last clause true. Note that stage 3 makes these models
*slower to notice*: once the fallback is also bff, a refused model no longer
announces itself by being on a visibly different code path, so the census is
the only thing that will say so.

---

## How to build and test

```bash
E=/Users/tpeulen/mambaforge/envs/arm64
B=/Users/tpeulen/dev/imp/cmake-build-arm64

ninja -C $B -j8                      # ONE at a time; ~1 min for bff-only
$E/bin/python -m pytest test/minimizer test/chi2 test/expression \
    test/spectrum test/sampler test/portnode test/factorgraph test/session -q
# 304 passed, 43 subtests passed

cd /Users/tpeulen/dev/chisurf
PYTHONPATH=. $E/bin/python -m pytest test/fitting -q
# 1000 passed, 2 known failures, ~3 min

$E/bin/python /Users/tpeulen/dev/imp.bff/test/minimizer/bench_fit.py
$E/bin/python /Users/tpeulen/dev/imp.bff/test/minimizer/census_models.py
```

ChiSurf is **not installed** in that env — it is imported from the sibling
checkout, so either `cd` there or set `PYTHONPATH`. The imp.bff tests that
compare against ChiSurf add the sibling path themselves.

---

## Why this is worth doing, in two sentences

The fit was made careful never to cross the boundary while it runs, and then
crosses it five times the moment it stops — so a third of the saving is
handed back, for a matrix that could be built where the model and the data
already are. And underneath that: there are **three** implementations of one
bounded Levenberg-Marquardt in this stack and **two** of its Jacobian, the
two Jacobians already use different step rules, and the last time two
implementations of one algorithm were finally run against each other here,
one of them turned out to have been wrong for months.
