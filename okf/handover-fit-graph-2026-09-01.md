# Handover: the fit graph — `GlobalFitModel` onto `JointChiSquared`

Written 2026-09-01. Everything below is measured on this machine (Apple
silicon, `arm64` conda env, Python 3.12) unless it says otherwise.

> **Done 2026-09-01.** The task described below -- the ChiSurf-side builder
> that maps a `GlobalFitModel` onto `JointChiSquared` -- landed; see
> `okf/log.md` **2026-09-01 (11)** for what it turned out to be and what the
> group measures now (**1.74x**, from 1.01x). This page is kept because the
> **Traps** section still applies to anything touching this code, and because
> the follow-on `T-20260901-03` starts from the same ground. Read the traps;
> treat "The task" as history.

**Ticket: `T-20260901-02`** (board: `okf/agent-board.md`, Resolved).
Follow-on: `T-20260901-03`.
Background: `okf/log.md` entries **2026-09-01 (9)**, **(10)** and **(11)**.

---

## The one-paragraph version

A plain ChiSurf `Fit` on a parse model now runs its whole optimisation in
C++ — parameters, model curve, residuals and Levenberg-Marquardt — and is
**2.9x** faster. A `FitGroup` is **not**, and a `FitGroup` is what the GUI
builds, so today's speed-up reaches almost nobody. The bff half of the fix
already exists (`JointChiSquared`, landed and tested); what is missing is the
ChiSurf-side builder that maps a `GlobalFitModel` onto it. That is this
ticket.

---

## What exists now

### In `imp.bff`

| File | What it is |
|---|---|
| `include/Minimizer.h`, `src/Minimizer.cpp` | MINPACK `lmdif` (`enorm`, `fdjac2`, `qrfac`, `qrsolv`, `lmpar`, `covar`) transcribed line by line, plus ChiSurf's `leastsqbound` bounds transform. Drives free-parameter `Port`s and a `Node` objective. Publishes a covariance from its final `R`. |
| `include/JointChiSquared.h`, `src/JointChiSquared.cpp` | **The grouping.** Members' residuals end to end; chi-square their sum. `add_member(node, "residuals")` creates an input port *linked* to that member's residual output, so `Node::update()` walks the whole tree from one call. |
| `include/ChiSquared.h` | Gained a `residuals` output port (written only when the node has one) and an `update()` override that clears sanitising on its transport. |
| `include/Port.h` | Gained `get_values_ref()` (no-copy read, C++ only), `set_values_array()` (numpy in), `set_sanitize()` (default unchanged). |
| `test/minimizer/test_minimizer.py` | 32 tests. Parity against scipy, bounds, cancellation, covariance. |
| `test/minimizer/test_joint.py` | 13 tests. Grouping, shared parameters, per-member windows, NaN. |
| `test/minimizer/bench_fit.py` | The A/B. Run it before and after. |

### In `chisurf`

| File | What it is |
|---|---|
| `chisurf/core/fitting/minimizer.py` | **New.** `minimize()` (drop-in for `leastsqbound` as `fit.run()` calls it), `graph_objective(fit, model)` (builds the C++ graph or returns `None`), `ResidualNode`, `ProgressObserver`. |
| `chisurf/core/fitting/fit.py` | `Fit.run()` and `FitGroup.run()` call `minimize`. `update_error_estimates` prefers the optimiser's covariance (`_optimiser_covariance`). |
| `chisurf/core/math/optimization/leastsqbound.py` | `_pivot_indices()` — fixes a permuted covariance. |
| `test/fitting/test_graph_fit.py` | 14 tests. What the graph must reproduce, and what it must refuse. |
| `test/fitting/test_leastsqbound_covariance.py` | 9 tests. |

### Numbers to beat

Plain `Fit`, 512 points, three free parameters, interleaved min-of-many
(`test/minimizer/bench_fit.py`):

```
scipy     1.324 ms
director  1.327 ms      <- the C++ optimiser driving a Python residual
graph     0.457 ms      <- 2.90x
```

`FitGroup` on the same data: **1.01x.** That is the gap this ticket closes.

---

## The task

`graph_objective` currently refuses a `GlobalFitModel`, so `FitGroup.run()`
falls back to scipy. Make it build a group graph.

**The mapping is exact, which is why this is tractable:**

| ChiSurf | bff |
|---|---|
| `GlobalFitModel.weighted_residuals` = concatenation of members' | `JointChiSquared` residuals = blocks end to end |
| `GlobalFitModel.parameters` = each member's free parameters, then `global_parameters` | the `Minimizer`'s parameter ports, in that order |
| a member parameter linked to a global one | `Port::set_link` between the private graph's own ports |
| each member's own `xmin`/`xmax`/mask/noise model | each member's own `ChiSquared` (already supported — see `test_members_keep_their_own_window`) |

Read `chisurf/core/models/global_model/globalfit.py`: `weighted_residuals`
is at line ~25, `parameters` at ~177.

### Shape of the change

Extend `graph_objective` in `chisurf/core/fitting/minimizer.py`. The
single-model path there is the template — copy its structure:

1. Detect a `GlobalFitModel` (it has `.fits`).
2. For each member fit, build that member's `Expression -> ChiSquared`
   exactly as the current code does for a plain `Fit`, using the member's own
   `fit.data`, `fit.xmin/xmax`, `fit.mask`, `fit.noise_model`. **Refuse the
   whole group if any member refuses** — a group half in C++ and half in
   Python is not a thing `JointChiSquared` can express as one objective.
3. Reproduce the links. A member parameter that is linked (to a global
   parameter or to another member's) is *not* in that member's
   `model.parameters`, so it will not get its own optimiser port; instead its
   private `Expression` port must `set_link` to the port of whatever it
   follows. Resolve through ChiSurf's `Parameter.link` / `is_linked`.
4. Build the `JointChiSquared`, `add_member` each member's `ChiSquared`.
5. Parameter ports in `GlobalFitModel.parameters` order.
6. Write back: the existing code already writes the answer through the
   ordinary setters and calls `model.update_model()`. For a group that means
   each member's model too — check `GlobalFitModel.update_model`.

### Done when

- A one-member and a two-member `FitGroup` over parse models reach the same
  parameters, `chi2r` and error estimates as the scipy path.
  `test/fitting/test_graph_fit.py` is the pattern to copy — in particular
  `test_the_answer_does_not_move` and `test_the_error_estimates_do_not_move`.
- A group with a **genuinely shared** parameter lands on the compromise, not
  on either member's own answer. `test_the_group_is_not_two_separate_fits` in
  `imp.bff/test/minimizer/test_joint.py` is that test at the C++ level; the
  ChiSurf one should mirror it.
- `chisurf test/fitting` still **932 passed** plus the two known
  pre-existing failures (`test_fit_state`, `test_pcf_experiment` — both fail
  on `HEAD` too; verify by swapping in `git show HEAD:<file>`, **not** by
  `git stash`, see Traps).
- `test/minimizer/bench_fit.py` extended with a `FitGroup` case showing a
  real speed-up rather than 1.01x.

---

## Traps

These each cost real time. None is obvious from the code.

**1. `Parameter.value` has a `_frozen_value` cache, so C++ must not write
ChiSurf's own ports.** Inside `factorgraph.frozen_structure()` — which every
fit runs in — a port written from C++ is *not seen* by the next Python read;
the model would report its pre-fit values, silently. This is why the graph is
**private to the fit**: built from the equation, the data and the current
values, run, then the answer written back through the ordinary setters. Keep
that. (`chisurf/core/parameter.py`, ~line 258.)

**2. Ports sanitise NaN by default, and for a fit that is backwards.**
chinet floors NaN to `np.finfo(float).tiny` so a stored value stays
serialisable. In a fit that reads as *zero* — a good fit for data near zero —
so the optimiser is attracted to the parameters that break the model. Fixed
via `Port::set_sanitize(false)` on the fit's own transport only. **Any new
port carrying a curve or a residual must opt out too.** Note `ChiSquared`
has to do it in an `update()` override, not `evaluate()`, because
`Node::update()` writes the model input first.

**3. The director path is slower than scipy. Do not use it as a fallback.**
1.33 ms against 1.32. Wrapping a Python callback in a C++ loop buys nothing
over wrapping it in a Fortran one. `ResidualNode` is kept because it is the
right tool if a *whole* loop ever needs to be C++ with a Python leaf, but
`minimize` falls back to `leastsqbound`.

**4. Two `ninja` runs in `/Users/tpeulen/dev/imp/cmake-build-arm64` corrupt
it, and the symptom lies.** You get `ld: library 'IMP.algebra-lib' not found`
and a half-linked `_IMP_bff.so` — nothing says "corrupt". Recovery:
`$E/bin/cmake .` then a single `ninja -j8`. A `ninja` that timed out into the
background is **still running**; wait for it. See the board's "Hazards in the
shared tree".

**5. Do not `git stash` in these checkouts.** Both have pre-existing
stashes and long-lived working-tree modifications that are not yours. To A/B
against `HEAD`, copy the file aside and `git show HEAD:<path> > <path>`, then
copy back.

**6. Never pin trajectory parity from a symmetric start.** A sum of
exponentials started at `(1, 1, 1, 1)` sits on its own permutation symmetry:
two *identical* pairs of Jacobian columns, so which of two equal norms the
pivoted QR takes is decided by the last bit. Both optimisers find the same
minimum and label its terms differently. Start off the ridge, e.g.
`(1.0, 0.5, 0.3, 2.0)`.

**7. Bounded parity cannot be exact, and that is not a bug.**
`leastsqbound` maps coordinates with numpy's `sin`/`arcsin`, the C++ with
libm's; `i2e(e2i(x))` round-trips to 2.2e-16 and twenty LM iterations
amplify one ulp to ~1e-7 — inside `xtol = 1.49012e-8` *relative*. scipy's own
answer moves by the same amount when the bounds are merely respelled.
Unbounded parity **is** exact (whole trajectory, same `nfev`).

**8. SWIG.** No `#` comments inside `%pythoncode` — SWIG reads them as its own
preprocessor directives. A director held *across* calls needs
`IMP_SWIG_DIRECTOR`, not a bare `%feature("director")`; the latter keeps only
a weak pointer to the Python proxy and segfaults when it is collected.
`std::vector<std::shared_ptr<Node> >` is not a template this module may name
and leaks if wrapped.

---

## How to build and test

```bash
E=/Users/tpeulen/mambaforge/envs/arm64
B=/Users/tpeulen/dev/imp/cmake-build-arm64

ninja -C $B -j8                      # ONE at a time; ~1 min for bff-only
$E/bin/python -m pytest test/minimizer test/chi2 test/expression \
    test/sampler test/portnode test/factorgraph test/session -q
# 284 passed

cd /Users/tpeulen/dev/chisurf
PYTHONPATH=. $E/bin/python -m pytest test/fitting -q
# 932 passed, 2 known failures, ~3 min

$E/bin/python /Users/tpeulen/dev/imp.bff/test/minimizer/bench_fit.py
```

ChiSurf is **not installed** in that env — it is imported from the sibling
checkout, so either `cd` there or set `PYTHONPATH`. The imp.bff tests that
compare against ChiSurf add the sibling path themselves.

---

## Why this is worth doing, in one measurement

Three attempts, of which two failed, are recorded in `okf/log.md`
2026-09-01 (9). The short version:

- Porting the optimiser alone: **1.02x**. MINPACK's arithmetic is ~4% of a
  fit; the callback is 57%.
- Optimiser + Python residual through a director: **0.99x** — a regression.
- The whole graph: **2.90x**, and reusing the optimiser's covariance instead
  of rebuilding the Jacobian removed a further 32%.

The lesson to carry into this ticket: **the crossing is the cost, and only
removing it entirely helps.** A group graph that still returns to Python once
per iteration for anything will measure like the middle row.
