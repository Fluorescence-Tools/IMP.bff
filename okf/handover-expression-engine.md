# Handover: the expression engine (bff + tttrlib)

Written 2026-08-31. Everything below is measured on this machine (Apple
silicon, arm64 conda env, Python 3.12) unless it says otherwise.

## What this is for

Three consumers want to evaluate a string over arrays, in C++:

1. **chisurf parse models** -- equations from `chisurf/core/models/**/models.yaml`
   (86 of them), evaluated per fit iteration over a curve of 100-4096 points.
2. **tttrlib gating** -- boolean queries over table columns of 1e5-1e7 rows,
   answering with a row mask. Replaces `DataFrame.eval`.
3. **bff fits** -- a whole fit as a C++ node graph (`Expression` ->
   `ChiSquared` -> `Sampler`), so a sampler move costs no interpreter.

## Current state

### Landed and green

- `imp.bff/include/Expression.h` + `src/standalone/Expression.cpp` (1491 lines).
  Built as its own translation unit via `IMP_bff_LIBRARY_EXTRA_SOURCES`, so it
  is not in bff's unity build (that mattered: 91s -> 73s incremental).
- `imp.bff/include/ChiSquared.h` -- the fit objective as a node. 1:1 with
  ChiSurf's `calculate_weighted_residuals`/`get_chi2`, both noise models.
- ~~`imp.bff/include/Table.h`~~ -- **retired 2026-08-31.** It had zero
  consumers outside its own test, `tttrlib::DataStore` already did everything
  it did and more, and it was the last thing keeping the vendored 1.6 MB
  ExprTk alive in this repository. `Expression::compute_pointers()` went with
  it -- Table was its only caller and it was `%ignore`d from Python.
- `tttrlib` `DataStore::select_expression` / `count_expression` -- **ported
  2026-08-31** onto `tttrlib/modules/core/ExpressionEngine.{h,cpp}`, which is
  this engine moved down a layer. Templated on the working type (float32
  columns run in float32) and packing `BitMask` words straight out of the
  block loop. Gating roughly halved: 1M rows 1.6-2.0 ms -> 0.76-0.91 ms.
  ExprTk survives there only as a fallback, and carries the same silent
  multi-argument bug as bff's (T-20260831-13).
- **imp.bff still has its own copy of the engine, so there are currently
  two.** T-20260831-12 deletes bff's and makes it consume tttrlib's; the
  feasibility question is settled (a module may carry its own
  `CMakeModules/Find*.cmake`, as `gsl`, `multifit` and `score_functor` do, so
  nothing needs committing to the `../imp` checkout).
- `Expression::compute_mask()` -- a gate answered as one byte per row,
  emitted from the block loop. Real as of 2026-08-31 (was a stub; see below).
- **1581 bff tests pass** (65 of them the expression suite), 4 xfailed.
  1,000,000 fuzz inputs, zero crashes.

### The engine, as it stands

Parsing is still **ExprTk** (vendored `include/internal/exprtk.h`, MIT), used
only to validate and to evaluate what the vector engine cannot represent.

Evaluation is **our own**: tokenizer -> shunting-yard -> RPN, evaluated
block-wise (512 doubles) with NEON intrinsics. Handles `+ - * / **`, unary
minus, `< <= > >= == !=`, `and or not`, and `abs exp sqrt log log10 sin cos
tan pow min max`.

Optimisations that are in and that earned their place:

- **Block-vectorised evaluation** -- one operation across a block, not a tree
  walk per element.
- **Reused block stack** -- was allocated per call; that malloc was most of
  the cost at 128-1024 points and its removal turned small-n losses into
  1.5-2x wins.
- **Scalar folding** -- a stack slot holding one repeated value stays scalar.
- **Direct-destination scalar kernels** (`simd_add_sv` etc.) -- the earlier
  form computed into the right-hand slot then memcpy'd it down, two wasted
  passes for `0.3+2.0*x`.
- **Constant folding**, and with it constant-power folding -- `**-1` ->
  reciprocal, `**0.5` -> sqrt, `**2` -> multiply. *Every* constant subtree
  is folded, not only a bare literal, which is what finally made the FCS
  exponent fold (see below). A scalar *parameter* used as an exponent takes
  the same kernels at run time.
- **Common subexpression elimination** -- a repeated subtree is computed
  once into a cache slot and loaded back, where recomputing it costs more
  than the two block copies that sharing costs (2026-08-31).
- **Typed boolean stack** -- comparisons write one BYTE per row, not a
  double. This was the last piece of ExprTk's everything-is-a-double model
  and it cost 8x the output traffic.

## Measurements

Curve lengths that actually occur (min-of-11, vs numpy `eval`).
Reproduce with `benchmark/expression_curves.py` -- these used to come from an
ad-hoc script, which is why the earlier column could not be re-checked:

| expression | 128 | 512 | 1024 | 2048 | 4096 |
|---|---|---|---|---|---|
| `exp(-x/1.5)+0.2*exp(-x/4.0)` | 1.80x | 1.43x | 1.21x | 1.09x | 1.01x |
| `sqrt(abs(x))+1/(x*x+1)` | 1.79x | 1.56x | 1.23x | 0.96x | 0.83x |
| `0.3+2.0*x` | 0.82x | 0.80x | 0.74x | 0.62x | 0.54x |
| FCS (below) | **3.2x** | **2.4x** | **1.75x** | **1.29x** | **1.01x** |

FCS = `0.3+1/1.7*(1+x/1.2)**(-1)/sqrt(1+1/2.1**2*x/1.2)`

Those are 2026-08-31 evening, min-of-11, three runs, load average 14. The FCS
row read 1.76x / 0.88x / 0.56x / 0.34x / 0.25x before that evening; the other
three rows are unchanged within run-to-run noise, and their compiled programs
are instruction for instruction the same. Absolute time for FCS at 4096 points
went 48.5 us -> 12.1-12.5 us.

**Independently reproduced** on a second run at load average 9.8, by the
session that opened the ticket rather than the one that did the work:
3.10x / 2.36x / 1.68x / 1.22x / **0.95x**, FCS at 4096 measuring 13.0 us
against the 46.5 us recorded earlier the same day. The first four cells
confirm. **The last one does not quite**: 0.95x here against 1.01x there, so
4096 points is *parity with numpy, not a win* -- the ticket asked for >= 1.0x
and that cell sits on the line rather than over it. Both sessions flagged the
same cell as the doubtful one, which is the useful part. Re-check on an idle
machine before quoting either number.

Three of those four rows reproduce the numbers this handover was written with
to within the noise. **`0.3+2.0*x` does not**: it was recorded at 0.98-1.00x
for 128-1024 and measures 0.72-0.80x here. That row is 1.4-2.0 us end to end,
so it is nearly all fixed cost -- the SWIG call and the malloc of the output
array -- rather than arithmetic, and it is the row most sensitive to how it
was called. Not chased further; the equation is not one anyone fits.

A/B'd against a build with today's type-reconciliation branches stripped out:
**identical within noise on every cell**, so the extra `is_bool` test per
operation costs nothing measurable.

Gating over large columns (vs pandas / vs numpy):

| rows | vs pandas | vs numpy |
|---|---|---|
| 100k | **4.6-4.9x** | 0.40x |
| 1M | **1.5-1.7x** | 0.33x |
| 5M | 0.93-1.11x | 0.33x |

tttrlib's `select_expression` (still ExprTk, writing a `BitMask`) measured
**4.5x pandas** at 200k rows: 0.32ms against 1.43ms.

The mask path, once real (`benchmark/expression_mask.py`, min-of-7, load
average 12 -- so read these as a floor with wide error bars):

| gate | rows | vs double path | vs numpy | vs pandas |
|---|---|---|---|---|
| `x > 0.3` | 100k / 1M / 5M | 1.18 / 1.72 / 2.17x | 0.50 / 0.48 / 0.57x | 9.1 / 3.3 / 2.6x |
| `x > 0.3 and x < 0.7` | 100k / 1M / 5M | 1.06 / 1.34 / 1.46x | 0.45 / 0.50 / 0.55x | 4.7 / 1.9 / 1.5x |
| three terms | 100k / 1M / 5M | 1.11 / 1.31 / 1.32x | 0.56 / 0.51 / 0.59x | 4.5 / 1.8 / 1.3x |
| `(x-.5)**2 + z*z < .25` | 100k / 1M / 5M | 1.07 / 1.31 / 1.29x | 0.85 / **1.37** / **1.56x** | 4.8 / 1.3 / 0.9x |

Two things to read out of that. The mask is **1.3-2.2x the double path**, and
the margin grows with row count, which is what a saving in output bytes should
do. And numpy wins every *pure comparison* and loses the *arithmetic* gate
from 1M rows up -- because `(x-0.5)*(x-0.5) + z*z` costs numpy three temporary
arrays of n doubles while the block evaluator keeps the intermediates in
cache. That is the same effect that makes FCS-shaped equations worth chasing
with common subexpression elimination, seen from the other side.

## Dead ends -- do not repeat these

Each was implemented, measured, and reverted. The reasoning is in the source
as comments.

1. **Fused element-major evaluation** (whole program per element, register
   stack, one pass over memory). Measured **0.05-0.29x of numpy** -- five to
   twenty times *worse*. The per-element switch costs far more than the
   memory traffic it saves. This is the same trap ExprTk falls into. Kept as
   `evaluate_fused()`, unused, with the finding recorded.
2. **Newton-Raphson reciprocal instead of hardware divide.** Accuracy was
   fine (2e-16, machine epsilon) but speed was **unchanged**: the per-lane
   guards needed for zero and non-finite inputs cost as much as the division.
3. **`__restrict` alone.** No measurable change. The loops needed explicit
   intrinsics, not a hint.
4. **Benchmarking with scalars materialised as full columns.** `np.full(n, v)`
   for a scalar parameter makes the engine memcpy a whole array per operand
   and is not what a real model does. It flattered numpy by ~30%.

## Done since: the mask path, and the bugs under it

`compute_mask()` was a placeholder -- it allocated a full
`std::vector<double>`, ran the ordinary double evaluator, and narrowed
afterwards, which is strictly *more* work than the evaluation it was meant to
beat. **Any mask measurement taken before 2026-08-31 measured that stub.**

The real version is what the handover predicted: the block loop already
carries a typed boolean stack, so only the per-block closing store differs.
`evaluate_vector_program` and `evaluate_mask_program` are now two names for
one `evaluate_program(..., double* out, unsigned char* out_mask)`, with
exactly one destination non-null.

Making that path live exposed **three real bugs in the typed stack**, all of
the same shape -- a slot's type was never reconciled at the boundaries:

1. **`is_bool` was never cleared when a push reused a slot.** `a>0 and b<1 or
   c>2` loads `c` into slot 1 while the slot still claimed to hold `b<1`.
   Latent before (nothing read the flag on that path), fatal once the mask
   writer trusts it. Now every `OP_CONST`/`OP_VAR` push clears it.
2. **`and`/`or`/`not` read the boolean stack unconditionally.** `(x>2) and y`
   has a column of doubles on the right, and the combiner read whatever bytes
   the boolean stack happened to hold. Now `booleanise()` casts a numeric slot
   to bytes first -- nonzero is true, as in numpy.
3. **A comparison used as a number read the wrong stack.** `(x>2)*3` computed
   `x*3`: the mask lived in the byte stack while the arithmetic read the
   double stack, which still held `x`. Now `numerify()` widens it back. This
   one was wrong in shipped behaviour, not merely latent.

Truthiness throughout is numpy's `arr.astype(bool)`: **anything not zero is
true**, so a negative operand and a NaN are both true. The stub's `> 0.5`
threshold got both wrong.

Tests: `test/expression/test_expression_mask.py`, 21 tests -- numpy agreement
over a length that is not a multiple of the 512-row block, every comparison
operator, the truthiness table including inf/NaN/denormal, comparisons used as
numbers, and a cross-check that the mask and double evaluators never disagree.

## Open questions

- **Is ExprTk worth keeping at all?** **Answered 2026-08-31: no — drop it,
  after four named gaps.** Full working in `okf/validation/exprtk_fate.md`.
  All 86 shipped equations (67 unique) take the vector path and *none* falls
  back. ExprTk is **45.3% of the shipped dylib's `__TEXT`** (3,186,952 of
  7,041,984 bytes) and takes `Expression.cpp` from 1.8-2.6 s to 22-24 s to
  compile; `Table.cpp` pays it a second time as `exprtk<float>`. Runtime is
  not an argument either way (the validation gate is 2.5-13.9 us of a 55-65 us
  `set_expression()`). The four gaps to close first: (1) a token-adjacency
  syntax check, or dropping the gate lets `1e+` evaluate as `1+e`; (2) `not`'s
  precedence without parentheses, currently unreachable only because ExprTk
  refuses that form; (3) a `program_depth_` cap, since the block stack is
  `depth x 512 x 8` bytes and ExprTk's nesting limit is the only thing
  bounding it; (4) 15 unary + 2 binary numpy-spelled functions.
  - **Correction to this document:** the claim that ExprTk still handles
    "comparisons, booleans" is **stale**, as is the comment at
    `Expression.cpp:1514`. The vector engine has handled both since the typed
    boolean stack landed.
  - **And the fallback is not merely redundant, it is wrong** — see
    T-20260831-10. ExprTk's multi-argument functions collapse to element 0 and
    broadcast, so `hypot(x,y)` returns a constant curve with no error.
    Independently reproduced. That is the strongest argument for removal.
- **`bff::Table` overlaps `tttrlib::DataStore`.** **Answered 2026-08-31:
  retire it.** Full working in `okf/validation/table_vs_datastore.md`. The
  decisive fact is that **Table has zero consumers** outside its own test —
  the migration it was written for already went to DataStore, and
  ndxplorer's `query_mask` calls `select_expression` today. The layering rule
  points to tttrlib on both halves (photon-derived columns in, chisurf gating
  out), so the "consumer wins" tie-break never fires. Deletion is confined to
  imp.bff. Note that `Expression::compute_pointers` exists only for Table and
  goes with it; it is ~30 lines and trivially reinstated.
  - Consequence worth keeping separate: Table's float64 path is currently the
    only in-tree caller of the fast block engine, while DataStore is still on
    ExprTk. imp.bff has **no** C++ dependency on tttrlib, so DataStore cannot
    call `bff::Expression` — which is what T-20260831-08 is for.
- **A stack slot cannot alias a column.** Every `OP_VAR` memcpy's a 4 kB
  block even when the value is only read, and in the FCS row two of about
  eleven block passes are that copy. Making a slot able to point at the
  caller's column, and materialising it only where a kernel writes in place,
  would take roughly a fifth off the arithmetic-light rows -- but it means
  every unary and scalar-binary kernel gaining a separate source and
  destination. Not attempted; the largest remaining lead.

## Done since: constant folding and CSE (2026-08-31 evening, T-20260831-05)

The ticket was written to add **common subexpression elimination**, on the
reading that FCS loses because it computes `x/1.2` twice. It does not. In
`0.3+1/1.7*(1+x/1.2)**(-1)/sqrt(1+1/2.1**2*x/1.2)` the second occurrence is
`(1/2.1**2*x)/1.2`, a different subtree, and the only thing genuinely shared
is the leaf `x`. **The whole of that row's loss was `**(-1)`.**

`x**(-1)` tokenises as a *negated* constant -- `OP_CONST 1`, `OP_FUN F_NEG` --
and the constant-power fold only looked for a bare `OP_CONST` under the power.
So the exponent stayed a buffer and every element paid a `std::pow()` call.
Measured in isolation: `x**(-1)` cost 39.7 us at 4096 points against 5.1 us
for `1/x`, ~9.7 ns a point of pure libm. numpy never pays it, because
`fast_scalar_power` in `number.c` rewrites `arr ** -1` into `np.reciprocal`
before the ufunc is reached -- so this was not numpy being fast, it was us
doing the thing numpy specifically avoids.

Two changes, both in `compile_vector_program`:

1. **Fold every constant subtree**, not only a bare literal, before the
   constant-power fold runs. Folding calls the same `apply_scalar_fun` /
   `apply_binary_scalar` the evaluator's own scalar path would have called,
   so it is bit-identical -- only *when* the arithmetic happens changes.
   `not` is deliberately excluded: the block evaluator gives it numpy's
   truthiness through `booleanise()` while `apply_scalar_fun` still carries
   ExprTk's `> 0.5` rule, and folding it would make `not 0.3` disagree with
   itself. Comparisons are excluded too -- folding one to a plain 1.0 would
   lose the fact that its slot is a mask.
2. **A run-time fast path for `vector ** scalar`**, so an exponent that
   arrives as a fit parameter rather than a literal takes the same
   reciprocal/sqrt/square kernels.

That alone took FCS from 0.25x to 1.00x at 4096 points.

**CSE was implemented anyway, and does pay -- just not on that row.** Half of
ChiSurf's FCS catalogue really does repeat vector subtrees: `4*D*x/w_r**2`
appears three times in the three-state FRET-FCS model, `exp(-(k_12+k_21)*x)`
three times. The pass hash-conses the RPN into a DAG (a node keyed by opcode
and child ids that has been seen before *is* the earlier node, the operands
being pure), then takes only the candidates that pay: a cached slot costs one
block copy to fill and one to read, so a subtree must cost more than that
before it is shared. `x/1.2` does not qualify; `exp(-x/tau)` does. Constants
are compared bitwise, so two subtrees are shared when they are the same text,
never when they merely compare equal.

A/B against a build with the cache capacity set to zero, same binary
otherwise, at 4096 points:

| equation | CSE off | CSE on |
|---|---|---|
| `exp(-x/1.5)*exp(-x/1.5)+exp(-x/1.5)` | 41.8 us | **17.5 us** |
| 3-state FRET-FCS (29 repeated subtrees) | 81.4 us | **46.6 us** |
| flow FCS (2 repeated) | 43.7 us | 42.5 us |
| 2-component FCS (6 repeated, all scalar) | 54.5 us | 53.7 us |
| the benchmark's FCS row | 12.4 us | 12.1 us |

So: **2.4x where subtrees really repeat, neutral where they do not**, and the
four benchmark rows are unmoved. Cost-blind CSE would have been a
pessimisation -- caching `x/1.2` replaces one multiply with two block copies.

CSE meets the typed stack at exactly the boundary the 2026-08-31 bugs were
about. `OP_SAVE` copies a slot *with* its type (byte mask, folded scalar or
doubles) and `OP_LOADC` restores all of it, because a load has to set the
flags a push would have set. `test_expression.py::SharedSubexpressionTests`
pins the three shapes that go wrong otherwise: a cached comparison combined
with `and`, a cached comparison used as a number, and a cached subtree under
`not`. 120,000 grammar-fuzz cases across three seeds agree with numpy and
with their own masks; 200,000 structural cases survive.

## Things that will bite you

- The **arm64 conda env's tttrlib is a scikit-build editable install**, which
  registers a `ScikitBuildRedirectingFinder` in `sys.meta_path`. That runs
  *before* `sys.path`, so `PYTHONPATH` cannot override it. To test a tttrlib
  change: `pip install -e /Users/tpeulen/dev/tttrlib --no-build-isolation
  --no-deps`. Copying build artifacts over the installed package breaks the
  import (SIGKILL) because the dylibs go out of sync.
- **tttrlib's SWIG target does not depend on `DataStore.h`.** A header change
  does not regenerate the wrapper, so new methods silently do not appear in
  Python. Delete `build_new/ext/CMakeFiles/tttrlib.dir/tttrlibPYTHON_wrap.cxx`
  to force it. Worth fixing properly.
- **`tttrlib/build/` had no `CMAKE_BUILD_TYPE`** -- everything built there was
  unoptimised. I set it to Release. Any benchmark taken in that directory
  before 2026-08-31 is meaningless. This one mistake cost several rounds of
  chasing phantom slowness.
- **This machine has been at load average 19-33** from a concurrent build for
  most of the session. Even min-of-11 timings swung 4x between runs. Re-measure
  on an idle machine before trusting any number here.
- `imp.bff`'s working tree carries unrelated WIP that another process edits
  live; builds failed several times mid-session on `AV.h`, `PathMap.h`,
  `ProbeParticle`. Retrying usually works.

## Test inventory

- `test/expression/test_expression.py` -- 35 tests: precedence (`2**3**2`=512,
  `-2**2`=-4), all 86 shipped equations compile and match numpy to 1e-10,
  compile-once, batch, shared subtrees keeping their slot type, and every
  spelling of a constant exponent.
- `test/expression/test_expression_mask.py` -- 21 tests: the gate path, its
  truthiness rule, and mask-vs-double agreement (see above).
- `test/expression/test_expression_robustness.py` -- 13 tests: malformed input,
  `is_supported` never raises, deep nesting, hostile floats (inf/nan/1e308/
  denormal), empty and mismatched columns, reuse after a failed parse, 20k
  in-suite fuzz.
- `test/chi2/`, `test/table/`, `test/sampler/`, `test/abtest/` (differential
  against the frozen Python chinet), `test/portnode/`, `test/session/`,
  `test/factorgraph/`.
- `test/expression/fuzz_expression.py` -- the fuzz harness, now in the repo
  (was `/tmp/bang.py`). Two modes. **structural** is the original: random
  characters at the parser, 200k per seed, and it has never crashed it -- but
  193,951 of 200,000 inputs are *refused*, so only 2.8% ever reach the
  evaluator. **grammar** builds a random expression tree and renders it twice,
  once for the engine and once for numpy, and demands agreement; it reaches
  the evaluator on 92% of cases and cross-checks the byte mask against the
  doubles every time. 75,000 valid expressions over five seeds, zero
  unexplained failures. Explainable disagreements are counted, not hidden:
  `known-divergence` is T-20260831-09 below, `ill-conditioned` is a case where
  one ulp of input moves numpy's own answer more than the tolerance.
  Deterministic in the seed, and it prints the seed to reproduce with.

**Known bug it found: `min`/`max` are not commutative under NaN.**
`min(y, nan)` is `y`, `min(nan, y)` is `nan` -- a ternary `(a > b) ? b : a`,
which is neither numpy's rule nor C's `fmin`, while this header promises
numpy's. Narrow (a NaN reaching `min` means the fit is already broken) but
real, since the parity suite advertises numpy agreement. Ticket
**T-20260831-09**; it collides with the CSE ticket, which owns the file.

## How to run things

```bash
# build (bff)
PATH=~/mambaforge/envs/arm64/bin:$PATH \
  ninja -C ~/dev/imp/cmake-build-arm64 IMP.bff-lib IMP.bff-python

# test (bff)
PYTHONPATH=~/dev/imp.bff/test:~/dev/imp/cmake-build-arm64/lib:~/dev/chisurf \
  ~/mambaforge/envs/arm64/bin/python -W ignore -m pytest ~/dev/imp.bff/test -q

# tttrlib
cd ~/dev/tttrlib/build_new && cmake --build . -j8
~/mambaforge/envs/arm64/bin/pip install -e ~/dev/tttrlib --no-build-isolation --no-deps
```
