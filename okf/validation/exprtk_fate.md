---
okf_version: "0.2"
title: "ExprTk's fate: drop it, after four named gaps"
description: "The vector parser already handles 86 of 86 shipped equations and every comparison and boolean operator. ExprTk is 45.3% of libimp_bff's code and ~45 s of compile time across two translation units, its fallback is silently wrong for every multi-argument function over a vector, and its parser refuses valid input the vector parser gets right. Four gaps stand between here and deleting it."
tags: [validation, imp.bff, expression, exprtk, build, performance, open]
---

# ExprTk's fate

**Date:** 2026-08-31 · **Ticket:** T-20260831-06 · **Verdict: drop it, after
the four gaps in [What has to close first](#what-has-to-close-first).**

Read-only investigation. Nothing was built and nothing in `src/` was touched —
another agent held the build lock. Compile-time and size figures come from
compiling the same translation units with the same flags outside the build
tree; everything else is measured against the module already in
`~/dev/imp/cmake-build-arm64/lib`.

## The short version

ExprTk is no longer the evaluator, and it is not the parser either — it is a
*gate* in front of the parser. It costs **45.3% of `libimp_bff`'s code**, about
**45 s of compile time** (it is instantiated twice: `<double>` in
`Expression.cpp`, `<float>` in `Table.cpp`), it **refuses valid input** the
vector parser handles correctly, and the fallback it exists to provide is
**silently wrong for every multi-argument function applied to a vector**.

What it does still buy is one thing: it is the only syntax check in the engine.
`compile_vector_program`'s sole well-formedness test is a stack-depth
invariant, which lets `1e+` through as `1 + e = 3.718`. That is a 30-line fix,
measured below.

## 1. Coverage: 86 of 86

Every equation ChiSurf ships takes the vector path. None falls back.

| | count |
|---|---|
| `equation:` entries under `chisurf/core/models/**/*.yaml` | **86** (67 unique) |
| compiled to a vector program | **86** |
| evaluated by ExprTk | **0** |
| refused by the gate | **0** |

Method: `has_compiled_plan()` read **immediately after `set_expression()` and
before any `compute()`**. That ordering is load-bearing. `compile()` sets
`impl_->ready` only when `compile_vector_program()` succeeds — but on the
ExprTk path the *first* `compute()` calls `build_plan()`, which also sets
`impl_->ready = true`. After an evaluation the flag no longer distinguishes the
two engines. (I got this wrong once mid-investigation and it inverted three
rows of a table.)

## 2. What actually falls back

I probed 78 constructs. The vector path takes: all four arithmetic operators,
`**` with right-associativity, unary minus at Python's precedence, `abs exp
sqrt log log10 sin cos tan pow min max`, `pi`, `e`, every numeric literal form
including `1e-3` and `.5`, **all six comparisons**, **`and` / `or` /
`not(...)`**, and the `&` `|` `~` rewrites.

> The handover and the comment at `Expression.cpp:1514` both say "ExprTk stays
> for everything else (comparisons, booleans)". **That is stale.** Comparisons
> and booleans are on the vector path — the typed boolean stack made them so.

The fallback set is exactly 27 ExprTk functions plus a little syntax:

| falls back to ExprTk | |
|---|---|
| functions (27) | `floor ceil round sgn erf erfc frac trunc clamp inrange root hypot logn log2 expm1 log1p deg2rad rad2deg asin acos atan atan2 sinh cosh tanh avg sum` |
| control flow | `if(c,a,b)`, `c ? a : b` |
| other | `%`, and the constants `inf epsilon true false` |

ExprTk *refuses* `xor(…)`, `nand(…)`, `nor(…)` (`ERR236 — Invalid use of
reserved symbol`) even though `is_reserved()` advertises them, so those are not
a fallback at all.

### The fallback is wrong wherever it is not redundant

This is the finding that decides the question. Evaluated over `x = [1,2,3,4]`
against numpy:

| equation | path | engine | numpy | |
|---|---|---|---|---|
| `floor(x)`, `ceil(x/2)`, `trunc(x)`, `tanh(x)`, `log2(x)`, `atan(x)`, `sinh(x)` | exprtk | — | — | ✅ agree |
| `hypot(x,y)` | exprtk | `[1.118 1.118 1.118 1.118]` | `[1.118 2.5 3.905 5.315]` | ❌ |
| `atan2(x,y)` | exprtk | `[0.464 0.464 0.464 0.464]` | `[0.464 0.644 0.695 0.719]` | ❌ |
| `clamp(1.5,x,3.0)` | exprtk | `[1.5 1.5 1.5 1.5]` | `[1.5 2 3 3]` | ❌ |
| `root(x,2)` | exprtk | `[1 1 1 1]` | `[1 1.414 1.732 2]` | ❌ |
| `if(x>2,1,0)` | exprtk | `[0 0 0 0]` | `[0 0 1 1]` | ❌ |
| `inrange(1.5,x,3.0)` | exprtk | `[0 0 0 0]` | `[0 1 1 0]` | ❌ |
| `round(x)` on `[0.5 1.5 2.5 3.5]` | exprtk | `[1 2 3 4]` | `[0 2 2 4]` | ❌ (half-away vs banker's) |

The pattern: `build_plan` binds each variable as an ExprTk *vector*, and
ExprTk's unary functions map over it while its **multi-argument functions
collapse to element 0 and broadcast the scalar answer**. Nothing tests this,
because no shipped equation reaches it.

### And most of the rest is unreachable by construction

`chisurf/core/models/parse/parse.py:143` evaluates the equation with
`from numpy import *` in scope. Of ExprTk's 27 extras, **9 do not exist in
numpy's namespace at all** (`sgn erf erfc frac clamp inrange root logn avg`),
so a user typing them would get a `NameError` from ChiSurf and a *number* from
this engine — a divergence, not a feature. `sum` exists but is a *reduction*
there, and the engine returns something else. (`min`/`max` are the same trap on
the vector path: `min(x,1)` gives elementwise `minimum`, while
`from numpy import *` makes `min` the reduction. Separate issue, not ExprTk's.)

**Net: the fallback delivers a correct, reachable answer for 15 unary
functions, and nothing else.**

## 3. The cost

Flags taken verbatim from `build.ninja:26452-26457` (`-w -O3 -DNDEBUG -arch
arm64 -fPIC`, same includes). The comparison is against a **30-line
API-compatible stub** standing in for `include/internal/exprtk.h` — the same
source, same flags, only the header swapped. Machine under load, so read the
absolute times as an upper band; the ratios are stable across runs.

| | with ExprTk | with the stub | ExprTk's share |
|---|---|---|---|
| `Expression.cpp` `-O3 -c` | **22.0 / 22.1 / 24.6 s** | 1.8 / 1.9 / 2.6 s | ~20–23 s (**91%**) |
| `Expression.cpp` object | 7,965,696 B | 99,888 B | **98.7%** |
| … `__TEXT` | 2,459,847 B | 65,693 B | 97.3% |
| … dead-stripped into a dylib alone | 6,267,744 B | 138,432 B | 45× |
| `Table.cpp` `-O3 -c` (`exprtk<float>`) | **23.2 / 24.3 s** | 1.1 / 1.2 s | ~22 s |
| `Table.cpp` object | 7,916,264 B | 49,776 B | 159× |

Where the compile time goes, for `Expression.cpp` with the real header:
`-fsyntax-only` **2.9–3.3 s**, `-O0 -c` **5.9 s**, `-O3 -c` **22–26 s**. So it
is not the 46,766 lines of header being *read* — it is `-O3` optimising the
~2.4 MB of template code the header instantiates. A `-O1` build of that one
file would recover most of it, if anyone wanted a stopgap.

In the shipped library, summing symbol sizes by address delta over all 52,751
`__TEXT` symbols of `libimp_bff.0.dylib`:

```
total __TEXT   7,041,984 B
exprtk         3,186,952 B   (45.3%)
```

909 of the 3,384 exported symbols are ExprTk's. The handover's "1.6 MB" is the
header on disk; the shipped cost is **twice that**, because it is instantiated
twice.

**Runtime is not an argument, and should not be used as one.** ExprTk's own
parse, benchmarked standalone (min of 200): FCS 13.9 µs, two-exponential
8.2 µs, a gate 5.3 µs, `0.3+2.0*x` 2.5 µs. A `set_expression()` from Python
costs 55–65 µs end to end, so the gate is 5–25% of a call that happens once per
equation. Removing it saves build time and binary size, not fit time.

## 4. The three consumers

None needs anything only ExprTk provides.

1. **chisurf parse models** — 86/86 on the vector path (§1). Nothing in
   `chisurf/` imports `bff.Expression` yet, so there is no installed base to
   break. Its reference semantics is numpy's, and numpy offers nothing the
   vector engine lacks that the fallback gets right (§2).
2. **tttrlib gating** — `DataStore.cpp` includes its **own** vendored copy
   (`tttrlib/thirdparty/exprtk/exprtk.hpp`), not bff's. bff's decision does not
   touch it. The three queries the design note benchmarks —
   `(g>2)&(r<10)`, `g>5|b<0.1`, `(g-b)/(r-b)>0.3` — are all vector-path here,
   so the engine tttrlib would adopt does not need ExprTk either.
3. **bff fits** — `Expression → ChiSquared → Sampler` is plain arithmetic.

The one internal consumer that *does* need work is **`bff::Table`**, which
instantiates `exprtk<float>` for its float32 column path
(`Table.cpp:98,134`). The handover already proposes retiring Table in favour
of `tttrlib::DataStore`; if that happens, half the ExprTk cost above leaves
with it.

## 5. What the gate costs in correctness — measured

ExprTk's refusals, with the errors it gives:

| input | accepted up to | error past it |
|---|---|---|
| `(((…1.0…)))` | **199** deep | `ERR000 — Current stack depth 401 exceeds maximum allowed stack depth of 400` |
| `1+(1+(1+…))` | **132** deep | same ERR000 |
| `abs(abs(…))` | **199** deep | same ERR000 |
| `x+x+…+x` | **10,002** terms | `ERR018 — Expression depth of 15998 exceeds maximum allowed expression depth of 10000` |
| `not x>2` | — | `ERR029 — Expected a '(' at start of function call to 'not'` |
| `----------x` | 3 signs | `ERR008 — Invalid token sequence: '+' and '+'` |

Because `compile()` runs ExprTk *before* `compile_vector_program()`, **these are
the engine's limits**, not ExprTk's.

The vector parser handles all of them, and handles them correctly. Verified by
evaluating the RPN it emits against Python:

| input | vector parser | Python |
|---|---|---|
| 200- and 500-deep parens | 1.0, 1.0 | 1.0, `SyntaxError` |
| 200- and 400-deep right-nesting | 201.0, 401.0 | 201.0, `SyntaxError` |
| 250-deep `abs(` | 1.0 | `SyntaxError` |
| 20,000- and 60,000-term sum (`x=3`) | 60000.0, 180000.0 | `RecursionError` |
| `not x>2` (`x=3`) | 0.0 | 0.0 |

### How I know that, given the gate is always in front

The vector parser's verdict on ExprTk-refused input is not observable from
Python. So I ported `tokenize` + `compile_vector_program` (with the
constant-power fold and the depth walk) to Python and validated the port
against the built module on every input where a comparison is possible — the 86
catalogue equations, the suite's 60-entry `MALFORMED` list, and 20,000 fuzz
strings from `test_expression_robustness.py`'s own alphabet:

```
comparable inputs 1,753 · agreements 1,753 · disagreements 0
```

(One class of disagreement had to be fixed in the port first, and it is a real
property of the C++: `strtod` accepts hex float literals, so `0x1F` tokenises
as the constant 31.)

Every "the vector parser would…" claim above is that port plus source reading,
not the shipped binary.

## 6. The honest counterweight: what the gate catches

Over 20,063 adversarial strings, **485 are accepted by the vector parser and
refused by ExprTk**. Some of those are the `not` forms ExprTk is simply wrong
about. Most are junk, and the worst are *silently misparsed* rather than merely
accepted:

- `1e+` → `strtod` takes the longest valid prefix `1`, then `e` is read as
  Euler's number and `+` as a binary operator: **`1 + e = 3.718`**.

`compile_vector_program`'s only well-formedness test is the stack-depth
invariant at its tail (`if (depth < 1) return false;` and
`if (depth != 1) return false;`). That catches `1 2`. It does not catch `1e+`,
because that *is* a well-formed program — of the wrong expression.

I prototyped the missing check: a token-adjacency pass with one bit of state
(*do I expect a value or an operator next?*), 30 lines. Measured:

| | |
|---|---|
| catalogue equations it accepts | **86 / 86** |
| of the 485 leaks, it removes | **433 (89.3%)** |
| of the 1,118 inputs both ExprTk and the vector parser accept, it wrongly rejects | **1** |
| the 52 leaks it still admits | all contain `~` / `not` — i.e. the form ExprTk is wrong about |

The single false rejection is `e^m82n|o0l9-bbrq//1.r`: two binary divides in a
row. ExprTk accepts that; refusing it is the better answer, since `//` means
floor division in the Python syntax these equations are written in.

## What has to close first

Four gaps. Together they are the ticket list for deleting ExprTk.

**Gap 1 — a well-formedness check in `compile_vector_program`.** Token
adjacency, one bool of state, ~30 lines before the shunting-yard loop.
*Measured*: closes 433 of 485 leaks, accepts all 86 catalogue equations, one
false rejection on 1,118 valid inputs (and that one is `//`, which should be
rejected). Without this, dropping the gate makes `1e+` evaluate to 3.718.

**Gap 2 — `not`'s precedence is wrong without parentheses.**
`not x>2 or y>0.1` gives **0.0**; Python gives **1.0**. The pending `F_NOT` is
only flushed when the token stream ends (the `if (stack.back().second == F_NOT)
break;` in the pending-unary loop), so it binds the entire expression instead
of the comparison to its right. `not(...)` and `~(...)` are correct — verified
against numpy on five shapes. This is **unreachable today** precisely because
ExprTk's `ERR029` refuses the unparenthesised form, and becomes reachable the
moment the gate goes.

**Gap 3 — cap the program depth.** `program_depth_` sizes the block stack at
`depth × 512 × 8` bytes with no bound (`Expression.cpp:1215`). ExprTk's
`ERR000` currently caps depth near 200, i.e. 819 KB. A 100,000-deep
right-nested expression would ask for **410 MB** of doubles plus 51 MB of the
boolean stack. Pick a limit, refuse past it, and put it in
`test_expression_robustness.py`.

**Gap 4 — 15 unary functions and 2 binary ones.** Exactly the subset the
fallback delivers correctly *and* numpy spells the same way:

> unary: `floor ceil round trunc log2 expm1 log1p deg2rad rad2deg sinh cosh
> tanh asin acos atan`
> binary: `hypot atan2` — these are the two the ExprTk path gets **wrong**
> today (§2), so implementing them is a bug fix, not a port.

One line each in `function_id()`, one case each in the block kernel. Note
`round` must be numpy's banker's rounding, not ExprTk's half-away-from-zero.

**Explicitly not gaps** — do not port these: `sgn clamp root frac inrange logn
avg sum erf erfc`, `if(…)`, `?:`, `%`, `inf epsilon true false`. None appears
in the 86; nine of them do not exist in numpy's namespace, so a ChiSurf user
cannot reach them without the engine and Python disagreeing; and `sum` means
something else there. Several are silently wrong on the fallback path today, so
removing them fixes more than it breaks.

## Then

Delete, in this order: the ExprTk branches in `compute`, `compute_columns`,
`compute_pointers` and `compute_mask`; `build_plan()`; the `table` / `compiled`
/ `storage` / `result` members of `Expression::Impl`; the gate in `compile()`;
`Table.cpp`'s float path (or `Table` itself, per the handover); the two
`#include <IMP/bff/internal/exprtk.h>` lines; and finally
`include/internal/exprtk.h`. `set_expression()` and `is_supported()` must keep
their contract — `std::domain_error` / `ValueError` on bad input — which is what
Gap 1 exists to preserve.

Expected recovery, from §3: **~45 s of compile time** and **~3.2 MB (45.3%) of
`libimp_bff`'s code**. Both unverified until someone actually removes it — I
could not build.

## Caveats

- **The build lock was held by another agent throughout.** Nothing here was
  compiled into the build tree, and the removal itself is unmeasured.
- **`Expression.cpp` changed under me**, 1,561 → 1,891 lines, gaining
  `eliminate_common_subexpressions()` and the `OP_SAVE`/`OP_LOADC` opcodes
  (the concurrent CSE ticket). That work is *post-parse*: the tokenizer, the
  function table and the shunting-yard loop are unchanged, and the depth walk
  simply learned two more opcodes. The accept/reject set — which is all this
  note is about — is unaffected. The Python module measured in §1, §2 and §5 is
  the build as of session start.
- Machine load was high (the handover warns about this; it was still true).
  Compile times are an upper band; the with/without ratios held across every
  repeat.
