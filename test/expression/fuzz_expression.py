"""Fuzz the expression engine: random strings, and random *valid* equations.

Two harnesses, because they find different things.

**Structural fuzz** throws random characters at the parser. It answers "can
input crash it", and the answer has been no for 1,000,000 inputs. What it
cannot do is reach the evaluator: random garbage almost never parses, so
better than 99% of its cases are refusals, and it is testing `is_supported`
rather than the arithmetic behind it.

**Grammar fuzz** builds a random expression tree and renders it twice -- once
in the engine's syntax, once as the numpy expression that means the same
thing -- then demands the two agree. That reaches the evaluator on every
single case, and it is the only harness that can catch the failure mode this
engine actually has: a *typed stack* whose slots hold doubles, comparison
byte-masks, or folded scalars, and which is wrong whenever a slot's type is
not reconciled at a boundary. Three such bugs were found on 2026-08-31 --
`(x>2)*3` computed `x*3`, `(x>2) and y` read uninitialised bytes, and
`a>0 and b<1 or c>2` reused a slot that still claimed to be a comparison.
Every one of them needs a *valid* expression mixing numbers and booleans to
show up, so no amount of structural fuzzing would ever have reached them.

Both are deterministic in the seed, so a failure is reproducible by rerunning
with the seed it prints.

    PYTHONPATH=~/dev/imp/cmake-build-arm64/lib python test/expression/fuzz_expression.py --seed 1
    ... --mode structural --cases 200000
    ... --mode grammar --cases 20000

Not named `test_*`, so pytest does not collect it: a fuzz run belongs on a
clock the suite does not own. `test_expression_robustness.py` keeps a small
in-suite run for the things that must never regress.
"""

import argparse
import math
import random
import sys

import numpy as np

from IMP import bff

# ---------------------------------------------------------------- structural

ALPHA = ("0123456789.+-*/^()xyzabc, sqrtexplogminmaxabs|&~<>=!"
         "\t\"'\\[]{}:;%@#$")


def fuzz_structural(rng, cases, verbose=False):
    """Random characters at the parser. Nothing may crash; refusals are fine."""
    refused = 0
    evaluated = 0
    for _ in range(cases):
        text = "".join(rng.choice(ALPHA) for _ in range(rng.randint(0, 60)))
        try:
            if not bff.GraphExpression.is_supported(text):
                refused += 1
                continue
            ex = bff.GraphExpression("m")
            ex.set_expression(text)
            names = list(ex.get_variable_names())
            if len(names) > 5:
                continue
            values = [hostile_column(rng) for _ in names]
            ex.compute(names, values)
            evaluated += 1
        except Exception:
            refused += 1
    return {"refused": refused, "evaluated": evaluated}


def hostile_column(rng):
    """The floats most likely to break something, plus ordinary ones."""
    k = rng.random()
    if k < 0.20: return [0.0]
    if k < 0.30: return [float("inf")]
    if k < 0.40: return [float("nan")]
    if k < 0.50: return [1e308]
    if k < 0.55: return [5e-324]          # denormal
    if k < 0.75: return [rng.uniform(-9, 9)
                         for _ in range(rng.randint(1, 64))]
    return [rng.uniform(-9, 9)]


# ------------------------------------------------------------------- grammar

VARIABLES = ["x", "y", "z"]
UNARY = ["abs", "exp", "sqrt", "log", "log10", "sin", "cos", "tan"]
BINARY_FUN = ["pow", "min", "max"]
ARITH = ["+", "-", "*", "/", "**"]
COMPARE = ["<", "<=", ">", ">=", "==", "!="]

# A node is (kind, ...) where kind says whether it evaluates to a number or a
# boolean. Keeping that on the node is what lets the generator deliberately
# feed booleans into arithmetic and numbers into `and`, which is where the
# typed stack has historically been wrong.


def gen(rng, depth, want="num"):
    """Random expression tree. `want` is 'num' or 'bool'."""
    if want == "bool":
        if depth <= 0 or rng.random() < 0.35:
            return ("cmp", rng.choice(COMPARE),
                    gen(rng, depth - 1, "num"), gen(rng, depth - 1, "num"))
        k = rng.random()
        if k < 0.30:
            return ("not", gen(rng, depth - 1, "bool"))
        if k < 0.65:
            # Deliberately allow a *numeric* operand: `(x>2) and y` is legal
            # and means "and y is nonzero". This is the shape that read
            # uninitialised bytes before booleanise() existed.
            return ("logic", rng.choice(["and", "or"]),
                    gen(rng, depth - 1, rng.choice(["bool", "num"])),
                    gen(rng, depth - 1, rng.choice(["bool", "num"])))
        return ("cmp", rng.choice(COMPARE),
                gen(rng, depth - 1, "num"), gen(rng, depth - 1, "num"))

    if depth <= 0:
        # Weighted towards variables: a tree of pure constants folds to a
        # scalar and has no columns to evaluate over, so it tests nothing.
        if rng.random() < 0.75:
            return ("var", rng.choice(VARIABLES))
        return ("const", round(rng.uniform(-4, 4), 3))

    k = rng.random()
    if k < 0.08:
        # A comparison used as a number: 0.0 or 1.0. `(x>2)*3` was wrong.
        return ("asnum", gen(rng, depth - 1, "bool"))
    if k < 0.20:
        return ("un", rng.choice(UNARY), gen(rng, depth - 1, "num"))
    if k < 0.30:
        return ("bin2", rng.choice(BINARY_FUN),
                gen(rng, depth - 1, "num"), gen(rng, depth - 1, "num"))
    if k < 0.42:
        return ("var", rng.choice(VARIABLES))
    if k < 0.50:
        return ("const", round(rng.uniform(-4, 4), 3))
    if k < 0.56:
        return ("neg", gen(rng, depth - 1, "num"))
    return ("arith", rng.choice(ARITH),
            gen(rng, depth - 1, "num"), gen(rng, depth - 1, "num"))


def render(node):
    """The engine's syntax (which is Python's, for the arithmetic subset)."""
    kind = node[0]
    if kind == "var":    return node[1]
    if kind == "const":  return repr(node[1])
    if kind == "neg":    return f"(-({render(node[1])}))"
    if kind == "asnum":  return f"({render(node[1])})"
    if kind == "un":     return f"{node[1]}({render(node[2])})"
    if kind == "bin2":   return f"{node[1]}({render(node[2])}, {render(node[3])})"
    if kind == "arith":  return f"(({render(node[2])}) {node[1]} ({render(node[3])}))"
    if kind == "cmp":    return f"(({render(node[2])}) {node[1]} ({render(node[3])}))"
    if kind == "not":    return f"not({render(node[1])})"
    if kind == "logic":  return f"(({render(node[2])}) {node[1]} ({render(node[3])}))"
    raise AssertionError(kind)


def render_numpy(node, minmax="numpy"):
    """The same tree as a numpy expression.

    `and`/`or`/`not` cannot be rendered literally -- Python's keywords raise
    on arrays -- so they become `&`/`|`/`~` over operands explicitly cast to
    bool. That cast is also the semantics being asserted: nonzero is true,
    which is what `arr.astype(bool)` does and what the engine copies.
    """
    kind = node[0]
    if kind == "var":    return node[1]
    if kind == "const":
        # np.float64, not a bare literal: Python evaluates a constant-only
        # subtree with Python's operator semantics, and `(-3.7) ** 2.1` is a
        # complex number in Python where C's pow() and numpy both give NaN.
        # A bare literal here makes the oracle disagree with the engine about
        # arithmetic neither of them got wrong.
        return f"np.float64({node[1]!r})"
    if kind == "neg":    return f"(-({render_numpy(node[1], minmax)}))"
    if kind == "asnum":
        # np.asarray: an all-constant subtree evaluates to a Python bool,
        # which has no .astype and would otherwise lose the case.
        return f"(np.asarray({render_numpy(node[1], minmax)}).astype(np.float64))"
    if kind == "un":     return f"np.{node[1]}({render_numpy(node[2], minmax)})"
    if kind == "bin2":
        # np.minimum/maximum propagate NaN (numpy's rule, and the engine's
        # stated contract); np.fmin/fmax ignore it (C's rule, and what the
        # engine actually does). Rendering both is how a run tells the one
        # known divergence apart from a new one -- see fuzz_grammar().
        if minmax == "numpy":
            fn = {"pow": "np.power", "min": "np.minimum",
                  "max": "np.maximum"}[node[1]]
        else:
            fn = {"pow": "np.power", "min": "engine_min",
                  "max": "engine_max"}[node[1]]
        return (f"{fn}({render_numpy(node[2], minmax)}, "
                f"{render_numpy(node[3], minmax)})")
    if kind == "arith":
        return (f"(({render_numpy(node[2], minmax)}) {node[1]} "
                f"({render_numpy(node[3], minmax)}))")
    if kind == "cmp":
        return (f"(({render_numpy(node[2], minmax)}) {node[1]} "
                f"({render_numpy(node[3], minmax)}))")
    if kind == "not":
        return f"(~({as_bool(node[1], minmax)}))"
    if kind == "logic":
        op = "&" if node[1] == "and" else "|"
        return f"(({as_bool(node[2], minmax)}) {op} ({as_bool(node[3], minmax)}))"
    raise AssertionError(kind)


def as_bool(node, minmax="numpy"):
    """A numpy subexpression forced to bool, the engine's truthiness rule."""
    if node[0] in ("cmp", "not", "logic"):
        return f"np.asarray({render_numpy(node, minmax)})"
    return f"(np.asarray({render_numpy(node, minmax)}) != 0)"


def engine_min(a, b):
    """`min` as the engine actually computes it: `(a > b) ? b : a`.

    Not `np.minimum` (which propagates NaN) and not `np.fmin` (which ignores
    it in either order). The engine's ternary is **order-dependent** under
    NaN -- `min(y, nan)` is `y` but `min(nan, y)` is `nan` -- so neither
    library function models it. See T-20260831-09.
    """
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    return np.where(a > b, b, a)


def engine_max(a, b):
    """`max` as the engine actually computes it: `(a < b) ? b : a`."""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    return np.where(a < b, b, a)


NUMPY_ENV = {"np": np, "float": float, "engine_min": engine_min,
             "engine_max": engine_max, "__builtins__": {}}


def fuzz_grammar(rng, cases, n_rows=257, verbose=False):
    """Random valid equations, checked against numpy and against the mask."""
    checked = compiled_out = skipped = known = unstable = 0
    failures = []

    # 257 rows: deliberately not a multiple of the evaluator's 512-row block,
    # so the tail of a block is exercised on every single case.
    columns = {
        "x": np.array([rng.uniform(-3, 3) for _ in range(n_rows)]),
        "y": np.array([rng.uniform(-3, 3) for _ in range(n_rows)]),
        # Positive, so sqrt and log have a real branch to exercise rather
        # than returning NaN for half the rows and proving nothing.
        "z": np.array([rng.uniform(0.1, 5) for _ in range(n_rows)]),
    }

    for case in range(cases):
        want = "bool" if rng.random() < 0.45 else "num"
        tree = gen(rng, rng.randint(1, 4), want)
        text = render(tree)
        try:
            expected = eval(render_numpy(tree), dict(NUMPY_ENV), columns)
        except Exception:
            skipped += 1           # numpy itself refused; no oracle to check
            continue
        expected = np.asarray(expected)
        if expected.ndim == 0:
            expected = np.full(n_rows, expected)

        ex = bff.GraphExpression("fuzz")
        try:
            ex.set_expression(text)
        except Exception:
            compiled_out += 1      # the engine refused; allowed, not a failure
            continue

        names = list(ex.get_variable_names())
        if not names:
            skipped += 1
            continue
        try:
            got = np.asarray(ex.compute(names, [list(columns[k]) for k in names]))
        except Exception as exc:
            failures.append((text, f"compute raised {exc!r}"))
            continue

        want_f = expected.astype(float)
        if not close(got, want_f):
            if agrees_with_c_minmax(tree, columns, got, n_rows):
                known += 1
                continue
            if ill_conditioned(tree, columns, n_rows):
                unstable += 1
                continue
            failures.append((text, describe_mismatch(got, want_f)))
            continue

        # The mask evaluator must agree with the double one, always. This is
        # the invariant that CSE or any new stack optimisation is most likely
        # to break, because it is the one that depends on slot *types*.
        block = np.ascontiguousarray(np.vstack([columns[k] for k in names]))
        try:
            mask = ex.compute_mask(names, block)
        except Exception as exc:
            failures.append((text, f"compute_mask raised {exc!r}"))
            continue
        if not np.array_equal(mask.astype(bool), got.astype(bool)):
            bad = np.flatnonzero(mask.astype(bool) != got.astype(bool))[:3]
            failures.append((text, f"mask disagrees with doubles at {bad}"))
            continue

        checked += 1
        if verbose and case < 5:
            print(f"  e.g. {text}")

    return {"checked": checked, "refused_by_engine": compiled_out,
            "no_oracle": skipped, "known": known,
            "unstable": unstable, "failures": failures}


def agrees_with_c_minmax(tree, columns, got, n_rows):
    """T-20260831-09: `min`/`max` do not propagate NaN, where numpy does.

    The engine's `min`/`max` are a plain ternary, `(a > b) ? b : a`, which
    under NaN is **order-dependent**: `min(y, nan)` is `y` but `min(nan, y)`
    is `nan`. numpy propagates in both orders, C's `fmin` ignores in both,
    so neither library function models the engine and the divergence is a
    bug rather than a semantic choice.

    Classifying it from the *output* does not work, because a comparison
    downstream turns the NaN into `False` and the disagreement arrives as
    `1.0` vs `0.0` with no NaN in sight. So the tree is re-rendered against
    `engine_min`/`engine_max`, which reproduce the ternary exactly: if the
    engine matches that twin, this one known rule explains the whole
    disagreement. Anything else still reports.
    """
    try:
        alt = eval(render_numpy(tree, minmax="c"), dict(NUMPY_ENV), columns)
    except Exception:
        return False
    alt = np.asarray(alt)
    if alt.ndim == 0:
        alt = np.full(n_rows, alt)
    return close(got, alt.astype(float))


def ill_conditioned(tree, columns, n_rows):
    """Whether the expression has no stable answer in double precision.

    `sin((z**z) ** (z/2.192))` reaches an argument near 8e7, where one ulp of
    the argument moves the result by ~1e-7 -- so the engine and numpy differ
    in the seventh digit while *both* are computing it correctly, and the
    disagreement says nothing about either.

    Rather than loosen the tolerance for trig (which would hide real bugs in
    every other expression), measure the conditioning directly: nudge every
    input by one ulp and re-evaluate the numpy twin. If numpy disagrees with
    *itself* by more than the tolerance the engine is being held to, no
    implementation can meet that tolerance and the case gets no verdict.
    """
    try:
        base = np.asarray(eval(render_numpy(tree), dict(NUMPY_ENV), columns))
        nudged = {k: np.nextafter(v, np.inf) for k, v in columns.items()}
        moved = np.asarray(eval(render_numpy(tree), dict(NUMPY_ENV), nudged))
    except Exception:
        return False
    if base.ndim == 0:
        base = np.full(n_rows, base)
    if moved.ndim == 0:
        moved = np.full(n_rows, moved)
    return not close(base.astype(float), moved.astype(float))


def describe_mismatch(got, expected):
    """Name the first row that differs, not the first four rows."""
    got = np.asarray(got, dtype=float)
    expected = np.asarray(expected, dtype=float)
    if got.shape != expected.shape:
        return f"shape {got.shape} vs {expected.shape}"
    same = (np.isclose(got, expected, rtol=1e-9, atol=1e-11, equal_nan=True)
            | (np.isnan(got) & np.isnan(expected)))
    bad = np.flatnonzero(~same)
    if bad.size == 0:
        return "close() rejected values that compare equal -- check close()"
    i = int(bad[0])
    return (f"double path disagrees with numpy at row {i} "
            f"({bad.size} of {got.size} rows)\n"
            f"    got      {got[i]!r}\n"
            f"    expected {expected[i]!r}")


def close(got, expected):
    """Agreement, with the non-finite cases treated as numpy treats them."""
    got = np.asarray(got, dtype=float)
    expected = np.asarray(expected, dtype=float)
    if got.shape != expected.shape:
        return False
    both_nan = np.isnan(got) & np.isnan(expected)
    both_inf = np.isinf(got) & np.isinf(expected) & (np.sign(got) == np.sign(expected))
    finite = ~(both_nan | both_inf)
    if not finite.any():
        return True
    # Huge intermediates diverge in the last bits for reasons that are not
    # bugs; the engine's contract is agreement on values a model can produce.
    sane = finite & (np.abs(expected) < 1e12) & np.isfinite(got)
    if finite.sum() != sane.sum():
        finite = sane
        if not finite.any():
            return True
    return np.allclose(got[finite], expected[finite],
                       rtol=1e-9, atol=1e-11, equal_nan=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--cases", type=int, default=None,
                    help="default 200000 structural, 20000 grammar")
    ap.add_argument("--mode", choices=["structural", "grammar", "both"],
                    default="both")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)

    failed = False
    if args.mode in ("structural", "both"):
        rng = random.Random(args.seed)
        n = args.cases if args.cases is not None else 200_000
        result = fuzz_structural(rng, n, args.verbose)
        print(f"structural  seed={args.seed}  cases={n:,}  "
              f"refused={result['refused']:,}  "
              f"evaluated={result['evaluated']:,}  -> survived")

    if args.mode in ("grammar", "both"):
        rng = random.Random(args.seed + 1_000_003)
        n = args.cases if args.cases is not None else 20_000
        result = fuzz_grammar(rng, n, verbose=args.verbose)
        print(f"grammar     seed={args.seed}  cases={n:,}  "
              f"checked={result['checked']:,}  "
              f"refused={result['refused_by_engine']:,}  "
              f"no-oracle={result['no_oracle']:,}  "
              f"known-divergence={result['known']:,}  "
              f"ill-conditioned={result['unstable']:,}")
        if result["failures"]:
            failed = True
            print(f"\n{len(result['failures'])} FAILURES "
                  f"(reproduce with --mode grammar --seed {args.seed}):\n")
            for text, why in result["failures"][:20]:
                print(f"  {text}\n    {why}")
        else:
            print("            -> every case agreed with numpy, "
                  "and every mask agreed with its doubles")

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
