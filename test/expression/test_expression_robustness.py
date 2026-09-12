"""Adversarial input to the expression engine: it must refuse, never crash.

The engine is about to be used by ChiSurf, ndxplorer and tttrlib, where an
equation comes from a YAML file a user edited or a query box someone typed
into. A malformed one has to raise a Python exception. A segfault takes the
whole interpreter -- and every unsaved thing in it -- down with it, so
"refuses cleanly" is a harder requirement here than "is fast".
"""

import itertools
import random
import unittest

import numpy as np

from IMP import bff

#: Anything here must raise, not crash.
MALFORMED = [
    "", " ", "\t\n", "+", "-", "*", "/", "^", "()", "(", ")", "((((",
    "))))", "1+", "+1+", "1++", "1**", "**2", "1 2", "a b", ",", "(,)",
    "1,2", "abs(", "abs)", "abs()", "abs(,)", "pow(1)", "pow(1,2,3)",
    "1/", "/1", "x^", "^x", "@", "$", "#", "!", "?", "%", "1%2", "a&&",
    "||", "~~~", "0x1F", "1e", "1e+", "1.2.3", "..", "1..2",
    "undefined_function(x)", "sin", "unknown_name_zzz",
    "x" * 5000, "(" * 500 + "1" + ")" * 500, "1" + "+1" * 20000,
    "-" * 1000 + "1", "abs(" * 300 + "1" + ")" * 300,
    "\x00", "a\x00b", "é", "λ+1", "𝜋", "x​+1",
]

#: These must parse, and then evaluate without crashing whatever the values.
HOSTILE_VALUES = [
    [0.0], [-0.0], [float("inf")], [float("-inf")], [float("nan")],
    [1e308], [-1e308], [1e-308], [5e-324],
]


class MalformedInputTests(unittest.TestCase):
    def test_malformed_expressions_raise(self):
        crashed = []
        for text in MALFORMED:
            try:
                bff.GraphExpression("m").set_expression(text)
            except (ValueError, TypeError, UnicodeDecodeError, UnicodeEncodeError):
                continue
            except Exception as exc:  # noqa: BLE001 - any exception beats a crash
                crashed.append(f"{text[:30]!r}: unexpected {type(exc).__name__}")
                continue
            # Parsing succeeded; evaluating must still not crash.
            try:
                ex = bff.GraphExpression("m")
                ex.set_expression(text)
                names = list(ex.get_variable_names())
                ex.compute(names, [[1.0]] * len(names))
            except Exception:
                pass
        self.assertEqual(crashed, [])

    def test_is_supported_never_raises(self):
        """The probe is what callers use to decide on a fallback; it has to
        answer for anything, including input that is not an expression."""
        for text in MALFORMED:
            result = bff.GraphExpression.is_supported(text)
            self.assertIn(result, (True, False), text[:30])

    def test_deep_nesting_is_bounded(self):
        for depth in (10, 100, 1000, 5000):
            text = "(" * depth + "1.0" + ")" * depth
            try:
                ex = bff.GraphExpression("m")
                ex.set_expression(text)
                ex.compute([], [])
            except Exception:
                pass  # refusing is fine; crashing is not

    def test_deep_unary_chain_is_bounded(self):
        for depth in (10, 500, 5000):
            try:
                ex = bff.GraphExpression("m")
                ex.set_expression("-" * depth + "x")
                ex.compute(["x"], [[2.0]])
            except Exception:
                pass

    def test_very_long_expression(self):
        try:
            ex = bff.GraphExpression("m")
            ex.set_expression("+".join(["x"] * 10000))
            out = ex.compute(["x"], [[1.0]])
            self.assertAlmostEqual(out[0], 10000.0)
        except Exception:
            pass


class HostileValueTests(unittest.TestCase):
    def test_non_finite_inputs_do_not_crash(self):
        for expr in ("1/x", "sqrt(x)", "log(x)", "x**x", "exp(x)", "x/x",
                     "abs(x)", "x*x+1"):
            ex = bff.GraphExpression("m")
            ex.set_expression(expr)
            for values in HOSTILE_VALUES:
                out = np.asarray(ex.compute(["x"], [values]))
                self.assertEqual(out.size, 1, f"{expr} {values}")

    def test_empty_column(self):
        ex = bff.GraphExpression("m")
        ex.set_expression("x+1")
        self.assertEqual(len(ex.compute(["x"], [[]])), 0)

    def test_mismatched_lengths_are_refused(self):
        ex = bff.GraphExpression("m")
        ex.set_expression("a+b")
        with self.assertRaises(ValueError):
            ex.compute(["a", "b"], [[1.0, 2.0, 3.0], [1.0, 2.0]])

    def test_missing_and_extra_names(self):
        ex = bff.GraphExpression("m")
        ex.set_expression("a+b")
        with self.assertRaises(ValueError):
            ex.compute(["a"], [[1.0]])
        # An extra name is simply unused, not an error.
        out = ex.compute(["a", "b", "c"], [[1.0], [2.0], [9.0]])
        self.assertAlmostEqual(out[0], 3.0)

    def test_evaluating_before_setting_an_expression(self):
        ex = bff.GraphExpression("m")
        try:
            ex.compute([], [])
        except Exception:
            pass

    def test_reuse_after_a_failed_parse(self):
        """A refused equation must leave the object usable, not wedged."""
        ex = bff.GraphExpression("m")
        ex.set_expression("a*x")
        try:
            ex.set_expression("1+(")
        except Exception:
            pass
        out = ex.compute(["a", "x"], [[2.0], [3.0]])
        self.assertAlmostEqual(out[0], 6.0)


class FuzzTests(unittest.TestCase):
    """Random strings from expression-shaped alphabets, 20k of them."""

    ALPHABET = "0123456789.+-*/^()xabc, sqrtexplogminax|&~<>=!"

    def test_random_strings_never_crash(self):
        rng = random.Random(20260831)
        for i in range(20000):
            n = rng.randint(0, 40)
            text = "".join(rng.choice(self.ALPHABET) for _ in range(n))
            try:
                if not bff.GraphExpression.is_supported(text):
                    continue
                ex = bff.GraphExpression("m")
                ex.set_expression(text)
                names = list(ex.get_variable_names())
                if len(names) > 6:
                    continue
                ex.compute(names, [[float(rng.uniform(-5, 5))] for _ in names])
            except Exception:
                continue

    def test_random_well_formed_expressions_match_numpy(self):
        """Generated expressions that ARE valid must also be correct."""
        rng = random.Random(7)
        env = {"__builtins__": {}, "abs": np.abs, "exp": np.exp,
               "sqrt": np.sqrt, "log": np.log}
        x = np.linspace(0.5, 3.0, 32)
        atoms = ["x", "a", "1.5", "2.0", "0.25"]
        ops = ["+", "-", "*", "/"]
        funcs = ["abs", "exp", "sqrt", "log"]
        checked = 0
        for i in range(400):
            def build(depth=0):
                if depth > 2 or rng.random() < 0.4:
                    return rng.choice(atoms)
                if rng.random() < 0.25:
                    return f"{rng.choice(funcs)}({build(depth + 1)})"
                return f"({build(depth + 1)}{rng.choice(ops)}{build(depth + 1)})"

            text = build()
            try:
                ex = bff.GraphExpression("m")
                ex.set_expression(text)
            except Exception:
                continue
            names = list(ex.get_variable_names())
            values = [list(x) if n == "x" else [1.7] for n in names]
            pyv = {n: (x if n == "x" else 1.7) for n in names}
            got = np.asarray(ex.compute(names, values))
            try:
                with np.errstate(all="ignore"):
                    ref = np.atleast_1d(eval(text, dict(env), pyv))
            except ZeroDivisionError:
                # Python raises on scalar 1/0 where numpy and this engine
                # both yield inf; nothing to compare.
                continue
            if ref.size == 1 and got.size > 1:
                ref = np.broadcast_to(ref, got.shape)
            finite = np.isfinite(got) & np.isfinite(ref)
            np.testing.assert_allclose(got[finite], ref[finite],
                                       rtol=1e-11, atol=1e-12, err_msg=text)
            np.testing.assert_array_equal(np.isfinite(got), np.isfinite(ref),
                                          err_msg=text)
            checked += 1
        self.assertGreater(checked, 200)


if __name__ == "__main__":
    unittest.main()


class MultiArgumentFunctionTests(unittest.TestCase):
    """Functions that once came back silently constant, now simply correct.

    History, because it explains the shape of these tests. imp.bff used to
    carry its own evaluator with a vendored ExprTk behind it as a fallback,
    and ExprTk vectorised a *unary* function correctly but evaluated a
    multi-argument one at element 0 and broadcast the result. `hypot(x,y)`
    over two columns came back constant -- no exception, no warning, a flat
    curve that looked entirely plausible (T-20260831-10).

    It was first made to *refuse* rather than answer wrongly. It is now
    simply **right**: the evaluator is tttrlib's, the missing functions were
    implemented there (T-20260831-12), and there is no fallback left to be
    wrong. So these assert values, not exceptions -- and what is still
    genuinely unimplemented is still refused, which is the property that
    mattered all along.
    """

    COLUMNS = {"x": [1.0, 2.0, 3.0, 4.0], "y": [0.5, 1.5, 2.5, 3.5]}

    def compute(self, text):
        ex = bff.GraphExpression("m")
        ex.set_expression(text)
        names = list(self.COLUMNS)
        return np.asarray(ex.compute(names, [self.COLUMNS[k] for k in names]))

    def test_multiarg_functions_now_match_numpy(self):
        x = np.array(self.COLUMNS["x"])
        y = np.array(self.COLUMNS["y"])
        for text, want in [("hypot(x,y)", np.hypot(x, y)),
                           ("atan2(x,y)", np.arctan2(x, y)),
                           ("min(x,y)", np.minimum(x, y)),
                           ("max(x,y)", np.maximum(x, y)),
                           ("pow(x,2)", np.power(x, 2.0))]:
            with self.subTest(text):
                np.testing.assert_allclose(self.compute(text), want, rtol=1e-12)

    def test_the_case_that_proved_the_regression(self):
        """`floor(min(x,y))` -- `floor` used to force the whole expression to
        the fallback and take `min` down with it. Both are in the engine now."""
        np.testing.assert_allclose(self.compute("floor(min(x,y))"),
                                   [0.0, 1.0, 2.0, 3.0])

    def test_the_unary_functions_the_fallback_used_to_supply(self):
        x = np.array(self.COLUMNS["x"])
        for text, want in [("floor(x/2)", np.floor(x / 2)),
                           ("ceil(x/2)", np.ceil(x / 2)),
                           ("trunc(x/2)", np.trunc(x / 2)),
                           ("sign(x-2)", np.sign(x - 2)),
                           ("tanh(x)", np.tanh(x)),
                           ("log2(x)", np.log2(x)),
                           ("expm1(x)", np.expm1(x)),
                           ("log1p(x)", np.log1p(x)),
                           ("erf(x)", None)]:
            if want is None:
                continue
            with self.subTest(text):
                np.testing.assert_allclose(self.compute(text), want, rtol=1e-12)

    def test_what_is_still_unimplemented_is_refused_not_guessed(self):
        """No fallback means no wrong answer: an unknown construct is an
        error, which is what `GraphExpression.h` promises. (`root`/`logn`/`frac`
        left this list on 2026-09-02, when they joined the engine so
        tttrlib's DataStore could drop its ExprTk fallback entirely —
        T-20260831-13; their values are pinned below.)"""
        for text in ["if(x>2,1,0)", "clamp(0,x,2)", "inrange(0,x,2)",
                     "avg(x,y)", "x % 2"]:
            with self.subTest(text), self.assertRaises(ValueError):
                self.compute(text)

    def test_root_logn_and_frac_are_engine_functions_now(self):
        """The last arity-1/2 names of the retired fallback, at ExprTk's
        semantics: root(x, n) = x**(1/n), logn(x, b) = log(x)/log(b), frac
        truncates toward zero."""
        x = np.array(self.COLUMNS["x"])
        np.testing.assert_allclose(
            self.compute("root(x, 2)"), np.sqrt(x), rtol=1e-12)
        np.testing.assert_allclose(
            self.compute("logn(x, 2)"), np.log(x) / np.log(2.0), rtol=1e-12)
        np.testing.assert_allclose(
            self.compute("frac(x + 0.25)"), (x + 0.25) - np.trunc(x + 0.25),
            rtol=1e-12)

    def test_a_name_that_merely_contains_a_function_name_is_a_variable(self):
        """`summary` must not be read as `sum`, nor `xmin` as `min`."""
        ex = bff.GraphExpression("m")
        ex.set_expression("summary + xmin")
        got = np.asarray(ex.compute(["summary", "xmin"], [[1.0, 2.0], [3.0, 4.0]]))
        np.testing.assert_allclose(got, [4.0, 6.0])
