"""bff.GraphExpression against Python's eval, on ChiSurf's own equation catalogue.

The compiler exists to replace ``eval`` in ChiSurf's parse models, so the
bar is not "plausible arithmetic" but "the same number numpy produces" --
for every equation actually shipped, not just hand-picked ones.
"""

import pathlib
import re
import unittest

import numpy as np

from IMP import bff

CHISURF_MODELS = pathlib.Path(
    "/Users/tpeulen/dev/chisurf/chisurf/core/models")
CATALOGUE = CHISURF_MODELS / "parse" / "models.yaml"


def all_catalogue_equations():
    """Every ``equation:`` in every model catalogue ChiSurf ships.

    Parsed as YAML rather than line by line: several equations wrap across
    lines, and a line-based reader silently truncates them into unbalanced
    expressions that then look like compiler gaps.
    """
    import yaml
    found = []

    def walk(node):
        if isinstance(node, dict):
            for key, value in node.items():
                if key == "equation" and isinstance(value, str):
                    found.append(" ".join(value.split()))
                else:
                    walk(value)
        elif isinstance(node, list):
            for item in node:
                walk(item)

    for path in sorted(CHISURF_MODELS.rglob("*.yaml")):
        try:
            walk(yaml.safe_load(path.read_text()))
        except Exception:
            continue
    return found

NUMPY_ENV = {
    "__builtins__": {},
    "abs": np.abs, "exp": np.exp, "sqrt": np.sqrt, "log": np.log,
    "log10": np.log10, "sin": np.sin, "cos": np.cos, "tan": np.tan,
    "pi": np.pi, "e": np.e,
    "pow": np.power, "min": np.minimum, "max": np.maximum,
}


def numpy_eval(expression, names, values):
    env = dict(zip(names, values))
    return np.atleast_1d(eval(expression, dict(NUMPY_ENV), env))


def bind(expression, x, seed=0):
    """Compile, then invent a plausible value for every free name."""
    rng = np.random.default_rng(seed)
    ex = bff.GraphExpression("m")
    ex.set_expression(expression)
    names = list(ex.get_variable_names())
    values, py_values = [], []
    for name in names:
        if name == "x":
            values.append(list(x))
            py_values.append(x)
        else:
            # Positive and away from zero: these are rates, times and
            # amplitudes, and a negative time would put sqrt() into NaN on
            # both sides and prove nothing.
            v = float(rng.uniform(0.4, 2.5))
            values.append([v])
            py_values.append(v)
    return ex, names, values, py_values


class ExpressionArithmeticTests(unittest.TestCase):
    def _check(self, expression, **env):
        ex = bff.GraphExpression("m")
        ex.set_expression(expression)
        names = list(ex.get_variable_names())
        values = [np.atleast_1d(np.asarray(env[n], dtype=float)).tolist()
                  for n in names]
        got = np.asarray(ex.compute(names, values))
        ref = numpy_eval(expression, names, [env[n] for n in names])
        if ref.size == 1 and got.size > 1:
            ref = np.broadcast_to(ref, got.shape)
        np.testing.assert_allclose(got, ref, rtol=1e-12, atol=1e-12,
                                   err_msg=expression)

    def test_precedence(self):
        for e in ("1+2*3", "2+3*4-5", "(1+2)*3", "2*3/4", "1-2-3"):
            self._check(e)

    def test_power_is_right_associative(self):
        """2**3**2 is 512, not 64 -- Python's rule, not left-to-right."""
        self._check("2**3**2")

    def test_unary_minus_binds_looser_than_power(self):
        """-2**2 is -4 in Python, and must be here too."""
        self._check("-2**2")
        self._check("2*-3")
        self._check("-(3+4)")

    def test_functions(self):
        for e in ("abs(-3)", "exp(1.5)", "sqrt(2)", "log(4)", "log10(100)",
                  "sin(1)+cos(1)", "pow(2,10)", "min(3,4)", "max(3,4)"):
            self._check(e)

    def test_broadcasting(self):
        x = np.linspace(0.5, 3.0, 9)
        self._check("a*x+b", x=x, a=2.0, b=1.0)
        self._check("x*x", x=x)

    def test_division_by_zero_is_inf_not_an_error(self):
        ex = bff.GraphExpression("m")
        ex.set_expression("1/x")
        self.assertEqual(ex.compute(["x"], [[0.0]])[0], float("inf"))

    def test_variables_are_reported_in_first_appearance_order(self):
        ex = bff.GraphExpression("m")
        ex.set_expression("c+a*x+b*x**2")
        self.assertEqual(list(ex.get_variable_names()), ["c", "a", "x", "b"])

    def test_constants_are_not_variables(self):
        ex = bff.GraphExpression("m")
        ex.set_expression("pi*x+e")
        self.assertEqual(list(ex.get_variable_names()), ["x"])

    def test_unsupported_input_is_refused(self):
        for bad in ("x @ y", "foo(x)", "1+", "(1+2", "1+2)", "x $ 2"):
            with self.assertRaises(ValueError, msg=bad):
                bff.GraphExpression("m").set_expression(bad)
            self.assertFalse(bff.GraphExpression.is_supported(bad), bad)

    def test_is_supported_accepts_the_real_thing(self):
        self.assertTrue(bff.GraphExpression.is_supported(
            "b+1/abs(N)*(1+x/td)**(-1)/sqrt(1+1/s**2*x/td)"))

    def test_missing_variable_is_refused(self):
        ex = bff.GraphExpression("m")
        ex.set_expression("a*x")
        with self.assertRaises(ValueError):
            ex.compute(["a"], [[1.0]])


class ExpressionAsANodeTests(unittest.TestCase):
    def test_evaluating_writes_the_curve_to_the_output_port(self):
        x = np.linspace(0.0, 1.0, 5)
        ex = bff.GraphExpression("model")
        ex.set_expression("a*x+b")
        ex.add_input_port("a", bff.GraphPort(2.0, name="a"))
        ex.add_input_port("x", bff.GraphPort(list(x), name="x"))
        ex.add_input_port("b", bff.GraphPort(1.0, name="b"))
        ex.add_output_port("model", bff.GraphPort(0.0, False, True))
        ex.evaluate()
        got = np.asarray(ex.get_output_port("model").value)
        np.testing.assert_allclose(got, 2.0 * x + 1.0, rtol=1e-12)

    def test_a_missing_input_port_is_refused(self):
        ex = bff.GraphExpression("model")
        ex.set_expression("a*x")
        ex.add_input_port("a", bff.GraphPort(2.0, name="a"))
        ex.add_output_port("model", bff.GraphPort(0.0, False, True))
        with self.assertRaises(ValueError):
            ex.evaluate()


@unittest.skipUnless(CATALOGUE.exists(), "ChiSurf equation catalogue not found")
class CatalogueParityTests(unittest.TestCase):
    """Every shipped equation must evaluate exactly as numpy does."""

    @classmethod
    def setUpClass(cls):
        # Parsed as YAML, not line by line: several equations wrap across
        # lines, and a line-based reader silently truncates them into
        # unbalanced expressions that then look like compiler gaps.
        import yaml
        doc = yaml.safe_load(CATALOGUE.read_text())
        cls.equations = []

        def walk(node):
            if isinstance(node, dict):
                for key, value in node.items():
                    if key == "equation" and isinstance(value, str):
                        cls.equations.append(" ".join(value.split()))
                    else:
                        walk(value)
            elif isinstance(node, list):
                for item in node:
                    walk(item)

        walk(doc)

    def test_the_catalogue_is_not_empty(self):
        self.assertGreaterEqual(len(self.equations), 8)

    def test_every_equation_matches_numpy(self):
        x = np.linspace(1e-3, 5.0, 32)
        compiled = skipped = 0
        for index, equation in enumerate(self.equations):
            if not bff.GraphExpression.is_supported(equation):
                skipped += 1
                continue
            compiled += 1
            ex, names, values, py_values = bind(equation, x, seed=index)
            got = np.asarray(ex.compute(names, values))
            ref = numpy_eval(equation, names, py_values)
            if ref.size == 1 and got.size > 1:
                ref = np.broadcast_to(ref, got.shape)
            finite = np.isfinite(got) & np.isfinite(ref)
            np.testing.assert_allclose(
                got[finite], ref[finite], rtol=1e-10, atol=1e-12,
                err_msg=f"equation {index}: {equation}")
            np.testing.assert_array_equal(
                np.isfinite(got), np.isfinite(ref),
                err_msg=f"equation {index} finiteness: {equation}")
        print(f"\ncatalogue: {compiled} equations matched numpy, "
              f"{skipped} unsupported")
        self.assertGreaterEqual(compiled, 6)


if __name__ == "__main__":
    unittest.main()


@unittest.skipUnless(CHISURF_MODELS.exists(), "ChiSurf model catalogues not found")
class EveryShippedEquationTests(unittest.TestCase):
    """Every equation ChiSurf ships must compile -- not just the parse models.

    A model whose equation the engine cannot compile silently falls back to
    Python, so a gap here is a performance cliff nobody would notice.
    """

    @classmethod
    def setUpClass(cls):
        cls.equations = all_catalogue_equations()

    def test_the_catalogues_were_found(self):
        self.assertGreater(len(self.equations), 50,
                           "expected the full set of shipped equations")

    def test_every_shipped_equation_compiles(self):
        failures = []
        for equation in self.equations:
            try:
                bff.GraphExpression("m").set_expression(equation)
            except Exception as exc:
                failures.append(f"{equation!r}: {exc}")
        print(f"\nshipped equations: {len(self.equations)} total, "
              f"{len(self.equations) - len(failures)} compile")
        self.assertEqual(failures, [], "equations that do not compile:\n" +
                         "\n".join(failures))

    def test_every_shipped_equation_matches_numpy(self):
        x = np.linspace(1e-3, 5.0, 16)
        checked = 0
        for index, equation in enumerate(self.equations):
            ex, names, values, py_values = bind(equation, x, seed=index)
            got = np.asarray(ex.compute(names, values))
            try:
                ref = numpy_eval(equation, names, py_values)
            except Exception:
                continue  # numpy itself cannot evaluate it; nothing to compare
            if ref.size == 1 and got.size > 1:
                ref = np.broadcast_to(ref, got.shape)
            finite = np.isfinite(got) & np.isfinite(ref)
            np.testing.assert_allclose(
                got[finite], ref[finite], rtol=1e-10, atol=1e-12,
                err_msg=f"equation {index}: {equation}")
            checked += 1
        print(f"shipped equations: {checked} verified against numpy")
        self.assertGreater(checked, 50)


class CompiledOncePlanTests(unittest.TestCase):
    """The equation is parsed once, not per evaluation.

    During a fit neither the equation nor the batch shape changes, so a
    sampler must only ever re-bind numbers. This is pinned by a test
    because re-parsing per step is invisible in results and ruinous in
    time -- it is exactly the bug this class was rewritten to fix.
    """

    def test_repeated_evaluation_never_reparses(self):
        ex = bff.GraphExpression("m")
        ex.set_expression("b+a*x+c*x**2")
        x = list(np.linspace(0.0, 1.0, 64))
        names = ["a", "b", "c", "x"]
        for step in range(500):
            ex.compute(names, [[1.0 + step * 1e-6], [0.5], [2.0], x])
        self.assertEqual(ex.get_number_of_compilations(), 1)
        self.assertTrue(ex.has_compiled_plan())

    def test_a_new_equation_recompiles(self):
        """A new equation is a new program, compiled once again."""
        ex = bff.GraphExpression("m")
        ex.set_expression("a*x")
        ex.compute(["a", "x"], [[2.0], [1.0, 2.0]])
        self.assertEqual(ex.get_number_of_compilations(), 1)
        ex.set_expression("a+x")
        self.assertEqual(ex.get_number_of_compilations(), 1)
        np.testing.assert_allclose(
            np.asarray(ex.compute(["a", "x"], [[2.0], [1.0, 2.0]])), [3.0, 4.0])

    def test_changing_the_batch_length_never_recompiles(self):
        """The block evaluator carries no per-length state, so a curve of a
        different length costs nothing extra."""
        ex = bff.GraphExpression("m")
        ex.set_expression("a*x")
        ex.compute(["a", "x"], [[2.0], [1.0, 2.0, 3.0]])
        for _ in range(20):
            ex.compute(["a", "x"], [[2.0], [1.0, 2.0, 3.0]])
        ex.compute(["a", "x"], [[2.0], [1.0, 2.0]])
        self.assertEqual(ex.get_number_of_compilations(), 1)

    def test_results_are_correct_across_a_shape_change(self):
        ex = bff.GraphExpression("m")
        ex.set_expression("a*x+1")
        np.testing.assert_allclose(
            np.asarray(ex.compute(["a", "x"], [[2.0], [1.0, 2.0, 3.0]])),
            [3.0, 5.0, 7.0])
        np.testing.assert_allclose(
            np.asarray(ex.compute(["a", "x"], [[3.0], [1.0]])), [4.0])
        np.testing.assert_allclose(
            np.asarray(ex.compute(["a", "x"], [[2.0], [1.0, 2.0, 3.0]])),
            [3.0, 5.0, 7.0])


class BatchEvaluationTests(unittest.TestCase):
    """Whole columns in one call -- what ndxplorer's per-row eval needs."""

    def test_a_column_is_one_call(self):
        ex = bff.GraphExpression("m")
        ex.set_expression("(g-b)/(r-b)")
        n = 5000
        rng = np.random.default_rng(3)
        g, r, b = rng.uniform(1, 9, n), rng.uniform(10, 20, n), rng.uniform(0, 1, n)
        got = np.asarray(ex.compute(["g", "r", "b"],
                                    [list(g), list(r), list(b)]))
        np.testing.assert_allclose(got, (g - b) / (r - b), rtol=1e-12)
        self.assertEqual(ex.get_number_of_compilations(), 1)

    def test_compute_batch_scores_one_row_at_a_time(self):
        ex = bff.GraphExpression("m")
        ex.set_expression("a*a+b")
        rows = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
        got = np.asarray(ex.compute_batch(["a", "b"], rows))
        np.testing.assert_allclose(got, [3.0, 13.0, 31.0])

    def test_compute_batch_refuses_a_ragged_row(self):
        ex = bff.GraphExpression("m")
        ex.set_expression("a+b")
        with self.assertRaises(ValueError):
            ex.compute_batch(["a", "b"], [[1.0, 2.0], [3.0]])


class SharedSubexpressionTests(unittest.TestCase):
    """A repeated subtree is computed once, and that must change nothing.

    Common subexpression elimination hands a finished subtree to a cache
    slot and reads it back at every later occurrence, which means a slot's
    *type* now travels: a cached comparison is a byte mask, a cached
    constant subtree is a folded scalar, and a load has to restore the flags
    a push would have set. The three shapes below are exactly the ones that
    go wrong when it does not.
    """

    def setUp(self):
        rng = np.random.default_rng(17)
        self.n = 2049  # not a multiple of the 512-row block
        self.x = rng.normal(0.0, 1.0, self.n)
        self.z = rng.normal(0.0, 1.0, self.n)
        self.block = np.ascontiguousarray(np.vstack([self.x, self.z]))

    def _got(self, expression):
        ex = bff.GraphExpression("m")
        ex.set_expression(expression)
        return np.asarray(ex.compute_columns(["x", "z"], self.block))

    def _check(self, expression, expected):
        np.testing.assert_allclose(self._got(expression), expected,
                                   rtol=1e-12, err_msg=expression)

    def test_a_shared_arithmetic_subtree(self):
        x = self.x
        self._check("exp(-x/1.5)*exp(-x/1.5)+exp(-x/1.5)",
                    np.exp(-x / 1.5) ** 2 + np.exp(-x / 1.5))

    def test_a_shared_subtree_of_two_columns(self):
        x, z = self.x, self.z
        self._check("sqrt(abs(x*z))+log(sqrt(abs(x*z))+2)",
                    np.sqrt(np.abs(x * z)) + np.log(np.sqrt(np.abs(x * z)) + 2))

    def test_a_shared_comparison_stays_a_mask(self):
        """The cached slot holds bytes, and `and` must still read bytes."""
        x = self.x
        self._check("(exp(x) > 1.5) and (exp(x) > 1.5)",
                    (np.exp(x) > 1.5).astype(float))

    def test_a_shared_comparison_used_as_a_number(self):
        """...and arithmetic on the same cached slot must widen it back."""
        x = self.x
        self._check("(exp(x) > 1.5)*3+(exp(x) > 1.5)",
                    (np.exp(x) > 1.5) * 3.0 + (np.exp(x) > 1.5))

    def test_a_shared_subtree_under_a_negation(self):
        x = self.x
        self._check("not (exp(x) > 1.5) or (exp(x) > 1.5)",
                    np.ones(self.n))

    def test_sharing_survives_a_change_of_batch_length(self):
        """The cache is per block, so a shorter curve must not see stale rows."""
        ex = bff.GraphExpression("m")
        ex.set_expression("exp(-x/1.5)*exp(-x/1.5)")
        for n in (7, 512, 513, 4096, 3):
            x = np.linspace(0.01, 4.0, n)
            got = np.asarray(ex.compute_columns(["x"],
                                                np.ascontiguousarray(x.reshape(1, n))))
            np.testing.assert_allclose(got, np.exp(-x / 1.5) ** 2, rtol=1e-12)


class ConstantExponentTests(unittest.TestCase):
    """``x**(-1)`` is a reciprocal, however the exponent was written.

    The exponent of an FCS model is a *parenthesised negative* constant,
    which parses as OP_CONST followed by a negation rather than as a bare
    constant. The compiler folded only the bare form, so the whole
    catalogue's `**(-1)` and `**(-0.5)` fell through to a ``std::pow()``
    call per element -- correct, and an order of magnitude slower than the
    reciprocal numpy rewrites the same expression into.
    """

    def setUp(self):
        rng = np.random.default_rng(23)
        self.n = 1024
        self.x = rng.uniform(0.5, 4.0, self.n)
        self.block = np.ascontiguousarray(self.x.reshape(1, self.n))

    def _check(self, expression, expected):
        ex = bff.GraphExpression("m")
        ex.set_expression(expression)
        got = np.asarray(ex.compute_columns(["x"], self.block))
        np.testing.assert_allclose(got, expected, rtol=1e-12,
                                   err_msg=expression)

    def test_every_spelling_of_a_negative_exponent(self):
        x = self.x
        for text in ("x**(-1)", "x**-1", "x**(-(1))"):
            with self.subTest(text):
                self._check(text, 1.0 / x)

    def test_fractional_and_square_exponents(self):
        x = self.x
        self._check("x**(-0.5)", x ** -0.5)
        self._check("x**(-2)", x ** -2.0)
        self._check("(1+x/1.2)**(-1)", (1 + x / 1.2) ** -1.0)

    def test_an_exponent_that_is_a_scalar_parameter(self):
        """A fit supplies the exponent as a scalar; it takes the same path."""
        ex = bff.GraphExpression("m")
        ex.set_expression("x**p")
        for p in (-1.0, -0.5, 0.5, 2.0, 1.7):
            with self.subTest(p=p):
                got = np.asarray(ex.compute(["x", "p"],
                                            [list(self.x), [p]]))
                np.testing.assert_allclose(got, self.x ** p, rtol=1e-12)

    def test_a_wholly_constant_expression_still_broadcasts(self):
        """Constant folding must not turn a curve into a single number."""
        ex = bff.GraphExpression("m")
        ex.set_expression("2**(-1)+sqrt(9)*x/x")
        got = np.asarray(ex.compute_columns(["x"], self.block))
        np.testing.assert_allclose(got, np.full(self.n, 3.5), rtol=1e-12)


class BoundParameterTests(unittest.TestCase):
    """`bind_parameters` / `compute_curve_bound`: no strings per iteration.

    A fit re-evaluates one equation thousands of times, changing only values.
    `compute_curve()` takes the names on every call, so SWIG rebuilds a
    ``std::vector<std::string>`` per step -- about 0.1 us a name, which is
    26-30% of a 512-point evaluation for a typical model. Binding resolves the
    mapping once; after that only numbers cross the boundary.

    Being an optimisation, the standard it is held to is that it changes
    nothing: same curve, to the last bit.
    """

    def test_the_bound_result_is_bit_identical_to_the_named_one(self):
        for text, names in [("b+a1*exp(-x/t1)", ["b", "a1", "t1"]),
                            ("b+1/N*(1+x/td)**(-1)/sqrt(1+1/s**2*x/td)",
                             ["b", "N", "td", "s"]),
                            ("(a-b)/(c-b)", ["a", "b", "c"])]:
            ex = bff.GraphExpression("m")
            ex.set_expression(text)
            ex.bind_parameters(names, "x")
            values = np.arange(1.0, len(names) + 1.0)
            x = np.linspace(0.01, 10.0, 337)   # not a multiple of the block
            with self.subTest(text):
                np.testing.assert_allclose(
                    ex.compute_curve_bound(values, x),
                    ex.compute_curve(names, values, "x", x), rtol=1e-15)

    def test_binding_survives_changing_the_values(self):
        """The whole point: rebind never, re-evaluate often."""
        ex = bff.GraphExpression("m")
        ex.set_expression("b+a1*exp(-x/t1)")
        ex.bind_parameters(["b", "a1", "t1"], "x")
        x = np.linspace(0.1, 5.0, 64)
        for a1 in (1.0, 2.0, 3.0):
            got = ex.compute_curve_bound(np.array([0.5, a1, 1.5]), x)
            np.testing.assert_allclose(got, 0.5 + a1 * np.exp(-x / 1.5),
                                       rtol=1e-12)

    def test_binding_order_is_the_callers_not_the_engines(self):
        """`variables()` is in first-appearance order; a caller's list is not.

        Binding by position against the engine's order would evaluate the
        right equation against the wrong values, silently.
        """
        ex = bff.GraphExpression("m")
        ex.set_expression("(g-b)/(r-b)")           # engine order: g, b, r
        ex.bind_parameters(["b", "g", "r"], "x")   # caller order: b, g, r
        x = np.linspace(0.1, 1.0, 8)
        got = ex.compute_curve_bound(np.array([1.0, 5.0, 9.0]), x)
        np.testing.assert_allclose(got, (5.0 - 1.0) / (9.0 - 1.0))

    def test_a_wrong_parameter_count_is_refused(self):
        ex = bff.GraphExpression("m")
        ex.set_expression("b+a1*exp(-x/t1)")
        ex.bind_parameters(["b", "a1", "t1"], "x")
        with self.assertRaises(ValueError):
            ex.compute_curve_bound(np.array([1.0, 2.0]), np.linspace(1, 2, 4))

    def test_evaluating_before_binding_is_refused(self):
        ex = bff.GraphExpression("m")
        ex.set_expression("b+a1*exp(-x/t1)")
        self.assertFalse(ex.has_parameter_binding())
        with self.assertRaises(ValueError):
            ex.compute_curve_bound(np.array([1.0, 2.0, 3.0]),
                                   np.linspace(1, 2, 4))

    def test_a_new_equation_drops_the_binding(self):
        """Keeping it would read the new equation from the old slots."""
        ex = bff.GraphExpression("m")
        ex.set_expression("b+a1*exp(-x/t1)")
        ex.bind_parameters(["b", "a1", "t1"], "x")
        self.assertTrue(ex.has_parameter_binding())
        ex.set_expression("p+q*x")
        self.assertFalse(ex.has_parameter_binding())

    def test_binding_a_name_the_equation_does_not_use_is_refused(self):
        ex = bff.GraphExpression("m")
        ex.set_expression("b+a1*exp(-x/t1)")
        with self.assertRaises(ValueError):
            ex.bind_parameters(["b", "a1"], "x")   # t1 unaccounted for
