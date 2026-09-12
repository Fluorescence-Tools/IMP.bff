"""``GraphExpression.compute_mask`` -- a gate answered in bytes, not doubles.

A selection asks a yes/no question of every row, and the answer to it is
one bit. Returning that as an array of doubles costs eight bytes a row to
carry one bit of information, so the mask path exists to hand back the
form a selection actually wants.

The bar here is pandas/numpy agreement: whatever ``DataFrame.eval`` or a
numpy comparison says a row is, the mask says the same, including where
the two disagree with naive C -- a NaN operand, a negative number used as
a truth value, a comparison result used as a number.
"""

import unittest

import numpy as np

from IMP import bff


def mask_of(expression, columns):
    """Compile ``expression`` and gate ``columns`` (a name -> array dict)."""
    ex = bff.GraphExpression("gate")
    ex.set_expression(expression)
    names = list(columns)
    block = np.ascontiguousarray(
        np.vstack([np.asarray(columns[k], dtype=float) for k in names]))
    return ex.compute_mask(names, block)


class TestMaskMatchesNumpy(unittest.TestCase):
    """The mask is numpy's answer to the same question, row for row."""

    def setUp(self):
        rng = np.random.default_rng(11)
        self.n = 10007  # deliberately not a multiple of the 512 block
        self.cols = {
            "x": rng.uniform(-5, 5, self.n),
            "y": rng.uniform(0, 100, self.n),
            "z": rng.normal(0, 1, self.n),
        }

    def check(self, expression, expected):
        got = mask_of(expression, self.cols)
        self.assertEqual(got.dtype, np.uint8)
        self.assertEqual(len(got), self.n)
        np.testing.assert_array_equal(got.astype(bool), expected)

    def test_a_single_comparison(self):
        self.check("x > 0", self.cols["x"] > 0)

    def test_every_comparison_operator(self):
        x, y = self.cols["x"], self.cols["y"]
        for text, expected in [
            ("x < 0", x < 0), ("x <= 0", x <= 0),
            ("x > 1.5", x > 1.5), ("x >= 1.5", x >= 1.5),
            ("y == 0", y == 0), ("y != 0", y != 0),
        ]:
            with self.subTest(text):
                self.check(text, expected)

    def test_and_or_and_not(self):
        x, y, z = self.cols["x"], self.cols["y"], self.cols["z"]
        self.check("x > 0 and y < 50", (x > 0) & (y < 50))
        self.check("x > 0 or y < 50", (x > 0) | (y < 50))
        self.check("not(x > 0)", ~(x > 0))
        self.check("x > 0 and y < 50 or z > 1",
                   ((x > 0) & (y < 50)) | (z > 1))

    def test_arithmetic_inside_the_gate(self):
        x, y, z = self.cols["x"], self.cols["y"], self.cols["z"]
        self.check("(x*x + z*z) < 4", (x * x + z * z) < 4)
        self.check("sqrt(abs(x)) > 1.5", np.sqrt(np.abs(x)) > 1.5)
        self.check("y/(abs(x)+1) > 20", y / (np.abs(x) + 1) > 20)

    def test_a_gate_that_spans_many_blocks(self):
        """The block loop writes 512 rows at a time; the seams must line up."""
        n = 5 * 512 + 137
        x = np.arange(n, dtype=float)
        got = mask_of("x > 1000 and x < 2000", {"x": x})
        np.testing.assert_array_equal(got.astype(bool), (x > 1000) & (x < 2000))


class TestMaskTruthiness(unittest.TestCase):
    """Where the answer is a number rather than a comparison.

    numpy's rule is ``arr.astype(bool)``: anything that is not zero is
    true. That keeps a negative operand true and a NaN true, which a naive
    ``> 0.5`` threshold gets wrong in both directions.
    """

    def test_a_numeric_column_is_true_where_nonzero(self):
        x = np.array([0.0, 1.0, -1.0, 0.5, -0.5, 1e-300, np.inf, np.nan])
        got = mask_of("x", {"x": x})
        np.testing.assert_array_equal(got.astype(bool), x.astype(bool))

    def test_a_negative_number_is_true(self):
        got = mask_of("x", {"x": np.full(600, -3.0)})
        self.assertTrue(got.all())

    def test_a_nan_is_true(self):
        got = mask_of("x", {"x": np.full(600, np.nan)})
        self.assertTrue(got.all())

    def test_a_numeric_operand_of_and(self):
        """`(x > 2) and y` -- the right operand is a column of doubles."""
        rng = np.random.default_rng(5)
        x = rng.uniform(0, 5, 3000)
        y = np.where(rng.uniform(0, 1, 3000) > 0.5, 0.0, 7.0)
        got = mask_of("x > 2 and y", {"x": x, "y": y})
        np.testing.assert_array_equal(got.astype(bool),
                                      (x > 2) & y.astype(bool))

    def test_a_numeric_operand_of_not(self):
        x = np.array([0.0, 1.0, -1.0, 0.0, 2.0] * 200)
        got = mask_of("not(x)", {"x": x})
        np.testing.assert_array_equal(got.astype(bool), ~x.astype(bool))

    def test_a_constant_gate_is_the_same_for_every_row(self):
        x = np.zeros(1500)
        np.testing.assert_array_equal(mask_of("1 > 0", {"x": x}), 1)
        np.testing.assert_array_equal(mask_of("1 < 0", {"x": x}), 0)


class TestComparisonUsedAsANumber(unittest.TestCase):
    """A mask flowing back into arithmetic, which numpy allows as 0 or 1."""

    def test_a_comparison_times_a_number(self):
        ex = bff.GraphExpression("m")
        ex.set_expression("(x > 2)*3")
        x = np.linspace(0, 5, 1000)
        got = np.asarray(ex.compute(["x"], [list(x)]))
        np.testing.assert_allclose(got, (x > 2) * 3.0)

    def test_a_comparison_added_to_a_column(self):
        ex = bff.GraphExpression("m")
        ex.set_expression("(x > 2) + x")
        x = np.linspace(0, 5, 1000)
        got = np.asarray(ex.compute(["x"], [list(x)]))
        np.testing.assert_allclose(got, (x > 2) + x)

    def test_a_comparison_under_a_function(self):
        ex = bff.GraphExpression("m")
        ex.set_expression("sqrt(x > 2)")
        x = np.linspace(0, 5, 1000)
        got = np.asarray(ex.compute(["x"], [list(x)]))
        np.testing.assert_allclose(got, np.sqrt((x > 2).astype(float)))


class TestMaskRefusals(unittest.TestCase):
    """What the mask path says no to, rather than answering wrongly."""

    def test_an_unnamed_variable(self):
        ex = bff.GraphExpression("gate")
        ex.set_expression("x > 0 and w < 1")
        block = np.ascontiguousarray(np.vstack([np.arange(4, dtype=float)]))
        with self.assertRaises(ValueError):
            ex.compute_mask(["x"], block)

    def test_an_empty_column_gives_an_empty_mask(self):
        got = mask_of("x > 0", {"x": np.zeros(0)})
        self.assertEqual(len(got), 0)


class TestMaskAgreesWithTheDoublePath(unittest.TestCase):
    """The two evaluators are one program; they may not drift apart."""

    def test_mask_equals_the_double_answer_cast_to_bool(self):
        rng = np.random.default_rng(29)
        n = 4321
        cols = {"a": rng.normal(0, 3, n), "b": rng.uniform(0, 10, n)}
        for text in [
            "a > 0",
            "a > 0 and b < 5",
            "not(a*a > 2)",
            "(a + b) > 5 or a < -2",
            "b/(abs(a)+0.1) > 3",
            "a",
            "a*b",
        ]:
            with self.subTest(text):
                ex = bff.GraphExpression("m")
                ex.set_expression(text)
                names = list(cols)
                doubles = np.asarray(
                    ex.compute(names, [list(cols[k]) for k in names]))
                block = np.ascontiguousarray(
                    np.vstack([cols[k] for k in names]))
                mask = ex.compute_mask(names, block)
                np.testing.assert_array_equal(mask.astype(bool),
                                              doubles.astype(bool))


if __name__ == "__main__":
    unittest.main()
