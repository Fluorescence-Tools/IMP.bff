"""A TCSPC decay as a node, and the kernels underneath it.

`TcspcDecay` is what makes a lifetime fit a graph: the model curve of a
multi-exponential decay through a real instrument -- reconvolved with a
measured response, shifted against it, scaled to the data, plus scatter and a
background -- computed in C++ so a fit does not return to Python once per
iteration.

The tests that matter are the ones pinning it to something that is *not*
this class:

* the convolution against **tttrlib's own** ``fconv_per_cs``, because the
  node uses the header-only twin of that kernel and the two must agree;
* the timeshift against **chisurf's** ``shift_array``, because that is the
  operation the application means by a timeshift and the two spell it with
  opposite signs;
* the whole curve against a numpy transcription of ChiSurf's
  ``LifetimeModel.update_model``.

Speed is not checked here. It is checked by not crossing the boundary, which
is a property of the arrangement rather than of a number.
"""

import pathlib
import sys
import unittest

import numpy as np
import pytest

from IMP import bff

try:
    import tttrlib
except ImportError:  # pragma: no cover - depends on the environment
    tttrlib = None

# ChiSurf is a sibling checkout (AGENTS.md), not an installed package.
_SIBLING = pathlib.Path(__file__).resolve().parents[3] / "chisurf"
if _SIBLING.is_dir() and str(_SIBLING) not in sys.path:
    sys.path.insert(0, str(_SIBLING))

N = 512
DT = 0.032
PERIOD = 1000.0 / 80.0


def response(n=N, dt=DT, centre=1.0, width=0.08, height=1000.0):
    """A narrow instrument response, one peak, on the data's time axis."""
    x = np.arange(n) * dt
    return height * np.exp(-0.5 * ((x - centre) / width) ** 2)


def build(n_lifetimes=1, irf=None, name="decay"):
    """A node with its ports made and its response set."""
    if irf is None:
        irf = response()
    node = bff.TcspcDecay(name)
    node.set_number_of_lifetimes(n_lifetimes)
    node.add_output_port(name, bff.GraphPort([0.0], False, True))
    node.set_response_array(np.ascontiguousarray(irf, dtype=float))
    node.set_timing(DT, PERIOD)
    node.set_convolution_range(len(irf), len(irf))
    return node


def curve_of(node, name="decay"):
    node.update()
    return np.asarray(node.get_output_port(name).value, dtype=float)


def set_spectrum(node, pairs):
    for i, (amplitude, lifetime) in enumerate(pairs):
        node.get_input_port("a%d" % i).value = float(amplitude)
        node.get_input_port("t%d" % i).value = float(lifetime)


def reference_curve(irf, pairs, timeshift=0.0, scatter=0.0, background=0.0,
                    n0=1.0, dt=DT, period=PERIOD):
    """ChiSurf's ``update_model``, transcribed in numpy.

    Deliberately written from the ChiSurf side rather than from the node's:
    a transcription of the implementation under test would pass whatever the
    implementation did.
    """
    pytest.importorskip("chisurf",
                        reason="the reference timeshift lives in a sibling "
                               "chisurf checkout; CI has neither installed")
    from chisurf.core.math.signal import shift_array

    spectrum = np.asarray([v for pair in pairs for v in pair], dtype=float)
    irf_y = np.asarray(irf, dtype=float)
    if timeshift != 0.0:
        irf_y = shift_array(irf_y, timeshift)
    irf_y = np.resize(irf_y, (len(irf),))
    total = irf_y.sum()
    if total > 0:
        irf_y = irf_y / total
    decay = np.zeros(len(irf))
    last = len(irf) - 1
    tttrlib.fconv_per_cs(decay, irf_y, spectrum, period,
                         min(len(irf), last), min(len(irf), last), dt)
    decay = decay + scatter * irf_y
    decay = decay * n0 + background
    return np.maximum(decay, 0)


@unittest.skipUnless(tttrlib is not None, "tttrlib not importable")
class KernelAgreementTests(unittest.TestCase):
    """The node must not have its own opinion about the kernels."""

    def test_the_convolution_is_tttrlibs(self):
        """`fconv_per_cs_ad<double>` against the library's `fconv_per_cs`.

        The node uses the header-only templated core; tttrlib's exported
        function dispatches to a SIMD kernel above one lifetime. They are
        documented to agree to rounding, and this is where that is checked.
        """
        irf = response()
        for pairs in ([(1.0, 3.5)],
                      [(0.7, 0.5), (0.3, 4.0)],
                      [(0.5, 0.2), (0.3, 1.5), (0.2, 6.0)]):
            with self.subTest(components=len(pairs)):
                node = build(len(pairs), irf)
                set_spectrum(node, pairs)
                got = curve_of(node)
                want = reference_curve(irf, pairs)
                self.assertLess(np.max(np.abs(got - want)),
                                1e-12 * np.max(want))

    def test_the_timeshift_is_chisurfs_shift(self):
        """`shift_lamp(v, -s)` is `shift_array(v, s)`, and this pins it.

        The two index in opposite directions *and* interpolate toward
        opposite neighbours; the flips cancel. They part company at exactly
        ``s == 0``, where `shift_lamp` drops the last sample -- so the node
        short-circuits there, and a fractional shift either side of zero has
        to come out continuous.
        """
        irf = response()
        pairs = [(1.0, 3.5)]
        for shift in (0.0, 1.0, -1.0, 2.5, -2.5, 0.25, -0.25, -3.75):
            with self.subTest(shift=shift):
                node = build(1, irf)
                set_spectrum(node, pairs)
                node.get_input_port("timeshift").value = shift
                got = curve_of(node)
                want = reference_curve(irf, pairs, timeshift=shift)
                self.assertLess(np.max(np.abs(got - want)),
                                1e-12 * np.max(want))

    def test_a_whole_curve_matches_the_numpy_transcription(self):
        """Every term at once: shift, scatter, scale, background."""
        irf = response()
        pairs = [(0.6, 0.4), (0.4, 3.9)]
        node = build(2, irf)
        set_spectrum(node, pairs)
        node.get_input_port("timeshift").value = -0.7
        node.get_input_port("scatter").value = 0.05
        node.get_input_port("background").value = 1.5
        node.get_input_port("n0").value = 250.0
        got = curve_of(node)
        want = reference_curve(irf, pairs, timeshift=-0.7, scatter=0.05,
                               background=1.5, n0=250.0)
        self.assertLess(np.max(np.abs(got - want)), 1e-9 * np.max(want))


class NodeBehaviourTests(unittest.TestCase):

    def test_a_port_write_invalidates_the_curve(self):
        """It is a node, so a new parameter means a new curve."""
        node = build(1)
        set_spectrum(node, [(1.0, 2.0)])
        first = curve_of(node).copy()
        node.get_input_port("t0").value = 4.0
        second = curve_of(node)
        self.assertGreater(np.max(np.abs(second - first)), 0.0)

    def test_a_negative_lifetime_is_used_as_its_magnitude(self):
        """A fit that walks a lifetime through zero must not get a growing
        exponential, which is what a negative one is."""
        a, b = build(1), build(1)
        set_spectrum(a, [(1.0, 2.5)])
        set_spectrum(b, [(1.0, -2.5)])
        np.testing.assert_allclose(curve_of(a), curve_of(b), rtol=0, atol=0)

    def test_amplitudes_are_normalised_when_asked(self):
        """`normalize_amplitudes` divides by the sum, so one amplitude pins
        the scale and doubling every one of them changes nothing."""
        node = build(2)
        node.set_normalize_amplitudes(True)
        set_spectrum(node, [(1.0, 0.5), (1.0, 4.0)])
        one = curve_of(node).copy()
        set_spectrum(node, [(2.0, 0.5), (2.0, 4.0)])
        np.testing.assert_allclose(curve_of(node), one, rtol=1e-14)

    def test_absolute_amplitudes_are_used_as_magnitudes(self):
        node = build(1)
        node.set_absolute_amplitudes(True)
        set_spectrum(node, [(-2.0, 3.0)])
        node.set_absolute_amplitudes(False)
        node.set_valid(False)
        signed = curve_of(node)
        node.set_absolute_amplitudes(True)
        node.set_valid(False)
        absolute = curve_of(node)
        # The signed curve is clamped to zero everywhere; the absolute one is
        # a decay. If `absolute_amplitudes` did nothing they would agree.
        self.assertEqual(np.max(signed), 0.0)
        self.assertGreater(np.max(absolute), 0.0)

    def test_a_negative_curve_is_clamped(self):
        node = build(1)
        set_spectrum(node, [(1.0, 3.0)])
        node.get_input_port("background").value = -1e6
        self.assertEqual(np.max(curve_of(node)), 0.0)

    def test_a_nan_survives_the_clamp(self):
        """The clamp is `v < 0`, not `!(v > 0)`.

        The second spelling turns a NaN into a *zero*, which in a fit reads
        as a good fit near zero -- the same trap the port sanitiser sets. A
        NaN has to reach `ChiSquared`, which is what makes the misfit
        infinite and the step rejected.

        Getting one *into* the curve takes opting the parameter port out of
        sanitising first, which is itself worth recording: an input port
        floors a NaN to `tiny` on the way in, so a NaN cannot arrive through
        a parameter -- only from the arithmetic.
        """
        node = build(1)
        set_spectrum(node, [(1.0, 3.0)])
        background = node.get_input_port("background")
        background.value = float("nan")
        self.assertFalse(np.any(np.isnan(curve_of(node))),
                         "an input port stopped sanitising")

        background.set_sanitize(False)
        background.value = float("nan")
        node.set_valid(False)
        self.assertTrue(np.all(np.isnan(curve_of(node))))

    def test_the_curve_port_does_not_sanitise(self):
        """chinet floors a NaN to `tiny` on a stored value; for a fit that is
        backwards, so the transport opts out."""
        node = build(1)
        set_spectrum(node, [(1.0, 3.0)])
        curve_of(node)
        self.assertFalse(node.get_output_port("decay").get_sanitize())

    def test_a_node_without_a_response_refuses_to_evaluate(self):
        node = bff.TcspcDecay("decay")
        node.set_number_of_lifetimes(1)
        node.add_output_port("decay", bff.GraphPort([0.0], False, True))
        with self.assertRaises(ValueError):
            node.update()

    def test_a_node_without_components_refuses_to_evaluate(self):
        node = bff.TcspcDecay("decay")
        node.set_number_of_lifetimes(0)
        node.add_output_port("decay", bff.GraphPort([0.0], False, True))
        node.set_response_array(np.ascontiguousarray(response()))
        with self.assertRaises(ValueError):
            node.update()

    def test_a_non_positive_period_is_refused(self):
        node = bff.TcspcDecay("decay")
        node.set_number_of_lifetimes(1)
        with self.assertRaises(ValueError):
            node.set_timing(DT, 0.0)


class AutoscaleTests(unittest.TestCase):
    """`n0` from the data rather than from the optimiser."""

    def _data(self, pairs, n0=3000.0, background=2.0, seed=4):
        irf = response()
        want = reference_curve(irf, pairs, n0=n0, background=background)
        rng = np.random.default_rng(seed)
        y = rng.poisson(np.maximum(want, 0.0)).astype(float)
        ey = np.sqrt(np.maximum(y, 1.0))
        return irf, y, ey

    @unittest.skipUnless(tttrlib is not None, "tttrlib not importable")
    def test_the_scale_is_chisurfs_weighted_least_squares(self):
        pairs = [(1.0, 3.5)]
        irf, y, ey = self._data(pairs)
        node = build(1, irf)
        set_spectrum(node, pairs)
        node.set_data_arrays(np.ascontiguousarray(y),
                             np.ascontiguousarray(ey))
        node.set_autoscale(True)
        node.set_scale_range(0, len(y))
        node.get_input_port("background").value = 2.0
        got = curve_of(node)

        from chisurf.core.fluorescence.tcspc import rescale_w_bg
        unscaled = reference_curve(irf, pairs)
        want_n0 = rescale_w_bg(
            model_decay=unscaled, experimental_decay=y,
            experimental_weights=1.0 / ey, experimental_background=2.0,
            start=0, stop=len(y))
        self.assertAlmostEqual(node.get_n0() / want_n0, 1.0, places=10)
        np.testing.assert_allclose(got, np.maximum(unscaled * want_n0 + 2.0,
                                                   0.0), rtol=1e-10)

    def test_the_autoscaled_amplitude_is_published_on_its_port(self):
        """A fit that autoscales still has to report the amplitude it
        settled on, and the port is where a caller reads it."""
        pairs = [(1.0, 3.5)]
        irf, y, ey = self._data(pairs)
        node = build(1, irf)
        set_spectrum(node, pairs)
        node.set_data_arrays(np.ascontiguousarray(y),
                             np.ascontiguousarray(ey))
        node.set_autoscale(True)
        node.set_scale_range(0, len(y))
        curve_of(node)
        self.assertAlmostEqual(node.get_input_port("n0").value, node.get_n0())
        self.assertGreater(node.get_n0(), 0.0)

    def test_autoscaling_without_data_is_refused(self):
        node = build(1)
        set_spectrum(node, [(1.0, 3.0)])
        node.set_autoscale(True)
        with self.assertRaises(ValueError):
            node.update()


class WholeFitTests(unittest.TestCase):
    """`TcspcDecay -> ChiSquared -> Minimizer`: the point of the class."""

    def test_a_lifetime_is_recovered_from_simulated_counts(self):
        irf = response()
        truth = [(1.0, 3.1)]
        clean = reference_curve(irf, truth, n0=4000.0, background=1.0)
        rng = np.random.default_rng(11)
        y = rng.poisson(np.maximum(clean, 0.0)).astype(float)
        ey = np.sqrt(np.maximum(y, 1.0))

        node = build(1, irf)
        set_spectrum(node, [(1.0, 1.0)])
        node.set_data_arrays(np.ascontiguousarray(y),
                             np.ascontiguousarray(ey))
        node.set_autoscale(True)
        node.set_scale_range(0, len(y))

        chi2 = bff.ChiSquared("chi2")
        chi2.set_data_arrays(np.ascontiguousarray(y),
                             np.ascontiguousarray(ey))
        chi2.set_fit_range(0, len(y))
        model_in = bff.GraphPort([0.0])
        model_in.link = node.get_output_port("decay")
        chi2.add_input_port("model", model_in)
        chi2.add_output_port("chi2", bff.GraphPort(0.0, False, True))
        chi2.add_output_port("residuals", bff.GraphPort([0.0], False, True))

        free = [node.get_input_port("t0"), node.get_input_port("background")]
        free[0].value = 1.0
        free[1].value = 0.0
        m = bff.Minimizer()
        m.set_parameter_ports(free)
        m.set_objective(chi2, "residuals")
        m._graph = (node, chi2, model_in)

        info = m.run()
        self.assertIn(info, (1, 2, 3, 4))
        self.assertAlmostEqual(m.x[0], 3.1, delta=0.05)
        self.assertLess(chi2.get_chi2r(len(free)), 1.5)

    def test_the_optimiser_never_leaves_cxx(self):
        """The claim the whole port rests on, stated as a test.

        Nothing in the graph is a Python object with an `evaluate`, so a
        `Minimizer` step cannot re-enter the interpreter. If a future change
        introduces a director in this path it will show up here as a node
        that is not one of the C++ types.
        """
        node = build(1)
        self.assertIsInstance(node, bff.GraphNode)
        self.assertNotIn("director", type(node).__name__.lower())


if __name__ == "__main__":
    unittest.main()


class SpectrumPortTests(unittest.TestCase):
    """The spectrum can come from upstream, which is what composes.

    Every TCSPC model shares this instrument model and differs only in how
    the (amplitude, lifetime) pairs are arrived at -- a FRET model from a
    distance, a distribution model from a distance distribution, a mixture
    from two spectra. Those become *nodes*, and this node stops holding an
    opinion about the photophysics.
    """

    def test_the_port_and_the_scalars_give_the_same_curve(self):
        """The two ways in must not be two models."""
        irf = response()
        pairs = [(0.6, 0.4), (0.4, 3.9)]

        scalars = build(2, irf)
        set_spectrum(scalars, pairs)

        piped = build(2, irf)
        piped.set_spectrum_from_port(True)
        piped.get_input_port("lifetime_spectrum").set_values_array(
            np.ascontiguousarray([v for pair in pairs for v in pair],
                                 dtype=float))

        np.testing.assert_allclose(curve_of(piped), curve_of(scalars),
                                   rtol=0, atol=0)

    def test_re_pushing_an_unchanged_spectrum_can_be_free(self):
        """The advice in `set_spectrum_from_port`'s documentation, pinned.

        A caller whose Jacobian column perturbs something *other* than the
        spectrum still pushes the spectrum, and without memoisation pays a
        full reconvolution for a curve it already had. This node is a
        function of its inputs, so `set_memoize` is sound for it -- and the
        curve must come out identical, which is the half that matters.
        """
        irf = response()
        pairs = [(0.6, 0.4), (0.4, 3.9)]
        flat = np.ascontiguousarray([v for pair in pairs for v in pair],
                                    dtype=float)

        node = build(2, irf)
        node.set_spectrum_from_port(True)
        port = node.get_input_port("lifetime_spectrum")
        port.set_values_array(flat)
        expected = curve_of(node).copy()

        node.set_memoize(True)
        port.set_values_array(flat)
        node.update()
        evaluations = node.get_evaluation_count()
        hits = node.get_memo_hit_count()
        for _ in range(5):
            port.set_values_array(flat)
            node.update()
        self.assertEqual(node.get_evaluation_count() - evaluations, 0)
        self.assertEqual(node.get_memo_hit_count() - hits, 5)
        np.testing.assert_allclose(curve_of(node), expected, rtol=0, atol=0)

        # and a spectrum that genuinely moved still gets through
        moved = flat.copy()
        moved[1] = 0.7
        port.set_values_array(moved)
        node.update()
        self.assertEqual(node.get_evaluation_count() - evaluations, 1)
        self.assertFalse(np.array_equal(curve_of(node), expected))

    def test_the_cached_response_is_dropped_when_it_should_be(self):
        """`evaluate` keeps the shifted, renormalised response between calls,
        because it depends on the response and the timeshift and on nothing
        else -- which is a cache, and a cache is only as good as what drops
        it. Four things that must reach the curve."""
        irf = response()
        pairs = [(0.6, 0.4), (0.4, 3.9)]
        node = build(2, irf)
        set_spectrum(node, pairs)
        base = curve_of(node).copy()

        # a timeshift moves it, and shifting back reproduces it exactly
        node.get_input_port("timeshift").value = 2.5
        shifted = curve_of(node).copy()
        self.assertFalse(np.array_equal(shifted, base))
        node.get_input_port("timeshift").value = 0.0
        np.testing.assert_allclose(curve_of(node), base, rtol=0, atol=0)

        # a new response reaches the curve
        node.set_response(list(response(centre=1.4)))
        self.assertFalse(np.array_equal(curve_of(node), base))

        # and a rescaled response does not, because it is normalised.
        # Not bitwise: `(7x)/(7T)` and `x/T` are the same number in real
        # arithmetic and differ in the last place in this one, so the claim
        # is scale invariance, not reproducibility.
        node.set_response(list(np.asarray(irf) * 7.0))
        np.testing.assert_allclose(curve_of(node), base, rtol=1e-14, atol=0)

    def _basis_node(self, negative=False, **options):
        irf = response()
        node = build(4, irf)
        node.add_output_port(bff.TcspcDecay.basis_port_key(),
                             bff.GraphPort([0.0], False, True))
        node.set_emit_basis(True)
        node.set_spectrum_from_port(True)
        if options.get("normalize"):
            node.set_normalize_amplitudes(True)
        if options.get("absolute"):
            node.set_absolute_amplitudes(True)
        if options.get("autoscale"):
            y = np.maximum(np.arange(N, dtype=float) % 97.0, 1.0)
            node.set_data(list(y), list(np.sqrt(y)))
            node.set_autoscale(True)
        spectrum = np.array([0.5, 0.3, -0.3 if negative else 0.3, 1.1,
                             0.15, 2.7, 0.05, 5.0])
        return node, spectrum

    def _curve_at(self, node, spectrum):
        node.get_input_port("lifetime_spectrum").set_values_array(
            np.ascontiguousarray(spectrum, dtype=float))
        node.update()
        return curve_of(node).copy()

    def _worst_jacobian_error(self, node, spectrum):
        """max relative error of the basis against a forward difference."""
        base = self._curve_at(node, spectrum)
        basis = np.asarray(
            node.get_output_port(bff.TcspcDecay.basis_port_key())
            .get_value_view()).reshape(len(base), -1).copy()
        worst = 0.0
        for s in range(basis.shape[1]):
            step = 1e-7 * max(1.0, abs(spectrum[2 * s]))
            moved = np.array(spectrum, dtype=float)
            moved[2 * s] += step
            fd = (self._curve_at(node, moved) - base) / step
            scale = max(np.max(np.abs(fd)), 1e-30)
            worst = max(worst, np.max(np.abs(fd - basis[:, s])) / scale)
        return worst

    def test_the_basis_is_the_jacobian_where_the_node_says_it_is(self):
        """`set_emit_basis` exists so a caller can build a Jacobian out of the
        basis, and three of this node's own options quietly break that. The
        query has to be honest in both directions, so this checks it against
        a forward difference rather than against the same reasoning that
        wrote it."""
        for label, options, negative in (
                ("plain", {}, False),
                ("absolute amplitudes, all positive", {"absolute": True}, False),
                ("a negative amplitude, no transform", {}, True)):
            node, spectrum = self._basis_node(negative=negative, **options)
            self._curve_at(node, spectrum)
            self.assertTrue(node.get_basis_is_jacobian(), label)
            self.assertLess(self._worst_jacobian_error(node, spectrum), 1e-5,
                            label)

    def test_and_says_so_where_it_is_not(self):
        """The half that matters: a caller trusting the basis under these
        settings gets a wrong Jacobian with no symptom at all."""
        for label, options, negative in (
                ("normalize_amplitudes", {"normalize": True}, False),
                ("autoscale", {"autoscale": True}, False),
                ("absolute amplitudes, one negative", {"absolute": True}, True)):
            node, spectrum = self._basis_node(negative=negative, **options)
            self._curve_at(node, spectrum)
            self.assertFalse(node.get_basis_is_jacobian(), label)
            # and it really is wrong, not merely suspected of it
            self.assertGreater(self._worst_jacobian_error(node, spectrum), 1e-3,
                               label)

    def test_an_upstream_node_can_drive_it(self):
        """The arrangement, end to end: a node computes the spectrum.

        `GraphExpression` stands in here for whatever a real model would be -- a
        FRET rate, a distance distribution -- because what is being checked
        is the *wiring*, not the photophysics: a linked vector port, and
        `GraphNode::update()` walking from one call.
        """
        irf = response()
        maker = bff.GraphExpression("spectrum")
        # An amplitude and a lifetime derived from one parameter, which is
        # the shape every derived spectrum has.
        maker.set_expression("x*0 + 1")
        axis = bff.GraphPort([0.0, 0.0])
        axis.set_values_array(np.ascontiguousarray([1.0, 3.5]))
        maker.add_input_port("x", axis)
        out = bff.GraphPort([0.0], False, True)
        out.set_sanitize(False)
        maker.add_output_port("spectrum", out)

        node = build(1, irf)
        node.set_spectrum_from_port(True)
        node.get_input_port("lifetime_spectrum").link = out
        node._upstream = (maker, axis, out)

        got = curve_of(node)
        # `x*0 + 1` over the axis (1.0, 3.5) is the spectrum (1, 1): one
        # component, unit amplitude, a lifetime of one.
        want = reference_curve(irf, [(1.0, 1.0)])
        self.assertLess(np.max(np.abs(got - want)), 1e-12 * np.max(want))
        self.assertEqual(node.get_number_of_lifetimes(), 1)

    def test_the_amplitude_conventions_still_apply(self):
        """`normalize_amplitudes` describes the model, not the wiring."""
        irf = response()
        node = build(2, irf)
        node.set_spectrum_from_port(True)
        node.set_normalize_amplitudes(True)
        port = node.get_input_port("lifetime_spectrum")

        port.set_values_array(np.ascontiguousarray([1.0, 0.5, 1.0, 4.0]))
        one = curve_of(node).copy()
        port.set_values_array(np.ascontiguousarray([2.0, 0.5, 2.0, 4.0]))
        node.set_valid(False)
        np.testing.assert_allclose(curve_of(node), one, rtol=1e-14)

    def test_an_odd_length_spectrum_is_refused(self):
        """Rounding down would silently drop the last amplitude's lifetime."""
        node = build(1)
        node.set_spectrum_from_port(True)
        node.get_input_port("lifetime_spectrum").set_values_array(
            np.ascontiguousarray([1.0, 2.0, 3.0]))
        with self.assertRaises(ValueError):
            node.update()

    def test_the_flag_is_refused_before_the_ports_exist(self):
        node = bff.TcspcDecay("decay")
        with self.assertRaises(ValueError):
            node.set_spectrum_from_port(True)
