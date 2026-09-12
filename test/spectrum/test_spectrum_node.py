"""The nodes that produce a lifetime spectrum, and the polarisation one.

`TcspcDecay` stopped holding an opinion about the photophysics when it grew a
`lifetime_spectrum` input port; these are the nodes on the other end of that
port. What has to be true of them is narrow and worth stating:

* the polarised spectrum is the one ChiSurf's `calculcate_spectrum` returns,
  **entry for entry and in the same order** -- a mixed channel is the union
  of two scaled spectra, so its length and its ordering are observable, and a
  fit built on a differently-ordered one is not wrong in a way any chi-square
  would show;
* `G` divides the perpendicular channel rather than multiplying it, which is
  a factor `G**2` this stack has already got wrong once;
* magic angle passes the spectrum through untouched, because a VM decay has
  no anisotropy in it at all;
* and the composition works: `LifetimeSpectrumNode -> AnisotropySpectrum ->
  TcspcDecay` is one graph and `GraphNode::update()` walks all of it.
"""

import math
import unittest

import numpy as np

from IMP import bff


def interleaved(amplitudes, times):
    out = np.empty(2 * len(amplitudes), dtype=float)
    out[0::2] = amplitudes
    out[1::2] = times
    return out


def make_source(spectrum, name="src"):
    """A `LifetimeSpectrumNode` holding *spectrum* on its scalar ports."""
    node = bff.PhotophysicsLifetimeSpectrumNode(name)
    node.set_number_of_lifetimes(len(spectrum) // 2)
    node.add_output_port(name, bff.GraphPort([0.0], False, True))
    for i in range(len(spectrum) // 2):
        node.get_input_port("a%d" % i).value = float(spectrum[2 * i])
        node.get_input_port("t%d" % i).value = float(spectrum[2 * i + 1])
    return node


def make_anisotropy(source, b, rho, r0=0.38, g=1.0, l1=0.0, l2=0.0,
                    polarization="vv", name="aniso"):
    node = bff.PhotophysicsAnisotropySpectrumNode(name)
    node.set_number_of_rotations(len(b))
    node.add_output_port(name, bff.GraphPort([0.0], False, True))
    node.set_polarization_name(polarization)
    node.get_input_port("lifetime_spectrum").link = source.get_output_port(
        source.get_name())
    node.get_input_port("r0").value = float(r0)
    node.get_input_port("g").value = float(g)
    node.get_input_port("l1").value = float(l1)
    node.get_input_port("l2").value = float(l2)
    for i in range(len(b)):
        node.get_input_port("b%d" % i).value = float(b[i])
        node.get_input_port("rho%d" % i).value = float(rho[i])
    return node


def evaluated(node):
    node.update()
    return np.asarray(node.get_output_port(node.get_name()).value, dtype=float)


class LifetimeSpectrumNodeTests(unittest.TestCase):

    def test_the_ports_are_the_spectrum(self):
        node = make_source([0.3, 1.5, 0.7, 4.0])
        np.testing.assert_allclose(evaluated(node), [0.3, 1.5, 0.7, 4.0])

    def test_a_lifetime_is_taken_absolute(self):
        """A lifetime walked through zero is a growing exponential.

        The optimiser has no way to know that; the model does, and says so
        here rather than letting a negative time constant into a convolution.
        """
        node = make_source([0.3, -1.5])
        np.testing.assert_allclose(evaluated(node), [0.3, 1.5])

    def test_normalising_the_amplitudes_is_opt_in(self):
        node = make_source([2.0, 1.0, 6.0, 4.0])
        node.set_normalize_amplitudes(True)
        np.testing.assert_allclose(evaluated(node), [0.25, 1.0, 0.75, 4.0])

    def test_absolute_amplitudes_is_opt_in(self):
        node = make_source([-2.0, 1.0])
        np.testing.assert_allclose(evaluated(node), [-2.0, 1.0])
        node.set_absolute_amplitudes(True)
        np.testing.assert_allclose(evaluated(node), [2.0, 1.0])

    def test_a_spectrum_with_no_species_is_refused(self):
        node = bff.PhotophysicsLifetimeSpectrumNode("empty")
        self.assertRaises(ValueError, node.set_number_of_lifetimes, 0)


class AnisotropySpectrumTests(unittest.TestCase):

    def test_magic_angle_passes_the_spectrum_through(self):
        """Not "returns something close": the *same* array.

        A VM model carrying an anisotropy that is merely small would fit
        almost as well and be a different model. It has none.
        """
        source = make_source([0.3, 1.5, 0.7, 4.0])
        node = make_anisotropy(source, [1.0], [2.0], polarization="vm")
        np.testing.assert_array_equal(evaluated(node), [0.3, 1.5, 0.7, 4.0])

    def test_the_unmixed_channels_start_where_the_anisotropy_says(self):
        """r(0) = r0, so VV starts at 1 + 2 r0 and VH at (1 - r0) / G.

        Checked on the decay at t = 0 rather than on the spectrum, because
        that is the statement anyone can verify against a textbook: the
        amplitudes of a spectrum sum to the initial intensity.
        """
        r0, g = 0.38, 1.0
        source = make_source([1.0, 4.0])
        for polarization, expected in (("vv", 1.0 + 2.0 * r0),
                                       ("vh", (1.0 - r0) / g)):
            spectrum = evaluated(make_anisotropy(
                source, [1.0], [10.0], r0=r0, g=g, polarization=polarization))
            self.assertAlmostEqual(float(spectrum[0::2].sum()), expected, 12)

    def test_g_divides_the_perpendicular_channel(self):
        """The sensitivity ratio, in the place that makes the pair invert.

        `G = S_par / S_perp`, so the perpendicular channel records `1/G` of
        what an equally sensitive one would. Multiplying instead -- which is
        the mistake available here -- puts two forward models a factor `G**2`
        apart, and both still fit.
        """
        source = make_source([1.0, 4.0])
        one = evaluated(make_anisotropy(source, [1.0], [10.0], g=1.0,
                                        polarization="vh"))
        two = evaluated(make_anisotropy(source, [1.0], [10.0], g=2.0,
                                        polarization="vh"))
        np.testing.assert_allclose(two[0::2], one[0::2] / 2.0, rtol=1e-14)
        np.testing.assert_allclose(two[1::2], one[1::2], rtol=1e-14)

    def test_the_amplitudes_are_normalised_to_r0(self):
        """`r0` alone sets r(0); the b_i only divide it up.

        Doubling every rotational amplitude is therefore *not* a different
        model -- which is what makes r0 identifiable at all.
        """
        source = make_source([1.0, 4.0])
        a = evaluated(make_anisotropy(source, [0.2, 0.3], [1.0, 10.0],
                                      r0=0.4, polarization="vv"))
        b = evaluated(make_anisotropy(source, [0.4, 0.6], [1.0, 10.0],
                                      r0=0.4, polarization="vv"))
        np.testing.assert_allclose(a, b, rtol=1e-14)

    def test_a_mixed_channel_is_a_union_and_not_a_sum(self):
        """l1 mixes VH into VV by *appending* its components, scaled.

        Adding element-wise -- the shape of the bug this guards -- would also
        add the time constants and so double the decay times, which reads as
        a slower dye rather than as an error.
        """
        source = make_source([1.0, 4.0])
        mixed = evaluated(make_anisotropy(source, [1.0], [10.0],
                                          polarization="vv", l1=0.1))
        # One species and one rotation make each unmixed channel two
        # components (the pass-through and the product), so the mixed one is
        # four -- their union, not their sum.
        self.assertEqual(mixed.size, 8)
        unmixed_vv = evaluated(make_anisotropy(source, [1.0], [10.0],
                                               polarization="vv", l1=0.0))[:4]
        np.testing.assert_allclose(mixed[:4:2], unmixed_vv[0::2] * 0.9,
                                   rtol=1e-14)
        # And the time constants are the unmixed ones. A sum would have
        # doubled them, which reads as a slower dye rather than as an error.
        np.testing.assert_allclose(mixed[1:4:2], unmixed_vv[1::2], rtol=1e-14)

    def test_vv_vh_is_the_two_channels_end_to_end(self):
        source = make_source([1.0, 4.0])
        both = evaluated(make_anisotropy(source, [1.0], [10.0],
                                         polarization="vv/vh", l1=0.1, l2=0.2))
        vv = evaluated(make_anisotropy(source, [1.0], [10.0],
                                       polarization="vv", l1=0.1, l2=0.2))
        vh = evaluated(make_anisotropy(source, [1.0], [10.0],
                                       polarization="vh", l1=0.1, l2=0.2))
        np.testing.assert_allclose(both, np.hstack([vv, vh]), rtol=1e-14)

    def test_the_product_is_the_harmonic_mean_of_the_time_constants(self):
        """Why a polarisation can be a *spectrum* transform at all.

        The measured decay is the fluorescence decay times the anisotropy
        decay, and the product of two exponentials is an exponential whose
        time constant is `1/(1/t + 1/rho)`. If that ever became something
        else, every polarised model would silently acquire a different
        rotational correlation time.
        """
        tau, rho = 4.0, 10.0
        source = make_source([1.0, tau])
        spectrum = evaluated(make_anisotropy(source, [1.0], [rho],
                                             polarization="vv"))
        # The pass-through term first, then the product term.
        self.assertAlmostEqual(float(spectrum[1]), tau, 12)
        self.assertAlmostEqual(float(spectrum[3]),
                               1.0 / (1.0 / tau + 1.0 / rho), 12)

    def test_an_unknown_polarisation_is_refused(self):
        """And is not quietly taken for magic angle.

        Defaulting would return the spectrum unchanged, which is a model with
        no anisotropy -- a wrong answer wearing the shape of a right one.
        """
        node = bff.PhotophysicsAnisotropySpectrumNode("a")
        node.set_number_of_rotations(1)
        self.assertRaises(ValueError, node.set_polarization_name, "vertical")

    def test_an_odd_spectrum_is_refused(self):
        node = bff.PhotophysicsAnisotropySpectrumNode("a")
        node.set_number_of_rotations(1)
        node.add_output_port("a", bff.GraphPort([0.0], False, True))
        node.set_polarization_name("vv")
        node.get_input_port("lifetime_spectrum").set_values_array(
            np.array([1.0, 4.0, 0.5]))
        self.assertRaises(ValueError, node.update)


class ChisurfParityTests(unittest.TestCase):
    """The same numbers as ChiSurf's own kernel, which is the contract.

    ChiSurf is a sibling checkout rather than an installed package, so this
    skips where it is absent; where it is present it is the test that decides
    whether a polarised fit may run on the graph at all.
    """

    def setUp(self):
        import sys
        from pathlib import Path
        sibling = Path(__file__).resolve().parents[3] / "chisurf"
        if sibling.is_dir() and str(sibling) not in sys.path:
            sys.path.insert(0, str(sibling))
        try:
            from chisurf.core.fluorescence.anisotropy.decay import (
                calculcate_spectrum)
        except Exception as exc:                       # pragma: no cover
            self.skipTest("ChiSurf is not importable: %s" % exc)
        self.reference = calculcate_spectrum

    def test_every_channel_matches_over_random_spectra(self):
        rng = np.random.default_rng(11)
        for _ in range(25):
            n_f = int(rng.integers(1, 4))
            n_b = int(rng.integers(1, 4))
            f = interleaved(rng.uniform(0.1, 1.0, n_f),
                            rng.uniform(0.2, 6.0, n_f))
            b = rng.uniform(0.05, 1.0, n_b)
            rho = rng.uniform(0.1, 20.0, n_b)
            r0 = float(rng.uniform(0.1, 0.4))
            g = float(rng.uniform(0.7, 1.8))
            l1 = float(rng.uniform(0.0, 0.1))
            l2 = float(rng.uniform(0.0, 0.1))
            # ChiSurf normalises the rotational amplitudes in the *getter*,
            # so its kernel is handed the normalised ones; the node does the
            # normalising itself, from the fitted values.
            normalised = np.abs(b) / np.abs(b).sum() * r0
            source = make_source(f)
            for polarization in ("vm", "vv", "vh", "vv/vh"):
                got = evaluated(make_anisotropy(
                    source, b, rho, r0=r0, g=g, l1=l1, l2=l2,
                    polarization=polarization))
                want = np.asarray(self.reference(
                    f, interleaved(normalised, rho), polarization, g, l1, l2),
                    dtype=float)
                self.assertEqual(got.shape, want.shape)
                np.testing.assert_allclose(got, want, rtol=1e-12, atol=1e-14)


class CompositionTests(unittest.TestCase):
    """One graph, walked by one `update()`, which is the whole point."""

    def test_the_decay_takes_the_polarised_spectrum(self):
        source = make_source([1.0, 4.0])
        aniso = make_anisotropy(source, [1.0], [10.0], polarization="vv")

        decay = bff.TCSPCDecay("decay")
        decay.set_number_of_lifetimes(1)
        decay.add_output_port("decay", bff.GraphPort([0.0], False, True))
        decay.set_spectrum_from_port(True)
        decay.get_input_port("lifetime_spectrum").link = \
            aniso.get_output_port("aniso")
        n = 64
        x = np.arange(n) * 0.1
        response = np.exp(-0.5 * ((x - 1.0) / 0.15) ** 2)
        decay.set_response_array(np.ascontiguousarray(response))
        decay.set_timing(0.1, 12.5)
        decay.update()

        # The decay used the polarised spectrum, not the source's, and it got
        # there without the caller evaluating anything in between.
        #
        # Compared against the *contributing* part of what the anisotropy
        # published, not against all of it: at `l1 = 0` a VV channel still
        # appends the VH components scaled by zero, and the decay drops
        # amplitudes that are exactly zero because a species multiplied by
        # zero costs a full recursion over every channel to add nothing. That
        # is bit-exact, which is why it is the default.
        upstream = evaluated(aniso)
        contributing = upstream.reshape(-1, 2)
        contributing = contributing[contributing[:, 0] != 0.0].ravel()
        np.testing.assert_allclose(
            np.asarray(decay.get_lifetime_spectrum(), dtype=float),
            contributing, rtol=1e-14)
        self.assertTrue(np.all(np.isfinite(decay.get_curve())))

    def test_moving_a_source_parameter_moves_the_decay(self):
        """Invalidation reaches the far end, which is what makes it a graph.

        A chain that had to be evaluated by hand in the right order would
        work exactly once and then quietly serve a stale curve to the
        optimiser.
        """
        source = make_source([1.0, 4.0])
        aniso = make_anisotropy(source, [1.0], [10.0], polarization="vv")
        before = evaluated(aniso).copy()
        source.get_input_port("t0").value = 2.0
        after = evaluated(aniso)
        self.assertFalse(np.allclose(before, after))
        self.assertAlmostEqual(float(after[1]), 2.0, 12)


if __name__ == "__main__":
    unittest.main()


class AmplitudeThresholdTests(unittest.TestCase):
    """Dropping species that cannot pay for themselves.

    A multi-exponential reconvolution costs one serial recursion over every
    channel *per species*, so it is linear in the species count -- and a
    model whose spectrum came from a *distribution* has as many species as
    that distribution has bins, whatever their weight. On the FRET decay this
    module builds, 97 species carry the weight of 53.
    """

    def make_decay(self, spectrum, threshold=0.0, n=64):
        decay = bff.TCSPCDecay("decay")
        decay.set_number_of_lifetimes(len(spectrum) // 2)
        decay.add_output_port("decay", bff.GraphPort([0.0], False, True))
        decay.set_spectrum_from_port(True)
        decay.get_input_port("lifetime_spectrum").set_values_array(
            np.ascontiguousarray(np.asarray(spectrum, dtype=float)))
        x = np.arange(n) * 0.1
        decay.set_response_array(
            np.ascontiguousarray(np.exp(-0.5 * ((x - 1.0) / 0.15) ** 2)))
        decay.set_timing(0.1, 12.5)
        decay.set_amplitude_threshold(threshold)
        return decay

    def test_the_default_is_exact(self):
        """Zero drops only exact zeros, so the curve cannot move.

        A library that silently drops terms is the failure this module spent
        2026-09-01 fixing one layer up; the caller that knows what its data
        are worth chooses the threshold, and the default chooses nothing.
        """
        decay = self.make_decay([1.0, 4.0, 0.5, 1.0])
        self.assertEqual(decay.get_amplitude_threshold(), 0.0)
        decay.update()
        self.assertEqual(decay.get_number_of_active_lifetimes(), 2)

    def test_an_exactly_zero_amplitude_is_dropped_and_the_curve_is_identical(self):
        kept = self.make_decay([1.0, 4.0, 0.5, 1.0])
        kept.update()
        padded = self.make_decay([1.0, 4.0, 0.0, 7.0, 0.5, 1.0, 0.0, 0.25])
        padded.update()
        self.assertEqual(padded.get_number_of_active_lifetimes(), 2)
        np.testing.assert_array_equal(np.asarray(padded.get_curve(), float),
                                      np.asarray(kept.get_curve(), float))

    def test_the_threshold_is_relative_to_the_largest_amplitude(self):
        decay = self.make_decay([1.0, 4.0, 1e-13, 1.0], threshold=1e-12)
        decay.update()
        self.assertEqual(decay.get_number_of_active_lifetimes(), 1)
        # Scale everything: the same species must still be dropped, because
        # the threshold is about proportion and not about magnitude.
        scaled = self.make_decay([1e6, 4.0, 1e-7, 1.0], threshold=1e-12)
        scaled.update()
        self.assertEqual(scaled.get_number_of_active_lifetimes(), 1)

    def test_pruning_does_not_consume_the_ports(self):
        """The port count must survive an evaluation that prunes.

        On the scalar path nothing re-reads the component count, so pruning
        it in place would make a component with a momentarily zero amplitude
        disappear *permanently* -- and an optimiser walking an amplitude
        through zero does exactly that.
        """
        decay = bff.TCSPCDecay("decay")
        decay.set_number_of_lifetimes(2)
        decay.add_output_port("decay", bff.GraphPort([0.0], False, True))
        n = 64
        x = np.arange(n) * 0.1
        decay.set_response_array(
            np.ascontiguousarray(np.exp(-0.5 * ((x - 1.0) / 0.15) ** 2)))
        decay.set_timing(0.1, 12.5)
        decay.get_input_port("a0").value = 1.0
        decay.get_input_port("t0").value = 4.0
        decay.get_input_port("a1").value = 0.0        # walks through zero
        decay.get_input_port("t1").value = 1.0
        decay.update()
        self.assertEqual(decay.get_number_of_active_lifetimes(), 1)
        self.assertEqual(decay.get_number_of_lifetimes(), 2)

        decay.get_input_port("a1").value = 0.5        # and comes back
        decay.update()
        self.assertEqual(decay.get_number_of_active_lifetimes(), 2)
        np.testing.assert_allclose(
            np.asarray(decay.get_lifetime_spectrum(), float),
            [1.0, 4.0, 0.5, 1.0], rtol=1e-14)

    def test_a_threshold_that_would_empty_the_spectrum_keeps_the_largest(self):
        """An empty spectrum is not something the kernel can be handed."""
        decay = self.make_decay([1.0, 4.0, 0.5, 1.0], threshold=2.0)
        decay.update()
        self.assertEqual(decay.get_number_of_active_lifetimes(), 1)
        np.testing.assert_allclose(
            np.asarray(decay.get_lifetime_spectrum(), float), [1.0, 4.0])

    def test_a_negative_threshold_is_refused(self):
        decay = bff.TCSPCDecay("decay")
        decay.set_number_of_lifetimes(1)
        self.assertRaises(ValueError, decay.set_amplitude_threshold, -1e-12)

    def test_the_threshold_invalidates(self):
        """It changes the value the node publishes, so it must be seen."""
        decay = self.make_decay([1.0, 4.0, 1e-13, 1.0])
        decay.update()
        self.assertEqual(decay.get_number_of_active_lifetimes(), 2)
        decay.set_amplitude_threshold(1e-12)
        decay.update()
        self.assertEqual(decay.get_number_of_active_lifetimes(), 1)
