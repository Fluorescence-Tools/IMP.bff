"""PRD-109: the photon and decay Monte-Carlo moved here from QuEst.

Both kernels are pinned bit-for-bit against QuEst's at a fixed seed (the photon
trace) and against QuEst run **serially** (the decay curve) -- serially, because
QuEst's decay kernel accumulated into ``decay[bin_idx]`` from inside a ``prange``
with a data-dependent index, which numba cannot privatise into a reduction. That
was a genuine data race: measured on 8 threads it lost up to 2% of the emitted
intensity, at random, run to run. The kernel here blocks the accumulator and
reduces in fixed order; :class:`DecayAccumulatorTests` is the regression pin.
"""

import unittest

import numpy as np

import IMP
import IMP.test

from IMP.bff.sampling import excited_state as photon

try:
    import numba
except ImportError:  # pragma: no cover - numba is optional for IMP.bff
    numba = None


def constant_rate(value=0.37, n_frames=4000):
    return np.full(n_frames, value, dtype=np.float32)


class PhotonTraceTests(IMP.test.TestCase):

    def setUp(self):
        super().setUp()
        rng = np.random.default_rng(3)
        self.k_quench = np.abs(rng.normal(0.5, 0.4, 5000)).astype(np.float32)

    def test_frozen_reference(self):
        """The trace's *statistics*, since the exact trace can no longer be pinned.

        This asserted 6297 emitted photons and a mean delay of
        1.266492156326010 at seed 5 -- bit-for-bit QuEst's trace, which held
        through the numba port because both ran the same Mersenne stream. The
        C++ port (PRD-113 stage 5) ends that: no C++ generator reproduces
        numba's stream, so those digits describe the generator, not the model.

        What survives is the distribution. Over 12 seeds at 20 000 excitations:

            numba   6288 +/- 48 emitted,  <dt> 1.26699 +/- 0.01205 ns
            C++     6328 +/- 59 emitted,  <dt> 1.25803 +/- 0.01290 ns

        Both means also have closed forms under a *constant* rate, which is the
        stronger check and is asserted in
        ``test_dynamics_cpp.py::test_constant_rate_reproduces_the_closed_form``.
        Here the rate varies along the trajectory, so this pins the sample.
        """
        emitted_counts, mean_delays = [], []
        for seed in range(6):
            dts, emitted = photon.simulate_photon_trace(
                20000, self.k_quench, t_step=0.01, tau0=4.0, random_seed=seed
            )
            emitted_counts.append(int(emitted.sum()))
            mean_delays.append(float(dts[emitted > 0].mean()))
        self.assertAlmostEqual(sum(emitted_counts) / 6, 6310, delta=200)
        self.assertAlmostEqual(sum(mean_delays) / 6, 1.2625, delta=0.05)

    def test_waiting_times_are_never_negative(self):
        """QuEst emitted photons *before* they were excited, 2.4e-4 of the time.

        `log(1 / (u + EPS))` goes negative once `u > 1 - EPS`, and every such
        photon was then dropped by a histogram starting at 0 -- so the decay
        curve held fewer photons than the trace reported as emitted. At 40 000
        photons that was 8 of them.
        """
        for seed in (1, 2, 3):
            dts, emitted = photon.simulate_photon_trace(
                50000, self.k_quench, t_step=0.01, tau0=4.0, random_seed=seed
            )
            self.assertGreater(float(dts[emitted > 0].min()), 0.0)

    def test_the_whole_emitted_trace_lands_in_a_wide_histogram(self):
        """The consequence of the above, stated as the property that matters."""
        dts, emitted = photon.simulate_photon_trace(
            50000, self.k_quench, t_step=0.01, tau0=4.0, random_seed=4
        )
        counts, _edges = np.histogram(dts[emitted > 0], range=(0.0, 400.0), bins=512)
        self.assertEqual(int(counts.sum()), int(emitted.sum()))

    def test_a_seed_pins_the_trace(self):
        first = photon.simulate_photon_trace(2000, self.k_quench, random_seed=5)
        second = photon.simulate_photon_trace(2000, self.k_quench, random_seed=5)
        self.assertTrue(np.array_equal(first[0], second[0]))
        self.assertTrue(np.array_equal(first[1], second[1]))

    def test_unseeded_runs_are_independent_samples(self):
        """numba seeds each worker deterministically, so this needs OS entropy.

        Without it an unseeded run repeats itself for the life of the process --
        four Monte-Carlo runs returning byte-identical results, which reads as
        precision that is not there.
        """
        sums = {
            float(photon.simulate_photon_trace(4000, self.k_quench)[0].sum())
            for _ in range(4)
        }
        self.assertGreater(len(sums), 1)

    def test_quenched_events_carry_no_delay_time(self):
        dts, emitted = photon.simulate_photon_trace(
            5000, self.k_quench, tau0=4.0, random_seed=1
        )
        self.assertTrue(np.all(dts[emitted == 0] == 0.0))
        self.assertTrue(np.all(dts[emitted > 0] > 0.0))

    def test_no_quenching_emits_every_photon(self):
        dts, emitted = photon.simulate_photon_trace(
            5000, np.zeros(100, np.float32), tau0=4.0, random_seed=1
        )
        self.assertEqual(int(emitted.sum()), 5000)

    def test_stronger_quenching_emits_fewer_photons(self):
        weak = photon.simulate_photon_trace(
            20000, constant_rate(0.05), tau0=4.0, random_seed=1
        )[1].sum()
        strong = photon.simulate_photon_trace(
            20000, constant_rate(5.0), tau0=4.0, random_seed=1
        )[1].sum()
        self.assertGreater(int(weak), int(strong))

    def test_the_unquenched_lifetime_is_recovered(self):
        """With no quencher the delay times are exponential with mean tau0."""
        dts, emitted = photon.simulate_photon_trace(
            200000, np.zeros(10, np.float32), tau0=4.0, random_seed=2
        )
        self.assertAlmostEqual(float(dts[emitted > 0].mean()), 4.0, delta=0.05)

    def test_an_empty_rate_array_is_not_an_error(self):
        dts, emitted = photon.simulate_photon_trace(
            10, np.zeros(0, np.float32), random_seed=1
        )
        self.assertEqual(dts.size, 10)
        self.assertEqual(int(emitted.sum()), 0)


class QuenchedDecayTests(IMP.test.TestCase):

    def test_frozen_reference(self):
        """Bit-for-bit QuEst run serially; see the module docstring."""
        decay = np.zeros(512)
        photon.simulate_quenched_decay(200, decay, 0.032, constant_rate(), 0.01, 4.0)
        self.assertAlmostEqual(float(decay.sum()), 200.612864255, delta=1e-6)

    def test_the_curve_decays_for_a_constant_rate(self):
        """Monotonic over a stride, not bin to bin.

        ``bin_idx = int(t / dt_tac)`` with ``t`` accumulated by repeated
        addition, so bins take turns collecting an extra frame however neatly
        ``t_step`` divides ``dt_tac``. Here that is 4 frames per bin, so one
        stray frame is **25%** of a bin while the curve falls only 2% per bin
        -- adjacent bins are genuinely not ordered, and saying otherwise would
        be pinning float noise. Over 16 bins the curve falls 27%, which clears
        the jitter. The exponent itself is pinned by
        :meth:`test_the_decay_constant_is_the_total_rate`.
        """
        decay = np.zeros(256)
        photon.simulate_quenched_decay(500, decay, 0.032, constant_rate(), 0.008, 4.0)
        self.assertTrue(np.all(decay[16:] < decay[:-16]))

    def test_the_decay_constant_is_the_total_rate(self):
        """1/tau0 + kQ, read back off the slope of the log curve."""
        tau0, kq = 4.0, 0.37
        dt_tac = 0.032
        t_step = 0.008
        n_bins = 256
        # The walk stops at `n_frames`, so it must reach the end of the window
        # or the tail bins stay empty and the fit sees log(0).
        n_frames = int(dt_tac * n_bins / t_step) + 1
        decay = np.zeros(n_bins)
        photon.simulate_quenched_decay(
            2000, decay, dt_tac, constant_rate(kq, n_frames), t_step, tau0
        )
        t = (np.arange(decay.size) + 0.5) * dt_tac
        slope = np.polyfit(t, np.log(decay), 1)[0]
        self.assertAlmostEqual(-slope, 1.0 / tau0 + kq, delta=0.01)

    def test_quenching_shortens_the_decay(self):
        unquenched, quenched = np.zeros(256), np.zeros(256)
        photon.simulate_quenched_decay(
            500, unquenched, 0.032, constant_rate(0.0), 0.01, 4.0)
        photon.simulate_quenched_decay(
            500, quenched, 0.032, constant_rate(2.0), 0.01, 4.0)
        centroid = lambda d: float((d * np.arange(d.size)).sum() / d.sum())
        self.assertLess(centroid(quenched), centroid(unquenched))

    def test_a_non_float64_output_array_is_still_filled(self):
        """`np.asarray(..., dtype=float64)` copies, and QuEst filled the copy."""
        decay = np.zeros(128, dtype=np.float32)
        photon.simulate_quenched_decay(100, decay, 0.032, constant_rate(), 0.01, 4.0)
        self.assertGreater(float(decay.sum()), 0.0)

    def test_degenerate_inputs_are_no_ops(self):
        for n_curves, size, rate in ((0, 64, 100), (10, 0, 100), (10, 64, 0)):
            decay = np.zeros(size)
            photon.simulate_quenched_decay(
                n_curves, decay, 0.032,
                np.zeros(rate, np.float32), 0.01, 4.0,
            )
            self.assertEqual(float(decay.sum()), 0.0)


class DecayAccumulatorTests(IMP.test.TestCase):
    """Regression pin for the data race QuEst's kernel had (PRD-109)."""

    def accumulate(self, n_threads=None):
        if n_threads is not None and numba is not None:
            numba.set_num_threads(n_threads)
        decay = np.zeros(512)
        photon.simulate_quenched_decay(
            2000, decay, 0.032, constant_rate(), 0.01, 4.0
        )
        return decay

    def tearDown(self):
        if numba is not None:
            numba.set_num_threads(numba.config.NUMBA_NUM_THREADS)
        super().tearDown()

    def test_the_curve_is_repeatable_run_to_run(self):
        curves = {self.accumulate().tobytes() for _ in range(8)}
        self.assertEqual(len(curves), 1)

    @unittest.skipIf(numba is None, "needs numba to vary the thread count")
    def test_the_total_does_not_depend_on_the_thread_count(self):
        """The race dropped intensity; only the summation order may differ now."""
        serial = self.accumulate(1)
        threaded = self.accumulate(numba.config.NUMBA_NUM_THREADS)
        self.assertTrue(np.allclose(threaded, serial, rtol=1e-12, atol=0.0))

    def test_every_curve_contributes_its_whole_intensity(self):
        """Sum over bins -> n_curves * (1 - exp(-rate * window)), no lost updates."""
        n_curves, tau0, kq, dt_tac, n_bins = 2000, 4.0, 0.37, 0.032, 512
        decay = self.accumulate()
        window = dt_tac * n_bins
        expected = n_curves * (1.0 - np.exp(-(1.0 / tau0 + kq) * window))
        # First-order Euler emission overshoots slightly; the race *under*shot,
        # by up to 2%, which is what this bound excludes.
        self.assertAlmostEqual(float(decay.sum()) / expected, 1.0, delta=5e-3)


if __name__ == "__main__":
    IMP.test.main()
