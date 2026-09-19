"""Scientific regression checks against independent analytic models and quadrature."""

import unittest

import numpy as np
from scipy.integrate import quad
from seispy.vsapp import VsAppError, gaussian_vertical, seispy_iter_vertical
from seispy.vsapp import _compute_vsapp_from_components as compute_vsapp


class CoreTests(unittest.TestCase):
    def setUp(self):
        self.time = np.arange(-1200, 1201) * 0.01
        self.z = gaussian_vertical(self.time, 3.5, normalization="area")
        self.periods = np.array([0.1, 0.3, 0.7, 1.6, 4.0, 10.0])

    def measure(self, r=None, z=None, times=None, periods=None):
        return compute_vsapp(
            self.time if times is None else times,
            0.35 * self.z if r is None else r,
            self.z if z is None else z,
            rayp=0.06,
            periods_s=self.periods if periods is None else periods,
        )

    def test_halfspace_velocity_recovered_at_every_period(self):
        # Free-surface half-space relation independently generates the R/Z ratio.
        for speed in (0.8, 2.4, 3.6, 4.5):
            x = 0.06 * speed
            ratio = 2 * x * np.sqrt(1 - x * x) / (1 - 2 * x * x)
            result = self.measure(r=ratio * self.z)
            np.testing.assert_allclose(result.vs_km_s, speed, atol=1e-12)
            self.assertAlmostEqual(result.vs0_km_s, speed, places=12)

    def test_delayed_phases_against_independent_adaptive_quadrature(self):
        a, p = 3.5, 0.06

        def vertical(t):
            return a / np.sqrt(np.pi) * np.exp(-((a * t) ** 2))

        def radial(t):
            return 0.3 * vertical(t) + 0.18 * vertical(t - 1.1) - 0.1 * vertical(t - 2.6)

        result = self.measure(r=radial(self.time))
        expected = []
        for period in self.periods:

            def weight(t, half_width=period):
                return np.cos(np.pi * t / (2 * half_width)) ** 2

            r = quad(lambda t: radial(t) * weight(t), -period, period, epsabs=1e-12)[0]
            z = quad(lambda t: vertical(t) * weight(t), -period, period, epsabs=1e-12)[0]
            expected.append(np.sin(np.arctan(r / z) / 2) / p)
        np.testing.assert_allclose(result.vs_km_s, expected, atol=2e-7, rtol=1e-6)

    def test_common_amplitude_scale_cancels(self):
        baseline = self.measure()
        for scale in (1e-12, 1e12):
            scaled = self.measure(r=0.35 * self.z * scale, z=self.z * scale)
            np.testing.assert_allclose(scaled.vs_km_s, baseline.vs_km_s, atol=1e-12)

    def test_area_and_peak_conventions_are_equivalent_if_both_components_scaled(self):
        peak = 3.5 / np.sqrt(np.pi)
        other_z = gaussian_vertical(self.time, 3.5, normalization="peak")
        result = self.measure(r=0.35 * self.z / peak, z=other_z)
        np.testing.assert_allclose(result.vs_km_s, self.measure().vs_km_s, atol=1e-12)

    def test_seispy_discrete_self_deconvolution_matches_resolved_analytic_gaussian(self):
        times = np.arange(12001) * 0.01 - 20
        z = seispy_iter_vertical(times, 3.5)
        analytic = gaussian_vertical(times, 3.5, normalization="area")
        np.testing.assert_allclose(z, analytic, atol=1e-11)

    def test_off_sample_zero_and_noninteger_kernel_boundaries(self):
        time = self.time + 0.003
        z = gaussian_vertical(time, 3.5, normalization="peak")
        result = self.measure(times=time, r=0.4 * z, z=z, periods=[0.105, 0.337])
        speed = np.sin(np.arctan(0.4) / 2) / 0.06
        np.testing.assert_allclose(result.vs_km_s, speed, atol=1e-12)
        self.assertAlmostEqual(result.vs0_km_s, speed, places=12)

    def test_negative_signal_is_flagged_not_flipped(self):
        result = self.measure(r=-self.z)
        self.assertTrue(np.isnan(result.vs_km_s).all())
        self.assertEqual(result.vs0_status, "nonpositive_radial")
        self.assertEqual(set(result.status), {"nonpositive_radial"})
        self.assertTrue(np.all(result.ratio < 0))

    def test_negative_vertical_is_not_valid(self):
        result = self.measure(z=-self.z)
        self.assertEqual(set(result.status), {"invalid_vertical"})

    def test_cancellation_in_vertical_integral_is_flagged(self):
        z = self.z - np.roll(self.z, 200) * 10
        result = self.measure(z=z)
        self.assertEqual(result.status[-1], "invalid_vertical")

    def test_invalid_arrays_periods_and_windows_raise(self):
        for kwargs in (
            {"times": self.time[::-1]},
            {"times": self.time + 50},
            {"r": self.z[:-1]},
            {"r": np.full_like(self.z, np.nan)},
            {"z": np.zeros_like(self.z)},
            {"periods": [13]},
            {"periods": [0.001]},
            {"periods": [0]},
            {"periods": [1, 0.5]},
            {"periods": [1, 1]},
            {"periods": []},
            {"r": [[1, 2, 3]]},
        ):
            with self.subTest(kwargs=list(kwargs)):
                with self.assertRaises(VsAppError):
                    self.measure(**kwargs)

    def test_invalid_ray_parameter_and_irregular_time_raise(self):
        for p in (0, -1, np.nan, np.inf):
            with self.assertRaises(VsAppError):
                compute_vsapp(self.time, self.z, self.z, rayp=p, periods_s=[1])
        times = self.time.copy()
        times[100] += 0.003
        with self.assertRaises(VsAppError):
            self.measure(times=times)

    def test_no_input_mutation_and_owned_readonly_results(self):
        before = self.z.copy()
        result = self.measure()
        np.testing.assert_array_equal(before, self.z)
        self.assertTrue(self.periods.flags.writeable)
        self.assertFalse(result.periods_s.flags.writeable)
        self.assertFalse(result.vs_km_s.flags.writeable)

    def test_nonfinite_derived_velocity_is_not_marked_ok(self):
        result = compute_vsapp(
            self.time, self.z, self.z, rayp=1e-320, periods_s=[1.0]
        )
        self.assertEqual(result.status, ("nonfinite_velocity",))
        self.assertTrue(np.isnan(result.vs0_km_s))


if __name__ == "__main__":
    unittest.main()
