"""Public measurements take only radial RFs and reconstruct their reference."""

import inspect
import unittest

import numpy as np
from seispy import RFStation
from seispy.decon import deconit
from seispy.vsapp import (
    VsAppError,
    _compute_vsapp_from_components,
    compute_station_vsapp,
    compute_vsapp,
    gaussian_vertical,
)


class TestRadialVsApp(unittest.TestCase):
    @staticmethod
    def native_pair(npts=512, dt=0.05, shift=5.0, f0=2.0):
        t = np.arange(npts) * dt
        source = np.exp(-((t - 3.0) / 0.4) ** 2)
        numerator = 0.3 * source + 0.1 * np.roll(source, 5)
        radial = deconit(numerator, source, dt, tshift=shift, f0=f0)[0]
        vertical = deconit(source, source, dt, tshift=shift, f0=f0)[0]
        return t - shift, radial, vertical

    def assert_curve_equal(self, actual, expected, rtol=3e-12, atol=3e-12):
        for field in ('periods_s', 'radial_integral', 'vertical_integral', 'ratio',
                      'angle_deg', 'vs_km_s', 'radial_zero', 'vertical_zero', 'vs0_km_s'):
            np.testing.assert_allclose(getattr(actual, field), getattr(expected, field),
                                       rtol=rtol, atol=atol)
        self.assertEqual(actual.status, expected.status)
        self.assertEqual(actual.vs0_status, expected.vs0_status)

    def test_default_iterative_reference_matches_actual_native_deconvolution(self):
        for npts, f0 in ((512, 2.0), (301, 3.5)):
            with self.subTest(npts=npts, f0=f0):
                times, radial, vertical = self.native_pair(npts=npts, f0=f0)
                actual = compute_vsapp(times, radial, rayp=0.06,
                                       periods_s=[0.1, 0.73, 2.0], f0=f0)
                expected = _compute_vsapp_from_components(
                    times, radial, vertical, rayp=0.06, periods_s=[0.1, 0.73, 2.0],
                )
                self.assert_curve_equal(actual, expected)

    def test_ideal_gaussian_references_preserve_their_amplitude_conventions(self):
        times = np.arange(301) * 0.05 - 5.0
        ratio, rayp = 0.4, 0.065
        expected_speed = np.sin(np.arctan(ratio) / 2) / rayp
        for convention in ('area', 'peak'):
            with self.subTest(convention=convention):
                radial = ratio * gaussian_vertical(times, 2.7, normalization=convention)
                actual = compute_vsapp(times, radial, rayp=rayp,
                                       periods_s=[0.1, 0.43, 3.0], f0=2.7,
                                       reference='gaussian-' + convention)
                np.testing.assert_allclose(actual.vs_km_s, expected_speed)
                self.assertAlmostEqual(actual.vs0_km_s, expected_speed)

    def test_radial_only_signatures_and_calls(self):
        for function in (compute_vsapp, compute_station_vsapp, RFStation.compute_vsapp):
            self.assertNotIn('vertical', inspect.signature(function).parameters)
        times, radial, vertical = self.native_pair()
        kwargs = dict(rayp=0.06, periods_s=[0.2, 1.0], f0=2.0)
        with self.assertRaises(TypeError):
            compute_vsapp(times, radial, vertical, **kwargs)
        with self.assertRaises(TypeError):
            compute_vsapp(times, radial, vertical=vertical, **kwargs)
        with self.assertRaises(TypeError):
            compute_vsapp(times, radial, rayp=0.06, periods_s=[0.2, 1.0])

    def test_gaussian_factor_validation(self):
        times, radial, _ = self.native_pair()
        for reference in ('seispy-iter', 'gaussian-area', 'gaussian-peak'):
            for factor in (0.0, -1.0, np.nan, np.inf, None, 'invalid', 20.0):
                with self.subTest(reference=reference, factor=factor):
                    with self.assertRaises(VsAppError):
                        compute_vsapp(times, radial, rayp=0.06, periods_s=[0.2],
                                      f0=factor, reference=reference)

    def test_reference_and_timing_option_conflicts(self):
        times, radial, _ = self.native_pair()
        cases = [
            {'reference': None}, {'reference': 'measured'}, {'reference': 'water'},
            {'reference': 'gaussian-area', 'deconvolution_npts': 1024},
            {'reference': 'gaussian-peak', 'deconvolution_shift_samples': 100},
            {'deconvolution_npts': 511}, {'deconvolution_npts': True},
            {'deconvolution_shift_samples': -1},
            {'deconvolution_shift_samples': np.bool_(True)},
        ]
        for options in cases:
            with self.subTest(options=options), self.assertRaises(VsAppError):
                compute_vsapp(times, radial, rayp=0.06, periods_s=[0.2], f0=2.0, **options)

    def test_tail_truncation_retains_original_deconvolution_length(self):
        # A broad pulse and short saved record make the FFT length observable.
        times, radial, vertical = self.native_pair(npts=128, shift=1.0, f0=0.5)
        times, radial, vertical = times[:48], radial[:48], vertical[:48]
        actual = compute_vsapp(times, radial, rayp=0.06, periods_s=[0.1, 0.6, 0.8],
                               f0=0.5, deconvolution_npts=128)
        expected = _compute_vsapp_from_components(
            times, radial, vertical, rayp=0.06, periods_s=[0.1, 0.6, 0.8],
        )
        self.assert_curve_equal(actual, expected)
        assumed_saved_length = compute_vsapp(
            times, radial, rayp=0.06, periods_s=[0.1, 0.6, 0.8], f0=0.5,
        )
        self.assertGreater(np.max(np.abs(assumed_saved_length.vs_km_s - actual.vs_km_s)),
                           1e-3)

    def test_explicit_original_shift_handles_native_truncation_and_non_grid_arrivals(self):
        for shift, index in ((0.3, 2), (0.35, 3)):
            with self.subTest(shift=shift):
                times, radial, vertical = self.native_pair(dt=0.1, shift=shift)
                actual = compute_vsapp(times, radial, rayp=0.06, periods_s=[0.22, 0.25],
                                       f0=2.0, deconvolution_shift_samples=index)
                expected = _compute_vsapp_from_components(
                    times, radial, vertical, rayp=0.06, periods_s=[0.22, 0.25],
                )
                self.assert_curve_equal(actual, expected)

    def test_default_shift_rounds_non_grid_arrivals_without_changing_input(self):
        # Test both sides of a sample and offsets exceeding the old 1e-5 tolerance.
        for shift, index in ((5.00002, 100), (5.019, 100), (5.031, 101)):
            with self.subTest(shift=shift):
                times, radial, _ = self.native_pair(dt=0.05, shift=shift)
                original_times, original_radial = times.copy(), radial.copy()
                actual = compute_vsapp(times, radial, rayp=0.06,
                                       periods_s=[0.11, 0.7, 2.0], f0=2.0)
                expected = compute_vsapp(times, radial, rayp=0.06,
                                         periods_s=[0.11, 0.7, 2.0], f0=2.0,
                                         deconvolution_shift_samples=index)
                self.assert_curve_equal(actual, expected)
                np.testing.assert_array_equal(times, original_times)
                np.testing.assert_array_equal(radial, original_radial)

    def test_inputs_are_not_mutated_and_result_arrays_are_owned(self):
        times, radial, _ = self.native_pair()
        periods = np.array([0.1, 0.73, 2.0])
        snapshots = [array.copy() for array in (times, radial, periods)]
        for array in (times, radial, periods):
            array.setflags(write=False)
        result = compute_vsapp(times, radial, rayp=0.06, periods_s=periods, f0=2.0)
        for array, snapshot in zip((times, radial, periods), snapshots, strict=True):
            np.testing.assert_array_equal(array, snapshot)
        self.assertFalse(np.shares_memory(result.periods_s, periods))
        for array in (result.periods_s, result.radial_integral, result.vertical_integral,
                      result.ratio, result.angle_deg, result.vs_km_s):
            self.assertFalse(array.flags.writeable)


if __name__ == '__main__':
    unittest.main()
