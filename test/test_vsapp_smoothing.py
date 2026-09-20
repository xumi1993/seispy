"""Depth-unit, boundary and compatibility checks for gradient smoothing."""
import unittest

import numpy as np

from seispy.inversion import smooth_gradient
from seispy.signal import smooth


class TestGaussianSmoothing(unittest.TestCase):
    def test_impulse_has_requested_standard_deviation_in_km(self):
        depths = np.arange(101) * 0.1
        impulse = np.zeros(101)
        impulse[50] = 1.
        result = smooth_gradient(impulse, depths, 0.3)
        offsets = np.arange(-12, 13) * 0.1
        expected = np.exp(-0.5 * (offsets / 0.3)**2)
        expected /= expected.sum()
        np.testing.assert_allclose(result[38:63], expected, atol=1e-15)
        np.testing.assert_allclose(result.sum(), 1.)
        self.assertLess(result[50], 1.)

    def test_fixed_entries_do_not_dilute_constant_gradient(self):
        depths = np.array([0., 0.1, 0.4, 1., 1.3])
        fixed = np.array([False, False, True, False, True])
        values = np.where(fixed, 1e9, 2.)
        result = smooth_gradient(values, depths, 0.4, fixed_mask=fixed)
        np.testing.assert_allclose(result[~fixed], 2.)
        np.testing.assert_array_equal(result[fixed], 0.)
        scaled = smooth_gradient(values, 10 * depths, 4., fixed_mask=fixed)
        np.testing.assert_allclose(result, scaled)
        np.testing.assert_array_equal(values[fixed], 1e9)

    def test_short_arrays_and_wide_gaussian_windows(self):
        for size in (1, 2, 4):
            np.testing.assert_allclose(smooth(np.ones(size), 20, 'gaussian', sigma=5.), 1.)
        with self.assertRaises(ValueError):
            smooth(np.ones(5), window='gaussian', sigma=0.)
        with self.assertRaises(ValueError):
            smooth(np.ones(5), window='gaussian')

    def test_existing_window_results_are_unchanged(self):
        x = np.arange(12.)**2
        half_len = 2
        padded = np.r_[x[4:0:-1], x, x[-1:-5:-1]]
        for window in ('flat', 'hanning', 'hamming', 'bartlett', 'blackman'):
            weights = np.ones(5) if window == 'flat' else getattr(np, window)(5)
            expected = np.convolve(weights / weights.sum(), padded, mode='valid')[2:-2]
            np.testing.assert_array_equal(smooth(x, half_len, window), expected)
