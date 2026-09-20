"""Adjacent-window misfit convergence, excluding the initial model and failed trials."""

import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import test_vsapp_sd_angle as sd_fixtures
from seispy.inversion import invert_vsapp
from seispy.inversion.forward import VsappForward
from seispy.inversion.optimize import windowed_misfit_change
from test_vsapp_inversion import PERIODS, make_model


def history(misfits):
    return [SimpleNamespace(misfit=value) for value in misfits]


class TestWindowedMisfit(unittest.TestCase):
    def test_requires_two_full_windows_of_updates_excluding_iteration_zero(self):
        self.assertIsNone(windowed_misfit_change(history([999., 1., 1., 1.]), 2))
        self.assertEqual(windowed_misfit_change(history([999., 1., 1., 1., 1.]), 2), 0.)

    def test_adjacent_nonoverlapping_window_means(self):
        # Earlier mean B=3, recent mean A=4.5. Old history must not affect the result.
        value = windowed_misfit_change(history([999., 1e8, 2., 4., 4., 5.]), 2)
        self.assertAlmostEqual(value, .5)

    def test_absolute_change_detects_increase_and_decrease_symmetrically(self):
        up = windowed_misfit_change(history([99., 2., 2., 3., 3.]), 2)
        down = windowed_misfit_change(history([99., 2., 2., 1., 1.]), 2)
        self.assertAlmostEqual(up, .5)
        self.assertAlmostEqual(down, .5)

    def test_zero_denominator_and_extreme_finite_misfits(self):
        self.assertEqual(windowed_misfit_change(history([99., 0., 0., 0., 0.]), 2), 0.)
        self.assertTrue(np.isinf(windowed_misfit_change(history([99., 0., 0., 1., 1.]), 2)))
        for scale in [1e-250, 1e250]:
            self.assertAlmostEqual(windowed_misfit_change(
                history([99.] + [2*scale, 2*scale, 3*scale, 3*scale]), 2), .5)


class TestWindowStopIntegration(unittest.TestCase):
    def run_misfits(self, misfits, **kwargs):
        fixture = sd_fixtures.TestSDStepControl()
        return fixture.run_sd([[1., 0., 0.]] * len(misfits), misfits=misfits,
                              max_step_lnvs=.001, **kwargs)

    def test_oscillating_misfits_stop_only_after_two_full_windows(self):
        result, _ = self.run_misfits([999., 1., 3., 3., 1., 8.], misfit_window=2, misfit_tol=.01)
        self.assertTrue(result.converged)
        self.assertEqual(result.n_iterations, 4)
        self.assertEqual(result.n_evaluations, 5)
        self.assertIn('Windowed misfit', result.message)
        np.testing.assert_array_equal(result.misfit_history, [999., 1., 3., 3., 1.])

    def test_equal_consecutive_misfits_do_not_stop_before_window_is_full(self):
        result, _ = self.run_misfits([999., 1., 1., 1.], misfit_window=2, misfit_tol=.01)
        self.assertFalse(result.converged)
        self.assertEqual(result.n_iterations, 3)

    def test_strict_tolerance_and_sliding_windows(self):
        # N=1: 2 -> 1 gives exactly 0.5, so it must not pass tol=0.5.
        result, _ = self.run_misfits([999., 2., 1., 1.1, 8.], misfit_window=1, misfit_tol=.5)
        self.assertEqual(result.n_iterations, 3)
        self.assertIn('Windowed misfit', result.message)

    def test_zero_tolerance_disables_window_stopping(self):
        result, _ = self.run_misfits([999., 1., 1., 1., 1.], misfit_window=1, misfit_tol=0.)
        self.assertFalse(result.converged)
        self.assertEqual(result.n_iterations, 4)

    def test_lbfgs_uses_the_same_accepted_iteration_windows(self):
        from test_vsapp_inversion import FORWARD, native_curve

        observed = native_curve(make_model((1.2, 2.5, 3.5)), 'water')
        result = invert_vsapp(make_model(), PERIODS, observed, rayp=.06, f0=3.5,
                              optimizer='lbfgs', method='water', kernel_kwargs=FORWARD,
                              misfit_window=2, misfit_tol=1., gradient_tol=0.)
        self.assertTrue(result.converged)
        self.assertEqual(result.n_iterations, 4)
        self.assertIn('Windowed misfit', result.message)
        self.assertLess(windowed_misfit_change(result.history, 2), 1.)

    def test_bad_window_rejected_before_forward(self):
        for window in [0, -1, 1.5, True]:
            with self.subTest(window=window), patch.object(VsappForward, 'evaluate') as forward:
                with self.assertRaises((ValueError, TypeError)):
                    invert_vsapp(make_model(), PERIODS, np.ones(PERIODS.size),
                                 rayp=.06, f0=3.5, misfit_window=window)
                forward.assert_not_called()
