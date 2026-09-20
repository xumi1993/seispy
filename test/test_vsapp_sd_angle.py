"""Gradient-angle step control for steepest descent in ln(Vs)."""

import unittest

import numpy as np
from seispy.inversion.optimize import gradient_angle


class TestGradientAngle(unittest.TestCase):
    def test_angles_and_magnitude_invariance(self):
        first = np.array([1., 0.])
        for angle in [0., 90., 119., 120., 121., 180.]:
            radians = np.radians(angle)
            second = np.array([np.cos(radians), np.sin(radians)])
            with self.subTest(angle=angle):
                self.assertAlmostEqual(gradient_angle(first, second), angle)
                self.assertAlmostEqual(gradient_angle(first * 1e-250, second * 1e250), angle)

    def test_missing_or_zero_gradient_has_no_angle(self):
        gradient = np.array([1., 2.])
        self.assertIsNone(gradient_angle(None, gradient))
        self.assertIsNone(gradient_angle(np.zeros(2), gradient))
        self.assertIsNone(gradient_angle(gradient, np.zeros(2)))

    def test_inputs_are_not_modified(self):
        previous, current = np.array([1., 2.]), np.array([2., -3.])
        gradient_angle(previous, current)
        np.testing.assert_array_equal(previous, [1., 2.])
        np.testing.assert_array_equal(current, [2., -3.])


class TestSDStepControl(unittest.TestCase):
    def run_sd(self, gradients, misfits=None, initial=None, **kwargs):
        from dataclasses import replace
        from unittest.mock import patch

        from seispy.inversion import invert_vsapp
        from seispy.inversion.forward import VsappForward
        from test_vsapp_inversion import FORWARD, PERIODS, make_model, native_curve

        if initial is None:
            initial = make_model()
        observed = native_curve(initial, 'water')
        original = VsappForward.evaluate
        evaluations = []

        def evaluate(forward, speed):
            state = original(forward, speed)
            index = len(evaluations)
            # Prescribed gradients isolate step control; each trial still runs a real kernel.
            gradient = np.array(gradients[index], dtype=float)
            misfit = 1. + index if misfits is None else misfits[index]
            state = replace(state, gradient=gradient, misfit=misfit)
            evaluations.append(state)
            return state

        options = dict(max_iterations=len(gradients)-1, max_step_lnvs=.1,
                       smooth_sigma_km=0., misfit_tol=.1, max_backtracks=0)
        options.update(kwargs)
        with patch.object(VsappForward, 'evaluate', evaluate):
            result = invert_vsapp(initial, PERIODS, observed, rayp=.06, f0=3.5,
                                  method='water', optimizer='gd', kernel_kwargs=FORWARD, **options)
        return result, evaluations

    @staticmethod
    def gradients(angles):
        return [[np.cos(np.radians(a)), np.sin(np.radians(a)), 0.] for a in angles]

    def test_persistent_halving_once_per_large_angle_and_no_reset(self):
        result, evaluations = self.run_sd(self.gradients([0., 121., 242., 242., 242.]))
        np.testing.assert_allclose([h.step_length_lnvs for h in result.history[1:]],
                                   [.1, .05, .025, .025])
        self.assertIsNone(result.history[1].gradient_angle_deg)
        np.testing.assert_allclose([h.gradient_angle_deg for h in result.history[2:]],
                                   [121., 121., 0.], atol=1e-12)
        self.assertEqual(result.n_evaluations, 5)
        self.assertEqual(len(evaluations), 5)
        for previous, current, record in zip(evaluations[:-1], evaluations[1:],
                                             result.history[1:], strict=True):
            expected = -record.step_length_lnvs * previous.gradient / max(abs(previous.gradient))
            np.testing.assert_allclose(np.log(current.model.vs / previous.model.vs),
                                       expected, atol=1e-14)
            self.assertEqual(current.model.vs[-1], evaluations[0].model.vs[-1])

    def test_threshold_is_strictly_greater_than_120_degrees(self):
        for angle in [90., 119., 120., 120.01, 180.]:
            with self.subTest(angle=angle):
                result, _ = self.run_sd(self.gradients([0., angle, angle]))
                expected = .05 if angle > 120 else .1
                self.assertEqual(result.history[2].step_length_lnvs, expected)

    def test_increasing_misfit_is_accepted_without_false_convergence(self):
        result, _ = self.run_sd(self.gradients([0., 0., 0., 0.]))
        np.testing.assert_array_equal(result.misfit_history, [1., 2., 3., 4.])
        self.assertEqual(result.n_iterations, 3)
        self.assertFalse(result.converged)
        self.assertIn('Maximum iterations', result.message)
        self.assertEqual(result.n_evaluations, 4)

    def test_same_fractional_step_for_different_layer_velocities(self):
        result, evaluations = self.run_sd([[1., 1., 0.], [1., 1., 0.]])
        ratio = result.model.vs / evaluations[0].model.vs
        np.testing.assert_allclose(ratio, [np.exp(-.1), np.exp(-.1), 1.])

    def test_no_velocity_bounds_are_reintroduced(self):
        from test_vsapp_inversion import make_model

        result, _ = self.run_sd([[1., -1., 0.], [1., -1., 0.]],
                                initial=make_model((.36, 4.4, 4.)))
        self.assertLess(result.model.vs[0], .35)
        self.assertGreater(result.model.vs[1], 4.5)
