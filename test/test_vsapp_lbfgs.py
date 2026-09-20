"""Independent L-BFGS recursion, strong Wolfe and coupled RF inversion checks."""

import logging
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import test_vsapp_inversion as fixtures
from seispy.inversion import invert_station_vsapp, invert_vsapp
from seispy.inversion.forward import VsappForward
from seispy.inversion.invpara import InvPara
from seispy.inversion.lbfgs import LBFGSHistory
from seispy.inversion.linesearch import strong_wolfe
from test_vsapp_inversion import FORWARD, PERIODS, make_model, native_curve


def state_at(x, gradient, misfit=0.):
    return SimpleNamespace(model=SimpleNamespace(vs=np.exp(x)),
                           gradient=np.array(gradient, dtype=float), misfit=misfit)


class TestStrongWolfe(unittest.TestCase):
    def setUp(self):
        self.problem = SimpleNamespace(fixed=np.array([False]), initial_vs=np.ones(1))

    def search_quadratic(self, target, direction, *, invalid_above=None, **options):
        calls = []

        def evaluate(speed):
            x = np.log(speed)
            calls.append(x.copy())
            if invalid_above is not None and x[0] > invalid_above:
                raise ValueError('Outside valid forward domain')
            return state_at(x, x - target, .5 * float(np.sum((x - target)**2)))

        start = evaluate(np.ones(1))
        calls.clear()
        controls = InvPara(**options)
        trial = strong_wolfe(self.problem, controls, SimpleNamespace(evaluate=evaluate),
                             start, np.array([direction]))
        if trial is not None:
            delta = np.log(trial.model.vs)
            alpha = float(delta[0] / direction)
            slope0 = start.gradient[0] * direction
            slope1 = trial.gradient[0] * direction
            self.assertLessEqual(trial.misfit, start.misfit + controls.wolfe_c1 * alpha * slope0)
            if abs(slope1) > controls.wolfe_c2 * abs(slope0):
                self.assertAlmostEqual(abs(delta[0]), controls.max_step_lnvs)
                self.assertLess(slope1, 0.)
            self.assertLessEqual(abs(delta[0]), controls.max_step_lnvs + 1e-14)
        return trial, calls

    def test_expands_step_when_curvature_requires_it(self):
        trial, calls = self.search_quadratic(.1, .01, max_step_lnvs=.5, wolfe_c2=.1)
        self.assertIsNotNone(trial)
        self.assertGreater(len(calls), 2)
        self.assertGreater(calls[1][0], calls[0][0])

    def test_zoom_recovers_from_overshoot(self):
        trial, calls = self.search_quadratic(.1, 1., max_step_lnvs=1., wolfe_c2=.1)
        self.assertIsNotNone(trial)
        self.assertLess(calls[1][0], calls[0][0])

    def test_invalid_trial_is_bracketed_and_recovered(self):
        trial, calls = self.search_quadratic(.1, 1., max_step_lnvs=1., invalid_above=.3)
        self.assertIsNotNone(trial)
        self.assertGreater(calls[0][0], .3)

    def test_cap_accepts_decrease_even_without_strong_curvature(self):
        trial, calls = self.search_quadratic(1., 1., max_step_lnvs=.05)
        self.assertIsNotNone(trial)
        self.assertEqual(len(calls), 1)
        self.assertAlmostEqual(np.log(trial.model.vs[0]), .05)
        self.assertGreater(abs(trial.gradient[0]), .9)

    def test_before_cap_budget_does_not_accept_armijo_only(self):
        trial, calls = self.search_quadratic(1., .01, max_step_lnvs=.05, max_backtracks=0)
        self.assertIsNone(trial)
        self.assertEqual(len(calls), 1)

    def test_cap_with_positive_slope_still_zooms(self):
        trial, calls = self.search_quadratic(.1, 1., max_step_lnvs=.15, wolfe_c2=.1)
        self.assertIsNotNone(trial)
        self.assertGreater(len(calls), 1)
        self.assertLess(np.log(trial.model.vs[0]), .15)

    def test_budget_exhaustion_returns_no_trial(self):
        trial, calls = self.search_quadratic(.1, 1., max_step_lnvs=1., max_backtracks=0)
        self.assertIsNone(trial)
        self.assertEqual(len(calls), 1)

    def test_ascent_direction_never_evaluates_a_trial(self):
        with patch.object(VsappForward, 'evaluate') as evaluator:
            trial = strong_wolfe(self.problem, InvPara(),
                                 SimpleNamespace(evaluate=evaluator),
                                 state_at(np.zeros(1), [1.]), np.ones(1))
        self.assertIsNone(trial)
        evaluator.assert_not_called()


class TestLBFGSHistory(unittest.TestCase):
    def test_two_loop_matches_independent_dense_inverse_bfgs(self):
        history = LBFGSHistory(2)
        hessian = np.diag([2., 5., 9.])
        x = np.zeros(3)
        for step in [np.array([.1, .2, .3]), np.array([.2, -.1, .1]),
                     np.array([-.15, .2, .05])]:
            gradient = np.einsum('ij,j->i', hessian, x)
            new_gradient = np.einsum('ij,j->i', hessian, x + step)
            self.assertTrue(history.update(state_at(x, gradient),
                                           state_at(x + step, new_gradient), np.zeros(3, bool)))
            x += step
        self.assertEqual(len(history.pairs), 2)
        last_s, last_y, _ = history.pairs[-1]
        dense = np.eye(3) * np.dot(last_s, last_y) / np.dot(last_y, last_y)
        for step, change, reciprocal in history.pairs:
            transform = np.eye(3) - reciprocal * np.outer(step, change)
            dense = np.einsum('ij,jk,lk->il', transform, dense, transform)
            dense += reciprocal * np.outer(step, step)
        gradient = np.array([1., -2., 3.])
        np.testing.assert_allclose(history.direction(gradient),
                                   -np.einsum('ij,j->i', dense, gradient), atol=1e-13)

    def test_bad_curvature_skipped_and_fixed_coordinate_excluded(self):
        history = LBFGSHistory(2)
        start = state_at(np.zeros(2), [0., 0.])
        self.assertFalse(history.update(start, state_at([.1, 0.], [-1., 0.]),
                                        np.zeros(2, bool)))
        self.assertFalse(history.update(start, start, np.zeros(2, bool)))
        self.assertTrue(history.update(start, state_at([.1, .3], [1., 1e9]),
                                       np.array([False, True])))
        self.assertEqual(history.pairs[0][0][-1], 0.)
        self.assertEqual(history.pairs[0][1][-1], 0.)
        history.clear()
        self.assertEqual(len(history.pairs), 0)


class TestLBFGSInversion(unittest.TestCase):
    def test_real_inversion_recomputes_kernels_and_satisfies_wolfe(self):
        original = VsappForward.evaluate
        for method in ('water', 'iter'):
            with self.subTest(method=method):
                evaluated = []

                def capture(forward, speed, evaluated=evaluated):
                    state = original(forward, speed)
                    evaluated.append(state)
                    return state

                initial = make_model()
                observed = native_curve(make_model((1.2, 2.5, 3.5)), method)
                with patch.object(VsappForward, 'evaluate', capture):
                    result = invert_vsapp(initial, PERIODS, observed, rayp=.06, f0=3.5,
                                          optimizer='lbfgs', method=method, max_iterations=30,
                                          kernel_kwargs=FORWARD)
                self.assertEqual(result.optimizer, 'lbfgs')
                self.assertEqual(result.n_evaluations, len(evaluated))
                self.assertGreater(result.n_iterations, 0, result.message)
                self.assertLess(result.misfit_history[-1], result.misfit_history[0] * .1)
                self.assertNotEqual(result.history[1].vs_km_s[-1], initial.vs[-1])
                self.assertNotEqual(result.model.vs[-1], initial.vs[-1])
                for state in evaluated:
                    np.testing.assert_array_equal(state.kernel.jacobian[:, -1], 0.)
                    self.assertEqual(state.gradient[-1], 0.)
                states = []
                for record in result.history:
                    states.append(next(s for s in evaluated
                                       if np.array_equal(s.model.vs, record.vs_km_s)))
                    self.assertLessEqual(record.max_update_lnvs, .05 + 1e-14)
                for previous, current in zip(states[:-1], states[1:], strict=True):
                    step = np.log(current.model.vs) - np.log(previous.model.vs)
                    slope0 = float(np.sum(previous.gradient * step))
                    slope1 = float(np.sum(current.gradient * step))
                    self.assertLess(slope0, 0.)
                    self.assertLessEqual(current.misfit, previous.misfit + 1e-4 * slope0 + 1e-14)
                    if abs(slope1) > .9 * abs(slope0) + 1e-12:
                        self.assertAlmostEqual(np.max(np.abs(step)), .05)
                        self.assertLess(slope1, 0.)
                np.testing.assert_allclose(
                    result.predicted_vsapp, native_curve(result.model, method),
                    rtol=1e-8, atol=1e-9,
                )
                np.testing.assert_array_equal(initial.vs, [1.3, 2.35, 3.5])

    def test_far_initial_model_uses_cap_and_reduces_iterative_vsapp_misfit(self):
        initial = make_model((2.8, 3.1, 3.5))
        observed = native_curve(make_model((1.2, 2.5, 3.5)), 'iter')
        with self.assertLogs('test.lbfgs.cap', level='INFO') as captured:
            # Allow the half-space to evolve along with the finite layers.
            result = invert_vsapp(initial, PERIODS, observed, rayp=.06, f0=3.5,
                                  optimizer='lbfgs', method='iter', max_iterations=30,
                                  kernel_kwargs=FORWARD, log=logging.getLogger('test.lbfgs.cap'))
        self.assertGreater(result.n_iterations, 0, result.message)
        self.assertAlmostEqual(result.history[1].max_update_lnvs, .05)
        self.assertIn('Accepted capped L-BFGS step', '\n'.join(captured.output))
        self.assertLess(result.misfit_history[-1], .1 * result.misfit_history[0])
        self.assertNotEqual(result.history[1].vs_km_s[-1], initial.vs[-1])
        self.assertNotEqual(result.model.vs[-1], initial.vs[-1])
        np.testing.assert_array_equal(result.kernel.jacobian[:, -1], 0.)
        for record in result.history:
            self.assertLessEqual(record.max_update_lnvs, .05 + 1e-14)
        np.testing.assert_array_equal(initial.vs, [2.8, 3.1, 3.5])

    def test_station_entry_point_accepts_lbfgs(self):
        fixture = fixtures.TestStationInversion()
        fixture.setUp()
        result = invert_station_vsapp(make_model(), fixture.measurements, f0=3.5,
                                       optimizer='lbfgs', max_iterations=0,
                                       kernel_kwargs=FORWARD)
        self.assertEqual(result.optimizer, 'lbfgs')
        np.testing.assert_array_equal(result.event_mask, [True, True, False])

    def test_bad_options_rejected_before_forward(self):
        for options in [dict(optimizer='bad'), dict(lbfgs_memory=0), dict(lbfgs_memory=1.5),
                        dict(wolfe_c1=0), dict(wolfe_c2=1), dict(wolfe_c1=.9, wolfe_c2=.1)]:
            with self.subTest(options=options), patch.object(VsappForward, 'evaluate') as forward:
                with self.assertRaises((ValueError, TypeError)):
                    invert_vsapp(make_model(), PERIODS, np.ones(PERIODS.size),
                                 rayp=.06, f0=3.5, **options)
                forward.assert_not_called()

    def test_wolfe_failure_keeps_last_accepted_model(self):
        observed = native_curve(make_model((1.2, 2.5, 3.5)), 'water')
        original = VsappForward.evaluate

        def reject_trials(forward, speed):
            if not np.array_equal(speed, forward.problem.initial_vs):
                raise ValueError('Invalid trial model')
            return original(forward, speed)

        with patch.object(VsappForward, 'evaluate', reject_trials):
            result = invert_vsapp(make_model(), PERIODS, observed, rayp=.06, f0=3.5,
                                  optimizer='lbfgs', method='water', kernel_kwargs=FORWARD)
        self.assertFalse(result.converged)
        self.assertIn('L-BFGS line search failed', result.message)
        self.assertEqual(result.n_iterations, 0)
        np.testing.assert_array_equal(result.model.vs, [1.3, 2.35, 3.5])
