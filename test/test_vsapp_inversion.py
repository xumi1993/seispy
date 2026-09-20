"""Independent derivative, optimization and station-average inversion checks."""

import unittest
from dataclasses import replace
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from seispy import DepModel, SynSeis
from seispy.inversion import (
    brocher_properties,
    invert_station_vsapp,
    invert_vsapp,
    smooth_gradient,
    vsapp_gradient,
)
from seispy.inversion.optimize import gradient_norm
from seispy.utils import vs2vprho
from seispy.vsapp import StationVsAppResult

PERIODS = np.array([.1, .3, .5, .8, 1.2, 2., 3.])
FORWARD = dict(dt=.05, npts=512, shift=5.)


def make_model(vs=(1.3, 2.35, 3.5)):
    speed = np.array(vs)
    vp, rho = vs2vprho(speed)
    return DepModel.read_layer_model(np.array([0., 1., 3.]), [1., 2., 0.], vp, speed, rho)


def native_curve(model, method='iter'):
    return SynSeis(model, .06, .05, npts=512).compute_vsapp(
        PERIODS, method=method, f0=3.5, shift=5.,
    )[0].vs_km_s


class TestInversionGradients(unittest.TestCase):
    def test_brocher_slopes_match_independent_central_differences(self):
        vs = np.array([.35, 1.2, 2.5, 3.5, 4.5])
        vp, rho, dvp, drho = brocher_properties(vs)
        expected = vs2vprho(vs)
        np.testing.assert_array_equal(vp, expected[0])
        np.testing.assert_array_equal(rho, expected[1])
        h = 1e-5
        plus, minus = vs2vprho(vs + h), vs2vprho(vs - h)
        np.testing.assert_allclose(dvp, (plus[0] - minus[0]) / (2*h), rtol=1e-8)
        np.testing.assert_allclose(drho, (plus[1] - minus[1]) / (2*h), rtol=1e-7)

    def test_coupled_weighted_gradient_matches_native_objective_difference(self):
        model = make_model()
        observed = native_curve(make_model((1.2, 2.5, 3.5)), method='water')
        sigma = np.linspace(.03, .1, PERIODS.size)
        _, _, dvp, drho = brocher_properties(model.vs)
        kernel = model.vsapp_kernel(.06, PERIODS, 3.5, method='water',
                                    vp_vs_derivative=dvp, rho_vs_derivative=drho, **FORWARD)
        misfit, gradient = vsapp_gradient(kernel, observed, vs_km_s=model.vs, sigma=sigma)
        np.testing.assert_allclose(misfit, .5 * np.mean(((native_curve(model, 'water')
                                                       - observed) / sigma)**2))
        result = invert_vsapp(model, PERIODS, observed, rayp=.06, f0=3.5,
                              method='water', sigma=sigma, zero_halfspace=False,
                              max_iterations=0, kernel_kwargs=FORWARD)
        np.testing.assert_allclose(result.gradient, gradient)
        for j in range(model.vs.size):
            h = 1e-5
            plus, minus = model.vs.copy(), model.vs.copy()
            plus[j] *= np.exp(h)
            minus[j] *= np.exp(-h)
            curve_plus = native_curve(make_model(plus), 'water')
            curve_minus = native_curve(make_model(minus), 'water')
            cost_plus = .5 * np.mean(((curve_plus - observed)/sigma)**2)
            cost_minus = .5 * np.mean(((curve_minus - observed)/sigma)**2)
            self.assertAlmostEqual(gradient[j], (cost_plus - cost_minus)/(2*h), delta=2e-6)

    def test_gaussian_uses_kilometres_and_excludes_fixed_halfspace(self):
        gradient = np.array([1., 0., 0., 1e9])
        depths = np.array([0., 1., 3., 4.])
        fixed = np.array([False, False, False, True])
        result = smooth_gradient(gradient, depths, 1., fixed_mask=fixed)
        self.assertEqual(result[-1], 0.)
        self.assertGreater(result[1], result[2])
        np.testing.assert_array_equal(gradient, [1., 0., 0., 1e9])
        np.testing.assert_array_equal(smooth_gradient(gradient, depths, 0., fixed_mask=fixed),
                                      [1., 0., 0., 0.])
        np.testing.assert_allclose(smooth_gradient(np.ones(4), depths, .8), 1.)


class TestLogParameterSearch(unittest.TestCase):
    def setUp(self):
        self.problem = SimpleNamespace(
            initial_vs=np.array([1., 3., 4.]), fixed=np.array([False, False, True]),
        )
        self.state = SimpleNamespace(model=SimpleNamespace(vs=self.problem.initial_vs.copy()),
                                     gradient=np.array([1., 1., 0.]), misfit=1.)

    def test_gradient_norm_uses_log_gradient(self):
        self.state.gradient = np.array([-3., 2., 0.])
        self.assertEqual(gradient_norm(self.state), 3.)



class TestVsappInversion(unittest.TestCase):
    def test_synthetic_inversion_decreases_misfit_and_updates_halfspace(self):
        for method in ('iter', 'water'):
            with self.subTest(method=method):
                initial = make_model()
                snapshots = {key: value.copy() for key, value in vars(initial).items()
                             if isinstance(value, np.ndarray)}
                observed = native_curve(make_model((1.2, 2.5, 3.5)), method)
                result = invert_vsapp(initial, PERIODS, observed, rayp=.06, f0=3.5,
                                       method=method, smooth_sigma_km=.5, max_iterations=80,
                                       kernel_kwargs=FORWARD)
                self.assertGreater(result.n_iterations, 0, result.message)
                self.assertLess(result.misfit_history[-1], result.misfit_history[0] * .05)
                self.assertEqual(result.n_evaluations, result.n_iterations + 1)
                self.assertNotEqual(result.history[1].vs_km_s[-1], initial.vs[-1])
                self.assertNotEqual(result.model.vs[-1], initial.vs[-1])
                for record in result.history:
                    self.assertLessEqual(record.max_update_lnvs, .05 + 1e-14)
                np.testing.assert_array_equal(result.kernel.jacobian[:, -1], 0.)
                self.assertEqual(result.gradient[-1], 0.)
                self.assertNotEqual(result.smoothed_gradient[-1], 0.)
                expected_vp, expected_rho = vs2vprho(result.model.vs)
                np.testing.assert_array_equal(result.model.vp, expected_vp)
                np.testing.assert_array_equal(result.model.rho, expected_rho)
                np.testing.assert_array_equal(result.model.model_array[:, 2], result.model.vs)
                np.testing.assert_allclose(
                    result.predicted_vsapp, native_curve(result.model, method),
                    rtol=1e-8, atol=1e-9,
                )
                for name, before in snapshots.items():
                    np.testing.assert_array_equal(getattr(initial, name), before)

    def test_every_updated_model_recomputes_coupled_kernel(self):
        calls = []
        original = DepModel.vsapp_kernel

        def capture(model, *args, **kwargs):
            vp, rho, dvp, drho = brocher_properties(model.vs)
            np.testing.assert_array_equal(model.vp, vp)
            np.testing.assert_array_equal(model.rho, rho)
            np.testing.assert_array_equal(kwargs['vp_vs_derivative'], dvp)
            np.testing.assert_array_equal(kwargs['rho_vs_derivative'], drho)
            calls.append(model.vs.copy())
            return original(model, *args, **kwargs)

        observed = native_curve(make_model((1.2, 2.5, 3.5)), 'water')
        with patch.object(DepModel, 'vsapp_kernel', capture):
            result = invert_vsapp(make_model(), PERIODS, observed, rayp=.06, f0=3.5,
                                   method='water', max_iterations=3, kernel_kwargs=FORWARD)
        self.assertEqual(result.n_iterations, 3)
        self.assertEqual(result.n_evaluations, len(calls))
        for record in result.history:
            self.assertTrue(any(np.array_equal(record.vs_km_s, speed) for speed in calls))
        self.assertGreater(len({tuple(speed) for speed in calls}), 1)

    def test_iteration_limit_has_explicit_status(self):
        observed = native_curve(make_model((1.2, 2.5, 3.5)), 'water')
        result = invert_vsapp(make_model(), PERIODS, observed, rayp=.06, f0=3.5,
                              method='water', max_iterations=1, kernel_kwargs=FORWARD)
        self.assertFalse(result.converged)
        self.assertIn('Maximum iterations', result.message)

    def test_initial_model_outside_former_bounds_is_not_clipped(self):
        model = make_model((1.3, 2.35, 4.6))
        observed = native_curve(model, 'water')
        result = invert_vsapp(model, PERIODS, observed, rayp=.06, f0=3.5,
                              method='water', max_iterations=0, kernel_kwargs=FORWARD)
        np.testing.assert_array_equal(result.model.vs, model.vs)

    def test_initial_exact_fit_needs_no_update(self):
        initial = make_model()
        observed = native_curve(initial)
        result = invert_vsapp(initial, PERIODS, observed, rayp=.06, f0=3.5,
                               kernel_kwargs=FORWARD)
        self.assertTrue(result.converged)
        self.assertEqual(result.n_iterations, 0)
        self.assertEqual(result.n_evaluations, 1)

    def test_initial_vp_and_density_are_rederived_without_mutating_input(self):
        initial = make_model()
        observed = native_curve(initial)
        initial.vp[:] = 8.
        initial.rho[:] = 3.8
        result = invert_vsapp(initial, PERIODS, observed, rayp=.06, f0=3.5,
                               max_iterations=0, kernel_kwargs=FORWARD)
        vp, rho = vs2vprho(initial.vs)
        np.testing.assert_array_equal(result.initial_model.vp, vp)
        np.testing.assert_array_equal(result.initial_model.rho, rho)
        np.testing.assert_array_equal(result.model.vp, vp)
        np.testing.assert_array_equal(initial.vp, 8.)
        np.testing.assert_array_equal(initial.rho, 3.8)

    def test_rejected_trials_do_not_overwrite_last_accepted_model(self):
        initial = make_model()
        observed = native_curve(make_model((1.2, 2.5, 3.5)))
        original = DepModel.vsapp_kernel
        counter = 0

        def fail_trials(model, *args, **kwargs):
            nonlocal counter
            counter += 1
            if counter > 1:
                raise ValueError('Invalid trial')
            return original(model, *args, **kwargs)

        with patch.object(DepModel, 'vsapp_kernel', fail_trials):
            result = invert_vsapp(initial, PERIODS, observed, rayp=.06, f0=3.5,
                                   max_backtracks=2, kernel_kwargs=FORWARD)
        self.assertFalse(result.converged)
        self.assertIn('SD forward calculation failed', result.message)
        self.assertEqual(result.n_evaluations, 2)
        np.testing.assert_array_equal(result.model.vs, initial.vs)
        self.assertEqual(result.n_iterations, 0)

    def test_single_halfspace_needs_nonzero_kernel_to_update(self):
        # A one-parameter half-space is also supported by the existing kernel.
        # With no neighbours, smoothing cannot propagate a nonzero update.
        model = make_model()
        model.vs = np.array([3.3])
        model.thickness = np.array([0.])
        model.depths_elev = np.array([0.])
        model.depths = np.array([0.])
        observed = np.full(PERIODS.size, 3.5)
        masked = invert_vsapp(model, PERIODS, observed, rayp=.06, f0=3.5,
                               kernel_kwargs=FORWARD)
        self.assertEqual(masked.model.vs[0], 3.3)
        free = invert_vsapp(model, PERIODS, observed, rayp=.06, f0=3.5,
                             zero_halfspace=False, kernel_kwargs=FORWARD)
        np.testing.assert_allclose(free.model.vs, [3.5], atol=1e-5)
        self.assertTrue(free.converged)

    def test_bad_options_fail_before_forward_evaluation(self):
        observed = np.ones(PERIODS.size) * 2.
        cases = [dict(sigma=0.), dict(sigma=[.1]), dict(smooth_sigma_km=-1.),
                 dict(max_step_km_s=.05), dict(max_iterations=1.5), dict(max_step_lnvs=0.),
                 dict(zero_halfspace='yes'),
                 dict(kernel_kwargs={'vp_vs_derivative': [0.]}),
                 dict(vs_bounds=(2., 4.)), dict(method='bad')]
        for options in cases:
            with self.subTest(options=options), patch.object(DepModel, 'vsapp_kernel') as forward:
                with self.assertRaises((TypeError, ValueError)):
                    invert_vsapp(make_model(), PERIODS, observed, rayp=.06, f0=3.5, **options)
                forward.assert_not_called()


class TestStationInversion(unittest.TestCase):
    def setUp(self):
        mean = native_curve(make_model())
        values = np.stack([mean-.01, mean+.01, mean+1.])
        statuses = np.full(values.shape, 'ok', dtype='<U24')
        statuses[2, 0] = 'nonpositive_radial'
        self.measurements = StationVsAppResult(
            curves=(), event=np.array(['A', 'B', 'C']), rayp=np.array([.05, .07, .09]),
            periods_s=PERIODS, vs_km_s=values, vs0_km_s=np.ones(3), status=statuses,
            vs0_status=('ok',)*3, reference='seispy-iter',
        )

    def test_station_uses_same_retained_events_for_all_periods_and_rayp(self):
        before = self.measurements.vs_km_s.copy()
        result = invert_station_vsapp(make_model(), self.measurements, f0=3.5,
                                       max_iterations=0, kernel_kwargs=FORWARD)
        np.testing.assert_array_equal(result.event_mask, [True, True, False])
        np.testing.assert_array_equal(result.observed_vsapp, before[:2].mean(axis=0))
        self.assertAlmostEqual(result.rayp, .06)
        np.testing.assert_array_equal(self.measurements.vs_km_s, before)

    def test_no_complete_event_is_rejected(self):
        bad = replace(self.measurements, status=np.full((3, PERIODS.size), 'invalid'))
        with self.assertRaisesRegex(ValueError, 'No event is valid'):
            invert_station_vsapp(make_model(), bad, f0=3.5)


if __name__ == '__main__':
    unittest.main()
